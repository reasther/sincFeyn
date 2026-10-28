/* checkpoint.c -- point lists, and save/restore of the refinement state so that a killed run
 * resumes where it stopped.  File: SINCQED_CKPT_DIR/sincqed_<graph>.ckpt (default: current
 * directory), written every SINCQED_CKPT seconds (default 600, 0 = off) between refinement blocks
 * and at the end of every level, replaced atomically, removed when the graph finishes.  A checkpoint
 * is used only if h, N, eps and top all match the new run.  It also keeps the predictor's table of 2h
 * values and the eps shadow sums, so a resumed run predicts and reports the shadow as an unbroken one would. */
#include "sq.h"

double CKPT_EVERY = 600;

void push(vec_t *v, const int *k){
    if(v->n == v->cap){ v->cap = v->cap ? 2*v->cap : 1024; v->p = (pt_t*)realloc(v->p, v->cap*sizeof(pt_t)); }
    pt_t *p = &v->p[v->n++]; memset(p, 0, sizeof *p); for(int i = 0; i < DIM && i < MAXD; i++) p->k[i] = (signed char)k[i];
}
/* not the PFLAG slot.  The order is that of memcmp on the little-endian ints the points used to be stored
 * as (before 2026-10-22), so the summation order, and so every sum, is unchanged. */
static int cmp_pt(const void *a, const void *b){
    const signed char *x = ((const pt_t*)a)->k, *y = ((const pt_t*)b)->k;
    for(int i = 0; i < DIM; i++) if(x[i] != y[i]){
        unsigned u = __builtin_bswap32((unsigned)(int)x[i]), w = __builtin_bswap32((unsigned)(int)y[i]);
        return u < w ? -1 : 1; }
    return 0;
}
/* qsort in parallel: P chunks sorted at once, then merged pairwise (same order as one qsort with cmp_pt) */
static void psort(pt_t *a, long n){
    int P = 1;
#ifdef _OPENMP
    while(P < omp_get_max_threads() && P < 64) P <<= 1;
#endif
    if(P == 1 || n < 65536){ qsort(a, n, sizeof(pt_t), cmp_pt); return; }
    long *o = (long*)malloc((P + 1)*sizeof(long)); for(int i = 0; i <= P; i++) o[i] = n*i/P;
    pt_t *t = (pt_t*)malloc(n*sizeof(pt_t)), *src = a, *dst = t;
    #pragma omp parallel for schedule(dynamic,1)
    for(int i = 0; i < P; i++) qsort(a + o[i], o[i+1] - o[i], sizeof(pt_t), cmp_pt);
    for(int w = 1; w < P; w <<= 1){
        #pragma omp parallel for schedule(dynamic,1)
        for(int i = 0; i < P; i += 2*w){
            long x = o[i], xe = o[i+w], y = xe, ye = o[i+2*w], z = o[i];
            while(x < xe && y < ye) dst[z++] = cmp_pt(&src[y], &src[x]) < 0 ? src[y++] : src[x++];
            while(x < xe) dst[z++] = src[x++];
            while(y < ye) dst[z++] = src[y++];
        }
        pt_t *u = src; src = dst; dst = u;
    }
    if(src != a) memcpy(a, src, n*sizeof(pt_t));
    free(t); free(o);
}
void uniq(vec_t *v){
    psort(v->p, v->n);
    long m = 0;
    for(long q = 0; q < v->n; q++){
        if(!m || cmp_pt(&v->p[q], &v->p[m-1])) v->p[m++] = v->p[q];
        else if(DIM < MAXD){ int a = v->p[m-1].k[PFLAG], b = v->p[q].k[PFLAG];   /* bit 0: both may skip; tier: the larger */
            v->p[m-1].k[PFLAG] = (a & b & 1) | ((a >> 1) > (b >> 1) ? a & ~1 : b & ~1); }
    }
    v->n = m;
}

typedef struct { char magic[8], name[32]; double h, eps; int N, top, dim, k; } ckhead_t;

static void head(ckhead_t *hd, const graph_t *g, double h, int N, double eps, int top){
    memset(hd, 0, sizeof *hd); memcpy(hd->magic, "SQCK4", 6); snprintf(hd->name, sizeof hd->name, "%s", g->name);
    hd->h = h; hd->eps = eps; hd->N = N; hd->top = top; hd->dim = DIM; hd->k = K;
}

void ckpt_name(char *buf, size_t n, const graph_t *g){
    const char *d = getenv("SINCQED_CKPT_DIR"); snprintf(buf, n, "%s/sincqed_%s.ckpt", d && *d ? d : ".", g->name);
}

int ckpt_save(const graph_t *g, double h, int N, double eps, int top, const ckpt_t *c){
    char fn[4096], tmp[4200]; ckpt_name(fn, sizeof fn, g); snprintf(tmp, sizeof tmp, "%s.tmp", fn);
    FILE *f = fopen(tmp, "wb"); if(!f) return -1;
    ckhead_t hd; head(&hd, g, h, N, eps, top);
    int ok = fwrite(&hd, sizeof hd, 1, f) == 1
          && fwrite(&c->lev, sizeof c->lev, 1, f) == 1 && fwrite(&c->next, sizeof c->next, 1, f) == 1
          && fwrite(&c->s, sizeof(double), 1, f) == 1 && fwrite(&c->a, sizeof(double), 1, f) == 1
          && fwrite(&c->drop, sizeof(double), 1, f) == 1 && fwrite(&c->u, sizeof(double), 1, f) == 1
          && fwrite(&c->ne, sizeof(long), 1, f) == 1 && fwrite(&c->q, sizeof(long), 1, f) == 1
          && fwrite(c->Slev, sizeof c->Slev, 1, f) == 1 && fwrite(&c->nl, sizeof c->nl, 1, f) == 1
          && fwrite(&c->R, sizeof c->R, 1, f) == 1
          && fwrite(&c->cur.n, sizeof(long), 1, f) == 1 && fwrite(c->cur.p, sizeof(pt_t), c->cur.n, f) == (size_t)c->cur.n
          && fwrite(&c->nxt.n, sizeof(long), 1, f) == 1 && fwrite(c->nxt.p, sizeof(pt_t), c->nxt.n, f) == (size_t)c->nxt.n
          && fwrite(c->es, sizeof c->es, 1, f) == 1 && fwrite(c->cnt, sizeof c->cnt, 1, f) == 1 && fwrite(&c->esok, sizeof c->esok, 1, f) == 1
          && fwrite(&c->vsz, sizeof c->vsz, 1, f) == 1 && fwrite(&c->vcap, sizeof c->vcap, 1, f) == 1 && fwrite(&c->vn, sizeof c->vn, 1, f) == 1
          && fwrite(c->vs, c->vsz, c->vcap, f) == (size_t)c->vcap;
    ok = (fclose(f) == 0) && ok;
    if(!ok || rename(tmp, fn)){ remove(tmp); say("# warning: could not write checkpoint %s\n", fn); return -1; }
    return 0;
}

static int readvec(FILE *f, vec_t *v){
    long n; if(fread(&n, sizeof n, 1, f) != 1 || n < 0) return 0;
    v->p = (pt_t*)malloc((n ? n : 1)*sizeof(pt_t)); v->cap = n ? n : 1; v->n = n;
    return fread(v->p, sizeof(pt_t), n, f) == (size_t)n;
}

int ckpt_load(const graph_t *g, double h, int N, double eps, int top, ckpt_t *c){
    if(CKPT_EVERY <= 0) return 0;
    char fn[4096]; ckpt_name(fn, sizeof fn, g);
    FILE *f = fopen(fn, "rb"); if(!f) return 0;
    ckhead_t want, got; head(&want, g, h, N, eps, top);
    int ok = fread(&got, sizeof got, 1, f) == 1 && !memcmp(&got, &want, sizeof got);
    if(!ok){ fclose(f); say("# checkpoint %s is for a different run (h, N, eps or top differ) or an older format: ignored\n", fn); return 0; }
    ok = fread(&c->lev, sizeof c->lev, 1, f) == 1 && fread(&c->next, sizeof c->next, 1, f) == 1
      && fread(&c->s, sizeof(double), 1, f) == 1 && fread(&c->a, sizeof(double), 1, f) == 1
      && fread(&c->drop, sizeof(double), 1, f) == 1 && fread(&c->u, sizeof(double), 1, f) == 1
      && fread(&c->ne, sizeof(long), 1, f) == 1 && fread(&c->q, sizeof(long), 1, f) == 1
      && fread(c->Slev, sizeof c->Slev, 1, f) == 1 && fread(&c->nl, sizeof c->nl, 1, f) == 1
      && fread(&c->R, sizeof c->R, 1, f) == 1 && readvec(f, &c->cur) && readvec(f, &c->nxt);
    int want_vsz = c->vsz; c->vs = NULL; c->vcap = c->vn = 0; c->esok = 0;
    if(ok){
        ok = fread(c->es, sizeof c->es, 1, f) == 1 && fread(c->cnt, sizeof c->cnt, 1, f) == 1 && fread(&c->esok, sizeof c->esok, 1, f) == 1
          && fread(&c->vsz, sizeof c->vsz, 1, f) == 1 && c->vsz == want_vsz && fread(&c->vcap, sizeof c->vcap, 1, f) == 1
          && fread(&c->vn, sizeof c->vn, 1, f) == 1 && c->vcap >= 0 && (c->vcap & (c->vcap - 1)) == 0;
        if(ok && c->vcap){ c->vs = malloc((size_t)c->vcap*c->vsz); ok = c->vs && fread(c->vs, c->vsz, c->vcap, f) == (size_t)c->vcap; }
    }
    fclose(f);
    if(!ok){ say("# checkpoint %s is damaged: ignored\n", fn); free(c->cur.p); free(c->nxt.p); free(c->vs); memset(c, 0, sizeof *c); c->vsz = want_vsz; return 0; }
    return 1;
}

void ckpt_remove(const graph_t *g){ char fn[4096]; ckpt_name(fn, sizeof fn, g); remove(fn); }
