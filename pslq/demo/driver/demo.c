/* demo.c -- "sincqed demo": the plain Sinc sum of the one- and two-loop graphs over all their free
 * Schwinger parameters, without any of the reductions the production sum uses.
 *
 *   sincqed demo h[,h,...] [eps|auto] [GRAPH...]       (default graphs: all1 all2)
 *
 * Every line has its own parameter s_i = e^{z_i}, z_i = k_i h, except the Cheng-Wu photon (s = 1):
 * 2 variables at one loop, 5 at two loops (the Pi2 graph splits its photon into its two segments).
 * There are no series classes, no rule for the external-vertex split, no coarse levels and no
 * predictor.  The only truncation is eps: the sum is a flood fill from z = 0 over the lattice, and the
 * neighbours of a node are visited when the node's integrand density in z, |f(s) prod s_i|, exceeds
 * eps (nodes below it are still summed: they are the boundary).  The integrand is evaluated after a
 * projective rescaling (term()), so large |z| cannot overflow, and the box |z_i| <= SINCQED_ZMAX
 * (default 300) is only a safety limit far outside the support; the log says if the fill reaches it.
 * Some directions decay slowly (e.g. a photon and an electron growing together, about e^-|z|/2 in
 * the largest z), so at eps = 1e-16 the fill reaches |z| ~ 70.
 *
 * eps = auto (the default) ties the truncation to the Sinc error of a single propagator,
 * eps(h) = SINCQED_EPSC * exp(-pi^2/h), floored at 1e-20.  The truncation error at two loops is about
 * 1e4-1e5 eps (measured at h = 1), and the Sinc error about 2e-3 exp(-pi^2/h), so the default
 * SINCQED_EPSC = 1e-10 keeps the truncation about 100 times below the discretization error.
 *
 * A point is evaluated in double precision and recomputed in QREAL when 16 eps_double sum|terms|
 * exceeds G2_ABSERR (1e-15), as in the production sum, or exceeds eps/10 in the density, so that the
 * fill is never steered by round-off.  Results go to stdout, the log and
 * SINCQED_DEMO_RESULTS (default sincqed_demo.txt), one line per graph and h; progress goes to stderr
 * (every SINCQED_PROGRESS seconds).
 */
#include "sq.h"
#include <float.h>
#include <stdint.h>

#define DMAX 6
static int DD_, KMAX, BITS, FIX;
static const graph_t *G;
static double NORM, HH, EPSV;
static uint64_t OFF;

static uint64_t pack(const int *k, int n){ uint64_t c = 0; for(int i = 0; i < n; i++) c = (c << BITS) | (uint64_t)(k[i] + (int64_t)OFF); return c; }
static void unpack(uint64_t c, int *k, int n){ for(int i = n-1; i >= 0; i--){ k[i] = (int)(c & ((1ull << BITS) - 1)) - (int)OFF; c >>= BITS; } }

/* the integrand at node k times prod(s_i h); *fb counts the high-precision recomputations.
 * The integrand times the measure is homogeneous of degree -n in all n line parameters (the Cheng-Wu
 * one included), so it is evaluated at s/lambda, lambda = max(1, max s_i), and multiplied by
 * lambda^-n prod s_i, all in logarithms: no parameter exceeds 1, and no |z| overflows. */
static double term(const int *k, double *dens, long *fb){
    double z[DMAX+1], ab, v, zl = 0, sz = 0;
    for(int i = 0; i < DD_; i++){ double zz = k[i]*HH; sz += zz; if(zz > zl) zl = zz; }
    for(int i = 0; i < DD_; i++) z[i] = exp(k[i]*HH - zl);
    const double cw = exp(-zl), fac = exp(sz - (DD_ + 1)*zl);
    if(G->ud){ z[DD_] = cw; v = G->ud(z, &ab); }
    else { double s[16]; int j = 0; for(int i = 0; i < G->n; i++) s[i] = i == FIX ? cw : z[j++];
           v = NORM*G->ptd(s, 0, &ab); ab *= fabs(NORM); }
    double hw = fac; for(int i = 0; i < DD_; i++) hw *= HH;
    /* recompute in QREAL when the rounding bound exceeds G2_ABSERR per term, or exceeds eps/10 in the
     * density: there the double value is round-off noise (about 1e-19 in slow directions such as an
     * electron and a photon shrinking together), and deciding the fill on it would follow the noise */
    int noisy = 0;
    if(!(ab*hw*16*DBL_EPSILON <= ABSERR0) || !(ab*fac*16*DBL_EPSILON <= 0.1*EPSV)){
        QREAL zq[DMAX+1], abq, vq;
        for(int i = 0; i <= DD_; i++) zq[i] = (QREAL)z[i];
        if(G->uq) vq = G->uq(zq, &abq);
        else { QREAL s[16]; int j = 0; for(int i = 0; i < G->n; i++) s[i] = i == FIX ? (QREAL)cw : zq[j++];
               vq = (QREAL)NORM*G->ptq(s, 0, &abq); }
        v = (double)vq; (*fb)++;
        noisy = !((double)abq*fac*16*QEPS <= 0.1*EPSV);      /* still noise: summed, not expanded */
    }
    if(!isfinite(v)) v = 0;
    *dens = noisy ? 0 : v*fac; return v*hw;
}

/* The flood fill is organised by columns: a column is a line of nodes along the last variable, keyed by
 * the other D-1 indices.  For each column the visited nodes are kept as a list of intervals, so memory
 * goes with the number of columns (about h^-(D-1)) rather than the number of nodes (h^-D), typically
 * 50-100 times less.  The visited set is the same as a node-by-node fill: every neighbour of an
 * expanded node is evaluated, along the column and in the 2(D-1) adjacent columns.  Each generation
 * takes the work items (column, run of rows) produced by the previous one, sorted, so the sum is
 * accumulated in a fixed order (compensated) whatever the number of threads. */
typedef struct { uint64_t key; int16_t lo, hi; } item_t;
typedef struct { int16_t *iv; int n, cap; } col_t;
static col_t *COLS; static long NCOL, CAPCOL;
static uint64_t *HK; static long *HV; static size_t HCAP, HN;      /* column key -> index */
static inline size_t hsh(uint64_t x){ x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33; return (size_t)x; }
static void hput(uint64_t k, long v){
    size_t i = hsh(k) & (HCAP - 1);
    while(HK[i]) i = (i + 1) & (HCAP - 1);
    HK[i] = k; HV[i] = v; HN++;
}
static long hget(uint64_t k){              /* index of column k, created if new (serial) */
    if(2*(HN + 1) > HCAP){
        size_t oc = HCAP; uint64_t *ok = HK; long *ov = HV;
        HCAP = HCAP ? 2*HCAP : (1u << 16); HN = 0;
        HK = (uint64_t*)calloc(HCAP, sizeof *HK); HV = (long*)malloc(HCAP*sizeof *HV);
        if(!HK || !HV){ fprintf(stderr, "demo: out of memory (%zu columns)\n", HCAP); exit(1); }
        for(size_t i = 0; i < oc; i++) if(ok[i]) hput(ok[i], ov[i]);
        free(ok); free(ov);
    }
    size_t i = hsh(k) & (HCAP - 1);
    while(HK[i]){ if(HK[i] == k) return HV[i]; i = (i + 1) & (HCAP - 1); }
    if(NCOL == CAPCOL){ CAPCOL = CAPCOL ? 2*CAPCOL : 4096; COLS = (col_t*)realloc(COLS, CAPCOL*sizeof *COLS); }
    COLS[NCOL].iv = NULL; COLS[NCOL].n = COLS[NCOL].cap = 0;
    HK[i] = k; HV[i] = NCOL; HN++;
    return NCOL++;
}
static int cmpi(const void *a, const void *b){
    const item_t *x = (const item_t*)a, *y = (const item_t*)b;
    if(x->key != y->key) return x->key < y->key ? -1 : 1;
    return x->lo != y->lo ? x->lo - y->lo : x->hi - y->hi;
}
/* compensated (Neumaier) accumulation */
static inline void nadd(double *s, double *c, double x){
    double t = *s + x; *c += fabs(*s) >= fabs(x) ? (*s - t) + x : (x - t) + *s; *s = t;
}

typedef struct { double S, A, edge, t; long nev, nfb, nexp, ncol; int zhit; double zreach; } dres_t;
typedef struct { double s, c, a, e; long nev, nfb, nexp; int zhit, maxk; } gres_t;

static dres_t demo_sum(const graph_t *g, double h, double eps, double zmax){
    dres_t R = {0}; double t0 = wtime();
    G = g; HH = h; FIX = g->fixed; EPSV = eps;
    DD_ = g->ud ? g->Du : g->n - 1;
    NORM = pow(-0.25, g->L)/16;
    int DC = DD_ - 1;                                   /* column key: the first DD_-1 indices */
    BITS = DC ? 64/DC : 16; if(BITS > 16) BITS = 16; OFF = 1ull << (BITS - 1);
    KMAX = (int)floor(zmax/h + 1e-9); if(KMAX > (int)OFF - 2) KMAX = (int)OFF - 2;
    NCOL = 0; HN = 0; if(HK) memset(HK, 0, HCAP*sizeof *HK);
    int k0[DMAX] = {0};
    item_t *W = (item_t*)malloc(sizeof *W); long nw = 1;
    W[0].key = DC ? pack(k0, DC) : 1; W[0].lo = W[0].hi = 0;
    double S = 0, Sc = 0, A = 0, edge = 0; long nfb = 0, nexp = 0, nev = 0; int zhit = 0, maxk = 0;
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
#endif
    int NR = 2*KMAX + 1;
    unsigned char *MK = (unsigned char*)calloc((size_t)nth*NR, 1);
    for(int gen = 0; nw; gen++){
        qsort(W, nw, sizeof *W, cmpi);
        long ng = 0; long *gs = (long*)malloc((nw + 1)*sizeof *gs), *gi;
        for(long i = 0; i < nw; i++) if(!i || W[i].key != W[i-1].key) gs[ng++] = i;
        gs[ng] = nw;
        gi = (long*)malloc((ng ? ng : 1)*sizeof *gi);
        for(long q = 0; q < ng; q++) gi[q] = hget(W[gs[q]].key);
        gres_t *GR = (gres_t*)calloc(ng ? ng : 1, sizeof *GR);
        item_t **ob = (item_t**)calloc(nth, sizeof *ob); long *on = (long*)calloc(nth, sizeof *on), *oc = (long*)calloc(nth, sizeof *oc);
        #pragma omp parallel
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            unsigned char *M = MK + (size_t)tid*NR;            /* 0 unvisited, 1 visited earlier, 2 new, 3 new and expanded */
            int *stk = (int*)malloc(NR*sizeof *stk);
            #pragma omp for schedule(dynamic, 16)
            for(long q = 0; q < ng; q++){
                col_t *C = &COLS[gi[q]]; gres_t *r = &GR[q];
                int k[DMAX];
                if(DC) unpack(W[gs[q]].key, k, DC);
                int tmin = NR, tmax = -1, ns = 0;
                for(int m = 0; m < C->n; m++) for(int x = C->iv[2*m]; x <= C->iv[2*m+1]; x++){
                    M[x + KMAX] = 1; if(x + KMAX < tmin) tmin = x + KMAX; if(x + KMAX > tmax) tmax = x + KMAX; }
                #define VISIT(x) do { int x_ = (x); k[DC] = x_; double d_, v_ = term(k, &d_, &r->nfb); \
                    nadd(&r->s, &r->c, v_); r->a += fabs(v_); r->nev++; \
                    int ex_ = fabs(d_) > eps; M[x_ + KMAX] = 2 + ex_; \
                    if(x_ + KMAX < tmin) tmin = x_ + KMAX; if(x_ + KMAX > tmax) tmax = x_ + KMAX; \
                    for(int j_ = 0; j_ <= DC; j_++){ int a_ = abs(k[j_]); if(a_ > r->maxk) r->maxk = a_; if(ex_ && a_ == KMAX) r->zhit = 1; } \
                    if(ex_){ r->nexp++; stk[ns++] = x_; } else r->e += v_; } while(0)
                for(long i = gs[q]; i < gs[q+1]; i++)
                    for(int x = W[i].lo; x <= W[i].hi; x++) if(!M[x + KMAX]) VISIT(x);
                while(ns){
                    int x = stk[--ns];
                    if(x - 1 >= -KMAX && !M[x - 1 + KMAX]) VISIT(x - 1);
                    if(x + 1 <=  KMAX && !M[x + 1 + KMAX]) VISIT(x + 1);
                }
                #undef VISIT
                /* new interval list, and the runs of new expanded rows sent to the adjacent columns */
                int n = 0;
                for(int x = tmin; x <= tmax; x++) if(M[x] && (x == tmin || !M[x-1])) n++;
                if(n > C->cap){ C->cap = n; C->iv = (int16_t*)realloc(C->iv, 2*n*sizeof *C->iv); }
                n = 0;
                for(int x = tmin; x <= tmax; x++){
                    if(M[x] && (x == tmin || !M[x-1])) C->iv[2*n] = (int16_t)(x - KMAX);
                    if(M[x] && (x == tmax || !M[x+1])){ C->iv[2*n+1] = (int16_t)(x - KMAX); n++; }
                }
                C->n = n;
                for(int x = tmin; x <= tmax; x++){
                    if(M[x] != 3 || (x > tmin && M[x-1] == 3)) continue;
                    int y = x; while(y + 1 <= tmax && M[y+1] == 3) y++;
                    for(int j = 0; j < DC; j++) for(int sg = -1; sg <= 1; sg += 2){
                        k[j] += sg;
                        if(abs(k[j]) <= KMAX){
                            if(on[tid] == oc[tid]){ oc[tid] = oc[tid] ? 2*oc[tid] : 4096; ob[tid] = (item_t*)realloc(ob[tid], oc[tid]*sizeof **ob); }
                            ob[tid][on[tid]].key = pack(k, DC);
                            ob[tid][on[tid]].lo = (int16_t)(x - KMAX); ob[tid][on[tid]].hi = (int16_t)(y - KMAX); on[tid]++;
                        }
                        k[j] -= sg;
                    }
                }
                memset(M + tmin, 0, tmax >= tmin ? tmax - tmin + 1 : 0);
            }
            free(stk);
        }
        for(long q = 0; q < ng; q++){                     /* fixed order */
            gres_t *r = &GR[q];
            nadd(&S, &Sc, r->s); nadd(&S, &Sc, r->c); A += r->a; edge += r->e;
            nev += r->nev; nfb += r->nfb; nexp += r->nexp; zhit |= r->zhit; if(r->maxk > maxk) maxk = r->maxk;
        }
        free(W); nw = 0; for(int t = 0; t < nth; t++) nw += on[t];
        W = (item_t*)malloc((nw ? nw : 1)*sizeof *W);
        { long at = 0; for(int t = 0; t < nth; t++){ if(on[t]) memcpy(W + at, ob[t], on[t]*sizeof *W); at += on[t]; free(ob[t]); } }
        free(ob); free(on); free(oc); free(GR); free(gs); free(gi);
        if(PROG > 0 && wtime() - TLAST > PROG){              /* progress on stderr only */
            TLAST = wtime();
            fprintf(stderr, "  [%6.0fs]  %s h=%g eps=%.0e generation %d: %ld nodes, %ld columns, next %ld items, sum %+.12f\n",
                    wtime() - T0, g->name, h, eps, gen, nev, NCOL, nw, (S + Sc)*g->weight);
        }
    }
    free(W); free(MK);
    R.ncol = NCOL;
    for(long i = 0; i < NCOL; i++) free(COLS[i].iv);
    free(COLS); COLS = NULL; NCOL = CAPCOL = 0; free(HK); free(HV); HK = NULL; HV = NULL; HCAP = HN = 0;
    R.S = S + Sc; R.A = A; R.edge = edge; R.nev = nev; R.nfb = nfb; R.nexp = nexp; R.zhit = zhit; R.zreach = maxk*h; R.t = wtime() - t0;
    return R;
}

/* exact values (coefficients of (alpha/pi)^L, mirror graph included) */
#define A1 0.5
#define A2 (-0.328478965579193)

int demo(int argc, char **argv){
    if(argc < 3){ fprintf(stderr, "demo h[,h,...] [eps|auto] [GRAPH...]   (default graphs: all1 all2)\n"); return 1; }
    double hs[64]; int nh = 0;
    for(char *t = strtok(argv[2], ","); t && nh < 64; t = strtok(NULL, ",")) hs[nh++] = atof(t);
    int a = 3; double eps = -1;
    if(a < argc && !find(argv[a]) && !isgroup(argv[a])){ if(strcmp(argv[a], "auto")) eps = atof(argv[a]); a++; }
    double epsc = getenv("SINCQED_EPSC") ? atof(getenv("SINCQED_EPSC")) : 1e-10;
    double zmax = getenv("SINCQED_ZMAX") ? atof(getenv("SINCQED_ZMAX")) : 300;
    const graph_t *sel[64]; int ns = 0;
    const char *defs[] = {"all1", "all2"};
    int nsrc = a < argc ? argc - a : 2;
    for(int b = 0; b < nsrc; b++){
        const char *s = a < argc ? argv[a + b] : defs[b]; int hit = 0;
        for(int i = 0; i < NALL && ns < 64; i++){
            if(!selected(ALL[i], s)) continue;
            hit = 1;
            if(ALL[i]->L > 2){ continue; }
            if(!ALL[i]->ud && !ALL[i]->ptd){ say("# %s: no unreduced integrand in this build, skipped\n", ALL[i]->name); continue; }
            int dup = 0; for(int j = 0; j < ns; j++) dup |= sel[j] == ALL[i];
            if(!dup) sel[ns++] = ALL[i];
        }
        if(!hit) fprintf(stderr, "# no graph matches %s in this build\n", s);
    }
    if(!ns){ fprintf(stderr, "demo: no one- or two-loop graphs selected (build with the g2 class)\n"); return 1; }
    const char *lf = getenv("SINCQED_LOG") ? getenv("SINCQED_LOG") : "sincqed.log";
    const char *rf = getenv("SINCQED_DEMO_RESULTS") ? getenv("SINCQED_DEMO_RESULTS") : "sincqed_demo.txt";
    if(*lf && strcmp(lf, "-")) LOGF = fopen(lf, "a");
    FILE *rs = NULL;
    if(*rf && strcmp(rf, "-")){
        FILE *t = fopen(rf, "r"); int fresh = !t; if(t) fclose(t);
        rs = fopen(rf, "a");
        if(rs && fresh) fprintf(rs, "# name    L D w      h       eps  value                   edge       sum|terms|  evals        fallback  zreach  zhit  time_s\n");
    }
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
#endif
    T0 = TLAST = wtime();
    say("# sincqed demo: plain Sinc sum over all free Schwinger parameters (no reductions), %d graphs, %d steps, eps %s, |z| <= %g, threads=%d\n",
        ns, nh, eps < 0 ? "auto" : argv[3], zmax, nth);
    if(eps < 0) say("# eps(h) = %g exp(-pi^2/h), floor 1e-20\n", epsc);
    for(int ih = 0; ih < nh; ih++){
        double h = hs[ih], e = eps >= 0 ? eps : fmax(epsc*exp(-M_PI*M_PI/h), 1e-20);
        double tot[3] = {0}; int have[3] = {0};
        for(int j = 0; j < ns; j++){
            const graph_t *g = sel[j];
            CURNAME = g->name;
            fprintf(stderr, "  [%6.0fs]  %s h=%g eps=%.0e started\n", wtime() - T0, g->name, h, e);
            dres_t R = demo_sum(g, h, e, zmax);
            double w = g->weight, v = w*R.S;
            say("%-7s %-26s L=%d D=%d w=%d  h=%.4f eps=%.1e  %+.15f  edge %+.1e  evals %ld (expanded %ld, fallback %ld, columns %ld)  |z|<=%.1f%s  %.1fs",
                g->name, g->volkov, g->L, g->ud ? g->Du : g->n - 1, g->weight, h, e, v, w*R.edge, R.nev, R.nexp, R.nfb, R.ncol, R.zreach,
                R.zhit ? " (box limit reached)" : "", R.t);
            const ref_t *r = findref(g->volkov);
            if(r) say("  ref %+.12f diff %+.2e", r->val, v - r->val);
            say("\n");
            if(rs){ fprintf(rs, "%-8s %d %d %d %7.4f %9.2e %+.16e %+.3e %.3e %12ld %9ld %6.1f %d %8.1f\n",
                            g->name, g->L, g->ud ? g->Du : g->n - 1, g->weight, h, e, v, w*R.edge, w*R.A, R.nev, R.nfb, R.zreach, R.zhit, R.t); fflush(rs); }
            if(g->L <= 2){ tot[g->L] += v; have[g->L]++; }
        }
        if(have[1]) say("# h=%.4f  one loop (%d graph):  %+.15f  exact %+.15f  diff %+.3e\n", h, have[1], tot[1], A1, tot[1] - A1);
        if(have[2]) say("# h=%.4f  two loops (%d graphs): %+.15f  exact %+.15f  diff %+.3e%s\n", h, have[2], tot[2], A2, tot[2] - A2,
                        have[2] == 5 ? "" : "  (not all five two-loop graphs)");
        if(rs){
            if(have[1]) fprintf(rs, "# SUM L=1 h=%.4f eps=%.2e %+.16e diff %+.3e\n", h, e, tot[1], tot[1] - A1);
            if(have[2]) fprintf(rs, "# SUM L=2 h=%.4f eps=%.2e %+.16e diff %+.3e%s\n", h, e, tot[2], tot[2] - A2, have[2] == 5 ? "" : " (partial)");
            fflush(rs);
        }
    }
    if(rs) fclose(rs);
    if(LOGF) fclose(LOGF);
    return 0;
}
