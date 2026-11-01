/* h15.c -- experimental two-level sum with the coarse level at step 1.5h (SINCQED_C15=r).
 *
 * The coarse lattice c in Z^D sits at z = 1.5 h c.  A support-following grow from c = 0 (as grow() in
 * sum.c: keep |F| >= eps*1e-3 * running sum|F|, expand axis neighbours of kept points), the final
 * threshold re-applied, then the kept set dilated to its axis neighbours.  The fine points (step h,
 * z = h k) evaluated are every k with |k_i - 1.5 c_i| <= r on every axis for some kept c, each counted
 * once (r = 0.5: the nearest coarse point, so every fine point has one parent within half a step;
 * r = 1: overlapping neighbourhoods, deduplicated).  The value is the sum over all those fine points.
 * No shadow tiers, symmetry fold, predictor or checkpoints. */
#include "sq.h"
#include <stdint.h>

typedef unsigned __int128 u128;
typedef struct { u128 *k; long cap, n; } hset_t;
static u128 key(const int *k){ u128 y = 1; for(int i = 0; i < DIM; i++) y = (y << 7) | (u128)(unsigned)(k[i] + 64); return y; }
static void unkey(u128 y, int *k){ for(int i = DIM - 1; i >= 0; i--){ k[i] = (int)(y & 127) - 64; y >>= 7; } }
static inline long hpos(u128 y, long mask){ uint64_t a = (uint64_t)y ^ (uint64_t)(y >> 64) * 0x9E3779B97F4A7C15ULL; a ^= a >> 29; a *= 0xBF58476D1CE4E5B9ULL; a ^= a >> 32; return (long)(a & mask); }
static int hins(hset_t *s, u128 y){           /* 1 if new */
    if(2*(s->n + 1) > s->cap){
        long oc = s->cap; u128 *o = s->k; s->cap = oc ? 2*oc : 1 << 16; s->k = (u128*)calloc(s->cap, sizeof(u128)); s->n = 0;
        for(long i = 0; i < oc; i++) if(o[i]) hins(s, o[i]);
        free(o);
    }
    long m = s->cap - 1, i = hpos(y, m);
    while(s->k[i]){ if(s->k[i] == y) return 0; i = (i + 1) & m; }
    s->k[i] = y; s->n++; return 1;
}
static int hhas(const hset_t *s, u128 y){
    if(!s->cap) return 0; long m = s->cap - 1, i = hpos(y, m);
    while(s->k[i]){ if(s->k[i] == y) return 1; i = (i + 1) & m; }
    return 0;
}
/* F*(1.5h)^D at coarse point c, with the double-double fallback as in eval.c */
static double coarse_val(const int *c, long *nq){
    double sf[MAXD], w = 1, ab;
    for(int i = 0; i < DIM; i++){ sf[i] = exp(1.5*H*c[i]); w *= sf[i]*1.5*H; }
    double v = FD(sf, &ab);
    if(!(ab*w*16*DBL_EPSILON <= ABSERR)){
        QREAL sq[MAXD], abq; for(int i = 0; i < DIM; i++) sq[i] = (QREAL)sf[i];
        QREAL vq = FQ(sq, &abq); (*nq)++; v = (double)vq;
        if(!(abq*QEPS <= RABS(vq))) v = 0;
    }
    v *= w; return isfinite(v) ? v : 0;
}

static res_t sumlist(const graph_t *g, const char *fn, double eps, double ref){
    res_t R; memset(&R, 0, sizeof R); double t0 = wtime();
    FILE *fp = fopen(fn, "rb"); if(!fp){ perror(fn); exit(1); }
    fseek(fp, 0, SEEK_END); long nb = ftell(fp); rewind(fp); long n = nb/(DIM*sizeof(int));
    int (*fk)[MAXD] = (int(*)[MAXD])calloc(n, sizeof *fk);
    for(long i = 0; i < n; i++) if(fread(fk[i], sizeof(int), DIM, fp) != (size_t)DIM){ perror("read"); exit(1); }
    fclose(fp);
    R.ref = ref; ABSERR = fmax(ABSERR0, ABREL*eps*ref);
    double s = 0, c = 0, a = 0, un = 0; long q2 = 0; const long B = 1 << 16; double *vb = (double*)malloc(B*sizeof(double));
    for(long b0 = 0; b0 < n; b0 += B){
        long b1 = b0 + B < n ? b0 + B : n;
        #pragma omp parallel for reduction(+:q2,un) schedule(dynamic,64)
        for(long i = b0; i < b1; i++) vb[i - b0] = leaf(fk[i], &q2, &un);
        for(long i = b0; i < b1; i++){ double x = vb[i - b0], t = s + x; if(fabs(s) >= fabs(x)) c += (s - t) + x; else c += (x - t) + s; s = t; a += fabs(x); }
    }
    free(vb); free(fk);
    R.S[0] = s + c; R.Ah = a; R.nev = n; R.nq = q2; R.unres = sqrt(un); R.S[1] = R.S[2] = NAN; R.t = wtime() - t0;
    say("# %s list %s: %ld points\n", g->name, fn, n);
    return R;
}

/* SINCQED_POLY=file: sum F over every lattice point k with L k <= c (file: "nf D" then nf rows of D ints and c).
 * Enumerated coordinate by coordinate: the range of k_m is cut by every facet with the later coordinates relaxed
 * to the bounding box, and the last coordinate's range is exact. */
typedef struct { int nf; int *L; long *c; int (*buf)[MAXD]; long nb, n; double s, cc, a, un; long q2; } poly_t;
static void pflush(poly_t *P){
    double *vb = (double*)malloc(P->nb*sizeof(double)); long q2 = 0; double un = 0;
    #pragma omp parallel for reduction(+:q2,un) schedule(dynamic,64)
    for(long i = 0; i < P->nb; i++) vb[i] = leaf(P->buf[i], &q2, &un);
    for(long i = 0; i < P->nb; i++){ double x = vb[i], t = P->s + x; if(fabs(P->s) >= fabs(x)) P->cc += (P->s - t) + x; else P->cc += (x - t) + P->s; P->s = t; P->a += fabs(x); }
    P->q2 += q2; P->un += un; P->n += P->nb; P->nb = 0; free(vb);
    progress("polytope", P->n, 0);
}
static long *RMIN;   /* RMIN[m*nf + f] = sum over i >= m of min(L_fi lo_i, L_fi hi_i) */
static int LO[MAXD], HI[MAXD];
static void pscan(poly_t *P, int *k, int m, const long *dot, long *work){
    int nf = P->nf; long lo = LO[m], hi = HI[m];
    for(int f = 0; f < nf; f++){
        int l = P->L[f*DIM + m]; if(!l) { if(dot[f] + RMIN[(m + 1)*nf + f] > P->c[f]) return; continue; }
        long r = P->c[f] - dot[f] - RMIN[(m + 1)*nf + f];      /* l k_m <= r */
        if(l > 0){ long u = r >= 0 ? r/l : -((-r + l - 1)/l); if(u < hi) hi = u; }
        else { long L2 = -l; long w = -r; long u = w >= 0 ? (w + L2 - 1)/L2 : -((-w)/L2); if(u > lo) lo = u; }   /* k_m >= ceil(-r/|l|) */
        if(lo > hi) return;
    }
    long *nd = work;
    for(long x = lo; x <= hi; x++){
        k[m] = (int)x;
        if(m == DIM - 1){ memcpy(P->buf[P->nb++], k, DIM*sizeof(int)); if(P->nb == (1 << 16)) pflush(P); continue; }
        for(int f = 0; f < nf; f++) nd[f] = dot[f] + (long)P->L[f*DIM + m]*x;
        pscan(P, k, m + 1, nd, work + nf);
    }
}

/* BFS version (SINCQED_POLYBFS=1): layers of the lattice graph from 0; a layer's neighbours lie in the
 * previous, current or next layer, so only two layers are kept in hash sets. */
static void pbfs(poly_t *P){
    int nf = P->nf;
    int (*cur)[MAXD] = (int(*)[MAXD])calloc(1, sizeof *cur); long nc = 1;
    hset_t hp = {0}, hc = {0}; { int z[MAXD] = {0}; hins(&hc, key(z)); }
    memcpy(P->buf[P->nb++], cur[0], DIM*sizeof(int));
    while(nc){
        long cap = 1024, nn = 0; int (*nx)[MAXD] = (int(*)[MAXD])malloc(cap*sizeof *nx); hset_t hn = {0};
        for(long i = 0; i < nc; i++) for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
            int q[MAXD]; memcpy(q, cur[i], sizeof q); q[d] += sg; if(abs(q[d]) > K) continue;
            u128 y = key(q); if(hhas(&hp, y) || hhas(&hc, y)) continue;
            if(!hins(&hn, y)) continue;
            if(nn == cap){ cap *= 2; nx = (int(*)[MAXD])realloc(nx, cap*sizeof *nx); } memcpy(nx[nn++], q, sizeof *nx);
        }
        char *ok = (char*)malloc(nn ? nn : 1);
        #pragma omp parallel for schedule(dynamic,256)
        for(long i = 0; i < nn; i++){
            int in = 1; for(int f = 0; f < nf && in; f++){ long t = 0; const int *l = P->L + f*DIM; for(int j = 0; j < DIM; j++) t += (long)l[j]*nx[i][j]; if(t > P->c[f]) in = 0; }
            ok[i] = (char)in;
        }
        long m = 0; hset_t hk = {0};
        for(long i = 0; i < nn; i++) if(ok[i]){ memcpy(nx[m++], nx[i], sizeof *nx); hins(&hk, key(nx[i])); }
        free(ok); free(hn.k);
        for(long i = 0; i < m; i++){ memcpy(P->buf[P->nb++], nx[i], DIM*sizeof(int)); if(P->nb == (1 << 16)) pflush(P); }
        free(hp.k); hp = hc; hc = hk; free(cur); cur = nx; nc = m;
    }
    free(cur); free(hp.k); free(hc.k);
}
static res_t sumpoly(const graph_t *g, const char *fn, double eps, double ref){
    res_t R; memset(&R, 0, sizeof R); double t0 = wtime();
    FILE *fp = fopen(fn, "r"); if(!fp){ perror(fn); exit(1); }
    poly_t P; memset(&P, 0, sizeof P); int d;
    if(fscanf(fp, "%d %d", &P.nf, &d) != 2 || d != DIM){ fprintf(stderr, "SINCQED_POLY: bad header\n"); exit(1); }
    P.L = (int*)malloc(P.nf*DIM*sizeof(int)); P.c = (long*)malloc(P.nf*sizeof(long));
    for(int f = 0; f < P.nf; f++){ for(int i = 0; i < DIM; i++) if(fscanf(fp, "%d", &P.L[f*DIM + i]) != 1) exit(1); if(fscanf(fp, "%ld", &P.c[f]) != 1) exit(1); }
    fclose(fp);
    R.ref = ref; ABSERR = fmax(ABSERR0, ABREL*eps*ref);
    P.buf = (int(*)[MAXD])malloc((1 << 16)*sizeof *P.buf);
    long *dot = (long*)calloc(P.nf, sizeof(long)), *work = (long*)malloc((size_t)P.nf*(DIM*(K + 1) + 2)*sizeof(long));
    /* box from the +-e_i facets */
    for(int i = 0; i < DIM; i++){ LO[i] = -K; HI[i] = K; }
    for(int f = 0; f < P.nf; f++){ int nz = 0, ii = -1; for(int i = 0; i < DIM; i++) if(P.L[f*DIM + i]){ nz++; ii = i; }
        if(nz == 1){ int l = P.L[f*DIM + ii]; if(l == 1 && P.c[f] < HI[ii]) HI[ii] = (int)P.c[f]; if(l == -1 && -P.c[f] > LO[ii]) LO[ii] = (int)-P.c[f]; } }
    RMIN = (long*)calloc((size_t)(DIM + 1)*P.nf, sizeof(long));
    for(int m = DIM - 1; m >= 0; m--) for(int f = 0; f < P.nf; f++){ long l = P.L[f*DIM + m]; RMIN[m*P.nf + f] = RMIN[(m + 1)*P.nf + f] + (l > 0 ? l*LO[m] : l*HI[m]); }
    int k[MAXD] = {0};
    if(getenv("SINCQED_POLYBFS")) pbfs(&P); else pscan(&P, k, 0, dot, work); if(P.nb) pflush(&P);
    R.S[0] = P.s + P.cc; R.Ah = P.a; R.nev = P.n; R.nq = P.q2; R.unres = sqrt(P.un); R.S[1] = R.S[2] = NAN; R.t = wtime() - t0;
    say("# %s polytope %s: %d facets, %ld points\n", g->name, fn, P.nf, P.n);
    free(P.buf); free(P.L); free(P.c); free(dot); free(work); free(RMIN);
    return R;
}
res_t hsum15(const graph_t *g, double h, int N, double eps, double r){
    { const char *fp = getenv("SINCQED_POLY"); if(fp && *fp) return sumpoly(g, fp, eps, r); }
    { const char *fl = getenv("SINCQED_PTS"); if(fl && *fl) return sumlist(g, fl, eps, r); }
    res_t R; memset(&R, 0, sizeof R); double t0 = wtime();
    int KC = (int)floor(K/1.5);
    ABSERR = ABSERR0;
    /* coarse grow */
    hset_t seen = {0}; long nkeep = 0, capk = 1024; int (*kp)[MAXD] = (int(*)[MAXD])malloc(capk*sizeof *kp); double *kv = (double*)malloc(capk*sizeof(double));
    long nf = 1, capf = 1024; int (*fr)[MAXD] = (int(*)[MAXD])calloc(capf, sizeof *fr); { int z[MAXD] = {0}; hins(&seen, key(z)); }
    double A = 0, S = 0; long nev = 0, nq = 0;
    while(nf){
        double *v = (double*)malloc(nf*sizeof(double));
        #pragma omp parallel for reduction(+:nq) schedule(dynamic,16)
        for(long i = 0; i < nf; i++) v[i] = coarse_val(fr[i], &nq);
        nev += nf;
        for(long i = 0; i < nf; i++){ A += fabs(v[i]); S += v[i]; }
        double thr = eps*1e-3*A;
        long nn = 0, capn = 1024; int (*nx)[MAXD] = (int(*)[MAXD])malloc(capn*sizeof *nx);
        for(long i = 0; i < nf; i++){
            if(!(fabs(v[i]) >= thr)) continue;
            if(nkeep == capk){ capk *= 2; kp = (int(*)[MAXD])realloc(kp, capk*sizeof *kp); kv = (double*)realloc(kv, capk*sizeof(double)); }
            memcpy(kp[nkeep], fr[i], sizeof *kp); kv[nkeep++] = fabs(v[i]);
            for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
                int c[MAXD]; memcpy(c, fr[i], sizeof c); c[d] += sg; if(abs(c[d]) > KC) continue;
                if(hins(&seen, key(c))){ if(nn == capn){ capn *= 2; nx = (int(*)[MAXD])realloc(nx, capn*sizeof *nx); } memcpy(nx[nn++], c, sizeof *nx); }
            }
        }
        free(v); free(fr); fr = nx; nf = nn;
    }
    free(fr);
    R.ref = A; R.S[1] = S;                       /* S_1.5h over the grow (not the whole coarse sum) */
    double thr = eps*1e-3*A; long m = 0;
    for(long i = 0; i < nkeep; i++) if(kv[i] >= thr){ memcpy(kp[m], kp[i], sizeof *kp); kv[m++] = kv[i]; }
    nkeep = m;
    /* dilate: kept plus axis neighbours */
    hset_t kept = {0}; long nk2 = 0, capk2 = 2*nkeep*(2*DIM+1) + 16; int (*kc)[MAXD] = (int(*)[MAXD])malloc(capk2*sizeof *kc);
    for(long i = 0; i < nkeep; i++){
        if(hins(&kept, key(kp[i]))) memcpy(kc[nk2++], kp[i], sizeof *kc);
        if(!getenv("SINCQED_C15NODIL")) for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
            int c[MAXD]; memcpy(c, kp[i], sizeof c); c[d] += sg; if(abs(c[d]) > KC) continue;
            if(hins(&kept, key(c))) memcpy(kc[nk2++], c, sizeof *kc);
        }
    }
    free(kp); free(kv); free(seen.k);
    R.kept_top = nk2;
    /* fine candidates, deduplicated */
    hset_t fine = {0}; long nfine = 0, capfi = 1 << 16; int (*fk)[MAXD] = (int(*)[MAXD])malloc(capfi*sizeof *fk);
    for(long i = 0; i < nk2; i++){
        int lo[MAXD], hi[MAXD], k[MAXD]; long tot = 1;
        for(int d = 0; d < DIM; d++){ double x = 1.5*kc[i][d]; lo[d] = (int)ceil(x - r - 1e-9); hi[d] = (int)floor(x + r + 1e-9);
            if(lo[d] < -K) lo[d] = -K; if(hi[d] > K) hi[d] = K; if(hi[d] < lo[d]){ tot = 0; break; } tot *= hi[d] - lo[d] + 1; }
        for(long j = 0; j < tot; j++){
            long q = j; for(int d = 0; d < DIM; d++){ int w = hi[d] - lo[d] + 1; k[d] = lo[d] + (int)(q % w); q /= w; }
            if(hins(&fine, key(k))){ if(nfine == capfi){ capfi *= 2; fk = (int(*)[MAXD])realloc(fk, capfi*sizeof *fk); } memcpy(fk[nfine++], k, sizeof *fk); }
        }
    }
    free(kc); free(kept.k); free(fine.k);
    ABSERR = fmax(ABSERR0, ABREL*eps*R.ref);
    /* evaluate and sum (compensated, in a fixed order) */
    double s = 0, c = 0, a = 0, un = 0; long q2 = 0;
    const long B = 1 << 16;
    double *vb = (double*)malloc(B*sizeof(double));
    for(long b0 = 0; b0 < nfine; b0 += B){
        long b1 = b0 + B < nfine ? b0 + B : nfine;
        #pragma omp parallel for reduction(+:q2,un) schedule(dynamic,64)
        for(long i = b0; i < b1; i++) vb[i - b0] = leaf(fk[i], &q2, &un);
        for(long i = b0; i < b1; i++){ double x = vb[i - b0], t = s + x; if(fabs(s) >= fabs(x)) c += (s - t) + x; else c += (x - t) + s; s = t; a += fabs(x); }
        progress("fine level (c15)", b1, nfine);
    }
    free(vb); free(fk);
    R.S[0] = s + c; R.Ah = a; R.nev = nev + nfine; R.nq = nq + q2; R.unres = sqrt(un); R.S[2] = NAN;
    R.t = wtime() - t0;
    say("# %s c15 r=%g: coarse points evaluated %ld, kept+dilated %ld, fine points %ld\n", g->name, r, nev, nk2, nfine);
    return R;
}
