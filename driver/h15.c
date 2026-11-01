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
res_t hsum15(const graph_t *g, double h, int N, double eps, double r){
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
