/* eval.c -- the integrand at lattice nodes.  Double precision first; a point whose terms cancel
 * beyond double precision (|f| sum of |terms| * 16 eps > ABSERR) is recomputed in the high-precision
 * type QREAL (double-double by default).  eval_batch is the unit a GPU back end would replace. */
#include "sq.h"
#include <float.h>

double ABSERR = 1e-15, ABSERR0 = 1e-15, ABREL = 1e-3;   /* ABSERR = max(ABSERR0, ABREL*eps*Ref) once Ref is known (sum.c) */
__thread int SQ_AB = 1;
/* skipping the |.| sums: a point is "near" the fallback when ab*w*16 eps > NEARF*ABSERR; the children
 * of cells that were not near are evaluated without the |.| sums (and so never fall back).
 * SINCQED_ABSKIP=0 turns this off; SINCQED_ABCHECK=1 evaluates them in full anyway and counts the
 * points that would have needed the fallback (ABMISS, printed in the log). */
/* DDSLACK: the test above bounds the rounding error by 16 eps sum|terms|, which overestimates the
 * actual double error by 1e3-1e5 on most points (measured 2026-10-10); the fallback fires when that
 * bound exceeds DDSLACK*ABSERR.  SINCQED_DDSLACK=100 removed 91-100% of the fallbacks with an added
 * error far below the pruning error; the default 1 keeps earlier runs reproducible. */
double DDSLACK = 1;
double NEARF = 1e-3; int ABSKIP = 1, ABCHECK = 0; long ABMISS = 0, ABSKIPPED = 0;

double leafx(const int *k, int skip, long *nq, double *unres, int *near){
    double sf[MAXD], w, ab;
    lattice_node(k, sf, &w);
    *near = 0;
    if(skip && !ABCHECK){ SQ_AB = 0; double v = FD(sf, &ab)*w; SQ_AB = 1; return isfinite(v) ? v : 0; }
    double v = FD(sf, &ab);
    double bnd = ab*w*16*DBL_EPSILON/DDSLACK;
    if(!(bnd <= NEARF*ABSERR)) *near = 1;
    if(!(bnd <= ABSERR)){
        if(skip){
            #pragma omp atomic
            ABMISS++;
        }
        QREAL sq[MAXD], abq; for(int i = 0; i < DIM; i++) sq[i] = (QREAL)sf[i];
        QREAL vq = FQ(sq, &abq); (*nq)++;
        double e = (double)abq*w*16*QEPS; if(!(e <= ABSERR)){ e /= 16; *unres += isfinite(e) ? e*e : 0; }   /* estimate eps*sum|terms| per point, added in quadrature */
        v = isfinite((double)vq) ? (double)vq : 0; return v*w;
    }
    v *= w; return isfinite(v) ? v : 0;
}
double leaf(const int *k, long *nq, double *unres){ int nr; return leafx(k, 0, nq, unres, &nr); }

/* n nodes; runs inside the caller's parallel region (orphaned omp for) or serially */
void eval_batch(long n, const int (*k)[MAXD], const unsigned char *skip, unsigned char *near, double *out, long *nq, double *unres){
    long q = 0; double u = 0;
    #pragma omp for schedule(dynamic,64)   /* q, u are per thread here */
    for(long i = 0; i < n; i++){ int nr; out[i] = leafx(k[i], skip && skip[i], &q, &u, &nr); if(near) near[i] = (unsigned char)nr; }
    #pragma omp atomic
    *nq += q;
    #pragma omp atomic
    *unres += u;
}

#ifdef SQ_GPU
#include "gpu.h"
int USEGPU = 0;
#ifdef SQ_DDV
const sq_ddv_t *FQV;                            /* this graph's SIMD double-double fallback, or NULL (lattice.c) */
#endif
/* eval_batch with the double half on the GPU, in two steps so the caller can overlap them:
 * eval_gpu_submit starts a batch in a slot (0 or 1) and returns; eval_gpu_finish waits for it and applies
 * leafx's tests on the host: the near flag, and the double-double fallback (in parallel on the CPU cores)
 * for the points whose terms cancel beyond double precision.  Same values as eval_batch up to the
 * rounding of the double evaluation.  Call both outside a parallel region; k[i] has stride MAXD. */
void eval_gpu_submit(int slot, long n, const int (*k)[MAXD], const unsigned char *skip){
    gpu_submit(slot, n, &k[0][0], MAXD, ABCHECK ? NULL : skip);
}
/* the double-double fallback for nodes FB[0..nf-1] of a batch, at coordinates KF[j]: results into out[FB[j]].
 * SIMD double-double (ddvk.h): FQV->w points per call, groups spread over the cores; else one at a time
 * (the fallback points are often clustered, so one point or group per thread at a time) */
static void fallback(long nf, const long *FB, const int (*KF)[MAXD], double *out, long *nq, double *unres){
    long q = 0; double u = 0;
#ifdef SQ_DDV
    if(nf && FQV){
        int W = FQV->w; long ng = (nf + W - 1)/W;
        #pragma omp parallel for schedule(dynamic,1) reduction(+:u)
        for(long gi = 0; gi < ng; gi++){
            double sf[16*MAXD], ww[16], vv[16], ab[16]; long j0 = gi*W; int m = nf - j0 < W ? (int)(nf - j0) : W;
            for(int l = 0; l < m; l++) lattice_node(KF[j0 + l], sf + l*DIM, &ww[l]);
            FQV->fv(m, sf, vv, ab);
            for(int l = 0; l < m; l++){
                double e = ab[l]*ww[l]*16*QEPS; if(!(e <= ABSERR)){ e /= 16; u += isfinite(e) ? e*e : 0; }
                double v = isfinite(vv[l]) ? vv[l] : 0; out[FB[j0 + l]] = v*ww[l];
            }
        }
        q = nf; nf = 0;
    }
#endif
    #pragma omp parallel for schedule(dynamic,1) reduction(+:q,u)
    for(long j = 0; j < nf; j++){
        double sf[MAXD], w; lattice_node(KF[j], sf, &w);
        QREAL sq[MAXD], abq; for(int d = 0; d < DIM; d++) sq[d] = (QREAL)sf[d];
        QREAL vq = FQ(sq, &abq); q++;
        double e = (double)abq*w*16*QEPS; if(!(e <= ABSERR)){ e /= 16; u += isfinite(e) ? e*e : 0; }
        double v = isfinite((double)vq) ? (double)vq : 0; out[FB[j]] = v*w;
    }
    *nq += q; *unres += u;
}
static long *FB, fcap; static int (*KF)[MAXD]; static long kfcap;
static void fbgrow(long n){ if(n > fcap){ free(FB); fcap = n; FB = (long*)malloc(fcap*sizeof(long)); } }
static void kfgrow(long n){ if(n > kfcap){ free(KF); kfcap = n; KF = (int(*)[MAXD])malloc(kfcap*sizeof *KF); } }

void eval_gpu_finish(int slot, long n, const int (*k)[MAXD], const unsigned char *skip, unsigned char *near, double *out, long *nq, double *unres){
    const double *FW, *AW; gpu_wait(slot, &FW, &AW);
    if(n <= 0) return;
    /* the double results and the tests in parallel; then the fallback for the marked points */
    fbgrow(n);
    long nf = 0, miss = 0;
    #pragma omp parallel for schedule(static) reduction(+:miss)
    for(long i = 0; i < n; i++){
        double v = FW[i], abw = AW[i]; int sk = skip && skip[i];
        if(near) near[i] = 0;
        if(sk && !ABCHECK){ out[i] = isfinite(v) ? v : 0; continue; }
        double bnd = abw*16*DBL_EPSILON/DDSLACK;                       /* as leafx() */
        if(near && !(bnd <= NEARF*ABSERR)) near[i] = 1;
        if(!(bnd <= ABSERR)){ if(sk) miss++; out[i] = NAN; }   /* marked: fallback */
        else out[i] = isfinite(v) ? v : 0;
    }
    for(long i = 0; i < n; i++) if(out[i] != out[i]) FB[nf++] = i;
    kfgrow(nf); for(long j = 0; j < nf; j++) memcpy(KF[j], k[FB[j]], sizeof KF[j]);
    fallback(nf, FB, (const int(*)[MAXD])KF, out, nq, unres);
    ABMISS += miss;
}

/* the same for the children of nc cells made on the GPU (gpu_submit_cells): out[c << DIM | b] is child b of
 * cell c (k = cell - st b), 0 for a child below the box; skip is per cell.  No near flags (last level). */
void eval_gpu_submit_cells(int slot, long nc, const int (*cells)[MAXD], const unsigned char *skip, int st){
    gpu_submit_cells(slot, nc, &cells[0][0], MAXD, ABCHECK ? NULL : skip, st);
}
void eval_gpu_finish_cells(int slot, long nc, const int (*cells)[MAXD], const unsigned char *skip, int st, double *out, long *nq, double *unres){
    const double *FW, *AW; gpu_wait(slot, &FW, &AW);
    long n = nc << DIM; if(nc <= 0) return;
    fbgrow(n);
    static unsigned char *HAS; static long hcap; if(nc > hcap){ free(HAS); hcap = nc; HAS = (unsigned char*)malloc(hcap); }
    long nf = 0, miss = 0;
    #pragma omp parallel for schedule(static) reduction(+:miss)
    for(long c = 0; c < nc; c++){
        int bad = 0; for(int d = 0; d < DIM; d++) if(cells[c][d] - st < -K) bad |= 1 << d;
        int sk = skip && skip[c], has = 0;
        for(long b = 0; b < (1L << DIM); b++){
            long i = (c << DIM) | b; double v = FW[i], abw = AW[i];
            if(b & bad){ out[i] = 0; continue; }
            if(sk && !ABCHECK){ out[i] = isfinite(v) ? v : 0; continue; }
            if(!(abw*16*DBL_EPSILON/DDSLACK <= ABSERR)){ if(sk) miss++; out[i] = NAN; has = 1; }   /* as leafx() */
            else out[i] = isfinite(v) ? v : 0;
        }
        HAS[c] = (unsigned char)has;
    }
    for(long c = 0; c < nc; c++) if(HAS[c]) for(long i = c << DIM; i < (c + 1) << DIM; i++) if(out[i] != out[i]) FB[nf++] = i;
    kfgrow(nf);
    for(long j = 0; j < nf; j++){ long c = FB[j] >> DIM, b = FB[j] & ((1L << DIM) - 1);
        for(int d = 0; d < DIM; d++) KF[j][d] = cells[c][d] - (((b >> d) & 1) ? st : 0); }
    fallback(nf, FB, (const int(*)[MAXD])KF, out, nq, unres);
    ABMISS += miss;
}
void eval_batch_gpu(long n, const int (*k)[MAXD], const unsigned char *skip, unsigned char *near, double *out, long *nq, double *unres){
    const long CH = 1L << 22;
    for(long i = 0; i < n; i += CH){
        long m = n - i < CH ? n - i : CH;
        eval_gpu_submit(0, m, k + i, skip ? skip + i : NULL);
        eval_gpu_finish(0, m, k + i, skip ? skip + i : NULL, near ? near + i : NULL, out + i, nq, unres);
    }
}
#endif
