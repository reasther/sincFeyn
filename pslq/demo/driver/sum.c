/* sum.c -- the hierarchical Sinc sum.
 *
 * Coarsest level (step 2^top h): a support-following scan from z = 0 (grow), or a scan of the whole
 * box.  Cells above eps*Ref*1e-3 are kept and dilated to their axis neighbours; each finer level
 * evaluates the 2^D children of the kept cells, keeping those above eps*Ref*1e-2.  S_h, S_2h, S_4h
 * are the sums at the last three levels.  The refinement runs in blocks of cells; each block's
 * nodes are evaluated as one batch and summed in a fixed order (so a run is reproducible), and the
 * state is checkpointed between blocks (checkpoint.c). */
#include "sq.h"

int GROW = 1;
static double TQ[4], TGROW;              /* SINCQED_TIMING: seconds in the refinement's children / pack / eval / sum steps */
static long DONE; static char STAGE[48];

/* SINCQED_HIST=prefix: per level, a histogram of the term sizes |f w| in half-decade bins (count,
 * signed and absolute contribution), written to <prefix>_<graph>.txt when the graph finishes */
#define HB 80                                   /* bins: log10|term| from -35 to +5 */
#define RB 64                                   /* radial bins: max_i |z_i| = h max|k_i| in unit bins */
static double HC[16][HB], HS_[16][HB], HA[16][HB], RC[16][RB], RS[16][RB], RA[16][RB]; static int HON = -1;
static inline void hadd(int lev, double t, const int *k){
    if(HON <= 0 || t == 0) return;
    int b = (int)floor(2*(log10(fabs(t)) + 35)); if(b < 0) b = 0; if(b >= HB) b = HB - 1;
    HC[lev][b] += 1; HS_[lev][b] += t; HA[lev][b] += fabs(t);
    int mx = 0; for(int d = 0; d < DIM; d++) if(abs(k[d]) > mx) mx = abs(k[d]);
    int r = (int)(mx*H); if(r >= RB) r = RB - 1;
    RC[lev][r] += 1; RS[lev][r] += t; RA[lev][r] += fabs(t);
}
/* SINCQED_PREDICT=alpha: at the last level, child k - h b of a kept 2h point k (the centre of the 2h cube with
 * corners k - 2h c, c a subset of b) is evaluated only if the mean of |f| over those corners is at least
 * alpha*thr; a corner never evaluated counts as thr.  The 2h values come from a hash of every point evaluated
 * at that level (lost on a checkpoint resume, which turns the prediction off).  SINCQED_PREDDIAG=1 evaluates
 * the points it would skip anyway, sums them as usual and reports their total. */
/* Shadow runs: every kept point also carries a tier t in {0,1,2} (bits 1-2 of its PFLAG slot): it would
 * still be kept, and its children evaluated, with eps raised by 10^t.  At the final level the points
 * of tier >= 1 and >= 2 are summed separately, so one run reports what eps*10 and eps*100 would give:
 * D1 = S(eps) - S(10 eps), D2 = S(10 eps) - S(100 eps), and from their ratio r an estimate D1/(r-1) of the
 * loss still left at eps.  Top level: tier from the point's own |f| (the grow's path is not tracked). */
static double ES1, ES2; static int ESOK;
#define TIER(f, t) (((f) >= 10*(t)) + ((f) >= 100*(t)))
static double PRED = -1; long PSKIP = 0; static int PDIAG = -1; static double PDS, PDA, PDN, *ST;
/* the 2h values: open addressing on 16-byte slots, so a lookup touches one cache line */
typedef struct { unsigned long long k; float v; unsigned char t; } vslot_t;   /* t: shadow tier up to which the point is evaluated */
static vslot_t *VS; static long VCAP, VN;
static unsigned long long vhash(const int *k){
    unsigned long long y = 0x9E3779B97F4A7C15ULL;
    for(int i = 0; i < DIM; i++){ y ^= (unsigned long long)(unsigned)k[i] + 0x9E3779B97F4A7C15ULL + (y << 6) + (y >> 2); y *= 0xBF58476D1CE4E5B9ULL; }
    return y | 1;
}
static inline long vpos(unsigned long long y){ return (long)(y >> 17) & (VCAP - 1); }
static void vput_h(unsigned long long y, float v, unsigned char t){
    long i = vpos(y);
    while(VS[i].k && VS[i].k != y) i = (i + 1) & (VCAP - 1);
    if(!VS[i].k){ VS[i].k = y; VN++; } VS[i].v = v; VS[i].t = t;
}
static void vput(const int *k, double a, int t){
    if(VN*2 >= VCAP){
        long oc = VCAP; vslot_t *os = VS;
        VCAP = oc ? 2*oc : 1 << 20; VS = (vslot_t*)calloc(VCAP, sizeof(vslot_t)); VN = 0;
        for(long i = 0; i < oc; i++) if(os[i].k) vput_h(os[i].k, os[i].v, os[i].t);
        free(os);
    }
    vput_h(vhash(k), (float)log(a > 1e-300 ? a : 1e-300), (unsigned char)t);
}
static inline int vget_h(unsigned long long y, double *lv){
    long i = vpos(y);
    while(VS[i].k){ if(VS[i].k == y){ *lv = VS[i].v; return 1 + VS[i].t; } i = (i + 1) & (VCAP - 1); }
    return 0;
}
static int vget(const int *k, double *lv){ return VCAP ? vget_h(vhash(k), lv) : 0; }
static void vraise(const int *k, int t){           /* the point is evaluated in the eps*10^t run too (thread-safe) */
    if(!VCAP) return;
    unsigned long long y = vhash(k); long i = vpos(y);
    while(VS[i].k){ if(VS[i].k == y){ unsigned char o = __atomic_load_n(&VS[i].t, __ATOMIC_RELAXED);
            while(o < t && !__atomic_compare_exchange_n(&VS[i].t, &o, (unsigned char)t, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) ;
            return; }
        i = (i + 1) & (VCAP - 1); }
}
/* SINCQED_DUMP=prefix: every point evaluated at the coarser levels, and a 1-in-SINCQED_DUMPSAMPLE (default 1000)
 * sample of the final level chosen by a hash of k (so the same points in every run), go to <prefix>_<graph>.pts.
 * With SINCQED_PREDICT, sampled points the prediction skips are evaluated for the dump but not summed.
 * Format: the text line "SQPT1 name D h N eps top K\n", then records of
 *   uint8 level | 0x80 if skipped by the prediction;  int16 k[D] (z_i = k_i h);  double f w (before the graph weight). */
static FILE *DUMPF; static long DUMPS = 1000;
static int sampled(const int *k){ return (vhash(k) >> 20) % (unsigned long long)DUMPS == 0; }
static void dump(int lev, int skipped, const int *k, double v){
    unsigned char l = (unsigned char)(lev | (skipped ? 0x80 : 0)); short kk[MAXD];
    for(int i = 0; i < DIM; i++) kk[i] = (short)k[i];
    fwrite(&l, 1, 1, DUMPF); fwrite(kk, sizeof(short), DIM, DUMPF); fwrite(&v, sizeof v, 1, DUMPF);
}
static void hwrite(const graph_t *g, int top){
    const char *pre = getenv("SINCQED_HIST"); if(!pre || !*pre) return;
    char fn[4096]; snprintf(fn, sizeof fn, "%s_%s.txt", pre, g->name); FILE *f = fopen(fn, "w"); if(!f) return;
    fprintf(f, "# %s term-size histogram: level (step 2^lev h), bin [lo,hi) of log10|f w|, count, sum f w, sum |f w| (before the graph weight w)\n", g->name);
    for(int l = top; l >= 0; l--) for(int b = 0; b < HB; b++) if(HC[l][b] > 0)
        fprintf(f, "%d %6.1f %6.1f %.0f %+.6e %.6e\n", l, b/2.0 - 35, b/2.0 - 34.5, HC[l][b], HS_[l][b], HA[l][b]);
    fprintf(f, "# %s radial histogram: R level, bin [r,r+1) of max_i |z_i|, count, sum f w, sum |f w|\n", g->name);
    for(int l = top; l >= 0; l--) for(int r = 0; r < RB; r++) if(RC[l][r] > 0)
        fprintf(f, "R %d %2d %2d %.0f %+.6e %.6e\n", l, r, r + 1, RC[l][r], RS[l][r], RA[l][r]);
    fclose(f);
}

/* full scan at level lev (step 2^lev): returns sum, sum|.|, outer-shell weights; optionally keeps */
static void scan(int lev, double thr, vec_t *keep, double *S, double *A, double *e1, double *e2, double *drop, long *nev, long *nq, double *un){
    int st = 1 << lev, W = 2*(K/st) + 1; long n = 1; for(int i = 0; i < DIM; i++) n *= W;
    double pw = pow(st, DIM), s = 0, a = 0, x1 = 0, x2 = 0, dr = 0; long ne = 0, q = 0; double u = 0;
    DONE = 0; snprintf(STAGE, sizeof STAGE, "scan at step %dh", st);
    #pragma omp parallel
    { vec_t loc = {0};
      #pragma omp for reduction(+:s,a,x1,x2,dr,ne,q,u) schedule(dynamic,256)
      for(long id = 0; id < n; id++){
          int k[MAXD], mx = 0; long r = id;
          for(int i = 0; i < DIM; i++){ k[i] = st*((int)(r % W) - K/st); r /= W; if(abs(k[i]) > mx) mx = abs(k[i]); }
          double v = leaf(k, &q, &u); ne++;
          { long dn;
            #pragma omp atomic capture
            dn = ++DONE;
            if(!(dn & 255)) progress(STAGE, dn, n); }
          s += v*pw; a += fabs(v)*pw;
          if(mx == K) x1 += fabs(v)*pw; else if(mx == K - st) x2 += fabs(v)*pw;
          if(keep){ if(fabs(v) >= thr) push(&loc, k); else dr += fabs(v)*pw; }
      }
      if(keep){
        #pragma omp critical
        { for(long i = 0; i < loc.n; i++) push(keep, loc.p[i].k); }
      }
      free(loc.p); }
    *S = s; *A = a; *e1 = x1; *e2 = x2; *drop = dr; *nev += ne; *nq += q; *un += u;
}

/* add the axis neighbours of every cell (step st) and remove duplicates; a point's flags merge as in
 * uniq() (bit 0: AND, tier: max).  The copies used to be pushed as 48-byte points and qsorted serially:
 * at h = 0.8 that was ~15 copies of ~40M cells (30 GB) and minutes on one core.  Now each copy is a
 * packed 64-bit key in a parallel lock-free hash set (~10 B per distinct point), and the distinct keys
 * are radix-sorted, so the list order (and so the summation order) is fixed. */
/* A point's packed key is the mixed-radix number 1 c_0 c_1 ... c_{D-1} in base W = 2*(K/st) + 1, with
 * c_i = k_i/st + K/st (the leading 1 keeps 0 free as the empty slot).  An axis step on dimension j adds
 * +-P[j] = +-W^(D-1-j).  Numeric order is lexicographic order in k, as with the earlier power-of-two
 * fields, but the key fits in 64 bits while W^D <= 2^63 rather than ceil(log2 W)*D <= 63: at top step 2h
 * D = 10 now reaches N = 153 (was 121) and D = 11 N = 97 (was 57).  keyradix() fills P and
 * the key's bit length and returns W, or 0 if the keys would not fit. */
static int keyradix(int st, unsigned long long *P, int *nb){
    int W = 2*(K/st) + 1; unsigned long long p = 1;
    for(int i = DIM - 1; i >= 0; i--){ P[i] = p; if(p > (1ULL << 63)/W) return 0; p *= W; }
    *nb = 65 - __builtin_clzll(p);              /* keys < 2 W^D */
    return W;
}
static void dilate_old(vec_t *cur, int st){
    long n0 = cur->n;
    for(long q = 0; q < n0; q++) for(int i = 0; i < DIM; i++) for(int sg = -1; sg <= 1; sg += 2){
        int k[MAXD]; memcpy(k, cur->p[q].k, sizeof(int)*DIM); k[i] += st*sg; int fl = cur->p[q].k[PFLAG];
        if(k[i] > K || k[i] < -K) continue; push(cur, k); if(DIM < MAXD) cur->p[cur->n-1].k[PFLAG] = fl; }
    uniq(cur);
}
static long hslot(unsigned long long *HK, long mask, unsigned long long y, long *cnt){
    long i = (long)((y*0x9E3779B97F4A7C15ULL) >> 21) & mask;
    for(;;){
        unsigned long long c = __atomic_load_n(&HK[i], __ATOMIC_RELAXED);
        if(c == y) return i;
        if(!c){ unsigned long long z = 0;
            if(__atomic_compare_exchange_n(&HK[i], &z, y, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)){ __atomic_fetch_add(cnt, 1, __ATOMIC_RELAXED); return i; }
            if(z == y) return i; }
        i = (i + 1) & mask;
    }
}
static long hfind(const unsigned long long *HK, long mask, unsigned long long y){
    long i = (long)((y*0x9E3779B97F4A7C15ULL) >> 21) & mask;
    while(HK[i] != y) i = (i + 1) & mask;
    return i;
}
static void dilate(vec_t *cur, int st){
    unsigned long long P[MAXD]; int nb, W = keyradix(st, P, &nb);
    long n0 = cur->n;
    if(!W || n0 == 0 || getenv("SINCQED_OLDDILATE")){ dilate_old(cur, st); return; }
    long cap = 1 << 20; while(cap < 4*n0) cap <<= 1;
    unsigned long long *HK; unsigned char *SK, *TR; long cnt;
    for(;;){
        HK = (unsigned long long*)calloc(cap, 8); SK = (unsigned char*)malloc(cap); TR = (unsigned char*)calloc(cap, 1);
        memset(SK, 1, cap); cnt = 0; int full = 0; long lim = cap/4*3, mask = cap - 1;
        #pragma omp parallel for schedule(dynamic,4096)
        for(long q = 0; q < n0; q++){
            if(__atomic_load_n(&full, __ATOMIC_RELAXED)) continue;
            const int *k0 = cur->p[q].k; int fl = DIM < MAXD ? k0[PFLAG] : 0;
            unsigned char sk = (unsigned char)(fl & 1), tr = (unsigned char)(fl >> 1);
            unsigned long long y0 = 1; for(int i = 0; i < DIM; i++) y0 = y0*W + (unsigned long long)(k0[i]/st + K/st);
            for(int j = -1; j < DIM; j++) for(int sg = -1; sg <= 1; sg += 2){
                if(j < 0 && sg > 0) continue;               /* j = -1: the cell itself */
                unsigned long long y = y0;
                if(j >= 0){ int kk = k0[j] + st*sg; if(kk > K || kk < -K) continue;
                    y = sg > 0 ? y + P[j] : y - P[j]; }
                long i = hslot(HK, mask, y, &cnt);
                if(!sk) __atomic_store_n(&SK[i], 0, __ATOMIC_RELAXED);
                unsigned char o = __atomic_load_n(&TR[i], __ATOMIC_RELAXED);
                while(o < tr && !__atomic_compare_exchange_n(&TR[i], &o, tr, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) ;
                if(__atomic_load_n(&cnt, __ATOMIC_RELAXED) > lim){ __atomic_store_n(&full, 1, __ATOMIC_RELAXED); break; }
            }
        }
        if(!full) break;
        free(HK); free(SK); free(TR); cap *= 2;      /* rare: more distinct points than 3/4 of the table */
    }
    /* distinct keys, radix-sorted (LSD, 16-bit digits) */
    long m = 0, mask = cap - 1;
    unsigned long long *U = (unsigned long long*)malloc(cnt*8), *T = (unsigned long long*)malloc(cnt*8);
    for(long i = 0; i < cap; i++) if(HK[i]) U[m++] = HK[i];
    long *cn = (long*)malloc((1L << 16)*sizeof(long));
    for(int sh = 0; sh < nb; sh += 16){
        memset(cn, 0, (1L << 16)*sizeof(long));
        for(long i = 0; i < m; i++) cn[(U[i] >> sh) & 0xffff]++;
        long sum = 0; for(long d = 0; d < (1L << 16); d++){ long c = cn[d]; cn[d] = sum; sum += c; }
        for(long i = 0; i < m; i++) T[cn[(U[i] >> sh) & 0xffff]++] = U[i];
        unsigned long long *t = U; U = T; T = t;
    }
    free(cn); free(T);
    free(cur->p); cur->p = (pt_t*)malloc((m ? m : 1)*sizeof(pt_t)); cur->n = cur->cap = m;
    #pragma omp parallel for schedule(static)
    for(long j = 0; j < m; j++){
        unsigned long long y = U[j]; pt_t *p = &cur->p[j]; memset(p, 0, sizeof *p);
        for(int i = DIM - 1; i >= 0; i--){ p->k[i] = ((int)(y % W) - K/st)*st; y /= W; }
        if(DIM < MAXD){ long s = hfind(HK, mask, U[j]); p->k[PFLAG] = SK[s] | (TR[s] << 1); }
    }
    free(U); free(HK); free(SK); free(TR);
}

/* support-following scan at level lev (step 2^lev): start at z = 0 and grow through axis neighbours
 * while |f| is above thr = eps_rel * (running sum|f|); points below it are evaluated (the boundary)
 * but not expanded.  Cost ~ support size x 2D, not (box width)^D.  Same outputs as scan(); *A is
 * also the reference scale.  Returns 0 if the lattice is too large to index in 64 bits (keyradix). */
static int grow(int lev, double eps_rel, vec_t *keep, double *S, double *A, double *e1, double *e2, double *drop, long *nev, long *nq, double *un){
    int st = 1 << lev, nb; unsigned long long P[MAXD]; int W = keyradix(st, P, &nb);
    if(!W) return 0;
    long cap = 1 << 20; unsigned long long *HS = (unsigned long long*)calloc(cap, 8); long nh = 0;
    #define KEY(k) ({ unsigned long long _y = 1; for(int _i = 0; _i < DIM; _i++) _y = _y*W + (unsigned long long)((k)[_i]/st + K/st); _y; })
    vec_t fr = {0}, nx = {0};
    double pw = pow(st, DIM), s = 0, a = 0, x1 = 0, x2 = 0, dr = 0, u = 0; long ne = 0, q = 0;
    int k0[MAXD] = {0}; push(&fr, k0);
    double *kv = NULL; long nkv = 0;   /* |f| of the kept points, for the final threshold */
    #define SEEN(key) ({ int _new = 0; unsigned long long _k = (key); long _h = (long)((_k*0x9E3779B97F4A7C15ULL) >> 20) & (cap-1); \
        while(HS[_h] && HS[_h] != _k) _h = (_h + 1) & (cap-1); if(!HS[_h]){ HS[_h] = _k; nh++; _new = 1; } _new; })
    SEEN(KEY(k0));
    DONE = 0; snprintf(STAGE, sizeof STAGE, "grow at step %dh", st);
    double tw0 = wtime();
    while(fr.n){
        double *v = (double*)malloc(fr.n*sizeof(double)); unsigned char *nr = (unsigned char*)malloc(fr.n);
        if(USEGPU){ eval_batch_gpu(fr.n, (const int(*)[MAXD])fr.p, NULL, nr, v, &q, &u); DONE += fr.n; progress(STAGE, DONE, DONE); }
        else
        #pragma omp parallel for reduction(+:q,u) schedule(dynamic,64)
        for(long i = 0; i < fr.n; i++){
            int nri; v[i] = leafx(fr.p[i].k, 0, &q, &u, &nri); nr[i] = (unsigned char)nri;
            { long dn;
              #pragma omp atomic capture
              dn = ++DONE;
              if(!(dn & 255)) progress(STAGE, dn, dn + fr.n - i); }
        }
        ne += fr.n;
        for(long i = 0; i < fr.n; i++){ s += v[i]*pw; a += fabs(v[i])*pw; hadd(lev, v[i]*pw, fr.p[i].k); if(PRED > 0) vput(fr.p[i].k, fabs(v[i]), 0); if(DUMPF) dump(lev, 0, fr.p[i].k, v[i]*pw); }
        double thr = eps_rel*a;                  /* per-point threshold from the running sum (as scan) */
        nx.n = 0;
        for(long i = 0; i < fr.n; i++){
            int *k = fr.p[i].k, mx = 0; for(int d = 0; d < DIM; d++) if(abs(k[d]) > mx) mx = abs(k[d]);
            if(mx == K) x1 += fabs(v[i])*pw; else if(mx == K - st) x2 += fabs(v[i])*pw;
            if(!(fabs(v[i]) >= thr)){ dr += fabs(v[i])*pw; continue; }
            if(keep->n >= nkv){ nkv = 2*nkv + 1024; kv = (double*)realloc(kv, nkv*sizeof(double)); }
            kv[keep->n] = fabs(v[i]); push(keep, k); if(ABSKIP && DIM < MAXD) keep->p[keep->n-1].k[PFLAG] = !nr[i];
            for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
                int kk[MAXD]; memcpy(kk, k, sizeof(int)*DIM); kk[d] += st*sg;
                if(kk[d] > K || kk[d] < -K) continue;
                if(2*(nh + 1) > cap){                                   /* rehash */
                    long oc = cap; unsigned long long *O = HS; cap *= 2; HS = (unsigned long long*)calloc(cap, 8); nh = 0;
                    for(long j = 0; j < oc; j++) if(O[j]) SEEN(O[j]);
                    free(O); }
                if(SEEN(KEY(kk))) push(&nx, kk);
            }
        }
        free(v); free(nr);
        vec_t t = fr; fr = nx; nx = t;
    }
    TGROW = wtime() - tw0;
    /* the running threshold was lower early on: re-apply the final one */
    { double thr = eps_rel*a; long m = 0;
      for(long i = 0; i < keep->n; i++){ if(kv[i] >= thr){ kv[m] = kv[i]; keep->p[m++] = keep->p[i]; } else dr += kv[i]*pw; }
      keep->n = m;
      /* shadow tiers: a run with eps*10^t keeps a point only if it reaches it from z = 0 through points with
       * |f| >= 10^t thr (the grow expands nothing else); flood-fill the kept set for t = 1, 2 */
      if(DIM < MAXD && m){
          /* index of the kept points by packed key (built in parallel), then a serial flood fill that
           * compares 64-bit keys only */
          unsigned long long *KY = (unsigned long long*)malloc(m*8);
          #pragma omp parallel for schedule(static)
          for(long i = 0; i < m; i++) KY[i] = KEY(keep->p[i].k);
          long hc = 1 << 20; while(hc < 2*m) hc <<= 1; long hm = hc - 1, nidx = 0;
          unsigned long long *HK2 = (unsigned long long*)calloc(hc, 8); long *HV = (long*)malloc(hc*sizeof(long));
          #pragma omp parallel for schedule(static)
          for(long i = 0; i < m; i++) HV[hslot(HK2, hm, KY[i], &nidx)] = i;
          unsigned char *tr = (unsigned char*)calloc(m, 1); long *stk = (long*)malloc(m*sizeof(long));
          for(int t = 1; t <= 2; t++){
              double lim = (t == 1 ? 10 : 100)*thr; long ns = 0;
              for(long i = 0; i < m; i++){ int z = 1; for(int d = 0; d < DIM; d++) if(keep->p[i].k[d]){ z = 0; break; }
                  if(z && kv[i] >= lim){ tr[i] = t; stk[ns++] = i; } }
              while(ns){
                  long i = stk[--ns];
                  for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
                      int kk = keep->p[i].k[d] + st*sg;
                      if(kk > K || kk < -K) continue;
                      unsigned long long y = sg > 0 ? KY[i] + P[d] : KY[i] - P[d];
                      long q = (long)((y*0x9E3779B97F4A7C15ULL) >> 21) & hm;
                      while(HK2[q]){
                          if(HK2[q] == y){ long j = HV[q]; if(tr[j] < t && kv[j] >= lim){ tr[j] = t; stk[ns++] = j; } break; }
                          q = (q + 1) & hm; }
                  }
              }
          }
          free(KY); free(HK2); free(HV);
          #pragma omp parallel for schedule(static)
          for(long i = 0; i < m; i++) keep->p[i].k[PFLAG] = (keep->p[i].k[PFLAG] & 1) | (tr[i] << 1);
          if(PRED > 0){                         /* what the eps*10^t grow evaluated: z = 0 and the neighbours of its expanded points */
              int k0[MAXD] = {0}; vraise(k0, 2);
              #pragma omp parallel for schedule(dynamic,4096)
              for(long i = 0; i < m; i++) if(tr[i]) for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
                  int kk[MAXD]; memcpy(kk, keep->p[i].k, sizeof(int)*DIM); kk[d] += st*sg;
                  if(kk[d] <= K && kk[d] >= -K) vraise(kk, tr[i]); }
          }
          free(tr); free(stk);
      }
      keep->n = m; }
    #undef KEY
    #undef SEEN
    free(kv); free(fr.p); free(nx.p); free(HS);
    *S = s; *A = a; *e1 = x1; *e2 = x2; *drop = dr; *nev += ne; *nq += q; *un += u;
    return 1;
}

/* an integrand that is exactly zero everywhere (e.g. a VP family that vanishes by construction, v3_p4Vc)
 * gives grow() a zero threshold, so every point counts as above it and the whole box is evaluated.
 * Probe z = 0 and ZPROBE fixed pseudo-random nodes of the step-2^lev lattice first: if every value is
 * exactly 0 the graph is reported as zero at once.  A nonzero integrand that vanishes exactly at all of
 * them is not a practical case. */
#define ZPROBE 4096
static int box_zero(int lev, long *nev){
    int st = 1 << lev, W = 2*(K/st) + 1, nz = 0;
    #pragma omp parallel for reduction(+:nz) schedule(static)
    for(long i = 0; i <= ZPROBE; i++){
        int k[MAXD] = {0}; unsigned long long y = 0x9E3779B97F4A7C15ULL*(unsigned long long)(i + 1);
        if(i) for(int d = 0; d < DIM; d++){ y ^= y >> 29; y *= 0xBF58476D1CE4E5B9ULL; y ^= y >> 32; k[d] = ((int)(y % (unsigned long long)W) - K/st)*st; }
        long q = 0; double u = 0; int nr;
        if(leafx(k, 0, &q, &u, &nr) != 0) nz++;
    }
    if(nz) return 0;
    *nev += ZPROBE + 1;                 /* counted only when they replace the scan */
    return 1;
}

/* a block of the refinement: the children of cells c0..c1-1 of the current level, packed in cell order */
typedef struct { long c0, c1, n; int (*k)[MAXD]; unsigned char *sk, *fl, *te, *nr; double *v; } blk_t;
typedef struct { int (*k)[MAXD]; unsigned char *sk, *fl, *te; long *cn, *of; } scratch_t;   /* 2^D slots per cell */
static void prep(blk_t *bk, scratch_t *S, const vec_t *cur, long c0, long B, int last, int st_, double thr){
    long c1 = c0 + B < cur->n ? c0 + B : cur->n, nsk = 0, n = 0;
    int (*pk)[MAXD] = S->k; unsigned char *psk = S->sk, *pfl = S->fl, *pte = S->te; long *pcn = S->cn, *pof = S->of;
    /* children of each cell go to its own 2^D slots in parallel (the predictor's corner lookups are
     * most of this work), then are packed in cell order: same points, same order as a serial pass */
    #pragma omp parallel reduction(+:nsk)
    { double *st = (double*)malloc(sizeof(double) << DIM), *st1 = (double*)malloc(sizeof(double) << DIM), *st2 = (double*)malloc(sizeof(double) << DIM);
      unsigned long long *hy = (unsigned long long*)malloc(sizeof(unsigned long long) << DIM);
      #pragma omp for schedule(dynamic,4)
      for(long c = c0; c < c1; c++){
      const int *kp = cur->p[c].k; long base = (c - c0) << DIM, m = base;
              /* prediction: child k - st b is the centre of the 2h cube with corners k - 2st c, c a subset of b;
               * st[b] = mean |f| over those corners (a corner never evaluated counts as thr), by a subset-sum transform */
              int pon = last && PRED > 0 && VN > 0, tp = DIM < MAXD ? kp[PFLAG] >> 1 : 0;
              if(pon){
                  int kn[MAXD]; double lv;
                  for(int cc = 0; cc < (1 << DIM); cc++){          /* hashes first, then the lookups with the slots prefetched ahead */
                      for(int i = 0; i < DIM; i++) kn[i] = kp[i] - (((cc >> i) & 1) ? 2*st_ : 0);
                      hy[cc] = vhash(kn); }
                  for(int cc = 0; cc < 16 && cc < (1 << DIM); cc++) __builtin_prefetch(&VS[vpos(hy[cc])]);
                  for(int cc = 0; cc < (1 << DIM); cc++){
                      if(cc + 16 < (1 << DIM)) __builtin_prefetch(&VS[vpos(hy[cc + 16])]);
                      int f = vget_h(hy[cc], &lv);
                      st[cc] = f ? exp(lv) : thr;
                      if(ESOK){ st1[cc] = f >= 2 ? st[cc] : 10*thr; st2[cc] = f >= 3 ? st[cc] : 100*thr; }   /* that run's unseen corners count as its thr */
                  }
                  for(int i = 0; i < DIM; i++) for(int cc = 0; cc < (1 << DIM); cc++) if((cc >> i) & 1) st[cc] += st[cc ^ (1 << i)];
                  for(int cc = 0; cc < (1 << DIM); cc++) st[cc] /= (double)(1 << __builtin_popcount(cc));
                  if(ESOK){
                      for(int i = 0; i < DIM; i++) for(int cc = 0; cc < (1 << DIM); cc++) if((cc >> i) & 1){ st1[cc] += st1[cc ^ (1 << i)]; st2[cc] += st2[cc ^ (1 << i)]; }
                      for(int cc = 0; cc < (1 << DIM); cc++){ double dv = (double)(1 << __builtin_popcount(cc)); st1[cc] /= dv; st2[cc] /= dv; }
                  }
              }
              for(int b = 0; b < (1 << DIM); b++){
                int ok = 1, sk = 0;
                for(int i = 0; i < DIM; i++){ pk[m][i] = kp[i] - (((b >> i) & 1) ? st_ : 0); if(pk[m][i] < -K) ok = 0; }
                if(!ok) continue;
                if(pon && st[b] < PRED*thr){
                    nsk++;
                    if(PDIAG) sk = 1;                                   /* evaluated and summed anyway; reported */
                    else if(DUMPF && sampled(pk[m])) sk = 2;            /* evaluated for the dump only, not summed */
                    else continue;
                }
                psk[m] = ABSKIP && (kp[PFLAG] & 1); pfl[m] = (unsigned char)sk;
                { int te = tp; if(pon && ESOK){ int tq = (st1[b] >= PRED*10*thr) + (st1[b] >= PRED*10*thr && st2[b] >= PRED*100*thr); if(tq < te) te = tq; }
                  pte[m] = (unsigned char)te; } m++;
              }
              pcn[c - c0] = m - base;
      }
      free(st); free(st1); free(st2); free(hy); }
    PSKIP += nsk;
    for(long j = 0; j < c1 - c0; j++){ pof[j] = n; n += pcn[j]; }
    #pragma omp parallel for schedule(dynamic,16)
    for(long j = 0; j < c1 - c0; j++){
        long base = j << DIM, m = pcn[j], o = pof[j];
        memcpy(bk->k[o], pk[base], m*sizeof *pk); memcpy(bk->sk + o, psk + base, m); memcpy(bk->fl + o, pfl + base, m); memcpy(bk->te + o, pte + base, m);
    }
    bk->c0 = c0; bk->c1 = c1; bk->n = n;
}

#ifdef SQ_GPU
/* The last level on the plain grid (no predictor, dump or histogram) with the GPU: a cell's children are
 * just its 2^D corners k - st b, so only the cells go to the GPU, which makes the children itself
 * (gpu_submit_cells), instead of 2^D 48-byte nodes per cell made, packed and copied on the host.  The
 * children come back in prep()'s order (b ascending; those below the box are skipped by the sum), so the
 * sums are the same.  prepcells() fills bk with the cells c0..: coordinates, skip flag, shadow tier. */
void eval_gpu_submit_cells(int slot, long nc, const int (*cells)[MAXD], const unsigned char *skip, int st);
void eval_gpu_finish_cells(int slot, long nc, const int (*cells)[MAXD], const unsigned char *skip, int st, double *out, long *nq, double *unres);
static void prepcells(blk_t *bk, const vec_t *cur, long c0, long B){
    long c1 = c0 + B < cur->n ? c0 + B : cur->n;
    #pragma omp parallel for schedule(static)
    for(long j = 0; j < c1 - c0; j++){
        const int *kp = cur->p[c0 + j].k; memcpy(bk->k[j], kp, sizeof bk->k[j]);
        bk->sk[j] = ABSKIP && (kp[PFLAG] & 1); bk->te[j] = (unsigned char)(DIM < MAXD ? kp[PFLAG] >> 1 : 0); bk->fl[j] = 0;
    }
    bk->c0 = c0; bk->c1 = c1; bk->n = c1 - c0;
}
#endif

res_t hsum(const graph_t *g, double h, int N, double eps, int top){
    ckpt_t C; memset(&C, 0, sizeof C);
    if(HON < 0){ const char *e = getenv("SINCQED_HIST"); HON = e && *e; }
    memset(HC, 0, sizeof HC); memset(HS_, 0, sizeof HS_); memset(HA, 0, sizeof HA); memset(RC, 0, sizeof RC); memset(RS, 0, sizeof RS); memset(RA, 0, sizeof RA);
    if(PRED < 0){ const char *e = getenv("SINCQED_PREDICT"); PRED = e ? atof(e) : 0; }
    if(VCAP){ memset(VS, 0, VCAP*sizeof(vslot_t)); VN = 0; } PSKIP = 0; PDS = PDA = PDN = 0;
    if(PDIAG < 0){ PDIAG = getenv("SINCQED_PREDDIAG") != NULL; ST = (double*)malloc(sizeof(double) << MAXD); }
    { const char *e = getenv("SINCQED_DUMP"), *d = getenv("SINCQED_DUMPSAMPLE"); if(d && atol(d) > 0) DUMPS = atol(d);
      if(e && *e){ char fn[4096]; snprintf(fn, sizeof fn, "%s_%s.pts", e, g->name); DUMPF = fopen(fn, "wb");
        if(DUMPF) fprintf(DUMPF, "SQPT1 %s %d %.17g %d %.17g %d %d\n", g->name, g->D, h, N, eps, top, K); else say("# could not open %s\n", fn); } }
    double t0 = wtime(), tck = t0;
    memset(TQ, 0, sizeof TQ); TGROW = 0;
    ABSERR = ABSERR0;
    int resumed = ckpt_load(g, h, N, eps, top, &C);
    ES1 = ES2 = 0; ESOK = !resumed && DIM < MAXD;
    if(resumed){
        say("# resuming %s from checkpoint: refine to step %dh, cell %ld of %ld\n", g->name, 1 << C.lev, C.next, C.cur.n);
        t0 -= C.R.t;                            /* R.t holds the time already spent */
    } else {
        res_t *R = &C.R; double S, A, e1, e2, dr;
        if(box_zero(top, &R->nev)){
            say("# %s: integrand is exactly zero at %d probe nodes; reported as zero without a scan\n", g->name, ZPROBE + 1);
            R->t = wtime() - t0; if(DUMPF){ fclose(DUMPF); DUMPF = NULL; } return *R; }
        double tg0 = wtime();
        if(!(GROW && grow(top, eps*1e-3, &C.cur, &S, &A, &e1, &e2, &dr, &R->nev, &R->nq, &R->unres))){
            /* reference scale from a coarser scan, then the full box at step 2^top h */
            if(GROW) say("# %s: lattice too wide for 64-bit grid keys (D=%d, N=%d, step %dh): full box scan instead of grow\n", g->name, DIM, 2*K+1, 1 << top);
            scan(top+1, 0, NULL, &S, &A, &e1, &e2, &dr, &R->nev, &R->nq, &R->unres);
            ESOK = 0;                           /* the box scan does not set tiers */
            R->ref = A;
            scan(top, eps*R->ref*1e-3, &C.cur, &S, &A, &e1, &e2, &dr, &R->nev, &R->nq, &R->unres);
        } else R->ref = A;
        R->prune = dr;
        if(A == 0){ free(C.cur.p); R->t = wtime() - t0; if(DUMPF){ fclose(DUMPF); DUMPF = NULL; } return *R; }     /* integrand identically zero (e.g. a vanishing VP family) */
        { double r = e2 > 0 ? e1/e2 : 0; R->trunc = (r > 0 && r < 1) ? e1*r/(1-r) : e1; }
        C.nl = 0; C.Slev[C.nl++] = S;
        double tg1 = wtime();
        dilate(&C.cur, 1 << top);                /* add the axis neighbours */
        if(getenv("SINCQED_TIMING")) say("# %s timing: top level %.1fs (grow loop %.1fs), dilate+sort %.1fs\n", g->name, tg1 - tg0, TGROW, wtime() - tg1);
        C.R.kept_top = C.cur.n;
        C.lev = top - 1; C.next = 0;
    }
    res_t *R = &C.R;
    const long MAXPTS = USEGPU ? 1L << 22 : 1L << 18;   /* nodes per batch (the GPU wants large ones; the sums do not depend on it) */
    scratch_t SC; blk_t BL[2];
    SC.k = (int(*)[MAXD])malloc(MAXPTS*sizeof *SC.k); SC.sk = (unsigned char*)malloc(MAXPTS); SC.fl = (unsigned char*)malloc(MAXPTS); SC.te = (unsigned char*)malloc(MAXPTS);
    SC.cn = (long*)malloc(MAXPTS*sizeof(long)); SC.of = (long*)malloc(MAXPTS*sizeof(long));
    for(int j = 0; j < 2; j++){ BL[j].k = (int(*)[MAXD])malloc(MAXPTS*sizeof *BL[j].k); BL[j].sk = (unsigned char*)malloc(MAXPTS); BL[j].fl = (unsigned char*)malloc(MAXPTS);
        BL[j].te = (unsigned char*)malloc(MAXPTS); BL[j].nr = (unsigned char*)malloc(MAXPTS); BL[j].v = (double*)malloc(MAXPTS*sizeof(double)); BL[j].n = 0; }
    if(DIM >= MAXD) ABSKIP = 0;
    /* a point's roundoff only needs to be small next to what the pruning drops anyway (cells below eps*Ref*1e-2) */
    ABSERR = fmax(ABSERR0, ABREL*eps*R->ref);
    for(; C.lev >= 0; C.lev--, C.next = 0){
        int lev = C.lev, st = 1 << lev, st_ = st, last = (lev == 0); double pw = pow(st, DIM);
        double thr = eps*R->ref*1e-2;
        long B = MAXPTS >> DIM; if(B < 1) B = 1;
        snprintf(STAGE, sizeof STAGE, "refine to step %dh", st);
        if(C.next == 0){ C.s = C.a = C.drop = C.u = 0; C.ne = C.q = 0; C.nxt.n = 0; }
        /* blocks of B cells, two in flight with the GPU: block j+1 is prepared and submitted before block j
         * is finished and summed, so that host work overlaps the GPU's evaluation of the other block.
         * Everything on the host stays in this thread, in cell order: same sums as a serial pass. */
        int sl = 0;
#ifdef SQ_GPU
        static int CELLS = -1; if(CELLS < 0){ const char *e = getenv("SINCQED_GPUCELLS"); CELLS = !(e && !strcmp(e, "0")); }
        int fastc = CELLS && USEGPU && last && PRED <= 0 && !DUMPF && HON <= 0;   /* children made on the GPU (prepcells); SINCQED_GPUCELLS=0: off */
        #define PREP(b_, c0_) do{ if(fastc) prepcells(b_, &C.cur, c0_, B); else prep(b_, &SC, &C.cur, c0_, B, last, st_, thr); }while(0)
        #define SUBMIT(s_, b_) do{ if(fastc) eval_gpu_submit_cells(s_, (b_)->n, (const int(*)[MAXD])(b_)->k, (b_)->sk, st_); \
                                   else eval_gpu_submit(s_, (b_)->n, (const int(*)[MAXD])(b_)->k, (b_)->sk); }while(0)
#else
        int fastc = 0;
        #define PREP(b_, c0_) prep(b_, &SC, &C.cur, c0_, B, last, st_, thr)
        #define SUBMIT(s_, b_) ((void)0)
#endif
        if(C.next < C.cur.n){ double tq6 = wtime(); PREP(&BL[sl], C.next); TQ[0] += wtime() - tq6;
                              if(USEGPU) SUBMIT(sl, &BL[sl]); }
        while(C.next < C.cur.n){
            blk_t *bk = &BL[sl], *bn = &BL[sl ^ 1]; long n = bk->n; double tq2 = wtime();
            if(USEGPU && bk->c1 < C.cur.n){ PREP(bn, bk->c1); SUBMIT(sl ^ 1, bn); }
            double tq3 = wtime(); TQ[0] += tq3 - tq2;
            int (*pk)[MAXD] = bk->k; unsigned char *psk = bk->sk, *pfl = bk->fl, *pte = bk->te, *pnr = bk->nr; double *pv = bk->v;
#ifdef SQ_GPU
            if(fastc){                                        /* n cells; their children in prep() order */
                eval_gpu_finish_cells(sl, n, (const int(*)[MAXD])pk, psk, st_, pv, &C.q, &C.u);
                double tq4 = wtime(); long nv = 0;
                for(long c = 0; c < n; c++){                  /* fixed order: the sums of the loop below */
                    int bad = 0; for(int d = 0; d < DIM; d++) if(pk[c][d] - st_ < -K) bad |= 1 << d;
                    int ps = psk[c], te = pte[c]; const double *vc = pv + (c << DIM);
                    for(long b = 0; b < (1L << DIM); b++){
                        if(b & bad) continue;
                        double v = vc[b];
                        C.s += v*pw; C.a += fabs(v); ABSKIPPED += ps;
                        if(te >= 1) ES1 += v*pw; if(te >= 2) ES2 += v*pw;
                        nv++;
                    }
                }
                n = nv; TQ[2] += tq4 - tq3; tq3 = tq4;
            } else
#endif
            {
            if(USEGPU) eval_gpu_finish(sl, n, (const int(*)[MAXD])pk, psk, pnr, pv, &C.q, &C.u);
            else
            #pragma omp parallel
            eval_batch(n, (const int(*)[MAXD])pk, psk, pnr, pv, &C.q, &C.u);
            double tq4 = wtime();
            TQ[2] += tq4 - tq3; tq3 = tq4;
            for(long i = 0; i < n; i++){                      /* fixed order: reproducible sums */
                double v = pv[i];
                if(DUMPF && (!last || sampled(pk[i]))) dump(lev, pfl[i] != 0, pk[i], v*pw);
                if(pfl[i]){ PDS += v; PDA += fabs(v); PDN += fabs(v) >= thr; if(pfl[i] == 2) continue; }
                C.s += v*pw; C.a += fabs(v); hadd(lev, v*pw, pk[i]); ABSKIPPED += psk[i];
                if(last){ if(pte[i] >= 1) ES1 += v*pw; if(pte[i] >= 2) ES2 += v*pw; }
                if(PRED > 0 && lev == 1) vput(pk[i], fabs(v), pte[i]);
                if(!last){ if(fabs(v) >= thr){ push(&C.nxt, pk[i]); if(DIM < MAXD){ int tv = TIER(fabs(v), thr); C.nxt.p[C.nxt.n-1].k[PFLAG] = (ABSKIP && (psk[i] || !pnr[i])) | ((pte[i] < tv ? pte[i] : tv) << 1); } } else C.drop += fabs(v)*pw; }
            }
            }
            C.ne += n; C.next = bk->c1; { double tq5 = wtime(); TQ[3] += tq5 - tq3; }
            progress(STAGE, C.next, C.cur.n);
            if(CKPT_EVERY > 0 && wtime() - tck >= CKPT_EVERY && C.next < C.cur.n){
                R->t = wtime() - t0; ckpt_save(g, h, N, eps, top, &C); tck = wtime(); }
            sl ^= 1;
            if(!USEGPU && C.next < C.cur.n){ double tq6 = wtime(); prep(&BL[sl], &SC, &C.cur, C.next, B, last, st_, thr); TQ[0] += wtime() - tq6; }
        }
        #undef PREP
        #undef SUBMIT
        R->prune += C.drop;                     /* weight dropped below the threshold at this level counts too */
        R->nev += C.ne; R->nq += C.q; if(last) R->unres = sqrt(C.u); C.Slev[C.nl++] = C.s;   /* only S_h's own points count */
        if(last) R->Ah = C.a;
        free(C.cur.p); C.cur = C.nxt; memset(&C.nxt, 0, sizeof C.nxt);
        if(CKPT_EVERY > 0 && !last){ R->t = wtime() - t0; C.lev--; C.next = 0; ckpt_save(g, h, N, eps, top, &C); C.lev++; tck = wtime(); }
    }
    free(SC.k); free(SC.sk); free(SC.fl); free(SC.te); free(SC.cn); free(SC.of);
    for(int j = 0; j < 2; j++){ free(BL[j].k); free(BL[j].sk); free(BL[j].fl); free(BL[j].te); free(BL[j].nr); free(BL[j].v); } free(C.cur.p); free(C.nxt.p);
    R->S[0] = C.Slev[C.nl-1]; R->S[1] = C.nl > 1 ? C.Slev[C.nl-2] : NAN; R->S[2] = C.nl > 2 ? C.Slev[C.nl-3] : NAN;
    R->t = wtime() - t0;
    ckpt_remove(g);
    if(getenv("SINCQED_TIMING")) say("# %s timing (refinement, host thread): prepare %.1fs  wait+check %.1fs  sum %.1fs\n", g->name, TQ[0], TQ[2], TQ[3]);
    if(PRED > 0) say("# %s predictor alpha %g: skipped %ld final-level points\n", g->name, PRED, PSKIP);
    if(ESOK){ double w = g->weight, d1 = w*(R->S[0] - ES1), d2 = w*(ES1 - ES2), r = d1 != 0 ? d2/d1 : 0;   /* weighted, like the results */
        say("# %s eps shadow: eps*10 would change the value by %+.3e, eps*100 by a further %+.3e (ratio %.2f)", g->name, -d1, -d2, r);
        if(r > 1.5) say("; pruning loss left at eps ~ %+.2e (add to the value)\n", d1/(r-1)); else say("\n"); }
    if(PRED > 0 && PDIAG) say("# %s predictor check: skipped points sum %+.3e, sum|.| %.3e, %.0f of them above thr\n", g->name, PDS, PDA, PDN);
    hwrite(g, top);
    if(DUMPF){ fclose(DUMPF); DUMPF = NULL; }
    return *R;
}
