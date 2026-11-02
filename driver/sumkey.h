/* sumkey.h -- the key-based parts of the hierarchical sum (grid keys, the dilation and the support-
 * following grow), included twice by sum.c: once with 64-bit keys (KT = unsigned long long, the
 * original code) and once with 128-bit keys for lattices with W^D > 2^63 (top = 0 at five loops, D = 13).
 * The including file defines KT (key type), KS (hash slot type), KTOP (largest allowed W^D), KBITLEN(p)
 * (bit length), KHASH(y), KEMPTY(H, i), KGET(H, i), KEQ(H, i, y), SFX(name), and the hash routines
 * SFX(hslot) / SFX(hfind). */
/* A point's packed key is the mixed-radix number 1 c_0 c_1 ... c_{D-1} in base W = 2*(K/st) + 1, with
 * c_i = k_i/st + K/st (the leading 1 keeps 0 free as the empty slot).  An axis step on dimension j adds
 * +-P[j] = +-W^(D-1-j).  Numeric order is lexicographic order in k, as with the earlier power-of-two
 * fields, but the key fits in 64 bits while W^D <= 2^63 rather than ceil(log2 W)*D <= 63: at top step 2h
 * D = 10 now reaches N = 153 (was 121) and D = 11 N = 97 (was 57).  keyradix() fills P and
 * the key's bit length and returns W, or 0 if the keys would not fit. */
static int SFX(keyradix)(int st, KT *P, int *nb){
    int W = 2*(K/st) + 1; KT p = 1;
    for(int i = DIM - 1; i >= 0; i--){ P[i] = p; if(p > KTOP/W) return 0; p *= W; }
    *nb = KBITLEN(p) + 1;                      /* keys < 2 W^D */
    return W;
}

static int SFX(dilate)(vec_t *cur, int st){
    KT P[MAXD]; int nb, W = SFX(keyradix)(st, P, &nb);
    long n0 = cur->n;
    if(!W) return 0;
    if(n0 == 0 || getenv("SINCQED_OLDDILATE")){ dilate_old(cur, st); return 1; }
    long cap = 1 << 20; while(cap < 4*n0) cap <<= 1;
    KS *HK; unsigned char *SK, *TR; long cnt;
    for(;;){
        HK = (KS*)calloc(cap, sizeof(KS)); SK = (unsigned char*)malloc(cap); TR = (unsigned char*)calloc(cap, 1);
        memset(SK, 1, cap); cnt = 0; int full = 0; long lim = cap/4*3, mask = cap - 1;
        #pragma omp parallel for schedule(dynamic,4096)
        for(long q = 0; q < n0; q++){
            if(__atomic_load_n(&full, __ATOMIC_RELAXED)) continue;
            int k0[MAXD]; ptget(&cur->p[q], k0); int fl = DIM < MAXD ? k0[PFLAG] : 0;
            unsigned char sk = (unsigned char)(fl & 1), tr = (unsigned char)(fl >> 1);
            KT y0 = 1; for(int i = 0; i < DIM; i++) y0 = y0*W + (KT)(k0[i]/st + K/st);
            for(int j = -1; j < DIM; j++) for(int sg = -1; sg <= 1; sg += 2){
                if(j < 0 && sg > 0) continue;               /* j = -1: the cell itself */
                KT y = y0;
                if(j >= 0){ int kk = k0[j] + st*sg; if(kk > K || kk < -K) continue;
                    y = sg > 0 ? y + P[j] : y - P[j]; }
                long i = SFX(hslot)(HK, mask, y, &cnt);
                if(!sk) __atomic_store_n(&SK[i], 0, __ATOMIC_RELAXED);
                unsigned char o = __atomic_load_n(&TR[i], __ATOMIC_RELAXED);
                while(o < tr && !__atomic_compare_exchange_n(&TR[i], &o, tr, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) ;
                if(__atomic_load_n(&cnt, __ATOMIC_RELAXED) > lim){ __atomic_store_n(&full, 1, __ATOMIC_RELAXED); break; }
            }
        }
        if(!full) break;
        free(HK); free(SK); free(TR); cap *= 2;      /* rare: more distinct points than 3/4 of the table */
    }
    /* distinct keys, radix-sorted (LSD, 16-bit digits); each thread takes a fixed slice, so the result
     * is the serial one */
    long m = 0, mask = cap - 1; int nt = omp_get_max_threads();
    KT *U = (KT*)malloc((cnt ? cnt : 1)*sizeof(KT)), *T = (KT*)malloc((cnt ? cnt : 1)*sizeof(KT));
    long *cn = (long*)malloc((size_t)nt*(1L << 16)*sizeof(long)), *tc = (long*)calloc(nt + 1, sizeof(long));
    #pragma omp parallel num_threads(nt)
    { int t = omp_get_thread_num(); long n = 0, *c = cn + ((long)t << 16);
      #pragma omp for schedule(static)
      for(long i = 0; i < cap; i++) if(!KEMPTY(HK, i)) n++;
      tc[t + 1] = n;
      #pragma omp barrier
      #pragma omp single
      for(int j = 0; j < nt; j++) tc[j + 1] += tc[j];
      long o = tc[t];
      #pragma omp for schedule(static)
      for(long i = 0; i < cap; i++) if(!KEMPTY(HK, i)) U[o++] = KGET(HK, i);
      long mm = tc[nt];
      for(int sh = 0; sh < nb; sh += 16){
          memset(c, 0, (1L << 16)*sizeof(long));
          #pragma omp for schedule(static)
          for(long i = 0; i < mm; i++) c[(U[i] >> sh) & 0xffff]++;
          #pragma omp single
          { long sum = 0; for(long d = 0; d < (1L << 16); d++) for(int j = 0; j < nt; j++){ long *cj = cn + ((long)j << 16), x = cj[d]; cj[d] = sum; sum += x; } }
          #pragma omp for schedule(static)
          for(long i = 0; i < mm; i++) T[c[(U[i] >> sh) & 0xffff]++] = U[i];
          #pragma omp single
          { KT *x = U; U = T; T = x; }
      }
    }
    m = tc[nt]; free(cn); free(tc); free(T);
    free(cur->p); cur->p = (pt_t*)malloc((m ? m : 1)*sizeof(pt_t)); cur->n = cur->cap = m;
    #pragma omp parallel for schedule(static)
    for(long j = 0; j < m; j++){
        KT y = U[j]; pt_t *p = &cur->p[j]; memset(p, 0, sizeof *p);
        for(int i = DIM - 1; i >= 0; i--){ p->k[i] = ((int)(y % W) - K/st)*st; y /= W; }
        if(DIM < MAXD){ long s = SFX(hfind)(HK, mask, U[j]); p->k[PFLAG] = SK[s] | (TR[s] << 1); }
    }
    free(U); free(HK); free(SK); free(TR);
    return 1;
}

/* every key of O (oc slots) into the empty table H (mask): blocks of 64 slots, their keys' new slots
 * prefetched before the inserts.  Returns the number of keys */
static long SFX(rehash)(KS *H, long mask, const KS *O, long oc){
    long n = 0;
    #pragma omp parallel for schedule(static) reduction(+:n)
    for(long b = 0; b < (oc + 63)/64; b++){
        KT yy[64]; int m = 0;
        for(long j = 64*b; j < oc && j < 64*b + 64; j++) if(!KEMPTY(O, j)){ yy[m] = KGET(O, j); __builtin_prefetch(&H[KHASH(yy[m]) & mask], 1); m++; }
        for(int t = 0; t < m; t++) SFX(hslot)(H, mask, yy[t], &n);
    }
    return n;
}
/* the same into the grow's interleaved table H (slot i at H[2i], its claim word at H[2i+1]), every claim
 * set to UINT32_MAX (seen) */
static long SFX(rehash2)(KS *H, long mask, const KS *O, long oc){
    long n = 0;
    #pragma omp parallel for schedule(static) reduction(+:n)
    for(long b = 0; b < (oc + 63)/64; b++){
        KT yy[64]; int m = 0;
        for(long j = 64*b; j < oc && j < 64*b + 64; j++) if(!KEMPTY(O, 2*j)){ yy[m] = KGET(O, 2*j); __builtin_prefetch(&H[2*(KHASH(yy[m]) & mask)], 1); m++; }
        for(int t = 0; t < m; t++) *(unsigned*)&H[2*SFX(hslots)(H, mask, yy[t], &n, 1) + 1] = UINT32_MAX;
    }
    return n;
}
/* is y in the table (unlike hfind, y may be absent) */
static inline int SFX(hhas)(const KS *H, long mask, KT y){
    long q = KHASH(y) & mask;
    while(!KEMPTY(H, q)){ if(KEQ(H, q, y)) return 1; q = (q + 1) & mask; }
    return 0;
}

/* support-following scan at level lev (step 2^lev): start at z = 0 and grow through axis neighbours
 * while |f| is above thr = eps_rel * (running sum|f|); points below it are evaluated (the boundary)
 * but not expanded.  Cost ~ support size x 2D, not (box width)^D.  Same outputs as scan(); *A is
 * also the reference scale.  Returns 0 if the lattice is too large to index with these keys (keyradix). */
static int SFX(grow)(int lev, double eps_rel, vec_t *keep, double *S, double *A, double *e1, double *e2, double *drop, long *nev, long *nq, double *un){
    int st = 1 << lev, nb; KT P[MAXD]; int W = SFX(keyradix)(st, P, &nb);
    if(!W) return 0;
    /* Which evaluated points must a new candidate be checked against?  The lattice is bipartite (a step
     * flips the parity of sum k/st) and frontier L holds the points first reached in L steps, so the
     * neighbours of a kept point of frontier L lie in frontiers L-1 and L+1.  A kept point of an earlier
     * frontier j would have put its neighbours into frontier j+1, so the only kept points a candidate
     * can repeat are those of frontier L-1.  Points evaluated but not kept (the boundary) are not expanded
     * and can be met again from a later frontier, but only from a nearby one: a step changes the norm
     * n = sum |k|/st by 1, so a point of frontier L has n <= L, and a candidate y can only repeat a boundary
     * point of a frontier L' with n(y) <= L' <= L-1 (same parity).  Nearly every kept point has n = L (the
     * grow from z = 0 is close to monotone), so its outward steps (n = L+1) need no boundary lookup at all and
     * its inward ones look only at frontier L-1's boundary.  So the boundary is one hash table per frontier
     * (BT[L], sized once, never rehashed), and HS is rebuilt each frontier from the previous frontier's kept
     * points (pk) plus this frontier's candidates.  HS interleaves each key with its claim word OW(slot) (one cache
     * line for both): UINT32_MAX for a point already seen, else the claim of the first
     * candidate (lowest id) that reached it.  Same frontiers, in the same order, as one table of every point. */
    long cap = 1 << 20; KS *HS = (KS*)aligned_alloc(64, 2*cap*sizeof(KS)); long nh = 0;
    memset(HS, 0, 2*cap*sizeof(KS));
    #define OW(q) (*(unsigned*)&HS[2*(q) + 1])
    KS **BT = NULL; long *BM = NULL; int nbt = 0;     /* boundary of frontier L: table BT[L], mask BM[L] (NULL if none) */
    vec_t pk = {0};
    int store = lev > 0;     /* top = 0: nothing refines the kept points, so only their |f| is kept (final threshold) */
    #define KEY(k) ({ KT _y = 1; for(int _i = 0; _i < DIM; _i++) _y = _y*W + (KT)((k)[_i]/st + K/st); _y; })
    vec_t fr = {0}, nx = {0};
    double pw = pow(st, DIM), a = 0, x1 = 0, x2 = 0, dr = 0, u = 0; long ne = 0, q = 0; ksum_t ks = {0, 0};
    GSHS = 0; GSHN = 0;
    int k0[MAXD] = {0}; push(&fr, k0);
    double *kv = NULL; long nkv = 0;   /* |f| of the kept points, for the final threshold */
    OW(SFX(hslots)(HS, cap - 1, KEY(k0), &nh, 1)) = UINT32_MAX;
    long *kj = NULL, nkj = 0; unsigned char *own = NULL; long nown = 0, *qs = NULL;      /* kept frontier points; per-candidate owner flags */
    DONE = 0; snprintf(STAGE, sizeof STAGE, "grow at step %dh", st);
    double tw0 = wtime(); TGE = TGS = 0; TGB[0] = TGB[1] = TGB[2] = 0; int lay = 0;
    while(fr.n){ double te0 = wtime(); double sk0 = kval(&ks);
        double *v = (double*)malloc(fr.n*sizeof(double)); unsigned char *nr = (unsigned char*)malloc(fr.n);
        /* gs > 0 (top = 0 only, where no refinement follows): a frontier point skips the |.| sums when the
         * point that reached it was checked and not near the fallback, or skipped fewer than gs generations
         * after such a point (the generation count is in PFLAG); 0 = every point checked, as before */
        int gs = lev == 0 && ABSKIP && DIM < MAXD && !USEGPU ? GROWSKIP : 0;
        unsigned char *sk = (unsigned char*)calloc(fr.n, 1);
        if(USEGPU){ int (*fk)[MAXD] = (int(*)[MAXD])malloc(fr.n*sizeof *fk); for(long i = 0; i < fr.n; i++) ptget(&fr.p[i], fk[i]);
            eval_batch_gpu(fr.n, (const int(*)[MAXD])fk, NULL, nr, v, &q, &u); free(fk); DONE += fr.n; grow_progress(lay, DONE, tw0, keep->n, fr.n, fr.n, sk0); }
        else if(GROWBATCH){      /* batches through eval_batch: the kernel's packed SIMD form, same values */
            long ns = 0, qq = 0; double uu = 0; const long CE = 1L << 20;
            int (*fk)[MAXD] = (int(*)[MAXD])malloc((fr.n < CE ? fr.n : CE)*sizeof *fk);   /* one chunk's coordinates at a time */
            for(long c0 = 0; c0 < fr.n; c0 += CE){
                long m = fr.n - c0 < CE ? fr.n - c0 : CE;
                #pragma omp parallel for schedule(static) reduction(+:ns)
                for(long i = 0; i < m; i++){ ptget(&fr.p[c0 + i], fk[i]); sk[c0 + i] = gs && fr.p[c0 + i].k[PFLAG] > 0; ns += sk[c0 + i]; }
                #pragma omp parallel
                eval_batch(m, (const int(*)[MAXD])fk, sk + c0, nr + c0, v + c0, &qq, &uu);
                DONE += m; grow_progress(lay, DONE, tw0, keep->n, c0 + m, fr.n, sk0);
            }
            ABSKIPPED += ns; q += qq; u += uu; free(fk);
        }
        else
        #pragma omp parallel for reduction(+:q,u) schedule(dynamic,64)
        for(long i = 0; i < fr.n; i++){
            int nri, k[MAXD]; ptget(&fr.p[i], k); v[i] = leafx(k, 0, &q, &u, &nri); nr[i] = (unsigned char)nri;
            { long dn;
              #pragma omp atomic capture
              dn = ++DONE;
              if(!(dn & 255)) grow_progress(lay, dn, tw0, keep->n, i, fr.n, sk0); }
        }
        ne += fr.n; double te1 = wtime(); TGE += te1 - te0;
        if(PRED > 0) vreserve(fr.n);
        if(PRED > 0 || HON > 0 || DUMPF){
            #pragma omp parallel for schedule(static) if(!HON && !DUMPF)
            for(long i = 0; i < fr.n; i++){ int k[MAXD]; ptget(&fr.p[i], k); if(PRED > 0) vput_par(k, fabs(v[i]), 0); hadd(lev, v[i]*pw, k); if(DUMPF) dump(lev, 0, k, v[i]*pw); }
        }
        double tb0 = wtime();
        for(long i = 0; i < fr.n; i++){ kadd(&ks, v[i]*pw); a += fabs(v[i])*pw; }
        double thr = eps_rel*a;                  /* per-point threshold from the running sum (as scan) */
        if(fr.n > nkj){ nkj = fr.n; kj = (long*)realloc(kj, nkj*sizeof(long)); }
        /* which points are kept, and the edge flags, in parallel; then the sums serially in point order (same
         * bits as one serial pass), and the kept points compacted in order */
        unsigned char *cl = (unsigned char*)malloc(fr.n);
        long nt = omp_get_max_threads(), *tc = (long*)calloc(nt + 1, sizeof(long));
        #pragma omp parallel num_threads(nt)
        { int tid = omp_get_thread_num(); long n = 0;
          #pragma omp for schedule(static)
          for(long i = 0; i < fr.n; i++){
              int mx = 0; for(int d = 0; d < DIM; d++){ int x = abs(fr.p[i].k[d]); if(x > mx) mx = x; }
              int kp = fabs(v[i]) >= thr; n += kp;
              cl[i] = (unsigned char)(kp | (mx == K ? 2 : mx == K - st ? 4 : 0));
          }
          tc[tid + 1] = n; }
        for(long i = 0; i < fr.n; i++){
            double w = fabs(v[i])*pw;
            if(cl[i] & 2) x1 += w; else if(cl[i] & 4) x2 += w;
            if(!(cl[i] & 1)){ dr += w; GSHS += v[i]*pw; GSHN++; }
        }
        for(int t = 0; t < nt; t++) tc[t + 1] += tc[t];
        long nk = tc[nt], kb = keep->n;
        if(kb + nk > nkv){ nkv = 2*(kb + nk) + 1024; kv = (double*)realloc(kv, nkv*sizeof(double)); }
        if(store && kb + nk > keep->cap){ keep->cap = 2*(kb + nk) + 1024; keep->p = (pt_t*)realloc(keep->p, keep->cap*sizeof(pt_t)); }
        #pragma omp parallel num_threads(nt)
        { long o = tc[omp_get_thread_num()];
          #pragma omp for schedule(static)
          for(long i = 0; i < fr.n; i++) if(cl[i] & 1){
              kv[kb + o] = fabs(v[i]); kj[o] = i;
              if(store){ pt_t *p = &keep->p[kb + o]; memset(p, 0, sizeof *p); for(int d = 0; d < DIM; d++) p->k[d] = fr.p[i].k[d];
                  if(ABSKIP && DIM < MAXD) p->k[PFLAG] = !nr[i]; }
              o++;
          } }
        keep->n = kb + nk; free(cl); free(tc);
        double tb1 = wtime(); TGB[0] += tb1 - tb0;
        { long nbd = fr.n - nk;
          if(lay >= nbt){ nbt = 2*lay + 64; BT = (KS**)realloc(BT, nbt*sizeof(KS*)); BM = (long*)realloc(BM, nbt*sizeof(long)); }
          BT[lay] = NULL; BM[lay] = 0;
          if(nbd > 0){
              long bc = 64; while(nbd > GROWLOAD*bc) bc *= 2;
              KS *B = (KS*)calloc(bc, sizeof(KS)); long nb_ = 0; BT[lay] = B; BM[lay] = bc - 1;
              #pragma omp parallel for schedule(static) reduction(+:nb_)
              for(long b = 0; b < (fr.n + 31)/32; b++){           /* blocks of 32: keys and prefetches first, then the inserts */
                  KT yy[32]; int m = 0;
                  for(long i = 32*b; i < fr.n && i < 32*b + 32; i++) if(!(fabs(v[i]) >= thr)){
                      int k[MAXD]; ptget(&fr.p[i], k); yy[m] = KEY(k); __builtin_prefetch(&B[KHASH(yy[m]) & (bc - 1)], 1); m++; }
                  for(int t = 0; t < m; t++) SFX(hslot)(B, bc - 1, yy[t], &nb_);
              }
          }
          long want = 1 << 20; while(pk.n + fr.n > GROWLOAD*want) want *= 2;   /* this frontier's size as a guess at the next */
          if(cap < want || cap > 4*want){ free(HS); cap = want; HS = (KS*)aligned_alloc(64, 2*cap*sizeof(KS)); }
          memset(HS, 0, 2*cap*sizeof(KS)); nh = 0;
          #pragma omp parallel for schedule(static) reduction(+:nh)
          for(long b = 0; b < (pk.n + 31)/32; b++){
              KT yy[32]; int m = 0;
              for(long i = 32*b; i < pk.n && i < 32*b + 32; i++){
                  int k[MAXD]; ptget(&pk.p[i], k); yy[m] = KEY(k); long hy = KHASH(yy[m]) & (cap - 1);
                  __builtin_prefetch(&HS[2*hy], 1); m++; }
              for(int t = 0; t < m; t++) OW(SFX(hslots)(HS, cap - 1, yy[t], &nh, 1)) = UINT32_MAX;
          }
        }
        /* the axis neighbours not seen before go to the next frontier, in the order of a serial pass over
         * (kept point, axis, sign): all candidates claim their slot in parallel (the lowest id wins), then
         * the winners are written in id order.  Chunks keep the ids below 2^32. */
        double tb2 = wtime(); TGB[1] += tb2 - tb1;
        nx.n = 0;
        /* chunks are sized so that even if every candidate were new the table stays below GROWLOAD full
         * (the frontier order does not depend on the chunking: a point claimed in an earlier chunk is skipped) */
        const long CH = (1L << 26)/(2*DIM);    /* also bounds qs (a slot index per candidate) to 512 MB */
        for(long c0 = 0, c1; c0 < nk; c0 = c1){
            long room;
            while((room = (long)(GROWLOAD*cap) - nh)/(2*DIM) < (nk - c0 < cap/(64*DIM) ? nk - c0 : cap/(64*DIM)) || room < 2*DIM){
                long oc = cap; KS *O = HS; cap *= 2; HS = (KS*)aligned_alloc(64, 2*cap*sizeof(KS));   /* rehash: every key in it is from an earlier chunk */
                memset(HS, 0, 2*cap*sizeof(KS)); nh = SFX(rehash2)(HS, cap - 1, O, oc); free(O); }
            c1 = c0 + room/(2*DIM); if(c1 > c0 + CH) c1 = c0 + CH; if(c1 > nk) c1 = nk;
            long nc = (c1 - c0)*2*DIM;
            if(nc > nown){ nown = nc; own = (unsigned char*)realloc(own, nown); qs = (long*)realloc(qs, nown*sizeof(long)); }
            long mask = cap - 1;
            /* claims, in blocks of QB kept points: first every candidate's key, with its table slots prefetched,
             * then the lookups (the hash probes are cache misses; this keeps many of them in flight) */
            enum { QB = 8 };
            long nbk = (c1 - c0 + QB - 1)/QB;
            #pragma omp parallel for schedule(static) reduction(+:nh)
            for(long bk = 0; bk < nbk; bk++){
                long j0 = c0 + bk*QB, j1 = j0 + QB < c1 ? j0 + QB : c1;
                KT yy[QB*2*MAXD]; int ny[QB*2*MAXD];
                for(long j = j0; j < j1; j++){
                    int k[MAXD], n = 0; ptget(&fr.p[kj[j]], k); KT y0 = KEY(k); unsigned char *ow = own + (j - c0)*2*DIM;
                    for(int d = 0; d < DIM; d++) n += abs(k[d]);
                    n /= st;
                    for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
                        int kk = k[d] + st*sg, c = 2*d + (sg > 0); ow[c] = 0; if(kk > K || kk < -K) continue;
                        KT y = sg > 0 ? y0 + P[d] : y0 - P[d]; long hy = KHASH(y);
                        int t = (j - j0)*2*DIM + c, n1 = abs(kk) > abs(k[d]) ? n + 1 : n - 1;
                        yy[t] = y; ny[t] = n1; ow[c] = 1;
                        if(n1 < lay && BT[n1]) __builtin_prefetch(&BT[n1][hy & BM[n1]]);
                        __builtin_prefetch(&HS[2*(hy & mask)], 1);
                    }
                }
                for(long j = j0; j < j1; j++){
                    unsigned char *ow = own + (j - c0)*2*DIM; long *qq = qs + (j - c0)*2*DIM;
                    for(int c = 0; c < 2*DIM; c++){
                        if(!ow[c]) continue;
                        int t = (j - j0)*2*DIM + c, hit = 0; KT y = yy[t];
                        for(int l = ny[t]; l < lay && !hit; l += 2) hit = BT[l] && SFX(hhas)(BT[l], BM[l], y);
                        if(hit){ ow[c] = 0; continue; }     /* a boundary point of an earlier frontier */
                        unsigned cl = UINT32_MAX - 1 - (unsigned)((j - c0)*2*DIM + c);
                        long q = SFX(hslots)(HS, mask, y, &nh, 1); qq[c] = q;
                        unsigned o = __atomic_load_n(&OW(q), __ATOMIC_RELAXED);
                        while(o < cl && !__atomic_compare_exchange_n(&OW(q), &o, cl, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) ;
                    }
                }
            }
            long nt = omp_get_max_threads(), *cnt = (long*)calloc(nt + 1, sizeof(long));
            #pragma omp parallel num_threads(nt)
            { int tid = omp_get_thread_num(); long n = 0;
              #pragma omp for schedule(static)
              for(long j = c0; j < c1; j++){
                  unsigned char *ow = own + (j - c0)*2*DIM; const long *qq = qs + (j - c0)*2*DIM;
                  for(int c = 0; c < 2*DIM; c++){
                      if(!ow[c]) continue;
                      unsigned cl = UINT32_MAX - 1 - (unsigned)((j - c0)*2*DIM + c);
                      if(OW(qq[c]) == cl) n++; else ow[c] = 0;
                  }
              }
              cnt[tid + 1] = n;
              #pragma omp barrier
              #pragma omp single
              { for(int t = 0; t < nt; t++) cnt[t + 1] += cnt[t];
                long tot = nx.n + cnt[nt]; if(tot > nx.cap){ nx.cap = tot + tot/2 + 1024; nx.p = (pt_t*)realloc(nx.p, nx.cap*sizeof(pt_t)); } }
              long o = nx.n + cnt[tid];
              #pragma omp for schedule(static)
              for(long j = c0; j < c1; j++){
                  int k[MAXD]; ptget(&fr.p[kj[j]], k); const unsigned char *ow = own + (j - c0)*2*DIM; const long *qq = qs + (j - c0)*2*DIM;
                  for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2) if(ow[2*d + (sg > 0)]){
                      pt_t *p = &nx.p[o++]; memset(p, 0, sizeof *p); for(int i = 0; i < DIM; i++) p->k[i] = (signed char)k[i];
                      p->k[d] = (signed char)(k[d] + st*sg);
                      if(gs){ long pi = kj[j]; int pg = fr.p[pi].k[PFLAG];
                          p->k[PFLAG] = (signed char)(!sk[pi] ? !nr[pi] : (pg < gs ? pg + 1 : 0)); }
                      OW(qq[2*d + (sg > 0)]) = UINT32_MAX;
                  }
              }
            }
            nx.n += cnt[nt]; free(cnt);
        }
        TGB[2] += wtime() - tb2;
        if(nk > pk.cap){ pk.cap = nk + nk/2 + 1024; pk.p = (pt_t*)realloc(pk.p, pk.cap*sizeof(pt_t)); }
        #pragma omp parallel for schedule(static)
        for(long j = 0; j < nk; j++) pk.p[j] = fr.p[kj[j]];
        pk.n = nk;
        free(v); free(nr); free(sk);
        vec_t t = fr; fr = nx; nx = t; TGS += wtime() - te1; lay++;
    }
    TGROW = wtime() - tw0; double tt0 = wtime();
    /* the running threshold was lower early on: re-apply the final one */
    { double thr = eps_rel*a; long m = 0;
      for(long i = 0; i < keep->n; i++){ if(kv[i] >= thr){ kv[m] = kv[i]; if(store) keep->p[m] = keep->p[i]; m++; } else dr += kv[i]*pw; }
      keep->n = m;
      /* shadow tiers: a run with eps*10^t keeps a point only if it reaches it from z = 0 through points with
       * |f| >= 10^t thr (the grow expands nothing else); flood-fill the kept set for t = 1, 2 */
      if(DIM < MAXD && m && store){       /* top = 0 has no refinement to use the tiers */
          /* index of the kept points by packed key (built in parallel), then a serial flood fill that
           * compares 64-bit keys only */
          KT *KY = (KT*)malloc(m*sizeof(KT));
          #pragma omp parallel for schedule(static)
          for(long i = 0; i < m; i++) KY[i] = KEY(keep->p[i].k);
          long hc = 1 << 20; while(hc < 2*m) hc <<= 1; long hm = hc - 1, nidx = 0;
          KS *HK2 = (KS*)calloc(hc, sizeof(KS)); long *HV = (long*)malloc(hc*sizeof(long));
          #pragma omp parallel for schedule(static) reduction(+:nidx)
          for(long i = 0; i < m; i++) HV[SFX(hslot)(HK2, hm, KY[i], &nidx)] = i;
          /* points reached from z = 0 through |f| >= lim, breadth first (the set does not depend on the order) */
          unsigned char *tr = (unsigned char*)calloc(m, 1); long *cu = (long*)malloc(m*sizeof(long)), *nw = (long*)malloc(m*sizeof(long));
          for(int t = 1; t <= 2; t++){
              double lim = (t == 1 ? 10 : 100)*thr; long nc = 0;
              for(long i = 0; i < m; i++){ int z = 1; for(int d = 0; d < DIM; d++) if(keep->p[i].k[d]){ z = 0; break; }
                  if(z && kv[i] >= lim){ tr[i] = t; cu[nc++] = i; } }
              while(nc){
                  long nn = 0;
                  #pragma omp parallel for schedule(dynamic,1024)
                  for(long c = 0; c < nc; c++){
                      long i = cu[c];
                      for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
                          int kk = keep->p[i].k[d] + st*sg;
                          if(kk > K || kk < -K) continue;
                          KT y = sg > 0 ? KY[i] + P[d] : KY[i] - P[d];
                          long q = KHASH(y) & hm;
                          while(!KEMPTY(HK2, q)){
                              if(KEQ(HK2, q, y)){ long j = HV[q];
                                  if(kv[j] >= lim){ unsigned char o = __atomic_load_n(&tr[j], __ATOMIC_RELAXED);
                                      while(o < t && !__atomic_compare_exchange_n(&tr[j], &o, (unsigned char)t, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) ;
                                      if(o < t) nw[__atomic_fetch_add(&nn, 1, __ATOMIC_RELAXED)] = j; }
                                  break; }
                              q = (q + 1) & hm; }
                      }
                  }
                  long *x = cu; cu = nw; nw = x; nc = nn;
              }
          }
          free(KY); free(HK2); free(HV);
          #pragma omp parallel for schedule(static)
          for(long i = 0; i < m; i++) keep->p[i].k[PFLAG] = (keep->p[i].k[PFLAG] & 1) | (tr[i] << 1);
          if(PRED > 0){                         /* what the eps*10^t grow evaluated: z = 0 and the neighbours of its expanded points */
              int k0[MAXD] = {0}; vraise(k0, 2);
              #pragma omp parallel for schedule(dynamic,4096)
              for(long i = 0; i < m; i++) if(tr[i]) for(int d = 0; d < DIM; d++) for(int sg = -1; sg <= 1; sg += 2){
                  int kk[MAXD]; ptget(&keep->p[i], kk); kk[d] += st*sg;
                  if(kk[d] <= K && kk[d] >= -K) vraise(kk, tr[i]); }
          }
          free(tr); free(cu); free(nw);
      }
      keep->n = m; }
    #undef KEY
    #undef SEEN
    free(kv); free(fr.p); free(nx.p); free(HS); for(int l = 0; l < lay && l < nbt; l++) free(BT[l]); free(BT); free(BM); free(pk.p); free(kj); free(own); free(qs); TGT = wtime() - tt0;
    #undef OW
    *S = kval(&ks); *A = a; *e1 = x1; *e2 = x2; *drop = dr; *nev += ne; *nq += q; *un += u;
    return 1;
}
