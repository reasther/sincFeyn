/* g2_01 (self-energy): reduced 3-D integrand summed in extended precision (quad-double, ~64 digits, or
 * double-quad), one-point rule in the external-vertex split (exact: see below).  usage: se_qd h eps [zmax] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "common.h"
#include "qio.h"
__thread int SQ_AB = 1;
extern const graph_t G_g2_01;
typedef QREAL Q;
static Q H, *EX; static double NZ[256]; static int K; static Q GX, GW;
static Q qabs(const Q &x){ return x < 0 ? -x : x; }
static Q fval(int a, int b, int c, double *noise){
    Q sig = EX[a+K], s2 = EX[b+K], s5 = EX[c+K], acc = 0, ab, aa = 0;
    {   Q u = GX, sb = sig*u, sa = sig*(1-u), s[6];
        s[4] = 1; s[2] = s2; s[5] = s5; s[0] = sb; s[1] = sa/2; s[3] = sa/2;
        acc = GW*sa*G_g2_01.ptq(s, 0, &ab); aa = GW*sa*ab;
    }
    *noise = hi_(aa*sig/256*sig*s2*s5)*16*QEPS;
    return acc*sig/256*sig*s2*s5;
}
int main(int argc, char **argv){
    double h = atof(argv[1]), eps = atof(argv[2]); double zmax = argc > 3 ? atof(argv[3]) : 200;
    H = Q(h); K = (int)(zmax/h); int N = 2*K+1;
    Q E = exp(H), Ei = 1/E;
    EX = new Q[N]; EX[K] = 1; for(int k = 1; k <= K; k++){ EX[K+k] = EX[K+k-1]*E; EX[K-k] = EX[K-k+1]*Ei; }
    // The split u of the external-vertex line enters linearly, f = f0 + f1 u, under the weight sa = sig (1-u),
    // so one point at the weighted mean u = 1/3 with weight 3/4 integrates it exactly (tools/glpts.py in
    // sincFeyn).  Both are formed in Q, so the rule is exact to the working precision.
    GX = Q(1)/Q(3); GW = Q(3)/Q(4);
    unsigned char *vis = (unsigned char*)calloc((size_t)N*N*N, 1);
    #define IDX(a,b,c) (((size_t)((a)+K)*N + ((b)+K))*N + ((c)+K))
    int (*cur)[3] = (int(*)[3])malloc(sizeof *cur); long nc = 1; cur[0][0] = cur[0][1] = cur[0][2] = 0; vis[IDX(0,0,0)] = 1;
    Q S = 0; long nev = 0; int hit = 0; double far = 0;
    while(nc){
        unsigned char *ex = (unsigned char*)malloc(nc); Q *v = new Q[nc];
        #pragma omp parallel for schedule(dynamic, 64)
        for(long i = 0; i < nc; i++){ double nz; v[i] = fval(cur[i][0], cur[i][1], cur[i][2], &nz); ex[i] = hi_(qabs(v[i])) > eps && nz < 0.1*eps; if(nz > NZ[omp_get_thread_num()]) NZ[omp_get_thread_num()] = nz; }
        long cap = 6*nc + 1, nn = 0; int (*nx)[3] = (int(*)[3])malloc(cap*sizeof *nx);
        for(long i = 0; i < nc; i++){
            S += v[i]; nev++;
            if(!ex[i]) continue;
            for(int j = 0; j < 3; j++){ double z = abs(cur[i][j])*h; if(z > far) far = z; }
            for(int j = 0; j < 3; j++) for(int sg = -1; sg <= 1; sg += 2){
                int k[3] = {cur[i][0], cur[i][1], cur[i][2]}; k[j] += sg;
                if(abs(k[j]) > K){ hit = 1; continue; }
                size_t id = IDX(k[0], k[1], k[2]); if(vis[id]) continue; vis[id] = 1;
                memcpy(nx[nn++], k, sizeof k);
            }
        }
        free(ex); delete[] v; free(cur); cur = nx; nc = nn;
    }
    S = S*2*H*H*H;
    char b1[400]; qout(S, b1);
    double mz = 0; for(int t = 0; t < 256; t++) if(NZ[t] > mz) mz = NZ[t];
    fprintf(stderr, "# %s, %d threads\n", QNAME, omp_get_max_threads());
    printf("g2_01 h=%.4f eps=%.0e  value = %s  evals %ld  reach %.1f  max noise %.1e%s\n", h, eps, b1, nev, far, mz, hit ? "  BOX HIT" : "");
}
