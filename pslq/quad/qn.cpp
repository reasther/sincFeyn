/* Sinc sum of one graph's reduced integrand f_q in extended precision (QREAL: quad-double qd.h, double-quad qq.h
 * or __float128), flood fill on |density| > eps.  Graph chosen at build time (-DGRAPH=G_g2_00).  usage: corner_qd h eps [zmax] */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include "common.h"
#include "qio.h"
__thread int SQ_AB = 1;
extern const graph_t GRAPH;
static const graph_t *G = &GRAPH;
typedef QREAL Q;
static Q *EX; static int K, D;
static uint64_t pack(const int *k){ uint64_t c = 0; for(int i = 0; i < D; i++) c = c << 16 | (uint16_t)(k[i] + 32768); return c; }
static void unpack(uint64_t c, int *k){ for(int i = D-1; i >= 0; i--){ k[i] = (int)(c & 0xffff) - 32768; c >>= 16; } }
int main(int argc, char **argv){
    double h = atof(argv[1]), eps = atof(argv[2]), zmax = argc > 3 ? atof(argv[3]) : 150;
    D = G->D; K = (int)(zmax/h); int N = 2*K+1;
    Q H = Q(h);
    Q E = qexp(H);
    Q Ei = 1/E;
    EX = new Q[N]; EX[K] = Q(1); for(int k = 1; k <= K; k++){ EX[K+k] = EX[K+k-1]*E; EX[K-k] = EX[K-k+1]*Ei; }
    std::unordered_set<uint64_t> vis; std::vector<uint64_t> cur; int k0[8] = {0}; cur.push_back(pack(k0)); vis.insert(cur[0]);
    Q S = Q(0); long nev = 0; int hit = 0; double far = 0, mz = 0;
    while(!cur.empty()){
        long nc = cur.size(); std::vector<Q> v(nc); std::vector<char> ex(nc); std::vector<double> nzv(nc);
        #pragma omp parallel for schedule(dynamic, 64)
        for(long i = 0; i < nc; i++){
            int k[8]; unpack(cur[i], k); Q sf[8], ab, pr = Q(1);
            for(int j = 0; j < D; j++){ sf[j] = EX[k[j]+K]; pr = pr*sf[j]; }
            Q f = G->fq(sf, &ab); v[i] = f*pr;
            double nz = hi_(ab*pr)*16*QEPS; nzv[i] = nz;
            ex[i] = fabs(hi_(v[i])) > eps && nz < 0.1*eps;
        }
        std::vector<uint64_t> nx;
        for(long i = 0; i < nc; i++){
            S = S + v[i]; nev++; if(nzv[i] > mz) mz = nzv[i];
            if(!ex[i]) continue;
            int k[8]; unpack(cur[i], k);
            for(int j = 0; j < D; j++){ double z = abs(k[j])*h; if(z > far) far = z; }
            for(int j = 0; j < D; j++) for(int sg = -1; sg <= 1; sg += 2){
                k[j] += sg;
                if(abs(k[j]) > K) hit = 1;
                else { uint64_t c = pack(k); if(vis.insert(c).second) nx.push_back(c); }
                k[j] -= sg;
            }
        }
        cur.swap(nx);
    }
    Q w = Q(G->weight); for(int j = 0; j < D; j++) w = w*H;
    S = S*w;
    char b1[400]; qout(S, b1);
    fprintf(stderr, "# %s, %d threads\n", QNAME, omp_get_max_threads());
    printf("%s h=%.4f eps=%.0e  value = %s  evals %ld  reach %.1f  max noise %.1e%s\n", G->name, h, eps, b1, nev, far, mz, hit ? "  BOX HIT" : "");
}
