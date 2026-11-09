/* lattice.c -- the Sinc lattice and its map to the integration variables z = ln s.
 * Node k (integers, |k_i| <= K) sits at z_i = k_i h; the quadrature weight is h^D prod_i s_i. */
#include "sq.h"
#ifdef SQ_GPU
#include "gpu.h"
#endif

double H;
int K, DIM;
fd_t FD; fq_t FQ;
static double *EX;                          /* e^{k h}, k = -K..K */

void lattice_setup(const graph_t *g, double h, int N, int top){
    H = h; DIM = g->D; FD = g->fd; FQ = g->fq;
    int m = 1 << (top+1); K = (N-1)/2; K = ((K + m - 1)/m)*m;   /* rounded UP to a multiple of 2^(top+1) */
    free(EX); EX = (double*)malloc((2*K+1)*sizeof(double));
    for(int k = -K; k <= K; k++) EX[k+K] = exp(k*H);
#ifdef SQ_GPU
    USEGPU = gpu_setup(g->name, g->D, EX, K, H);
#ifdef SQ_DDV
    { const char *e = getenv("SINCQED_DDV"); FQV = NULL;
      if(!(e && !strcmp(e, "0"))) for(int i = 0; i < DDV_N; i++) if(DDV_ALL[i] && !strcmp(DDV_ALL[i]->name, g->name) && DDV_ALL[i]->D == g->D){ FQV = DDV_ALL[i]; break; } }
#endif
#endif
}

/* s = e^z and the weight h^D prod(s) for node k */
void lattice_node(const int *k, double *sf, double *w){
    double ww = 1;
    for(int i = 0; i < DIM; i++){ sf[i] = EX[k[i]+K]; ww *= sf[i]*H; }
    *w = ww;
}
