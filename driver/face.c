/* face.c -- experimental: sum the tail of the lattice along a "face" direction in closed form.
 *
 * SINCQED_FACE="GRAPH:u0,u1,...,u_{D-1}:a[:m[:W[:M]]]" registers two extra graphs next to GRAPH:
 *   GRAPH_fd  (D dims)    F(k) - T(k)
 *   GRAPH_ff  (D-1 dims)  sum over n of T(p + n u), done in closed form
 * whose lattice sums add up exactly to the lattice sum of GRAPH (same h, any eps).  u has entries in
 * {-1,0,1}; along it the integrand decays as exp(-a t) (a an integer: measured with `sincqed probe`).
 * Every lattice point is k = p + n u with p_i = 0 at the pivot i (first nonzero u_i), n = u_i k_i.
 * On the line through p the asymptote is T = G(p) r^n, r = e^{-a h}, with G(p) read off one evaluation far
 * out on the line (n = n0 + M, M = ceil(20/h) by default).  T is switched on by a ramp chi that rises from 0
 * at n = n0(p) to 1 at n0 + W, n0(p) = max_j |p_j| + m.  The split is an exact rearrangement of the lattice
 * sum for any n0, m, W, M; they only decide how small F - T is, and so how much the pruning can drop.
 * Several splits separated by ';' are applied in order, and GRAPH may be a graph made by an earlier split
 * (e.g. "g3_00:...;g3_00_fd:...", giving g3_00_fd_fd and g3_00_fd_ff).  Up to NSLOT splits. */
#include "sq.h"

#define NSLOT 8
typedef struct { const graph_t *bg; graph_t gd, gf; char nd[96], nf[96]; int u[MAXD], piv, m, w, M; double a; } face_t;
static face_t FS[NSLOT];

static inline int n0of(const face_t *F, const int *p){ int mx = 0; for(int j = 0; j < F->bg->D; j++) if(j != F->piv && abs(p[j]) > mx) mx = abs(p[j]); return mx + F->m; }
static inline double chi(const face_t *F, int n, int n0){ double x = (double)(n - n0)/F->w; return x <= 0 ? 0 : x >= 1 ? 1 : x; }
static inline int farN(const face_t *F, int n0){ return n0 + (F->M > 0 ? F->M : (int)ceil(20.0/H)); }
static inline double usum(const face_t *F){ double s = 0; for(int j = 0; j < F->bg->D; j++) s += F->u[j]; return s; }
static double linesum(const face_t *F, int n0, int N){
    double r = exp(-F->a*H), s = 0;
    for(int j = 1; j < F->w; j++) s += chi(F, n0 + j, n0)*exp(-F->a*H*(n0 + j - N));
    return s + exp(-F->a*H*(n0 + F->w - N))/(1 - r);
}
/* _fd: f - t, t/prod(s) at k = f(far) prod(s_far)/prod(s_k) r^{n-N} = f(far) e^{(N-n) h (sum u + a)} */
#define FD_BODY(T_, FN) \
    int D = F->bg->D, k[MAXD], p[MAXD]; \
    for(int j = 0; j < D; j++) k[j] = (int)lrint(log((double)sf[j])/H); \
    T_ v = F->bg->FN(sf, ab); \
    int n = F->u[F->piv]*k[F->piv]; for(int j = 0; j < D; j++) p[j] = k[j] - n*F->u[j]; \
    int n0 = n0of(F, p); double c = chi(F, n, n0); if(c == 0) return v; \
    int N = farN(F, n0); T_ sg[MAXD], abf; \
    for(int j = 0; j < D; j++) sg[j] = (T_)exp((p[j] + N*F->u[j])*H); \
    T_ vf = F->bg->FN(sg, &abf); \
    double fac = c*exp((N - n)*H*(usum(F) + F->a)); \
    *ab = *ab + abf*(T_)fac; return v - vf*(T_)fac;
/* _ff: value before the driver's h^{D-1} prod s: sum_n chi r^n G(p) / (h^{D-1} prod_{j != piv} s_j(p)),
 * G(p) = F(far) h^D r^{-N}, F(far) = f(far) prod s_far, prod s_far / prod_{j != piv} s_j(p) = e^{N h sum u} */
#define FF_BODY(T_, FN) \
    int D = F->bg->D, p[MAXD]; T_ sg[MAXD]; \
    for(int j = 0, q = 0; j < D; j++) p[j] = j == F->piv ? 0 : (int)lrint(log((double)sf[q++])/H); \
    int n0 = n0of(F, p), N = farN(F, n0); \
    for(int j = 0; j < D; j++) sg[j] = (T_)exp((p[j] + N*F->u[j])*H); \
    double fac = H*exp(N*H*usum(F))*linesum(F, n0, N); \
    T_ abf; T_ vf = F->bg->FN(sg, &abf); *ab = abf*(T_)fac; return vf*(T_)fac;
#define SLOT(i) \
    static double fdd##i(const double *sf, double *ab){ const face_t *F = &FS[i]; FD_BODY(double, fd) } \
    static QREAL  fdq##i(const QREAL *sf, QREAL *ab){ const face_t *F = &FS[i]; FD_BODY(QREAL, fq) } \
    static double ffd##i(const double *sf, double *ab){ const face_t *F = &FS[i]; FF_BODY(double, fd) } \
    static QREAL  ffq##i(const QREAL *sf, QREAL *ab){ const face_t *F = &FS[i]; FF_BODY(QREAL, fq) }
SLOT(0) SLOT(1) SLOT(2) SLOT(3) SLOT(4) SLOT(5) SLOT(6) SLOT(7)
static const fd_t FDD[NSLOT] = {fdd0, fdd1, fdd2, fdd3, fdd4, fdd5, fdd6, fdd7}, FFD[NSLOT] = {ffd0, ffd1, ffd2, ffd3, ffd4, ffd5, ffd6, ffd7};
static const fq_t FDQ[NSLOT] = {fdq0, fdq1, fdq2, fdq3, fdq4, fdq5, fdq6, fdq7}, FFQ[NSLOT] = {ffq0, ffq1, ffq2, ffq3, ffq4, ffq5, ffq6, ffq7};

void face_register(void){
    const char *e = getenv("SINCQED_FACE"); if(!e || !*e) return;
    char all[2048]; snprintf(all, sizeof all, "%s", e);
    int ns = 0; char *save = NULL;
    for(char *spec = strtok_r(all, ";", &save); spec; spec = strtok_r(NULL, ";", &save)){
        if(ns == NSLOT){ fprintf(stderr, "SINCQED_FACE: at most %d splits\n", NSLOT); exit(1); }
        face_t *F = &FS[ns]; char *s2 = NULL;
        char *gname = strtok_r(spec, ":", &s2), *us = strtok_r(NULL, ":", &s2), *as = strtok_r(NULL, ":", &s2),
             *ms = strtok_r(NULL, ":", &s2), *ws = strtok_r(NULL, ":", &s2), *Ms = strtok_r(NULL, ":", &s2);
        if(!gname || !us || !as){ fprintf(stderr, "SINCQED_FACE=GRAPH:u0,u1,...:a[:m[:W[:M]]][;...]\n"); exit(1); }
        F->bg = find(gname); if(!F->bg){ fprintf(stderr, "SINCQED_FACE: no graph %s\n", gname); exit(1); }
        int d = 0; char *s3 = NULL; for(char *t = strtok_r(us, ",", &s3); t; t = strtok_r(NULL, ",", &s3)) if(d < MAXD) F->u[d++] = atoi(t);
        if(d != F->bg->D){ fprintf(stderr, "SINCQED_FACE: u has %d entries, %s has D = %d\n", d, F->bg->name, F->bg->D); exit(1); }
        F->piv = -1; for(int j = 0; j < d; j++) if(F->u[j]){ F->piv = j; break; }
        if(F->piv < 0){ fprintf(stderr, "SINCQED_FACE: u = 0\n"); exit(1); }
        F->a = atof(as); F->m = ms ? atoi(ms) : 0; F->w = ws ? atoi(ws) : 4; if(F->w < 1) F->w = 1; F->M = Ms ? atoi(Ms) : 0;
        snprintf(F->nd, sizeof F->nd, "%s_fd", F->bg->name); snprintf(F->nf, sizeof F->nf, "%s_ff", F->bg->name);
        F->gd = *F->bg; F->gd.name = F->nd; F->gd.volkov = F->nd; F->gd.fd = FDD[ns]; F->gd.fq = FDQ[ns]; F->gd.fv = NULL; F->gd.ptd = NULL; F->gd.ptq = NULL;
        F->gf = *F->bg; F->gf.name = F->nf; F->gf.volkov = F->nf; F->gf.D = F->bg->D - 1; F->gf.fd = FFD[ns]; F->gf.fq = FFQ[ns]; F->gf.fv = NULL; F->gf.ptd = NULL; F->gf.ptq = NULL;
        add_graph(&F->gd); add_graph(&F->gf); ns++;
    }
}
