/* gpu.cu -- host side of the CUDA back end.  The double half of a batch of lattice nodes runs on the
 * device; eval.c does the |.|-sum test and the double-double fallback on the host.  Two slots, each with
 * its own stream and buffers, so the driver can prepare and finish one batch while the other computes:
 *   gpu_submit(slot, ...)   pack the nodes (int16), copy, launch, copy back: returns at once
 *   gpu_wait(slot, ...)     wait for that slot; f w and sum|terms| w in pinned host buffers
 * SINCQED_GPU=0 turns the GPU off (CPU path), SINCQED_GPU=d picks device d-1 (default: device 0). */
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "gpu.h"

extern "C" const sq_gpu_t *GPU_ALL[];          /* build/gpu/registry.c */

#define CK(x) do { cudaError_t e_ = (x); if(e_ != cudaSuccess){ fprintf(stderr, "CUDA error %s at %s:%d: %s\n", #x, __FILE__, __LINE__, cudaGetErrorString(e_)); exit(3); } } while(0)

static int ON = -1, DEV = 0, D_, K_; static double H_;
static const sq_gpu_t *CUR;
static double *dEX; static long nEX;
typedef struct {
    cudaStream_t st; cudaEvent_t ev; long cap, n;
    short *hK, *dK; unsigned char *hS, *dS; double *hF, *hA, *dF, *dA;
} slot_t;
static slot_t SL[2];
static char INFO[256];

static int init(void){
    if(ON >= 0) return ON;
    const char *e = getenv("SINCQED_GPU"); ON = 0;
    if(e && !strcmp(e, "0")) return 0;
    int nd = 0; if(cudaGetDeviceCount(&nd) != cudaSuccess || nd == 0) return 0;
    DEV = e && atoi(e) > 0 ? atoi(e) - 1 : 0; if(DEV >= nd) DEV = 0;
    CK(cudaSetDevice(DEV));
    cudaDeviceProp p; CK(cudaGetDeviceProperties(&p, DEV));
    snprintf(INFO, sizeof INFO, "%s (sm_%d%d, %d SMs)", p.name, p.major, p.minor, p.multiProcessorCount);
    for(int b = 0; b < 2; b++){ CK(cudaStreamCreateWithFlags(&SL[b].st, cudaStreamNonBlocking)); CK(cudaEventCreateWithFlags(&SL[b].ev, cudaEventDisableTiming)); }
    return ON = 1;
}

static void reserve(slot_t *s, long n){
    if(n <= s->cap) return;
    if(s->cap){ CK(cudaFreeHost(s->hK)); CK(cudaFree(s->dK)); CK(cudaFreeHost(s->hS)); CK(cudaFree(s->dS));
                CK(cudaFreeHost(s->hF)); CK(cudaFree(s->dF)); CK(cudaFreeHost(s->hA)); CK(cudaFree(s->dA)); }
    long c = 1L << 16; while(c < n) c <<= 1; s->cap = c;
    CK(cudaMallocHost(&s->hK, c*16*sizeof(short))); CK(cudaMalloc(&s->dK, c*16*sizeof(short)));
    CK(cudaMallocHost(&s->hS, c)); CK(cudaMalloc(&s->dS, c));
    CK(cudaMallocHost(&s->hF, c*sizeof(double))); CK(cudaMalloc(&s->dF, c*sizeof(double)));
    CK(cudaMallocHost(&s->hA, c*sizeof(double))); CK(cudaMalloc(&s->dA, c*sizeof(double)));
}

const char *gpu_info(void){ return ON == 1 ? INFO : "off"; }

int gpu_setup(const char *name, int D, const double *EX, int K, double H){
    CUR = NULL;
    if(!init()) return 0;
    for(int i = 0; GPU_ALL[i]; i++) if(!strcmp(GPU_ALL[i]->name, name)){ CUR = GPU_ALL[i]; break; }
    if(!CUR) return 0;
    if(CUR->D != D){ fprintf(stderr, "gpu: %s has D=%d, expected %d\n", name, CUR->D, D); exit(3); }
    D_ = D; K_ = K; H_ = H;
    CK(cudaDeviceSynchronize());
    if(nEX < 2*K + 1){ if(dEX) CK(cudaFree(dEX)); nEX = 2*K + 1; CK(cudaMalloc(&dEX, nEX*sizeof(double))); }
    CK(cudaMemcpy(dEX, EX, (2*K + 1)*sizeof(double), cudaMemcpyHostToDevice));
    return 1;
}

/* n nodes k[i*stride + d]; skip may be NULL (all |.| sums) */
void gpu_submit(int slot, long n, const int *k, int stride, const unsigned char *skip){
    slot_t *s = &SL[slot]; reserve(s, n); s->n = n;
    if(n <= 0) return;
    short *hk = s->hK; int D = D_;
    #pragma omp parallel for schedule(static)
    for(long i = 0; i < n; i++) for(int d = 0; d < D; d++) hk[i*D + d] = (short)k[i*stride + d];
    CK(cudaMemcpyAsync(s->dK, hk, n*D*sizeof(short), cudaMemcpyHostToDevice, s->st));
    if(skip){ memcpy(s->hS, skip, n); CK(cudaMemcpyAsync(s->dS, s->hS, n, cudaMemcpyHostToDevice, s->st)); }
    CUR->launch(s->st, n, s->dK, skip ? s->dS : NULL, dEX, K_, H_, s->dF, s->dA, 0, 0);
    CK(cudaGetLastError());
    CK(cudaMemcpyAsync(s->hF, s->dF, n*sizeof(double), cudaMemcpyDeviceToHost, s->st));
    CK(cudaMemcpyAsync(s->hA, s->dA, n*sizeof(double), cudaMemcpyDeviceToHost, s->st));
    CK(cudaEventRecord(s->ev, s->st));
}

void gpu_submit_cells(int slot, long nc, const int *cells, int stride, const unsigned char *skip, int st){
    slot_t *s = &SL[slot]; long n = nc << D_; reserve(s, n); s->n = n;
    if(nc <= 0) return;
    short *hk = s->hK; int D = D_;
    #pragma omp parallel for schedule(static)
    for(long c = 0; c < nc; c++) for(int d = 0; d < D; d++) hk[c*D + d] = (short)cells[c*stride + d];
    CK(cudaMemcpyAsync(s->dK, hk, nc*D*sizeof(short), cudaMemcpyHostToDevice, s->st));
    if(skip){ memcpy(s->hS, skip, nc); CK(cudaMemcpyAsync(s->dS, s->hS, nc, cudaMemcpyHostToDevice, s->st)); }
    CUR->launch(s->st, n, s->dK, skip ? s->dS : NULL, dEX, K_, H_, s->dF, s->dA, 1, st);
    CK(cudaGetLastError());
    CK(cudaMemcpyAsync(s->hF, s->dF, n*sizeof(double), cudaMemcpyDeviceToHost, s->st));
    CK(cudaMemcpyAsync(s->hA, s->dA, n*sizeof(double), cudaMemcpyDeviceToHost, s->st));
    CK(cudaEventRecord(s->ev, s->st));
}

/* fw = f w, abw = sum|terms| w (0 where skipped); valid until the slot's next submit */
void gpu_wait(int slot, const double **fw, const double **abw){
    slot_t *s = &SL[slot];
    if(s->n > 0) CK(cudaEventSynchronize(s->ev));
    *fw = s->hF; *abw = s->hA;
}
