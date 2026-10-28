/* gpu_kernel.h -- included by the CUDA kernels that scripts/gpugen.py makes from kernels/<class>/*.c.
 * One thread per lattice node: s = e^{k h} from the host's table, w = prod(s h) in the same order as
 * lattice_node(), then f(s) in double with the |.| sums only where skip[i] == 0. */
#pragma once
#include <cuda_runtime.h>
#include "gpu.h"
#undef RABS
#define RABS(x) fabs(x)
#ifndef SQ_CHUNK
#define SQ_CHUNK __device__ __noinline__       /* the chunk functions: -DSQ_CHUNK="__device__ __forceinline__" inlines them */
#endif
#ifndef SQ_MINB
#define SQ_MINB 1
#endif
#ifndef SQ_BLOCK
#define SQ_BLOCK 128
#endif

/* Two modes.  Nodes (cmode = 0): k[i*D + d] is node i.  Cells (cmode = 1): thread i is child b = i mod 2^D
 * of cell c = i / 2^D, k = cell - st b (bit d of b: coordinate d), exactly as prep() in sum.c makes them;
 * a child below the box (k_d < -K) is not a node and gives 0, 0; skip[c] is the cell's flag. */
#define SQ_GPU_KERNEL(NAME, D)                                                                         \
__global__ void __launch_bounds__(SQ_BLOCK, SQ_MINB) gk_##NAME(long n, const short *k, const unsigned char *skip, \
        const double *EX, int K, double H, double *of, double *oab, int cmode, int st){                  \
    long i = blockIdx.x*(long)blockDim.x + threadIdx.x; if(i >= n) return;                              \
    double sf[D], w = 1, ab = 0; int sk;                                                                 \
    if(cmode){ long c = i >> D; int b = (int)(i & ((1 << D) - 1));                                       \
        for(int d = 0; d < D; d++){ int kd = k[c*D + d] - (((b >> d) & 1) ? st : 0);                      \
            if(kd < -K){ of[i] = 0; oab[i] = 0; return; }                                                \
            double s = EX[kd + K]; sf[d] = s; w *= s*H; }                                                \
        sk = skip ? skip[c] : 0; }                                                                       \
    else { for(int d = 0; d < D; d++){ double s = EX[k[i*D + d] + K]; sf[d] = s; w *= s*H; }               \
        sk = skip ? skip[i] : 0; }                                                                       \
    double f = f_##NAME##_g(sf, &ab, !sk);                                                               \
    of[i] = f*w; oab[i] = ab*w;                                                                          \
}                                                                                                        \
static void gl_##NAME(cudaStream_t st, long n, const short *k, const unsigned char *skip,                \
        const double *EX, int K, double H, double *of, double *oab, int cmode, int cst){                 \
    if(n > 0) gk_##NAME<<<(unsigned)((n + SQ_BLOCK - 1)/SQ_BLOCK), SQ_BLOCK, 0, st>>>(n, k, skip, EX, K, H, of, oab, cmode, cst); \
}                                                                                                        \
extern "C" const sq_gpu_t GG_##NAME = {#NAME, D, gl_##NAME};
