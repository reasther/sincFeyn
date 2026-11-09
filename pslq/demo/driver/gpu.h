/* gpu.h -- the CUDA back end (make GPU=1): the double half of eval_batch on the GPU. */
#pragma once
#ifdef __cplusplus
extern "C" {
#endif
#ifdef __CUDACC__
typedef cudaStream_t sq_stream_t;
#else
typedef void *sq_stream_t;
#endif
/* launcher of one graph's batch kernel (gpu_kernel.h), on the device buffers */
typedef void (*sq_gpu_launch_t)(sq_stream_t st, long n, const short *k, const unsigned char *skip,
                                const double *EX, int K, double H, double *of, double *oab, int cmode, int cst);
typedef struct { const char *name; int D; sq_gpu_launch_t launch; } sq_gpu_t;

/* gpu.cu */
int  gpu_setup(const char *name, int D, const double *EX, int K, double H);   /* 1 if this graph has a GPU kernel */
void gpu_submit(int slot, long n, const int *k, int stride, const unsigned char *skip);
void gpu_wait(int slot, const double **fw, const double **abw);
/* the children of nc cells (cells[c*stride + d], skip per cell) at step st: nc << D results, in prep() order */
void gpu_submit_cells(int slot, long nc, const int *cells, int stride, const unsigned char *skip, int st);
const char *gpu_info(void);
#ifdef __cplusplus
}
#endif
