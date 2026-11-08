/* ddvk.h -- the SIMD double-double fallback kernels (tools/ddvgen.py, make DDV=1, the default with GPU=1):
 * fv(m, sf, v, ab) evaluates m <= w points (sf: m x D Sinc variables) in double-double, one per lane, and
 * returns v = (double) f and ab = (double) sum|terms| per point, as the scalar FQ does. */
#pragma once
#ifdef __cplusplus
extern "C" {
#endif
typedef void (*sq_ddv_fn_t)(int m, const double *sf, double *v, double *ab);
typedef struct { const char *name; int D, w; sq_ddv_fn_t fv; } sq_ddv_t;
#ifdef __cplusplus
}
#endif
