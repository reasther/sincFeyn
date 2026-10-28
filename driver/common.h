/* Shared definitions for the generated graph kernels and the driver. */
#pragma once
#include <math.h>
#include <float.h>
#if defined(QREAL_DD)              /* double-double (dd.h, C++): make DD=1, recommended on Apple silicon */
#include "dd.h"
typedef dd QREAL;
#define QEPS 4.93e-32
#define QNAME "double-double"
#elif defined(QREAL_DOUBLE)          /* test only: no extended precision, as on Apple silicon */
typedef double QREAL;
#define QEPS DBL_EPSILON
#define QNAME "double (QREAL_DOUBLE)"
#elif defined(__SIZEOF_FLOAT128__) && !defined(NO_FLOAT128)
typedef __float128 QREAL;
#define QEPS 1.93e-34
#define QNAME "__float128"
#else
typedef long double QREAL;
#define QEPS ((double)LDBL_EPSILON)
#define QNAME "long double"
#endif
#define RABS(x) ((x)<0?-(x):(x))
/* SQ_AB = 0 lets a kernel skip the |.| companion sums of its double half (make inserts the test);
 * the driver clears it for points whose neighbourhood is far from needing the high-precision fallback */
extern __thread int SQ_AB;
typedef double (*fd_t)(const double*, double*);
typedef QREAL  (*fq_t)(const QREAL*, QREAL*);
/* f*: integrand in the D Sinc variables sf[] (= class sums, weight prod(sf h) applied by the driver)
 * pt*: raw point function of all n Schwinger parameters (sum over subtraction terms, before
 *      normalization; Vx is added to every term's exponent, used for VP insertions), used by the
 *      scaling tests. */
typedef double (*ptd_t)(const double*, const double, double*);
typedef QREAL  (*ptq_t)(const QREAL*, const QREAL, QREAL*);
/* fv (tools/packpts.py; NULL in kernels without it): f* in double for n points at once, point j at
 * sf[j*ld ...], the points in the SIMD lanes; out[j], ab[j] as fd's value and |.| sum */
typedef void (*fv_t)(long n, const double *sf, long ld, double *out, double *ab);
typedef struct { const char *name, *volkov; int L, D, weight, n, fixed; fd_t fd; fq_t fq; ptd_t ptd; ptq_t ptq; fv_t fv; } graph_t;
