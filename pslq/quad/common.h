/* Shared definitions for the generated graph kernels and the driver. */
#pragma once
#include <math.h>
#include <float.h>
#if defined(QREAL_QD)                /* quad-double from plain doubles (qd.h), portable */
#include "qd.h"
typedef qd QREAL;
#define QEPS 2.5e-63
#define QNAME "quad-double"
#elif defined(QREAL_QQ)              /* double-quad, two __float128 (x86-64 with libquadmath only) */
#include "qq.h"
typedef qq QREAL;
#define QEPS 1e-66
#define QNAME "double-quad"
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
 *      scaling tests.
 * ud/uq: optional integrand in all Du free Schwinger parameters, normalized (sincqed demo); graphs without
 *      lepton loops leave it NULL and the driver builds it from pt*. */
typedef double (*ptd_t)(const double*, const double, double*);
typedef QREAL  (*ptq_t)(const QREAL*, const QREAL, QREAL*);
typedef struct { const char *name, *volkov; int L, D, weight, n, fixed; fd_t fd; fq_t fq; ptd_t ptd; ptq_t ptq; fd_t ud; fq_t uq; int Du; } graph_t;
