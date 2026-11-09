/* Precision-independent helpers for the summers: leading double, decimal output. */
#pragma once
#include <stdio.h>
#if defined(QREAL_QD)
static inline double hi_(const QREAL &x){ return x.x[0]; }
static inline void qout(const QREAL &x, char *b){ char s[160], p[100]; qd_str(x, 40, s); qd_parts(x, p); sprintf(b, "%s  [parts %s]", s, p); }
static inline QREAL qexp(const QREAL &x){ return exp(x); }
#elif defined(QREAL_QQ)
static inline double hi_(const QREAL &x){ return (double)x.hi; }
static inline void qout(const QREAL &x, char *b){ char s1[64], s2[64]; quadmath_snprintf(s1, 64, "%+.36Qe", x.hi); quadmath_snprintf(s2, 64, "%+.36Qe", x.lo); sprintf(b, "%s %s", s1, s2); }
static inline QREAL qexp(const QREAL &x){ return exp(x); }
#else
#include <quadmath.h>
static inline double hi_(const QREAL &x){ return (double)x; }
static inline void qout(const QREAL &x, char *b){ quadmath_snprintf(b, 64, "%+.34Qe", x); }
static inline QREAL qexp(const QREAL &x){ return expq(x); }
#endif
