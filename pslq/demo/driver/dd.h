/* Double-double arithmetic (about 32 significant digits) for the high-precision fallback.
 * A value is hi + lo with |lo| <= ulp(hi)/2.  Only + - * / and comparisons are needed: the
 * generated kernels are rational functions of the Schwinger parameters.  Plain FP64 operations
 * (with fma), so results are the same on x86, Apple silicon and GPUs.  C++ only (operator
 * overloading); built with `make DD=1` (recommended on Apple silicon).
 * Algorithms: Dekker/Knuth two-sum and fma two-prod, as in the QD library (Hida, Li, Bailey). */
#pragma once
#include <math.h>

struct dd {
    double hi, lo;
    dd() : hi(0), lo(0) {}
    dd(double x) : hi(x), lo(0) {}
    dd(double h, double l) : hi(h), lo(l) {}
    explicit operator double() const { return hi + lo; }
};

static inline dd dd_two_sum(double a, double b){ double s = a + b, bb = s - a; return dd(s, (a - (s - bb)) + (b - bb)); }
static inline dd dd_quick_two_sum(double a, double b){ double s = a + b; return dd(s, b - (s - a)); }
static inline dd dd_two_prod(double a, double b){ double p = a*b; return dd(p, fma(a, b, -p)); }

#ifdef DD_IEEE_ADD
static inline dd operator+(const dd &a, const dd &b){
    dd s = dd_two_sum(a.hi, b.hi), t = dd_two_sum(a.lo, b.lo);
    s.lo += t.hi; s = dd_quick_two_sum(s.hi, s.lo); s.lo += t.lo; return dd_quick_two_sum(s.hi, s.lo);
}
#else
/* "sloppy" add (QD's default): absolute error ~ 2^-104 (|a| + |b|), which is the error model the
 * fallback test assumes (QEPS * sum of |terms|); the IEEE add above also bounds the relative error. */
static inline dd operator+(const dd &a, const dd &b){
    dd s = dd_two_sum(a.hi, b.hi); s.lo += a.lo + b.lo; return dd_quick_two_sum(s.hi, s.lo);
}
#endif
/* dd op double: cheaper than promoting the double to dd */
static inline dd dd_add_d(const dd &a, double b){ dd s = dd_two_sum(a.hi, b); s.lo += a.lo; return dd_quick_two_sum(s.hi, s.lo); }
static inline dd dd_mul_d(const dd &a, double b){ dd p = dd_two_prod(a.hi, b); p.lo += a.lo*b; return dd_quick_two_sum(p.hi, p.lo); }
static inline dd operator-(const dd &a){ return dd(-a.hi, -a.lo); }
static inline dd operator-(const dd &a, const dd &b){ return a + (-b); }
static inline dd operator*(const dd &a, const dd &b){
    dd p = dd_two_prod(a.hi, b.hi); p.lo += a.hi*b.lo + a.lo*b.hi; return dd_quick_two_sum(p.hi, p.lo);
}
static inline dd operator/(const dd &a, const dd &b){
    double q1 = a.hi/b.hi; dd r = a - b*dd(q1);
    double q2 = r.hi/b.hi; r = r - b*dd(q2);
    double q3 = r.hi/b.hi;
    dd q = dd_quick_two_sum(q1, q2); return q + dd(q3);
}
static inline dd &operator+=(dd &a, const dd &b){ a = a + b; return a; }
static inline dd &operator-=(dd &a, const dd &b){ a = a - b; return a; }
static inline dd &operator*=(dd &a, const dd &b){ a = a * b; return a; }
static inline dd &operator/=(dd &a, const dd &b){ a = a / b; return a; }
static inline bool operator<(const dd &a, const dd &b){ return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
static inline bool operator>(const dd &a, const dd &b){ return b < a; }
static inline bool operator<=(const dd &a, const dd &b){ return !(b < a); }
static inline bool operator>=(const dd &a, const dd &b){ return !(a < b); }
static inline bool operator==(const dd &a, const dd &b){ return a.hi == b.hi && a.lo == b.lo; }
static inline bool operator!=(const dd &a, const dd &b){ return !(a == b); }
/* mixed forms with int literals (the kernels contain e.g. 128+x, (REAL)(-2)) */
static inline dd operator+(int a, const dd &b){ return dd_add_d(b, (double)a); }
static inline dd operator+(const dd &a, int b){ return dd_add_d(a, (double)b); }
static inline dd operator-(int a, const dd &b){ return dd_add_d(-b, (double)a); }
static inline dd operator-(const dd &a, int b){ return dd_add_d(a, -(double)b); }
static inline dd operator*(int a, const dd &b){ return dd_mul_d(b, (double)a); }
static inline dd operator*(const dd &a, int b){ return dd_mul_d(a, (double)b); }
static inline dd operator/(int a, const dd &b){ return dd((double)a) / b; }
static inline dd operator/(const dd &a, int b){ return a / dd((double)b); }
static inline dd operator+(double a, const dd &b){ return dd_add_d(b, a); }
static inline dd operator+(const dd &a, double b){ return dd_add_d(a, b); }
static inline dd operator-(double a, const dd &b){ return dd_add_d(-b, a); }
static inline dd operator-(const dd &a, double b){ return dd_add_d(a, -b); }
static inline dd operator*(double a, const dd &b){ return dd_mul_d(b, a); }
static inline dd operator*(const dd &a, double b){ return dd_mul_d(a, b); }
static inline bool operator<(const dd &a, int b){ return a < dd((double)b); }
static inline bool operator>(const dd &a, int b){ return a > dd((double)b); }
static inline dd exp(const dd &a){ return dd(exp(a.hi)); }   /* only used by the scaling test, where double is enough */
