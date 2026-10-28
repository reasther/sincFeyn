/* ddv.h -- double-double arithmetic on SIMD vectors: VW points at once, one per lane (C++).
 * The operations of dd.h, lane by lane: Dekker/Knuth two-sum and fma two-prod (as in the QD library), the
 * "sloppy" add unless DD_IEEE_ADD, and dd op double without promoting the double, so each lane gives what
 * dd.h gives for that point.  Used by the vectorised fallback kernels made by
 * scripts/ddvgen.py: the double-double half of a graph has no branches on values (only on SQ_AB, and
 * RABS, done here with a lane-wise select), so VW fallback points can run in lockstep.
 * VW = 8 with AVX-512, 4 with AVX2/FMA, else 2 (fma per lane). */
#pragma once
#include <math.h>
#if defined(__AVX512F__) || defined(__FMA__)
#include <immintrin.h>
#endif

#if defined(__AVX512F__)
#define SQ_DDVW 8
#elif defined(__FMA__)
#define SQ_DDVW 4
#else
#define SQ_DDVW 2
#endif
typedef double vd __attribute__((vector_size(8*SQ_DDVW)));
typedef long long vl __attribute__((vector_size(8*SQ_DDVW)));

static inline vd vfma(vd a, vd b, vd c){
#if defined(__AVX512F__)
    return (vd)_mm512_fmadd_pd((__m512d)a, (__m512d)b, (__m512d)c);
#elif defined(__FMA__)
    return (vd)_mm256_fmadd_pd((__m256d)a, (__m256d)b, (__m256d)c);
#else
    vd r; for(int i = 0; i < SQ_DDVW; i++) r[i] = fma(a[i], b[i], c[i]); return r;
#endif
}

struct ddv {
    vd hi, lo;
    ddv() : hi(vd{} ), lo(vd{}) {}
    ddv(double x) : hi(vd{} + x), lo(vd{}) {}
    ddv(vd h, vd l) : hi(h), lo(l) {}
};

static inline ddv ddv_two_sum(vd a, vd b){ vd s = a + b, bb = s - a; return ddv(s, (a - (s - bb)) + (b - bb)); }
static inline ddv ddv_quick_two_sum(vd a, vd b){ vd s = a + b; return ddv(s, b - (s - a)); }
static inline ddv ddv_two_prod(vd a, vd b){ vd p = a*b; return ddv(p, vfma(a, b, -p)); }

#ifdef DD_IEEE_ADD
static inline ddv operator+(const ddv &a, const ddv &b){
    ddv s = ddv_two_sum(a.hi, b.hi), t = ddv_two_sum(a.lo, b.lo);
    s.lo += t.hi; s = ddv_quick_two_sum(s.hi, s.lo); s.lo += t.lo; return ddv_quick_two_sum(s.hi, s.lo);
}
#else
/* the "sloppy" add of dd.h (QD's default) */
static inline ddv operator+(const ddv &a, const ddv &b){
    ddv s = ddv_two_sum(a.hi, b.hi); s.lo += a.lo + b.lo; return ddv_quick_two_sum(s.hi, s.lo);
}
#endif
/* dd op double of dd.h: the double is not promoted */
static inline ddv ddv_add_d(const ddv &a, double b){ ddv s = ddv_two_sum(a.hi, vd{} + b); s.lo += a.lo; return ddv_quick_two_sum(s.hi, s.lo); }
static inline ddv ddv_mul_d(const ddv &a, double b){ ddv p = ddv_two_prod(a.hi, vd{} + b); p.lo += a.lo*b; return ddv_quick_two_sum(p.hi, p.lo); }
static inline ddv operator-(const ddv &a){ return ddv(-a.hi, -a.lo); }
static inline ddv operator-(const ddv &a, const ddv &b){ return a + (-b); }
static inline ddv operator*(const ddv &a, const ddv &b){
    ddv p = ddv_two_prod(a.hi, b.hi); p.lo += a.hi*b.lo + a.lo*b.hi; return ddv_quick_two_sum(p.hi, p.lo);
}
static inline ddv operator/(const ddv &a, const ddv &b){
    vd q1 = a.hi/b.hi; ddv r = a - b*ddv(q1, vd{});
    vd q2 = r.hi/b.hi; r = r - b*ddv(q2, vd{});
    vd q3 = r.hi/b.hi;
    ddv q = ddv_quick_two_sum(q1, q2); return q + ddv(q3, vd{});
}
static inline ddv &operator+=(ddv &a, const ddv &b){ a = a + b; return a; }
static inline ddv &operator-=(ddv &a, const ddv &b){ a = a - b; return a; }
static inline ddv &operator*=(ddv &a, const ddv &b){ a = a * b; return a; }
static inline ddv &operator/=(ddv &a, const ddv &b){ a = a / b; return a; }
/* mixed forms with int and double operands, as in dd.h (the kernels contain e.g. 128+x, (REAL)(-2)) */
static inline ddv operator+(int a, const ddv &b){ return ddv_add_d(b, (double)a); }
static inline ddv operator+(const ddv &a, int b){ return ddv_add_d(a, (double)b); }
static inline ddv operator-(int a, const ddv &b){ return ddv_add_d(-b, (double)a); }
static inline ddv operator-(const ddv &a, int b){ return ddv_add_d(a, -(double)b); }
static inline ddv operator*(int a, const ddv &b){ return ddv_mul_d(b, (double)a); }
static inline ddv operator*(const ddv &a, int b){ return ddv_mul_d(a, (double)b); }
static inline ddv operator/(int a, const ddv &b){ return ddv((double)a) / b; }
static inline ddv operator/(const ddv &a, int b){ return a / ddv((double)b); }
static inline ddv operator+(double a, const ddv &b){ return ddv_add_d(b, a); }
static inline ddv operator+(const ddv &a, double b){ return ddv_add_d(a, b); }
static inline ddv operator-(double a, const ddv &b){ return ddv_add_d(-b, a); }
static inline ddv operator-(const ddv &a, double b){ return ddv_add_d(a, -b); }
static inline ddv operator*(double a, const ddv &b){ return ddv_mul_d(b, a); }
static inline ddv operator*(const ddv &a, double b){ return ddv_mul_d(a, b); }
/* RABS(x) = x < 0 ? -x : x of dd.h (x < 0: hi < 0, or hi == 0 and lo < 0), lane by lane */
static inline ddv ddv_abs(const ddv &x){
    vl m = (x.hi < 0) | ((x.hi == 0) & (x.lo < 0));
    return ddv(m ? -x.hi : x.hi, m ? -x.lo : x.lo);
}
