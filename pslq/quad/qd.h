/* qd.h: quad-double arithmetic (four non-overlapping doubles, ~64 significant digits) from plain
 * IEEE double operations and fma, so it builds anywhere, Apple silicon included (no __float128).
 *
 * Every operation is done as an exact floating-point expansion (Shewchuk's two_sum and
 * grow-expansion, products split exactly with fma), then compressed and cut to four doubles.
 * So each +, -, *, / is correct to about 2^-208 relative; there are no sloppy shortcuts.
 * Compile with -ffp-contract=off (the Makefile does) so the compiler does not fuse the
 * error-free transformations into fma.  It is slower than a tuned quad-double library but simple
 * enough to check by eye.  Interface: the same as qq.h, so the generated kernels compile unchanged. */
#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>

struct qd {
    double x[4];
    qd(){ x[0] = x[1] = x[2] = x[3] = 0; }
    qd(double a){ x[0] = a; x[1] = x[2] = x[3] = 0; }
    qd(int a){ x[0] = a; x[1] = x[2] = x[3] = 0; }
    qd(long a){ x[0] = (double)a; x[1] = (double)(a - (long)x[0]); x[2] = x[3] = 0; }
    explicit operator double() const { return x[0] + (x[1] + (x[2] + x[3])); }
};

namespace qdimpl {
static inline void two_sum(double a, double b, double &s, double &e){ s = a + b; double bb = s - a; e = (a - (s - bb)) + (b - bb); }
static inline void fast_two_sum(double a, double b, double &s, double &e){ s = a + b; e = b - (s - a); }   /* |a| >= |b| */
/* exact sum of n doubles t[] -> qd.  The expansion e[] is kept non-overlapping in increasing magnitude
 * (grow-expansion with zero elimination), then compressed (Shewchuk 1997, sec. 2.8) and the top four
 * components taken; anything below them is added into the fourth. */
static inline qd sum(const double *t, int n){
    double e[64], f[64]; int m = 0;
    for(int k = 0; k < n; k++){
        double Q = t[k]; if(Q == 0) continue;
        int mm = 0;
        for(int i = 0; i < m; i++){ double s, h; two_sum(Q, e[i], s, h); if(h != 0) f[mm++] = h; Q = s; }
        f[mm++] = Q;
        memcpy(e, f, mm*sizeof(double)); m = mm;
    }
    qd r; if(m == 0) return r;
    /* compress: top-down then bottom-up pass */
    double g[64]; int bottom = m - 1; double Q = e[m-1];
    for(int i = m - 2; i >= 0; i--){ double s, q; fast_two_sum(Q, e[i], s, q); if(q != 0){ g[bottom--] = s; Q = q; } else Q = s; }
    g[bottom] = Q;
    double hh[64]; int top = 0; Q = g[bottom];
    for(int i = bottom + 1; i < m; i++){ double s, q; fast_two_sum(g[i], Q, s, q); if(q != 0) hh[top++] = q; Q = s; }
    hh[top++] = Q;                     /* hh[0..top-1], increasing magnitude, hh[top-1] ~ the sum */
    for(int j = 0; j < 4 && j < top; j++) r.x[j] = hh[top-1-j];
    if(top > 4){ double rest = 0; for(int i = 0; i < top - 4; i++) rest += hh[i]; r.x[3] += rest; }
    return r;
}
}

static inline qd operator+(const qd &a, const qd &b){ double t[8] = {a.x[3], b.x[3], a.x[2], b.x[2], a.x[1], b.x[1], a.x[0], b.x[0]}; return qdimpl::sum(t, 8); }
static inline qd operator-(const qd &a){ qd r; for(int i = 0; i < 4; i++) r.x[i] = -a.x[i]; return r; }
static inline qd operator-(const qd &a, const qd &b){ return a + (-b); }
static inline qd operator*(const qd &a, const qd &b){
    double t[26]; int n = 0;
    /* products a_i b_j with i+j <= 3 exactly (value + fma error), i+j = 4 rounded; the rest is below 2^-260 */
    for(int s = 4; s >= 0; s--) for(int i = 0; i <= s && i < 4; i++){ int j = s - i; if(j > 3) continue;
        double p = a.x[i]*b.x[j];
        if(s <= 3){ t[n++] = fma(a.x[i], b.x[j], -p); }
        t[n++] = p; }
    return qdimpl::sum(t, n);
}
static inline qd mul_d(const qd &a, double b){
    double t[8]; int n = 0;
    for(int i = 3; i >= 0; i--){ double p = a.x[i]*b; t[n++] = fma(a.x[i], b, -p); t[n++] = p; }
    return qdimpl::sum(t, n);
}
/* long division: five quotient digits, each from the leading double of the remainder */
static inline qd operator/(const qd &a, const qd &b){
    double q[5]; qd r = a;
    for(int k = 0; k < 5; k++){ q[k] = r.x[0]/b.x[0]; r = r - mul_d(b, q[k]); }
    double t[5] = {q[4], q[3], q[2], q[1], q[0]}; return qdimpl::sum(t, 5);
}
#define QDOPS(T) \
static inline qd operator+(T a,const qd&b){return qd(a)+b;} static inline qd operator+(const qd&a,T b){return a+qd(b);} \
static inline qd operator-(T a,const qd&b){return qd(a)-b;} static inline qd operator-(const qd&a,T b){return a-qd(b);} \
static inline qd operator*(T a,const qd&b){return mul_d(b,(double)a);} static inline qd operator*(const qd&a,T b){return mul_d(a,(double)b);} \
static inline qd operator/(T a,const qd&b){return qd(a)/b;} static inline qd operator/(const qd&a,T b){return a/qd(b);}
QDOPS(int) QDOPS(double)
static inline qd &operator+=(qd &a, const qd &b){ a = a + b; return a; } static inline qd &operator-=(qd &a, const qd &b){ a = a - b; return a; }
static inline qd &operator*=(qd &a, const qd &b){ a = a * b; return a; } static inline qd &operator/=(qd &a, const qd &b){ a = a / b; return a; }
static inline int qd_cmp(const qd &a, const qd &b){ qd d = a - b; return d.x[0] < 0 ? -1 : d.x[0] > 0 ? 1 : 0; }
static inline bool operator<(const qd &a, const qd &b){ return qd_cmp(a, b) < 0; }
static inline bool operator>(const qd &a, const qd &b){ return qd_cmp(a, b) > 0; }
static inline bool operator<=(const qd &a, const qd &b){ return qd_cmp(a, b) <= 0; }
static inline bool operator>=(const qd &a, const qd &b){ return qd_cmp(a, b) >= 0; }
/* comparisons with a plain number: the leading component decides unless it equals the number */
static inline bool operator<(const qd &a, int b){ return a.x[0] < b || (a.x[0] == b && a.x[1] < 0); }
static inline bool operator>(const qd &a, int b){ return a.x[0] > b || (a.x[0] == b && a.x[1] > 0); }
static inline bool operator<(const qd &a, double b){ return a.x[0] < b || (a.x[0] == b && a.x[1] < 0); }
static inline bool operator>(const qd &a, double b){ return a.x[0] > b || (a.x[0] == b && a.x[1] > 0); }

static inline qd sqrt(const qd &a){ if(a.x[0] <= 0) return qd(0); qd x = qd(sqrt(a.x[0])); for(int i = 0; i < 3; i++) x = (x + a/x)*0.5; return x; }
/* exp by argument halving: e^a = (e^{a/2^20})^{2^20}, Taylor series for the small argument */
static inline qd exp(const qd &a){ qd r = a*(1.0/1048576), t = qd(1), s = qd(1); for(int n = 1; n < 30; n++){ t = t*r/n; s += t; } for(int i = 0; i < 20; i++) s = s*s; return s; }

/* decimal string with nd significant digits, e.g. "-5.8997802228274214549e-01" */
static inline void qd_str(const qd &v, int nd, char *out){
    qd a = v; char *o = out;
    if(a.x[0] == 0){ strcpy(out, "+0"); return; }
    *o++ = a.x[0] < 0 ? '-' : '+'; if(a.x[0] < 0) a = -a;
    int e = (int)floor(log10(a.x[0]));
    qd p = qd(1), ten = qd(10); int ae = e < 0 ? -e : e;
    for(int i = 0; i < ae; i++) p = p*ten;
    a = e >= 0 ? a/p : a*p;
    if(a.x[0] >= 10){ a = a/ten; e++; } else if(a.x[0] < 1){ a = a*ten; e--; }
    int d[128]; if(nd > 120) nd = 120;
    for(int i = 0; i <= nd; i++){
        int k = (int)floor(a.x[0]); a = a - qd(k);
        if(a.x[0] < 0){ k--; a = a + qd(1); } else if(a.x[0] >= 1){ k++; a = a - qd(1); }
        d[i] = k; a = mul_d(a, 10.0);
    }
    if(d[nd] >= 5) d[nd-1]++;           /* round, then carry */
    for(int i = nd - 1; i > 0; i--) if(d[i] >= 10){ d[i] -= 10; d[i-1]++; } else if(d[i] < 0){ d[i] += 10; d[i-1]--; }
    if(d[0] >= 10){ d[0] -= 10; for(int i = nd - 1; i > 0; i--) d[i] = d[i-1]; d[1] = 0; d[0] = 1; e++; }
    *o++ = '0' + d[0]; *o++ = '.';
    for(int i = 1; i < nd; i++) *o++ = '0' + d[i];
    sprintf(o, "e%+03d", e);
}
/* the four components, round-trip exact (for exact post-processing, e.g. in mpmath) */
static inline void qd_parts(const qd &v, char *out){ sprintf(out, "%.17e %.17e %.17e %.17e", v.x[0], v.x[1], v.x[2], v.x[3]); }
