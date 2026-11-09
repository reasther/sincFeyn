/* double-quad: hi + lo with __float128 parts (~68 digits). Same algorithms as dd.h. */
#pragma once
#include <quadmath.h>
typedef __float128 F;
struct qq { F hi, lo; qq():hi(0),lo(0){} qq(F x):hi(x),lo(0){} qq(double x):hi(x),lo(0){} qq(int x):hi(x),lo(0){} qq(F h,F l):hi(h),lo(l){} };
static inline qq two_sum(F a, F b){ F s=a+b, bb=s-a; return qq(s,(a-(s-bb))+(b-bb)); }
static inline qq qtwo_sum(F a, F b){ F s=a+b; return qq(s,b-(s-a)); }
static inline qq two_prod(F a, F b){ F p=a*b; return qq(p,fmaq(a,b,-p)); }
static inline qq operator+(const qq&a,const qq&b){ qq s=two_sum(a.hi,b.hi), t=two_sum(a.lo,b.lo); s.lo+=t.hi; s=qtwo_sum(s.hi,s.lo); s.lo+=t.lo; return qtwo_sum(s.hi,s.lo); }
static inline qq operator-(const qq&a){ return qq(-a.hi,-a.lo); }
static inline qq operator-(const qq&a,const qq&b){ return a+(-b); }
static inline qq operator*(const qq&a,const qq&b){ qq p=two_prod(a.hi,b.hi); p.lo+=a.hi*b.lo+a.lo*b.hi; return qtwo_sum(p.hi,p.lo); }
static inline qq operator/(const qq&a,const qq&b){ F q1=a.hi/b.hi; qq r=a-b*qq(q1); F q2=r.hi/b.hi; r=r-b*qq(q2); F q3=r.hi/b.hi; return qtwo_sum(q1,q2)+qq(q3); }
#define QOPS(T) \
static inline qq operator+(T a,const qq&b){return qq(a)+b;} static inline qq operator+(const qq&a,T b){return a+qq(b);} \
static inline qq operator-(T a,const qq&b){return qq(a)-b;} static inline qq operator-(const qq&a,T b){return a-qq(b);} \
static inline qq operator*(T a,const qq&b){return qq(a)*b;} static inline qq operator*(const qq&a,T b){return a*qq(b);} \
static inline qq operator/(T a,const qq&b){return qq(a)/b;} static inline qq operator/(const qq&a,T b){return a/qq(b);}
QOPS(int) QOPS(double)
static inline qq&operator+=(qq&a,const qq&b){a=a+b;return a;} static inline qq&operator-=(qq&a,const qq&b){a=a-b;return a;}
static inline qq&operator*=(qq&a,const qq&b){a=a*b;return a;} static inline qq&operator/=(qq&a,const qq&b){a=a/b;return a;}
static inline bool operator<(const qq&a,const qq&b){return a.hi<b.hi||(a.hi==b.hi&&a.lo<b.lo);}
static inline bool operator>(const qq&a,const qq&b){return b<a;}
static inline bool operator<=(const qq&a,const qq&b){return !(b<a);} static inline bool operator>=(const qq&a,const qq&b){return !(a<b);}
static inline bool operator<(const qq&a,int b){return a<qq(b);} static inline bool operator>(const qq&a,int b){return a>qq(b);}
static inline qq exp(const qq&a){ qq r=a/qq((F)1048576), t=1, s=1; for(int n=1;n<40;n++){ t=t*r/n; s+=t; } for(int i=0;i<20;i++) s=s*s; return s; }
static inline qq sqrt(const qq&a){ qq x=sqrtq(a.hi); for(int i=0;i<3;i++) x=(x+a/x)/2; return x; }
