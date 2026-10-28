/* DeltaGL terms and full-box sums |k_i| <= K of the C kernels kernels/g2/*.c (as sincqed sums them, without pruning).
   From the repo top:  gcc -O2 -Idriver -I. -DQREAL_DOUBLE -o boxsum_gl <this file> -lm;  ./boxsum_gl 1.0 35 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sq.h"
__thread int SQ_AB=1;
#include "kernels/g2/g2_00.c"
#include "kernels/g2/g2_01.c"
#include "kernels/g2/g2_02.c"
#include "kernels/g2/g2_03.c"
#include "kernels/g2/v2_00.c"
typedef double (*F)(const double*, double*);
static double term(F f, int D, const int *k, double h){ double s[8], z=0, ab; for(int i=0;i<D;i++){ s[i]=exp(h*k[i]); z+=h*k[i]; } return pow(h,D)*exp(z)*f(s,&ab); }
static double box(F f, int D, double h, int K){ int k[8]; long n=1; for(int i=0;i<D;i++) n*=2*K+1; double tot=0, c=0;
  for(long j=0;j<n;j++){ long r=j; for(int i=0;i<D;i++){ k[i]=r%(2*K+1)-K; r/=2*K+1; } double y=term(f,D,k,h)-c, t=tot+y; c=(t-tot)-y; tot=t; } return tot; }
int main(int argc, char **argv){
  F fs[5]={f_g2_00_d,f_g2_01_d,f_g2_02_d,f_g2_03_d,f_v2_00_d}; int Ds[5]={4,3,4,3,3}; const char *nm[5]={"g2_00","g2_01","g2_02","g2_03","v2_00"};
  int z[4]={0,0,0,0}, k1[4]={1,-1,2,0};
  for(int g=0;g<5;g++) printf("%s term0 %.17g term1 %.17g\n", nm[g], term(fs[g],Ds[g],z,1.0), term(fs[g],Ds[g],k1,1.0));
  if(argc>1) for(int g=0;g<5;g++) printf("%s box h=%s K=%s %.10f\n", nm[g], argv[1], argv[2], box(fs[g],Ds[g],atof(argv[1]),atoi(argv[2])));
}
