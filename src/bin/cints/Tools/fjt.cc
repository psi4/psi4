/*! \file fjt.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"fjt.h"


namespace psi { namespace cints {

/*!-------------------------
  Compute F_m, 0 <= m <= n
 -------------------------*/
void calc_f(double *F, int n, double t)
{
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  static double K = 1.0/M_2_SQRTPI;
  double et;


  if (t>20.0){
    t2 = 2*t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K*erf(t)/t;
    for(m=0; m<=n-1; m++){
      F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
      }
    }
  /*else {
    for(m=0; m<=n; m++){
      m2 = 2*m;
      num = df[m2];
      i = 0;
      sum = 1.0/(m2+1);
      do{
        i++;
        num = num*t2;
        term1 = num/df[m2+2*i+2];
        sum += term1;
        } while (fabs(term1) > EPS && i < MAXFACT);
      F[m] = sum*exp(-t);
      }
    }*/
  else {
    et = exp(-t);
    t2 = 2*t;
    m2 = 2*n;
    num = df[m2];
    i=0;
    sum = 1.0/(m2+1);
    do{
      i++;
      num = num*t2;
      term1 = num/df[m2+2*i+2];
      sum += term1;
      } while (fabs(term1) > EPS && i < MAXFACT);
    F[n] = sum*et;
    for(m=n-1;m>=0;m--){
      F[m] = (t2*F[m+1] + et)/(2*m+1);
      }
    }
}
}}
