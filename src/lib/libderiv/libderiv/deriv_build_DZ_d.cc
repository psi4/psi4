#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DZ_d(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[2] ;
    *(vp++) = twotzeta*I0[4] ;
    *(vp++) = twotzeta*I0[5] - 1.000000*I1[0];
    *(vp++) = twotzeta*I0[7] ;
    *(vp++) = twotzeta*I0[8] - 1.000000*I1[1];
    *(vp++) = twotzeta*I0[9] - 2.000000*I1[2];
  I0 += 10;  I1 += 3;
  }
}
