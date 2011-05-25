#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DX_f(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[0] - 3.000000*I1[0];
    *(vp++) = twotzeta*I0[1] - 2.000000*I1[1];
    *(vp++) = twotzeta*I0[2] - 2.000000*I1[2];
    *(vp++) = twotzeta*I0[3] - 1.000000*I1[3];
    *(vp++) = twotzeta*I0[4] - 1.000000*I1[4];
    *(vp++) = twotzeta*I0[5] - 1.000000*I1[5];
    *(vp++) = twotzeta*I0[6] ;
    *(vp++) = twotzeta*I0[7] ;
    *(vp++) = twotzeta*I0[8] ;
    *(vp++) = twotzeta*I0[9] ;
  I0 += 15;  I1 += 6;
  }
}
