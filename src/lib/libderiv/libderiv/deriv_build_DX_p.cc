#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DX_p(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[0] - 1.000000*I1[0];
    *(vp++) = twotzeta*I0[1] ;
    *(vp++) = twotzeta*I0[2] ;
  I0 += 6;  I1 += 1;
  }
}
