#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_CY_p(prim_data *Data, const int ab_num, const int d_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_c;
  const double *i0, *i1;
  int ab,d;

  for(ab=0;ab<ab_num;ab++) {
  i0 = I0 + 1*d_num;
  for(d=0;d<d_num;d++)
    *(vp++) = twotzeta*(*(i0++)) ;
  i0 = I0 + 3*d_num;
  i1 = I1;
  for(d=0;d<d_num;d++)
    *(vp++) = twotzeta*(*(i0++)) - 1.000000*(*(i1++));
  i0 = I0 + 4*d_num;
  for(d=0;d<d_num;d++)
    *(vp++) = twotzeta*(*(i0++)) ;
  I0 += 6*d_num;  I1 += 1*d_num;
  }
}
