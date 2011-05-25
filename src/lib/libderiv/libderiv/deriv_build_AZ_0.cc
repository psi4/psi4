#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_AZ_0(prim_data *Data, const int bcd_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_a;
  const double *i0, *i1;
  int bcd;

  i0 = I0 + 2*bcd_num;
  for(bcd=0;bcd<bcd_num;bcd++)
    *(vp++) = twotzeta*(*(i0++)) ;
}
