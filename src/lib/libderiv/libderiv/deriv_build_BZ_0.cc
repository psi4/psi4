#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_BZ_0(prim_data *Data, const int a_num, const int cd_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_b;
  const double *i0, *i1;
  int a,cd;

  for(a=0;a<a_num;a++) {
  i0 = I0 + 2*cd_num;
  for(cd=0;cd<cd_num;cd++)
    *(vp++) = twotzeta*(*(i0++)) ;
  I0 += 3*cd_num;  I1 += 0*cd_num;
  }
}
