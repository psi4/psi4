  /* These machine-generated functions compute a quartet of |0p) integrals */

#include "libint.h"

void hrr3_build_0p(const REALTYPE *CD, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1, int ab_num)
{
  const REALTYPE CD0 = CD[0];
  const REALTYPE CD1 = CD[1];
  const REALTYPE CD2 = CD[2];
  int ab;

  for(ab=0;ab<ab_num;ab++) {
    *(vp++) = I0[0] + CD0*I1[0];
    *(vp++) = I0[1] + CD1*I1[0];
    *(vp++) = I0[2] + CD2*I1[0];
    I0 += 3;
    I1 += 1;
  }
}
/* Total number of FLOPs = 6 * ab_num */
