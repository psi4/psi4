  /* These machine-generated functions compute a quartet of |pp) integrals */

#include "libint.h"

void hrr3_build_pp(const REALTYPE *CD, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1, int ab_num)
{
  const REALTYPE CD0 = CD[0];
  const REALTYPE CD1 = CD[1];
  const REALTYPE CD2 = CD[2];
  int ab;

  for(ab=0;ab<ab_num;ab++) {
    *(vp++) = I0[0] + CD0*I1[0];
    *(vp++) = I0[1] + CD1*I1[0];
    *(vp++) = I0[2] + CD2*I1[0];
    *(vp++) = I0[1] + CD0*I1[1];
    *(vp++) = I0[3] + CD1*I1[1];
    *(vp++) = I0[4] + CD2*I1[1];
    *(vp++) = I0[2] + CD0*I1[2];
    *(vp++) = I0[4] + CD1*I1[2];
    *(vp++) = I0[5] + CD2*I1[2];
    I0 += 6;
    I1 += 3;
  }
}
/* Total number of FLOPs = 18 * ab_num */
