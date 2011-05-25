  /* These machine-generated functions compute a quartet of |dd) integrals */

#include "libint.h"

void hrr3_build_dd(const REALTYPE *CD, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1, int ab_num)
{
  const REALTYPE CD0 = CD[0];
  const REALTYPE CD1 = CD[1];
  const REALTYPE CD2 = CD[2];
  int ab;

  for(ab=0;ab<ab_num;ab++) {
    *(vp++) = I0[0] + CD0*I1[0];
    *(vp++) = I0[1] + CD0*I1[1];
    *(vp++) = I0[2] + CD0*I1[2];
    *(vp++) = I0[4] + CD1*I1[1];
    *(vp++) = I0[5] + CD1*I1[2];
    *(vp++) = I0[8] + CD2*I1[2];
    *(vp++) = I0[3] + CD0*I1[3];
    *(vp++) = I0[4] + CD0*I1[4];
    *(vp++) = I0[5] + CD0*I1[5];
    *(vp++) = I0[10] + CD1*I1[4];
    *(vp++) = I0[11] + CD1*I1[5];
    *(vp++) = I0[14] + CD2*I1[5];
    *(vp++) = I0[6] + CD0*I1[6];
    *(vp++) = I0[7] + CD0*I1[7];
    *(vp++) = I0[8] + CD0*I1[8];
    *(vp++) = I0[13] + CD1*I1[7];
    *(vp++) = I0[14] + CD1*I1[8];
    *(vp++) = I0[17] + CD2*I1[8];
    *(vp++) = I0[9] + CD0*I1[9];
    *(vp++) = I0[10] + CD0*I1[10];
    *(vp++) = I0[11] + CD0*I1[11];
    *(vp++) = I0[19] + CD1*I1[10];
    *(vp++) = I0[20] + CD1*I1[11];
    *(vp++) = I0[23] + CD2*I1[11];
    *(vp++) = I0[12] + CD0*I1[12];
    *(vp++) = I0[13] + CD0*I1[13];
    *(vp++) = I0[14] + CD0*I1[14];
    *(vp++) = I0[22] + CD1*I1[13];
    *(vp++) = I0[23] + CD1*I1[14];
    *(vp++) = I0[26] + CD2*I1[14];
    *(vp++) = I0[15] + CD0*I1[15];
    *(vp++) = I0[16] + CD0*I1[16];
    *(vp++) = I0[17] + CD0*I1[17];
    *(vp++) = I0[25] + CD1*I1[16];
    *(vp++) = I0[26] + CD1*I1[17];
    *(vp++) = I0[29] + CD2*I1[17];
    I0 += 30;
    I1 += 18;
  }
}
/* Total number of FLOPs = 72 * ab_num */
