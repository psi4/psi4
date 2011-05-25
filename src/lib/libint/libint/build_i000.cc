  /* These machine-generated functions compute a quartet of (is|0s) integrals */

#include "libint.h"

void _build_i000(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
{
  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;
  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;
  REALTYPE lpoz = Data->poz;
  REALTYPE lpon = Data->pon;
  REALTYPE oneo2z;
  REALTYPE twoo2z;
  REALTYPE threeo2z;
  REALTYPE fouro2z;
  REALTYPE fiveo2z;
  oneo2z = 1.0*Data->oo2z;
  twoo2z = 2.0*Data->oo2z;
  threeo2z = 3.0*Data->oo2z;
  fouro2z = 4.0*Data->oo2z;
  fiveo2z = 5.0*Data->oo2z;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (fiveo2z)*(I2[0] - (lpoz)*I3[0]);
*(vp++) = U00*I0[1] + U40*I1[1]
           + (fouro2z)*(I2[1] - (lpoz)*I3[1]);
*(vp++) = U00*I0[2] + U40*I1[2]
           + (fouro2z)*(I2[2] - (lpoz)*I3[2]);
*(vp++) = U00*I0[3] + U40*I1[3]
           + (threeo2z)*(I2[3] - (lpoz)*I3[3]);
*(vp++) = U00*I0[4] + U40*I1[4]
           + (threeo2z)*(I2[4] - (lpoz)*I3[4]);
*(vp++) = U00*I0[5] + U40*I1[5]
           + (threeo2z)*(I2[5] - (lpoz)*I3[5]);
*(vp++) = U00*I0[6] + U40*I1[6]
           + (twoo2z)*(I2[6] - (lpoz)*I3[6]);
*(vp++) = U00*I0[7] + U40*I1[7]
           + (twoo2z)*(I2[7] - (lpoz)*I3[7]);
*(vp++) = U00*I0[8] + U40*I1[8]
           + (twoo2z)*(I2[8] - (lpoz)*I3[8]);
*(vp++) = U00*I0[9] + U40*I1[9]
           + (twoo2z)*(I2[9] - (lpoz)*I3[9]);
*(vp++) = U00*I0[10] + U40*I1[10]
           + (oneo2z)*(I2[10] - (lpoz)*I3[10]);
*(vp++) = U00*I0[11] + U40*I1[11]
           + (oneo2z)*(I2[11] - (lpoz)*I3[11]);
*(vp++) = U00*I0[12] + U40*I1[12]
           + (oneo2z)*(I2[12] - (lpoz)*I3[12]);
*(vp++) = U00*I0[13] + U40*I1[13]
           + (oneo2z)*(I2[13] - (lpoz)*I3[13]);
*(vp++) = U00*I0[14] + U40*I1[14]
           + (oneo2z)*(I2[14] - (lpoz)*I3[14]);
*(vp++) = U00*I0[15] + U40*I1[15];
*(vp++) = U00*I0[16] + U40*I1[16];
*(vp++) = U00*I0[17] + U40*I1[17];
*(vp++) = U00*I0[18] + U40*I1[18];
*(vp++) = U00*I0[19] + U40*I1[19];
*(vp++) = U00*I0[20] + U40*I1[20];
*(vp++) = U01*I0[15] + U41*I1[15]
           + (fiveo2z)*(I2[10] - (lpoz)*I3[10]);
*(vp++) = U01*I0[16] + U41*I1[16]
           + (fouro2z)*(I2[11] - (lpoz)*I3[11]);
*(vp++) = U01*I0[17] + U41*I1[17]
           + (threeo2z)*(I2[12] - (lpoz)*I3[12]);
*(vp++) = U01*I0[18] + U41*I1[18]
           + (twoo2z)*(I2[13] - (lpoz)*I3[13]);
*(vp++) = U01*I0[19] + U41*I1[19]
           + (oneo2z)*(I2[14] - (lpoz)*I3[14]);
*(vp++) = U01*I0[20] + U41*I1[20];
*(vp++) = U02*I0[20] + U42*I1[20]
           + (fiveo2z)*(I2[14] - (lpoz)*I3[14]);

}
/* Total number of FLOPs = 168 */
