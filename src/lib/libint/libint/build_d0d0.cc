  /* These machine-generated functions compute a quartet of (ds|ds) integrals */

#include "libint.h"

void _build_d0d0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
{
  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;
  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;
  REALTYPE lpoz = Data->poz;
  REALTYPE lpon = Data->pon;
  REALTYPE oneo2zn;
  REALTYPE twoo2zn;
  REALTYPE oneo2z;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  oneo2z = 1.0*Data->oo2z;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (oneo2z)*(I2[0] - (lpoz)*I3[0])
           + (twoo2zn)*I4[0];
*(vp++) = U00*I0[1] + U40*I1[1]
           + (oneo2z)*(I2[1] - (lpoz)*I3[1])
           + (oneo2zn)*I4[1];
*(vp++) = U00*I0[2] + U40*I1[2]
           + (oneo2z)*(I2[2] - (lpoz)*I3[2])
           + (oneo2zn)*I4[2];
*(vp++) = U00*I0[3] + U40*I1[3]
           + (oneo2z)*(I2[3] - (lpoz)*I3[3]);
*(vp++) = U00*I0[4] + U40*I1[4]
           + (oneo2z)*(I2[4] - (lpoz)*I3[4]);
*(vp++) = U00*I0[5] + U40*I1[5]
           + (oneo2z)*(I2[5] - (lpoz)*I3[5]);
*(vp++) = U00*I0[6] + U40*I1[6]
           + (twoo2zn)*I4[3];
*(vp++) = U00*I0[7] + U40*I1[7]
           + (oneo2zn)*I4[4];
*(vp++) = U00*I0[8] + U40*I1[8]
           + (oneo2zn)*I4[5];
*(vp++) = U00*I0[9] + U40*I1[9];
*(vp++) = U00*I0[10] + U40*I1[10];
*(vp++) = U00*I0[11] + U40*I1[11];
*(vp++) = U00*I0[12] + U40*I1[12]
           + (twoo2zn)*I4[6];
*(vp++) = U00*I0[13] + U40*I1[13]
           + (oneo2zn)*I4[7];
*(vp++) = U00*I0[14] + U40*I1[14]
           + (oneo2zn)*I4[8];
*(vp++) = U00*I0[15] + U40*I1[15];
*(vp++) = U00*I0[16] + U40*I1[16];
*(vp++) = U00*I0[17] + U40*I1[17];
*(vp++) = U01*I0[6] + U41*I1[6]
           + (oneo2z)*(I2[0] - (lpoz)*I3[0]);
*(vp++) = U01*I0[7] + U41*I1[7]
           + (oneo2z)*(I2[1] - (lpoz)*I3[1])
           + (oneo2zn)*I4[3];
*(vp++) = U01*I0[8] + U41*I1[8]
           + (oneo2z)*(I2[2] - (lpoz)*I3[2]);
*(vp++) = U01*I0[9] + U41*I1[9]
           + (oneo2z)*(I2[3] - (lpoz)*I3[3])
           + (twoo2zn)*I4[4];
*(vp++) = U01*I0[10] + U41*I1[10]
           + (oneo2z)*(I2[4] - (lpoz)*I3[4])
           + (oneo2zn)*I4[5];
*(vp++) = U01*I0[11] + U41*I1[11]
           + (oneo2z)*(I2[5] - (lpoz)*I3[5]);
*(vp++) = U01*I0[12] + U41*I1[12];
*(vp++) = U01*I0[13] + U41*I1[13]
           + (oneo2zn)*I4[6];
*(vp++) = U01*I0[14] + U41*I1[14];
*(vp++) = U01*I0[15] + U41*I1[15]
           + (twoo2zn)*I4[7];
*(vp++) = U01*I0[16] + U41*I1[16]
           + (oneo2zn)*I4[8];
*(vp++) = U01*I0[17] + U41*I1[17];
*(vp++) = U02*I0[12] + U42*I1[12]
           + (oneo2z)*(I2[0] - (lpoz)*I3[0]);
*(vp++) = U02*I0[13] + U42*I1[13]
           + (oneo2z)*(I2[1] - (lpoz)*I3[1]);
*(vp++) = U02*I0[14] + U42*I1[14]
           + (oneo2z)*(I2[2] - (lpoz)*I3[2])
           + (oneo2zn)*I4[6];
*(vp++) = U02*I0[15] + U42*I1[15]
           + (oneo2z)*(I2[3] - (lpoz)*I3[3]);
*(vp++) = U02*I0[16] + U42*I1[16]
           + (oneo2z)*(I2[4] - (lpoz)*I3[4])
           + (oneo2zn)*I4[7];
*(vp++) = U02*I0[17] + U42*I1[17]
           + (oneo2z)*(I2[5] - (lpoz)*I3[5])
           + (twoo2zn)*I4[8];

}
/* Total number of FLOPs = 216 */
