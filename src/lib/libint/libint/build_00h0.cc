  /* These machine-generated functions compute a quartet of (0s|hs) integrals */

#include "libint.h"

void _build_00h0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
{
  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;
  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;
  REALTYPE lpoz = Data->poz;
  REALTYPE lpon = Data->pon;
  REALTYPE oneo2n;
  REALTYPE twoo2n;
  REALTYPE threeo2n;
  REALTYPE fouro2n;
  oneo2n = 1.0*Data->oo2n;
  twoo2n = 2.0*Data->oo2n;
  threeo2n = 3.0*Data->oo2n;
  fouro2n = 4.0*Data->oo2n;
  U20 = Data->U[2][0];
  U21 = Data->U[2][1];
  U22 = Data->U[2][2];
  U50 = Data->U[5][0];
  U51 = Data->U[5][1];
  U52 = Data->U[5][2];


*(vp++) = U20*I0[0] + U50*I1[0]
           + (fouro2n)*(I2[0] - (lpon)*I3[0]);
*(vp++) = U20*I0[1] + U50*I1[1]
           + (threeo2n)*(I2[1] - (lpon)*I3[1]);
*(vp++) = U20*I0[2] + U50*I1[2]
           + (threeo2n)*(I2[2] - (lpon)*I3[2]);
*(vp++) = U20*I0[3] + U50*I1[3]
           + (twoo2n)*(I2[3] - (lpon)*I3[3]);
*(vp++) = U20*I0[4] + U50*I1[4]
           + (twoo2n)*(I2[4] - (lpon)*I3[4]);
*(vp++) = U20*I0[5] + U50*I1[5]
           + (twoo2n)*(I2[5] - (lpon)*I3[5]);
*(vp++) = U20*I0[6] + U50*I1[6]
           + (oneo2n)*(I2[6] - (lpon)*I3[6]);
*(vp++) = U20*I0[7] + U50*I1[7]
           + (oneo2n)*(I2[7] - (lpon)*I3[7]);
*(vp++) = U20*I0[8] + U50*I1[8]
           + (oneo2n)*(I2[8] - (lpon)*I3[8]);
*(vp++) = U20*I0[9] + U50*I1[9]
           + (oneo2n)*(I2[9] - (lpon)*I3[9]);
*(vp++) = U20*I0[10] + U50*I1[10];
*(vp++) = U20*I0[11] + U50*I1[11];
*(vp++) = U20*I0[12] + U50*I1[12];
*(vp++) = U20*I0[13] + U50*I1[13];
*(vp++) = U20*I0[14] + U50*I1[14];
*(vp++) = U21*I0[10] + U51*I1[10]
           + (fouro2n)*(I2[6] - (lpon)*I3[6]);
*(vp++) = U21*I0[11] + U51*I1[11]
           + (threeo2n)*(I2[7] - (lpon)*I3[7]);
*(vp++) = U21*I0[12] + U51*I1[12]
           + (twoo2n)*(I2[8] - (lpon)*I3[8]);
*(vp++) = U21*I0[13] + U51*I1[13]
           + (oneo2n)*(I2[9] - (lpon)*I3[9]);
*(vp++) = U21*I0[14] + U51*I1[14];
*(vp++) = U22*I0[14] + U52*I1[14]
           + (fouro2n)*(I2[9] - (lpon)*I3[9]);

}
/* Total number of FLOPs = 123 */
