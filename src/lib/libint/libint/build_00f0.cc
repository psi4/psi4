  /* These machine-generated functions compute a quartet of (0s|fs) integrals */

#include "libint.h"

void _build_00f0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
{
  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;
  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;
  REALTYPE lpoz = Data->poz;
  REALTYPE lpon = Data->pon;
  REALTYPE oneo2n;
  REALTYPE twoo2n;
  oneo2n = 1.0*Data->oo2n;
  twoo2n = 2.0*Data->oo2n;
  U20 = Data->U[2][0];
  U21 = Data->U[2][1];
  U22 = Data->U[2][2];
  U50 = Data->U[5][0];
  U51 = Data->U[5][1];
  U52 = Data->U[5][2];


*(vp++) = U20*I0[0] + U50*I1[0]
           + (twoo2n)*(I2[0] - (lpon)*I3[0]);
*(vp++) = U20*I0[1] + U50*I1[1]
           + (oneo2n)*(I2[1] - (lpon)*I3[1]);
*(vp++) = U20*I0[2] + U50*I1[2]
           + (oneo2n)*(I2[2] - (lpon)*I3[2]);
*(vp++) = U20*I0[3] + U50*I1[3];
*(vp++) = U20*I0[4] + U50*I1[4];
*(vp++) = U20*I0[5] + U50*I1[5];
*(vp++) = U21*I0[3] + U51*I1[3]
           + (twoo2n)*(I2[1] - (lpon)*I3[1]);
*(vp++) = U21*I0[4] + U51*I1[4]
           + (oneo2n)*(I2[2] - (lpon)*I3[2]);
*(vp++) = U21*I0[5] + U51*I1[5];
*(vp++) = U22*I0[5] + U52*I1[5]
           + (twoo2n)*(I2[2] - (lpon)*I3[2]);

}
/* Total number of FLOPs = 54 */
