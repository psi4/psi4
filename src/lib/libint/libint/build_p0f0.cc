  /* These machine-generated functions compute a quartet of (ps|fs) integrals */

#include "libint.h"

void _build_p0f0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
{
  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;
  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;
  REALTYPE lpoz = Data->poz;
  REALTYPE lpon = Data->pon;
  REALTYPE oneo2zn;
  REALTYPE twoo2zn;
  REALTYPE threeo2zn;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  threeo2zn = 3.0*Data->oo2zn;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (threeo2zn)*I4[0];
*(vp++) = U00*I0[1] + U40*I1[1]
           + (twoo2zn)*I4[1];
*(vp++) = U00*I0[2] + U40*I1[2]
           + (twoo2zn)*I4[2];
*(vp++) = U00*I0[3] + U40*I1[3]
           + (oneo2zn)*I4[3];
*(vp++) = U00*I0[4] + U40*I1[4]
           + (oneo2zn)*I4[4];
*(vp++) = U00*I0[5] + U40*I1[5]
           + (oneo2zn)*I4[5];
*(vp++) = U00*I0[6] + U40*I1[6];
*(vp++) = U00*I0[7] + U40*I1[7];
*(vp++) = U00*I0[8] + U40*I1[8];
*(vp++) = U00*I0[9] + U40*I1[9];
*(vp++) = U01*I0[0] + U41*I1[0];
*(vp++) = U01*I0[1] + U41*I1[1]
           + (oneo2zn)*I4[0];
*(vp++) = U01*I0[2] + U41*I1[2];
*(vp++) = U01*I0[3] + U41*I1[3]
           + (twoo2zn)*I4[1];
*(vp++) = U01*I0[4] + U41*I1[4]
           + (oneo2zn)*I4[2];
*(vp++) = U01*I0[5] + U41*I1[5];
*(vp++) = U01*I0[6] + U41*I1[6]
           + (threeo2zn)*I4[3];
*(vp++) = U01*I0[7] + U41*I1[7]
           + (twoo2zn)*I4[4];
*(vp++) = U01*I0[8] + U41*I1[8]
           + (oneo2zn)*I4[5];
*(vp++) = U01*I0[9] + U41*I1[9];
*(vp++) = U02*I0[0] + U42*I1[0];
*(vp++) = U02*I0[1] + U42*I1[1];
*(vp++) = U02*I0[2] + U42*I1[2]
           + (oneo2zn)*I4[0];
*(vp++) = U02*I0[3] + U42*I1[3];
*(vp++) = U02*I0[4] + U42*I1[4]
           + (oneo2zn)*I4[1];
*(vp++) = U02*I0[5] + U42*I1[5]
           + (twoo2zn)*I4[2];
*(vp++) = U02*I0[6] + U42*I1[6];
*(vp++) = U02*I0[7] + U42*I1[7]
           + (oneo2zn)*I4[3];
*(vp++) = U02*I0[8] + U42*I1[8]
           + (twoo2zn)*I4[4];
*(vp++) = U02*I0[9] + U42*I1[9]
           + (threeo2zn)*I4[5];

}
/* Total number of FLOPs = 126 */
