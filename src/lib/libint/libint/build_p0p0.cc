  /* These machine-generated functions compute a quartet of (ps|ps) integrals */

#include "libint.h"

void _build_p0p0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
{
  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;
  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;
  REALTYPE lpoz = Data->poz;
  REALTYPE lpon = Data->pon;
  REALTYPE oneo2zn;
  oneo2zn = 1.0*Data->oo2zn;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (oneo2zn)*I4[0];
*(vp++) = U00*I0[1] + U40*I1[1];
*(vp++) = U00*I0[2] + U40*I1[2];
*(vp++) = U01*I0[0] + U41*I1[0];
*(vp++) = U01*I0[1] + U41*I1[1]
           + (oneo2zn)*I4[0];
*(vp++) = U01*I0[2] + U41*I1[2];
*(vp++) = U02*I0[0] + U42*I1[0];
*(vp++) = U02*I0[1] + U42*I1[1];
*(vp++) = U02*I0[2] + U42*I1[2]
           + (oneo2zn)*I4[0];

}
/* Total number of FLOPs = 33 */
