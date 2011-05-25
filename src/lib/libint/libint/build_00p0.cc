  /* These machine-generated functions compute a quartet of (0s|ps) integrals */

#include "libint.h"

void _build_00p0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
{
  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;
  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;
  REALTYPE lpoz = Data->poz;
  REALTYPE lpon = Data->pon;
  U20 = Data->U[2][0];
  U21 = Data->U[2][1];
  U22 = Data->U[2][2];
  U50 = Data->U[5][0];
  U51 = Data->U[5][1];
  U52 = Data->U[5][2];


*(vp++) = U20*I0[0] + U50*I1[0];
*(vp++) = U21*I0[0] + U51*I1[0];
*(vp++) = U22*I0[0] + U52*I1[0];

}
/* Total number of FLOPs = 9 */
