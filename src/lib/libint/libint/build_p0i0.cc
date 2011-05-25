  /* These machine-generated functions compute a quartet of (ps|is) integrals */

#include "libint.h"

void _build_p0i0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
{
  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;
  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;
  REALTYPE lpoz = Data->poz;
  REALTYPE lpon = Data->pon;
  REALTYPE oneo2zn;
  REALTYPE twoo2zn;
  REALTYPE threeo2zn;
  REALTYPE fouro2zn;
  REALTYPE fiveo2zn;
  REALTYPE sixo2zn;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  threeo2zn = 3.0*Data->oo2zn;
  fouro2zn = 4.0*Data->oo2zn;
  fiveo2zn = 5.0*Data->oo2zn;
  sixo2zn = 6.0*Data->oo2zn;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (sixo2zn)*I4[0];
*(vp++) = U00*I0[1] + U40*I1[1]
           + (fiveo2zn)*I4[1];
*(vp++) = U00*I0[2] + U40*I1[2]
           + (fiveo2zn)*I4[2];
*(vp++) = U00*I0[3] + U40*I1[3]
           + (fouro2zn)*I4[3];
*(vp++) = U00*I0[4] + U40*I1[4]
           + (fouro2zn)*I4[4];
*(vp++) = U00*I0[5] + U40*I1[5]
           + (fouro2zn)*I4[5];
*(vp++) = U00*I0[6] + U40*I1[6]
           + (threeo2zn)*I4[6];
*(vp++) = U00*I0[7] + U40*I1[7]
           + (threeo2zn)*I4[7];
*(vp++) = U00*I0[8] + U40*I1[8]
           + (threeo2zn)*I4[8];
*(vp++) = U00*I0[9] + U40*I1[9]
           + (threeo2zn)*I4[9];
*(vp++) = U00*I0[10] + U40*I1[10]
           + (twoo2zn)*I4[10];
*(vp++) = U00*I0[11] + U40*I1[11]
           + (twoo2zn)*I4[11];
*(vp++) = U00*I0[12] + U40*I1[12]
           + (twoo2zn)*I4[12];
*(vp++) = U00*I0[13] + U40*I1[13]
           + (twoo2zn)*I4[13];
*(vp++) = U00*I0[14] + U40*I1[14]
           + (twoo2zn)*I4[14];
*(vp++) = U00*I0[15] + U40*I1[15]
           + (oneo2zn)*I4[15];
*(vp++) = U00*I0[16] + U40*I1[16]
           + (oneo2zn)*I4[16];
*(vp++) = U00*I0[17] + U40*I1[17]
           + (oneo2zn)*I4[17];
*(vp++) = U00*I0[18] + U40*I1[18]
           + (oneo2zn)*I4[18];
*(vp++) = U00*I0[19] + U40*I1[19]
           + (oneo2zn)*I4[19];
*(vp++) = U00*I0[20] + U40*I1[20]
           + (oneo2zn)*I4[20];
*(vp++) = U00*I0[21] + U40*I1[21];
*(vp++) = U00*I0[22] + U40*I1[22];
*(vp++) = U00*I0[23] + U40*I1[23];
*(vp++) = U00*I0[24] + U40*I1[24];
*(vp++) = U00*I0[25] + U40*I1[25];
*(vp++) = U00*I0[26] + U40*I1[26];
*(vp++) = U00*I0[27] + U40*I1[27];
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
*(vp++) = U01*I0[10] + U41*I1[10]
           + (fouro2zn)*I4[6];
*(vp++) = U01*I0[11] + U41*I1[11]
           + (threeo2zn)*I4[7];
*(vp++) = U01*I0[12] + U41*I1[12]
           + (twoo2zn)*I4[8];
*(vp++) = U01*I0[13] + U41*I1[13]
           + (oneo2zn)*I4[9];
*(vp++) = U01*I0[14] + U41*I1[14];
*(vp++) = U01*I0[15] + U41*I1[15]
           + (fiveo2zn)*I4[10];
*(vp++) = U01*I0[16] + U41*I1[16]
           + (fouro2zn)*I4[11];
*(vp++) = U01*I0[17] + U41*I1[17]
           + (threeo2zn)*I4[12];
*(vp++) = U01*I0[18] + U41*I1[18]
           + (twoo2zn)*I4[13];
*(vp++) = U01*I0[19] + U41*I1[19]
           + (oneo2zn)*I4[14];
*(vp++) = U01*I0[20] + U41*I1[20];
*(vp++) = U01*I0[21] + U41*I1[21]
           + (sixo2zn)*I4[15];
*(vp++) = U01*I0[22] + U41*I1[22]
           + (fiveo2zn)*I4[16];
*(vp++) = U01*I0[23] + U41*I1[23]
           + (fouro2zn)*I4[17];
*(vp++) = U01*I0[24] + U41*I1[24]
           + (threeo2zn)*I4[18];
*(vp++) = U01*I0[25] + U41*I1[25]
           + (twoo2zn)*I4[19];
*(vp++) = U01*I0[26] + U41*I1[26]
           + (oneo2zn)*I4[20];
*(vp++) = U01*I0[27] + U41*I1[27];
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
*(vp++) = U02*I0[10] + U42*I1[10];
*(vp++) = U02*I0[11] + U42*I1[11]
           + (oneo2zn)*I4[6];
*(vp++) = U02*I0[12] + U42*I1[12]
           + (twoo2zn)*I4[7];
*(vp++) = U02*I0[13] + U42*I1[13]
           + (threeo2zn)*I4[8];
*(vp++) = U02*I0[14] + U42*I1[14]
           + (fouro2zn)*I4[9];
*(vp++) = U02*I0[15] + U42*I1[15];
*(vp++) = U02*I0[16] + U42*I1[16]
           + (oneo2zn)*I4[10];
*(vp++) = U02*I0[17] + U42*I1[17]
           + (twoo2zn)*I4[11];
*(vp++) = U02*I0[18] + U42*I1[18]
           + (threeo2zn)*I4[12];
*(vp++) = U02*I0[19] + U42*I1[19]
           + (fouro2zn)*I4[13];
*(vp++) = U02*I0[20] + U42*I1[20]
           + (fiveo2zn)*I4[14];
*(vp++) = U02*I0[21] + U42*I1[21];
*(vp++) = U02*I0[22] + U42*I1[22]
           + (oneo2zn)*I4[15];
*(vp++) = U02*I0[23] + U42*I1[23]
           + (twoo2zn)*I4[16];
*(vp++) = U02*I0[24] + U42*I1[24]
           + (threeo2zn)*I4[17];
*(vp++) = U02*I0[25] + U42*I1[25]
           + (fouro2zn)*I4[18];
*(vp++) = U02*I0[26] + U42*I1[26]
           + (fiveo2zn)*I4[19];
*(vp++) = U02*I0[27] + U42*I1[27]
           + (sixo2zn)*I4[20];

}
/* Total number of FLOPs = 378 */
