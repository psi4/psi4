  /* These machine-generated functions compute a quartet of (ds|is) integrals */

#include "libint.h"

void _build_d0i0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
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
  REALTYPE oneo2z;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  threeo2zn = 3.0*Data->oo2zn;
  fouro2zn = 4.0*Data->oo2zn;
  fiveo2zn = 5.0*Data->oo2zn;
  sixo2zn = 6.0*Data->oo2zn;
  oneo2z = 1.0*Data->oo2z;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (oneo2z)*(I2[0] - (lpoz)*I3[0])
           + (sixo2zn)*I4[0];
*(vp++) = U00*I0[1] + U40*I1[1]
           + (oneo2z)*(I2[1] - (lpoz)*I3[1])
           + (fiveo2zn)*I4[1];
*(vp++) = U00*I0[2] + U40*I1[2]
           + (oneo2z)*(I2[2] - (lpoz)*I3[2])
           + (fiveo2zn)*I4[2];
*(vp++) = U00*I0[3] + U40*I1[3]
           + (oneo2z)*(I2[3] - (lpoz)*I3[3])
           + (fouro2zn)*I4[3];
*(vp++) = U00*I0[4] + U40*I1[4]
           + (oneo2z)*(I2[4] - (lpoz)*I3[4])
           + (fouro2zn)*I4[4];
*(vp++) = U00*I0[5] + U40*I1[5]
           + (oneo2z)*(I2[5] - (lpoz)*I3[5])
           + (fouro2zn)*I4[5];
*(vp++) = U00*I0[6] + U40*I1[6]
           + (oneo2z)*(I2[6] - (lpoz)*I3[6])
           + (threeo2zn)*I4[6];
*(vp++) = U00*I0[7] + U40*I1[7]
           + (oneo2z)*(I2[7] - (lpoz)*I3[7])
           + (threeo2zn)*I4[7];
*(vp++) = U00*I0[8] + U40*I1[8]
           + (oneo2z)*(I2[8] - (lpoz)*I3[8])
           + (threeo2zn)*I4[8];
*(vp++) = U00*I0[9] + U40*I1[9]
           + (oneo2z)*(I2[9] - (lpoz)*I3[9])
           + (threeo2zn)*I4[9];
*(vp++) = U00*I0[10] + U40*I1[10]
           + (oneo2z)*(I2[10] - (lpoz)*I3[10])
           + (twoo2zn)*I4[10];
*(vp++) = U00*I0[11] + U40*I1[11]
           + (oneo2z)*(I2[11] - (lpoz)*I3[11])
           + (twoo2zn)*I4[11];
*(vp++) = U00*I0[12] + U40*I1[12]
           + (oneo2z)*(I2[12] - (lpoz)*I3[12])
           + (twoo2zn)*I4[12];
*(vp++) = U00*I0[13] + U40*I1[13]
           + (oneo2z)*(I2[13] - (lpoz)*I3[13])
           + (twoo2zn)*I4[13];
*(vp++) = U00*I0[14] + U40*I1[14]
           + (oneo2z)*(I2[14] - (lpoz)*I3[14])
           + (twoo2zn)*I4[14];
*(vp++) = U00*I0[15] + U40*I1[15]
           + (oneo2z)*(I2[15] - (lpoz)*I3[15])
           + (oneo2zn)*I4[15];
*(vp++) = U00*I0[16] + U40*I1[16]
           + (oneo2z)*(I2[16] - (lpoz)*I3[16])
           + (oneo2zn)*I4[16];
*(vp++) = U00*I0[17] + U40*I1[17]
           + (oneo2z)*(I2[17] - (lpoz)*I3[17])
           + (oneo2zn)*I4[17];
*(vp++) = U00*I0[18] + U40*I1[18]
           + (oneo2z)*(I2[18] - (lpoz)*I3[18])
           + (oneo2zn)*I4[18];
*(vp++) = U00*I0[19] + U40*I1[19]
           + (oneo2z)*(I2[19] - (lpoz)*I3[19])
           + (oneo2zn)*I4[19];
*(vp++) = U00*I0[20] + U40*I1[20]
           + (oneo2z)*(I2[20] - (lpoz)*I3[20])
           + (oneo2zn)*I4[20];
*(vp++) = U00*I0[21] + U40*I1[21]
           + (oneo2z)*(I2[21] - (lpoz)*I3[21]);
*(vp++) = U00*I0[22] + U40*I1[22]
           + (oneo2z)*(I2[22] - (lpoz)*I3[22]);
*(vp++) = U00*I0[23] + U40*I1[23]
           + (oneo2z)*(I2[23] - (lpoz)*I3[23]);
*(vp++) = U00*I0[24] + U40*I1[24]
           + (oneo2z)*(I2[24] - (lpoz)*I3[24]);
*(vp++) = U00*I0[25] + U40*I1[25]
           + (oneo2z)*(I2[25] - (lpoz)*I3[25]);
*(vp++) = U00*I0[26] + U40*I1[26]
           + (oneo2z)*(I2[26] - (lpoz)*I3[26]);
*(vp++) = U00*I0[27] + U40*I1[27]
           + (oneo2z)*(I2[27] - (lpoz)*I3[27]);
*(vp++) = U00*I0[28] + U40*I1[28]
           + (sixo2zn)*I4[21];
*(vp++) = U00*I0[29] + U40*I1[29]
           + (fiveo2zn)*I4[22];
*(vp++) = U00*I0[30] + U40*I1[30]
           + (fiveo2zn)*I4[23];
*(vp++) = U00*I0[31] + U40*I1[31]
           + (fouro2zn)*I4[24];
*(vp++) = U00*I0[32] + U40*I1[32]
           + (fouro2zn)*I4[25];
*(vp++) = U00*I0[33] + U40*I1[33]
           + (fouro2zn)*I4[26];
*(vp++) = U00*I0[34] + U40*I1[34]
           + (threeo2zn)*I4[27];
*(vp++) = U00*I0[35] + U40*I1[35]
           + (threeo2zn)*I4[28];
*(vp++) = U00*I0[36] + U40*I1[36]
           + (threeo2zn)*I4[29];
*(vp++) = U00*I0[37] + U40*I1[37]
           + (threeo2zn)*I4[30];
*(vp++) = U00*I0[38] + U40*I1[38]
           + (twoo2zn)*I4[31];
*(vp++) = U00*I0[39] + U40*I1[39]
           + (twoo2zn)*I4[32];
*(vp++) = U00*I0[40] + U40*I1[40]
           + (twoo2zn)*I4[33];
*(vp++) = U00*I0[41] + U40*I1[41]
           + (twoo2zn)*I4[34];
*(vp++) = U00*I0[42] + U40*I1[42]
           + (twoo2zn)*I4[35];
*(vp++) = U00*I0[43] + U40*I1[43]
           + (oneo2zn)*I4[36];
*(vp++) = U00*I0[44] + U40*I1[44]
           + (oneo2zn)*I4[37];
*(vp++) = U00*I0[45] + U40*I1[45]
           + (oneo2zn)*I4[38];
*(vp++) = U00*I0[46] + U40*I1[46]
           + (oneo2zn)*I4[39];
*(vp++) = U00*I0[47] + U40*I1[47]
           + (oneo2zn)*I4[40];
*(vp++) = U00*I0[48] + U40*I1[48]
           + (oneo2zn)*I4[41];
*(vp++) = U00*I0[49] + U40*I1[49];
*(vp++) = U00*I0[50] + U40*I1[50];
*(vp++) = U00*I0[51] + U40*I1[51];
*(vp++) = U00*I0[52] + U40*I1[52];
*(vp++) = U00*I0[53] + U40*I1[53];
*(vp++) = U00*I0[54] + U40*I1[54];
*(vp++) = U00*I0[55] + U40*I1[55];
*(vp++) = U00*I0[56] + U40*I1[56]
           + (sixo2zn)*I4[42];
*(vp++) = U00*I0[57] + U40*I1[57]
           + (fiveo2zn)*I4[43];
*(vp++) = U00*I0[58] + U40*I1[58]
           + (fiveo2zn)*I4[44];
*(vp++) = U00*I0[59] + U40*I1[59]
           + (fouro2zn)*I4[45];
*(vp++) = U00*I0[60] + U40*I1[60]
           + (fouro2zn)*I4[46];
*(vp++) = U00*I0[61] + U40*I1[61]
           + (fouro2zn)*I4[47];
*(vp++) = U00*I0[62] + U40*I1[62]
           + (threeo2zn)*I4[48];
*(vp++) = U00*I0[63] + U40*I1[63]
           + (threeo2zn)*I4[49];
*(vp++) = U00*I0[64] + U40*I1[64]
           + (threeo2zn)*I4[50];
*(vp++) = U00*I0[65] + U40*I1[65]
           + (threeo2zn)*I4[51];
*(vp++) = U00*I0[66] + U40*I1[66]
           + (twoo2zn)*I4[52];
*(vp++) = U00*I0[67] + U40*I1[67]
           + (twoo2zn)*I4[53];
*(vp++) = U00*I0[68] + U40*I1[68]
           + (twoo2zn)*I4[54];
*(vp++) = U00*I0[69] + U40*I1[69]
           + (twoo2zn)*I4[55];
*(vp++) = U00*I0[70] + U40*I1[70]
           + (twoo2zn)*I4[56];
*(vp++) = U00*I0[71] + U40*I1[71]
           + (oneo2zn)*I4[57];
*(vp++) = U00*I0[72] + U40*I1[72]
           + (oneo2zn)*I4[58];
*(vp++) = U00*I0[73] + U40*I1[73]
           + (oneo2zn)*I4[59];
*(vp++) = U00*I0[74] + U40*I1[74]
           + (oneo2zn)*I4[60];
*(vp++) = U00*I0[75] + U40*I1[75]
           + (oneo2zn)*I4[61];
*(vp++) = U00*I0[76] + U40*I1[76]
           + (oneo2zn)*I4[62];
*(vp++) = U00*I0[77] + U40*I1[77];
*(vp++) = U00*I0[78] + U40*I1[78];
*(vp++) = U00*I0[79] + U40*I1[79];
*(vp++) = U00*I0[80] + U40*I1[80];
*(vp++) = U00*I0[81] + U40*I1[81];
*(vp++) = U00*I0[82] + U40*I1[82];
*(vp++) = U00*I0[83] + U40*I1[83];
*(vp++) = U01*I0[28] + U41*I1[28]
           + (oneo2z)*(I2[0] - (lpoz)*I3[0]);
*(vp++) = U01*I0[29] + U41*I1[29]
           + (oneo2z)*(I2[1] - (lpoz)*I3[1])
           + (oneo2zn)*I4[21];
*(vp++) = U01*I0[30] + U41*I1[30]
           + (oneo2z)*(I2[2] - (lpoz)*I3[2]);
*(vp++) = U01*I0[31] + U41*I1[31]
           + (oneo2z)*(I2[3] - (lpoz)*I3[3])
           + (twoo2zn)*I4[22];
*(vp++) = U01*I0[32] + U41*I1[32]
           + (oneo2z)*(I2[4] - (lpoz)*I3[4])
           + (oneo2zn)*I4[23];
*(vp++) = U01*I0[33] + U41*I1[33]
           + (oneo2z)*(I2[5] - (lpoz)*I3[5]);
*(vp++) = U01*I0[34] + U41*I1[34]
           + (oneo2z)*(I2[6] - (lpoz)*I3[6])
           + (threeo2zn)*I4[24];
*(vp++) = U01*I0[35] + U41*I1[35]
           + (oneo2z)*(I2[7] - (lpoz)*I3[7])
           + (twoo2zn)*I4[25];
*(vp++) = U01*I0[36] + U41*I1[36]
           + (oneo2z)*(I2[8] - (lpoz)*I3[8])
           + (oneo2zn)*I4[26];
*(vp++) = U01*I0[37] + U41*I1[37]
           + (oneo2z)*(I2[9] - (lpoz)*I3[9]);
*(vp++) = U01*I0[38] + U41*I1[38]
           + (oneo2z)*(I2[10] - (lpoz)*I3[10])
           + (fouro2zn)*I4[27];
*(vp++) = U01*I0[39] + U41*I1[39]
           + (oneo2z)*(I2[11] - (lpoz)*I3[11])
           + (threeo2zn)*I4[28];
*(vp++) = U01*I0[40] + U41*I1[40]
           + (oneo2z)*(I2[12] - (lpoz)*I3[12])
           + (twoo2zn)*I4[29];
*(vp++) = U01*I0[41] + U41*I1[41]
           + (oneo2z)*(I2[13] - (lpoz)*I3[13])
           + (oneo2zn)*I4[30];
*(vp++) = U01*I0[42] + U41*I1[42]
           + (oneo2z)*(I2[14] - (lpoz)*I3[14]);
*(vp++) = U01*I0[43] + U41*I1[43]
           + (oneo2z)*(I2[15] - (lpoz)*I3[15])
           + (fiveo2zn)*I4[31];
*(vp++) = U01*I0[44] + U41*I1[44]
           + (oneo2z)*(I2[16] - (lpoz)*I3[16])
           + (fouro2zn)*I4[32];
*(vp++) = U01*I0[45] + U41*I1[45]
           + (oneo2z)*(I2[17] - (lpoz)*I3[17])
           + (threeo2zn)*I4[33];
*(vp++) = U01*I0[46] + U41*I1[46]
           + (oneo2z)*(I2[18] - (lpoz)*I3[18])
           + (twoo2zn)*I4[34];
*(vp++) = U01*I0[47] + U41*I1[47]
           + (oneo2z)*(I2[19] - (lpoz)*I3[19])
           + (oneo2zn)*I4[35];
*(vp++) = U01*I0[48] + U41*I1[48]
           + (oneo2z)*(I2[20] - (lpoz)*I3[20]);
*(vp++) = U01*I0[49] + U41*I1[49]
           + (oneo2z)*(I2[21] - (lpoz)*I3[21])
           + (sixo2zn)*I4[36];
*(vp++) = U01*I0[50] + U41*I1[50]
           + (oneo2z)*(I2[22] - (lpoz)*I3[22])
           + (fiveo2zn)*I4[37];
*(vp++) = U01*I0[51] + U41*I1[51]
           + (oneo2z)*(I2[23] - (lpoz)*I3[23])
           + (fouro2zn)*I4[38];
*(vp++) = U01*I0[52] + U41*I1[52]
           + (oneo2z)*(I2[24] - (lpoz)*I3[24])
           + (threeo2zn)*I4[39];
*(vp++) = U01*I0[53] + U41*I1[53]
           + (oneo2z)*(I2[25] - (lpoz)*I3[25])
           + (twoo2zn)*I4[40];
*(vp++) = U01*I0[54] + U41*I1[54]
           + (oneo2z)*(I2[26] - (lpoz)*I3[26])
           + (oneo2zn)*I4[41];
*(vp++) = U01*I0[55] + U41*I1[55]
           + (oneo2z)*(I2[27] - (lpoz)*I3[27]);
*(vp++) = U01*I0[56] + U41*I1[56];
*(vp++) = U01*I0[57] + U41*I1[57]
           + (oneo2zn)*I4[42];
*(vp++) = U01*I0[58] + U41*I1[58];
*(vp++) = U01*I0[59] + U41*I1[59]
           + (twoo2zn)*I4[43];
*(vp++) = U01*I0[60] + U41*I1[60]
           + (oneo2zn)*I4[44];
*(vp++) = U01*I0[61] + U41*I1[61];
*(vp++) = U01*I0[62] + U41*I1[62]
           + (threeo2zn)*I4[45];
*(vp++) = U01*I0[63] + U41*I1[63]
           + (twoo2zn)*I4[46];
*(vp++) = U01*I0[64] + U41*I1[64]
           + (oneo2zn)*I4[47];
*(vp++) = U01*I0[65] + U41*I1[65];
*(vp++) = U01*I0[66] + U41*I1[66]
           + (fouro2zn)*I4[48];
*(vp++) = U01*I0[67] + U41*I1[67]
           + (threeo2zn)*I4[49];
*(vp++) = U01*I0[68] + U41*I1[68]
           + (twoo2zn)*I4[50];
*(vp++) = U01*I0[69] + U41*I1[69]
           + (oneo2zn)*I4[51];
*(vp++) = U01*I0[70] + U41*I1[70];
*(vp++) = U01*I0[71] + U41*I1[71]
           + (fiveo2zn)*I4[52];
*(vp++) = U01*I0[72] + U41*I1[72]
           + (fouro2zn)*I4[53];
*(vp++) = U01*I0[73] + U41*I1[73]
           + (threeo2zn)*I4[54];
*(vp++) = U01*I0[74] + U41*I1[74]
           + (twoo2zn)*I4[55];
*(vp++) = U01*I0[75] + U41*I1[75]
           + (oneo2zn)*I4[56];
*(vp++) = U01*I0[76] + U41*I1[76];
*(vp++) = U01*I0[77] + U41*I1[77]
           + (sixo2zn)*I4[57];
*(vp++) = U01*I0[78] + U41*I1[78]
           + (fiveo2zn)*I4[58];
*(vp++) = U01*I0[79] + U41*I1[79]
           + (fouro2zn)*I4[59];
*(vp++) = U01*I0[80] + U41*I1[80]
           + (threeo2zn)*I4[60];
*(vp++) = U01*I0[81] + U41*I1[81]
           + (twoo2zn)*I4[61];
*(vp++) = U01*I0[82] + U41*I1[82]
           + (oneo2zn)*I4[62];
*(vp++) = U01*I0[83] + U41*I1[83];
*(vp++) = U02*I0[56] + U42*I1[56]
           + (oneo2z)*(I2[0] - (lpoz)*I3[0]);
*(vp++) = U02*I0[57] + U42*I1[57]
           + (oneo2z)*(I2[1] - (lpoz)*I3[1]);
*(vp++) = U02*I0[58] + U42*I1[58]
           + (oneo2z)*(I2[2] - (lpoz)*I3[2])
           + (oneo2zn)*I4[42];
*(vp++) = U02*I0[59] + U42*I1[59]
           + (oneo2z)*(I2[3] - (lpoz)*I3[3]);
*(vp++) = U02*I0[60] + U42*I1[60]
           + (oneo2z)*(I2[4] - (lpoz)*I3[4])
           + (oneo2zn)*I4[43];
*(vp++) = U02*I0[61] + U42*I1[61]
           + (oneo2z)*(I2[5] - (lpoz)*I3[5])
           + (twoo2zn)*I4[44];
*(vp++) = U02*I0[62] + U42*I1[62]
           + (oneo2z)*(I2[6] - (lpoz)*I3[6]);
*(vp++) = U02*I0[63] + U42*I1[63]
           + (oneo2z)*(I2[7] - (lpoz)*I3[7])
           + (oneo2zn)*I4[45];
*(vp++) = U02*I0[64] + U42*I1[64]
           + (oneo2z)*(I2[8] - (lpoz)*I3[8])
           + (twoo2zn)*I4[46];
*(vp++) = U02*I0[65] + U42*I1[65]
           + (oneo2z)*(I2[9] - (lpoz)*I3[9])
           + (threeo2zn)*I4[47];
*(vp++) = U02*I0[66] + U42*I1[66]
           + (oneo2z)*(I2[10] - (lpoz)*I3[10]);
*(vp++) = U02*I0[67] + U42*I1[67]
           + (oneo2z)*(I2[11] - (lpoz)*I3[11])
           + (oneo2zn)*I4[48];
*(vp++) = U02*I0[68] + U42*I1[68]
           + (oneo2z)*(I2[12] - (lpoz)*I3[12])
           + (twoo2zn)*I4[49];
*(vp++) = U02*I0[69] + U42*I1[69]
           + (oneo2z)*(I2[13] - (lpoz)*I3[13])
           + (threeo2zn)*I4[50];
*(vp++) = U02*I0[70] + U42*I1[70]
           + (oneo2z)*(I2[14] - (lpoz)*I3[14])
           + (fouro2zn)*I4[51];
*(vp++) = U02*I0[71] + U42*I1[71]
           + (oneo2z)*(I2[15] - (lpoz)*I3[15]);
*(vp++) = U02*I0[72] + U42*I1[72]
           + (oneo2z)*(I2[16] - (lpoz)*I3[16])
           + (oneo2zn)*I4[52];
*(vp++) = U02*I0[73] + U42*I1[73]
           + (oneo2z)*(I2[17] - (lpoz)*I3[17])
           + (twoo2zn)*I4[53];
*(vp++) = U02*I0[74] + U42*I1[74]
           + (oneo2z)*(I2[18] - (lpoz)*I3[18])
           + (threeo2zn)*I4[54];
*(vp++) = U02*I0[75] + U42*I1[75]
           + (oneo2z)*(I2[19] - (lpoz)*I3[19])
           + (fouro2zn)*I4[55];
*(vp++) = U02*I0[76] + U42*I1[76]
           + (oneo2z)*(I2[20] - (lpoz)*I3[20])
           + (fiveo2zn)*I4[56];
*(vp++) = U02*I0[77] + U42*I1[77]
           + (oneo2z)*(I2[21] - (lpoz)*I3[21]);
*(vp++) = U02*I0[78] + U42*I1[78]
           + (oneo2z)*(I2[22] - (lpoz)*I3[22])
           + (oneo2zn)*I4[57];
*(vp++) = U02*I0[79] + U42*I1[79]
           + (oneo2z)*(I2[23] - (lpoz)*I3[23])
           + (twoo2zn)*I4[58];
*(vp++) = U02*I0[80] + U42*I1[80]
           + (oneo2z)*(I2[24] - (lpoz)*I3[24])
           + (threeo2zn)*I4[59];
*(vp++) = U02*I0[81] + U42*I1[81]
           + (oneo2z)*(I2[25] - (lpoz)*I3[25])
           + (fouro2zn)*I4[60];
*(vp++) = U02*I0[82] + U42*I1[82]
           + (oneo2z)*(I2[26] - (lpoz)*I3[26])
           + (fiveo2zn)*I4[61];
*(vp++) = U02*I0[83] + U42*I1[83]
           + (oneo2z)*(I2[27] - (lpoz)*I3[27])
           + (sixo2zn)*I4[62];

}
/* Total number of FLOPs = 1092 */
