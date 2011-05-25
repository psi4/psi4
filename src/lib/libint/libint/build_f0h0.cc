  /* These machine-generated functions compute a quartet of (fs|hs) integrals */

#include "libint.h"

void _build_f0h0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
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
  REALTYPE oneo2z;
  REALTYPE twoo2z;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  threeo2zn = 3.0*Data->oo2zn;
  fouro2zn = 4.0*Data->oo2zn;
  fiveo2zn = 5.0*Data->oo2zn;
  oneo2z = 1.0*Data->oo2z;
  twoo2z = 2.0*Data->oo2z;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (twoo2z)*(I2[0] - (lpoz)*I3[0])
           + (fiveo2zn)*I4[0];
*(vp++) = U00*I0[1] + U40*I1[1]
           + (twoo2z)*(I2[1] - (lpoz)*I3[1])
           + (fouro2zn)*I4[1];
*(vp++) = U00*I0[2] + U40*I1[2]
           + (twoo2z)*(I2[2] - (lpoz)*I3[2])
           + (fouro2zn)*I4[2];
*(vp++) = U00*I0[3] + U40*I1[3]
           + (twoo2z)*(I2[3] - (lpoz)*I3[3])
           + (threeo2zn)*I4[3];
*(vp++) = U00*I0[4] + U40*I1[4]
           + (twoo2z)*(I2[4] - (lpoz)*I3[4])
           + (threeo2zn)*I4[4];
*(vp++) = U00*I0[5] + U40*I1[5]
           + (twoo2z)*(I2[5] - (lpoz)*I3[5])
           + (threeo2zn)*I4[5];
*(vp++) = U00*I0[6] + U40*I1[6]
           + (twoo2z)*(I2[6] - (lpoz)*I3[6])
           + (twoo2zn)*I4[6];
*(vp++) = U00*I0[7] + U40*I1[7]
           + (twoo2z)*(I2[7] - (lpoz)*I3[7])
           + (twoo2zn)*I4[7];
*(vp++) = U00*I0[8] + U40*I1[8]
           + (twoo2z)*(I2[8] - (lpoz)*I3[8])
           + (twoo2zn)*I4[8];
*(vp++) = U00*I0[9] + U40*I1[9]
           + (twoo2z)*(I2[9] - (lpoz)*I3[9])
           + (twoo2zn)*I4[9];
*(vp++) = U00*I0[10] + U40*I1[10]
           + (twoo2z)*(I2[10] - (lpoz)*I3[10])
           + (oneo2zn)*I4[10];
*(vp++) = U00*I0[11] + U40*I1[11]
           + (twoo2z)*(I2[11] - (lpoz)*I3[11])
           + (oneo2zn)*I4[11];
*(vp++) = U00*I0[12] + U40*I1[12]
           + (twoo2z)*(I2[12] - (lpoz)*I3[12])
           + (oneo2zn)*I4[12];
*(vp++) = U00*I0[13] + U40*I1[13]
           + (twoo2z)*(I2[13] - (lpoz)*I3[13])
           + (oneo2zn)*I4[13];
*(vp++) = U00*I0[14] + U40*I1[14]
           + (twoo2z)*(I2[14] - (lpoz)*I3[14])
           + (oneo2zn)*I4[14];
*(vp++) = U00*I0[15] + U40*I1[15]
           + (twoo2z)*(I2[15] - (lpoz)*I3[15]);
*(vp++) = U00*I0[16] + U40*I1[16]
           + (twoo2z)*(I2[16] - (lpoz)*I3[16]);
*(vp++) = U00*I0[17] + U40*I1[17]
           + (twoo2z)*(I2[17] - (lpoz)*I3[17]);
*(vp++) = U00*I0[18] + U40*I1[18]
           + (twoo2z)*(I2[18] - (lpoz)*I3[18]);
*(vp++) = U00*I0[19] + U40*I1[19]
           + (twoo2z)*(I2[19] - (lpoz)*I3[19]);
*(vp++) = U00*I0[20] + U40*I1[20]
           + (twoo2z)*(I2[20] - (lpoz)*I3[20]);
*(vp++) = U00*I0[21] + U40*I1[21]
           + (oneo2z)*(I2[21] - (lpoz)*I3[21])
           + (fiveo2zn)*I4[15];
*(vp++) = U00*I0[22] + U40*I1[22]
           + (oneo2z)*(I2[22] - (lpoz)*I3[22])
           + (fouro2zn)*I4[16];
*(vp++) = U00*I0[23] + U40*I1[23]
           + (oneo2z)*(I2[23] - (lpoz)*I3[23])
           + (fouro2zn)*I4[17];
*(vp++) = U00*I0[24] + U40*I1[24]
           + (oneo2z)*(I2[24] - (lpoz)*I3[24])
           + (threeo2zn)*I4[18];
*(vp++) = U00*I0[25] + U40*I1[25]
           + (oneo2z)*(I2[25] - (lpoz)*I3[25])
           + (threeo2zn)*I4[19];
*(vp++) = U00*I0[26] + U40*I1[26]
           + (oneo2z)*(I2[26] - (lpoz)*I3[26])
           + (threeo2zn)*I4[20];
*(vp++) = U00*I0[27] + U40*I1[27]
           + (oneo2z)*(I2[27] - (lpoz)*I3[27])
           + (twoo2zn)*I4[21];
*(vp++) = U00*I0[28] + U40*I1[28]
           + (oneo2z)*(I2[28] - (lpoz)*I3[28])
           + (twoo2zn)*I4[22];
*(vp++) = U00*I0[29] + U40*I1[29]
           + (oneo2z)*(I2[29] - (lpoz)*I3[29])
           + (twoo2zn)*I4[23];
*(vp++) = U00*I0[30] + U40*I1[30]
           + (oneo2z)*(I2[30] - (lpoz)*I3[30])
           + (twoo2zn)*I4[24];
*(vp++) = U00*I0[31] + U40*I1[31]
           + (oneo2z)*(I2[31] - (lpoz)*I3[31])
           + (oneo2zn)*I4[25];
*(vp++) = U00*I0[32] + U40*I1[32]
           + (oneo2z)*(I2[32] - (lpoz)*I3[32])
           + (oneo2zn)*I4[26];
*(vp++) = U00*I0[33] + U40*I1[33]
           + (oneo2z)*(I2[33] - (lpoz)*I3[33])
           + (oneo2zn)*I4[27];
*(vp++) = U00*I0[34] + U40*I1[34]
           + (oneo2z)*(I2[34] - (lpoz)*I3[34])
           + (oneo2zn)*I4[28];
*(vp++) = U00*I0[35] + U40*I1[35]
           + (oneo2z)*(I2[35] - (lpoz)*I3[35])
           + (oneo2zn)*I4[29];
*(vp++) = U00*I0[36] + U40*I1[36]
           + (oneo2z)*(I2[36] - (lpoz)*I3[36]);
*(vp++) = U00*I0[37] + U40*I1[37]
           + (oneo2z)*(I2[37] - (lpoz)*I3[37]);
*(vp++) = U00*I0[38] + U40*I1[38]
           + (oneo2z)*(I2[38] - (lpoz)*I3[38]);
*(vp++) = U00*I0[39] + U40*I1[39]
           + (oneo2z)*(I2[39] - (lpoz)*I3[39]);
*(vp++) = U00*I0[40] + U40*I1[40]
           + (oneo2z)*(I2[40] - (lpoz)*I3[40]);
*(vp++) = U00*I0[41] + U40*I1[41]
           + (oneo2z)*(I2[41] - (lpoz)*I3[41]);
*(vp++) = U00*I0[42] + U40*I1[42]
           + (oneo2z)*(I2[42] - (lpoz)*I3[42])
           + (fiveo2zn)*I4[30];
*(vp++) = U00*I0[43] + U40*I1[43]
           + (oneo2z)*(I2[43] - (lpoz)*I3[43])
           + (fouro2zn)*I4[31];
*(vp++) = U00*I0[44] + U40*I1[44]
           + (oneo2z)*(I2[44] - (lpoz)*I3[44])
           + (fouro2zn)*I4[32];
*(vp++) = U00*I0[45] + U40*I1[45]
           + (oneo2z)*(I2[45] - (lpoz)*I3[45])
           + (threeo2zn)*I4[33];
*(vp++) = U00*I0[46] + U40*I1[46]
           + (oneo2z)*(I2[46] - (lpoz)*I3[46])
           + (threeo2zn)*I4[34];
*(vp++) = U00*I0[47] + U40*I1[47]
           + (oneo2z)*(I2[47] - (lpoz)*I3[47])
           + (threeo2zn)*I4[35];
*(vp++) = U00*I0[48] + U40*I1[48]
           + (oneo2z)*(I2[48] - (lpoz)*I3[48])
           + (twoo2zn)*I4[36];
*(vp++) = U00*I0[49] + U40*I1[49]
           + (oneo2z)*(I2[49] - (lpoz)*I3[49])
           + (twoo2zn)*I4[37];
*(vp++) = U00*I0[50] + U40*I1[50]
           + (oneo2z)*(I2[50] - (lpoz)*I3[50])
           + (twoo2zn)*I4[38];
*(vp++) = U00*I0[51] + U40*I1[51]
           + (oneo2z)*(I2[51] - (lpoz)*I3[51])
           + (twoo2zn)*I4[39];
*(vp++) = U00*I0[52] + U40*I1[52]
           + (oneo2z)*(I2[52] - (lpoz)*I3[52])
           + (oneo2zn)*I4[40];
*(vp++) = U00*I0[53] + U40*I1[53]
           + (oneo2z)*(I2[53] - (lpoz)*I3[53])
           + (oneo2zn)*I4[41];
*(vp++) = U00*I0[54] + U40*I1[54]
           + (oneo2z)*(I2[54] - (lpoz)*I3[54])
           + (oneo2zn)*I4[42];
*(vp++) = U00*I0[55] + U40*I1[55]
           + (oneo2z)*(I2[55] - (lpoz)*I3[55])
           + (oneo2zn)*I4[43];
*(vp++) = U00*I0[56] + U40*I1[56]
           + (oneo2z)*(I2[56] - (lpoz)*I3[56])
           + (oneo2zn)*I4[44];
*(vp++) = U00*I0[57] + U40*I1[57]
           + (oneo2z)*(I2[57] - (lpoz)*I3[57]);
*(vp++) = U00*I0[58] + U40*I1[58]
           + (oneo2z)*(I2[58] - (lpoz)*I3[58]);
*(vp++) = U00*I0[59] + U40*I1[59]
           + (oneo2z)*(I2[59] - (lpoz)*I3[59]);
*(vp++) = U00*I0[60] + U40*I1[60]
           + (oneo2z)*(I2[60] - (lpoz)*I3[60]);
*(vp++) = U00*I0[61] + U40*I1[61]
           + (oneo2z)*(I2[61] - (lpoz)*I3[61]);
*(vp++) = U00*I0[62] + U40*I1[62]
           + (oneo2z)*(I2[62] - (lpoz)*I3[62]);
*(vp++) = U00*I0[63] + U40*I1[63]
           + (fiveo2zn)*I4[45];
*(vp++) = U00*I0[64] + U40*I1[64]
           + (fouro2zn)*I4[46];
*(vp++) = U00*I0[65] + U40*I1[65]
           + (fouro2zn)*I4[47];
*(vp++) = U00*I0[66] + U40*I1[66]
           + (threeo2zn)*I4[48];
*(vp++) = U00*I0[67] + U40*I1[67]
           + (threeo2zn)*I4[49];
*(vp++) = U00*I0[68] + U40*I1[68]
           + (threeo2zn)*I4[50];
*(vp++) = U00*I0[69] + U40*I1[69]
           + (twoo2zn)*I4[51];
*(vp++) = U00*I0[70] + U40*I1[70]
           + (twoo2zn)*I4[52];
*(vp++) = U00*I0[71] + U40*I1[71]
           + (twoo2zn)*I4[53];
*(vp++) = U00*I0[72] + U40*I1[72]
           + (twoo2zn)*I4[54];
*(vp++) = U00*I0[73] + U40*I1[73]
           + (oneo2zn)*I4[55];
*(vp++) = U00*I0[74] + U40*I1[74]
           + (oneo2zn)*I4[56];
*(vp++) = U00*I0[75] + U40*I1[75]
           + (oneo2zn)*I4[57];
*(vp++) = U00*I0[76] + U40*I1[76]
           + (oneo2zn)*I4[58];
*(vp++) = U00*I0[77] + U40*I1[77]
           + (oneo2zn)*I4[59];
*(vp++) = U00*I0[78] + U40*I1[78];
*(vp++) = U00*I0[79] + U40*I1[79];
*(vp++) = U00*I0[80] + U40*I1[80];
*(vp++) = U00*I0[81] + U40*I1[81];
*(vp++) = U00*I0[82] + U40*I1[82];
*(vp++) = U00*I0[83] + U40*I1[83];
*(vp++) = U00*I0[84] + U40*I1[84]
           + (fiveo2zn)*I4[60];
*(vp++) = U00*I0[85] + U40*I1[85]
           + (fouro2zn)*I4[61];
*(vp++) = U00*I0[86] + U40*I1[86]
           + (fouro2zn)*I4[62];
*(vp++) = U00*I0[87] + U40*I1[87]
           + (threeo2zn)*I4[63];
*(vp++) = U00*I0[88] + U40*I1[88]
           + (threeo2zn)*I4[64];
*(vp++) = U00*I0[89] + U40*I1[89]
           + (threeo2zn)*I4[65];
*(vp++) = U00*I0[90] + U40*I1[90]
           + (twoo2zn)*I4[66];
*(vp++) = U00*I0[91] + U40*I1[91]
           + (twoo2zn)*I4[67];
*(vp++) = U00*I0[92] + U40*I1[92]
           + (twoo2zn)*I4[68];
*(vp++) = U00*I0[93] + U40*I1[93]
           + (twoo2zn)*I4[69];
*(vp++) = U00*I0[94] + U40*I1[94]
           + (oneo2zn)*I4[70];
*(vp++) = U00*I0[95] + U40*I1[95]
           + (oneo2zn)*I4[71];
*(vp++) = U00*I0[96] + U40*I1[96]
           + (oneo2zn)*I4[72];
*(vp++) = U00*I0[97] + U40*I1[97]
           + (oneo2zn)*I4[73];
*(vp++) = U00*I0[98] + U40*I1[98]
           + (oneo2zn)*I4[74];
*(vp++) = U00*I0[99] + U40*I1[99];
*(vp++) = U00*I0[100] + U40*I1[100];
*(vp++) = U00*I0[101] + U40*I1[101];
*(vp++) = U00*I0[102] + U40*I1[102];
*(vp++) = U00*I0[103] + U40*I1[103];
*(vp++) = U00*I0[104] + U40*I1[104];
*(vp++) = U00*I0[105] + U40*I1[105]
           + (fiveo2zn)*I4[75];
*(vp++) = U00*I0[106] + U40*I1[106]
           + (fouro2zn)*I4[76];
*(vp++) = U00*I0[107] + U40*I1[107]
           + (fouro2zn)*I4[77];
*(vp++) = U00*I0[108] + U40*I1[108]
           + (threeo2zn)*I4[78];
*(vp++) = U00*I0[109] + U40*I1[109]
           + (threeo2zn)*I4[79];
*(vp++) = U00*I0[110] + U40*I1[110]
           + (threeo2zn)*I4[80];
*(vp++) = U00*I0[111] + U40*I1[111]
           + (twoo2zn)*I4[81];
*(vp++) = U00*I0[112] + U40*I1[112]
           + (twoo2zn)*I4[82];
*(vp++) = U00*I0[113] + U40*I1[113]
           + (twoo2zn)*I4[83];
*(vp++) = U00*I0[114] + U40*I1[114]
           + (twoo2zn)*I4[84];
*(vp++) = U00*I0[115] + U40*I1[115]
           + (oneo2zn)*I4[85];
*(vp++) = U00*I0[116] + U40*I1[116]
           + (oneo2zn)*I4[86];
*(vp++) = U00*I0[117] + U40*I1[117]
           + (oneo2zn)*I4[87];
*(vp++) = U00*I0[118] + U40*I1[118]
           + (oneo2zn)*I4[88];
*(vp++) = U00*I0[119] + U40*I1[119]
           + (oneo2zn)*I4[89];
*(vp++) = U00*I0[120] + U40*I1[120];
*(vp++) = U00*I0[121] + U40*I1[121];
*(vp++) = U00*I0[122] + U40*I1[122];
*(vp++) = U00*I0[123] + U40*I1[123];
*(vp++) = U00*I0[124] + U40*I1[124];
*(vp++) = U00*I0[125] + U40*I1[125];
*(vp++) = U01*I0[63] + U41*I1[63]
           + (twoo2z)*(I2[21] - (lpoz)*I3[21]);
*(vp++) = U01*I0[64] + U41*I1[64]
           + (twoo2z)*(I2[22] - (lpoz)*I3[22])
           + (oneo2zn)*I4[45];
*(vp++) = U01*I0[65] + U41*I1[65]
           + (twoo2z)*(I2[23] - (lpoz)*I3[23]);
*(vp++) = U01*I0[66] + U41*I1[66]
           + (twoo2z)*(I2[24] - (lpoz)*I3[24])
           + (twoo2zn)*I4[46];
*(vp++) = U01*I0[67] + U41*I1[67]
           + (twoo2z)*(I2[25] - (lpoz)*I3[25])
           + (oneo2zn)*I4[47];
*(vp++) = U01*I0[68] + U41*I1[68]
           + (twoo2z)*(I2[26] - (lpoz)*I3[26]);
*(vp++) = U01*I0[69] + U41*I1[69]
           + (twoo2z)*(I2[27] - (lpoz)*I3[27])
           + (threeo2zn)*I4[48];
*(vp++) = U01*I0[70] + U41*I1[70]
           + (twoo2z)*(I2[28] - (lpoz)*I3[28])
           + (twoo2zn)*I4[49];
*(vp++) = U01*I0[71] + U41*I1[71]
           + (twoo2z)*(I2[29] - (lpoz)*I3[29])
           + (oneo2zn)*I4[50];
*(vp++) = U01*I0[72] + U41*I1[72]
           + (twoo2z)*(I2[30] - (lpoz)*I3[30]);
*(vp++) = U01*I0[73] + U41*I1[73]
           + (twoo2z)*(I2[31] - (lpoz)*I3[31])
           + (fouro2zn)*I4[51];
*(vp++) = U01*I0[74] + U41*I1[74]
           + (twoo2z)*(I2[32] - (lpoz)*I3[32])
           + (threeo2zn)*I4[52];
*(vp++) = U01*I0[75] + U41*I1[75]
           + (twoo2z)*(I2[33] - (lpoz)*I3[33])
           + (twoo2zn)*I4[53];
*(vp++) = U01*I0[76] + U41*I1[76]
           + (twoo2z)*(I2[34] - (lpoz)*I3[34])
           + (oneo2zn)*I4[54];
*(vp++) = U01*I0[77] + U41*I1[77]
           + (twoo2z)*(I2[35] - (lpoz)*I3[35]);
*(vp++) = U01*I0[78] + U41*I1[78]
           + (twoo2z)*(I2[36] - (lpoz)*I3[36])
           + (fiveo2zn)*I4[55];
*(vp++) = U01*I0[79] + U41*I1[79]
           + (twoo2z)*(I2[37] - (lpoz)*I3[37])
           + (fouro2zn)*I4[56];
*(vp++) = U01*I0[80] + U41*I1[80]
           + (twoo2z)*(I2[38] - (lpoz)*I3[38])
           + (threeo2zn)*I4[57];
*(vp++) = U01*I0[81] + U41*I1[81]
           + (twoo2z)*(I2[39] - (lpoz)*I3[39])
           + (twoo2zn)*I4[58];
*(vp++) = U01*I0[82] + U41*I1[82]
           + (twoo2z)*(I2[40] - (lpoz)*I3[40])
           + (oneo2zn)*I4[59];
*(vp++) = U01*I0[83] + U41*I1[83]
           + (twoo2z)*(I2[41] - (lpoz)*I3[41]);
*(vp++) = U01*I0[84] + U41*I1[84]
           + (oneo2z)*(I2[42] - (lpoz)*I3[42]);
*(vp++) = U01*I0[85] + U41*I1[85]
           + (oneo2z)*(I2[43] - (lpoz)*I3[43])
           + (oneo2zn)*I4[60];
*(vp++) = U01*I0[86] + U41*I1[86]
           + (oneo2z)*(I2[44] - (lpoz)*I3[44]);
*(vp++) = U01*I0[87] + U41*I1[87]
           + (oneo2z)*(I2[45] - (lpoz)*I3[45])
           + (twoo2zn)*I4[61];
*(vp++) = U01*I0[88] + U41*I1[88]
           + (oneo2z)*(I2[46] - (lpoz)*I3[46])
           + (oneo2zn)*I4[62];
*(vp++) = U01*I0[89] + U41*I1[89]
           + (oneo2z)*(I2[47] - (lpoz)*I3[47]);
*(vp++) = U01*I0[90] + U41*I1[90]
           + (oneo2z)*(I2[48] - (lpoz)*I3[48])
           + (threeo2zn)*I4[63];
*(vp++) = U01*I0[91] + U41*I1[91]
           + (oneo2z)*(I2[49] - (lpoz)*I3[49])
           + (twoo2zn)*I4[64];
*(vp++) = U01*I0[92] + U41*I1[92]
           + (oneo2z)*(I2[50] - (lpoz)*I3[50])
           + (oneo2zn)*I4[65];
*(vp++) = U01*I0[93] + U41*I1[93]
           + (oneo2z)*(I2[51] - (lpoz)*I3[51]);
*(vp++) = U01*I0[94] + U41*I1[94]
           + (oneo2z)*(I2[52] - (lpoz)*I3[52])
           + (fouro2zn)*I4[66];
*(vp++) = U01*I0[95] + U41*I1[95]
           + (oneo2z)*(I2[53] - (lpoz)*I3[53])
           + (threeo2zn)*I4[67];
*(vp++) = U01*I0[96] + U41*I1[96]
           + (oneo2z)*(I2[54] - (lpoz)*I3[54])
           + (twoo2zn)*I4[68];
*(vp++) = U01*I0[97] + U41*I1[97]
           + (oneo2z)*(I2[55] - (lpoz)*I3[55])
           + (oneo2zn)*I4[69];
*(vp++) = U01*I0[98] + U41*I1[98]
           + (oneo2z)*(I2[56] - (lpoz)*I3[56]);
*(vp++) = U01*I0[99] + U41*I1[99]
           + (oneo2z)*(I2[57] - (lpoz)*I3[57])
           + (fiveo2zn)*I4[70];
*(vp++) = U01*I0[100] + U41*I1[100]
           + (oneo2z)*(I2[58] - (lpoz)*I3[58])
           + (fouro2zn)*I4[71];
*(vp++) = U01*I0[101] + U41*I1[101]
           + (oneo2z)*(I2[59] - (lpoz)*I3[59])
           + (threeo2zn)*I4[72];
*(vp++) = U01*I0[102] + U41*I1[102]
           + (oneo2z)*(I2[60] - (lpoz)*I3[60])
           + (twoo2zn)*I4[73];
*(vp++) = U01*I0[103] + U41*I1[103]
           + (oneo2z)*(I2[61] - (lpoz)*I3[61])
           + (oneo2zn)*I4[74];
*(vp++) = U01*I0[104] + U41*I1[104]
           + (oneo2z)*(I2[62] - (lpoz)*I3[62]);
*(vp++) = U01*I0[105] + U41*I1[105];
*(vp++) = U01*I0[106] + U41*I1[106]
           + (oneo2zn)*I4[75];
*(vp++) = U01*I0[107] + U41*I1[107];
*(vp++) = U01*I0[108] + U41*I1[108]
           + (twoo2zn)*I4[76];
*(vp++) = U01*I0[109] + U41*I1[109]
           + (oneo2zn)*I4[77];
*(vp++) = U01*I0[110] + U41*I1[110];
*(vp++) = U01*I0[111] + U41*I1[111]
           + (threeo2zn)*I4[78];
*(vp++) = U01*I0[112] + U41*I1[112]
           + (twoo2zn)*I4[79];
*(vp++) = U01*I0[113] + U41*I1[113]
           + (oneo2zn)*I4[80];
*(vp++) = U01*I0[114] + U41*I1[114];
*(vp++) = U01*I0[115] + U41*I1[115]
           + (fouro2zn)*I4[81];
*(vp++) = U01*I0[116] + U41*I1[116]
           + (threeo2zn)*I4[82];
*(vp++) = U01*I0[117] + U41*I1[117]
           + (twoo2zn)*I4[83];
*(vp++) = U01*I0[118] + U41*I1[118]
           + (oneo2zn)*I4[84];
*(vp++) = U01*I0[119] + U41*I1[119];
*(vp++) = U01*I0[120] + U41*I1[120]
           + (fiveo2zn)*I4[85];
*(vp++) = U01*I0[121] + U41*I1[121]
           + (fouro2zn)*I4[86];
*(vp++) = U01*I0[122] + U41*I1[122]
           + (threeo2zn)*I4[87];
*(vp++) = U01*I0[123] + U41*I1[123]
           + (twoo2zn)*I4[88];
*(vp++) = U01*I0[124] + U41*I1[124]
           + (oneo2zn)*I4[89];
*(vp++) = U01*I0[125] + U41*I1[125];
*(vp++) = U02*I0[105] + U42*I1[105]
           + (twoo2z)*(I2[42] - (lpoz)*I3[42]);
*(vp++) = U02*I0[106] + U42*I1[106]
           + (twoo2z)*(I2[43] - (lpoz)*I3[43]);
*(vp++) = U02*I0[107] + U42*I1[107]
           + (twoo2z)*(I2[44] - (lpoz)*I3[44])
           + (oneo2zn)*I4[75];
*(vp++) = U02*I0[108] + U42*I1[108]
           + (twoo2z)*(I2[45] - (lpoz)*I3[45]);
*(vp++) = U02*I0[109] + U42*I1[109]
           + (twoo2z)*(I2[46] - (lpoz)*I3[46])
           + (oneo2zn)*I4[76];
*(vp++) = U02*I0[110] + U42*I1[110]
           + (twoo2z)*(I2[47] - (lpoz)*I3[47])
           + (twoo2zn)*I4[77];
*(vp++) = U02*I0[111] + U42*I1[111]
           + (twoo2z)*(I2[48] - (lpoz)*I3[48]);
*(vp++) = U02*I0[112] + U42*I1[112]
           + (twoo2z)*(I2[49] - (lpoz)*I3[49])
           + (oneo2zn)*I4[78];
*(vp++) = U02*I0[113] + U42*I1[113]
           + (twoo2z)*(I2[50] - (lpoz)*I3[50])
           + (twoo2zn)*I4[79];
*(vp++) = U02*I0[114] + U42*I1[114]
           + (twoo2z)*(I2[51] - (lpoz)*I3[51])
           + (threeo2zn)*I4[80];
*(vp++) = U02*I0[115] + U42*I1[115]
           + (twoo2z)*(I2[52] - (lpoz)*I3[52]);
*(vp++) = U02*I0[116] + U42*I1[116]
           + (twoo2z)*(I2[53] - (lpoz)*I3[53])
           + (oneo2zn)*I4[81];
*(vp++) = U02*I0[117] + U42*I1[117]
           + (twoo2z)*(I2[54] - (lpoz)*I3[54])
           + (twoo2zn)*I4[82];
*(vp++) = U02*I0[118] + U42*I1[118]
           + (twoo2z)*(I2[55] - (lpoz)*I3[55])
           + (threeo2zn)*I4[83];
*(vp++) = U02*I0[119] + U42*I1[119]
           + (twoo2z)*(I2[56] - (lpoz)*I3[56])
           + (fouro2zn)*I4[84];
*(vp++) = U02*I0[120] + U42*I1[120]
           + (twoo2z)*(I2[57] - (lpoz)*I3[57]);
*(vp++) = U02*I0[121] + U42*I1[121]
           + (twoo2z)*(I2[58] - (lpoz)*I3[58])
           + (oneo2zn)*I4[85];
*(vp++) = U02*I0[122] + U42*I1[122]
           + (twoo2z)*(I2[59] - (lpoz)*I3[59])
           + (twoo2zn)*I4[86];
*(vp++) = U02*I0[123] + U42*I1[123]
           + (twoo2z)*(I2[60] - (lpoz)*I3[60])
           + (threeo2zn)*I4[87];
*(vp++) = U02*I0[124] + U42*I1[124]
           + (twoo2z)*(I2[61] - (lpoz)*I3[61])
           + (fouro2zn)*I4[88];
*(vp++) = U02*I0[125] + U42*I1[125]
           + (twoo2z)*(I2[62] - (lpoz)*I3[62])
           + (fiveo2zn)*I4[89];

}
/* Total number of FLOPs = 1434 */
