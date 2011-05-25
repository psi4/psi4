  /* These machine-generated functions compute a quartet of (gs|is) integrals */

#include "libint.h"

REALTYPE *_build_g0i0_0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
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
  REALTYPE twoo2z;
  REALTYPE threeo2z;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  threeo2zn = 3.0*Data->oo2zn;
  fouro2zn = 4.0*Data->oo2zn;
  fiveo2zn = 5.0*Data->oo2zn;
  sixo2zn = 6.0*Data->oo2zn;
  oneo2z = 1.0*Data->oo2z;
  twoo2z = 2.0*Data->oo2z;
  threeo2z = 3.0*Data->oo2z;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (threeo2z)*(I2[0] - (lpoz)*I3[0])
           + (sixo2zn)*I4[0];
*(vp++) = U00*I0[1] + U40*I1[1]
           + (threeo2z)*(I2[1] - (lpoz)*I3[1])
           + (fiveo2zn)*I4[1];
*(vp++) = U00*I0[2] + U40*I1[2]
           + (threeo2z)*(I2[2] - (lpoz)*I3[2])
           + (fiveo2zn)*I4[2];
*(vp++) = U00*I0[3] + U40*I1[3]
           + (threeo2z)*(I2[3] - (lpoz)*I3[3])
           + (fouro2zn)*I4[3];
*(vp++) = U00*I0[4] + U40*I1[4]
           + (threeo2z)*(I2[4] - (lpoz)*I3[4])
           + (fouro2zn)*I4[4];
*(vp++) = U00*I0[5] + U40*I1[5]
           + (threeo2z)*(I2[5] - (lpoz)*I3[5])
           + (fouro2zn)*I4[5];
*(vp++) = U00*I0[6] + U40*I1[6]
           + (threeo2z)*(I2[6] - (lpoz)*I3[6])
           + (threeo2zn)*I4[6];
*(vp++) = U00*I0[7] + U40*I1[7]
           + (threeo2z)*(I2[7] - (lpoz)*I3[7])
           + (threeo2zn)*I4[7];
*(vp++) = U00*I0[8] + U40*I1[8]
           + (threeo2z)*(I2[8] - (lpoz)*I3[8])
           + (threeo2zn)*I4[8];
*(vp++) = U00*I0[9] + U40*I1[9]
           + (threeo2z)*(I2[9] - (lpoz)*I3[9])
           + (threeo2zn)*I4[9];
*(vp++) = U00*I0[10] + U40*I1[10]
           + (threeo2z)*(I2[10] - (lpoz)*I3[10])
           + (twoo2zn)*I4[10];
*(vp++) = U00*I0[11] + U40*I1[11]
           + (threeo2z)*(I2[11] - (lpoz)*I3[11])
           + (twoo2zn)*I4[11];
*(vp++) = U00*I0[12] + U40*I1[12]
           + (threeo2z)*(I2[12] - (lpoz)*I3[12])
           + (twoo2zn)*I4[12];
*(vp++) = U00*I0[13] + U40*I1[13]
           + (threeo2z)*(I2[13] - (lpoz)*I3[13])
           + (twoo2zn)*I4[13];
*(vp++) = U00*I0[14] + U40*I1[14]
           + (threeo2z)*(I2[14] - (lpoz)*I3[14])
           + (twoo2zn)*I4[14];
*(vp++) = U00*I0[15] + U40*I1[15]
           + (threeo2z)*(I2[15] - (lpoz)*I3[15])
           + (oneo2zn)*I4[15];
*(vp++) = U00*I0[16] + U40*I1[16]
           + (threeo2z)*(I2[16] - (lpoz)*I3[16])
           + (oneo2zn)*I4[16];
*(vp++) = U00*I0[17] + U40*I1[17]
           + (threeo2z)*(I2[17] - (lpoz)*I3[17])
           + (oneo2zn)*I4[17];
*(vp++) = U00*I0[18] + U40*I1[18]
           + (threeo2z)*(I2[18] - (lpoz)*I3[18])
           + (oneo2zn)*I4[18];
*(vp++) = U00*I0[19] + U40*I1[19]
           + (threeo2z)*(I2[19] - (lpoz)*I3[19])
           + (oneo2zn)*I4[19];
*(vp++) = U00*I0[20] + U40*I1[20]
           + (threeo2z)*(I2[20] - (lpoz)*I3[20])
           + (oneo2zn)*I4[20];
*(vp++) = U00*I0[21] + U40*I1[21]
           + (threeo2z)*(I2[21] - (lpoz)*I3[21]);
*(vp++) = U00*I0[22] + U40*I1[22]
           + (threeo2z)*(I2[22] - (lpoz)*I3[22]);
*(vp++) = U00*I0[23] + U40*I1[23]
           + (threeo2z)*(I2[23] - (lpoz)*I3[23]);
*(vp++) = U00*I0[24] + U40*I1[24]
           + (threeo2z)*(I2[24] - (lpoz)*I3[24]);
*(vp++) = U00*I0[25] + U40*I1[25]
           + (threeo2z)*(I2[25] - (lpoz)*I3[25]);
*(vp++) = U00*I0[26] + U40*I1[26]
           + (threeo2z)*(I2[26] - (lpoz)*I3[26]);
*(vp++) = U00*I0[27] + U40*I1[27]
           + (threeo2z)*(I2[27] - (lpoz)*I3[27]);
*(vp++) = U00*I0[28] + U40*I1[28]
           + (twoo2z)*(I2[28] - (lpoz)*I3[28])
           + (sixo2zn)*I4[21];
*(vp++) = U00*I0[29] + U40*I1[29]
           + (twoo2z)*(I2[29] - (lpoz)*I3[29])
           + (fiveo2zn)*I4[22];
*(vp++) = U00*I0[30] + U40*I1[30]
           + (twoo2z)*(I2[30] - (lpoz)*I3[30])
           + (fiveo2zn)*I4[23];
*(vp++) = U00*I0[31] + U40*I1[31]
           + (twoo2z)*(I2[31] - (lpoz)*I3[31])
           + (fouro2zn)*I4[24];
*(vp++) = U00*I0[32] + U40*I1[32]
           + (twoo2z)*(I2[32] - (lpoz)*I3[32])
           + (fouro2zn)*I4[25];
*(vp++) = U00*I0[33] + U40*I1[33]
           + (twoo2z)*(I2[33] - (lpoz)*I3[33])
           + (fouro2zn)*I4[26];
*(vp++) = U00*I0[34] + U40*I1[34]
           + (twoo2z)*(I2[34] - (lpoz)*I3[34])
           + (threeo2zn)*I4[27];
*(vp++) = U00*I0[35] + U40*I1[35]
           + (twoo2z)*(I2[35] - (lpoz)*I3[35])
           + (threeo2zn)*I4[28];
*(vp++) = U00*I0[36] + U40*I1[36]
           + (twoo2z)*(I2[36] - (lpoz)*I3[36])
           + (threeo2zn)*I4[29];
*(vp++) = U00*I0[37] + U40*I1[37]
           + (twoo2z)*(I2[37] - (lpoz)*I3[37])
           + (threeo2zn)*I4[30];
*(vp++) = U00*I0[38] + U40*I1[38]
           + (twoo2z)*(I2[38] - (lpoz)*I3[38])
           + (twoo2zn)*I4[31];
*(vp++) = U00*I0[39] + U40*I1[39]
           + (twoo2z)*(I2[39] - (lpoz)*I3[39])
           + (twoo2zn)*I4[32];
*(vp++) = U00*I0[40] + U40*I1[40]
           + (twoo2z)*(I2[40] - (lpoz)*I3[40])
           + (twoo2zn)*I4[33];
*(vp++) = U00*I0[41] + U40*I1[41]
           + (twoo2z)*(I2[41] - (lpoz)*I3[41])
           + (twoo2zn)*I4[34];
*(vp++) = U00*I0[42] + U40*I1[42]
           + (twoo2z)*(I2[42] - (lpoz)*I3[42])
           + (twoo2zn)*I4[35];
*(vp++) = U00*I0[43] + U40*I1[43]
           + (twoo2z)*(I2[43] - (lpoz)*I3[43])
           + (oneo2zn)*I4[36];
*(vp++) = U00*I0[44] + U40*I1[44]
           + (twoo2z)*(I2[44] - (lpoz)*I3[44])
           + (oneo2zn)*I4[37];
*(vp++) = U00*I0[45] + U40*I1[45]
           + (twoo2z)*(I2[45] - (lpoz)*I3[45])
           + (oneo2zn)*I4[38];
*(vp++) = U00*I0[46] + U40*I1[46]
           + (twoo2z)*(I2[46] - (lpoz)*I3[46])
           + (oneo2zn)*I4[39];
*(vp++) = U00*I0[47] + U40*I1[47]
           + (twoo2z)*(I2[47] - (lpoz)*I3[47])
           + (oneo2zn)*I4[40];
*(vp++) = U00*I0[48] + U40*I1[48]
           + (twoo2z)*(I2[48] - (lpoz)*I3[48])
           + (oneo2zn)*I4[41];
*(vp++) = U00*I0[49] + U40*I1[49]
           + (twoo2z)*(I2[49] - (lpoz)*I3[49]);
*(vp++) = U00*I0[50] + U40*I1[50]
           + (twoo2z)*(I2[50] - (lpoz)*I3[50]);
*(vp++) = U00*I0[51] + U40*I1[51]
           + (twoo2z)*(I2[51] - (lpoz)*I3[51]);
*(vp++) = U00*I0[52] + U40*I1[52]
           + (twoo2z)*(I2[52] - (lpoz)*I3[52]);
*(vp++) = U00*I0[53] + U40*I1[53]
           + (twoo2z)*(I2[53] - (lpoz)*I3[53]);
*(vp++) = U00*I0[54] + U40*I1[54]
           + (twoo2z)*(I2[54] - (lpoz)*I3[54]);
*(vp++) = U00*I0[55] + U40*I1[55]
           + (twoo2z)*(I2[55] - (lpoz)*I3[55]);
*(vp++) = U00*I0[56] + U40*I1[56]
           + (twoo2z)*(I2[56] - (lpoz)*I3[56])
           + (sixo2zn)*I4[42];
*(vp++) = U00*I0[57] + U40*I1[57]
           + (twoo2z)*(I2[57] - (lpoz)*I3[57])
           + (fiveo2zn)*I4[43];
*(vp++) = U00*I0[58] + U40*I1[58]
           + (twoo2z)*(I2[58] - (lpoz)*I3[58])
           + (fiveo2zn)*I4[44];
*(vp++) = U00*I0[59] + U40*I1[59]
           + (twoo2z)*(I2[59] - (lpoz)*I3[59])
           + (fouro2zn)*I4[45];
*(vp++) = U00*I0[60] + U40*I1[60]
           + (twoo2z)*(I2[60] - (lpoz)*I3[60])
           + (fouro2zn)*I4[46];
*(vp++) = U00*I0[61] + U40*I1[61]
           + (twoo2z)*(I2[61] - (lpoz)*I3[61])
           + (fouro2zn)*I4[47];
*(vp++) = U00*I0[62] + U40*I1[62]
           + (twoo2z)*(I2[62] - (lpoz)*I3[62])
           + (threeo2zn)*I4[48];
*(vp++) = U00*I0[63] + U40*I1[63]
           + (twoo2z)*(I2[63] - (lpoz)*I3[63])
           + (threeo2zn)*I4[49];
*(vp++) = U00*I0[64] + U40*I1[64]
           + (twoo2z)*(I2[64] - (lpoz)*I3[64])
           + (threeo2zn)*I4[50];
*(vp++) = U00*I0[65] + U40*I1[65]
           + (twoo2z)*(I2[65] - (lpoz)*I3[65])
           + (threeo2zn)*I4[51];
*(vp++) = U00*I0[66] + U40*I1[66]
           + (twoo2z)*(I2[66] - (lpoz)*I3[66])
           + (twoo2zn)*I4[52];
*(vp++) = U00*I0[67] + U40*I1[67]
           + (twoo2z)*(I2[67] - (lpoz)*I3[67])
           + (twoo2zn)*I4[53];
*(vp++) = U00*I0[68] + U40*I1[68]
           + (twoo2z)*(I2[68] - (lpoz)*I3[68])
           + (twoo2zn)*I4[54];
*(vp++) = U00*I0[69] + U40*I1[69]
           + (twoo2z)*(I2[69] - (lpoz)*I3[69])
           + (twoo2zn)*I4[55];
*(vp++) = U00*I0[70] + U40*I1[70]
           + (twoo2z)*(I2[70] - (lpoz)*I3[70])
           + (twoo2zn)*I4[56];
*(vp++) = U00*I0[71] + U40*I1[71]
           + (twoo2z)*(I2[71] - (lpoz)*I3[71])
           + (oneo2zn)*I4[57];
*(vp++) = U00*I0[72] + U40*I1[72]
           + (twoo2z)*(I2[72] - (lpoz)*I3[72])
           + (oneo2zn)*I4[58];
*(vp++) = U00*I0[73] + U40*I1[73]
           + (twoo2z)*(I2[73] - (lpoz)*I3[73])
           + (oneo2zn)*I4[59];
*(vp++) = U00*I0[74] + U40*I1[74]
           + (twoo2z)*(I2[74] - (lpoz)*I3[74])
           + (oneo2zn)*I4[60];
*(vp++) = U00*I0[75] + U40*I1[75]
           + (twoo2z)*(I2[75] - (lpoz)*I3[75])
           + (oneo2zn)*I4[61];
*(vp++) = U00*I0[76] + U40*I1[76]
           + (twoo2z)*(I2[76] - (lpoz)*I3[76])
           + (oneo2zn)*I4[62];
*(vp++) = U00*I0[77] + U40*I1[77]
           + (twoo2z)*(I2[77] - (lpoz)*I3[77]);
*(vp++) = U00*I0[78] + U40*I1[78]
           + (twoo2z)*(I2[78] - (lpoz)*I3[78]);
*(vp++) = U00*I0[79] + U40*I1[79]
           + (twoo2z)*(I2[79] - (lpoz)*I3[79]);
*(vp++) = U00*I0[80] + U40*I1[80]
           + (twoo2z)*(I2[80] - (lpoz)*I3[80]);
*(vp++) = U00*I0[81] + U40*I1[81]
           + (twoo2z)*(I2[81] - (lpoz)*I3[81]);
*(vp++) = U00*I0[82] + U40*I1[82]
           + (twoo2z)*(I2[82] - (lpoz)*I3[82]);
*(vp++) = U00*I0[83] + U40*I1[83]
           + (twoo2z)*(I2[83] - (lpoz)*I3[83]);
*(vp++) = U00*I0[84] + U40*I1[84]
           + (oneo2z)*(I2[84] - (lpoz)*I3[84])
           + (sixo2zn)*I4[63];
*(vp++) = U00*I0[85] + U40*I1[85]
           + (oneo2z)*(I2[85] - (lpoz)*I3[85])
           + (fiveo2zn)*I4[64];
*(vp++) = U00*I0[86] + U40*I1[86]
           + (oneo2z)*(I2[86] - (lpoz)*I3[86])
           + (fiveo2zn)*I4[65];
*(vp++) = U00*I0[87] + U40*I1[87]
           + (oneo2z)*(I2[87] - (lpoz)*I3[87])
           + (fouro2zn)*I4[66];
*(vp++) = U00*I0[88] + U40*I1[88]
           + (oneo2z)*(I2[88] - (lpoz)*I3[88])
           + (fouro2zn)*I4[67];
*(vp++) = U00*I0[89] + U40*I1[89]
           + (oneo2z)*(I2[89] - (lpoz)*I3[89])
           + (fouro2zn)*I4[68];
*(vp++) = U00*I0[90] + U40*I1[90]
           + (oneo2z)*(I2[90] - (lpoz)*I3[90])
           + (threeo2zn)*I4[69];
*(vp++) = U00*I0[91] + U40*I1[91]
           + (oneo2z)*(I2[91] - (lpoz)*I3[91])
           + (threeo2zn)*I4[70];
*(vp++) = U00*I0[92] + U40*I1[92]
           + (oneo2z)*(I2[92] - (lpoz)*I3[92])
           + (threeo2zn)*I4[71];
*(vp++) = U00*I0[93] + U40*I1[93]
           + (oneo2z)*(I2[93] - (lpoz)*I3[93])
           + (threeo2zn)*I4[72];
*(vp++) = U00*I0[94] + U40*I1[94]
           + (oneo2z)*(I2[94] - (lpoz)*I3[94])
           + (twoo2zn)*I4[73];
*(vp++) = U00*I0[95] + U40*I1[95]
           + (oneo2z)*(I2[95] - (lpoz)*I3[95])
           + (twoo2zn)*I4[74];
*(vp++) = U00*I0[96] + U40*I1[96]
           + (oneo2z)*(I2[96] - (lpoz)*I3[96])
           + (twoo2zn)*I4[75];
*(vp++) = U00*I0[97] + U40*I1[97]
           + (oneo2z)*(I2[97] - (lpoz)*I3[97])
           + (twoo2zn)*I4[76];
*(vp++) = U00*I0[98] + U40*I1[98]
           + (oneo2z)*(I2[98] - (lpoz)*I3[98])
           + (twoo2zn)*I4[77];
*(vp++) = U00*I0[99] + U40*I1[99]
           + (oneo2z)*(I2[99] - (lpoz)*I3[99])
           + (oneo2zn)*I4[78];
*(vp++) = U00*I0[100] + U40*I1[100]
           + (oneo2z)*(I2[100] - (lpoz)*I3[100])
           + (oneo2zn)*I4[79];
*(vp++) = U00*I0[101] + U40*I1[101]
           + (oneo2z)*(I2[101] - (lpoz)*I3[101])
           + (oneo2zn)*I4[80];
*(vp++) = U00*I0[102] + U40*I1[102]
           + (oneo2z)*(I2[102] - (lpoz)*I3[102])
           + (oneo2zn)*I4[81];
*(vp++) = U00*I0[103] + U40*I1[103]
           + (oneo2z)*(I2[103] - (lpoz)*I3[103])
           + (oneo2zn)*I4[82];
*(vp++) = U00*I0[104] + U40*I1[104]
           + (oneo2z)*(I2[104] - (lpoz)*I3[104])
           + (oneo2zn)*I4[83];
*(vp++) = U00*I0[105] + U40*I1[105]
           + (oneo2z)*(I2[105] - (lpoz)*I3[105]);
*(vp++) = U00*I0[106] + U40*I1[106]
           + (oneo2z)*(I2[106] - (lpoz)*I3[106]);
*(vp++) = U00*I0[107] + U40*I1[107]
           + (oneo2z)*(I2[107] - (lpoz)*I3[107]);
*(vp++) = U00*I0[108] + U40*I1[108]
           + (oneo2z)*(I2[108] - (lpoz)*I3[108]);
*(vp++) = U00*I0[109] + U40*I1[109]
           + (oneo2z)*(I2[109] - (lpoz)*I3[109]);
*(vp++) = U00*I0[110] + U40*I1[110]
           + (oneo2z)*(I2[110] - (lpoz)*I3[110]);
*(vp++) = U00*I0[111] + U40*I1[111]
           + (oneo2z)*(I2[111] - (lpoz)*I3[111]);
*(vp++) = U00*I0[112] + U40*I1[112]
           + (oneo2z)*(I2[112] - (lpoz)*I3[112])
           + (sixo2zn)*I4[84];
*(vp++) = U00*I0[113] + U40*I1[113]
           + (oneo2z)*(I2[113] - (lpoz)*I3[113])
           + (fiveo2zn)*I4[85];
*(vp++) = U00*I0[114] + U40*I1[114]
           + (oneo2z)*(I2[114] - (lpoz)*I3[114])
           + (fiveo2zn)*I4[86];
*(vp++) = U00*I0[115] + U40*I1[115]
           + (oneo2z)*(I2[115] - (lpoz)*I3[115])
           + (fouro2zn)*I4[87];
*(vp++) = U00*I0[116] + U40*I1[116]
           + (oneo2z)*(I2[116] - (lpoz)*I3[116])
           + (fouro2zn)*I4[88];
*(vp++) = U00*I0[117] + U40*I1[117]
           + (oneo2z)*(I2[117] - (lpoz)*I3[117])
           + (fouro2zn)*I4[89];
*(vp++) = U00*I0[118] + U40*I1[118]
           + (oneo2z)*(I2[118] - (lpoz)*I3[118])
           + (threeo2zn)*I4[90];
*(vp++) = U00*I0[119] + U40*I1[119]
           + (oneo2z)*(I2[119] - (lpoz)*I3[119])
           + (threeo2zn)*I4[91];
*(vp++) = U00*I0[120] + U40*I1[120]
           + (oneo2z)*(I2[120] - (lpoz)*I3[120])
           + (threeo2zn)*I4[92];
*(vp++) = U00*I0[121] + U40*I1[121]
           + (oneo2z)*(I2[121] - (lpoz)*I3[121])
           + (threeo2zn)*I4[93];
*(vp++) = U00*I0[122] + U40*I1[122]
           + (oneo2z)*(I2[122] - (lpoz)*I3[122])
           + (twoo2zn)*I4[94];
*(vp++) = U00*I0[123] + U40*I1[123]
           + (oneo2z)*(I2[123] - (lpoz)*I3[123])
           + (twoo2zn)*I4[95];
*(vp++) = U00*I0[124] + U40*I1[124]
           + (oneo2z)*(I2[124] - (lpoz)*I3[124])
           + (twoo2zn)*I4[96];
*(vp++) = U00*I0[125] + U40*I1[125]
           + (oneo2z)*(I2[125] - (lpoz)*I3[125])
           + (twoo2zn)*I4[97];
*(vp++) = U00*I0[126] + U40*I1[126]
           + (oneo2z)*(I2[126] - (lpoz)*I3[126])
           + (twoo2zn)*I4[98];
*(vp++) = U00*I0[127] + U40*I1[127]
           + (oneo2z)*(I2[127] - (lpoz)*I3[127])
           + (oneo2zn)*I4[99];
*(vp++) = U00*I0[128] + U40*I1[128]
           + (oneo2z)*(I2[128] - (lpoz)*I3[128])
           + (oneo2zn)*I4[100];
*(vp++) = U00*I0[129] + U40*I1[129]
           + (oneo2z)*(I2[129] - (lpoz)*I3[129])
           + (oneo2zn)*I4[101];
*(vp++) = U00*I0[130] + U40*I1[130]
           + (oneo2z)*(I2[130] - (lpoz)*I3[130])
           + (oneo2zn)*I4[102];
*(vp++) = U00*I0[131] + U40*I1[131]
           + (oneo2z)*(I2[131] - (lpoz)*I3[131])
           + (oneo2zn)*I4[103];
*(vp++) = U00*I0[132] + U40*I1[132]
           + (oneo2z)*(I2[132] - (lpoz)*I3[132])
           + (oneo2zn)*I4[104];
*(vp++) = U00*I0[133] + U40*I1[133]
           + (oneo2z)*(I2[133] - (lpoz)*I3[133]);
*(vp++) = U00*I0[134] + U40*I1[134]
           + (oneo2z)*(I2[134] - (lpoz)*I3[134]);
*(vp++) = U00*I0[135] + U40*I1[135]
           + (oneo2z)*(I2[135] - (lpoz)*I3[135]);
*(vp++) = U00*I0[136] + U40*I1[136]
           + (oneo2z)*(I2[136] - (lpoz)*I3[136]);
*(vp++) = U00*I0[137] + U40*I1[137]
           + (oneo2z)*(I2[137] - (lpoz)*I3[137]);
*(vp++) = U00*I0[138] + U40*I1[138]
           + (oneo2z)*(I2[138] - (lpoz)*I3[138]);
*(vp++) = U00*I0[139] + U40*I1[139]
           + (oneo2z)*(I2[139] - (lpoz)*I3[139]);
*(vp++) = U00*I0[140] + U40*I1[140]
           + (oneo2z)*(I2[140] - (lpoz)*I3[140])
           + (sixo2zn)*I4[105];
*(vp++) = U00*I0[141] + U40*I1[141]
           + (oneo2z)*(I2[141] - (lpoz)*I3[141])
           + (fiveo2zn)*I4[106];
*(vp++) = U00*I0[142] + U40*I1[142]
           + (oneo2z)*(I2[142] - (lpoz)*I3[142])
           + (fiveo2zn)*I4[107];
*(vp++) = U00*I0[143] + U40*I1[143]
           + (oneo2z)*(I2[143] - (lpoz)*I3[143])
           + (fouro2zn)*I4[108];
*(vp++) = U00*I0[144] + U40*I1[144]
           + (oneo2z)*(I2[144] - (lpoz)*I3[144])
           + (fouro2zn)*I4[109];
*(vp++) = U00*I0[145] + U40*I1[145]
           + (oneo2z)*(I2[145] - (lpoz)*I3[145])
           + (fouro2zn)*I4[110];
*(vp++) = U00*I0[146] + U40*I1[146]
           + (oneo2z)*(I2[146] - (lpoz)*I3[146])
           + (threeo2zn)*I4[111];
*(vp++) = U00*I0[147] + U40*I1[147]
           + (oneo2z)*(I2[147] - (lpoz)*I3[147])
           + (threeo2zn)*I4[112];
*(vp++) = U00*I0[148] + U40*I1[148]
           + (oneo2z)*(I2[148] - (lpoz)*I3[148])
           + (threeo2zn)*I4[113];
*(vp++) = U00*I0[149] + U40*I1[149]
           + (oneo2z)*(I2[149] - (lpoz)*I3[149])
           + (threeo2zn)*I4[114];
*(vp++) = U00*I0[150] + U40*I1[150]
           + (oneo2z)*(I2[150] - (lpoz)*I3[150])
           + (twoo2zn)*I4[115];
*(vp++) = U00*I0[151] + U40*I1[151]
           + (oneo2z)*(I2[151] - (lpoz)*I3[151])
           + (twoo2zn)*I4[116];
*(vp++) = U00*I0[152] + U40*I1[152]
           + (oneo2z)*(I2[152] - (lpoz)*I3[152])
           + (twoo2zn)*I4[117];
*(vp++) = U00*I0[153] + U40*I1[153]
           + (oneo2z)*(I2[153] - (lpoz)*I3[153])
           + (twoo2zn)*I4[118];
*(vp++) = U00*I0[154] + U40*I1[154]
           + (oneo2z)*(I2[154] - (lpoz)*I3[154])
           + (twoo2zn)*I4[119];
*(vp++) = U00*I0[155] + U40*I1[155]
           + (oneo2z)*(I2[155] - (lpoz)*I3[155])
           + (oneo2zn)*I4[120];
*(vp++) = U00*I0[156] + U40*I1[156]
           + (oneo2z)*(I2[156] - (lpoz)*I3[156])
           + (oneo2zn)*I4[121];
*(vp++) = U00*I0[157] + U40*I1[157]
           + (oneo2z)*(I2[157] - (lpoz)*I3[157])
           + (oneo2zn)*I4[122];
*(vp++) = U00*I0[158] + U40*I1[158]
           + (oneo2z)*(I2[158] - (lpoz)*I3[158])
           + (oneo2zn)*I4[123];
*(vp++) = U00*I0[159] + U40*I1[159]
           + (oneo2z)*(I2[159] - (lpoz)*I3[159])
           + (oneo2zn)*I4[124];
*(vp++) = U00*I0[160] + U40*I1[160]
           + (oneo2z)*(I2[160] - (lpoz)*I3[160])
           + (oneo2zn)*I4[125];
*(vp++) = U00*I0[161] + U40*I1[161]
           + (oneo2z)*(I2[161] - (lpoz)*I3[161]);
*(vp++) = U00*I0[162] + U40*I1[162]
           + (oneo2z)*(I2[162] - (lpoz)*I3[162]);
*(vp++) = U00*I0[163] + U40*I1[163]
           + (oneo2z)*(I2[163] - (lpoz)*I3[163]);
*(vp++) = U00*I0[164] + U40*I1[164]
           + (oneo2z)*(I2[164] - (lpoz)*I3[164]);
*(vp++) = U00*I0[165] + U40*I1[165]
           + (oneo2z)*(I2[165] - (lpoz)*I3[165]);
*(vp++) = U00*I0[166] + U40*I1[166]
           + (oneo2z)*(I2[166] - (lpoz)*I3[166]);
*(vp++) = U00*I0[167] + U40*I1[167]
           + (oneo2z)*(I2[167] - (lpoz)*I3[167]);
*(vp++) = U00*I0[168] + U40*I1[168]
           + (sixo2zn)*I4[126];
*(vp++) = U00*I0[169] + U40*I1[169]
           + (fiveo2zn)*I4[127];
*(vp++) = U00*I0[170] + U40*I1[170]
           + (fiveo2zn)*I4[128];
*(vp++) = U00*I0[171] + U40*I1[171]
           + (fouro2zn)*I4[129];
*(vp++) = U00*I0[172] + U40*I1[172]
           + (fouro2zn)*I4[130];
*(vp++) = U00*I0[173] + U40*I1[173]
           + (fouro2zn)*I4[131];
*(vp++) = U00*I0[174] + U40*I1[174]
           + (threeo2zn)*I4[132];
*(vp++) = U00*I0[175] + U40*I1[175]
           + (threeo2zn)*I4[133];
*(vp++) = U00*I0[176] + U40*I1[176]
           + (threeo2zn)*I4[134];
*(vp++) = U00*I0[177] + U40*I1[177]
           + (threeo2zn)*I4[135];
*(vp++) = U00*I0[178] + U40*I1[178]
           + (twoo2zn)*I4[136];
*(vp++) = U00*I0[179] + U40*I1[179]
           + (twoo2zn)*I4[137];
*(vp++) = U00*I0[180] + U40*I1[180]
           + (twoo2zn)*I4[138];
*(vp++) = U00*I0[181] + U40*I1[181]
           + (twoo2zn)*I4[139];
*(vp++) = U00*I0[182] + U40*I1[182]
           + (twoo2zn)*I4[140];
*(vp++) = U00*I0[183] + U40*I1[183]
           + (oneo2zn)*I4[141];
*(vp++) = U00*I0[184] + U40*I1[184]
           + (oneo2zn)*I4[142];
*(vp++) = U00*I0[185] + U40*I1[185]
           + (oneo2zn)*I4[143];
*(vp++) = U00*I0[186] + U40*I1[186]
           + (oneo2zn)*I4[144];
*(vp++) = U00*I0[187] + U40*I1[187]
           + (oneo2zn)*I4[145];
*(vp++) = U00*I0[188] + U40*I1[188]
           + (oneo2zn)*I4[146];
*(vp++) = U00*I0[189] + U40*I1[189];
*(vp++) = U00*I0[190] + U40*I1[190];
*(vp++) = U00*I0[191] + U40*I1[191];
*(vp++) = U00*I0[192] + U40*I1[192];
*(vp++) = U00*I0[193] + U40*I1[193];
*(vp++) = U00*I0[194] + U40*I1[194];
*(vp++) = U00*I0[195] + U40*I1[195];
*(vp++) = U00*I0[196] + U40*I1[196]
           + (sixo2zn)*I4[147];
*(vp++) = U00*I0[197] + U40*I1[197]
           + (fiveo2zn)*I4[148];
*(vp++) = U00*I0[198] + U40*I1[198]
           + (fiveo2zn)*I4[149];
*(vp++) = U00*I0[199] + U40*I1[199]
           + (fouro2zn)*I4[150];
*(vp++) = U00*I0[200] + U40*I1[200]
           + (fouro2zn)*I4[151];
*(vp++) = U00*I0[201] + U40*I1[201]
           + (fouro2zn)*I4[152];
*(vp++) = U00*I0[202] + U40*I1[202]
           + (threeo2zn)*I4[153];
*(vp++) = U00*I0[203] + U40*I1[203]
           + (threeo2zn)*I4[154];
*(vp++) = U00*I0[204] + U40*I1[204]
           + (threeo2zn)*I4[155];
*(vp++) = U00*I0[205] + U40*I1[205]
           + (threeo2zn)*I4[156];
*(vp++) = U00*I0[206] + U40*I1[206]
           + (twoo2zn)*I4[157];
*(vp++) = U00*I0[207] + U40*I1[207]
           + (twoo2zn)*I4[158];
*(vp++) = U00*I0[208] + U40*I1[208]
           + (twoo2zn)*I4[159];
*(vp++) = U00*I0[209] + U40*I1[209]
           + (twoo2zn)*I4[160];
*(vp++) = U00*I0[210] + U40*I1[210]
           + (twoo2zn)*I4[161];
return vp;
}

REALTYPE *_build_g0i0_1(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
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
  REALTYPE twoo2z;
  REALTYPE threeo2z;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  threeo2zn = 3.0*Data->oo2zn;
  fouro2zn = 4.0*Data->oo2zn;
  fiveo2zn = 5.0*Data->oo2zn;
  sixo2zn = 6.0*Data->oo2zn;
  oneo2z = 1.0*Data->oo2z;
  twoo2z = 2.0*Data->oo2z;
  threeo2z = 3.0*Data->oo2z;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[211] + U40*I1[211]
           + (oneo2zn)*I4[162];
*(vp++) = U00*I0[212] + U40*I1[212]
           + (oneo2zn)*I4[163];
*(vp++) = U00*I0[213] + U40*I1[213]
           + (oneo2zn)*I4[164];
*(vp++) = U00*I0[214] + U40*I1[214]
           + (oneo2zn)*I4[165];
*(vp++) = U00*I0[215] + U40*I1[215]
           + (oneo2zn)*I4[166];
*(vp++) = U00*I0[216] + U40*I1[216]
           + (oneo2zn)*I4[167];
*(vp++) = U00*I0[217] + U40*I1[217];
*(vp++) = U00*I0[218] + U40*I1[218];
*(vp++) = U00*I0[219] + U40*I1[219];
*(vp++) = U00*I0[220] + U40*I1[220];
*(vp++) = U00*I0[221] + U40*I1[221];
*(vp++) = U00*I0[222] + U40*I1[222];
*(vp++) = U00*I0[223] + U40*I1[223];
*(vp++) = U00*I0[224] + U40*I1[224]
           + (sixo2zn)*I4[168];
*(vp++) = U00*I0[225] + U40*I1[225]
           + (fiveo2zn)*I4[169];
*(vp++) = U00*I0[226] + U40*I1[226]
           + (fiveo2zn)*I4[170];
*(vp++) = U00*I0[227] + U40*I1[227]
           + (fouro2zn)*I4[171];
*(vp++) = U00*I0[228] + U40*I1[228]
           + (fouro2zn)*I4[172];
*(vp++) = U00*I0[229] + U40*I1[229]
           + (fouro2zn)*I4[173];
*(vp++) = U00*I0[230] + U40*I1[230]
           + (threeo2zn)*I4[174];
*(vp++) = U00*I0[231] + U40*I1[231]
           + (threeo2zn)*I4[175];
*(vp++) = U00*I0[232] + U40*I1[232]
           + (threeo2zn)*I4[176];
*(vp++) = U00*I0[233] + U40*I1[233]
           + (threeo2zn)*I4[177];
*(vp++) = U00*I0[234] + U40*I1[234]
           + (twoo2zn)*I4[178];
*(vp++) = U00*I0[235] + U40*I1[235]
           + (twoo2zn)*I4[179];
*(vp++) = U00*I0[236] + U40*I1[236]
           + (twoo2zn)*I4[180];
*(vp++) = U00*I0[237] + U40*I1[237]
           + (twoo2zn)*I4[181];
*(vp++) = U00*I0[238] + U40*I1[238]
           + (twoo2zn)*I4[182];
*(vp++) = U00*I0[239] + U40*I1[239]
           + (oneo2zn)*I4[183];
*(vp++) = U00*I0[240] + U40*I1[240]
           + (oneo2zn)*I4[184];
*(vp++) = U00*I0[241] + U40*I1[241]
           + (oneo2zn)*I4[185];
*(vp++) = U00*I0[242] + U40*I1[242]
           + (oneo2zn)*I4[186];
*(vp++) = U00*I0[243] + U40*I1[243]
           + (oneo2zn)*I4[187];
*(vp++) = U00*I0[244] + U40*I1[244]
           + (oneo2zn)*I4[188];
*(vp++) = U00*I0[245] + U40*I1[245];
*(vp++) = U00*I0[246] + U40*I1[246];
*(vp++) = U00*I0[247] + U40*I1[247];
*(vp++) = U00*I0[248] + U40*I1[248];
*(vp++) = U00*I0[249] + U40*I1[249];
*(vp++) = U00*I0[250] + U40*I1[250];
*(vp++) = U00*I0[251] + U40*I1[251];
*(vp++) = U00*I0[252] + U40*I1[252]
           + (sixo2zn)*I4[189];
*(vp++) = U00*I0[253] + U40*I1[253]
           + (fiveo2zn)*I4[190];
*(vp++) = U00*I0[254] + U40*I1[254]
           + (fiveo2zn)*I4[191];
*(vp++) = U00*I0[255] + U40*I1[255]
           + (fouro2zn)*I4[192];
*(vp++) = U00*I0[256] + U40*I1[256]
           + (fouro2zn)*I4[193];
*(vp++) = U00*I0[257] + U40*I1[257]
           + (fouro2zn)*I4[194];
*(vp++) = U00*I0[258] + U40*I1[258]
           + (threeo2zn)*I4[195];
*(vp++) = U00*I0[259] + U40*I1[259]
           + (threeo2zn)*I4[196];
*(vp++) = U00*I0[260] + U40*I1[260]
           + (threeo2zn)*I4[197];
*(vp++) = U00*I0[261] + U40*I1[261]
           + (threeo2zn)*I4[198];
*(vp++) = U00*I0[262] + U40*I1[262]
           + (twoo2zn)*I4[199];
*(vp++) = U00*I0[263] + U40*I1[263]
           + (twoo2zn)*I4[200];
*(vp++) = U00*I0[264] + U40*I1[264]
           + (twoo2zn)*I4[201];
*(vp++) = U00*I0[265] + U40*I1[265]
           + (twoo2zn)*I4[202];
*(vp++) = U00*I0[266] + U40*I1[266]
           + (twoo2zn)*I4[203];
*(vp++) = U00*I0[267] + U40*I1[267]
           + (oneo2zn)*I4[204];
*(vp++) = U00*I0[268] + U40*I1[268]
           + (oneo2zn)*I4[205];
*(vp++) = U00*I0[269] + U40*I1[269]
           + (oneo2zn)*I4[206];
*(vp++) = U00*I0[270] + U40*I1[270]
           + (oneo2zn)*I4[207];
*(vp++) = U00*I0[271] + U40*I1[271]
           + (oneo2zn)*I4[208];
*(vp++) = U00*I0[272] + U40*I1[272]
           + (oneo2zn)*I4[209];
*(vp++) = U00*I0[273] + U40*I1[273];
*(vp++) = U00*I0[274] + U40*I1[274];
*(vp++) = U00*I0[275] + U40*I1[275];
*(vp++) = U00*I0[276] + U40*I1[276];
*(vp++) = U00*I0[277] + U40*I1[277];
*(vp++) = U00*I0[278] + U40*I1[278];
*(vp++) = U00*I0[279] + U40*I1[279];
*(vp++) = U01*I0[168] + U41*I1[168]
           + (threeo2z)*(I2[84] - (lpoz)*I3[84]);
*(vp++) = U01*I0[169] + U41*I1[169]
           + (threeo2z)*(I2[85] - (lpoz)*I3[85])
           + (oneo2zn)*I4[126];
*(vp++) = U01*I0[170] + U41*I1[170]
           + (threeo2z)*(I2[86] - (lpoz)*I3[86]);
*(vp++) = U01*I0[171] + U41*I1[171]
           + (threeo2z)*(I2[87] - (lpoz)*I3[87])
           + (twoo2zn)*I4[127];
*(vp++) = U01*I0[172] + U41*I1[172]
           + (threeo2z)*(I2[88] - (lpoz)*I3[88])
           + (oneo2zn)*I4[128];
*(vp++) = U01*I0[173] + U41*I1[173]
           + (threeo2z)*(I2[89] - (lpoz)*I3[89]);
*(vp++) = U01*I0[174] + U41*I1[174]
           + (threeo2z)*(I2[90] - (lpoz)*I3[90])
           + (threeo2zn)*I4[129];
*(vp++) = U01*I0[175] + U41*I1[175]
           + (threeo2z)*(I2[91] - (lpoz)*I3[91])
           + (twoo2zn)*I4[130];
*(vp++) = U01*I0[176] + U41*I1[176]
           + (threeo2z)*(I2[92] - (lpoz)*I3[92])
           + (oneo2zn)*I4[131];
*(vp++) = U01*I0[177] + U41*I1[177]
           + (threeo2z)*(I2[93] - (lpoz)*I3[93]);
*(vp++) = U01*I0[178] + U41*I1[178]
           + (threeo2z)*(I2[94] - (lpoz)*I3[94])
           + (fouro2zn)*I4[132];
*(vp++) = U01*I0[179] + U41*I1[179]
           + (threeo2z)*(I2[95] - (lpoz)*I3[95])
           + (threeo2zn)*I4[133];
*(vp++) = U01*I0[180] + U41*I1[180]
           + (threeo2z)*(I2[96] - (lpoz)*I3[96])
           + (twoo2zn)*I4[134];
*(vp++) = U01*I0[181] + U41*I1[181]
           + (threeo2z)*(I2[97] - (lpoz)*I3[97])
           + (oneo2zn)*I4[135];
*(vp++) = U01*I0[182] + U41*I1[182]
           + (threeo2z)*(I2[98] - (lpoz)*I3[98]);
*(vp++) = U01*I0[183] + U41*I1[183]
           + (threeo2z)*(I2[99] - (lpoz)*I3[99])
           + (fiveo2zn)*I4[136];
*(vp++) = U01*I0[184] + U41*I1[184]
           + (threeo2z)*(I2[100] - (lpoz)*I3[100])
           + (fouro2zn)*I4[137];
*(vp++) = U01*I0[185] + U41*I1[185]
           + (threeo2z)*(I2[101] - (lpoz)*I3[101])
           + (threeo2zn)*I4[138];
*(vp++) = U01*I0[186] + U41*I1[186]
           + (threeo2z)*(I2[102] - (lpoz)*I3[102])
           + (twoo2zn)*I4[139];
*(vp++) = U01*I0[187] + U41*I1[187]
           + (threeo2z)*(I2[103] - (lpoz)*I3[103])
           + (oneo2zn)*I4[140];
*(vp++) = U01*I0[188] + U41*I1[188]
           + (threeo2z)*(I2[104] - (lpoz)*I3[104]);
*(vp++) = U01*I0[189] + U41*I1[189]
           + (threeo2z)*(I2[105] - (lpoz)*I3[105])
           + (sixo2zn)*I4[141];
*(vp++) = U01*I0[190] + U41*I1[190]
           + (threeo2z)*(I2[106] - (lpoz)*I3[106])
           + (fiveo2zn)*I4[142];
*(vp++) = U01*I0[191] + U41*I1[191]
           + (threeo2z)*(I2[107] - (lpoz)*I3[107])
           + (fouro2zn)*I4[143];
*(vp++) = U01*I0[192] + U41*I1[192]
           + (threeo2z)*(I2[108] - (lpoz)*I3[108])
           + (threeo2zn)*I4[144];
*(vp++) = U01*I0[193] + U41*I1[193]
           + (threeo2z)*(I2[109] - (lpoz)*I3[109])
           + (twoo2zn)*I4[145];
*(vp++) = U01*I0[194] + U41*I1[194]
           + (threeo2z)*(I2[110] - (lpoz)*I3[110])
           + (oneo2zn)*I4[146];
*(vp++) = U01*I0[195] + U41*I1[195]
           + (threeo2z)*(I2[111] - (lpoz)*I3[111]);
*(vp++) = U01*I0[196] + U41*I1[196]
           + (twoo2z)*(I2[112] - (lpoz)*I3[112]);
*(vp++) = U01*I0[197] + U41*I1[197]
           + (twoo2z)*(I2[113] - (lpoz)*I3[113])
           + (oneo2zn)*I4[147];
*(vp++) = U01*I0[198] + U41*I1[198]
           + (twoo2z)*(I2[114] - (lpoz)*I3[114]);
*(vp++) = U01*I0[199] + U41*I1[199]
           + (twoo2z)*(I2[115] - (lpoz)*I3[115])
           + (twoo2zn)*I4[148];
*(vp++) = U01*I0[200] + U41*I1[200]
           + (twoo2z)*(I2[116] - (lpoz)*I3[116])
           + (oneo2zn)*I4[149];
*(vp++) = U01*I0[201] + U41*I1[201]
           + (twoo2z)*(I2[117] - (lpoz)*I3[117]);
*(vp++) = U01*I0[202] + U41*I1[202]
           + (twoo2z)*(I2[118] - (lpoz)*I3[118])
           + (threeo2zn)*I4[150];
*(vp++) = U01*I0[203] + U41*I1[203]
           + (twoo2z)*(I2[119] - (lpoz)*I3[119])
           + (twoo2zn)*I4[151];
*(vp++) = U01*I0[204] + U41*I1[204]
           + (twoo2z)*(I2[120] - (lpoz)*I3[120])
           + (oneo2zn)*I4[152];
*(vp++) = U01*I0[205] + U41*I1[205]
           + (twoo2z)*(I2[121] - (lpoz)*I3[121]);
*(vp++) = U01*I0[206] + U41*I1[206]
           + (twoo2z)*(I2[122] - (lpoz)*I3[122])
           + (fouro2zn)*I4[153];
*(vp++) = U01*I0[207] + U41*I1[207]
           + (twoo2z)*(I2[123] - (lpoz)*I3[123])
           + (threeo2zn)*I4[154];
*(vp++) = U01*I0[208] + U41*I1[208]
           + (twoo2z)*(I2[124] - (lpoz)*I3[124])
           + (twoo2zn)*I4[155];
*(vp++) = U01*I0[209] + U41*I1[209]
           + (twoo2z)*(I2[125] - (lpoz)*I3[125])
           + (oneo2zn)*I4[156];
*(vp++) = U01*I0[210] + U41*I1[210]
           + (twoo2z)*(I2[126] - (lpoz)*I3[126]);
*(vp++) = U01*I0[211] + U41*I1[211]
           + (twoo2z)*(I2[127] - (lpoz)*I3[127])
           + (fiveo2zn)*I4[157];
*(vp++) = U01*I0[212] + U41*I1[212]
           + (twoo2z)*(I2[128] - (lpoz)*I3[128])
           + (fouro2zn)*I4[158];
*(vp++) = U01*I0[213] + U41*I1[213]
           + (twoo2z)*(I2[129] - (lpoz)*I3[129])
           + (threeo2zn)*I4[159];
*(vp++) = U01*I0[214] + U41*I1[214]
           + (twoo2z)*(I2[130] - (lpoz)*I3[130])
           + (twoo2zn)*I4[160];
*(vp++) = U01*I0[215] + U41*I1[215]
           + (twoo2z)*(I2[131] - (lpoz)*I3[131])
           + (oneo2zn)*I4[161];
*(vp++) = U01*I0[216] + U41*I1[216]
           + (twoo2z)*(I2[132] - (lpoz)*I3[132]);
*(vp++) = U01*I0[217] + U41*I1[217]
           + (twoo2z)*(I2[133] - (lpoz)*I3[133])
           + (sixo2zn)*I4[162];
*(vp++) = U01*I0[218] + U41*I1[218]
           + (twoo2z)*(I2[134] - (lpoz)*I3[134])
           + (fiveo2zn)*I4[163];
*(vp++) = U01*I0[219] + U41*I1[219]
           + (twoo2z)*(I2[135] - (lpoz)*I3[135])
           + (fouro2zn)*I4[164];
*(vp++) = U01*I0[220] + U41*I1[220]
           + (twoo2z)*(I2[136] - (lpoz)*I3[136])
           + (threeo2zn)*I4[165];
*(vp++) = U01*I0[221] + U41*I1[221]
           + (twoo2z)*(I2[137] - (lpoz)*I3[137])
           + (twoo2zn)*I4[166];
*(vp++) = U01*I0[222] + U41*I1[222]
           + (twoo2z)*(I2[138] - (lpoz)*I3[138])
           + (oneo2zn)*I4[167];
*(vp++) = U01*I0[223] + U41*I1[223]
           + (twoo2z)*(I2[139] - (lpoz)*I3[139]);
*(vp++) = U01*I0[224] + U41*I1[224]
           + (oneo2z)*(I2[140] - (lpoz)*I3[140]);
*(vp++) = U01*I0[225] + U41*I1[225]
           + (oneo2z)*(I2[141] - (lpoz)*I3[141])
           + (oneo2zn)*I4[168];
*(vp++) = U01*I0[226] + U41*I1[226]
           + (oneo2z)*(I2[142] - (lpoz)*I3[142]);
*(vp++) = U01*I0[227] + U41*I1[227]
           + (oneo2z)*(I2[143] - (lpoz)*I3[143])
           + (twoo2zn)*I4[169];
*(vp++) = U01*I0[228] + U41*I1[228]
           + (oneo2z)*(I2[144] - (lpoz)*I3[144])
           + (oneo2zn)*I4[170];
*(vp++) = U01*I0[229] + U41*I1[229]
           + (oneo2z)*(I2[145] - (lpoz)*I3[145]);
*(vp++) = U01*I0[230] + U41*I1[230]
           + (oneo2z)*(I2[146] - (lpoz)*I3[146])
           + (threeo2zn)*I4[171];
*(vp++) = U01*I0[231] + U41*I1[231]
           + (oneo2z)*(I2[147] - (lpoz)*I3[147])
           + (twoo2zn)*I4[172];
*(vp++) = U01*I0[232] + U41*I1[232]
           + (oneo2z)*(I2[148] - (lpoz)*I3[148])
           + (oneo2zn)*I4[173];
*(vp++) = U01*I0[233] + U41*I1[233]
           + (oneo2z)*(I2[149] - (lpoz)*I3[149]);
*(vp++) = U01*I0[234] + U41*I1[234]
           + (oneo2z)*(I2[150] - (lpoz)*I3[150])
           + (fouro2zn)*I4[174];
*(vp++) = U01*I0[235] + U41*I1[235]
           + (oneo2z)*(I2[151] - (lpoz)*I3[151])
           + (threeo2zn)*I4[175];
*(vp++) = U01*I0[236] + U41*I1[236]
           + (oneo2z)*(I2[152] - (lpoz)*I3[152])
           + (twoo2zn)*I4[176];
*(vp++) = U01*I0[237] + U41*I1[237]
           + (oneo2z)*(I2[153] - (lpoz)*I3[153])
           + (oneo2zn)*I4[177];
*(vp++) = U01*I0[238] + U41*I1[238]
           + (oneo2z)*(I2[154] - (lpoz)*I3[154]);
*(vp++) = U01*I0[239] + U41*I1[239]
           + (oneo2z)*(I2[155] - (lpoz)*I3[155])
           + (fiveo2zn)*I4[178];
*(vp++) = U01*I0[240] + U41*I1[240]
           + (oneo2z)*(I2[156] - (lpoz)*I3[156])
           + (fouro2zn)*I4[179];
*(vp++) = U01*I0[241] + U41*I1[241]
           + (oneo2z)*(I2[157] - (lpoz)*I3[157])
           + (threeo2zn)*I4[180];
*(vp++) = U01*I0[242] + U41*I1[242]
           + (oneo2z)*(I2[158] - (lpoz)*I3[158])
           + (twoo2zn)*I4[181];
*(vp++) = U01*I0[243] + U41*I1[243]
           + (oneo2z)*(I2[159] - (lpoz)*I3[159])
           + (oneo2zn)*I4[182];
*(vp++) = U01*I0[244] + U41*I1[244]
           + (oneo2z)*(I2[160] - (lpoz)*I3[160]);
*(vp++) = U01*I0[245] + U41*I1[245]
           + (oneo2z)*(I2[161] - (lpoz)*I3[161])
           + (sixo2zn)*I4[183];
*(vp++) = U01*I0[246] + U41*I1[246]
           + (oneo2z)*(I2[162] - (lpoz)*I3[162])
           + (fiveo2zn)*I4[184];
*(vp++) = U01*I0[247] + U41*I1[247]
           + (oneo2z)*(I2[163] - (lpoz)*I3[163])
           + (fouro2zn)*I4[185];
*(vp++) = U01*I0[248] + U41*I1[248]
           + (oneo2z)*(I2[164] - (lpoz)*I3[164])
           + (threeo2zn)*I4[186];
*(vp++) = U01*I0[249] + U41*I1[249]
           + (oneo2z)*(I2[165] - (lpoz)*I3[165])
           + (twoo2zn)*I4[187];
*(vp++) = U01*I0[250] + U41*I1[250]
           + (oneo2z)*(I2[166] - (lpoz)*I3[166])
           + (oneo2zn)*I4[188];
*(vp++) = U01*I0[251] + U41*I1[251]
           + (oneo2z)*(I2[167] - (lpoz)*I3[167]);
*(vp++) = U01*I0[252] + U41*I1[252];
*(vp++) = U01*I0[253] + U41*I1[253]
           + (oneo2zn)*I4[189];
*(vp++) = U01*I0[254] + U41*I1[254];
*(vp++) = U01*I0[255] + U41*I1[255]
           + (twoo2zn)*I4[190];
*(vp++) = U01*I0[256] + U41*I1[256]
           + (oneo2zn)*I4[191];
*(vp++) = U01*I0[257] + U41*I1[257];
*(vp++) = U01*I0[258] + U41*I1[258]
           + (threeo2zn)*I4[192];
*(vp++) = U01*I0[259] + U41*I1[259]
           + (twoo2zn)*I4[193];
*(vp++) = U01*I0[260] + U41*I1[260]
           + (oneo2zn)*I4[194];
*(vp++) = U01*I0[261] + U41*I1[261];
*(vp++) = U01*I0[262] + U41*I1[262]
           + (fouro2zn)*I4[195];
*(vp++) = U01*I0[263] + U41*I1[263]
           + (threeo2zn)*I4[196];
*(vp++) = U01*I0[264] + U41*I1[264]
           + (twoo2zn)*I4[197];
*(vp++) = U01*I0[265] + U41*I1[265]
           + (oneo2zn)*I4[198];
*(vp++) = U01*I0[266] + U41*I1[266];
*(vp++) = U01*I0[267] + U41*I1[267]
           + (fiveo2zn)*I4[199];
*(vp++) = U01*I0[268] + U41*I1[268]
           + (fouro2zn)*I4[200];
*(vp++) = U01*I0[269] + U41*I1[269]
           + (threeo2zn)*I4[201];
*(vp++) = U01*I0[270] + U41*I1[270]
           + (twoo2zn)*I4[202];
*(vp++) = U01*I0[271] + U41*I1[271]
           + (oneo2zn)*I4[203];
*(vp++) = U01*I0[272] + U41*I1[272];
*(vp++) = U01*I0[273] + U41*I1[273]
           + (sixo2zn)*I4[204];
*(vp++) = U01*I0[274] + U41*I1[274]
           + (fiveo2zn)*I4[205];
*(vp++) = U01*I0[275] + U41*I1[275]
           + (fouro2zn)*I4[206];
*(vp++) = U01*I0[276] + U41*I1[276]
           + (threeo2zn)*I4[207];
*(vp++) = U01*I0[277] + U41*I1[277]
           + (twoo2zn)*I4[208];
*(vp++) = U01*I0[278] + U41*I1[278]
           + (oneo2zn)*I4[209];
*(vp++) = U01*I0[279] + U41*I1[279];
*(vp++) = U02*I0[252] + U42*I1[252]
           + (threeo2z)*(I2[140] - (lpoz)*I3[140]);
*(vp++) = U02*I0[253] + U42*I1[253]
           + (threeo2z)*(I2[141] - (lpoz)*I3[141]);
*(vp++) = U02*I0[254] + U42*I1[254]
           + (threeo2z)*(I2[142] - (lpoz)*I3[142])
           + (oneo2zn)*I4[189];
*(vp++) = U02*I0[255] + U42*I1[255]
           + (threeo2z)*(I2[143] - (lpoz)*I3[143]);
*(vp++) = U02*I0[256] + U42*I1[256]
           + (threeo2z)*(I2[144] - (lpoz)*I3[144])
           + (oneo2zn)*I4[190];
*(vp++) = U02*I0[257] + U42*I1[257]
           + (threeo2z)*(I2[145] - (lpoz)*I3[145])
           + (twoo2zn)*I4[191];
*(vp++) = U02*I0[258] + U42*I1[258]
           + (threeo2z)*(I2[146] - (lpoz)*I3[146]);
*(vp++) = U02*I0[259] + U42*I1[259]
           + (threeo2z)*(I2[147] - (lpoz)*I3[147])
           + (oneo2zn)*I4[192];
*(vp++) = U02*I0[260] + U42*I1[260]
           + (threeo2z)*(I2[148] - (lpoz)*I3[148])
           + (twoo2zn)*I4[193];
*(vp++) = U02*I0[261] + U42*I1[261]
           + (threeo2z)*(I2[149] - (lpoz)*I3[149])
           + (threeo2zn)*I4[194];
*(vp++) = U02*I0[262] + U42*I1[262]
           + (threeo2z)*(I2[150] - (lpoz)*I3[150]);
*(vp++) = U02*I0[263] + U42*I1[263]
           + (threeo2z)*(I2[151] - (lpoz)*I3[151])
           + (oneo2zn)*I4[195];
*(vp++) = U02*I0[264] + U42*I1[264]
           + (threeo2z)*(I2[152] - (lpoz)*I3[152])
           + (twoo2zn)*I4[196];
*(vp++) = U02*I0[265] + U42*I1[265]
           + (threeo2z)*(I2[153] - (lpoz)*I3[153])
           + (threeo2zn)*I4[197];
*(vp++) = U02*I0[266] + U42*I1[266]
           + (threeo2z)*(I2[154] - (lpoz)*I3[154])
           + (fouro2zn)*I4[198];
*(vp++) = U02*I0[267] + U42*I1[267]
           + (threeo2z)*(I2[155] - (lpoz)*I3[155]);
*(vp++) = U02*I0[268] + U42*I1[268]
           + (threeo2z)*(I2[156] - (lpoz)*I3[156])
           + (oneo2zn)*I4[199];
*(vp++) = U02*I0[269] + U42*I1[269]
           + (threeo2z)*(I2[157] - (lpoz)*I3[157])
           + (twoo2zn)*I4[200];
*(vp++) = U02*I0[270] + U42*I1[270]
           + (threeo2z)*(I2[158] - (lpoz)*I3[158])
           + (threeo2zn)*I4[201];
*(vp++) = U02*I0[271] + U42*I1[271]
           + (threeo2z)*(I2[159] - (lpoz)*I3[159])
           + (fouro2zn)*I4[202];
*(vp++) = U02*I0[272] + U42*I1[272]
           + (threeo2z)*(I2[160] - (lpoz)*I3[160])
           + (fiveo2zn)*I4[203];
*(vp++) = U02*I0[273] + U42*I1[273]
           + (threeo2z)*(I2[161] - (lpoz)*I3[161]);
*(vp++) = U02*I0[274] + U42*I1[274]
           + (threeo2z)*(I2[162] - (lpoz)*I3[162])
           + (oneo2zn)*I4[204];
*(vp++) = U02*I0[275] + U42*I1[275]
           + (threeo2z)*(I2[163] - (lpoz)*I3[163])
           + (twoo2zn)*I4[205];
*(vp++) = U02*I0[276] + U42*I1[276]
           + (threeo2z)*(I2[164] - (lpoz)*I3[164])
           + (threeo2zn)*I4[206];
*(vp++) = U02*I0[277] + U42*I1[277]
           + (threeo2z)*(I2[165] - (lpoz)*I3[165])
           + (fouro2zn)*I4[207];
*(vp++) = U02*I0[278] + U42*I1[278]
           + (threeo2z)*(I2[166] - (lpoz)*I3[166])
           + (fiveo2zn)*I4[208];
*(vp++) = U02*I0[279] + U42*I1[279]
           + (threeo2z)*(I2[167] - (lpoz)*I3[167])
           + (sixo2zn)*I4[209];
return vp;
}
/* Total number of FLOPs = 3010 */
