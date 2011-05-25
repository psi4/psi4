  /* These machine-generated functions compute a quartet of (hs|hs) integrals */

#include "libint.h"

REALTYPE *_build_h0h0_0(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
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
  REALTYPE threeo2z;
  REALTYPE fouro2z;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  threeo2zn = 3.0*Data->oo2zn;
  fouro2zn = 4.0*Data->oo2zn;
  fiveo2zn = 5.0*Data->oo2zn;
  oneo2z = 1.0*Data->oo2z;
  twoo2z = 2.0*Data->oo2z;
  threeo2z = 3.0*Data->oo2z;
  fouro2z = 4.0*Data->oo2z;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[0] + U40*I1[0]
           + (fouro2z)*(I2[0] - (lpoz)*I3[0])
           + (fiveo2zn)*I4[0];
*(vp++) = U00*I0[1] + U40*I1[1]
           + (fouro2z)*(I2[1] - (lpoz)*I3[1])
           + (fouro2zn)*I4[1];
*(vp++) = U00*I0[2] + U40*I1[2]
           + (fouro2z)*(I2[2] - (lpoz)*I3[2])
           + (fouro2zn)*I4[2];
*(vp++) = U00*I0[3] + U40*I1[3]
           + (fouro2z)*(I2[3] - (lpoz)*I3[3])
           + (threeo2zn)*I4[3];
*(vp++) = U00*I0[4] + U40*I1[4]
           + (fouro2z)*(I2[4] - (lpoz)*I3[4])
           + (threeo2zn)*I4[4];
*(vp++) = U00*I0[5] + U40*I1[5]
           + (fouro2z)*(I2[5] - (lpoz)*I3[5])
           + (threeo2zn)*I4[5];
*(vp++) = U00*I0[6] + U40*I1[6]
           + (fouro2z)*(I2[6] - (lpoz)*I3[6])
           + (twoo2zn)*I4[6];
*(vp++) = U00*I0[7] + U40*I1[7]
           + (fouro2z)*(I2[7] - (lpoz)*I3[7])
           + (twoo2zn)*I4[7];
*(vp++) = U00*I0[8] + U40*I1[8]
           + (fouro2z)*(I2[8] - (lpoz)*I3[8])
           + (twoo2zn)*I4[8];
*(vp++) = U00*I0[9] + U40*I1[9]
           + (fouro2z)*(I2[9] - (lpoz)*I3[9])
           + (twoo2zn)*I4[9];
*(vp++) = U00*I0[10] + U40*I1[10]
           + (fouro2z)*(I2[10] - (lpoz)*I3[10])
           + (oneo2zn)*I4[10];
*(vp++) = U00*I0[11] + U40*I1[11]
           + (fouro2z)*(I2[11] - (lpoz)*I3[11])
           + (oneo2zn)*I4[11];
*(vp++) = U00*I0[12] + U40*I1[12]
           + (fouro2z)*(I2[12] - (lpoz)*I3[12])
           + (oneo2zn)*I4[12];
*(vp++) = U00*I0[13] + U40*I1[13]
           + (fouro2z)*(I2[13] - (lpoz)*I3[13])
           + (oneo2zn)*I4[13];
*(vp++) = U00*I0[14] + U40*I1[14]
           + (fouro2z)*(I2[14] - (lpoz)*I3[14])
           + (oneo2zn)*I4[14];
*(vp++) = U00*I0[15] + U40*I1[15]
           + (fouro2z)*(I2[15] - (lpoz)*I3[15]);
*(vp++) = U00*I0[16] + U40*I1[16]
           + (fouro2z)*(I2[16] - (lpoz)*I3[16]);
*(vp++) = U00*I0[17] + U40*I1[17]
           + (fouro2z)*(I2[17] - (lpoz)*I3[17]);
*(vp++) = U00*I0[18] + U40*I1[18]
           + (fouro2z)*(I2[18] - (lpoz)*I3[18]);
*(vp++) = U00*I0[19] + U40*I1[19]
           + (fouro2z)*(I2[19] - (lpoz)*I3[19]);
*(vp++) = U00*I0[20] + U40*I1[20]
           + (fouro2z)*(I2[20] - (lpoz)*I3[20]);
*(vp++) = U00*I0[21] + U40*I1[21]
           + (threeo2z)*(I2[21] - (lpoz)*I3[21])
           + (fiveo2zn)*I4[15];
*(vp++) = U00*I0[22] + U40*I1[22]
           + (threeo2z)*(I2[22] - (lpoz)*I3[22])
           + (fouro2zn)*I4[16];
*(vp++) = U00*I0[23] + U40*I1[23]
           + (threeo2z)*(I2[23] - (lpoz)*I3[23])
           + (fouro2zn)*I4[17];
*(vp++) = U00*I0[24] + U40*I1[24]
           + (threeo2z)*(I2[24] - (lpoz)*I3[24])
           + (threeo2zn)*I4[18];
*(vp++) = U00*I0[25] + U40*I1[25]
           + (threeo2z)*(I2[25] - (lpoz)*I3[25])
           + (threeo2zn)*I4[19];
*(vp++) = U00*I0[26] + U40*I1[26]
           + (threeo2z)*(I2[26] - (lpoz)*I3[26])
           + (threeo2zn)*I4[20];
*(vp++) = U00*I0[27] + U40*I1[27]
           + (threeo2z)*(I2[27] - (lpoz)*I3[27])
           + (twoo2zn)*I4[21];
*(vp++) = U00*I0[28] + U40*I1[28]
           + (threeo2z)*(I2[28] - (lpoz)*I3[28])
           + (twoo2zn)*I4[22];
*(vp++) = U00*I0[29] + U40*I1[29]
           + (threeo2z)*(I2[29] - (lpoz)*I3[29])
           + (twoo2zn)*I4[23];
*(vp++) = U00*I0[30] + U40*I1[30]
           + (threeo2z)*(I2[30] - (lpoz)*I3[30])
           + (twoo2zn)*I4[24];
*(vp++) = U00*I0[31] + U40*I1[31]
           + (threeo2z)*(I2[31] - (lpoz)*I3[31])
           + (oneo2zn)*I4[25];
*(vp++) = U00*I0[32] + U40*I1[32]
           + (threeo2z)*(I2[32] - (lpoz)*I3[32])
           + (oneo2zn)*I4[26];
*(vp++) = U00*I0[33] + U40*I1[33]
           + (threeo2z)*(I2[33] - (lpoz)*I3[33])
           + (oneo2zn)*I4[27];
*(vp++) = U00*I0[34] + U40*I1[34]
           + (threeo2z)*(I2[34] - (lpoz)*I3[34])
           + (oneo2zn)*I4[28];
*(vp++) = U00*I0[35] + U40*I1[35]
           + (threeo2z)*(I2[35] - (lpoz)*I3[35])
           + (oneo2zn)*I4[29];
*(vp++) = U00*I0[36] + U40*I1[36]
           + (threeo2z)*(I2[36] - (lpoz)*I3[36]);
*(vp++) = U00*I0[37] + U40*I1[37]
           + (threeo2z)*(I2[37] - (lpoz)*I3[37]);
*(vp++) = U00*I0[38] + U40*I1[38]
           + (threeo2z)*(I2[38] - (lpoz)*I3[38]);
*(vp++) = U00*I0[39] + U40*I1[39]
           + (threeo2z)*(I2[39] - (lpoz)*I3[39]);
*(vp++) = U00*I0[40] + U40*I1[40]
           + (threeo2z)*(I2[40] - (lpoz)*I3[40]);
*(vp++) = U00*I0[41] + U40*I1[41]
           + (threeo2z)*(I2[41] - (lpoz)*I3[41]);
*(vp++) = U00*I0[42] + U40*I1[42]
           + (threeo2z)*(I2[42] - (lpoz)*I3[42])
           + (fiveo2zn)*I4[30];
*(vp++) = U00*I0[43] + U40*I1[43]
           + (threeo2z)*(I2[43] - (lpoz)*I3[43])
           + (fouro2zn)*I4[31];
*(vp++) = U00*I0[44] + U40*I1[44]
           + (threeo2z)*(I2[44] - (lpoz)*I3[44])
           + (fouro2zn)*I4[32];
*(vp++) = U00*I0[45] + U40*I1[45]
           + (threeo2z)*(I2[45] - (lpoz)*I3[45])
           + (threeo2zn)*I4[33];
*(vp++) = U00*I0[46] + U40*I1[46]
           + (threeo2z)*(I2[46] - (lpoz)*I3[46])
           + (threeo2zn)*I4[34];
*(vp++) = U00*I0[47] + U40*I1[47]
           + (threeo2z)*(I2[47] - (lpoz)*I3[47])
           + (threeo2zn)*I4[35];
*(vp++) = U00*I0[48] + U40*I1[48]
           + (threeo2z)*(I2[48] - (lpoz)*I3[48])
           + (twoo2zn)*I4[36];
*(vp++) = U00*I0[49] + U40*I1[49]
           + (threeo2z)*(I2[49] - (lpoz)*I3[49])
           + (twoo2zn)*I4[37];
*(vp++) = U00*I0[50] + U40*I1[50]
           + (threeo2z)*(I2[50] - (lpoz)*I3[50])
           + (twoo2zn)*I4[38];
*(vp++) = U00*I0[51] + U40*I1[51]
           + (threeo2z)*(I2[51] - (lpoz)*I3[51])
           + (twoo2zn)*I4[39];
*(vp++) = U00*I0[52] + U40*I1[52]
           + (threeo2z)*(I2[52] - (lpoz)*I3[52])
           + (oneo2zn)*I4[40];
*(vp++) = U00*I0[53] + U40*I1[53]
           + (threeo2z)*(I2[53] - (lpoz)*I3[53])
           + (oneo2zn)*I4[41];
*(vp++) = U00*I0[54] + U40*I1[54]
           + (threeo2z)*(I2[54] - (lpoz)*I3[54])
           + (oneo2zn)*I4[42];
*(vp++) = U00*I0[55] + U40*I1[55]
           + (threeo2z)*(I2[55] - (lpoz)*I3[55])
           + (oneo2zn)*I4[43];
*(vp++) = U00*I0[56] + U40*I1[56]
           + (threeo2z)*(I2[56] - (lpoz)*I3[56])
           + (oneo2zn)*I4[44];
*(vp++) = U00*I0[57] + U40*I1[57]
           + (threeo2z)*(I2[57] - (lpoz)*I3[57]);
*(vp++) = U00*I0[58] + U40*I1[58]
           + (threeo2z)*(I2[58] - (lpoz)*I3[58]);
*(vp++) = U00*I0[59] + U40*I1[59]
           + (threeo2z)*(I2[59] - (lpoz)*I3[59]);
*(vp++) = U00*I0[60] + U40*I1[60]
           + (threeo2z)*(I2[60] - (lpoz)*I3[60]);
*(vp++) = U00*I0[61] + U40*I1[61]
           + (threeo2z)*(I2[61] - (lpoz)*I3[61]);
*(vp++) = U00*I0[62] + U40*I1[62]
           + (threeo2z)*(I2[62] - (lpoz)*I3[62]);
*(vp++) = U00*I0[63] + U40*I1[63]
           + (twoo2z)*(I2[63] - (lpoz)*I3[63])
           + (fiveo2zn)*I4[45];
*(vp++) = U00*I0[64] + U40*I1[64]
           + (twoo2z)*(I2[64] - (lpoz)*I3[64])
           + (fouro2zn)*I4[46];
*(vp++) = U00*I0[65] + U40*I1[65]
           + (twoo2z)*(I2[65] - (lpoz)*I3[65])
           + (fouro2zn)*I4[47];
*(vp++) = U00*I0[66] + U40*I1[66]
           + (twoo2z)*(I2[66] - (lpoz)*I3[66])
           + (threeo2zn)*I4[48];
*(vp++) = U00*I0[67] + U40*I1[67]
           + (twoo2z)*(I2[67] - (lpoz)*I3[67])
           + (threeo2zn)*I4[49];
*(vp++) = U00*I0[68] + U40*I1[68]
           + (twoo2z)*(I2[68] - (lpoz)*I3[68])
           + (threeo2zn)*I4[50];
*(vp++) = U00*I0[69] + U40*I1[69]
           + (twoo2z)*(I2[69] - (lpoz)*I3[69])
           + (twoo2zn)*I4[51];
*(vp++) = U00*I0[70] + U40*I1[70]
           + (twoo2z)*(I2[70] - (lpoz)*I3[70])
           + (twoo2zn)*I4[52];
*(vp++) = U00*I0[71] + U40*I1[71]
           + (twoo2z)*(I2[71] - (lpoz)*I3[71])
           + (twoo2zn)*I4[53];
*(vp++) = U00*I0[72] + U40*I1[72]
           + (twoo2z)*(I2[72] - (lpoz)*I3[72])
           + (twoo2zn)*I4[54];
*(vp++) = U00*I0[73] + U40*I1[73]
           + (twoo2z)*(I2[73] - (lpoz)*I3[73])
           + (oneo2zn)*I4[55];
*(vp++) = U00*I0[74] + U40*I1[74]
           + (twoo2z)*(I2[74] - (lpoz)*I3[74])
           + (oneo2zn)*I4[56];
*(vp++) = U00*I0[75] + U40*I1[75]
           + (twoo2z)*(I2[75] - (lpoz)*I3[75])
           + (oneo2zn)*I4[57];
*(vp++) = U00*I0[76] + U40*I1[76]
           + (twoo2z)*(I2[76] - (lpoz)*I3[76])
           + (oneo2zn)*I4[58];
*(vp++) = U00*I0[77] + U40*I1[77]
           + (twoo2z)*(I2[77] - (lpoz)*I3[77])
           + (oneo2zn)*I4[59];
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
           + (twoo2z)*(I2[84] - (lpoz)*I3[84])
           + (fiveo2zn)*I4[60];
*(vp++) = U00*I0[85] + U40*I1[85]
           + (twoo2z)*(I2[85] - (lpoz)*I3[85])
           + (fouro2zn)*I4[61];
*(vp++) = U00*I0[86] + U40*I1[86]
           + (twoo2z)*(I2[86] - (lpoz)*I3[86])
           + (fouro2zn)*I4[62];
*(vp++) = U00*I0[87] + U40*I1[87]
           + (twoo2z)*(I2[87] - (lpoz)*I3[87])
           + (threeo2zn)*I4[63];
*(vp++) = U00*I0[88] + U40*I1[88]
           + (twoo2z)*(I2[88] - (lpoz)*I3[88])
           + (threeo2zn)*I4[64];
*(vp++) = U00*I0[89] + U40*I1[89]
           + (twoo2z)*(I2[89] - (lpoz)*I3[89])
           + (threeo2zn)*I4[65];
*(vp++) = U00*I0[90] + U40*I1[90]
           + (twoo2z)*(I2[90] - (lpoz)*I3[90])
           + (twoo2zn)*I4[66];
*(vp++) = U00*I0[91] + U40*I1[91]
           + (twoo2z)*(I2[91] - (lpoz)*I3[91])
           + (twoo2zn)*I4[67];
*(vp++) = U00*I0[92] + U40*I1[92]
           + (twoo2z)*(I2[92] - (lpoz)*I3[92])
           + (twoo2zn)*I4[68];
*(vp++) = U00*I0[93] + U40*I1[93]
           + (twoo2z)*(I2[93] - (lpoz)*I3[93])
           + (twoo2zn)*I4[69];
*(vp++) = U00*I0[94] + U40*I1[94]
           + (twoo2z)*(I2[94] - (lpoz)*I3[94])
           + (oneo2zn)*I4[70];
*(vp++) = U00*I0[95] + U40*I1[95]
           + (twoo2z)*(I2[95] - (lpoz)*I3[95])
           + (oneo2zn)*I4[71];
*(vp++) = U00*I0[96] + U40*I1[96]
           + (twoo2z)*(I2[96] - (lpoz)*I3[96])
           + (oneo2zn)*I4[72];
*(vp++) = U00*I0[97] + U40*I1[97]
           + (twoo2z)*(I2[97] - (lpoz)*I3[97])
           + (oneo2zn)*I4[73];
*(vp++) = U00*I0[98] + U40*I1[98]
           + (twoo2z)*(I2[98] - (lpoz)*I3[98])
           + (oneo2zn)*I4[74];
*(vp++) = U00*I0[99] + U40*I1[99]
           + (twoo2z)*(I2[99] - (lpoz)*I3[99]);
*(vp++) = U00*I0[100] + U40*I1[100]
           + (twoo2z)*(I2[100] - (lpoz)*I3[100]);
*(vp++) = U00*I0[101] + U40*I1[101]
           + (twoo2z)*(I2[101] - (lpoz)*I3[101]);
*(vp++) = U00*I0[102] + U40*I1[102]
           + (twoo2z)*(I2[102] - (lpoz)*I3[102]);
*(vp++) = U00*I0[103] + U40*I1[103]
           + (twoo2z)*(I2[103] - (lpoz)*I3[103]);
*(vp++) = U00*I0[104] + U40*I1[104]
           + (twoo2z)*(I2[104] - (lpoz)*I3[104]);
*(vp++) = U00*I0[105] + U40*I1[105]
           + (twoo2z)*(I2[105] - (lpoz)*I3[105])
           + (fiveo2zn)*I4[75];
*(vp++) = U00*I0[106] + U40*I1[106]
           + (twoo2z)*(I2[106] - (lpoz)*I3[106])
           + (fouro2zn)*I4[76];
*(vp++) = U00*I0[107] + U40*I1[107]
           + (twoo2z)*(I2[107] - (lpoz)*I3[107])
           + (fouro2zn)*I4[77];
*(vp++) = U00*I0[108] + U40*I1[108]
           + (twoo2z)*(I2[108] - (lpoz)*I3[108])
           + (threeo2zn)*I4[78];
*(vp++) = U00*I0[109] + U40*I1[109]
           + (twoo2z)*(I2[109] - (lpoz)*I3[109])
           + (threeo2zn)*I4[79];
*(vp++) = U00*I0[110] + U40*I1[110]
           + (twoo2z)*(I2[110] - (lpoz)*I3[110])
           + (threeo2zn)*I4[80];
*(vp++) = U00*I0[111] + U40*I1[111]
           + (twoo2z)*(I2[111] - (lpoz)*I3[111])
           + (twoo2zn)*I4[81];
*(vp++) = U00*I0[112] + U40*I1[112]
           + (twoo2z)*(I2[112] - (lpoz)*I3[112])
           + (twoo2zn)*I4[82];
*(vp++) = U00*I0[113] + U40*I1[113]
           + (twoo2z)*(I2[113] - (lpoz)*I3[113])
           + (twoo2zn)*I4[83];
*(vp++) = U00*I0[114] + U40*I1[114]
           + (twoo2z)*(I2[114] - (lpoz)*I3[114])
           + (twoo2zn)*I4[84];
*(vp++) = U00*I0[115] + U40*I1[115]
           + (twoo2z)*(I2[115] - (lpoz)*I3[115])
           + (oneo2zn)*I4[85];
*(vp++) = U00*I0[116] + U40*I1[116]
           + (twoo2z)*(I2[116] - (lpoz)*I3[116])
           + (oneo2zn)*I4[86];
*(vp++) = U00*I0[117] + U40*I1[117]
           + (twoo2z)*(I2[117] - (lpoz)*I3[117])
           + (oneo2zn)*I4[87];
*(vp++) = U00*I0[118] + U40*I1[118]
           + (twoo2z)*(I2[118] - (lpoz)*I3[118])
           + (oneo2zn)*I4[88];
*(vp++) = U00*I0[119] + U40*I1[119]
           + (twoo2z)*(I2[119] - (lpoz)*I3[119])
           + (oneo2zn)*I4[89];
*(vp++) = U00*I0[120] + U40*I1[120]
           + (twoo2z)*(I2[120] - (lpoz)*I3[120]);
*(vp++) = U00*I0[121] + U40*I1[121]
           + (twoo2z)*(I2[121] - (lpoz)*I3[121]);
*(vp++) = U00*I0[122] + U40*I1[122]
           + (twoo2z)*(I2[122] - (lpoz)*I3[122]);
*(vp++) = U00*I0[123] + U40*I1[123]
           + (twoo2z)*(I2[123] - (lpoz)*I3[123]);
*(vp++) = U00*I0[124] + U40*I1[124]
           + (twoo2z)*(I2[124] - (lpoz)*I3[124]);
*(vp++) = U00*I0[125] + U40*I1[125]
           + (twoo2z)*(I2[125] - (lpoz)*I3[125]);
*(vp++) = U00*I0[126] + U40*I1[126]
           + (oneo2z)*(I2[126] - (lpoz)*I3[126])
           + (fiveo2zn)*I4[90];
*(vp++) = U00*I0[127] + U40*I1[127]
           + (oneo2z)*(I2[127] - (lpoz)*I3[127])
           + (fouro2zn)*I4[91];
*(vp++) = U00*I0[128] + U40*I1[128]
           + (oneo2z)*(I2[128] - (lpoz)*I3[128])
           + (fouro2zn)*I4[92];
*(vp++) = U00*I0[129] + U40*I1[129]
           + (oneo2z)*(I2[129] - (lpoz)*I3[129])
           + (threeo2zn)*I4[93];
*(vp++) = U00*I0[130] + U40*I1[130]
           + (oneo2z)*(I2[130] - (lpoz)*I3[130])
           + (threeo2zn)*I4[94];
*(vp++) = U00*I0[131] + U40*I1[131]
           + (oneo2z)*(I2[131] - (lpoz)*I3[131])
           + (threeo2zn)*I4[95];
*(vp++) = U00*I0[132] + U40*I1[132]
           + (oneo2z)*(I2[132] - (lpoz)*I3[132])
           + (twoo2zn)*I4[96];
*(vp++) = U00*I0[133] + U40*I1[133]
           + (oneo2z)*(I2[133] - (lpoz)*I3[133])
           + (twoo2zn)*I4[97];
*(vp++) = U00*I0[134] + U40*I1[134]
           + (oneo2z)*(I2[134] - (lpoz)*I3[134])
           + (twoo2zn)*I4[98];
*(vp++) = U00*I0[135] + U40*I1[135]
           + (oneo2z)*(I2[135] - (lpoz)*I3[135])
           + (twoo2zn)*I4[99];
*(vp++) = U00*I0[136] + U40*I1[136]
           + (oneo2z)*(I2[136] - (lpoz)*I3[136])
           + (oneo2zn)*I4[100];
*(vp++) = U00*I0[137] + U40*I1[137]
           + (oneo2z)*(I2[137] - (lpoz)*I3[137])
           + (oneo2zn)*I4[101];
*(vp++) = U00*I0[138] + U40*I1[138]
           + (oneo2z)*(I2[138] - (lpoz)*I3[138])
           + (oneo2zn)*I4[102];
*(vp++) = U00*I0[139] + U40*I1[139]
           + (oneo2z)*(I2[139] - (lpoz)*I3[139])
           + (oneo2zn)*I4[103];
*(vp++) = U00*I0[140] + U40*I1[140]
           + (oneo2z)*(I2[140] - (lpoz)*I3[140])
           + (oneo2zn)*I4[104];
*(vp++) = U00*I0[141] + U40*I1[141]
           + (oneo2z)*(I2[141] - (lpoz)*I3[141]);
*(vp++) = U00*I0[142] + U40*I1[142]
           + (oneo2z)*(I2[142] - (lpoz)*I3[142]);
*(vp++) = U00*I0[143] + U40*I1[143]
           + (oneo2z)*(I2[143] - (lpoz)*I3[143]);
*(vp++) = U00*I0[144] + U40*I1[144]
           + (oneo2z)*(I2[144] - (lpoz)*I3[144]);
*(vp++) = U00*I0[145] + U40*I1[145]
           + (oneo2z)*(I2[145] - (lpoz)*I3[145]);
*(vp++) = U00*I0[146] + U40*I1[146]
           + (oneo2z)*(I2[146] - (lpoz)*I3[146]);
*(vp++) = U00*I0[147] + U40*I1[147]
           + (oneo2z)*(I2[147] - (lpoz)*I3[147])
           + (fiveo2zn)*I4[105];
*(vp++) = U00*I0[148] + U40*I1[148]
           + (oneo2z)*(I2[148] - (lpoz)*I3[148])
           + (fouro2zn)*I4[106];
*(vp++) = U00*I0[149] + U40*I1[149]
           + (oneo2z)*(I2[149] - (lpoz)*I3[149])
           + (fouro2zn)*I4[107];
*(vp++) = U00*I0[150] + U40*I1[150]
           + (oneo2z)*(I2[150] - (lpoz)*I3[150])
           + (threeo2zn)*I4[108];
*(vp++) = U00*I0[151] + U40*I1[151]
           + (oneo2z)*(I2[151] - (lpoz)*I3[151])
           + (threeo2zn)*I4[109];
*(vp++) = U00*I0[152] + U40*I1[152]
           + (oneo2z)*(I2[152] - (lpoz)*I3[152])
           + (threeo2zn)*I4[110];
*(vp++) = U00*I0[153] + U40*I1[153]
           + (oneo2z)*(I2[153] - (lpoz)*I3[153])
           + (twoo2zn)*I4[111];
*(vp++) = U00*I0[154] + U40*I1[154]
           + (oneo2z)*(I2[154] - (lpoz)*I3[154])
           + (twoo2zn)*I4[112];
*(vp++) = U00*I0[155] + U40*I1[155]
           + (oneo2z)*(I2[155] - (lpoz)*I3[155])
           + (twoo2zn)*I4[113];
*(vp++) = U00*I0[156] + U40*I1[156]
           + (oneo2z)*(I2[156] - (lpoz)*I3[156])
           + (twoo2zn)*I4[114];
*(vp++) = U00*I0[157] + U40*I1[157]
           + (oneo2z)*(I2[157] - (lpoz)*I3[157])
           + (oneo2zn)*I4[115];
*(vp++) = U00*I0[158] + U40*I1[158]
           + (oneo2z)*(I2[158] - (lpoz)*I3[158])
           + (oneo2zn)*I4[116];
*(vp++) = U00*I0[159] + U40*I1[159]
           + (oneo2z)*(I2[159] - (lpoz)*I3[159])
           + (oneo2zn)*I4[117];
*(vp++) = U00*I0[160] + U40*I1[160]
           + (oneo2z)*(I2[160] - (lpoz)*I3[160])
           + (oneo2zn)*I4[118];
*(vp++) = U00*I0[161] + U40*I1[161]
           + (oneo2z)*(I2[161] - (lpoz)*I3[161])
           + (oneo2zn)*I4[119];
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
           + (oneo2z)*(I2[168] - (lpoz)*I3[168])
           + (fiveo2zn)*I4[120];
*(vp++) = U00*I0[169] + U40*I1[169]
           + (oneo2z)*(I2[169] - (lpoz)*I3[169])
           + (fouro2zn)*I4[121];
*(vp++) = U00*I0[170] + U40*I1[170]
           + (oneo2z)*(I2[170] - (lpoz)*I3[170])
           + (fouro2zn)*I4[122];
*(vp++) = U00*I0[171] + U40*I1[171]
           + (oneo2z)*(I2[171] - (lpoz)*I3[171])
           + (threeo2zn)*I4[123];
*(vp++) = U00*I0[172] + U40*I1[172]
           + (oneo2z)*(I2[172] - (lpoz)*I3[172])
           + (threeo2zn)*I4[124];
*(vp++) = U00*I0[173] + U40*I1[173]
           + (oneo2z)*(I2[173] - (lpoz)*I3[173])
           + (threeo2zn)*I4[125];
*(vp++) = U00*I0[174] + U40*I1[174]
           + (oneo2z)*(I2[174] - (lpoz)*I3[174])
           + (twoo2zn)*I4[126];
*(vp++) = U00*I0[175] + U40*I1[175]
           + (oneo2z)*(I2[175] - (lpoz)*I3[175])
           + (twoo2zn)*I4[127];
*(vp++) = U00*I0[176] + U40*I1[176]
           + (oneo2z)*(I2[176] - (lpoz)*I3[176])
           + (twoo2zn)*I4[128];
*(vp++) = U00*I0[177] + U40*I1[177]
           + (oneo2z)*(I2[177] - (lpoz)*I3[177])
           + (twoo2zn)*I4[129];
*(vp++) = U00*I0[178] + U40*I1[178]
           + (oneo2z)*(I2[178] - (lpoz)*I3[178])
           + (oneo2zn)*I4[130];
*(vp++) = U00*I0[179] + U40*I1[179]
           + (oneo2z)*(I2[179] - (lpoz)*I3[179])
           + (oneo2zn)*I4[131];
*(vp++) = U00*I0[180] + U40*I1[180]
           + (oneo2z)*(I2[180] - (lpoz)*I3[180])
           + (oneo2zn)*I4[132];
*(vp++) = U00*I0[181] + U40*I1[181]
           + (oneo2z)*(I2[181] - (lpoz)*I3[181])
           + (oneo2zn)*I4[133];
*(vp++) = U00*I0[182] + U40*I1[182]
           + (oneo2z)*(I2[182] - (lpoz)*I3[182])
           + (oneo2zn)*I4[134];
*(vp++) = U00*I0[183] + U40*I1[183]
           + (oneo2z)*(I2[183] - (lpoz)*I3[183]);
*(vp++) = U00*I0[184] + U40*I1[184]
           + (oneo2z)*(I2[184] - (lpoz)*I3[184]);
*(vp++) = U00*I0[185] + U40*I1[185]
           + (oneo2z)*(I2[185] - (lpoz)*I3[185]);
*(vp++) = U00*I0[186] + U40*I1[186]
           + (oneo2z)*(I2[186] - (lpoz)*I3[186]);
*(vp++) = U00*I0[187] + U40*I1[187]
           + (oneo2z)*(I2[187] - (lpoz)*I3[187]);
*(vp++) = U00*I0[188] + U40*I1[188]
           + (oneo2z)*(I2[188] - (lpoz)*I3[188]);
*(vp++) = U00*I0[189] + U40*I1[189]
           + (oneo2z)*(I2[189] - (lpoz)*I3[189])
           + (fiveo2zn)*I4[135];
*(vp++) = U00*I0[190] + U40*I1[190]
           + (oneo2z)*(I2[190] - (lpoz)*I3[190])
           + (fouro2zn)*I4[136];
*(vp++) = U00*I0[191] + U40*I1[191]
           + (oneo2z)*(I2[191] - (lpoz)*I3[191])
           + (fouro2zn)*I4[137];
*(vp++) = U00*I0[192] + U40*I1[192]
           + (oneo2z)*(I2[192] - (lpoz)*I3[192])
           + (threeo2zn)*I4[138];
*(vp++) = U00*I0[193] + U40*I1[193]
           + (oneo2z)*(I2[193] - (lpoz)*I3[193])
           + (threeo2zn)*I4[139];
*(vp++) = U00*I0[194] + U40*I1[194]
           + (oneo2z)*(I2[194] - (lpoz)*I3[194])
           + (threeo2zn)*I4[140];
*(vp++) = U00*I0[195] + U40*I1[195]
           + (oneo2z)*(I2[195] - (lpoz)*I3[195])
           + (twoo2zn)*I4[141];
*(vp++) = U00*I0[196] + U40*I1[196]
           + (oneo2z)*(I2[196] - (lpoz)*I3[196])
           + (twoo2zn)*I4[142];
*(vp++) = U00*I0[197] + U40*I1[197]
           + (oneo2z)*(I2[197] - (lpoz)*I3[197])
           + (twoo2zn)*I4[143];
*(vp++) = U00*I0[198] + U40*I1[198]
           + (oneo2z)*(I2[198] - (lpoz)*I3[198])
           + (twoo2zn)*I4[144];
*(vp++) = U00*I0[199] + U40*I1[199]
           + (oneo2z)*(I2[199] - (lpoz)*I3[199])
           + (oneo2zn)*I4[145];
*(vp++) = U00*I0[200] + U40*I1[200]
           + (oneo2z)*(I2[200] - (lpoz)*I3[200])
           + (oneo2zn)*I4[146];
*(vp++) = U00*I0[201] + U40*I1[201]
           + (oneo2z)*(I2[201] - (lpoz)*I3[201])
           + (oneo2zn)*I4[147];
*(vp++) = U00*I0[202] + U40*I1[202]
           + (oneo2z)*(I2[202] - (lpoz)*I3[202])
           + (oneo2zn)*I4[148];
*(vp++) = U00*I0[203] + U40*I1[203]
           + (oneo2z)*(I2[203] - (lpoz)*I3[203])
           + (oneo2zn)*I4[149];
*(vp++) = U00*I0[204] + U40*I1[204]
           + (oneo2z)*(I2[204] - (lpoz)*I3[204]);
*(vp++) = U00*I0[205] + U40*I1[205]
           + (oneo2z)*(I2[205] - (lpoz)*I3[205]);
*(vp++) = U00*I0[206] + U40*I1[206]
           + (oneo2z)*(I2[206] - (lpoz)*I3[206]);
*(vp++) = U00*I0[207] + U40*I1[207]
           + (oneo2z)*(I2[207] - (lpoz)*I3[207]);
*(vp++) = U00*I0[208] + U40*I1[208]
           + (oneo2z)*(I2[208] - (lpoz)*I3[208]);
*(vp++) = U00*I0[209] + U40*I1[209]
           + (oneo2z)*(I2[209] - (lpoz)*I3[209]);
*(vp++) = U00*I0[210] + U40*I1[210]
           + (fiveo2zn)*I4[150];
*(vp++) = U00*I0[211] + U40*I1[211]
           + (fouro2zn)*I4[151];
*(vp++) = U00*I0[212] + U40*I1[212]
           + (fouro2zn)*I4[152];
*(vp++) = U00*I0[213] + U40*I1[213]
           + (threeo2zn)*I4[153];
*(vp++) = U00*I0[214] + U40*I1[214]
           + (threeo2zn)*I4[154];
*(vp++) = U00*I0[215] + U40*I1[215]
           + (threeo2zn)*I4[155];
*(vp++) = U00*I0[216] + U40*I1[216]
           + (twoo2zn)*I4[156];
*(vp++) = U00*I0[217] + U40*I1[217]
           + (twoo2zn)*I4[157];
*(vp++) = U00*I0[218] + U40*I1[218]
           + (twoo2zn)*I4[158];
*(vp++) = U00*I0[219] + U40*I1[219]
           + (twoo2zn)*I4[159];
*(vp++) = U00*I0[220] + U40*I1[220]
           + (oneo2zn)*I4[160];
return vp;
}

REALTYPE *_build_h0h0_1(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)
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
  REALTYPE threeo2z;
  REALTYPE fouro2z;
  oneo2zn = 1.0*Data->oo2zn;
  twoo2zn = 2.0*Data->oo2zn;
  threeo2zn = 3.0*Data->oo2zn;
  fouro2zn = 4.0*Data->oo2zn;
  fiveo2zn = 5.0*Data->oo2zn;
  oneo2z = 1.0*Data->oo2z;
  twoo2z = 2.0*Data->oo2z;
  threeo2z = 3.0*Data->oo2z;
  fouro2z = 4.0*Data->oo2z;
  U00 = Data->U[0][0];
  U01 = Data->U[0][1];
  U02 = Data->U[0][2];
  U40 = Data->U[4][0];
  U41 = Data->U[4][1];
  U42 = Data->U[4][2];


*(vp++) = U00*I0[221] + U40*I1[221]
           + (oneo2zn)*I4[161];
*(vp++) = U00*I0[222] + U40*I1[222]
           + (oneo2zn)*I4[162];
*(vp++) = U00*I0[223] + U40*I1[223]
           + (oneo2zn)*I4[163];
*(vp++) = U00*I0[224] + U40*I1[224]
           + (oneo2zn)*I4[164];
*(vp++) = U00*I0[225] + U40*I1[225];
*(vp++) = U00*I0[226] + U40*I1[226];
*(vp++) = U00*I0[227] + U40*I1[227];
*(vp++) = U00*I0[228] + U40*I1[228];
*(vp++) = U00*I0[229] + U40*I1[229];
*(vp++) = U00*I0[230] + U40*I1[230];
*(vp++) = U00*I0[231] + U40*I1[231]
           + (fiveo2zn)*I4[165];
*(vp++) = U00*I0[232] + U40*I1[232]
           + (fouro2zn)*I4[166];
*(vp++) = U00*I0[233] + U40*I1[233]
           + (fouro2zn)*I4[167];
*(vp++) = U00*I0[234] + U40*I1[234]
           + (threeo2zn)*I4[168];
*(vp++) = U00*I0[235] + U40*I1[235]
           + (threeo2zn)*I4[169];
*(vp++) = U00*I0[236] + U40*I1[236]
           + (threeo2zn)*I4[170];
*(vp++) = U00*I0[237] + U40*I1[237]
           + (twoo2zn)*I4[171];
*(vp++) = U00*I0[238] + U40*I1[238]
           + (twoo2zn)*I4[172];
*(vp++) = U00*I0[239] + U40*I1[239]
           + (twoo2zn)*I4[173];
*(vp++) = U00*I0[240] + U40*I1[240]
           + (twoo2zn)*I4[174];
*(vp++) = U00*I0[241] + U40*I1[241]
           + (oneo2zn)*I4[175];
*(vp++) = U00*I0[242] + U40*I1[242]
           + (oneo2zn)*I4[176];
*(vp++) = U00*I0[243] + U40*I1[243]
           + (oneo2zn)*I4[177];
*(vp++) = U00*I0[244] + U40*I1[244]
           + (oneo2zn)*I4[178];
*(vp++) = U00*I0[245] + U40*I1[245]
           + (oneo2zn)*I4[179];
*(vp++) = U00*I0[246] + U40*I1[246];
*(vp++) = U00*I0[247] + U40*I1[247];
*(vp++) = U00*I0[248] + U40*I1[248];
*(vp++) = U00*I0[249] + U40*I1[249];
*(vp++) = U00*I0[250] + U40*I1[250];
*(vp++) = U00*I0[251] + U40*I1[251];
*(vp++) = U00*I0[252] + U40*I1[252]
           + (fiveo2zn)*I4[180];
*(vp++) = U00*I0[253] + U40*I1[253]
           + (fouro2zn)*I4[181];
*(vp++) = U00*I0[254] + U40*I1[254]
           + (fouro2zn)*I4[182];
*(vp++) = U00*I0[255] + U40*I1[255]
           + (threeo2zn)*I4[183];
*(vp++) = U00*I0[256] + U40*I1[256]
           + (threeo2zn)*I4[184];
*(vp++) = U00*I0[257] + U40*I1[257]
           + (threeo2zn)*I4[185];
*(vp++) = U00*I0[258] + U40*I1[258]
           + (twoo2zn)*I4[186];
*(vp++) = U00*I0[259] + U40*I1[259]
           + (twoo2zn)*I4[187];
*(vp++) = U00*I0[260] + U40*I1[260]
           + (twoo2zn)*I4[188];
*(vp++) = U00*I0[261] + U40*I1[261]
           + (twoo2zn)*I4[189];
*(vp++) = U00*I0[262] + U40*I1[262]
           + (oneo2zn)*I4[190];
*(vp++) = U00*I0[263] + U40*I1[263]
           + (oneo2zn)*I4[191];
*(vp++) = U00*I0[264] + U40*I1[264]
           + (oneo2zn)*I4[192];
*(vp++) = U00*I0[265] + U40*I1[265]
           + (oneo2zn)*I4[193];
*(vp++) = U00*I0[266] + U40*I1[266]
           + (oneo2zn)*I4[194];
*(vp++) = U00*I0[267] + U40*I1[267];
*(vp++) = U00*I0[268] + U40*I1[268];
*(vp++) = U00*I0[269] + U40*I1[269];
*(vp++) = U00*I0[270] + U40*I1[270];
*(vp++) = U00*I0[271] + U40*I1[271];
*(vp++) = U00*I0[272] + U40*I1[272];
*(vp++) = U00*I0[273] + U40*I1[273]
           + (fiveo2zn)*I4[195];
*(vp++) = U00*I0[274] + U40*I1[274]
           + (fouro2zn)*I4[196];
*(vp++) = U00*I0[275] + U40*I1[275]
           + (fouro2zn)*I4[197];
*(vp++) = U00*I0[276] + U40*I1[276]
           + (threeo2zn)*I4[198];
*(vp++) = U00*I0[277] + U40*I1[277]
           + (threeo2zn)*I4[199];
*(vp++) = U00*I0[278] + U40*I1[278]
           + (threeo2zn)*I4[200];
*(vp++) = U00*I0[279] + U40*I1[279]
           + (twoo2zn)*I4[201];
*(vp++) = U00*I0[280] + U40*I1[280]
           + (twoo2zn)*I4[202];
*(vp++) = U00*I0[281] + U40*I1[281]
           + (twoo2zn)*I4[203];
*(vp++) = U00*I0[282] + U40*I1[282]
           + (twoo2zn)*I4[204];
*(vp++) = U00*I0[283] + U40*I1[283]
           + (oneo2zn)*I4[205];
*(vp++) = U00*I0[284] + U40*I1[284]
           + (oneo2zn)*I4[206];
*(vp++) = U00*I0[285] + U40*I1[285]
           + (oneo2zn)*I4[207];
*(vp++) = U00*I0[286] + U40*I1[286]
           + (oneo2zn)*I4[208];
*(vp++) = U00*I0[287] + U40*I1[287]
           + (oneo2zn)*I4[209];
*(vp++) = U00*I0[288] + U40*I1[288];
*(vp++) = U00*I0[289] + U40*I1[289];
*(vp++) = U00*I0[290] + U40*I1[290];
*(vp++) = U00*I0[291] + U40*I1[291];
*(vp++) = U00*I0[292] + U40*I1[292];
*(vp++) = U00*I0[293] + U40*I1[293];
*(vp++) = U00*I0[294] + U40*I1[294]
           + (fiveo2zn)*I4[210];
*(vp++) = U00*I0[295] + U40*I1[295]
           + (fouro2zn)*I4[211];
*(vp++) = U00*I0[296] + U40*I1[296]
           + (fouro2zn)*I4[212];
*(vp++) = U00*I0[297] + U40*I1[297]
           + (threeo2zn)*I4[213];
*(vp++) = U00*I0[298] + U40*I1[298]
           + (threeo2zn)*I4[214];
*(vp++) = U00*I0[299] + U40*I1[299]
           + (threeo2zn)*I4[215];
*(vp++) = U00*I0[300] + U40*I1[300]
           + (twoo2zn)*I4[216];
*(vp++) = U00*I0[301] + U40*I1[301]
           + (twoo2zn)*I4[217];
*(vp++) = U00*I0[302] + U40*I1[302]
           + (twoo2zn)*I4[218];
*(vp++) = U00*I0[303] + U40*I1[303]
           + (twoo2zn)*I4[219];
*(vp++) = U00*I0[304] + U40*I1[304]
           + (oneo2zn)*I4[220];
*(vp++) = U00*I0[305] + U40*I1[305]
           + (oneo2zn)*I4[221];
*(vp++) = U00*I0[306] + U40*I1[306]
           + (oneo2zn)*I4[222];
*(vp++) = U00*I0[307] + U40*I1[307]
           + (oneo2zn)*I4[223];
*(vp++) = U00*I0[308] + U40*I1[308]
           + (oneo2zn)*I4[224];
*(vp++) = U00*I0[309] + U40*I1[309];
*(vp++) = U00*I0[310] + U40*I1[310];
*(vp++) = U00*I0[311] + U40*I1[311];
*(vp++) = U00*I0[312] + U40*I1[312];
*(vp++) = U00*I0[313] + U40*I1[313];
*(vp++) = U00*I0[314] + U40*I1[314];
*(vp++) = U01*I0[210] + U41*I1[210]
           + (fouro2z)*(I2[126] - (lpoz)*I3[126]);
*(vp++) = U01*I0[211] + U41*I1[211]
           + (fouro2z)*(I2[127] - (lpoz)*I3[127])
           + (oneo2zn)*I4[150];
*(vp++) = U01*I0[212] + U41*I1[212]
           + (fouro2z)*(I2[128] - (lpoz)*I3[128]);
*(vp++) = U01*I0[213] + U41*I1[213]
           + (fouro2z)*(I2[129] - (lpoz)*I3[129])
           + (twoo2zn)*I4[151];
*(vp++) = U01*I0[214] + U41*I1[214]
           + (fouro2z)*(I2[130] - (lpoz)*I3[130])
           + (oneo2zn)*I4[152];
*(vp++) = U01*I0[215] + U41*I1[215]
           + (fouro2z)*(I2[131] - (lpoz)*I3[131]);
*(vp++) = U01*I0[216] + U41*I1[216]
           + (fouro2z)*(I2[132] - (lpoz)*I3[132])
           + (threeo2zn)*I4[153];
*(vp++) = U01*I0[217] + U41*I1[217]
           + (fouro2z)*(I2[133] - (lpoz)*I3[133])
           + (twoo2zn)*I4[154];
*(vp++) = U01*I0[218] + U41*I1[218]
           + (fouro2z)*(I2[134] - (lpoz)*I3[134])
           + (oneo2zn)*I4[155];
*(vp++) = U01*I0[219] + U41*I1[219]
           + (fouro2z)*(I2[135] - (lpoz)*I3[135]);
*(vp++) = U01*I0[220] + U41*I1[220]
           + (fouro2z)*(I2[136] - (lpoz)*I3[136])
           + (fouro2zn)*I4[156];
*(vp++) = U01*I0[221] + U41*I1[221]
           + (fouro2z)*(I2[137] - (lpoz)*I3[137])
           + (threeo2zn)*I4[157];
*(vp++) = U01*I0[222] + U41*I1[222]
           + (fouro2z)*(I2[138] - (lpoz)*I3[138])
           + (twoo2zn)*I4[158];
*(vp++) = U01*I0[223] + U41*I1[223]
           + (fouro2z)*(I2[139] - (lpoz)*I3[139])
           + (oneo2zn)*I4[159];
*(vp++) = U01*I0[224] + U41*I1[224]
           + (fouro2z)*(I2[140] - (lpoz)*I3[140]);
*(vp++) = U01*I0[225] + U41*I1[225]
           + (fouro2z)*(I2[141] - (lpoz)*I3[141])
           + (fiveo2zn)*I4[160];
*(vp++) = U01*I0[226] + U41*I1[226]
           + (fouro2z)*(I2[142] - (lpoz)*I3[142])
           + (fouro2zn)*I4[161];
*(vp++) = U01*I0[227] + U41*I1[227]
           + (fouro2z)*(I2[143] - (lpoz)*I3[143])
           + (threeo2zn)*I4[162];
*(vp++) = U01*I0[228] + U41*I1[228]
           + (fouro2z)*(I2[144] - (lpoz)*I3[144])
           + (twoo2zn)*I4[163];
*(vp++) = U01*I0[229] + U41*I1[229]
           + (fouro2z)*(I2[145] - (lpoz)*I3[145])
           + (oneo2zn)*I4[164];
*(vp++) = U01*I0[230] + U41*I1[230]
           + (fouro2z)*(I2[146] - (lpoz)*I3[146]);
*(vp++) = U01*I0[231] + U41*I1[231]
           + (threeo2z)*(I2[147] - (lpoz)*I3[147]);
*(vp++) = U01*I0[232] + U41*I1[232]
           + (threeo2z)*(I2[148] - (lpoz)*I3[148])
           + (oneo2zn)*I4[165];
*(vp++) = U01*I0[233] + U41*I1[233]
           + (threeo2z)*(I2[149] - (lpoz)*I3[149]);
*(vp++) = U01*I0[234] + U41*I1[234]
           + (threeo2z)*(I2[150] - (lpoz)*I3[150])
           + (twoo2zn)*I4[166];
*(vp++) = U01*I0[235] + U41*I1[235]
           + (threeo2z)*(I2[151] - (lpoz)*I3[151])
           + (oneo2zn)*I4[167];
*(vp++) = U01*I0[236] + U41*I1[236]
           + (threeo2z)*(I2[152] - (lpoz)*I3[152]);
*(vp++) = U01*I0[237] + U41*I1[237]
           + (threeo2z)*(I2[153] - (lpoz)*I3[153])
           + (threeo2zn)*I4[168];
*(vp++) = U01*I0[238] + U41*I1[238]
           + (threeo2z)*(I2[154] - (lpoz)*I3[154])
           + (twoo2zn)*I4[169];
*(vp++) = U01*I0[239] + U41*I1[239]
           + (threeo2z)*(I2[155] - (lpoz)*I3[155])
           + (oneo2zn)*I4[170];
*(vp++) = U01*I0[240] + U41*I1[240]
           + (threeo2z)*(I2[156] - (lpoz)*I3[156]);
*(vp++) = U01*I0[241] + U41*I1[241]
           + (threeo2z)*(I2[157] - (lpoz)*I3[157])
           + (fouro2zn)*I4[171];
*(vp++) = U01*I0[242] + U41*I1[242]
           + (threeo2z)*(I2[158] - (lpoz)*I3[158])
           + (threeo2zn)*I4[172];
*(vp++) = U01*I0[243] + U41*I1[243]
           + (threeo2z)*(I2[159] - (lpoz)*I3[159])
           + (twoo2zn)*I4[173];
*(vp++) = U01*I0[244] + U41*I1[244]
           + (threeo2z)*(I2[160] - (lpoz)*I3[160])
           + (oneo2zn)*I4[174];
*(vp++) = U01*I0[245] + U41*I1[245]
           + (threeo2z)*(I2[161] - (lpoz)*I3[161]);
*(vp++) = U01*I0[246] + U41*I1[246]
           + (threeo2z)*(I2[162] - (lpoz)*I3[162])
           + (fiveo2zn)*I4[175];
*(vp++) = U01*I0[247] + U41*I1[247]
           + (threeo2z)*(I2[163] - (lpoz)*I3[163])
           + (fouro2zn)*I4[176];
*(vp++) = U01*I0[248] + U41*I1[248]
           + (threeo2z)*(I2[164] - (lpoz)*I3[164])
           + (threeo2zn)*I4[177];
*(vp++) = U01*I0[249] + U41*I1[249]
           + (threeo2z)*(I2[165] - (lpoz)*I3[165])
           + (twoo2zn)*I4[178];
*(vp++) = U01*I0[250] + U41*I1[250]
           + (threeo2z)*(I2[166] - (lpoz)*I3[166])
           + (oneo2zn)*I4[179];
*(vp++) = U01*I0[251] + U41*I1[251]
           + (threeo2z)*(I2[167] - (lpoz)*I3[167]);
*(vp++) = U01*I0[252] + U41*I1[252]
           + (twoo2z)*(I2[168] - (lpoz)*I3[168]);
*(vp++) = U01*I0[253] + U41*I1[253]
           + (twoo2z)*(I2[169] - (lpoz)*I3[169])
           + (oneo2zn)*I4[180];
*(vp++) = U01*I0[254] + U41*I1[254]
           + (twoo2z)*(I2[170] - (lpoz)*I3[170]);
*(vp++) = U01*I0[255] + U41*I1[255]
           + (twoo2z)*(I2[171] - (lpoz)*I3[171])
           + (twoo2zn)*I4[181];
*(vp++) = U01*I0[256] + U41*I1[256]
           + (twoo2z)*(I2[172] - (lpoz)*I3[172])
           + (oneo2zn)*I4[182];
*(vp++) = U01*I0[257] + U41*I1[257]
           + (twoo2z)*(I2[173] - (lpoz)*I3[173]);
*(vp++) = U01*I0[258] + U41*I1[258]
           + (twoo2z)*(I2[174] - (lpoz)*I3[174])
           + (threeo2zn)*I4[183];
*(vp++) = U01*I0[259] + U41*I1[259]
           + (twoo2z)*(I2[175] - (lpoz)*I3[175])
           + (twoo2zn)*I4[184];
*(vp++) = U01*I0[260] + U41*I1[260]
           + (twoo2z)*(I2[176] - (lpoz)*I3[176])
           + (oneo2zn)*I4[185];
*(vp++) = U01*I0[261] + U41*I1[261]
           + (twoo2z)*(I2[177] - (lpoz)*I3[177]);
*(vp++) = U01*I0[262] + U41*I1[262]
           + (twoo2z)*(I2[178] - (lpoz)*I3[178])
           + (fouro2zn)*I4[186];
*(vp++) = U01*I0[263] + U41*I1[263]
           + (twoo2z)*(I2[179] - (lpoz)*I3[179])
           + (threeo2zn)*I4[187];
*(vp++) = U01*I0[264] + U41*I1[264]
           + (twoo2z)*(I2[180] - (lpoz)*I3[180])
           + (twoo2zn)*I4[188];
*(vp++) = U01*I0[265] + U41*I1[265]
           + (twoo2z)*(I2[181] - (lpoz)*I3[181])
           + (oneo2zn)*I4[189];
*(vp++) = U01*I0[266] + U41*I1[266]
           + (twoo2z)*(I2[182] - (lpoz)*I3[182]);
*(vp++) = U01*I0[267] + U41*I1[267]
           + (twoo2z)*(I2[183] - (lpoz)*I3[183])
           + (fiveo2zn)*I4[190];
*(vp++) = U01*I0[268] + U41*I1[268]
           + (twoo2z)*(I2[184] - (lpoz)*I3[184])
           + (fouro2zn)*I4[191];
*(vp++) = U01*I0[269] + U41*I1[269]
           + (twoo2z)*(I2[185] - (lpoz)*I3[185])
           + (threeo2zn)*I4[192];
*(vp++) = U01*I0[270] + U41*I1[270]
           + (twoo2z)*(I2[186] - (lpoz)*I3[186])
           + (twoo2zn)*I4[193];
*(vp++) = U01*I0[271] + U41*I1[271]
           + (twoo2z)*(I2[187] - (lpoz)*I3[187])
           + (oneo2zn)*I4[194];
*(vp++) = U01*I0[272] + U41*I1[272]
           + (twoo2z)*(I2[188] - (lpoz)*I3[188]);
*(vp++) = U01*I0[273] + U41*I1[273]
           + (oneo2z)*(I2[189] - (lpoz)*I3[189]);
*(vp++) = U01*I0[274] + U41*I1[274]
           + (oneo2z)*(I2[190] - (lpoz)*I3[190])
           + (oneo2zn)*I4[195];
*(vp++) = U01*I0[275] + U41*I1[275]
           + (oneo2z)*(I2[191] - (lpoz)*I3[191]);
*(vp++) = U01*I0[276] + U41*I1[276]
           + (oneo2z)*(I2[192] - (lpoz)*I3[192])
           + (twoo2zn)*I4[196];
*(vp++) = U01*I0[277] + U41*I1[277]
           + (oneo2z)*(I2[193] - (lpoz)*I3[193])
           + (oneo2zn)*I4[197];
*(vp++) = U01*I0[278] + U41*I1[278]
           + (oneo2z)*(I2[194] - (lpoz)*I3[194]);
*(vp++) = U01*I0[279] + U41*I1[279]
           + (oneo2z)*(I2[195] - (lpoz)*I3[195])
           + (threeo2zn)*I4[198];
*(vp++) = U01*I0[280] + U41*I1[280]
           + (oneo2z)*(I2[196] - (lpoz)*I3[196])
           + (twoo2zn)*I4[199];
*(vp++) = U01*I0[281] + U41*I1[281]
           + (oneo2z)*(I2[197] - (lpoz)*I3[197])
           + (oneo2zn)*I4[200];
*(vp++) = U01*I0[282] + U41*I1[282]
           + (oneo2z)*(I2[198] - (lpoz)*I3[198]);
*(vp++) = U01*I0[283] + U41*I1[283]
           + (oneo2z)*(I2[199] - (lpoz)*I3[199])
           + (fouro2zn)*I4[201];
*(vp++) = U01*I0[284] + U41*I1[284]
           + (oneo2z)*(I2[200] - (lpoz)*I3[200])
           + (threeo2zn)*I4[202];
*(vp++) = U01*I0[285] + U41*I1[285]
           + (oneo2z)*(I2[201] - (lpoz)*I3[201])
           + (twoo2zn)*I4[203];
*(vp++) = U01*I0[286] + U41*I1[286]
           + (oneo2z)*(I2[202] - (lpoz)*I3[202])
           + (oneo2zn)*I4[204];
*(vp++) = U01*I0[287] + U41*I1[287]
           + (oneo2z)*(I2[203] - (lpoz)*I3[203]);
*(vp++) = U01*I0[288] + U41*I1[288]
           + (oneo2z)*(I2[204] - (lpoz)*I3[204])
           + (fiveo2zn)*I4[205];
*(vp++) = U01*I0[289] + U41*I1[289]
           + (oneo2z)*(I2[205] - (lpoz)*I3[205])
           + (fouro2zn)*I4[206];
*(vp++) = U01*I0[290] + U41*I1[290]
           + (oneo2z)*(I2[206] - (lpoz)*I3[206])
           + (threeo2zn)*I4[207];
*(vp++) = U01*I0[291] + U41*I1[291]
           + (oneo2z)*(I2[207] - (lpoz)*I3[207])
           + (twoo2zn)*I4[208];
*(vp++) = U01*I0[292] + U41*I1[292]
           + (oneo2z)*(I2[208] - (lpoz)*I3[208])
           + (oneo2zn)*I4[209];
*(vp++) = U01*I0[293] + U41*I1[293]
           + (oneo2z)*(I2[209] - (lpoz)*I3[209]);
*(vp++) = U01*I0[294] + U41*I1[294];
*(vp++) = U01*I0[295] + U41*I1[295]
           + (oneo2zn)*I4[210];
*(vp++) = U01*I0[296] + U41*I1[296];
*(vp++) = U01*I0[297] + U41*I1[297]
           + (twoo2zn)*I4[211];
*(vp++) = U01*I0[298] + U41*I1[298]
           + (oneo2zn)*I4[212];
*(vp++) = U01*I0[299] + U41*I1[299];
*(vp++) = U01*I0[300] + U41*I1[300]
           + (threeo2zn)*I4[213];
*(vp++) = U01*I0[301] + U41*I1[301]
           + (twoo2zn)*I4[214];
*(vp++) = U01*I0[302] + U41*I1[302]
           + (oneo2zn)*I4[215];
*(vp++) = U01*I0[303] + U41*I1[303];
*(vp++) = U01*I0[304] + U41*I1[304]
           + (fouro2zn)*I4[216];
*(vp++) = U01*I0[305] + U41*I1[305]
           + (threeo2zn)*I4[217];
*(vp++) = U01*I0[306] + U41*I1[306]
           + (twoo2zn)*I4[218];
*(vp++) = U01*I0[307] + U41*I1[307]
           + (oneo2zn)*I4[219];
*(vp++) = U01*I0[308] + U41*I1[308];
*(vp++) = U01*I0[309] + U41*I1[309]
           + (fiveo2zn)*I4[220];
*(vp++) = U01*I0[310] + U41*I1[310]
           + (fouro2zn)*I4[221];
*(vp++) = U01*I0[311] + U41*I1[311]
           + (threeo2zn)*I4[222];
*(vp++) = U01*I0[312] + U41*I1[312]
           + (twoo2zn)*I4[223];
*(vp++) = U01*I0[313] + U41*I1[313]
           + (oneo2zn)*I4[224];
*(vp++) = U01*I0[314] + U41*I1[314];
*(vp++) = U02*I0[294] + U42*I1[294]
           + (fouro2z)*(I2[189] - (lpoz)*I3[189]);
*(vp++) = U02*I0[295] + U42*I1[295]
           + (fouro2z)*(I2[190] - (lpoz)*I3[190]);
*(vp++) = U02*I0[296] + U42*I1[296]
           + (fouro2z)*(I2[191] - (lpoz)*I3[191])
           + (oneo2zn)*I4[210];
*(vp++) = U02*I0[297] + U42*I1[297]
           + (fouro2z)*(I2[192] - (lpoz)*I3[192]);
*(vp++) = U02*I0[298] + U42*I1[298]
           + (fouro2z)*(I2[193] - (lpoz)*I3[193])
           + (oneo2zn)*I4[211];
*(vp++) = U02*I0[299] + U42*I1[299]
           + (fouro2z)*(I2[194] - (lpoz)*I3[194])
           + (twoo2zn)*I4[212];
*(vp++) = U02*I0[300] + U42*I1[300]
           + (fouro2z)*(I2[195] - (lpoz)*I3[195]);
*(vp++) = U02*I0[301] + U42*I1[301]
           + (fouro2z)*(I2[196] - (lpoz)*I3[196])
           + (oneo2zn)*I4[213];
*(vp++) = U02*I0[302] + U42*I1[302]
           + (fouro2z)*(I2[197] - (lpoz)*I3[197])
           + (twoo2zn)*I4[214];
*(vp++) = U02*I0[303] + U42*I1[303]
           + (fouro2z)*(I2[198] - (lpoz)*I3[198])
           + (threeo2zn)*I4[215];
*(vp++) = U02*I0[304] + U42*I1[304]
           + (fouro2z)*(I2[199] - (lpoz)*I3[199]);
*(vp++) = U02*I0[305] + U42*I1[305]
           + (fouro2z)*(I2[200] - (lpoz)*I3[200])
           + (oneo2zn)*I4[216];
*(vp++) = U02*I0[306] + U42*I1[306]
           + (fouro2z)*(I2[201] - (lpoz)*I3[201])
           + (twoo2zn)*I4[217];
*(vp++) = U02*I0[307] + U42*I1[307]
           + (fouro2z)*(I2[202] - (lpoz)*I3[202])
           + (threeo2zn)*I4[218];
*(vp++) = U02*I0[308] + U42*I1[308]
           + (fouro2z)*(I2[203] - (lpoz)*I3[203])
           + (fouro2zn)*I4[219];
*(vp++) = U02*I0[309] + U42*I1[309]
           + (fouro2z)*(I2[204] - (lpoz)*I3[204]);
*(vp++) = U02*I0[310] + U42*I1[310]
           + (fouro2z)*(I2[205] - (lpoz)*I3[205])
           + (oneo2zn)*I4[220];
*(vp++) = U02*I0[311] + U42*I1[311]
           + (fouro2z)*(I2[206] - (lpoz)*I3[206])
           + (twoo2zn)*I4[221];
*(vp++) = U02*I0[312] + U42*I1[312]
           + (fouro2z)*(I2[207] - (lpoz)*I3[207])
           + (threeo2zn)*I4[222];
*(vp++) = U02*I0[313] + U42*I1[313]
           + (fouro2z)*(I2[208] - (lpoz)*I3[208])
           + (fouro2zn)*I4[223];
*(vp++) = U02*I0[314] + U42*I1[314]
           + (fouro2z)*(I2[209] - (lpoz)*I3[209])
           + (fiveo2zn)*I4[224];
return vp;
}
/* Total number of FLOPs = 3213 */
