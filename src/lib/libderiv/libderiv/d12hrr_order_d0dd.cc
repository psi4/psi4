#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_d0dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|dd) integrals */

void d12hrr_order_d0dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][10] = int_stack + 90;
 Libderiv->deriv_classes[2][4][9] = int_stack + 180;
 Libderiv->deriv_classes[2][4][8] = int_stack + 270;
 Libderiv->deriv_classes[2][4][7] = int_stack + 360;
 Libderiv->dvrr_classes[2][3] = int_stack + 450;
 Libderiv->deriv_classes[2][4][6] = int_stack + 510;
 Libderiv->deriv_classes[2][4][2] = int_stack + 600;
 Libderiv->deriv_classes[2][4][1] = int_stack + 690;
 Libderiv->deriv_classes[2][4][0] = int_stack + 780;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 870;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 906;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 966;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 1056;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1092;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 1152;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 1242;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 1278;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 1338;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 1428;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 1464;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 1524;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 1614;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 1650;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 1710;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 1800;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 1836;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 1896;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 1986;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 2022;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 2082;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 2172;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 2208;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 2268;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 2358;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 2394;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 2454;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 2544;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 2580;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 2640;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 2730;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 2766;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 2826;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 2916;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 2952;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 3012;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 3102;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 3138;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 3198;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 3288;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 3324;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 3384;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 3474;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 3510;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 3570;
 Libderiv->deriv_classes[2][2][11] = int_stack + 3660;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 3696;
 Libderiv->deriv_classes[2][3][11] = int_stack + 3732;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 3792;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 3852;
 Libderiv->deriv_classes[2][2][10] = int_stack + 3942;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 3978;
 Libderiv->deriv_classes[2][3][10] = int_stack + 4014;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 4074;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 4134;
 Libderiv->deriv_classes[2][2][9] = int_stack + 4224;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 4260;
 Libderiv->deriv_classes[2][3][9] = int_stack + 4296;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 4356;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 4416;
 Libderiv->deriv_classes[2][2][8] = int_stack + 4506;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 4542;
 Libderiv->deriv_classes[2][3][8] = int_stack + 4578;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 4638;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 4698;
 Libderiv->deriv_classes[2][2][7] = int_stack + 4788;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 4824;
 Libderiv->deriv_classes[2][3][7] = int_stack + 4860;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 4920;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 4980;
 Libderiv->dvrr_classes[2][2] = int_stack + 5070;
 Libderiv->deriv_classes[2][2][6] = int_stack + 5106;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 5142;
 Libderiv->deriv_classes[2][3][6] = int_stack + 5178;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 5238;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 5298;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 5388;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 5424;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 5484;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 5574;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 5610;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 5670;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 5760;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 5796;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 5856;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 5946;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 5982;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 6042;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 6132;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 6168;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 6228;
 Libderiv->deriv_classes[2][2][2] = int_stack + 6318;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 6354;
 Libderiv->deriv_classes[2][3][2] = int_stack + 6390;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 6450;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 6510;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 6600;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 6636;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 6696;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 6786;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 6822;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 6882;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 6972;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 7008;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 7068;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 7158;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 7194;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 7254;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 7344;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 7380;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 7440;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 7530;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 7566;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 7626;
 Libderiv->deriv_classes[2][2][1] = int_stack + 7716;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 7752;
 Libderiv->deriv_classes[2][3][1] = int_stack + 7788;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 7848;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 7908;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 7998;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 8034;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 8094;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 8184;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 8220;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 8280;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 8370;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 8406;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 8466;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 8556;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 8592;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 8652;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 8742;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 8778;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 8838;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 8928;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 8964;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 9024;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 9114;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 9150;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 9210;
 Libderiv->deriv_classes[2][2][0] = int_stack + 9300;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 9336;
 Libderiv->deriv_classes[2][3][0] = int_stack + 9372;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 9432;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 9492;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 9582;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 9618;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 9678;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 9768;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 9804;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 9864;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 9954;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 9990;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 10050;
 memset(int_stack,0,81120);

 Libderiv->dvrr_stack = int_stack + 16512;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_d0dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+10140,int_stack+450,int_stack+5070,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10248,int_stack+3732,int_stack+3660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5070,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10356,int_stack+0,int_stack+3732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10536,int_stack+4014,int_stack+3942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5070, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10644,int_stack+90,int_stack+4014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+4296,int_stack+4224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5070, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10824,int_stack+180,int_stack+4296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+4578,int_stack+4506, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11004,int_stack+270,int_stack+4578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+4860,int_stack+4788, 0.0, zero_stack, 1.0, int_stack+5070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11184,int_stack+360,int_stack+4860, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+5178,int_stack+5106, 1.0, int_stack+5070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11364,int_stack+510,int_stack+5178, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+6390,int_stack+6318,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11544,int_stack+600,int_stack+6390,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+7788,int_stack+7716,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11724,int_stack+690,int_stack+7788,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+9372,int_stack+9300,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11904,int_stack+780,int_stack+9372,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+756,int_stack+906,int_stack+870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3660,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12084,int_stack+966,int_stack+906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3732,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+864,int_stack+1092,int_stack+1056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3660, 1.0, int_stack+3942,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12264,int_stack+1152,int_stack+1092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3732, 1.0, int_stack+4014,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+972,int_stack+1278,int_stack+1242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3942, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+1338,int_stack+1278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4014, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1260,int_stack+1464,int_stack+1428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3660, 0.0, zero_stack, 1.0, int_stack+4224,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12444,int_stack+1524,int_stack+1464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3732, 0.0, zero_stack, 1.0, int_stack+4296,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1368,int_stack+1650,int_stack+1614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3942, 1.0, int_stack+4224, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12624,int_stack+1710,int_stack+1650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4014, 1.0, int_stack+4296, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1476,int_stack+1836,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4224, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1584,int_stack+1896,int_stack+1836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4296, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1764,int_stack+2022,int_stack+1986, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4506,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12804,int_stack+2082,int_stack+2022, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3732, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4578,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1872,int_stack+2208,int_stack+2172, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3942, 0.0, zero_stack, 1.0, int_stack+4506, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1980,int_stack+2268,int_stack+2208, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4014, 0.0, zero_stack, 1.0, int_stack+4578, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2160,int_stack+2394,int_stack+2358, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4224, 1.0, int_stack+4506, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12984,int_stack+2454,int_stack+2394, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4296, 1.0, int_stack+4578, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2268,int_stack+2580,int_stack+2544, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2376,int_stack+2640,int_stack+2580, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2556,int_stack+2766,int_stack+2730, 0.0, zero_stack, 1.0, int_stack+3660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4788,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13164,int_stack+2826,int_stack+2766, 0.0, zero_stack, 1.0, int_stack+3732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4860,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2664,int_stack+2952,int_stack+2916, 0.0, zero_stack, 1.0, int_stack+3942, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4788, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2772,int_stack+3012,int_stack+2952, 0.0, zero_stack, 1.0, int_stack+4014, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4860, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2952,int_stack+3138,int_stack+3102, 0.0, zero_stack, 1.0, int_stack+4224, 0.0, zero_stack, 1.0, int_stack+4788, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13344,int_stack+3198,int_stack+3138, 0.0, zero_stack, 1.0, int_stack+4296, 0.0, zero_stack, 1.0, int_stack+4860, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3060,int_stack+3324,int_stack+3288, 0.0, zero_stack, 1.0, int_stack+4506, 1.0, int_stack+4788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13524,int_stack+3384,int_stack+3324, 0.0, zero_stack, 1.0, int_stack+4578, 1.0, int_stack+4860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3168,int_stack+3510,int_stack+3474, 0.0, zero_stack, 2.0, int_stack+4788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3276,int_stack+3570,int_stack+3510, 0.0, zero_stack, 2.0, int_stack+4860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3456,int_stack+3792,int_stack+3696, 1.0, int_stack+3660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5106,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13704,int_stack+3852,int_stack+3792, 1.0, int_stack+3732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5178,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3564,int_stack+4074,int_stack+3978, 1.0, int_stack+3942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5106, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3672,int_stack+4134,int_stack+4074, 1.0, int_stack+4014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5178, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3852,int_stack+4356,int_stack+4260, 1.0, int_stack+4224, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5106, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3960,int_stack+4416,int_stack+4356, 1.0, int_stack+4296, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5178, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4140,int_stack+4638,int_stack+4542, 1.0, int_stack+4506, 0.0, zero_stack, 1.0, int_stack+5106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4248,int_stack+4698,int_stack+4638, 1.0, int_stack+4578, 0.0, zero_stack, 1.0, int_stack+5178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4428,int_stack+4920,int_stack+4824, 1.0, int_stack+4788, 1.0, int_stack+5106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4536,int_stack+4980,int_stack+4920, 1.0, int_stack+4860, 1.0, int_stack+5178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4716,int_stack+5238,int_stack+5142, 2.0, int_stack+5106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4824,int_stack+5298,int_stack+5238, 2.0, int_stack+5178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5004,int_stack+5424,int_stack+5388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6318,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5112,int_stack+5484,int_stack+5424, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6390,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5292,int_stack+5610,int_stack+5574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6318, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5400,int_stack+5670,int_stack+5610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6390, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5580,int_stack+5796,int_stack+5760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6318, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13884,int_stack+5856,int_stack+5796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6390, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5688,int_stack+5982,int_stack+5946, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5796,int_stack+6042,int_stack+5982, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5976,int_stack+6168,int_stack+6132, 0.0, zero_stack, 1.0, int_stack+6318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14064,int_stack+6228,int_stack+6168, 0.0, zero_stack, 1.0, int_stack+6390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6084,int_stack+6450,int_stack+6354, 1.0, int_stack+6318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6192,int_stack+6510,int_stack+6450, 1.0, int_stack+6390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+6372,int_stack+6636,int_stack+6600,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14244,int_stack+6696,int_stack+6636,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6480,int_stack+6822,int_stack+6786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7716,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6588,int_stack+6882,int_stack+6822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7788,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6768,int_stack+7008,int_stack+6972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7716, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14424,int_stack+7068,int_stack+7008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7788, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6876,int_stack+7194,int_stack+7158, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7716, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6984,int_stack+7254,int_stack+7194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7788, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7164,int_stack+7380,int_stack+7344, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14604,int_stack+7440,int_stack+7380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7272,int_stack+7566,int_stack+7530, 0.0, zero_stack, 1.0, int_stack+7716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7380,int_stack+7626,int_stack+7566, 0.0, zero_stack, 1.0, int_stack+7788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7560,int_stack+7848,int_stack+7752, 1.0, int_stack+7716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14784,int_stack+7908,int_stack+7848, 1.0, int_stack+7788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7668,int_stack+8034,int_stack+7998,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7776,int_stack+8094,int_stack+8034,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7956,int_stack+8220,int_stack+8184,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14964,int_stack+8280,int_stack+8220,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8064,int_stack+8406,int_stack+8370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9300,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8172,int_stack+8466,int_stack+8406, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9372,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8352,int_stack+8592,int_stack+8556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9300, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15144,int_stack+8652,int_stack+8592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9372, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8460,int_stack+8778,int_stack+8742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9300, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8568,int_stack+8838,int_stack+8778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9372, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8748,int_stack+8964,int_stack+8928, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15324,int_stack+9024,int_stack+8964, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8856,int_stack+9150,int_stack+9114, 0.0, zero_stack, 1.0, int_stack+9300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8964,int_stack+9210,int_stack+9150, 0.0, zero_stack, 1.0, int_stack+9372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9144,int_stack+9432,int_stack+9336, 1.0, int_stack+9300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15504,int_stack+9492,int_stack+9432, 1.0, int_stack+9372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+9252,int_stack+9618,int_stack+9582,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+9360,int_stack+9678,int_stack+9618,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+9540,int_stack+9804,int_stack+9768,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15684,int_stack+9864,int_stack+9804,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+9648,int_stack+9990,int_stack+9954,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+9756,int_stack+10050,int_stack+9990,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15864,int_stack+10356,int_stack+10248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10140,6);
     Libderiv->ABCD[11] = int_stack + 15864;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16080,int_stack+10644,int_stack+10536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10140, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 16080;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16296,int_stack+10824,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10140, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 16296;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10644,int_stack+11004,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 10644;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10860,int_stack+11184,int_stack+216, 0.0, zero_stack, 1.0, int_stack+10140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 10860;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11076,int_stack+11364,int_stack+324, 1.0, int_stack+10140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 11076;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+11292,int_stack+11544,int_stack+432,6);
     Libderiv->ABCD[2] = int_stack + 11292;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+11508,int_stack+11724,int_stack+540,6);
     Libderiv->ABCD[1] = int_stack + 11508;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+9936,int_stack+11904,int_stack+648,6);
     Libderiv->ABCD[0] = int_stack + 9936;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11724,int_stack+12084,int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10248,6);
     Libderiv->ABCD[155] = int_stack + 11724;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11940,int_stack+12264,int_stack+864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10248, 1.0, int_stack+10536,6);
     Libderiv->ABCD[143] = int_stack + 11940;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+756,int_stack+1080,int_stack+972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10536, 0.0, zero_stack,6);
     Libderiv->ABCD[142] = int_stack + 756;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+972,int_stack+12444,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10248, 0.0, zero_stack, 1.0, int_stack+0,6);
     Libderiv->ABCD[131] = int_stack + 972;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12156,int_stack+12624,int_stack+1368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10536, 1.0, int_stack+0, 0.0, zero_stack,6);
     Libderiv->ABCD[130] = int_stack + 12156;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1188,int_stack+1584,int_stack+1476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[129] = int_stack + 1188;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1404,int_stack+12804,int_stack+1764, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10248, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108,6);
     Libderiv->ABCD[119] = int_stack + 1404;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1620,int_stack+1980,int_stack+1872, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10536, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack,6);
     Libderiv->ABCD[118] = int_stack + 1620;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1836,int_stack+12984,int_stack+2160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[117] = int_stack + 1836;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2052,int_stack+2376,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[116] = int_stack + 2052;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2268,int_stack+13164,int_stack+2556, 0.0, zero_stack, 1.0, int_stack+10248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216,6);
     Libderiv->ABCD[107] = int_stack + 2268;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12372,int_stack+2772,int_stack+2664, 0.0, zero_stack, 1.0, int_stack+10536, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack,6);
     Libderiv->ABCD[106] = int_stack + 12372;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2484,int_stack+13344,int_stack+2952, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[105] = int_stack + 2484;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2700,int_stack+13524,int_stack+3060, 0.0, zero_stack, 1.0, int_stack+108, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[104] = int_stack + 2700;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2916,int_stack+3276,int_stack+3168, 0.0, zero_stack, 2.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[103] = int_stack + 2916;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3132,int_stack+13704,int_stack+3456, 1.0, int_stack+10248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324,6);
     Libderiv->ABCD[95] = int_stack + 3132;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3348,int_stack+3672,int_stack+3564, 1.0, int_stack+10536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack,6);
     Libderiv->ABCD[94] = int_stack + 3348;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3564,int_stack+3960,int_stack+3852, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[93] = int_stack + 3564;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3780,int_stack+4248,int_stack+4140, 1.0, int_stack+108, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[92] = int_stack + 3780;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+4536,int_stack+4428, 1.0, int_stack+216, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[91] = int_stack + 0;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3996,int_stack+4824,int_stack+4716, 2.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[90] = int_stack + 3996;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+216,int_stack+5112,int_stack+5004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432,6);
     Libderiv->ABCD[47] = int_stack + 216;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4212,int_stack+5400,int_stack+5292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack,6);
     Libderiv->ABCD[46] = int_stack + 4212;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4428,int_stack+13884,int_stack+5580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[45] = int_stack + 4428;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4644,int_stack+5796,int_stack+5688, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[44] = int_stack + 4644;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4860,int_stack+14064,int_stack+5976, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[43] = int_stack + 4860;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5076,int_stack+6192,int_stack+6084, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[42] = int_stack + 5076;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+5292,int_stack+14244,int_stack+6372,6);
     Libderiv->ABCD[38] = int_stack + 5292;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5508,int_stack+6588,int_stack+6480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540,6);
     Libderiv->ABCD[35] = int_stack + 5508;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5724,int_stack+14424,int_stack+6768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack,6);
     Libderiv->ABCD[34] = int_stack + 5724;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5940,int_stack+6984,int_stack+6876, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[33] = int_stack + 5940;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6156,int_stack+14604,int_stack+7164, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[32] = int_stack + 6156;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6372,int_stack+7380,int_stack+7272, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[31] = int_stack + 6372;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6588,int_stack+14784,int_stack+7560, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[30] = int_stack + 6588;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+432,int_stack+7776,int_stack+7668,6);
     Libderiv->ABCD[26] = int_stack + 432;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6804,int_stack+14964,int_stack+7956,6);
     Libderiv->ABCD[25] = int_stack + 6804;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7020,int_stack+8172,int_stack+8064, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648,6);
     Libderiv->ABCD[23] = int_stack + 7020;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7236,int_stack+15144,int_stack+8352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack,6);
     Libderiv->ABCD[22] = int_stack + 7236;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7452,int_stack+8568,int_stack+8460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[21] = int_stack + 7452;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7668,int_stack+15324,int_stack+8748, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[20] = int_stack + 7668;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7884,int_stack+8964,int_stack+8856, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[19] = int_stack + 7884;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8100,int_stack+15504,int_stack+9144, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[18] = int_stack + 8100;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+8316,int_stack+9360,int_stack+9252,6);
     Libderiv->ABCD[14] = int_stack + 8316;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+8532,int_stack+15684,int_stack+9540,6);
     Libderiv->ABCD[13] = int_stack + 8532;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+8748,int_stack+9756,int_stack+9648,6);
     Libderiv->ABCD[12] = int_stack + 8748;

}
