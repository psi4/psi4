#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0dp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|dp) integrals */

void d12hrr_order_p0dp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[1][3][10] = int_stack + 30;
 Libderiv->deriv_classes[1][3][9] = int_stack + 60;
 Libderiv->deriv_classes[1][3][8] = int_stack + 90;
 Libderiv->deriv_classes[1][3][7] = int_stack + 120;
 Libderiv->dvrr_classes[1][2] = int_stack + 150;
 Libderiv->deriv_classes[1][3][6] = int_stack + 168;
 Libderiv->deriv_classes[1][3][2] = int_stack + 198;
 Libderiv->deriv_classes[1][3][1] = int_stack + 228;
 Libderiv->deriv_classes[1][3][0] = int_stack + 258;
 Libderiv->deriv2_classes[1][2][143] = int_stack + 288;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 306;
 Libderiv->deriv2_classes[1][2][131] = int_stack + 336;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 354;
 Libderiv->deriv2_classes[1][2][130] = int_stack + 384;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 402;
 Libderiv->deriv2_classes[1][2][119] = int_stack + 432;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 450;
 Libderiv->deriv2_classes[1][2][118] = int_stack + 480;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 498;
 Libderiv->deriv2_classes[1][2][117] = int_stack + 528;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 546;
 Libderiv->deriv2_classes[1][2][107] = int_stack + 576;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 594;
 Libderiv->deriv2_classes[1][2][106] = int_stack + 624;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 642;
 Libderiv->deriv2_classes[1][2][105] = int_stack + 672;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 690;
 Libderiv->deriv2_classes[1][2][104] = int_stack + 720;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 738;
 Libderiv->deriv2_classes[1][2][95] = int_stack + 768;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 786;
 Libderiv->deriv2_classes[1][2][94] = int_stack + 816;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 834;
 Libderiv->deriv2_classes[1][2][93] = int_stack + 864;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 882;
 Libderiv->deriv2_classes[1][2][92] = int_stack + 912;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 930;
 Libderiv->deriv2_classes[1][2][91] = int_stack + 960;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 978;
 Libderiv->deriv_classes[1][2][11] = int_stack + 1008;
 Libderiv->deriv2_classes[1][2][83] = int_stack + 1026;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 1044;
 Libderiv->deriv_classes[1][2][10] = int_stack + 1074;
 Libderiv->deriv2_classes[1][2][82] = int_stack + 1092;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 1110;
 Libderiv->deriv_classes[1][2][9] = int_stack + 1140;
 Libderiv->deriv2_classes[1][2][81] = int_stack + 1158;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 1176;
 Libderiv->deriv_classes[1][2][8] = int_stack + 1206;
 Libderiv->deriv2_classes[1][2][80] = int_stack + 1224;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 1242;
 Libderiv->deriv_classes[1][2][7] = int_stack + 1272;
 Libderiv->deriv2_classes[1][2][79] = int_stack + 1290;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 1308;
 Libderiv->deriv_classes[1][2][6] = int_stack + 1338;
 Libderiv->deriv2_classes[1][2][78] = int_stack + 1356;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 1374;
 Libderiv->deriv2_classes[1][2][35] = int_stack + 1404;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 1422;
 Libderiv->deriv2_classes[1][2][34] = int_stack + 1452;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 1470;
 Libderiv->deriv2_classes[1][2][33] = int_stack + 1500;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 1518;
 Libderiv->deriv2_classes[1][2][32] = int_stack + 1548;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 1566;
 Libderiv->deriv2_classes[1][2][31] = int_stack + 1596;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 1614;
 Libderiv->deriv_classes[1][2][2] = int_stack + 1644;
 Libderiv->deriv2_classes[1][2][30] = int_stack + 1662;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 1680;
 Libderiv->deriv2_classes[1][2][26] = int_stack + 1710;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 1728;
 Libderiv->deriv2_classes[1][2][23] = int_stack + 1758;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 1776;
 Libderiv->deriv2_classes[1][2][22] = int_stack + 1806;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 1824;
 Libderiv->deriv2_classes[1][2][21] = int_stack + 1854;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 1872;
 Libderiv->deriv2_classes[1][2][20] = int_stack + 1902;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 1920;
 Libderiv->deriv2_classes[1][2][19] = int_stack + 1950;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 1968;
 Libderiv->deriv_classes[1][2][1] = int_stack + 1998;
 Libderiv->deriv2_classes[1][2][18] = int_stack + 2016;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 2034;
 Libderiv->deriv2_classes[1][2][14] = int_stack + 2064;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 2082;
 Libderiv->deriv2_classes[1][2][13] = int_stack + 2112;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 2130;
 Libderiv->deriv2_classes[1][2][11] = int_stack + 2160;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 2178;
 Libderiv->deriv2_classes[1][2][10] = int_stack + 2208;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 2226;
 Libderiv->deriv2_classes[1][2][9] = int_stack + 2256;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 2274;
 Libderiv->deriv2_classes[1][2][8] = int_stack + 2304;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 2322;
 Libderiv->deriv2_classes[1][2][7] = int_stack + 2352;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 2370;
 Libderiv->deriv_classes[1][2][0] = int_stack + 2400;
 Libderiv->deriv2_classes[1][2][6] = int_stack + 2418;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 2436;
 Libderiv->deriv2_classes[1][2][2] = int_stack + 2466;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 2484;
 Libderiv->deriv2_classes[1][2][1] = int_stack + 2514;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 2532;
 Libderiv->deriv2_classes[1][2][0] = int_stack + 2562;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 2580;
 memset(int_stack,0,20880);

 Libderiv->dvrr_stack = int_stack + 2988;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0dp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2610,int_stack+0,int_stack+1008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150,3);
     Libderiv->ABCD[11] = int_stack + 2610;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2664,int_stack+30,int_stack+1074, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 2664;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+60,int_stack+1140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2718,int_stack+90,int_stack+1206, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 2718;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+120,int_stack+1272, 0.0, zero_stack, 1.0, int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 54;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2772,int_stack+168,int_stack+1338, 1.0, int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 2772;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+198,int_stack+1644,3);
     Libderiv->ABCD[2] = int_stack + 108;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+162,int_stack+228,int_stack+1998,3);
     Libderiv->ABCD[1] = int_stack + 162;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2826,int_stack+258,int_stack+2400,3);
     Libderiv->ABCD[0] = int_stack + 2826;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+306,int_stack+288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1008,3);
     Libderiv->ABCD[155] = int_stack + 216;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+354,int_stack+336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1008, 1.0, int_stack+1074,3);
     Libderiv->ABCD[143] = int_stack + 270;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+402,int_stack+384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1074, 0.0, zero_stack,3);
     Libderiv->ABCD[142] = int_stack + 324;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+378,int_stack+450,int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1008, 0.0, zero_stack, 1.0, int_stack+1140,3);
     Libderiv->ABCD[131] = int_stack + 378;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2880,int_stack+498,int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1074, 1.0, int_stack+1140, 0.0, zero_stack,3);
     Libderiv->ABCD[130] = int_stack + 2880;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+546,int_stack+528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1140, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[129] = int_stack + 432;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+486,int_stack+594,int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1008, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1206,3);
     Libderiv->ABCD[119] = int_stack + 486;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+642,int_stack+624, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1074, 0.0, zero_stack, 1.0, int_stack+1206, 0.0, zero_stack,3);
     Libderiv->ABCD[118] = int_stack + 540;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+594,int_stack+690,int_stack+672, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1140, 1.0, int_stack+1206, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[117] = int_stack + 594;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+738,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[116] = int_stack + 648;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+702,int_stack+786,int_stack+768, 0.0, zero_stack, 1.0, int_stack+1008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1272,3);
     Libderiv->ABCD[107] = int_stack + 702;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+756,int_stack+834,int_stack+816, 0.0, zero_stack, 1.0, int_stack+1074, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1272, 0.0, zero_stack,3);
     Libderiv->ABCD[106] = int_stack + 756;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+810,int_stack+882,int_stack+864, 0.0, zero_stack, 1.0, int_stack+1140, 0.0, zero_stack, 1.0, int_stack+1272, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[105] = int_stack + 810;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2934,int_stack+930,int_stack+912, 0.0, zero_stack, 1.0, int_stack+1206, 1.0, int_stack+1272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[104] = int_stack + 2934;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+864,int_stack+978,int_stack+960, 0.0, zero_stack, 2.0, int_stack+1272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[103] = int_stack + 864;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+918,int_stack+1044,int_stack+1026, 1.0, int_stack+1008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1338,3);
     Libderiv->ABCD[95] = int_stack + 918;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+972,int_stack+1110,int_stack+1092, 1.0, int_stack+1074, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1338, 0.0, zero_stack,3);
     Libderiv->ABCD[94] = int_stack + 972;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1026,int_stack+1176,int_stack+1158, 1.0, int_stack+1140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1338, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[93] = int_stack + 1026;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1242,int_stack+1224, 1.0, int_stack+1206, 0.0, zero_stack, 1.0, int_stack+1338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[92] = int_stack + 1080;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1134,int_stack+1308,int_stack+1290, 1.0, int_stack+1272, 1.0, int_stack+1338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[91] = int_stack + 1134;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1188,int_stack+1374,int_stack+1356, 2.0, int_stack+1338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[90] = int_stack + 1188;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1242,int_stack+1422,int_stack+1404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1644,3);
     Libderiv->ABCD[47] = int_stack + 1242;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1296,int_stack+1470,int_stack+1452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1644, 0.0, zero_stack,3);
     Libderiv->ABCD[46] = int_stack + 1296;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1350,int_stack+1518,int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1644, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[45] = int_stack + 1350;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1404,int_stack+1566,int_stack+1548, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[44] = int_stack + 1404;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1458,int_stack+1614,int_stack+1596, 0.0, zero_stack, 1.0, int_stack+1644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[43] = int_stack + 1458;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1512,int_stack+1680,int_stack+1662, 1.0, int_stack+1644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[42] = int_stack + 1512;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1566,int_stack+1728,int_stack+1710,3);
     Libderiv->ABCD[38] = int_stack + 1566;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1620,int_stack+1776,int_stack+1758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1998,3);
     Libderiv->ABCD[35] = int_stack + 1620;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1674,int_stack+1824,int_stack+1806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1998, 0.0, zero_stack,3);
     Libderiv->ABCD[34] = int_stack + 1674;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1728,int_stack+1872,int_stack+1854, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1998, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[33] = int_stack + 1728;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1782,int_stack+1920,int_stack+1902, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[32] = int_stack + 1782;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1836,int_stack+1968,int_stack+1950, 0.0, zero_stack, 1.0, int_stack+1998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[31] = int_stack + 1836;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1890,int_stack+2034,int_stack+2016, 1.0, int_stack+1998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[30] = int_stack + 1890;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1944,int_stack+2082,int_stack+2064,3);
     Libderiv->ABCD[26] = int_stack + 1944;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1998,int_stack+2130,int_stack+2112,3);
     Libderiv->ABCD[25] = int_stack + 1998;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2052,int_stack+2178,int_stack+2160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2400,3);
     Libderiv->ABCD[23] = int_stack + 2052;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2106,int_stack+2226,int_stack+2208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2400, 0.0, zero_stack,3);
     Libderiv->ABCD[22] = int_stack + 2106;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2160,int_stack+2274,int_stack+2256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2400, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[21] = int_stack + 2160;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2214,int_stack+2322,int_stack+2304, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[20] = int_stack + 2214;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2268,int_stack+2370,int_stack+2352, 0.0, zero_stack, 1.0, int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[19] = int_stack + 2268;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2322,int_stack+2436,int_stack+2418, 1.0, int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[18] = int_stack + 2322;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2376,int_stack+2484,int_stack+2466,3);
     Libderiv->ABCD[14] = int_stack + 2376;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2430,int_stack+2532,int_stack+2514,3);
     Libderiv->ABCD[13] = int_stack + 2430;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2484,int_stack+2580,int_stack+2562,3);
     Libderiv->ABCD[12] = int_stack + 2484;

}
