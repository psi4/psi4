#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|dd) integrals */

void d12hrr_order_p0dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[1][4][10] = int_stack + 45;
 Libderiv->deriv_classes[1][4][9] = int_stack + 90;
 Libderiv->deriv_classes[1][4][8] = int_stack + 135;
 Libderiv->deriv_classes[1][4][7] = int_stack + 180;
 Libderiv->dvrr_classes[1][3] = int_stack + 225;
 Libderiv->deriv_classes[1][4][6] = int_stack + 255;
 Libderiv->deriv_classes[1][4][2] = int_stack + 300;
 Libderiv->deriv_classes[1][4][1] = int_stack + 345;
 Libderiv->deriv_classes[1][4][0] = int_stack + 390;
 Libderiv->deriv2_classes[1][2][143] = int_stack + 435;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 453;
 Libderiv->deriv2_classes[1][4][143] = int_stack + 483;
 Libderiv->deriv2_classes[1][2][131] = int_stack + 528;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 546;
 Libderiv->deriv2_classes[1][4][131] = int_stack + 576;
 Libderiv->deriv2_classes[1][2][130] = int_stack + 621;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 639;
 Libderiv->deriv2_classes[1][4][130] = int_stack + 669;
 Libderiv->deriv2_classes[1][2][119] = int_stack + 714;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 732;
 Libderiv->deriv2_classes[1][4][119] = int_stack + 762;
 Libderiv->deriv2_classes[1][2][118] = int_stack + 807;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 825;
 Libderiv->deriv2_classes[1][4][118] = int_stack + 855;
 Libderiv->deriv2_classes[1][2][117] = int_stack + 900;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 918;
 Libderiv->deriv2_classes[1][4][117] = int_stack + 948;
 Libderiv->deriv2_classes[1][2][107] = int_stack + 993;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 1011;
 Libderiv->deriv2_classes[1][4][107] = int_stack + 1041;
 Libderiv->deriv2_classes[1][2][106] = int_stack + 1086;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 1104;
 Libderiv->deriv2_classes[1][4][106] = int_stack + 1134;
 Libderiv->deriv2_classes[1][2][105] = int_stack + 1179;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 1197;
 Libderiv->deriv2_classes[1][4][105] = int_stack + 1227;
 Libderiv->deriv2_classes[1][2][104] = int_stack + 1272;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 1290;
 Libderiv->deriv2_classes[1][4][104] = int_stack + 1320;
 Libderiv->deriv2_classes[1][2][95] = int_stack + 1365;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 1383;
 Libderiv->deriv2_classes[1][4][95] = int_stack + 1413;
 Libderiv->deriv2_classes[1][2][94] = int_stack + 1458;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 1476;
 Libderiv->deriv2_classes[1][4][94] = int_stack + 1506;
 Libderiv->deriv2_classes[1][2][93] = int_stack + 1551;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 1569;
 Libderiv->deriv2_classes[1][4][93] = int_stack + 1599;
 Libderiv->deriv2_classes[1][2][92] = int_stack + 1644;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 1662;
 Libderiv->deriv2_classes[1][4][92] = int_stack + 1692;
 Libderiv->deriv2_classes[1][2][91] = int_stack + 1737;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 1755;
 Libderiv->deriv2_classes[1][4][91] = int_stack + 1785;
 Libderiv->deriv_classes[1][2][11] = int_stack + 1830;
 Libderiv->deriv2_classes[1][2][83] = int_stack + 1848;
 Libderiv->deriv_classes[1][3][11] = int_stack + 1866;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 1896;
 Libderiv->deriv2_classes[1][4][83] = int_stack + 1926;
 Libderiv->deriv_classes[1][2][10] = int_stack + 1971;
 Libderiv->deriv2_classes[1][2][82] = int_stack + 1989;
 Libderiv->deriv_classes[1][3][10] = int_stack + 2007;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 2037;
 Libderiv->deriv2_classes[1][4][82] = int_stack + 2067;
 Libderiv->deriv_classes[1][2][9] = int_stack + 2112;
 Libderiv->deriv2_classes[1][2][81] = int_stack + 2130;
 Libderiv->deriv_classes[1][3][9] = int_stack + 2148;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 2178;
 Libderiv->deriv2_classes[1][4][81] = int_stack + 2208;
 Libderiv->deriv_classes[1][2][8] = int_stack + 2253;
 Libderiv->deriv2_classes[1][2][80] = int_stack + 2271;
 Libderiv->deriv_classes[1][3][8] = int_stack + 2289;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 2319;
 Libderiv->deriv2_classes[1][4][80] = int_stack + 2349;
 Libderiv->deriv_classes[1][2][7] = int_stack + 2394;
 Libderiv->deriv2_classes[1][2][79] = int_stack + 2412;
 Libderiv->deriv_classes[1][3][7] = int_stack + 2430;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 2460;
 Libderiv->deriv2_classes[1][4][79] = int_stack + 2490;
 Libderiv->dvrr_classes[1][2] = int_stack + 2535;
 Libderiv->deriv_classes[1][2][6] = int_stack + 2553;
 Libderiv->deriv2_classes[1][2][78] = int_stack + 2571;
 Libderiv->deriv_classes[1][3][6] = int_stack + 2589;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 2619;
 Libderiv->deriv2_classes[1][4][78] = int_stack + 2649;
 Libderiv->deriv2_classes[1][2][35] = int_stack + 2694;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 2712;
 Libderiv->deriv2_classes[1][4][35] = int_stack + 2742;
 Libderiv->deriv2_classes[1][2][34] = int_stack + 2787;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 2805;
 Libderiv->deriv2_classes[1][4][34] = int_stack + 2835;
 Libderiv->deriv2_classes[1][2][33] = int_stack + 2880;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 2898;
 Libderiv->deriv2_classes[1][4][33] = int_stack + 2928;
 Libderiv->deriv2_classes[1][2][32] = int_stack + 2973;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 2991;
 Libderiv->deriv2_classes[1][4][32] = int_stack + 3021;
 Libderiv->deriv2_classes[1][2][31] = int_stack + 3066;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 3084;
 Libderiv->deriv2_classes[1][4][31] = int_stack + 3114;
 Libderiv->deriv_classes[1][2][2] = int_stack + 3159;
 Libderiv->deriv2_classes[1][2][30] = int_stack + 3177;
 Libderiv->deriv_classes[1][3][2] = int_stack + 3195;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 3225;
 Libderiv->deriv2_classes[1][4][30] = int_stack + 3255;
 Libderiv->deriv2_classes[1][2][26] = int_stack + 3300;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 3318;
 Libderiv->deriv2_classes[1][4][26] = int_stack + 3348;
 Libderiv->deriv2_classes[1][2][23] = int_stack + 3393;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 3411;
 Libderiv->deriv2_classes[1][4][23] = int_stack + 3441;
 Libderiv->deriv2_classes[1][2][22] = int_stack + 3486;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 3504;
 Libderiv->deriv2_classes[1][4][22] = int_stack + 3534;
 Libderiv->deriv2_classes[1][2][21] = int_stack + 3579;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 3597;
 Libderiv->deriv2_classes[1][4][21] = int_stack + 3627;
 Libderiv->deriv2_classes[1][2][20] = int_stack + 3672;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 3690;
 Libderiv->deriv2_classes[1][4][20] = int_stack + 3720;
 Libderiv->deriv2_classes[1][2][19] = int_stack + 3765;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 3783;
 Libderiv->deriv2_classes[1][4][19] = int_stack + 3813;
 Libderiv->deriv_classes[1][2][1] = int_stack + 3858;
 Libderiv->deriv2_classes[1][2][18] = int_stack + 3876;
 Libderiv->deriv_classes[1][3][1] = int_stack + 3894;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 3924;
 Libderiv->deriv2_classes[1][4][18] = int_stack + 3954;
 Libderiv->deriv2_classes[1][2][14] = int_stack + 3999;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 4017;
 Libderiv->deriv2_classes[1][4][14] = int_stack + 4047;
 Libderiv->deriv2_classes[1][2][13] = int_stack + 4092;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 4110;
 Libderiv->deriv2_classes[1][4][13] = int_stack + 4140;
 Libderiv->deriv2_classes[1][2][11] = int_stack + 4185;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 4203;
 Libderiv->deriv2_classes[1][4][11] = int_stack + 4233;
 Libderiv->deriv2_classes[1][2][10] = int_stack + 4278;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 4296;
 Libderiv->deriv2_classes[1][4][10] = int_stack + 4326;
 Libderiv->deriv2_classes[1][2][9] = int_stack + 4371;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 4389;
 Libderiv->deriv2_classes[1][4][9] = int_stack + 4419;
 Libderiv->deriv2_classes[1][2][8] = int_stack + 4464;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 4482;
 Libderiv->deriv2_classes[1][4][8] = int_stack + 4512;
 Libderiv->deriv2_classes[1][2][7] = int_stack + 4557;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 4575;
 Libderiv->deriv2_classes[1][4][7] = int_stack + 4605;
 Libderiv->deriv_classes[1][2][0] = int_stack + 4650;
 Libderiv->deriv2_classes[1][2][6] = int_stack + 4668;
 Libderiv->deriv_classes[1][3][0] = int_stack + 4686;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 4716;
 Libderiv->deriv2_classes[1][4][6] = int_stack + 4746;
 Libderiv->deriv2_classes[1][2][2] = int_stack + 4791;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 4809;
 Libderiv->deriv2_classes[1][4][2] = int_stack + 4839;
 Libderiv->deriv2_classes[1][2][1] = int_stack + 4884;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 4902;
 Libderiv->deriv2_classes[1][4][1] = int_stack + 4932;
 Libderiv->deriv2_classes[1][2][0] = int_stack + 4977;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 4995;
 Libderiv->deriv2_classes[1][4][0] = int_stack + 5025;
 memset(int_stack,0,40560);

 Libderiv->dvrr_stack = int_stack + 8256;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+5070,int_stack+225,int_stack+2535,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5124,int_stack+1866,int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2535,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5178,int_stack+0,int_stack+1866, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5268,int_stack+2007,int_stack+1971, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2535, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5322,int_stack+45,int_stack+2007, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+2148,int_stack+2112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2535, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5412,int_stack+90,int_stack+2148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+2289,int_stack+2253, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5502,int_stack+135,int_stack+2289, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+2430,int_stack+2394, 0.0, zero_stack, 1.0, int_stack+2535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5592,int_stack+180,int_stack+2430, 0.0, zero_stack, 1.0, int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+162,int_stack+2589,int_stack+2553, 1.0, int_stack+2535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5682,int_stack+255,int_stack+2589, 1.0, int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+3195,int_stack+3159,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5772,int_stack+300,int_stack+3195,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+3894,int_stack+3858,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5862,int_stack+345,int_stack+3894,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+4686,int_stack+4650,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5952,int_stack+390,int_stack+4686,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+378,int_stack+453,int_stack+435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1830,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6042,int_stack+483,int_stack+453, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1866,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+546,int_stack+528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1830, 1.0, int_stack+1971,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6132,int_stack+576,int_stack+546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1866, 1.0, int_stack+2007,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+486,int_stack+639,int_stack+621, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1971, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+669,int_stack+639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2007, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+630,int_stack+732,int_stack+714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1830, 0.0, zero_stack, 1.0, int_stack+2112,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6222,int_stack+762,int_stack+732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1866, 0.0, zero_stack, 1.0, int_stack+2148,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+684,int_stack+825,int_stack+807, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1971, 1.0, int_stack+2112, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6312,int_stack+855,int_stack+825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2007, 1.0, int_stack+2148, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+738,int_stack+918,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2112, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+792,int_stack+948,int_stack+918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2148, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+882,int_stack+1011,int_stack+993, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2253,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6402,int_stack+1041,int_stack+1011, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1866, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2289,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+936,int_stack+1104,int_stack+1086, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1971, 0.0, zero_stack, 1.0, int_stack+2253, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+990,int_stack+1134,int_stack+1104, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2007, 0.0, zero_stack, 1.0, int_stack+2289, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1197,int_stack+1179, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2112, 1.0, int_stack+2253, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6492,int_stack+1227,int_stack+1197, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2148, 1.0, int_stack+2289, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1134,int_stack+1290,int_stack+1272, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2253, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1188,int_stack+1320,int_stack+1290, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2289, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1278,int_stack+1383,int_stack+1365, 0.0, zero_stack, 1.0, int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2394,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6582,int_stack+1413,int_stack+1383, 0.0, zero_stack, 1.0, int_stack+1866, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2430,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1332,int_stack+1476,int_stack+1458, 0.0, zero_stack, 1.0, int_stack+1971, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2394, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1386,int_stack+1506,int_stack+1476, 0.0, zero_stack, 1.0, int_stack+2007, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2430, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1476,int_stack+1569,int_stack+1551, 0.0, zero_stack, 1.0, int_stack+2112, 0.0, zero_stack, 1.0, int_stack+2394, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6672,int_stack+1599,int_stack+1569, 0.0, zero_stack, 1.0, int_stack+2148, 0.0, zero_stack, 1.0, int_stack+2430, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1530,int_stack+1662,int_stack+1644, 0.0, zero_stack, 1.0, int_stack+2253, 1.0, int_stack+2394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6762,int_stack+1692,int_stack+1662, 0.0, zero_stack, 1.0, int_stack+2289, 1.0, int_stack+2430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1584,int_stack+1755,int_stack+1737, 0.0, zero_stack, 2.0, int_stack+2394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1638,int_stack+1785,int_stack+1755, 0.0, zero_stack, 2.0, int_stack+2430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1728,int_stack+1896,int_stack+1848, 1.0, int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2553,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6852,int_stack+1926,int_stack+1896, 1.0, int_stack+1866, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2589,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1782,int_stack+2037,int_stack+1989, 1.0, int_stack+1971, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2553, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1836,int_stack+2067,int_stack+2037, 1.0, int_stack+2007, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2589, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1926,int_stack+2178,int_stack+2130, 1.0, int_stack+2112, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2553, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1980,int_stack+2208,int_stack+2178, 1.0, int_stack+2148, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2589, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2070,int_stack+2319,int_stack+2271, 1.0, int_stack+2253, 0.0, zero_stack, 1.0, int_stack+2553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2124,int_stack+2349,int_stack+2319, 1.0, int_stack+2289, 0.0, zero_stack, 1.0, int_stack+2589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2214,int_stack+2460,int_stack+2412, 1.0, int_stack+2394, 1.0, int_stack+2553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2268,int_stack+2490,int_stack+2460, 1.0, int_stack+2430, 1.0, int_stack+2589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2358,int_stack+2619,int_stack+2571, 2.0, int_stack+2553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2412,int_stack+2649,int_stack+2619, 2.0, int_stack+2589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2502,int_stack+2712,int_stack+2694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3159,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2556,int_stack+2742,int_stack+2712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3195,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2646,int_stack+2805,int_stack+2787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3159, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+2835,int_stack+2805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3195, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2790,int_stack+2898,int_stack+2880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3159, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6942,int_stack+2928,int_stack+2898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3195, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2844,int_stack+2991,int_stack+2973, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3159, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2898,int_stack+3021,int_stack+2991, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2988,int_stack+3084,int_stack+3066, 0.0, zero_stack, 1.0, int_stack+3159, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7032,int_stack+3114,int_stack+3084, 0.0, zero_stack, 1.0, int_stack+3195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3042,int_stack+3225,int_stack+3177, 1.0, int_stack+3159, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3096,int_stack+3255,int_stack+3225, 1.0, int_stack+3195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+3186,int_stack+3318,int_stack+3300,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7122,int_stack+3348,int_stack+3318,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3240,int_stack+3411,int_stack+3393, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3858,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3294,int_stack+3441,int_stack+3411, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3894,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3384,int_stack+3504,int_stack+3486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3858, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7212,int_stack+3534,int_stack+3504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3894, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3438,int_stack+3597,int_stack+3579, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3858, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3492,int_stack+3627,int_stack+3597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3894, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3582,int_stack+3690,int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7302,int_stack+3720,int_stack+3690, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3636,int_stack+3783,int_stack+3765, 0.0, zero_stack, 1.0, int_stack+3858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3690,int_stack+3813,int_stack+3783, 0.0, zero_stack, 1.0, int_stack+3894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3780,int_stack+3924,int_stack+3876, 1.0, int_stack+3858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7392,int_stack+3954,int_stack+3924, 1.0, int_stack+3894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+3834,int_stack+4017,int_stack+3999,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3888,int_stack+4047,int_stack+4017,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+3978,int_stack+4110,int_stack+4092,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7482,int_stack+4140,int_stack+4110,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4032,int_stack+4203,int_stack+4185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4650,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4086,int_stack+4233,int_stack+4203, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4686,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4176,int_stack+4296,int_stack+4278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4650, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7572,int_stack+4326,int_stack+4296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4686, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4230,int_stack+4389,int_stack+4371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4650, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4284,int_stack+4419,int_stack+4389, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4686, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4374,int_stack+4482,int_stack+4464, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7662,int_stack+4512,int_stack+4482, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4686, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4428,int_stack+4575,int_stack+4557, 0.0, zero_stack, 1.0, int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4482,int_stack+4605,int_stack+4575, 0.0, zero_stack, 1.0, int_stack+4686, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4572,int_stack+4716,int_stack+4668, 1.0, int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7752,int_stack+4746,int_stack+4716, 1.0, int_stack+4686, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4626,int_stack+4809,int_stack+4791,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4680,int_stack+4839,int_stack+4809,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4770,int_stack+4902,int_stack+4884,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7842,int_stack+4932,int_stack+4902,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4824,int_stack+4995,int_stack+4977,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4878,int_stack+5025,int_stack+4995,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7932,int_stack+5178,int_stack+5124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5070,3);
     Libderiv->ABCD[11] = int_stack + 7932;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8040,int_stack+5322,int_stack+5268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5070, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 8040;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8148,int_stack+5412,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5070, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 8148;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5322,int_stack+5502,int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 5322;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5430,int_stack+5592,int_stack+108, 0.0, zero_stack, 1.0, int_stack+5070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 5430;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5538,int_stack+5682,int_stack+162, 1.0, int_stack+5070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 5538;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+5646,int_stack+5772,int_stack+216,3);
     Libderiv->ABCD[2] = int_stack + 5646;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+5754,int_stack+5862,int_stack+270,3);
     Libderiv->ABCD[1] = int_stack + 5754;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+4968,int_stack+5952,int_stack+324,3);
     Libderiv->ABCD[0] = int_stack + 4968;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5862,int_stack+6042,int_stack+378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5124,3);
     Libderiv->ABCD[155] = int_stack + 5862;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5970,int_stack+6132,int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5124, 1.0, int_stack+5268,3);
     Libderiv->ABCD[143] = int_stack + 5970;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+378,int_stack+540,int_stack+486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5268, 0.0, zero_stack,3);
     Libderiv->ABCD[142] = int_stack + 378;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+486,int_stack+6222,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5124, 0.0, zero_stack, 1.0, int_stack+0,3);
     Libderiv->ABCD[131] = int_stack + 486;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6078,int_stack+6312,int_stack+684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5268, 1.0, int_stack+0, 0.0, zero_stack,3);
     Libderiv->ABCD[130] = int_stack + 6078;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+594,int_stack+792,int_stack+738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[129] = int_stack + 594;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+702,int_stack+6402,int_stack+882, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5124, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54,3);
     Libderiv->ABCD[119] = int_stack + 702;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+810,int_stack+990,int_stack+936, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5268, 0.0, zero_stack, 1.0, int_stack+54, 0.0, zero_stack,3);
     Libderiv->ABCD[118] = int_stack + 810;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+918,int_stack+6492,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[117] = int_stack + 918;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1026,int_stack+1188,int_stack+1134, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[116] = int_stack + 1026;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1134,int_stack+6582,int_stack+1278, 0.0, zero_stack, 1.0, int_stack+5124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108,3);
     Libderiv->ABCD[107] = int_stack + 1134;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6186,int_stack+1386,int_stack+1332, 0.0, zero_stack, 1.0, int_stack+5268, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack,3);
     Libderiv->ABCD[106] = int_stack + 6186;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1242,int_stack+6672,int_stack+1476, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[105] = int_stack + 1242;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1350,int_stack+6762,int_stack+1530, 0.0, zero_stack, 1.0, int_stack+54, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[104] = int_stack + 1350;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1458,int_stack+1638,int_stack+1584, 0.0, zero_stack, 2.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[103] = int_stack + 1458;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1566,int_stack+6852,int_stack+1728, 1.0, int_stack+5124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162,3);
     Libderiv->ABCD[95] = int_stack + 1566;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1674,int_stack+1836,int_stack+1782, 1.0, int_stack+5268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162, 0.0, zero_stack,3);
     Libderiv->ABCD[94] = int_stack + 1674;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1782,int_stack+1980,int_stack+1926, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[93] = int_stack + 1782;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1890,int_stack+2124,int_stack+2070, 1.0, int_stack+54, 0.0, zero_stack, 1.0, int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[92] = int_stack + 1890;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+2268,int_stack+2214, 1.0, int_stack+108, 1.0, int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[91] = int_stack + 0;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1998,int_stack+2412,int_stack+2358, 2.0, int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[90] = int_stack + 1998;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+108,int_stack+2556,int_stack+2502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216,3);
     Libderiv->ABCD[47] = int_stack + 108;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2106,int_stack+2700,int_stack+2646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack,3);
     Libderiv->ABCD[46] = int_stack + 2106;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2214,int_stack+6942,int_stack+2790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[45] = int_stack + 2214;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2322,int_stack+2898,int_stack+2844, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[44] = int_stack + 2322;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2430,int_stack+7032,int_stack+2988, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[43] = int_stack + 2430;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2538,int_stack+3096,int_stack+3042, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[42] = int_stack + 2538;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2646,int_stack+7122,int_stack+3186,3);
     Libderiv->ABCD[38] = int_stack + 2646;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2754,int_stack+3294,int_stack+3240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270,3);
     Libderiv->ABCD[35] = int_stack + 2754;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2862,int_stack+7212,int_stack+3384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack,3);
     Libderiv->ABCD[34] = int_stack + 2862;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2970,int_stack+3492,int_stack+3438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[33] = int_stack + 2970;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3078,int_stack+7302,int_stack+3582, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[32] = int_stack + 3078;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3186,int_stack+3690,int_stack+3636, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[31] = int_stack + 3186;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3294,int_stack+7392,int_stack+3780, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[30] = int_stack + 3294;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+216,int_stack+3888,int_stack+3834,3);
     Libderiv->ABCD[26] = int_stack + 216;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+3402,int_stack+7482,int_stack+3978,3);
     Libderiv->ABCD[25] = int_stack + 3402;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3510,int_stack+4086,int_stack+4032, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324,3);
     Libderiv->ABCD[23] = int_stack + 3510;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3618,int_stack+7572,int_stack+4176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack,3);
     Libderiv->ABCD[22] = int_stack + 3618;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3726,int_stack+4284,int_stack+4230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[21] = int_stack + 3726;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3834,int_stack+7662,int_stack+4374, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[20] = int_stack + 3834;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3942,int_stack+4482,int_stack+4428, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[19] = int_stack + 3942;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4050,int_stack+7752,int_stack+4572, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[18] = int_stack + 4050;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+4158,int_stack+4680,int_stack+4626,3);
     Libderiv->ABCD[14] = int_stack + 4158;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+4266,int_stack+7842,int_stack+4770,3);
     Libderiv->ABCD[13] = int_stack + 4266;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+4374,int_stack+4878,int_stack+4824,3);
     Libderiv->ABCD[12] = int_stack + 4374;

}
