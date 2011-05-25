#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppdd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|dd) integrals */

void d1hrr_order_ppdd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][2][11] = int_stack + 0;
 Libderiv->deriv_classes[1][3][11] = int_stack + 18;
 Libderiv->deriv_classes[1][4][11] = int_stack + 48;
 Libderiv->deriv_classes[2][2][11] = int_stack + 93;
 Libderiv->deriv_classes[2][3][11] = int_stack + 129;
 Libderiv->deriv_classes[2][4][11] = int_stack + 189;
 Libderiv->deriv_classes[1][2][10] = int_stack + 279;
 Libderiv->deriv_classes[1][3][10] = int_stack + 297;
 Libderiv->deriv_classes[1][4][10] = int_stack + 327;
 Libderiv->deriv_classes[2][2][10] = int_stack + 372;
 Libderiv->deriv_classes[2][3][10] = int_stack + 408;
 Libderiv->deriv_classes[2][4][10] = int_stack + 468;
 Libderiv->deriv_classes[1][2][9] = int_stack + 558;
 Libderiv->deriv_classes[1][3][9] = int_stack + 576;
 Libderiv->deriv_classes[1][4][9] = int_stack + 606;
 Libderiv->deriv_classes[2][2][9] = int_stack + 651;
 Libderiv->deriv_classes[2][3][9] = int_stack + 687;
 Libderiv->deriv_classes[2][4][9] = int_stack + 747;
 Libderiv->deriv_classes[1][2][8] = int_stack + 837;
 Libderiv->deriv_classes[1][3][8] = int_stack + 855;
 Libderiv->deriv_classes[1][4][8] = int_stack + 885;
 Libderiv->deriv_classes[2][2][8] = int_stack + 930;
 Libderiv->deriv_classes[2][3][8] = int_stack + 966;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1026;
 Libderiv->deriv_classes[1][2][7] = int_stack + 1116;
 Libderiv->deriv_classes[1][3][7] = int_stack + 1134;
 Libderiv->deriv_classes[1][4][7] = int_stack + 1164;
 Libderiv->deriv_classes[2][2][7] = int_stack + 1209;
 Libderiv->deriv_classes[2][3][7] = int_stack + 1245;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1305;
 Libderiv->deriv_classes[1][2][6] = int_stack + 1395;
 Libderiv->deriv_classes[1][3][6] = int_stack + 1413;
 Libderiv->deriv_classes[1][4][6] = int_stack + 1443;
 Libderiv->dvrr_classes[2][2] = int_stack + 1488;
 Libderiv->deriv_classes[2][2][6] = int_stack + 1524;
 Libderiv->dvrr_classes[2][3] = int_stack + 1560;
 Libderiv->deriv_classes[2][3][6] = int_stack + 1620;
 Libderiv->deriv_classes[2][4][6] = int_stack + 1680;
 Libderiv->deriv_classes[1][2][2] = int_stack + 1770;
 Libderiv->deriv_classes[1][3][2] = int_stack + 1788;
 Libderiv->deriv_classes[1][4][2] = int_stack + 1818;
 Libderiv->deriv_classes[2][2][2] = int_stack + 1863;
 Libderiv->deriv_classes[2][3][2] = int_stack + 1899;
 Libderiv->deriv_classes[2][4][2] = int_stack + 1959;
 Libderiv->deriv_classes[1][2][1] = int_stack + 2049;
 Libderiv->deriv_classes[1][3][1] = int_stack + 2067;
 Libderiv->deriv_classes[1][4][1] = int_stack + 2097;
 Libderiv->deriv_classes[2][2][1] = int_stack + 2142;
 Libderiv->deriv_classes[2][3][1] = int_stack + 2178;
 Libderiv->deriv_classes[2][4][1] = int_stack + 2238;
 Libderiv->dvrr_classes[1][2] = int_stack + 2328;
 Libderiv->dvrr_classes[1][3] = int_stack + 2346;
 Libderiv->dvrr_classes[1][4] = int_stack + 2376;
 Libderiv->deriv_classes[1][2][0] = int_stack + 2421;
 Libderiv->deriv_classes[1][3][0] = int_stack + 2439;
 Libderiv->deriv_classes[1][4][0] = int_stack + 2469;
 Libderiv->deriv_classes[2][2][0] = int_stack + 2514;
 Libderiv->deriv_classes[2][3][0] = int_stack + 2550;
 Libderiv->deriv_classes[2][4][0] = int_stack + 2610;
 memset(int_stack,0,21600);

 Libderiv->dvrr_stack = int_stack + 4842;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppdd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+2346,int_stack+2328,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2754,int_stack+18,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2328,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2808,int_stack+48,int_stack+18, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2346,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2898,int_stack+2808,int_stack+2754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2754,int_stack+1560,int_stack+1488,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+129,int_stack+93, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1488,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3114,int_stack+189,int_stack+129, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+3114,int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2754,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+297,int_stack+279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2328, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3060,int_stack+327,int_stack+297, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2346, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3150,int_stack+3060,int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+408,int_stack+372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1488, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3258,int_stack+468,int_stack+408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3438,int_stack+3258,int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2754, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+576,int_stack+558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2328, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3060,int_stack+606,int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2346, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3258,int_stack+3060,int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+687,int_stack+651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1488, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3654,int_stack+747,int_stack+687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+216,int_stack+3654,int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2754, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+855,int_stack+837, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3060,int_stack+885,int_stack+855, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3654,int_stack+3060,int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+966,int_stack+930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3762,int_stack+1026,int_stack+966, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+432,int_stack+3762,int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+1134,int_stack+1116, 0.0, zero_stack, 1.0, int_stack+2328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3060,int_stack+1164,int_stack+1134, 0.0, zero_stack, 1.0, int_stack+2346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3762,int_stack+3060,int_stack+3006, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+1245,int_stack+1209, 0.0, zero_stack, 1.0, int_stack+1488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+648,int_stack+1305,int_stack+1245, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+828,int_stack+648,int_stack+3006, 0.0, zero_stack, 1.0, int_stack+2754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+1413,int_stack+1395, 1.0, int_stack+2328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3060,int_stack+1443,int_stack+1413, 1.0, int_stack+2346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+648,int_stack+3060,int_stack+3006, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3006,int_stack+1620,int_stack+1524, 1.0, int_stack+1488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1044,int_stack+1680,int_stack+1620, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1224,int_stack+1044,int_stack+3006, 1.0, int_stack+2754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2754,int_stack+2376,int_stack+2346,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+3006,int_stack+2754,int_stack+2700,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+1788,int_stack+1770,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2754,int_stack+1818,int_stack+1788,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1044,int_stack+2754,int_stack+2700,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+1899,int_stack+1863,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1440,int_stack+1959,int_stack+1899,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1620,int_stack+1440,int_stack+2700,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+2067,int_stack+2049,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2754,int_stack+2097,int_stack+2067,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1440,int_stack+2754,int_stack+2700,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+2178,int_stack+2142,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1836,int_stack+2238,int_stack+2178,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2016,int_stack+1836,int_stack+2700,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+2439,int_stack+2421,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2754,int_stack+2469,int_stack+2439,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1836,int_stack+2754,int_stack+2700,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+2550,int_stack+2514,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2232,int_stack+2610,int_stack+2550,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2412,int_stack+2232,int_stack+2700,6);
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3870,int_stack+0,int_stack+2898,36);
     Libderiv->ABCD[11] = int_stack + 3870;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4194,int_stack+3438,int_stack+3150,36);
     Libderiv->ABCD[10] = int_stack + 4194;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4518,int_stack+216,int_stack+3258,36);
     Libderiv->ABCD[9] = int_stack + 4518;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+432,int_stack+3654,36);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+324,int_stack+828,int_stack+3762,36);
     Libderiv->ABCD[7] = int_stack + 324;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3114,int_stack+1224,int_stack+648,36);
     Libderiv->ABCD[6] = int_stack + 3114;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+648,int_stack+1620,int_stack+1044, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[2] = int_stack + 648;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+972,int_stack+2016,int_stack+1440, 0.0, zero_stack, 1.0, int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[1] = int_stack + 972;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1296,int_stack+2412,int_stack+1836, 1.0, int_stack+3006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[0] = int_stack + 1296;

}
