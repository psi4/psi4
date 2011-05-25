#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpgd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|gd) integrals */

void d1hrr_order_dpgd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[2][6][11] = int_stack + 216;
 Libderiv->deriv_classes[3][4][11] = int_stack + 384;
 Libderiv->deriv_classes[3][5][11] = int_stack + 534;
 Libderiv->deriv_classes[3][6][11] = int_stack + 744;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1024;
 Libderiv->deriv_classes[2][5][10] = int_stack + 1114;
 Libderiv->deriv_classes[2][6][10] = int_stack + 1240;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1408;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1558;
 Libderiv->deriv_classes[3][6][10] = int_stack + 1768;
 Libderiv->deriv_classes[2][4][9] = int_stack + 2048;
 Libderiv->deriv_classes[2][5][9] = int_stack + 2138;
 Libderiv->deriv_classes[2][6][9] = int_stack + 2264;
 Libderiv->deriv_classes[3][4][9] = int_stack + 2432;
 Libderiv->deriv_classes[3][5][9] = int_stack + 2582;
 Libderiv->deriv_classes[3][6][9] = int_stack + 2792;
 Libderiv->deriv_classes[2][4][8] = int_stack + 3072;
 Libderiv->deriv_classes[2][5][8] = int_stack + 3162;
 Libderiv->deriv_classes[2][6][8] = int_stack + 3288;
 Libderiv->deriv_classes[3][4][8] = int_stack + 3456;
 Libderiv->deriv_classes[3][5][8] = int_stack + 3606;
 Libderiv->deriv_classes[3][6][8] = int_stack + 3816;
 Libderiv->deriv_classes[2][4][7] = int_stack + 4096;
 Libderiv->deriv_classes[2][5][7] = int_stack + 4186;
 Libderiv->deriv_classes[2][6][7] = int_stack + 4312;
 Libderiv->deriv_classes[3][4][7] = int_stack + 4480;
 Libderiv->deriv_classes[3][5][7] = int_stack + 4630;
 Libderiv->deriv_classes[3][6][7] = int_stack + 4840;
 Libderiv->deriv_classes[2][4][6] = int_stack + 5120;
 Libderiv->deriv_classes[2][5][6] = int_stack + 5210;
 Libderiv->deriv_classes[2][6][6] = int_stack + 5336;
 Libderiv->dvrr_classes[3][4] = int_stack + 5504;
 Libderiv->deriv_classes[3][4][6] = int_stack + 5654;
 Libderiv->dvrr_classes[3][5] = int_stack + 5804;
 Libderiv->deriv_classes[3][5][6] = int_stack + 6014;
 Libderiv->deriv_classes[3][6][6] = int_stack + 6224;
 Libderiv->deriv_classes[2][4][2] = int_stack + 6504;
 Libderiv->deriv_classes[2][5][2] = int_stack + 6594;
 Libderiv->deriv_classes[2][6][2] = int_stack + 6720;
 Libderiv->deriv_classes[3][4][2] = int_stack + 6888;
 Libderiv->deriv_classes[3][5][2] = int_stack + 7038;
 Libderiv->deriv_classes[3][6][2] = int_stack + 7248;
 Libderiv->deriv_classes[2][4][1] = int_stack + 7528;
 Libderiv->deriv_classes[2][5][1] = int_stack + 7618;
 Libderiv->deriv_classes[2][6][1] = int_stack + 7744;
 Libderiv->deriv_classes[3][4][1] = int_stack + 7912;
 Libderiv->deriv_classes[3][5][1] = int_stack + 8062;
 Libderiv->deriv_classes[3][6][1] = int_stack + 8272;
 Libderiv->dvrr_classes[2][4] = int_stack + 8552;
 Libderiv->dvrr_classes[2][5] = int_stack + 8642;
 Libderiv->dvrr_classes[2][6] = int_stack + 8768;
 Libderiv->deriv_classes[2][4][0] = int_stack + 8936;
 Libderiv->deriv_classes[2][5][0] = int_stack + 9026;
 Libderiv->deriv_classes[2][6][0] = int_stack + 9152;
 Libderiv->deriv_classes[3][4][0] = int_stack + 9320;
 Libderiv->deriv_classes[3][5][0] = int_stack + 9470;
 Libderiv->deriv_classes[3][6][0] = int_stack + 9680;
 memset(int_stack,0,79680);

 Libderiv->dvrr_stack = int_stack + 21066;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpgd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9960,int_stack+8642,int_stack+8552,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10230,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8552,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10500,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8642,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10878,int_stack+10500,int_stack+10230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9960,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10230,int_stack+5804,int_stack+5504,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+534,int_stack+384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5504,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11868,int_stack+744,int_stack+534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5804,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+11868,int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10230,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+1114,int_stack+1024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8552, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11688,int_stack+1240,int_stack+1114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8642, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+12066,int_stack+11688,int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9960, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+1558,int_stack+1408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5504, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+900,int_stack+1768,int_stack+1558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5804, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+12606,int_stack+900,int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10230, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+2138,int_stack+2048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8552, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11688,int_stack+2264,int_stack+2138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8642, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+900,int_stack+11688,int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9960, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+2582,int_stack+2432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5504, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1440,int_stack+2792,int_stack+2582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5804, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2070,int_stack+1440,int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10230, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+3162,int_stack+3072, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11688,int_stack+3288,int_stack+3162, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1440,int_stack+11688,int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+3606,int_stack+3456, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2970,int_stack+3816,int_stack+3606, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5804, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13506,int_stack+2970,int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+4186,int_stack+4096, 0.0, zero_stack, 1.0, int_stack+8552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11688,int_stack+4312,int_stack+4186, 0.0, zero_stack, 1.0, int_stack+8642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2970,int_stack+11688,int_stack+11418, 0.0, zero_stack, 1.0, int_stack+9960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+4630,int_stack+4480, 0.0, zero_stack, 1.0, int_stack+5504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3510,int_stack+4840,int_stack+4630, 0.0, zero_stack, 1.0, int_stack+5804, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4140,int_stack+3510,int_stack+11418, 0.0, zero_stack, 1.0, int_stack+10230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+5210,int_stack+5120, 1.0, int_stack+8552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11688,int_stack+5336,int_stack+5210, 1.0, int_stack+8642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3510,int_stack+11688,int_stack+11418, 1.0, int_stack+9960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11418,int_stack+6014,int_stack+5654, 1.0, int_stack+5504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5040,int_stack+6224,int_stack+6014, 1.0, int_stack+5804, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14406,int_stack+5040,int_stack+11418, 1.0, int_stack+10230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10230,int_stack+8768,int_stack+8642,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11418,int_stack+10230,int_stack+9960,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9960,int_stack+6594,int_stack+6504,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10230,int_stack+6720,int_stack+6594,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5040,int_stack+10230,int_stack+9960,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9960,int_stack+7038,int_stack+6888,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5580,int_stack+7248,int_stack+7038,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6210,int_stack+5580,int_stack+9960,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9960,int_stack+7618,int_stack+7528,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10230,int_stack+7744,int_stack+7618,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5580,int_stack+10230,int_stack+9960,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9960,int_stack+8062,int_stack+7912,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7110,int_stack+8272,int_stack+8062,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+7740,int_stack+7110,int_stack+9960,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9960,int_stack+9026,int_stack+8936,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10230,int_stack+9152,int_stack+9026,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+7110,int_stack+10230,int_stack+9960,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9960,int_stack+9470,int_stack+9320,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8640,int_stack+9680,int_stack+9470,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15306,int_stack+8640,int_stack+9960,10);
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8640,int_stack+0,int_stack+10878,90);
     Libderiv->ABCD[11] = int_stack + 8640;
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+16206,int_stack+12606,int_stack+12066,90);
     Libderiv->ABCD[10] = int_stack + 16206;
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+17826,int_stack+2070,int_stack+900,90);
     Libderiv->ABCD[9] = int_stack + 17826;
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+19446,int_stack+13506,int_stack+1440,90);
     Libderiv->ABCD[8] = int_stack + 19446;
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+4140,int_stack+2970,90);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1620,int_stack+14406,int_stack+3510,90);
     Libderiv->ABCD[6] = int_stack + 1620;
 /*--- compute (dp|gd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+3240,int_stack+6210,int_stack+5040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[2] = int_stack + 3240;
 /*--- compute (dp|gd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+11958,int_stack+7740,int_stack+5580, 0.0, zero_stack, 1.0, int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[1] = int_stack + 11958;
 /*--- compute (dp|gd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+4860,int_stack+15306,int_stack+7110, 1.0, int_stack+11418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[0] = int_stack + 4860;

}
