#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpdd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|dd) integrals */

void d1hrr_order_fpdd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][2][11] = int_stack + 0;
 Libderiv->deriv_classes[3][3][11] = int_stack + 60;
 Libderiv->deriv_classes[3][4][11] = int_stack + 160;
 Libderiv->deriv_classes[4][2][11] = int_stack + 310;
 Libderiv->deriv_classes[4][3][11] = int_stack + 400;
 Libderiv->deriv_classes[4][4][11] = int_stack + 550;
 Libderiv->deriv_classes[3][2][10] = int_stack + 775;
 Libderiv->deriv_classes[3][3][10] = int_stack + 835;
 Libderiv->deriv_classes[3][4][10] = int_stack + 935;
 Libderiv->deriv_classes[4][2][10] = int_stack + 1085;
 Libderiv->deriv_classes[4][3][10] = int_stack + 1175;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1325;
 Libderiv->deriv_classes[3][2][9] = int_stack + 1550;
 Libderiv->deriv_classes[3][3][9] = int_stack + 1610;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1710;
 Libderiv->deriv_classes[4][2][9] = int_stack + 1860;
 Libderiv->deriv_classes[4][3][9] = int_stack + 1950;
 Libderiv->deriv_classes[4][4][9] = int_stack + 2100;
 Libderiv->deriv_classes[3][2][8] = int_stack + 2325;
 Libderiv->deriv_classes[3][3][8] = int_stack + 2385;
 Libderiv->deriv_classes[3][4][8] = int_stack + 2485;
 Libderiv->deriv_classes[4][2][8] = int_stack + 2635;
 Libderiv->deriv_classes[4][3][8] = int_stack + 2725;
 Libderiv->deriv_classes[4][4][8] = int_stack + 2875;
 Libderiv->deriv_classes[3][2][7] = int_stack + 3100;
 Libderiv->deriv_classes[3][3][7] = int_stack + 3160;
 Libderiv->deriv_classes[3][4][7] = int_stack + 3260;
 Libderiv->deriv_classes[4][2][7] = int_stack + 3410;
 Libderiv->deriv_classes[4][3][7] = int_stack + 3500;
 Libderiv->deriv_classes[4][4][7] = int_stack + 3650;
 Libderiv->deriv_classes[3][2][6] = int_stack + 3875;
 Libderiv->deriv_classes[3][3][6] = int_stack + 3935;
 Libderiv->deriv_classes[3][4][6] = int_stack + 4035;
 Libderiv->dvrr_classes[4][2] = int_stack + 4185;
 Libderiv->deriv_classes[4][2][6] = int_stack + 4275;
 Libderiv->dvrr_classes[4][3] = int_stack + 4365;
 Libderiv->deriv_classes[4][3][6] = int_stack + 4515;
 Libderiv->deriv_classes[4][4][6] = int_stack + 4665;
 Libderiv->deriv_classes[3][2][2] = int_stack + 4890;
 Libderiv->deriv_classes[3][3][2] = int_stack + 4950;
 Libderiv->deriv_classes[3][4][2] = int_stack + 5050;
 Libderiv->deriv_classes[4][2][2] = int_stack + 5200;
 Libderiv->deriv_classes[4][3][2] = int_stack + 5290;
 Libderiv->deriv_classes[4][4][2] = int_stack + 5440;
 Libderiv->deriv_classes[3][2][1] = int_stack + 5665;
 Libderiv->deriv_classes[3][3][1] = int_stack + 5725;
 Libderiv->deriv_classes[3][4][1] = int_stack + 5825;
 Libderiv->deriv_classes[4][2][1] = int_stack + 5975;
 Libderiv->deriv_classes[4][3][1] = int_stack + 6065;
 Libderiv->deriv_classes[4][4][1] = int_stack + 6215;
 Libderiv->dvrr_classes[3][2] = int_stack + 6440;
 Libderiv->dvrr_classes[3][3] = int_stack + 6500;
 Libderiv->dvrr_classes[3][4] = int_stack + 6600;
 Libderiv->deriv_classes[3][2][0] = int_stack + 6750;
 Libderiv->deriv_classes[3][3][0] = int_stack + 6810;
 Libderiv->deriv_classes[3][4][0] = int_stack + 6910;
 Libderiv->deriv_classes[4][2][0] = int_stack + 7060;
 Libderiv->deriv_classes[4][3][0] = int_stack + 7150;
 Libderiv->deriv_classes[4][4][0] = int_stack + 7300;
 memset(int_stack,0,60200);

 Libderiv->dvrr_stack = int_stack + 13315;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpdd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7525,int_stack+6500,int_stack+6440,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7705,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6440,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7885,int_stack+160,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6500,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8185,int_stack+7885,int_stack+7705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7525,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7705,int_stack+4365,int_stack+4185,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+400,int_stack+310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4185,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+550,int_stack+400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4365,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8995,int_stack+8545,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7705,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+835,int_stack+775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6440, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+935,int_stack+835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6500, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+480,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7525, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+1175,int_stack+1085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4185, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+1325,int_stack+1175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4365, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+840,int_stack+8545,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7705, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+1610,int_stack+1550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6440, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+1710,int_stack+1610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8545,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7525, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+1950,int_stack+1860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4185, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1380,int_stack+2100,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4365, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9535,int_stack+1380,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7705, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+2385,int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+2485,int_stack+2385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1380,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+2725,int_stack+2635, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1740,int_stack+2875,int_stack+2725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2190,int_stack+1740,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+3160,int_stack+3100, 0.0, zero_stack, 1.0, int_stack+6440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+3260,int_stack+3160, 0.0, zero_stack, 1.0, int_stack+6500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1740,int_stack+180,int_stack+0, 0.0, zero_stack, 1.0, int_stack+7525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+3500,int_stack+3410, 0.0, zero_stack, 1.0, int_stack+4185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2730,int_stack+3650,int_stack+3500, 0.0, zero_stack, 1.0, int_stack+4365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3180,int_stack+2730,int_stack+0, 0.0, zero_stack, 1.0, int_stack+7705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+3935,int_stack+3875, 1.0, int_stack+6440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+4035,int_stack+3935, 1.0, int_stack+6500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2730,int_stack+180,int_stack+0, 1.0, int_stack+7525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+4515,int_stack+4275, 1.0, int_stack+4185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3720,int_stack+4665,int_stack+4515, 1.0, int_stack+4365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4170,int_stack+3720,int_stack+0, 1.0, int_stack+7705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7705,int_stack+6600,int_stack+6500,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+7705,int_stack+7525,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4710,int_stack+4950,int_stack+4890,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+5050,int_stack+4950,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+7825,int_stack+7525,int_stack+4710,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4710,int_stack+5290,int_stack+5200,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3720,int_stack+5440,int_stack+5290,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+4980,int_stack+3720,int_stack+4710,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4710,int_stack+5725,int_stack+5665,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+5825,int_stack+5725,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+3720,int_stack+7525,int_stack+4710,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4710,int_stack+6065,int_stack+5975,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5520,int_stack+6215,int_stack+6065,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+5970,int_stack+5520,int_stack+4710,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4710,int_stack+6810,int_stack+6750,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+6910,int_stack+6810,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+5520,int_stack+7525,int_stack+4710,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4710,int_stack+7150,int_stack+7060,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6510,int_stack+7300,int_stack+7150,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6960,int_stack+6510,int_stack+4710,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+10075,int_stack+8995,int_stack+8185,36);
     Libderiv->ABCD[11] = int_stack + 10075;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11155,int_stack+840,int_stack+480,36);
     Libderiv->ABCD[10] = int_stack + 11155;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12235,int_stack+9535,int_stack+8545,36);
     Libderiv->ABCD[9] = int_stack + 12235;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8185,int_stack+2190,int_stack+1380,36);
     Libderiv->ABCD[8] = int_stack + 8185;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+360,int_stack+3180,int_stack+1740,36);
     Libderiv->ABCD[7] = int_stack + 360;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1440,int_stack+4170,int_stack+2730,36);
     Libderiv->ABCD[6] = int_stack + 1440;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+2520,int_stack+4980,int_stack+7825, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[2] = int_stack + 2520;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+4080,int_stack+5970,int_stack+3720, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[1] = int_stack + 4080;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+5880,int_stack+6960,int_stack+5520, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[0] = int_stack + 5880;

}
