#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ddfp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|fp) integrals */

void d1hrr_order_ddfp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[3][3][11] = int_stack + 150;
 Libderiv->deriv_classes[3][4][11] = int_stack + 250;
 Libderiv->deriv_classes[4][3][11] = int_stack + 400;
 Libderiv->deriv_classes[4][4][11] = int_stack + 550;
 Libderiv->deriv_classes[2][3][10] = int_stack + 775;
 Libderiv->deriv_classes[2][4][10] = int_stack + 835;
 Libderiv->deriv_classes[3][3][10] = int_stack + 925;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1025;
 Libderiv->deriv_classes[4][3][10] = int_stack + 1175;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1325;
 Libderiv->deriv_classes[2][3][9] = int_stack + 1550;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1610;
 Libderiv->deriv_classes[3][3][9] = int_stack + 1700;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1800;
 Libderiv->deriv_classes[4][3][9] = int_stack + 1950;
 Libderiv->deriv_classes[4][4][9] = int_stack + 2100;
 Libderiv->deriv_classes[2][3][8] = int_stack + 2325;
 Libderiv->deriv_classes[2][4][8] = int_stack + 2385;
 Libderiv->deriv_classes[3][3][8] = int_stack + 2475;
 Libderiv->deriv_classes[3][4][8] = int_stack + 2575;
 Libderiv->deriv_classes[4][3][8] = int_stack + 2725;
 Libderiv->deriv_classes[4][4][8] = int_stack + 2875;
 Libderiv->deriv_classes[2][3][7] = int_stack + 3100;
 Libderiv->deriv_classes[2][4][7] = int_stack + 3160;
 Libderiv->deriv_classes[3][3][7] = int_stack + 3250;
 Libderiv->deriv_classes[3][4][7] = int_stack + 3350;
 Libderiv->deriv_classes[4][3][7] = int_stack + 3500;
 Libderiv->deriv_classes[4][4][7] = int_stack + 3650;
 Libderiv->deriv_classes[2][3][6] = int_stack + 3875;
 Libderiv->deriv_classes[2][4][6] = int_stack + 3935;
 Libderiv->deriv_classes[3][3][6] = int_stack + 4025;
 Libderiv->deriv_classes[3][4][6] = int_stack + 4125;
 Libderiv->dvrr_classes[4][3] = int_stack + 4275;
 Libderiv->deriv_classes[4][3][6] = int_stack + 4425;
 Libderiv->deriv_classes[4][4][6] = int_stack + 4575;
 Libderiv->deriv_classes[2][3][2] = int_stack + 4800;
 Libderiv->deriv_classes[2][4][2] = int_stack + 4860;
 Libderiv->deriv_classes[3][3][2] = int_stack + 4950;
 Libderiv->deriv_classes[3][4][2] = int_stack + 5050;
 Libderiv->deriv_classes[4][3][2] = int_stack + 5200;
 Libderiv->deriv_classes[4][4][2] = int_stack + 5350;
 Libderiv->deriv_classes[2][3][1] = int_stack + 5575;
 Libderiv->deriv_classes[2][4][1] = int_stack + 5635;
 Libderiv->deriv_classes[3][3][1] = int_stack + 5725;
 Libderiv->deriv_classes[3][4][1] = int_stack + 5825;
 Libderiv->deriv_classes[4][3][1] = int_stack + 5975;
 Libderiv->deriv_classes[4][4][1] = int_stack + 6125;
 Libderiv->dvrr_classes[2][3] = int_stack + 6350;
 Libderiv->dvrr_classes[2][4] = int_stack + 6410;
 Libderiv->deriv_classes[2][3][0] = int_stack + 6500;
 Libderiv->deriv_classes[2][4][0] = int_stack + 6560;
 Libderiv->dvrr_classes[3][3] = int_stack + 6650;
 Libderiv->dvrr_classes[3][4] = int_stack + 6750;
 Libderiv->deriv_classes[3][3][0] = int_stack + 6900;
 Libderiv->deriv_classes[3][4][0] = int_stack + 7000;
 Libderiv->deriv_classes[4][3][0] = int_stack + 7150;
 Libderiv->deriv_classes[4][4][0] = int_stack + 7300;
 memset(int_stack,0,60200);

 Libderiv->dvrr_stack = int_stack + 16465;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ddfp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6350,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7705,int_stack+250,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6650,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8005,int_stack+7705,int_stack+7525,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+550,int_stack+400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4275,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8995,int_stack+8545,int_stack+7705,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+835,int_stack+775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6350, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+1025,int_stack+925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6650, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+7525,int_stack+8545,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+1325,int_stack+1175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4275, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+540,int_stack+8545,int_stack+7525,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+1610,int_stack+1550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6350, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7705,int_stack+1800,int_stack+1700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6650, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9895,int_stack+7705,int_stack+7525,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+2100,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4275, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+10435,int_stack+8545,int_stack+7705,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+2385,int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+2575,int_stack+2475, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1440,int_stack+7525,int_stack+8545,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+2875,int_stack+2725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1980,int_stack+8545,int_stack+7525,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+3160,int_stack+3100, 0.0, zero_stack, 1.0, int_stack+6350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7705,int_stack+3350,int_stack+3250, 0.0, zero_stack, 1.0, int_stack+6650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2880,int_stack+7705,int_stack+7525,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+3650,int_stack+3500, 0.0, zero_stack, 1.0, int_stack+4275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11335,int_stack+8545,int_stack+7705,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+3935,int_stack+3875, 1.0, int_stack+6350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+4125,int_stack+4025, 1.0, int_stack+6650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3420,int_stack+7525,int_stack+8545,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+4575,int_stack+4425, 1.0, int_stack+4275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12235,int_stack+8545,int_stack+7525,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7525,int_stack+6410,int_stack+6350,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7705,int_stack+6750,int_stack+6650,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3960,int_stack+7705,int_stack+7525,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6650,int_stack+4860,int_stack+4800,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+5050,int_stack+4950,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+4500,int_stack+8545,int_stack+6650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13135,int_stack+5350,int_stack+5200,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+13585,int_stack+13135,int_stack+8545, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+5635,int_stack+5575,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13135,int_stack+5825,int_stack+5725,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5040,int_stack+13135,int_stack+8545, 0.0, zero_stack, 1.0, int_stack+7525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+6125,int_stack+5975,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+5580,int_stack+8545,int_stack+13135, 0.0, zero_stack, 1.0, int_stack+7705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13135,int_stack+6560,int_stack+6500,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8545,int_stack+7000,int_stack+6900,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6480,int_stack+8545,int_stack+13135, 1.0, int_stack+7525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13135,int_stack+7300,int_stack+7150,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+14485,int_stack+13135,int_stack+8545, 1.0, int_stack+7705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+15385,int_stack+8995,int_stack+8005,30);
     Libderiv->ABCD[11] = int_stack + 15385;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+7020,int_stack+540,int_stack+0,30);
     Libderiv->ABCD[10] = int_stack + 7020;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+10435,int_stack+9895,30);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+8100,int_stack+1980,int_stack+1440,30);
     Libderiv->ABCD[8] = int_stack + 8100;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+1080,int_stack+11335,int_stack+2880,30);
     Libderiv->ABCD[7] = int_stack + 1080;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+2160,int_stack+12235,int_stack+3420,30);
     Libderiv->ABCD[6] = int_stack + 2160;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+9180,int_stack+13585,int_stack+4500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[2] = int_stack + 9180;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+10260,int_stack+5580,int_stack+5040, 0.0, zero_stack, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[1] = int_stack + 10260;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+4500,int_stack+14485,int_stack+6480, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[0] = int_stack + 4500;

}
