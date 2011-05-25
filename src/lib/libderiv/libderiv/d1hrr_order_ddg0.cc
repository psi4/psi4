#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ddg0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|g0) integrals */

void d1hrr_order_ddg0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 90;
 Libderiv->deriv_classes[4][4][11] = int_stack + 240;
 Libderiv->deriv_classes[2][4][10] = int_stack + 465;
 Libderiv->deriv_classes[3][4][10] = int_stack + 555;
 Libderiv->deriv_classes[4][4][10] = int_stack + 705;
 Libderiv->deriv_classes[2][4][9] = int_stack + 930;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1020;
 Libderiv->deriv_classes[4][4][9] = int_stack + 1170;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1395;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1485;
 Libderiv->deriv_classes[4][4][8] = int_stack + 1635;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1860;
 Libderiv->deriv_classes[3][4][7] = int_stack + 1950;
 Libderiv->deriv_classes[4][4][7] = int_stack + 2100;
 Libderiv->deriv_classes[2][4][6] = int_stack + 2325;
 Libderiv->deriv_classes[3][4][6] = int_stack + 2415;
 Libderiv->deriv_classes[4][4][6] = int_stack + 2565;
 Libderiv->deriv_classes[2][4][2] = int_stack + 2790;
 Libderiv->deriv_classes[3][4][2] = int_stack + 2880;
 Libderiv->deriv_classes[4][4][2] = int_stack + 3030;
 Libderiv->deriv_classes[2][4][1] = int_stack + 3255;
 Libderiv->deriv_classes[3][4][1] = int_stack + 3345;
 Libderiv->deriv_classes[4][4][1] = int_stack + 3495;
 Libderiv->dvrr_classes[2][4] = int_stack + 3720;
 Libderiv->deriv_classes[2][4][0] = int_stack + 3810;
 Libderiv->dvrr_classes[3][4] = int_stack + 3900;
 Libderiv->deriv_classes[3][4][0] = int_stack + 4050;
 Libderiv->deriv_classes[4][4][0] = int_stack + 4200;
 memset(int_stack,0,35400);

 Libderiv->dvrr_stack = int_stack + 7395;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ddg0(Libderiv, Data);
   Data++;
 }

 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4425,int_stack+90,int_stack+0,15);
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+4695,int_stack+240,int_stack+90,15);
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+555,int_stack+465,15);
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5145,int_stack+705,int_stack+555,15);
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+270,int_stack+1020,int_stack+930,15);
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+540,int_stack+1170,int_stack+1020,15);
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+990,int_stack+1485,int_stack+1395,15);
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5595,int_stack+1635,int_stack+1485,15);
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1260,int_stack+1950,int_stack+1860,15);
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+6045,int_stack+2100,int_stack+1950,15);
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1530,int_stack+2415,int_stack+2325,15);
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1800,int_stack+2565,int_stack+2415,15);
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2250,int_stack+3900,int_stack+3720,15);
 /*--- compute (dp|g0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+2520,int_stack+2880,int_stack+2790, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|g0) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+6495,int_stack+3030,int_stack+2880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (dp|g0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+2790,int_stack+3345,int_stack+3255, 0.0, zero_stack, 1.0, int_stack+3720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|g0) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+6945,int_stack+3495,int_stack+3345, 0.0, zero_stack, 1.0, int_stack+3900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (dp|g0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+3060,int_stack+4050,int_stack+3810, 1.0, int_stack+3720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|g0) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+3330,int_stack+4200,int_stack+4050, 1.0, int_stack+3900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (dd|g0) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+3780,int_stack+4695,int_stack+4425,15);
     Libderiv->ABCD[11] = int_stack + 3780;
 /*--- compute (dd|g0) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+4320,int_stack+5145,int_stack+0,15);
     Libderiv->ABCD[10] = int_stack + 4320;
 /*--- compute (dd|g0) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+4860,int_stack+540,int_stack+270,15);
     Libderiv->ABCD[9] = int_stack + 4860;
 /*--- compute (dd|g0) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+5595,int_stack+990,15);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (dd|g0) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+540,int_stack+6045,int_stack+1260,15);
     Libderiv->ABCD[7] = int_stack + 540;
 /*--- compute (dd|g0) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+5400,int_stack+1800,int_stack+1530,15);
     Libderiv->ABCD[6] = int_stack + 5400;
 /*--- compute (dd|g0) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+1080,int_stack+6495,int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[2] = int_stack + 1080;
 /*--- compute (dd|g0) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+1620,int_stack+6945,int_stack+2790, 0.0, zero_stack, 1.0, int_stack+2250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[1] = int_stack + 1620;
 /*--- compute (dd|g0) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+2520,int_stack+3330,int_stack+3060, 1.0, int_stack+2250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[0] = int_stack + 2520;

}
