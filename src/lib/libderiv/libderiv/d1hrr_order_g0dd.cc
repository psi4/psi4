#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_g0dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (g0|dd) integrals */

void d1hrr_order_g0dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][2][11] = int_stack + 0;
 Libderiv->deriv_classes[4][3][11] = int_stack + 90;
 Libderiv->deriv_classes[4][4][11] = int_stack + 240;
 Libderiv->deriv_classes[4][2][10] = int_stack + 465;
 Libderiv->deriv_classes[4][3][10] = int_stack + 555;
 Libderiv->deriv_classes[4][4][10] = int_stack + 705;
 Libderiv->deriv_classes[4][2][9] = int_stack + 930;
 Libderiv->deriv_classes[4][3][9] = int_stack + 1020;
 Libderiv->deriv_classes[4][4][9] = int_stack + 1170;
 Libderiv->deriv_classes[4][2][8] = int_stack + 1395;
 Libderiv->deriv_classes[4][3][8] = int_stack + 1485;
 Libderiv->deriv_classes[4][4][8] = int_stack + 1635;
 Libderiv->deriv_classes[4][2][7] = int_stack + 1860;
 Libderiv->deriv_classes[4][3][7] = int_stack + 1950;
 Libderiv->deriv_classes[4][4][7] = int_stack + 2100;
 Libderiv->dvrr_classes[4][2] = int_stack + 2325;
 Libderiv->deriv_classes[4][2][6] = int_stack + 2415;
 Libderiv->dvrr_classes[4][3] = int_stack + 2505;
 Libderiv->deriv_classes[4][3][6] = int_stack + 2655;
 Libderiv->deriv_classes[4][4][6] = int_stack + 2805;
 Libderiv->deriv_classes[4][2][2] = int_stack + 3030;
 Libderiv->deriv_classes[4][3][2] = int_stack + 3120;
 Libderiv->deriv_classes[4][4][2] = int_stack + 3270;
 Libderiv->deriv_classes[4][2][1] = int_stack + 3495;
 Libderiv->deriv_classes[4][3][1] = int_stack + 3585;
 Libderiv->deriv_classes[4][4][1] = int_stack + 3735;
 Libderiv->deriv_classes[4][2][0] = int_stack + 3960;
 Libderiv->deriv_classes[4][3][0] = int_stack + 4050;
 Libderiv->deriv_classes[4][4][0] = int_stack + 4200;
 memset(int_stack,0,35400);

 Libderiv->dvrr_stack = int_stack + 7755;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_g0dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4425,int_stack+2505,int_stack+2325,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4695,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2325,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4965,int_stack+240,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2505,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+555,int_stack+465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2325, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5415,int_stack+705,int_stack+555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2505, 0.0, zero_stack,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+1020,int_stack+930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2325, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+1170,int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2505, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+990,int_stack+1485,int_stack+1395, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5865,int_stack+1635,int_stack+1485, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1260,int_stack+1950,int_stack+1860, 0.0, zero_stack, 1.0, int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6315,int_stack+2100,int_stack+1950, 0.0, zero_stack, 1.0, int_stack+2505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1530,int_stack+2655,int_stack+2415, 1.0, int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2805,int_stack+2655, 1.0, int_stack+2505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2250,int_stack+3120,int_stack+3030,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2520,int_stack+3270,int_stack+3120,15);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2970,int_stack+3585,int_stack+3495,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6765,int_stack+3735,int_stack+3585,15);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+3240,int_stack+4050,int_stack+3960,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3510,int_stack+4200,int_stack+4050,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7215,int_stack+4965,int_stack+4695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4425,15);
     Libderiv->ABCD[11] = int_stack + 7215;
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4695,int_stack+5415,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4425, 0.0, zero_stack,15);
     Libderiv->ABCD[10] = int_stack + 4695;
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5235,int_stack+540,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4425, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[9] = int_stack + 5235;
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+5865,int_stack+990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5775,int_stack+6315,int_stack+1260, 0.0, zero_stack, 1.0, int_stack+4425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[7] = int_stack + 5775;
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+540,int_stack+1800,int_stack+1530, 1.0, int_stack+4425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[6] = int_stack + 540;
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1080,int_stack+2520,int_stack+2250,15);
     Libderiv->ABCD[2] = int_stack + 1080;
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1620,int_stack+6765,int_stack+2970,15);
     Libderiv->ABCD[1] = int_stack + 1620;
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6315,int_stack+3510,int_stack+3240,15);
     Libderiv->ABCD[0] = int_stack + 6315;

}
