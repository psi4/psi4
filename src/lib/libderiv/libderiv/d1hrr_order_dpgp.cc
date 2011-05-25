#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpgp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|gp) integrals */

void d1hrr_order_dpgp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[3][4][11] = int_stack + 216;
 Libderiv->deriv_classes[3][5][11] = int_stack + 366;
 Libderiv->deriv_classes[2][4][10] = int_stack + 576;
 Libderiv->deriv_classes[2][5][10] = int_stack + 666;
 Libderiv->deriv_classes[3][4][10] = int_stack + 792;
 Libderiv->deriv_classes[3][5][10] = int_stack + 942;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1152;
 Libderiv->deriv_classes[2][5][9] = int_stack + 1242;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1368;
 Libderiv->deriv_classes[3][5][9] = int_stack + 1518;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1728;
 Libderiv->deriv_classes[2][5][8] = int_stack + 1818;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1944;
 Libderiv->deriv_classes[3][5][8] = int_stack + 2094;
 Libderiv->deriv_classes[2][4][7] = int_stack + 2304;
 Libderiv->deriv_classes[2][5][7] = int_stack + 2394;
 Libderiv->deriv_classes[3][4][7] = int_stack + 2520;
 Libderiv->deriv_classes[3][5][7] = int_stack + 2670;
 Libderiv->deriv_classes[2][4][6] = int_stack + 2880;
 Libderiv->deriv_classes[2][5][6] = int_stack + 2970;
 Libderiv->dvrr_classes[3][4] = int_stack + 3096;
 Libderiv->deriv_classes[3][4][6] = int_stack + 3246;
 Libderiv->deriv_classes[3][5][6] = int_stack + 3396;
 Libderiv->deriv_classes[2][4][2] = int_stack + 3606;
 Libderiv->deriv_classes[2][5][2] = int_stack + 3696;
 Libderiv->deriv_classes[3][4][2] = int_stack + 3822;
 Libderiv->deriv_classes[3][5][2] = int_stack + 3972;
 Libderiv->deriv_classes[2][4][1] = int_stack + 4182;
 Libderiv->deriv_classes[2][5][1] = int_stack + 4272;
 Libderiv->deriv_classes[3][4][1] = int_stack + 4398;
 Libderiv->deriv_classes[3][5][1] = int_stack + 4548;
 Libderiv->dvrr_classes[2][4] = int_stack + 4758;
 Libderiv->dvrr_classes[2][5] = int_stack + 4848;
 Libderiv->deriv_classes[2][4][0] = int_stack + 4974;
 Libderiv->deriv_classes[2][5][0] = int_stack + 5064;
 Libderiv->deriv_classes[3][4][0] = int_stack + 5190;
 Libderiv->deriv_classes[3][5][0] = int_stack + 5340;
 memset(int_stack,0,44400);

 Libderiv->dvrr_stack = int_stack + 9240;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpgp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5550,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4758,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5820,int_stack+366,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3096,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+666,int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4758, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+942,int_stack+792, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3096, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+720,int_stack+1242,int_stack+1152, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4758, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6270,int_stack+1518,int_stack+1368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3096, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+990,int_stack+1818,int_stack+1728, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+2094,int_stack+1944, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1710,int_stack+2394,int_stack+2304, 0.0, zero_stack, 1.0, int_stack+4758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1980,int_stack+2670,int_stack+2520, 0.0, zero_stack, 1.0, int_stack+3096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2430,int_stack+2970,int_stack+2880, 1.0, int_stack+4758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6720,int_stack+3396,int_stack+3246, 1.0, int_stack+3096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2700,int_stack+4848,int_stack+4758,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2970,int_stack+3696,int_stack+3606,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3240,int_stack+3972,int_stack+3822,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3690,int_stack+4272,int_stack+4182,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7170,int_stack+4548,int_stack+4398,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3960,int_stack+5064,int_stack+4974,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4230,int_stack+5340,int_stack+5190,10);
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4680,int_stack+5820,int_stack+5550,45);
     Libderiv->ABCD[11] = int_stack + 4680;
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7620,int_stack+270,int_stack+0,45);
     Libderiv->ABCD[10] = int_stack + 7620;
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8430,int_stack+6270,int_stack+720,45);
     Libderiv->ABCD[9] = int_stack + 8430;
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+1260,int_stack+990,45);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+810,int_stack+1980,int_stack+1710,45);
     Libderiv->ABCD[7] = int_stack + 810;
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1620,int_stack+6720,int_stack+2430,45);
     Libderiv->ABCD[6] = int_stack + 1620;
 /*--- compute (dp|gp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5490,int_stack+3240,int_stack+2970, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[2] = int_stack + 5490;
 /*--- compute (dp|gp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6300,int_stack+7170,int_stack+3690, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[1] = int_stack + 6300;
 /*--- compute (dp|gp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+2970,int_stack+4230,int_stack+3960, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[0] = int_stack + 2970;

}
