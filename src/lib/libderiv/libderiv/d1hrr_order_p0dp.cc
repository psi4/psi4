#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0dp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|dp) integrals */

void d1hrr_order_p0dp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][2][11] = int_stack + 0;
 Libderiv->deriv_classes[1][3][11] = int_stack + 18;
 Libderiv->deriv_classes[1][2][10] = int_stack + 48;
 Libderiv->deriv_classes[1][3][10] = int_stack + 66;
 Libderiv->deriv_classes[1][2][9] = int_stack + 96;
 Libderiv->deriv_classes[1][3][9] = int_stack + 114;
 Libderiv->deriv_classes[1][2][8] = int_stack + 144;
 Libderiv->deriv_classes[1][3][8] = int_stack + 162;
 Libderiv->deriv_classes[1][2][7] = int_stack + 192;
 Libderiv->deriv_classes[1][3][7] = int_stack + 210;
 Libderiv->dvrr_classes[1][2] = int_stack + 240;
 Libderiv->deriv_classes[1][2][6] = int_stack + 258;
 Libderiv->deriv_classes[1][3][6] = int_stack + 276;
 Libderiv->deriv_classes[1][2][2] = int_stack + 306;
 Libderiv->deriv_classes[1][3][2] = int_stack + 324;
 Libderiv->deriv_classes[1][2][1] = int_stack + 354;
 Libderiv->deriv_classes[1][3][1] = int_stack + 372;
 Libderiv->deriv_classes[1][2][0] = int_stack + 402;
 Libderiv->deriv_classes[1][3][0] = int_stack + 420;
 memset(int_stack,0,3600);

 Libderiv->dvrr_stack = int_stack + 558;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0dp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+450,int_stack+18,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+240,3);
     Libderiv->ABCD[11] = int_stack + 450;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+504,int_stack+66,int_stack+48, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+240, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 504;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+114,int_stack+96, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+240, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+162,int_stack+144, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 54;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+210,int_stack+192, 0.0, zero_stack, 1.0, int_stack+240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 108;
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+162,int_stack+276,int_stack+258, 1.0, int_stack+240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 162;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+324,int_stack+306,3);
     Libderiv->ABCD[2] = int_stack + 216;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+372,int_stack+354,3);
     Libderiv->ABCD[1] = int_stack + 270;
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+420,int_stack+402,3);
     Libderiv->ABCD[0] = int_stack + 324;

}
