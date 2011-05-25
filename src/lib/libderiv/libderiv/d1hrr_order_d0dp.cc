#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0dp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|dp) integrals */

void d1hrr_order_d0dp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][2][11] = int_stack + 0;
 Libderiv->deriv_classes[2][3][11] = int_stack + 36;
 Libderiv->deriv_classes[2][2][10] = int_stack + 96;
 Libderiv->deriv_classes[2][3][10] = int_stack + 132;
 Libderiv->deriv_classes[2][2][9] = int_stack + 192;
 Libderiv->deriv_classes[2][3][9] = int_stack + 228;
 Libderiv->deriv_classes[2][2][8] = int_stack + 288;
 Libderiv->deriv_classes[2][3][8] = int_stack + 324;
 Libderiv->deriv_classes[2][2][7] = int_stack + 384;
 Libderiv->deriv_classes[2][3][7] = int_stack + 420;
 Libderiv->dvrr_classes[2][2] = int_stack + 480;
 Libderiv->deriv_classes[2][2][6] = int_stack + 516;
 Libderiv->deriv_classes[2][3][6] = int_stack + 552;
 Libderiv->deriv_classes[2][2][2] = int_stack + 612;
 Libderiv->deriv_classes[2][3][2] = int_stack + 648;
 Libderiv->deriv_classes[2][2][1] = int_stack + 708;
 Libderiv->deriv_classes[2][3][1] = int_stack + 744;
 Libderiv->deriv_classes[2][2][0] = int_stack + 804;
 Libderiv->deriv_classes[2][3][0] = int_stack + 840;
 memset(int_stack,0,7200);

 Libderiv->dvrr_stack = int_stack + 1116;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0dp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+900,int_stack+36,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+480,6);
     Libderiv->ABCD[11] = int_stack + 900;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1008,int_stack+132,int_stack+96, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+480, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 1008;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+228,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+480, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+324,int_stack+288, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 108;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+420,int_stack+384, 0.0, zero_stack, 1.0, int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 216;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+552,int_stack+516, 1.0, int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 324;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+648,int_stack+612,6);
     Libderiv->ABCD[2] = int_stack + 432;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+744,int_stack+708,6);
     Libderiv->ABCD[1] = int_stack + 540;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+840,int_stack+804,6);
     Libderiv->ABCD[0] = int_stack + 648;

}
