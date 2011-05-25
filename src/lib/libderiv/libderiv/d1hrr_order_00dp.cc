#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00dp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|dp) integrals */

void d1hrr_order_00dp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][2][11] = int_stack + 0;
 Libderiv->deriv_classes[0][3][11] = int_stack + 6;
 Libderiv->deriv_classes[0][2][10] = int_stack + 16;
 Libderiv->deriv_classes[0][3][10] = int_stack + 22;
 Libderiv->deriv_classes[0][2][9] = int_stack + 32;
 Libderiv->deriv_classes[0][3][9] = int_stack + 38;
 Libderiv->deriv_classes[0][2][8] = int_stack + 48;
 Libderiv->deriv_classes[0][3][8] = int_stack + 54;
 Libderiv->deriv_classes[0][2][7] = int_stack + 64;
 Libderiv->deriv_classes[0][3][7] = int_stack + 70;
 Libderiv->dvrr_classes[0][2] = int_stack + 80;
 Libderiv->deriv_classes[0][2][6] = int_stack + 86;
 Libderiv->deriv_classes[0][3][6] = int_stack + 92;
 Libderiv->deriv_classes[0][2][2] = int_stack + 102;
 Libderiv->deriv_classes[0][3][2] = int_stack + 108;
 Libderiv->deriv_classes[0][2][1] = int_stack + 118;
 Libderiv->deriv_classes[0][3][1] = int_stack + 124;
 Libderiv->deriv_classes[0][2][0] = int_stack + 134;
 Libderiv->deriv_classes[0][3][0] = int_stack + 140;
 memset(int_stack,0,1200);

 Libderiv->dvrr_stack = int_stack + 186;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00dp(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+150,int_stack+6,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80,1);
     Libderiv->ABCD[11] = int_stack + 150;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+168,int_stack+22,int_stack+16, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 168;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+38,int_stack+32, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18,int_stack+54,int_stack+48, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 18;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+36,int_stack+70,int_stack+64, 0.0, zero_stack, 1.0, int_stack+80, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 36;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+92,int_stack+86, 1.0, int_stack+80, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 54;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+72,int_stack+108,int_stack+102,1);
     Libderiv->ABCD[2] = int_stack + 72;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+90,int_stack+124,int_stack+118,1);
     Libderiv->ABCD[1] = int_stack + 90;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+140,int_stack+134,1);
     Libderiv->ABCD[0] = int_stack + 108;

}
