#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00gp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|gp) integrals */

void d1hrr_order_00gp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][4][11] = int_stack + 0;
 Libderiv->deriv_classes[0][5][11] = int_stack + 15;
 Libderiv->deriv_classes[0][4][10] = int_stack + 36;
 Libderiv->deriv_classes[0][5][10] = int_stack + 51;
 Libderiv->deriv_classes[0][4][9] = int_stack + 72;
 Libderiv->deriv_classes[0][5][9] = int_stack + 87;
 Libderiv->deriv_classes[0][4][8] = int_stack + 108;
 Libderiv->deriv_classes[0][5][8] = int_stack + 123;
 Libderiv->deriv_classes[0][4][7] = int_stack + 144;
 Libderiv->deriv_classes[0][5][7] = int_stack + 159;
 Libderiv->dvrr_classes[0][4] = int_stack + 180;
 Libderiv->deriv_classes[0][4][6] = int_stack + 195;
 Libderiv->deriv_classes[0][5][6] = int_stack + 210;
 Libderiv->deriv_classes[0][4][2] = int_stack + 231;
 Libderiv->deriv_classes[0][5][2] = int_stack + 246;
 Libderiv->deriv_classes[0][4][1] = int_stack + 267;
 Libderiv->deriv_classes[0][5][1] = int_stack + 282;
 Libderiv->deriv_classes[0][4][0] = int_stack + 303;
 Libderiv->deriv_classes[0][5][0] = int_stack + 318;
 memset(int_stack,0,2712);

 Libderiv->dvrr_stack = int_stack + 474;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00gp(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+339,int_stack+15,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180,1);
     Libderiv->ABCD[11] = int_stack + 339;
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+384,int_stack+51,int_stack+36, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 384;
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+87,int_stack+72, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45,int_stack+123,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 45;
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+159,int_stack+144, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 90;
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+135,int_stack+210,int_stack+195, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 135;
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+180,int_stack+246,int_stack+231,1);
     Libderiv->ABCD[2] = int_stack + 180;
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+429,int_stack+282,int_stack+267,1);
     Libderiv->ABCD[1] = int_stack + 429;
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+225,int_stack+318,int_stack+303,1);
     Libderiv->ABCD[0] = int_stack + 225;

}
