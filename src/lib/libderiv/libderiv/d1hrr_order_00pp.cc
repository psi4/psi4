#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00pp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|pp) integrals */

void d1hrr_order_00pp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][1][11] = int_stack + 0;
 Libderiv->deriv_classes[0][2][11] = int_stack + 3;
 Libderiv->deriv_classes[0][1][10] = int_stack + 9;
 Libderiv->deriv_classes[0][2][10] = int_stack + 12;
 Libderiv->deriv_classes[0][1][9] = int_stack + 18;
 Libderiv->deriv_classes[0][2][9] = int_stack + 21;
 Libderiv->deriv_classes[0][1][8] = int_stack + 27;
 Libderiv->deriv_classes[0][2][8] = int_stack + 30;
 Libderiv->deriv_classes[0][1][7] = int_stack + 36;
 Libderiv->deriv_classes[0][2][7] = int_stack + 39;
 Libderiv->dvrr_classes[0][1] = int_stack + 45;
 Libderiv->deriv_classes[0][1][6] = int_stack + 48;
 Libderiv->deriv_classes[0][2][6] = int_stack + 51;
 Libderiv->deriv_classes[0][1][2] = int_stack + 57;
 Libderiv->deriv_classes[0][2][2] = int_stack + 60;
 Libderiv->deriv_classes[0][1][1] = int_stack + 66;
 Libderiv->deriv_classes[0][2][1] = int_stack + 69;
 Libderiv->deriv_classes[0][1][0] = int_stack + 75;
 Libderiv->deriv_classes[0][2][0] = int_stack + 78;
 memset(int_stack,0,672);

 Libderiv->dvrr_stack = int_stack + 93;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00pp(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+84,int_stack+3,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45,1);
     Libderiv->ABCD[11] = int_stack + 84;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+0,int_stack+12,int_stack+9, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 0;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+9,int_stack+21,int_stack+18, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 9;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+18,int_stack+30,int_stack+27, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 18;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+27,int_stack+39,int_stack+36, 0.0, zero_stack, 1.0, int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 27;
 /*--- compute (00|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+36,int_stack+51,int_stack+48, 1.0, int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 36;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+45,int_stack+60,int_stack+57,1);
     Libderiv->ABCD[2] = int_stack + 45;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+54,int_stack+69,int_stack+66,1);
     Libderiv->ABCD[1] = int_stack + 54;
 /*--- compute (00|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+63,int_stack+78,int_stack+75,1);
     Libderiv->ABCD[0] = int_stack + 63;

}
