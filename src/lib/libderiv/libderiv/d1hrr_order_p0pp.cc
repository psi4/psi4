#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0pp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|pp) integrals */

void d1hrr_order_p0pp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][1][11] = int_stack + 0;
 Libderiv->deriv_classes[1][2][11] = int_stack + 9;
 Libderiv->deriv_classes[1][1][10] = int_stack + 27;
 Libderiv->deriv_classes[1][2][10] = int_stack + 36;
 Libderiv->deriv_classes[1][1][9] = int_stack + 54;
 Libderiv->deriv_classes[1][2][9] = int_stack + 63;
 Libderiv->deriv_classes[1][1][8] = int_stack + 81;
 Libderiv->deriv_classes[1][2][8] = int_stack + 90;
 Libderiv->deriv_classes[1][1][7] = int_stack + 108;
 Libderiv->deriv_classes[1][2][7] = int_stack + 117;
 Libderiv->dvrr_classes[1][1] = int_stack + 135;
 Libderiv->deriv_classes[1][1][6] = int_stack + 144;
 Libderiv->deriv_classes[1][2][6] = int_stack + 153;
 Libderiv->deriv_classes[1][1][2] = int_stack + 171;
 Libderiv->deriv_classes[1][2][2] = int_stack + 180;
 Libderiv->deriv_classes[1][1][1] = int_stack + 198;
 Libderiv->deriv_classes[1][2][1] = int_stack + 207;
 Libderiv->deriv_classes[1][1][0] = int_stack + 225;
 Libderiv->deriv_classes[1][2][0] = int_stack + 234;
 memset(int_stack,0,2016);

 Libderiv->dvrr_stack = int_stack + 279;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0pp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+252,int_stack+9,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135,3);
     Libderiv->ABCD[11] = int_stack + 252;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+0,int_stack+36,int_stack+27, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 0;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+27,int_stack+63,int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 27;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+54,int_stack+90,int_stack+81, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 54;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+81,int_stack+117,int_stack+108, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 81;
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+108,int_stack+153,int_stack+144, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 108;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+135,int_stack+180,int_stack+171,3);
     Libderiv->ABCD[2] = int_stack + 135;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+162,int_stack+207,int_stack+198,3);
     Libderiv->ABCD[1] = int_stack + 162;
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+189,int_stack+234,int_stack+225,3);
     Libderiv->ABCD[0] = int_stack + 189;

}
