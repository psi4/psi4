#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0pp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|pp) integrals */

void d1hrr_order_d0pp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][1][11] = int_stack + 0;
 Libderiv->deriv_classes[2][2][11] = int_stack + 18;
 Libderiv->deriv_classes[2][1][10] = int_stack + 54;
 Libderiv->deriv_classes[2][2][10] = int_stack + 72;
 Libderiv->deriv_classes[2][1][9] = int_stack + 108;
 Libderiv->deriv_classes[2][2][9] = int_stack + 126;
 Libderiv->deriv_classes[2][1][8] = int_stack + 162;
 Libderiv->deriv_classes[2][2][8] = int_stack + 180;
 Libderiv->deriv_classes[2][1][7] = int_stack + 216;
 Libderiv->deriv_classes[2][2][7] = int_stack + 234;
 Libderiv->dvrr_classes[2][1] = int_stack + 270;
 Libderiv->deriv_classes[2][1][6] = int_stack + 288;
 Libderiv->deriv_classes[2][2][6] = int_stack + 306;
 Libderiv->deriv_classes[2][1][2] = int_stack + 342;
 Libderiv->deriv_classes[2][2][2] = int_stack + 360;
 Libderiv->deriv_classes[2][1][1] = int_stack + 396;
 Libderiv->deriv_classes[2][2][1] = int_stack + 414;
 Libderiv->deriv_classes[2][1][0] = int_stack + 450;
 Libderiv->deriv_classes[2][2][0] = int_stack + 468;
 memset(int_stack,0,4032);

 Libderiv->dvrr_stack = int_stack + 558;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0pp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+504,int_stack+18,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270,6);
     Libderiv->ABCD[11] = int_stack + 504;
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+0,int_stack+72,int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 0;
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+54,int_stack+126,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 54;
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+108,int_stack+180,int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 108;
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+162,int_stack+234,int_stack+216, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 162;
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+216,int_stack+306,int_stack+288, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 216;
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+270,int_stack+360,int_stack+342,6);
     Libderiv->ABCD[2] = int_stack + 270;
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+324,int_stack+414,int_stack+396,6);
     Libderiv->ABCD[1] = int_stack + 324;
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+378,int_stack+468,int_stack+450,6);
     Libderiv->ABCD[0] = int_stack + 378;

}
