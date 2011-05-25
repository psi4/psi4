#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppd0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|d0) integrals */

void d1hrr_order_ppd0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][2][11] = int_stack + 0;
 Libderiv->deriv_classes[2][2][11] = int_stack + 18;
 Libderiv->deriv_classes[1][2][10] = int_stack + 54;
 Libderiv->deriv_classes[2][2][10] = int_stack + 72;
 Libderiv->deriv_classes[1][2][9] = int_stack + 108;
 Libderiv->deriv_classes[2][2][9] = int_stack + 126;
 Libderiv->deriv_classes[1][2][8] = int_stack + 162;
 Libderiv->deriv_classes[2][2][8] = int_stack + 180;
 Libderiv->deriv_classes[1][2][7] = int_stack + 216;
 Libderiv->deriv_classes[2][2][7] = int_stack + 234;
 Libderiv->deriv_classes[1][2][6] = int_stack + 270;
 Libderiv->deriv_classes[2][2][6] = int_stack + 288;
 Libderiv->deriv_classes[1][2][2] = int_stack + 324;
 Libderiv->deriv_classes[2][2][2] = int_stack + 342;
 Libderiv->deriv_classes[1][2][1] = int_stack + 378;
 Libderiv->deriv_classes[2][2][1] = int_stack + 396;
 Libderiv->dvrr_classes[1][2] = int_stack + 432;
 Libderiv->deriv_classes[1][2][0] = int_stack + 450;
 Libderiv->deriv_classes[2][2][0] = int_stack + 468;
 memset(int_stack,0,4032);

 Libderiv->dvrr_stack = int_stack + 558;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppd0(Libderiv, Data);
   Data++;
 }

 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+504,int_stack+18,int_stack+0,6);
     Libderiv->ABCD[11] = int_stack + 504;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+72,int_stack+54,6);
     Libderiv->ABCD[10] = int_stack + 0;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+54,int_stack+126,int_stack+108,6);
     Libderiv->ABCD[9] = int_stack + 54;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+108,int_stack+180,int_stack+162,6);
     Libderiv->ABCD[8] = int_stack + 108;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+162,int_stack+234,int_stack+216,6);
     Libderiv->ABCD[7] = int_stack + 162;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+216,int_stack+288,int_stack+270,6);
     Libderiv->ABCD[6] = int_stack + 216;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+270,int_stack+342,int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[2] = int_stack + 270;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+324,int_stack+396,int_stack+378, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[1] = int_stack + 324;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+378,int_stack+468,int_stack+450, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[0] = int_stack + 378;

}
