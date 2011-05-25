#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppf0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|f0) integrals */

void d1hrr_order_ppf0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][3][11] = int_stack + 30;
 Libderiv->deriv_classes[1][3][10] = int_stack + 90;
 Libderiv->deriv_classes[2][3][10] = int_stack + 120;
 Libderiv->deriv_classes[1][3][9] = int_stack + 180;
 Libderiv->deriv_classes[2][3][9] = int_stack + 210;
 Libderiv->deriv_classes[1][3][8] = int_stack + 270;
 Libderiv->deriv_classes[2][3][8] = int_stack + 300;
 Libderiv->deriv_classes[1][3][7] = int_stack + 360;
 Libderiv->deriv_classes[2][3][7] = int_stack + 390;
 Libderiv->deriv_classes[1][3][6] = int_stack + 450;
 Libderiv->deriv_classes[2][3][6] = int_stack + 480;
 Libderiv->deriv_classes[1][3][2] = int_stack + 540;
 Libderiv->deriv_classes[2][3][2] = int_stack + 570;
 Libderiv->deriv_classes[1][3][1] = int_stack + 630;
 Libderiv->deriv_classes[2][3][1] = int_stack + 660;
 Libderiv->dvrr_classes[1][3] = int_stack + 720;
 Libderiv->deriv_classes[1][3][0] = int_stack + 750;
 Libderiv->deriv_classes[2][3][0] = int_stack + 780;
 memset(int_stack,0,6720);

 Libderiv->dvrr_stack = int_stack + 930;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppf0(Libderiv, Data);
   Data++;
 }

 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+840,int_stack+30,int_stack+0,10);
     Libderiv->ABCD[11] = int_stack + 840;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+120,int_stack+90,10);
     Libderiv->ABCD[10] = int_stack + 0;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+90,int_stack+210,int_stack+180,10);
     Libderiv->ABCD[9] = int_stack + 90;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+180,int_stack+300,int_stack+270,10);
     Libderiv->ABCD[8] = int_stack + 180;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+270,int_stack+390,int_stack+360,10);
     Libderiv->ABCD[7] = int_stack + 270;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+360,int_stack+480,int_stack+450,10);
     Libderiv->ABCD[6] = int_stack + 360;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+450,int_stack+570,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[2] = int_stack + 450;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+540,int_stack+660,int_stack+630, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[1] = int_stack + 540;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+630,int_stack+780,int_stack+750, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[0] = int_stack + 630;

}
