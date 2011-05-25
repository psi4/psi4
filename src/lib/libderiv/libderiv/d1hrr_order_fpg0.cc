#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpg0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|g0) integrals */

void d1hrr_order_fpg0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][4][11] = int_stack + 150;
 Libderiv->deriv_classes[3][4][10] = int_stack + 375;
 Libderiv->deriv_classes[4][4][10] = int_stack + 525;
 Libderiv->deriv_classes[3][4][9] = int_stack + 750;
 Libderiv->deriv_classes[4][4][9] = int_stack + 900;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1125;
 Libderiv->deriv_classes[4][4][8] = int_stack + 1275;
 Libderiv->deriv_classes[3][4][7] = int_stack + 1500;
 Libderiv->deriv_classes[4][4][7] = int_stack + 1650;
 Libderiv->deriv_classes[3][4][6] = int_stack + 1875;
 Libderiv->deriv_classes[4][4][6] = int_stack + 2025;
 Libderiv->deriv_classes[3][4][2] = int_stack + 2250;
 Libderiv->deriv_classes[4][4][2] = int_stack + 2400;
 Libderiv->deriv_classes[3][4][1] = int_stack + 2625;
 Libderiv->deriv_classes[4][4][1] = int_stack + 2775;
 Libderiv->dvrr_classes[3][4] = int_stack + 3000;
 Libderiv->deriv_classes[3][4][0] = int_stack + 3150;
 Libderiv->deriv_classes[4][4][0] = int_stack + 3300;
 memset(int_stack,0,28200);

 Libderiv->dvrr_stack = int_stack + 4875;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpg0(Libderiv, Data);
   Data++;
 }

 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3525,int_stack+150,int_stack+0,15);
     Libderiv->ABCD[11] = int_stack + 3525;
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3975,int_stack+525,int_stack+375,15);
     Libderiv->ABCD[10] = int_stack + 3975;
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+900,int_stack+750,15);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+450,int_stack+1275,int_stack+1125,15);
     Libderiv->ABCD[8] = int_stack + 450;
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+900,int_stack+1650,int_stack+1500,15);
     Libderiv->ABCD[7] = int_stack + 900;
 /*--- compute (fp|g0) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1350,int_stack+2025,int_stack+1875,15);
     Libderiv->ABCD[6] = int_stack + 1350;
 /*--- compute (fp|g0) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+1800,int_stack+2400,int_stack+2250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[2] = int_stack + 1800;
 /*--- compute (fp|g0) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+4425,int_stack+2775,int_stack+2625, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[1] = int_stack + 4425;
 /*--- compute (fp|g0) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+2250,int_stack+3300,int_stack+3150, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[0] = int_stack + 2250;

}
