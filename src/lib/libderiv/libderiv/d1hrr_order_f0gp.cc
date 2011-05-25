#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0gp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|gp) integrals */

void d1hrr_order_f0gp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][4][10] = int_stack + 360;
 Libderiv->deriv_classes[3][5][10] = int_stack + 510;
 Libderiv->deriv_classes[3][4][9] = int_stack + 720;
 Libderiv->deriv_classes[3][5][9] = int_stack + 870;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1080;
 Libderiv->deriv_classes[3][5][8] = int_stack + 1230;
 Libderiv->deriv_classes[3][4][7] = int_stack + 1440;
 Libderiv->deriv_classes[3][5][7] = int_stack + 1590;
 Libderiv->dvrr_classes[3][4] = int_stack + 1800;
 Libderiv->deriv_classes[3][4][6] = int_stack + 1950;
 Libderiv->deriv_classes[3][5][6] = int_stack + 2100;
 Libderiv->deriv_classes[3][4][2] = int_stack + 2310;
 Libderiv->deriv_classes[3][5][2] = int_stack + 2460;
 Libderiv->deriv_classes[3][4][1] = int_stack + 2670;
 Libderiv->deriv_classes[3][5][1] = int_stack + 2820;
 Libderiv->deriv_classes[3][4][0] = int_stack + 3030;
 Libderiv->deriv_classes[3][5][0] = int_stack + 3180;
 memset(int_stack,0,27120);

 Libderiv->dvrr_stack = int_stack + 4740;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0gp(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3390,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1800,10);
     Libderiv->ABCD[11] = int_stack + 3390;
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3840,int_stack+510,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1800, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 3840;
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+870,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1800, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+1230,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 450;
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+1590,int_stack+1440, 0.0, zero_stack, 1.0, int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 900;
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1350,int_stack+2100,int_stack+1950, 1.0, int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 1350;
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1800,int_stack+2460,int_stack+2310,10);
     Libderiv->ABCD[2] = int_stack + 1800;
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4290,int_stack+2820,int_stack+2670,10);
     Libderiv->ABCD[1] = int_stack + 4290;
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2250,int_stack+3180,int_stack+3030,10);
     Libderiv->ABCD[0] = int_stack + 2250;

}
