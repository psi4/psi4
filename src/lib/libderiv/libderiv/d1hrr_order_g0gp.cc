#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_g0gp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (g0|gp) integrals */

void d1hrr_order_g0gp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][4][10] = int_stack + 540;
 Libderiv->deriv_classes[4][5][10] = int_stack + 765;
 Libderiv->deriv_classes[4][4][9] = int_stack + 1080;
 Libderiv->deriv_classes[4][5][9] = int_stack + 1305;
 Libderiv->deriv_classes[4][4][8] = int_stack + 1620;
 Libderiv->deriv_classes[4][5][8] = int_stack + 1845;
 Libderiv->deriv_classes[4][4][7] = int_stack + 2160;
 Libderiv->deriv_classes[4][5][7] = int_stack + 2385;
 Libderiv->dvrr_classes[4][4] = int_stack + 2700;
 Libderiv->deriv_classes[4][4][6] = int_stack + 2925;
 Libderiv->deriv_classes[4][5][6] = int_stack + 3150;
 Libderiv->deriv_classes[4][4][2] = int_stack + 3465;
 Libderiv->deriv_classes[4][5][2] = int_stack + 3690;
 Libderiv->deriv_classes[4][4][1] = int_stack + 4005;
 Libderiv->deriv_classes[4][5][1] = int_stack + 4230;
 Libderiv->deriv_classes[4][4][0] = int_stack + 4545;
 Libderiv->deriv_classes[4][5][0] = int_stack + 4770;
 memset(int_stack,0,40680);

 Libderiv->dvrr_stack = int_stack + 7110;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_g0gp(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5085,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700,15);
     Libderiv->ABCD[11] = int_stack + 5085;
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5760,int_stack+765,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack,15);
     Libderiv->ABCD[10] = int_stack + 5760;
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1305,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+675,int_stack+1845,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[8] = int_stack + 675;
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1350,int_stack+2385,int_stack+2160, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[7] = int_stack + 1350;
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2025,int_stack+3150,int_stack+2925, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[6] = int_stack + 2025;
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2700,int_stack+3690,int_stack+3465,15);
     Libderiv->ABCD[2] = int_stack + 2700;
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6435,int_stack+4230,int_stack+4005,15);
     Libderiv->ABCD[1] = int_stack + 6435;
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3375,int_stack+4770,int_stack+4545,15);
     Libderiv->ABCD[0] = int_stack + 3375;

}
