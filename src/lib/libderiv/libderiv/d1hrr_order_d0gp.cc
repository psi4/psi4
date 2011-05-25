#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0gp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|gp) integrals */

void d1hrr_order_d0gp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[2][4][10] = int_stack + 216;
 Libderiv->deriv_classes[2][5][10] = int_stack + 306;
 Libderiv->deriv_classes[2][4][9] = int_stack + 432;
 Libderiv->deriv_classes[2][5][9] = int_stack + 522;
 Libderiv->deriv_classes[2][4][8] = int_stack + 648;
 Libderiv->deriv_classes[2][5][8] = int_stack + 738;
 Libderiv->deriv_classes[2][4][7] = int_stack + 864;
 Libderiv->deriv_classes[2][5][7] = int_stack + 954;
 Libderiv->dvrr_classes[2][4] = int_stack + 1080;
 Libderiv->deriv_classes[2][4][6] = int_stack + 1170;
 Libderiv->deriv_classes[2][5][6] = int_stack + 1260;
 Libderiv->deriv_classes[2][4][2] = int_stack + 1386;
 Libderiv->deriv_classes[2][5][2] = int_stack + 1476;
 Libderiv->deriv_classes[2][4][1] = int_stack + 1602;
 Libderiv->deriv_classes[2][5][1] = int_stack + 1692;
 Libderiv->deriv_classes[2][4][0] = int_stack + 1818;
 Libderiv->deriv_classes[2][5][0] = int_stack + 1908;
 memset(int_stack,0,16272);

 Libderiv->dvrr_stack = int_stack + 2844;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0gp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2034,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080,6);
     Libderiv->ABCD[11] = int_stack + 2034;
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2304,int_stack+306,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 2304;
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+522,int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+738,int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 270;
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+954,int_stack+864, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 540;
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+810,int_stack+1260,int_stack+1170, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 810;
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1080,int_stack+1476,int_stack+1386,6);
     Libderiv->ABCD[2] = int_stack + 1080;
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2574,int_stack+1692,int_stack+1602,6);
     Libderiv->ABCD[1] = int_stack + 2574;
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1350,int_stack+1908,int_stack+1818,6);
     Libderiv->ABCD[0] = int_stack + 1350;

}
