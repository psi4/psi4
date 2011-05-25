#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0gp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|gp) integrals */

void d1hrr_order_p0gp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[1][5][11] = int_stack + 45;
 Libderiv->deriv_classes[1][4][10] = int_stack + 108;
 Libderiv->deriv_classes[1][5][10] = int_stack + 153;
 Libderiv->deriv_classes[1][4][9] = int_stack + 216;
 Libderiv->deriv_classes[1][5][9] = int_stack + 261;
 Libderiv->deriv_classes[1][4][8] = int_stack + 324;
 Libderiv->deriv_classes[1][5][8] = int_stack + 369;
 Libderiv->deriv_classes[1][4][7] = int_stack + 432;
 Libderiv->deriv_classes[1][5][7] = int_stack + 477;
 Libderiv->dvrr_classes[1][4] = int_stack + 540;
 Libderiv->deriv_classes[1][4][6] = int_stack + 585;
 Libderiv->deriv_classes[1][5][6] = int_stack + 630;
 Libderiv->deriv_classes[1][4][2] = int_stack + 693;
 Libderiv->deriv_classes[1][5][2] = int_stack + 738;
 Libderiv->deriv_classes[1][4][1] = int_stack + 801;
 Libderiv->deriv_classes[1][5][1] = int_stack + 846;
 Libderiv->deriv_classes[1][4][0] = int_stack + 909;
 Libderiv->deriv_classes[1][5][0] = int_stack + 954;
 memset(int_stack,0,8136);

 Libderiv->dvrr_stack = int_stack + 1422;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0gp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1017,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540,3);
     Libderiv->ABCD[11] = int_stack + 1017;
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1152,int_stack+153,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 1152;
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+261,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+135,int_stack+369,int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 135;
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+477,int_stack+432, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 270;
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+405,int_stack+630,int_stack+585, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 405;
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+738,int_stack+693,3);
     Libderiv->ABCD[2] = int_stack + 540;
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1287,int_stack+846,int_stack+801,3);
     Libderiv->ABCD[1] = int_stack + 1287;
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+675,int_stack+954,int_stack+909,3);
     Libderiv->ABCD[0] = int_stack + 675;

}
