#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0fp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|fp) integrals */

void d1hrr_order_p0fp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[1][4][11] = int_stack + 30;
 Libderiv->deriv_classes[1][3][10] = int_stack + 75;
 Libderiv->deriv_classes[1][4][10] = int_stack + 105;
 Libderiv->deriv_classes[1][3][9] = int_stack + 150;
 Libderiv->deriv_classes[1][4][9] = int_stack + 180;
 Libderiv->deriv_classes[1][3][8] = int_stack + 225;
 Libderiv->deriv_classes[1][4][8] = int_stack + 255;
 Libderiv->deriv_classes[1][3][7] = int_stack + 300;
 Libderiv->deriv_classes[1][4][7] = int_stack + 330;
 Libderiv->dvrr_classes[1][3] = int_stack + 375;
 Libderiv->deriv_classes[1][3][6] = int_stack + 405;
 Libderiv->deriv_classes[1][4][6] = int_stack + 435;
 Libderiv->deriv_classes[1][3][2] = int_stack + 480;
 Libderiv->deriv_classes[1][4][2] = int_stack + 510;
 Libderiv->deriv_classes[1][3][1] = int_stack + 555;
 Libderiv->deriv_classes[1][4][1] = int_stack + 585;
 Libderiv->deriv_classes[1][3][0] = int_stack + 630;
 Libderiv->deriv_classes[1][4][0] = int_stack + 660;
 memset(int_stack,0,5640);

 Libderiv->dvrr_stack = int_stack + 885;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0fp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+705,int_stack+30,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+375,3);
     Libderiv->ABCD[11] = int_stack + 705;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+795,int_stack+105,int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+375, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 795;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+180,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+375, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+255,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 90;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+330,int_stack+300, 0.0, zero_stack, 1.0, int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 180;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+270,int_stack+435,int_stack+405, 1.0, int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 270;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+510,int_stack+480,3);
     Libderiv->ABCD[2] = int_stack + 360;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+585,int_stack+555,3);
     Libderiv->ABCD[1] = int_stack + 450;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+660,int_stack+630,3);
     Libderiv->ABCD[0] = int_stack + 540;

}
