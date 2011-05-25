#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0fp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|fp) integrals */

void d1hrr_order_d0fp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[2][3][10] = int_stack + 150;
 Libderiv->deriv_classes[2][4][10] = int_stack + 210;
 Libderiv->deriv_classes[2][3][9] = int_stack + 300;
 Libderiv->deriv_classes[2][4][9] = int_stack + 360;
 Libderiv->deriv_classes[2][3][8] = int_stack + 450;
 Libderiv->deriv_classes[2][4][8] = int_stack + 510;
 Libderiv->deriv_classes[2][3][7] = int_stack + 600;
 Libderiv->deriv_classes[2][4][7] = int_stack + 660;
 Libderiv->dvrr_classes[2][3] = int_stack + 750;
 Libderiv->deriv_classes[2][3][6] = int_stack + 810;
 Libderiv->deriv_classes[2][4][6] = int_stack + 870;
 Libderiv->deriv_classes[2][3][2] = int_stack + 960;
 Libderiv->deriv_classes[2][4][2] = int_stack + 1020;
 Libderiv->deriv_classes[2][3][1] = int_stack + 1110;
 Libderiv->deriv_classes[2][4][1] = int_stack + 1170;
 Libderiv->deriv_classes[2][3][0] = int_stack + 1260;
 Libderiv->deriv_classes[2][4][0] = int_stack + 1320;
 memset(int_stack,0,11280);

 Libderiv->dvrr_stack = int_stack + 1770;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0fp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1410,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750,6);
     Libderiv->ABCD[11] = int_stack + 1410;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1590,int_stack+210,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 1590;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+360,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+510,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 180;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+660,int_stack+600, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 360;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+870,int_stack+810, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 540;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+1020,int_stack+960,6);
     Libderiv->ABCD[2] = int_stack + 720;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1170,int_stack+1110,6);
     Libderiv->ABCD[1] = int_stack + 900;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+1320,int_stack+1260,6);
     Libderiv->ABCD[0] = int_stack + 1080;

}
