#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpfp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|fp) integrals */

void d1hrr_order_fpfp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 100;
 Libderiv->deriv_classes[4][3][11] = int_stack + 250;
 Libderiv->deriv_classes[4][4][11] = int_stack + 400;
 Libderiv->deriv_classes[3][3][10] = int_stack + 625;
 Libderiv->deriv_classes[3][4][10] = int_stack + 725;
 Libderiv->deriv_classes[4][3][10] = int_stack + 875;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1025;
 Libderiv->deriv_classes[3][3][9] = int_stack + 1250;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1350;
 Libderiv->deriv_classes[4][3][9] = int_stack + 1500;
 Libderiv->deriv_classes[4][4][9] = int_stack + 1650;
 Libderiv->deriv_classes[3][3][8] = int_stack + 1875;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1975;
 Libderiv->deriv_classes[4][3][8] = int_stack + 2125;
 Libderiv->deriv_classes[4][4][8] = int_stack + 2275;
 Libderiv->deriv_classes[3][3][7] = int_stack + 2500;
 Libderiv->deriv_classes[3][4][7] = int_stack + 2600;
 Libderiv->deriv_classes[4][3][7] = int_stack + 2750;
 Libderiv->deriv_classes[4][4][7] = int_stack + 2900;
 Libderiv->deriv_classes[3][3][6] = int_stack + 3125;
 Libderiv->deriv_classes[3][4][6] = int_stack + 3225;
 Libderiv->dvrr_classes[4][3] = int_stack + 3375;
 Libderiv->deriv_classes[4][3][6] = int_stack + 3525;
 Libderiv->deriv_classes[4][4][6] = int_stack + 3675;
 Libderiv->deriv_classes[3][3][2] = int_stack + 3900;
 Libderiv->deriv_classes[3][4][2] = int_stack + 4000;
 Libderiv->deriv_classes[4][3][2] = int_stack + 4150;
 Libderiv->deriv_classes[4][4][2] = int_stack + 4300;
 Libderiv->deriv_classes[3][3][1] = int_stack + 4525;
 Libderiv->deriv_classes[3][4][1] = int_stack + 4625;
 Libderiv->deriv_classes[4][3][1] = int_stack + 4775;
 Libderiv->deriv_classes[4][4][1] = int_stack + 4925;
 Libderiv->dvrr_classes[3][3] = int_stack + 5150;
 Libderiv->dvrr_classes[3][4] = int_stack + 5250;
 Libderiv->deriv_classes[3][3][0] = int_stack + 5400;
 Libderiv->deriv_classes[3][4][0] = int_stack + 5500;
 Libderiv->deriv_classes[4][3][0] = int_stack + 5650;
 Libderiv->deriv_classes[4][4][0] = int_stack + 5800;
 memset(int_stack,0,48200);

 Libderiv->dvrr_stack = int_stack + 10375;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpfp(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6025,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5150,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6325,int_stack+400,int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3375,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+725,int_stack+625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5150, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+1025,int_stack+875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3375, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+750,int_stack+1350,int_stack+1250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5150, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1050,int_stack+1650,int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3375, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1500,int_stack+1975,int_stack+1875, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6775,int_stack+2275,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2600,int_stack+2500, 0.0, zero_stack, 1.0, int_stack+5150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+2900,int_stack+2750, 0.0, zero_stack, 1.0, int_stack+3375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+3225,int_stack+3125, 1.0, int_stack+5150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+3675,int_stack+3525, 1.0, int_stack+3375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3300,int_stack+5250,int_stack+5150,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4000,int_stack+3900,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7225,int_stack+4300,int_stack+4150,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3900,int_stack+4625,int_stack+4525,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4200,int_stack+4925,int_stack+4775,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4650,int_stack+5500,int_stack+5400,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4950,int_stack+5800,int_stack+5650,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+7675,int_stack+6325,int_stack+6025,30);
     Libderiv->ABCD[11] = int_stack + 7675;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5400,int_stack+300,int_stack+0,30);
     Libderiv->ABCD[10] = int_stack + 5400;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8575,int_stack+1050,int_stack+750,30);
     Libderiv->ABCD[9] = int_stack + 8575;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+6775,int_stack+1500,30);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+900,int_stack+2100,int_stack+1800,30);
     Libderiv->ABCD[7] = int_stack + 900;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+6300,int_stack+2850,int_stack+2550,30);
     Libderiv->ABCD[6] = int_stack + 6300;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+1800,int_stack+7225,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[2] = int_stack + 1800;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+9475,int_stack+4200,int_stack+3900, 0.0, zero_stack, 1.0, int_stack+3300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[1] = int_stack + 9475;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+3600,int_stack+4950,int_stack+4650, 1.0, int_stack+3300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[0] = int_stack + 3600;

}
