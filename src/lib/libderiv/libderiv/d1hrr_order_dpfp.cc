#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpfp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|fp) integrals */

void d1hrr_order_dpfp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[3][3][11] = int_stack + 150;
 Libderiv->deriv_classes[3][4][11] = int_stack + 250;
 Libderiv->deriv_classes[2][3][10] = int_stack + 400;
 Libderiv->deriv_classes[2][4][10] = int_stack + 460;
 Libderiv->deriv_classes[3][3][10] = int_stack + 550;
 Libderiv->deriv_classes[3][4][10] = int_stack + 650;
 Libderiv->deriv_classes[2][3][9] = int_stack + 800;
 Libderiv->deriv_classes[2][4][9] = int_stack + 860;
 Libderiv->deriv_classes[3][3][9] = int_stack + 950;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1050;
 Libderiv->deriv_classes[2][3][8] = int_stack + 1200;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1260;
 Libderiv->deriv_classes[3][3][8] = int_stack + 1350;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1450;
 Libderiv->deriv_classes[2][3][7] = int_stack + 1600;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1660;
 Libderiv->deriv_classes[3][3][7] = int_stack + 1750;
 Libderiv->deriv_classes[3][4][7] = int_stack + 1850;
 Libderiv->deriv_classes[2][3][6] = int_stack + 2000;
 Libderiv->deriv_classes[2][4][6] = int_stack + 2060;
 Libderiv->dvrr_classes[3][3] = int_stack + 2150;
 Libderiv->deriv_classes[3][3][6] = int_stack + 2250;
 Libderiv->deriv_classes[3][4][6] = int_stack + 2350;
 Libderiv->deriv_classes[2][3][2] = int_stack + 2500;
 Libderiv->deriv_classes[2][4][2] = int_stack + 2560;
 Libderiv->deriv_classes[3][3][2] = int_stack + 2650;
 Libderiv->deriv_classes[3][4][2] = int_stack + 2750;
 Libderiv->deriv_classes[2][3][1] = int_stack + 2900;
 Libderiv->deriv_classes[2][4][1] = int_stack + 2960;
 Libderiv->deriv_classes[3][3][1] = int_stack + 3050;
 Libderiv->deriv_classes[3][4][1] = int_stack + 3150;
 Libderiv->dvrr_classes[2][3] = int_stack + 3300;
 Libderiv->dvrr_classes[2][4] = int_stack + 3360;
 Libderiv->deriv_classes[2][3][0] = int_stack + 3450;
 Libderiv->deriv_classes[2][4][0] = int_stack + 3510;
 Libderiv->deriv_classes[3][3][0] = int_stack + 3600;
 Libderiv->deriv_classes[3][4][0] = int_stack + 3700;
 memset(int_stack,0,30800);

 Libderiv->dvrr_stack = int_stack + 6550;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpfp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3850,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3300,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4030,int_stack+250,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2150,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+460,int_stack+400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3300, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+650,int_stack+550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2150, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+480,int_stack+860,int_stack+800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3300, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4330,int_stack+1050,int_stack+950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2150, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+660,int_stack+1260,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+840,int_stack+1450,int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1140,int_stack+1660,int_stack+1600, 0.0, zero_stack, 1.0, int_stack+3300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1320,int_stack+1850,int_stack+1750, 0.0, zero_stack, 1.0, int_stack+2150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1620,int_stack+2060,int_stack+2000, 1.0, int_stack+3300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2350,int_stack+2250, 1.0, int_stack+2150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+3360,int_stack+3300,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2280,int_stack+2560,int_stack+2500,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4630,int_stack+2750,int_stack+2650,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2460,int_stack+2960,int_stack+2900,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2640,int_stack+3150,int_stack+3050,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2940,int_stack+3510,int_stack+3450,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3120,int_stack+3700,int_stack+3600,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4930,int_stack+4030,int_stack+3850,30);
     Libderiv->ABCD[11] = int_stack + 4930;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3420,int_stack+180,int_stack+0,30);
     Libderiv->ABCD[10] = int_stack + 3420;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5470,int_stack+4330,int_stack+480,30);
     Libderiv->ABCD[9] = int_stack + 5470;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+840,int_stack+660,30);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+540,int_stack+1320,int_stack+1140,30);
     Libderiv->ABCD[7] = int_stack + 540;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1080,int_stack+1800,int_stack+1620,30);
     Libderiv->ABCD[6] = int_stack + 1080;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+3960,int_stack+4630,int_stack+2280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[2] = int_stack + 3960;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6010,int_stack+2640,int_stack+2460, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[1] = int_stack + 6010;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+2280,int_stack+3120,int_stack+2940, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[0] = int_stack + 2280;

}
