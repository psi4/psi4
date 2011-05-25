#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpgd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|gd) integrals */

void d1hrr_order_fpgd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][6][11] = int_stack + 360;
 Libderiv->deriv_classes[4][4][11] = int_stack + 640;
 Libderiv->deriv_classes[4][5][11] = int_stack + 865;
 Libderiv->deriv_classes[4][6][11] = int_stack + 1180;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1600;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1750;
 Libderiv->deriv_classes[3][6][10] = int_stack + 1960;
 Libderiv->deriv_classes[4][4][10] = int_stack + 2240;
 Libderiv->deriv_classes[4][5][10] = int_stack + 2465;
 Libderiv->deriv_classes[4][6][10] = int_stack + 2780;
 Libderiv->deriv_classes[3][4][9] = int_stack + 3200;
 Libderiv->deriv_classes[3][5][9] = int_stack + 3350;
 Libderiv->deriv_classes[3][6][9] = int_stack + 3560;
 Libderiv->deriv_classes[4][4][9] = int_stack + 3840;
 Libderiv->deriv_classes[4][5][9] = int_stack + 4065;
 Libderiv->deriv_classes[4][6][9] = int_stack + 4380;
 Libderiv->deriv_classes[3][4][8] = int_stack + 4800;
 Libderiv->deriv_classes[3][5][8] = int_stack + 4950;
 Libderiv->deriv_classes[3][6][8] = int_stack + 5160;
 Libderiv->deriv_classes[4][4][8] = int_stack + 5440;
 Libderiv->deriv_classes[4][5][8] = int_stack + 5665;
 Libderiv->deriv_classes[4][6][8] = int_stack + 5980;
 Libderiv->deriv_classes[3][4][7] = int_stack + 6400;
 Libderiv->deriv_classes[3][5][7] = int_stack + 6550;
 Libderiv->deriv_classes[3][6][7] = int_stack + 6760;
 Libderiv->deriv_classes[4][4][7] = int_stack + 7040;
 Libderiv->deriv_classes[4][5][7] = int_stack + 7265;
 Libderiv->deriv_classes[4][6][7] = int_stack + 7580;
 Libderiv->deriv_classes[3][4][6] = int_stack + 8000;
 Libderiv->deriv_classes[3][5][6] = int_stack + 8150;
 Libderiv->deriv_classes[3][6][6] = int_stack + 8360;
 Libderiv->dvrr_classes[4][4] = int_stack + 8640;
 Libderiv->deriv_classes[4][4][6] = int_stack + 8865;
 Libderiv->dvrr_classes[4][5] = int_stack + 9090;
 Libderiv->deriv_classes[4][5][6] = int_stack + 9405;
 Libderiv->deriv_classes[4][6][6] = int_stack + 9720;
 Libderiv->deriv_classes[3][4][2] = int_stack + 10140;
 Libderiv->deriv_classes[3][5][2] = int_stack + 10290;
 Libderiv->deriv_classes[3][6][2] = int_stack + 10500;
 Libderiv->deriv_classes[4][4][2] = int_stack + 10780;
 Libderiv->deriv_classes[4][5][2] = int_stack + 11005;
 Libderiv->deriv_classes[4][6][2] = int_stack + 11320;
 Libderiv->deriv_classes[3][4][1] = int_stack + 11740;
 Libderiv->deriv_classes[3][5][1] = int_stack + 11890;
 Libderiv->deriv_classes[3][6][1] = int_stack + 12100;
 Libderiv->deriv_classes[4][4][1] = int_stack + 12380;
 Libderiv->deriv_classes[4][5][1] = int_stack + 12605;
 Libderiv->deriv_classes[4][6][1] = int_stack + 12920;
 Libderiv->dvrr_classes[3][4] = int_stack + 13340;
 Libderiv->dvrr_classes[3][5] = int_stack + 13490;
 Libderiv->dvrr_classes[3][6] = int_stack + 13700;
 Libderiv->deriv_classes[3][4][0] = int_stack + 13980;
 Libderiv->deriv_classes[3][5][0] = int_stack + 14130;
 Libderiv->deriv_classes[3][6][0] = int_stack + 14340;
 Libderiv->deriv_classes[4][4][0] = int_stack + 14620;
 Libderiv->deriv_classes[4][5][0] = int_stack + 14845;
 Libderiv->deriv_classes[4][6][0] = int_stack + 15160;
 memset(int_stack,0,124640);

 Libderiv->dvrr_stack = int_stack + 34840;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpgd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15580,int_stack+13490,int_stack+13340,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16030,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13340,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16480,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13490,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17110,int_stack+16480,int_stack+16030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15580,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16030,int_stack+9090,int_stack+8640,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+865,int_stack+640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8640,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18685,int_stack+1180,int_stack+865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9090,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+18685,int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16030,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+1750,int_stack+1600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13340, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18460,int_stack+1960,int_stack+1750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13490, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19090,int_stack+18460,int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15580, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+2465,int_stack+2240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8640, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1350,int_stack+2780,int_stack+2465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9090, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19990,int_stack+1350,int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16030, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+3350,int_stack+3200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13340, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18460,int_stack+3560,int_stack+3350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13490, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1350,int_stack+18460,int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15580, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+4065,int_stack+3840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8640, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2250,int_stack+4380,int_stack+4065, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9090, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3195,int_stack+2250,int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16030, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+4950,int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18460,int_stack+5160,int_stack+4950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2250,int_stack+18460,int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+5665,int_stack+5440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4545,int_stack+5980,int_stack+5665, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21340,int_stack+4545,int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+6550,int_stack+6400, 0.0, zero_stack, 1.0, int_stack+13340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18460,int_stack+6760,int_stack+6550, 0.0, zero_stack, 1.0, int_stack+13490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4545,int_stack+18460,int_stack+18010, 0.0, zero_stack, 1.0, int_stack+15580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+7265,int_stack+7040, 0.0, zero_stack, 1.0, int_stack+8640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5445,int_stack+7580,int_stack+7265, 0.0, zero_stack, 1.0, int_stack+9090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6390,int_stack+5445,int_stack+18010, 0.0, zero_stack, 1.0, int_stack+16030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+8150,int_stack+8000, 1.0, int_stack+13340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18460,int_stack+8360,int_stack+8150, 1.0, int_stack+13490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7740,int_stack+18460,int_stack+18010, 1.0, int_stack+15580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18010,int_stack+9405,int_stack+8865, 1.0, int_stack+8640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5445,int_stack+9720,int_stack+9405, 1.0, int_stack+9090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8640,int_stack+5445,int_stack+18010, 1.0, int_stack+16030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16030,int_stack+13700,int_stack+13490,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18010,int_stack+16030,int_stack+15580,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15580,int_stack+10290,int_stack+10140,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16030,int_stack+10500,int_stack+10290,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5445,int_stack+16030,int_stack+15580,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15580,int_stack+11005,int_stack+10780,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+9990,int_stack+11320,int_stack+11005,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+22690,int_stack+9990,int_stack+15580,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15580,int_stack+11890,int_stack+11740,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16030,int_stack+12100,int_stack+11890,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9990,int_stack+16030,int_stack+15580,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15580,int_stack+12605,int_stack+12380,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10890,int_stack+12920,int_stack+12605,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11835,int_stack+10890,int_stack+15580,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15580,int_stack+14130,int_stack+13980,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16030,int_stack+14340,int_stack+14130,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10890,int_stack+16030,int_stack+15580,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15580,int_stack+14845,int_stack+14620,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+13185,int_stack+15160,int_stack+14845,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+14130,int_stack+13185,int_stack+15580,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+24040,int_stack+0,int_stack+17110,90);
     Libderiv->ABCD[11] = int_stack + 24040;
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+26740,int_stack+19990,int_stack+19090,90);
     Libderiv->ABCD[10] = int_stack + 26740;
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+29440,int_stack+3195,int_stack+1350,90);
     Libderiv->ABCD[9] = int_stack + 29440;
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+32140,int_stack+21340,int_stack+2250,90);
     Libderiv->ABCD[8] = int_stack + 32140;
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+6390,int_stack+4545,90);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+2700,int_stack+8640,int_stack+7740,90);
     Libderiv->ABCD[6] = int_stack + 2700;
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+18910,int_stack+22690,int_stack+5445, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[2] = int_stack + 18910;
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+5400,int_stack+11835,int_stack+9990, 0.0, zero_stack, 1.0, int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[1] = int_stack + 5400;
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+8100,int_stack+14130,int_stack+10890, 1.0, int_stack+18010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[0] = int_stack + 8100;

}
