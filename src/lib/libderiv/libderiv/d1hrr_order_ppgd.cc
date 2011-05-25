#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppgd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|gd) integrals */

void d1hrr_order_ppgd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[1][5][11] = int_stack + 45;
 Libderiv->deriv_classes[1][6][11] = int_stack + 108;
 Libderiv->deriv_classes[2][4][11] = int_stack + 192;
 Libderiv->deriv_classes[2][5][11] = int_stack + 282;
 Libderiv->deriv_classes[2][6][11] = int_stack + 408;
 Libderiv->deriv_classes[1][4][10] = int_stack + 576;
 Libderiv->deriv_classes[1][5][10] = int_stack + 621;
 Libderiv->deriv_classes[1][6][10] = int_stack + 684;
 Libderiv->deriv_classes[2][4][10] = int_stack + 768;
 Libderiv->deriv_classes[2][5][10] = int_stack + 858;
 Libderiv->deriv_classes[2][6][10] = int_stack + 984;
 Libderiv->deriv_classes[1][4][9] = int_stack + 1152;
 Libderiv->deriv_classes[1][5][9] = int_stack + 1197;
 Libderiv->deriv_classes[1][6][9] = int_stack + 1260;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1344;
 Libderiv->deriv_classes[2][5][9] = int_stack + 1434;
 Libderiv->deriv_classes[2][6][9] = int_stack + 1560;
 Libderiv->deriv_classes[1][4][8] = int_stack + 1728;
 Libderiv->deriv_classes[1][5][8] = int_stack + 1773;
 Libderiv->deriv_classes[1][6][8] = int_stack + 1836;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1920;
 Libderiv->deriv_classes[2][5][8] = int_stack + 2010;
 Libderiv->deriv_classes[2][6][8] = int_stack + 2136;
 Libderiv->deriv_classes[1][4][7] = int_stack + 2304;
 Libderiv->deriv_classes[1][5][7] = int_stack + 2349;
 Libderiv->deriv_classes[1][6][7] = int_stack + 2412;
 Libderiv->deriv_classes[2][4][7] = int_stack + 2496;
 Libderiv->deriv_classes[2][5][7] = int_stack + 2586;
 Libderiv->deriv_classes[2][6][7] = int_stack + 2712;
 Libderiv->deriv_classes[1][4][6] = int_stack + 2880;
 Libderiv->deriv_classes[1][5][6] = int_stack + 2925;
 Libderiv->deriv_classes[1][6][6] = int_stack + 2988;
 Libderiv->dvrr_classes[2][4] = int_stack + 3072;
 Libderiv->deriv_classes[2][4][6] = int_stack + 3162;
 Libderiv->dvrr_classes[2][5] = int_stack + 3252;
 Libderiv->deriv_classes[2][5][6] = int_stack + 3378;
 Libderiv->deriv_classes[2][6][6] = int_stack + 3504;
 Libderiv->deriv_classes[1][4][2] = int_stack + 3672;
 Libderiv->deriv_classes[1][5][2] = int_stack + 3717;
 Libderiv->deriv_classes[1][6][2] = int_stack + 3780;
 Libderiv->deriv_classes[2][4][2] = int_stack + 3864;
 Libderiv->deriv_classes[2][5][2] = int_stack + 3954;
 Libderiv->deriv_classes[2][6][2] = int_stack + 4080;
 Libderiv->deriv_classes[1][4][1] = int_stack + 4248;
 Libderiv->deriv_classes[1][5][1] = int_stack + 4293;
 Libderiv->deriv_classes[1][6][1] = int_stack + 4356;
 Libderiv->deriv_classes[2][4][1] = int_stack + 4440;
 Libderiv->deriv_classes[2][5][1] = int_stack + 4530;
 Libderiv->deriv_classes[2][6][1] = int_stack + 4656;
 Libderiv->dvrr_classes[1][4] = int_stack + 4824;
 Libderiv->dvrr_classes[1][5] = int_stack + 4869;
 Libderiv->dvrr_classes[1][6] = int_stack + 4932;
 Libderiv->deriv_classes[1][4][0] = int_stack + 5016;
 Libderiv->deriv_classes[1][5][0] = int_stack + 5061;
 Libderiv->deriv_classes[1][6][0] = int_stack + 5124;
 Libderiv->deriv_classes[2][4][0] = int_stack + 5208;
 Libderiv->deriv_classes[2][5][0] = int_stack + 5298;
 Libderiv->deriv_classes[2][6][0] = int_stack + 5424;
 memset(int_stack,0,44736);

 Libderiv->dvrr_stack = int_stack + 11073;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppgd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5592,int_stack+4869,int_stack+4824,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5727,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4824,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5862,int_stack+108,int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4869,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6051,int_stack+5862,int_stack+5727, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5592,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5727,int_stack+3252,int_stack+3072,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+282,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3072,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6591,int_stack+408,int_stack+282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3252,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+6591,int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5727,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+621,int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4824, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6456,int_stack+684,int_stack+621, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4869, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6645,int_stack+6456,int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5592, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+858,int_stack+768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3072, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6915,int_stack+984,int_stack+858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3252, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7293,int_stack+6915,int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5727, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+1197,int_stack+1152, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4824, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6456,int_stack+1260,int_stack+1197, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4869, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6915,int_stack+6456,int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5592, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+1434,int_stack+1344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3072, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+540,int_stack+1560,int_stack+1434, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3252, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+918,int_stack+540,int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5727, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+1773,int_stack+1728, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6456,int_stack+1836,int_stack+1773, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4869, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+540,int_stack+6456,int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+2010,int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1458,int_stack+2136,int_stack+2010, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7833,int_stack+1458,int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5727, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+2349,int_stack+2304, 0.0, zero_stack, 1.0, int_stack+4824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6456,int_stack+2412,int_stack+2349, 0.0, zero_stack, 1.0, int_stack+4869, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1458,int_stack+6456,int_stack+6321, 0.0, zero_stack, 1.0, int_stack+5592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+2586,int_stack+2496, 0.0, zero_stack, 1.0, int_stack+3072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1728,int_stack+2712,int_stack+2586, 0.0, zero_stack, 1.0, int_stack+3252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2106,int_stack+1728,int_stack+6321, 0.0, zero_stack, 1.0, int_stack+5727, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+2925,int_stack+2880, 1.0, int_stack+4824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6456,int_stack+2988,int_stack+2925, 1.0, int_stack+4869, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1728,int_stack+6456,int_stack+6321, 1.0, int_stack+5592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6321,int_stack+3378,int_stack+3162, 1.0, int_stack+3072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2646,int_stack+3504,int_stack+3378, 1.0, int_stack+3252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3024,int_stack+2646,int_stack+6321, 1.0, int_stack+5727, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5727,int_stack+4932,int_stack+4869,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6321,int_stack+5727,int_stack+5592,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5592,int_stack+3717,int_stack+3672,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5727,int_stack+3780,int_stack+3717,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2646,int_stack+5727,int_stack+5592,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5592,int_stack+3954,int_stack+3864,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3564,int_stack+4080,int_stack+3954,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8373,int_stack+3564,int_stack+5592,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5592,int_stack+4293,int_stack+4248,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5727,int_stack+4356,int_stack+4293,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3564,int_stack+5727,int_stack+5592,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5592,int_stack+4530,int_stack+4440,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3834,int_stack+4656,int_stack+4530,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4212,int_stack+3834,int_stack+5592,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5592,int_stack+5061,int_stack+5016,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5727,int_stack+5124,int_stack+5061,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3834,int_stack+5727,int_stack+5592,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5592,int_stack+5298,int_stack+5208,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4752,int_stack+5424,int_stack+5298,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8913,int_stack+4752,int_stack+5592,6);
 /*--- compute (pp|gd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4752,int_stack+0,int_stack+6051,90);
     Libderiv->ABCD[11] = int_stack + 4752;
 /*--- compute (pp|gd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+9453,int_stack+7293,int_stack+6645,90);
     Libderiv->ABCD[10] = int_stack + 9453;
 /*--- compute (pp|gd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+10263,int_stack+918,int_stack+6915,90);
     Libderiv->ABCD[9] = int_stack + 10263;
 /*--- compute (pp|gd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+6591,int_stack+7833,int_stack+540,90);
     Libderiv->ABCD[8] = int_stack + 6591;
 /*--- compute (pp|gd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+2106,int_stack+1458,90);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (pp|gd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+810,int_stack+3024,int_stack+1728,90);
     Libderiv->ABCD[6] = int_stack + 810;
 /*--- compute (pp|gd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1620,int_stack+8373,int_stack+2646, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[2] = int_stack + 1620;
 /*--- compute (pp|gd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2430,int_stack+4212,int_stack+3564, 0.0, zero_stack, 1.0, int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[1] = int_stack + 2430;
 /*--- compute (pp|gd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7401,int_stack+8913,int_stack+3834, 1.0, int_stack+6321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[0] = int_stack + 7401;

}
