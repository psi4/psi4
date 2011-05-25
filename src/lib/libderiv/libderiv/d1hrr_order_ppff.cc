#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|ff) integrals */

void d1hrr_order_ppff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[1][4][11] = int_stack + 30;
 Libderiv->deriv_classes[1][5][11] = int_stack + 75;
 Libderiv->deriv_classes[1][6][11] = int_stack + 138;
 Libderiv->deriv_classes[2][3][11] = int_stack + 222;
 Libderiv->deriv_classes[2][4][11] = int_stack + 282;
 Libderiv->deriv_classes[2][5][11] = int_stack + 372;
 Libderiv->deriv_classes[2][6][11] = int_stack + 498;
 Libderiv->deriv_classes[1][3][10] = int_stack + 666;
 Libderiv->deriv_classes[1][4][10] = int_stack + 696;
 Libderiv->deriv_classes[1][5][10] = int_stack + 741;
 Libderiv->deriv_classes[1][6][10] = int_stack + 804;
 Libderiv->deriv_classes[2][3][10] = int_stack + 888;
 Libderiv->deriv_classes[2][4][10] = int_stack + 948;
 Libderiv->deriv_classes[2][5][10] = int_stack + 1038;
 Libderiv->deriv_classes[2][6][10] = int_stack + 1164;
 Libderiv->deriv_classes[1][3][9] = int_stack + 1332;
 Libderiv->deriv_classes[1][4][9] = int_stack + 1362;
 Libderiv->deriv_classes[1][5][9] = int_stack + 1407;
 Libderiv->deriv_classes[1][6][9] = int_stack + 1470;
 Libderiv->deriv_classes[2][3][9] = int_stack + 1554;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1614;
 Libderiv->deriv_classes[2][5][9] = int_stack + 1704;
 Libderiv->deriv_classes[2][6][9] = int_stack + 1830;
 Libderiv->deriv_classes[1][3][8] = int_stack + 1998;
 Libderiv->deriv_classes[1][4][8] = int_stack + 2028;
 Libderiv->deriv_classes[1][5][8] = int_stack + 2073;
 Libderiv->deriv_classes[1][6][8] = int_stack + 2136;
 Libderiv->deriv_classes[2][3][8] = int_stack + 2220;
 Libderiv->deriv_classes[2][4][8] = int_stack + 2280;
 Libderiv->deriv_classes[2][5][8] = int_stack + 2370;
 Libderiv->deriv_classes[2][6][8] = int_stack + 2496;
 Libderiv->deriv_classes[1][3][7] = int_stack + 2664;
 Libderiv->deriv_classes[1][4][7] = int_stack + 2694;
 Libderiv->deriv_classes[1][5][7] = int_stack + 2739;
 Libderiv->deriv_classes[1][6][7] = int_stack + 2802;
 Libderiv->deriv_classes[2][3][7] = int_stack + 2886;
 Libderiv->deriv_classes[2][4][7] = int_stack + 2946;
 Libderiv->deriv_classes[2][5][7] = int_stack + 3036;
 Libderiv->deriv_classes[2][6][7] = int_stack + 3162;
 Libderiv->deriv_classes[1][3][6] = int_stack + 3330;
 Libderiv->deriv_classes[1][4][6] = int_stack + 3360;
 Libderiv->deriv_classes[1][5][6] = int_stack + 3405;
 Libderiv->deriv_classes[1][6][6] = int_stack + 3468;
 Libderiv->dvrr_classes[2][3] = int_stack + 3552;
 Libderiv->deriv_classes[2][3][6] = int_stack + 3612;
 Libderiv->dvrr_classes[2][4] = int_stack + 3672;
 Libderiv->deriv_classes[2][4][6] = int_stack + 3762;
 Libderiv->dvrr_classes[2][5] = int_stack + 3852;
 Libderiv->deriv_classes[2][5][6] = int_stack + 3978;
 Libderiv->deriv_classes[2][6][6] = int_stack + 4104;
 Libderiv->deriv_classes[1][3][2] = int_stack + 4272;
 Libderiv->deriv_classes[1][4][2] = int_stack + 4302;
 Libderiv->deriv_classes[1][5][2] = int_stack + 4347;
 Libderiv->deriv_classes[1][6][2] = int_stack + 4410;
 Libderiv->deriv_classes[2][3][2] = int_stack + 4494;
 Libderiv->deriv_classes[2][4][2] = int_stack + 4554;
 Libderiv->deriv_classes[2][5][2] = int_stack + 4644;
 Libderiv->deriv_classes[2][6][2] = int_stack + 4770;
 Libderiv->deriv_classes[1][3][1] = int_stack + 4938;
 Libderiv->deriv_classes[1][4][1] = int_stack + 4968;
 Libderiv->deriv_classes[1][5][1] = int_stack + 5013;
 Libderiv->deriv_classes[1][6][1] = int_stack + 5076;
 Libderiv->deriv_classes[2][3][1] = int_stack + 5160;
 Libderiv->deriv_classes[2][4][1] = int_stack + 5220;
 Libderiv->deriv_classes[2][5][1] = int_stack + 5310;
 Libderiv->deriv_classes[2][6][1] = int_stack + 5436;
 Libderiv->dvrr_classes[1][3] = int_stack + 5604;
 Libderiv->dvrr_classes[1][4] = int_stack + 5634;
 Libderiv->dvrr_classes[1][5] = int_stack + 5679;
 Libderiv->dvrr_classes[1][6] = int_stack + 5742;
 Libderiv->deriv_classes[1][3][0] = int_stack + 5826;
 Libderiv->deriv_classes[1][4][0] = int_stack + 5856;
 Libderiv->deriv_classes[1][5][0] = int_stack + 5901;
 Libderiv->deriv_classes[1][6][0] = int_stack + 5964;
 Libderiv->deriv_classes[2][3][0] = int_stack + 6048;
 Libderiv->deriv_classes[2][4][0] = int_stack + 6108;
 Libderiv->deriv_classes[2][5][0] = int_stack + 6198;
 Libderiv->deriv_classes[2][6][0] = int_stack + 6324;
 memset(int_stack,0,51936);

 Libderiv->dvrr_stack = int_stack + 12768;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppff(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6492,int_stack+5634,int_stack+5604,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6582,int_stack+5679,int_stack+5634,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+6717,int_stack+6582,int_stack+6492,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6897,int_stack+30,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5604,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6987,int_stack+75,int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5634,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7122,int_stack+6987,int_stack+6897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6492,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7302,int_stack+138,int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5679,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7491,int_stack+7302,int_stack+6987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6582,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7761,int_stack+7491,int_stack+7122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6717,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6897,int_stack+3672,int_stack+3552,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7077,int_stack+3852,int_stack+3672,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7347,int_stack+7077,int_stack+6897,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+282,int_stack+222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3552,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8061,int_stack+372,int_stack+282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3672,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8331,int_stack+8061,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6897,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8691,int_stack+498,int_stack+372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3852,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+8691,int_stack+8061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7077,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8691,int_stack+0,int_stack+8331, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7347,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+696,int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5604, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+741,int_stack+696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5634, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+225,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6492, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+405,int_stack+804,int_stack+741, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5679, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+594,int_stack+405,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6582, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8061,int_stack+594,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6717, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+948,int_stack+888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3552, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180,int_stack+1038,int_stack+948, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3672, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+450,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6897, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9291,int_stack+1164,int_stack+1038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3852, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9669,int_stack+9291,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7077, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10209,int_stack+9669,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7347, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9291,int_stack+1362,int_stack+1332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5604, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9381,int_stack+1407,int_stack+1362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5634, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9516,int_stack+9381,int_stack+9291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6492, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9696,int_stack+1470,int_stack+1407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5679, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9885,int_stack+9696,int_stack+9381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6582, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+9885,int_stack+9516, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6717, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9291,int_stack+1614,int_stack+1554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3552, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9471,int_stack+1704,int_stack+1614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3672, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9741,int_stack+9471,int_stack+9291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6897, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+300,int_stack+1830,int_stack+1704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3852, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+678,int_stack+300,int_stack+9471, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7077, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1218,int_stack+678,int_stack+9741, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7347, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+2028,int_stack+1998, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+390,int_stack+2073,int_stack+2028, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+525,int_stack+390,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+705,int_stack+2136,int_stack+2073, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5679, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+894,int_stack+705,int_stack+390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9291,int_stack+894,int_stack+525, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6717, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+2280,int_stack+2220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+480,int_stack+2370,int_stack+2280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+750,int_stack+480,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9591,int_stack+2496,int_stack+2370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1818,int_stack+9591,int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7077, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9591,int_stack+1818,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1818,int_stack+2694,int_stack+2664, 0.0, zero_stack, 1.0, int_stack+5604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1908,int_stack+2739,int_stack+2694, 0.0, zero_stack, 1.0, int_stack+5634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2043,int_stack+1908,int_stack+1818, 0.0, zero_stack, 1.0, int_stack+6492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2223,int_stack+2802,int_stack+2739, 0.0, zero_stack, 1.0, int_stack+5679, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2412,int_stack+2223,int_stack+1908, 0.0, zero_stack, 1.0, int_stack+6582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+300,int_stack+2412,int_stack+2043, 0.0, zero_stack, 1.0, int_stack+6717, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1818,int_stack+2946,int_stack+2886, 0.0, zero_stack, 1.0, int_stack+3552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1998,int_stack+3036,int_stack+2946, 0.0, zero_stack, 1.0, int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2268,int_stack+1998,int_stack+1818, 0.0, zero_stack, 1.0, int_stack+6897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2628,int_stack+3162,int_stack+3036, 0.0, zero_stack, 1.0, int_stack+3852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+600,int_stack+2628,int_stack+1998, 0.0, zero_stack, 1.0, int_stack+7077, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2628,int_stack+600,int_stack+2268, 0.0, zero_stack, 1.0, int_stack+7347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+3360,int_stack+3330, 1.0, int_stack+5604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+690,int_stack+3405,int_stack+3360, 1.0, int_stack+5634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+825,int_stack+690,int_stack+600, 1.0, int_stack+6492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1005,int_stack+3468,int_stack+3405, 1.0, int_stack+5679, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3228,int_stack+1005,int_stack+690, 1.0, int_stack+6582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1818,int_stack+3228,int_stack+825, 1.0, int_stack+6717, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3228,int_stack+3762,int_stack+3612, 1.0, int_stack+3552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+600,int_stack+3978,int_stack+3762, 1.0, int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3408,int_stack+600,int_stack+3228, 1.0, int_stack+6897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2118,int_stack+4104,int_stack+3978, 1.0, int_stack+3852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10809,int_stack+2118,int_stack+600, 1.0, int_stack+7077, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+600,int_stack+10809,int_stack+3408, 1.0, int_stack+7347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10809,int_stack+5742,int_stack+5679,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10998,int_stack+10809,int_stack+6582,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+11268,int_stack+10998,int_stack+6717,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10809,int_stack+4302,int_stack+4272,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10899,int_stack+4347,int_stack+4302,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11034,int_stack+10899,int_stack+10809,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11568,int_stack+4410,int_stack+4347,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2118,int_stack+11568,int_stack+10899,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+11568,int_stack+2118,int_stack+11034,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2118,int_stack+4554,int_stack+4494,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2298,int_stack+4644,int_stack+4554,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10809,int_stack+2298,int_stack+2118,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3228,int_stack+4770,int_stack+4644,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3606,int_stack+3228,int_stack+2298,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+4146,int_stack+3606,int_stack+10809,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10809,int_stack+4968,int_stack+4938,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10899,int_stack+5013,int_stack+4968,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11034,int_stack+10899,int_stack+10809,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3228,int_stack+5076,int_stack+5013,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3417,int_stack+3228,int_stack+10899,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+3687,int_stack+3417,int_stack+11034,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3228,int_stack+5220,int_stack+5160,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3408,int_stack+5310,int_stack+5220,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10809,int_stack+3408,int_stack+3228,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2118,int_stack+5436,int_stack+5310,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4746,int_stack+2118,int_stack+3408,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+6492,int_stack+4746,int_stack+10809,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10809,int_stack+5856,int_stack+5826,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10899,int_stack+5901,int_stack+5856,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11034,int_stack+10899,int_stack+10809,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4746,int_stack+5964,int_stack+5901,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4935,int_stack+4746,int_stack+10899,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+5205,int_stack+4935,int_stack+11034,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4746,int_stack+6108,int_stack+6048,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4926,int_stack+6198,int_stack+6108,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10809,int_stack+4926,int_stack+4746,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5505,int_stack+6324,int_stack+6198,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5883,int_stack+5505,int_stack+4926,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+7092,int_stack+5883,int_stack+10809,6);
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+5505,int_stack+8691,int_stack+7761,100);
     Libderiv->ABCD[11] = int_stack + 5505;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+8361,int_stack+10209,int_stack+8061,100);
     Libderiv->ABCD[10] = int_stack + 8361;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+10191,int_stack+1218,int_stack+0,100);
     Libderiv->ABCD[9] = int_stack + 10191;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+11868,int_stack+9591,int_stack+9291,100);
     Libderiv->ABCD[8] = int_stack + 11868;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+9261,int_stack+2628,int_stack+300,100);
     Libderiv->ABCD[7] = int_stack + 9261;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2118,int_stack+600,int_stack+1818,100);
     Libderiv->ABCD[6] = int_stack + 2118;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+4146,int_stack+11568, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 0;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+900,int_stack+6492,int_stack+3687, 0.0, zero_stack, 1.0, int_stack+11268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 900;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3018,int_stack+7092,int_stack+5205, 1.0, int_stack+11268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 3018;

}
