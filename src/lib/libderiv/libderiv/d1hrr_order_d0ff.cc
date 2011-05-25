#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|ff) integrals */

void d1hrr_order_d0ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[2][5][11] = int_stack + 150;
 Libderiv->deriv_classes[2][6][11] = int_stack + 276;
 Libderiv->deriv_classes[2][3][10] = int_stack + 444;
 Libderiv->deriv_classes[2][4][10] = int_stack + 504;
 Libderiv->deriv_classes[2][5][10] = int_stack + 594;
 Libderiv->deriv_classes[2][6][10] = int_stack + 720;
 Libderiv->deriv_classes[2][3][9] = int_stack + 888;
 Libderiv->deriv_classes[2][4][9] = int_stack + 948;
 Libderiv->deriv_classes[2][5][9] = int_stack + 1038;
 Libderiv->deriv_classes[2][6][9] = int_stack + 1164;
 Libderiv->deriv_classes[2][3][8] = int_stack + 1332;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1392;
 Libderiv->deriv_classes[2][5][8] = int_stack + 1482;
 Libderiv->deriv_classes[2][6][8] = int_stack + 1608;
 Libderiv->deriv_classes[2][3][7] = int_stack + 1776;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1836;
 Libderiv->deriv_classes[2][5][7] = int_stack + 1926;
 Libderiv->deriv_classes[2][6][7] = int_stack + 2052;
 Libderiv->dvrr_classes[2][3] = int_stack + 2220;
 Libderiv->deriv_classes[2][3][6] = int_stack + 2280;
 Libderiv->dvrr_classes[2][4] = int_stack + 2340;
 Libderiv->deriv_classes[2][4][6] = int_stack + 2430;
 Libderiv->dvrr_classes[2][5] = int_stack + 2520;
 Libderiv->deriv_classes[2][5][6] = int_stack + 2646;
 Libderiv->deriv_classes[2][6][6] = int_stack + 2772;
 Libderiv->deriv_classes[2][3][2] = int_stack + 2940;
 Libderiv->deriv_classes[2][4][2] = int_stack + 3000;
 Libderiv->deriv_classes[2][5][2] = int_stack + 3090;
 Libderiv->deriv_classes[2][6][2] = int_stack + 3216;
 Libderiv->deriv_classes[2][3][1] = int_stack + 3384;
 Libderiv->deriv_classes[2][4][1] = int_stack + 3444;
 Libderiv->deriv_classes[2][5][1] = int_stack + 3534;
 Libderiv->deriv_classes[2][6][1] = int_stack + 3660;
 Libderiv->deriv_classes[2][3][0] = int_stack + 3828;
 Libderiv->deriv_classes[2][4][0] = int_stack + 3888;
 Libderiv->deriv_classes[2][5][0] = int_stack + 3978;
 Libderiv->deriv_classes[2][6][0] = int_stack + 4104;
 memset(int_stack,0,34176);

 Libderiv->dvrr_stack = int_stack + 10668;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4272,int_stack+2340,int_stack+2220,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4452,int_stack+2520,int_stack+2340,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4722,int_stack+4452,int_stack+4272,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5082,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2220,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5262,int_stack+150,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2340,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5532,int_stack+5262,int_stack+5082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4272,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5892,int_stack+276,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2520,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6270,int_stack+5892,int_stack+5262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4452,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5892,int_stack+504,int_stack+444, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2220, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5082,int_stack+594,int_stack+504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2340, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+5082,int_stack+5892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4272, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5892,int_stack+720,int_stack+594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2520, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6810,int_stack+5892,int_stack+5082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4452, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5082,int_stack+948,int_stack+888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2220, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5262,int_stack+1038,int_stack+948, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2340, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5892,int_stack+5262,int_stack+5082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4272, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+360,int_stack+1164,int_stack+1038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2520, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+738,int_stack+360,int_stack+5262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4452, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+1392,int_stack+1332, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5082,int_stack+1482,int_stack+1392, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7350,int_stack+5082,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+360,int_stack+1608,int_stack+1482, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7710,int_stack+360,int_stack+5082, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5082,int_stack+1836,int_stack+1776, 0.0, zero_stack, 1.0, int_stack+2220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5262,int_stack+1926,int_stack+1836, 0.0, zero_stack, 1.0, int_stack+2340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+360,int_stack+5262,int_stack+5082, 0.0, zero_stack, 1.0, int_stack+4272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1278,int_stack+2052,int_stack+1926, 0.0, zero_stack, 1.0, int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1656,int_stack+1278,int_stack+5262, 0.0, zero_stack, 1.0, int_stack+4452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1278,int_stack+2430,int_stack+2280, 1.0, int_stack+2220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5082,int_stack+2646,int_stack+2430, 1.0, int_stack+2340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8250,int_stack+5082,int_stack+1278, 1.0, int_stack+4272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1278,int_stack+2772,int_stack+2646, 1.0, int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2196,int_stack+1278,int_stack+5082, 1.0, int_stack+4452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5082,int_stack+3000,int_stack+2940,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5262,int_stack+3090,int_stack+3000,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1278,int_stack+5262,int_stack+5082,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4272,int_stack+3216,int_stack+3090,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2736,int_stack+4272,int_stack+5262,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4272,int_stack+3444,int_stack+3384,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4452,int_stack+3534,int_stack+3444,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5082,int_stack+4452,int_stack+4272,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8610,int_stack+3660,int_stack+3534,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3276,int_stack+8610,int_stack+4452,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8610,int_stack+3888,int_stack+3828,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4272,int_stack+3978,int_stack+3888,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8790,int_stack+4272,int_stack+8610,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+9150,int_stack+4104,int_stack+3978,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9528,int_stack+9150,int_stack+4272,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10068,int_stack+6270,int_stack+5532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4722,6);
     Libderiv->ABCD[11] = int_stack + 10068;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3816,int_stack+6810,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4722, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 3816;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6252,int_stack+738,int_stack+5892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4722, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 6252;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5442,int_stack+7710,int_stack+7350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 5442;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6852,int_stack+1656,int_stack+360, 0.0, zero_stack, 1.0, int_stack+4722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 6852;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+2196,int_stack+8250, 1.0, int_stack+4722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 0;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+600,int_stack+2736,int_stack+1278,6);
     Libderiv->ABCD[2] = int_stack + 600;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1200,int_stack+3276,int_stack+5082,6);
     Libderiv->ABCD[1] = int_stack + 1200;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1800,int_stack+9528,int_stack+8790,6);
     Libderiv->ABCD[0] = int_stack + 1800;

}
