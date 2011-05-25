#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpdd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|dd) integrals */

void d1hrr_order_dpdd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][2][11] = int_stack + 0;
 Libderiv->deriv_classes[2][3][11] = int_stack + 36;
 Libderiv->deriv_classes[2][4][11] = int_stack + 96;
 Libderiv->deriv_classes[3][2][11] = int_stack + 186;
 Libderiv->deriv_classes[3][3][11] = int_stack + 246;
 Libderiv->deriv_classes[3][4][11] = int_stack + 346;
 Libderiv->deriv_classes[2][2][10] = int_stack + 496;
 Libderiv->deriv_classes[2][3][10] = int_stack + 532;
 Libderiv->deriv_classes[2][4][10] = int_stack + 592;
 Libderiv->deriv_classes[3][2][10] = int_stack + 682;
 Libderiv->deriv_classes[3][3][10] = int_stack + 742;
 Libderiv->deriv_classes[3][4][10] = int_stack + 842;
 Libderiv->deriv_classes[2][2][9] = int_stack + 992;
 Libderiv->deriv_classes[2][3][9] = int_stack + 1028;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1088;
 Libderiv->deriv_classes[3][2][9] = int_stack + 1178;
 Libderiv->deriv_classes[3][3][9] = int_stack + 1238;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1338;
 Libderiv->deriv_classes[2][2][8] = int_stack + 1488;
 Libderiv->deriv_classes[2][3][8] = int_stack + 1524;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1584;
 Libderiv->deriv_classes[3][2][8] = int_stack + 1674;
 Libderiv->deriv_classes[3][3][8] = int_stack + 1734;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1834;
 Libderiv->deriv_classes[2][2][7] = int_stack + 1984;
 Libderiv->deriv_classes[2][3][7] = int_stack + 2020;
 Libderiv->deriv_classes[2][4][7] = int_stack + 2080;
 Libderiv->deriv_classes[3][2][7] = int_stack + 2170;
 Libderiv->deriv_classes[3][3][7] = int_stack + 2230;
 Libderiv->deriv_classes[3][4][7] = int_stack + 2330;
 Libderiv->deriv_classes[2][2][6] = int_stack + 2480;
 Libderiv->deriv_classes[2][3][6] = int_stack + 2516;
 Libderiv->deriv_classes[2][4][6] = int_stack + 2576;
 Libderiv->dvrr_classes[3][2] = int_stack + 2666;
 Libderiv->deriv_classes[3][2][6] = int_stack + 2726;
 Libderiv->dvrr_classes[3][3] = int_stack + 2786;
 Libderiv->deriv_classes[3][3][6] = int_stack + 2886;
 Libderiv->deriv_classes[3][4][6] = int_stack + 2986;
 Libderiv->deriv_classes[2][2][2] = int_stack + 3136;
 Libderiv->deriv_classes[2][3][2] = int_stack + 3172;
 Libderiv->deriv_classes[2][4][2] = int_stack + 3232;
 Libderiv->deriv_classes[3][2][2] = int_stack + 3322;
 Libderiv->deriv_classes[3][3][2] = int_stack + 3382;
 Libderiv->deriv_classes[3][4][2] = int_stack + 3482;
 Libderiv->deriv_classes[2][2][1] = int_stack + 3632;
 Libderiv->deriv_classes[2][3][1] = int_stack + 3668;
 Libderiv->deriv_classes[2][4][1] = int_stack + 3728;
 Libderiv->deriv_classes[3][2][1] = int_stack + 3818;
 Libderiv->deriv_classes[3][3][1] = int_stack + 3878;
 Libderiv->deriv_classes[3][4][1] = int_stack + 3978;
 Libderiv->dvrr_classes[2][2] = int_stack + 4128;
 Libderiv->dvrr_classes[2][3] = int_stack + 4164;
 Libderiv->dvrr_classes[2][4] = int_stack + 4224;
 Libderiv->deriv_classes[2][2][0] = int_stack + 4314;
 Libderiv->deriv_classes[2][3][0] = int_stack + 4350;
 Libderiv->deriv_classes[2][4][0] = int_stack + 4410;
 Libderiv->deriv_classes[3][2][0] = int_stack + 4500;
 Libderiv->deriv_classes[3][3][0] = int_stack + 4560;
 Libderiv->deriv_classes[3][4][0] = int_stack + 4660;
 memset(int_stack,0,38480);

 Libderiv->dvrr_stack = int_stack + 7954;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpdd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4810,int_stack+4164,int_stack+4128,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4918,int_stack+36,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4128,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5026,int_stack+96,int_stack+36, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4164,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5206,int_stack+5026,int_stack+4918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4810,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4918,int_stack+2786,int_stack+2666,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+246,int_stack+186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2666,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5422,int_stack+346,int_stack+246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2786,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5722,int_stack+5422,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4918,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5098,int_stack+532,int_stack+496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4128, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+592,int_stack+532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4164, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+180,int_stack+0,int_stack+5098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4810, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+742,int_stack+682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2666, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5422,int_stack+842,int_stack+742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2786, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+396,int_stack+5422,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4918, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5098,int_stack+1028,int_stack+992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4128, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+1088,int_stack+1028, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4164, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5422,int_stack+0,int_stack+5098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4810, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+1238,int_stack+1178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2666, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+756,int_stack+1338,int_stack+1238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2786, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1056,int_stack+756,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4918, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5098,int_stack+1524,int_stack+1488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+1584,int_stack+1524, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4164, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+756,int_stack+0,int_stack+5098, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+1734,int_stack+1674, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1416,int_stack+1834,int_stack+1734, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6082,int_stack+1416,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5098,int_stack+2020,int_stack+1984, 0.0, zero_stack, 1.0, int_stack+4128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+2080,int_stack+2020, 0.0, zero_stack, 1.0, int_stack+4164, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1416,int_stack+0,int_stack+5098, 0.0, zero_stack, 1.0, int_stack+4810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+2230,int_stack+2170, 0.0, zero_stack, 1.0, int_stack+2666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1632,int_stack+2330,int_stack+2230, 0.0, zero_stack, 1.0, int_stack+2786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1932,int_stack+1632,int_stack+0, 0.0, zero_stack, 1.0, int_stack+4918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5098,int_stack+2516,int_stack+2480, 1.0, int_stack+4128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+2576,int_stack+2516, 1.0, int_stack+4164, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1632,int_stack+0,int_stack+5098, 1.0, int_stack+4810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+2886,int_stack+2726, 1.0, int_stack+2666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2292,int_stack+2986,int_stack+2886, 1.0, int_stack+2786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2592,int_stack+2292,int_stack+0, 1.0, int_stack+4918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+4224,int_stack+4164,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+4918,int_stack+0,int_stack+4810,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4810,int_stack+3172,int_stack+3136,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+3232,int_stack+3172,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2292,int_stack+0,int_stack+4810,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+3382,int_stack+3322,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2952,int_stack+3482,int_stack+3382,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+3252,int_stack+2952,int_stack+0,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4810,int_stack+3668,int_stack+3632,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+3728,int_stack+3668,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2952,int_stack+0,int_stack+4810,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+3878,int_stack+3818,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6442,int_stack+3978,int_stack+3878,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+3612,int_stack+6442,int_stack+0,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4810,int_stack+4350,int_stack+4314,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+4410,int_stack+4350,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6442,int_stack+0,int_stack+4810,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+4560,int_stack+4500,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6658,int_stack+4660,int_stack+4560,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+3972,int_stack+6658,int_stack+0,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6658,int_stack+5722,int_stack+5206,36);
     Libderiv->ABCD[11] = int_stack + 6658;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7306,int_stack+396,int_stack+180,36);
     Libderiv->ABCD[10] = int_stack + 7306;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+1056,int_stack+5422,36);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5134,int_stack+6082,int_stack+756,36);
     Libderiv->ABCD[8] = int_stack + 5134;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+648,int_stack+1932,int_stack+1416,36);
     Libderiv->ABCD[7] = int_stack + 648;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5782,int_stack+2592,int_stack+1632,36);
     Libderiv->ABCD[6] = int_stack + 5782;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1296,int_stack+3252,int_stack+2292, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[2] = int_stack + 1296;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1944,int_stack+3612,int_stack+2952, 0.0, zero_stack, 1.0, int_stack+4918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[1] = int_stack + 1944;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+2592,int_stack+3972,int_stack+6442, 1.0, int_stack+4918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[0] = int_stack + 2592;

}
