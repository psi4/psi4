#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dddd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|dd) integrals */

void d1hrr_order_dddd(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[4][2][11] = int_stack + 496;
 Libderiv->deriv_classes[4][3][11] = int_stack + 586;
 Libderiv->deriv_classes[4][4][11] = int_stack + 736;
 Libderiv->deriv_classes[2][2][10] = int_stack + 961;
 Libderiv->deriv_classes[2][3][10] = int_stack + 997;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1057;
 Libderiv->deriv_classes[3][2][10] = int_stack + 1147;
 Libderiv->deriv_classes[3][3][10] = int_stack + 1207;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1307;
 Libderiv->deriv_classes[4][2][10] = int_stack + 1457;
 Libderiv->deriv_classes[4][3][10] = int_stack + 1547;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1697;
 Libderiv->deriv_classes[2][2][9] = int_stack + 1922;
 Libderiv->deriv_classes[2][3][9] = int_stack + 1958;
 Libderiv->deriv_classes[2][4][9] = int_stack + 2018;
 Libderiv->deriv_classes[3][2][9] = int_stack + 2108;
 Libderiv->deriv_classes[3][3][9] = int_stack + 2168;
 Libderiv->deriv_classes[3][4][9] = int_stack + 2268;
 Libderiv->deriv_classes[4][2][9] = int_stack + 2418;
 Libderiv->deriv_classes[4][3][9] = int_stack + 2508;
 Libderiv->deriv_classes[4][4][9] = int_stack + 2658;
 Libderiv->deriv_classes[2][2][8] = int_stack + 2883;
 Libderiv->deriv_classes[2][3][8] = int_stack + 2919;
 Libderiv->deriv_classes[2][4][8] = int_stack + 2979;
 Libderiv->deriv_classes[3][2][8] = int_stack + 3069;
 Libderiv->deriv_classes[3][3][8] = int_stack + 3129;
 Libderiv->deriv_classes[3][4][8] = int_stack + 3229;
 Libderiv->deriv_classes[4][2][8] = int_stack + 3379;
 Libderiv->deriv_classes[4][3][8] = int_stack + 3469;
 Libderiv->deriv_classes[4][4][8] = int_stack + 3619;
 Libderiv->deriv_classes[2][2][7] = int_stack + 3844;
 Libderiv->deriv_classes[2][3][7] = int_stack + 3880;
 Libderiv->deriv_classes[2][4][7] = int_stack + 3940;
 Libderiv->deriv_classes[3][2][7] = int_stack + 4030;
 Libderiv->deriv_classes[3][3][7] = int_stack + 4090;
 Libderiv->deriv_classes[3][4][7] = int_stack + 4190;
 Libderiv->deriv_classes[4][2][7] = int_stack + 4340;
 Libderiv->deriv_classes[4][3][7] = int_stack + 4430;
 Libderiv->deriv_classes[4][4][7] = int_stack + 4580;
 Libderiv->deriv_classes[2][2][6] = int_stack + 4805;
 Libderiv->deriv_classes[2][3][6] = int_stack + 4841;
 Libderiv->deriv_classes[2][4][6] = int_stack + 4901;
 Libderiv->deriv_classes[3][2][6] = int_stack + 4991;
 Libderiv->deriv_classes[3][3][6] = int_stack + 5051;
 Libderiv->deriv_classes[3][4][6] = int_stack + 5151;
 Libderiv->dvrr_classes[4][2] = int_stack + 5301;
 Libderiv->deriv_classes[4][2][6] = int_stack + 5391;
 Libderiv->dvrr_classes[4][3] = int_stack + 5481;
 Libderiv->deriv_classes[4][3][6] = int_stack + 5631;
 Libderiv->deriv_classes[4][4][6] = int_stack + 5781;
 Libderiv->deriv_classes[2][2][2] = int_stack + 6006;
 Libderiv->deriv_classes[2][3][2] = int_stack + 6042;
 Libderiv->deriv_classes[2][4][2] = int_stack + 6102;
 Libderiv->deriv_classes[3][2][2] = int_stack + 6192;
 Libderiv->deriv_classes[3][3][2] = int_stack + 6252;
 Libderiv->deriv_classes[3][4][2] = int_stack + 6352;
 Libderiv->deriv_classes[4][2][2] = int_stack + 6502;
 Libderiv->deriv_classes[4][3][2] = int_stack + 6592;
 Libderiv->deriv_classes[4][4][2] = int_stack + 6742;
 Libderiv->deriv_classes[2][2][1] = int_stack + 6967;
 Libderiv->deriv_classes[2][3][1] = int_stack + 7003;
 Libderiv->deriv_classes[2][4][1] = int_stack + 7063;
 Libderiv->deriv_classes[3][2][1] = int_stack + 7153;
 Libderiv->deriv_classes[3][3][1] = int_stack + 7213;
 Libderiv->deriv_classes[3][4][1] = int_stack + 7313;
 Libderiv->deriv_classes[4][2][1] = int_stack + 7463;
 Libderiv->deriv_classes[4][3][1] = int_stack + 7553;
 Libderiv->deriv_classes[4][4][1] = int_stack + 7703;
 Libderiv->dvrr_classes[2][2] = int_stack + 7928;
 Libderiv->dvrr_classes[2][3] = int_stack + 7964;
 Libderiv->dvrr_classes[2][4] = int_stack + 8024;
 Libderiv->deriv_classes[2][2][0] = int_stack + 8114;
 Libderiv->deriv_classes[2][3][0] = int_stack + 8150;
 Libderiv->deriv_classes[2][4][0] = int_stack + 8210;
 Libderiv->dvrr_classes[3][2] = int_stack + 8300;
 Libderiv->dvrr_classes[3][3] = int_stack + 8360;
 Libderiv->dvrr_classes[3][4] = int_stack + 8460;
 Libderiv->deriv_classes[3][2][0] = int_stack + 8610;
 Libderiv->deriv_classes[3][3][0] = int_stack + 8670;
 Libderiv->deriv_classes[3][4][0] = int_stack + 8770;
 Libderiv->deriv_classes[4][2][0] = int_stack + 8920;
 Libderiv->deriv_classes[4][3][0] = int_stack + 9010;
 Libderiv->deriv_classes[4][4][0] = int_stack + 9160;
 memset(int_stack,0,75080);

 Libderiv->dvrr_stack = int_stack + 20395;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dddd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+9385,int_stack+7964,int_stack+7928,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9493,int_stack+36,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7928,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9601,int_stack+96,int_stack+36, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7964,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9781,int_stack+9601,int_stack+9493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9385,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+9493,int_stack+8360,int_stack+8300,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+246,int_stack+186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8300,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9997,int_stack+346,int_stack+246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8360,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10297,int_stack+9997,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9493,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+10657,int_stack+10297,int_stack+9781,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+5481,int_stack+5301,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9673,int_stack+586,int_stack+496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5301,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11305,int_stack+736,int_stack+586, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5481,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+270,int_stack+11305,int_stack+9673, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11305,int_stack+270,int_stack+10297,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+997,int_stack+961, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7928, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+378,int_stack+1057,int_stack+997, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7964, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+558,int_stack+378,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9385, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+1207,int_stack+1147, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8300, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+774,int_stack+1307,int_stack+1207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8360, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1074,int_stack+774,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9493, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9673,int_stack+1074,int_stack+558,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+1547,int_stack+1457, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5301, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+1697,int_stack+1547, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5481, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12385,int_stack+540,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12925,int_stack+12385,int_stack+1074,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+1958,int_stack+1922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7928, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12493,int_stack+2018,int_stack+1958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7964, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12673,int_stack+12493,int_stack+12385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9385, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+2168,int_stack+2108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8300, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+270,int_stack+2268,int_stack+2168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8360, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+570,int_stack+270,int_stack+12385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9493, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+930,int_stack+570,int_stack+12673,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+2508,int_stack+2418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5301, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1578,int_stack+2658,int_stack+2508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5481, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2028,int_stack+1578,int_stack+12385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+14005,int_stack+2028,int_stack+570,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+2919,int_stack+2883, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12493,int_stack+2979,int_stack+2919, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12673,int_stack+12493,int_stack+12385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+3129,int_stack+3069, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1578,int_stack+3229,int_stack+3129, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1878,int_stack+1578,int_stack+12385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2238,int_stack+1878,int_stack+12673,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+3469,int_stack+3379, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5301, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2886,int_stack+3619,int_stack+3469, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5481, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+270,int_stack+2886,int_stack+12385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15085,int_stack+270,int_stack+1878,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+3880,int_stack+3844, 0.0, zero_stack, 1.0, int_stack+7928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+378,int_stack+3940,int_stack+3880, 0.0, zero_stack, 1.0, int_stack+7964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+558,int_stack+378,int_stack+270, 0.0, zero_stack, 1.0, int_stack+9385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+4090,int_stack+4030, 0.0, zero_stack, 1.0, int_stack+8300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12385,int_stack+4190,int_stack+4090, 0.0, zero_stack, 1.0, int_stack+8360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2886,int_stack+12385,int_stack+270, 0.0, zero_stack, 1.0, int_stack+9493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3246,int_stack+2886,int_stack+558,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+4430,int_stack+4340, 0.0, zero_stack, 1.0, int_stack+5301, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12385,int_stack+4580,int_stack+4430, 0.0, zero_stack, 1.0, int_stack+5481, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3894,int_stack+12385,int_stack+270, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+16165,int_stack+3894,int_stack+2886,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2886,int_stack+4841,int_stack+4805, 1.0, int_stack+7928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2994,int_stack+4901,int_stack+4841, 1.0, int_stack+7964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3894,int_stack+2994,int_stack+2886, 1.0, int_stack+9385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2886,int_stack+5051,int_stack+4991, 1.0, int_stack+8300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4110,int_stack+5151,int_stack+5051, 1.0, int_stack+8360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4410,int_stack+4110,int_stack+2886, 1.0, int_stack+9493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+270,int_stack+4410,int_stack+3894,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3894,int_stack+5631,int_stack+5391, 1.0, int_stack+5301, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4770,int_stack+5781,int_stack+5631, 1.0, int_stack+5481, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12385,int_stack+4770,int_stack+3894, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+4770,int_stack+12385,int_stack+4410,36);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12385,int_stack+8024,int_stack+7964,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+12565,int_stack+12385,int_stack+9385,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3894,int_stack+8460,int_stack+8360,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2886,int_stack+3894,int_stack+9493,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3894,int_stack+2886,int_stack+12565,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4542,int_stack+6042,int_stack+6006,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12385,int_stack+6102,int_stack+6042,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+9385,int_stack+12385,int_stack+4542,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+6252,int_stack+6192,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5850,int_stack+6352,int_stack+6252,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1578,int_stack+5850,int_stack+12385,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5850,int_stack+1578,int_stack+9385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+6592,int_stack+6502,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+17245,int_stack+6742,int_stack+6592,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+17695,int_stack+17245,int_stack+0,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+18235,int_stack+17695,int_stack+1578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1578,int_stack+7003,int_stack+6967,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12385,int_stack+7063,int_stack+7003,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1686,int_stack+12385,int_stack+1578,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+7213,int_stack+7153,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1902,int_stack+7313,int_stack+7213,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+17245,int_stack+1902,int_stack+12385,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6498,int_stack+17245,int_stack+1686, 0.0, zero_stack, 1.0, int_stack+12565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+7553,int_stack+7463,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1578,int_stack+7703,int_stack+7553,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+17605,int_stack+1578,int_stack+0,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+19315,int_stack+17605,int_stack+17245, 0.0, zero_stack, 1.0, int_stack+2886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+17245,int_stack+8150,int_stack+8114,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12385,int_stack+8210,int_stack+8150,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+17353,int_stack+12385,int_stack+17245,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12385,int_stack+8670,int_stack+8610,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+17569,int_stack+8770,int_stack+8670,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+17869,int_stack+17569,int_stack+12385,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1578,int_stack+17869,int_stack+17353, 1.0, int_stack+12565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+9010,int_stack+8920,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12385,int_stack+9160,int_stack+9010,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+17245,int_stack+12385,int_stack+0,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+7146,int_stack+17245,int_stack+17869, 1.0, int_stack+2886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+8226,int_stack+11305,int_stack+10657,36);
     Libderiv->ABCD[11] = int_stack + 8226;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+10321,int_stack+12925,int_stack+9673,36);
     Libderiv->ABCD[10] = int_stack + 10321;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+11617,int_stack+14005,int_stack+930,36);
     Libderiv->ABCD[9] = int_stack + 11617;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+12913,int_stack+15085,int_stack+2238,36);
     Libderiv->ABCD[8] = int_stack + 12913;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+14209,int_stack+16165,int_stack+3246,36);
     Libderiv->ABCD[7] = int_stack + 14209;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+2226,int_stack+4770,int_stack+270,36);
     Libderiv->ABCD[6] = int_stack + 2226;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+18235,int_stack+5850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[2] = int_stack + 0;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+4542,int_stack+19315,int_stack+6498, 0.0, zero_stack, 1.0, int_stack+3894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[1] = int_stack + 4542;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+5838,int_stack+7146,int_stack+1578, 1.0, int_stack+3894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[0] = int_stack + 5838;

}
