#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ddgp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|gp) integrals */

void d1hrr_order_ddgp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[3][4][11] = int_stack + 216;
 Libderiv->deriv_classes[3][5][11] = int_stack + 366;
 Libderiv->deriv_classes[4][4][11] = int_stack + 576;
 Libderiv->deriv_classes[4][5][11] = int_stack + 801;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1116;
 Libderiv->deriv_classes[2][5][10] = int_stack + 1206;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1332;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1482;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1692;
 Libderiv->deriv_classes[4][5][10] = int_stack + 1917;
 Libderiv->deriv_classes[2][4][9] = int_stack + 2232;
 Libderiv->deriv_classes[2][5][9] = int_stack + 2322;
 Libderiv->deriv_classes[3][4][9] = int_stack + 2448;
 Libderiv->deriv_classes[3][5][9] = int_stack + 2598;
 Libderiv->deriv_classes[4][4][9] = int_stack + 2808;
 Libderiv->deriv_classes[4][5][9] = int_stack + 3033;
 Libderiv->deriv_classes[2][4][8] = int_stack + 3348;
 Libderiv->deriv_classes[2][5][8] = int_stack + 3438;
 Libderiv->deriv_classes[3][4][8] = int_stack + 3564;
 Libderiv->deriv_classes[3][5][8] = int_stack + 3714;
 Libderiv->deriv_classes[4][4][8] = int_stack + 3924;
 Libderiv->deriv_classes[4][5][8] = int_stack + 4149;
 Libderiv->deriv_classes[2][4][7] = int_stack + 4464;
 Libderiv->deriv_classes[2][5][7] = int_stack + 4554;
 Libderiv->deriv_classes[3][4][7] = int_stack + 4680;
 Libderiv->deriv_classes[3][5][7] = int_stack + 4830;
 Libderiv->deriv_classes[4][4][7] = int_stack + 5040;
 Libderiv->deriv_classes[4][5][7] = int_stack + 5265;
 Libderiv->deriv_classes[2][4][6] = int_stack + 5580;
 Libderiv->deriv_classes[2][5][6] = int_stack + 5670;
 Libderiv->deriv_classes[3][4][6] = int_stack + 5796;
 Libderiv->deriv_classes[3][5][6] = int_stack + 5946;
 Libderiv->dvrr_classes[4][4] = int_stack + 6156;
 Libderiv->deriv_classes[4][4][6] = int_stack + 6381;
 Libderiv->deriv_classes[4][5][6] = int_stack + 6606;
 Libderiv->deriv_classes[2][4][2] = int_stack + 6921;
 Libderiv->deriv_classes[2][5][2] = int_stack + 7011;
 Libderiv->deriv_classes[3][4][2] = int_stack + 7137;
 Libderiv->deriv_classes[3][5][2] = int_stack + 7287;
 Libderiv->deriv_classes[4][4][2] = int_stack + 7497;
 Libderiv->deriv_classes[4][5][2] = int_stack + 7722;
 Libderiv->deriv_classes[2][4][1] = int_stack + 8037;
 Libderiv->deriv_classes[2][5][1] = int_stack + 8127;
 Libderiv->deriv_classes[3][4][1] = int_stack + 8253;
 Libderiv->deriv_classes[3][5][1] = int_stack + 8403;
 Libderiv->deriv_classes[4][4][1] = int_stack + 8613;
 Libderiv->deriv_classes[4][5][1] = int_stack + 8838;
 Libderiv->dvrr_classes[2][4] = int_stack + 9153;
 Libderiv->dvrr_classes[2][5] = int_stack + 9243;
 Libderiv->deriv_classes[2][4][0] = int_stack + 9369;
 Libderiv->deriv_classes[2][5][0] = int_stack + 9459;
 Libderiv->dvrr_classes[3][4] = int_stack + 9585;
 Libderiv->dvrr_classes[3][5] = int_stack + 9735;
 Libderiv->deriv_classes[3][4][0] = int_stack + 9945;
 Libderiv->deriv_classes[3][5][0] = int_stack + 10095;
 Libderiv->deriv_classes[4][4][0] = int_stack + 10305;
 Libderiv->deriv_classes[4][5][0] = int_stack + 10530;
 memset(int_stack,0,86760);

 Libderiv->dvrr_stack = int_stack + 24525;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ddgp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10845,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9153,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11115,int_stack+366,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9585,10);
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+11565,int_stack+11115,int_stack+10845,45);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+801,int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6156,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+13050,int_stack+12375,int_stack+11115,45);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+1206,int_stack+1116, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9153, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10845,int_stack+1482,int_stack+1332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9585, 0.0, zero_stack,10);
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+10845,int_stack+12375,45);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+1917,int_stack+1692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6156, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+810,int_stack+12375,int_stack+10845,45);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10845,int_stack+2322,int_stack+2232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9153, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11115,int_stack+2598,int_stack+2448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9585, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14400,int_stack+11115,int_stack+10845,45);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+3033,int_stack+2808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6156, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15210,int_stack+12375,int_stack+11115,45);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+3438,int_stack+3348, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10845,int_stack+3714,int_stack+3564, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2160,int_stack+10845,int_stack+12375,45);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+4149,int_stack+3924, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+2970,int_stack+12375,int_stack+10845,45);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10845,int_stack+4554,int_stack+4464, 0.0, zero_stack, 1.0, int_stack+9153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11115,int_stack+4830,int_stack+4680, 0.0, zero_stack, 1.0, int_stack+9585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+16560,int_stack+11115,int_stack+10845,45);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+5265,int_stack+5040, 0.0, zero_stack, 1.0, int_stack+6156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+17370,int_stack+12375,int_stack+11115,45);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+5670,int_stack+5580, 1.0, int_stack+9153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10845,int_stack+5946,int_stack+5796, 1.0, int_stack+9585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4320,int_stack+10845,int_stack+12375,45);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+6606,int_stack+6381, 1.0, int_stack+6156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5130,int_stack+12375,int_stack+10845,45);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10845,int_stack+9243,int_stack+9153,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11115,int_stack+9735,int_stack+9585,10);
 /*--- compute (dp|gp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+18720,int_stack+11115,int_stack+10845,45);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9585,int_stack+7011,int_stack+6921,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+7287,int_stack+7137,10);
 /*--- compute (dp|gp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6480,int_stack+12375,int_stack+9585, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19530,int_stack+7722,int_stack+7497,15);
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+20205,int_stack+19530,int_stack+12375, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+8127,int_stack+8037,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19530,int_stack+8403,int_stack+8253,10);
 /*--- compute (dp|gp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7290,int_stack+19530,int_stack+12375, 0.0, zero_stack, 1.0, int_stack+10845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+8838,int_stack+8613,15);
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+21555,int_stack+12375,int_stack+19530, 0.0, zero_stack, 1.0, int_stack+11115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19530,int_stack+9459,int_stack+9369,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12375,int_stack+10095,int_stack+9945,10);
 /*--- compute (dp|gp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+8100,int_stack+12375,int_stack+19530, 1.0, int_stack+10845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19530,int_stack+10530,int_stack+10305,15);
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+8910,int_stack+19530,int_stack+12375, 1.0, int_stack+11115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
 /*--- compute (dd|gp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+22905,int_stack+13050,int_stack+11565,45);
     Libderiv->ABCD[11] = int_stack + 22905;
 /*--- compute (dd|gp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+10260,int_stack+810,int_stack+0,45);
     Libderiv->ABCD[10] = int_stack + 10260;
 /*--- compute (dd|gp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+15210,int_stack+14400,45);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (dd|gp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+11880,int_stack+2970,int_stack+2160,45);
     Libderiv->ABCD[8] = int_stack + 11880;
 /*--- compute (dd|gp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+1620,int_stack+17370,int_stack+16560,45);
     Libderiv->ABCD[7] = int_stack + 1620;
 /*--- compute (dd|gp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+13500,int_stack+5130,int_stack+4320,45);
     Libderiv->ABCD[6] = int_stack + 13500;
 /*--- compute (dd|gp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+3240,int_stack+20205,int_stack+6480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[2] = int_stack + 3240;
 /*--- compute (dd|gp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+4860,int_stack+21555,int_stack+7290, 0.0, zero_stack, 1.0, int_stack+18720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[1] = int_stack + 4860;
 /*--- compute (dd|gp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+6480,int_stack+8910,int_stack+8100, 1.0, int_stack+18720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[0] = int_stack + 6480;

}
