#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|fd) integrals */

void d1hrr_order_dpfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[2][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][3][11] = int_stack + 276;
 Libderiv->deriv_classes[3][4][11] = int_stack + 376;
 Libderiv->deriv_classes[3][5][11] = int_stack + 526;
 Libderiv->deriv_classes[2][3][10] = int_stack + 736;
 Libderiv->deriv_classes[2][4][10] = int_stack + 796;
 Libderiv->deriv_classes[2][5][10] = int_stack + 886;
 Libderiv->deriv_classes[3][3][10] = int_stack + 1012;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1112;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1262;
 Libderiv->deriv_classes[2][3][9] = int_stack + 1472;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1532;
 Libderiv->deriv_classes[2][5][9] = int_stack + 1622;
 Libderiv->deriv_classes[3][3][9] = int_stack + 1748;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1848;
 Libderiv->deriv_classes[3][5][9] = int_stack + 1998;
 Libderiv->deriv_classes[2][3][8] = int_stack + 2208;
 Libderiv->deriv_classes[2][4][8] = int_stack + 2268;
 Libderiv->deriv_classes[2][5][8] = int_stack + 2358;
 Libderiv->deriv_classes[3][3][8] = int_stack + 2484;
 Libderiv->deriv_classes[3][4][8] = int_stack + 2584;
 Libderiv->deriv_classes[3][5][8] = int_stack + 2734;
 Libderiv->deriv_classes[2][3][7] = int_stack + 2944;
 Libderiv->deriv_classes[2][4][7] = int_stack + 3004;
 Libderiv->deriv_classes[2][5][7] = int_stack + 3094;
 Libderiv->deriv_classes[3][3][7] = int_stack + 3220;
 Libderiv->deriv_classes[3][4][7] = int_stack + 3320;
 Libderiv->deriv_classes[3][5][7] = int_stack + 3470;
 Libderiv->deriv_classes[2][3][6] = int_stack + 3680;
 Libderiv->deriv_classes[2][4][6] = int_stack + 3740;
 Libderiv->deriv_classes[2][5][6] = int_stack + 3830;
 Libderiv->dvrr_classes[3][3] = int_stack + 3956;
 Libderiv->deriv_classes[3][3][6] = int_stack + 4056;
 Libderiv->dvrr_classes[3][4] = int_stack + 4156;
 Libderiv->deriv_classes[3][4][6] = int_stack + 4306;
 Libderiv->deriv_classes[3][5][6] = int_stack + 4456;
 Libderiv->deriv_classes[2][3][2] = int_stack + 4666;
 Libderiv->deriv_classes[2][4][2] = int_stack + 4726;
 Libderiv->deriv_classes[2][5][2] = int_stack + 4816;
 Libderiv->deriv_classes[3][3][2] = int_stack + 4942;
 Libderiv->deriv_classes[3][4][2] = int_stack + 5042;
 Libderiv->deriv_classes[3][5][2] = int_stack + 5192;
 Libderiv->deriv_classes[2][3][1] = int_stack + 5402;
 Libderiv->deriv_classes[2][4][1] = int_stack + 5462;
 Libderiv->deriv_classes[2][5][1] = int_stack + 5552;
 Libderiv->deriv_classes[3][3][1] = int_stack + 5678;
 Libderiv->deriv_classes[3][4][1] = int_stack + 5778;
 Libderiv->deriv_classes[3][5][1] = int_stack + 5928;
 Libderiv->dvrr_classes[2][3] = int_stack + 6138;
 Libderiv->dvrr_classes[2][4] = int_stack + 6198;
 Libderiv->dvrr_classes[2][5] = int_stack + 6288;
 Libderiv->deriv_classes[2][3][0] = int_stack + 6414;
 Libderiv->deriv_classes[2][4][0] = int_stack + 6474;
 Libderiv->deriv_classes[2][5][0] = int_stack + 6564;
 Libderiv->deriv_classes[3][3][0] = int_stack + 6690;
 Libderiv->deriv_classes[3][4][0] = int_stack + 6790;
 Libderiv->deriv_classes[3][5][0] = int_stack + 6940;
 memset(int_stack,0,57200);

 Libderiv->dvrr_stack = int_stack + 13990;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7150,int_stack+6198,int_stack+6138,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7330,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6138,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7510,int_stack+150,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6198,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7780,int_stack+7510,int_stack+7330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7150,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7330,int_stack+4156,int_stack+3956,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+376,int_stack+276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3956,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8440,int_stack+526,int_stack+376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4156,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+8440,int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7330,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+796,int_stack+736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6138, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8320,int_stack+886,int_stack+796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6198, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8590,int_stack+8320,int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7150, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+1112,int_stack+1012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3956, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+600,int_stack+1262,int_stack+1112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4156, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8950,int_stack+600,int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7330, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+1532,int_stack+1472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6138, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8320,int_stack+1622,int_stack+1532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6198, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+600,int_stack+8320,int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7150, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+1848,int_stack+1748, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3956, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+960,int_stack+1998,int_stack+1848, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4156, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1410,int_stack+960,int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7330, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+2268,int_stack+2208, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8320,int_stack+2358,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+960,int_stack+8320,int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+2584,int_stack+2484, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2010,int_stack+2734,int_stack+2584, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9550,int_stack+2010,int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+3004,int_stack+2944, 0.0, zero_stack, 1.0, int_stack+6138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8320,int_stack+3094,int_stack+3004, 0.0, zero_stack, 1.0, int_stack+6198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2010,int_stack+8320,int_stack+8140, 0.0, zero_stack, 1.0, int_stack+7150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+3320,int_stack+3220, 0.0, zero_stack, 1.0, int_stack+3956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2370,int_stack+3470,int_stack+3320, 0.0, zero_stack, 1.0, int_stack+4156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2820,int_stack+2370,int_stack+8140, 0.0, zero_stack, 1.0, int_stack+7330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+3740,int_stack+3680, 1.0, int_stack+6138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8320,int_stack+3830,int_stack+3740, 1.0, int_stack+6198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2370,int_stack+8320,int_stack+8140, 1.0, int_stack+7150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8140,int_stack+4306,int_stack+4056, 1.0, int_stack+3956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3420,int_stack+4456,int_stack+4306, 1.0, int_stack+4156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3870,int_stack+3420,int_stack+8140, 1.0, int_stack+7330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7330,int_stack+6288,int_stack+6198,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8140,int_stack+7330,int_stack+7150,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7150,int_stack+4726,int_stack+4666,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7330,int_stack+4816,int_stack+4726,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3420,int_stack+7330,int_stack+7150,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7150,int_stack+5042,int_stack+4942,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4470,int_stack+5192,int_stack+5042,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10150,int_stack+4470,int_stack+7150,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7150,int_stack+5462,int_stack+5402,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7330,int_stack+5552,int_stack+5462,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4470,int_stack+7330,int_stack+7150,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7150,int_stack+5778,int_stack+5678,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4830,int_stack+5928,int_stack+5778,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5280,int_stack+4830,int_stack+7150,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7150,int_stack+6474,int_stack+6414,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7330,int_stack+6564,int_stack+6474,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4830,int_stack+7330,int_stack+7150,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7150,int_stack+6790,int_stack+6690,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5880,int_stack+6940,int_stack+6790,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+6330,int_stack+5880,int_stack+7150,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+10750,int_stack+0,int_stack+7780,60);
     Libderiv->ABCD[11] = int_stack + 10750;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6930,int_stack+8950,int_stack+8590,60);
     Libderiv->ABCD[10] = int_stack + 6930;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+11830,int_stack+1410,int_stack+600,60);
     Libderiv->ABCD[9] = int_stack + 11830;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+12910,int_stack+9550,int_stack+960,60);
     Libderiv->ABCD[8] = int_stack + 12910;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+2820,int_stack+2010,60);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1080,int_stack+3870,int_stack+2370,60);
     Libderiv->ABCD[6] = int_stack + 1080;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+2160,int_stack+10150,int_stack+3420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 2160;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+3240,int_stack+5280,int_stack+4470, 0.0, zero_stack, 1.0, int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 3240;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5190,int_stack+6330,int_stack+4830, 1.0, int_stack+8140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 5190;

}
