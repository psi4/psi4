#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_d0ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|ff) integrals */

void d12hrr_order_d0ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][6][11] = int_stack + 0;
 Libderiv->deriv_classes[2][6][10] = int_stack + 168;
 Libderiv->deriv_classes[2][6][9] = int_stack + 336;
 Libderiv->deriv_classes[2][6][8] = int_stack + 504;
 Libderiv->deriv_classes[2][6][7] = int_stack + 672;
 Libderiv->dvrr_classes[2][5] = int_stack + 840;
 Libderiv->deriv_classes[2][6][6] = int_stack + 966;
 Libderiv->deriv_classes[2][6][2] = int_stack + 1134;
 Libderiv->deriv_classes[2][6][1] = int_stack + 1302;
 Libderiv->deriv_classes[2][6][0] = int_stack + 1470;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 1638;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 1698;
 Libderiv->deriv2_classes[2][5][143] = int_stack + 1788;
 Libderiv->deriv2_classes[2][6][143] = int_stack + 1914;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 2082;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 2142;
 Libderiv->deriv2_classes[2][5][131] = int_stack + 2232;
 Libderiv->deriv2_classes[2][6][131] = int_stack + 2358;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 2526;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 2586;
 Libderiv->deriv2_classes[2][5][130] = int_stack + 2676;
 Libderiv->deriv2_classes[2][6][130] = int_stack + 2802;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 2970;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 3030;
 Libderiv->deriv2_classes[2][5][119] = int_stack + 3120;
 Libderiv->deriv2_classes[2][6][119] = int_stack + 3246;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 3414;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 3474;
 Libderiv->deriv2_classes[2][5][118] = int_stack + 3564;
 Libderiv->deriv2_classes[2][6][118] = int_stack + 3690;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 3858;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 3918;
 Libderiv->deriv2_classes[2][5][117] = int_stack + 4008;
 Libderiv->deriv2_classes[2][6][117] = int_stack + 4134;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 4302;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 4362;
 Libderiv->deriv2_classes[2][5][107] = int_stack + 4452;
 Libderiv->deriv2_classes[2][6][107] = int_stack + 4578;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 4746;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 4806;
 Libderiv->deriv2_classes[2][5][106] = int_stack + 4896;
 Libderiv->deriv2_classes[2][6][106] = int_stack + 5022;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 5190;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 5250;
 Libderiv->deriv2_classes[2][5][105] = int_stack + 5340;
 Libderiv->deriv2_classes[2][6][105] = int_stack + 5466;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 5634;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 5694;
 Libderiv->deriv2_classes[2][5][104] = int_stack + 5784;
 Libderiv->deriv2_classes[2][6][104] = int_stack + 5910;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 6078;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 6138;
 Libderiv->deriv2_classes[2][5][95] = int_stack + 6228;
 Libderiv->deriv2_classes[2][6][95] = int_stack + 6354;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 6522;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 6582;
 Libderiv->deriv2_classes[2][5][94] = int_stack + 6672;
 Libderiv->deriv2_classes[2][6][94] = int_stack + 6798;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 6966;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 7026;
 Libderiv->deriv2_classes[2][5][93] = int_stack + 7116;
 Libderiv->deriv2_classes[2][6][93] = int_stack + 7242;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 7410;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 7470;
 Libderiv->deriv2_classes[2][5][92] = int_stack + 7560;
 Libderiv->deriv2_classes[2][6][92] = int_stack + 7686;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 7854;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 7914;
 Libderiv->deriv2_classes[2][5][91] = int_stack + 8004;
 Libderiv->deriv2_classes[2][6][91] = int_stack + 8130;
 Libderiv->deriv_classes[2][3][11] = int_stack + 8298;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 8358;
 Libderiv->deriv_classes[2][4][11] = int_stack + 8418;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 8508;
 Libderiv->deriv_classes[2][5][11] = int_stack + 8598;
 Libderiv->deriv2_classes[2][5][83] = int_stack + 8724;
 Libderiv->deriv2_classes[2][6][83] = int_stack + 8850;
 Libderiv->deriv_classes[2][3][10] = int_stack + 9018;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 9078;
 Libderiv->deriv_classes[2][4][10] = int_stack + 9138;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 9228;
 Libderiv->deriv_classes[2][5][10] = int_stack + 9318;
 Libderiv->deriv2_classes[2][5][82] = int_stack + 9444;
 Libderiv->deriv2_classes[2][6][82] = int_stack + 9570;
 Libderiv->deriv_classes[2][3][9] = int_stack + 9738;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 9798;
 Libderiv->deriv_classes[2][4][9] = int_stack + 9858;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 9948;
 Libderiv->deriv_classes[2][5][9] = int_stack + 10038;
 Libderiv->deriv2_classes[2][5][81] = int_stack + 10164;
 Libderiv->deriv2_classes[2][6][81] = int_stack + 10290;
 Libderiv->deriv_classes[2][3][8] = int_stack + 10458;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 10518;
 Libderiv->deriv_classes[2][4][8] = int_stack + 10578;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 10668;
 Libderiv->deriv_classes[2][5][8] = int_stack + 10758;
 Libderiv->deriv2_classes[2][5][80] = int_stack + 10884;
 Libderiv->deriv2_classes[2][6][80] = int_stack + 11010;
 Libderiv->deriv_classes[2][3][7] = int_stack + 11178;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 11238;
 Libderiv->deriv_classes[2][4][7] = int_stack + 11298;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 11388;
 Libderiv->deriv_classes[2][5][7] = int_stack + 11478;
 Libderiv->deriv2_classes[2][5][79] = int_stack + 11604;
 Libderiv->deriv2_classes[2][6][79] = int_stack + 11730;
 Libderiv->dvrr_classes[2][3] = int_stack + 11898;
 Libderiv->deriv_classes[2][3][6] = int_stack + 11958;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 12018;
 Libderiv->dvrr_classes[2][4] = int_stack + 12078;
 Libderiv->deriv_classes[2][4][6] = int_stack + 12168;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 12258;
 Libderiv->deriv_classes[2][5][6] = int_stack + 12348;
 Libderiv->deriv2_classes[2][5][78] = int_stack + 12474;
 Libderiv->deriv2_classes[2][6][78] = int_stack + 12600;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 12768;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 12828;
 Libderiv->deriv2_classes[2][5][35] = int_stack + 12918;
 Libderiv->deriv2_classes[2][6][35] = int_stack + 13044;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 13212;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 13272;
 Libderiv->deriv2_classes[2][5][34] = int_stack + 13362;
 Libderiv->deriv2_classes[2][6][34] = int_stack + 13488;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 13656;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 13716;
 Libderiv->deriv2_classes[2][5][33] = int_stack + 13806;
 Libderiv->deriv2_classes[2][6][33] = int_stack + 13932;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 14100;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 14160;
 Libderiv->deriv2_classes[2][5][32] = int_stack + 14250;
 Libderiv->deriv2_classes[2][6][32] = int_stack + 14376;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 14544;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 14604;
 Libderiv->deriv2_classes[2][5][31] = int_stack + 14694;
 Libderiv->deriv2_classes[2][6][31] = int_stack + 14820;
 Libderiv->deriv_classes[2][3][2] = int_stack + 14988;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 15048;
 Libderiv->deriv_classes[2][4][2] = int_stack + 15108;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 15198;
 Libderiv->deriv_classes[2][5][2] = int_stack + 15288;
 Libderiv->deriv2_classes[2][5][30] = int_stack + 15414;
 Libderiv->deriv2_classes[2][6][30] = int_stack + 15540;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 15708;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 15768;
 Libderiv->deriv2_classes[2][5][26] = int_stack + 15858;
 Libderiv->deriv2_classes[2][6][26] = int_stack + 15984;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 16152;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 16212;
 Libderiv->deriv2_classes[2][5][23] = int_stack + 16302;
 Libderiv->deriv2_classes[2][6][23] = int_stack + 16428;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 16596;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 16656;
 Libderiv->deriv2_classes[2][5][22] = int_stack + 16746;
 Libderiv->deriv2_classes[2][6][22] = int_stack + 16872;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 17040;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 17100;
 Libderiv->deriv2_classes[2][5][21] = int_stack + 17190;
 Libderiv->deriv2_classes[2][6][21] = int_stack + 17316;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 17484;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 17544;
 Libderiv->deriv2_classes[2][5][20] = int_stack + 17634;
 Libderiv->deriv2_classes[2][6][20] = int_stack + 17760;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 17928;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 17988;
 Libderiv->deriv2_classes[2][5][19] = int_stack + 18078;
 Libderiv->deriv2_classes[2][6][19] = int_stack + 18204;
 Libderiv->deriv_classes[2][3][1] = int_stack + 18372;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 18432;
 Libderiv->deriv_classes[2][4][1] = int_stack + 18492;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 18582;
 Libderiv->deriv_classes[2][5][1] = int_stack + 18672;
 Libderiv->deriv2_classes[2][5][18] = int_stack + 18798;
 Libderiv->deriv2_classes[2][6][18] = int_stack + 18924;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 19092;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 19152;
 Libderiv->deriv2_classes[2][5][14] = int_stack + 19242;
 Libderiv->deriv2_classes[2][6][14] = int_stack + 19368;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 19536;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 19596;
 Libderiv->deriv2_classes[2][5][13] = int_stack + 19686;
 Libderiv->deriv2_classes[2][6][13] = int_stack + 19812;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 19980;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 20040;
 Libderiv->deriv2_classes[2][5][11] = int_stack + 20130;
 Libderiv->deriv2_classes[2][6][11] = int_stack + 20256;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 20424;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 20484;
 Libderiv->deriv2_classes[2][5][10] = int_stack + 20574;
 Libderiv->deriv2_classes[2][6][10] = int_stack + 20700;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 20868;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 20928;
 Libderiv->deriv2_classes[2][5][9] = int_stack + 21018;
 Libderiv->deriv2_classes[2][6][9] = int_stack + 21144;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 21312;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 21372;
 Libderiv->deriv2_classes[2][5][8] = int_stack + 21462;
 Libderiv->deriv2_classes[2][6][8] = int_stack + 21588;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 21756;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 21816;
 Libderiv->deriv2_classes[2][5][7] = int_stack + 21906;
 Libderiv->deriv2_classes[2][6][7] = int_stack + 22032;
 Libderiv->deriv_classes[2][3][0] = int_stack + 22200;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 22260;
 Libderiv->deriv_classes[2][4][0] = int_stack + 22320;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 22410;
 Libderiv->deriv_classes[2][5][0] = int_stack + 22500;
 Libderiv->deriv2_classes[2][5][6] = int_stack + 22626;
 Libderiv->deriv2_classes[2][6][6] = int_stack + 22752;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 22920;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 22980;
 Libderiv->deriv2_classes[2][5][2] = int_stack + 23070;
 Libderiv->deriv2_classes[2][6][2] = int_stack + 23196;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 23364;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 23424;
 Libderiv->deriv2_classes[2][5][1] = int_stack + 23514;
 Libderiv->deriv2_classes[2][6][1] = int_stack + 23640;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 23808;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 23868;
 Libderiv->deriv2_classes[2][5][0] = int_stack + 23958;
 Libderiv->deriv2_classes[2][6][0] = int_stack + 24084;
 memset(int_stack,0,194016);

 Libderiv->dvrr_stack = int_stack + 55626;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_d0ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24252,int_stack+12078,int_stack+11898,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24432,int_stack+840,int_stack+12078,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+24702,int_stack+24432,int_stack+24252,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25062,int_stack+8418,int_stack+8298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11898,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25242,int_stack+8598,int_stack+8418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12078,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25512,int_stack+25242,int_stack+25062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24252,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+25872,int_stack+0,int_stack+8598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+26250,int_stack+25872,int_stack+25242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24432,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25872,int_stack+9138,int_stack+9018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11898, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26790,int_stack+9318,int_stack+9138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12078, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+27060,int_stack+26790,int_stack+25872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24252, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+27420,int_stack+168,int_stack+9318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27798,int_stack+27420,int_stack+26790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24432, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27420,int_stack+9858,int_stack+9738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11898, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+10038,int_stack+9858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12078, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+28338,int_stack+0,int_stack+27420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24252, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+28698,int_stack+336,int_stack+10038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29076,int_stack+28698,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24432, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28698,int_stack+10578,int_stack+10458, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29616,int_stack+10758,int_stack+10578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29886,int_stack+29616,int_stack+28698, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30246,int_stack+504,int_stack+10758, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+30624,int_stack+30246,int_stack+29616, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30246,int_stack+11298,int_stack+11178, 0.0, zero_stack, 1.0, int_stack+11898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+11478,int_stack+11298, 0.0, zero_stack, 1.0, int_stack+12078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31164,int_stack+270,int_stack+30246, 0.0, zero_stack, 1.0, int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+31524,int_stack+672,int_stack+11478, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31902,int_stack+31524,int_stack+270, 0.0, zero_stack, 1.0, int_stack+24432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+31524,int_stack+12168,int_stack+11958, 1.0, int_stack+11898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+12348,int_stack+12168, 1.0, int_stack+12078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32442,int_stack+540,int_stack+31524, 1.0, int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+32802,int_stack+966,int_stack+12348, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33180,int_stack+32802,int_stack+540, 1.0, int_stack+24432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+32802,int_stack+15108,int_stack+14988,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24252,int_stack+15288,int_stack+15108,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+33720,int_stack+24252,int_stack+32802,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+34080,int_stack+1134,int_stack+15288,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+34458,int_stack+34080,int_stack+24252,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24522,int_stack+18492,int_stack+18372,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+34080,int_stack+18672,int_stack+18492,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+810,int_stack+34080,int_stack+24522,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+34998,int_stack+1302,int_stack+18672,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+35376,int_stack+34998,int_stack+34080,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+34998,int_stack+22320,int_stack+22200,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1170,int_stack+22500,int_stack+22320,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+35916,int_stack+1170,int_stack+34998,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+36276,int_stack+1470,int_stack+22500,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+36654,int_stack+36276,int_stack+1170,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36276,int_stack+1698,int_stack+1638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8298,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+1788,int_stack+1698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8418,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+37464,int_stack+37194,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25062,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36276,int_stack+1914,int_stack+1788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8598,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1440,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25242,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37194,int_stack+2142,int_stack+2082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8298, 1.0, int_stack+9018,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36276,int_stack+2232,int_stack+2142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8418, 1.0, int_stack+9138,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+37824,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25062, 1.0, int_stack+25872,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38184,int_stack+2358,int_stack+2232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8598, 1.0, int_stack+9318,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1980,int_stack+38184,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25242, 1.0, int_stack+26790,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36276,int_stack+2586,int_stack+2526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9018, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+2676,int_stack+2586, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9138, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+38184,int_stack+37194,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25872, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36276,int_stack+2802,int_stack+2676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9318, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+38544,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+26790, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37194,int_stack+3030,int_stack+2970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8298, 0.0, zero_stack, 1.0, int_stack+9738,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36276,int_stack+3120,int_stack+3030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8418, 0.0, zero_stack, 1.0, int_stack+9858,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+39084,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25062, 0.0, zero_stack, 1.0, int_stack+27420,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+39444,int_stack+3246,int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8598, 0.0, zero_stack, 1.0, int_stack+10038,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2520,int_stack+39444,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25242, 0.0, zero_stack, 1.0, int_stack+0,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36276,int_stack+3474,int_stack+3414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9018, 1.0, int_stack+9738, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+3564,int_stack+3474, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9138, 1.0, int_stack+9858, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+39444,int_stack+37194,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25872, 1.0, int_stack+27420, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36276,int_stack+3690,int_stack+3564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9318, 1.0, int_stack+10038, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3060,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26790, 1.0, int_stack+0, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37194,int_stack+3918,int_stack+3858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9738, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36276,int_stack+4008,int_stack+3918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9858, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3600,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27420, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+39804,int_stack+4134,int_stack+4008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10038, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+40182,int_stack+39804,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36276,int_stack+4362,int_stack+4302, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8298, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10458,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+4452,int_stack+4362, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8418, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10578,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+39804,int_stack+37194,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25062, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28698,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36276,int_stack+4578,int_stack+4452, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8598, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10758,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3960,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25242, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29616,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37194,int_stack+4806,int_stack+4746, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9018, 0.0, zero_stack, 1.0, int_stack+10458, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36276,int_stack+4896,int_stack+4806, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9138, 0.0, zero_stack, 1.0, int_stack+10578, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4500,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25872, 0.0, zero_stack, 1.0, int_stack+28698, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+40722,int_stack+5022,int_stack+4896, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9318, 0.0, zero_stack, 1.0, int_stack+10758, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41100,int_stack+40722,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26790, 0.0, zero_stack, 1.0, int_stack+29616, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36276,int_stack+5250,int_stack+5190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9738, 1.0, int_stack+10458, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+5340,int_stack+5250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9858, 1.0, int_stack+10578, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+40722,int_stack+37194,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27420, 1.0, int_stack+28698, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36276,int_stack+5466,int_stack+5340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10038, 1.0, int_stack+10758, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4860,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+29616, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37194,int_stack+5694,int_stack+5634, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36276,int_stack+5784,int_stack+5694, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41640,int_stack+36276,int_stack+37194, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5400,int_stack+5910,int_stack+5784, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42000,int_stack+5400,int_stack+36276, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36276,int_stack+6138,int_stack+6078, 0.0, zero_stack, 1.0, int_stack+8298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11178,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+6228,int_stack+6138, 0.0, zero_stack, 1.0, int_stack+8418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11298,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5400,int_stack+37194,int_stack+36276, 0.0, zero_stack, 1.0, int_stack+25062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30246,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36276,int_stack+6354,int_stack+6228, 0.0, zero_stack, 1.0, int_stack+8598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11478,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5760,int_stack+36276,int_stack+37194, 0.0, zero_stack, 1.0, int_stack+25242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37194,int_stack+6582,int_stack+6522, 0.0, zero_stack, 1.0, int_stack+9018, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11178, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36276,int_stack+6672,int_stack+6582, 0.0, zero_stack, 1.0, int_stack+9138, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11298, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6300,int_stack+36276,int_stack+37194, 0.0, zero_stack, 1.0, int_stack+25872, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30246, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42540,int_stack+6798,int_stack+6672, 0.0, zero_stack, 1.0, int_stack+9318, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11478, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42918,int_stack+42540,int_stack+36276, 0.0, zero_stack, 1.0, int_stack+26790, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36276,int_stack+7026,int_stack+6966, 0.0, zero_stack, 1.0, int_stack+9738, 0.0, zero_stack, 1.0, int_stack+11178, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+7116,int_stack+7026, 0.0, zero_stack, 1.0, int_stack+9858, 0.0, zero_stack, 1.0, int_stack+11298, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42540,int_stack+37194,int_stack+36276, 0.0, zero_stack, 1.0, int_stack+27420, 0.0, zero_stack, 1.0, int_stack+30246, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36276,int_stack+7242,int_stack+7116, 0.0, zero_stack, 1.0, int_stack+10038, 0.0, zero_stack, 1.0, int_stack+11478, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6660,int_stack+36276,int_stack+37194, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37194,int_stack+7470,int_stack+7410, 0.0, zero_stack, 1.0, int_stack+10458, 1.0, int_stack+11178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7200,int_stack+7560,int_stack+7470, 0.0, zero_stack, 1.0, int_stack+10578, 1.0, int_stack+11298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+36276,int_stack+7200,int_stack+37194, 0.0, zero_stack, 1.0, int_stack+28698, 1.0, int_stack+30246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43458,int_stack+7686,int_stack+7560, 0.0, zero_stack, 1.0, int_stack+10758, 1.0, int_stack+11478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+43836,int_stack+43458,int_stack+7200, 0.0, zero_stack, 1.0, int_stack+29616, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+7914,int_stack+7854, 0.0, zero_stack, 2.0, int_stack+11178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+8004,int_stack+7914, 0.0, zero_stack, 2.0, int_stack+11298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7380,int_stack+37194,int_stack+7200, 0.0, zero_stack, 2.0, int_stack+30246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43458,int_stack+8130,int_stack+8004, 0.0, zero_stack, 2.0, int_stack+11478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7740,int_stack+43458,int_stack+37194, 0.0, zero_stack, 2.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+8508,int_stack+8358, 1.0, int_stack+8298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11958,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+8724,int_stack+8508, 1.0, int_stack+8418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12168,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+43458,int_stack+37194,int_stack+7200, 1.0, int_stack+25062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31524,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+44376,int_stack+8850,int_stack+8724, 1.0, int_stack+8598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12348,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8280,int_stack+44376,int_stack+37194, 1.0, int_stack+25242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+9228,int_stack+9078, 1.0, int_stack+9018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11958, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37194,int_stack+9444,int_stack+9228, 1.0, int_stack+9138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12168, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+44376,int_stack+37194,int_stack+7200, 1.0, int_stack+25872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31524, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+25872,int_stack+9570,int_stack+9444, 1.0, int_stack+9318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12348, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8820,int_stack+25872,int_stack+37194, 1.0, int_stack+26790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+9948,int_stack+9798, 1.0, int_stack+9738, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11958, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26790,int_stack+10164,int_stack+9948, 1.0, int_stack+9858, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12168, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25872,int_stack+26790,int_stack+7200, 1.0, int_stack+27420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31524, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+27420,int_stack+10290,int_stack+10164, 1.0, int_stack+10038, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12348, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9360,int_stack+27420,int_stack+26790, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+10668,int_stack+10518, 1.0, int_stack+10458, 0.0, zero_stack, 1.0, int_stack+11958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+10884,int_stack+10668, 1.0, int_stack+10578, 0.0, zero_stack, 1.0, int_stack+12168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+27420,int_stack+0,int_stack+7200, 1.0, int_stack+28698, 0.0, zero_stack, 1.0, int_stack+31524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+28698,int_stack+11010,int_stack+10884, 1.0, int_stack+10758, 0.0, zero_stack, 1.0, int_stack+12348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9900,int_stack+28698,int_stack+0, 1.0, int_stack+29616, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+11388,int_stack+11238, 1.0, int_stack+11178, 1.0, int_stack+11958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29616,int_stack+11604,int_stack+11388, 1.0, int_stack+11298, 1.0, int_stack+12168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+28698,int_stack+29616,int_stack+7200, 1.0, int_stack+30246, 1.0, int_stack+31524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30246,int_stack+11730,int_stack+11604, 1.0, int_stack+11478, 1.0, int_stack+12348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10440,int_stack+30246,int_stack+29616, 1.0, int_stack+270, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+12258,int_stack+12018, 2.0, int_stack+11958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29616,int_stack+12474,int_stack+12258, 2.0, int_stack+12168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30246,int_stack+29616,int_stack+7200, 2.0, int_stack+31524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+31524,int_stack+12600,int_stack+12474, 2.0, int_stack+12348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+31524,int_stack+29616, 2.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+12828,int_stack+12768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14988,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+12918,int_stack+12828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15108,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31524,int_stack+540,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32802,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+25062,int_stack+13044,int_stack+12918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15288,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10980,int_stack+25062,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24252,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+13272,int_stack+13212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14988, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+13362,int_stack+13272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15108, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25062,int_stack+540,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32802, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11520,int_stack+13488,int_stack+13362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15288, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11898,int_stack+11520,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24252, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+13716,int_stack+13656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14988, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+13806,int_stack+13716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15108, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11520,int_stack+540,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32802, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12438,int_stack+13932,int_stack+13806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15288, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+12816,int_stack+12438,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24252, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+14160,int_stack+14100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+14250,int_stack+14160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12438,int_stack+540,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13356,int_stack+14376,int_stack+14250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13734,int_stack+13356,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+14604,int_stack+14544, 0.0, zero_stack, 1.0, int_stack+14988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+14694,int_stack+14604, 0.0, zero_stack, 1.0, int_stack+15108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13356,int_stack+540,int_stack+7200, 0.0, zero_stack, 1.0, int_stack+32802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14274,int_stack+14820,int_stack+14694, 0.0, zero_stack, 1.0, int_stack+15288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+44736,int_stack+14274,int_stack+540, 0.0, zero_stack, 1.0, int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+15198,int_stack+15048, 1.0, int_stack+14988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+15414,int_stack+15198, 1.0, int_stack+15108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14274,int_stack+540,int_stack+7200, 1.0, int_stack+32802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+32802,int_stack+15540,int_stack+15414, 1.0, int_stack+15288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14634,int_stack+32802,int_stack+540, 1.0, int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+15768,int_stack+15708,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24252,int_stack+15858,int_stack+15768,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+32802,int_stack+24252,int_stack+7200,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15174,int_stack+15984,int_stack+15858,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15552,int_stack+15174,int_stack+24252,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+16212,int_stack+16152, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18372,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24252,int_stack+16302,int_stack+16212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18492,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15174,int_stack+24252,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24522,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45276,int_stack+16428,int_stack+16302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18672,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45654,int_stack+45276,int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34080,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+16656,int_stack+16596, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18372, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24252,int_stack+16746,int_stack+16656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18492, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45276,int_stack+24252,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24522, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16092,int_stack+16872,int_stack+16746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18672, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16470,int_stack+16092,int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34080, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+17100,int_stack+17040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18372, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24252,int_stack+17190,int_stack+17100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18492, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16092,int_stack+24252,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24522, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+46194,int_stack+17316,int_stack+17190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18672, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+46572,int_stack+46194,int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34080, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+17544,int_stack+17484, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24252,int_stack+17634,int_stack+17544, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+46194,int_stack+24252,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17010,int_stack+17760,int_stack+17634, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17388,int_stack+17010,int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+17988,int_stack+17928, 0.0, zero_stack, 1.0, int_stack+18372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24252,int_stack+18078,int_stack+17988, 0.0, zero_stack, 1.0, int_stack+18492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17010,int_stack+24252,int_stack+7200, 0.0, zero_stack, 1.0, int_stack+24522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47112,int_stack+18204,int_stack+18078, 0.0, zero_stack, 1.0, int_stack+18672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+47490,int_stack+47112,int_stack+24252, 0.0, zero_stack, 1.0, int_stack+34080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+18582,int_stack+18432, 1.0, int_stack+18372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24252,int_stack+18798,int_stack+18582, 1.0, int_stack+18492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+47112,int_stack+24252,int_stack+7200, 1.0, int_stack+24522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17928,int_stack+18924,int_stack+18798, 1.0, int_stack+18672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18306,int_stack+17928,int_stack+24252, 1.0, int_stack+34080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+19152,int_stack+19092,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+19242,int_stack+19152,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+34080,int_stack+540,int_stack+7200,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+17928,int_stack+19368,int_stack+19242,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18846,int_stack+17928,int_stack+540,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+19596,int_stack+19536,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+19686,int_stack+19596,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17928,int_stack+540,int_stack+7200,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24252,int_stack+19812,int_stack+19686,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+19386,int_stack+24252,int_stack+540,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+20040,int_stack+19980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22200,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+20130,int_stack+20040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22320,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24252,int_stack+540,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34998,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+48030,int_stack+20256,int_stack+20130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22500,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48408,int_stack+48030,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+20484,int_stack+20424, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22200, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+20574,int_stack+20484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22320, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+48030,int_stack+540,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34998, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+19926,int_stack+20700,int_stack+20574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22500, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20304,int_stack+19926,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+20928,int_stack+20868, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22200, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+21018,int_stack+20928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22320, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19926,int_stack+540,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34998, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+48948,int_stack+21144,int_stack+21018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22500, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+49326,int_stack+48948,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+21372,int_stack+21312, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+21462,int_stack+21372, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+48948,int_stack+540,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20844,int_stack+21588,int_stack+21462, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+49866,int_stack+20844,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+21816,int_stack+21756, 0.0, zero_stack, 1.0, int_stack+22200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+21906,int_stack+21816, 0.0, zero_stack, 1.0, int_stack+22320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20844,int_stack+540,int_stack+7200, 0.0, zero_stack, 1.0, int_stack+34998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21204,int_stack+22032,int_stack+21906, 0.0, zero_stack, 1.0, int_stack+22500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21582,int_stack+21204,int_stack+540, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+22410,int_stack+22260, 1.0, int_stack+22200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+22626,int_stack+22410, 1.0, int_stack+22320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21204,int_stack+540,int_stack+7200, 1.0, int_stack+34998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+34998,int_stack+22752,int_stack+22626, 1.0, int_stack+22500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22122,int_stack+34998,int_stack+540, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+22980,int_stack+22920,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1170,int_stack+23070,int_stack+22980,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+34998,int_stack+1170,int_stack+7200,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+22662,int_stack+23196,int_stack+23070,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+50406,int_stack+22662,int_stack+1170,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+23424,int_stack+23364,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1170,int_stack+23514,int_stack+23424,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+22662,int_stack+1170,int_stack+7200,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+23022,int_stack+23640,int_stack+23514,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+50946,int_stack+23022,int_stack+1170,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+23868,int_stack+23808,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1170,int_stack+23958,int_stack+23868,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+23022,int_stack+1170,int_stack+7200,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+23382,int_stack+24084,int_stack+23958,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+51486,int_stack+23382,int_stack+1170,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+23382,int_stack+26250,int_stack+25512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24702,6);
     Libderiv->ABCD[11] = int_stack + 23382;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+26232,int_stack+27798,int_stack+27060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24702, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 26232;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+52026,int_stack+29076,int_stack+28338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24702, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 52026;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+29058,int_stack+30624,int_stack+29886, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 29058;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+52626,int_stack+31902,int_stack+31164, 0.0, zero_stack, 1.0, int_stack+24702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 52626;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+53226,int_stack+33180,int_stack+32442, 1.0, int_stack+24702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 53226;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+53826,int_stack+34458,int_stack+33720,6);
     Libderiv->ABCD[2] = int_stack + 53826;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+54426,int_stack+35376,int_stack+810,6);
     Libderiv->ABCD[1] = int_stack + 54426;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+55026,int_stack+36654,int_stack+35916,6);
     Libderiv->ABCD[0] = int_stack + 55026;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+36636,int_stack+1440,int_stack+37464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25512,6);
     Libderiv->ABCD[155] = int_stack + 36636;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1170,int_stack+1980,int_stack+37824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25512, 1.0, int_stack+27060,6);
     Libderiv->ABCD[143] = int_stack + 1170;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1770,int_stack+38544,int_stack+38184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27060, 0.0, zero_stack,6);
     Libderiv->ABCD[142] = int_stack + 1770;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+37236,int_stack+2520,int_stack+39084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25512, 0.0, zero_stack, 1.0, int_stack+28338,6);
     Libderiv->ABCD[131] = int_stack + 37236;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2370,int_stack+3060,int_stack+39444, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27060, 1.0, int_stack+28338, 0.0, zero_stack,6);
     Libderiv->ABCD[130] = int_stack + 2370;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2970,int_stack+40182,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28338, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[129] = int_stack + 2970;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+37836,int_stack+3960,int_stack+39804, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25512, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29886,6);
     Libderiv->ABCD[119] = int_stack + 37836;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3570,int_stack+41100,int_stack+4500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27060, 0.0, zero_stack, 1.0, int_stack+29886, 0.0, zero_stack,6);
     Libderiv->ABCD[118] = int_stack + 3570;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4170,int_stack+4860,int_stack+40722, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28338, 1.0, int_stack+29886, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[117] = int_stack + 4170;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4770,int_stack+42000,int_stack+41640, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[116] = int_stack + 4770;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+38436,int_stack+5760,int_stack+5400, 0.0, zero_stack, 1.0, int_stack+25512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31164,6);
     Libderiv->ABCD[107] = int_stack + 38436;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5370,int_stack+42918,int_stack+6300, 0.0, zero_stack, 1.0, int_stack+27060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31164, 0.0, zero_stack,6);
     Libderiv->ABCD[106] = int_stack + 5370;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5970,int_stack+6660,int_stack+42540, 0.0, zero_stack, 1.0, int_stack+28338, 0.0, zero_stack, 1.0, int_stack+31164, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[105] = int_stack + 5970;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6570,int_stack+43836,int_stack+36276, 0.0, zero_stack, 1.0, int_stack+29886, 1.0, int_stack+31164, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[104] = int_stack + 6570;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+39036,int_stack+7740,int_stack+7380, 0.0, zero_stack, 2.0, int_stack+31164, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[103] = int_stack + 39036;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7170,int_stack+8280,int_stack+43458, 1.0, int_stack+25512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32442,6);
     Libderiv->ABCD[95] = int_stack + 7170;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7770,int_stack+8820,int_stack+44376, 1.0, int_stack+27060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32442, 0.0, zero_stack,6);
     Libderiv->ABCD[94] = int_stack + 7770;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8370,int_stack+9360,int_stack+25872, 1.0, int_stack+28338, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32442, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[93] = int_stack + 8370;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8970,int_stack+9900,int_stack+27420, 1.0, int_stack+29886, 0.0, zero_stack, 1.0, int_stack+32442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[92] = int_stack + 8970;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9570,int_stack+10440,int_stack+28698, 1.0, int_stack+31164, 1.0, int_stack+32442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[91] = int_stack + 9570;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10170,int_stack+0,int_stack+30246, 2.0, int_stack+32442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[90] = int_stack + 10170;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+10980,int_stack+31524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33720,6);
     Libderiv->ABCD[47] = int_stack + 0;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10770,int_stack+11898,int_stack+25062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33720, 0.0, zero_stack,6);
     Libderiv->ABCD[46] = int_stack + 10770;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+24612,int_stack+12816,int_stack+11520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33720, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[45] = int_stack + 24612;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11370,int_stack+13734,int_stack+12438, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[44] = int_stack + 11370;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11970,int_stack+44736,int_stack+13356, 0.0, zero_stack, 1.0, int_stack+33720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[43] = int_stack + 11970;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12570,int_stack+14634,int_stack+14274, 1.0, int_stack+33720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[42] = int_stack + 12570;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+13170,int_stack+15552,int_stack+32802,6);
     Libderiv->ABCD[38] = int_stack + 13170;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13770,int_stack+45654,int_stack+15174, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+810,6);
     Libderiv->ABCD[35] = int_stack + 13770;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14370,int_stack+16470,int_stack+45276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+810, 0.0, zero_stack,6);
     Libderiv->ABCD[34] = int_stack + 14370;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14970,int_stack+46572,int_stack+16092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+810, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[33] = int_stack + 14970;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+15570,int_stack+17388,int_stack+46194, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[32] = int_stack + 15570;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16170,int_stack+47490,int_stack+17010, 0.0, zero_stack, 1.0, int_stack+810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[31] = int_stack + 16170;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16770,int_stack+18306,int_stack+47112, 1.0, int_stack+810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[30] = int_stack + 16770;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+25212,int_stack+18846,int_stack+34080,6);
     Libderiv->ABCD[26] = int_stack + 25212;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+18288,int_stack+19386,int_stack+17928,6);
     Libderiv->ABCD[25] = int_stack + 18288;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+17370,int_stack+48408,int_stack+24252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35916,6);
     Libderiv->ABCD[23] = int_stack + 17370;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+23982,int_stack+20304,int_stack+48030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35916, 0.0, zero_stack,6);
     Libderiv->ABCD[22] = int_stack + 23982;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18888,int_stack+49326,int_stack+19926, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35916, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[21] = int_stack + 18888;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+19488,int_stack+49866,int_stack+48948, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[20] = int_stack + 19488;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+20088,int_stack+21582,int_stack+20844, 0.0, zero_stack, 1.0, int_stack+35916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[19] = int_stack + 20088;
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+29658,int_stack+22122,int_stack+21204, 1.0, int_stack+35916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[18] = int_stack + 29658;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+35358,int_stack+50406,int_stack+34998,6);
     Libderiv->ABCD[14] = int_stack + 35358;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+35958,int_stack+50946,int_stack+22662,6);
     Libderiv->ABCD[13] = int_stack + 35958;
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+20688,int_stack+51486,int_stack+23022,6);
     Libderiv->ABCD[12] = int_stack + 20688;

}
