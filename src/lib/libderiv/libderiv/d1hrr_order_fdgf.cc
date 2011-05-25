#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fdgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fd|gf) integrals */

void d1hrr_order_fdgf(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][6][11] = int_stack + 360;
 Libderiv->deriv_classes[3][7][11] = int_stack + 640;
 Libderiv->deriv_classes[4][4][11] = int_stack + 1000;
 Libderiv->deriv_classes[4][5][11] = int_stack + 1225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 1540;
 Libderiv->deriv_classes[4][7][11] = int_stack + 1960;
 Libderiv->deriv_classes[5][4][11] = int_stack + 2500;
 Libderiv->deriv_classes[5][5][11] = int_stack + 2815;
 Libderiv->deriv_classes[5][6][11] = int_stack + 3256;
 Libderiv->deriv_classes[5][7][11] = int_stack + 3844;
 Libderiv->deriv_classes[3][4][10] = int_stack + 4600;
 Libderiv->deriv_classes[3][5][10] = int_stack + 4750;
 Libderiv->deriv_classes[3][6][10] = int_stack + 4960;
 Libderiv->deriv_classes[3][7][10] = int_stack + 5240;
 Libderiv->deriv_classes[4][4][10] = int_stack + 5600;
 Libderiv->deriv_classes[4][5][10] = int_stack + 5825;
 Libderiv->deriv_classes[4][6][10] = int_stack + 6140;
 Libderiv->deriv_classes[4][7][10] = int_stack + 6560;
 Libderiv->deriv_classes[5][4][10] = int_stack + 7100;
 Libderiv->deriv_classes[5][5][10] = int_stack + 7415;
 Libderiv->deriv_classes[5][6][10] = int_stack + 7856;
 Libderiv->deriv_classes[5][7][10] = int_stack + 8444;
 Libderiv->deriv_classes[3][4][9] = int_stack + 9200;
 Libderiv->deriv_classes[3][5][9] = int_stack + 9350;
 Libderiv->deriv_classes[3][6][9] = int_stack + 9560;
 Libderiv->deriv_classes[3][7][9] = int_stack + 9840;
 Libderiv->deriv_classes[4][4][9] = int_stack + 10200;
 Libderiv->deriv_classes[4][5][9] = int_stack + 10425;
 Libderiv->deriv_classes[4][6][9] = int_stack + 10740;
 Libderiv->deriv_classes[4][7][9] = int_stack + 11160;
 Libderiv->deriv_classes[5][4][9] = int_stack + 11700;
 Libderiv->deriv_classes[5][5][9] = int_stack + 12015;
 Libderiv->deriv_classes[5][6][9] = int_stack + 12456;
 Libderiv->deriv_classes[5][7][9] = int_stack + 13044;
 Libderiv->deriv_classes[3][4][8] = int_stack + 13800;
 Libderiv->deriv_classes[3][5][8] = int_stack + 13950;
 Libderiv->deriv_classes[3][6][8] = int_stack + 14160;
 Libderiv->deriv_classes[3][7][8] = int_stack + 14440;
 Libderiv->deriv_classes[4][4][8] = int_stack + 14800;
 Libderiv->deriv_classes[4][5][8] = int_stack + 15025;
 Libderiv->deriv_classes[4][6][8] = int_stack + 15340;
 Libderiv->deriv_classes[4][7][8] = int_stack + 15760;
 Libderiv->deriv_classes[5][4][8] = int_stack + 16300;
 Libderiv->deriv_classes[5][5][8] = int_stack + 16615;
 Libderiv->deriv_classes[5][6][8] = int_stack + 17056;
 Libderiv->deriv_classes[5][7][8] = int_stack + 17644;
 Libderiv->deriv_classes[3][4][7] = int_stack + 18400;
 Libderiv->deriv_classes[3][5][7] = int_stack + 18550;
 Libderiv->deriv_classes[3][6][7] = int_stack + 18760;
 Libderiv->deriv_classes[3][7][7] = int_stack + 19040;
 Libderiv->deriv_classes[4][4][7] = int_stack + 19400;
 Libderiv->deriv_classes[4][5][7] = int_stack + 19625;
 Libderiv->deriv_classes[4][6][7] = int_stack + 19940;
 Libderiv->deriv_classes[4][7][7] = int_stack + 20360;
 Libderiv->deriv_classes[5][4][7] = int_stack + 20900;
 Libderiv->deriv_classes[5][5][7] = int_stack + 21215;
 Libderiv->deriv_classes[5][6][7] = int_stack + 21656;
 Libderiv->deriv_classes[5][7][7] = int_stack + 22244;
 Libderiv->deriv_classes[3][4][6] = int_stack + 23000;
 Libderiv->deriv_classes[3][5][6] = int_stack + 23150;
 Libderiv->deriv_classes[3][6][6] = int_stack + 23360;
 Libderiv->deriv_classes[3][7][6] = int_stack + 23640;
 Libderiv->deriv_classes[4][4][6] = int_stack + 24000;
 Libderiv->deriv_classes[4][5][6] = int_stack + 24225;
 Libderiv->deriv_classes[4][6][6] = int_stack + 24540;
 Libderiv->deriv_classes[4][7][6] = int_stack + 24960;
 Libderiv->dvrr_classes[5][4] = int_stack + 25500;
 Libderiv->deriv_classes[5][4][6] = int_stack + 25815;
 Libderiv->dvrr_classes[5][5] = int_stack + 26130;
 Libderiv->deriv_classes[5][5][6] = int_stack + 26571;
 Libderiv->dvrr_classes[5][6] = int_stack + 27012;
 Libderiv->deriv_classes[5][6][6] = int_stack + 27600;
 Libderiv->deriv_classes[5][7][6] = int_stack + 28188;
 Libderiv->deriv_classes[3][4][2] = int_stack + 28944;
 Libderiv->deriv_classes[3][5][2] = int_stack + 29094;
 Libderiv->deriv_classes[3][6][2] = int_stack + 29304;
 Libderiv->deriv_classes[3][7][2] = int_stack + 29584;
 Libderiv->deriv_classes[4][4][2] = int_stack + 29944;
 Libderiv->deriv_classes[4][5][2] = int_stack + 30169;
 Libderiv->deriv_classes[4][6][2] = int_stack + 30484;
 Libderiv->deriv_classes[4][7][2] = int_stack + 30904;
 Libderiv->deriv_classes[5][4][2] = int_stack + 31444;
 Libderiv->deriv_classes[5][5][2] = int_stack + 31759;
 Libderiv->deriv_classes[5][6][2] = int_stack + 32200;
 Libderiv->deriv_classes[5][7][2] = int_stack + 32788;
 Libderiv->deriv_classes[3][4][1] = int_stack + 33544;
 Libderiv->deriv_classes[3][5][1] = int_stack + 33694;
 Libderiv->deriv_classes[3][6][1] = int_stack + 33904;
 Libderiv->deriv_classes[3][7][1] = int_stack + 34184;
 Libderiv->deriv_classes[4][4][1] = int_stack + 34544;
 Libderiv->deriv_classes[4][5][1] = int_stack + 34769;
 Libderiv->deriv_classes[4][6][1] = int_stack + 35084;
 Libderiv->deriv_classes[4][7][1] = int_stack + 35504;
 Libderiv->deriv_classes[5][4][1] = int_stack + 36044;
 Libderiv->deriv_classes[5][5][1] = int_stack + 36359;
 Libderiv->deriv_classes[5][6][1] = int_stack + 36800;
 Libderiv->deriv_classes[5][7][1] = int_stack + 37388;
 Libderiv->dvrr_classes[3][4] = int_stack + 38144;
 Libderiv->dvrr_classes[3][5] = int_stack + 38294;
 Libderiv->dvrr_classes[3][6] = int_stack + 38504;
 Libderiv->dvrr_classes[3][7] = int_stack + 38784;
 Libderiv->deriv_classes[3][4][0] = int_stack + 39144;
 Libderiv->deriv_classes[3][5][0] = int_stack + 39294;
 Libderiv->deriv_classes[3][6][0] = int_stack + 39504;
 Libderiv->deriv_classes[3][7][0] = int_stack + 39784;
 Libderiv->dvrr_classes[4][4] = int_stack + 40144;
 Libderiv->dvrr_classes[4][5] = int_stack + 40369;
 Libderiv->dvrr_classes[4][6] = int_stack + 40684;
 Libderiv->dvrr_classes[4][7] = int_stack + 41104;
 Libderiv->deriv_classes[4][4][0] = int_stack + 41644;
 Libderiv->deriv_classes[4][5][0] = int_stack + 41869;
 Libderiv->deriv_classes[4][6][0] = int_stack + 42184;
 Libderiv->deriv_classes[4][7][0] = int_stack + 42604;
 Libderiv->deriv_classes[5][4][0] = int_stack + 43144;
 Libderiv->deriv_classes[5][5][0] = int_stack + 43459;
 Libderiv->deriv_classes[5][6][0] = int_stack + 43900;
 Libderiv->deriv_classes[5][7][0] = int_stack + 44488;
 memset(int_stack,0,361952);

 Libderiv->dvrr_stack = int_stack + 138724;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fdgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+45244,int_stack+38294,int_stack+38144,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+45694,int_stack+38504,int_stack+38294,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+46324,int_stack+45694,int_stack+45244,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47224,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38144,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47674,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38294,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48304,int_stack+47674,int_stack+47224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45244,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+49204,int_stack+640,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38504,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+50044,int_stack+49204,int_stack+47674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45694,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+51304,int_stack+50044,int_stack+48304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46324,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+47224,int_stack+40369,int_stack+40144,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+47899,int_stack+40684,int_stack+40369,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+48844,int_stack+47899,int_stack+47224,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+50194,int_stack+1225,int_stack+1000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40144,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+1540,int_stack+1225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40369,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+52804,int_stack+0,int_stack+50194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47224,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+54154,int_stack+1960,int_stack+1540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40684,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+55414,int_stack+54154,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47899,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+55414,int_stack+52804, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48844,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+52804,int_stack+0,int_stack+51304,150);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+50194,int_stack+26130,int_stack+25500,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+51139,int_stack+27012,int_stack+26130,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+57304,int_stack+51139,int_stack+50194,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59194,int_stack+2815,int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25500,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60139,int_stack+3256,int_stack+2815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26130,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+61462,int_stack+60139,int_stack+59194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50194,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+63352,int_stack+3844,int_stack+3256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27012,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+65116,int_stack+63352,int_stack+60139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51139,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+67762,int_stack+65116,int_stack+61462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57304,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+59194,int_stack+67762,int_stack+0,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+4750,int_stack+4600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38144, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+450,int_stack+4960,int_stack+4750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38294, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45244, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1980,int_stack+5240,int_stack+4960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38504, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2820,int_stack+1980,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45694, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4080,int_stack+2820,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46324, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+5825,int_stack+5600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40144, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+675,int_stack+6140,int_stack+5825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40369, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1620,int_stack+675,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47224, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+65944,int_stack+6560,int_stack+6140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40684, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+67204,int_stack+65944,int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47899, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+69094,int_stack+67204,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48844, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+71344,int_stack+69094,int_stack+4080,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65944,int_stack+7415,int_stack+7100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25500, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66889,int_stack+7856,int_stack+7415, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26130, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+66889,int_stack+65944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50194, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1890,int_stack+8444,int_stack+7856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27012, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3654,int_stack+1890,int_stack+66889, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51139, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+65944,int_stack+3654,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57304, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+65944,int_stack+69094,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65944,int_stack+9350,int_stack+9200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38144, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66394,int_stack+9560,int_stack+9350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38294, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67024,int_stack+66394,int_stack+65944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45244, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+67924,int_stack+9840,int_stack+9560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38504, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+68764,int_stack+67924,int_stack+66394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45694, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6750,int_stack+68764,int_stack+67024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46324, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65944,int_stack+10425,int_stack+10200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40144, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66619,int_stack+10740,int_stack+10425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40369, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67564,int_stack+66619,int_stack+65944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47224, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+68914,int_stack+11160,int_stack+10740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40684, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8250,int_stack+68914,int_stack+66619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47899, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+68914,int_stack+8250,int_stack+67564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48844, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+75844,int_stack+68914,int_stack+6750,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6750,int_stack+12015,int_stack+11700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25500, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7695,int_stack+12456,int_stack+12015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26130, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9018,int_stack+7695,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50194, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+65944,int_stack+13044,int_stack+12456, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27012, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10908,int_stack+65944,int_stack+7695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51139, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+80344,int_stack+10908,int_stack+9018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57304, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+6750,int_stack+80344,int_stack+68914,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80344,int_stack+13950,int_stack+13800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80794,int_stack+14160,int_stack+13950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81424,int_stack+80794,int_stack+80344, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+82324,int_stack+14440,int_stack+14160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+65944,int_stack+82324,int_stack+80794, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+82324,int_stack+65944,int_stack+81424, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65944,int_stack+15025,int_stack+14800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66619,int_stack+15340,int_stack+15025, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67564,int_stack+66619,int_stack+65944, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+68914,int_stack+15760,int_stack+15340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+80344,int_stack+68914,int_stack+66619, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47899, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+68914,int_stack+80344,int_stack+67564, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+83824,int_stack+68914,int_stack+82324,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80344,int_stack+16615,int_stack+16300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81289,int_stack+17056,int_stack+16615, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+65944,int_stack+81289,int_stack+80344, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+13500,int_stack+17644,int_stack+17056, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+15264,int_stack+13500,int_stack+81289, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+80344,int_stack+15264,int_stack+65944, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+88324,int_stack+80344,int_stack+68914,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80344,int_stack+18550,int_stack+18400, 0.0, zero_stack, 1.0, int_stack+38144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80794,int_stack+18760,int_stack+18550, 0.0, zero_stack, 1.0, int_stack+38294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81424,int_stack+80794,int_stack+80344, 0.0, zero_stack, 1.0, int_stack+45244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+82324,int_stack+19040,int_stack+18760, 0.0, zero_stack, 1.0, int_stack+38504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+65944,int_stack+82324,int_stack+80794, 0.0, zero_stack, 1.0, int_stack+45694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+82324,int_stack+65944,int_stack+81424, 0.0, zero_stack, 1.0, int_stack+46324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65944,int_stack+19625,int_stack+19400, 0.0, zero_stack, 1.0, int_stack+40144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66619,int_stack+19940,int_stack+19625, 0.0, zero_stack, 1.0, int_stack+40369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67564,int_stack+66619,int_stack+65944, 0.0, zero_stack, 1.0, int_stack+47224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+68914,int_stack+20360,int_stack+19940, 0.0, zero_stack, 1.0, int_stack+40684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+80344,int_stack+68914,int_stack+66619, 0.0, zero_stack, 1.0, int_stack+47899, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+68914,int_stack+80344,int_stack+67564, 0.0, zero_stack, 1.0, int_stack+48844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+13500,int_stack+68914,int_stack+82324,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80344,int_stack+21215,int_stack+20900, 0.0, zero_stack, 1.0, int_stack+25500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81289,int_stack+21656,int_stack+21215, 0.0, zero_stack, 1.0, int_stack+26130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+65944,int_stack+81289,int_stack+80344, 0.0, zero_stack, 1.0, int_stack+50194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+18000,int_stack+22244,int_stack+21656, 0.0, zero_stack, 1.0, int_stack+27012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+19764,int_stack+18000,int_stack+81289, 0.0, zero_stack, 1.0, int_stack+51139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+80344,int_stack+19764,int_stack+65944, 0.0, zero_stack, 1.0, int_stack+57304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+95074,int_stack+80344,int_stack+68914,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80344,int_stack+23150,int_stack+23000, 1.0, int_stack+38144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80794,int_stack+23360,int_stack+23150, 1.0, int_stack+38294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81424,int_stack+80794,int_stack+80344, 1.0, int_stack+45244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+82324,int_stack+23640,int_stack+23360, 1.0, int_stack+38504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+65944,int_stack+82324,int_stack+80794, 1.0, int_stack+45694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+82324,int_stack+65944,int_stack+81424, 1.0, int_stack+46324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65944,int_stack+24225,int_stack+24000, 1.0, int_stack+40144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66619,int_stack+24540,int_stack+24225, 1.0, int_stack+40369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67564,int_stack+66619,int_stack+65944, 1.0, int_stack+47224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+68914,int_stack+24960,int_stack+24540, 1.0, int_stack+40684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+80344,int_stack+68914,int_stack+66619, 1.0, int_stack+47899, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+68914,int_stack+80344,int_stack+67564, 1.0, int_stack+48844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+18000,int_stack+68914,int_stack+82324,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80344,int_stack+26571,int_stack+25815, 1.0, int_stack+25500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81289,int_stack+27600,int_stack+26571, 1.0, int_stack+26130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+65944,int_stack+81289,int_stack+80344, 1.0, int_stack+50194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+22500,int_stack+28188,int_stack+27600, 1.0, int_stack+27012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+24264,int_stack+22500,int_stack+81289, 1.0, int_stack+51139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+80344,int_stack+24264,int_stack+65944, 1.0, int_stack+57304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+101824,int_stack+80344,int_stack+68914,150);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+80344,int_stack+38784,int_stack+38504,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+81184,int_stack+80344,int_stack+45694,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+57304,int_stack+81184,int_stack+46324,10);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+80344,int_stack+41104,int_stack+40684,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+81604,int_stack+80344,int_stack+47899,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+65944,int_stack+81604,int_stack+48844,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+22500,int_stack+65944,int_stack+57304,150);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+80344,int_stack+29094,int_stack+28944,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+80794,int_stack+29304,int_stack+29094,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81424,int_stack+80794,int_stack+80344,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+82324,int_stack+29584,int_stack+29304,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+68194,int_stack+82324,int_stack+80794,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+82324,int_stack+68194,int_stack+81424,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+68194,int_stack+30169,int_stack+29944,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+68869,int_stack+30484,int_stack+30169,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+69814,int_stack+68869,int_stack+68194,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+80344,int_stack+30904,int_stack+30484,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+27000,int_stack+80344,int_stack+68869,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+28890,int_stack+27000,int_stack+69814,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+45244,int_stack+28890,int_stack+82324, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27000,int_stack+31759,int_stack+31444,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+80344,int_stack+32200,int_stack+31759,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81667,int_stack+80344,int_stack+27000,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+27000,int_stack+32788,int_stack+32200,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+68194,int_stack+27000,int_stack+80344,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+108574,int_stack+68194,int_stack+81667,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+111724,int_stack+108574,int_stack+28890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+108574,int_stack+33694,int_stack+33544,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+109024,int_stack+33904,int_stack+33694,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+109654,int_stack+109024,int_stack+108574,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+110554,int_stack+34184,int_stack+33904,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+68194,int_stack+110554,int_stack+109024,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+40144,int_stack+68194,int_stack+109654,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+68194,int_stack+34769,int_stack+34544,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+68869,int_stack+35084,int_stack+34769,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+69814,int_stack+68869,int_stack+68194,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+108574,int_stack+35504,int_stack+35084,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+109834,int_stack+108574,int_stack+68869,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+80344,int_stack+109834,int_stack+69814,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+27000,int_stack+80344,int_stack+40144, 0.0, zero_stack, 1.0, int_stack+57304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+40144,int_stack+36359,int_stack+36044,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+108574,int_stack+36800,int_stack+36359,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+68194,int_stack+108574,int_stack+40144,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+109897,int_stack+37388,int_stack+36800,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+31500,int_stack+109897,int_stack+108574,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+108574,int_stack+31500,int_stack+68194,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+31500,int_stack+108574,int_stack+80344, 0.0, zero_stack, 1.0, int_stack+65944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+80344,int_stack+39294,int_stack+39144,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+80794,int_stack+39504,int_stack+39294,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81424,int_stack+80794,int_stack+80344,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+82324,int_stack+39784,int_stack+39504,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+108574,int_stack+82324,int_stack+80794,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+82324,int_stack+108574,int_stack+81424,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+108574,int_stack+41869,int_stack+41644,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+109249,int_stack+42184,int_stack+41869,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+110194,int_stack+109249,int_stack+108574,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+80344,int_stack+42604,int_stack+42184,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+68194,int_stack+80344,int_stack+109249,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+38250,int_stack+68194,int_stack+110194,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+118474,int_stack+38250,int_stack+82324, 1.0, int_stack+57304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+57304,int_stack+43459,int_stack+43144,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+68194,int_stack+43900,int_stack+43459,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+80344,int_stack+68194,int_stack+57304,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+57304,int_stack+44488,int_stack+43900,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+108574,int_stack+57304,int_stack+68194,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+68194,int_stack+108574,int_stack+80344,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+122974,int_stack+68194,int_stack+38250, 1.0, int_stack+65944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+129724,int_stack+59194,int_stack+52804,150);
     Libderiv->ABCD[11] = int_stack + 129724;
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+49744,int_stack+0,int_stack+71344,150);
     Libderiv->ABCD[10] = int_stack + 49744;
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+58744,int_stack+6750,int_stack+75844,150);
     Libderiv->ABCD[9] = int_stack + 58744;
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+0,int_stack+88324,int_stack+83824,150);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+67744,int_stack+95074,int_stack+13500,150);
     Libderiv->ABCD[7] = int_stack + 67744;
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+9000,int_stack+101824,int_stack+18000,150);
     Libderiv->ABCD[6] = int_stack + 9000;
 /*--- compute (fd|gf) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+76744,int_stack+111724,int_stack+45244, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 76744;
 /*--- compute (fd|gf) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+38250,int_stack+31500,int_stack+27000, 0.0, zero_stack, 1.0, int_stack+22500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 38250;
 /*--- compute (fd|gf) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+27000,int_stack+122974,int_stack+118474, 1.0, int_stack+22500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 27000;

}
