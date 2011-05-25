#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|gf) integrals */

void d1hrr_order_dpgf(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[2][6][11] = int_stack + 216;
 Libderiv->deriv_classes[2][7][11] = int_stack + 384;
 Libderiv->deriv_classes[3][4][11] = int_stack + 600;
 Libderiv->deriv_classes[3][5][11] = int_stack + 750;
 Libderiv->deriv_classes[3][6][11] = int_stack + 960;
 Libderiv->deriv_classes[3][7][11] = int_stack + 1240;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1600;
 Libderiv->deriv_classes[2][5][10] = int_stack + 1690;
 Libderiv->deriv_classes[2][6][10] = int_stack + 1816;
 Libderiv->deriv_classes[2][7][10] = int_stack + 1984;
 Libderiv->deriv_classes[3][4][10] = int_stack + 2200;
 Libderiv->deriv_classes[3][5][10] = int_stack + 2350;
 Libderiv->deriv_classes[3][6][10] = int_stack + 2560;
 Libderiv->deriv_classes[3][7][10] = int_stack + 2840;
 Libderiv->deriv_classes[2][4][9] = int_stack + 3200;
 Libderiv->deriv_classes[2][5][9] = int_stack + 3290;
 Libderiv->deriv_classes[2][6][9] = int_stack + 3416;
 Libderiv->deriv_classes[2][7][9] = int_stack + 3584;
 Libderiv->deriv_classes[3][4][9] = int_stack + 3800;
 Libderiv->deriv_classes[3][5][9] = int_stack + 3950;
 Libderiv->deriv_classes[3][6][9] = int_stack + 4160;
 Libderiv->deriv_classes[3][7][9] = int_stack + 4440;
 Libderiv->deriv_classes[2][4][8] = int_stack + 4800;
 Libderiv->deriv_classes[2][5][8] = int_stack + 4890;
 Libderiv->deriv_classes[2][6][8] = int_stack + 5016;
 Libderiv->deriv_classes[2][7][8] = int_stack + 5184;
 Libderiv->deriv_classes[3][4][8] = int_stack + 5400;
 Libderiv->deriv_classes[3][5][8] = int_stack + 5550;
 Libderiv->deriv_classes[3][6][8] = int_stack + 5760;
 Libderiv->deriv_classes[3][7][8] = int_stack + 6040;
 Libderiv->deriv_classes[2][4][7] = int_stack + 6400;
 Libderiv->deriv_classes[2][5][7] = int_stack + 6490;
 Libderiv->deriv_classes[2][6][7] = int_stack + 6616;
 Libderiv->deriv_classes[2][7][7] = int_stack + 6784;
 Libderiv->deriv_classes[3][4][7] = int_stack + 7000;
 Libderiv->deriv_classes[3][5][7] = int_stack + 7150;
 Libderiv->deriv_classes[3][6][7] = int_stack + 7360;
 Libderiv->deriv_classes[3][7][7] = int_stack + 7640;
 Libderiv->deriv_classes[2][4][6] = int_stack + 8000;
 Libderiv->deriv_classes[2][5][6] = int_stack + 8090;
 Libderiv->deriv_classes[2][6][6] = int_stack + 8216;
 Libderiv->deriv_classes[2][7][6] = int_stack + 8384;
 Libderiv->dvrr_classes[3][4] = int_stack + 8600;
 Libderiv->deriv_classes[3][4][6] = int_stack + 8750;
 Libderiv->dvrr_classes[3][5] = int_stack + 8900;
 Libderiv->deriv_classes[3][5][6] = int_stack + 9110;
 Libderiv->dvrr_classes[3][6] = int_stack + 9320;
 Libderiv->deriv_classes[3][6][6] = int_stack + 9600;
 Libderiv->deriv_classes[3][7][6] = int_stack + 9880;
 Libderiv->deriv_classes[2][4][2] = int_stack + 10240;
 Libderiv->deriv_classes[2][5][2] = int_stack + 10330;
 Libderiv->deriv_classes[2][6][2] = int_stack + 10456;
 Libderiv->deriv_classes[2][7][2] = int_stack + 10624;
 Libderiv->deriv_classes[3][4][2] = int_stack + 10840;
 Libderiv->deriv_classes[3][5][2] = int_stack + 10990;
 Libderiv->deriv_classes[3][6][2] = int_stack + 11200;
 Libderiv->deriv_classes[3][7][2] = int_stack + 11480;
 Libderiv->deriv_classes[2][4][1] = int_stack + 11840;
 Libderiv->deriv_classes[2][5][1] = int_stack + 11930;
 Libderiv->deriv_classes[2][6][1] = int_stack + 12056;
 Libderiv->deriv_classes[2][7][1] = int_stack + 12224;
 Libderiv->deriv_classes[3][4][1] = int_stack + 12440;
 Libderiv->deriv_classes[3][5][1] = int_stack + 12590;
 Libderiv->deriv_classes[3][6][1] = int_stack + 12800;
 Libderiv->deriv_classes[3][7][1] = int_stack + 13080;
 Libderiv->dvrr_classes[2][4] = int_stack + 13440;
 Libderiv->dvrr_classes[2][5] = int_stack + 13530;
 Libderiv->dvrr_classes[2][6] = int_stack + 13656;
 Libderiv->dvrr_classes[2][7] = int_stack + 13824;
 Libderiv->deriv_classes[2][4][0] = int_stack + 14040;
 Libderiv->deriv_classes[2][5][0] = int_stack + 14130;
 Libderiv->deriv_classes[2][6][0] = int_stack + 14256;
 Libderiv->deriv_classes[2][7][0] = int_stack + 14424;
 Libderiv->deriv_classes[3][4][0] = int_stack + 14640;
 Libderiv->deriv_classes[3][5][0] = int_stack + 14790;
 Libderiv->deriv_classes[3][6][0] = int_stack + 15000;
 Libderiv->deriv_classes[3][7][0] = int_stack + 15280;
 memset(int_stack,0,125120);

 Libderiv->dvrr_stack = int_stack + 36472;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15640,int_stack+13530,int_stack+13440,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15910,int_stack+13656,int_stack+13530,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16288,int_stack+15910,int_stack+15640,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16828,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13440,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17098,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13530,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17476,int_stack+17098,int_stack+16828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15640,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+18016,int_stack+384,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13656,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+18520,int_stack+18016,int_stack+17098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15910,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+19276,int_stack+18520,int_stack+17476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16288,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16828,int_stack+8900,int_stack+8600,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+17278,int_stack+9320,int_stack+8900,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17908,int_stack+17278,int_stack+16828,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18808,int_stack+750,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8600,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+960,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8900,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20176,int_stack+0,int_stack+18808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16828,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+21076,int_stack+1240,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9320,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21916,int_stack+21076,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17278,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+21916,int_stack+20176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17908,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20176,int_stack+1690,int_stack+1600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13440, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20446,int_stack+1816,int_stack+1690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13530, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20824,int_stack+20446,int_stack+20176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15640, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+21364,int_stack+1984,int_stack+1816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13656, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21868,int_stack+21364,int_stack+20446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15910, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+22624,int_stack+21868,int_stack+20824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16288, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20176,int_stack+2350,int_stack+2200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8600, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20626,int_stack+2560,int_stack+2350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8900, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21256,int_stack+20626,int_stack+20176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16828, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1500,int_stack+2840,int_stack+2560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9320, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+23524,int_stack+1500,int_stack+20626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17278, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1500,int_stack+23524,int_stack+21256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17908, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23524,int_stack+3290,int_stack+3200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13440, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23794,int_stack+3416,int_stack+3290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13530, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24172,int_stack+23794,int_stack+23524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15640, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20176,int_stack+3584,int_stack+3416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13656, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3000,int_stack+20176,int_stack+23794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15910, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+20176,int_stack+3000,int_stack+24172, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16288, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3000,int_stack+3950,int_stack+3800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8600, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21076,int_stack+4160,int_stack+3950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8900, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21706,int_stack+21076,int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16828, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3000,int_stack+4440,int_stack+4160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9320, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+23524,int_stack+3000,int_stack+21076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17278, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3000,int_stack+23524,int_stack+21706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17908, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23524,int_stack+4890,int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23794,int_stack+5016,int_stack+4890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24172,int_stack+23794,int_stack+23524, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+21076,int_stack+5184,int_stack+5016, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21580,int_stack+21076,int_stack+23794, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4500,int_stack+21580,int_stack+24172, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21076,int_stack+5550,int_stack+5400, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21526,int_stack+5760,int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23524,int_stack+21526,int_stack+21076, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+24424,int_stack+6040,int_stack+5760, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+25264,int_stack+24424,int_stack+21526, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+21076,int_stack+25264,int_stack+23524, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23524,int_stack+6490,int_stack+6400, 0.0, zero_stack, 1.0, int_stack+13440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23794,int_stack+6616,int_stack+6490, 0.0, zero_stack, 1.0, int_stack+13530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24172,int_stack+23794,int_stack+23524, 0.0, zero_stack, 1.0, int_stack+15640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+24712,int_stack+6784,int_stack+6616, 0.0, zero_stack, 1.0, int_stack+13656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+25216,int_stack+24712,int_stack+23794, 0.0, zero_stack, 1.0, int_stack+15910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+25972,int_stack+25216,int_stack+24172, 0.0, zero_stack, 1.0, int_stack+16288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23524,int_stack+7150,int_stack+7000, 0.0, zero_stack, 1.0, int_stack+8600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23974,int_stack+7360,int_stack+7150, 0.0, zero_stack, 1.0, int_stack+8900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24604,int_stack+23974,int_stack+23524, 0.0, zero_stack, 1.0, int_stack+16828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+5400,int_stack+7640,int_stack+7360, 0.0, zero_stack, 1.0, int_stack+9320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+6240,int_stack+5400,int_stack+23974, 0.0, zero_stack, 1.0, int_stack+17278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+26872,int_stack+6240,int_stack+24604, 0.0, zero_stack, 1.0, int_stack+17908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5400,int_stack+8090,int_stack+8000, 1.0, int_stack+13440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5670,int_stack+8216,int_stack+8090, 1.0, int_stack+13530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6048,int_stack+5670,int_stack+5400, 1.0, int_stack+15640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+6588,int_stack+8384,int_stack+8216, 1.0, int_stack+13656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7092,int_stack+6588,int_stack+5670, 1.0, int_stack+15910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+23524,int_stack+7092,int_stack+6048, 1.0, int_stack+16288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5400,int_stack+9110,int_stack+8750, 1.0, int_stack+8600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5850,int_stack+9600,int_stack+9110, 1.0, int_stack+8900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6480,int_stack+5850,int_stack+5400, 1.0, int_stack+16828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+7380,int_stack+9880,int_stack+9600, 1.0, int_stack+9320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8220,int_stack+7380,int_stack+5850, 1.0, int_stack+17278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+24424,int_stack+8220,int_stack+6480, 1.0, int_stack+17908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+16828,int_stack+13824,int_stack+13656,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+17332,int_stack+16828,int_stack+15910,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+18088,int_stack+17332,int_stack+16288,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18988,int_stack+10330,int_stack+10240,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5400,int_stack+10456,int_stack+10330,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5778,int_stack+5400,int_stack+18988,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+6318,int_stack+10624,int_stack+10456,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+6822,int_stack+6318,int_stack+5400,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+7578,int_stack+6822,int_stack+5778,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5400,int_stack+10990,int_stack+10840,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5850,int_stack+11200,int_stack+10990,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6480,int_stack+5850,int_stack+5400,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8478,int_stack+11480,int_stack+11200,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+9318,int_stack+8478,int_stack+5850,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+15640,int_stack+9318,int_stack+6480,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8478,int_stack+11930,int_stack+11840,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8748,int_stack+12056,int_stack+11930,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9126,int_stack+8748,int_stack+8478,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+9666,int_stack+12224,int_stack+12056,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+10170,int_stack+9666,int_stack+8748,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+10926,int_stack+10170,int_stack+9126,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8478,int_stack+12590,int_stack+12440,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8928,int_stack+12800,int_stack+12590,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9558,int_stack+8928,int_stack+8478,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+11826,int_stack+13080,int_stack+12800,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+12666,int_stack+11826,int_stack+8928,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+5400,int_stack+12666,int_stack+9558,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11826,int_stack+14130,int_stack+14040,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12096,int_stack+14256,int_stack+14130,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12474,int_stack+12096,int_stack+11826,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+13014,int_stack+14424,int_stack+14256,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+13518,int_stack+13014,int_stack+12096,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+8478,int_stack+13518,int_stack+12474,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11826,int_stack+14790,int_stack+14640,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12276,int_stack+15000,int_stack+14790,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12906,int_stack+12276,int_stack+11826,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+13806,int_stack+15280,int_stack+15000,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+9378,int_stack+13806,int_stack+12276,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+13806,int_stack+9378,int_stack+12906,10);
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+28372,int_stack+0,int_stack+19276,150);
     Libderiv->ABCD[11] = int_stack + 28372;
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+31072,int_stack+1500,int_stack+22624,150);
     Libderiv->ABCD[10] = int_stack + 31072;
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+3000,int_stack+20176,150);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+33772,int_stack+21076,int_stack+4500,150);
     Libderiv->ABCD[8] = int_stack + 33772;
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2700,int_stack+26872,int_stack+25972,150);
     Libderiv->ABCD[7] = int_stack + 2700;
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+18988,int_stack+24424,int_stack+23524,150);
     Libderiv->ABCD[6] = int_stack + 18988;
 /*--- compute (dp|gf) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+21688,int_stack+15640,int_stack+7578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 21688;
 /*--- compute (dp|gf) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+15306,int_stack+5400,int_stack+10926, 0.0, zero_stack, 1.0, int_stack+18088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 15306;
 /*--- compute (dp|gf) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5400,int_stack+13806,int_stack+8478, 1.0, int_stack+18088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 5400;

}
