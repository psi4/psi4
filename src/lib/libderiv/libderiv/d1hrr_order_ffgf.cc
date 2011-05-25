#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ffgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (ff|gf) integrals */

void d1hrr_order_ffgf(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[6][4][11] = int_stack + 4600;
 Libderiv->deriv_classes[6][5][11] = int_stack + 5020;
 Libderiv->deriv_classes[6][6][11] = int_stack + 5608;
 Libderiv->deriv_classes[6][7][11] = int_stack + 6392;
 Libderiv->deriv_classes[3][4][10] = int_stack + 7400;
 Libderiv->deriv_classes[3][5][10] = int_stack + 7550;
 Libderiv->deriv_classes[3][6][10] = int_stack + 7760;
 Libderiv->deriv_classes[3][7][10] = int_stack + 8040;
 Libderiv->deriv_classes[4][4][10] = int_stack + 8400;
 Libderiv->deriv_classes[4][5][10] = int_stack + 8625;
 Libderiv->deriv_classes[4][6][10] = int_stack + 8940;
 Libderiv->deriv_classes[4][7][10] = int_stack + 9360;
 Libderiv->deriv_classes[5][4][10] = int_stack + 9900;
 Libderiv->deriv_classes[5][5][10] = int_stack + 10215;
 Libderiv->deriv_classes[5][6][10] = int_stack + 10656;
 Libderiv->deriv_classes[5][7][10] = int_stack + 11244;
 Libderiv->deriv_classes[6][4][10] = int_stack + 12000;
 Libderiv->deriv_classes[6][5][10] = int_stack + 12420;
 Libderiv->deriv_classes[6][6][10] = int_stack + 13008;
 Libderiv->deriv_classes[6][7][10] = int_stack + 13792;
 Libderiv->deriv_classes[3][4][9] = int_stack + 14800;
 Libderiv->deriv_classes[3][5][9] = int_stack + 14950;
 Libderiv->deriv_classes[3][6][9] = int_stack + 15160;
 Libderiv->deriv_classes[3][7][9] = int_stack + 15440;
 Libderiv->deriv_classes[4][4][9] = int_stack + 15800;
 Libderiv->deriv_classes[4][5][9] = int_stack + 16025;
 Libderiv->deriv_classes[4][6][9] = int_stack + 16340;
 Libderiv->deriv_classes[4][7][9] = int_stack + 16760;
 Libderiv->deriv_classes[5][4][9] = int_stack + 17300;
 Libderiv->deriv_classes[5][5][9] = int_stack + 17615;
 Libderiv->deriv_classes[5][6][9] = int_stack + 18056;
 Libderiv->deriv_classes[5][7][9] = int_stack + 18644;
 Libderiv->deriv_classes[6][4][9] = int_stack + 19400;
 Libderiv->deriv_classes[6][5][9] = int_stack + 19820;
 Libderiv->deriv_classes[6][6][9] = int_stack + 20408;
 Libderiv->deriv_classes[6][7][9] = int_stack + 21192;
 Libderiv->deriv_classes[3][4][8] = int_stack + 22200;
 Libderiv->deriv_classes[3][5][8] = int_stack + 22350;
 Libderiv->deriv_classes[3][6][8] = int_stack + 22560;
 Libderiv->deriv_classes[3][7][8] = int_stack + 22840;
 Libderiv->deriv_classes[4][4][8] = int_stack + 23200;
 Libderiv->deriv_classes[4][5][8] = int_stack + 23425;
 Libderiv->deriv_classes[4][6][8] = int_stack + 23740;
 Libderiv->deriv_classes[4][7][8] = int_stack + 24160;
 Libderiv->deriv_classes[5][4][8] = int_stack + 24700;
 Libderiv->deriv_classes[5][5][8] = int_stack + 25015;
 Libderiv->deriv_classes[5][6][8] = int_stack + 25456;
 Libderiv->deriv_classes[5][7][8] = int_stack + 26044;
 Libderiv->deriv_classes[6][4][8] = int_stack + 26800;
 Libderiv->deriv_classes[6][5][8] = int_stack + 27220;
 Libderiv->deriv_classes[6][6][8] = int_stack + 27808;
 Libderiv->deriv_classes[6][7][8] = int_stack + 28592;
 Libderiv->deriv_classes[3][4][7] = int_stack + 29600;
 Libderiv->deriv_classes[3][5][7] = int_stack + 29750;
 Libderiv->deriv_classes[3][6][7] = int_stack + 29960;
 Libderiv->deriv_classes[3][7][7] = int_stack + 30240;
 Libderiv->deriv_classes[4][4][7] = int_stack + 30600;
 Libderiv->deriv_classes[4][5][7] = int_stack + 30825;
 Libderiv->deriv_classes[4][6][7] = int_stack + 31140;
 Libderiv->deriv_classes[4][7][7] = int_stack + 31560;
 Libderiv->deriv_classes[5][4][7] = int_stack + 32100;
 Libderiv->deriv_classes[5][5][7] = int_stack + 32415;
 Libderiv->deriv_classes[5][6][7] = int_stack + 32856;
 Libderiv->deriv_classes[5][7][7] = int_stack + 33444;
 Libderiv->deriv_classes[6][4][7] = int_stack + 34200;
 Libderiv->deriv_classes[6][5][7] = int_stack + 34620;
 Libderiv->deriv_classes[6][6][7] = int_stack + 35208;
 Libderiv->deriv_classes[6][7][7] = int_stack + 35992;
 Libderiv->deriv_classes[3][4][6] = int_stack + 37000;
 Libderiv->deriv_classes[3][5][6] = int_stack + 37150;
 Libderiv->deriv_classes[3][6][6] = int_stack + 37360;
 Libderiv->deriv_classes[3][7][6] = int_stack + 37640;
 Libderiv->deriv_classes[4][4][6] = int_stack + 38000;
 Libderiv->deriv_classes[4][5][6] = int_stack + 38225;
 Libderiv->deriv_classes[4][6][6] = int_stack + 38540;
 Libderiv->deriv_classes[4][7][6] = int_stack + 38960;
 Libderiv->deriv_classes[5][4][6] = int_stack + 39500;
 Libderiv->deriv_classes[5][5][6] = int_stack + 39815;
 Libderiv->deriv_classes[5][6][6] = int_stack + 40256;
 Libderiv->deriv_classes[5][7][6] = int_stack + 40844;
 Libderiv->dvrr_classes[6][4] = int_stack + 41600;
 Libderiv->deriv_classes[6][4][6] = int_stack + 42020;
 Libderiv->dvrr_classes[6][5] = int_stack + 42440;
 Libderiv->deriv_classes[6][5][6] = int_stack + 43028;
 Libderiv->dvrr_classes[6][6] = int_stack + 43616;
 Libderiv->deriv_classes[6][6][6] = int_stack + 44400;
 Libderiv->deriv_classes[6][7][6] = int_stack + 45184;
 Libderiv->deriv_classes[3][4][2] = int_stack + 46192;
 Libderiv->deriv_classes[3][5][2] = int_stack + 46342;
 Libderiv->deriv_classes[3][6][2] = int_stack + 46552;
 Libderiv->deriv_classes[3][7][2] = int_stack + 46832;
 Libderiv->deriv_classes[4][4][2] = int_stack + 47192;
 Libderiv->deriv_classes[4][5][2] = int_stack + 47417;
 Libderiv->deriv_classes[4][6][2] = int_stack + 47732;
 Libderiv->deriv_classes[4][7][2] = int_stack + 48152;
 Libderiv->deriv_classes[5][4][2] = int_stack + 48692;
 Libderiv->deriv_classes[5][5][2] = int_stack + 49007;
 Libderiv->deriv_classes[5][6][2] = int_stack + 49448;
 Libderiv->deriv_classes[5][7][2] = int_stack + 50036;
 Libderiv->deriv_classes[6][4][2] = int_stack + 50792;
 Libderiv->deriv_classes[6][5][2] = int_stack + 51212;
 Libderiv->deriv_classes[6][6][2] = int_stack + 51800;
 Libderiv->deriv_classes[6][7][2] = int_stack + 52584;
 Libderiv->deriv_classes[3][4][1] = int_stack + 53592;
 Libderiv->deriv_classes[3][5][1] = int_stack + 53742;
 Libderiv->deriv_classes[3][6][1] = int_stack + 53952;
 Libderiv->deriv_classes[3][7][1] = int_stack + 54232;
 Libderiv->deriv_classes[4][4][1] = int_stack + 54592;
 Libderiv->deriv_classes[4][5][1] = int_stack + 54817;
 Libderiv->deriv_classes[4][6][1] = int_stack + 55132;
 Libderiv->deriv_classes[4][7][1] = int_stack + 55552;
 Libderiv->deriv_classes[5][4][1] = int_stack + 56092;
 Libderiv->deriv_classes[5][5][1] = int_stack + 56407;
 Libderiv->deriv_classes[5][6][1] = int_stack + 56848;
 Libderiv->deriv_classes[5][7][1] = int_stack + 57436;
 Libderiv->deriv_classes[6][4][1] = int_stack + 58192;
 Libderiv->deriv_classes[6][5][1] = int_stack + 58612;
 Libderiv->deriv_classes[6][6][1] = int_stack + 59200;
 Libderiv->deriv_classes[6][7][1] = int_stack + 59984;
 Libderiv->dvrr_classes[3][4] = int_stack + 60992;
 Libderiv->dvrr_classes[3][5] = int_stack + 61142;
 Libderiv->dvrr_classes[3][6] = int_stack + 61352;
 Libderiv->dvrr_classes[3][7] = int_stack + 61632;
 Libderiv->deriv_classes[3][4][0] = int_stack + 61992;
 Libderiv->deriv_classes[3][5][0] = int_stack + 62142;
 Libderiv->deriv_classes[3][6][0] = int_stack + 62352;
 Libderiv->deriv_classes[3][7][0] = int_stack + 62632;
 Libderiv->dvrr_classes[4][4] = int_stack + 62992;
 Libderiv->dvrr_classes[4][5] = int_stack + 63217;
 Libderiv->dvrr_classes[4][6] = int_stack + 63532;
 Libderiv->dvrr_classes[4][7] = int_stack + 63952;
 Libderiv->deriv_classes[4][4][0] = int_stack + 64492;
 Libderiv->deriv_classes[4][5][0] = int_stack + 64717;
 Libderiv->deriv_classes[4][6][0] = int_stack + 65032;
 Libderiv->deriv_classes[4][7][0] = int_stack + 65452;
 Libderiv->dvrr_classes[5][4] = int_stack + 65992;
 Libderiv->dvrr_classes[5][5] = int_stack + 66307;
 Libderiv->dvrr_classes[5][6] = int_stack + 66748;
 Libderiv->dvrr_classes[5][7] = int_stack + 67336;
 Libderiv->deriv_classes[5][4][0] = int_stack + 68092;
 Libderiv->deriv_classes[5][5][0] = int_stack + 68407;
 Libderiv->deriv_classes[5][6][0] = int_stack + 68848;
 Libderiv->deriv_classes[5][7][0] = int_stack + 69436;
 Libderiv->deriv_classes[6][4][0] = int_stack + 70192;
 Libderiv->deriv_classes[6][5][0] = int_stack + 70612;
 Libderiv->deriv_classes[6][6][0] = int_stack + 71200;
 Libderiv->deriv_classes[6][7][0] = int_stack + 71984;
 memset(int_stack,0,583936);

 Libderiv->dvrr_stack = int_stack + 278654;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ffgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+72992,int_stack+61142,int_stack+60992,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+73442,int_stack+61352,int_stack+61142,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+74072,int_stack+73442,int_stack+72992,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+74972,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60992,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+75422,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61142,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+76052,int_stack+75422,int_stack+74972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72992,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+76952,int_stack+640,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61352,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+77792,int_stack+76952,int_stack+75422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73442,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+79052,int_stack+77792,int_stack+76052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74072,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+74972,int_stack+63217,int_stack+62992,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+75647,int_stack+63532,int_stack+63217,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+76592,int_stack+75647,int_stack+74972,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+77942,int_stack+1225,int_stack+1000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+1540,int_stack+1225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63217,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+80552,int_stack+0,int_stack+77942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74972,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+81902,int_stack+1960,int_stack+1540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63532,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+83162,int_stack+81902,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75647,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+83162,int_stack+80552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76592,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+80552,int_stack+0,int_stack+79052,150);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+85052,int_stack+66307,int_stack+65992,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+77942,int_stack+66748,int_stack+66307,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+85997,int_stack+77942,int_stack+85052,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79265,int_stack+2815,int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65992,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+87887,int_stack+3256,int_stack+2815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66307,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+89210,int_stack+87887,int_stack+79265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85052,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+91100,int_stack+3844,int_stack+3256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66748,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+92864,int_stack+91100,int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77942,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+95510,int_stack+92864,int_stack+89210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85997,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+87887,int_stack+95510,int_stack+0,150);
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+98660,int_stack+87887,int_stack+80552,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+42440,int_stack+41600,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1260,int_stack+43616,int_stack+42440,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+79265,int_stack+1260,int_stack+0,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+5020,int_stack+4600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41600,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81785,int_stack+5608,int_stack+5020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42440,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+107660,int_stack+81785,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3024,int_stack+6392,int_stack+5608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43616,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+110180,int_stack+3024,int_stack+81785, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3024,int_stack+110180,int_stack+107660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79265,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+107660,int_stack+3024,int_stack+95510,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+117110,int_stack+107660,int_stack+87887,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87887,int_stack+7550,int_stack+7400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60992, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88337,int_stack+7760,int_stack+7550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61142, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+88967,int_stack+88337,int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72992, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+89867,int_stack+8040,int_stack+7760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61352, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+90707,int_stack+89867,int_stack+88337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73442, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+91967,int_stack+90707,int_stack+88967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74072, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87887,int_stack+8625,int_stack+8400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88562,int_stack+8940,int_stack+8625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63217, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+89507,int_stack+88562,int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74972, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+93467,int_stack+9360,int_stack+8940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63532, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+94727,int_stack+93467,int_stack+88562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75647, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+107660,int_stack+94727,int_stack+89507, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76592, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+93467,int_stack+107660,int_stack+91967,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87887,int_stack+10215,int_stack+9900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65992, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88832,int_stack+10656,int_stack+10215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66307, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+90155,int_stack+88832,int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85052, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+109910,int_stack+11244,int_stack+10656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66748, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+111674,int_stack+109910,int_stack+88832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77942, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3024,int_stack+111674,int_stack+90155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85997, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+109910,int_stack+3024,int_stack+107660,150);
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+130610,int_stack+109910,int_stack+93467,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+12420,int_stack+12000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41600, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+87887,int_stack+13008,int_stack+12420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42440, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+89651,int_stack+87887,int_stack+107660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+92171,int_stack+13792,int_stack+13008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43616, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+94523,int_stack+92171,int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6174,int_stack+94523,int_stack+89651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79265, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+87887,int_stack+6174,int_stack+3024,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+139610,int_stack+87887,int_stack+109910,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87887,int_stack+14950,int_stack+14800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60992, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88337,int_stack+15160,int_stack+14950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61142, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+88967,int_stack+88337,int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72992, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+89867,int_stack+15440,int_stack+15160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61352, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+90707,int_stack+89867,int_stack+88337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73442, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+91967,int_stack+90707,int_stack+88967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74072, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87887,int_stack+16025,int_stack+15800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88562,int_stack+16340,int_stack+16025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63217, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+89507,int_stack+88562,int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74972, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+93467,int_stack+16760,int_stack+16340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63532, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+94727,int_stack+93467,int_stack+88562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75647, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3024,int_stack+94727,int_stack+89507, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76592, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+93467,int_stack+3024,int_stack+91967,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87887,int_stack+17615,int_stack+17300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65992, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88832,int_stack+18056,int_stack+17615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66307, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+90155,int_stack+88832,int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85052, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+5274,int_stack+18644,int_stack+18056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66748, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7038,int_stack+5274,int_stack+88832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77942, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+9684,int_stack+7038,int_stack+90155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85997, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+107660,int_stack+9684,int_stack+3024,150);
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+153110,int_stack+107660,int_stack+93467,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+19820,int_stack+19400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41600, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4284,int_stack+20408,int_stack+19820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42440, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6048,int_stack+4284,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+87887,int_stack+21192,int_stack+20408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43616, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+90239,int_stack+87887,int_stack+4284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+93767,int_stack+90239,int_stack+6048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79265, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+162110,int_stack+93767,int_stack+9684,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+3024,int_stack+162110,int_stack+107660,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+22350,int_stack+22200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108110,int_stack+22560,int_stack+22350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+108740,int_stack+108110,int_stack+107660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+109640,int_stack+22840,int_stack+22560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+110480,int_stack+109640,int_stack+108110, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+111740,int_stack+110480,int_stack+108740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+23425,int_stack+23200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108335,int_stack+23740,int_stack+23425, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109280,int_stack+108335,int_stack+107660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+113240,int_stack+24160,int_stack+23740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+114500,int_stack+113240,int_stack+108335, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75647, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+162110,int_stack+114500,int_stack+109280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+164360,int_stack+162110,int_stack+111740,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+25015,int_stack+24700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108605,int_stack+25456,int_stack+25015, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109928,int_stack+108605,int_stack+107660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+111818,int_stack+26044,int_stack+25456, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66748, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+113582,int_stack+111818,int_stack+108605, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+87887,int_stack+113582,int_stack+109928, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85997, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+107660,int_stack+87887,int_stack+162110,150);
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+16524,int_stack+107660,int_stack+164360,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+162110,int_stack+27220,int_stack+26800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+163370,int_stack+27808,int_stack+27220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+165134,int_stack+163370,int_stack+162110, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+167654,int_stack+28592,int_stack+27808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+91037,int_stack+167654,int_stack+163370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+167654,int_stack+91037,int_stack+165134, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+171854,int_stack+167654,int_stack+87887,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+181304,int_stack+171854,int_stack+107660,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+29750,int_stack+29600, 0.0, zero_stack, 1.0, int_stack+60992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108110,int_stack+29960,int_stack+29750, 0.0, zero_stack, 1.0, int_stack+61142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+108740,int_stack+108110,int_stack+107660, 0.0, zero_stack, 1.0, int_stack+72992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+109640,int_stack+30240,int_stack+29960, 0.0, zero_stack, 1.0, int_stack+61352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+110480,int_stack+109640,int_stack+108110, 0.0, zero_stack, 1.0, int_stack+73442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+111740,int_stack+110480,int_stack+108740, 0.0, zero_stack, 1.0, int_stack+74072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+30825,int_stack+30600, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108335,int_stack+31140,int_stack+30825, 0.0, zero_stack, 1.0, int_stack+63217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109280,int_stack+108335,int_stack+107660, 0.0, zero_stack, 1.0, int_stack+74972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+113240,int_stack+31560,int_stack+31140, 0.0, zero_stack, 1.0, int_stack+63532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+114500,int_stack+113240,int_stack+108335, 0.0, zero_stack, 1.0, int_stack+75647, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+87887,int_stack+114500,int_stack+109280, 0.0, zero_stack, 1.0, int_stack+76592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+90137,int_stack+87887,int_stack+111740,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+32415,int_stack+32100, 0.0, zero_stack, 1.0, int_stack+65992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108605,int_stack+32856,int_stack+32415, 0.0, zero_stack, 1.0, int_stack+66307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109928,int_stack+108605,int_stack+107660, 0.0, zero_stack, 1.0, int_stack+85052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+111818,int_stack+33444,int_stack+32856, 0.0, zero_stack, 1.0, int_stack+66748, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+113582,int_stack+111818,int_stack+108605, 0.0, zero_stack, 1.0, int_stack+77942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+94637,int_stack+113582,int_stack+109928, 0.0, zero_stack, 1.0, int_stack+85997, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+107660,int_stack+94637,int_stack+87887,150);
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+162110,int_stack+107660,int_stack+90137,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87887,int_stack+34620,int_stack+34200, 0.0, zero_stack, 1.0, int_stack+41600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+89147,int_stack+35208,int_stack+34620, 0.0, zero_stack, 1.0, int_stack+42440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+90911,int_stack+89147,int_stack+87887, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+114410,int_stack+35992,int_stack+35208, 0.0, zero_stack, 1.0, int_stack+43616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+171110,int_stack+114410,int_stack+89147, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+174638,int_stack+171110,int_stack+90911, 0.0, zero_stack, 1.0, int_stack+79265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+25524,int_stack+174638,int_stack+94637,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+194804,int_stack+25524,int_stack+107660,150);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+37150,int_stack+37000, 1.0, int_stack+60992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108110,int_stack+37360,int_stack+37150, 1.0, int_stack+61142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+108740,int_stack+108110,int_stack+107660, 1.0, int_stack+72992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+109640,int_stack+37640,int_stack+37360, 1.0, int_stack+61352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+110480,int_stack+109640,int_stack+108110, 1.0, int_stack+73442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+111740,int_stack+110480,int_stack+108740, 1.0, int_stack+74072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+38225,int_stack+38000, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108335,int_stack+38540,int_stack+38225, 1.0, int_stack+63217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109280,int_stack+108335,int_stack+107660, 1.0, int_stack+74972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+113240,int_stack+38960,int_stack+38540, 1.0, int_stack+63532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+114500,int_stack+113240,int_stack+108335, 1.0, int_stack+75647, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+25524,int_stack+114500,int_stack+109280, 1.0, int_stack+76592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+27774,int_stack+25524,int_stack+111740,150);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+39815,int_stack+39500, 1.0, int_stack+65992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108605,int_stack+40256,int_stack+39815, 1.0, int_stack+66307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109928,int_stack+108605,int_stack+107660, 1.0, int_stack+85052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+111818,int_stack+40844,int_stack+40256, 1.0, int_stack+66748, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+113582,int_stack+111818,int_stack+108605, 1.0, int_stack+77942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+32274,int_stack+113582,int_stack+109928, 1.0, int_stack+85997, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+107660,int_stack+32274,int_stack+25524,150);
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+171110,int_stack+107660,int_stack+27774,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25524,int_stack+43028,int_stack+42020, 1.0, int_stack+41600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26784,int_stack+44400,int_stack+43028, 1.0, int_stack+42440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+28548,int_stack+26784,int_stack+25524, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+114410,int_stack+45184,int_stack+44400, 1.0, int_stack+43616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+35424,int_stack+114410,int_stack+26784, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+38952,int_stack+35424,int_stack+28548, 1.0, int_stack+79265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+87887,int_stack+38952,int_stack+32274,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+25524,int_stack+87887,int_stack+107660,150);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+107660,int_stack+61632,int_stack+61352,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+108500,int_stack+107660,int_stack+73442,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+109760,int_stack+108500,int_stack+74072,10);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+107660,int_stack+63952,int_stack+63532,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+111260,int_stack+107660,int_stack+75647,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+113150,int_stack+111260,int_stack+76592,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+87887,int_stack+113150,int_stack+109760,150);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+111260,int_stack+67336,int_stack+66748,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+92387,int_stack+111260,int_stack+77942,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+95033,int_stack+92387,int_stack+85997,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+39024,int_stack+95033,int_stack+113150,150);
 /*--- compute (fd|gf) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+72992,int_stack+39024,int_stack+87887,150);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+92387,int_stack+46342,int_stack+46192,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+92837,int_stack+46552,int_stack+46342,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+93467,int_stack+92837,int_stack+92387,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+111260,int_stack+46832,int_stack+46552,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+107660,int_stack+111260,int_stack+92837,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+62992,int_stack+107660,int_stack+93467,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+107660,int_stack+47417,int_stack+47192,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+108335,int_stack+47732,int_stack+47417,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+111260,int_stack+108335,int_stack+107660,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+92387,int_stack+48152,int_stack+47732,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+0,int_stack+92387,int_stack+108335,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+92387,int_stack+0,int_stack+111260,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+81992,int_stack+92387,int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+62992,int_stack+49007,int_stack+48692,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+111260,int_stack+49448,int_stack+49007,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+111260,int_stack+62992,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+107660,int_stack+50036,int_stack+49448,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+45774,int_stack+107660,int_stack+111260,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+208304,int_stack+45774,int_stack+0,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+211454,int_stack+208304,int_stack+92387, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (fd|gf) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+218204,int_stack+211454,int_stack+81992, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+81992,int_stack+51212,int_stack+50792,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+83252,int_stack+51800,int_stack+51212,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+85016,int_stack+83252,int_stack+81992,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+92387,int_stack+52584,int_stack+51800,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+45774,int_stack+92387,int_stack+83252,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+49302,int_stack+45774,int_stack+85016,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+227204,int_stack+49302,int_stack+208304, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95033, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+236654,int_stack+227204,int_stack+211454, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+227204,int_stack+53742,int_stack+53592,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+227654,int_stack+53952,int_stack+53742,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+228284,int_stack+227654,int_stack+227204,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+229184,int_stack+54232,int_stack+53952,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+230024,int_stack+229184,int_stack+227654,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+62992,int_stack+230024,int_stack+228284,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+227204,int_stack+54817,int_stack+54592,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+227879,int_stack+55132,int_stack+54817,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+228824,int_stack+227879,int_stack+227204,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+230174,int_stack+55552,int_stack+55132,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+111260,int_stack+230174,int_stack+227879,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+230174,int_stack+111260,int_stack+228824,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+208304,int_stack+230174,int_stack+62992, 0.0, zero_stack, 1.0, int_stack+109760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+62992,int_stack+56407,int_stack+56092,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+111260,int_stack+56848,int_stack+56407,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+232424,int_stack+111260,int_stack+62992,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+234314,int_stack+57436,int_stack+56848,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+92387,int_stack+234314,int_stack+111260,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+212804,int_stack+92387,int_stack+232424,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+45774,int_stack+212804,int_stack+230174, 0.0, zero_stack, 1.0, int_stack+113150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (fd|gf) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+227204,int_stack+45774,int_stack+208304, 0.0, zero_stack, 1.0, int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+208304,int_stack+58612,int_stack+58192,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+209564,int_stack+59200,int_stack+58612,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+92387,int_stack+209564,int_stack+208304,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+52524,int_stack+59984,int_stack+59200,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+54876,int_stack+52524,int_stack+209564,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+208304,int_stack+54876,int_stack+92387,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+52524,int_stack+208304,int_stack+212804, 0.0, zero_stack, 1.0, int_stack+95033, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+250154,int_stack+52524,int_stack+45774, 0.0, zero_stack, 1.0, int_stack+39024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+236204,int_stack+62142,int_stack+61992,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+45774,int_stack+62352,int_stack+62142,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+46404,int_stack+45774,int_stack+236204,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+47304,int_stack+62632,int_stack+62352,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+48144,int_stack+47304,int_stack+45774,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+49404,int_stack+48144,int_stack+46404,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+45774,int_stack+64717,int_stack+64492,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+46449,int_stack+65032,int_stack+64717,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+47394,int_stack+46449,int_stack+45774,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+50904,int_stack+65452,int_stack+65032,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+111260,int_stack+50904,int_stack+46449,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+50904,int_stack+111260,int_stack+47394,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+53154,int_stack+50904,int_stack+49404, 1.0, int_stack+109760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+57654,int_stack+68407,int_stack+68092,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+58599,int_stack+68848,int_stack+68407,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+59922,int_stack+58599,int_stack+57654,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+61812,int_stack+69436,int_stack+68848,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+92387,int_stack+61812,int_stack+58599,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+61812,int_stack+92387,int_stack+59922,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+208304,int_stack+61812,int_stack+50904, 1.0, int_stack+113150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (fd|gf) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+107660,int_stack+208304,int_stack+53154, 1.0, int_stack+87887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+64962,int_stack+70612,int_stack+70192,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+66222,int_stack+71200,int_stack+70612,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+67986,int_stack+66222,int_stack+64962,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+45774,int_stack+71984,int_stack+71200,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+48126,int_stack+45774,int_stack+66222,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+51654,int_stack+48126,int_stack+67986,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+81992,int_stack+51654,int_stack+61812, 1.0, int_stack+95033, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+45774,int_stack+81992,int_stack+208304, 1.0, int_stack+39024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (ff|gf) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+81992,int_stack+117110,int_stack+98660,150);
     Libderiv->ABCD[11] = int_stack + 81992;
 /*--- compute (ff|gf) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+263654,int_stack+139610,int_stack+130610,150);
     Libderiv->ABCD[10] = int_stack + 263654;
 /*--- compute (ff|gf) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+116660,int_stack+3024,int_stack+153110,150);
     Libderiv->ABCD[9] = int_stack + 116660;
 /*--- compute (ff|gf) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+0,int_stack+181304,int_stack+16524,150);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (ff|gf) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+131660,int_stack+194804,int_stack+162110,150);
     Libderiv->ABCD[7] = int_stack + 131660;
 /*--- compute (ff|gf) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+146660,int_stack+25524,int_stack+171110,150);
     Libderiv->ABCD[6] = int_stack + 146660;
 /*--- compute (ff|gf) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+15000,int_stack+236654,int_stack+218204, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 15000;
 /*--- compute (ff|gf) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+30000,int_stack+250154,int_stack+227204, 0.0, zero_stack, 1.0, int_stack+72992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 30000;
 /*--- compute (ff|gf) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+161660,int_stack+45774,int_stack+107660, 1.0, int_stack+72992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 161660;

}
