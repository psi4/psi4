#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gpgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gp|gg) integrals */

void d1hrr_order_gpgg(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 540;
 Libderiv->deriv_classes[4][7][11] = int_stack + 960;
 Libderiv->deriv_classes[4][8][11] = int_stack + 1500;
 Libderiv->deriv_classes[5][4][11] = int_stack + 2175;
 Libderiv->deriv_classes[5][5][11] = int_stack + 2490;
 Libderiv->deriv_classes[5][6][11] = int_stack + 2931;
 Libderiv->deriv_classes[5][7][11] = int_stack + 3519;
 Libderiv->deriv_classes[5][8][11] = int_stack + 4275;
 Libderiv->deriv_classes[4][4][10] = int_stack + 5220;
 Libderiv->deriv_classes[4][5][10] = int_stack + 5445;
 Libderiv->deriv_classes[4][6][10] = int_stack + 5760;
 Libderiv->deriv_classes[4][7][10] = int_stack + 6180;
 Libderiv->deriv_classes[4][8][10] = int_stack + 6720;
 Libderiv->deriv_classes[5][4][10] = int_stack + 7395;
 Libderiv->deriv_classes[5][5][10] = int_stack + 7710;
 Libderiv->deriv_classes[5][6][10] = int_stack + 8151;
 Libderiv->deriv_classes[5][7][10] = int_stack + 8739;
 Libderiv->deriv_classes[5][8][10] = int_stack + 9495;
 Libderiv->deriv_classes[4][4][9] = int_stack + 10440;
 Libderiv->deriv_classes[4][5][9] = int_stack + 10665;
 Libderiv->deriv_classes[4][6][9] = int_stack + 10980;
 Libderiv->deriv_classes[4][7][9] = int_stack + 11400;
 Libderiv->deriv_classes[4][8][9] = int_stack + 11940;
 Libderiv->deriv_classes[5][4][9] = int_stack + 12615;
 Libderiv->deriv_classes[5][5][9] = int_stack + 12930;
 Libderiv->deriv_classes[5][6][9] = int_stack + 13371;
 Libderiv->deriv_classes[5][7][9] = int_stack + 13959;
 Libderiv->deriv_classes[5][8][9] = int_stack + 14715;
 Libderiv->deriv_classes[4][4][8] = int_stack + 15660;
 Libderiv->deriv_classes[4][5][8] = int_stack + 15885;
 Libderiv->deriv_classes[4][6][8] = int_stack + 16200;
 Libderiv->deriv_classes[4][7][8] = int_stack + 16620;
 Libderiv->deriv_classes[4][8][8] = int_stack + 17160;
 Libderiv->deriv_classes[5][4][8] = int_stack + 17835;
 Libderiv->deriv_classes[5][5][8] = int_stack + 18150;
 Libderiv->deriv_classes[5][6][8] = int_stack + 18591;
 Libderiv->deriv_classes[5][7][8] = int_stack + 19179;
 Libderiv->deriv_classes[5][8][8] = int_stack + 19935;
 Libderiv->deriv_classes[4][4][7] = int_stack + 20880;
 Libderiv->deriv_classes[4][5][7] = int_stack + 21105;
 Libderiv->deriv_classes[4][6][7] = int_stack + 21420;
 Libderiv->deriv_classes[4][7][7] = int_stack + 21840;
 Libderiv->deriv_classes[4][8][7] = int_stack + 22380;
 Libderiv->deriv_classes[5][4][7] = int_stack + 23055;
 Libderiv->deriv_classes[5][5][7] = int_stack + 23370;
 Libderiv->deriv_classes[5][6][7] = int_stack + 23811;
 Libderiv->deriv_classes[5][7][7] = int_stack + 24399;
 Libderiv->deriv_classes[5][8][7] = int_stack + 25155;
 Libderiv->deriv_classes[4][4][6] = int_stack + 26100;
 Libderiv->deriv_classes[4][5][6] = int_stack + 26325;
 Libderiv->deriv_classes[4][6][6] = int_stack + 26640;
 Libderiv->deriv_classes[4][7][6] = int_stack + 27060;
 Libderiv->deriv_classes[4][8][6] = int_stack + 27600;
 Libderiv->dvrr_classes[5][4] = int_stack + 28275;
 Libderiv->deriv_classes[5][4][6] = int_stack + 28590;
 Libderiv->dvrr_classes[5][5] = int_stack + 28905;
 Libderiv->deriv_classes[5][5][6] = int_stack + 29346;
 Libderiv->dvrr_classes[5][6] = int_stack + 29787;
 Libderiv->deriv_classes[5][6][6] = int_stack + 30375;
 Libderiv->dvrr_classes[5][7] = int_stack + 30963;
 Libderiv->deriv_classes[5][7][6] = int_stack + 31719;
 Libderiv->deriv_classes[5][8][6] = int_stack + 32475;
 Libderiv->deriv_classes[4][4][2] = int_stack + 33420;
 Libderiv->deriv_classes[4][5][2] = int_stack + 33645;
 Libderiv->deriv_classes[4][6][2] = int_stack + 33960;
 Libderiv->deriv_classes[4][7][2] = int_stack + 34380;
 Libderiv->deriv_classes[4][8][2] = int_stack + 34920;
 Libderiv->deriv_classes[5][4][2] = int_stack + 35595;
 Libderiv->deriv_classes[5][5][2] = int_stack + 35910;
 Libderiv->deriv_classes[5][6][2] = int_stack + 36351;
 Libderiv->deriv_classes[5][7][2] = int_stack + 36939;
 Libderiv->deriv_classes[5][8][2] = int_stack + 37695;
 Libderiv->deriv_classes[4][4][1] = int_stack + 38640;
 Libderiv->deriv_classes[4][5][1] = int_stack + 38865;
 Libderiv->deriv_classes[4][6][1] = int_stack + 39180;
 Libderiv->deriv_classes[4][7][1] = int_stack + 39600;
 Libderiv->deriv_classes[4][8][1] = int_stack + 40140;
 Libderiv->deriv_classes[5][4][1] = int_stack + 40815;
 Libderiv->deriv_classes[5][5][1] = int_stack + 41130;
 Libderiv->deriv_classes[5][6][1] = int_stack + 41571;
 Libderiv->deriv_classes[5][7][1] = int_stack + 42159;
 Libderiv->deriv_classes[5][8][1] = int_stack + 42915;
 Libderiv->dvrr_classes[4][4] = int_stack + 43860;
 Libderiv->dvrr_classes[4][5] = int_stack + 44085;
 Libderiv->dvrr_classes[4][6] = int_stack + 44400;
 Libderiv->dvrr_classes[4][7] = int_stack + 44820;
 Libderiv->dvrr_classes[4][8] = int_stack + 45360;
 Libderiv->deriv_classes[4][4][0] = int_stack + 46035;
 Libderiv->deriv_classes[4][5][0] = int_stack + 46260;
 Libderiv->deriv_classes[4][6][0] = int_stack + 46575;
 Libderiv->deriv_classes[4][7][0] = int_stack + 46995;
 Libderiv->deriv_classes[4][8][0] = int_stack + 47535;
 Libderiv->deriv_classes[5][4][0] = int_stack + 48210;
 Libderiv->deriv_classes[5][5][0] = int_stack + 48525;
 Libderiv->deriv_classes[5][6][0] = int_stack + 48966;
 Libderiv->deriv_classes[5][7][0] = int_stack + 49554;
 Libderiv->deriv_classes[5][8][0] = int_stack + 50310;
 memset(int_stack,0,410040);

 Libderiv->dvrr_stack = int_stack + 122526;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gpgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51255,int_stack+44085,int_stack+43860,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+51930,int_stack+44400,int_stack+44085,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+52875,int_stack+51930,int_stack+51255,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+54225,int_stack+44820,int_stack+44400,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+55485,int_stack+54225,int_stack+51930,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+57375,int_stack+55485,int_stack+52875,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+59625,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43860,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60300,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44085,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+61245,int_stack+60300,int_stack+59625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51255,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+62595,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44400,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+63855,int_stack+62595,int_stack+60300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51930,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+65745,int_stack+63855,int_stack+61245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52875,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+59625,int_stack+1500,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44820,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+67995,int_stack+59625,int_stack+62595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54225,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+59625,int_stack+67995,int_stack+63855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55485,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+67995,int_stack+59625,int_stack+65745, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57375,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+59625,int_stack+28905,int_stack+28275,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+60570,int_stack+29787,int_stack+28905,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+61893,int_stack+60570,int_stack+59625,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+63783,int_stack+30963,int_stack+29787,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+71370,int_stack+63783,int_stack+60570,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+74016,int_stack+71370,int_stack+61893,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65547,int_stack+2490,int_stack+2175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28275,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66492,int_stack+2931,int_stack+2490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28905,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+66492,int_stack+65547, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59625,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+77166,int_stack+3519,int_stack+2931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29787,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+78930,int_stack+77166,int_stack+66492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60570,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+81576,int_stack+78930,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61893,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+4275,int_stack+3519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30963,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+84726,int_stack+0,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63783,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+84726,int_stack+78930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+84726,int_stack+0,int_stack+81576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74016,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+5445,int_stack+5220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43860, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+675,int_stack+5760,int_stack+5445, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44085, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1620,int_stack+675,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51255, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2970,int_stack+6180,int_stack+5760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44400, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+4230,int_stack+2970,int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51930, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+77166,int_stack+4230,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52875, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+6720,int_stack+6180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44820, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+79416,int_stack+0,int_stack+2970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54225, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+79416,int_stack+4230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55485, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+79416,int_stack+0,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57375, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+77166,int_stack+7710,int_stack+7395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28275, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+8151,int_stack+7710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28905, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1323,int_stack+0,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59625, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+77166,int_stack+8739,int_stack+8151, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29787, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3213,int_stack+77166,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60570, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+89451,int_stack+3213,int_stack+1323, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61893, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+9495,int_stack+8739, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30963, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+5859,int_stack+0,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63783, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+92601,int_stack+5859,int_stack+3213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+0,int_stack+92601,int_stack+89451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74016, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+89451,int_stack+10665,int_stack+10440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43860, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+90126,int_stack+10980,int_stack+10665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44085, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+91071,int_stack+90126,int_stack+89451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51255, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+92421,int_stack+11400,int_stack+10980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44400, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+93681,int_stack+92421,int_stack+90126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51930, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+77166,int_stack+93681,int_stack+91071, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52875, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+89451,int_stack+11940,int_stack+11400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44820, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+4725,int_stack+89451,int_stack+92421, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54225, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+89451,int_stack+4725,int_stack+93681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55485, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+4725,int_stack+89451,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57375, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+77166,int_stack+12930,int_stack+12615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28275, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+89451,int_stack+13371,int_stack+12930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28905, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+90774,int_stack+89451,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59625, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+77166,int_stack+13959,int_stack+13371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29787, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+92664,int_stack+77166,int_stack+89451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60570, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8100,int_stack+92664,int_stack+90774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61893, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+89451,int_stack+14715,int_stack+13959, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30963, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+11250,int_stack+89451,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63783, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+95310,int_stack+11250,int_stack+92664, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+89451,int_stack+95310,int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74016, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+15885,int_stack+15660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8775,int_stack+16200,int_stack+15885, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9720,int_stack+8775,int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+11070,int_stack+16620,int_stack+16200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+12330,int_stack+11070,int_stack+8775, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+77166,int_stack+12330,int_stack+9720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+8100,int_stack+17160,int_stack+16620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+14220,int_stack+8100,int_stack+11070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+8100,int_stack+14220,int_stack+12330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+11250,int_stack+8100,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+77166,int_stack+18150,int_stack+17835, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8100,int_stack+18591,int_stack+18150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14625,int_stack+8100,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+77166,int_stack+19179,int_stack+18591, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+16515,int_stack+77166,int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8100,int_stack+16515,int_stack+14625, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61893, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+94176,int_stack+19935,int_stack+19179, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+96444,int_stack+94176,int_stack+77166, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63783, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+99972,int_stack+96444,int_stack+16515, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+94176,int_stack+99972,int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+21105,int_stack+20880, 0.0, zero_stack, 1.0, int_stack+43860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8775,int_stack+21420,int_stack+21105, 0.0, zero_stack, 1.0, int_stack+44085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9720,int_stack+8775,int_stack+8100, 0.0, zero_stack, 1.0, int_stack+51255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+77166,int_stack+21840,int_stack+21420, 0.0, zero_stack, 1.0, int_stack+44400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+98901,int_stack+77166,int_stack+8775, 0.0, zero_stack, 1.0, int_stack+51930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+100791,int_stack+98901,int_stack+9720, 0.0, zero_stack, 1.0, int_stack+52875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+8100,int_stack+22380,int_stack+21840, 0.0, zero_stack, 1.0, int_stack+44820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+14625,int_stack+8100,int_stack+77166, 0.0, zero_stack, 1.0, int_stack+54225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+8100,int_stack+14625,int_stack+98901, 0.0, zero_stack, 1.0, int_stack+55485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+14625,int_stack+8100,int_stack+100791, 0.0, zero_stack, 1.0, int_stack+57375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+23370,int_stack+23055, 0.0, zero_stack, 1.0, int_stack+28275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9045,int_stack+23811,int_stack+23370, 0.0, zero_stack, 1.0, int_stack+28905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+98901,int_stack+9045,int_stack+8100, 0.0, zero_stack, 1.0, int_stack+59625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+100791,int_stack+24399,int_stack+23811, 0.0, zero_stack, 1.0, int_stack+29787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+18000,int_stack+100791,int_stack+9045, 0.0, zero_stack, 1.0, int_stack+60570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8100,int_stack+18000,int_stack+98901, 0.0, zero_stack, 1.0, int_stack+61893, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+102555,int_stack+25155,int_stack+24399, 0.0, zero_stack, 1.0, int_stack+30963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+20646,int_stack+102555,int_stack+100791, 0.0, zero_stack, 1.0, int_stack+63783, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+98901,int_stack+20646,int_stack+18000, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+18000,int_stack+98901,int_stack+8100, 0.0, zero_stack, 1.0, int_stack+74016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+26325,int_stack+26100, 1.0, int_stack+43860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8775,int_stack+26640,int_stack+26325, 1.0, int_stack+44085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9720,int_stack+8775,int_stack+8100, 1.0, int_stack+51255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+98901,int_stack+27060,int_stack+26640, 1.0, int_stack+44400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+100161,int_stack+98901,int_stack+8775, 1.0, int_stack+51930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+77166,int_stack+100161,int_stack+9720, 1.0, int_stack+52875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+51255,int_stack+27600,int_stack+27060, 1.0, int_stack+44820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+8100,int_stack+51255,int_stack+98901, 1.0, int_stack+54225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+22725,int_stack+8100,int_stack+100161, 1.0, int_stack+55485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+98901,int_stack+22725,int_stack+77166, 1.0, int_stack+57375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+77166,int_stack+29346,int_stack+28590, 1.0, int_stack+28275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22725,int_stack+30375,int_stack+29346, 1.0, int_stack+28905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24048,int_stack+22725,int_stack+77166, 1.0, int_stack+59625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+77166,int_stack+31719,int_stack+30375, 1.0, int_stack+29787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+25938,int_stack+77166,int_stack+22725, 1.0, int_stack+60570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8100,int_stack+25938,int_stack+24048, 1.0, int_stack+61893, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+22725,int_stack+32475,int_stack+31719, 1.0, int_stack+30963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+59625,int_stack+22725,int_stack+77166, 1.0, int_stack+63783, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+63153,int_stack+59625,int_stack+25938, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+22725,int_stack+63153,int_stack+8100, 1.0, int_stack+74016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+8100,int_stack+45360,int_stack+44820,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+71370,int_stack+8100,int_stack+54225,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+8100,int_stack+71370,int_stack+55485,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+71370,int_stack+8100,int_stack+57375,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+33645,int_stack+33420,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8775,int_stack+33960,int_stack+33645,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9720,int_stack+8775,int_stack+8100,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+74745,int_stack+34380,int_stack+33960,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+76005,int_stack+74745,int_stack+8775,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+27450,int_stack+76005,int_stack+9720,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+8100,int_stack+34920,int_stack+34380,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+29700,int_stack+8100,int_stack+74745,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+8100,int_stack+29700,int_stack+76005,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+29700,int_stack+8100,int_stack+27450,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27450,int_stack+35910,int_stack+35595,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8100,int_stack+36351,int_stack+35910,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+33075,int_stack+8100,int_stack+27450,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+27450,int_stack+36939,int_stack+36351,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+74745,int_stack+27450,int_stack+8100,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+8100,int_stack+74745,int_stack+33075,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+33075,int_stack+37695,int_stack+36939,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+51255,int_stack+33075,int_stack+27450,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+33075,int_stack+51255,int_stack+74745,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+51255,int_stack+33075,int_stack+8100,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+38865,int_stack+38640,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8775,int_stack+39180,int_stack+38865,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9720,int_stack+8775,int_stack+8100,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+33075,int_stack+39600,int_stack+39180,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+34335,int_stack+33075,int_stack+8775,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+27450,int_stack+34335,int_stack+9720,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+8100,int_stack+40140,int_stack+39600,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+36225,int_stack+8100,int_stack+33075,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+8100,int_stack+36225,int_stack+34335,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+33075,int_stack+8100,int_stack+27450,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27450,int_stack+41130,int_stack+40815,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8100,int_stack+41571,int_stack+41130,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+36450,int_stack+8100,int_stack+27450,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+27450,int_stack+42159,int_stack+41571,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+38340,int_stack+27450,int_stack+8100,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+8100,int_stack+38340,int_stack+36450,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+74745,int_stack+42915,int_stack+42159,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+40986,int_stack+74745,int_stack+27450,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+74745,int_stack+40986,int_stack+38340,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+36450,int_stack+74745,int_stack+8100,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+46260,int_stack+46035,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8775,int_stack+46575,int_stack+46260,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9720,int_stack+8775,int_stack+8100,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+74745,int_stack+46995,int_stack+46575,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+76005,int_stack+74745,int_stack+8775,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+27450,int_stack+76005,int_stack+9720,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+8100,int_stack+47535,int_stack+46995,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+41175,int_stack+8100,int_stack+74745,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+8100,int_stack+41175,int_stack+76005,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+41175,int_stack+8100,int_stack+27450,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27450,int_stack+48525,int_stack+48210,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8100,int_stack+48966,int_stack+48525,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+44550,int_stack+8100,int_stack+27450,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+27450,int_stack+49554,int_stack+48966,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+46440,int_stack+27450,int_stack+8100,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+8100,int_stack+46440,int_stack+44550,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+74745,int_stack+50310,int_stack+49554,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+55980,int_stack+74745,int_stack+27450,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+74745,int_stack+55980,int_stack+46440,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+55980,int_stack+74745,int_stack+8100,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+102276,int_stack+84726,int_stack+67995,225);
     Libderiv->ABCD[11] = int_stack + 102276;
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+60705,int_stack+0,int_stack+79416,225);
     Libderiv->ABCD[10] = int_stack + 60705;
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+74745,int_stack+89451,int_stack+4725,225);
     Libderiv->ABCD[9] = int_stack + 74745;
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+94176,int_stack+11250,225);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+84870,int_stack+18000,int_stack+14625,225);
     Libderiv->ABCD[7] = int_stack + 84870;
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+10125,int_stack+22725,int_stack+98901,225);
     Libderiv->ABCD[6] = int_stack + 10125;
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+112401,int_stack+51255,int_stack+29700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 112401;
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+20250,int_stack+36450,int_stack+33075, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 20250;
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+30375,int_stack+55980,int_stack+41175, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 30375;

}
