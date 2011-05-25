#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0gg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|gg) integrals */

void d1hrr_order_f0gg(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[3][8][11] = int_stack + 1000;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1450;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1600;
 Libderiv->deriv_classes[3][6][10] = int_stack + 1810;
 Libderiv->deriv_classes[3][7][10] = int_stack + 2090;
 Libderiv->deriv_classes[3][8][10] = int_stack + 2450;
 Libderiv->deriv_classes[3][4][9] = int_stack + 2900;
 Libderiv->deriv_classes[3][5][9] = int_stack + 3050;
 Libderiv->deriv_classes[3][6][9] = int_stack + 3260;
 Libderiv->deriv_classes[3][7][9] = int_stack + 3540;
 Libderiv->deriv_classes[3][8][9] = int_stack + 3900;
 Libderiv->deriv_classes[3][4][8] = int_stack + 4350;
 Libderiv->deriv_classes[3][5][8] = int_stack + 4500;
 Libderiv->deriv_classes[3][6][8] = int_stack + 4710;
 Libderiv->deriv_classes[3][7][8] = int_stack + 4990;
 Libderiv->deriv_classes[3][8][8] = int_stack + 5350;
 Libderiv->deriv_classes[3][4][7] = int_stack + 5800;
 Libderiv->deriv_classes[3][5][7] = int_stack + 5950;
 Libderiv->deriv_classes[3][6][7] = int_stack + 6160;
 Libderiv->deriv_classes[3][7][7] = int_stack + 6440;
 Libderiv->deriv_classes[3][8][7] = int_stack + 6800;
 Libderiv->dvrr_classes[3][4] = int_stack + 7250;
 Libderiv->deriv_classes[3][4][6] = int_stack + 7400;
 Libderiv->dvrr_classes[3][5] = int_stack + 7550;
 Libderiv->deriv_classes[3][5][6] = int_stack + 7760;
 Libderiv->dvrr_classes[3][6] = int_stack + 7970;
 Libderiv->deriv_classes[3][6][6] = int_stack + 8250;
 Libderiv->dvrr_classes[3][7] = int_stack + 8530;
 Libderiv->deriv_classes[3][7][6] = int_stack + 8890;
 Libderiv->deriv_classes[3][8][6] = int_stack + 9250;
 Libderiv->deriv_classes[3][4][2] = int_stack + 9700;
 Libderiv->deriv_classes[3][5][2] = int_stack + 9850;
 Libderiv->deriv_classes[3][6][2] = int_stack + 10060;
 Libderiv->deriv_classes[3][7][2] = int_stack + 10340;
 Libderiv->deriv_classes[3][8][2] = int_stack + 10700;
 Libderiv->deriv_classes[3][4][1] = int_stack + 11150;
 Libderiv->deriv_classes[3][5][1] = int_stack + 11300;
 Libderiv->deriv_classes[3][6][1] = int_stack + 11510;
 Libderiv->deriv_classes[3][7][1] = int_stack + 11790;
 Libderiv->deriv_classes[3][8][1] = int_stack + 12150;
 Libderiv->deriv_classes[3][4][0] = int_stack + 12600;
 Libderiv->deriv_classes[3][5][0] = int_stack + 12750;
 Libderiv->deriv_classes[3][6][0] = int_stack + 12960;
 Libderiv->deriv_classes[3][7][0] = int_stack + 13240;
 Libderiv->deriv_classes[3][8][0] = int_stack + 13600;
 memset(int_stack,0,112400);

 Libderiv->dvrr_stack = int_stack + 47230;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0gg(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14050,int_stack+7550,int_stack+7250,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+14500,int_stack+7970,int_stack+7550,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15130,int_stack+14500,int_stack+14050,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+16030,int_stack+8530,int_stack+7970,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+16870,int_stack+16030,int_stack+14500,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+18130,int_stack+16870,int_stack+15130,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19630,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7250,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20080,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7550,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20710,int_stack+20080,int_stack+19630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14050,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+21610,int_stack+640,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7970,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+22450,int_stack+21610,int_stack+20080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14500,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+23710,int_stack+22450,int_stack+20710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15130,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+19630,int_stack+1000,int_stack+640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8530,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+25210,int_stack+19630,int_stack+21610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16030,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+19630,int_stack+25210,int_stack+22450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16870,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25210,int_stack+1600,int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7250, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+25660,int_stack+1810,int_stack+1600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7550, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21730,int_stack+25660,int_stack+25210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14050, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+22630,int_stack+2090,int_stack+1810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7970, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+0,int_stack+22630,int_stack+25660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14500, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+25210,int_stack+0,int_stack+21730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15130, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+26710,int_stack+2450,int_stack+2090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8530, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+27790,int_stack+26710,int_stack+22630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16030, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+29470,int_stack+27790,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16870, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3050,int_stack+2900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7250, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+450,int_stack+3260,int_stack+3050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7550, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14050, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1980,int_stack+3540,int_stack+3260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7970, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+26710,int_stack+1980,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+27970,int_stack+26710,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15130, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+3900,int_stack+3540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8530, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+21730,int_stack+0,int_stack+1980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16030, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+21730,int_stack+26710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16870, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26710,int_stack+4500,int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+27160,int_stack+4710,int_stack+4500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21730,int_stack+27160,int_stack+26710, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+22630,int_stack+4990,int_stack+4710, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2100,int_stack+22630,int_stack+27160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3360,int_stack+2100,int_stack+21730, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+26710,int_stack+5350,int_stack+4990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+31570,int_stack+26710,int_stack+22630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+33250,int_stack+31570,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2100,int_stack+5950,int_stack+5800, 0.0, zero_stack, 1.0, int_stack+7250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2550,int_stack+6160,int_stack+5950, 0.0, zero_stack, 1.0, int_stack+7550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31570,int_stack+2550,int_stack+2100, 0.0, zero_stack, 1.0, int_stack+14050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+26710,int_stack+6440,int_stack+6160, 0.0, zero_stack, 1.0, int_stack+7970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21730,int_stack+26710,int_stack+2550, 0.0, zero_stack, 1.0, int_stack+14500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4860,int_stack+21730,int_stack+31570, 0.0, zero_stack, 1.0, int_stack+15130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+31570,int_stack+6800,int_stack+6440, 0.0, zero_stack, 1.0, int_stack+8530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+35350,int_stack+31570,int_stack+26710, 0.0, zero_stack, 1.0, int_stack+16030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+37030,int_stack+35350,int_stack+21730, 0.0, zero_stack, 1.0, int_stack+16870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21730,int_stack+7760,int_stack+7400, 1.0, int_stack+7250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22180,int_stack+8250,int_stack+7760, 1.0, int_stack+7550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22810,int_stack+22180,int_stack+21730, 1.0, int_stack+14050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+35350,int_stack+8890,int_stack+8250, 1.0, int_stack+7970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+26710,int_stack+35350,int_stack+22180, 1.0, int_stack+14500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+31570,int_stack+26710,int_stack+22810, 1.0, int_stack+15130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+14050,int_stack+9250,int_stack+8890, 1.0, int_stack+8530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+21730,int_stack+14050,int_stack+35350, 1.0, int_stack+16030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+14050,int_stack+21730,int_stack+26710, 1.0, int_stack+16870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+26710,int_stack+9850,int_stack+9700,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+27160,int_stack+10060,int_stack+9850,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+21730,int_stack+27160,int_stack+26710,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+22630,int_stack+10340,int_stack+10060,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+2100,int_stack+22630,int_stack+27160,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+35350,int_stack+2100,int_stack+21730,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+26710,int_stack+10700,int_stack+10340,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+16150,int_stack+26710,int_stack+22630,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+6360,int_stack+16150,int_stack+2100,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2100,int_stack+11300,int_stack+11150,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2550,int_stack+11510,int_stack+11300,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16150,int_stack+2550,int_stack+2100,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+17050,int_stack+11790,int_stack+11510,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+26710,int_stack+17050,int_stack+2550,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+21730,int_stack+26710,int_stack+16150,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+2100,int_stack+12150,int_stack+11790,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+8460,int_stack+2100,int_stack+17050,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+10140,int_stack+8460,int_stack+26710,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+26710,int_stack+12750,int_stack+12600,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+27160,int_stack+12960,int_stack+12750,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8460,int_stack+27160,int_stack+26710,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+2100,int_stack+13240,int_stack+12960,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+16150,int_stack+2100,int_stack+27160,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+39130,int_stack+16150,int_stack+8460,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+8460,int_stack+13600,int_stack+13240,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+12240,int_stack+8460,int_stack+2100,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+40630,int_stack+12240,int_stack+16150,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+42730,int_stack+19630,int_stack+23710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18130,10);
     Libderiv->ABCD[11] = int_stack + 42730;
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+44980,int_stack+29470,int_stack+25210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18130, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 44980;
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+23230,int_stack+0,int_stack+27970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18130, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 23230;
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+0,int_stack+33250,int_stack+3360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+2250,int_stack+37030,int_stack+4860, 0.0, zero_stack, 1.0, int_stack+18130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 2250;
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+25480,int_stack+14050,int_stack+31570, 1.0, int_stack+18130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 25480;
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+12240,int_stack+6360,int_stack+35350,10);
     Libderiv->ABCD[2] = int_stack + 12240;
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+4500,int_stack+10140,int_stack+21730,10);
     Libderiv->ABCD[1] = int_stack + 4500;
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+6750,int_stack+40630,int_stack+39130,10);
     Libderiv->ABCD[0] = int_stack + 6750;

}
