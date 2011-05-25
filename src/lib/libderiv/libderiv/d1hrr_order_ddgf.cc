#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ddgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|gf) integrals */

void d1hrr_order_ddgf(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[4][4][11] = int_stack + 1600;
 Libderiv->deriv_classes[4][5][11] = int_stack + 1825;
 Libderiv->deriv_classes[4][6][11] = int_stack + 2140;
 Libderiv->deriv_classes[4][7][11] = int_stack + 2560;
 Libderiv->deriv_classes[2][4][10] = int_stack + 3100;
 Libderiv->deriv_classes[2][5][10] = int_stack + 3190;
 Libderiv->deriv_classes[2][6][10] = int_stack + 3316;
 Libderiv->deriv_classes[2][7][10] = int_stack + 3484;
 Libderiv->deriv_classes[3][4][10] = int_stack + 3700;
 Libderiv->deriv_classes[3][5][10] = int_stack + 3850;
 Libderiv->deriv_classes[3][6][10] = int_stack + 4060;
 Libderiv->deriv_classes[3][7][10] = int_stack + 4340;
 Libderiv->deriv_classes[4][4][10] = int_stack + 4700;
 Libderiv->deriv_classes[4][5][10] = int_stack + 4925;
 Libderiv->deriv_classes[4][6][10] = int_stack + 5240;
 Libderiv->deriv_classes[4][7][10] = int_stack + 5660;
 Libderiv->deriv_classes[2][4][9] = int_stack + 6200;
 Libderiv->deriv_classes[2][5][9] = int_stack + 6290;
 Libderiv->deriv_classes[2][6][9] = int_stack + 6416;
 Libderiv->deriv_classes[2][7][9] = int_stack + 6584;
 Libderiv->deriv_classes[3][4][9] = int_stack + 6800;
 Libderiv->deriv_classes[3][5][9] = int_stack + 6950;
 Libderiv->deriv_classes[3][6][9] = int_stack + 7160;
 Libderiv->deriv_classes[3][7][9] = int_stack + 7440;
 Libderiv->deriv_classes[4][4][9] = int_stack + 7800;
 Libderiv->deriv_classes[4][5][9] = int_stack + 8025;
 Libderiv->deriv_classes[4][6][9] = int_stack + 8340;
 Libderiv->deriv_classes[4][7][9] = int_stack + 8760;
 Libderiv->deriv_classes[2][4][8] = int_stack + 9300;
 Libderiv->deriv_classes[2][5][8] = int_stack + 9390;
 Libderiv->deriv_classes[2][6][8] = int_stack + 9516;
 Libderiv->deriv_classes[2][7][8] = int_stack + 9684;
 Libderiv->deriv_classes[3][4][8] = int_stack + 9900;
 Libderiv->deriv_classes[3][5][8] = int_stack + 10050;
 Libderiv->deriv_classes[3][6][8] = int_stack + 10260;
 Libderiv->deriv_classes[3][7][8] = int_stack + 10540;
 Libderiv->deriv_classes[4][4][8] = int_stack + 10900;
 Libderiv->deriv_classes[4][5][8] = int_stack + 11125;
 Libderiv->deriv_classes[4][6][8] = int_stack + 11440;
 Libderiv->deriv_classes[4][7][8] = int_stack + 11860;
 Libderiv->deriv_classes[2][4][7] = int_stack + 12400;
 Libderiv->deriv_classes[2][5][7] = int_stack + 12490;
 Libderiv->deriv_classes[2][6][7] = int_stack + 12616;
 Libderiv->deriv_classes[2][7][7] = int_stack + 12784;
 Libderiv->deriv_classes[3][4][7] = int_stack + 13000;
 Libderiv->deriv_classes[3][5][7] = int_stack + 13150;
 Libderiv->deriv_classes[3][6][7] = int_stack + 13360;
 Libderiv->deriv_classes[3][7][7] = int_stack + 13640;
 Libderiv->deriv_classes[4][4][7] = int_stack + 14000;
 Libderiv->deriv_classes[4][5][7] = int_stack + 14225;
 Libderiv->deriv_classes[4][6][7] = int_stack + 14540;
 Libderiv->deriv_classes[4][7][7] = int_stack + 14960;
 Libderiv->deriv_classes[2][4][6] = int_stack + 15500;
 Libderiv->deriv_classes[2][5][6] = int_stack + 15590;
 Libderiv->deriv_classes[2][6][6] = int_stack + 15716;
 Libderiv->deriv_classes[2][7][6] = int_stack + 15884;
 Libderiv->deriv_classes[3][4][6] = int_stack + 16100;
 Libderiv->deriv_classes[3][5][6] = int_stack + 16250;
 Libderiv->deriv_classes[3][6][6] = int_stack + 16460;
 Libderiv->deriv_classes[3][7][6] = int_stack + 16740;
 Libderiv->dvrr_classes[4][4] = int_stack + 17100;
 Libderiv->deriv_classes[4][4][6] = int_stack + 17325;
 Libderiv->dvrr_classes[4][5] = int_stack + 17550;
 Libderiv->deriv_classes[4][5][6] = int_stack + 17865;
 Libderiv->dvrr_classes[4][6] = int_stack + 18180;
 Libderiv->deriv_classes[4][6][6] = int_stack + 18600;
 Libderiv->deriv_classes[4][7][6] = int_stack + 19020;
 Libderiv->deriv_classes[2][4][2] = int_stack + 19560;
 Libderiv->deriv_classes[2][5][2] = int_stack + 19650;
 Libderiv->deriv_classes[2][6][2] = int_stack + 19776;
 Libderiv->deriv_classes[2][7][2] = int_stack + 19944;
 Libderiv->deriv_classes[3][4][2] = int_stack + 20160;
 Libderiv->deriv_classes[3][5][2] = int_stack + 20310;
 Libderiv->deriv_classes[3][6][2] = int_stack + 20520;
 Libderiv->deriv_classes[3][7][2] = int_stack + 20800;
 Libderiv->deriv_classes[4][4][2] = int_stack + 21160;
 Libderiv->deriv_classes[4][5][2] = int_stack + 21385;
 Libderiv->deriv_classes[4][6][2] = int_stack + 21700;
 Libderiv->deriv_classes[4][7][2] = int_stack + 22120;
 Libderiv->deriv_classes[2][4][1] = int_stack + 22660;
 Libderiv->deriv_classes[2][5][1] = int_stack + 22750;
 Libderiv->deriv_classes[2][6][1] = int_stack + 22876;
 Libderiv->deriv_classes[2][7][1] = int_stack + 23044;
 Libderiv->deriv_classes[3][4][1] = int_stack + 23260;
 Libderiv->deriv_classes[3][5][1] = int_stack + 23410;
 Libderiv->deriv_classes[3][6][1] = int_stack + 23620;
 Libderiv->deriv_classes[3][7][1] = int_stack + 23900;
 Libderiv->deriv_classes[4][4][1] = int_stack + 24260;
 Libderiv->deriv_classes[4][5][1] = int_stack + 24485;
 Libderiv->deriv_classes[4][6][1] = int_stack + 24800;
 Libderiv->deriv_classes[4][7][1] = int_stack + 25220;
 Libderiv->dvrr_classes[2][4] = int_stack + 25760;
 Libderiv->dvrr_classes[2][5] = int_stack + 25850;
 Libderiv->dvrr_classes[2][6] = int_stack + 25976;
 Libderiv->dvrr_classes[2][7] = int_stack + 26144;
 Libderiv->deriv_classes[2][4][0] = int_stack + 26360;
 Libderiv->deriv_classes[2][5][0] = int_stack + 26450;
 Libderiv->deriv_classes[2][6][0] = int_stack + 26576;
 Libderiv->deriv_classes[2][7][0] = int_stack + 26744;
 Libderiv->dvrr_classes[3][4] = int_stack + 26960;
 Libderiv->dvrr_classes[3][5] = int_stack + 27110;
 Libderiv->dvrr_classes[3][6] = int_stack + 27320;
 Libderiv->dvrr_classes[3][7] = int_stack + 27600;
 Libderiv->deriv_classes[3][4][0] = int_stack + 27960;
 Libderiv->deriv_classes[3][5][0] = int_stack + 28110;
 Libderiv->deriv_classes[3][6][0] = int_stack + 28320;
 Libderiv->deriv_classes[3][7][0] = int_stack + 28600;
 Libderiv->deriv_classes[4][4][0] = int_stack + 28960;
 Libderiv->deriv_classes[4][5][0] = int_stack + 29185;
 Libderiv->deriv_classes[4][6][0] = int_stack + 29500;
 Libderiv->deriv_classes[4][7][0] = int_stack + 29920;
 memset(int_stack,0,243680);

 Libderiv->dvrr_stack = int_stack + 90601;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ddgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+30460,int_stack+25850,int_stack+25760,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+30730,int_stack+25976,int_stack+25850,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+31108,int_stack+30730,int_stack+30460,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31648,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25760,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+31918,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25850,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+32296,int_stack+31918,int_stack+31648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30460,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+32836,int_stack+384,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25976,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+33340,int_stack+32836,int_stack+31918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30730,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+34096,int_stack+33340,int_stack+32296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31108,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+31648,int_stack+27110,int_stack+26960,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+32098,int_stack+27320,int_stack+27110,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+32728,int_stack+32098,int_stack+31648,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33628,int_stack+750,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26960,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+960,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27110,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+34996,int_stack+0,int_stack+33628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31648,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+35896,int_stack+1240,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27320,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+36736,int_stack+35896,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32098,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+36736,int_stack+34996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32728,10);
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+34996,int_stack+0,int_stack+34096,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+33628,int_stack+17550,int_stack+17100,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+37696,int_stack+18180,int_stack+17550,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+38641,int_stack+37696,int_stack+33628,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34303,int_stack+1825,int_stack+1600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17100,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+39991,int_stack+2140,int_stack+1825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17550,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+40936,int_stack+39991,int_stack+34303, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33628,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+42286,int_stack+2560,int_stack+2140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18180,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+43546,int_stack+42286,int_stack+39991, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37696,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+45436,int_stack+43546,int_stack+40936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38641,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+39991,int_stack+45436,int_stack+0,150);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3190,int_stack+3100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25760, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+270,int_stack+3316,int_stack+3190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25850, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+648,int_stack+270,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30460, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1188,int_stack+3484,int_stack+3316, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25976, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1692,int_stack+1188,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30730, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2448,int_stack+1692,int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31108, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3850,int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26960, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+450,int_stack+4060,int_stack+3850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27110, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31648, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+44491,int_stack+4340,int_stack+4060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27320, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3348,int_stack+44491,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32098, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+44491,int_stack+3348,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32728, 0.0, zero_stack,10);
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+45991,int_stack+44491,int_stack+2448,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+4925,int_stack+4700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17100, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+675,int_stack+5240,int_stack+4925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17550, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1620,int_stack+675,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33628, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2970,int_stack+5660,int_stack+5240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18180, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+4230,int_stack+2970,int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37696, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+48691,int_stack+4230,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38641, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+48691,int_stack+44491,150);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+44491,int_stack+6290,int_stack+6200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25760, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+44761,int_stack+6416,int_stack+6290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25850, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45139,int_stack+44761,int_stack+44491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30460, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+48691,int_stack+6584,int_stack+6416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25976, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+49195,int_stack+48691,int_stack+44761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30730, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+49951,int_stack+49195,int_stack+45139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31108, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48691,int_stack+6950,int_stack+6800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26960, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49141,int_stack+7160,int_stack+6950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27110, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+44491,int_stack+49141,int_stack+48691, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31648, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4500,int_stack+7440,int_stack+7160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27320, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+5340,int_stack+4500,int_stack+49141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32098, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+50851,int_stack+5340,int_stack+44491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32728, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4500,int_stack+50851,int_stack+49951,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+44491,int_stack+8025,int_stack+7800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17100, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+48691,int_stack+8340,int_stack+8025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17550, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+52351,int_stack+48691,int_stack+44491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33628, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+44491,int_stack+8760,int_stack+8340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18180, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7200,int_stack+44491,int_stack+48691, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37696, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+53701,int_stack+7200,int_stack+52351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38641, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+55951,int_stack+53701,int_stack+50851,150);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7200,int_stack+9390,int_stack+9300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7470,int_stack+9516,int_stack+9390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7848,int_stack+7470,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8388,int_stack+9684,int_stack+9516, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8892,int_stack+8388,int_stack+7470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+48691,int_stack+8892,int_stack+7848, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7200,int_stack+10050,int_stack+9900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7650,int_stack+10260,int_stack+10050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8280,int_stack+7650,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9180,int_stack+10540,int_stack+10260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+49591,int_stack+9180,int_stack+7650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+44491,int_stack+49591,int_stack+8280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+49591,int_stack+44491,int_stack+48691,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48691,int_stack+11125,int_stack+10900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52291,int_stack+11440,int_stack+11125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+53236,int_stack+52291,int_stack+48691, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+54586,int_stack+11860,int_stack+11440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7200,int_stack+54586,int_stack+52291, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+9090,int_stack+7200,int_stack+53236, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+60451,int_stack+9090,int_stack+44491,150);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+44491,int_stack+12490,int_stack+12400, 0.0, zero_stack, 1.0, int_stack+25760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+44761,int_stack+12616,int_stack+12490, 0.0, zero_stack, 1.0, int_stack+25850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45139,int_stack+44761,int_stack+44491, 0.0, zero_stack, 1.0, int_stack+30460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+7200,int_stack+12784,int_stack+12616, 0.0, zero_stack, 1.0, int_stack+25976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7704,int_stack+7200,int_stack+44761, 0.0, zero_stack, 1.0, int_stack+30730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+48691,int_stack+7704,int_stack+45139, 0.0, zero_stack, 1.0, int_stack+31108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7200,int_stack+13150,int_stack+13000, 0.0, zero_stack, 1.0, int_stack+26960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7650,int_stack+13360,int_stack+13150, 0.0, zero_stack, 1.0, int_stack+27110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8280,int_stack+7650,int_stack+7200, 0.0, zero_stack, 1.0, int_stack+31648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9180,int_stack+13640,int_stack+13360, 0.0, zero_stack, 1.0, int_stack+27320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10020,int_stack+9180,int_stack+7650, 0.0, zero_stack, 1.0, int_stack+32098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+44491,int_stack+10020,int_stack+8280, 0.0, zero_stack, 1.0, int_stack+32728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7200,int_stack+44491,int_stack+48691,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48691,int_stack+14225,int_stack+14000, 0.0, zero_stack, 1.0, int_stack+17100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9900,int_stack+14540,int_stack+14225, 0.0, zero_stack, 1.0, int_stack+17550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10845,int_stack+9900,int_stack+48691, 0.0, zero_stack, 1.0, int_stack+33628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+12195,int_stack+14960,int_stack+14540, 0.0, zero_stack, 1.0, int_stack+18180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13455,int_stack+12195,int_stack+9900, 0.0, zero_stack, 1.0, int_stack+37696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+52291,int_stack+13455,int_stack+10845, 0.0, zero_stack, 1.0, int_stack+38641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+9900,int_stack+52291,int_stack+44491,150);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+44491,int_stack+15590,int_stack+15500, 1.0, int_stack+25760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+44761,int_stack+15716,int_stack+15590, 1.0, int_stack+25850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45139,int_stack+44761,int_stack+44491, 1.0, int_stack+30460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+52291,int_stack+15884,int_stack+15716, 1.0, int_stack+25976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+52795,int_stack+52291,int_stack+44761, 1.0, int_stack+30730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+48691,int_stack+52795,int_stack+45139, 1.0, int_stack+31108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52291,int_stack+16250,int_stack+16100, 1.0, int_stack+26960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52741,int_stack+16460,int_stack+16250, 1.0, int_stack+27110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+53371,int_stack+52741,int_stack+52291, 1.0, int_stack+31648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+54271,int_stack+16740,int_stack+16460, 1.0, int_stack+27320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+44491,int_stack+54271,int_stack+52741, 1.0, int_stack+32098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+54271,int_stack+44491,int_stack+53371, 1.0, int_stack+32728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14400,int_stack+54271,int_stack+48691,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48691,int_stack+17865,int_stack+17325, 1.0, int_stack+17100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+44491,int_stack+18600,int_stack+17865, 1.0, int_stack+17550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+52291,int_stack+44491,int_stack+48691, 1.0, int_stack+33628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+33628,int_stack+19020,int_stack+18600, 1.0, int_stack+18180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+17100,int_stack+33628,int_stack+44491, 1.0, int_stack+37696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+64951,int_stack+17100,int_stack+52291, 1.0, int_stack+38641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+67201,int_stack+64951,int_stack+54271,150);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+64951,int_stack+26144,int_stack+25976,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+65455,int_stack+64951,int_stack+30730,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+48691,int_stack+65455,int_stack+31108,6);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+64951,int_stack+27600,int_stack+27320,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+65791,int_stack+64951,int_stack+32098,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+44491,int_stack+65791,int_stack+32728,10);
 /*--- compute (dp|gf) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+52291,int_stack+44491,int_stack+48691,150);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+64951,int_stack+19650,int_stack+19560,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+65221,int_stack+19776,int_stack+19650,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+65599,int_stack+65221,int_stack+64951,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+66139,int_stack+19944,int_stack+19776,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+54991,int_stack+66139,int_stack+65221,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+66139,int_stack+54991,int_stack+65599,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+54991,int_stack+20310,int_stack+20160,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+64951,int_stack+20520,int_stack+20310,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17100,int_stack+64951,int_stack+54991,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+54991,int_stack+20800,int_stack+20520,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+18000,int_stack+54991,int_stack+64951,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+19260,int_stack+18000,int_stack+17100,10);
 /*--- compute (dp|gf) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+30460,int_stack+19260,int_stack+66139, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48691, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17100,int_stack+21385,int_stack+21160,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+17775,int_stack+21700,int_stack+21385,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+64951,int_stack+17775,int_stack+17100,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+37696,int_stack+22120,int_stack+21700,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+20760,int_stack+37696,int_stack+17775,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+37696,int_stack+20760,int_stack+64951,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+71701,int_stack+37696,int_stack+19260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+37696,int_stack+22750,int_stack+22660,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+37966,int_stack+22876,int_stack+22750,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+38344,int_stack+37966,int_stack+37696,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+38884,int_stack+23044,int_stack+22876,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+64951,int_stack+38884,int_stack+37966,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+38884,int_stack+64951,int_stack+38344,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+64951,int_stack+23410,int_stack+23260,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+65401,int_stack+23620,int_stack+23410,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+66031,int_stack+65401,int_stack+64951,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+37696,int_stack+23900,int_stack+23620,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+17100,int_stack+37696,int_stack+65401,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+18360,int_stack+17100,int_stack+66031,10);
 /*--- compute (dp|gf) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+19860,int_stack+18360,int_stack+38884, 0.0, zero_stack, 1.0, int_stack+48691, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17100,int_stack+24485,int_stack+24260,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+37696,int_stack+24800,int_stack+24485,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+38641,int_stack+37696,int_stack+17100,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+17100,int_stack+25220,int_stack+24800,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+64951,int_stack+17100,int_stack+37696,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+22560,int_stack+64951,int_stack+38641,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+76201,int_stack+22560,int_stack+18360, 0.0, zero_stack, 1.0, int_stack+44491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22560,int_stack+26450,int_stack+26360,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+22830,int_stack+26576,int_stack+26450,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23208,int_stack+22830,int_stack+22560,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+23748,int_stack+26744,int_stack+26576,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+24252,int_stack+23748,int_stack+22830,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+25008,int_stack+24252,int_stack+23208,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22560,int_stack+28110,int_stack+27960,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+23010,int_stack+28320,int_stack+28110,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23640,int_stack+23010,int_stack+22560,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+25908,int_stack+28600,int_stack+28320,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+26748,int_stack+25908,int_stack+23010,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+64951,int_stack+26748,int_stack+23640,10);
 /*--- compute (dp|gf) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+25908,int_stack+64951,int_stack+25008, 1.0, int_stack+48691, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+48691,int_stack+29185,int_stack+28960,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+22560,int_stack+29500,int_stack+29185,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23505,int_stack+22560,int_stack+48691,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+37696,int_stack+29920,int_stack+29500,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+17100,int_stack+37696,int_stack+22560,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+37696,int_stack+17100,int_stack+23505,15);
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+80701,int_stack+37696,int_stack+64951, 1.0, int_stack+44491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (dd|gf) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+85201,int_stack+39991,int_stack+34996,150);
     Libderiv->ABCD[11] = int_stack + 85201;
 /*--- compute (dd|gf) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+33160,int_stack+0,int_stack+45991,150);
     Libderiv->ABCD[10] = int_stack + 33160;
 /*--- compute (dd|gf) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+38560,int_stack+55951,int_stack+4500,150);
     Libderiv->ABCD[9] = int_stack + 38560;
 /*--- compute (dd|gf) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+60451,int_stack+49591,150);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (dd|gf) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+54991,int_stack+9900,int_stack+7200,150);
     Libderiv->ABCD[7] = int_stack + 54991;
 /*--- compute (dd|gf) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+5400,int_stack+67201,int_stack+14400,150);
     Libderiv->ABCD[6] = int_stack + 5400;
 /*--- compute (dd|gf) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+10800,int_stack+71701,int_stack+30460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 10800;
 /*--- compute (dd|gf) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+60391,int_stack+76201,int_stack+19860, 0.0, zero_stack, 1.0, int_stack+52291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 60391;
 /*--- compute (dd|gf) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+16200,int_stack+80701,int_stack+25908, 1.0, int_stack+52291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 16200;

}
