#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|gf) integrals */

void d1hrr_order_fpgf(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[3][4][10] = int_stack + 2500;
 Libderiv->deriv_classes[3][5][10] = int_stack + 2650;
 Libderiv->deriv_classes[3][6][10] = int_stack + 2860;
 Libderiv->deriv_classes[3][7][10] = int_stack + 3140;
 Libderiv->deriv_classes[4][4][10] = int_stack + 3500;
 Libderiv->deriv_classes[4][5][10] = int_stack + 3725;
 Libderiv->deriv_classes[4][6][10] = int_stack + 4040;
 Libderiv->deriv_classes[4][7][10] = int_stack + 4460;
 Libderiv->deriv_classes[3][4][9] = int_stack + 5000;
 Libderiv->deriv_classes[3][5][9] = int_stack + 5150;
 Libderiv->deriv_classes[3][6][9] = int_stack + 5360;
 Libderiv->deriv_classes[3][7][9] = int_stack + 5640;
 Libderiv->deriv_classes[4][4][9] = int_stack + 6000;
 Libderiv->deriv_classes[4][5][9] = int_stack + 6225;
 Libderiv->deriv_classes[4][6][9] = int_stack + 6540;
 Libderiv->deriv_classes[4][7][9] = int_stack + 6960;
 Libderiv->deriv_classes[3][4][8] = int_stack + 7500;
 Libderiv->deriv_classes[3][5][8] = int_stack + 7650;
 Libderiv->deriv_classes[3][6][8] = int_stack + 7860;
 Libderiv->deriv_classes[3][7][8] = int_stack + 8140;
 Libderiv->deriv_classes[4][4][8] = int_stack + 8500;
 Libderiv->deriv_classes[4][5][8] = int_stack + 8725;
 Libderiv->deriv_classes[4][6][8] = int_stack + 9040;
 Libderiv->deriv_classes[4][7][8] = int_stack + 9460;
 Libderiv->deriv_classes[3][4][7] = int_stack + 10000;
 Libderiv->deriv_classes[3][5][7] = int_stack + 10150;
 Libderiv->deriv_classes[3][6][7] = int_stack + 10360;
 Libderiv->deriv_classes[3][7][7] = int_stack + 10640;
 Libderiv->deriv_classes[4][4][7] = int_stack + 11000;
 Libderiv->deriv_classes[4][5][7] = int_stack + 11225;
 Libderiv->deriv_classes[4][6][7] = int_stack + 11540;
 Libderiv->deriv_classes[4][7][7] = int_stack + 11960;
 Libderiv->deriv_classes[3][4][6] = int_stack + 12500;
 Libderiv->deriv_classes[3][5][6] = int_stack + 12650;
 Libderiv->deriv_classes[3][6][6] = int_stack + 12860;
 Libderiv->deriv_classes[3][7][6] = int_stack + 13140;
 Libderiv->dvrr_classes[4][4] = int_stack + 13500;
 Libderiv->deriv_classes[4][4][6] = int_stack + 13725;
 Libderiv->dvrr_classes[4][5] = int_stack + 13950;
 Libderiv->deriv_classes[4][5][6] = int_stack + 14265;
 Libderiv->dvrr_classes[4][6] = int_stack + 14580;
 Libderiv->deriv_classes[4][6][6] = int_stack + 15000;
 Libderiv->deriv_classes[4][7][6] = int_stack + 15420;
 Libderiv->deriv_classes[3][4][2] = int_stack + 15960;
 Libderiv->deriv_classes[3][5][2] = int_stack + 16110;
 Libderiv->deriv_classes[3][6][2] = int_stack + 16320;
 Libderiv->deriv_classes[3][7][2] = int_stack + 16600;
 Libderiv->deriv_classes[4][4][2] = int_stack + 16960;
 Libderiv->deriv_classes[4][5][2] = int_stack + 17185;
 Libderiv->deriv_classes[4][6][2] = int_stack + 17500;
 Libderiv->deriv_classes[4][7][2] = int_stack + 17920;
 Libderiv->deriv_classes[3][4][1] = int_stack + 18460;
 Libderiv->deriv_classes[3][5][1] = int_stack + 18610;
 Libderiv->deriv_classes[3][6][1] = int_stack + 18820;
 Libderiv->deriv_classes[3][7][1] = int_stack + 19100;
 Libderiv->deriv_classes[4][4][1] = int_stack + 19460;
 Libderiv->deriv_classes[4][5][1] = int_stack + 19685;
 Libderiv->deriv_classes[4][6][1] = int_stack + 20000;
 Libderiv->deriv_classes[4][7][1] = int_stack + 20420;
 Libderiv->dvrr_classes[3][4] = int_stack + 20960;
 Libderiv->dvrr_classes[3][5] = int_stack + 21110;
 Libderiv->dvrr_classes[3][6] = int_stack + 21320;
 Libderiv->dvrr_classes[3][7] = int_stack + 21600;
 Libderiv->deriv_classes[3][4][0] = int_stack + 21960;
 Libderiv->deriv_classes[3][5][0] = int_stack + 22110;
 Libderiv->deriv_classes[3][6][0] = int_stack + 22320;
 Libderiv->deriv_classes[3][7][0] = int_stack + 22600;
 Libderiv->deriv_classes[4][4][0] = int_stack + 22960;
 Libderiv->deriv_classes[4][5][0] = int_stack + 23185;
 Libderiv->deriv_classes[4][6][0] = int_stack + 23500;
 Libderiv->deriv_classes[4][7][0] = int_stack + 23920;
 memset(int_stack,0,195680);

 Libderiv->dvrr_stack = int_stack + 54250;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24460,int_stack+21110,int_stack+20960,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24910,int_stack+21320,int_stack+21110,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+25540,int_stack+24910,int_stack+24460,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26440,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20960,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26890,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21110,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27520,int_stack+26890,int_stack+26440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24460,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+28420,int_stack+640,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21320,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+29260,int_stack+28420,int_stack+26890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24910,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+30520,int_stack+29260,int_stack+27520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25540,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+26440,int_stack+13950,int_stack+13500,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+27115,int_stack+14580,int_stack+13950,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+28060,int_stack+27115,int_stack+26440,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29410,int_stack+1225,int_stack+1000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13500,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+1540,int_stack+1225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13950,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+32020,int_stack+0,int_stack+29410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26440,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+33370,int_stack+1960,int_stack+1540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14580,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+34630,int_stack+33370,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27115,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+34630,int_stack+32020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28060,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32020,int_stack+2650,int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20960, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+32470,int_stack+2860,int_stack+2650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21110, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33100,int_stack+32470,int_stack+32020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24460, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+34000,int_stack+3140,int_stack+2860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21320, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+34840,int_stack+34000,int_stack+32470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24910, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+36100,int_stack+34840,int_stack+33100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25540, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32020,int_stack+3725,int_stack+3500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13500, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+32695,int_stack+4040,int_stack+3725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13950, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33640,int_stack+32695,int_stack+32020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26440, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2250,int_stack+4460,int_stack+4040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14580, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+37600,int_stack+2250,int_stack+32695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27115, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2250,int_stack+37600,int_stack+33640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28060, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37600,int_stack+5150,int_stack+5000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20960, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38050,int_stack+5360,int_stack+5150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21110, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+38680,int_stack+38050,int_stack+37600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24460, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4500,int_stack+5640,int_stack+5360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21320, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+32020,int_stack+4500,int_stack+38050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24910, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4500,int_stack+32020,int_stack+38680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25540, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32020,int_stack+6225,int_stack+6000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+32695,int_stack+6540,int_stack+6225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13950, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33640,int_stack+32695,int_stack+32020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26440, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+37600,int_stack+6960,int_stack+6540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14580, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+38860,int_stack+37600,int_stack+32695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27115, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+40750,int_stack+38860,int_stack+33640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28060, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37600,int_stack+7650,int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38050,int_stack+7860,int_stack+7650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+38680,int_stack+38050,int_stack+37600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+39580,int_stack+8140,int_stack+7860, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+32020,int_stack+39580,int_stack+38050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+33280,int_stack+32020,int_stack+38680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32020,int_stack+8725,int_stack+8500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+37600,int_stack+9040,int_stack+8725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+38545,int_stack+37600,int_stack+32020, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+32020,int_stack+9460,int_stack+9040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+6000,int_stack+32020,int_stack+37600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+43000,int_stack+6000,int_stack+38545, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6000,int_stack+10150,int_stack+10000, 0.0, zero_stack, 1.0, int_stack+20960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6450,int_stack+10360,int_stack+10150, 0.0, zero_stack, 1.0, int_stack+21110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7080,int_stack+6450,int_stack+6000, 0.0, zero_stack, 1.0, int_stack+24460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+7980,int_stack+10640,int_stack+10360, 0.0, zero_stack, 1.0, int_stack+21320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+32020,int_stack+7980,int_stack+6450, 0.0, zero_stack, 1.0, int_stack+24910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+7980,int_stack+32020,int_stack+7080, 0.0, zero_stack, 1.0, int_stack+25540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32020,int_stack+11225,int_stack+11000, 0.0, zero_stack, 1.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9480,int_stack+11540,int_stack+11225, 0.0, zero_stack, 1.0, int_stack+13950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6000,int_stack+9480,int_stack+32020, 0.0, zero_stack, 1.0, int_stack+26440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+32020,int_stack+11960,int_stack+11540, 0.0, zero_stack, 1.0, int_stack+14580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10425,int_stack+32020,int_stack+9480, 0.0, zero_stack, 1.0, int_stack+27115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+37600,int_stack+10425,int_stack+6000, 0.0, zero_stack, 1.0, int_stack+28060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6000,int_stack+12650,int_stack+12500, 1.0, int_stack+20960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6450,int_stack+12860,int_stack+12650, 1.0, int_stack+21110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7080,int_stack+6450,int_stack+6000, 1.0, int_stack+24460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9480,int_stack+13140,int_stack+12860, 1.0, int_stack+21320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+32020,int_stack+9480,int_stack+6450, 1.0, int_stack+24910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+9480,int_stack+32020,int_stack+7080, 1.0, int_stack+25540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32020,int_stack+14265,int_stack+13725, 1.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10980,int_stack+15000,int_stack+14265, 1.0, int_stack+13950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11925,int_stack+10980,int_stack+32020, 1.0, int_stack+26440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+32020,int_stack+15420,int_stack+15000, 1.0, int_stack+14580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13275,int_stack+32020,int_stack+10980, 1.0, int_stack+27115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+45250,int_stack+13275,int_stack+11925, 1.0, int_stack+28060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+10980,int_stack+21600,int_stack+21320,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+32020,int_stack+10980,int_stack+24910,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+10980,int_stack+32020,int_stack+25540,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+32020,int_stack+16110,int_stack+15960,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+32470,int_stack+16320,int_stack+16110,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+39850,int_stack+32470,int_stack+32020,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+12480,int_stack+16600,int_stack+16320,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+13320,int_stack+12480,int_stack+32470,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+14580,int_stack+13320,int_stack+39850,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+39850,int_stack+17185,int_stack+16960,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12480,int_stack+17500,int_stack+17185,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16080,int_stack+12480,int_stack+39850,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+32020,int_stack+17920,int_stack+17500,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+24460,int_stack+32020,int_stack+12480,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+26350,int_stack+24460,int_stack+16080,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16080,int_stack+18610,int_stack+18460,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16530,int_stack+18820,int_stack+18610,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+39850,int_stack+16530,int_stack+16080,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+17160,int_stack+19100,int_stack+18820,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+32020,int_stack+17160,int_stack+16530,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+16080,int_stack+32020,int_stack+39850,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+39850,int_stack+19685,int_stack+19460,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+32020,int_stack+20000,int_stack+19685,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17580,int_stack+32020,int_stack+39850,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+24460,int_stack+20420,int_stack+20000,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+18930,int_stack+24460,int_stack+32020,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+47500,int_stack+18930,int_stack+17580,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17580,int_stack+22110,int_stack+21960,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+18030,int_stack+22320,int_stack+22110,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+39850,int_stack+18030,int_stack+17580,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+18660,int_stack+22600,int_stack+22320,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+32020,int_stack+18660,int_stack+18030,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+17580,int_stack+32020,int_stack+39850,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+39850,int_stack+23185,int_stack+22960,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+32020,int_stack+23500,int_stack+23185,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+19080,int_stack+32020,int_stack+39850,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+20430,int_stack+23920,int_stack+23500,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+21690,int_stack+20430,int_stack+32020,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+23580,int_stack+21690,int_stack+19080,15);
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+19080,int_stack+0,int_stack+30520,150);
     Libderiv->ABCD[11] = int_stack + 19080;
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28600,int_stack+2250,int_stack+36100,150);
     Libderiv->ABCD[10] = int_stack + 28600;
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+40750,int_stack+4500,150);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+49750,int_stack+43000,int_stack+33280,150);
     Libderiv->ABCD[8] = int_stack + 49750;
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+33100,int_stack+37600,int_stack+7980,150);
     Libderiv->ABCD[7] = int_stack + 33100;
 /*--- compute (fp|gf) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+37600,int_stack+45250,int_stack+9480,150);
     Libderiv->ABCD[6] = int_stack + 37600;
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+42100,int_stack+26350,int_stack+14580, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 42100;
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+4500,int_stack+47500,int_stack+16080, 0.0, zero_stack, 1.0, int_stack+10980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 4500;
 /*--- compute (fp|gf) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+12480,int_stack+23580,int_stack+17580, 1.0, int_stack+10980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 12480;

}
