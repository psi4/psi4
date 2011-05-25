#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|gg) integrals */

void d1hrr_order_fpgg(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[4][4][11] = int_stack + 1450;
 Libderiv->deriv_classes[4][5][11] = int_stack + 1675;
 Libderiv->deriv_classes[4][6][11] = int_stack + 1990;
 Libderiv->deriv_classes[4][7][11] = int_stack + 2410;
 Libderiv->deriv_classes[4][8][11] = int_stack + 2950;
 Libderiv->deriv_classes[3][4][10] = int_stack + 3625;
 Libderiv->deriv_classes[3][5][10] = int_stack + 3775;
 Libderiv->deriv_classes[3][6][10] = int_stack + 3985;
 Libderiv->deriv_classes[3][7][10] = int_stack + 4265;
 Libderiv->deriv_classes[3][8][10] = int_stack + 4625;
 Libderiv->deriv_classes[4][4][10] = int_stack + 5075;
 Libderiv->deriv_classes[4][5][10] = int_stack + 5300;
 Libderiv->deriv_classes[4][6][10] = int_stack + 5615;
 Libderiv->deriv_classes[4][7][10] = int_stack + 6035;
 Libderiv->deriv_classes[4][8][10] = int_stack + 6575;
 Libderiv->deriv_classes[3][4][9] = int_stack + 7250;
 Libderiv->deriv_classes[3][5][9] = int_stack + 7400;
 Libderiv->deriv_classes[3][6][9] = int_stack + 7610;
 Libderiv->deriv_classes[3][7][9] = int_stack + 7890;
 Libderiv->deriv_classes[3][8][9] = int_stack + 8250;
 Libderiv->deriv_classes[4][4][9] = int_stack + 8700;
 Libderiv->deriv_classes[4][5][9] = int_stack + 8925;
 Libderiv->deriv_classes[4][6][9] = int_stack + 9240;
 Libderiv->deriv_classes[4][7][9] = int_stack + 9660;
 Libderiv->deriv_classes[4][8][9] = int_stack + 10200;
 Libderiv->deriv_classes[3][4][8] = int_stack + 10875;
 Libderiv->deriv_classes[3][5][8] = int_stack + 11025;
 Libderiv->deriv_classes[3][6][8] = int_stack + 11235;
 Libderiv->deriv_classes[3][7][8] = int_stack + 11515;
 Libderiv->deriv_classes[3][8][8] = int_stack + 11875;
 Libderiv->deriv_classes[4][4][8] = int_stack + 12325;
 Libderiv->deriv_classes[4][5][8] = int_stack + 12550;
 Libderiv->deriv_classes[4][6][8] = int_stack + 12865;
 Libderiv->deriv_classes[4][7][8] = int_stack + 13285;
 Libderiv->deriv_classes[4][8][8] = int_stack + 13825;
 Libderiv->deriv_classes[3][4][7] = int_stack + 14500;
 Libderiv->deriv_classes[3][5][7] = int_stack + 14650;
 Libderiv->deriv_classes[3][6][7] = int_stack + 14860;
 Libderiv->deriv_classes[3][7][7] = int_stack + 15140;
 Libderiv->deriv_classes[3][8][7] = int_stack + 15500;
 Libderiv->deriv_classes[4][4][7] = int_stack + 15950;
 Libderiv->deriv_classes[4][5][7] = int_stack + 16175;
 Libderiv->deriv_classes[4][6][7] = int_stack + 16490;
 Libderiv->deriv_classes[4][7][7] = int_stack + 16910;
 Libderiv->deriv_classes[4][8][7] = int_stack + 17450;
 Libderiv->deriv_classes[3][4][6] = int_stack + 18125;
 Libderiv->deriv_classes[3][5][6] = int_stack + 18275;
 Libderiv->deriv_classes[3][6][6] = int_stack + 18485;
 Libderiv->deriv_classes[3][7][6] = int_stack + 18765;
 Libderiv->deriv_classes[3][8][6] = int_stack + 19125;
 Libderiv->dvrr_classes[4][4] = int_stack + 19575;
 Libderiv->deriv_classes[4][4][6] = int_stack + 19800;
 Libderiv->dvrr_classes[4][5] = int_stack + 20025;
 Libderiv->deriv_classes[4][5][6] = int_stack + 20340;
 Libderiv->dvrr_classes[4][6] = int_stack + 20655;
 Libderiv->deriv_classes[4][6][6] = int_stack + 21075;
 Libderiv->dvrr_classes[4][7] = int_stack + 21495;
 Libderiv->deriv_classes[4][7][6] = int_stack + 22035;
 Libderiv->deriv_classes[4][8][6] = int_stack + 22575;
 Libderiv->deriv_classes[3][4][2] = int_stack + 23250;
 Libderiv->deriv_classes[3][5][2] = int_stack + 23400;
 Libderiv->deriv_classes[3][6][2] = int_stack + 23610;
 Libderiv->deriv_classes[3][7][2] = int_stack + 23890;
 Libderiv->deriv_classes[3][8][2] = int_stack + 24250;
 Libderiv->deriv_classes[4][4][2] = int_stack + 24700;
 Libderiv->deriv_classes[4][5][2] = int_stack + 24925;
 Libderiv->deriv_classes[4][6][2] = int_stack + 25240;
 Libderiv->deriv_classes[4][7][2] = int_stack + 25660;
 Libderiv->deriv_classes[4][8][2] = int_stack + 26200;
 Libderiv->deriv_classes[3][4][1] = int_stack + 26875;
 Libderiv->deriv_classes[3][5][1] = int_stack + 27025;
 Libderiv->deriv_classes[3][6][1] = int_stack + 27235;
 Libderiv->deriv_classes[3][7][1] = int_stack + 27515;
 Libderiv->deriv_classes[3][8][1] = int_stack + 27875;
 Libderiv->deriv_classes[4][4][1] = int_stack + 28325;
 Libderiv->deriv_classes[4][5][1] = int_stack + 28550;
 Libderiv->deriv_classes[4][6][1] = int_stack + 28865;
 Libderiv->deriv_classes[4][7][1] = int_stack + 29285;
 Libderiv->deriv_classes[4][8][1] = int_stack + 29825;
 Libderiv->dvrr_classes[3][4] = int_stack + 30500;
 Libderiv->dvrr_classes[3][5] = int_stack + 30650;
 Libderiv->dvrr_classes[3][6] = int_stack + 30860;
 Libderiv->dvrr_classes[3][7] = int_stack + 31140;
 Libderiv->dvrr_classes[3][8] = int_stack + 31500;
 Libderiv->deriv_classes[3][4][0] = int_stack + 31950;
 Libderiv->deriv_classes[3][5][0] = int_stack + 32100;
 Libderiv->deriv_classes[3][6][0] = int_stack + 32310;
 Libderiv->deriv_classes[3][7][0] = int_stack + 32590;
 Libderiv->deriv_classes[3][8][0] = int_stack + 32950;
 Libderiv->deriv_classes[4][4][0] = int_stack + 33400;
 Libderiv->deriv_classes[4][5][0] = int_stack + 33625;
 Libderiv->deriv_classes[4][6][0] = int_stack + 33940;
 Libderiv->deriv_classes[4][7][0] = int_stack + 34360;
 Libderiv->deriv_classes[4][8][0] = int_stack + 34900;
 memset(int_stack,0,284600);

 Libderiv->dvrr_stack = int_stack + 86665;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35575,int_stack+30650,int_stack+30500,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+36025,int_stack+30860,int_stack+30650,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+36655,int_stack+36025,int_stack+35575,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+37555,int_stack+31140,int_stack+30860,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+38395,int_stack+37555,int_stack+36025,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+39655,int_stack+38395,int_stack+36655,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41155,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30500,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41605,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30650,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42235,int_stack+41605,int_stack+41155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35575,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+43135,int_stack+640,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30860,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+43975,int_stack+43135,int_stack+41605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36025,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+45235,int_stack+43975,int_stack+42235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36655,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+41155,int_stack+1000,int_stack+640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31140,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+46735,int_stack+41155,int_stack+43135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37555,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+41155,int_stack+46735,int_stack+43975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38395,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+46735,int_stack+41155,int_stack+45235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39655,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+41155,int_stack+20025,int_stack+19575,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+41830,int_stack+20655,int_stack+20025,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+42775,int_stack+41830,int_stack+41155,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+44125,int_stack+21495,int_stack+20655,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+48985,int_stack+44125,int_stack+41830,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+50875,int_stack+48985,int_stack+42775,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45385,int_stack+1675,int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19575,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+1990,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20025,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+53125,int_stack+0,int_stack+45385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41155,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+45385,int_stack+2410,int_stack+1990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20655,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+54475,int_stack+45385,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41830,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+54475,int_stack+53125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42775,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+56365,int_stack+2950,int_stack+2410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21495,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+57985,int_stack+56365,int_stack+45385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44125,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+60505,int_stack+57985,int_stack+54475, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48985,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+53125,int_stack+60505,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50875,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3775,int_stack+3625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30500, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+450,int_stack+3985,int_stack+3775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30650, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35575, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1980,int_stack+4265,int_stack+3985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30860, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2820,int_stack+1980,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36025, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+56500,int_stack+2820,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36655, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+4625,int_stack+4265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31140, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+58000,int_stack+0,int_stack+1980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37555, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+58000,int_stack+2820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38395, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+58000,int_stack+0,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39655, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+5300,int_stack+5075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19575, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+5615,int_stack+5300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20025, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45385,int_stack+0,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41155, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+56500,int_stack+6035,int_stack+5615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20655, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+945,int_stack+56500,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41830, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2835,int_stack+945,int_stack+45385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42775, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+60250,int_stack+6575,int_stack+6035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21495, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+61870,int_stack+60250,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44125, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+64390,int_stack+61870,int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48985, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+60250,int_stack+64390,int_stack+2835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50875, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+7400,int_stack+7250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+56950,int_stack+7610,int_stack+7400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30650, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+63625,int_stack+56950,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35575, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+64525,int_stack+7890,int_stack+7610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30860, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+65365,int_stack+64525,int_stack+56950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36025, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+56500,int_stack+65365,int_stack+63625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36655, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+66625,int_stack+8250,int_stack+7890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31140, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+0,int_stack+66625,int_stack+64525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37555, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+1680,int_stack+0,int_stack+65365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38395, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+63625,int_stack+1680,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39655, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+8925,int_stack+8700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19575, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+9240,int_stack+8925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20025, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45385,int_stack+0,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41155, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+56500,int_stack+9660,int_stack+9240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20655, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+945,int_stack+56500,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41830, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2835,int_stack+945,int_stack+45385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42775, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+5085,int_stack+10200,int_stack+9660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21495, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+6705,int_stack+5085,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44125, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+65875,int_stack+6705,int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48985, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+5085,int_stack+65875,int_stack+2835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50875, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65875,int_stack+11025,int_stack+10875, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66325,int_stack+11235,int_stack+11025, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+66955,int_stack+66325,int_stack+65875, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+67855,int_stack+11515,int_stack+11235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+68695,int_stack+67855,int_stack+66325, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+56500,int_stack+68695,int_stack+66955, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+65875,int_stack+11875,int_stack+11515, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+8460,int_stack+65875,int_stack+67855, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+65875,int_stack+8460,int_stack+68695, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+8460,int_stack+65875,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+12550,int_stack+12325, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+65875,int_stack+12865,int_stack+12550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45385,int_stack+65875,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+56500,int_stack+13285,int_stack+12865, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+66820,int_stack+56500,int_stack+65875, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+10710,int_stack+66820,int_stack+45385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+13825,int_stack+13285, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1620,int_stack+0,int_stack+56500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+68710,int_stack+1620,int_stack+66820, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+0,int_stack+68710,int_stack+10710, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10710,int_stack+14650,int_stack+14500, 0.0, zero_stack, 1.0, int_stack+30500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11160,int_stack+14860,int_stack+14650, 0.0, zero_stack, 1.0, int_stack+30650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11790,int_stack+11160,int_stack+10710, 0.0, zero_stack, 1.0, int_stack+35575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+12690,int_stack+15140,int_stack+14860, 0.0, zero_stack, 1.0, int_stack+30860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13530,int_stack+12690,int_stack+11160, 0.0, zero_stack, 1.0, int_stack+36025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+56500,int_stack+13530,int_stack+11790, 0.0, zero_stack, 1.0, int_stack+36655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+10710,int_stack+15500,int_stack+15140, 0.0, zero_stack, 1.0, int_stack+31140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3375,int_stack+10710,int_stack+12690, 0.0, zero_stack, 1.0, int_stack+37555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+10710,int_stack+3375,int_stack+13530, 0.0, zero_stack, 1.0, int_stack+38395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+12810,int_stack+10710,int_stack+56500, 0.0, zero_stack, 1.0, int_stack+39655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+16175,int_stack+15950, 0.0, zero_stack, 1.0, int_stack+19575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10710,int_stack+16490,int_stack+16175, 0.0, zero_stack, 1.0, int_stack+20025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45385,int_stack+10710,int_stack+56500, 0.0, zero_stack, 1.0, int_stack+41155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+56500,int_stack+16910,int_stack+16490, 0.0, zero_stack, 1.0, int_stack+20655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+65875,int_stack+56500,int_stack+10710, 0.0, zero_stack, 1.0, int_stack+41830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+67765,int_stack+65875,int_stack+45385, 0.0, zero_stack, 1.0, int_stack+42775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+10710,int_stack+17450,int_stack+16910, 0.0, zero_stack, 1.0, int_stack+21495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+15060,int_stack+10710,int_stack+56500, 0.0, zero_stack, 1.0, int_stack+44125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+70015,int_stack+15060,int_stack+65875, 0.0, zero_stack, 1.0, int_stack+48985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+73165,int_stack+70015,int_stack+67765, 0.0, zero_stack, 1.0, int_stack+50875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65875,int_stack+18275,int_stack+18125, 1.0, int_stack+30500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66325,int_stack+18485,int_stack+18275, 1.0, int_stack+30650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+66955,int_stack+66325,int_stack+65875, 1.0, int_stack+35575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+67855,int_stack+18765,int_stack+18485, 1.0, int_stack+30860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+68695,int_stack+67855,int_stack+66325, 1.0, int_stack+36025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+56500,int_stack+68695,int_stack+66955, 1.0, int_stack+36655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+35575,int_stack+19125,int_stack+18765, 1.0, int_stack+31140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+65875,int_stack+35575,int_stack+67855, 1.0, int_stack+37555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+10710,int_stack+65875,int_stack+68695, 1.0, int_stack+38395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+65875,int_stack+10710,int_stack+56500, 1.0, int_stack+39655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+20340,int_stack+19800, 1.0, int_stack+19575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10710,int_stack+21075,int_stack+20340, 1.0, int_stack+20025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45385,int_stack+10710,int_stack+56500, 1.0, int_stack+41155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+56500,int_stack+22035,int_stack+21075, 1.0, int_stack+20655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+68125,int_stack+56500,int_stack+10710, 1.0, int_stack+41830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+70015,int_stack+68125,int_stack+45385, 1.0, int_stack+42775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+10710,int_stack+22575,int_stack+22035, 1.0, int_stack+21495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+41155,int_stack+10710,int_stack+56500, 1.0, int_stack+44125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+15060,int_stack+41155,int_stack+68125, 1.0, int_stack+48985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+41155,int_stack+15060,int_stack+70015, 1.0, int_stack+50875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+15060,int_stack+31500,int_stack+31140,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+16140,int_stack+15060,int_stack+37555,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+10710,int_stack+16140,int_stack+38395,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+15060,int_stack+10710,int_stack+39655,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10710,int_stack+23400,int_stack+23250,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11160,int_stack+23610,int_stack+23400,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11790,int_stack+11160,int_stack+10710,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+17310,int_stack+23890,int_stack+23610,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+18150,int_stack+17310,int_stack+11160,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+56500,int_stack+18150,int_stack+11790,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+10710,int_stack+24250,int_stack+23890,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+19410,int_stack+10710,int_stack+17310,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+10710,int_stack+19410,int_stack+18150,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+17310,int_stack+10710,int_stack+56500,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+24925,int_stack+24700,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10710,int_stack+25240,int_stack+24925,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+19560,int_stack+10710,int_stack+56500,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+56500,int_stack+25660,int_stack+25240,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+20910,int_stack+56500,int_stack+10710,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+22800,int_stack+20910,int_stack+19560,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+10710,int_stack+26200,int_stack+25660,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+48985,int_stack+10710,int_stack+56500,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+68125,int_stack+48985,int_stack+20910,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+48985,int_stack+68125,int_stack+22800,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+68125,int_stack+27025,int_stack+26875,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+68575,int_stack+27235,int_stack+27025,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+69205,int_stack+68575,int_stack+68125,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+70105,int_stack+27515,int_stack+27235,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+70945,int_stack+70105,int_stack+68575,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+56500,int_stack+70945,int_stack+69205,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+68125,int_stack+27875,int_stack+27515,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+10710,int_stack+68125,int_stack+70105,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+68125,int_stack+10710,int_stack+70945,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+70225,int_stack+68125,int_stack+56500,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+28550,int_stack+28325,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+68125,int_stack+28865,int_stack+28550,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10710,int_stack+68125,int_stack+56500,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+56500,int_stack+29285,int_stack+28865,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+19560,int_stack+56500,int_stack+68125,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+21450,int_stack+19560,int_stack+10710,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+10710,int_stack+29825,int_stack+29285,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+23700,int_stack+10710,int_stack+56500,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+26220,int_stack+23700,int_stack+19560,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+35575,int_stack+26220,int_stack+21450,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19560,int_stack+32100,int_stack+31950,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+20010,int_stack+32310,int_stack+32100,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+20640,int_stack+20010,int_stack+19560,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+21540,int_stack+32590,int_stack+32310,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+22380,int_stack+21540,int_stack+20010,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+56500,int_stack+22380,int_stack+20640,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+19560,int_stack+32950,int_stack+32590,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+23640,int_stack+19560,int_stack+21540,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+10710,int_stack+23640,int_stack+22380,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+19560,int_stack+10710,int_stack+56500,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+56500,int_stack+33625,int_stack+33400,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10710,int_stack+33940,int_stack+33625,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+21810,int_stack+10710,int_stack+56500,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+56500,int_stack+34360,int_stack+33940,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+23160,int_stack+56500,int_stack+10710,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+25050,int_stack+23160,int_stack+21810,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+10710,int_stack+34900,int_stack+34360,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+27300,int_stack+10710,int_stack+56500,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+29820,int_stack+27300,int_stack+23160,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+76540,int_stack+29820,int_stack+25050,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+21810,int_stack+53125,int_stack+46735,225);
     Libderiv->ABCD[11] = int_stack + 21810;
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28560,int_stack+60250,int_stack+58000,225);
     Libderiv->ABCD[10] = int_stack + 28560;
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+52360,int_stack+5085,int_stack+63625,225);
     Libderiv->ABCD[9] = int_stack + 52360;
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+59110,int_stack+0,int_stack+8460,225);
     Libderiv->ABCD[8] = int_stack + 59110;
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+73165,int_stack+12810,225);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+6750,int_stack+41155,int_stack+65875,225);
     Libderiv->ABCD[6] = int_stack + 6750;
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+38950,int_stack+48985,int_stack+17310, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 38950;
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+79915,int_stack+35575,int_stack+70225, 0.0, zero_stack, 1.0, int_stack+15060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 79915;
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+65860,int_stack+76540,int_stack+19560, 1.0, int_stack+15060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 65860;

}
