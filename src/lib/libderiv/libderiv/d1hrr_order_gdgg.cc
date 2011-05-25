#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gdgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gd|gg) integrals */

void d1hrr_order_gdgg(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[6][4][11] = int_stack + 5220;
 Libderiv->deriv_classes[6][5][11] = int_stack + 5640;
 Libderiv->deriv_classes[6][6][11] = int_stack + 6228;
 Libderiv->deriv_classes[6][7][11] = int_stack + 7012;
 Libderiv->deriv_classes[6][8][11] = int_stack + 8020;
 Libderiv->deriv_classes[4][4][10] = int_stack + 9280;
 Libderiv->deriv_classes[4][5][10] = int_stack + 9505;
 Libderiv->deriv_classes[4][6][10] = int_stack + 9820;
 Libderiv->deriv_classes[4][7][10] = int_stack + 10240;
 Libderiv->deriv_classes[4][8][10] = int_stack + 10780;
 Libderiv->deriv_classes[5][4][10] = int_stack + 11455;
 Libderiv->deriv_classes[5][5][10] = int_stack + 11770;
 Libderiv->deriv_classes[5][6][10] = int_stack + 12211;
 Libderiv->deriv_classes[5][7][10] = int_stack + 12799;
 Libderiv->deriv_classes[5][8][10] = int_stack + 13555;
 Libderiv->deriv_classes[6][4][10] = int_stack + 14500;
 Libderiv->deriv_classes[6][5][10] = int_stack + 14920;
 Libderiv->deriv_classes[6][6][10] = int_stack + 15508;
 Libderiv->deriv_classes[6][7][10] = int_stack + 16292;
 Libderiv->deriv_classes[6][8][10] = int_stack + 17300;
 Libderiv->deriv_classes[4][4][9] = int_stack + 18560;
 Libderiv->deriv_classes[4][5][9] = int_stack + 18785;
 Libderiv->deriv_classes[4][6][9] = int_stack + 19100;
 Libderiv->deriv_classes[4][7][9] = int_stack + 19520;
 Libderiv->deriv_classes[4][8][9] = int_stack + 20060;
 Libderiv->deriv_classes[5][4][9] = int_stack + 20735;
 Libderiv->deriv_classes[5][5][9] = int_stack + 21050;
 Libderiv->deriv_classes[5][6][9] = int_stack + 21491;
 Libderiv->deriv_classes[5][7][9] = int_stack + 22079;
 Libderiv->deriv_classes[5][8][9] = int_stack + 22835;
 Libderiv->deriv_classes[6][4][9] = int_stack + 23780;
 Libderiv->deriv_classes[6][5][9] = int_stack + 24200;
 Libderiv->deriv_classes[6][6][9] = int_stack + 24788;
 Libderiv->deriv_classes[6][7][9] = int_stack + 25572;
 Libderiv->deriv_classes[6][8][9] = int_stack + 26580;
 Libderiv->deriv_classes[4][4][8] = int_stack + 27840;
 Libderiv->deriv_classes[4][5][8] = int_stack + 28065;
 Libderiv->deriv_classes[4][6][8] = int_stack + 28380;
 Libderiv->deriv_classes[4][7][8] = int_stack + 28800;
 Libderiv->deriv_classes[4][8][8] = int_stack + 29340;
 Libderiv->deriv_classes[5][4][8] = int_stack + 30015;
 Libderiv->deriv_classes[5][5][8] = int_stack + 30330;
 Libderiv->deriv_classes[5][6][8] = int_stack + 30771;
 Libderiv->deriv_classes[5][7][8] = int_stack + 31359;
 Libderiv->deriv_classes[5][8][8] = int_stack + 32115;
 Libderiv->deriv_classes[6][4][8] = int_stack + 33060;
 Libderiv->deriv_classes[6][5][8] = int_stack + 33480;
 Libderiv->deriv_classes[6][6][8] = int_stack + 34068;
 Libderiv->deriv_classes[6][7][8] = int_stack + 34852;
 Libderiv->deriv_classes[6][8][8] = int_stack + 35860;
 Libderiv->deriv_classes[4][4][7] = int_stack + 37120;
 Libderiv->deriv_classes[4][5][7] = int_stack + 37345;
 Libderiv->deriv_classes[4][6][7] = int_stack + 37660;
 Libderiv->deriv_classes[4][7][7] = int_stack + 38080;
 Libderiv->deriv_classes[4][8][7] = int_stack + 38620;
 Libderiv->deriv_classes[5][4][7] = int_stack + 39295;
 Libderiv->deriv_classes[5][5][7] = int_stack + 39610;
 Libderiv->deriv_classes[5][6][7] = int_stack + 40051;
 Libderiv->deriv_classes[5][7][7] = int_stack + 40639;
 Libderiv->deriv_classes[5][8][7] = int_stack + 41395;
 Libderiv->deriv_classes[6][4][7] = int_stack + 42340;
 Libderiv->deriv_classes[6][5][7] = int_stack + 42760;
 Libderiv->deriv_classes[6][6][7] = int_stack + 43348;
 Libderiv->deriv_classes[6][7][7] = int_stack + 44132;
 Libderiv->deriv_classes[6][8][7] = int_stack + 45140;
 Libderiv->deriv_classes[4][4][6] = int_stack + 46400;
 Libderiv->deriv_classes[4][5][6] = int_stack + 46625;
 Libderiv->deriv_classes[4][6][6] = int_stack + 46940;
 Libderiv->deriv_classes[4][7][6] = int_stack + 47360;
 Libderiv->deriv_classes[4][8][6] = int_stack + 47900;
 Libderiv->deriv_classes[5][4][6] = int_stack + 48575;
 Libderiv->deriv_classes[5][5][6] = int_stack + 48890;
 Libderiv->deriv_classes[5][6][6] = int_stack + 49331;
 Libderiv->deriv_classes[5][7][6] = int_stack + 49919;
 Libderiv->deriv_classes[5][8][6] = int_stack + 50675;
 Libderiv->dvrr_classes[6][4] = int_stack + 51620;
 Libderiv->deriv_classes[6][4][6] = int_stack + 52040;
 Libderiv->dvrr_classes[6][5] = int_stack + 52460;
 Libderiv->deriv_classes[6][5][6] = int_stack + 53048;
 Libderiv->dvrr_classes[6][6] = int_stack + 53636;
 Libderiv->deriv_classes[6][6][6] = int_stack + 54420;
 Libderiv->dvrr_classes[6][7] = int_stack + 55204;
 Libderiv->deriv_classes[6][7][6] = int_stack + 56212;
 Libderiv->deriv_classes[6][8][6] = int_stack + 57220;
 Libderiv->deriv_classes[4][4][2] = int_stack + 58480;
 Libderiv->deriv_classes[4][5][2] = int_stack + 58705;
 Libderiv->deriv_classes[4][6][2] = int_stack + 59020;
 Libderiv->deriv_classes[4][7][2] = int_stack + 59440;
 Libderiv->deriv_classes[4][8][2] = int_stack + 59980;
 Libderiv->deriv_classes[5][4][2] = int_stack + 60655;
 Libderiv->deriv_classes[5][5][2] = int_stack + 60970;
 Libderiv->deriv_classes[5][6][2] = int_stack + 61411;
 Libderiv->deriv_classes[5][7][2] = int_stack + 61999;
 Libderiv->deriv_classes[5][8][2] = int_stack + 62755;
 Libderiv->deriv_classes[6][4][2] = int_stack + 63700;
 Libderiv->deriv_classes[6][5][2] = int_stack + 64120;
 Libderiv->deriv_classes[6][6][2] = int_stack + 64708;
 Libderiv->deriv_classes[6][7][2] = int_stack + 65492;
 Libderiv->deriv_classes[6][8][2] = int_stack + 66500;
 Libderiv->deriv_classes[4][4][1] = int_stack + 67760;
 Libderiv->deriv_classes[4][5][1] = int_stack + 67985;
 Libderiv->deriv_classes[4][6][1] = int_stack + 68300;
 Libderiv->deriv_classes[4][7][1] = int_stack + 68720;
 Libderiv->deriv_classes[4][8][1] = int_stack + 69260;
 Libderiv->deriv_classes[5][4][1] = int_stack + 69935;
 Libderiv->deriv_classes[5][5][1] = int_stack + 70250;
 Libderiv->deriv_classes[5][6][1] = int_stack + 70691;
 Libderiv->deriv_classes[5][7][1] = int_stack + 71279;
 Libderiv->deriv_classes[5][8][1] = int_stack + 72035;
 Libderiv->deriv_classes[6][4][1] = int_stack + 72980;
 Libderiv->deriv_classes[6][5][1] = int_stack + 73400;
 Libderiv->deriv_classes[6][6][1] = int_stack + 73988;
 Libderiv->deriv_classes[6][7][1] = int_stack + 74772;
 Libderiv->deriv_classes[6][8][1] = int_stack + 75780;
 Libderiv->dvrr_classes[4][4] = int_stack + 77040;
 Libderiv->dvrr_classes[4][5] = int_stack + 77265;
 Libderiv->dvrr_classes[4][6] = int_stack + 77580;
 Libderiv->dvrr_classes[4][7] = int_stack + 78000;
 Libderiv->dvrr_classes[4][8] = int_stack + 78540;
 Libderiv->deriv_classes[4][4][0] = int_stack + 79215;
 Libderiv->deriv_classes[4][5][0] = int_stack + 79440;
 Libderiv->deriv_classes[4][6][0] = int_stack + 79755;
 Libderiv->deriv_classes[4][7][0] = int_stack + 80175;
 Libderiv->deriv_classes[4][8][0] = int_stack + 80715;
 Libderiv->dvrr_classes[5][4] = int_stack + 81390;
 Libderiv->dvrr_classes[5][5] = int_stack + 81705;
 Libderiv->dvrr_classes[5][6] = int_stack + 82146;
 Libderiv->dvrr_classes[5][7] = int_stack + 82734;
 Libderiv->dvrr_classes[5][8] = int_stack + 83490;
 Libderiv->deriv_classes[5][4][0] = int_stack + 84435;
 Libderiv->deriv_classes[5][5][0] = int_stack + 84750;
 Libderiv->deriv_classes[5][6][0] = int_stack + 85191;
 Libderiv->deriv_classes[5][7][0] = int_stack + 85779;
 Libderiv->deriv_classes[5][8][0] = int_stack + 86535;
 Libderiv->deriv_classes[6][4][0] = int_stack + 87480;
 Libderiv->deriv_classes[6][5][0] = int_stack + 87900;
 Libderiv->deriv_classes[6][6][0] = int_stack + 88488;
 Libderiv->deriv_classes[6][7][0] = int_stack + 89272;
 Libderiv->deriv_classes[6][8][0] = int_stack + 90280;
 memset(int_stack,0,732320);

 Libderiv->dvrr_stack = int_stack + 308464;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gdgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+91540,int_stack+77265,int_stack+77040,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+92215,int_stack+77580,int_stack+77265,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+93160,int_stack+92215,int_stack+91540,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+94510,int_stack+78000,int_stack+77580,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+95770,int_stack+94510,int_stack+92215,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+97660,int_stack+95770,int_stack+93160,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+99910,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77040,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+100585,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77265,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+101530,int_stack+100585,int_stack+99910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+102880,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77580,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+104140,int_stack+102880,int_stack+100585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92215,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+106030,int_stack+104140,int_stack+101530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93160,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+99910,int_stack+1500,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78000,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+108280,int_stack+99910,int_stack+102880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+94510,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+99910,int_stack+108280,int_stack+104140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95770,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+108280,int_stack+99910,int_stack+106030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97660,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+99910,int_stack+81705,int_stack+81390,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+100855,int_stack+82146,int_stack+81705,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+102178,int_stack+100855,int_stack+99910,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+104068,int_stack+82734,int_stack+82146,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+111655,int_stack+104068,int_stack+100855,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+114301,int_stack+111655,int_stack+102178,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+105832,int_stack+2490,int_stack+2175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81390,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+106777,int_stack+2931,int_stack+2490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81705,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+106777,int_stack+105832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99910,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+117451,int_stack+3519,int_stack+2931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82146,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+119215,int_stack+117451,int_stack+106777, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100855,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+121861,int_stack+119215,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102178,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+4275,int_stack+3519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82734,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+125011,int_stack+0,int_stack+117451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104068,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+125011,int_stack+119215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111655,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+125011,int_stack+0,int_stack+121861, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114301,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+129736,int_stack+125011,int_stack+108280,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+52460,int_stack+51620,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1260,int_stack+53636,int_stack+52460,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+117451,int_stack+1260,int_stack+0,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+119971,int_stack+55204,int_stack+53636,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+105832,int_stack+119971,int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+139861,int_stack+105832,int_stack+117451,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+5640,int_stack+5220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51620,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+122323,int_stack+6228,int_stack+5640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52460,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+144061,int_stack+122323,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3024,int_stack+7012,int_stack+6228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53636,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+146581,int_stack+3024,int_stack+122323, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+150109,int_stack+146581,int_stack+144061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+117451,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+154309,int_stack+8020,int_stack+7012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55204,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+157333,int_stack+154309,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119971,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3024,int_stack+157333,int_stack+146581, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105832,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+154309,int_stack+3024,int_stack+150109, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139861,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+160609,int_stack+154309,int_stack+125011,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+9505,int_stack+9280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77040, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3699,int_stack+9820,int_stack+9505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77265, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4644,int_stack+3699,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+5994,int_stack+10240,int_stack+9820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77580, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7254,int_stack+5994,int_stack+3699, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92215, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+144061,int_stack+7254,int_stack+4644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93160, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+3024,int_stack+10780,int_stack+10240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78000, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+146311,int_stack+3024,int_stack+5994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+94510, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3024,int_stack+146311,int_stack+7254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95770, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+146311,int_stack+3024,int_stack+144061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97660, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+144061,int_stack+11770,int_stack+11455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81390, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3024,int_stack+12211,int_stack+11770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81705, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4347,int_stack+3024,int_stack+144061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99910, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+144061,int_stack+12799,int_stack+12211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82146, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+6237,int_stack+144061,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100855, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8883,int_stack+6237,int_stack+4347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102178, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+3024,int_stack+13555,int_stack+12799, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82734, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+149686,int_stack+3024,int_stack+144061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104068, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+153214,int_stack+149686,int_stack+6237, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111655, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+3024,int_stack+153214,int_stack+8883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114301, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+149686,int_stack+3024,int_stack+146311,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+144061,int_stack+14920,int_stack+14500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51620, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+145321,int_stack+15508,int_stack+14920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52460, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+147085,int_stack+145321,int_stack+144061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+7749,int_stack+16292,int_stack+15508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53636, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10101,int_stack+7749,int_stack+145321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+122323,int_stack+10101,int_stack+147085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+117451, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+144061,int_stack+17300,int_stack+16292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55204, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+13629,int_stack+144061,int_stack+7749, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119971, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+174784,int_stack+13629,int_stack+10101, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105832, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+7749,int_stack+174784,int_stack+122323, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139861, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+174784,int_stack+7749,int_stack+3024,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+18785,int_stack+18560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77040, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3699,int_stack+19100,int_stack+18785, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77265, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4644,int_stack+3699,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+5994,int_stack+19520,int_stack+19100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77580, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7254,int_stack+5994,int_stack+3699, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92215, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+9144,int_stack+7254,int_stack+4644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93160, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+3024,int_stack+20060,int_stack+19520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78000, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+11394,int_stack+3024,int_stack+5994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+94510, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3024,int_stack+11394,int_stack+7254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95770, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+11394,int_stack+3024,int_stack+9144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97660, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+21050,int_stack+20735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81390, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3969,int_stack+21491,int_stack+21050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81705, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5292,int_stack+3969,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99910, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+7182,int_stack+22079,int_stack+21491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82146, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+14769,int_stack+7182,int_stack+3969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100855, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+17415,int_stack+14769,int_stack+5292, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102178, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+3024,int_stack+22835,int_stack+22079, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82734, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+122323,int_stack+3024,int_stack+7182, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104068, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3024,int_stack+122323,int_stack+14769, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111655, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+122323,int_stack+3024,int_stack+17415, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114301, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+188959,int_stack+122323,int_stack+11394,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+24200,int_stack+23780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51620, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4284,int_stack+24788,int_stack+24200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52460, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6048,int_stack+4284,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8568,int_stack+25572,int_stack+24788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53636, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10920,int_stack+8568,int_stack+4284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+14448,int_stack+10920,int_stack+6048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+117451, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+3024,int_stack+26580,int_stack+25572, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55204, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+18648,int_stack+3024,int_stack+8568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119971, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3024,int_stack+18648,int_stack+10920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105832, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+18648,int_stack+3024,int_stack+14448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139861, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+3024,int_stack+18648,int_stack+122323,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+122323,int_stack+28065,int_stack+27840, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+122998,int_stack+28380,int_stack+28065, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+123943,int_stack+122998,int_stack+122323, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+125293,int_stack+28800,int_stack+28380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+126553,int_stack+125293,int_stack+122998, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+17199,int_stack+126553,int_stack+123943, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+122323,int_stack+29340,int_stack+28800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+19449,int_stack+122323,int_stack+125293, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+94510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+122323,int_stack+19449,int_stack+126553, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+19449,int_stack+122323,int_stack+17199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17199,int_stack+30330,int_stack+30015, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+122323,int_stack+30771,int_stack+30330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+123646,int_stack+122323,int_stack+17199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+17199,int_stack+31359,int_stack+30771, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+125536,int_stack+17199,int_stack+122323, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+22824,int_stack+125536,int_stack+123646, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+122323,int_stack+32115,int_stack+31359, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+25974,int_stack+122323,int_stack+17199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+144061,int_stack+25974,int_stack+125536, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+25974,int_stack+144061,int_stack+22824, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114301, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+199084,int_stack+25974,int_stack+19449,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+144061,int_stack+33480,int_stack+33060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+145321,int_stack+34068,int_stack+33480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+147085,int_stack+145321,int_stack+144061, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+30699,int_stack+34852,int_stack+34068, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+17199,int_stack+30699,int_stack+145321, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+20727,int_stack+17199,int_stack+147085, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+117451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+144061,int_stack+35860,int_stack+34852, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+122323,int_stack+144061,int_stack+30699, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119971, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+30699,int_stack+122323,int_stack+17199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+122323,int_stack+30699,int_stack+20727, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139861, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+209209,int_stack+122323,int_stack+25974,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+122323,int_stack+37345,int_stack+37120, 0.0, zero_stack, 1.0, int_stack+77040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+122998,int_stack+37660,int_stack+37345, 0.0, zero_stack, 1.0, int_stack+77265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+123943,int_stack+122998,int_stack+122323, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+125293,int_stack+38080,int_stack+37660, 0.0, zero_stack, 1.0, int_stack+77580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+126553,int_stack+125293,int_stack+122998, 0.0, zero_stack, 1.0, int_stack+92215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+17199,int_stack+126553,int_stack+123943, 0.0, zero_stack, 1.0, int_stack+93160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+122323,int_stack+38620,int_stack+38080, 0.0, zero_stack, 1.0, int_stack+78000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+19449,int_stack+122323,int_stack+125293, 0.0, zero_stack, 1.0, int_stack+94510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+122323,int_stack+19449,int_stack+126553, 0.0, zero_stack, 1.0, int_stack+95770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+19449,int_stack+122323,int_stack+17199, 0.0, zero_stack, 1.0, int_stack+97660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17199,int_stack+39610,int_stack+39295, 0.0, zero_stack, 1.0, int_stack+81390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+122323,int_stack+40051,int_stack+39610, 0.0, zero_stack, 1.0, int_stack+81705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+123646,int_stack+122323,int_stack+17199, 0.0, zero_stack, 1.0, int_stack+99910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+17199,int_stack+40639,int_stack+40051, 0.0, zero_stack, 1.0, int_stack+82146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+125536,int_stack+17199,int_stack+122323, 0.0, zero_stack, 1.0, int_stack+100855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+22824,int_stack+125536,int_stack+123646, 0.0, zero_stack, 1.0, int_stack+102178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+122323,int_stack+41395,int_stack+40639, 0.0, zero_stack, 1.0, int_stack+82734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+25974,int_stack+122323,int_stack+17199, 0.0, zero_stack, 1.0, int_stack+104068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+29502,int_stack+25974,int_stack+125536, 0.0, zero_stack, 1.0, int_stack+111655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+122323,int_stack+29502,int_stack+22824, 0.0, zero_stack, 1.0, int_stack+114301, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+22824,int_stack+122323,int_stack+19449,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32949,int_stack+42760,int_stack+42340, 0.0, zero_stack, 1.0, int_stack+51620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+34209,int_stack+43348,int_stack+42760, 0.0, zero_stack, 1.0, int_stack+52460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+35973,int_stack+34209,int_stack+32949, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+38493,int_stack+44132,int_stack+43348, 0.0, zero_stack, 1.0, int_stack+53636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+17199,int_stack+38493,int_stack+34209, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+144061,int_stack+17199,int_stack+35973, 0.0, zero_stack, 1.0, int_stack+117451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+32949,int_stack+45140,int_stack+44132, 0.0, zero_stack, 1.0, int_stack+55204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+40845,int_stack+32949,int_stack+38493, 0.0, zero_stack, 1.0, int_stack+119971, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+32949,int_stack+40845,int_stack+17199, 0.0, zero_stack, 1.0, int_stack+105832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+38829,int_stack+32949,int_stack+144061, 0.0, zero_stack, 1.0, int_stack+139861, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+223384,int_stack+38829,int_stack+122323,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+122323,int_stack+46625,int_stack+46400, 1.0, int_stack+77040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+122998,int_stack+46940,int_stack+46625, 1.0, int_stack+77265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+123943,int_stack+122998,int_stack+122323, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+125293,int_stack+47360,int_stack+46940, 1.0, int_stack+77580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+126553,int_stack+125293,int_stack+122998, 1.0, int_stack+92215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+144061,int_stack+126553,int_stack+123943, 1.0, int_stack+93160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+91540,int_stack+47900,int_stack+47360, 1.0, int_stack+78000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+122323,int_stack+91540,int_stack+125293, 1.0, int_stack+94510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+146311,int_stack+122323,int_stack+126553, 1.0, int_stack+95770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+122323,int_stack+146311,int_stack+144061, 1.0, int_stack+97660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+144061,int_stack+48890,int_stack+48575, 1.0, int_stack+81390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+145006,int_stack+49331,int_stack+48890, 1.0, int_stack+81705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+146329,int_stack+145006,int_stack+144061, 1.0, int_stack+99910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+125698,int_stack+49919,int_stack+49331, 1.0, int_stack+82146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+91540,int_stack+125698,int_stack+145006, 1.0, int_stack+100855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+32949,int_stack+91540,int_stack+146329, 1.0, int_stack+102178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+99910,int_stack+50675,int_stack+49919, 1.0, int_stack+82734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+144061,int_stack+99910,int_stack+125698, 1.0, int_stack+104068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+36099,int_stack+144061,int_stack+91540, 1.0, int_stack+111655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+144061,int_stack+36099,int_stack+32949, 1.0, int_stack+114301, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+32949,int_stack+144061,int_stack+122323,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+122323,int_stack+53048,int_stack+52040, 1.0, int_stack+51620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+123583,int_stack+54420,int_stack+53048, 1.0, int_stack+52460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+125347,int_stack+123583,int_stack+122323, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+43074,int_stack+56212,int_stack+54420, 1.0, int_stack+53636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+45426,int_stack+43074,int_stack+123583, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+48954,int_stack+45426,int_stack+125347, 1.0, int_stack+117451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+57220,int_stack+56212, 1.0, int_stack+55204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+122323,int_stack+0,int_stack+43074, 1.0, int_stack+119971, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+237559,int_stack+122323,int_stack+45426, 1.0, int_stack+105832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+117451,int_stack+237559,int_stack+48954, 1.0, int_stack+139861, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+43074,int_stack+117451,int_stack+144061,225);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+117451,int_stack+78540,int_stack+78000,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+119071,int_stack+117451,int_stack+94510,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+121591,int_stack+119071,int_stack+95770,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+117451,int_stack+121591,int_stack+97660,15);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+120826,int_stack+83490,int_stack+82734,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+123094,int_stack+120826,int_stack+104068,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+139861,int_stack+123094,int_stack+111655,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+120826,int_stack+139861,int_stack+114301,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+91540,int_stack+120826,int_stack+117451,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+139861,int_stack+58705,int_stack+58480,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+140536,int_stack+59020,int_stack+58705,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+141481,int_stack+140536,int_stack+139861,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+142831,int_stack+59440,int_stack+59020,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+144091,int_stack+142831,int_stack+140536,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+145981,int_stack+144091,int_stack+141481,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+139861,int_stack+59980,int_stack+59440,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+125551,int_stack+139861,int_stack+142831,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+139861,int_stack+125551,int_stack+144091,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+125551,int_stack+139861,int_stack+145981,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+139861,int_stack+60970,int_stack+60655,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+140806,int_stack+61411,int_stack+60970,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+142129,int_stack+140806,int_stack+139861,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+144019,int_stack+61999,int_stack+61411,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+145783,int_stack+144019,int_stack+140806,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+237559,int_stack+145783,int_stack+142129,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+139861,int_stack+62755,int_stack+61999,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+57249,int_stack+139861,int_stack+144019,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+139861,int_stack+57249,int_stack+145783,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+57249,int_stack+139861,int_stack+237559,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+101665,int_stack+57249,int_stack+125551, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+117451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+125551,int_stack+64120,int_stack+63700,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+126811,int_stack+64708,int_stack+64120,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+237559,int_stack+126811,int_stack+125551,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+240079,int_stack+65492,int_stack+64708,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+139861,int_stack+240079,int_stack+126811,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+143389,int_stack+139861,int_stack+237559,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+66500,int_stack+65492,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+61974,int_stack+0,int_stack+240079,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+237559,int_stack+61974,int_stack+139861,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+243439,int_stack+237559,int_stack+143389,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+249739,int_stack+243439,int_stack+57249, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+57249,int_stack+67985,int_stack+67760,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+57924,int_stack+68300,int_stack+67985,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+58869,int_stack+57924,int_stack+57249,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+60219,int_stack+68720,int_stack+68300,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+61479,int_stack+60219,int_stack+57924,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+63369,int_stack+61479,int_stack+58869,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+57249,int_stack+69260,int_stack+68720,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+65619,int_stack+57249,int_stack+60219,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+57249,int_stack+65619,int_stack+61479,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+65619,int_stack+57249,int_stack+63369,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+57249,int_stack+70250,int_stack+69935,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+58194,int_stack+70691,int_stack+70250,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+59517,int_stack+58194,int_stack+57249,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+61407,int_stack+71279,int_stack+70691,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+237559,int_stack+61407,int_stack+58194,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+240205,int_stack+237559,int_stack+59517,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+57249,int_stack+72035,int_stack+71279,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+68994,int_stack+57249,int_stack+61407,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+57249,int_stack+68994,int_stack+237559,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+243355,int_stack+57249,int_stack+240205,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+263914,int_stack+243355,int_stack+65619, 0.0, zero_stack, 1.0, int_stack+117451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+57249,int_stack+73400,int_stack+72980,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+58509,int_stack+73988,int_stack+73400,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+60273,int_stack+58509,int_stack+57249,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+62793,int_stack+74772,int_stack+73988,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+65145,int_stack+62793,int_stack+58509,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+68673,int_stack+65145,int_stack+60273,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+75780,int_stack+74772,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+57249,int_stack+0,int_stack+62793,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+72873,int_stack+57249,int_stack+65145,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+57249,int_stack+72873,int_stack+68673,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+63549,int_stack+57249,int_stack+243355, 0.0, zero_stack, 1.0, int_stack+120826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+57249,int_stack+79440,int_stack+79215,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+57924,int_stack+79755,int_stack+79440,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+58869,int_stack+57924,int_stack+57249,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+60219,int_stack+80175,int_stack+79755,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+61479,int_stack+60219,int_stack+57924,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+77724,int_stack+61479,int_stack+58869,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+57249,int_stack+80715,int_stack+80175,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+79974,int_stack+57249,int_stack+60219,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+57249,int_stack+79974,int_stack+61479,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+79974,int_stack+57249,int_stack+77724,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+77724,int_stack+84750,int_stack+84435,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+57249,int_stack+85191,int_stack+84750,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+58572,int_stack+57249,int_stack+77724,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+77724,int_stack+85779,int_stack+85191,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+60462,int_stack+77724,int_stack+57249,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+237559,int_stack+60462,int_stack+58572,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+57249,int_stack+86535,int_stack+85779,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+83349,int_stack+57249,int_stack+77724,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+240709,int_stack+83349,int_stack+60462,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+57249,int_stack+240709,int_stack+237559,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+237559,int_stack+57249,int_stack+79974, 1.0, int_stack+117451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+247684,int_stack+87900,int_stack+87480,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+77724,int_stack+88488,int_stack+87900,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+79488,int_stack+77724,int_stack+247684,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+82008,int_stack+89272,int_stack+88488,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+84360,int_stack+82008,int_stack+77724,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+139861,int_stack+84360,int_stack+79488,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+90280,int_stack+89272,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+144061,int_stack+0,int_stack+82008,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+77724,int_stack+144061,int_stack+84360,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+83604,int_stack+77724,int_stack+139861,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+274039,int_stack+83604,int_stack+57249, 1.0, int_stack+120826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+288214,int_stack+160609,int_stack+129736,225);
     Libderiv->ABCD[11] = int_stack + 288214;
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+111790,int_stack+174784,int_stack+149686,225);
     Libderiv->ABCD[10] = int_stack + 111790;
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+132040,int_stack+3024,int_stack+188959,225);
     Libderiv->ABCD[9] = int_stack + 132040;
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+0,int_stack+209209,int_stack+199084,225);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+152290,int_stack+223384,int_stack+22824,225);
     Libderiv->ABCD[7] = int_stack + 152290;
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+172540,int_stack+43074,int_stack+32949,225);
     Libderiv->ABCD[6] = int_stack + 172540;
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+20250,int_stack+249739,int_stack+101665, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 20250;
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+40500,int_stack+63549,int_stack+263914, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 40500;
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+60750,int_stack+274039,int_stack+237559, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 60750;

}
