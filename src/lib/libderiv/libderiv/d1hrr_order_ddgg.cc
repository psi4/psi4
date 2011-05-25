#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ddgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|gg) integrals */

void d1hrr_order_ddgg(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[2][8][11] = int_stack + 600;
 Libderiv->deriv_classes[3][4][11] = int_stack + 870;
 Libderiv->deriv_classes[3][5][11] = int_stack + 1020;
 Libderiv->deriv_classes[3][6][11] = int_stack + 1230;
 Libderiv->deriv_classes[3][7][11] = int_stack + 1510;
 Libderiv->deriv_classes[3][8][11] = int_stack + 1870;
 Libderiv->deriv_classes[4][4][11] = int_stack + 2320;
 Libderiv->deriv_classes[4][5][11] = int_stack + 2545;
 Libderiv->deriv_classes[4][6][11] = int_stack + 2860;
 Libderiv->deriv_classes[4][7][11] = int_stack + 3280;
 Libderiv->deriv_classes[4][8][11] = int_stack + 3820;
 Libderiv->deriv_classes[2][4][10] = int_stack + 4495;
 Libderiv->deriv_classes[2][5][10] = int_stack + 4585;
 Libderiv->deriv_classes[2][6][10] = int_stack + 4711;
 Libderiv->deriv_classes[2][7][10] = int_stack + 4879;
 Libderiv->deriv_classes[2][8][10] = int_stack + 5095;
 Libderiv->deriv_classes[3][4][10] = int_stack + 5365;
 Libderiv->deriv_classes[3][5][10] = int_stack + 5515;
 Libderiv->deriv_classes[3][6][10] = int_stack + 5725;
 Libderiv->deriv_classes[3][7][10] = int_stack + 6005;
 Libderiv->deriv_classes[3][8][10] = int_stack + 6365;
 Libderiv->deriv_classes[4][4][10] = int_stack + 6815;
 Libderiv->deriv_classes[4][5][10] = int_stack + 7040;
 Libderiv->deriv_classes[4][6][10] = int_stack + 7355;
 Libderiv->deriv_classes[4][7][10] = int_stack + 7775;
 Libderiv->deriv_classes[4][8][10] = int_stack + 8315;
 Libderiv->deriv_classes[2][4][9] = int_stack + 8990;
 Libderiv->deriv_classes[2][5][9] = int_stack + 9080;
 Libderiv->deriv_classes[2][6][9] = int_stack + 9206;
 Libderiv->deriv_classes[2][7][9] = int_stack + 9374;
 Libderiv->deriv_classes[2][8][9] = int_stack + 9590;
 Libderiv->deriv_classes[3][4][9] = int_stack + 9860;
 Libderiv->deriv_classes[3][5][9] = int_stack + 10010;
 Libderiv->deriv_classes[3][6][9] = int_stack + 10220;
 Libderiv->deriv_classes[3][7][9] = int_stack + 10500;
 Libderiv->deriv_classes[3][8][9] = int_stack + 10860;
 Libderiv->deriv_classes[4][4][9] = int_stack + 11310;
 Libderiv->deriv_classes[4][5][9] = int_stack + 11535;
 Libderiv->deriv_classes[4][6][9] = int_stack + 11850;
 Libderiv->deriv_classes[4][7][9] = int_stack + 12270;
 Libderiv->deriv_classes[4][8][9] = int_stack + 12810;
 Libderiv->deriv_classes[2][4][8] = int_stack + 13485;
 Libderiv->deriv_classes[2][5][8] = int_stack + 13575;
 Libderiv->deriv_classes[2][6][8] = int_stack + 13701;
 Libderiv->deriv_classes[2][7][8] = int_stack + 13869;
 Libderiv->deriv_classes[2][8][8] = int_stack + 14085;
 Libderiv->deriv_classes[3][4][8] = int_stack + 14355;
 Libderiv->deriv_classes[3][5][8] = int_stack + 14505;
 Libderiv->deriv_classes[3][6][8] = int_stack + 14715;
 Libderiv->deriv_classes[3][7][8] = int_stack + 14995;
 Libderiv->deriv_classes[3][8][8] = int_stack + 15355;
 Libderiv->deriv_classes[4][4][8] = int_stack + 15805;
 Libderiv->deriv_classes[4][5][8] = int_stack + 16030;
 Libderiv->deriv_classes[4][6][8] = int_stack + 16345;
 Libderiv->deriv_classes[4][7][8] = int_stack + 16765;
 Libderiv->deriv_classes[4][8][8] = int_stack + 17305;
 Libderiv->deriv_classes[2][4][7] = int_stack + 17980;
 Libderiv->deriv_classes[2][5][7] = int_stack + 18070;
 Libderiv->deriv_classes[2][6][7] = int_stack + 18196;
 Libderiv->deriv_classes[2][7][7] = int_stack + 18364;
 Libderiv->deriv_classes[2][8][7] = int_stack + 18580;
 Libderiv->deriv_classes[3][4][7] = int_stack + 18850;
 Libderiv->deriv_classes[3][5][7] = int_stack + 19000;
 Libderiv->deriv_classes[3][6][7] = int_stack + 19210;
 Libderiv->deriv_classes[3][7][7] = int_stack + 19490;
 Libderiv->deriv_classes[3][8][7] = int_stack + 19850;
 Libderiv->deriv_classes[4][4][7] = int_stack + 20300;
 Libderiv->deriv_classes[4][5][7] = int_stack + 20525;
 Libderiv->deriv_classes[4][6][7] = int_stack + 20840;
 Libderiv->deriv_classes[4][7][7] = int_stack + 21260;
 Libderiv->deriv_classes[4][8][7] = int_stack + 21800;
 Libderiv->deriv_classes[2][4][6] = int_stack + 22475;
 Libderiv->deriv_classes[2][5][6] = int_stack + 22565;
 Libderiv->deriv_classes[2][6][6] = int_stack + 22691;
 Libderiv->deriv_classes[2][7][6] = int_stack + 22859;
 Libderiv->deriv_classes[2][8][6] = int_stack + 23075;
 Libderiv->deriv_classes[3][4][6] = int_stack + 23345;
 Libderiv->deriv_classes[3][5][6] = int_stack + 23495;
 Libderiv->deriv_classes[3][6][6] = int_stack + 23705;
 Libderiv->deriv_classes[3][7][6] = int_stack + 23985;
 Libderiv->deriv_classes[3][8][6] = int_stack + 24345;
 Libderiv->dvrr_classes[4][4] = int_stack + 24795;
 Libderiv->deriv_classes[4][4][6] = int_stack + 25020;
 Libderiv->dvrr_classes[4][5] = int_stack + 25245;
 Libderiv->deriv_classes[4][5][6] = int_stack + 25560;
 Libderiv->dvrr_classes[4][6] = int_stack + 25875;
 Libderiv->deriv_classes[4][6][6] = int_stack + 26295;
 Libderiv->dvrr_classes[4][7] = int_stack + 26715;
 Libderiv->deriv_classes[4][7][6] = int_stack + 27255;
 Libderiv->deriv_classes[4][8][6] = int_stack + 27795;
 Libderiv->deriv_classes[2][4][2] = int_stack + 28470;
 Libderiv->deriv_classes[2][5][2] = int_stack + 28560;
 Libderiv->deriv_classes[2][6][2] = int_stack + 28686;
 Libderiv->deriv_classes[2][7][2] = int_stack + 28854;
 Libderiv->deriv_classes[2][8][2] = int_stack + 29070;
 Libderiv->deriv_classes[3][4][2] = int_stack + 29340;
 Libderiv->deriv_classes[3][5][2] = int_stack + 29490;
 Libderiv->deriv_classes[3][6][2] = int_stack + 29700;
 Libderiv->deriv_classes[3][7][2] = int_stack + 29980;
 Libderiv->deriv_classes[3][8][2] = int_stack + 30340;
 Libderiv->deriv_classes[4][4][2] = int_stack + 30790;
 Libderiv->deriv_classes[4][5][2] = int_stack + 31015;
 Libderiv->deriv_classes[4][6][2] = int_stack + 31330;
 Libderiv->deriv_classes[4][7][2] = int_stack + 31750;
 Libderiv->deriv_classes[4][8][2] = int_stack + 32290;
 Libderiv->deriv_classes[2][4][1] = int_stack + 32965;
 Libderiv->deriv_classes[2][5][1] = int_stack + 33055;
 Libderiv->deriv_classes[2][6][1] = int_stack + 33181;
 Libderiv->deriv_classes[2][7][1] = int_stack + 33349;
 Libderiv->deriv_classes[2][8][1] = int_stack + 33565;
 Libderiv->deriv_classes[3][4][1] = int_stack + 33835;
 Libderiv->deriv_classes[3][5][1] = int_stack + 33985;
 Libderiv->deriv_classes[3][6][1] = int_stack + 34195;
 Libderiv->deriv_classes[3][7][1] = int_stack + 34475;
 Libderiv->deriv_classes[3][8][1] = int_stack + 34835;
 Libderiv->deriv_classes[4][4][1] = int_stack + 35285;
 Libderiv->deriv_classes[4][5][1] = int_stack + 35510;
 Libderiv->deriv_classes[4][6][1] = int_stack + 35825;
 Libderiv->deriv_classes[4][7][1] = int_stack + 36245;
 Libderiv->deriv_classes[4][8][1] = int_stack + 36785;
 Libderiv->dvrr_classes[2][4] = int_stack + 37460;
 Libderiv->dvrr_classes[2][5] = int_stack + 37550;
 Libderiv->dvrr_classes[2][6] = int_stack + 37676;
 Libderiv->dvrr_classes[2][7] = int_stack + 37844;
 Libderiv->dvrr_classes[2][8] = int_stack + 38060;
 Libderiv->deriv_classes[2][4][0] = int_stack + 38330;
 Libderiv->deriv_classes[2][5][0] = int_stack + 38420;
 Libderiv->deriv_classes[2][6][0] = int_stack + 38546;
 Libderiv->deriv_classes[2][7][0] = int_stack + 38714;
 Libderiv->deriv_classes[2][8][0] = int_stack + 38930;
 Libderiv->dvrr_classes[3][4] = int_stack + 39200;
 Libderiv->dvrr_classes[3][5] = int_stack + 39350;
 Libderiv->dvrr_classes[3][6] = int_stack + 39560;
 Libderiv->dvrr_classes[3][7] = int_stack + 39840;
 Libderiv->dvrr_classes[3][8] = int_stack + 40200;
 Libderiv->deriv_classes[3][4][0] = int_stack + 40650;
 Libderiv->deriv_classes[3][5][0] = int_stack + 40800;
 Libderiv->deriv_classes[3][6][0] = int_stack + 41010;
 Libderiv->deriv_classes[3][7][0] = int_stack + 41290;
 Libderiv->deriv_classes[3][8][0] = int_stack + 41650;
 Libderiv->deriv_classes[4][4][0] = int_stack + 42100;
 Libderiv->deriv_classes[4][5][0] = int_stack + 42325;
 Libderiv->deriv_classes[4][6][0] = int_stack + 42640;
 Libderiv->deriv_classes[4][7][0] = int_stack + 43060;
 Libderiv->deriv_classes[4][8][0] = int_stack + 43600;
 memset(int_stack,0,354200);

 Libderiv->dvrr_stack = int_stack + 135487;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ddgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+44275,int_stack+37550,int_stack+37460,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+44545,int_stack+37676,int_stack+37550,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+44923,int_stack+44545,int_stack+44275,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+45463,int_stack+37844,int_stack+37676,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+45967,int_stack+45463,int_stack+44545,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+46723,int_stack+45967,int_stack+44923,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47623,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37460,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47893,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37550,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48271,int_stack+47893,int_stack+47623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44275,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+48811,int_stack+384,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37676,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+49315,int_stack+48811,int_stack+47893, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44545,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+50071,int_stack+49315,int_stack+48271, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44923,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+47623,int_stack+600,int_stack+384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37844,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+50971,int_stack+47623,int_stack+48811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45463,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+47623,int_stack+50971,int_stack+49315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45967,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+50971,int_stack+47623,int_stack+50071, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46723,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+47623,int_stack+39350,int_stack+39200,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+48073,int_stack+39560,int_stack+39350,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+48703,int_stack+48073,int_stack+47623,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+49603,int_stack+39840,int_stack+39560,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+52321,int_stack+49603,int_stack+48073,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+53581,int_stack+52321,int_stack+48703,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+50443,int_stack+1020,int_stack+870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39200,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55081,int_stack+1230,int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39350,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+55081,int_stack+50443, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47623,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+55711,int_stack+1510,int_stack+1230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39560,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+56551,int_stack+55711,int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48073,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+57811,int_stack+56551,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48703,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+1870,int_stack+1510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39840,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+59311,int_stack+0,int_stack+55711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49603,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+59311,int_stack+56551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52321,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+55081,int_stack+0,int_stack+57811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53581,10);
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+57331,int_stack+55081,int_stack+50971,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+25245,int_stack+24795,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+675,int_stack+25875,int_stack+25245,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+50443,int_stack+675,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+61381,int_stack+26715,int_stack+25875,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+62641,int_stack+61381,int_stack+675,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+64531,int_stack+62641,int_stack+50443,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1620,int_stack+2545,int_stack+2320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24795,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66781,int_stack+2860,int_stack+2545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25245,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67726,int_stack+66781,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+69076,int_stack+3280,int_stack+2860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25875,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+70336,int_stack+69076,int_stack+66781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+675,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+72226,int_stack+70336,int_stack+67726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50443,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+66781,int_stack+3820,int_stack+3280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26715,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1620,int_stack+66781,int_stack+69076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61381,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+66781,int_stack+1620,int_stack+70336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62641,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+74476,int_stack+66781,int_stack+72226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64531,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+66781,int_stack+74476,int_stack+55081,225);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+55081,int_stack+4585,int_stack+4495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37460, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55351,int_stack+4711,int_stack+4585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37550, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55729,int_stack+55351,int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44275, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+56269,int_stack+4879,int_stack+4711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37676, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+73531,int_stack+56269,int_stack+55351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44545, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+74287,int_stack+73531,int_stack+55729, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44923, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+55081,int_stack+5095,int_stack+4879, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37844, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+75187,int_stack+55081,int_stack+56269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45463, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+55081,int_stack+75187,int_stack+73531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45967, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+75187,int_stack+55081,int_stack+74287, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46723, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+55081,int_stack+5515,int_stack+5365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39200, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55531,int_stack+5725,int_stack+5515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39350, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56161,int_stack+55531,int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47623, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+73531,int_stack+6005,int_stack+5725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39560, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+76537,int_stack+73531,int_stack+55531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48073, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1620,int_stack+76537,int_stack+56161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48703, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+55081,int_stack+6365,int_stack+6005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39840, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3120,int_stack+55081,int_stack+73531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49603, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+55081,int_stack+3120,int_stack+76537, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52321, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+3120,int_stack+55081,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53581, 0.0, zero_stack,10);
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+76537,int_stack+3120,int_stack+75187,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1620,int_stack+7040,int_stack+6815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24795, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55081,int_stack+7355,int_stack+7040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25245, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5370,int_stack+55081,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1620,int_stack+7775,int_stack+7355, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25875, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+73531,int_stack+1620,int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+675, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+55081,int_stack+73531,int_stack+5370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50443, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+5370,int_stack+8315,int_stack+7775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26715, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+80587,int_stack+5370,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61381, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+5370,int_stack+80587,int_stack+73531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62641, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+80587,int_stack+5370,int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64531, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+83962,int_stack+80587,int_stack+3120,225);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80587,int_stack+9080,int_stack+8990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37460, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80857,int_stack+9206,int_stack+9080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37550, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81235,int_stack+80857,int_stack+80587, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44275, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+81775,int_stack+9374,int_stack+9206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37676, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+82279,int_stack+81775,int_stack+80857, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44545, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+83035,int_stack+82279,int_stack+81235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44923, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+80587,int_stack+9590,int_stack+9374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37844, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+55081,int_stack+80587,int_stack+81775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45463, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+80587,int_stack+55081,int_stack+82279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45967, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+55081,int_stack+80587,int_stack+83035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46723, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80587,int_stack+10010,int_stack+9860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39200, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81037,int_stack+10220,int_stack+10010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39350, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56431,int_stack+81037,int_stack+80587, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47623, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+81667,int_stack+10500,int_stack+10220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39560, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+82507,int_stack+81667,int_stack+81037, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48073, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+73531,int_stack+82507,int_stack+56431, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48703, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+80587,int_stack+10860,int_stack+10500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39840, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1620,int_stack+80587,int_stack+81667, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49603, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3300,int_stack+1620,int_stack+82507, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52321, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+80587,int_stack+3300,int_stack+73531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53581, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1620,int_stack+80587,int_stack+55081,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+55081,int_stack+11535,int_stack+11310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24795, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55756,int_stack+11850,int_stack+11535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25245, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+73531,int_stack+55756,int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+74881,int_stack+12270,int_stack+11850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25875, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+5670,int_stack+74881,int_stack+55756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+675, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+55081,int_stack+5670,int_stack+73531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50443, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+7560,int_stack+12810,int_stack+12270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26715, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+9180,int_stack+7560,int_stack+74881, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61381, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+90712,int_stack+9180,int_stack+5670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62641, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+5670,int_stack+90712,int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64531, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+90712,int_stack+5670,int_stack+80587,225);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80587,int_stack+13575,int_stack+13485, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80857,int_stack+13701,int_stack+13575, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81235,int_stack+80857,int_stack+80587, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+81775,int_stack+13869,int_stack+13701, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+82279,int_stack+81775,int_stack+80857, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+83035,int_stack+82279,int_stack+81235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+80587,int_stack+14085,int_stack+13869, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+5670,int_stack+80587,int_stack+81775, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45463, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+80587,int_stack+5670,int_stack+82279, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+5670,int_stack+80587,int_stack+83035, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46723, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80587,int_stack+14505,int_stack+14355, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81037,int_stack+14715,int_stack+14505, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81667,int_stack+81037,int_stack+80587, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+82567,int_stack+14995,int_stack+14715, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7020,int_stack+82567,int_stack+81037, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8280,int_stack+7020,int_stack+81667, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48703, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+80587,int_stack+15355,int_stack+14995, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+9780,int_stack+80587,int_stack+82567, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+80587,int_stack+9780,int_stack+7020, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+55081,int_stack+80587,int_stack+8280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53581, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7020,int_stack+55081,int_stack+5670,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5670,int_stack+16030,int_stack+15805, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80587,int_stack+16345,int_stack+16030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81532,int_stack+80587,int_stack+5670, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+5670,int_stack+16765,int_stack+16345, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+11070,int_stack+5670,int_stack+80587, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+12960,int_stack+11070,int_stack+81532, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50443, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+80587,int_stack+17305,int_stack+16765, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+15210,int_stack+80587,int_stack+5670, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+80587,int_stack+15210,int_stack+11070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+97462,int_stack+80587,int_stack+12960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11070,int_stack+97462,int_stack+55081,225);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+55081,int_stack+18070,int_stack+17980, 0.0, zero_stack, 1.0, int_stack+37460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55351,int_stack+18196,int_stack+18070, 0.0, zero_stack, 1.0, int_stack+37550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55729,int_stack+55351,int_stack+55081, 0.0, zero_stack, 1.0, int_stack+44275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+56269,int_stack+18364,int_stack+18196, 0.0, zero_stack, 1.0, int_stack+37676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+97462,int_stack+56269,int_stack+55351, 0.0, zero_stack, 1.0, int_stack+44545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+98218,int_stack+97462,int_stack+55729, 0.0, zero_stack, 1.0, int_stack+44923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+55081,int_stack+18580,int_stack+18364, 0.0, zero_stack, 1.0, int_stack+37844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+99118,int_stack+55081,int_stack+56269, 0.0, zero_stack, 1.0, int_stack+45463, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+55081,int_stack+99118,int_stack+97462, 0.0, zero_stack, 1.0, int_stack+45967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+5670,int_stack+55081,int_stack+98218, 0.0, zero_stack, 1.0, int_stack+46723, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+55081,int_stack+19000,int_stack+18850, 0.0, zero_stack, 1.0, int_stack+39200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55531,int_stack+19210,int_stack+19000, 0.0, zero_stack, 1.0, int_stack+39350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56161,int_stack+55531,int_stack+55081, 0.0, zero_stack, 1.0, int_stack+47623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+97462,int_stack+19490,int_stack+19210, 0.0, zero_stack, 1.0, int_stack+39560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+98302,int_stack+97462,int_stack+55531, 0.0, zero_stack, 1.0, int_stack+48073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+80587,int_stack+98302,int_stack+56161, 0.0, zero_stack, 1.0, int_stack+48703, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+55081,int_stack+19850,int_stack+19490, 0.0, zero_stack, 1.0, int_stack+39840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+82087,int_stack+55081,int_stack+97462, 0.0, zero_stack, 1.0, int_stack+49603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+55081,int_stack+82087,int_stack+98302, 0.0, zero_stack, 1.0, int_stack+52321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+97462,int_stack+55081,int_stack+80587, 0.0, zero_stack, 1.0, int_stack+53581, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+99712,int_stack+97462,int_stack+5670,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5670,int_stack+20525,int_stack+20300, 0.0, zero_stack, 1.0, int_stack+24795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80587,int_stack+20840,int_stack+20525, 0.0, zero_stack, 1.0, int_stack+25245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81532,int_stack+80587,int_stack+5670, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+5670,int_stack+21260,int_stack+20840, 0.0, zero_stack, 1.0, int_stack+25875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+55081,int_stack+5670,int_stack+80587, 0.0, zero_stack, 1.0, int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+17820,int_stack+55081,int_stack+81532, 0.0, zero_stack, 1.0, int_stack+50443, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+80587,int_stack+21800,int_stack+21260, 0.0, zero_stack, 1.0, int_stack+26715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+73531,int_stack+80587,int_stack+5670, 0.0, zero_stack, 1.0, int_stack+61381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+80587,int_stack+73531,int_stack+55081, 0.0, zero_stack, 1.0, int_stack+62641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+103762,int_stack+80587,int_stack+17820, 0.0, zero_stack, 1.0, int_stack+64531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+107137,int_stack+103762,int_stack+97462,225);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97462,int_stack+22565,int_stack+22475, 1.0, int_stack+37460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+97732,int_stack+22691,int_stack+22565, 1.0, int_stack+37550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+98110,int_stack+97732,int_stack+97462, 1.0, int_stack+44275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+98650,int_stack+22859,int_stack+22691, 1.0, int_stack+37676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+103762,int_stack+98650,int_stack+97732, 1.0, int_stack+44545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+104518,int_stack+103762,int_stack+98110, 1.0, int_stack+44923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+44275,int_stack+23075,int_stack+22859, 1.0, int_stack+37844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+97462,int_stack+44275,int_stack+98650, 1.0, int_stack+45463, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+105418,int_stack+97462,int_stack+103762, 1.0, int_stack+45967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+5670,int_stack+105418,int_stack+104518, 1.0, int_stack+46723, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+103762,int_stack+23495,int_stack+23345, 1.0, int_stack+39200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+104212,int_stack+23705,int_stack+23495, 1.0, int_stack+39350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+104842,int_stack+104212,int_stack+103762, 1.0, int_stack+47623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+105742,int_stack+23985,int_stack+23705, 1.0, int_stack+39560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+97462,int_stack+105742,int_stack+104212, 1.0, int_stack+48073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+17820,int_stack+97462,int_stack+104842, 1.0, int_stack+48703, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+47623,int_stack+24345,int_stack+23985, 1.0, int_stack+39840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+103762,int_stack+47623,int_stack+105742, 1.0, int_stack+49603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+19320,int_stack+103762,int_stack+97462, 1.0, int_stack+52321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+97462,int_stack+19320,int_stack+17820, 1.0, int_stack+53581, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+17820,int_stack+97462,int_stack+5670,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5670,int_stack+25560,int_stack+25020, 1.0, int_stack+24795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21870,int_stack+26295,int_stack+25560, 1.0, int_stack+25245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22815,int_stack+21870,int_stack+5670, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+5670,int_stack+27255,int_stack+26295, 1.0, int_stack+25875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+24165,int_stack+5670,int_stack+21870, 1.0, int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+55081,int_stack+24165,int_stack+22815, 1.0, int_stack+50443, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+27795,int_stack+27255, 1.0, int_stack+26715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+103762,int_stack+0,int_stack+5670, 1.0, int_stack+61381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+80587,int_stack+103762,int_stack+24165, 1.0, int_stack+62641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+103762,int_stack+80587,int_stack+55081, 1.0, int_stack+64531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+113887,int_stack+103762,int_stack+97462,225);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+97462,int_stack+38060,int_stack+37844,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+98110,int_stack+97462,int_stack+45463,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+103762,int_stack+98110,int_stack+45967,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+5670,int_stack+103762,int_stack+46723,6);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+103762,int_stack+40200,int_stack+39840,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+104842,int_stack+103762,int_stack+49603,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+97462,int_stack+104842,int_stack+52321,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+55081,int_stack+97462,int_stack+53581,10);
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+61381,int_stack+55081,int_stack+5670,225);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+97462,int_stack+28560,int_stack+28470,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+97732,int_stack+28686,int_stack+28560,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+98110,int_stack+97732,int_stack+97462,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+98650,int_stack+28854,int_stack+28686,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+103762,int_stack+98650,int_stack+97732,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+104518,int_stack+103762,int_stack+98110,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+97462,int_stack+29070,int_stack+28854,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+105418,int_stack+97462,int_stack+98650,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+97462,int_stack+105418,int_stack+103762,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+65431,int_stack+97462,int_stack+104518,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+97462,int_stack+29490,int_stack+29340,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+97912,int_stack+29700,int_stack+29490,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+98542,int_stack+97912,int_stack+97462,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+103762,int_stack+29980,int_stack+29700,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+104602,int_stack+103762,int_stack+97912,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+80587,int_stack+104602,int_stack+98542,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+97462,int_stack+30340,int_stack+29980,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+82087,int_stack+97462,int_stack+103762,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+97462,int_stack+82087,int_stack+104602,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+103762,int_stack+97462,int_stack+80587,10);
 /*--- compute (dp|gg) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+21870,int_stack+103762,int_stack+65431, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+65431,int_stack+31015,int_stack+30790,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+80587,int_stack+31330,int_stack+31015,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81532,int_stack+80587,int_stack+65431,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+65431,int_stack+31750,int_stack+31330,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+97462,int_stack+65431,int_stack+80587,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+25920,int_stack+97462,int_stack+81532,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+32290,int_stack+31750,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+80587,int_stack+0,int_stack+65431,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+28170,int_stack+80587,int_stack+97462,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+80587,int_stack+28170,int_stack+25920,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+25920,int_stack+80587,int_stack+103762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+103762,int_stack+33055,int_stack+32965,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+104032,int_stack+33181,int_stack+33055,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+104410,int_stack+104032,int_stack+103762,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+104950,int_stack+33349,int_stack+33181,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+105454,int_stack+104950,int_stack+104032,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+106210,int_stack+105454,int_stack+104410,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+103762,int_stack+33565,int_stack+33349,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+80587,int_stack+103762,int_stack+104950,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+103762,int_stack+80587,int_stack+105454,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+65431,int_stack+103762,int_stack+106210,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+103762,int_stack+33985,int_stack+33835,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+104212,int_stack+34195,int_stack+33985,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+104842,int_stack+104212,int_stack+103762,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+105742,int_stack+34475,int_stack+34195,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+80587,int_stack+105742,int_stack+104212,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+81847,int_stack+80587,int_stack+104842,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+103762,int_stack+34835,int_stack+34475,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+32670,int_stack+103762,int_stack+105742,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+103762,int_stack+32670,int_stack+80587,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+97462,int_stack+103762,int_stack+81847,10);
 /*--- compute (dp|gg) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+44275,int_stack+97462,int_stack+65431, 0.0, zero_stack, 1.0, int_stack+5670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+65431,int_stack+35510,int_stack+35285,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+103762,int_stack+35825,int_stack+35510,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+104707,int_stack+103762,int_stack+65431,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+65431,int_stack+36245,int_stack+35825,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+80587,int_stack+65431,int_stack+103762,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+32670,int_stack+80587,int_stack+104707,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+36785,int_stack+36245,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+103762,int_stack+0,int_stack+65431,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+34920,int_stack+103762,int_stack+80587,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+80587,int_stack+34920,int_stack+32670,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+48325,int_stack+80587,int_stack+97462, 0.0, zero_stack, 1.0, int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+97462,int_stack+38420,int_stack+38330,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+97732,int_stack+38546,int_stack+38420,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+98110,int_stack+97732,int_stack+97462,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+98650,int_stack+38714,int_stack+38546,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+80587,int_stack+98650,int_stack+97732,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+81343,int_stack+80587,int_stack+98110,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+97462,int_stack+38930,int_stack+38714,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+82243,int_stack+97462,int_stack+98650,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+97462,int_stack+82243,int_stack+80587,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+65431,int_stack+97462,int_stack+81343,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+97462,int_stack+40800,int_stack+40650,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+97912,int_stack+41010,int_stack+40800,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+98542,int_stack+97912,int_stack+97462,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+80587,int_stack+41290,int_stack+41010,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+81427,int_stack+80587,int_stack+97912,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+32670,int_stack+81427,int_stack+98542,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+97462,int_stack+41650,int_stack+41290,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+34170,int_stack+97462,int_stack+80587,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+97462,int_stack+34170,int_stack+81427,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+34170,int_stack+97462,int_stack+32670,10);
 /*--- compute (dp|gg) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+36420,int_stack+34170,int_stack+65431, 1.0, int_stack+5670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5670,int_stack+42325,int_stack+42100,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+65431,int_stack+42640,int_stack+42325,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+32670,int_stack+65431,int_stack+5670,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+5670,int_stack+43060,int_stack+42640,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+97462,int_stack+5670,int_stack+65431,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+40470,int_stack+97462,int_stack+32670,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+43600,int_stack+43060,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+80587,int_stack+0,int_stack+5670,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+103762,int_stack+80587,int_stack+97462,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+80587,int_stack+103762,int_stack+40470,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+120637,int_stack+80587,int_stack+34170, 1.0, int_stack+55081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (dd|gg) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+127387,int_stack+66781,int_stack+57331,225);
     Libderiv->ABCD[11] = int_stack + 127387;
 /*--- compute (dd|gg) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+65431,int_stack+83962,int_stack+76537,225);
     Libderiv->ABCD[10] = int_stack + 65431;
 /*--- compute (dd|gg) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+73531,int_stack+90712,int_stack+1620,225);
     Libderiv->ABCD[9] = int_stack + 73531;
 /*--- compute (dd|gg) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+81631,int_stack+11070,int_stack+7020,225);
     Libderiv->ABCD[8] = int_stack + 81631;
 /*--- compute (dd|gg) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+107137,int_stack+99712,225);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (dd|gg) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+8100,int_stack+113887,int_stack+17820,225);
     Libderiv->ABCD[6] = int_stack + 8100;
 /*--- compute (dd|gg) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+89731,int_stack+25920,int_stack+21870, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 89731;
 /*--- compute (dd|gg) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+16200,int_stack+48325,int_stack+44275, 0.0, zero_stack, 1.0, int_stack+61381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 16200;
 /*--- compute (dd|gg) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+40470,int_stack+120637,int_stack+36420, 1.0, int_stack+61381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 40470;

}
