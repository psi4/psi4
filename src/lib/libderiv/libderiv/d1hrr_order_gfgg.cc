#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gfgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gf|gg) integrals */

void d1hrr_order_gfgg(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[7][4][11] = int_stack + 9280;
 Libderiv->deriv_classes[7][5][11] = int_stack + 9820;
 Libderiv->deriv_classes[7][6][11] = int_stack + 10576;
 Libderiv->deriv_classes[7][7][11] = int_stack + 11584;
 Libderiv->deriv_classes[7][8][11] = int_stack + 12880;
 Libderiv->deriv_classes[4][4][10] = int_stack + 14500;
 Libderiv->deriv_classes[4][5][10] = int_stack + 14725;
 Libderiv->deriv_classes[4][6][10] = int_stack + 15040;
 Libderiv->deriv_classes[4][7][10] = int_stack + 15460;
 Libderiv->deriv_classes[4][8][10] = int_stack + 16000;
 Libderiv->deriv_classes[5][4][10] = int_stack + 16675;
 Libderiv->deriv_classes[5][5][10] = int_stack + 16990;
 Libderiv->deriv_classes[5][6][10] = int_stack + 17431;
 Libderiv->deriv_classes[5][7][10] = int_stack + 18019;
 Libderiv->deriv_classes[5][8][10] = int_stack + 18775;
 Libderiv->deriv_classes[6][4][10] = int_stack + 19720;
 Libderiv->deriv_classes[6][5][10] = int_stack + 20140;
 Libderiv->deriv_classes[6][6][10] = int_stack + 20728;
 Libderiv->deriv_classes[6][7][10] = int_stack + 21512;
 Libderiv->deriv_classes[6][8][10] = int_stack + 22520;
 Libderiv->deriv_classes[7][4][10] = int_stack + 23780;
 Libderiv->deriv_classes[7][5][10] = int_stack + 24320;
 Libderiv->deriv_classes[7][6][10] = int_stack + 25076;
 Libderiv->deriv_classes[7][7][10] = int_stack + 26084;
 Libderiv->deriv_classes[7][8][10] = int_stack + 27380;
 Libderiv->deriv_classes[4][4][9] = int_stack + 29000;
 Libderiv->deriv_classes[4][5][9] = int_stack + 29225;
 Libderiv->deriv_classes[4][6][9] = int_stack + 29540;
 Libderiv->deriv_classes[4][7][9] = int_stack + 29960;
 Libderiv->deriv_classes[4][8][9] = int_stack + 30500;
 Libderiv->deriv_classes[5][4][9] = int_stack + 31175;
 Libderiv->deriv_classes[5][5][9] = int_stack + 31490;
 Libderiv->deriv_classes[5][6][9] = int_stack + 31931;
 Libderiv->deriv_classes[5][7][9] = int_stack + 32519;
 Libderiv->deriv_classes[5][8][9] = int_stack + 33275;
 Libderiv->deriv_classes[6][4][9] = int_stack + 34220;
 Libderiv->deriv_classes[6][5][9] = int_stack + 34640;
 Libderiv->deriv_classes[6][6][9] = int_stack + 35228;
 Libderiv->deriv_classes[6][7][9] = int_stack + 36012;
 Libderiv->deriv_classes[6][8][9] = int_stack + 37020;
 Libderiv->deriv_classes[7][4][9] = int_stack + 38280;
 Libderiv->deriv_classes[7][5][9] = int_stack + 38820;
 Libderiv->deriv_classes[7][6][9] = int_stack + 39576;
 Libderiv->deriv_classes[7][7][9] = int_stack + 40584;
 Libderiv->deriv_classes[7][8][9] = int_stack + 41880;
 Libderiv->deriv_classes[4][4][8] = int_stack + 43500;
 Libderiv->deriv_classes[4][5][8] = int_stack + 43725;
 Libderiv->deriv_classes[4][6][8] = int_stack + 44040;
 Libderiv->deriv_classes[4][7][8] = int_stack + 44460;
 Libderiv->deriv_classes[4][8][8] = int_stack + 45000;
 Libderiv->deriv_classes[5][4][8] = int_stack + 45675;
 Libderiv->deriv_classes[5][5][8] = int_stack + 45990;
 Libderiv->deriv_classes[5][6][8] = int_stack + 46431;
 Libderiv->deriv_classes[5][7][8] = int_stack + 47019;
 Libderiv->deriv_classes[5][8][8] = int_stack + 47775;
 Libderiv->deriv_classes[6][4][8] = int_stack + 48720;
 Libderiv->deriv_classes[6][5][8] = int_stack + 49140;
 Libderiv->deriv_classes[6][6][8] = int_stack + 49728;
 Libderiv->deriv_classes[6][7][8] = int_stack + 50512;
 Libderiv->deriv_classes[6][8][8] = int_stack + 51520;
 Libderiv->deriv_classes[7][4][8] = int_stack + 52780;
 Libderiv->deriv_classes[7][5][8] = int_stack + 53320;
 Libderiv->deriv_classes[7][6][8] = int_stack + 54076;
 Libderiv->deriv_classes[7][7][8] = int_stack + 55084;
 Libderiv->deriv_classes[7][8][8] = int_stack + 56380;
 Libderiv->deriv_classes[4][4][7] = int_stack + 58000;
 Libderiv->deriv_classes[4][5][7] = int_stack + 58225;
 Libderiv->deriv_classes[4][6][7] = int_stack + 58540;
 Libderiv->deriv_classes[4][7][7] = int_stack + 58960;
 Libderiv->deriv_classes[4][8][7] = int_stack + 59500;
 Libderiv->deriv_classes[5][4][7] = int_stack + 60175;
 Libderiv->deriv_classes[5][5][7] = int_stack + 60490;
 Libderiv->deriv_classes[5][6][7] = int_stack + 60931;
 Libderiv->deriv_classes[5][7][7] = int_stack + 61519;
 Libderiv->deriv_classes[5][8][7] = int_stack + 62275;
 Libderiv->deriv_classes[6][4][7] = int_stack + 63220;
 Libderiv->deriv_classes[6][5][7] = int_stack + 63640;
 Libderiv->deriv_classes[6][6][7] = int_stack + 64228;
 Libderiv->deriv_classes[6][7][7] = int_stack + 65012;
 Libderiv->deriv_classes[6][8][7] = int_stack + 66020;
 Libderiv->deriv_classes[7][4][7] = int_stack + 67280;
 Libderiv->deriv_classes[7][5][7] = int_stack + 67820;
 Libderiv->deriv_classes[7][6][7] = int_stack + 68576;
 Libderiv->deriv_classes[7][7][7] = int_stack + 69584;
 Libderiv->deriv_classes[7][8][7] = int_stack + 70880;
 Libderiv->deriv_classes[4][4][6] = int_stack + 72500;
 Libderiv->deriv_classes[4][5][6] = int_stack + 72725;
 Libderiv->deriv_classes[4][6][6] = int_stack + 73040;
 Libderiv->deriv_classes[4][7][6] = int_stack + 73460;
 Libderiv->deriv_classes[4][8][6] = int_stack + 74000;
 Libderiv->deriv_classes[5][4][6] = int_stack + 74675;
 Libderiv->deriv_classes[5][5][6] = int_stack + 74990;
 Libderiv->deriv_classes[5][6][6] = int_stack + 75431;
 Libderiv->deriv_classes[5][7][6] = int_stack + 76019;
 Libderiv->deriv_classes[5][8][6] = int_stack + 76775;
 Libderiv->deriv_classes[6][4][6] = int_stack + 77720;
 Libderiv->deriv_classes[6][5][6] = int_stack + 78140;
 Libderiv->deriv_classes[6][6][6] = int_stack + 78728;
 Libderiv->deriv_classes[6][7][6] = int_stack + 79512;
 Libderiv->deriv_classes[6][8][6] = int_stack + 80520;
 Libderiv->dvrr_classes[7][4] = int_stack + 81780;
 Libderiv->deriv_classes[7][4][6] = int_stack + 82320;
 Libderiv->dvrr_classes[7][5] = int_stack + 82860;
 Libderiv->deriv_classes[7][5][6] = int_stack + 83616;
 Libderiv->dvrr_classes[7][6] = int_stack + 84372;
 Libderiv->deriv_classes[7][6][6] = int_stack + 85380;
 Libderiv->dvrr_classes[7][7] = int_stack + 86388;
 Libderiv->deriv_classes[7][7][6] = int_stack + 87684;
 Libderiv->deriv_classes[7][8][6] = int_stack + 88980;
 Libderiv->deriv_classes[4][4][2] = int_stack + 90600;
 Libderiv->deriv_classes[4][5][2] = int_stack + 90825;
 Libderiv->deriv_classes[4][6][2] = int_stack + 91140;
 Libderiv->deriv_classes[4][7][2] = int_stack + 91560;
 Libderiv->deriv_classes[4][8][2] = int_stack + 92100;
 Libderiv->deriv_classes[5][4][2] = int_stack + 92775;
 Libderiv->deriv_classes[5][5][2] = int_stack + 93090;
 Libderiv->deriv_classes[5][6][2] = int_stack + 93531;
 Libderiv->deriv_classes[5][7][2] = int_stack + 94119;
 Libderiv->deriv_classes[5][8][2] = int_stack + 94875;
 Libderiv->deriv_classes[6][4][2] = int_stack + 95820;
 Libderiv->deriv_classes[6][5][2] = int_stack + 96240;
 Libderiv->deriv_classes[6][6][2] = int_stack + 96828;
 Libderiv->deriv_classes[6][7][2] = int_stack + 97612;
 Libderiv->deriv_classes[6][8][2] = int_stack + 98620;
 Libderiv->deriv_classes[7][4][2] = int_stack + 99880;
 Libderiv->deriv_classes[7][5][2] = int_stack + 100420;
 Libderiv->deriv_classes[7][6][2] = int_stack + 101176;
 Libderiv->deriv_classes[7][7][2] = int_stack + 102184;
 Libderiv->deriv_classes[7][8][2] = int_stack + 103480;
 Libderiv->deriv_classes[4][4][1] = int_stack + 105100;
 Libderiv->deriv_classes[4][5][1] = int_stack + 105325;
 Libderiv->deriv_classes[4][6][1] = int_stack + 105640;
 Libderiv->deriv_classes[4][7][1] = int_stack + 106060;
 Libderiv->deriv_classes[4][8][1] = int_stack + 106600;
 Libderiv->deriv_classes[5][4][1] = int_stack + 107275;
 Libderiv->deriv_classes[5][5][1] = int_stack + 107590;
 Libderiv->deriv_classes[5][6][1] = int_stack + 108031;
 Libderiv->deriv_classes[5][7][1] = int_stack + 108619;
 Libderiv->deriv_classes[5][8][1] = int_stack + 109375;
 Libderiv->deriv_classes[6][4][1] = int_stack + 110320;
 Libderiv->deriv_classes[6][5][1] = int_stack + 110740;
 Libderiv->deriv_classes[6][6][1] = int_stack + 111328;
 Libderiv->deriv_classes[6][7][1] = int_stack + 112112;
 Libderiv->deriv_classes[6][8][1] = int_stack + 113120;
 Libderiv->deriv_classes[7][4][1] = int_stack + 114380;
 Libderiv->deriv_classes[7][5][1] = int_stack + 114920;
 Libderiv->deriv_classes[7][6][1] = int_stack + 115676;
 Libderiv->deriv_classes[7][7][1] = int_stack + 116684;
 Libderiv->deriv_classes[7][8][1] = int_stack + 117980;
 Libderiv->dvrr_classes[4][4] = int_stack + 119600;
 Libderiv->dvrr_classes[4][5] = int_stack + 119825;
 Libderiv->dvrr_classes[4][6] = int_stack + 120140;
 Libderiv->dvrr_classes[4][7] = int_stack + 120560;
 Libderiv->dvrr_classes[4][8] = int_stack + 121100;
 Libderiv->deriv_classes[4][4][0] = int_stack + 121775;
 Libderiv->deriv_classes[4][5][0] = int_stack + 122000;
 Libderiv->deriv_classes[4][6][0] = int_stack + 122315;
 Libderiv->deriv_classes[4][7][0] = int_stack + 122735;
 Libderiv->deriv_classes[4][8][0] = int_stack + 123275;
 Libderiv->dvrr_classes[5][4] = int_stack + 123950;
 Libderiv->dvrr_classes[5][5] = int_stack + 124265;
 Libderiv->dvrr_classes[5][6] = int_stack + 124706;
 Libderiv->dvrr_classes[5][7] = int_stack + 125294;
 Libderiv->dvrr_classes[5][8] = int_stack + 126050;
 Libderiv->deriv_classes[5][4][0] = int_stack + 126995;
 Libderiv->deriv_classes[5][5][0] = int_stack + 127310;
 Libderiv->deriv_classes[5][6][0] = int_stack + 127751;
 Libderiv->deriv_classes[5][7][0] = int_stack + 128339;
 Libderiv->deriv_classes[5][8][0] = int_stack + 129095;
 Libderiv->dvrr_classes[6][4] = int_stack + 130040;
 Libderiv->dvrr_classes[6][5] = int_stack + 130460;
 Libderiv->dvrr_classes[6][6] = int_stack + 131048;
 Libderiv->dvrr_classes[6][7] = int_stack + 131832;
 Libderiv->dvrr_classes[6][8] = int_stack + 132840;
 Libderiv->deriv_classes[6][4][0] = int_stack + 134100;
 Libderiv->deriv_classes[6][5][0] = int_stack + 134520;
 Libderiv->deriv_classes[6][6][0] = int_stack + 135108;
 Libderiv->deriv_classes[6][7][0] = int_stack + 135892;
 Libderiv->deriv_classes[6][8][0] = int_stack + 136900;
 Libderiv->deriv_classes[7][4][0] = int_stack + 138160;
 Libderiv->deriv_classes[7][5][0] = int_stack + 138700;
 Libderiv->deriv_classes[7][6][0] = int_stack + 139456;
 Libderiv->deriv_classes[7][7][0] = int_stack + 140464;
 Libderiv->deriv_classes[7][8][0] = int_stack + 141760;
 memset(int_stack,0,1147040);

 Libderiv->dvrr_stack = int_stack + 647770;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gfgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+143380,int_stack+119825,int_stack+119600,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+144055,int_stack+120140,int_stack+119825,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+145000,int_stack+144055,int_stack+143380,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+146350,int_stack+120560,int_stack+120140,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+147610,int_stack+146350,int_stack+144055,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+149500,int_stack+147610,int_stack+145000,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+151750,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119600,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+152425,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119825,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+153370,int_stack+152425,int_stack+151750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143380,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+154720,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120140,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+155980,int_stack+154720,int_stack+152425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144055,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+157870,int_stack+155980,int_stack+153370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145000,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+151750,int_stack+1500,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120560,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+160120,int_stack+151750,int_stack+154720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146350,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+151750,int_stack+160120,int_stack+155980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+147610,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+160120,int_stack+151750,int_stack+157870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149500,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+151750,int_stack+124265,int_stack+123950,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+152695,int_stack+124706,int_stack+124265,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+154018,int_stack+152695,int_stack+151750,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+155908,int_stack+125294,int_stack+124706,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+163495,int_stack+155908,int_stack+152695,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+166141,int_stack+163495,int_stack+154018,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157672,int_stack+2490,int_stack+2175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123950,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+158617,int_stack+2931,int_stack+2490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124265,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+158617,int_stack+157672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151750,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+169291,int_stack+3519,int_stack+2931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124706,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+171055,int_stack+169291,int_stack+158617, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152695,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+173701,int_stack+171055,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+154018,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+4275,int_stack+3519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125294,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+176851,int_stack+0,int_stack+169291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155908,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+176851,int_stack+171055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163495,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+176851,int_stack+0,int_stack+173701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+166141,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+181576,int_stack+176851,int_stack+160120,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+130460,int_stack+130040,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1260,int_stack+131048,int_stack+130460,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+169291,int_stack+1260,int_stack+0,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+171811,int_stack+131832,int_stack+131048,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+157672,int_stack+171811,int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+191701,int_stack+157672,int_stack+169291,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+5640,int_stack+5220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130040,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+174163,int_stack+6228,int_stack+5640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130460,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+195901,int_stack+174163,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3024,int_stack+7012,int_stack+6228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131048,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+198421,int_stack+3024,int_stack+174163, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+201949,int_stack+198421,int_stack+195901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+169291,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+206149,int_stack+8020,int_stack+7012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131832,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+209173,int_stack+206149,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171811,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3024,int_stack+209173,int_stack+198421, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+157672,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+206149,int_stack+3024,int_stack+201949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191701,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+212449,int_stack+206149,int_stack+176851,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+226624,int_stack+212449,int_stack+181576,225);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3024,int_stack+82860,int_stack+81780,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4644,int_stack+84372,int_stack+82860,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+195901,int_stack+4644,int_stack+3024,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+199141,int_stack+86388,int_stack+84372,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+174163,int_stack+199141,int_stack+4644,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+178699,int_stack+174163,int_stack+195901,36);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6912,int_stack+9820,int_stack+9280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81780,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202165,int_stack+10576,int_stack+9820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82860,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+184099,int_stack+202165,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3024,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+6912,int_stack+11584,int_stack+10576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84372,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+246874,int_stack+6912,int_stack+202165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4644,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+251410,int_stack+246874,int_stack+184099, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195901,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+184099,int_stack+12880,int_stack+11584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86388,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+256810,int_stack+184099,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+199141,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+6912,int_stack+256810,int_stack+246874, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+174163,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+256810,int_stack+6912,int_stack+251410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178699,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+264910,int_stack+256810,int_stack+206149,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+283810,int_stack+264910,int_stack+212449,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6912,int_stack+14725,int_stack+14500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119600, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7587,int_stack+15040,int_stack+14725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119825, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8532,int_stack+7587,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143380, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9882,int_stack+15460,int_stack+15040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120140, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+11142,int_stack+9882,int_stack+7587, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144055, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+13032,int_stack+11142,int_stack+8532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145000, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+6912,int_stack+16000,int_stack+15460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120560, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+246874,int_stack+6912,int_stack+9882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146350, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+6912,int_stack+246874,int_stack+11142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+147610, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+246874,int_stack+6912,int_stack+13032, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149500, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6912,int_stack+16990,int_stack+16675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123950, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7857,int_stack+17431,int_stack+16990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124265, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9180,int_stack+7857,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151750, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+11070,int_stack+18019,int_stack+17431, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124706, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+12834,int_stack+11070,int_stack+7857, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152695, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+250249,int_stack+12834,int_stack+9180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+154018, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+6912,int_stack+18775,int_stack+18019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125294, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+15480,int_stack+6912,int_stack+11070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155908, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+6912,int_stack+15480,int_stack+12834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163495, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+11322,int_stack+6912,int_stack+250249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+166141, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+250249,int_stack+11322,int_stack+246874,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+246874,int_stack+20140,int_stack+19720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130040, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+248134,int_stack+20728,int_stack+20140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130460, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+260374,int_stack+248134,int_stack+246874, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+262894,int_stack+21512,int_stack+20728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131048, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+265246,int_stack+262894,int_stack+248134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+268774,int_stack+265246,int_stack+260374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+169291, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+246874,int_stack+22520,int_stack+21512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131832, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+272974,int_stack+246874,int_stack+262894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171811, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+277678,int_stack+272974,int_stack+265246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+157672, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+260374,int_stack+277678,int_stack+268774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191701, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+266674,int_stack+260374,int_stack+11322,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+202165,int_stack+266674,int_stack+250249,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+246874,int_stack+24320,int_stack+23780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81780, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+248494,int_stack+25076,int_stack+24320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82860, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+250762,int_stack+248494,int_stack+246874, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3024, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+254002,int_stack+26084,int_stack+25076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84372, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+6912,int_stack+254002,int_stack+248494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4644, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+11448,int_stack+6912,int_stack+250762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195901, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+246874,int_stack+27380,int_stack+26084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86388, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+16848,int_stack+246874,int_stack+254002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+199141, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+246874,int_stack+16848,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+174163, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+16848,int_stack+246874,int_stack+11448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178699, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+312160,int_stack+16848,int_stack+260374,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+331060,int_stack+312160,int_stack+266674,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+29225,int_stack+29000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119600, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+312835,int_stack+29540,int_stack+29225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119825, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+313780,int_stack+312835,int_stack+312160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143380, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+315130,int_stack+29960,int_stack+29540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120140, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+316390,int_stack+315130,int_stack+312835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144055, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+318280,int_stack+316390,int_stack+313780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145000, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+30500,int_stack+29960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120560, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+320530,int_stack+312160,int_stack+315130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146350, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+320530,int_stack+316390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+147610, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+320530,int_stack+312160,int_stack+318280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149500, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+31490,int_stack+31175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123950, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+313105,int_stack+31931,int_stack+31490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124265, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+314428,int_stack+313105,int_stack+312160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151750, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+316318,int_stack+32519,int_stack+31931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124706, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+323905,int_stack+316318,int_stack+313105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152695, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+326551,int_stack+323905,int_stack+314428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+154018, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+33275,int_stack+32519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125294, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+246874,int_stack+312160,int_stack+316318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155908, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+246874,int_stack+323905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163495, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+246874,int_stack+312160,int_stack+326551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+166141, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+251599,int_stack+246874,int_stack+320530,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+34640,int_stack+34220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130040, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+313420,int_stack+35228,int_stack+34640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130460, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+315184,int_stack+313420,int_stack+312160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+317704,int_stack+36012,int_stack+35228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131048, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+320056,int_stack+317704,int_stack+313420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+323584,int_stack+320056,int_stack+315184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+169291, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+37020,int_stack+36012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131832, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+261724,int_stack+312160,int_stack+317704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171811, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+261724,int_stack+320056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+157672, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+261724,int_stack+312160,int_stack+323584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191701, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+312160,int_stack+261724,int_stack+246874,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+6912,int_stack+312160,int_stack+251599,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+246874,int_stack+38820,int_stack+38280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81780, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+248494,int_stack+39576,int_stack+38820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82860, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+250762,int_stack+248494,int_stack+246874, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3024, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+254002,int_stack+40584,int_stack+39576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84372, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+257026,int_stack+254002,int_stack+248494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4644, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+268024,int_stack+257026,int_stack+250762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195901, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+246874,int_stack+41880,int_stack+40584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86388, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+273424,int_stack+246874,int_stack+254002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+199141, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+246874,int_stack+273424,int_stack+257026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+174163, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+273424,int_stack+246874,int_stack+268024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178699, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+359410,int_stack+273424,int_stack+261724,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+246874,int_stack+359410,int_stack+312160,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+43725,int_stack+43500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+312835,int_stack+44040,int_stack+43725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+313780,int_stack+312835,int_stack+312160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+315130,int_stack+44460,int_stack+44040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+316390,int_stack+315130,int_stack+312835, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+318280,int_stack+316390,int_stack+313780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+45000,int_stack+44460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+320530,int_stack+312160,int_stack+315130, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+320530,int_stack+316390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+147610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+320530,int_stack+312160,int_stack+318280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+45990,int_stack+45675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+313105,int_stack+46431,int_stack+45990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+314428,int_stack+313105,int_stack+312160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+316318,int_stack+47019,int_stack+46431, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+323905,int_stack+316318,int_stack+313105, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+326551,int_stack+323905,int_stack+314428, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+154018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+47775,int_stack+47019, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+359410,int_stack+312160,int_stack+316318, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+359410,int_stack+323905, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+359410,int_stack+312160,int_stack+326551, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+166141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+364135,int_stack+359410,int_stack+320530,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+49140,int_stack+48720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+313420,int_stack+49728,int_stack+49140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+315184,int_stack+313420,int_stack+312160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+317704,int_stack+50512,int_stack+49728, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+320056,int_stack+317704,int_stack+313420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+323584,int_stack+320056,int_stack+315184, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+169291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+51520,int_stack+50512, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+374260,int_stack+312160,int_stack+317704, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+374260,int_stack+320056, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+157672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+275224,int_stack+312160,int_stack+323584, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+191701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+312160,int_stack+275224,int_stack+359410,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+27162,int_stack+312160,int_stack+364135,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+359410,int_stack+53320,int_stack+52780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+361030,int_stack+54076,int_stack+53320, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+363298,int_stack+361030,int_stack+359410, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+366538,int_stack+55084,int_stack+54076, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+369562,int_stack+366538,int_stack+361030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+47412,int_stack+369562,int_stack+363298, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+359410,int_stack+56380,int_stack+55084, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+184099,int_stack+359410,int_stack+366538, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+199141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+359410,int_stack+184099,int_stack+369562, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+174163, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+366970,int_stack+359410,int_stack+47412, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178699, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+375070,int_stack+366970,int_stack+275224,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+393970,int_stack+375070,int_stack+312160,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+58225,int_stack+58000, 0.0, zero_stack, 1.0, int_stack+119600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+312835,int_stack+58540,int_stack+58225, 0.0, zero_stack, 1.0, int_stack+119825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+313780,int_stack+312835,int_stack+312160, 0.0, zero_stack, 1.0, int_stack+143380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+315130,int_stack+58960,int_stack+58540, 0.0, zero_stack, 1.0, int_stack+120140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+316390,int_stack+315130,int_stack+312835, 0.0, zero_stack, 1.0, int_stack+144055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+318280,int_stack+316390,int_stack+313780, 0.0, zero_stack, 1.0, int_stack+145000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+59500,int_stack+58960, 0.0, zero_stack, 1.0, int_stack+120560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+320530,int_stack+312160,int_stack+315130, 0.0, zero_stack, 1.0, int_stack+146350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+320530,int_stack+316390, 0.0, zero_stack, 1.0, int_stack+147610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+320530,int_stack+312160,int_stack+318280, 0.0, zero_stack, 1.0, int_stack+149500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+60490,int_stack+60175, 0.0, zero_stack, 1.0, int_stack+123950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+313105,int_stack+60931,int_stack+60490, 0.0, zero_stack, 1.0, int_stack+124265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+314428,int_stack+313105,int_stack+312160, 0.0, zero_stack, 1.0, int_stack+151750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+316318,int_stack+61519,int_stack+60931, 0.0, zero_stack, 1.0, int_stack+124706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+323905,int_stack+316318,int_stack+313105, 0.0, zero_stack, 1.0, int_stack+152695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+326551,int_stack+323905,int_stack+314428, 0.0, zero_stack, 1.0, int_stack+154018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+62275,int_stack+61519, 0.0, zero_stack, 1.0, int_stack+125294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+275224,int_stack+312160,int_stack+316318, 0.0, zero_stack, 1.0, int_stack+155908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+275224,int_stack+323905, 0.0, zero_stack, 1.0, int_stack+163495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+275224,int_stack+312160,int_stack+326551, 0.0, zero_stack, 1.0, int_stack+166141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+47412,int_stack+275224,int_stack+320530,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+63640,int_stack+63220, 0.0, zero_stack, 1.0, int_stack+130040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+313420,int_stack+64228,int_stack+63640, 0.0, zero_stack, 1.0, int_stack+130460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+315184,int_stack+313420,int_stack+312160, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+317704,int_stack+65012,int_stack+64228, 0.0, zero_stack, 1.0, int_stack+131048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+320056,int_stack+317704,int_stack+313420, 0.0, zero_stack, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+323584,int_stack+320056,int_stack+315184, 0.0, zero_stack, 1.0, int_stack+169291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+66020,int_stack+65012, 0.0, zero_stack, 1.0, int_stack+131832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+57537,int_stack+312160,int_stack+317704, 0.0, zero_stack, 1.0, int_stack+171811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+57537,int_stack+320056, 0.0, zero_stack, 1.0, int_stack+157672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+57537,int_stack+312160,int_stack+323584, 0.0, zero_stack, 1.0, int_stack+191701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+312160,int_stack+57537,int_stack+275224,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+359410,int_stack+312160,int_stack+47412,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47412,int_stack+67820,int_stack+67280, 0.0, zero_stack, 1.0, int_stack+81780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49032,int_stack+68576,int_stack+67820, 0.0, zero_stack, 1.0, int_stack+82860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+51300,int_stack+49032,int_stack+47412, 0.0, zero_stack, 1.0, int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+275224,int_stack+69584,int_stack+68576, 0.0, zero_stack, 1.0, int_stack+84372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+278248,int_stack+275224,int_stack+49032, 0.0, zero_stack, 1.0, int_stack+4644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+63837,int_stack+278248,int_stack+51300, 0.0, zero_stack, 1.0, int_stack+195901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+47412,int_stack+70880,int_stack+69584, 0.0, zero_stack, 1.0, int_stack+86388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+51300,int_stack+47412,int_stack+275224, 0.0, zero_stack, 1.0, int_stack+199141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+379660,int_stack+51300,int_stack+278248, 0.0, zero_stack, 1.0, int_stack+174163, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+275224,int_stack+379660,int_stack+63837, 0.0, zero_stack, 1.0, int_stack+178699, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+422320,int_stack+275224,int_stack+57537,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+441220,int_stack+422320,int_stack+312160,225);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+72725,int_stack+72500, 1.0, int_stack+119600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+312835,int_stack+73040,int_stack+72725, 1.0, int_stack+119825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+313780,int_stack+312835,int_stack+312160, 1.0, int_stack+143380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+315130,int_stack+73460,int_stack+73040, 1.0, int_stack+120140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+316390,int_stack+315130,int_stack+312835, 1.0, int_stack+144055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+318280,int_stack+316390,int_stack+313780, 1.0, int_stack+145000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+143380,int_stack+74000,int_stack+73460, 1.0, int_stack+120560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+312160,int_stack+143380,int_stack+315130, 1.0, int_stack+146350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+320530,int_stack+312160,int_stack+316390, 1.0, int_stack+147610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+312160,int_stack+320530,int_stack+318280, 1.0, int_stack+149500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+315535,int_stack+74990,int_stack+74675, 1.0, int_stack+123950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+316480,int_stack+75431,int_stack+74990, 1.0, int_stack+124265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+317803,int_stack+316480,int_stack+315535, 1.0, int_stack+151750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+319693,int_stack+76019,int_stack+75431, 1.0, int_stack+124706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+321457,int_stack+319693,int_stack+316480, 1.0, int_stack+152695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+324103,int_stack+321457,int_stack+317803, 1.0, int_stack+154018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+151750,int_stack+76775,int_stack+76019, 1.0, int_stack+125294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+315535,int_stack+151750,int_stack+319693, 1.0, int_stack+155908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+422320,int_stack+315535,int_stack+321457, 1.0, int_stack+163495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+315535,int_stack+422320,int_stack+324103, 1.0, int_stack+166141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+422320,int_stack+315535,int_stack+312160,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+78140,int_stack+77720, 1.0, int_stack+130040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+313420,int_stack+78728,int_stack+78140, 1.0, int_stack+130460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+432445,int_stack+313420,int_stack+312160, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+434965,int_stack+79512,int_stack+78728, 1.0, int_stack+131048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+437317,int_stack+434965,int_stack+313420, 1.0, int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+320260,int_stack+437317,int_stack+432445, 1.0, int_stack+169291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+80520,int_stack+79512, 1.0, int_stack+131832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+324460,int_stack+0,int_stack+434965, 1.0, int_stack+171811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+275224,int_stack+324460,int_stack+437317, 1.0, int_stack+157672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+324460,int_stack+275224,int_stack+320260, 1.0, int_stack+191701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+379660,int_stack+324460,int_stack+315535,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+47412,int_stack+379660,int_stack+422320,225);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+422320,int_stack+83616,int_stack+82320, 1.0, int_stack+81780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+423940,int_stack+85380,int_stack+83616, 1.0, int_stack+82860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+426208,int_stack+423940,int_stack+422320, 1.0, int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+429448,int_stack+87684,int_stack+85380, 1.0, int_stack+84372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+432472,int_stack+429448,int_stack+423940, 1.0, int_stack+4644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+275224,int_stack+432472,int_stack+426208, 1.0, int_stack+195901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+422320,int_stack+88980,int_stack+87684, 1.0, int_stack+86388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+0,int_stack+422320,int_stack+429448, 1.0, int_stack+199141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+422320,int_stack+0,int_stack+432472, 1.0, int_stack+174163, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+429880,int_stack+422320,int_stack+275224, 1.0, int_stack+178699, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gg) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+67662,int_stack+429880,int_stack+324460,225);
 /*--- compute (hd|gg) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+469570,int_stack+67662,int_stack+379660,225);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+379660,int_stack+121100,int_stack+120560,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+169291,int_stack+379660,int_stack+146350,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+379660,int_stack+169291,int_stack+147610,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+382810,int_stack+379660,int_stack+149500,15);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+379660,int_stack+126050,int_stack+125294,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+386185,int_stack+379660,int_stack+155908,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+67662,int_stack+386185,int_stack+163495,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+386185,int_stack+67662,int_stack+166141,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+67662,int_stack+386185,int_stack+382810,225);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+77787,int_stack+132840,int_stack+131832,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+80811,int_stack+77787,int_stack+171811,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+275224,int_stack+80811,int_stack+157672,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+77787,int_stack+275224,int_stack+191701,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+422320,int_stack+77787,int_stack+386185,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+143380,int_stack+422320,int_stack+67662,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+275224,int_stack+90825,int_stack+90600,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+275899,int_stack+91140,int_stack+90825,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+276844,int_stack+275899,int_stack+275224,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+278194,int_stack+91560,int_stack+91140,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+279454,int_stack+278194,int_stack+275899,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+281344,int_stack+279454,int_stack+276844,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+275224,int_stack+92100,int_stack+91560,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+84087,int_stack+275224,int_stack+278194,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+379660,int_stack+84087,int_stack+279454,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+84087,int_stack+379660,int_stack+281344,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+379660,int_stack+93090,int_stack+92775,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+380605,int_stack+93531,int_stack+93090,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+87462,int_stack+380605,int_stack+379660,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+89352,int_stack+94119,int_stack+93531,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+91116,int_stack+89352,int_stack+380605,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+379660,int_stack+91116,int_stack+87462,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+275224,int_stack+94875,int_stack+94119,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+277492,int_stack+275224,int_stack+89352,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+436495,int_stack+277492,int_stack+91116,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+275224,int_stack+436495,int_stack+379660,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+312160,int_stack+275224,int_stack+84087, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+382810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+84087,int_stack+96240,int_stack+95820,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+85347,int_stack+96828,int_stack+96240,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+87111,int_stack+85347,int_stack+84087,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+89631,int_stack+97612,int_stack+96828,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+91983,int_stack+89631,int_stack+85347,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+436495,int_stack+91983,int_stack+87111,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+84087,int_stack+98620,int_stack+97612,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+0,int_stack+84087,int_stack+89631,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+84087,int_stack+0,int_stack+91983,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+0,int_stack+84087,int_stack+436495,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+84087,int_stack+0,int_stack+275224, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+386185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+163630,int_stack+84087,int_stack+312160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67662, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+312160,int_stack+100420,int_stack+99880,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+313780,int_stack+101176,int_stack+100420,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+316048,int_stack+313780,int_stack+312160,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+319288,int_stack+102184,int_stack+101176,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+322312,int_stack+319288,int_stack+313780,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+275224,int_stack+322312,int_stack+316048,36);
 /*--- compute (k0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+312160,int_stack+103480,int_stack+102184,36);
 /*--- compute (k0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+98262,int_stack+312160,int_stack+319288,36);
 /*--- compute (k0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+98262,int_stack+322312,36);
 /*--- compute (k0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+319720,int_stack+312160,int_stack+275224,36);
 /*--- compute (ip|gg) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+497920,int_stack+319720,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hd|gg) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+516820,int_stack+497920,int_stack+84087, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+422320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+84087,int_stack+105325,int_stack+105100,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+84762,int_stack+105640,int_stack+105325,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+85707,int_stack+84762,int_stack+84087,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+87057,int_stack+106060,int_stack+105640,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+88317,int_stack+87057,int_stack+84762,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+90207,int_stack+88317,int_stack+85707,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+84087,int_stack+106600,int_stack+106060,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+92457,int_stack+84087,int_stack+87057,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+379660,int_stack+92457,int_stack+88317,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+92457,int_stack+379660,int_stack+90207,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+379660,int_stack+107590,int_stack+107275,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+380605,int_stack+108031,int_stack+107590,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+95832,int_stack+380605,int_stack+379660,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+97722,int_stack+108619,int_stack+108031,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+99486,int_stack+97722,int_stack+380605,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+379660,int_stack+99486,int_stack+95832,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+102132,int_stack+109375,int_stack+108619,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+104400,int_stack+102132,int_stack+97722,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+84087,int_stack+104400,int_stack+99486,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+436495,int_stack+84087,int_stack+379660,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+95832,int_stack+436495,int_stack+92457, 0.0, zero_stack, 1.0, int_stack+382810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+379660,int_stack+110740,int_stack+110320,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+380920,int_stack+111328,int_stack+110740,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+84087,int_stack+380920,int_stack+379660,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+86607,int_stack+112112,int_stack+111328,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+88959,int_stack+86607,int_stack+380920,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+105957,int_stack+88959,int_stack+84087,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+379660,int_stack+113120,int_stack+112112,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+497920,int_stack+379660,int_stack+86607,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+502624,int_stack+497920,int_stack+88959,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+84087,int_stack+502624,int_stack+105957,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+497920,int_stack+84087,int_stack+436495, 0.0, zero_stack, 1.0, int_stack+386185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+545170,int_stack+497920,int_stack+95832, 0.0, zero_stack, 1.0, int_stack+67662, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+436495,int_stack+114920,int_stack+114380,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+438115,int_stack+115676,int_stack+114920,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+512095,int_stack+438115,int_stack+436495,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+379660,int_stack+116684,int_stack+115676,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+90387,int_stack+379660,int_stack+438115,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+94923,int_stack+90387,int_stack+512095,36);
 /*--- compute (k0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+512095,int_stack+117980,int_stack+116684,36);
 /*--- compute (k0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+100323,int_stack+512095,int_stack+379660,36);
 /*--- compute (k0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+106371,int_stack+100323,int_stack+90387,36);
 /*--- compute (k0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+275224,int_stack+106371,int_stack+94923,36);
 /*--- compute (ip|gg) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+312160,int_stack+275224,int_stack+84087, 0.0, zero_stack, 1.0, int_stack+77787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hd|gg) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+84087,int_stack+312160,int_stack+497920, 0.0, zero_stack, 1.0, int_stack+422320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+497920,int_stack+122000,int_stack+121775,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+498595,int_stack+122315,int_stack+122000,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+499540,int_stack+498595,int_stack+497920,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+500890,int_stack+122735,int_stack+122315,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+502150,int_stack+500890,int_stack+498595,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+504040,int_stack+502150,int_stack+499540,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+497920,int_stack+123275,int_stack+122735,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+506290,int_stack+497920,int_stack+500890,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+379660,int_stack+506290,int_stack+502150,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+506290,int_stack+379660,int_stack+504040,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+379660,int_stack+127310,int_stack+126995,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+380605,int_stack+127751,int_stack+127310,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+509665,int_stack+380605,int_stack+379660,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+511555,int_stack+128339,int_stack+127751,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+513319,int_stack+511555,int_stack+380605,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+379660,int_stack+513319,int_stack+509665,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+497920,int_stack+129095,int_stack+128339,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+500188,int_stack+497920,int_stack+511555,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+312160,int_stack+500188,int_stack+513319,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+436495,int_stack+312160,int_stack+379660,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+312160,int_stack+436495,int_stack+506290, 1.0, int_stack+382810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+379660,int_stack+134520,int_stack+134100,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+380920,int_stack+135108,int_stack+134520,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+382684,int_stack+380920,int_stack+379660,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+322285,int_stack+135892,int_stack+135108,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+324637,int_stack+322285,int_stack+380920,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+497920,int_stack+324637,int_stack+382684,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+379660,int_stack+136900,int_stack+135892,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+502120,int_stack+379660,int_stack+322285,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+379660,int_stack+502120,int_stack+324637,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+502120,int_stack+379660,int_stack+497920,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+112437,int_stack+502120,int_stack+436495, 1.0, int_stack+386185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+565420,int_stack+112437,int_stack+312160, 1.0, int_stack+67662, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+67662,int_stack+138700,int_stack+138160,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+69282,int_stack+139456,int_stack+138700,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+71550,int_stack+69282,int_stack+67662,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+312160,int_stack+140464,int_stack+139456,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+315184,int_stack+312160,int_stack+69282,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+319720,int_stack+315184,int_stack+71550,36);
 /*--- compute (k0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+67662,int_stack+141760,int_stack+140464,36);
 /*--- compute (k0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+71550,int_stack+67662,int_stack+312160,36);
 /*--- compute (k0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+379660,int_stack+71550,int_stack+315184,36);
 /*--- compute (k0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+67662,int_stack+379660,int_stack+319720,36);
 /*--- compute (ip|gg) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+312160,int_stack+67662,int_stack+502120, 1.0, int_stack+77787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (hd|gg) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+585670,int_stack+312160,int_stack+112437, 1.0, int_stack+422320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+614020,int_stack+283810,int_stack+226624,225);
     Libderiv->ABCD[11] = int_stack + 614020;
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+275224,int_stack+331060,int_stack+202165,225);
     Libderiv->ABCD[10] = int_stack + 275224;
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+308974,int_stack+246874,int_stack+6912,225);
     Libderiv->ABCD[9] = int_stack + 308974;
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+183880,int_stack+393970,int_stack+27162,225);
     Libderiv->ABCD[8] = int_stack + 183880;
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+379660,int_stack+441220,int_stack+359410,225);
     Libderiv->ABCD[7] = int_stack + 379660;
 /*--- compute (gf|gg) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+413410,int_stack+469570,int_stack+47412,225);
     Libderiv->ABCD[6] = int_stack + 413410;
 /*--- compute (gf|gg) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+447160,int_stack+516820,int_stack+163630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 447160;
 /*--- compute (gf|gg) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+480910,int_stack+84087,int_stack+545170, 0.0, zero_stack, 1.0, int_stack+143380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 480910;
 /*--- compute (gf|gg) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+514660,int_stack+585670,int_stack+565420, 1.0, int_stack+143380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 514660;

}
