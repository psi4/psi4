#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fdgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fd|gg) integrals */

void d1hrr_order_fdgg(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[5][4][11] = int_stack + 3625;
 Libderiv->deriv_classes[5][5][11] = int_stack + 3940;
 Libderiv->deriv_classes[5][6][11] = int_stack + 4381;
 Libderiv->deriv_classes[5][7][11] = int_stack + 4969;
 Libderiv->deriv_classes[5][8][11] = int_stack + 5725;
 Libderiv->deriv_classes[3][4][10] = int_stack + 6670;
 Libderiv->deriv_classes[3][5][10] = int_stack + 6820;
 Libderiv->deriv_classes[3][6][10] = int_stack + 7030;
 Libderiv->deriv_classes[3][7][10] = int_stack + 7310;
 Libderiv->deriv_classes[3][8][10] = int_stack + 7670;
 Libderiv->deriv_classes[4][4][10] = int_stack + 8120;
 Libderiv->deriv_classes[4][5][10] = int_stack + 8345;
 Libderiv->deriv_classes[4][6][10] = int_stack + 8660;
 Libderiv->deriv_classes[4][7][10] = int_stack + 9080;
 Libderiv->deriv_classes[4][8][10] = int_stack + 9620;
 Libderiv->deriv_classes[5][4][10] = int_stack + 10295;
 Libderiv->deriv_classes[5][5][10] = int_stack + 10610;
 Libderiv->deriv_classes[5][6][10] = int_stack + 11051;
 Libderiv->deriv_classes[5][7][10] = int_stack + 11639;
 Libderiv->deriv_classes[5][8][10] = int_stack + 12395;
 Libderiv->deriv_classes[3][4][9] = int_stack + 13340;
 Libderiv->deriv_classes[3][5][9] = int_stack + 13490;
 Libderiv->deriv_classes[3][6][9] = int_stack + 13700;
 Libderiv->deriv_classes[3][7][9] = int_stack + 13980;
 Libderiv->deriv_classes[3][8][9] = int_stack + 14340;
 Libderiv->deriv_classes[4][4][9] = int_stack + 14790;
 Libderiv->deriv_classes[4][5][9] = int_stack + 15015;
 Libderiv->deriv_classes[4][6][9] = int_stack + 15330;
 Libderiv->deriv_classes[4][7][9] = int_stack + 15750;
 Libderiv->deriv_classes[4][8][9] = int_stack + 16290;
 Libderiv->deriv_classes[5][4][9] = int_stack + 16965;
 Libderiv->deriv_classes[5][5][9] = int_stack + 17280;
 Libderiv->deriv_classes[5][6][9] = int_stack + 17721;
 Libderiv->deriv_classes[5][7][9] = int_stack + 18309;
 Libderiv->deriv_classes[5][8][9] = int_stack + 19065;
 Libderiv->deriv_classes[3][4][8] = int_stack + 20010;
 Libderiv->deriv_classes[3][5][8] = int_stack + 20160;
 Libderiv->deriv_classes[3][6][8] = int_stack + 20370;
 Libderiv->deriv_classes[3][7][8] = int_stack + 20650;
 Libderiv->deriv_classes[3][8][8] = int_stack + 21010;
 Libderiv->deriv_classes[4][4][8] = int_stack + 21460;
 Libderiv->deriv_classes[4][5][8] = int_stack + 21685;
 Libderiv->deriv_classes[4][6][8] = int_stack + 22000;
 Libderiv->deriv_classes[4][7][8] = int_stack + 22420;
 Libderiv->deriv_classes[4][8][8] = int_stack + 22960;
 Libderiv->deriv_classes[5][4][8] = int_stack + 23635;
 Libderiv->deriv_classes[5][5][8] = int_stack + 23950;
 Libderiv->deriv_classes[5][6][8] = int_stack + 24391;
 Libderiv->deriv_classes[5][7][8] = int_stack + 24979;
 Libderiv->deriv_classes[5][8][8] = int_stack + 25735;
 Libderiv->deriv_classes[3][4][7] = int_stack + 26680;
 Libderiv->deriv_classes[3][5][7] = int_stack + 26830;
 Libderiv->deriv_classes[3][6][7] = int_stack + 27040;
 Libderiv->deriv_classes[3][7][7] = int_stack + 27320;
 Libderiv->deriv_classes[3][8][7] = int_stack + 27680;
 Libderiv->deriv_classes[4][4][7] = int_stack + 28130;
 Libderiv->deriv_classes[4][5][7] = int_stack + 28355;
 Libderiv->deriv_classes[4][6][7] = int_stack + 28670;
 Libderiv->deriv_classes[4][7][7] = int_stack + 29090;
 Libderiv->deriv_classes[4][8][7] = int_stack + 29630;
 Libderiv->deriv_classes[5][4][7] = int_stack + 30305;
 Libderiv->deriv_classes[5][5][7] = int_stack + 30620;
 Libderiv->deriv_classes[5][6][7] = int_stack + 31061;
 Libderiv->deriv_classes[5][7][7] = int_stack + 31649;
 Libderiv->deriv_classes[5][8][7] = int_stack + 32405;
 Libderiv->deriv_classes[3][4][6] = int_stack + 33350;
 Libderiv->deriv_classes[3][5][6] = int_stack + 33500;
 Libderiv->deriv_classes[3][6][6] = int_stack + 33710;
 Libderiv->deriv_classes[3][7][6] = int_stack + 33990;
 Libderiv->deriv_classes[3][8][6] = int_stack + 34350;
 Libderiv->deriv_classes[4][4][6] = int_stack + 34800;
 Libderiv->deriv_classes[4][5][6] = int_stack + 35025;
 Libderiv->deriv_classes[4][6][6] = int_stack + 35340;
 Libderiv->deriv_classes[4][7][6] = int_stack + 35760;
 Libderiv->deriv_classes[4][8][6] = int_stack + 36300;
 Libderiv->dvrr_classes[5][4] = int_stack + 36975;
 Libderiv->deriv_classes[5][4][6] = int_stack + 37290;
 Libderiv->dvrr_classes[5][5] = int_stack + 37605;
 Libderiv->deriv_classes[5][5][6] = int_stack + 38046;
 Libderiv->dvrr_classes[5][6] = int_stack + 38487;
 Libderiv->deriv_classes[5][6][6] = int_stack + 39075;
 Libderiv->dvrr_classes[5][7] = int_stack + 39663;
 Libderiv->deriv_classes[5][7][6] = int_stack + 40419;
 Libderiv->deriv_classes[5][8][6] = int_stack + 41175;
 Libderiv->deriv_classes[3][4][2] = int_stack + 42120;
 Libderiv->deriv_classes[3][5][2] = int_stack + 42270;
 Libderiv->deriv_classes[3][6][2] = int_stack + 42480;
 Libderiv->deriv_classes[3][7][2] = int_stack + 42760;
 Libderiv->deriv_classes[3][8][2] = int_stack + 43120;
 Libderiv->deriv_classes[4][4][2] = int_stack + 43570;
 Libderiv->deriv_classes[4][5][2] = int_stack + 43795;
 Libderiv->deriv_classes[4][6][2] = int_stack + 44110;
 Libderiv->deriv_classes[4][7][2] = int_stack + 44530;
 Libderiv->deriv_classes[4][8][2] = int_stack + 45070;
 Libderiv->deriv_classes[5][4][2] = int_stack + 45745;
 Libderiv->deriv_classes[5][5][2] = int_stack + 46060;
 Libderiv->deriv_classes[5][6][2] = int_stack + 46501;
 Libderiv->deriv_classes[5][7][2] = int_stack + 47089;
 Libderiv->deriv_classes[5][8][2] = int_stack + 47845;
 Libderiv->deriv_classes[3][4][1] = int_stack + 48790;
 Libderiv->deriv_classes[3][5][1] = int_stack + 48940;
 Libderiv->deriv_classes[3][6][1] = int_stack + 49150;
 Libderiv->deriv_classes[3][7][1] = int_stack + 49430;
 Libderiv->deriv_classes[3][8][1] = int_stack + 49790;
 Libderiv->deriv_classes[4][4][1] = int_stack + 50240;
 Libderiv->deriv_classes[4][5][1] = int_stack + 50465;
 Libderiv->deriv_classes[4][6][1] = int_stack + 50780;
 Libderiv->deriv_classes[4][7][1] = int_stack + 51200;
 Libderiv->deriv_classes[4][8][1] = int_stack + 51740;
 Libderiv->deriv_classes[5][4][1] = int_stack + 52415;
 Libderiv->deriv_classes[5][5][1] = int_stack + 52730;
 Libderiv->deriv_classes[5][6][1] = int_stack + 53171;
 Libderiv->deriv_classes[5][7][1] = int_stack + 53759;
 Libderiv->deriv_classes[5][8][1] = int_stack + 54515;
 Libderiv->dvrr_classes[3][4] = int_stack + 55460;
 Libderiv->dvrr_classes[3][5] = int_stack + 55610;
 Libderiv->dvrr_classes[3][6] = int_stack + 55820;
 Libderiv->dvrr_classes[3][7] = int_stack + 56100;
 Libderiv->dvrr_classes[3][8] = int_stack + 56460;
 Libderiv->deriv_classes[3][4][0] = int_stack + 56910;
 Libderiv->deriv_classes[3][5][0] = int_stack + 57060;
 Libderiv->deriv_classes[3][6][0] = int_stack + 57270;
 Libderiv->deriv_classes[3][7][0] = int_stack + 57550;
 Libderiv->deriv_classes[3][8][0] = int_stack + 57910;
 Libderiv->dvrr_classes[4][4] = int_stack + 58360;
 Libderiv->dvrr_classes[4][5] = int_stack + 58585;
 Libderiv->dvrr_classes[4][6] = int_stack + 58900;
 Libderiv->dvrr_classes[4][7] = int_stack + 59320;
 Libderiv->dvrr_classes[4][8] = int_stack + 59860;
 Libderiv->deriv_classes[4][4][0] = int_stack + 60535;
 Libderiv->deriv_classes[4][5][0] = int_stack + 60760;
 Libderiv->deriv_classes[4][6][0] = int_stack + 61075;
 Libderiv->deriv_classes[4][7][0] = int_stack + 61495;
 Libderiv->deriv_classes[4][8][0] = int_stack + 62035;
 Libderiv->deriv_classes[5][4][0] = int_stack + 62710;
 Libderiv->deriv_classes[5][5][0] = int_stack + 63025;
 Libderiv->deriv_classes[5][6][0] = int_stack + 63466;
 Libderiv->deriv_classes[5][7][0] = int_stack + 64054;
 Libderiv->deriv_classes[5][8][0] = int_stack + 64810;
 memset(int_stack,0,526040);

 Libderiv->dvrr_stack = int_stack + 222751;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fdgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+65755,int_stack+55610,int_stack+55460,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+66205,int_stack+55820,int_stack+55610,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+66835,int_stack+66205,int_stack+65755,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+67735,int_stack+56100,int_stack+55820,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+68575,int_stack+67735,int_stack+66205,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+69835,int_stack+68575,int_stack+66835,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+71335,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55460,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+71785,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55610,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72415,int_stack+71785,int_stack+71335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65755,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+73315,int_stack+640,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55820,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+74155,int_stack+73315,int_stack+71785, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66205,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+75415,int_stack+74155,int_stack+72415, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66835,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+71335,int_stack+1000,int_stack+640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56100,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+76915,int_stack+71335,int_stack+73315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67735,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+71335,int_stack+76915,int_stack+74155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68575,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+76915,int_stack+71335,int_stack+75415, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69835,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+71335,int_stack+58585,int_stack+58360,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+72010,int_stack+58900,int_stack+58585,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+72955,int_stack+72010,int_stack+71335,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+74305,int_stack+59320,int_stack+58900,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+79165,int_stack+74305,int_stack+72010,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+81055,int_stack+79165,int_stack+72955,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+75565,int_stack+1675,int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58360,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+1990,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58585,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83305,int_stack+0,int_stack+75565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71335,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+75565,int_stack+2410,int_stack+1990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58900,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+84655,int_stack+75565,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72010,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+84655,int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72955,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+86545,int_stack+2950,int_stack+2410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59320,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+88165,int_stack+86545,int_stack+75565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74305,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+90685,int_stack+88165,int_stack+84655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79165,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+83305,int_stack+90685,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81055,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+86680,int_stack+83305,int_stack+76915,225);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+37605,int_stack+36975,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+945,int_stack+38487,int_stack+37605,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+75565,int_stack+945,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+93430,int_stack+39663,int_stack+38487,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+95194,int_stack+93430,int_stack+945,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+97840,int_stack+95194,int_stack+75565,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2268,int_stack+3940,int_stack+3625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36975,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+77455,int_stack+4381,int_stack+3940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37605,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+100990,int_stack+77455,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2268,int_stack+4969,int_stack+4381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38487,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+102880,int_stack+2268,int_stack+77455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+105526,int_stack+102880,int_stack+100990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75565,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+108676,int_stack+5725,int_stack+4969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39663,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+110944,int_stack+108676,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93430,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+114472,int_stack+110944,int_stack+102880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95194,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+108676,int_stack+114472,int_stack+105526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97840,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+113401,int_stack+108676,int_stack+83305,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83305,int_stack+6820,int_stack+6670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55460, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83755,int_stack+7030,int_stack+6820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55610, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84385,int_stack+83755,int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65755, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+85285,int_stack+7310,int_stack+7030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55820, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2268,int_stack+85285,int_stack+83755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66205, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3528,int_stack+2268,int_stack+84385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66835, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+83305,int_stack+7670,int_stack+7310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56100, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+5028,int_stack+83305,int_stack+85285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67735, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+83305,int_stack+5028,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68575, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+5028,int_stack+83305,int_stack+3528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69835, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83305,int_stack+8345,int_stack+8120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58360, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83980,int_stack+8660,int_stack+8345, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58585, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84925,int_stack+83980,int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71335, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2268,int_stack+9080,int_stack+8660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58900, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+100990,int_stack+2268,int_stack+83980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72010, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+102880,int_stack+100990,int_stack+84925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72955, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+83305,int_stack+9620,int_stack+9080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59320, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+7278,int_stack+83305,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74305, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+83305,int_stack+7278,int_stack+100990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79165, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+105130,int_stack+83305,int_stack+102880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81055, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+123526,int_stack+105130,int_stack+5028,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83305,int_stack+10610,int_stack+10295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36975, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84250,int_stack+11051,int_stack+10610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37605, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+100990,int_stack+84250,int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+102880,int_stack+11639,int_stack+11051, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38487, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2268,int_stack+102880,int_stack+84250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+83305,int_stack+2268,int_stack+100990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75565, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+4914,int_stack+12395,int_stack+11639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39663, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+7182,int_stack+4914,int_stack+102880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93430, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+108505,int_stack+7182,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95194, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+2268,int_stack+108505,int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97840, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+130276,int_stack+2268,int_stack+105130,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2268,int_stack+13490,int_stack+13340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55460, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2718,int_stack+13700,int_stack+13490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55610, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3348,int_stack+2718,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65755, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4248,int_stack+13980,int_stack+13700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55820, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+5088,int_stack+4248,int_stack+2718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66205, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6348,int_stack+5088,int_stack+3348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66835, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+2268,int_stack+14340,int_stack+13980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56100, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+7848,int_stack+2268,int_stack+4248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67735, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+2268,int_stack+7848,int_stack+5088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68575, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+7848,int_stack+2268,int_stack+6348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69835, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2268,int_stack+15015,int_stack+14790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58360, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2943,int_stack+15330,int_stack+15015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58585, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3888,int_stack+2943,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71335, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+5238,int_stack+15750,int_stack+15330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58900, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10098,int_stack+5238,int_stack+2943, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72010, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+11988,int_stack+10098,int_stack+3888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72955, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+2268,int_stack+16290,int_stack+15750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59320, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+14238,int_stack+2268,int_stack+5238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74305, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+2268,int_stack+14238,int_stack+10098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79165, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+83305,int_stack+2268,int_stack+11988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81055, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+10098,int_stack+83305,int_stack+7848,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2268,int_stack+17280,int_stack+16965, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36975, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3213,int_stack+17721,int_stack+17280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37605, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4536,int_stack+3213,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+6426,int_stack+18309,int_stack+17721, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38487, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+100990,int_stack+6426,int_stack+3213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+103636,int_stack+100990,int_stack+4536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75565, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+2268,int_stack+19065,int_stack+18309, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39663, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+106786,int_stack+2268,int_stack+6426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93430, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+2268,int_stack+106786,int_stack+100990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95194, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+106786,int_stack+2268,int_stack+103636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97840, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+140401,int_stack+106786,int_stack+83305,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83305,int_stack+20160,int_stack+20010, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83755,int_stack+20370,int_stack+20160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84385,int_stack+83755,int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+85285,int_stack+20650,int_stack+20370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2268,int_stack+85285,int_stack+83755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3528,int_stack+2268,int_stack+84385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+83305,int_stack+21010,int_stack+20650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+5028,int_stack+83305,int_stack+85285, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+83305,int_stack+5028,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+5028,int_stack+83305,int_stack+3528, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83305,int_stack+21685,int_stack+21460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83980,int_stack+22000,int_stack+21685, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84925,int_stack+83980,int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2268,int_stack+22420,int_stack+22000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7278,int_stack+2268,int_stack+83980, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+100990,int_stack+7278,int_stack+84925, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+83305,int_stack+22960,int_stack+22420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+103240,int_stack+83305,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+83305,int_stack+103240,int_stack+7278, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+103240,int_stack+83305,int_stack+100990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+106615,int_stack+103240,int_stack+5028,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+100990,int_stack+23950,int_stack+23635, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83305,int_stack+24391,int_stack+23950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84628,int_stack+83305,int_stack+100990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+100990,int_stack+24979,int_stack+24391, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38487, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2268,int_stack+100990,int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4914,int_stack+2268,int_stack+84628, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+83305,int_stack+25735,int_stack+24979, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39663, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+16848,int_stack+83305,int_stack+100990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+20376,int_stack+16848,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+150526,int_stack+20376,int_stack+4914, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+155251,int_stack+150526,int_stack+103240,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+150526,int_stack+26830,int_stack+26680, 0.0, zero_stack, 1.0, int_stack+55460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+150976,int_stack+27040,int_stack+26830, 0.0, zero_stack, 1.0, int_stack+55610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+151606,int_stack+150976,int_stack+150526, 0.0, zero_stack, 1.0, int_stack+65755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+152506,int_stack+27320,int_stack+27040, 0.0, zero_stack, 1.0, int_stack+55820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+153346,int_stack+152506,int_stack+150976, 0.0, zero_stack, 1.0, int_stack+66205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2268,int_stack+153346,int_stack+151606, 0.0, zero_stack, 1.0, int_stack+66835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+150526,int_stack+27680,int_stack+27320, 0.0, zero_stack, 1.0, int_stack+56100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3768,int_stack+150526,int_stack+152506, 0.0, zero_stack, 1.0, int_stack+67735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+150526,int_stack+3768,int_stack+153346, 0.0, zero_stack, 1.0, int_stack+68575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+3768,int_stack+150526,int_stack+2268, 0.0, zero_stack, 1.0, int_stack+69835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2268,int_stack+28355,int_stack+28130, 0.0, zero_stack, 1.0, int_stack+58360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+150526,int_stack+28670,int_stack+28355, 0.0, zero_stack, 1.0, int_stack+58585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+151471,int_stack+150526,int_stack+2268, 0.0, zero_stack, 1.0, int_stack+71335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2268,int_stack+29090,int_stack+28670, 0.0, zero_stack, 1.0, int_stack+58900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+152821,int_stack+2268,int_stack+150526, 0.0, zero_stack, 1.0, int_stack+72010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6018,int_stack+152821,int_stack+151471, 0.0, zero_stack, 1.0, int_stack+72955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+150526,int_stack+29630,int_stack+29090, 0.0, zero_stack, 1.0, int_stack+59320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+16848,int_stack+150526,int_stack+2268, 0.0, zero_stack, 1.0, int_stack+74305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+19368,int_stack+16848,int_stack+152821, 0.0, zero_stack, 1.0, int_stack+79165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+83305,int_stack+19368,int_stack+6018, 0.0, zero_stack, 1.0, int_stack+81055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+16848,int_stack+83305,int_stack+3768,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23598,int_stack+30620,int_stack+30305, 0.0, zero_stack, 1.0, int_stack+36975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24543,int_stack+31061,int_stack+30620, 0.0, zero_stack, 1.0, int_stack+37605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+25866,int_stack+24543,int_stack+23598, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+27756,int_stack+31649,int_stack+31061, 0.0, zero_stack, 1.0, int_stack+38487, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2268,int_stack+27756,int_stack+24543, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4914,int_stack+2268,int_stack+25866, 0.0, zero_stack, 1.0, int_stack+75565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+23598,int_stack+32405,int_stack+31649, 0.0, zero_stack, 1.0, int_stack+39663, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+29520,int_stack+23598,int_stack+27756, 0.0, zero_stack, 1.0, int_stack+93430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+23598,int_stack+29520,int_stack+2268, 0.0, zero_stack, 1.0, int_stack+95194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+150526,int_stack+23598,int_stack+4914, 0.0, zero_stack, 1.0, int_stack+97840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+165376,int_stack+150526,int_stack+83305,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83305,int_stack+33500,int_stack+33350, 1.0, int_stack+55460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83755,int_stack+33710,int_stack+33500, 1.0, int_stack+55610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84385,int_stack+83755,int_stack+83305, 1.0, int_stack+65755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+85285,int_stack+33990,int_stack+33710, 1.0, int_stack+55820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+150526,int_stack+85285,int_stack+83755, 1.0, int_stack+66205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+151786,int_stack+150526,int_stack+84385, 1.0, int_stack+66835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+65755,int_stack+34350,int_stack+33990, 1.0, int_stack+56100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+83305,int_stack+65755,int_stack+85285, 1.0, int_stack+67735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+23598,int_stack+83305,int_stack+150526, 1.0, int_stack+68575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+83305,int_stack+23598,int_stack+151786, 1.0, int_stack+69835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23598,int_stack+35025,int_stack+34800, 1.0, int_stack+58360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24273,int_stack+35340,int_stack+35025, 1.0, int_stack+58585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+25218,int_stack+24273,int_stack+23598, 1.0, int_stack+71335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+26568,int_stack+35760,int_stack+35340, 1.0, int_stack+58900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+27828,int_stack+26568,int_stack+24273, 1.0, int_stack+72010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+29718,int_stack+27828,int_stack+25218, 1.0, int_stack+72955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+71335,int_stack+36300,int_stack+35760, 1.0, int_stack+59320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+23598,int_stack+71335,int_stack+26568, 1.0, int_stack+74305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+31968,int_stack+23598,int_stack+27828, 1.0, int_stack+79165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+23598,int_stack+31968,int_stack+29718, 1.0, int_stack+81055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+26973,int_stack+23598,int_stack+83305,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83305,int_stack+38046,int_stack+37290, 1.0, int_stack+36975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84250,int_stack+39075,int_stack+38046, 1.0, int_stack+37605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33723,int_stack+84250,int_stack+83305, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+35613,int_stack+40419,int_stack+39075, 1.0, int_stack+38487, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+71335,int_stack+35613,int_stack+84250, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+71335,int_stack+33723, 1.0, int_stack+75565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+75565,int_stack+41175,int_stack+40419, 1.0, int_stack+39663, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3150,int_stack+75565,int_stack+35613, 1.0, int_stack+93430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+33723,int_stack+3150,int_stack+71335, 1.0, int_stack+95194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+150526,int_stack+33723,int_stack+0, 1.0, int_stack+97840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+93430,int_stack+150526,int_stack+23598,225);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+23598,int_stack+56460,int_stack+56100,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+24678,int_stack+23598,int_stack+67735,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+150526,int_stack+24678,int_stack+68575,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+23598,int_stack+150526,int_stack+69835,10);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+150526,int_stack+59860,int_stack+59320,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+152146,int_stack+150526,int_stack+74305,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+152146,int_stack+79165,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+83305,int_stack+0,int_stack+81055,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+83305,int_stack+23598,225);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6750,int_stack+42270,int_stack+42120,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7200,int_stack+42480,int_stack+42270,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+7830,int_stack+7200,int_stack+6750,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8730,int_stack+42760,int_stack+42480,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+150526,int_stack+8730,int_stack+7200,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+151786,int_stack+150526,int_stack+7830,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+6750,int_stack+43120,int_stack+42760,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+153286,int_stack+6750,int_stack+8730,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+6750,int_stack+153286,int_stack+150526,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+33723,int_stack+6750,int_stack+151786,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6750,int_stack+43795,int_stack+43570,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7425,int_stack+44110,int_stack+43795,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8370,int_stack+7425,int_stack+6750,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+150526,int_stack+44530,int_stack+44110,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+151786,int_stack+150526,int_stack+7425,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+35973,int_stack+151786,int_stack+8370,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+6750,int_stack+45070,int_stack+44530,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+38223,int_stack+6750,int_stack+150526,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+6750,int_stack+38223,int_stack+151786,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+38223,int_stack+6750,int_stack+35973,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+65755,int_stack+38223,int_stack+33723, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+33723,int_stack+46060,int_stack+45745,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+34668,int_stack+46501,int_stack+46060,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+35991,int_stack+34668,int_stack+33723,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+6750,int_stack+47089,int_stack+46501,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+41598,int_stack+6750,int_stack+34668,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+150526,int_stack+41598,int_stack+35991,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+33723,int_stack+47845,int_stack+47089,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+44244,int_stack+33723,int_stack+6750,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+33723,int_stack+44244,int_stack+41598,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+41598,int_stack+33723,int_stack+150526,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+72505,int_stack+41598,int_stack+38223, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+150526,int_stack+48940,int_stack+48790,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+150976,int_stack+49150,int_stack+48940,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+151606,int_stack+150976,int_stack+150526,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+152506,int_stack+49430,int_stack+49150,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+153346,int_stack+152506,int_stack+150976,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+33723,int_stack+153346,int_stack+151606,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+150526,int_stack+49790,int_stack+49430,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+35223,int_stack+150526,int_stack+152506,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+150526,int_stack+35223,int_stack+153346,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+35223,int_stack+150526,int_stack+33723,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+82630,int_stack+50465,int_stack+50240,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+33723,int_stack+50780,int_stack+50465,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+150526,int_stack+33723,int_stack+82630,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+151876,int_stack+51200,int_stack+50780,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+153136,int_stack+151876,int_stack+33723,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+37473,int_stack+153136,int_stack+150526,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+39723,int_stack+51740,int_stack+51200,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+41343,int_stack+39723,int_stack+151876,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+43863,int_stack+41343,int_stack+153136,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+39723,int_stack+43863,int_stack+37473,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+43098,int_stack+39723,int_stack+35223, 0.0, zero_stack, 1.0, int_stack+23598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+49848,int_stack+52730,int_stack+52415,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+50793,int_stack+53171,int_stack+52730,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+150526,int_stack+50793,int_stack+49848,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+152416,int_stack+53759,int_stack+53171,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+33723,int_stack+152416,int_stack+50793,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+49848,int_stack+33723,int_stack+150526,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+36369,int_stack+54515,int_stack+53759,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+52998,int_stack+36369,int_stack+152416,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+150526,int_stack+52998,int_stack+33723,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+33723,int_stack+150526,int_stack+49848,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+175501,int_stack+33723,int_stack+39723, 0.0, zero_stack, 1.0, int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+33723,int_stack+57060,int_stack+56910,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+34173,int_stack+57270,int_stack+57060,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+34803,int_stack+34173,int_stack+33723,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+35703,int_stack+57550,int_stack+57270,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+36543,int_stack+35703,int_stack+34173,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+37803,int_stack+36543,int_stack+34803,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+33723,int_stack+57910,int_stack+57550,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+39303,int_stack+33723,int_stack+35703,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+33723,int_stack+39303,int_stack+36543,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+39303,int_stack+33723,int_stack+37803,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+82630,int_stack+60760,int_stack+60535,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+33723,int_stack+61075,int_stack+60760,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+34668,int_stack+33723,int_stack+82630,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+36018,int_stack+61495,int_stack+61075,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+37278,int_stack+36018,int_stack+33723,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+49848,int_stack+37278,int_stack+34668,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+33723,int_stack+62035,int_stack+61495,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+52098,int_stack+33723,int_stack+36018,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+33723,int_stack+52098,int_stack+37278,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+52098,int_stack+33723,int_stack+49848,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+55473,int_stack+52098,int_stack+39303, 1.0, int_stack+23598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23598,int_stack+63025,int_stack+62710,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24543,int_stack+63466,int_stack+63025,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+49848,int_stack+24543,int_stack+23598,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+33723,int_stack+64054,int_stack+63466,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+35487,int_stack+33723,int_stack+24543,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+23598,int_stack+35487,int_stack+49848,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+38133,int_stack+64810,int_stack+64054,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+62223,int_stack+38133,int_stack+33723,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+38133,int_stack+62223,int_stack+35487,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+150526,int_stack+38133,int_stack+23598,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+185626,int_stack+150526,int_stack+52098, 1.0, int_stack+83305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+195751,int_stack+113401,int_stack+86680,225);
     Libderiv->ABCD[11] = int_stack + 195751;
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+209251,int_stack+130276,int_stack+123526,225);
     Libderiv->ABCD[10] = int_stack + 209251;
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+113365,int_stack+140401,int_stack+10098,225);
     Libderiv->ABCD[9] = int_stack + 113365;
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+126865,int_stack+155251,int_stack+106615,225);
     Libderiv->ABCD[8] = int_stack + 126865;
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+140365,int_stack+165376,int_stack+16848,225);
     Libderiv->ABCD[7] = int_stack + 140365;
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+6750,int_stack+93430,int_stack+26973,225);
     Libderiv->ABCD[6] = int_stack + 6750;
 /*--- compute (fd|gg) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+82630,int_stack+72505,int_stack+65755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 82630;
 /*--- compute (fd|gg) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+62223,int_stack+175501,int_stack+43098, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 62223;
 /*--- compute (fd|gg) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+96130,int_stack+185626,int_stack+55473, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 96130;

}
