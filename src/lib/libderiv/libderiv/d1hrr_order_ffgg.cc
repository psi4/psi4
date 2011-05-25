#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ffgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (ff|gg) integrals */

void d1hrr_order_ffgg(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[6][4][11] = int_stack + 6670;
 Libderiv->deriv_classes[6][5][11] = int_stack + 7090;
 Libderiv->deriv_classes[6][6][11] = int_stack + 7678;
 Libderiv->deriv_classes[6][7][11] = int_stack + 8462;
 Libderiv->deriv_classes[6][8][11] = int_stack + 9470;
 Libderiv->deriv_classes[3][4][10] = int_stack + 10730;
 Libderiv->deriv_classes[3][5][10] = int_stack + 10880;
 Libderiv->deriv_classes[3][6][10] = int_stack + 11090;
 Libderiv->deriv_classes[3][7][10] = int_stack + 11370;
 Libderiv->deriv_classes[3][8][10] = int_stack + 11730;
 Libderiv->deriv_classes[4][4][10] = int_stack + 12180;
 Libderiv->deriv_classes[4][5][10] = int_stack + 12405;
 Libderiv->deriv_classes[4][6][10] = int_stack + 12720;
 Libderiv->deriv_classes[4][7][10] = int_stack + 13140;
 Libderiv->deriv_classes[4][8][10] = int_stack + 13680;
 Libderiv->deriv_classes[5][4][10] = int_stack + 14355;
 Libderiv->deriv_classes[5][5][10] = int_stack + 14670;
 Libderiv->deriv_classes[5][6][10] = int_stack + 15111;
 Libderiv->deriv_classes[5][7][10] = int_stack + 15699;
 Libderiv->deriv_classes[5][8][10] = int_stack + 16455;
 Libderiv->deriv_classes[6][4][10] = int_stack + 17400;
 Libderiv->deriv_classes[6][5][10] = int_stack + 17820;
 Libderiv->deriv_classes[6][6][10] = int_stack + 18408;
 Libderiv->deriv_classes[6][7][10] = int_stack + 19192;
 Libderiv->deriv_classes[6][8][10] = int_stack + 20200;
 Libderiv->deriv_classes[3][4][9] = int_stack + 21460;
 Libderiv->deriv_classes[3][5][9] = int_stack + 21610;
 Libderiv->deriv_classes[3][6][9] = int_stack + 21820;
 Libderiv->deriv_classes[3][7][9] = int_stack + 22100;
 Libderiv->deriv_classes[3][8][9] = int_stack + 22460;
 Libderiv->deriv_classes[4][4][9] = int_stack + 22910;
 Libderiv->deriv_classes[4][5][9] = int_stack + 23135;
 Libderiv->deriv_classes[4][6][9] = int_stack + 23450;
 Libderiv->deriv_classes[4][7][9] = int_stack + 23870;
 Libderiv->deriv_classes[4][8][9] = int_stack + 24410;
 Libderiv->deriv_classes[5][4][9] = int_stack + 25085;
 Libderiv->deriv_classes[5][5][9] = int_stack + 25400;
 Libderiv->deriv_classes[5][6][9] = int_stack + 25841;
 Libderiv->deriv_classes[5][7][9] = int_stack + 26429;
 Libderiv->deriv_classes[5][8][9] = int_stack + 27185;
 Libderiv->deriv_classes[6][4][9] = int_stack + 28130;
 Libderiv->deriv_classes[6][5][9] = int_stack + 28550;
 Libderiv->deriv_classes[6][6][9] = int_stack + 29138;
 Libderiv->deriv_classes[6][7][9] = int_stack + 29922;
 Libderiv->deriv_classes[6][8][9] = int_stack + 30930;
 Libderiv->deriv_classes[3][4][8] = int_stack + 32190;
 Libderiv->deriv_classes[3][5][8] = int_stack + 32340;
 Libderiv->deriv_classes[3][6][8] = int_stack + 32550;
 Libderiv->deriv_classes[3][7][8] = int_stack + 32830;
 Libderiv->deriv_classes[3][8][8] = int_stack + 33190;
 Libderiv->deriv_classes[4][4][8] = int_stack + 33640;
 Libderiv->deriv_classes[4][5][8] = int_stack + 33865;
 Libderiv->deriv_classes[4][6][8] = int_stack + 34180;
 Libderiv->deriv_classes[4][7][8] = int_stack + 34600;
 Libderiv->deriv_classes[4][8][8] = int_stack + 35140;
 Libderiv->deriv_classes[5][4][8] = int_stack + 35815;
 Libderiv->deriv_classes[5][5][8] = int_stack + 36130;
 Libderiv->deriv_classes[5][6][8] = int_stack + 36571;
 Libderiv->deriv_classes[5][7][8] = int_stack + 37159;
 Libderiv->deriv_classes[5][8][8] = int_stack + 37915;
 Libderiv->deriv_classes[6][4][8] = int_stack + 38860;
 Libderiv->deriv_classes[6][5][8] = int_stack + 39280;
 Libderiv->deriv_classes[6][6][8] = int_stack + 39868;
 Libderiv->deriv_classes[6][7][8] = int_stack + 40652;
 Libderiv->deriv_classes[6][8][8] = int_stack + 41660;
 Libderiv->deriv_classes[3][4][7] = int_stack + 42920;
 Libderiv->deriv_classes[3][5][7] = int_stack + 43070;
 Libderiv->deriv_classes[3][6][7] = int_stack + 43280;
 Libderiv->deriv_classes[3][7][7] = int_stack + 43560;
 Libderiv->deriv_classes[3][8][7] = int_stack + 43920;
 Libderiv->deriv_classes[4][4][7] = int_stack + 44370;
 Libderiv->deriv_classes[4][5][7] = int_stack + 44595;
 Libderiv->deriv_classes[4][6][7] = int_stack + 44910;
 Libderiv->deriv_classes[4][7][7] = int_stack + 45330;
 Libderiv->deriv_classes[4][8][7] = int_stack + 45870;
 Libderiv->deriv_classes[5][4][7] = int_stack + 46545;
 Libderiv->deriv_classes[5][5][7] = int_stack + 46860;
 Libderiv->deriv_classes[5][6][7] = int_stack + 47301;
 Libderiv->deriv_classes[5][7][7] = int_stack + 47889;
 Libderiv->deriv_classes[5][8][7] = int_stack + 48645;
 Libderiv->deriv_classes[6][4][7] = int_stack + 49590;
 Libderiv->deriv_classes[6][5][7] = int_stack + 50010;
 Libderiv->deriv_classes[6][6][7] = int_stack + 50598;
 Libderiv->deriv_classes[6][7][7] = int_stack + 51382;
 Libderiv->deriv_classes[6][8][7] = int_stack + 52390;
 Libderiv->deriv_classes[3][4][6] = int_stack + 53650;
 Libderiv->deriv_classes[3][5][6] = int_stack + 53800;
 Libderiv->deriv_classes[3][6][6] = int_stack + 54010;
 Libderiv->deriv_classes[3][7][6] = int_stack + 54290;
 Libderiv->deriv_classes[3][8][6] = int_stack + 54650;
 Libderiv->deriv_classes[4][4][6] = int_stack + 55100;
 Libderiv->deriv_classes[4][5][6] = int_stack + 55325;
 Libderiv->deriv_classes[4][6][6] = int_stack + 55640;
 Libderiv->deriv_classes[4][7][6] = int_stack + 56060;
 Libderiv->deriv_classes[4][8][6] = int_stack + 56600;
 Libderiv->deriv_classes[5][4][6] = int_stack + 57275;
 Libderiv->deriv_classes[5][5][6] = int_stack + 57590;
 Libderiv->deriv_classes[5][6][6] = int_stack + 58031;
 Libderiv->deriv_classes[5][7][6] = int_stack + 58619;
 Libderiv->deriv_classes[5][8][6] = int_stack + 59375;
 Libderiv->dvrr_classes[6][4] = int_stack + 60320;
 Libderiv->deriv_classes[6][4][6] = int_stack + 60740;
 Libderiv->dvrr_classes[6][5] = int_stack + 61160;
 Libderiv->deriv_classes[6][5][6] = int_stack + 61748;
 Libderiv->dvrr_classes[6][6] = int_stack + 62336;
 Libderiv->deriv_classes[6][6][6] = int_stack + 63120;
 Libderiv->dvrr_classes[6][7] = int_stack + 63904;
 Libderiv->deriv_classes[6][7][6] = int_stack + 64912;
 Libderiv->deriv_classes[6][8][6] = int_stack + 65920;
 Libderiv->deriv_classes[3][4][2] = int_stack + 67180;
 Libderiv->deriv_classes[3][5][2] = int_stack + 67330;
 Libderiv->deriv_classes[3][6][2] = int_stack + 67540;
 Libderiv->deriv_classes[3][7][2] = int_stack + 67820;
 Libderiv->deriv_classes[3][8][2] = int_stack + 68180;
 Libderiv->deriv_classes[4][4][2] = int_stack + 68630;
 Libderiv->deriv_classes[4][5][2] = int_stack + 68855;
 Libderiv->deriv_classes[4][6][2] = int_stack + 69170;
 Libderiv->deriv_classes[4][7][2] = int_stack + 69590;
 Libderiv->deriv_classes[4][8][2] = int_stack + 70130;
 Libderiv->deriv_classes[5][4][2] = int_stack + 70805;
 Libderiv->deriv_classes[5][5][2] = int_stack + 71120;
 Libderiv->deriv_classes[5][6][2] = int_stack + 71561;
 Libderiv->deriv_classes[5][7][2] = int_stack + 72149;
 Libderiv->deriv_classes[5][8][2] = int_stack + 72905;
 Libderiv->deriv_classes[6][4][2] = int_stack + 73850;
 Libderiv->deriv_classes[6][5][2] = int_stack + 74270;
 Libderiv->deriv_classes[6][6][2] = int_stack + 74858;
 Libderiv->deriv_classes[6][7][2] = int_stack + 75642;
 Libderiv->deriv_classes[6][8][2] = int_stack + 76650;
 Libderiv->deriv_classes[3][4][1] = int_stack + 77910;
 Libderiv->deriv_classes[3][5][1] = int_stack + 78060;
 Libderiv->deriv_classes[3][6][1] = int_stack + 78270;
 Libderiv->deriv_classes[3][7][1] = int_stack + 78550;
 Libderiv->deriv_classes[3][8][1] = int_stack + 78910;
 Libderiv->deriv_classes[4][4][1] = int_stack + 79360;
 Libderiv->deriv_classes[4][5][1] = int_stack + 79585;
 Libderiv->deriv_classes[4][6][1] = int_stack + 79900;
 Libderiv->deriv_classes[4][7][1] = int_stack + 80320;
 Libderiv->deriv_classes[4][8][1] = int_stack + 80860;
 Libderiv->deriv_classes[5][4][1] = int_stack + 81535;
 Libderiv->deriv_classes[5][5][1] = int_stack + 81850;
 Libderiv->deriv_classes[5][6][1] = int_stack + 82291;
 Libderiv->deriv_classes[5][7][1] = int_stack + 82879;
 Libderiv->deriv_classes[5][8][1] = int_stack + 83635;
 Libderiv->deriv_classes[6][4][1] = int_stack + 84580;
 Libderiv->deriv_classes[6][5][1] = int_stack + 85000;
 Libderiv->deriv_classes[6][6][1] = int_stack + 85588;
 Libderiv->deriv_classes[6][7][1] = int_stack + 86372;
 Libderiv->deriv_classes[6][8][1] = int_stack + 87380;
 Libderiv->dvrr_classes[3][4] = int_stack + 88640;
 Libderiv->dvrr_classes[3][5] = int_stack + 88790;
 Libderiv->dvrr_classes[3][6] = int_stack + 89000;
 Libderiv->dvrr_classes[3][7] = int_stack + 89280;
 Libderiv->dvrr_classes[3][8] = int_stack + 89640;
 Libderiv->deriv_classes[3][4][0] = int_stack + 90090;
 Libderiv->deriv_classes[3][5][0] = int_stack + 90240;
 Libderiv->deriv_classes[3][6][0] = int_stack + 90450;
 Libderiv->deriv_classes[3][7][0] = int_stack + 90730;
 Libderiv->deriv_classes[3][8][0] = int_stack + 91090;
 Libderiv->dvrr_classes[4][4] = int_stack + 91540;
 Libderiv->dvrr_classes[4][5] = int_stack + 91765;
 Libderiv->dvrr_classes[4][6] = int_stack + 92080;
 Libderiv->dvrr_classes[4][7] = int_stack + 92500;
 Libderiv->dvrr_classes[4][8] = int_stack + 93040;
 Libderiv->deriv_classes[4][4][0] = int_stack + 93715;
 Libderiv->deriv_classes[4][5][0] = int_stack + 93940;
 Libderiv->deriv_classes[4][6][0] = int_stack + 94255;
 Libderiv->deriv_classes[4][7][0] = int_stack + 94675;
 Libderiv->deriv_classes[4][8][0] = int_stack + 95215;
 Libderiv->dvrr_classes[5][4] = int_stack + 95890;
 Libderiv->dvrr_classes[5][5] = int_stack + 96205;
 Libderiv->dvrr_classes[5][6] = int_stack + 96646;
 Libderiv->dvrr_classes[5][7] = int_stack + 97234;
 Libderiv->dvrr_classes[5][8] = int_stack + 97990;
 Libderiv->deriv_classes[5][4][0] = int_stack + 98935;
 Libderiv->deriv_classes[5][5][0] = int_stack + 99250;
 Libderiv->deriv_classes[5][6][0] = int_stack + 99691;
 Libderiv->deriv_classes[5][7][0] = int_stack + 100279;
 Libderiv->deriv_classes[5][8][0] = int_stack + 101035;
 Libderiv->deriv_classes[6][4][0] = int_stack + 101980;
 Libderiv->deriv_classes[6][5][0] = int_stack + 102400;
 Libderiv->deriv_classes[6][6][0] = int_stack + 102988;
 Libderiv->deriv_classes[6][7][0] = int_stack + 103772;
 Libderiv->deriv_classes[6][8][0] = int_stack + 104780;
 memset(int_stack,0,848320);

 Libderiv->dvrr_stack = int_stack + 443711;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ffgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+106040,int_stack+88790,int_stack+88640,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+106490,int_stack+89000,int_stack+88790,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+107120,int_stack+106490,int_stack+106040,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+108020,int_stack+89280,int_stack+89000,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+108860,int_stack+108020,int_stack+106490,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+110120,int_stack+108860,int_stack+107120,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+111620,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88640,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+112070,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88790,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+112700,int_stack+112070,int_stack+111620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106040,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+113600,int_stack+640,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89000,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+114440,int_stack+113600,int_stack+112070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106490,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+115700,int_stack+114440,int_stack+112700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107120,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+111620,int_stack+1000,int_stack+640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89280,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+117200,int_stack+111620,int_stack+113600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108020,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+111620,int_stack+117200,int_stack+114440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108860,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+117200,int_stack+111620,int_stack+115700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110120,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+111620,int_stack+91765,int_stack+91540,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+112295,int_stack+92080,int_stack+91765,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+113240,int_stack+112295,int_stack+111620,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+114590,int_stack+92500,int_stack+92080,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+119450,int_stack+114590,int_stack+112295,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+121340,int_stack+119450,int_stack+113240,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+115850,int_stack+1675,int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+1990,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91765,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+123590,int_stack+0,int_stack+115850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111620,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+115850,int_stack+2410,int_stack+1990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92080,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+124940,int_stack+115850,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112295,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+124940,int_stack+123590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113240,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+126830,int_stack+2950,int_stack+2410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92500,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+128450,int_stack+126830,int_stack+115850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114590,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+130970,int_stack+128450,int_stack+124940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119450,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+123590,int_stack+130970,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121340,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+126965,int_stack+123590,int_stack+117200,225);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+96205,int_stack+95890,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+945,int_stack+96646,int_stack+96205,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+115850,int_stack+945,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+133715,int_stack+97234,int_stack+96646,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+135479,int_stack+133715,int_stack+945,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+138125,int_stack+135479,int_stack+115850,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2268,int_stack+3940,int_stack+3625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95890,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+117740,int_stack+4381,int_stack+3940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96205,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+141275,int_stack+117740,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2268,int_stack+4969,int_stack+4381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96646,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+143165,int_stack+2268,int_stack+117740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+145811,int_stack+143165,int_stack+141275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115850,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+148961,int_stack+5725,int_stack+4969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97234,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+151229,int_stack+148961,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133715,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+154757,int_stack+151229,int_stack+143165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135479,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+148961,int_stack+154757,int_stack+145811, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138125,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+153686,int_stack+148961,int_stack+123590,225);
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+163811,int_stack+153686,int_stack+126965,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+123590,int_stack+61160,int_stack+60320,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+124850,int_stack+62336,int_stack+61160,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+126614,int_stack+124850,int_stack+123590,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+129134,int_stack+63904,int_stack+62336,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+2268,int_stack+129134,int_stack+124850,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+141275,int_stack+2268,int_stack+126614,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+131486,int_stack+7090,int_stack+6670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60320,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+145475,int_stack+7678,int_stack+7090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61160,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+177311,int_stack+145475,int_stack+131486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123590,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+179831,int_stack+8462,int_stack+7678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62336,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+182183,int_stack+179831,int_stack+145475, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124850,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+185711,int_stack+182183,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126614,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+145475,int_stack+9470,int_stack+8462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63904,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+5796,int_stack+145475,int_stack+179831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129134,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+189911,int_stack+5796,int_stack+182183, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2268,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+177311,int_stack+189911,int_stack+185711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141275,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+183611,int_stack+177311,int_stack+148961,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+197786,int_stack+183611,int_stack+153686,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+10880,int_stack+10730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88640, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177761,int_stack+11090,int_stack+10880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88790, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178391,int_stack+177761,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106040, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+179291,int_stack+11370,int_stack+11090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89000, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+180131,int_stack+179291,int_stack+177761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106490, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+181391,int_stack+180131,int_stack+178391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107120, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+11730,int_stack+11370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89280, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+179291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108020, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+182891,int_stack+180131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108860, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+181391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110120, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+12405,int_stack+12180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177986,int_stack+12720,int_stack+12405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91765, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178931,int_stack+177986,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111620, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+180281,int_stack+13140,int_stack+12720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92080, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+185141,int_stack+180281,int_stack+177986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112295, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+187031,int_stack+185141,int_stack+178931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113240, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+13680,int_stack+13140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92500, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+189281,int_stack+177311,int_stack+180281, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114590, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+189281,int_stack+185141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119450, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+189281,int_stack+177311,int_stack+187031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121340, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5796,int_stack+189281,int_stack+182891,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+14670,int_stack+14355, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95890, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+178256,int_stack+15111,int_stack+14670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96205, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+179579,int_stack+178256,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+181469,int_stack+15699,int_stack+15111, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96646, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+183233,int_stack+181469,int_stack+178256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+185879,int_stack+183233,int_stack+179579, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115850, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+16455,int_stack+15699, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97234, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+192656,int_stack+177311,int_stack+181469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133715, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+192656,int_stack+183233, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135479, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+192656,int_stack+177311,int_stack+185879, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138125, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+177311,int_stack+192656,int_stack+189281,225);
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+145475,int_stack+177311,int_stack+5796,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5796,int_stack+17820,int_stack+17400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60320, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7056,int_stack+18408,int_stack+17820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61160, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8820,int_stack+7056,int_stack+5796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123590, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+11340,int_stack+19192,int_stack+18408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62336, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13692,int_stack+11340,int_stack+7056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124850, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+187436,int_stack+13692,int_stack+8820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126614, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+5796,int_stack+20200,int_stack+19192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63904, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+158975,int_stack+5796,int_stack+11340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129134, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+5796,int_stack+158975,int_stack+13692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2268, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+11676,int_stack+5796,int_stack+187436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141275, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+218036,int_stack+11676,int_stack+192656,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+232211,int_stack+218036,int_stack+177311,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+21610,int_stack+21460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88640, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177761,int_stack+21820,int_stack+21610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88790, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178391,int_stack+177761,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106040, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+179291,int_stack+22100,int_stack+21820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89000, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+180131,int_stack+179291,int_stack+177761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106490, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+181391,int_stack+180131,int_stack+178391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107120, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+22460,int_stack+22100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89280, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+179291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108020, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+182891,int_stack+180131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108860, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+181391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110120, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+23135,int_stack+22910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177986,int_stack+23450,int_stack+23135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91765, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178931,int_stack+177986,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111620, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+180281,int_stack+23870,int_stack+23450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92080, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+185141,int_stack+180281,int_stack+177986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112295, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+187031,int_stack+185141,int_stack+178931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113240, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+24410,int_stack+23870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92500, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+189281,int_stack+177311,int_stack+180281, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114590, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+189281,int_stack+185141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119450, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+189281,int_stack+177311,int_stack+187031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121340, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+218036,int_stack+189281,int_stack+182891,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+25400,int_stack+25085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95890, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+178256,int_stack+25841,int_stack+25400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96205, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+179579,int_stack+178256,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+181469,int_stack+26429,int_stack+25841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96646, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+183233,int_stack+181469,int_stack+178256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+185879,int_stack+183233,int_stack+179579, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115850, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+27185,int_stack+26429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97234, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+192656,int_stack+177311,int_stack+181469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133715, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+192656,int_stack+183233, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135479, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+192656,int_stack+177311,int_stack+185879, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138125, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+177311,int_stack+192656,int_stack+189281,225);
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+5796,int_stack+177311,int_stack+218036,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+218036,int_stack+28550,int_stack+28130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60320, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+219296,int_stack+29138,int_stack+28550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61160, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+221060,int_stack+219296,int_stack+218036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123590, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+223580,int_stack+29922,int_stack+29138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62336, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+225932,int_stack+223580,int_stack+219296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124850, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+187436,int_stack+225932,int_stack+221060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126614, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+218036,int_stack+30930,int_stack+29922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63904, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+19296,int_stack+218036,int_stack+223580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129134, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+218036,int_stack+19296,int_stack+225932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2268, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+19296,int_stack+218036,int_stack+187436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141275, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+218036,int_stack+19296,int_stack+192656,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+252461,int_stack+218036,int_stack+177311,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+32340,int_stack+32190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177761,int_stack+32550,int_stack+32340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178391,int_stack+177761,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+179291,int_stack+32830,int_stack+32550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+180131,int_stack+179291,int_stack+177761, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+181391,int_stack+180131,int_stack+178391, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+33190,int_stack+32830, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+179291, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+182891,int_stack+180131, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+181391, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+33865,int_stack+33640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177986,int_stack+34180,int_stack+33865, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178931,int_stack+177986,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+180281,int_stack+34600,int_stack+34180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+185141,int_stack+180281,int_stack+177986, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+187031,int_stack+185141,int_stack+178931, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+35140,int_stack+34600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+189281,int_stack+177311,int_stack+180281, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+189281,int_stack+185141, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+189281,int_stack+177311,int_stack+187031, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+218036,int_stack+189281,int_stack+182891,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+36130,int_stack+35815, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+178256,int_stack+36571,int_stack+36130, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+179579,int_stack+178256,int_stack+177311, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+181469,int_stack+37159,int_stack+36571, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+183233,int_stack+181469,int_stack+178256, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+185879,int_stack+183233,int_stack+179579, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+37915,int_stack+37159, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+192656,int_stack+177311,int_stack+181469, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+192656,int_stack+183233, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135479, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+192656,int_stack+177311,int_stack+185879, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+177311,int_stack+192656,int_stack+189281,225);
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+19296,int_stack+177311,int_stack+218036,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+218036,int_stack+39280,int_stack+38860, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+219296,int_stack+39868,int_stack+39280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+221060,int_stack+219296,int_stack+218036, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+223580,int_stack+40652,int_stack+39868, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+225932,int_stack+223580,int_stack+219296, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+124850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+187436,int_stack+225932,int_stack+221060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+218036,int_stack+41660,int_stack+40652, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+32796,int_stack+218036,int_stack+223580, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+218036,int_stack+32796,int_stack+225932, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+32796,int_stack+218036,int_stack+187436, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+218036,int_stack+32796,int_stack+192656,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+272711,int_stack+218036,int_stack+177311,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+43070,int_stack+42920, 0.0, zero_stack, 1.0, int_stack+88640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177761,int_stack+43280,int_stack+43070, 0.0, zero_stack, 1.0, int_stack+88790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178391,int_stack+177761,int_stack+177311, 0.0, zero_stack, 1.0, int_stack+106040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+179291,int_stack+43560,int_stack+43280, 0.0, zero_stack, 1.0, int_stack+89000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+180131,int_stack+179291,int_stack+177761, 0.0, zero_stack, 1.0, int_stack+106490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+181391,int_stack+180131,int_stack+178391, 0.0, zero_stack, 1.0, int_stack+107120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+43920,int_stack+43560, 0.0, zero_stack, 1.0, int_stack+89280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+179291, 0.0, zero_stack, 1.0, int_stack+108020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+182891,int_stack+180131, 0.0, zero_stack, 1.0, int_stack+108860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+181391, 0.0, zero_stack, 1.0, int_stack+110120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+44595,int_stack+44370, 0.0, zero_stack, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177986,int_stack+44910,int_stack+44595, 0.0, zero_stack, 1.0, int_stack+91765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178931,int_stack+177986,int_stack+177311, 0.0, zero_stack, 1.0, int_stack+111620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+180281,int_stack+45330,int_stack+44910, 0.0, zero_stack, 1.0, int_stack+92080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+185141,int_stack+180281,int_stack+177986, 0.0, zero_stack, 1.0, int_stack+112295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+187031,int_stack+185141,int_stack+178931, 0.0, zero_stack, 1.0, int_stack+113240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+45870,int_stack+45330, 0.0, zero_stack, 1.0, int_stack+92500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+189281,int_stack+177311,int_stack+180281, 0.0, zero_stack, 1.0, int_stack+114590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+189281,int_stack+185141, 0.0, zero_stack, 1.0, int_stack+119450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+189281,int_stack+177311,int_stack+187031, 0.0, zero_stack, 1.0, int_stack+121340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+218036,int_stack+189281,int_stack+182891,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+46860,int_stack+46545, 0.0, zero_stack, 1.0, int_stack+95890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+178256,int_stack+47301,int_stack+46860, 0.0, zero_stack, 1.0, int_stack+96205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+179579,int_stack+178256,int_stack+177311, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+181469,int_stack+47889,int_stack+47301, 0.0, zero_stack, 1.0, int_stack+96646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+183233,int_stack+181469,int_stack+178256, 0.0, zero_stack, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+185879,int_stack+183233,int_stack+179579, 0.0, zero_stack, 1.0, int_stack+115850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+177311,int_stack+48645,int_stack+47889, 0.0, zero_stack, 1.0, int_stack+97234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+192656,int_stack+177311,int_stack+181469, 0.0, zero_stack, 1.0, int_stack+133715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+177311,int_stack+192656,int_stack+183233, 0.0, zero_stack, 1.0, int_stack+135479, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+192656,int_stack+177311,int_stack+185879, 0.0, zero_stack, 1.0, int_stack+138125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+177311,int_stack+192656,int_stack+189281,225);
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+32796,int_stack+177311,int_stack+218036,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+218036,int_stack+50010,int_stack+49590, 0.0, zero_stack, 1.0, int_stack+60320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+219296,int_stack+50598,int_stack+50010, 0.0, zero_stack, 1.0, int_stack+61160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+221060,int_stack+219296,int_stack+218036, 0.0, zero_stack, 1.0, int_stack+123590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+223580,int_stack+51382,int_stack+50598, 0.0, zero_stack, 1.0, int_stack+62336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+225932,int_stack+223580,int_stack+219296, 0.0, zero_stack, 1.0, int_stack+124850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+187436,int_stack+225932,int_stack+221060, 0.0, zero_stack, 1.0, int_stack+126614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+218036,int_stack+52390,int_stack+51382, 0.0, zero_stack, 1.0, int_stack+63904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+46296,int_stack+218036,int_stack+223580, 0.0, zero_stack, 1.0, int_stack+129134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+218036,int_stack+46296,int_stack+225932, 0.0, zero_stack, 1.0, int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+46296,int_stack+218036,int_stack+187436, 0.0, zero_stack, 1.0, int_stack+141275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+218036,int_stack+46296,int_stack+192656,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+292961,int_stack+218036,int_stack+177311,225);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+53800,int_stack+53650, 1.0, int_stack+88640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177761,int_stack+54010,int_stack+53800, 1.0, int_stack+88790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+178391,int_stack+177761,int_stack+177311, 1.0, int_stack+106040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+179291,int_stack+54290,int_stack+54010, 1.0, int_stack+89000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+180131,int_stack+179291,int_stack+177761, 1.0, int_stack+106490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+181391,int_stack+180131,int_stack+178391, 1.0, int_stack+107120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+106040,int_stack+54650,int_stack+54290, 1.0, int_stack+89280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+177311,int_stack+106040,int_stack+179291, 1.0, int_stack+108020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+182891,int_stack+177311,int_stack+180131, 1.0, int_stack+108860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+177311,int_stack+182891,int_stack+181391, 1.0, int_stack+110120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+179561,int_stack+55325,int_stack+55100, 1.0, int_stack+91540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+180236,int_stack+55640,int_stack+55325, 1.0, int_stack+91765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+181181,int_stack+180236,int_stack+179561, 1.0, int_stack+111620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+182531,int_stack+56060,int_stack+55640, 1.0, int_stack+92080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+183791,int_stack+182531,int_stack+180236, 1.0, int_stack+112295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+185681,int_stack+183791,int_stack+181181, 1.0, int_stack+113240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+111620,int_stack+56600,int_stack+56060, 1.0, int_stack+92500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+179561,int_stack+111620,int_stack+182531, 1.0, int_stack+114590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+187931,int_stack+179561,int_stack+183791, 1.0, int_stack+119450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+179561,int_stack+187931,int_stack+185681, 1.0, int_stack+121340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+182936,int_stack+179561,int_stack+177311,225);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+177311,int_stack+57590,int_stack+57275, 1.0, int_stack+95890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+189686,int_stack+58031,int_stack+57590, 1.0, int_stack+96205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+191009,int_stack+189686,int_stack+177311, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+177311,int_stack+58619,int_stack+58031, 1.0, int_stack+96646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+192899,int_stack+177311,int_stack+189686, 1.0, int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+218036,int_stack+192899,int_stack+191009, 1.0, int_stack+115850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+59375,int_stack+58619, 1.0, int_stack+97234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+115850,int_stack+0,int_stack+177311, 1.0, int_stack+133715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+221186,int_stack+115850,int_stack+192899, 1.0, int_stack+135479, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+189686,int_stack+221186,int_stack+218036, 1.0, int_stack+138125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+218036,int_stack+189686,int_stack+179561,225);
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+46296,int_stack+218036,int_stack+182936,225);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+228161,int_stack+61748,int_stack+60740, 1.0, int_stack+60320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+229421,int_stack+63120,int_stack+61748, 1.0, int_stack+61160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+115850,int_stack+229421,int_stack+228161, 1.0, int_stack+123590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+177311,int_stack+64912,int_stack+63120, 1.0, int_stack+62336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+179663,int_stack+177311,int_stack+229421, 1.0, int_stack+124850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+183191,int_stack+179663,int_stack+115850, 1.0, int_stack+126614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+115850,int_stack+65920,int_stack+64912, 1.0, int_stack+63904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+123590,int_stack+115850,int_stack+177311, 1.0, int_stack+129134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+59796,int_stack+123590,int_stack+179663, 1.0, int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+123590,int_stack+59796,int_stack+183191, 1.0, int_stack+141275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gg) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+313211,int_stack+123590,int_stack+189686,225);
 /*--- compute (gd|gg) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+177311,int_stack+313211,int_stack+218036,225);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+218036,int_stack+89640,int_stack+89280,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+219116,int_stack+218036,int_stack+108020,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+220796,int_stack+219116,int_stack+108860,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+218036,int_stack+220796,int_stack+110120,10);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+220286,int_stack+93040,int_stack+92500,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+221906,int_stack+220286,int_stack+114590,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+224426,int_stack+221906,int_stack+119450,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+220286,int_stack+224426,int_stack+121340,15);
 /*--- compute (fp|gg) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+223661,int_stack+220286,int_stack+218036,225);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+313211,int_stack+97990,int_stack+97234,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+315479,int_stack+313211,int_stack+133715,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+319007,int_stack+315479,int_stack+135479,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+313211,int_stack+319007,int_stack+138125,21);
 /*--- compute (gp|gg) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+106040,int_stack+313211,int_stack+220286,225);
 /*--- compute (fd|gg) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+116165,int_stack+106040,int_stack+223661,225);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+317936,int_stack+67330,int_stack+67180,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+318386,int_stack+67540,int_stack+67330,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+319016,int_stack+318386,int_stack+317936,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+319916,int_stack+67820,int_stack+67540,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+320756,int_stack+319916,int_stack+318386,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+322016,int_stack+320756,int_stack+319016,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+317936,int_stack+68180,int_stack+67820,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+323516,int_stack+317936,int_stack+319916,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+317936,int_stack+323516,int_stack+320756,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+323516,int_stack+317936,int_stack+322016,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+317936,int_stack+68855,int_stack+68630,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+318611,int_stack+69170,int_stack+68855,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+319556,int_stack+318611,int_stack+317936,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+320906,int_stack+69590,int_stack+69170,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+325766,int_stack+320906,int_stack+318611,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+59796,int_stack+325766,int_stack+319556,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+317936,int_stack+70130,int_stack+69590,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+62046,int_stack+317936,int_stack+320906,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+317936,int_stack+62046,int_stack+325766,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+62046,int_stack+317936,int_stack+59796,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+129665,int_stack+62046,int_stack+323516, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+218036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+59796,int_stack+71120,int_stack+70805,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+317936,int_stack+71561,int_stack+71120,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+319259,int_stack+317936,int_stack+59796,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+59796,int_stack+72149,int_stack+71561,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+321149,int_stack+59796,int_stack+317936,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+323795,int_stack+321149,int_stack+319259,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+317936,int_stack+72905,int_stack+72149,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+65421,int_stack+317936,int_stack+59796,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+68949,int_stack+65421,int_stack+321149,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+317936,int_stack+68949,int_stack+323795,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+322661,int_stack+317936,int_stack+62046, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (fd|gg) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+59796,int_stack+322661,int_stack+129665, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+223661, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+129665,int_stack+74270,int_stack+73850,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130925,int_stack+74858,int_stack+74270,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+132689,int_stack+130925,int_stack+129665,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+135209,int_stack+75642,int_stack+74858,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+137561,int_stack+135209,int_stack+130925,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+141089,int_stack+137561,int_stack+132689,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+129665,int_stack+76650,int_stack+75642,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+0,int_stack+129665,int_stack+135209,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+129665,int_stack+0,int_stack+137561,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+332786,int_stack+129665,int_stack+141089,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+129665,int_stack+332786,int_stack+317936, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+313211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+332786,int_stack+129665,int_stack+322661, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+129665,int_stack+78060,int_stack+77910,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130115,int_stack+78270,int_stack+78060,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+130745,int_stack+130115,int_stack+129665,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+131645,int_stack+78550,int_stack+78270,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+132485,int_stack+131645,int_stack+130115,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+133745,int_stack+132485,int_stack+130745,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+129665,int_stack+78910,int_stack+78550,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+135245,int_stack+129665,int_stack+131645,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+129665,int_stack+135245,int_stack+132485,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+135245,int_stack+129665,int_stack+133745,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+129665,int_stack+79585,int_stack+79360,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130340,int_stack+79900,int_stack+79585,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+131285,int_stack+130340,int_stack+129665,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+132635,int_stack+80320,int_stack+79900,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+137495,int_stack+132635,int_stack+130340,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+139385,int_stack+137495,int_stack+131285,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+129665,int_stack+80860,int_stack+80320,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+141635,int_stack+129665,int_stack+132635,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+129665,int_stack+141635,int_stack+137495,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+141635,int_stack+129665,int_stack+139385,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+317936,int_stack+141635,int_stack+135245, 0.0, zero_stack, 1.0, int_stack+218036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+129665,int_stack+81850,int_stack+81535,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130610,int_stack+82291,int_stack+81850,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+131933,int_stack+130610,int_stack+129665,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+133823,int_stack+82879,int_stack+82291,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+135587,int_stack+133823,int_stack+130610,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+138233,int_stack+135587,int_stack+131933,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+129665,int_stack+83635,int_stack+82879,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+324686,int_stack+129665,int_stack+133823,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+129665,int_stack+324686,int_stack+135587,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+324686,int_stack+129665,int_stack+138233,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+129665,int_stack+324686,int_stack+141635, 0.0, zero_stack, 1.0, int_stack+220286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (fd|gg) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+353036,int_stack+129665,int_stack+317936, 0.0, zero_stack, 1.0, int_stack+223661, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+317936,int_stack+85000,int_stack+84580,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+319196,int_stack+85588,int_stack+85000,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+320960,int_stack+319196,int_stack+317936,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+139790,int_stack+86372,int_stack+85588,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+0,int_stack+139790,int_stack+319196,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+73296,int_stack+0,int_stack+320960,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+317936,int_stack+87380,int_stack+86372,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+77496,int_stack+317936,int_stack+139790,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+317936,int_stack+77496,int_stack+0,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+77496,int_stack+317936,int_stack+73296,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+366536,int_stack+77496,int_stack+324686, 0.0, zero_stack, 1.0, int_stack+313211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+380711,int_stack+366536,int_stack+129665, 0.0, zero_stack, 1.0, int_stack+106040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+129665,int_stack+90240,int_stack+90090,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130115,int_stack+90450,int_stack+90240,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+130745,int_stack+130115,int_stack+129665,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+131645,int_stack+90730,int_stack+90450,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+132485,int_stack+131645,int_stack+130115,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+133745,int_stack+132485,int_stack+130745,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+129665,int_stack+91090,int_stack+90730,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+135245,int_stack+129665,int_stack+131645,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+129665,int_stack+135245,int_stack+132485,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+135245,int_stack+129665,int_stack+133745,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+129665,int_stack+93940,int_stack+93715,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130340,int_stack+94255,int_stack+93940,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+131285,int_stack+130340,int_stack+129665,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+132635,int_stack+94675,int_stack+94255,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+137495,int_stack+132635,int_stack+130340,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+139385,int_stack+137495,int_stack+131285,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+129665,int_stack+95215,int_stack+94675,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+141635,int_stack+129665,int_stack+132635,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+129665,int_stack+141635,int_stack+137495,15);
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+141635,int_stack+129665,int_stack+139385,15);
 /*--- compute (fp|gg) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+366536,int_stack+141635,int_stack+135245, 1.0, int_stack+218036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+218036,int_stack+99250,int_stack+98935,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+129665,int_stack+99691,int_stack+99250,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+130988,int_stack+129665,int_stack+218036,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+218036,int_stack+100279,int_stack+99691,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+132878,int_stack+218036,int_stack+129665,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+135524,int_stack+132878,int_stack+130988,21);
 /*--- compute (h0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+129665,int_stack+101035,int_stack+100279,21);
 /*--- compute (h0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+373286,int_stack+129665,int_stack+218036,21);
 /*--- compute (h0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+73296,int_stack+373286,int_stack+132878,21);
 /*--- compute (h0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+373286,int_stack+73296,int_stack+135524,21);
 /*--- compute (gp|gg) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+73296,int_stack+373286,int_stack+141635, 1.0, int_stack+220286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (fd|gg) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+83421,int_stack+73296,int_stack+366536, 1.0, int_stack+223661, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+366536,int_stack+102400,int_stack+101980,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+367796,int_stack+102988,int_stack+102400,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+369560,int_stack+367796,int_stack+366536,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+96921,int_stack+103772,int_stack+102988,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+99273,int_stack+96921,int_stack+367796,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+218036,int_stack+99273,int_stack+369560,28);
 /*--- compute (i0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+366536,int_stack+104780,int_stack+103772,28);
 /*--- compute (i0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+222236,int_stack+366536,int_stack+96921,28);
 /*--- compute (i0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+366536,int_stack+222236,int_stack+99273,28);
 /*--- compute (i0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+222236,int_stack+366536,int_stack+218036,28);
 /*--- compute (hp|gg) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+129665,int_stack+222236,int_stack+373286, 1.0, int_stack+313211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (gd|gg) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+400961,int_stack+129665,int_stack+73296, 1.0, int_stack+106040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
 /*--- compute (ff|gg) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+421211,int_stack+197786,int_stack+163811,225);
     Libderiv->ABCD[11] = int_stack + 421211;
 /*--- compute (ff|gg) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+197561,int_stack+232211,int_stack+145475,225);
     Libderiv->ABCD[10] = int_stack + 197561;
 /*--- compute (ff|gg) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+129665,int_stack+252461,int_stack+5796,225);
     Libderiv->ABCD[9] = int_stack + 129665;
 /*--- compute (ff|gg) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+152165,int_stack+272711,int_stack+19296,225);
     Libderiv->ABCD[8] = int_stack + 152165;
 /*--- compute (ff|gg) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+0,int_stack+292961,int_stack+32796,225);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (ff|gg) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+22500,int_stack+177311,int_stack+46296,225);
     Libderiv->ABCD[6] = int_stack + 22500;
 /*--- compute (ff|gg) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+174665,int_stack+332786,int_stack+59796, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 174665;
 /*--- compute (ff|gg) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+45000,int_stack+380711,int_stack+353036, 0.0, zero_stack, 1.0, int_stack+116165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 45000;
 /*--- compute (ff|gg) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+220061,int_stack+400961,int_stack+83421, 1.0, int_stack+116165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 220061;

}
