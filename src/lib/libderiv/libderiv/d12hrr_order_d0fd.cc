#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_d0fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|fd) integrals */

void d12hrr_order_d0fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][5][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][10] = int_stack + 126;
 Libderiv->deriv_classes[2][5][9] = int_stack + 252;
 Libderiv->deriv_classes[2][5][8] = int_stack + 378;
 Libderiv->deriv_classes[2][5][7] = int_stack + 504;
 Libderiv->dvrr_classes[2][4] = int_stack + 630;
 Libderiv->deriv_classes[2][5][6] = int_stack + 720;
 Libderiv->deriv_classes[2][5][2] = int_stack + 846;
 Libderiv->deriv_classes[2][5][1] = int_stack + 972;
 Libderiv->deriv_classes[2][5][0] = int_stack + 1098;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 1224;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 1284;
 Libderiv->deriv2_classes[2][5][143] = int_stack + 1374;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1500;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 1560;
 Libderiv->deriv2_classes[2][5][131] = int_stack + 1650;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 1776;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 1836;
 Libderiv->deriv2_classes[2][5][130] = int_stack + 1926;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 2052;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 2112;
 Libderiv->deriv2_classes[2][5][119] = int_stack + 2202;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 2328;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 2388;
 Libderiv->deriv2_classes[2][5][118] = int_stack + 2478;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 2604;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 2664;
 Libderiv->deriv2_classes[2][5][117] = int_stack + 2754;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 2880;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 2940;
 Libderiv->deriv2_classes[2][5][107] = int_stack + 3030;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 3156;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 3216;
 Libderiv->deriv2_classes[2][5][106] = int_stack + 3306;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 3432;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 3492;
 Libderiv->deriv2_classes[2][5][105] = int_stack + 3582;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 3708;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 3768;
 Libderiv->deriv2_classes[2][5][104] = int_stack + 3858;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 3984;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 4044;
 Libderiv->deriv2_classes[2][5][95] = int_stack + 4134;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 4260;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 4320;
 Libderiv->deriv2_classes[2][5][94] = int_stack + 4410;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 4536;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 4596;
 Libderiv->deriv2_classes[2][5][93] = int_stack + 4686;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 4812;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 4872;
 Libderiv->deriv2_classes[2][5][92] = int_stack + 4962;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 5088;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 5148;
 Libderiv->deriv2_classes[2][5][91] = int_stack + 5238;
 Libderiv->deriv_classes[2][3][11] = int_stack + 5364;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 5424;
 Libderiv->deriv_classes[2][4][11] = int_stack + 5484;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 5574;
 Libderiv->deriv2_classes[2][5][83] = int_stack + 5664;
 Libderiv->deriv_classes[2][3][10] = int_stack + 5790;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 5850;
 Libderiv->deriv_classes[2][4][10] = int_stack + 5910;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 6000;
 Libderiv->deriv2_classes[2][5][82] = int_stack + 6090;
 Libderiv->deriv_classes[2][3][9] = int_stack + 6216;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 6276;
 Libderiv->deriv_classes[2][4][9] = int_stack + 6336;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 6426;
 Libderiv->deriv2_classes[2][5][81] = int_stack + 6516;
 Libderiv->deriv_classes[2][3][8] = int_stack + 6642;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 6702;
 Libderiv->deriv_classes[2][4][8] = int_stack + 6762;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 6852;
 Libderiv->deriv2_classes[2][5][80] = int_stack + 6942;
 Libderiv->deriv_classes[2][3][7] = int_stack + 7068;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 7128;
 Libderiv->deriv_classes[2][4][7] = int_stack + 7188;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 7278;
 Libderiv->deriv2_classes[2][5][79] = int_stack + 7368;
 Libderiv->dvrr_classes[2][3] = int_stack + 7494;
 Libderiv->deriv_classes[2][3][6] = int_stack + 7554;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 7614;
 Libderiv->deriv_classes[2][4][6] = int_stack + 7674;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 7764;
 Libderiv->deriv2_classes[2][5][78] = int_stack + 7854;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 7980;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 8040;
 Libderiv->deriv2_classes[2][5][35] = int_stack + 8130;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 8256;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 8316;
 Libderiv->deriv2_classes[2][5][34] = int_stack + 8406;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 8532;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 8592;
 Libderiv->deriv2_classes[2][5][33] = int_stack + 8682;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 8808;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 8868;
 Libderiv->deriv2_classes[2][5][32] = int_stack + 8958;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 9084;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 9144;
 Libderiv->deriv2_classes[2][5][31] = int_stack + 9234;
 Libderiv->deriv_classes[2][3][2] = int_stack + 9360;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 9420;
 Libderiv->deriv_classes[2][4][2] = int_stack + 9480;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 9570;
 Libderiv->deriv2_classes[2][5][30] = int_stack + 9660;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 9786;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 9846;
 Libderiv->deriv2_classes[2][5][26] = int_stack + 9936;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 10062;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 10122;
 Libderiv->deriv2_classes[2][5][23] = int_stack + 10212;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 10338;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 10398;
 Libderiv->deriv2_classes[2][5][22] = int_stack + 10488;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 10614;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 10674;
 Libderiv->deriv2_classes[2][5][21] = int_stack + 10764;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 10890;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 10950;
 Libderiv->deriv2_classes[2][5][20] = int_stack + 11040;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 11166;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 11226;
 Libderiv->deriv2_classes[2][5][19] = int_stack + 11316;
 Libderiv->deriv_classes[2][3][1] = int_stack + 11442;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 11502;
 Libderiv->deriv_classes[2][4][1] = int_stack + 11562;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 11652;
 Libderiv->deriv2_classes[2][5][18] = int_stack + 11742;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 11868;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 11928;
 Libderiv->deriv2_classes[2][5][14] = int_stack + 12018;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 12144;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 12204;
 Libderiv->deriv2_classes[2][5][13] = int_stack + 12294;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 12420;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 12480;
 Libderiv->deriv2_classes[2][5][11] = int_stack + 12570;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 12696;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 12756;
 Libderiv->deriv2_classes[2][5][10] = int_stack + 12846;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 12972;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 13032;
 Libderiv->deriv2_classes[2][5][9] = int_stack + 13122;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 13248;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 13308;
 Libderiv->deriv2_classes[2][5][8] = int_stack + 13398;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 13524;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 13584;
 Libderiv->deriv2_classes[2][5][7] = int_stack + 13674;
 Libderiv->deriv_classes[2][3][0] = int_stack + 13800;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 13860;
 Libderiv->deriv_classes[2][4][0] = int_stack + 13920;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 14010;
 Libderiv->deriv2_classes[2][5][6] = int_stack + 14100;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 14226;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 14286;
 Libderiv->deriv2_classes[2][5][2] = int_stack + 14376;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 14502;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 14562;
 Libderiv->deriv2_classes[2][5][1] = int_stack + 14652;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 14778;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 14838;
 Libderiv->deriv2_classes[2][5][0] = int_stack + 14928;
 memset(int_stack,0,120432);

 Libderiv->dvrr_stack = int_stack + 26304;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_d0fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15054,int_stack+630,int_stack+7494,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15234,int_stack+5484,int_stack+5364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7494,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15414,int_stack+0,int_stack+5484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+630,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15684,int_stack+5910,int_stack+5790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7494, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15864,int_stack+126,int_stack+5910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+6336,int_stack+6216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7494, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16134,int_stack+252,int_stack+6336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+6762,int_stack+6642, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16404,int_stack+378,int_stack+6762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16674,int_stack+7188,int_stack+7068, 0.0, zero_stack, 1.0, int_stack+7494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16854,int_stack+504,int_stack+7188, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+7674,int_stack+7554, 1.0, int_stack+7494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17124,int_stack+720,int_stack+7674, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+9480,int_stack+9360,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17394,int_stack+846,int_stack+9480,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+11562,int_stack+11442,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17664,int_stack+972,int_stack+11562,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+13920,int_stack+13800,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17934,int_stack+1098,int_stack+13920,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18204,int_stack+1284,int_stack+1224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5364,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18384,int_stack+1374,int_stack+1284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5484,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+1560,int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5364, 1.0, int_stack+5790,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+1650,int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5484, 1.0, int_stack+5910,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1530,int_stack+1836,int_stack+1776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5790, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18654,int_stack+1926,int_stack+1836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5910, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1710,int_stack+2112,int_stack+2052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5364, 0.0, zero_stack, 1.0, int_stack+6216,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18924,int_stack+2202,int_stack+2112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5484, 0.0, zero_stack, 1.0, int_stack+6336,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1890,int_stack+2388,int_stack+2328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5790, 1.0, int_stack+6216, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2070,int_stack+2478,int_stack+2388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5910, 1.0, int_stack+6336, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2340,int_stack+2664,int_stack+2604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6216, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19194,int_stack+2754,int_stack+2664, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6336, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2520,int_stack+2940,int_stack+2880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5364, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6642,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19464,int_stack+3030,int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5484, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6762,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+3216,int_stack+3156, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5790, 0.0, zero_stack, 1.0, int_stack+6642, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2880,int_stack+3306,int_stack+3216, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5910, 0.0, zero_stack, 1.0, int_stack+6762, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3150,int_stack+3492,int_stack+3432, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6216, 1.0, int_stack+6642, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19734,int_stack+3582,int_stack+3492, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6336, 1.0, int_stack+6762, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3330,int_stack+3768,int_stack+3708, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20004,int_stack+3858,int_stack+3768, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3510,int_stack+4044,int_stack+3984, 0.0, zero_stack, 1.0, int_stack+5364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7068,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3690,int_stack+4134,int_stack+4044, 0.0, zero_stack, 1.0, int_stack+5484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7188,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3960,int_stack+4320,int_stack+4260, 0.0, zero_stack, 1.0, int_stack+5790, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7068, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20274,int_stack+4410,int_stack+4320, 0.0, zero_stack, 1.0, int_stack+5910, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7188, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4140,int_stack+4596,int_stack+4536, 0.0, zero_stack, 1.0, int_stack+6216, 0.0, zero_stack, 1.0, int_stack+7068, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4320,int_stack+4686,int_stack+4596, 0.0, zero_stack, 1.0, int_stack+6336, 0.0, zero_stack, 1.0, int_stack+7188, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4590,int_stack+4872,int_stack+4812, 0.0, zero_stack, 1.0, int_stack+6642, 1.0, int_stack+7068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20544,int_stack+4962,int_stack+4872, 0.0, zero_stack, 1.0, int_stack+6762, 1.0, int_stack+7188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4770,int_stack+5148,int_stack+5088, 0.0, zero_stack, 2.0, int_stack+7068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20814,int_stack+5238,int_stack+5148, 0.0, zero_stack, 2.0, int_stack+7188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4950,int_stack+5574,int_stack+5424, 1.0, int_stack+5364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7554,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5130,int_stack+5664,int_stack+5574, 1.0, int_stack+5484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7674,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5400,int_stack+6000,int_stack+5850, 1.0, int_stack+5790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7554, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5580,int_stack+6090,int_stack+6000, 1.0, int_stack+5910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7674, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5850,int_stack+6426,int_stack+6276, 1.0, int_stack+6216, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7554, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6030,int_stack+6516,int_stack+6426, 1.0, int_stack+6336, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7674, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6300,int_stack+6852,int_stack+6702, 1.0, int_stack+6642, 0.0, zero_stack, 1.0, int_stack+7554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6480,int_stack+6942,int_stack+6852, 1.0, int_stack+6762, 0.0, zero_stack, 1.0, int_stack+7674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6750,int_stack+7278,int_stack+7128, 1.0, int_stack+7068, 1.0, int_stack+7554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21084,int_stack+7368,int_stack+7278, 1.0, int_stack+7188, 1.0, int_stack+7674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6930,int_stack+7764,int_stack+7614, 2.0, int_stack+7554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7110,int_stack+7854,int_stack+7764, 2.0, int_stack+7674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7380,int_stack+8040,int_stack+7980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9360,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7560,int_stack+8130,int_stack+8040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9480,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7830,int_stack+8316,int_stack+8256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9360, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8010,int_stack+8406,int_stack+8316, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9480, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8280,int_stack+8592,int_stack+8532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9360, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21354,int_stack+8682,int_stack+8592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9480, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8460,int_stack+8868,int_stack+8808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21624,int_stack+8958,int_stack+8868, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8640,int_stack+9144,int_stack+9084, 0.0, zero_stack, 1.0, int_stack+9360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8820,int_stack+9234,int_stack+9144, 0.0, zero_stack, 1.0, int_stack+9480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9090,int_stack+9570,int_stack+9420, 1.0, int_stack+9360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21894,int_stack+9660,int_stack+9570, 1.0, int_stack+9480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+9270,int_stack+9846,int_stack+9786,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9450,int_stack+9936,int_stack+9846,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9720,int_stack+10122,int_stack+10062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11442,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22164,int_stack+10212,int_stack+10122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11562,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9900,int_stack+10398,int_stack+10338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11442, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10080,int_stack+10488,int_stack+10398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11562, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10350,int_stack+10674,int_stack+10614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11442, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22434,int_stack+10764,int_stack+10674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11562, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10530,int_stack+10950,int_stack+10890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22704,int_stack+11040,int_stack+10950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10710,int_stack+11226,int_stack+11166, 0.0, zero_stack, 1.0, int_stack+11442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10890,int_stack+11316,int_stack+11226, 0.0, zero_stack, 1.0, int_stack+11562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11160,int_stack+11652,int_stack+11502, 1.0, int_stack+11442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22974,int_stack+11742,int_stack+11652, 1.0, int_stack+11562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11340,int_stack+11928,int_stack+11868,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11520,int_stack+12018,int_stack+11928,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11790,int_stack+12204,int_stack+12144,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23244,int_stack+12294,int_stack+12204,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11970,int_stack+12480,int_stack+12420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13800,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12150,int_stack+12570,int_stack+12480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13920,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12420,int_stack+12756,int_stack+12696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13800, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23514,int_stack+12846,int_stack+12756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13920, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12600,int_stack+13032,int_stack+12972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13800, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23784,int_stack+13122,int_stack+13032, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13920, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12780,int_stack+13308,int_stack+13248, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12960,int_stack+13398,int_stack+13308, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+13584,int_stack+13524, 0.0, zero_stack, 1.0, int_stack+13800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24054,int_stack+13674,int_stack+13584, 0.0, zero_stack, 1.0, int_stack+13920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13410,int_stack+14010,int_stack+13860, 1.0, int_stack+13800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13590,int_stack+14100,int_stack+14010, 1.0, int_stack+13920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13860,int_stack+14286,int_stack+14226,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24324,int_stack+14376,int_stack+14286,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14040,int_stack+14562,int_stack+14502,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14220,int_stack+14652,int_stack+14562,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14490,int_stack+14838,int_stack+14778,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24594,int_stack+14928,int_stack+14838,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14670,int_stack+15414,int_stack+15234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15054,6);
     Libderiv->ABCD[11] = int_stack + 14670;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24864,int_stack+15864,int_stack+15684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15054, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 24864;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25224,int_stack+16134,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15054, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 25224;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15864,int_stack+16404,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 15864;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16224,int_stack+16854,int_stack+16674, 0.0, zero_stack, 1.0, int_stack+15054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 16224;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25584,int_stack+17124,int_stack+360, 1.0, int_stack+15054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 25584;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+16854,int_stack+17394,int_stack+540,6);
     Libderiv->ABCD[2] = int_stack + 16854;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17214,int_stack+17664,int_stack+720,6);
     Libderiv->ABCD[1] = int_stack + 17214;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17574,int_stack+17934,int_stack+900,6);
     Libderiv->ABCD[0] = int_stack + 17574;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25944,int_stack+18384,int_stack+18204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15234,6);
     Libderiv->ABCD[155] = int_stack + 25944;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17934,int_stack+1260,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15234, 1.0, int_stack+15684,6);
     Libderiv->ABCD[143] = int_stack + 17934;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18294,int_stack+18654,int_stack+1530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15684, 0.0, zero_stack,6);
     Libderiv->ABCD[142] = int_stack + 18294;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1080,int_stack+18924,int_stack+1710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15234, 0.0, zero_stack, 1.0, int_stack+0,6);
     Libderiv->ABCD[131] = int_stack + 1080;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18654,int_stack+2070,int_stack+1890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15684, 1.0, int_stack+0, 0.0, zero_stack,6);
     Libderiv->ABCD[130] = int_stack + 18654;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1440,int_stack+19194,int_stack+2340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[129] = int_stack + 1440;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19014,int_stack+19464,int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15234, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180,6);
     Libderiv->ABCD[119] = int_stack + 19014;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19374,int_stack+2880,int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15684, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack,6);
     Libderiv->ABCD[118] = int_stack + 19374;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1800,int_stack+19734,int_stack+3150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[117] = int_stack + 1800;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2160,int_stack+20004,int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[116] = int_stack + 2160;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19734,int_stack+3690,int_stack+3510, 0.0, zero_stack, 1.0, int_stack+15234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16674,6);
     Libderiv->ABCD[107] = int_stack + 19734;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2520,int_stack+20274,int_stack+3960, 0.0, zero_stack, 1.0, int_stack+15684, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16674, 0.0, zero_stack,6);
     Libderiv->ABCD[106] = int_stack + 2520;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20094,int_stack+4320,int_stack+4140, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+16674, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[105] = int_stack + 20094;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2880,int_stack+20544,int_stack+4590, 0.0, zero_stack, 1.0, int_stack+180, 1.0, int_stack+16674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[104] = int_stack + 2880;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20454,int_stack+20814,int_stack+4770, 0.0, zero_stack, 2.0, int_stack+16674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[103] = int_stack + 20454;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3240,int_stack+5130,int_stack+4950, 1.0, int_stack+15234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360,6);
     Libderiv->ABCD[95] = int_stack + 3240;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3600,int_stack+5580,int_stack+5400, 1.0, int_stack+15684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack,6);
     Libderiv->ABCD[94] = int_stack + 3600;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3960,int_stack+6030,int_stack+5850, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[93] = int_stack + 3960;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4320,int_stack+6480,int_stack+6300, 1.0, int_stack+180, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[92] = int_stack + 4320;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+21084,int_stack+6750, 1.0, int_stack+16674, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[91] = int_stack + 0;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20814,int_stack+7110,int_stack+6930, 2.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[90] = int_stack + 20814;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4680,int_stack+7560,int_stack+7380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540,6);
     Libderiv->ABCD[47] = int_stack + 4680;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5040,int_stack+8010,int_stack+7830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack,6);
     Libderiv->ABCD[46] = int_stack + 5040;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5400,int_stack+21354,int_stack+8280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[45] = int_stack + 5400;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21174,int_stack+21624,int_stack+8460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[44] = int_stack + 21174;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21534,int_stack+8820,int_stack+8640, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[43] = int_stack + 21534;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5760,int_stack+21894,int_stack+9090, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[42] = int_stack + 5760;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+360,int_stack+9450,int_stack+9270,6);
     Libderiv->ABCD[38] = int_stack + 360;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6120,int_stack+22164,int_stack+9720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720,6);
     Libderiv->ABCD[35] = int_stack + 6120;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21894,int_stack+10080,int_stack+9900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack,6);
     Libderiv->ABCD[34] = int_stack + 21894;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6480,int_stack+22434,int_stack+10350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[33] = int_stack + 6480;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22254,int_stack+22704,int_stack+10530, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[32] = int_stack + 22254;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22614,int_stack+10890,int_stack+10710, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[31] = int_stack + 22614;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6840,int_stack+22974,int_stack+11160, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[30] = int_stack + 6840;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7200,int_stack+11520,int_stack+11340,6);
     Libderiv->ABCD[26] = int_stack + 7200;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7560,int_stack+23244,int_stack+11790,6);
     Libderiv->ABCD[25] = int_stack + 7560;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22974,int_stack+12150,int_stack+11970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900,6);
     Libderiv->ABCD[23] = int_stack + 22974;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7920,int_stack+23514,int_stack+12420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack,6);
     Libderiv->ABCD[22] = int_stack + 7920;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23334,int_stack+23784,int_stack+12600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[21] = int_stack + 23334;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23694,int_stack+12960,int_stack+12780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[20] = int_stack + 23694;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8280,int_stack+24054,int_stack+13230, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[19] = int_stack + 8280;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8640,int_stack+13590,int_stack+13410, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[18] = int_stack + 8640;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+720,int_stack+24324,int_stack+13860,6);
     Libderiv->ABCD[14] = int_stack + 720;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+24054,int_stack+14220,int_stack+14040,6);
     Libderiv->ABCD[13] = int_stack + 24054;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+9000,int_stack+24594,int_stack+14490,6);
     Libderiv->ABCD[12] = int_stack + 9000;

}
