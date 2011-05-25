#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_f0dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|dd) integrals */

void d12hrr_order_f0dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][10] = int_stack + 150;
 Libderiv->deriv_classes[3][4][9] = int_stack + 300;
 Libderiv->deriv_classes[3][4][8] = int_stack + 450;
 Libderiv->deriv_classes[3][4][7] = int_stack + 600;
 Libderiv->dvrr_classes[3][3] = int_stack + 750;
 Libderiv->deriv_classes[3][4][6] = int_stack + 850;
 Libderiv->deriv_classes[3][4][2] = int_stack + 1000;
 Libderiv->deriv_classes[3][4][1] = int_stack + 1150;
 Libderiv->deriv_classes[3][4][0] = int_stack + 1300;
 Libderiv->deriv2_classes[3][2][143] = int_stack + 1450;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 1510;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 1610;
 Libderiv->deriv2_classes[3][2][131] = int_stack + 1760;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 1820;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 1920;
 Libderiv->deriv2_classes[3][2][130] = int_stack + 2070;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 2130;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 2230;
 Libderiv->deriv2_classes[3][2][119] = int_stack + 2380;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 2440;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 2540;
 Libderiv->deriv2_classes[3][2][118] = int_stack + 2690;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 2750;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 2850;
 Libderiv->deriv2_classes[3][2][117] = int_stack + 3000;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 3060;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 3160;
 Libderiv->deriv2_classes[3][2][107] = int_stack + 3310;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 3370;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 3470;
 Libderiv->deriv2_classes[3][2][106] = int_stack + 3620;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 3680;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 3780;
 Libderiv->deriv2_classes[3][2][105] = int_stack + 3930;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 3990;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 4090;
 Libderiv->deriv2_classes[3][2][104] = int_stack + 4240;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 4300;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 4400;
 Libderiv->deriv2_classes[3][2][95] = int_stack + 4550;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 4610;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 4710;
 Libderiv->deriv2_classes[3][2][94] = int_stack + 4860;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 4920;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 5020;
 Libderiv->deriv2_classes[3][2][93] = int_stack + 5170;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 5230;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 5330;
 Libderiv->deriv2_classes[3][2][92] = int_stack + 5480;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 5540;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 5640;
 Libderiv->deriv2_classes[3][2][91] = int_stack + 5790;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 5850;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 5950;
 Libderiv->deriv_classes[3][2][11] = int_stack + 6100;
 Libderiv->deriv2_classes[3][2][83] = int_stack + 6160;
 Libderiv->deriv_classes[3][3][11] = int_stack + 6220;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 6320;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 6420;
 Libderiv->deriv_classes[3][2][10] = int_stack + 6570;
 Libderiv->deriv2_classes[3][2][82] = int_stack + 6630;
 Libderiv->deriv_classes[3][3][10] = int_stack + 6690;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 6790;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 6890;
 Libderiv->deriv_classes[3][2][9] = int_stack + 7040;
 Libderiv->deriv2_classes[3][2][81] = int_stack + 7100;
 Libderiv->deriv_classes[3][3][9] = int_stack + 7160;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 7260;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 7360;
 Libderiv->deriv_classes[3][2][8] = int_stack + 7510;
 Libderiv->deriv2_classes[3][2][80] = int_stack + 7570;
 Libderiv->deriv_classes[3][3][8] = int_stack + 7630;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 7730;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 7830;
 Libderiv->deriv_classes[3][2][7] = int_stack + 7980;
 Libderiv->deriv2_classes[3][2][79] = int_stack + 8040;
 Libderiv->deriv_classes[3][3][7] = int_stack + 8100;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 8200;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 8300;
 Libderiv->dvrr_classes[3][2] = int_stack + 8450;
 Libderiv->deriv_classes[3][2][6] = int_stack + 8510;
 Libderiv->deriv2_classes[3][2][78] = int_stack + 8570;
 Libderiv->deriv_classes[3][3][6] = int_stack + 8630;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 8730;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 8830;
 Libderiv->deriv2_classes[3][2][35] = int_stack + 8980;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 9040;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 9140;
 Libderiv->deriv2_classes[3][2][34] = int_stack + 9290;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 9350;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 9450;
 Libderiv->deriv2_classes[3][2][33] = int_stack + 9600;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 9660;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 9760;
 Libderiv->deriv2_classes[3][2][32] = int_stack + 9910;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 9970;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 10070;
 Libderiv->deriv2_classes[3][2][31] = int_stack + 10220;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 10280;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 10380;
 Libderiv->deriv_classes[3][2][2] = int_stack + 10530;
 Libderiv->deriv2_classes[3][2][30] = int_stack + 10590;
 Libderiv->deriv_classes[3][3][2] = int_stack + 10650;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 10750;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 10850;
 Libderiv->deriv2_classes[3][2][26] = int_stack + 11000;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 11060;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 11160;
 Libderiv->deriv2_classes[3][2][23] = int_stack + 11310;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 11370;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 11470;
 Libderiv->deriv2_classes[3][2][22] = int_stack + 11620;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 11680;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 11780;
 Libderiv->deriv2_classes[3][2][21] = int_stack + 11930;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 11990;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 12090;
 Libderiv->deriv2_classes[3][2][20] = int_stack + 12240;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 12300;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 12400;
 Libderiv->deriv2_classes[3][2][19] = int_stack + 12550;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 12610;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 12710;
 Libderiv->deriv_classes[3][2][1] = int_stack + 12860;
 Libderiv->deriv2_classes[3][2][18] = int_stack + 12920;
 Libderiv->deriv_classes[3][3][1] = int_stack + 12980;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 13080;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 13180;
 Libderiv->deriv2_classes[3][2][14] = int_stack + 13330;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 13390;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 13490;
 Libderiv->deriv2_classes[3][2][13] = int_stack + 13640;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 13700;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 13800;
 Libderiv->deriv2_classes[3][2][11] = int_stack + 13950;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 14010;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 14110;
 Libderiv->deriv2_classes[3][2][10] = int_stack + 14260;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 14320;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 14420;
 Libderiv->deriv2_classes[3][2][9] = int_stack + 14570;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 14630;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 14730;
 Libderiv->deriv2_classes[3][2][8] = int_stack + 14880;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 14940;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 15040;
 Libderiv->deriv2_classes[3][2][7] = int_stack + 15190;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 15250;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 15350;
 Libderiv->deriv_classes[3][2][0] = int_stack + 15500;
 Libderiv->deriv2_classes[3][2][6] = int_stack + 15560;
 Libderiv->deriv_classes[3][3][0] = int_stack + 15620;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 15720;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 15820;
 Libderiv->deriv2_classes[3][2][2] = int_stack + 15970;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 16030;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 16130;
 Libderiv->deriv2_classes[3][2][1] = int_stack + 16280;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 16340;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 16440;
 Libderiv->deriv2_classes[3][2][0] = int_stack + 16590;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 16650;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 16750;
 memset(int_stack,0,135200);

 Libderiv->dvrr_stack = int_stack + 27520;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_f0dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+16900,int_stack+750,int_stack+8450,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17080,int_stack+6220,int_stack+6100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8450,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17260,int_stack+0,int_stack+6220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17560,int_stack+6690,int_stack+6570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8450, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17740,int_stack+150,int_stack+6690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+7160,int_stack+7040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8450, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18040,int_stack+300,int_stack+7160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+180,int_stack+7630,int_stack+7510, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18340,int_stack+450,int_stack+7630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+360,int_stack+8100,int_stack+7980, 0.0, zero_stack, 1.0, int_stack+8450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18640,int_stack+600,int_stack+8100, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+8630,int_stack+8510, 1.0, int_stack+8450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18940,int_stack+850,int_stack+8630, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+720,int_stack+10650,int_stack+10530,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19240,int_stack+1000,int_stack+10650,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+900,int_stack+12980,int_stack+12860,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19540,int_stack+1150,int_stack+12980,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+15620,int_stack+15500,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19840,int_stack+1300,int_stack+15620,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1260,int_stack+1510,int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6100,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20140,int_stack+1610,int_stack+1510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6220,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1440,int_stack+1820,int_stack+1760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6100, 1.0, int_stack+6570,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20440,int_stack+1920,int_stack+1820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6220, 1.0, int_stack+6690,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1620,int_stack+2130,int_stack+2070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6570, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2230,int_stack+2130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6690, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2100,int_stack+2440,int_stack+2380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6100, 0.0, zero_stack, 1.0, int_stack+7040,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20740,int_stack+2540,int_stack+2440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6220, 0.0, zero_stack, 1.0, int_stack+7160,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2280,int_stack+2750,int_stack+2690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6570, 1.0, int_stack+7040, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21040,int_stack+2850,int_stack+2750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6690, 1.0, int_stack+7160, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2460,int_stack+3060,int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7040, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2640,int_stack+3160,int_stack+3060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7160, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2940,int_stack+3370,int_stack+3310, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7510,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21340,int_stack+3470,int_stack+3370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7630,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3120,int_stack+3680,int_stack+3620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6570, 0.0, zero_stack, 1.0, int_stack+7510, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3300,int_stack+3780,int_stack+3680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6690, 0.0, zero_stack, 1.0, int_stack+7630, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3600,int_stack+3990,int_stack+3930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7040, 1.0, int_stack+7510, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21640,int_stack+4090,int_stack+3990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7160, 1.0, int_stack+7630, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3780,int_stack+4300,int_stack+4240, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3960,int_stack+4400,int_stack+4300, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4260,int_stack+4610,int_stack+4550, 0.0, zero_stack, 1.0, int_stack+6100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7980,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21940,int_stack+4710,int_stack+4610, 0.0, zero_stack, 1.0, int_stack+6220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8100,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4440,int_stack+4920,int_stack+4860, 0.0, zero_stack, 1.0, int_stack+6570, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7980, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4620,int_stack+5020,int_stack+4920, 0.0, zero_stack, 1.0, int_stack+6690, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8100, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4920,int_stack+5230,int_stack+5170, 0.0, zero_stack, 1.0, int_stack+7040, 0.0, zero_stack, 1.0, int_stack+7980, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22240,int_stack+5330,int_stack+5230, 0.0, zero_stack, 1.0, int_stack+7160, 0.0, zero_stack, 1.0, int_stack+8100, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5100,int_stack+5540,int_stack+5480, 0.0, zero_stack, 1.0, int_stack+7510, 1.0, int_stack+7980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22540,int_stack+5640,int_stack+5540, 0.0, zero_stack, 1.0, int_stack+7630, 1.0, int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5280,int_stack+5850,int_stack+5790, 0.0, zero_stack, 2.0, int_stack+7980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5460,int_stack+5950,int_stack+5850, 0.0, zero_stack, 2.0, int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5760,int_stack+6320,int_stack+6160, 1.0, int_stack+6100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8510,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22840,int_stack+6420,int_stack+6320, 1.0, int_stack+6220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8630,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5940,int_stack+6790,int_stack+6630, 1.0, int_stack+6570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8510, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6120,int_stack+6890,int_stack+6790, 1.0, int_stack+6690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8630, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6420,int_stack+7260,int_stack+7100, 1.0, int_stack+7040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8510, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6600,int_stack+7360,int_stack+7260, 1.0, int_stack+7160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8630, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6900,int_stack+7730,int_stack+7570, 1.0, int_stack+7510, 0.0, zero_stack, 1.0, int_stack+8510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7080,int_stack+7830,int_stack+7730, 1.0, int_stack+7630, 0.0, zero_stack, 1.0, int_stack+8630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7380,int_stack+8200,int_stack+8040, 1.0, int_stack+7980, 1.0, int_stack+8510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7560,int_stack+8300,int_stack+8200, 1.0, int_stack+8100, 1.0, int_stack+8630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7860,int_stack+8730,int_stack+8570, 2.0, int_stack+8510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8040,int_stack+8830,int_stack+8730, 2.0, int_stack+8630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8340,int_stack+9040,int_stack+8980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10530,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8520,int_stack+9140,int_stack+9040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10650,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8820,int_stack+9350,int_stack+9290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10530, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9000,int_stack+9450,int_stack+9350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10650, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9300,int_stack+9660,int_stack+9600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10530, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23140,int_stack+9760,int_stack+9660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10650, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9480,int_stack+9970,int_stack+9910, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9660,int_stack+10070,int_stack+9970, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9960,int_stack+10280,int_stack+10220, 0.0, zero_stack, 1.0, int_stack+10530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23440,int_stack+10380,int_stack+10280, 0.0, zero_stack, 1.0, int_stack+10650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10140,int_stack+10750,int_stack+10590, 1.0, int_stack+10530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10320,int_stack+10850,int_stack+10750, 1.0, int_stack+10650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+10620,int_stack+11060,int_stack+11000,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+23740,int_stack+11160,int_stack+11060,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10800,int_stack+11370,int_stack+11310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12860,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10980,int_stack+11470,int_stack+11370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12980,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11280,int_stack+11680,int_stack+11620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12860, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24040,int_stack+11780,int_stack+11680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12980, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11460,int_stack+11990,int_stack+11930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12860, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11640,int_stack+12090,int_stack+11990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12980, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11940,int_stack+12300,int_stack+12240, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24340,int_stack+12400,int_stack+12300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12120,int_stack+12610,int_stack+12550, 0.0, zero_stack, 1.0, int_stack+12860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12300,int_stack+12710,int_stack+12610, 0.0, zero_stack, 1.0, int_stack+12980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12600,int_stack+13080,int_stack+12920, 1.0, int_stack+12860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24640,int_stack+13180,int_stack+13080, 1.0, int_stack+12980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12780,int_stack+13390,int_stack+13330,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12960,int_stack+13490,int_stack+13390,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+13260,int_stack+13700,int_stack+13640,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24940,int_stack+13800,int_stack+13700,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+13440,int_stack+14010,int_stack+13950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15500,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13620,int_stack+14110,int_stack+14010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15620,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+13920,int_stack+14320,int_stack+14260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15500, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25240,int_stack+14420,int_stack+14320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15620, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14100,int_stack+14630,int_stack+14570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14280,int_stack+14730,int_stack+14630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15620, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14580,int_stack+14940,int_stack+14880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25540,int_stack+15040,int_stack+14940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14760,int_stack+15250,int_stack+15190, 0.0, zero_stack, 1.0, int_stack+15500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14940,int_stack+15350,int_stack+15250, 0.0, zero_stack, 1.0, int_stack+15620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15240,int_stack+15720,int_stack+15560, 1.0, int_stack+15500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25840,int_stack+15820,int_stack+15720, 1.0, int_stack+15620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15420,int_stack+16030,int_stack+15970,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15600,int_stack+16130,int_stack+16030,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15900,int_stack+16340,int_stack+16280,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+26140,int_stack+16440,int_stack+16340,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+16080,int_stack+16650,int_stack+16590,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16260,int_stack+16750,int_stack+16650,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+26440,int_stack+17260,int_stack+17080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16900,10);
     Libderiv->ABCD[11] = int_stack + 26440;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+26800,int_stack+17740,int_stack+17560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16900, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 26800;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+27160,int_stack+18040,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16900, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 27160;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17740,int_stack+18340,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 17740;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+18100,int_stack+18640,int_stack+360, 0.0, zero_stack, 1.0, int_stack+16900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 18100;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+18460,int_stack+18940,int_stack+540, 1.0, int_stack+16900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 18460;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+18820,int_stack+19240,int_stack+720,10);
     Libderiv->ABCD[2] = int_stack + 18820;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+19180,int_stack+19540,int_stack+900,10);
     Libderiv->ABCD[1] = int_stack + 19180;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+16560,int_stack+19840,int_stack+1080,10);
     Libderiv->ABCD[0] = int_stack + 16560;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19540,int_stack+20140,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17080,10);
     Libderiv->ABCD[155] = int_stack + 19540;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19900,int_stack+20440,int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17080, 1.0, int_stack+17560,10);
     Libderiv->ABCD[143] = int_stack + 19900;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1260,int_stack+1800,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17560, 0.0, zero_stack,10);
     Libderiv->ABCD[142] = int_stack + 1260;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1620,int_stack+20740,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17080, 0.0, zero_stack, 1.0, int_stack+0,10);
     Libderiv->ABCD[131] = int_stack + 1620;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20260,int_stack+21040,int_stack+2280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17560, 1.0, int_stack+0, 0.0, zero_stack,10);
     Libderiv->ABCD[130] = int_stack + 20260;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1980,int_stack+2640,int_stack+2460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[129] = int_stack + 1980;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2340,int_stack+21340,int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180,10);
     Libderiv->ABCD[119] = int_stack + 2340;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2700,int_stack+3300,int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17560, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack,10);
     Libderiv->ABCD[118] = int_stack + 2700;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3060,int_stack+21640,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[117] = int_stack + 3060;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3420,int_stack+3960,int_stack+3780, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[116] = int_stack + 3420;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3780,int_stack+21940,int_stack+4260, 0.0, zero_stack, 1.0, int_stack+17080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360,10);
     Libderiv->ABCD[107] = int_stack + 3780;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20620,int_stack+4620,int_stack+4440, 0.0, zero_stack, 1.0, int_stack+17560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack,10);
     Libderiv->ABCD[106] = int_stack + 20620;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4140,int_stack+22240,int_stack+4920, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[105] = int_stack + 4140;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4500,int_stack+22540,int_stack+5100, 0.0, zero_stack, 1.0, int_stack+180, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[104] = int_stack + 4500;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4860,int_stack+5460,int_stack+5280, 0.0, zero_stack, 2.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[103] = int_stack + 4860;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5220,int_stack+22840,int_stack+5760, 1.0, int_stack+17080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540,10);
     Libderiv->ABCD[95] = int_stack + 5220;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5580,int_stack+6120,int_stack+5940, 1.0, int_stack+17560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack,10);
     Libderiv->ABCD[94] = int_stack + 5580;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5940,int_stack+6600,int_stack+6420, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[93] = int_stack + 5940;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6300,int_stack+7080,int_stack+6900, 1.0, int_stack+180, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[92] = int_stack + 6300;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+7560,int_stack+7380, 1.0, int_stack+360, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[91] = int_stack + 0;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6660,int_stack+8040,int_stack+7860, 2.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[90] = int_stack + 6660;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+360,int_stack+8520,int_stack+8340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720,10);
     Libderiv->ABCD[47] = int_stack + 360;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7020,int_stack+9000,int_stack+8820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack,10);
     Libderiv->ABCD[46] = int_stack + 7020;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7380,int_stack+23140,int_stack+9300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[45] = int_stack + 7380;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7740,int_stack+9660,int_stack+9480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[44] = int_stack + 7740;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8100,int_stack+23440,int_stack+9960, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[43] = int_stack + 8100;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8460,int_stack+10320,int_stack+10140, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[42] = int_stack + 8460;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+8820,int_stack+23740,int_stack+10620,10);
     Libderiv->ABCD[38] = int_stack + 8820;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9180,int_stack+10980,int_stack+10800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900,10);
     Libderiv->ABCD[35] = int_stack + 9180;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9540,int_stack+24040,int_stack+11280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack,10);
     Libderiv->ABCD[34] = int_stack + 9540;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9900,int_stack+11640,int_stack+11460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[33] = int_stack + 9900;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10260,int_stack+24340,int_stack+11940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[32] = int_stack + 10260;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10620,int_stack+12300,int_stack+12120, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[31] = int_stack + 10620;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10980,int_stack+24640,int_stack+12600, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[30] = int_stack + 10980;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+720,int_stack+12960,int_stack+12780,10);
     Libderiv->ABCD[26] = int_stack + 720;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+11340,int_stack+24940,int_stack+13260,10);
     Libderiv->ABCD[25] = int_stack + 11340;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11700,int_stack+13620,int_stack+13440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080,10);
     Libderiv->ABCD[23] = int_stack + 11700;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12060,int_stack+25240,int_stack+13920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack,10);
     Libderiv->ABCD[22] = int_stack + 12060;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12420,int_stack+14280,int_stack+14100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[21] = int_stack + 12420;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12780,int_stack+25540,int_stack+14580, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[20] = int_stack + 12780;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+13140,int_stack+14940,int_stack+14760, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[19] = int_stack + 13140;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+13500,int_stack+25840,int_stack+15240, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[18] = int_stack + 13500;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+13860,int_stack+15600,int_stack+15420,10);
     Libderiv->ABCD[14] = int_stack + 13860;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+14220,int_stack+26140,int_stack+15900,10);
     Libderiv->ABCD[13] = int_stack + 14220;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+14580,int_stack+16260,int_stack+16080,10);
     Libderiv->ABCD[12] = int_stack + 14580;

}
