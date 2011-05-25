#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00dp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|dp) integrals */

void d12hrr_order_00dp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][3][11] = int_stack + 0;
 Libderiv->deriv_classes[0][3][10] = int_stack + 10;
 Libderiv->deriv_classes[0][3][9] = int_stack + 20;
 Libderiv->deriv_classes[0][3][8] = int_stack + 30;
 Libderiv->deriv_classes[0][3][7] = int_stack + 40;
 Libderiv->dvrr_classes[0][2] = int_stack + 50;
 Libderiv->deriv_classes[0][3][6] = int_stack + 56;
 Libderiv->deriv_classes[0][3][2] = int_stack + 66;
 Libderiv->deriv_classes[0][3][1] = int_stack + 76;
 Libderiv->deriv_classes[0][3][0] = int_stack + 86;
 Libderiv->deriv2_classes[0][2][143] = int_stack + 96;
 Libderiv->deriv2_classes[0][3][143] = int_stack + 102;
 Libderiv->deriv2_classes[0][2][131] = int_stack + 112;
 Libderiv->deriv2_classes[0][3][131] = int_stack + 118;
 Libderiv->deriv2_classes[0][2][130] = int_stack + 128;
 Libderiv->deriv2_classes[0][3][130] = int_stack + 134;
 Libderiv->deriv2_classes[0][2][119] = int_stack + 144;
 Libderiv->deriv2_classes[0][3][119] = int_stack + 150;
 Libderiv->deriv2_classes[0][2][118] = int_stack + 160;
 Libderiv->deriv2_classes[0][3][118] = int_stack + 166;
 Libderiv->deriv2_classes[0][2][117] = int_stack + 176;
 Libderiv->deriv2_classes[0][3][117] = int_stack + 182;
 Libderiv->deriv2_classes[0][2][107] = int_stack + 192;
 Libderiv->deriv2_classes[0][3][107] = int_stack + 198;
 Libderiv->deriv2_classes[0][2][106] = int_stack + 208;
 Libderiv->deriv2_classes[0][3][106] = int_stack + 214;
 Libderiv->deriv2_classes[0][2][105] = int_stack + 224;
 Libderiv->deriv2_classes[0][3][105] = int_stack + 230;
 Libderiv->deriv2_classes[0][2][104] = int_stack + 240;
 Libderiv->deriv2_classes[0][3][104] = int_stack + 246;
 Libderiv->deriv2_classes[0][2][95] = int_stack + 256;
 Libderiv->deriv2_classes[0][3][95] = int_stack + 262;
 Libderiv->deriv2_classes[0][2][94] = int_stack + 272;
 Libderiv->deriv2_classes[0][3][94] = int_stack + 278;
 Libderiv->deriv2_classes[0][2][93] = int_stack + 288;
 Libderiv->deriv2_classes[0][3][93] = int_stack + 294;
 Libderiv->deriv2_classes[0][2][92] = int_stack + 304;
 Libderiv->deriv2_classes[0][3][92] = int_stack + 310;
 Libderiv->deriv2_classes[0][2][91] = int_stack + 320;
 Libderiv->deriv2_classes[0][3][91] = int_stack + 326;
 Libderiv->deriv_classes[0][2][11] = int_stack + 336;
 Libderiv->deriv2_classes[0][2][83] = int_stack + 342;
 Libderiv->deriv2_classes[0][3][83] = int_stack + 348;
 Libderiv->deriv_classes[0][2][10] = int_stack + 358;
 Libderiv->deriv2_classes[0][2][82] = int_stack + 364;
 Libderiv->deriv2_classes[0][3][82] = int_stack + 370;
 Libderiv->deriv_classes[0][2][9] = int_stack + 380;
 Libderiv->deriv2_classes[0][2][81] = int_stack + 386;
 Libderiv->deriv2_classes[0][3][81] = int_stack + 392;
 Libderiv->deriv_classes[0][2][8] = int_stack + 402;
 Libderiv->deriv2_classes[0][2][80] = int_stack + 408;
 Libderiv->deriv2_classes[0][3][80] = int_stack + 414;
 Libderiv->deriv_classes[0][2][7] = int_stack + 424;
 Libderiv->deriv2_classes[0][2][79] = int_stack + 430;
 Libderiv->deriv2_classes[0][3][79] = int_stack + 436;
 Libderiv->deriv_classes[0][2][6] = int_stack + 446;
 Libderiv->deriv2_classes[0][2][78] = int_stack + 452;
 Libderiv->deriv2_classes[0][3][78] = int_stack + 458;
 Libderiv->deriv2_classes[0][2][35] = int_stack + 468;
 Libderiv->deriv2_classes[0][3][35] = int_stack + 474;
 Libderiv->deriv2_classes[0][2][34] = int_stack + 484;
 Libderiv->deriv2_classes[0][3][34] = int_stack + 490;
 Libderiv->deriv2_classes[0][2][33] = int_stack + 500;
 Libderiv->deriv2_classes[0][3][33] = int_stack + 506;
 Libderiv->deriv2_classes[0][2][32] = int_stack + 516;
 Libderiv->deriv2_classes[0][3][32] = int_stack + 522;
 Libderiv->deriv2_classes[0][2][31] = int_stack + 532;
 Libderiv->deriv2_classes[0][3][31] = int_stack + 538;
 Libderiv->deriv_classes[0][2][2] = int_stack + 548;
 Libderiv->deriv2_classes[0][2][30] = int_stack + 554;
 Libderiv->deriv2_classes[0][3][30] = int_stack + 560;
 Libderiv->deriv2_classes[0][2][26] = int_stack + 570;
 Libderiv->deriv2_classes[0][3][26] = int_stack + 576;
 Libderiv->deriv2_classes[0][2][23] = int_stack + 586;
 Libderiv->deriv2_classes[0][3][23] = int_stack + 592;
 Libderiv->deriv2_classes[0][2][22] = int_stack + 602;
 Libderiv->deriv2_classes[0][3][22] = int_stack + 608;
 Libderiv->deriv2_classes[0][2][21] = int_stack + 618;
 Libderiv->deriv2_classes[0][3][21] = int_stack + 624;
 Libderiv->deriv2_classes[0][2][20] = int_stack + 634;
 Libderiv->deriv2_classes[0][3][20] = int_stack + 640;
 Libderiv->deriv2_classes[0][2][19] = int_stack + 650;
 Libderiv->deriv2_classes[0][3][19] = int_stack + 656;
 Libderiv->deriv_classes[0][2][1] = int_stack + 666;
 Libderiv->deriv2_classes[0][2][18] = int_stack + 672;
 Libderiv->deriv2_classes[0][3][18] = int_stack + 678;
 Libderiv->deriv2_classes[0][2][14] = int_stack + 688;
 Libderiv->deriv2_classes[0][3][14] = int_stack + 694;
 Libderiv->deriv2_classes[0][2][13] = int_stack + 704;
 Libderiv->deriv2_classes[0][3][13] = int_stack + 710;
 Libderiv->deriv2_classes[0][2][11] = int_stack + 720;
 Libderiv->deriv2_classes[0][3][11] = int_stack + 726;
 Libderiv->deriv2_classes[0][2][10] = int_stack + 736;
 Libderiv->deriv2_classes[0][3][10] = int_stack + 742;
 Libderiv->deriv2_classes[0][2][9] = int_stack + 752;
 Libderiv->deriv2_classes[0][3][9] = int_stack + 758;
 Libderiv->deriv2_classes[0][2][8] = int_stack + 768;
 Libderiv->deriv2_classes[0][3][8] = int_stack + 774;
 Libderiv->deriv2_classes[0][2][7] = int_stack + 784;
 Libderiv->deriv2_classes[0][3][7] = int_stack + 790;
 Libderiv->deriv_classes[0][2][0] = int_stack + 800;
 Libderiv->deriv2_classes[0][2][6] = int_stack + 806;
 Libderiv->deriv2_classes[0][3][6] = int_stack + 812;
 Libderiv->deriv2_classes[0][2][2] = int_stack + 822;
 Libderiv->deriv2_classes[0][3][2] = int_stack + 828;
 Libderiv->deriv2_classes[0][2][1] = int_stack + 838;
 Libderiv->deriv2_classes[0][3][1] = int_stack + 844;
 Libderiv->deriv2_classes[0][2][0] = int_stack + 854;
 Libderiv->deriv2_classes[0][3][0] = int_stack + 860;
 memset(int_stack,0,6960);

 Libderiv->dvrr_stack = int_stack + 996;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00dp(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+870,int_stack+0,int_stack+336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50,1);
     Libderiv->ABCD[11] = int_stack + 870;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+888,int_stack+10,int_stack+358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 888;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+20,int_stack+380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+906,int_stack+30,int_stack+402, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 906;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18,int_stack+40,int_stack+424, 0.0, zero_stack, 1.0, int_stack+50, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 18;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+924,int_stack+56,int_stack+446, 1.0, int_stack+50, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 924;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+36,int_stack+66,int_stack+548,1);
     Libderiv->ABCD[2] = int_stack + 36;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+76,int_stack+666,1);
     Libderiv->ABCD[1] = int_stack + 54;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+942,int_stack+86,int_stack+800,1);
     Libderiv->ABCD[0] = int_stack + 942;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72,int_stack+102,int_stack+96, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+336,1);
     Libderiv->ABCD[155] = int_stack + 72;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+90,int_stack+118,int_stack+112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+336, 1.0, int_stack+358,1);
     Libderiv->ABCD[143] = int_stack + 90;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+134,int_stack+128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+358, 0.0, zero_stack,1);
     Libderiv->ABCD[142] = int_stack + 108;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+126,int_stack+150,int_stack+144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+336, 0.0, zero_stack, 1.0, int_stack+380,1);
     Libderiv->ABCD[131] = int_stack + 126;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+960,int_stack+166,int_stack+160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+358, 1.0, int_stack+380, 0.0, zero_stack,1);
     Libderiv->ABCD[130] = int_stack + 960;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+144,int_stack+182,int_stack+176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+380, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[129] = int_stack + 144;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+162,int_stack+198,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+336, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+402,1);
     Libderiv->ABCD[119] = int_stack + 162;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+180,int_stack+214,int_stack+208, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+358, 0.0, zero_stack, 1.0, int_stack+402, 0.0, zero_stack,1);
     Libderiv->ABCD[118] = int_stack + 180;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+198,int_stack+230,int_stack+224, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+380, 1.0, int_stack+402, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[117] = int_stack + 198;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+246,int_stack+240, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+402, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[116] = int_stack + 216;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+234,int_stack+262,int_stack+256, 0.0, zero_stack, 1.0, int_stack+336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+424,1);
     Libderiv->ABCD[107] = int_stack + 234;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+252,int_stack+278,int_stack+272, 0.0, zero_stack, 1.0, int_stack+358, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+424, 0.0, zero_stack,1);
     Libderiv->ABCD[106] = int_stack + 252;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+294,int_stack+288, 0.0, zero_stack, 1.0, int_stack+380, 0.0, zero_stack, 1.0, int_stack+424, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[105] = int_stack + 270;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+978,int_stack+310,int_stack+304, 0.0, zero_stack, 1.0, int_stack+402, 1.0, int_stack+424, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[104] = int_stack + 978;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+288,int_stack+326,int_stack+320, 0.0, zero_stack, 2.0, int_stack+424, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[103] = int_stack + 288;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+306,int_stack+348,int_stack+342, 1.0, int_stack+336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+446,1);
     Libderiv->ABCD[95] = int_stack + 306;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+370,int_stack+364, 1.0, int_stack+358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+446, 0.0, zero_stack,1);
     Libderiv->ABCD[94] = int_stack + 324;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+342,int_stack+392,int_stack+386, 1.0, int_stack+380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+446, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[93] = int_stack + 342;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+360,int_stack+414,int_stack+408, 1.0, int_stack+402, 0.0, zero_stack, 1.0, int_stack+446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[92] = int_stack + 360;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+378,int_stack+436,int_stack+430, 1.0, int_stack+424, 1.0, int_stack+446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[91] = int_stack + 378;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+396,int_stack+458,int_stack+452, 2.0, int_stack+446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[90] = int_stack + 396;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+414,int_stack+474,int_stack+468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+548,1);
     Libderiv->ABCD[47] = int_stack + 414;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+490,int_stack+484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+548, 0.0, zero_stack,1);
     Libderiv->ABCD[46] = int_stack + 432;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+450,int_stack+506,int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+548, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[45] = int_stack + 450;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+468,int_stack+522,int_stack+516, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[44] = int_stack + 468;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+486,int_stack+538,int_stack+532, 0.0, zero_stack, 1.0, int_stack+548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[43] = int_stack + 486;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+504,int_stack+560,int_stack+554, 1.0, int_stack+548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[42] = int_stack + 504;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+522,int_stack+576,int_stack+570,1);
     Libderiv->ABCD[38] = int_stack + 522;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+592,int_stack+586, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+666,1);
     Libderiv->ABCD[35] = int_stack + 540;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+558,int_stack+608,int_stack+602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+666, 0.0, zero_stack,1);
     Libderiv->ABCD[34] = int_stack + 558;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+576,int_stack+624,int_stack+618, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+666, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[33] = int_stack + 576;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+594,int_stack+640,int_stack+634, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[32] = int_stack + 594;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+612,int_stack+656,int_stack+650, 0.0, zero_stack, 1.0, int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[31] = int_stack + 612;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+630,int_stack+678,int_stack+672, 1.0, int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[30] = int_stack + 630;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+694,int_stack+688,1);
     Libderiv->ABCD[26] = int_stack + 648;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+666,int_stack+710,int_stack+704,1);
     Libderiv->ABCD[25] = int_stack + 666;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+684,int_stack+726,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+800,1);
     Libderiv->ABCD[23] = int_stack + 684;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+702,int_stack+742,int_stack+736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+800, 0.0, zero_stack,1);
     Libderiv->ABCD[22] = int_stack + 702;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+720,int_stack+758,int_stack+752, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+800, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[21] = int_stack + 720;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+738,int_stack+774,int_stack+768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[20] = int_stack + 738;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+756,int_stack+790,int_stack+784, 0.0, zero_stack, 1.0, int_stack+800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[19] = int_stack + 756;
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+774,int_stack+812,int_stack+806, 1.0, int_stack+800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[18] = int_stack + 774;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+792,int_stack+828,int_stack+822,1);
     Libderiv->ABCD[14] = int_stack + 792;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+810,int_stack+844,int_stack+838,1);
     Libderiv->ABCD[13] = int_stack + 810;
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+828,int_stack+860,int_stack+854,1);
     Libderiv->ABCD[12] = int_stack + 828;

}
