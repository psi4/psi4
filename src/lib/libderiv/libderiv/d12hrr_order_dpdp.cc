#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_dpdp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|dp) integrals */

void d12hrr_order_dpdp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][3][10] = int_stack + 100;
 Libderiv->deriv_classes[3][3][9] = int_stack + 200;
 Libderiv->deriv_classes[3][3][8] = int_stack + 300;
 Libderiv->deriv_classes[3][3][7] = int_stack + 400;
 Libderiv->dvrr_classes[3][2] = int_stack + 500;
 Libderiv->deriv_classes[3][3][6] = int_stack + 560;
 Libderiv->deriv_classes[3][3][2] = int_stack + 660;
 Libderiv->deriv_classes[3][3][1] = int_stack + 760;
 Libderiv->dvrr_classes[2][3] = int_stack + 860;
 Libderiv->deriv_classes[3][3][0] = int_stack + 920;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 1020;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 1056;
 Libderiv->deriv2_classes[3][2][143] = int_stack + 1116;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 1176;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 1276;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1312;
 Libderiv->deriv2_classes[3][2][131] = int_stack + 1372;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 1432;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 1532;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 1568;
 Libderiv->deriv2_classes[3][2][130] = int_stack + 1628;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 1688;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 1788;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 1824;
 Libderiv->deriv2_classes[3][2][119] = int_stack + 1884;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 1944;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 2044;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 2080;
 Libderiv->deriv2_classes[3][2][118] = int_stack + 2140;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 2200;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 2300;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 2336;
 Libderiv->deriv2_classes[3][2][117] = int_stack + 2396;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 2456;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 2556;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 2592;
 Libderiv->deriv2_classes[3][2][107] = int_stack + 2652;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 2712;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 2812;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 2848;
 Libderiv->deriv2_classes[3][2][106] = int_stack + 2908;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 2968;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 3068;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 3104;
 Libderiv->deriv2_classes[3][2][105] = int_stack + 3164;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 3224;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 3324;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 3360;
 Libderiv->deriv2_classes[3][2][104] = int_stack + 3420;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 3480;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 3580;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 3616;
 Libderiv->deriv2_classes[3][2][95] = int_stack + 3676;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 3736;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 3836;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 3872;
 Libderiv->deriv2_classes[3][2][94] = int_stack + 3932;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 3992;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 4092;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 4128;
 Libderiv->deriv2_classes[3][2][93] = int_stack + 4188;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 4248;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 4348;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 4384;
 Libderiv->deriv2_classes[3][2][92] = int_stack + 4444;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 4504;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 4604;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 4640;
 Libderiv->deriv2_classes[3][2][91] = int_stack + 4700;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 4760;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 4860;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 4896;
 Libderiv->deriv_classes[3][2][11] = int_stack + 4956;
 Libderiv->deriv2_classes[3][2][83] = int_stack + 5016;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 5076;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 5176;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 5212;
 Libderiv->deriv_classes[3][2][10] = int_stack + 5272;
 Libderiv->deriv2_classes[3][2][82] = int_stack + 5332;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 5392;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 5492;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 5528;
 Libderiv->deriv_classes[3][2][9] = int_stack + 5588;
 Libderiv->deriv2_classes[3][2][81] = int_stack + 5648;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 5708;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 5808;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 5844;
 Libderiv->deriv_classes[3][2][8] = int_stack + 5904;
 Libderiv->deriv2_classes[3][2][80] = int_stack + 5964;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 6024;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 6124;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 6160;
 Libderiv->deriv_classes[3][2][7] = int_stack + 6220;
 Libderiv->deriv2_classes[3][2][79] = int_stack + 6280;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 6340;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 6440;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 6476;
 Libderiv->deriv_classes[3][2][6] = int_stack + 6536;
 Libderiv->deriv2_classes[3][2][78] = int_stack + 6596;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 6656;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 6756;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 6792;
 Libderiv->deriv2_classes[3][2][35] = int_stack + 6852;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 6912;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 7012;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 7048;
 Libderiv->deriv2_classes[3][2][34] = int_stack + 7108;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 7168;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 7268;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 7304;
 Libderiv->deriv2_classes[3][2][33] = int_stack + 7364;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 7424;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 7524;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 7560;
 Libderiv->deriv2_classes[3][2][32] = int_stack + 7620;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 7680;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 7780;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 7816;
 Libderiv->deriv2_classes[3][2][31] = int_stack + 7876;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 7936;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 8036;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 8072;
 Libderiv->deriv_classes[3][2][2] = int_stack + 8132;
 Libderiv->deriv2_classes[3][2][30] = int_stack + 8192;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 8252;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 8352;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 8388;
 Libderiv->deriv2_classes[3][2][26] = int_stack + 8448;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 8508;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 8608;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 8644;
 Libderiv->deriv2_classes[3][2][23] = int_stack + 8704;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 8764;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 8864;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 8900;
 Libderiv->deriv2_classes[3][2][22] = int_stack + 8960;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 9020;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 9120;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 9156;
 Libderiv->deriv2_classes[3][2][21] = int_stack + 9216;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 9276;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 9376;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 9412;
 Libderiv->deriv2_classes[3][2][20] = int_stack + 9472;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 9532;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 9632;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 9668;
 Libderiv->deriv2_classes[3][2][19] = int_stack + 9728;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 9788;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 9888;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 9924;
 Libderiv->deriv_classes[3][2][1] = int_stack + 9984;
 Libderiv->deriv2_classes[3][2][18] = int_stack + 10044;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 10104;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 10204;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 10240;
 Libderiv->deriv2_classes[3][2][14] = int_stack + 10300;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 10360;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 10460;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 10496;
 Libderiv->deriv2_classes[3][2][13] = int_stack + 10556;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 10616;
 Libderiv->deriv_classes[2][2][11] = int_stack + 10716;
 Libderiv->deriv_classes[2][3][11] = int_stack + 10752;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 10812;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 10848;
 Libderiv->deriv2_classes[3][2][11] = int_stack + 10908;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 10968;
 Libderiv->deriv_classes[2][2][10] = int_stack + 11068;
 Libderiv->deriv_classes[2][3][10] = int_stack + 11104;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 11164;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 11200;
 Libderiv->deriv2_classes[3][2][10] = int_stack + 11260;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 11320;
 Libderiv->deriv_classes[2][2][9] = int_stack + 11420;
 Libderiv->deriv_classes[2][3][9] = int_stack + 11456;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 11516;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 11552;
 Libderiv->deriv2_classes[3][2][9] = int_stack + 11612;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 11672;
 Libderiv->deriv_classes[2][2][8] = int_stack + 11772;
 Libderiv->deriv_classes[2][3][8] = int_stack + 11808;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 11868;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 11904;
 Libderiv->deriv2_classes[3][2][8] = int_stack + 11964;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 12024;
 Libderiv->deriv_classes[2][2][7] = int_stack + 12124;
 Libderiv->deriv_classes[2][3][7] = int_stack + 12160;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 12220;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 12256;
 Libderiv->deriv2_classes[3][2][7] = int_stack + 12316;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 12376;
 Libderiv->dvrr_classes[2][2] = int_stack + 12476;
 Libderiv->deriv_classes[2][2][6] = int_stack + 12512;
 Libderiv->deriv_classes[2][3][6] = int_stack + 12548;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 12608;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 12644;
 Libderiv->deriv_classes[3][2][0] = int_stack + 12704;
 Libderiv->deriv2_classes[3][2][6] = int_stack + 12764;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 12824;
 Libderiv->deriv_classes[2][2][2] = int_stack + 12924;
 Libderiv->deriv_classes[2][3][2] = int_stack + 12960;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 13020;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 13056;
 Libderiv->deriv2_classes[3][2][2] = int_stack + 13116;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 13176;
 Libderiv->deriv_classes[2][2][1] = int_stack + 13276;
 Libderiv->deriv_classes[2][3][1] = int_stack + 13312;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 13372;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 13408;
 Libderiv->deriv2_classes[3][2][1] = int_stack + 13468;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 13528;
 Libderiv->deriv_classes[2][2][0] = int_stack + 13628;
 Libderiv->deriv_classes[2][3][0] = int_stack + 13664;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 13724;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 13760;
 Libderiv->deriv2_classes[3][2][0] = int_stack + 13820;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 13880;
 memset(int_stack,0,111840);

 Libderiv->dvrr_stack = int_stack + 19452;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_dpdp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+13980,int_stack+10752,int_stack+10716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12476,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14088,int_stack+0,int_stack+4956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14268,int_stack+11104,int_stack+11068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12476, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14376,int_stack+100,int_stack+5272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+11456,int_stack+11420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12476, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14556,int_stack+200,int_stack+5588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+11808,int_stack+11772, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14736,int_stack+300,int_stack+5904, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+12160,int_stack+12124, 0.0, zero_stack, 1.0, int_stack+12476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+14916,int_stack+400,int_stack+6220, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+12548,int_stack+12512, 1.0, int_stack+12476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15096,int_stack+560,int_stack+6536, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+860,int_stack+12476,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+12960,int_stack+12924,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15276,int_stack+660,int_stack+8132,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+13312,int_stack+13276,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15456,int_stack+760,int_stack+9984,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+756,int_stack+13664,int_stack+13628,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15636,int_stack+920,int_stack+12704,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+864,int_stack+1056,int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10716,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15816,int_stack+1176,int_stack+1116, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4956,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+972,int_stack+1312,int_stack+1276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10716, 1.0, int_stack+11068,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1432,int_stack+1372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4956, 1.0, int_stack+5272,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1260,int_stack+1568,int_stack+1532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11068, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1368,int_stack+1688,int_stack+1628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5272, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1548,int_stack+1824,int_stack+1788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10716, 0.0, zero_stack, 1.0, int_stack+11420,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1656,int_stack+1944,int_stack+1884, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4956, 0.0, zero_stack, 1.0, int_stack+5588,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1836,int_stack+2080,int_stack+2044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11068, 1.0, int_stack+11420, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1944,int_stack+2200,int_stack+2140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5272, 1.0, int_stack+5588, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2124,int_stack+2336,int_stack+2300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11420, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15996,int_stack+2456,int_stack+2396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5588, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2232,int_stack+2592,int_stack+2556, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10716, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11772,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2340,int_stack+2712,int_stack+2652, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4956, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5904,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2520,int_stack+2848,int_stack+2812, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11068, 0.0, zero_stack, 1.0, int_stack+11772, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2628,int_stack+2968,int_stack+2908, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5272, 0.0, zero_stack, 1.0, int_stack+5904, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2808,int_stack+3104,int_stack+3068, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11420, 1.0, int_stack+11772, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2916,int_stack+3224,int_stack+3164, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5588, 1.0, int_stack+5904, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3096,int_stack+3360,int_stack+3324, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3204,int_stack+3480,int_stack+3420, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3384,int_stack+3616,int_stack+3580, 0.0, zero_stack, 1.0, int_stack+10716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12124,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3492,int_stack+3736,int_stack+3676, 0.0, zero_stack, 1.0, int_stack+4956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6220,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3672,int_stack+3872,int_stack+3836, 0.0, zero_stack, 1.0, int_stack+11068, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12124, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+16176,int_stack+3992,int_stack+3932, 0.0, zero_stack, 1.0, int_stack+5272, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6220, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3780,int_stack+4128,int_stack+4092, 0.0, zero_stack, 1.0, int_stack+11420, 0.0, zero_stack, 1.0, int_stack+12124, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3888,int_stack+4248,int_stack+4188, 0.0, zero_stack, 1.0, int_stack+5588, 0.0, zero_stack, 1.0, int_stack+6220, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4068,int_stack+4384,int_stack+4348, 0.0, zero_stack, 1.0, int_stack+11772, 1.0, int_stack+12124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4176,int_stack+4504,int_stack+4444, 0.0, zero_stack, 1.0, int_stack+5904, 1.0, int_stack+6220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4356,int_stack+4640,int_stack+4604, 0.0, zero_stack, 2.0, int_stack+12124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4464,int_stack+4760,int_stack+4700, 0.0, zero_stack, 2.0, int_stack+6220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4644,int_stack+4896,int_stack+4860, 1.0, int_stack+10716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12512,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4752,int_stack+5076,int_stack+5016, 1.0, int_stack+4956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6536,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4932,int_stack+5212,int_stack+5176, 1.0, int_stack+11068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12512, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5040,int_stack+5392,int_stack+5332, 1.0, int_stack+5272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6536, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5220,int_stack+5528,int_stack+5492, 1.0, int_stack+11420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12512, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5328,int_stack+5708,int_stack+5648, 1.0, int_stack+5588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6536, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5508,int_stack+5844,int_stack+5808, 1.0, int_stack+11772, 0.0, zero_stack, 1.0, int_stack+12512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5616,int_stack+6024,int_stack+5964, 1.0, int_stack+5904, 0.0, zero_stack, 1.0, int_stack+6536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5796,int_stack+6160,int_stack+6124, 1.0, int_stack+12124, 1.0, int_stack+12512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5904,int_stack+6340,int_stack+6280, 1.0, int_stack+6220, 1.0, int_stack+6536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6084,int_stack+6476,int_stack+6440, 2.0, int_stack+12512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6192,int_stack+6656,int_stack+6596, 2.0, int_stack+6536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12476,int_stack+6792,int_stack+6756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12924,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6372,int_stack+6912,int_stack+6852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8132,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6552,int_stack+7048,int_stack+7012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12924, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6660,int_stack+7168,int_stack+7108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8132, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6840,int_stack+7304,int_stack+7268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12924, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6948,int_stack+7424,int_stack+7364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8132, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7128,int_stack+7560,int_stack+7524, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12924, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7236,int_stack+7680,int_stack+7620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7416,int_stack+7816,int_stack+7780, 0.0, zero_stack, 1.0, int_stack+12924, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7524,int_stack+7936,int_stack+7876, 0.0, zero_stack, 1.0, int_stack+8132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7704,int_stack+8072,int_stack+8036, 1.0, int_stack+12924, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7812,int_stack+8252,int_stack+8192, 1.0, int_stack+8132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7992,int_stack+8388,int_stack+8352,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+8100,int_stack+8508,int_stack+8448,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8280,int_stack+8644,int_stack+8608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13276,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8388,int_stack+8764,int_stack+8704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9984,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8568,int_stack+8900,int_stack+8864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13276, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8676,int_stack+9020,int_stack+8960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9984, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8856,int_stack+9156,int_stack+9120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13276, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8964,int_stack+9276,int_stack+9216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9984, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9144,int_stack+9412,int_stack+9376, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9252,int_stack+9532,int_stack+9472, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9432,int_stack+9668,int_stack+9632, 0.0, zero_stack, 1.0, int_stack+13276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9540,int_stack+9788,int_stack+9728, 0.0, zero_stack, 1.0, int_stack+9984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9720,int_stack+9924,int_stack+9888, 1.0, int_stack+13276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+16356,int_stack+10104,int_stack+10044, 1.0, int_stack+9984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+9828,int_stack+10240,int_stack+10204,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+9936,int_stack+10360,int_stack+10300,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+10116,int_stack+10496,int_stack+10460,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+10224,int_stack+10616,int_stack+10556,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10404,int_stack+10848,int_stack+10812, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13628,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10512,int_stack+10968,int_stack+10908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12704,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10692,int_stack+11200,int_stack+11164, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13628, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10800,int_stack+11320,int_stack+11260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12704, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10980,int_stack+11552,int_stack+11516, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13628, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11088,int_stack+11672,int_stack+11612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12704, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11268,int_stack+11904,int_stack+11868, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11376,int_stack+12024,int_stack+11964, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11556,int_stack+12256,int_stack+12220, 0.0, zero_stack, 1.0, int_stack+13628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11664,int_stack+12376,int_stack+12316, 0.0, zero_stack, 1.0, int_stack+12704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11844,int_stack+12644,int_stack+12608, 1.0, int_stack+13628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11952,int_stack+12824,int_stack+12764, 1.0, int_stack+12704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12584,int_stack+13056,int_stack+13020,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12692,int_stack+13176,int_stack+13116,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12872,int_stack+13408,int_stack+13372,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12980,int_stack+13528,int_stack+13468,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+13160,int_stack+13760,int_stack+13724,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+13268,int_stack+13880,int_stack+13820,10);
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+13448,int_stack+14088,int_stack+13980,18);
     Libderiv->ABCD[11] = int_stack + 13448;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+12132,int_stack+14376,int_stack+14268,18);
     Libderiv->ABCD[10] = int_stack + 12132;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+16536,int_stack+14556,int_stack+0,18);
     Libderiv->ABCD[9] = int_stack + 16536;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14376,int_stack+14736,int_stack+108,18);
     Libderiv->ABCD[8] = int_stack + 14376;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+16860,int_stack+14916,int_stack+216,18);
     Libderiv->ABCD[7] = int_stack + 16860;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14700,int_stack+15096,int_stack+324,18);
     Libderiv->ABCD[6] = int_stack + 14700;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+17184,int_stack+15276,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[2] = int_stack + 17184;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+15024,int_stack+15456,int_stack+648, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[1] = int_stack + 15024;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+17508,int_stack+15636,int_stack+756, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[0] = int_stack + 17508;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15348,int_stack+15816,int_stack+864,18);
     Libderiv->ABCD[155] = int_stack + 15348;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15672,int_stack+1080,int_stack+972,18);
     Libderiv->ABCD[143] = int_stack + 15672;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+864,int_stack+1368,int_stack+1260,18);
     Libderiv->ABCD[142] = int_stack + 864;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1188,int_stack+1656,int_stack+1548,18);
     Libderiv->ABCD[131] = int_stack + 1188;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1512,int_stack+1944,int_stack+1836,18);
     Libderiv->ABCD[130] = int_stack + 1512;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+17832,int_stack+15996,int_stack+2124,18);
     Libderiv->ABCD[129] = int_stack + 17832;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1836,int_stack+2340,int_stack+2232,18);
     Libderiv->ABCD[119] = int_stack + 1836;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2160,int_stack+2628,int_stack+2520,18);
     Libderiv->ABCD[118] = int_stack + 2160;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2484,int_stack+2916,int_stack+2808,18);
     Libderiv->ABCD[117] = int_stack + 2484;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+18156,int_stack+3204,int_stack+3096,18);
     Libderiv->ABCD[116] = int_stack + 18156;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2808,int_stack+3492,int_stack+3384,18);
     Libderiv->ABCD[107] = int_stack + 2808;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3132,int_stack+16176,int_stack+3672,18);
     Libderiv->ABCD[106] = int_stack + 3132;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3456,int_stack+3888,int_stack+3780,18);
     Libderiv->ABCD[105] = int_stack + 3456;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15996,int_stack+4176,int_stack+4068,18);
     Libderiv->ABCD[104] = int_stack + 15996;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3780,int_stack+4464,int_stack+4356,18);
     Libderiv->ABCD[103] = int_stack + 3780;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4104,int_stack+4752,int_stack+4644,18);
     Libderiv->ABCD[95] = int_stack + 4104;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4428,int_stack+5040,int_stack+4932,18);
     Libderiv->ABCD[94] = int_stack + 4428;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4752,int_stack+5328,int_stack+5220,18);
     Libderiv->ABCD[93] = int_stack + 4752;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5076,int_stack+5616,int_stack+5508,18);
     Libderiv->ABCD[92] = int_stack + 5076;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5400,int_stack+5904,int_stack+5796,18);
     Libderiv->ABCD[91] = int_stack + 5400;
 /*--- compute (dp|dp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5724,int_stack+6192,int_stack+6084,18);
     Libderiv->ABCD[90] = int_stack + 5724;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6048,int_stack+6372,int_stack+12476, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[47] = int_stack + 6048;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+18480,int_stack+6660,int_stack+6552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[46] = int_stack + 18480;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6372,int_stack+6948,int_stack+6840, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[45] = int_stack + 6372;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6696,int_stack+7236,int_stack+7128, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[44] = int_stack + 6696;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7020,int_stack+7524,int_stack+7416, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[43] = int_stack + 7020;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7344,int_stack+7812,int_stack+7704, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[42] = int_stack + 7344;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7668,int_stack+8100,int_stack+7992, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[38] = int_stack + 7668;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+18804,int_stack+8388,int_stack+8280, 0.0, zero_stack, 1.0, int_stack+13980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[35] = int_stack + 18804;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7992,int_stack+8676,int_stack+8568, 0.0, zero_stack, 1.0, int_stack+14268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[34] = int_stack + 7992;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+8316,int_stack+8964,int_stack+8856, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[33] = int_stack + 8316;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+8640,int_stack+9252,int_stack+9144, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[32] = int_stack + 8640;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+8964,int_stack+9540,int_stack+9432, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[31] = int_stack + 8964;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9288,int_stack+16356,int_stack+9720, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[30] = int_stack + 9288;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+19128,int_stack+9936,int_stack+9828, 0.0, zero_stack, 1.0, int_stack+540, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[26] = int_stack + 19128;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9612,int_stack+10224,int_stack+10116, 0.0, zero_stack, 2.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[25] = int_stack + 9612;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9936,int_stack+10512,int_stack+10404, 1.0, int_stack+13980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[23] = int_stack + 9936;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+10260,int_stack+10800,int_stack+10692, 1.0, int_stack+14268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[22] = int_stack + 10260;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+10584,int_stack+11088,int_stack+10980, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[21] = int_stack + 10584;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+10908,int_stack+11376,int_stack+11268, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[20] = int_stack + 10908;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+11232,int_stack+11664,int_stack+11556, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[19] = int_stack + 11232;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+11952,int_stack+11844, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[18] = int_stack + 0;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+11556,int_stack+12692,int_stack+12584, 1.0, int_stack+540, 0.0, zero_stack, 1.0, int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[14] = int_stack + 11556;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+324,int_stack+12980,int_stack+12872, 1.0, int_stack+648, 1.0, int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[13] = int_stack + 324;
 /*--- compute (dp|dp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+13772,int_stack+13268,int_stack+13160, 2.0, int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[12] = int_stack + 13772;

}
