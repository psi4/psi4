#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ppdp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|dp) integrals */

void d12hrr_order_ppdp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][3][10] = int_stack + 60;
 Libderiv->deriv_classes[2][3][9] = int_stack + 120;
 Libderiv->deriv_classes[2][3][8] = int_stack + 180;
 Libderiv->deriv_classes[2][3][7] = int_stack + 240;
 Libderiv->dvrr_classes[2][2] = int_stack + 300;
 Libderiv->deriv_classes[2][3][6] = int_stack + 336;
 Libderiv->deriv_classes[2][3][2] = int_stack + 396;
 Libderiv->deriv_classes[2][3][1] = int_stack + 456;
 Libderiv->dvrr_classes[1][3] = int_stack + 516;
 Libderiv->deriv_classes[2][3][0] = int_stack + 546;
 Libderiv->deriv2_classes[1][2][143] = int_stack + 606;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 624;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 654;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 690;
 Libderiv->deriv2_classes[1][2][131] = int_stack + 750;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 768;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 798;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 834;
 Libderiv->deriv2_classes[1][2][130] = int_stack + 894;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 912;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 942;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 978;
 Libderiv->deriv2_classes[1][2][119] = int_stack + 1038;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 1056;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 1086;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 1122;
 Libderiv->deriv2_classes[1][2][118] = int_stack + 1182;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 1200;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 1230;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 1266;
 Libderiv->deriv2_classes[1][2][117] = int_stack + 1326;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 1344;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 1374;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 1410;
 Libderiv->deriv2_classes[1][2][107] = int_stack + 1470;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 1488;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 1518;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 1554;
 Libderiv->deriv2_classes[1][2][106] = int_stack + 1614;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 1632;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 1662;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 1698;
 Libderiv->deriv2_classes[1][2][105] = int_stack + 1758;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 1776;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 1806;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 1842;
 Libderiv->deriv2_classes[1][2][104] = int_stack + 1902;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 1920;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 1950;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 1986;
 Libderiv->deriv2_classes[1][2][95] = int_stack + 2046;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 2064;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 2094;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 2130;
 Libderiv->deriv2_classes[1][2][94] = int_stack + 2190;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 2208;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 2238;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 2274;
 Libderiv->deriv2_classes[1][2][93] = int_stack + 2334;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 2352;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 2382;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 2418;
 Libderiv->deriv2_classes[1][2][92] = int_stack + 2478;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 2496;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 2526;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 2562;
 Libderiv->deriv2_classes[1][2][91] = int_stack + 2622;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 2640;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 2670;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 2706;
 Libderiv->deriv2_classes[1][2][83] = int_stack + 2766;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 2784;
 Libderiv->deriv_classes[2][2][11] = int_stack + 2814;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 2850;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 2886;
 Libderiv->deriv2_classes[1][2][82] = int_stack + 2946;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 2964;
 Libderiv->deriv_classes[2][2][10] = int_stack + 2994;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 3030;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 3066;
 Libderiv->deriv2_classes[1][2][81] = int_stack + 3126;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 3144;
 Libderiv->deriv_classes[2][2][9] = int_stack + 3174;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 3210;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 3246;
 Libderiv->deriv2_classes[1][2][80] = int_stack + 3306;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 3324;
 Libderiv->deriv_classes[2][2][8] = int_stack + 3354;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 3390;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 3426;
 Libderiv->deriv2_classes[1][2][79] = int_stack + 3486;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 3504;
 Libderiv->deriv_classes[2][2][7] = int_stack + 3534;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 3570;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 3606;
 Libderiv->deriv2_classes[1][2][78] = int_stack + 3666;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 3684;
 Libderiv->deriv_classes[2][2][6] = int_stack + 3714;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 3750;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 3786;
 Libderiv->deriv2_classes[1][2][35] = int_stack + 3846;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 3864;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 3894;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 3930;
 Libderiv->deriv2_classes[1][2][34] = int_stack + 3990;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 4008;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 4038;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 4074;
 Libderiv->deriv2_classes[1][2][33] = int_stack + 4134;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 4152;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 4182;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 4218;
 Libderiv->deriv2_classes[1][2][32] = int_stack + 4278;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 4296;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 4326;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 4362;
 Libderiv->deriv2_classes[1][2][31] = int_stack + 4422;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 4440;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 4470;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 4506;
 Libderiv->deriv2_classes[1][2][30] = int_stack + 4566;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 4584;
 Libderiv->deriv_classes[2][2][2] = int_stack + 4614;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 4650;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 4686;
 Libderiv->deriv2_classes[1][2][26] = int_stack + 4746;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 4764;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 4794;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 4830;
 Libderiv->deriv2_classes[1][2][23] = int_stack + 4890;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 4908;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 4938;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 4974;
 Libderiv->deriv2_classes[1][2][22] = int_stack + 5034;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 5052;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 5082;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 5118;
 Libderiv->deriv2_classes[1][2][21] = int_stack + 5178;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 5196;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 5226;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 5262;
 Libderiv->deriv2_classes[1][2][20] = int_stack + 5322;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 5340;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 5370;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 5406;
 Libderiv->deriv2_classes[1][2][19] = int_stack + 5466;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 5484;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 5514;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 5550;
 Libderiv->deriv2_classes[1][2][18] = int_stack + 5610;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 5628;
 Libderiv->deriv_classes[2][2][1] = int_stack + 5658;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 5694;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 5730;
 Libderiv->deriv2_classes[1][2][14] = int_stack + 5790;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 5808;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 5838;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 5874;
 Libderiv->deriv2_classes[1][2][13] = int_stack + 5934;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 5952;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 5982;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 6018;
 Libderiv->deriv_classes[1][2][11] = int_stack + 6078;
 Libderiv->deriv_classes[1][3][11] = int_stack + 6096;
 Libderiv->deriv2_classes[1][2][11] = int_stack + 6126;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 6144;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 6174;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 6210;
 Libderiv->deriv_classes[1][2][10] = int_stack + 6270;
 Libderiv->deriv_classes[1][3][10] = int_stack + 6288;
 Libderiv->deriv2_classes[1][2][10] = int_stack + 6318;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 6336;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 6366;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 6402;
 Libderiv->deriv_classes[1][2][9] = int_stack + 6462;
 Libderiv->deriv_classes[1][3][9] = int_stack + 6480;
 Libderiv->deriv2_classes[1][2][9] = int_stack + 6510;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 6528;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 6558;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 6594;
 Libderiv->deriv_classes[1][2][8] = int_stack + 6654;
 Libderiv->deriv_classes[1][3][8] = int_stack + 6672;
 Libderiv->deriv2_classes[1][2][8] = int_stack + 6702;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 6720;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 6750;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 6786;
 Libderiv->deriv_classes[1][2][7] = int_stack + 6846;
 Libderiv->deriv_classes[1][3][7] = int_stack + 6864;
 Libderiv->deriv2_classes[1][2][7] = int_stack + 6894;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 6912;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 6942;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 6978;
 Libderiv->dvrr_classes[1][2] = int_stack + 7038;
 Libderiv->deriv_classes[1][2][6] = int_stack + 7056;
 Libderiv->deriv_classes[1][3][6] = int_stack + 7074;
 Libderiv->deriv2_classes[1][2][6] = int_stack + 7104;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 7122;
 Libderiv->deriv_classes[2][2][0] = int_stack + 7152;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 7188;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 7224;
 Libderiv->deriv_classes[1][2][2] = int_stack + 7284;
 Libderiv->deriv_classes[1][3][2] = int_stack + 7302;
 Libderiv->deriv2_classes[1][2][2] = int_stack + 7332;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 7350;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 7380;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 7416;
 Libderiv->deriv_classes[1][2][1] = int_stack + 7476;
 Libderiv->deriv_classes[1][3][1] = int_stack + 7494;
 Libderiv->deriv2_classes[1][2][1] = int_stack + 7524;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 7542;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 7572;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 7608;
 Libderiv->deriv_classes[1][2][0] = int_stack + 7668;
 Libderiv->deriv_classes[1][3][0] = int_stack + 7686;
 Libderiv->deriv2_classes[1][2][0] = int_stack + 7716;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 7734;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 7764;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 7800;
 memset(int_stack,0,62880);

 Libderiv->dvrr_stack = int_stack + 9858;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ppdp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7860,int_stack+6096,int_stack+6078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7038,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7914,int_stack+0,int_stack+2814, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+6288,int_stack+6270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7038, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8022,int_stack+60,int_stack+2994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+6480,int_stack+6462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7038, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8130,int_stack+120,int_stack+3174, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+6672,int_stack+6654, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8238,int_stack+180,int_stack+3354, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+162,int_stack+6864,int_stack+6846, 0.0, zero_stack, 1.0, int_stack+7038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8346,int_stack+240,int_stack+3534, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+7074,int_stack+7056, 1.0, int_stack+7038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8454,int_stack+336,int_stack+3714, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+270,int_stack+516,int_stack+7038,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+7302,int_stack+7284,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+8562,int_stack+396,int_stack+4614,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+378,int_stack+7494,int_stack+7476,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+8670,int_stack+456,int_stack+5658,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+7686,int_stack+7668,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+8778,int_stack+546,int_stack+7152,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+486,int_stack+624,int_stack+606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6078,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+690,int_stack+654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2814,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+768,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6078, 1.0, int_stack+6270,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8886,int_stack+834,int_stack+798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2814, 1.0, int_stack+2994,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+702,int_stack+912,int_stack+894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6270, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+756,int_stack+978,int_stack+942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2994, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+864,int_stack+1056,int_stack+1038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6078, 0.0, zero_stack, 1.0, int_stack+6462,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+918,int_stack+1122,int_stack+1086, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2814, 0.0, zero_stack, 1.0, int_stack+3174,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1026,int_stack+1200,int_stack+1182, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6270, 1.0, int_stack+6462, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1266,int_stack+1230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2994, 1.0, int_stack+3174, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1188,int_stack+1344,int_stack+1326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6462, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1242,int_stack+1410,int_stack+1374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3174, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1350,int_stack+1488,int_stack+1470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6078, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6654,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1404,int_stack+1554,int_stack+1518, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2814, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3354,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1512,int_stack+1632,int_stack+1614, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6270, 0.0, zero_stack, 1.0, int_stack+6654, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8994,int_stack+1698,int_stack+1662, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2994, 0.0, zero_stack, 1.0, int_stack+3354, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1566,int_stack+1776,int_stack+1758, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6462, 1.0, int_stack+6654, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1620,int_stack+1842,int_stack+1806, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3174, 1.0, int_stack+3354, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1728,int_stack+1920,int_stack+1902, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1782,int_stack+1986,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3354, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1890,int_stack+2064,int_stack+2046, 0.0, zero_stack, 1.0, int_stack+6078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6846,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1944,int_stack+2130,int_stack+2094, 0.0, zero_stack, 1.0, int_stack+2814, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3534,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2052,int_stack+2208,int_stack+2190, 0.0, zero_stack, 1.0, int_stack+6270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6846, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2106,int_stack+2274,int_stack+2238, 0.0, zero_stack, 1.0, int_stack+2994, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3534, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2214,int_stack+2352,int_stack+2334, 0.0, zero_stack, 1.0, int_stack+6462, 0.0, zero_stack, 1.0, int_stack+6846, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2268,int_stack+2418,int_stack+2382, 0.0, zero_stack, 1.0, int_stack+3174, 0.0, zero_stack, 1.0, int_stack+3534, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2376,int_stack+2496,int_stack+2478, 0.0, zero_stack, 1.0, int_stack+6654, 1.0, int_stack+6846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9102,int_stack+2562,int_stack+2526, 0.0, zero_stack, 1.0, int_stack+3354, 1.0, int_stack+3534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2430,int_stack+2640,int_stack+2622, 0.0, zero_stack, 2.0, int_stack+6846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2484,int_stack+2706,int_stack+2670, 0.0, zero_stack, 2.0, int_stack+3534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2592,int_stack+2784,int_stack+2766, 1.0, int_stack+6078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7056,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2646,int_stack+2886,int_stack+2850, 1.0, int_stack+2814, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3714,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2754,int_stack+2964,int_stack+2946, 1.0, int_stack+6270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7056, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2808,int_stack+3066,int_stack+3030, 1.0, int_stack+2994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3714, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2916,int_stack+3144,int_stack+3126, 1.0, int_stack+6462, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7056, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2970,int_stack+3246,int_stack+3210, 1.0, int_stack+3174, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3714, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3078,int_stack+3324,int_stack+3306, 1.0, int_stack+6654, 0.0, zero_stack, 1.0, int_stack+7056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3132,int_stack+3426,int_stack+3390, 1.0, int_stack+3354, 0.0, zero_stack, 1.0, int_stack+3714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3240,int_stack+3504,int_stack+3486, 1.0, int_stack+6846, 1.0, int_stack+7056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3294,int_stack+3606,int_stack+3570, 1.0, int_stack+3534, 1.0, int_stack+3714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3402,int_stack+3684,int_stack+3666, 2.0, int_stack+7056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3456,int_stack+3786,int_stack+3750, 2.0, int_stack+3714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7038,int_stack+3864,int_stack+3846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7284,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3564,int_stack+3930,int_stack+3894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4614,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3672,int_stack+4008,int_stack+3990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7284, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3726,int_stack+4074,int_stack+4038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4614, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3834,int_stack+4152,int_stack+4134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7284, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3888,int_stack+4218,int_stack+4182, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4614, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3996,int_stack+4296,int_stack+4278, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4050,int_stack+4362,int_stack+4326, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4158,int_stack+4440,int_stack+4422, 0.0, zero_stack, 1.0, int_stack+7284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4212,int_stack+4506,int_stack+4470, 0.0, zero_stack, 1.0, int_stack+4614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4320,int_stack+4584,int_stack+4566, 1.0, int_stack+7284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4374,int_stack+4686,int_stack+4650, 1.0, int_stack+4614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4482,int_stack+4764,int_stack+4746,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4536,int_stack+4830,int_stack+4794,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4644,int_stack+4908,int_stack+4890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7476,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4698,int_stack+4974,int_stack+4938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5658,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4806,int_stack+5052,int_stack+5034, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7476, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4860,int_stack+5118,int_stack+5082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5658, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4968,int_stack+5196,int_stack+5178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7476, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5022,int_stack+5262,int_stack+5226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5658, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5130,int_stack+5340,int_stack+5322, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5184,int_stack+5406,int_stack+5370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5658, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5292,int_stack+5484,int_stack+5466, 0.0, zero_stack, 1.0, int_stack+7476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5346,int_stack+5550,int_stack+5514, 0.0, zero_stack, 1.0, int_stack+5658, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5454,int_stack+5628,int_stack+5610, 1.0, int_stack+7476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5508,int_stack+5730,int_stack+5694, 1.0, int_stack+5658, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+5616,int_stack+5808,int_stack+5790,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+5670,int_stack+5874,int_stack+5838,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+5778,int_stack+5952,int_stack+5934,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+5832,int_stack+6018,int_stack+5982,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5940,int_stack+6144,int_stack+6126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7668,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5994,int_stack+6210,int_stack+6174, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7152,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6102,int_stack+6336,int_stack+6318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7668, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6156,int_stack+6402,int_stack+6366, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7152, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6264,int_stack+6528,int_stack+6510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7668, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6318,int_stack+6594,int_stack+6558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7152, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6426,int_stack+6720,int_stack+6702, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6480,int_stack+6786,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7152, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6588,int_stack+6912,int_stack+6894, 0.0, zero_stack, 1.0, int_stack+7668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6642,int_stack+6978,int_stack+6942, 0.0, zero_stack, 1.0, int_stack+7152, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6750,int_stack+7122,int_stack+7104, 1.0, int_stack+7668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6804,int_stack+7224,int_stack+7188, 1.0, int_stack+7152, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7092,int_stack+7350,int_stack+7332,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7146,int_stack+7416,int_stack+7380,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7254,int_stack+7542,int_stack+7524,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7308,int_stack+7608,int_stack+7572,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7416,int_stack+7734,int_stack+7716,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7470,int_stack+7800,int_stack+7764,6);
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+7578,int_stack+7914,int_stack+7860,18);
     Libderiv->ABCD[11] = int_stack + 7578;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+9210,int_stack+8022,int_stack+0,18);
     Libderiv->ABCD[10] = int_stack + 9210;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+7914,int_stack+8130,int_stack+54,18);
     Libderiv->ABCD[9] = int_stack + 7914;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+8076,int_stack+8238,int_stack+108,18);
     Libderiv->ABCD[8] = int_stack + 8076;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+9372,int_stack+8346,int_stack+162,18);
     Libderiv->ABCD[7] = int_stack + 9372;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+8238,int_stack+8454,int_stack+216,18);
     Libderiv->ABCD[6] = int_stack + 8238;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8400,int_stack+8562,int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[2] = int_stack + 8400;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9534,int_stack+8670,int_stack+378, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[1] = int_stack + 9534;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8562,int_stack+8778,int_stack+432, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[0] = int_stack + 8562;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+8724,int_stack+540,int_stack+486,18);
     Libderiv->ABCD[155] = int_stack + 8724;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+486,int_stack+8886,int_stack+648,18);
     Libderiv->ABCD[143] = int_stack + 486;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+9696,int_stack+756,int_stack+702,18);
     Libderiv->ABCD[142] = int_stack + 9696;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+648,int_stack+918,int_stack+864,18);
     Libderiv->ABCD[131] = int_stack + 648;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+810,int_stack+1080,int_stack+1026,18);
     Libderiv->ABCD[130] = int_stack + 810;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+972,int_stack+1242,int_stack+1188,18);
     Libderiv->ABCD[129] = int_stack + 972;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1134,int_stack+1404,int_stack+1350,18);
     Libderiv->ABCD[119] = int_stack + 1134;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1296,int_stack+8994,int_stack+1512,18);
     Libderiv->ABCD[118] = int_stack + 1296;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+8886,int_stack+1620,int_stack+1566,18);
     Libderiv->ABCD[117] = int_stack + 8886;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1458,int_stack+1782,int_stack+1728,18);
     Libderiv->ABCD[116] = int_stack + 1458;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1620,int_stack+1944,int_stack+1890,18);
     Libderiv->ABCD[107] = int_stack + 1620;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1782,int_stack+2106,int_stack+2052,18);
     Libderiv->ABCD[106] = int_stack + 1782;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1944,int_stack+2268,int_stack+2214,18);
     Libderiv->ABCD[105] = int_stack + 1944;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2106,int_stack+9102,int_stack+2376,18);
     Libderiv->ABCD[104] = int_stack + 2106;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2268,int_stack+2484,int_stack+2430,18);
     Libderiv->ABCD[103] = int_stack + 2268;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2430,int_stack+2646,int_stack+2592,18);
     Libderiv->ABCD[95] = int_stack + 2430;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2592,int_stack+2808,int_stack+2754,18);
     Libderiv->ABCD[94] = int_stack + 2592;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2754,int_stack+2970,int_stack+2916,18);
     Libderiv->ABCD[93] = int_stack + 2754;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2916,int_stack+3132,int_stack+3078,18);
     Libderiv->ABCD[92] = int_stack + 2916;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3078,int_stack+3294,int_stack+3240,18);
     Libderiv->ABCD[91] = int_stack + 3078;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3240,int_stack+3456,int_stack+3402,18);
     Libderiv->ABCD[90] = int_stack + 3240;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3402,int_stack+3564,int_stack+7038, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[47] = int_stack + 3402;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9048,int_stack+3726,int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[46] = int_stack + 9048;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3564,int_stack+3888,int_stack+3834, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[45] = int_stack + 3564;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3726,int_stack+4050,int_stack+3996, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[44] = int_stack + 3726;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3888,int_stack+4212,int_stack+4158, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[43] = int_stack + 3888;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4050,int_stack+4374,int_stack+4320, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[42] = int_stack + 4050;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4212,int_stack+4536,int_stack+4482, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[38] = int_stack + 4212;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4374,int_stack+4698,int_stack+4644, 0.0, zero_stack, 1.0, int_stack+7860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[35] = int_stack + 4374;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4536,int_stack+4860,int_stack+4806, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[34] = int_stack + 4536;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4698,int_stack+5022,int_stack+4968, 0.0, zero_stack, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[33] = int_stack + 4698;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4860,int_stack+5184,int_stack+5130, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[32] = int_stack + 4860;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5022,int_stack+5346,int_stack+5292, 0.0, zero_stack, 1.0, int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[31] = int_stack + 5022;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5184,int_stack+5508,int_stack+5454, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[30] = int_stack + 5184;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5346,int_stack+5670,int_stack+5616, 0.0, zero_stack, 1.0, int_stack+324, 1.0, int_stack+378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[26] = int_stack + 5346;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5508,int_stack+5832,int_stack+5778, 0.0, zero_stack, 2.0, int_stack+378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[25] = int_stack + 5508;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5670,int_stack+5994,int_stack+5940, 1.0, int_stack+7860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[23] = int_stack + 5670;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5832,int_stack+6156,int_stack+6102, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[22] = int_stack + 5832;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5994,int_stack+6318,int_stack+6264, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[21] = int_stack + 5994;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6156,int_stack+6480,int_stack+6426, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[20] = int_stack + 6156;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+6642,int_stack+6588, 1.0, int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[19] = int_stack + 0;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6318,int_stack+6804,int_stack+6750, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[18] = int_stack + 6318;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+162,int_stack+7146,int_stack+7092, 1.0, int_stack+324, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[14] = int_stack + 162;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6480,int_stack+7308,int_stack+7254, 1.0, int_stack+378, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[13] = int_stack + 6480;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6642,int_stack+7470,int_stack+7416, 2.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[12] = int_stack + 6642;

}
