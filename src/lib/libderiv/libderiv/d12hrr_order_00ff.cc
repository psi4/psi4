#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|ff) integrals */

void d12hrr_order_00ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][6][11] = int_stack + 0;
 Libderiv->deriv_classes[0][6][10] = int_stack + 28;
 Libderiv->deriv_classes[0][6][9] = int_stack + 56;
 Libderiv->deriv_classes[0][6][8] = int_stack + 84;
 Libderiv->deriv_classes[0][6][7] = int_stack + 112;
 Libderiv->dvrr_classes[0][5] = int_stack + 140;
 Libderiv->deriv_classes[0][6][6] = int_stack + 161;
 Libderiv->deriv_classes[0][6][2] = int_stack + 189;
 Libderiv->deriv_classes[0][6][1] = int_stack + 217;
 Libderiv->deriv_classes[0][6][0] = int_stack + 245;
 Libderiv->deriv2_classes[0][3][143] = int_stack + 273;
 Libderiv->deriv2_classes[0][4][143] = int_stack + 283;
 Libderiv->deriv2_classes[0][5][143] = int_stack + 298;
 Libderiv->deriv2_classes[0][6][143] = int_stack + 319;
 Libderiv->deriv2_classes[0][3][131] = int_stack + 347;
 Libderiv->deriv2_classes[0][4][131] = int_stack + 357;
 Libderiv->deriv2_classes[0][5][131] = int_stack + 372;
 Libderiv->deriv2_classes[0][6][131] = int_stack + 393;
 Libderiv->deriv2_classes[0][3][130] = int_stack + 421;
 Libderiv->deriv2_classes[0][4][130] = int_stack + 431;
 Libderiv->deriv2_classes[0][5][130] = int_stack + 446;
 Libderiv->deriv2_classes[0][6][130] = int_stack + 467;
 Libderiv->deriv2_classes[0][3][119] = int_stack + 495;
 Libderiv->deriv2_classes[0][4][119] = int_stack + 505;
 Libderiv->deriv2_classes[0][5][119] = int_stack + 520;
 Libderiv->deriv2_classes[0][6][119] = int_stack + 541;
 Libderiv->deriv2_classes[0][3][118] = int_stack + 569;
 Libderiv->deriv2_classes[0][4][118] = int_stack + 579;
 Libderiv->deriv2_classes[0][5][118] = int_stack + 594;
 Libderiv->deriv2_classes[0][6][118] = int_stack + 615;
 Libderiv->deriv2_classes[0][3][117] = int_stack + 643;
 Libderiv->deriv2_classes[0][4][117] = int_stack + 653;
 Libderiv->deriv2_classes[0][5][117] = int_stack + 668;
 Libderiv->deriv2_classes[0][6][117] = int_stack + 689;
 Libderiv->deriv2_classes[0][3][107] = int_stack + 717;
 Libderiv->deriv2_classes[0][4][107] = int_stack + 727;
 Libderiv->deriv2_classes[0][5][107] = int_stack + 742;
 Libderiv->deriv2_classes[0][6][107] = int_stack + 763;
 Libderiv->deriv2_classes[0][3][106] = int_stack + 791;
 Libderiv->deriv2_classes[0][4][106] = int_stack + 801;
 Libderiv->deriv2_classes[0][5][106] = int_stack + 816;
 Libderiv->deriv2_classes[0][6][106] = int_stack + 837;
 Libderiv->deriv2_classes[0][3][105] = int_stack + 865;
 Libderiv->deriv2_classes[0][4][105] = int_stack + 875;
 Libderiv->deriv2_classes[0][5][105] = int_stack + 890;
 Libderiv->deriv2_classes[0][6][105] = int_stack + 911;
 Libderiv->deriv2_classes[0][3][104] = int_stack + 939;
 Libderiv->deriv2_classes[0][4][104] = int_stack + 949;
 Libderiv->deriv2_classes[0][5][104] = int_stack + 964;
 Libderiv->deriv2_classes[0][6][104] = int_stack + 985;
 Libderiv->deriv2_classes[0][3][95] = int_stack + 1013;
 Libderiv->deriv2_classes[0][4][95] = int_stack + 1023;
 Libderiv->deriv2_classes[0][5][95] = int_stack + 1038;
 Libderiv->deriv2_classes[0][6][95] = int_stack + 1059;
 Libderiv->deriv2_classes[0][3][94] = int_stack + 1087;
 Libderiv->deriv2_classes[0][4][94] = int_stack + 1097;
 Libderiv->deriv2_classes[0][5][94] = int_stack + 1112;
 Libderiv->deriv2_classes[0][6][94] = int_stack + 1133;
 Libderiv->deriv2_classes[0][3][93] = int_stack + 1161;
 Libderiv->deriv2_classes[0][4][93] = int_stack + 1171;
 Libderiv->deriv2_classes[0][5][93] = int_stack + 1186;
 Libderiv->deriv2_classes[0][6][93] = int_stack + 1207;
 Libderiv->deriv2_classes[0][3][92] = int_stack + 1235;
 Libderiv->deriv2_classes[0][4][92] = int_stack + 1245;
 Libderiv->deriv2_classes[0][5][92] = int_stack + 1260;
 Libderiv->deriv2_classes[0][6][92] = int_stack + 1281;
 Libderiv->deriv2_classes[0][3][91] = int_stack + 1309;
 Libderiv->deriv2_classes[0][4][91] = int_stack + 1319;
 Libderiv->deriv2_classes[0][5][91] = int_stack + 1334;
 Libderiv->deriv2_classes[0][6][91] = int_stack + 1355;
 Libderiv->deriv_classes[0][3][11] = int_stack + 1383;
 Libderiv->deriv2_classes[0][3][83] = int_stack + 1393;
 Libderiv->deriv_classes[0][4][11] = int_stack + 1403;
 Libderiv->deriv2_classes[0][4][83] = int_stack + 1418;
 Libderiv->deriv_classes[0][5][11] = int_stack + 1433;
 Libderiv->deriv2_classes[0][5][83] = int_stack + 1454;
 Libderiv->deriv2_classes[0][6][83] = int_stack + 1475;
 Libderiv->deriv_classes[0][3][10] = int_stack + 1503;
 Libderiv->deriv2_classes[0][3][82] = int_stack + 1513;
 Libderiv->deriv_classes[0][4][10] = int_stack + 1523;
 Libderiv->deriv2_classes[0][4][82] = int_stack + 1538;
 Libderiv->deriv_classes[0][5][10] = int_stack + 1553;
 Libderiv->deriv2_classes[0][5][82] = int_stack + 1574;
 Libderiv->deriv2_classes[0][6][82] = int_stack + 1595;
 Libderiv->deriv_classes[0][3][9] = int_stack + 1623;
 Libderiv->deriv2_classes[0][3][81] = int_stack + 1633;
 Libderiv->deriv_classes[0][4][9] = int_stack + 1643;
 Libderiv->deriv2_classes[0][4][81] = int_stack + 1658;
 Libderiv->deriv_classes[0][5][9] = int_stack + 1673;
 Libderiv->deriv2_classes[0][5][81] = int_stack + 1694;
 Libderiv->deriv2_classes[0][6][81] = int_stack + 1715;
 Libderiv->deriv_classes[0][3][8] = int_stack + 1743;
 Libderiv->deriv2_classes[0][3][80] = int_stack + 1753;
 Libderiv->deriv_classes[0][4][8] = int_stack + 1763;
 Libderiv->deriv2_classes[0][4][80] = int_stack + 1778;
 Libderiv->deriv_classes[0][5][8] = int_stack + 1793;
 Libderiv->deriv2_classes[0][5][80] = int_stack + 1814;
 Libderiv->deriv2_classes[0][6][80] = int_stack + 1835;
 Libderiv->deriv_classes[0][3][7] = int_stack + 1863;
 Libderiv->deriv2_classes[0][3][79] = int_stack + 1873;
 Libderiv->deriv_classes[0][4][7] = int_stack + 1883;
 Libderiv->deriv2_classes[0][4][79] = int_stack + 1898;
 Libderiv->deriv_classes[0][5][7] = int_stack + 1913;
 Libderiv->deriv2_classes[0][5][79] = int_stack + 1934;
 Libderiv->deriv2_classes[0][6][79] = int_stack + 1955;
 Libderiv->dvrr_classes[0][3] = int_stack + 1983;
 Libderiv->deriv_classes[0][3][6] = int_stack + 1993;
 Libderiv->deriv2_classes[0][3][78] = int_stack + 2003;
 Libderiv->dvrr_classes[0][4] = int_stack + 2013;
 Libderiv->deriv_classes[0][4][6] = int_stack + 2028;
 Libderiv->deriv2_classes[0][4][78] = int_stack + 2043;
 Libderiv->deriv_classes[0][5][6] = int_stack + 2058;
 Libderiv->deriv2_classes[0][5][78] = int_stack + 2079;
 Libderiv->deriv2_classes[0][6][78] = int_stack + 2100;
 Libderiv->deriv2_classes[0][3][35] = int_stack + 2128;
 Libderiv->deriv2_classes[0][4][35] = int_stack + 2138;
 Libderiv->deriv2_classes[0][5][35] = int_stack + 2153;
 Libderiv->deriv2_classes[0][6][35] = int_stack + 2174;
 Libderiv->deriv2_classes[0][3][34] = int_stack + 2202;
 Libderiv->deriv2_classes[0][4][34] = int_stack + 2212;
 Libderiv->deriv2_classes[0][5][34] = int_stack + 2227;
 Libderiv->deriv2_classes[0][6][34] = int_stack + 2248;
 Libderiv->deriv2_classes[0][3][33] = int_stack + 2276;
 Libderiv->deriv2_classes[0][4][33] = int_stack + 2286;
 Libderiv->deriv2_classes[0][5][33] = int_stack + 2301;
 Libderiv->deriv2_classes[0][6][33] = int_stack + 2322;
 Libderiv->deriv2_classes[0][3][32] = int_stack + 2350;
 Libderiv->deriv2_classes[0][4][32] = int_stack + 2360;
 Libderiv->deriv2_classes[0][5][32] = int_stack + 2375;
 Libderiv->deriv2_classes[0][6][32] = int_stack + 2396;
 Libderiv->deriv2_classes[0][3][31] = int_stack + 2424;
 Libderiv->deriv2_classes[0][4][31] = int_stack + 2434;
 Libderiv->deriv2_classes[0][5][31] = int_stack + 2449;
 Libderiv->deriv2_classes[0][6][31] = int_stack + 2470;
 Libderiv->deriv_classes[0][3][2] = int_stack + 2498;
 Libderiv->deriv2_classes[0][3][30] = int_stack + 2508;
 Libderiv->deriv_classes[0][4][2] = int_stack + 2518;
 Libderiv->deriv2_classes[0][4][30] = int_stack + 2533;
 Libderiv->deriv_classes[0][5][2] = int_stack + 2548;
 Libderiv->deriv2_classes[0][5][30] = int_stack + 2569;
 Libderiv->deriv2_classes[0][6][30] = int_stack + 2590;
 Libderiv->deriv2_classes[0][3][26] = int_stack + 2618;
 Libderiv->deriv2_classes[0][4][26] = int_stack + 2628;
 Libderiv->deriv2_classes[0][5][26] = int_stack + 2643;
 Libderiv->deriv2_classes[0][6][26] = int_stack + 2664;
 Libderiv->deriv2_classes[0][3][23] = int_stack + 2692;
 Libderiv->deriv2_classes[0][4][23] = int_stack + 2702;
 Libderiv->deriv2_classes[0][5][23] = int_stack + 2717;
 Libderiv->deriv2_classes[0][6][23] = int_stack + 2738;
 Libderiv->deriv2_classes[0][3][22] = int_stack + 2766;
 Libderiv->deriv2_classes[0][4][22] = int_stack + 2776;
 Libderiv->deriv2_classes[0][5][22] = int_stack + 2791;
 Libderiv->deriv2_classes[0][6][22] = int_stack + 2812;
 Libderiv->deriv2_classes[0][3][21] = int_stack + 2840;
 Libderiv->deriv2_classes[0][4][21] = int_stack + 2850;
 Libderiv->deriv2_classes[0][5][21] = int_stack + 2865;
 Libderiv->deriv2_classes[0][6][21] = int_stack + 2886;
 Libderiv->deriv2_classes[0][3][20] = int_stack + 2914;
 Libderiv->deriv2_classes[0][4][20] = int_stack + 2924;
 Libderiv->deriv2_classes[0][5][20] = int_stack + 2939;
 Libderiv->deriv2_classes[0][6][20] = int_stack + 2960;
 Libderiv->deriv2_classes[0][3][19] = int_stack + 2988;
 Libderiv->deriv2_classes[0][4][19] = int_stack + 2998;
 Libderiv->deriv2_classes[0][5][19] = int_stack + 3013;
 Libderiv->deriv2_classes[0][6][19] = int_stack + 3034;
 Libderiv->deriv_classes[0][3][1] = int_stack + 3062;
 Libderiv->deriv2_classes[0][3][18] = int_stack + 3072;
 Libderiv->deriv_classes[0][4][1] = int_stack + 3082;
 Libderiv->deriv2_classes[0][4][18] = int_stack + 3097;
 Libderiv->deriv_classes[0][5][1] = int_stack + 3112;
 Libderiv->deriv2_classes[0][5][18] = int_stack + 3133;
 Libderiv->deriv2_classes[0][6][18] = int_stack + 3154;
 Libderiv->deriv2_classes[0][3][14] = int_stack + 3182;
 Libderiv->deriv2_classes[0][4][14] = int_stack + 3192;
 Libderiv->deriv2_classes[0][5][14] = int_stack + 3207;
 Libderiv->deriv2_classes[0][6][14] = int_stack + 3228;
 Libderiv->deriv2_classes[0][3][13] = int_stack + 3256;
 Libderiv->deriv2_classes[0][4][13] = int_stack + 3266;
 Libderiv->deriv2_classes[0][5][13] = int_stack + 3281;
 Libderiv->deriv2_classes[0][6][13] = int_stack + 3302;
 Libderiv->deriv2_classes[0][3][11] = int_stack + 3330;
 Libderiv->deriv2_classes[0][4][11] = int_stack + 3340;
 Libderiv->deriv2_classes[0][5][11] = int_stack + 3355;
 Libderiv->deriv2_classes[0][6][11] = int_stack + 3376;
 Libderiv->deriv2_classes[0][3][10] = int_stack + 3404;
 Libderiv->deriv2_classes[0][4][10] = int_stack + 3414;
 Libderiv->deriv2_classes[0][5][10] = int_stack + 3429;
 Libderiv->deriv2_classes[0][6][10] = int_stack + 3450;
 Libderiv->deriv2_classes[0][3][9] = int_stack + 3478;
 Libderiv->deriv2_classes[0][4][9] = int_stack + 3488;
 Libderiv->deriv2_classes[0][5][9] = int_stack + 3503;
 Libderiv->deriv2_classes[0][6][9] = int_stack + 3524;
 Libderiv->deriv2_classes[0][3][8] = int_stack + 3552;
 Libderiv->deriv2_classes[0][4][8] = int_stack + 3562;
 Libderiv->deriv2_classes[0][5][8] = int_stack + 3577;
 Libderiv->deriv2_classes[0][6][8] = int_stack + 3598;
 Libderiv->deriv2_classes[0][3][7] = int_stack + 3626;
 Libderiv->deriv2_classes[0][4][7] = int_stack + 3636;
 Libderiv->deriv2_classes[0][5][7] = int_stack + 3651;
 Libderiv->deriv2_classes[0][6][7] = int_stack + 3672;
 Libderiv->deriv_classes[0][3][0] = int_stack + 3700;
 Libderiv->deriv2_classes[0][3][6] = int_stack + 3710;
 Libderiv->deriv_classes[0][4][0] = int_stack + 3720;
 Libderiv->deriv2_classes[0][4][6] = int_stack + 3735;
 Libderiv->deriv_classes[0][5][0] = int_stack + 3750;
 Libderiv->deriv2_classes[0][5][6] = int_stack + 3771;
 Libderiv->deriv2_classes[0][6][6] = int_stack + 3792;
 Libderiv->deriv2_classes[0][3][2] = int_stack + 3820;
 Libderiv->deriv2_classes[0][4][2] = int_stack + 3830;
 Libderiv->deriv2_classes[0][5][2] = int_stack + 3845;
 Libderiv->deriv2_classes[0][6][2] = int_stack + 3866;
 Libderiv->deriv2_classes[0][3][1] = int_stack + 3894;
 Libderiv->deriv2_classes[0][4][1] = int_stack + 3904;
 Libderiv->deriv2_classes[0][5][1] = int_stack + 3919;
 Libderiv->deriv2_classes[0][6][1] = int_stack + 3940;
 Libderiv->deriv2_classes[0][3][0] = int_stack + 3968;
 Libderiv->deriv2_classes[0][4][0] = int_stack + 3978;
 Libderiv->deriv2_classes[0][5][0] = int_stack + 3993;
 Libderiv->deriv2_classes[0][6][0] = int_stack + 4014;
 memset(int_stack,0,32336);

 Libderiv->dvrr_stack = int_stack + 9928;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4042,int_stack+2013,int_stack+1983,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4072,int_stack+140,int_stack+2013,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4117,int_stack+4072,int_stack+4042,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4177,int_stack+1403,int_stack+1383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1983,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4207,int_stack+1433,int_stack+1403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2013,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4252,int_stack+4207,int_stack+4177, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4042,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4312,int_stack+0,int_stack+1433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4375,int_stack+4312,int_stack+4207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4072,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4312,int_stack+1523,int_stack+1503, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1983, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4465,int_stack+1553,int_stack+1523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2013, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4510,int_stack+4465,int_stack+4312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4042, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4570,int_stack+28,int_stack+1553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4633,int_stack+4570,int_stack+4465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4072, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4570,int_stack+1643,int_stack+1623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1983, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1673,int_stack+1643, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2013, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4723,int_stack+0,int_stack+4570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4042, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4783,int_stack+56,int_stack+1673, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4846,int_stack+4783,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4072, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4783,int_stack+1763,int_stack+1743, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1983, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+1793,int_stack+1763, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2013, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4981,int_stack+4936,int_stack+4783, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5041,int_stack+84,int_stack+1793, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5104,int_stack+5041,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5041,int_stack+1883,int_stack+1863, 0.0, zero_stack, 1.0, int_stack+1983, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45,int_stack+1913,int_stack+1883, 0.0, zero_stack, 1.0, int_stack+2013, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5194,int_stack+45,int_stack+5041, 0.0, zero_stack, 1.0, int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5254,int_stack+112,int_stack+1913, 0.0, zero_stack, 1.0, int_stack+140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5317,int_stack+5254,int_stack+45, 0.0, zero_stack, 1.0, int_stack+4072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5254,int_stack+2028,int_stack+1993, 1.0, int_stack+1983, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+2058,int_stack+2028, 1.0, int_stack+2013, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5407,int_stack+90,int_stack+5254, 1.0, int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5467,int_stack+161,int_stack+2058, 1.0, int_stack+140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5530,int_stack+5467,int_stack+90, 1.0, int_stack+4072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5467,int_stack+2518,int_stack+2498,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4042,int_stack+2548,int_stack+2518,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5620,int_stack+4042,int_stack+5467,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5680,int_stack+189,int_stack+2548,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5743,int_stack+5680,int_stack+4042,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4087,int_stack+3082,int_stack+3062,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5680,int_stack+3112,int_stack+3082,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+135,int_stack+5680,int_stack+4087,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5833,int_stack+217,int_stack+3112,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5896,int_stack+5833,int_stack+5680,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5833,int_stack+3720,int_stack+3700,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+195,int_stack+3750,int_stack+3720,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5986,int_stack+195,int_stack+5833,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+245,int_stack+3750,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6109,int_stack+6046,int_stack+195,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+283,int_stack+273, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1383,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+298,int_stack+283, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1403,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6244,int_stack+6199,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4177,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+319,int_stack+298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1433,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+240,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4207,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+357,int_stack+347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1383, 1.0, int_stack+1503,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6046,int_stack+372,int_stack+357, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1403, 1.0, int_stack+1523,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6304,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4177, 1.0, int_stack+4312,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6364,int_stack+393,int_stack+372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1433, 1.0, int_stack+1553,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+330,int_stack+6364,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4207, 1.0, int_stack+4465,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+431,int_stack+421, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1503, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+446,int_stack+431, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1523, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6364,int_stack+6199,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4312, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+467,int_stack+446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1553, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6424,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4465, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+505,int_stack+495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1383, 0.0, zero_stack, 1.0, int_stack+1623,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6046,int_stack+520,int_stack+505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1403, 0.0, zero_stack, 1.0, int_stack+1643,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6514,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4177, 0.0, zero_stack, 1.0, int_stack+4570,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6574,int_stack+541,int_stack+520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1433, 0.0, zero_stack, 1.0, int_stack+1673,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6637,int_stack+6574,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4207, 0.0, zero_stack, 1.0, int_stack+0,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+579,int_stack+569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1503, 1.0, int_stack+1623, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+594,int_stack+579, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1523, 1.0, int_stack+1643, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6574,int_stack+6199,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4312, 1.0, int_stack+4570, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+615,int_stack+594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1553, 1.0, int_stack+1673, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6727,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4465, 1.0, int_stack+0, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+653,int_stack+643, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1623, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6046,int_stack+668,int_stack+653, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1643, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6817,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4570, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6877,int_stack+689,int_stack+668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1673, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+420,int_stack+6877,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+727,int_stack+717, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1383, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1743,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+742,int_stack+727, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1403, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1763,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6877,int_stack+6199,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4177, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4783,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+763,int_stack+742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1433, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1793,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+510,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4207, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4936,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+801,int_stack+791, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1503, 0.0, zero_stack, 1.0, int_stack+1743, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6046,int_stack+816,int_stack+801, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1523, 0.0, zero_stack, 1.0, int_stack+1763, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6937,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4312, 0.0, zero_stack, 1.0, int_stack+4783, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+600,int_stack+837,int_stack+816, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1553, 0.0, zero_stack, 1.0, int_stack+1793, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+663,int_stack+600,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4465, 0.0, zero_stack, 1.0, int_stack+4936, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+875,int_stack+865, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1623, 1.0, int_stack+1743, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+890,int_stack+875, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1643, 1.0, int_stack+1763, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+600,int_stack+6199,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4570, 1.0, int_stack+4783, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+911,int_stack+890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1673, 1.0, int_stack+1793, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+753,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+4936, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+949,int_stack+939, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1743, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6046,int_stack+964,int_stack+949, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1763, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+843,int_stack+6046,int_stack+6199, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4783, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6997,int_stack+985,int_stack+964, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1793, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7060,int_stack+6997,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+1023,int_stack+1013, 0.0, zero_stack, 1.0, int_stack+1383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1863,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+1038,int_stack+1023, 0.0, zero_stack, 1.0, int_stack+1403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1883,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6997,int_stack+6199,int_stack+6046, 0.0, zero_stack, 1.0, int_stack+4177, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5041,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+1059,int_stack+1038, 0.0, zero_stack, 1.0, int_stack+1433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1913,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7150,int_stack+6046,int_stack+6199, 0.0, zero_stack, 1.0, int_stack+4207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+1097,int_stack+1087, 0.0, zero_stack, 1.0, int_stack+1503, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1863, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6046,int_stack+1112,int_stack+1097, 0.0, zero_stack, 1.0, int_stack+1523, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1883, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7240,int_stack+6046,int_stack+6199, 0.0, zero_stack, 1.0, int_stack+4312, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5041, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7300,int_stack+1133,int_stack+1112, 0.0, zero_stack, 1.0, int_stack+1553, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1913, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7363,int_stack+7300,int_stack+6046, 0.0, zero_stack, 1.0, int_stack+4465, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+1171,int_stack+1161, 0.0, zero_stack, 1.0, int_stack+1623, 0.0, zero_stack, 1.0, int_stack+1863, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+1186,int_stack+1171, 0.0, zero_stack, 1.0, int_stack+1643, 0.0, zero_stack, 1.0, int_stack+1883, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7300,int_stack+6199,int_stack+6046, 0.0, zero_stack, 1.0, int_stack+4570, 0.0, zero_stack, 1.0, int_stack+5041, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+1207,int_stack+1186, 0.0, zero_stack, 1.0, int_stack+1673, 0.0, zero_stack, 1.0, int_stack+1913, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7453,int_stack+6046,int_stack+6199, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+45, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+1245,int_stack+1235, 0.0, zero_stack, 1.0, int_stack+1743, 1.0, int_stack+1863, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6046,int_stack+1260,int_stack+1245, 0.0, zero_stack, 1.0, int_stack+1763, 1.0, int_stack+1883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7543,int_stack+6046,int_stack+6199, 0.0, zero_stack, 1.0, int_stack+4783, 1.0, int_stack+5041, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7603,int_stack+1281,int_stack+1260, 0.0, zero_stack, 1.0, int_stack+1793, 1.0, int_stack+1913, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7666,int_stack+7603,int_stack+6046, 0.0, zero_stack, 1.0, int_stack+4936, 1.0, int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+1319,int_stack+1309, 0.0, zero_stack, 2.0, int_stack+1863, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+1334,int_stack+1319, 0.0, zero_stack, 2.0, int_stack+1883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7603,int_stack+6199,int_stack+6046, 0.0, zero_stack, 2.0, int_stack+5041, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+1355,int_stack+1334, 0.0, zero_stack, 2.0, int_stack+1913, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7756,int_stack+6046,int_stack+6199, 0.0, zero_stack, 2.0, int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+1418,int_stack+1393, 1.0, int_stack+1383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1993,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6046,int_stack+1454,int_stack+1418, 1.0, int_stack+1403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2028,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7846,int_stack+6046,int_stack+6199, 1.0, int_stack+4177, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5254,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7906,int_stack+1475,int_stack+1454, 1.0, int_stack+1433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2058,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+903,int_stack+7906,int_stack+6046, 1.0, int_stack+4207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6046,int_stack+1538,int_stack+1513, 1.0, int_stack+1503, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1993, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+1574,int_stack+1538, 1.0, int_stack+1523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2028, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7906,int_stack+6199,int_stack+6046, 1.0, int_stack+4312, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5254, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4312,int_stack+1595,int_stack+1574, 1.0, int_stack+1553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2058, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+993,int_stack+4312,int_stack+6199, 1.0, int_stack+4465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4465,int_stack+1658,int_stack+1633, 1.0, int_stack+1623, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1993, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+1694,int_stack+1658, 1.0, int_stack+1643, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2028, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4312,int_stack+6199,int_stack+4465, 1.0, int_stack+4570, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5254, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4570,int_stack+1715,int_stack+1694, 1.0, int_stack+1673, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2058, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1083,int_stack+4570,int_stack+6199, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+1778,int_stack+1753, 1.0, int_stack+1743, 0.0, zero_stack, 1.0, int_stack+1993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+1814,int_stack+1778, 1.0, int_stack+1763, 0.0, zero_stack, 1.0, int_stack+2028, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4570,int_stack+6199,int_stack+0, 1.0, int_stack+4783, 0.0, zero_stack, 1.0, int_stack+5254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4783,int_stack+1835,int_stack+1814, 1.0, int_stack+1793, 0.0, zero_stack, 1.0, int_stack+2058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1173,int_stack+4783,int_stack+6199, 1.0, int_stack+4936, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4936,int_stack+1898,int_stack+1873, 1.0, int_stack+1863, 1.0, int_stack+1993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6199,int_stack+1934,int_stack+1898, 1.0, int_stack+1883, 1.0, int_stack+2028, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4783,int_stack+6199,int_stack+4936, 1.0, int_stack+5041, 1.0, int_stack+5254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5041,int_stack+1955,int_stack+1934, 1.0, int_stack+1913, 1.0, int_stack+2058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1263,int_stack+5041,int_stack+6199, 1.0, int_stack+45, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6199,int_stack+2043,int_stack+2003, 2.0, int_stack+1993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+2079,int_stack+2043, 2.0, int_stack+2028, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5041,int_stack+4936,int_stack+6199, 2.0, int_stack+5254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5254,int_stack+2100,int_stack+2079, 2.0, int_stack+2058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+5254,int_stack+4936, 2.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+2138,int_stack+2128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2498,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+2153,int_stack+2138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2518,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5254,int_stack+4936,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5467,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6046,int_stack+2174,int_stack+2153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2548,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1353,int_stack+6046,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4042,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4936,int_stack+2212,int_stack+2202, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2498, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+2227,int_stack+2212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2518, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6046,int_stack+90,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5467, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4177,int_stack+2248,int_stack+2227, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2548, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1443,int_stack+4177,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4042, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+2286,int_stack+2276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2498, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+2301,int_stack+2286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2518, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4177,int_stack+4936,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5467, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1533,int_stack+2322,int_stack+2301, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2548, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1596,int_stack+1533,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4042, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4936,int_stack+2360,int_stack+2350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2498, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+2375,int_stack+2360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1533,int_stack+90,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5467, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1686,int_stack+2396,int_stack+2375, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1749,int_stack+1686,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+2434,int_stack+2424, 0.0, zero_stack, 1.0, int_stack+2498, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+2449,int_stack+2434, 0.0, zero_stack, 1.0, int_stack+2518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1686,int_stack+4936,int_stack+90, 0.0, zero_stack, 1.0, int_stack+5467, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1839,int_stack+2470,int_stack+2449, 0.0, zero_stack, 1.0, int_stack+2548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1902,int_stack+1839,int_stack+4936, 0.0, zero_stack, 1.0, int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4936,int_stack+2533,int_stack+2508, 1.0, int_stack+2498, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+2569,int_stack+2533, 1.0, int_stack+2518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1839,int_stack+90,int_stack+4936, 1.0, int_stack+5467, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5467,int_stack+2590,int_stack+2569, 1.0, int_stack+2548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1992,int_stack+5467,int_stack+90, 1.0, int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4042,int_stack+2628,int_stack+2618,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+2643,int_stack+2628,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5467,int_stack+90,int_stack+4042,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2082,int_stack+2664,int_stack+2643,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2145,int_stack+2082,int_stack+90,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+2702,int_stack+2692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3062,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4042,int_stack+2717,int_stack+2702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3082,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2082,int_stack+4042,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4087,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2235,int_stack+2738,int_stack+2717, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3112,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2298,int_stack+2235,int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5680,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4042,int_stack+2776,int_stack+2766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3062, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+2791,int_stack+2776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3082, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2235,int_stack+90,int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4087, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2388,int_stack+2812,int_stack+2791, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3112, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2451,int_stack+2388,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5680, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+2850,int_stack+2840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3062, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4042,int_stack+2865,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3082, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2388,int_stack+4042,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4087, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2541,int_stack+2886,int_stack+2865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3112, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2604,int_stack+2541,int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5680, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4042,int_stack+2924,int_stack+2914, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+2939,int_stack+2924, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2541,int_stack+90,int_stack+4042, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4087, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2694,int_stack+2960,int_stack+2939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2757,int_stack+2694,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+2998,int_stack+2988, 0.0, zero_stack, 1.0, int_stack+3062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4042,int_stack+3013,int_stack+2998, 0.0, zero_stack, 1.0, int_stack+3082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2694,int_stack+4042,int_stack+90, 0.0, zero_stack, 1.0, int_stack+4087, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2847,int_stack+3034,int_stack+3013, 0.0, zero_stack, 1.0, int_stack+3112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2910,int_stack+2847,int_stack+4042, 0.0, zero_stack, 1.0, int_stack+5680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4042,int_stack+3097,int_stack+3072, 1.0, int_stack+3062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+3133,int_stack+3097, 1.0, int_stack+3082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2847,int_stack+90,int_stack+4042, 1.0, int_stack+4087, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4042,int_stack+3154,int_stack+3133, 1.0, int_stack+3112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3000,int_stack+4042,int_stack+90, 1.0, int_stack+5680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5680,int_stack+3192,int_stack+3182,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+3207,int_stack+3192,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4042,int_stack+90,int_stack+5680,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5680,int_stack+3228,int_stack+3207,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3090,int_stack+5680,int_stack+90,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+3266,int_stack+3256,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+3281,int_stack+3266,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5680,int_stack+4936,int_stack+90,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3180,int_stack+3302,int_stack+3281,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+7966,int_stack+3180,int_stack+4936,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4936,int_stack+3340,int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3700,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+3355,int_stack+3340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3720,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3180,int_stack+90,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5833,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3240,int_stack+3376,int_stack+3355, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3303,int_stack+3240,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+3414,int_stack+3404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3700, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+3429,int_stack+3414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3720, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3240,int_stack+4936,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5833, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8056,int_stack+3450,int_stack+3429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8119,int_stack+8056,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4936,int_stack+3488,int_stack+3478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3700, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+3503,int_stack+3488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3720, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8056,int_stack+90,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5833, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3393,int_stack+3524,int_stack+3503, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3456,int_stack+3393,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+3562,int_stack+3552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+3577,int_stack+3562, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3393,int_stack+4936,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8209,int_stack+3598,int_stack+3577, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8272,int_stack+8209,int_stack+4936, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4936,int_stack+3636,int_stack+3626, 0.0, zero_stack, 1.0, int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+90,int_stack+3651,int_stack+3636, 0.0, zero_stack, 1.0, int_stack+3720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8209,int_stack+90,int_stack+4936, 0.0, zero_stack, 1.0, int_stack+5833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3546,int_stack+3672,int_stack+3651, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3609,int_stack+3546,int_stack+90, 0.0, zero_stack, 1.0, int_stack+195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+3735,int_stack+3710, 1.0, int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+3771,int_stack+3735, 1.0, int_stack+3720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3546,int_stack+4936,int_stack+90, 1.0, int_stack+5833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5833,int_stack+3792,int_stack+3771, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3699,int_stack+5833,int_stack+4936, 1.0, int_stack+195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+195,int_stack+3830,int_stack+3820,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+3845,int_stack+3830,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+5833,int_stack+4936,int_stack+195,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8362,int_stack+3866,int_stack+3845,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+3789,int_stack+8362,int_stack+4936,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4936,int_stack+3904,int_stack+3894,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+195,int_stack+3919,int_stack+3904,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8362,int_stack+195,int_stack+4936,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8422,int_stack+3940,int_stack+3919,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8485,int_stack+8422,int_stack+195,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+195,int_stack+3978,int_stack+3968,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4936,int_stack+3993,int_stack+3978,1);
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8422,int_stack+4936,int_stack+195,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8575,int_stack+4014,int_stack+3993,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8638,int_stack+8575,int_stack+4936,1);
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8728,int_stack+4375,int_stack+4252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4117,1);
     Libderiv->ABCD[11] = int_stack + 8728;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8828,int_stack+4633,int_stack+4510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4117, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 8828;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3879,int_stack+4846,int_stack+4723, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4117, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 3879;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4843,int_stack+5104,int_stack+4981, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4117, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 4843;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4372,int_stack+5317,int_stack+5194, 0.0, zero_stack, 1.0, int_stack+4117, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 4372;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8928,int_stack+5530,int_stack+5407, 1.0, int_stack+4117, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 8928;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+9028,int_stack+5743,int_stack+5620,1);
     Libderiv->ABCD[2] = int_stack + 9028;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+9128,int_stack+5896,int_stack+135,1);
     Libderiv->ABCD[1] = int_stack + 9128;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+9228,int_stack+6109,int_stack+5986,1);
     Libderiv->ABCD[0] = int_stack + 9228;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9328,int_stack+240,int_stack+6244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4252,1);
     Libderiv->ABCD[155] = int_stack + 9328;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+195,int_stack+330,int_stack+6304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4252, 1.0, int_stack+4510,1);
     Libderiv->ABCD[143] = int_stack + 195;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+295,int_stack+6424,int_stack+6364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4510, 0.0, zero_stack,1);
     Libderiv->ABCD[142] = int_stack + 295;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9428,int_stack+6637,int_stack+6514, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4252, 0.0, zero_stack, 1.0, int_stack+4723,1);
     Libderiv->ABCD[131] = int_stack + 9428;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9528,int_stack+6727,int_stack+6574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4510, 1.0, int_stack+4723, 0.0, zero_stack,1);
     Libderiv->ABCD[130] = int_stack + 9528;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9628,int_stack+420,int_stack+6817, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4723, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[129] = int_stack + 9628;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+395,int_stack+510,int_stack+6877, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4252, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4981,1);
     Libderiv->ABCD[119] = int_stack + 395;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+495,int_stack+663,int_stack+6937, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4510, 0.0, zero_stack, 1.0, int_stack+4981, 0.0, zero_stack,1);
     Libderiv->ABCD[118] = int_stack + 495;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9728,int_stack+753,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4723, 1.0, int_stack+4981, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[117] = int_stack + 9728;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+595,int_stack+7060,int_stack+843, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4981, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[116] = int_stack + 595;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+695,int_stack+7150,int_stack+6997, 0.0, zero_stack, 1.0, int_stack+4252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5194,1);
     Libderiv->ABCD[107] = int_stack + 695;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+795,int_stack+7363,int_stack+7240, 0.0, zero_stack, 1.0, int_stack+4510, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5194, 0.0, zero_stack,1);
     Libderiv->ABCD[106] = int_stack + 795;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9828,int_stack+7453,int_stack+7300, 0.0, zero_stack, 1.0, int_stack+4723, 0.0, zero_stack, 1.0, int_stack+5194, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[105] = int_stack + 9828;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6106,int_stack+7666,int_stack+7543, 0.0, zero_stack, 1.0, int_stack+4981, 1.0, int_stack+5194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[104] = int_stack + 6106;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6206,int_stack+7756,int_stack+7603, 0.0, zero_stack, 2.0, int_stack+5194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[103] = int_stack + 6206;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6306,int_stack+903,int_stack+7846, 1.0, int_stack+4252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5407,1);
     Libderiv->ABCD[95] = int_stack + 6306;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6406,int_stack+993,int_stack+7906, 1.0, int_stack+4510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5407, 0.0, zero_stack,1);
     Libderiv->ABCD[94] = int_stack + 6406;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+895,int_stack+1083,int_stack+4312, 1.0, int_stack+4723, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5407, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[93] = int_stack + 895;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+995,int_stack+1173,int_stack+4570, 1.0, int_stack+4981, 0.0, zero_stack, 1.0, int_stack+5407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[92] = int_stack + 995;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1095,int_stack+1263,int_stack+4783, 1.0, int_stack+5194, 1.0, int_stack+5407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[91] = int_stack + 1095;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1195,int_stack+0,int_stack+5041, 2.0, int_stack+5407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[90] = int_stack + 1195;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+1353,int_stack+5254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5620,1);
     Libderiv->ABCD[47] = int_stack + 0;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1295,int_stack+1443,int_stack+6046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5620, 0.0, zero_stack,1);
     Libderiv->ABCD[46] = int_stack + 1295;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1395,int_stack+1596,int_stack+4177, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5620, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[45] = int_stack + 1395;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4102,int_stack+1749,int_stack+1533, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[44] = int_stack + 4102;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1495,int_stack+1902,int_stack+1686, 0.0, zero_stack, 1.0, int_stack+5620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[43] = int_stack + 1495;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1595,int_stack+1992,int_stack+1839, 1.0, int_stack+5620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[42] = int_stack + 1595;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1695,int_stack+2145,int_stack+5467,1);
     Libderiv->ABCD[38] = int_stack + 1695;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1795,int_stack+2298,int_stack+2082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135,1);
     Libderiv->ABCD[35] = int_stack + 1795;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1895,int_stack+2451,int_stack+2235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack,1);
     Libderiv->ABCD[34] = int_stack + 1895;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1995,int_stack+2604,int_stack+2388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[33] = int_stack + 1995;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2095,int_stack+2757,int_stack+2541, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[32] = int_stack + 2095;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2195,int_stack+2910,int_stack+2694, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[31] = int_stack + 2195;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2295,int_stack+3000,int_stack+2847, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[30] = int_stack + 2295;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+2395,int_stack+3090,int_stack+4042,1);
     Libderiv->ABCD[26] = int_stack + 2395;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+2495,int_stack+7966,int_stack+5680,1);
     Libderiv->ABCD[25] = int_stack + 2495;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2595,int_stack+3303,int_stack+3180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5986,1);
     Libderiv->ABCD[23] = int_stack + 2595;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2695,int_stack+8119,int_stack+3240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5986, 0.0, zero_stack,1);
     Libderiv->ABCD[22] = int_stack + 2695;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2795,int_stack+3456,int_stack+8056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5986, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[21] = int_stack + 2795;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2895,int_stack+8272,int_stack+3393, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[20] = int_stack + 2895;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2995,int_stack+3609,int_stack+8209, 0.0, zero_stack, 1.0, int_stack+5986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[19] = int_stack + 2995;
 /*--- compute (00|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3095,int_stack+3699,int_stack+3546, 1.0, int_stack+5986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[18] = int_stack + 3095;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+3195,int_stack+3789,int_stack+5833,1);
     Libderiv->ABCD[14] = int_stack + 3195;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+3295,int_stack+8485,int_stack+8362,1);
     Libderiv->ABCD[13] = int_stack + 3295;
 /*--- compute (00|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+3395,int_stack+8638,int_stack+8422,1);
     Libderiv->ABCD[12] = int_stack + 3395;

}
