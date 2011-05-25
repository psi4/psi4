#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|fd) integrals */

void d12hrr_order_00fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][5][11] = int_stack + 0;
 Libderiv->deriv_classes[0][5][10] = int_stack + 21;
 Libderiv->deriv_classes[0][5][9] = int_stack + 42;
 Libderiv->deriv_classes[0][5][8] = int_stack + 63;
 Libderiv->deriv_classes[0][5][7] = int_stack + 84;
 Libderiv->dvrr_classes[0][4] = int_stack + 105;
 Libderiv->deriv_classes[0][5][6] = int_stack + 120;
 Libderiv->deriv_classes[0][5][2] = int_stack + 141;
 Libderiv->deriv_classes[0][5][1] = int_stack + 162;
 Libderiv->deriv_classes[0][5][0] = int_stack + 183;
 Libderiv->deriv2_classes[0][3][143] = int_stack + 204;
 Libderiv->deriv2_classes[0][4][143] = int_stack + 214;
 Libderiv->deriv2_classes[0][5][143] = int_stack + 229;
 Libderiv->deriv2_classes[0][3][131] = int_stack + 250;
 Libderiv->deriv2_classes[0][4][131] = int_stack + 260;
 Libderiv->deriv2_classes[0][5][131] = int_stack + 275;
 Libderiv->deriv2_classes[0][3][130] = int_stack + 296;
 Libderiv->deriv2_classes[0][4][130] = int_stack + 306;
 Libderiv->deriv2_classes[0][5][130] = int_stack + 321;
 Libderiv->deriv2_classes[0][3][119] = int_stack + 342;
 Libderiv->deriv2_classes[0][4][119] = int_stack + 352;
 Libderiv->deriv2_classes[0][5][119] = int_stack + 367;
 Libderiv->deriv2_classes[0][3][118] = int_stack + 388;
 Libderiv->deriv2_classes[0][4][118] = int_stack + 398;
 Libderiv->deriv2_classes[0][5][118] = int_stack + 413;
 Libderiv->deriv2_classes[0][3][117] = int_stack + 434;
 Libderiv->deriv2_classes[0][4][117] = int_stack + 444;
 Libderiv->deriv2_classes[0][5][117] = int_stack + 459;
 Libderiv->deriv2_classes[0][3][107] = int_stack + 480;
 Libderiv->deriv2_classes[0][4][107] = int_stack + 490;
 Libderiv->deriv2_classes[0][5][107] = int_stack + 505;
 Libderiv->deriv2_classes[0][3][106] = int_stack + 526;
 Libderiv->deriv2_classes[0][4][106] = int_stack + 536;
 Libderiv->deriv2_classes[0][5][106] = int_stack + 551;
 Libderiv->deriv2_classes[0][3][105] = int_stack + 572;
 Libderiv->deriv2_classes[0][4][105] = int_stack + 582;
 Libderiv->deriv2_classes[0][5][105] = int_stack + 597;
 Libderiv->deriv2_classes[0][3][104] = int_stack + 618;
 Libderiv->deriv2_classes[0][4][104] = int_stack + 628;
 Libderiv->deriv2_classes[0][5][104] = int_stack + 643;
 Libderiv->deriv2_classes[0][3][95] = int_stack + 664;
 Libderiv->deriv2_classes[0][4][95] = int_stack + 674;
 Libderiv->deriv2_classes[0][5][95] = int_stack + 689;
 Libderiv->deriv2_classes[0][3][94] = int_stack + 710;
 Libderiv->deriv2_classes[0][4][94] = int_stack + 720;
 Libderiv->deriv2_classes[0][5][94] = int_stack + 735;
 Libderiv->deriv2_classes[0][3][93] = int_stack + 756;
 Libderiv->deriv2_classes[0][4][93] = int_stack + 766;
 Libderiv->deriv2_classes[0][5][93] = int_stack + 781;
 Libderiv->deriv2_classes[0][3][92] = int_stack + 802;
 Libderiv->deriv2_classes[0][4][92] = int_stack + 812;
 Libderiv->deriv2_classes[0][5][92] = int_stack + 827;
 Libderiv->deriv2_classes[0][3][91] = int_stack + 848;
 Libderiv->deriv2_classes[0][4][91] = int_stack + 858;
 Libderiv->deriv2_classes[0][5][91] = int_stack + 873;
 Libderiv->deriv_classes[0][3][11] = int_stack + 894;
 Libderiv->deriv2_classes[0][3][83] = int_stack + 904;
 Libderiv->deriv_classes[0][4][11] = int_stack + 914;
 Libderiv->deriv2_classes[0][4][83] = int_stack + 929;
 Libderiv->deriv2_classes[0][5][83] = int_stack + 944;
 Libderiv->deriv_classes[0][3][10] = int_stack + 965;
 Libderiv->deriv2_classes[0][3][82] = int_stack + 975;
 Libderiv->deriv_classes[0][4][10] = int_stack + 985;
 Libderiv->deriv2_classes[0][4][82] = int_stack + 1000;
 Libderiv->deriv2_classes[0][5][82] = int_stack + 1015;
 Libderiv->deriv_classes[0][3][9] = int_stack + 1036;
 Libderiv->deriv2_classes[0][3][81] = int_stack + 1046;
 Libderiv->deriv_classes[0][4][9] = int_stack + 1056;
 Libderiv->deriv2_classes[0][4][81] = int_stack + 1071;
 Libderiv->deriv2_classes[0][5][81] = int_stack + 1086;
 Libderiv->deriv_classes[0][3][8] = int_stack + 1107;
 Libderiv->deriv2_classes[0][3][80] = int_stack + 1117;
 Libderiv->deriv_classes[0][4][8] = int_stack + 1127;
 Libderiv->deriv2_classes[0][4][80] = int_stack + 1142;
 Libderiv->deriv2_classes[0][5][80] = int_stack + 1157;
 Libderiv->deriv_classes[0][3][7] = int_stack + 1178;
 Libderiv->deriv2_classes[0][3][79] = int_stack + 1188;
 Libderiv->deriv_classes[0][4][7] = int_stack + 1198;
 Libderiv->deriv2_classes[0][4][79] = int_stack + 1213;
 Libderiv->deriv2_classes[0][5][79] = int_stack + 1228;
 Libderiv->dvrr_classes[0][3] = int_stack + 1249;
 Libderiv->deriv_classes[0][3][6] = int_stack + 1259;
 Libderiv->deriv2_classes[0][3][78] = int_stack + 1269;
 Libderiv->deriv_classes[0][4][6] = int_stack + 1279;
 Libderiv->deriv2_classes[0][4][78] = int_stack + 1294;
 Libderiv->deriv2_classes[0][5][78] = int_stack + 1309;
 Libderiv->deriv2_classes[0][3][35] = int_stack + 1330;
 Libderiv->deriv2_classes[0][4][35] = int_stack + 1340;
 Libderiv->deriv2_classes[0][5][35] = int_stack + 1355;
 Libderiv->deriv2_classes[0][3][34] = int_stack + 1376;
 Libderiv->deriv2_classes[0][4][34] = int_stack + 1386;
 Libderiv->deriv2_classes[0][5][34] = int_stack + 1401;
 Libderiv->deriv2_classes[0][3][33] = int_stack + 1422;
 Libderiv->deriv2_classes[0][4][33] = int_stack + 1432;
 Libderiv->deriv2_classes[0][5][33] = int_stack + 1447;
 Libderiv->deriv2_classes[0][3][32] = int_stack + 1468;
 Libderiv->deriv2_classes[0][4][32] = int_stack + 1478;
 Libderiv->deriv2_classes[0][5][32] = int_stack + 1493;
 Libderiv->deriv2_classes[0][3][31] = int_stack + 1514;
 Libderiv->deriv2_classes[0][4][31] = int_stack + 1524;
 Libderiv->deriv2_classes[0][5][31] = int_stack + 1539;
 Libderiv->deriv_classes[0][3][2] = int_stack + 1560;
 Libderiv->deriv2_classes[0][3][30] = int_stack + 1570;
 Libderiv->deriv_classes[0][4][2] = int_stack + 1580;
 Libderiv->deriv2_classes[0][4][30] = int_stack + 1595;
 Libderiv->deriv2_classes[0][5][30] = int_stack + 1610;
 Libderiv->deriv2_classes[0][3][26] = int_stack + 1631;
 Libderiv->deriv2_classes[0][4][26] = int_stack + 1641;
 Libderiv->deriv2_classes[0][5][26] = int_stack + 1656;
 Libderiv->deriv2_classes[0][3][23] = int_stack + 1677;
 Libderiv->deriv2_classes[0][4][23] = int_stack + 1687;
 Libderiv->deriv2_classes[0][5][23] = int_stack + 1702;
 Libderiv->deriv2_classes[0][3][22] = int_stack + 1723;
 Libderiv->deriv2_classes[0][4][22] = int_stack + 1733;
 Libderiv->deriv2_classes[0][5][22] = int_stack + 1748;
 Libderiv->deriv2_classes[0][3][21] = int_stack + 1769;
 Libderiv->deriv2_classes[0][4][21] = int_stack + 1779;
 Libderiv->deriv2_classes[0][5][21] = int_stack + 1794;
 Libderiv->deriv2_classes[0][3][20] = int_stack + 1815;
 Libderiv->deriv2_classes[0][4][20] = int_stack + 1825;
 Libderiv->deriv2_classes[0][5][20] = int_stack + 1840;
 Libderiv->deriv2_classes[0][3][19] = int_stack + 1861;
 Libderiv->deriv2_classes[0][4][19] = int_stack + 1871;
 Libderiv->deriv2_classes[0][5][19] = int_stack + 1886;
 Libderiv->deriv_classes[0][3][1] = int_stack + 1907;
 Libderiv->deriv2_classes[0][3][18] = int_stack + 1917;
 Libderiv->deriv_classes[0][4][1] = int_stack + 1927;
 Libderiv->deriv2_classes[0][4][18] = int_stack + 1942;
 Libderiv->deriv2_classes[0][5][18] = int_stack + 1957;
 Libderiv->deriv2_classes[0][3][14] = int_stack + 1978;
 Libderiv->deriv2_classes[0][4][14] = int_stack + 1988;
 Libderiv->deriv2_classes[0][5][14] = int_stack + 2003;
 Libderiv->deriv2_classes[0][3][13] = int_stack + 2024;
 Libderiv->deriv2_classes[0][4][13] = int_stack + 2034;
 Libderiv->deriv2_classes[0][5][13] = int_stack + 2049;
 Libderiv->deriv2_classes[0][3][11] = int_stack + 2070;
 Libderiv->deriv2_classes[0][4][11] = int_stack + 2080;
 Libderiv->deriv2_classes[0][5][11] = int_stack + 2095;
 Libderiv->deriv2_classes[0][3][10] = int_stack + 2116;
 Libderiv->deriv2_classes[0][4][10] = int_stack + 2126;
 Libderiv->deriv2_classes[0][5][10] = int_stack + 2141;
 Libderiv->deriv2_classes[0][3][9] = int_stack + 2162;
 Libderiv->deriv2_classes[0][4][9] = int_stack + 2172;
 Libderiv->deriv2_classes[0][5][9] = int_stack + 2187;
 Libderiv->deriv2_classes[0][3][8] = int_stack + 2208;
 Libderiv->deriv2_classes[0][4][8] = int_stack + 2218;
 Libderiv->deriv2_classes[0][5][8] = int_stack + 2233;
 Libderiv->deriv2_classes[0][3][7] = int_stack + 2254;
 Libderiv->deriv2_classes[0][4][7] = int_stack + 2264;
 Libderiv->deriv2_classes[0][5][7] = int_stack + 2279;
 Libderiv->deriv_classes[0][3][0] = int_stack + 2300;
 Libderiv->deriv2_classes[0][3][6] = int_stack + 2310;
 Libderiv->deriv_classes[0][4][0] = int_stack + 2320;
 Libderiv->deriv2_classes[0][4][6] = int_stack + 2335;
 Libderiv->deriv2_classes[0][5][6] = int_stack + 2350;
 Libderiv->deriv2_classes[0][3][2] = int_stack + 2371;
 Libderiv->deriv2_classes[0][4][2] = int_stack + 2381;
 Libderiv->deriv2_classes[0][5][2] = int_stack + 2396;
 Libderiv->deriv2_classes[0][3][1] = int_stack + 2417;
 Libderiv->deriv2_classes[0][4][1] = int_stack + 2427;
 Libderiv->deriv2_classes[0][5][1] = int_stack + 2442;
 Libderiv->deriv2_classes[0][3][0] = int_stack + 2463;
 Libderiv->deriv2_classes[0][4][0] = int_stack + 2473;
 Libderiv->deriv2_classes[0][5][0] = int_stack + 2488;
 memset(int_stack,0,20072);

 Libderiv->dvrr_stack = int_stack + 4384;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2509,int_stack+105,int_stack+1249,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2539,int_stack+914,int_stack+894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1249,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2569,int_stack+0,int_stack+914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2614,int_stack+985,int_stack+965, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1249, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2644,int_stack+21,int_stack+985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+1056,int_stack+1036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1249, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2689,int_stack+42,int_stack+1056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30,int_stack+1127,int_stack+1107, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2734,int_stack+63,int_stack+1127, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2779,int_stack+1198,int_stack+1178, 0.0, zero_stack, 1.0, int_stack+1249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2809,int_stack+84,int_stack+1198, 0.0, zero_stack, 1.0, int_stack+105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60,int_stack+1279,int_stack+1259, 1.0, int_stack+1249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2854,int_stack+120,int_stack+1279, 1.0, int_stack+105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+1580,int_stack+1560,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2899,int_stack+141,int_stack+1580,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+120,int_stack+1927,int_stack+1907,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2944,int_stack+162,int_stack+1927,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+150,int_stack+2320,int_stack+2300,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2989,int_stack+183,int_stack+2320,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3034,int_stack+214,int_stack+204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+894,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3064,int_stack+229,int_stack+214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+914,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+260,int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+894, 1.0, int_stack+965,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+210,int_stack+275,int_stack+260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+914, 1.0, int_stack+985,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+255,int_stack+306,int_stack+296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+965, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3109,int_stack+321,int_stack+306, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+985, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+285,int_stack+352,int_stack+342, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+894, 0.0, zero_stack, 1.0, int_stack+1036,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3154,int_stack+367,int_stack+352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+914, 0.0, zero_stack, 1.0, int_stack+1056,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+315,int_stack+398,int_stack+388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+965, 1.0, int_stack+1036, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+345,int_stack+413,int_stack+398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+985, 1.0, int_stack+1056, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+390,int_stack+444,int_stack+434, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1036, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3199,int_stack+459,int_stack+444, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1056, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+420,int_stack+490,int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+894, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1107,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3244,int_stack+505,int_stack+490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+914, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1127,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+536,int_stack+526, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+965, 0.0, zero_stack, 1.0, int_stack+1107, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+480,int_stack+551,int_stack+536, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+985, 0.0, zero_stack, 1.0, int_stack+1127, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+525,int_stack+582,int_stack+572, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1036, 1.0, int_stack+1107, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3289,int_stack+597,int_stack+582, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1056, 1.0, int_stack+1127, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+555,int_stack+628,int_stack+618, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1107, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3334,int_stack+643,int_stack+628, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+585,int_stack+674,int_stack+664, 0.0, zero_stack, 1.0, int_stack+894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1178,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+615,int_stack+689,int_stack+674, 0.0, zero_stack, 1.0, int_stack+914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1198,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+660,int_stack+720,int_stack+710, 0.0, zero_stack, 1.0, int_stack+965, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1178, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3379,int_stack+735,int_stack+720, 0.0, zero_stack, 1.0, int_stack+985, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1198, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+690,int_stack+766,int_stack+756, 0.0, zero_stack, 1.0, int_stack+1036, 0.0, zero_stack, 1.0, int_stack+1178, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+720,int_stack+781,int_stack+766, 0.0, zero_stack, 1.0, int_stack+1056, 0.0, zero_stack, 1.0, int_stack+1198, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+765,int_stack+812,int_stack+802, 0.0, zero_stack, 1.0, int_stack+1107, 1.0, int_stack+1178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3424,int_stack+827,int_stack+812, 0.0, zero_stack, 1.0, int_stack+1127, 1.0, int_stack+1198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+795,int_stack+858,int_stack+848, 0.0, zero_stack, 2.0, int_stack+1178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3469,int_stack+873,int_stack+858, 0.0, zero_stack, 2.0, int_stack+1198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+825,int_stack+929,int_stack+904, 1.0, int_stack+894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1259,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+855,int_stack+944,int_stack+929, 1.0, int_stack+914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1279,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1000,int_stack+975, 1.0, int_stack+965, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1259, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+930,int_stack+1015,int_stack+1000, 1.0, int_stack+985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1279, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+975,int_stack+1071,int_stack+1046, 1.0, int_stack+1036, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1259, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1005,int_stack+1086,int_stack+1071, 1.0, int_stack+1056, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1279, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1050,int_stack+1142,int_stack+1117, 1.0, int_stack+1107, 0.0, zero_stack, 1.0, int_stack+1259, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1080,int_stack+1157,int_stack+1142, 1.0, int_stack+1127, 0.0, zero_stack, 1.0, int_stack+1279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1125,int_stack+1213,int_stack+1188, 1.0, int_stack+1178, 1.0, int_stack+1259, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3514,int_stack+1228,int_stack+1213, 1.0, int_stack+1198, 1.0, int_stack+1279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1155,int_stack+1294,int_stack+1269, 2.0, int_stack+1259, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1185,int_stack+1309,int_stack+1294, 2.0, int_stack+1279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1230,int_stack+1340,int_stack+1330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+1355,int_stack+1340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1580,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1305,int_stack+1386,int_stack+1376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1335,int_stack+1401,int_stack+1386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1580, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1380,int_stack+1432,int_stack+1422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3559,int_stack+1447,int_stack+1432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1580, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1410,int_stack+1478,int_stack+1468, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3604,int_stack+1493,int_stack+1478, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1440,int_stack+1524,int_stack+1514, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1470,int_stack+1539,int_stack+1524, 0.0, zero_stack, 1.0, int_stack+1580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1515,int_stack+1595,int_stack+1570, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3649,int_stack+1610,int_stack+1595, 1.0, int_stack+1580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1545,int_stack+1641,int_stack+1631,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1575,int_stack+1656,int_stack+1641,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1620,int_stack+1687,int_stack+1677, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1907,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3694,int_stack+1702,int_stack+1687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1927,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1650,int_stack+1733,int_stack+1723, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1907, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1680,int_stack+1748,int_stack+1733, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1927, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1725,int_stack+1779,int_stack+1769, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1907, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3739,int_stack+1794,int_stack+1779, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1927, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1755,int_stack+1825,int_stack+1815, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3784,int_stack+1840,int_stack+1825, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1785,int_stack+1871,int_stack+1861, 0.0, zero_stack, 1.0, int_stack+1907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1815,int_stack+1886,int_stack+1871, 0.0, zero_stack, 1.0, int_stack+1927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1860,int_stack+1942,int_stack+1917, 1.0, int_stack+1907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3829,int_stack+1957,int_stack+1942, 1.0, int_stack+1927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1890,int_stack+1988,int_stack+1978,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1920,int_stack+2003,int_stack+1988,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1965,int_stack+2034,int_stack+2024,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3874,int_stack+2049,int_stack+2034,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1995,int_stack+2080,int_stack+2070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2300,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2025,int_stack+2095,int_stack+2080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2320,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2070,int_stack+2126,int_stack+2116, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2300, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3919,int_stack+2141,int_stack+2126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2320, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+2172,int_stack+2162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2300, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3964,int_stack+2187,int_stack+2172, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2320, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2130,int_stack+2218,int_stack+2208, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2160,int_stack+2233,int_stack+2218, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2205,int_stack+2264,int_stack+2254, 0.0, zero_stack, 1.0, int_stack+2300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4009,int_stack+2279,int_stack+2264, 0.0, zero_stack, 1.0, int_stack+2320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2235,int_stack+2335,int_stack+2310, 1.0, int_stack+2300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2265,int_stack+2350,int_stack+2335, 1.0, int_stack+2320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2310,int_stack+2381,int_stack+2371,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4054,int_stack+2396,int_stack+2381,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2340,int_stack+2427,int_stack+2417,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2370,int_stack+2442,int_stack+2427,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2415,int_stack+2473,int_stack+2463,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4099,int_stack+2488,int_stack+2473,1);
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2445,int_stack+2569,int_stack+2539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2509,1);
     Libderiv->ABCD[11] = int_stack + 2445;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4144,int_stack+2644,int_stack+2614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2509, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 4144;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4204,int_stack+2689,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2509, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 4204;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2644,int_stack+2734,int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 2644;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2704,int_stack+2809,int_stack+2779, 0.0, zero_stack, 1.0, int_stack+2509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 2704;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4264,int_stack+2854,int_stack+60, 1.0, int_stack+2509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 4264;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2809,int_stack+2899,int_stack+90,1);
     Libderiv->ABCD[2] = int_stack + 2809;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2869,int_stack+2944,int_stack+120,1);
     Libderiv->ABCD[1] = int_stack + 2869;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2929,int_stack+2989,int_stack+150,1);
     Libderiv->ABCD[0] = int_stack + 2929;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4324,int_stack+3064,int_stack+3034, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2539,1);
     Libderiv->ABCD[155] = int_stack + 4324;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2989,int_stack+210,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2539, 1.0, int_stack+2614,1);
     Libderiv->ABCD[143] = int_stack + 2989;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3049,int_stack+3109,int_stack+255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2614, 0.0, zero_stack,1);
     Libderiv->ABCD[142] = int_stack + 3049;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+180,int_stack+3154,int_stack+285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2539, 0.0, zero_stack, 1.0, int_stack+0,1);
     Libderiv->ABCD[131] = int_stack + 180;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3109,int_stack+345,int_stack+315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2614, 1.0, int_stack+0, 0.0, zero_stack,1);
     Libderiv->ABCD[130] = int_stack + 3109;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+240,int_stack+3199,int_stack+390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[129] = int_stack + 240;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3169,int_stack+3244,int_stack+420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2539, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30,1);
     Libderiv->ABCD[119] = int_stack + 3169;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3229,int_stack+480,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2614, 0.0, zero_stack, 1.0, int_stack+30, 0.0, zero_stack,1);
     Libderiv->ABCD[118] = int_stack + 3229;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+300,int_stack+3289,int_stack+525, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+30, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[117] = int_stack + 300;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+360,int_stack+3334,int_stack+555, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[116] = int_stack + 360;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3289,int_stack+615,int_stack+585, 0.0, zero_stack, 1.0, int_stack+2539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2779,1);
     Libderiv->ABCD[107] = int_stack + 3289;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+420,int_stack+3379,int_stack+660, 0.0, zero_stack, 1.0, int_stack+2614, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2779, 0.0, zero_stack,1);
     Libderiv->ABCD[106] = int_stack + 420;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3349,int_stack+720,int_stack+690, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+2779, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[105] = int_stack + 3349;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+480,int_stack+3424,int_stack+765, 0.0, zero_stack, 1.0, int_stack+30, 1.0, int_stack+2779, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[104] = int_stack + 480;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3409,int_stack+3469,int_stack+795, 0.0, zero_stack, 2.0, int_stack+2779, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[103] = int_stack + 3409;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+540,int_stack+855,int_stack+825, 1.0, int_stack+2539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60,1);
     Libderiv->ABCD[95] = int_stack + 540;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+600,int_stack+930,int_stack+900, 1.0, int_stack+2614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60, 0.0, zero_stack,1);
     Libderiv->ABCD[94] = int_stack + 600;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+660,int_stack+1005,int_stack+975, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[93] = int_stack + 660;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+720,int_stack+1080,int_stack+1050, 1.0, int_stack+30, 0.0, zero_stack, 1.0, int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[92] = int_stack + 720;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+3514,int_stack+1125, 1.0, int_stack+2779, 1.0, int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[91] = int_stack + 0;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3469,int_stack+1185,int_stack+1155, 2.0, int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[90] = int_stack + 3469;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+780,int_stack+1260,int_stack+1230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90,1);
     Libderiv->ABCD[47] = int_stack + 780;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+840,int_stack+1335,int_stack+1305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack,1);
     Libderiv->ABCD[46] = int_stack + 840;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+900,int_stack+3559,int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[45] = int_stack + 900;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3529,int_stack+3604,int_stack+1410, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[44] = int_stack + 3529;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3589,int_stack+1470,int_stack+1440, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[43] = int_stack + 3589;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+960,int_stack+3649,int_stack+1515, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[42] = int_stack + 960;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+60,int_stack+1575,int_stack+1545,1);
     Libderiv->ABCD[38] = int_stack + 60;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1020,int_stack+3694,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120,1);
     Libderiv->ABCD[35] = int_stack + 1020;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3649,int_stack+1680,int_stack+1650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120, 0.0, zero_stack,1);
     Libderiv->ABCD[34] = int_stack + 3649;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1080,int_stack+3739,int_stack+1725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[33] = int_stack + 1080;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3709,int_stack+3784,int_stack+1755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[32] = int_stack + 3709;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3769,int_stack+1815,int_stack+1785, 0.0, zero_stack, 1.0, int_stack+120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[31] = int_stack + 3769;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1140,int_stack+3829,int_stack+1860, 1.0, int_stack+120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[30] = int_stack + 1140;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+1920,int_stack+1890,1);
     Libderiv->ABCD[26] = int_stack + 1200;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1260,int_stack+3874,int_stack+1965,1);
     Libderiv->ABCD[25] = int_stack + 1260;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3829,int_stack+2025,int_stack+1995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150,1);
     Libderiv->ABCD[23] = int_stack + 3829;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1320,int_stack+3919,int_stack+2070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150, 0.0, zero_stack,1);
     Libderiv->ABCD[22] = int_stack + 1320;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3889,int_stack+3964,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[21] = int_stack + 3889;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3949,int_stack+2160,int_stack+2130, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[20] = int_stack + 3949;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1380,int_stack+4009,int_stack+2205, 0.0, zero_stack, 1.0, int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[19] = int_stack + 1380;
 /*--- compute (00|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1440,int_stack+2265,int_stack+2235, 1.0, int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[18] = int_stack + 1440;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+120,int_stack+4054,int_stack+2310,1);
     Libderiv->ABCD[14] = int_stack + 120;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4009,int_stack+2370,int_stack+2340,1);
     Libderiv->ABCD[13] = int_stack + 4009;
 /*--- compute (00|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1500,int_stack+4099,int_stack+2415,1);
     Libderiv->ABCD[12] = int_stack + 1500;

}
