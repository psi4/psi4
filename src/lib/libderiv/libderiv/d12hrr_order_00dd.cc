#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|dd) integrals */

void d12hrr_order_00dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][4][11] = int_stack + 0;
 Libderiv->deriv_classes[0][4][10] = int_stack + 15;
 Libderiv->deriv_classes[0][4][9] = int_stack + 30;
 Libderiv->deriv_classes[0][4][8] = int_stack + 45;
 Libderiv->deriv_classes[0][4][7] = int_stack + 60;
 Libderiv->dvrr_classes[0][3] = int_stack + 75;
 Libderiv->deriv_classes[0][4][6] = int_stack + 85;
 Libderiv->deriv_classes[0][4][2] = int_stack + 100;
 Libderiv->deriv_classes[0][4][1] = int_stack + 115;
 Libderiv->deriv_classes[0][4][0] = int_stack + 130;
 Libderiv->deriv2_classes[0][2][143] = int_stack + 145;
 Libderiv->deriv2_classes[0][3][143] = int_stack + 151;
 Libderiv->deriv2_classes[0][4][143] = int_stack + 161;
 Libderiv->deriv2_classes[0][2][131] = int_stack + 176;
 Libderiv->deriv2_classes[0][3][131] = int_stack + 182;
 Libderiv->deriv2_classes[0][4][131] = int_stack + 192;
 Libderiv->deriv2_classes[0][2][130] = int_stack + 207;
 Libderiv->deriv2_classes[0][3][130] = int_stack + 213;
 Libderiv->deriv2_classes[0][4][130] = int_stack + 223;
 Libderiv->deriv2_classes[0][2][119] = int_stack + 238;
 Libderiv->deriv2_classes[0][3][119] = int_stack + 244;
 Libderiv->deriv2_classes[0][4][119] = int_stack + 254;
 Libderiv->deriv2_classes[0][2][118] = int_stack + 269;
 Libderiv->deriv2_classes[0][3][118] = int_stack + 275;
 Libderiv->deriv2_classes[0][4][118] = int_stack + 285;
 Libderiv->deriv2_classes[0][2][117] = int_stack + 300;
 Libderiv->deriv2_classes[0][3][117] = int_stack + 306;
 Libderiv->deriv2_classes[0][4][117] = int_stack + 316;
 Libderiv->deriv2_classes[0][2][107] = int_stack + 331;
 Libderiv->deriv2_classes[0][3][107] = int_stack + 337;
 Libderiv->deriv2_classes[0][4][107] = int_stack + 347;
 Libderiv->deriv2_classes[0][2][106] = int_stack + 362;
 Libderiv->deriv2_classes[0][3][106] = int_stack + 368;
 Libderiv->deriv2_classes[0][4][106] = int_stack + 378;
 Libderiv->deriv2_classes[0][2][105] = int_stack + 393;
 Libderiv->deriv2_classes[0][3][105] = int_stack + 399;
 Libderiv->deriv2_classes[0][4][105] = int_stack + 409;
 Libderiv->deriv2_classes[0][2][104] = int_stack + 424;
 Libderiv->deriv2_classes[0][3][104] = int_stack + 430;
 Libderiv->deriv2_classes[0][4][104] = int_stack + 440;
 Libderiv->deriv2_classes[0][2][95] = int_stack + 455;
 Libderiv->deriv2_classes[0][3][95] = int_stack + 461;
 Libderiv->deriv2_classes[0][4][95] = int_stack + 471;
 Libderiv->deriv2_classes[0][2][94] = int_stack + 486;
 Libderiv->deriv2_classes[0][3][94] = int_stack + 492;
 Libderiv->deriv2_classes[0][4][94] = int_stack + 502;
 Libderiv->deriv2_classes[0][2][93] = int_stack + 517;
 Libderiv->deriv2_classes[0][3][93] = int_stack + 523;
 Libderiv->deriv2_classes[0][4][93] = int_stack + 533;
 Libderiv->deriv2_classes[0][2][92] = int_stack + 548;
 Libderiv->deriv2_classes[0][3][92] = int_stack + 554;
 Libderiv->deriv2_classes[0][4][92] = int_stack + 564;
 Libderiv->deriv2_classes[0][2][91] = int_stack + 579;
 Libderiv->deriv2_classes[0][3][91] = int_stack + 585;
 Libderiv->deriv2_classes[0][4][91] = int_stack + 595;
 Libderiv->deriv_classes[0][2][11] = int_stack + 610;
 Libderiv->deriv2_classes[0][2][83] = int_stack + 616;
 Libderiv->deriv_classes[0][3][11] = int_stack + 622;
 Libderiv->deriv2_classes[0][3][83] = int_stack + 632;
 Libderiv->deriv2_classes[0][4][83] = int_stack + 642;
 Libderiv->deriv_classes[0][2][10] = int_stack + 657;
 Libderiv->deriv2_classes[0][2][82] = int_stack + 663;
 Libderiv->deriv_classes[0][3][10] = int_stack + 669;
 Libderiv->deriv2_classes[0][3][82] = int_stack + 679;
 Libderiv->deriv2_classes[0][4][82] = int_stack + 689;
 Libderiv->deriv_classes[0][2][9] = int_stack + 704;
 Libderiv->deriv2_classes[0][2][81] = int_stack + 710;
 Libderiv->deriv_classes[0][3][9] = int_stack + 716;
 Libderiv->deriv2_classes[0][3][81] = int_stack + 726;
 Libderiv->deriv2_classes[0][4][81] = int_stack + 736;
 Libderiv->deriv_classes[0][2][8] = int_stack + 751;
 Libderiv->deriv2_classes[0][2][80] = int_stack + 757;
 Libderiv->deriv_classes[0][3][8] = int_stack + 763;
 Libderiv->deriv2_classes[0][3][80] = int_stack + 773;
 Libderiv->deriv2_classes[0][4][80] = int_stack + 783;
 Libderiv->deriv_classes[0][2][7] = int_stack + 798;
 Libderiv->deriv2_classes[0][2][79] = int_stack + 804;
 Libderiv->deriv_classes[0][3][7] = int_stack + 810;
 Libderiv->deriv2_classes[0][3][79] = int_stack + 820;
 Libderiv->deriv2_classes[0][4][79] = int_stack + 830;
 Libderiv->dvrr_classes[0][2] = int_stack + 845;
 Libderiv->deriv_classes[0][2][6] = int_stack + 851;
 Libderiv->deriv2_classes[0][2][78] = int_stack + 857;
 Libderiv->deriv_classes[0][3][6] = int_stack + 863;
 Libderiv->deriv2_classes[0][3][78] = int_stack + 873;
 Libderiv->deriv2_classes[0][4][78] = int_stack + 883;
 Libderiv->deriv2_classes[0][2][35] = int_stack + 898;
 Libderiv->deriv2_classes[0][3][35] = int_stack + 904;
 Libderiv->deriv2_classes[0][4][35] = int_stack + 914;
 Libderiv->deriv2_classes[0][2][34] = int_stack + 929;
 Libderiv->deriv2_classes[0][3][34] = int_stack + 935;
 Libderiv->deriv2_classes[0][4][34] = int_stack + 945;
 Libderiv->deriv2_classes[0][2][33] = int_stack + 960;
 Libderiv->deriv2_classes[0][3][33] = int_stack + 966;
 Libderiv->deriv2_classes[0][4][33] = int_stack + 976;
 Libderiv->deriv2_classes[0][2][32] = int_stack + 991;
 Libderiv->deriv2_classes[0][3][32] = int_stack + 997;
 Libderiv->deriv2_classes[0][4][32] = int_stack + 1007;
 Libderiv->deriv2_classes[0][2][31] = int_stack + 1022;
 Libderiv->deriv2_classes[0][3][31] = int_stack + 1028;
 Libderiv->deriv2_classes[0][4][31] = int_stack + 1038;
 Libderiv->deriv_classes[0][2][2] = int_stack + 1053;
 Libderiv->deriv2_classes[0][2][30] = int_stack + 1059;
 Libderiv->deriv_classes[0][3][2] = int_stack + 1065;
 Libderiv->deriv2_classes[0][3][30] = int_stack + 1075;
 Libderiv->deriv2_classes[0][4][30] = int_stack + 1085;
 Libderiv->deriv2_classes[0][2][26] = int_stack + 1100;
 Libderiv->deriv2_classes[0][3][26] = int_stack + 1106;
 Libderiv->deriv2_classes[0][4][26] = int_stack + 1116;
 Libderiv->deriv2_classes[0][2][23] = int_stack + 1131;
 Libderiv->deriv2_classes[0][3][23] = int_stack + 1137;
 Libderiv->deriv2_classes[0][4][23] = int_stack + 1147;
 Libderiv->deriv2_classes[0][2][22] = int_stack + 1162;
 Libderiv->deriv2_classes[0][3][22] = int_stack + 1168;
 Libderiv->deriv2_classes[0][4][22] = int_stack + 1178;
 Libderiv->deriv2_classes[0][2][21] = int_stack + 1193;
 Libderiv->deriv2_classes[0][3][21] = int_stack + 1199;
 Libderiv->deriv2_classes[0][4][21] = int_stack + 1209;
 Libderiv->deriv2_classes[0][2][20] = int_stack + 1224;
 Libderiv->deriv2_classes[0][3][20] = int_stack + 1230;
 Libderiv->deriv2_classes[0][4][20] = int_stack + 1240;
 Libderiv->deriv2_classes[0][2][19] = int_stack + 1255;
 Libderiv->deriv2_classes[0][3][19] = int_stack + 1261;
 Libderiv->deriv2_classes[0][4][19] = int_stack + 1271;
 Libderiv->deriv_classes[0][2][1] = int_stack + 1286;
 Libderiv->deriv2_classes[0][2][18] = int_stack + 1292;
 Libderiv->deriv_classes[0][3][1] = int_stack + 1298;
 Libderiv->deriv2_classes[0][3][18] = int_stack + 1308;
 Libderiv->deriv2_classes[0][4][18] = int_stack + 1318;
 Libderiv->deriv2_classes[0][2][14] = int_stack + 1333;
 Libderiv->deriv2_classes[0][3][14] = int_stack + 1339;
 Libderiv->deriv2_classes[0][4][14] = int_stack + 1349;
 Libderiv->deriv2_classes[0][2][13] = int_stack + 1364;
 Libderiv->deriv2_classes[0][3][13] = int_stack + 1370;
 Libderiv->deriv2_classes[0][4][13] = int_stack + 1380;
 Libderiv->deriv2_classes[0][2][11] = int_stack + 1395;
 Libderiv->deriv2_classes[0][3][11] = int_stack + 1401;
 Libderiv->deriv2_classes[0][4][11] = int_stack + 1411;
 Libderiv->deriv2_classes[0][2][10] = int_stack + 1426;
 Libderiv->deriv2_classes[0][3][10] = int_stack + 1432;
 Libderiv->deriv2_classes[0][4][10] = int_stack + 1442;
 Libderiv->deriv2_classes[0][2][9] = int_stack + 1457;
 Libderiv->deriv2_classes[0][3][9] = int_stack + 1463;
 Libderiv->deriv2_classes[0][4][9] = int_stack + 1473;
 Libderiv->deriv2_classes[0][2][8] = int_stack + 1488;
 Libderiv->deriv2_classes[0][3][8] = int_stack + 1494;
 Libderiv->deriv2_classes[0][4][8] = int_stack + 1504;
 Libderiv->deriv2_classes[0][2][7] = int_stack + 1519;
 Libderiv->deriv2_classes[0][3][7] = int_stack + 1525;
 Libderiv->deriv2_classes[0][4][7] = int_stack + 1535;
 Libderiv->deriv_classes[0][2][0] = int_stack + 1550;
 Libderiv->deriv2_classes[0][2][6] = int_stack + 1556;
 Libderiv->deriv_classes[0][3][0] = int_stack + 1562;
 Libderiv->deriv2_classes[0][3][6] = int_stack + 1572;
 Libderiv->deriv2_classes[0][4][6] = int_stack + 1582;
 Libderiv->deriv2_classes[0][2][2] = int_stack + 1597;
 Libderiv->deriv2_classes[0][3][2] = int_stack + 1603;
 Libderiv->deriv2_classes[0][4][2] = int_stack + 1613;
 Libderiv->deriv2_classes[0][2][1] = int_stack + 1628;
 Libderiv->deriv2_classes[0][3][1] = int_stack + 1634;
 Libderiv->deriv2_classes[0][4][1] = int_stack + 1644;
 Libderiv->deriv2_classes[0][2][0] = int_stack + 1659;
 Libderiv->deriv2_classes[0][3][0] = int_stack + 1665;
 Libderiv->deriv2_classes[0][4][0] = int_stack + 1675;
 memset(int_stack,0,13520);

 Libderiv->dvrr_stack = int_stack + 2752;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1690,int_stack+75,int_stack+845,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1708,int_stack+622,int_stack+610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+845,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1726,int_stack+0,int_stack+622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1756,int_stack+669,int_stack+657, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+845, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1774,int_stack+15,int_stack+669, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+716,int_stack+704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+845, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1804,int_stack+30,int_stack+716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18,int_stack+763,int_stack+751, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1834,int_stack+45,int_stack+763, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+36,int_stack+810,int_stack+798, 0.0, zero_stack, 1.0, int_stack+845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1864,int_stack+60,int_stack+810, 0.0, zero_stack, 1.0, int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+863,int_stack+851, 1.0, int_stack+845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1894,int_stack+85,int_stack+863, 1.0, int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+72,int_stack+1065,int_stack+1053,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1924,int_stack+100,int_stack+1065,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+90,int_stack+1298,int_stack+1286,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1954,int_stack+115,int_stack+1298,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+1562,int_stack+1550,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1984,int_stack+130,int_stack+1562,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+126,int_stack+151,int_stack+145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+610,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2014,int_stack+161,int_stack+151, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+622,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+144,int_stack+182,int_stack+176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+610, 1.0, int_stack+657,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2044,int_stack+192,int_stack+182, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+622, 1.0, int_stack+669,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+162,int_stack+213,int_stack+207, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+657, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+223,int_stack+213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+669, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+210,int_stack+244,int_stack+238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+610, 0.0, zero_stack, 1.0, int_stack+704,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2074,int_stack+254,int_stack+244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+622, 0.0, zero_stack, 1.0, int_stack+716,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+228,int_stack+275,int_stack+269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+657, 1.0, int_stack+704, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2104,int_stack+285,int_stack+275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+669, 1.0, int_stack+716, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+246,int_stack+306,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+704, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+264,int_stack+316,int_stack+306, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+716, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+294,int_stack+337,int_stack+331, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+610, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+751,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2134,int_stack+347,int_stack+337, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+622, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+763,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+312,int_stack+368,int_stack+362, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+657, 0.0, zero_stack, 1.0, int_stack+751, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+330,int_stack+378,int_stack+368, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+669, 0.0, zero_stack, 1.0, int_stack+763, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+360,int_stack+399,int_stack+393, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+704, 1.0, int_stack+751, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2164,int_stack+409,int_stack+399, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+716, 1.0, int_stack+763, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+378,int_stack+430,int_stack+424, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+751, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+396,int_stack+440,int_stack+430, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+763, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+426,int_stack+461,int_stack+455, 0.0, zero_stack, 1.0, int_stack+610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+798,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2194,int_stack+471,int_stack+461, 0.0, zero_stack, 1.0, int_stack+622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+810,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+444,int_stack+492,int_stack+486, 0.0, zero_stack, 1.0, int_stack+657, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+798, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+462,int_stack+502,int_stack+492, 0.0, zero_stack, 1.0, int_stack+669, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+810, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+492,int_stack+523,int_stack+517, 0.0, zero_stack, 1.0, int_stack+704, 0.0, zero_stack, 1.0, int_stack+798, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2224,int_stack+533,int_stack+523, 0.0, zero_stack, 1.0, int_stack+716, 0.0, zero_stack, 1.0, int_stack+810, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+510,int_stack+554,int_stack+548, 0.0, zero_stack, 1.0, int_stack+751, 1.0, int_stack+798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2254,int_stack+564,int_stack+554, 0.0, zero_stack, 1.0, int_stack+763, 1.0, int_stack+810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+528,int_stack+585,int_stack+579, 0.0, zero_stack, 2.0, int_stack+798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+546,int_stack+595,int_stack+585, 0.0, zero_stack, 2.0, int_stack+810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+576,int_stack+632,int_stack+616, 1.0, int_stack+610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+851,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2284,int_stack+642,int_stack+632, 1.0, int_stack+622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+863,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+594,int_stack+679,int_stack+663, 1.0, int_stack+657, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+851, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+612,int_stack+689,int_stack+679, 1.0, int_stack+669, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+863, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+642,int_stack+726,int_stack+710, 1.0, int_stack+704, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+851, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+660,int_stack+736,int_stack+726, 1.0, int_stack+716, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+863, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+690,int_stack+773,int_stack+757, 1.0, int_stack+751, 0.0, zero_stack, 1.0, int_stack+851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+708,int_stack+783,int_stack+773, 1.0, int_stack+763, 0.0, zero_stack, 1.0, int_stack+863, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+738,int_stack+820,int_stack+804, 1.0, int_stack+798, 1.0, int_stack+851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+756,int_stack+830,int_stack+820, 1.0, int_stack+810, 1.0, int_stack+863, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+786,int_stack+873,int_stack+857, 2.0, int_stack+851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+804,int_stack+883,int_stack+873, 2.0, int_stack+863, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+834,int_stack+904,int_stack+898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1053,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+852,int_stack+914,int_stack+904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1065,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+882,int_stack+935,int_stack+929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1053, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+945,int_stack+935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1065, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+930,int_stack+966,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1053, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2314,int_stack+976,int_stack+966, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1065, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+948,int_stack+997,int_stack+991, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+966,int_stack+1007,int_stack+997, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1065, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+996,int_stack+1028,int_stack+1022, 0.0, zero_stack, 1.0, int_stack+1053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2344,int_stack+1038,int_stack+1028, 0.0, zero_stack, 1.0, int_stack+1065, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1014,int_stack+1075,int_stack+1059, 1.0, int_stack+1053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1032,int_stack+1085,int_stack+1075, 1.0, int_stack+1065, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1062,int_stack+1106,int_stack+1100,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2374,int_stack+1116,int_stack+1106,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1137,int_stack+1131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1286,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1098,int_stack+1147,int_stack+1137, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1298,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1128,int_stack+1168,int_stack+1162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1286, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2404,int_stack+1178,int_stack+1168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1298, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1146,int_stack+1199,int_stack+1193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1286, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1164,int_stack+1209,int_stack+1199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1298, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1194,int_stack+1230,int_stack+1224, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2434,int_stack+1240,int_stack+1230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1212,int_stack+1261,int_stack+1255, 0.0, zero_stack, 1.0, int_stack+1286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1230,int_stack+1271,int_stack+1261, 0.0, zero_stack, 1.0, int_stack+1298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1260,int_stack+1308,int_stack+1292, 1.0, int_stack+1286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2464,int_stack+1318,int_stack+1308, 1.0, int_stack+1298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1278,int_stack+1339,int_stack+1333,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1296,int_stack+1349,int_stack+1339,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1326,int_stack+1370,int_stack+1364,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2494,int_stack+1380,int_stack+1370,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1344,int_stack+1401,int_stack+1395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1550,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1362,int_stack+1411,int_stack+1401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1562,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1392,int_stack+1432,int_stack+1426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1550, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2524,int_stack+1442,int_stack+1432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1562, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1410,int_stack+1463,int_stack+1457, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1550, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1428,int_stack+1473,int_stack+1463, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1562, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1458,int_stack+1494,int_stack+1488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2554,int_stack+1504,int_stack+1494, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1476,int_stack+1525,int_stack+1519, 0.0, zero_stack, 1.0, int_stack+1550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1494,int_stack+1535,int_stack+1525, 0.0, zero_stack, 1.0, int_stack+1562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1524,int_stack+1572,int_stack+1556, 1.0, int_stack+1550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2584,int_stack+1582,int_stack+1572, 1.0, int_stack+1562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1542,int_stack+1603,int_stack+1597,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1560,int_stack+1613,int_stack+1603,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1590,int_stack+1634,int_stack+1628,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2614,int_stack+1644,int_stack+1634,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1608,int_stack+1665,int_stack+1659,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1626,int_stack+1675,int_stack+1665,1);
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2644,int_stack+1726,int_stack+1708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1690,1);
     Libderiv->ABCD[11] = int_stack + 2644;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2680,int_stack+1774,int_stack+1756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1690, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 2680;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2716,int_stack+1804,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1690, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 2716;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1774,int_stack+1834,int_stack+18, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 1774;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1810,int_stack+1864,int_stack+36, 0.0, zero_stack, 1.0, int_stack+1690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 1810;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1846,int_stack+1894,int_stack+54, 1.0, int_stack+1690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 1846;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1882,int_stack+1924,int_stack+72,1);
     Libderiv->ABCD[2] = int_stack + 1882;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1918,int_stack+1954,int_stack+90,1);
     Libderiv->ABCD[1] = int_stack + 1918;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1656,int_stack+1984,int_stack+108,1);
     Libderiv->ABCD[0] = int_stack + 1656;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1954,int_stack+2014,int_stack+126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1708,1);
     Libderiv->ABCD[155] = int_stack + 1954;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1990,int_stack+2044,int_stack+144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1708, 1.0, int_stack+1756,1);
     Libderiv->ABCD[143] = int_stack + 1990;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+126,int_stack+180,int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1756, 0.0, zero_stack,1);
     Libderiv->ABCD[142] = int_stack + 126;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+162,int_stack+2074,int_stack+210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1708, 0.0, zero_stack, 1.0, int_stack+0,1);
     Libderiv->ABCD[131] = int_stack + 162;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2026,int_stack+2104,int_stack+228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1756, 1.0, int_stack+0, 0.0, zero_stack,1);
     Libderiv->ABCD[130] = int_stack + 2026;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+198,int_stack+264,int_stack+246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[129] = int_stack + 198;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+234,int_stack+2134,int_stack+294, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1708, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18,1);
     Libderiv->ABCD[119] = int_stack + 234;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+270,int_stack+330,int_stack+312, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1756, 0.0, zero_stack, 1.0, int_stack+18, 0.0, zero_stack,1);
     Libderiv->ABCD[118] = int_stack + 270;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+306,int_stack+2164,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+18, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[117] = int_stack + 306;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+342,int_stack+396,int_stack+378, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+18, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[116] = int_stack + 342;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+378,int_stack+2194,int_stack+426, 0.0, zero_stack, 1.0, int_stack+1708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36,1);
     Libderiv->ABCD[107] = int_stack + 378;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2062,int_stack+462,int_stack+444, 0.0, zero_stack, 1.0, int_stack+1756, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36, 0.0, zero_stack,1);
     Libderiv->ABCD[106] = int_stack + 2062;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+414,int_stack+2224,int_stack+492, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+36, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[105] = int_stack + 414;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+450,int_stack+2254,int_stack+510, 0.0, zero_stack, 1.0, int_stack+18, 1.0, int_stack+36, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[104] = int_stack + 450;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+486,int_stack+546,int_stack+528, 0.0, zero_stack, 2.0, int_stack+36, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[103] = int_stack + 486;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+522,int_stack+2284,int_stack+576, 1.0, int_stack+1708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54,1);
     Libderiv->ABCD[95] = int_stack + 522;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+558,int_stack+612,int_stack+594, 1.0, int_stack+1756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54, 0.0, zero_stack,1);
     Libderiv->ABCD[94] = int_stack + 558;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+594,int_stack+660,int_stack+642, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[93] = int_stack + 594;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+630,int_stack+708,int_stack+690, 1.0, int_stack+18, 0.0, zero_stack, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[92] = int_stack + 630;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+756,int_stack+738, 1.0, int_stack+36, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[91] = int_stack + 0;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+666,int_stack+804,int_stack+786, 2.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[90] = int_stack + 666;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+36,int_stack+852,int_stack+834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72,1);
     Libderiv->ABCD[47] = int_stack + 36;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+702,int_stack+900,int_stack+882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72, 0.0, zero_stack,1);
     Libderiv->ABCD[46] = int_stack + 702;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+738,int_stack+2314,int_stack+930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[45] = int_stack + 738;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+774,int_stack+966,int_stack+948, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[44] = int_stack + 774;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+810,int_stack+2344,int_stack+996, 0.0, zero_stack, 1.0, int_stack+72, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[43] = int_stack + 810;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+846,int_stack+1032,int_stack+1014, 1.0, int_stack+72, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[42] = int_stack + 846;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+882,int_stack+2374,int_stack+1062,1);
     Libderiv->ABCD[38] = int_stack + 882;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+918,int_stack+1098,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90,1);
     Libderiv->ABCD[35] = int_stack + 918;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+954,int_stack+2404,int_stack+1128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack,1);
     Libderiv->ABCD[34] = int_stack + 954;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+990,int_stack+1164,int_stack+1146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[33] = int_stack + 990;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1026,int_stack+2434,int_stack+1194, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[32] = int_stack + 1026;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1062,int_stack+1230,int_stack+1212, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[31] = int_stack + 1062;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1098,int_stack+2464,int_stack+1260, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[30] = int_stack + 1098;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+72,int_stack+1296,int_stack+1278,1);
     Libderiv->ABCD[26] = int_stack + 72;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1134,int_stack+2494,int_stack+1326,1);
     Libderiv->ABCD[25] = int_stack + 1134;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1170,int_stack+1362,int_stack+1344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108,1);
     Libderiv->ABCD[23] = int_stack + 1170;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1206,int_stack+2524,int_stack+1392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack,1);
     Libderiv->ABCD[22] = int_stack + 1206;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1242,int_stack+1428,int_stack+1410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[21] = int_stack + 1242;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1278,int_stack+2554,int_stack+1458, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[20] = int_stack + 1278;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1314,int_stack+1494,int_stack+1476, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[19] = int_stack + 1314;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1350,int_stack+2584,int_stack+1524, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[18] = int_stack + 1350;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1386,int_stack+1560,int_stack+1542,1);
     Libderiv->ABCD[14] = int_stack + 1386;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1422,int_stack+2614,int_stack+1590,1);
     Libderiv->ABCD[13] = int_stack + 1422;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1458,int_stack+1626,int_stack+1608,1);
     Libderiv->ABCD[12] = int_stack + 1458;

}
