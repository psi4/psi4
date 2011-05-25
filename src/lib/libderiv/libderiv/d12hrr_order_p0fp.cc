#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0fp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|fp) integrals */

void d12hrr_order_p0fp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[1][4][10] = int_stack + 45;
 Libderiv->deriv_classes[1][4][9] = int_stack + 90;
 Libderiv->deriv_classes[1][4][8] = int_stack + 135;
 Libderiv->deriv_classes[1][4][7] = int_stack + 180;
 Libderiv->dvrr_classes[1][3] = int_stack + 225;
 Libderiv->deriv_classes[1][4][6] = int_stack + 255;
 Libderiv->deriv_classes[1][4][2] = int_stack + 300;
 Libderiv->deriv_classes[1][4][1] = int_stack + 345;
 Libderiv->deriv_classes[1][4][0] = int_stack + 390;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 435;
 Libderiv->deriv2_classes[1][4][143] = int_stack + 465;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 510;
 Libderiv->deriv2_classes[1][4][131] = int_stack + 540;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 585;
 Libderiv->deriv2_classes[1][4][130] = int_stack + 615;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 660;
 Libderiv->deriv2_classes[1][4][119] = int_stack + 690;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 735;
 Libderiv->deriv2_classes[1][4][118] = int_stack + 765;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 810;
 Libderiv->deriv2_classes[1][4][117] = int_stack + 840;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 885;
 Libderiv->deriv2_classes[1][4][107] = int_stack + 915;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 960;
 Libderiv->deriv2_classes[1][4][106] = int_stack + 990;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 1035;
 Libderiv->deriv2_classes[1][4][105] = int_stack + 1065;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 1110;
 Libderiv->deriv2_classes[1][4][104] = int_stack + 1140;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 1185;
 Libderiv->deriv2_classes[1][4][95] = int_stack + 1215;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 1260;
 Libderiv->deriv2_classes[1][4][94] = int_stack + 1290;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 1335;
 Libderiv->deriv2_classes[1][4][93] = int_stack + 1365;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 1410;
 Libderiv->deriv2_classes[1][4][92] = int_stack + 1440;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 1485;
 Libderiv->deriv2_classes[1][4][91] = int_stack + 1515;
 Libderiv->deriv_classes[1][3][11] = int_stack + 1560;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 1590;
 Libderiv->deriv2_classes[1][4][83] = int_stack + 1620;
 Libderiv->deriv_classes[1][3][10] = int_stack + 1665;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 1695;
 Libderiv->deriv2_classes[1][4][82] = int_stack + 1725;
 Libderiv->deriv_classes[1][3][9] = int_stack + 1770;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 1800;
 Libderiv->deriv2_classes[1][4][81] = int_stack + 1830;
 Libderiv->deriv_classes[1][3][8] = int_stack + 1875;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 1905;
 Libderiv->deriv2_classes[1][4][80] = int_stack + 1935;
 Libderiv->deriv_classes[1][3][7] = int_stack + 1980;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 2010;
 Libderiv->deriv2_classes[1][4][79] = int_stack + 2040;
 Libderiv->deriv_classes[1][3][6] = int_stack + 2085;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 2115;
 Libderiv->deriv2_classes[1][4][78] = int_stack + 2145;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 2190;
 Libderiv->deriv2_classes[1][4][35] = int_stack + 2220;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 2265;
 Libderiv->deriv2_classes[1][4][34] = int_stack + 2295;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 2340;
 Libderiv->deriv2_classes[1][4][33] = int_stack + 2370;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 2415;
 Libderiv->deriv2_classes[1][4][32] = int_stack + 2445;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 2490;
 Libderiv->deriv2_classes[1][4][31] = int_stack + 2520;
 Libderiv->deriv_classes[1][3][2] = int_stack + 2565;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 2595;
 Libderiv->deriv2_classes[1][4][30] = int_stack + 2625;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 2670;
 Libderiv->deriv2_classes[1][4][26] = int_stack + 2700;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 2745;
 Libderiv->deriv2_classes[1][4][23] = int_stack + 2775;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 2820;
 Libderiv->deriv2_classes[1][4][22] = int_stack + 2850;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 2895;
 Libderiv->deriv2_classes[1][4][21] = int_stack + 2925;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 2970;
 Libderiv->deriv2_classes[1][4][20] = int_stack + 3000;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 3045;
 Libderiv->deriv2_classes[1][4][19] = int_stack + 3075;
 Libderiv->deriv_classes[1][3][1] = int_stack + 3120;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 3150;
 Libderiv->deriv2_classes[1][4][18] = int_stack + 3180;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 3225;
 Libderiv->deriv2_classes[1][4][14] = int_stack + 3255;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 3300;
 Libderiv->deriv2_classes[1][4][13] = int_stack + 3330;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 3375;
 Libderiv->deriv2_classes[1][4][11] = int_stack + 3405;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 3450;
 Libderiv->deriv2_classes[1][4][10] = int_stack + 3480;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 3525;
 Libderiv->deriv2_classes[1][4][9] = int_stack + 3555;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 3600;
 Libderiv->deriv2_classes[1][4][8] = int_stack + 3630;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 3675;
 Libderiv->deriv2_classes[1][4][7] = int_stack + 3705;
 Libderiv->deriv_classes[1][3][0] = int_stack + 3750;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 3780;
 Libderiv->deriv2_classes[1][4][6] = int_stack + 3810;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 3855;
 Libderiv->deriv2_classes[1][4][2] = int_stack + 3885;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 3930;
 Libderiv->deriv2_classes[1][4][1] = int_stack + 3960;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 4005;
 Libderiv->deriv2_classes[1][4][0] = int_stack + 4035;
 memset(int_stack,0,32640);

 Libderiv->dvrr_stack = int_stack + 4980;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0fp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4080,int_stack+0,int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225,3);
     Libderiv->ABCD[11] = int_stack + 4080;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4170,int_stack+45,int_stack+1665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 4170;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+90,int_stack+1770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4260,int_stack+135,int_stack+1875, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 4260;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+180,int_stack+1980, 0.0, zero_stack, 1.0, int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 90;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4350,int_stack+255,int_stack+2085, 1.0, int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 4350;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+300,int_stack+2565,3);
     Libderiv->ABCD[2] = int_stack + 180;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4440,int_stack+345,int_stack+3120,3);
     Libderiv->ABCD[1] = int_stack + 4440;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+270,int_stack+390,int_stack+3750,3);
     Libderiv->ABCD[0] = int_stack + 270;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4530,int_stack+465,int_stack+435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1560,3);
     Libderiv->ABCD[155] = int_stack + 4530;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+540,int_stack+510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 1.0, int_stack+1665,3);
     Libderiv->ABCD[143] = int_stack + 360;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+615,int_stack+585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1665, 0.0, zero_stack,3);
     Libderiv->ABCD[142] = int_stack + 450;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+690,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 1.0, int_stack+1770,3);
     Libderiv->ABCD[131] = int_stack + 540;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+765,int_stack+735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1665, 1.0, int_stack+1770, 0.0, zero_stack,3);
     Libderiv->ABCD[130] = int_stack + 630;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+840,int_stack+810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1770, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[129] = int_stack + 720;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4620,int_stack+915,int_stack+885, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1875,3);
     Libderiv->ABCD[119] = int_stack + 4620;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+810,int_stack+990,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1665, 0.0, zero_stack, 1.0, int_stack+1875, 0.0, zero_stack,3);
     Libderiv->ABCD[118] = int_stack + 810;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1065,int_stack+1035, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1770, 1.0, int_stack+1875, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[117] = int_stack + 900;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+990,int_stack+1140,int_stack+1110, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[116] = int_stack + 990;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+1215,int_stack+1185, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1980,3);
     Libderiv->ABCD[107] = int_stack + 1080;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1170,int_stack+1290,int_stack+1260, 0.0, zero_stack, 1.0, int_stack+1665, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1980, 0.0, zero_stack,3);
     Libderiv->ABCD[106] = int_stack + 1170;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4710,int_stack+1365,int_stack+1335, 0.0, zero_stack, 1.0, int_stack+1770, 0.0, zero_stack, 1.0, int_stack+1980, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[105] = int_stack + 4710;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1260,int_stack+1440,int_stack+1410, 0.0, zero_stack, 1.0, int_stack+1875, 1.0, int_stack+1980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[104] = int_stack + 1260;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1350,int_stack+1515,int_stack+1485, 0.0, zero_stack, 2.0, int_stack+1980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[103] = int_stack + 1350;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1440,int_stack+1620,int_stack+1590, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2085,3);
     Libderiv->ABCD[95] = int_stack + 1440;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1530,int_stack+1725,int_stack+1695, 1.0, int_stack+1665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2085, 0.0, zero_stack,3);
     Libderiv->ABCD[94] = int_stack + 1530;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1620,int_stack+1830,int_stack+1800, 1.0, int_stack+1770, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2085, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[93] = int_stack + 1620;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1710,int_stack+1935,int_stack+1905, 1.0, int_stack+1875, 0.0, zero_stack, 1.0, int_stack+2085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[92] = int_stack + 1710;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2040,int_stack+2010, 1.0, int_stack+1980, 1.0, int_stack+2085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[91] = int_stack + 1800;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1890,int_stack+2145,int_stack+2115, 2.0, int_stack+2085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[90] = int_stack + 1890;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1980,int_stack+2220,int_stack+2190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2565,3);
     Libderiv->ABCD[47] = int_stack + 1980;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2070,int_stack+2295,int_stack+2265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2565, 0.0, zero_stack,3);
     Libderiv->ABCD[46] = int_stack + 2070;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2160,int_stack+2370,int_stack+2340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2565, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[45] = int_stack + 2160;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2250,int_stack+2445,int_stack+2415, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[44] = int_stack + 2250;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2340,int_stack+2520,int_stack+2490, 0.0, zero_stack, 1.0, int_stack+2565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[43] = int_stack + 2340;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2430,int_stack+2625,int_stack+2595, 1.0, int_stack+2565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[42] = int_stack + 2430;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2520,int_stack+2700,int_stack+2670,3);
     Libderiv->ABCD[38] = int_stack + 2520;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2610,int_stack+2775,int_stack+2745, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3120,3);
     Libderiv->ABCD[35] = int_stack + 2610;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+2850,int_stack+2820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3120, 0.0, zero_stack,3);
     Libderiv->ABCD[34] = int_stack + 2700;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2790,int_stack+2925,int_stack+2895, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3120, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[33] = int_stack + 2790;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2880,int_stack+3000,int_stack+2970, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[32] = int_stack + 2880;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4800,int_stack+3075,int_stack+3045, 0.0, zero_stack, 1.0, int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[31] = int_stack + 4800;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2970,int_stack+3180,int_stack+3150, 1.0, int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[30] = int_stack + 2970;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3060,int_stack+3255,int_stack+3225,3);
     Libderiv->ABCD[26] = int_stack + 3060;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3150,int_stack+3330,int_stack+3300,3);
     Libderiv->ABCD[25] = int_stack + 3150;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3240,int_stack+3405,int_stack+3375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750,3);
     Libderiv->ABCD[23] = int_stack + 3240;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3330,int_stack+3480,int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack,3);
     Libderiv->ABCD[22] = int_stack + 3330;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3420,int_stack+3555,int_stack+3525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[21] = int_stack + 3420;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3510,int_stack+3630,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[20] = int_stack + 3510;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4890,int_stack+3705,int_stack+3675, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[19] = int_stack + 4890;
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+3810,int_stack+3780, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[18] = int_stack + 3600;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3690,int_stack+3885,int_stack+3855,3);
     Libderiv->ABCD[14] = int_stack + 3690;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3780,int_stack+3960,int_stack+3930,3);
     Libderiv->ABCD[13] = int_stack + 3780;
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3870,int_stack+4035,int_stack+4005,3);
     Libderiv->ABCD[12] = int_stack + 3870;

}
