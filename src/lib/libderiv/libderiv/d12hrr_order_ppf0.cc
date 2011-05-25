#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ppf0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|f0) integrals */

void d12hrr_order_ppf0(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[2][3][6] = int_stack + 300;
 Libderiv->deriv_classes[2][3][2] = int_stack + 360;
 Libderiv->deriv_classes[2][3][1] = int_stack + 420;
 Libderiv->dvrr_classes[1][3] = int_stack + 480;
 Libderiv->deriv_classes[2][3][0] = int_stack + 510;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 570;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 600;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 660;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 690;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 750;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 780;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 840;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 870;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 930;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 960;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 1020;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 1050;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 1110;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 1140;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 1200;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 1230;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 1290;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 1320;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 1380;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 1410;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 1470;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 1500;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 1560;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 1590;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 1650;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 1680;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 1740;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 1770;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 1830;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 1860;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 1920;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 1950;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 2010;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 2040;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 2100;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 2130;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 2190;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 2220;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 2280;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 2310;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 2370;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 2400;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 2460;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 2490;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 2550;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 2580;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 2640;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 2670;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 2730;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 2760;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 2820;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 2850;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 2910;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 2940;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 3000;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 3030;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 3090;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 3120;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 3180;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 3210;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 3270;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 3300;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 3360;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 3390;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 3450;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 3480;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 3540;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 3570;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 3630;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 3660;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 3720;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 3750;
 Libderiv->deriv_classes[1][3][11] = int_stack + 3810;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 3840;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 3870;
 Libderiv->deriv_classes[1][3][10] = int_stack + 3930;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 3960;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 3990;
 Libderiv->deriv_classes[1][3][9] = int_stack + 4050;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 4080;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 4110;
 Libderiv->deriv_classes[1][3][8] = int_stack + 4170;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 4200;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 4230;
 Libderiv->deriv_classes[1][3][7] = int_stack + 4290;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 4320;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 4350;
 Libderiv->deriv_classes[1][3][6] = int_stack + 4410;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 4440;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 4470;
 Libderiv->deriv_classes[1][3][2] = int_stack + 4530;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 4560;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 4590;
 Libderiv->deriv_classes[1][3][1] = int_stack + 4650;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 4680;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 4710;
 Libderiv->deriv_classes[1][3][0] = int_stack + 4770;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 4800;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 4830;
 memset(int_stack,0,39120);

 Libderiv->dvrr_stack = int_stack + 5250;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ppf0(Libderiv, Data);
   Data++;
 }

 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4890,int_stack+0,int_stack+3810,10);
     Libderiv->ABCD[11] = int_stack + 4890;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4980,int_stack+60,int_stack+3930,10);
     Libderiv->ABCD[10] = int_stack + 4980;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+120,int_stack+4050,10);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+90,int_stack+180,int_stack+4170,10);
     Libderiv->ABCD[8] = int_stack + 90;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+5070,int_stack+240,int_stack+4290,10);
     Libderiv->ABCD[7] = int_stack + 5070;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+180,int_stack+300,int_stack+4410,10);
     Libderiv->ABCD[6] = int_stack + 180;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+270,int_stack+360,int_stack+4530, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[2] = int_stack + 270;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5160,int_stack+420,int_stack+4650, 0.0, zero_stack, 1.0, int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[1] = int_stack + 5160;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+360,int_stack+510,int_stack+4770, 1.0, int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[0] = int_stack + 360;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+450,int_stack+600,int_stack+570,10);
     Libderiv->ABCD[155] = int_stack + 450;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+540,int_stack+690,int_stack+660,10);
     Libderiv->ABCD[143] = int_stack + 540;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+630,int_stack+780,int_stack+750,10);
     Libderiv->ABCD[142] = int_stack + 630;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+720,int_stack+870,int_stack+840,10);
     Libderiv->ABCD[131] = int_stack + 720;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+810,int_stack+960,int_stack+930,10);
     Libderiv->ABCD[130] = int_stack + 810;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+900,int_stack+1050,int_stack+1020,10);
     Libderiv->ABCD[129] = int_stack + 900;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+990,int_stack+1140,int_stack+1110,10);
     Libderiv->ABCD[119] = int_stack + 990;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1080,int_stack+1230,int_stack+1200,10);
     Libderiv->ABCD[118] = int_stack + 1080;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1170,int_stack+1320,int_stack+1290,10);
     Libderiv->ABCD[117] = int_stack + 1170;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1260,int_stack+1410,int_stack+1380,10);
     Libderiv->ABCD[116] = int_stack + 1260;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1350,int_stack+1500,int_stack+1470,10);
     Libderiv->ABCD[107] = int_stack + 1350;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1440,int_stack+1590,int_stack+1560,10);
     Libderiv->ABCD[106] = int_stack + 1440;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1530,int_stack+1680,int_stack+1650,10);
     Libderiv->ABCD[105] = int_stack + 1530;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1620,int_stack+1770,int_stack+1740,10);
     Libderiv->ABCD[104] = int_stack + 1620;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1710,int_stack+1860,int_stack+1830,10);
     Libderiv->ABCD[103] = int_stack + 1710;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1800,int_stack+1950,int_stack+1920,10);
     Libderiv->ABCD[95] = int_stack + 1800;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1890,int_stack+2040,int_stack+2010,10);
     Libderiv->ABCD[94] = int_stack + 1890;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1980,int_stack+2130,int_stack+2100,10);
     Libderiv->ABCD[93] = int_stack + 1980;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2070,int_stack+2220,int_stack+2190,10);
     Libderiv->ABCD[92] = int_stack + 2070;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2160,int_stack+2310,int_stack+2280,10);
     Libderiv->ABCD[91] = int_stack + 2160;
 /*--- compute (pp|f0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2250,int_stack+2400,int_stack+2370,10);
     Libderiv->ABCD[90] = int_stack + 2250;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2340,int_stack+2490,int_stack+2460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[47] = int_stack + 2340;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2430,int_stack+2580,int_stack+2550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[46] = int_stack + 2430;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2520,int_stack+2670,int_stack+2640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[45] = int_stack + 2520;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2610,int_stack+2760,int_stack+2730, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[44] = int_stack + 2610;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2700,int_stack+2850,int_stack+2820, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[43] = int_stack + 2700;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2790,int_stack+2940,int_stack+2910, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[42] = int_stack + 2790;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2880,int_stack+3030,int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[38] = int_stack + 2880;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2970,int_stack+3120,int_stack+3090, 0.0, zero_stack, 1.0, int_stack+3810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[35] = int_stack + 2970;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3060,int_stack+3210,int_stack+3180, 0.0, zero_stack, 1.0, int_stack+3930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[34] = int_stack + 3060;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3150,int_stack+3300,int_stack+3270, 0.0, zero_stack, 1.0, int_stack+4050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[33] = int_stack + 3150;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3240,int_stack+3390,int_stack+3360, 0.0, zero_stack, 1.0, int_stack+4170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[32] = int_stack + 3240;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3330,int_stack+3480,int_stack+3450, 0.0, zero_stack, 1.0, int_stack+4290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[31] = int_stack + 3330;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3420,int_stack+3570,int_stack+3540, 0.0, zero_stack, 1.0, int_stack+4410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[30] = int_stack + 3420;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3510,int_stack+3660,int_stack+3630, 0.0, zero_stack, 1.0, int_stack+4530, 1.0, int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[26] = int_stack + 3510;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3600,int_stack+3750,int_stack+3720, 0.0, zero_stack, 2.0, int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[25] = int_stack + 3600;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3690,int_stack+3870,int_stack+3840, 1.0, int_stack+3810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[23] = int_stack + 3690;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3780,int_stack+3990,int_stack+3960, 1.0, int_stack+3930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[22] = int_stack + 3780;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3870,int_stack+4110,int_stack+4080, 1.0, int_stack+4050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[21] = int_stack + 3870;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3960,int_stack+4230,int_stack+4200, 1.0, int_stack+4170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[20] = int_stack + 3960;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4050,int_stack+4350,int_stack+4320, 1.0, int_stack+4290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[19] = int_stack + 4050;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4140,int_stack+4470,int_stack+4440, 1.0, int_stack+4410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[18] = int_stack + 4140;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4230,int_stack+4590,int_stack+4560, 1.0, int_stack+4530, 0.0, zero_stack, 1.0, int_stack+4770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[14] = int_stack + 4230;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4320,int_stack+4710,int_stack+4680, 1.0, int_stack+4650, 1.0, int_stack+4770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[13] = int_stack + 4320;
 /*--- compute (pp|f0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4410,int_stack+4830,int_stack+4800, 2.0, int_stack+4770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[12] = int_stack + 4410;

}
