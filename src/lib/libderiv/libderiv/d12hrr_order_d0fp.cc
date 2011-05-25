#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_d0fp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|fp) integrals */

void d12hrr_order_d0fp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][10] = int_stack + 90;
 Libderiv->deriv_classes[2][4][9] = int_stack + 180;
 Libderiv->deriv_classes[2][4][8] = int_stack + 270;
 Libderiv->deriv_classes[2][4][7] = int_stack + 360;
 Libderiv->dvrr_classes[2][3] = int_stack + 450;
 Libderiv->deriv_classes[2][4][6] = int_stack + 510;
 Libderiv->deriv_classes[2][4][2] = int_stack + 600;
 Libderiv->deriv_classes[2][4][1] = int_stack + 690;
 Libderiv->deriv_classes[2][4][0] = int_stack + 780;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 870;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 930;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1020;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 1080;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 1170;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 1230;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 1320;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 1380;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 1470;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 1530;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 1620;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 1680;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 1770;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 1830;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 1920;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 1980;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 2070;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 2130;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 2220;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 2280;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 2370;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 2430;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 2520;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 2580;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 2670;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 2730;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 2820;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 2880;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 2970;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 3030;
 Libderiv->deriv_classes[2][3][11] = int_stack + 3120;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 3180;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 3240;
 Libderiv->deriv_classes[2][3][10] = int_stack + 3330;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 3390;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 3450;
 Libderiv->deriv_classes[2][3][9] = int_stack + 3540;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 3600;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 3660;
 Libderiv->deriv_classes[2][3][8] = int_stack + 3750;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 3810;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 3870;
 Libderiv->deriv_classes[2][3][7] = int_stack + 3960;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 4020;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 4080;
 Libderiv->deriv_classes[2][3][6] = int_stack + 4170;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 4230;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 4290;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 4380;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 4440;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 4530;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 4590;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 4680;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 4740;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 4830;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 4890;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 4980;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 5040;
 Libderiv->deriv_classes[2][3][2] = int_stack + 5130;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 5190;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 5250;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 5340;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 5400;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 5490;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 5550;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 5640;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 5700;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 5790;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 5850;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 5940;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 6000;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 6090;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 6150;
 Libderiv->deriv_classes[2][3][1] = int_stack + 6240;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 6300;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 6360;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 6450;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 6510;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 6600;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 6660;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 6750;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 6810;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 6900;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 6960;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 7050;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 7110;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 7200;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 7260;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 7350;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 7410;
 Libderiv->deriv_classes[2][3][0] = int_stack + 7500;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 7560;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 7620;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 7710;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 7770;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 7860;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 7920;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 8010;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 8070;
 memset(int_stack,0,65280);

 Libderiv->dvrr_stack = int_stack + 9960;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_d0fp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8160,int_stack+0,int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,6);
     Libderiv->ABCD[11] = int_stack + 8160;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8340,int_stack+90,int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 8340;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+180,int_stack+3540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8520,int_stack+270,int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 8520;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+360,int_stack+3960, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 180;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8700,int_stack+510,int_stack+4170, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 8700;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+600,int_stack+5130,6);
     Libderiv->ABCD[2] = int_stack + 360;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8880,int_stack+690,int_stack+6240,6);
     Libderiv->ABCD[1] = int_stack + 8880;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+780,int_stack+7500,6);
     Libderiv->ABCD[0] = int_stack + 540;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9060,int_stack+930,int_stack+870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3120,6);
     Libderiv->ABCD[155] = int_stack + 9060;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+1080,int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3120, 1.0, int_stack+3330,6);
     Libderiv->ABCD[143] = int_stack + 720;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1230,int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3330, 0.0, zero_stack,6);
     Libderiv->ABCD[142] = int_stack + 900;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+1380,int_stack+1320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3120, 0.0, zero_stack, 1.0, int_stack+3540,6);
     Libderiv->ABCD[131] = int_stack + 1080;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1260,int_stack+1530,int_stack+1470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3330, 1.0, int_stack+3540, 0.0, zero_stack,6);
     Libderiv->ABCD[130] = int_stack + 1260;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1440,int_stack+1680,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3540, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[129] = int_stack + 1440;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9240,int_stack+1830,int_stack+1770, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750,6);
     Libderiv->ABCD[119] = int_stack + 9240;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1620,int_stack+1980,int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3330, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack,6);
     Libderiv->ABCD[118] = int_stack + 1620;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2130,int_stack+2070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3540, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[117] = int_stack + 1800;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1980,int_stack+2280,int_stack+2220, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[116] = int_stack + 1980;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2160,int_stack+2430,int_stack+2370, 0.0, zero_stack, 1.0, int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3960,6);
     Libderiv->ABCD[107] = int_stack + 2160;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2340,int_stack+2580,int_stack+2520, 0.0, zero_stack, 1.0, int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3960, 0.0, zero_stack,6);
     Libderiv->ABCD[106] = int_stack + 2340;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9420,int_stack+2730,int_stack+2670, 0.0, zero_stack, 1.0, int_stack+3540, 0.0, zero_stack, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[105] = int_stack + 9420;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2520,int_stack+2880,int_stack+2820, 0.0, zero_stack, 1.0, int_stack+3750, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[104] = int_stack + 2520;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+3030,int_stack+2970, 0.0, zero_stack, 2.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[103] = int_stack + 2700;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2880,int_stack+3240,int_stack+3180, 1.0, int_stack+3120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4170,6);
     Libderiv->ABCD[95] = int_stack + 2880;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3060,int_stack+3450,int_stack+3390, 1.0, int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4170, 0.0, zero_stack,6);
     Libderiv->ABCD[94] = int_stack + 3060;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3240,int_stack+3660,int_stack+3600, 1.0, int_stack+3540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4170, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[93] = int_stack + 3240;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3420,int_stack+3870,int_stack+3810, 1.0, int_stack+3750, 0.0, zero_stack, 1.0, int_stack+4170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[92] = int_stack + 3420;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4080,int_stack+4020, 1.0, int_stack+3960, 1.0, int_stack+4170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[91] = int_stack + 3600;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3780,int_stack+4290,int_stack+4230, 2.0, int_stack+4170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[90] = int_stack + 3780;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3960,int_stack+4440,int_stack+4380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5130,6);
     Libderiv->ABCD[47] = int_stack + 3960;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4140,int_stack+4590,int_stack+4530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5130, 0.0, zero_stack,6);
     Libderiv->ABCD[46] = int_stack + 4140;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4320,int_stack+4740,int_stack+4680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5130, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[45] = int_stack + 4320;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+4890,int_stack+4830, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[44] = int_stack + 4500;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4680,int_stack+5040,int_stack+4980, 0.0, zero_stack, 1.0, int_stack+5130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[43] = int_stack + 4680;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4860,int_stack+5250,int_stack+5190, 1.0, int_stack+5130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[42] = int_stack + 4860;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5040,int_stack+5400,int_stack+5340,6);
     Libderiv->ABCD[38] = int_stack + 5040;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5220,int_stack+5550,int_stack+5490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6240,6);
     Libderiv->ABCD[35] = int_stack + 5220;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5400,int_stack+5700,int_stack+5640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6240, 0.0, zero_stack,6);
     Libderiv->ABCD[34] = int_stack + 5400;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5580,int_stack+5850,int_stack+5790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6240, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[33] = int_stack + 5580;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5760,int_stack+6000,int_stack+5940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[32] = int_stack + 5760;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9600,int_stack+6150,int_stack+6090, 0.0, zero_stack, 1.0, int_stack+6240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[31] = int_stack + 9600;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5940,int_stack+6360,int_stack+6300, 1.0, int_stack+6240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[30] = int_stack + 5940;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6120,int_stack+6510,int_stack+6450,6);
     Libderiv->ABCD[26] = int_stack + 6120;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6300,int_stack+6660,int_stack+6600,6);
     Libderiv->ABCD[25] = int_stack + 6300;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6480,int_stack+6810,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7500,6);
     Libderiv->ABCD[23] = int_stack + 6480;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6660,int_stack+6960,int_stack+6900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7500, 0.0, zero_stack,6);
     Libderiv->ABCD[22] = int_stack + 6660;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6840,int_stack+7110,int_stack+7050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7500, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[21] = int_stack + 6840;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7020,int_stack+7260,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[20] = int_stack + 7020;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9780,int_stack+7410,int_stack+7350, 0.0, zero_stack, 1.0, int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[19] = int_stack + 9780;
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+7620,int_stack+7560, 1.0, int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[18] = int_stack + 7200;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7380,int_stack+7770,int_stack+7710,6);
     Libderiv->ABCD[14] = int_stack + 7380;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7560,int_stack+7920,int_stack+7860,6);
     Libderiv->ABCD[13] = int_stack + 7560;
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7740,int_stack+8070,int_stack+8010,6);
     Libderiv->ABCD[12] = int_stack + 7740;

}
