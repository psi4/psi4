#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_d0dp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|dp) integrals */

void d12hrr_order_d0dp(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[2][3][0] = int_stack + 516;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 576;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 612;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 672;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 708;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 768;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 804;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 864;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 900;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 960;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 996;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 1056;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 1092;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 1152;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 1188;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 1248;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 1284;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 1344;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 1380;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 1440;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 1476;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 1536;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 1572;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 1632;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 1668;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 1728;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 1764;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 1824;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 1860;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 1920;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 1956;
 Libderiv->deriv_classes[2][2][11] = int_stack + 2016;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 2052;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 2088;
 Libderiv->deriv_classes[2][2][10] = int_stack + 2148;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 2184;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 2220;
 Libderiv->deriv_classes[2][2][9] = int_stack + 2280;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 2316;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 2352;
 Libderiv->deriv_classes[2][2][8] = int_stack + 2412;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 2448;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 2484;
 Libderiv->deriv_classes[2][2][7] = int_stack + 2544;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 2580;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 2616;
 Libderiv->deriv_classes[2][2][6] = int_stack + 2676;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 2712;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 2748;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 2808;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 2844;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 2904;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 2940;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 3000;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 3036;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 3096;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 3132;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 3192;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 3228;
 Libderiv->deriv_classes[2][2][2] = int_stack + 3288;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 3324;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 3360;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 3420;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 3456;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 3516;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 3552;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 3612;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 3648;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 3708;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 3744;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 3804;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 3840;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 3900;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 3936;
 Libderiv->deriv_classes[2][2][1] = int_stack + 3996;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 4032;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 4068;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 4128;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 4164;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 4224;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 4260;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 4320;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 4356;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 4416;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 4452;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 4512;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 4548;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 4608;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 4644;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 4704;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 4740;
 Libderiv->deriv_classes[2][2][0] = int_stack + 4800;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 4836;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 4872;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 4932;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 4968;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 5028;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 5064;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 5124;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 5160;
 memset(int_stack,0,41760);

 Libderiv->dvrr_stack = int_stack + 5976;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_d0dp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5220,int_stack+0,int_stack+2016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300,6);
     Libderiv->ABCD[11] = int_stack + 5220;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5328,int_stack+60,int_stack+2148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 5328;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+120,int_stack+2280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5436,int_stack+180,int_stack+2412, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 5436;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+240,int_stack+2544, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 108;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5544,int_stack+336,int_stack+2676, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 5544;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+396,int_stack+3288,6);
     Libderiv->ABCD[2] = int_stack + 216;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+456,int_stack+3996,6);
     Libderiv->ABCD[1] = int_stack + 324;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+5652,int_stack+516,int_stack+4800,6);
     Libderiv->ABCD[0] = int_stack + 5652;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+612,int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2016,6);
     Libderiv->ABCD[155] = int_stack + 432;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+708,int_stack+672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2016, 1.0, int_stack+2148,6);
     Libderiv->ABCD[143] = int_stack + 540;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+804,int_stack+768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2148, 0.0, zero_stack,6);
     Libderiv->ABCD[142] = int_stack + 648;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+756,int_stack+900,int_stack+864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2016, 0.0, zero_stack, 1.0, int_stack+2280,6);
     Libderiv->ABCD[131] = int_stack + 756;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5760,int_stack+996,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2148, 1.0, int_stack+2280, 0.0, zero_stack,6);
     Libderiv->ABCD[130] = int_stack + 5760;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+864,int_stack+1092,int_stack+1056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2280, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[129] = int_stack + 864;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+972,int_stack+1188,int_stack+1152, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2016, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2412,6);
     Libderiv->ABCD[119] = int_stack + 972;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1284,int_stack+1248, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2148, 0.0, zero_stack, 1.0, int_stack+2412, 0.0, zero_stack,6);
     Libderiv->ABCD[118] = int_stack + 1080;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1188,int_stack+1380,int_stack+1344, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2280, 1.0, int_stack+2412, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[117] = int_stack + 1188;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1296,int_stack+1476,int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2412, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[116] = int_stack + 1296;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1404,int_stack+1572,int_stack+1536, 0.0, zero_stack, 1.0, int_stack+2016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2544,6);
     Libderiv->ABCD[107] = int_stack + 1404;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1512,int_stack+1668,int_stack+1632, 0.0, zero_stack, 1.0, int_stack+2148, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2544, 0.0, zero_stack,6);
     Libderiv->ABCD[106] = int_stack + 1512;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1620,int_stack+1764,int_stack+1728, 0.0, zero_stack, 1.0, int_stack+2280, 0.0, zero_stack, 1.0, int_stack+2544, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[105] = int_stack + 1620;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5868,int_stack+1860,int_stack+1824, 0.0, zero_stack, 1.0, int_stack+2412, 1.0, int_stack+2544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[104] = int_stack + 5868;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1728,int_stack+1956,int_stack+1920, 0.0, zero_stack, 2.0, int_stack+2544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[103] = int_stack + 1728;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1836,int_stack+2088,int_stack+2052, 1.0, int_stack+2016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2676,6);
     Libderiv->ABCD[95] = int_stack + 1836;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1944,int_stack+2220,int_stack+2184, 1.0, int_stack+2148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2676, 0.0, zero_stack,6);
     Libderiv->ABCD[94] = int_stack + 1944;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2052,int_stack+2352,int_stack+2316, 1.0, int_stack+2280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2676, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[93] = int_stack + 2052;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2160,int_stack+2484,int_stack+2448, 1.0, int_stack+2412, 0.0, zero_stack, 1.0, int_stack+2676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[92] = int_stack + 2160;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2268,int_stack+2616,int_stack+2580, 1.0, int_stack+2544, 1.0, int_stack+2676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[91] = int_stack + 2268;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2376,int_stack+2748,int_stack+2712, 2.0, int_stack+2676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[90] = int_stack + 2376;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2484,int_stack+2844,int_stack+2808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3288,6);
     Libderiv->ABCD[47] = int_stack + 2484;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2592,int_stack+2940,int_stack+2904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3288, 0.0, zero_stack,6);
     Libderiv->ABCD[46] = int_stack + 2592;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+3036,int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3288, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[45] = int_stack + 2700;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2808,int_stack+3132,int_stack+3096, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[44] = int_stack + 2808;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2916,int_stack+3228,int_stack+3192, 0.0, zero_stack, 1.0, int_stack+3288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[43] = int_stack + 2916;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3024,int_stack+3360,int_stack+3324, 1.0, int_stack+3288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[42] = int_stack + 3024;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+3132,int_stack+3456,int_stack+3420,6);
     Libderiv->ABCD[38] = int_stack + 3132;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3240,int_stack+3552,int_stack+3516, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3996,6);
     Libderiv->ABCD[35] = int_stack + 3240;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3348,int_stack+3648,int_stack+3612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3996, 0.0, zero_stack,6);
     Libderiv->ABCD[34] = int_stack + 3348;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3456,int_stack+3744,int_stack+3708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3996, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[33] = int_stack + 3456;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3564,int_stack+3840,int_stack+3804, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[32] = int_stack + 3564;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3672,int_stack+3936,int_stack+3900, 0.0, zero_stack, 1.0, int_stack+3996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[31] = int_stack + 3672;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3780,int_stack+4068,int_stack+4032, 1.0, int_stack+3996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[30] = int_stack + 3780;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+3888,int_stack+4164,int_stack+4128,6);
     Libderiv->ABCD[26] = int_stack + 3888;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+3996,int_stack+4260,int_stack+4224,6);
     Libderiv->ABCD[25] = int_stack + 3996;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4104,int_stack+4356,int_stack+4320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800,6);
     Libderiv->ABCD[23] = int_stack + 4104;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4212,int_stack+4452,int_stack+4416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack,6);
     Libderiv->ABCD[22] = int_stack + 4212;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4320,int_stack+4548,int_stack+4512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[21] = int_stack + 4320;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4428,int_stack+4644,int_stack+4608, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[20] = int_stack + 4428;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4536,int_stack+4740,int_stack+4704, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[19] = int_stack + 4536;
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4644,int_stack+4872,int_stack+4836, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[18] = int_stack + 4644;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4752,int_stack+4968,int_stack+4932,6);
     Libderiv->ABCD[14] = int_stack + 4752;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4860,int_stack+5064,int_stack+5028,6);
     Libderiv->ABCD[13] = int_stack + 4860;
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+4968,int_stack+5160,int_stack+5124,6);
     Libderiv->ABCD[12] = int_stack + 4968;

}
