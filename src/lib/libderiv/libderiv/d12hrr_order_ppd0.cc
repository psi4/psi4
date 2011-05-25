#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ppd0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|d0) integrals */

void d12hrr_order_ppd0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][2][11] = int_stack + 0;
 Libderiv->deriv_classes[2][2][10] = int_stack + 36;
 Libderiv->deriv_classes[2][2][9] = int_stack + 72;
 Libderiv->deriv_classes[2][2][8] = int_stack + 108;
 Libderiv->deriv_classes[2][2][7] = int_stack + 144;
 Libderiv->deriv_classes[2][2][6] = int_stack + 180;
 Libderiv->deriv_classes[2][2][2] = int_stack + 216;
 Libderiv->deriv_classes[2][2][1] = int_stack + 252;
 Libderiv->dvrr_classes[1][2] = int_stack + 288;
 Libderiv->deriv_classes[2][2][0] = int_stack + 306;
 Libderiv->deriv2_classes[1][2][143] = int_stack + 342;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 360;
 Libderiv->deriv2_classes[1][2][131] = int_stack + 396;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 414;
 Libderiv->deriv2_classes[1][2][130] = int_stack + 450;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 468;
 Libderiv->deriv2_classes[1][2][119] = int_stack + 504;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 522;
 Libderiv->deriv2_classes[1][2][118] = int_stack + 558;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 576;
 Libderiv->deriv2_classes[1][2][117] = int_stack + 612;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 630;
 Libderiv->deriv2_classes[1][2][107] = int_stack + 666;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 684;
 Libderiv->deriv2_classes[1][2][106] = int_stack + 720;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 738;
 Libderiv->deriv2_classes[1][2][105] = int_stack + 774;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 792;
 Libderiv->deriv2_classes[1][2][104] = int_stack + 828;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 846;
 Libderiv->deriv2_classes[1][2][95] = int_stack + 882;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 900;
 Libderiv->deriv2_classes[1][2][94] = int_stack + 936;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 954;
 Libderiv->deriv2_classes[1][2][93] = int_stack + 990;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 1008;
 Libderiv->deriv2_classes[1][2][92] = int_stack + 1044;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 1062;
 Libderiv->deriv2_classes[1][2][91] = int_stack + 1098;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 1116;
 Libderiv->deriv2_classes[1][2][83] = int_stack + 1152;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 1170;
 Libderiv->deriv2_classes[1][2][82] = int_stack + 1206;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 1224;
 Libderiv->deriv2_classes[1][2][81] = int_stack + 1260;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 1278;
 Libderiv->deriv2_classes[1][2][80] = int_stack + 1314;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 1332;
 Libderiv->deriv2_classes[1][2][79] = int_stack + 1368;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 1386;
 Libderiv->deriv2_classes[1][2][78] = int_stack + 1422;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 1440;
 Libderiv->deriv2_classes[1][2][35] = int_stack + 1476;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 1494;
 Libderiv->deriv2_classes[1][2][34] = int_stack + 1530;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 1548;
 Libderiv->deriv2_classes[1][2][33] = int_stack + 1584;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 1602;
 Libderiv->deriv2_classes[1][2][32] = int_stack + 1638;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 1656;
 Libderiv->deriv2_classes[1][2][31] = int_stack + 1692;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 1710;
 Libderiv->deriv2_classes[1][2][30] = int_stack + 1746;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 1764;
 Libderiv->deriv2_classes[1][2][26] = int_stack + 1800;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 1818;
 Libderiv->deriv2_classes[1][2][23] = int_stack + 1854;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 1872;
 Libderiv->deriv2_classes[1][2][22] = int_stack + 1908;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 1926;
 Libderiv->deriv2_classes[1][2][21] = int_stack + 1962;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 1980;
 Libderiv->deriv2_classes[1][2][20] = int_stack + 2016;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 2034;
 Libderiv->deriv2_classes[1][2][19] = int_stack + 2070;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 2088;
 Libderiv->deriv2_classes[1][2][18] = int_stack + 2124;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 2142;
 Libderiv->deriv2_classes[1][2][14] = int_stack + 2178;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 2196;
 Libderiv->deriv2_classes[1][2][13] = int_stack + 2232;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 2250;
 Libderiv->deriv_classes[1][2][11] = int_stack + 2286;
 Libderiv->deriv2_classes[1][2][11] = int_stack + 2304;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 2322;
 Libderiv->deriv_classes[1][2][10] = int_stack + 2358;
 Libderiv->deriv2_classes[1][2][10] = int_stack + 2376;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 2394;
 Libderiv->deriv_classes[1][2][9] = int_stack + 2430;
 Libderiv->deriv2_classes[1][2][9] = int_stack + 2448;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 2466;
 Libderiv->deriv_classes[1][2][8] = int_stack + 2502;
 Libderiv->deriv2_classes[1][2][8] = int_stack + 2520;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 2538;
 Libderiv->deriv_classes[1][2][7] = int_stack + 2574;
 Libderiv->deriv2_classes[1][2][7] = int_stack + 2592;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 2610;
 Libderiv->deriv_classes[1][2][6] = int_stack + 2646;
 Libderiv->deriv2_classes[1][2][6] = int_stack + 2664;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 2682;
 Libderiv->deriv_classes[1][2][2] = int_stack + 2718;
 Libderiv->deriv2_classes[1][2][2] = int_stack + 2736;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 2754;
 Libderiv->deriv_classes[1][2][1] = int_stack + 2790;
 Libderiv->deriv2_classes[1][2][1] = int_stack + 2808;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 2826;
 Libderiv->deriv_classes[1][2][0] = int_stack + 2862;
 Libderiv->deriv2_classes[1][2][0] = int_stack + 2880;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 2898;
 memset(int_stack,0,23472);

 Libderiv->dvrr_stack = int_stack + 3150;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ppd0(Libderiv, Data);
   Data++;
 }

 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2934,int_stack+0,int_stack+2286,6);
     Libderiv->ABCD[11] = int_stack + 2934;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2988,int_stack+36,int_stack+2358,6);
     Libderiv->ABCD[10] = int_stack + 2988;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+72,int_stack+2430,6);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+54,int_stack+108,int_stack+2502,6);
     Libderiv->ABCD[8] = int_stack + 54;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3042,int_stack+144,int_stack+2574,6);
     Libderiv->ABCD[7] = int_stack + 3042;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+108,int_stack+180,int_stack+2646,6);
     Libderiv->ABCD[6] = int_stack + 108;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+162,int_stack+216,int_stack+2718, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[2] = int_stack + 162;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3096,int_stack+252,int_stack+2790, 0.0, zero_stack, 1.0, int_stack+288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[1] = int_stack + 3096;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+216,int_stack+306,int_stack+2862, 1.0, int_stack+288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[0] = int_stack + 216;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+270,int_stack+360,int_stack+342,6);
     Libderiv->ABCD[155] = int_stack + 270;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+324,int_stack+414,int_stack+396,6);
     Libderiv->ABCD[143] = int_stack + 324;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+378,int_stack+468,int_stack+450,6);
     Libderiv->ABCD[142] = int_stack + 378;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+432,int_stack+522,int_stack+504,6);
     Libderiv->ABCD[131] = int_stack + 432;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+486,int_stack+576,int_stack+558,6);
     Libderiv->ABCD[130] = int_stack + 486;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+540,int_stack+630,int_stack+612,6);
     Libderiv->ABCD[129] = int_stack + 540;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+594,int_stack+684,int_stack+666,6);
     Libderiv->ABCD[119] = int_stack + 594;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+648,int_stack+738,int_stack+720,6);
     Libderiv->ABCD[118] = int_stack + 648;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+702,int_stack+792,int_stack+774,6);
     Libderiv->ABCD[117] = int_stack + 702;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+756,int_stack+846,int_stack+828,6);
     Libderiv->ABCD[116] = int_stack + 756;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+810,int_stack+900,int_stack+882,6);
     Libderiv->ABCD[107] = int_stack + 810;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+864,int_stack+954,int_stack+936,6);
     Libderiv->ABCD[106] = int_stack + 864;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+918,int_stack+1008,int_stack+990,6);
     Libderiv->ABCD[105] = int_stack + 918;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+972,int_stack+1062,int_stack+1044,6);
     Libderiv->ABCD[104] = int_stack + 972;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1026,int_stack+1116,int_stack+1098,6);
     Libderiv->ABCD[103] = int_stack + 1026;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1080,int_stack+1170,int_stack+1152,6);
     Libderiv->ABCD[95] = int_stack + 1080;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1134,int_stack+1224,int_stack+1206,6);
     Libderiv->ABCD[94] = int_stack + 1134;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1188,int_stack+1278,int_stack+1260,6);
     Libderiv->ABCD[93] = int_stack + 1188;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1242,int_stack+1332,int_stack+1314,6);
     Libderiv->ABCD[92] = int_stack + 1242;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1296,int_stack+1386,int_stack+1368,6);
     Libderiv->ABCD[91] = int_stack + 1296;
 /*--- compute (pp|d0) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1350,int_stack+1440,int_stack+1422,6);
     Libderiv->ABCD[90] = int_stack + 1350;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1404,int_stack+1494,int_stack+1476, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[47] = int_stack + 1404;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1458,int_stack+1548,int_stack+1530, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[46] = int_stack + 1458;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1512,int_stack+1602,int_stack+1584, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[45] = int_stack + 1512;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1566,int_stack+1656,int_stack+1638, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[44] = int_stack + 1566;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1620,int_stack+1710,int_stack+1692, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[43] = int_stack + 1620;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1674,int_stack+1764,int_stack+1746, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[42] = int_stack + 1674;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1728,int_stack+1818,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[38] = int_stack + 1728;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1782,int_stack+1872,int_stack+1854, 0.0, zero_stack, 1.0, int_stack+2286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[35] = int_stack + 1782;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1836,int_stack+1926,int_stack+1908, 0.0, zero_stack, 1.0, int_stack+2358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[34] = int_stack + 1836;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1890,int_stack+1980,int_stack+1962, 0.0, zero_stack, 1.0, int_stack+2430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[33] = int_stack + 1890;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1944,int_stack+2034,int_stack+2016, 0.0, zero_stack, 1.0, int_stack+2502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[32] = int_stack + 1944;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1998,int_stack+2088,int_stack+2070, 0.0, zero_stack, 1.0, int_stack+2574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[31] = int_stack + 1998;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2052,int_stack+2142,int_stack+2124, 0.0, zero_stack, 1.0, int_stack+2646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[30] = int_stack + 2052;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2106,int_stack+2196,int_stack+2178, 0.0, zero_stack, 1.0, int_stack+2718, 1.0, int_stack+2790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[26] = int_stack + 2106;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2160,int_stack+2250,int_stack+2232, 0.0, zero_stack, 2.0, int_stack+2790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[25] = int_stack + 2160;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2214,int_stack+2322,int_stack+2304, 1.0, int_stack+2286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[23] = int_stack + 2214;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2268,int_stack+2394,int_stack+2376, 1.0, int_stack+2358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[22] = int_stack + 2268;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2322,int_stack+2466,int_stack+2448, 1.0, int_stack+2430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[21] = int_stack + 2322;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2376,int_stack+2538,int_stack+2520, 1.0, int_stack+2502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[20] = int_stack + 2376;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2430,int_stack+2610,int_stack+2592, 1.0, int_stack+2574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[19] = int_stack + 2430;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2484,int_stack+2682,int_stack+2664, 1.0, int_stack+2646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[18] = int_stack + 2484;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2538,int_stack+2754,int_stack+2736, 1.0, int_stack+2718, 0.0, zero_stack, 1.0, int_stack+2862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[14] = int_stack + 2538;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2592,int_stack+2826,int_stack+2808, 1.0, int_stack+2790, 1.0, int_stack+2862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[13] = int_stack + 2592;
 /*--- compute (pp|d0) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2646,int_stack+2898,int_stack+2880, 2.0, int_stack+2862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[12] = int_stack + 2646;

}
