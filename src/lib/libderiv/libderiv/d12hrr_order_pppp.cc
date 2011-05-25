#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_pppp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|pp) integrals */

void d12hrr_order_pppp(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->dvrr_classes[2][1] = int_stack + 180;
 Libderiv->deriv_classes[2][2][6] = int_stack + 198;
 Libderiv->deriv_classes[2][2][2] = int_stack + 234;
 Libderiv->deriv_classes[2][2][1] = int_stack + 270;
 Libderiv->dvrr_classes[1][2] = int_stack + 306;
 Libderiv->deriv_classes[2][2][0] = int_stack + 324;
 Libderiv->deriv2_classes[1][1][143] = int_stack + 360;
 Libderiv->deriv2_classes[1][2][143] = int_stack + 369;
 Libderiv->deriv2_classes[2][1][143] = int_stack + 387;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 405;
 Libderiv->deriv2_classes[1][1][131] = int_stack + 441;
 Libderiv->deriv2_classes[1][2][131] = int_stack + 450;
 Libderiv->deriv2_classes[2][1][131] = int_stack + 468;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 486;
 Libderiv->deriv2_classes[1][1][130] = int_stack + 522;
 Libderiv->deriv2_classes[1][2][130] = int_stack + 531;
 Libderiv->deriv2_classes[2][1][130] = int_stack + 549;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 567;
 Libderiv->deriv2_classes[1][1][119] = int_stack + 603;
 Libderiv->deriv2_classes[1][2][119] = int_stack + 612;
 Libderiv->deriv2_classes[2][1][119] = int_stack + 630;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 648;
 Libderiv->deriv2_classes[1][1][118] = int_stack + 684;
 Libderiv->deriv2_classes[1][2][118] = int_stack + 693;
 Libderiv->deriv2_classes[2][1][118] = int_stack + 711;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 729;
 Libderiv->deriv2_classes[1][1][117] = int_stack + 765;
 Libderiv->deriv2_classes[1][2][117] = int_stack + 774;
 Libderiv->deriv2_classes[2][1][117] = int_stack + 792;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 810;
 Libderiv->deriv2_classes[1][1][107] = int_stack + 846;
 Libderiv->deriv2_classes[1][2][107] = int_stack + 855;
 Libderiv->deriv2_classes[2][1][107] = int_stack + 873;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 891;
 Libderiv->deriv2_classes[1][1][106] = int_stack + 927;
 Libderiv->deriv2_classes[1][2][106] = int_stack + 936;
 Libderiv->deriv2_classes[2][1][106] = int_stack + 954;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 972;
 Libderiv->deriv2_classes[1][1][105] = int_stack + 1008;
 Libderiv->deriv2_classes[1][2][105] = int_stack + 1017;
 Libderiv->deriv2_classes[2][1][105] = int_stack + 1035;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 1053;
 Libderiv->deriv2_classes[1][1][104] = int_stack + 1089;
 Libderiv->deriv2_classes[1][2][104] = int_stack + 1098;
 Libderiv->deriv2_classes[2][1][104] = int_stack + 1116;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 1134;
 Libderiv->deriv2_classes[1][1][95] = int_stack + 1170;
 Libderiv->deriv2_classes[1][2][95] = int_stack + 1179;
 Libderiv->deriv2_classes[2][1][95] = int_stack + 1197;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 1215;
 Libderiv->deriv2_classes[1][1][94] = int_stack + 1251;
 Libderiv->deriv2_classes[1][2][94] = int_stack + 1260;
 Libderiv->deriv2_classes[2][1][94] = int_stack + 1278;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 1296;
 Libderiv->deriv2_classes[1][1][93] = int_stack + 1332;
 Libderiv->deriv2_classes[1][2][93] = int_stack + 1341;
 Libderiv->deriv2_classes[2][1][93] = int_stack + 1359;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 1377;
 Libderiv->deriv2_classes[1][1][92] = int_stack + 1413;
 Libderiv->deriv2_classes[1][2][92] = int_stack + 1422;
 Libderiv->deriv2_classes[2][1][92] = int_stack + 1440;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 1458;
 Libderiv->deriv2_classes[1][1][91] = int_stack + 1494;
 Libderiv->deriv2_classes[1][2][91] = int_stack + 1503;
 Libderiv->deriv2_classes[2][1][91] = int_stack + 1521;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 1539;
 Libderiv->deriv2_classes[1][1][83] = int_stack + 1575;
 Libderiv->deriv2_classes[1][2][83] = int_stack + 1584;
 Libderiv->deriv_classes[2][1][11] = int_stack + 1602;
 Libderiv->deriv2_classes[2][1][83] = int_stack + 1620;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 1638;
 Libderiv->deriv2_classes[1][1][82] = int_stack + 1674;
 Libderiv->deriv2_classes[1][2][82] = int_stack + 1683;
 Libderiv->deriv_classes[2][1][10] = int_stack + 1701;
 Libderiv->deriv2_classes[2][1][82] = int_stack + 1719;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 1737;
 Libderiv->deriv2_classes[1][1][81] = int_stack + 1773;
 Libderiv->deriv2_classes[1][2][81] = int_stack + 1782;
 Libderiv->deriv_classes[2][1][9] = int_stack + 1800;
 Libderiv->deriv2_classes[2][1][81] = int_stack + 1818;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 1836;
 Libderiv->deriv2_classes[1][1][80] = int_stack + 1872;
 Libderiv->deriv2_classes[1][2][80] = int_stack + 1881;
 Libderiv->deriv_classes[2][1][8] = int_stack + 1899;
 Libderiv->deriv2_classes[2][1][80] = int_stack + 1917;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 1935;
 Libderiv->deriv2_classes[1][1][79] = int_stack + 1971;
 Libderiv->deriv2_classes[1][2][79] = int_stack + 1980;
 Libderiv->deriv_classes[2][1][7] = int_stack + 1998;
 Libderiv->deriv2_classes[2][1][79] = int_stack + 2016;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 2034;
 Libderiv->deriv2_classes[1][1][78] = int_stack + 2070;
 Libderiv->deriv2_classes[1][2][78] = int_stack + 2079;
 Libderiv->deriv_classes[2][1][6] = int_stack + 2097;
 Libderiv->deriv2_classes[2][1][78] = int_stack + 2115;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 2133;
 Libderiv->deriv2_classes[1][1][35] = int_stack + 2169;
 Libderiv->deriv2_classes[1][2][35] = int_stack + 2178;
 Libderiv->deriv2_classes[2][1][35] = int_stack + 2196;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 2214;
 Libderiv->deriv2_classes[1][1][34] = int_stack + 2250;
 Libderiv->deriv2_classes[1][2][34] = int_stack + 2259;
 Libderiv->deriv2_classes[2][1][34] = int_stack + 2277;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 2295;
 Libderiv->deriv2_classes[1][1][33] = int_stack + 2331;
 Libderiv->deriv2_classes[1][2][33] = int_stack + 2340;
 Libderiv->deriv2_classes[2][1][33] = int_stack + 2358;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 2376;
 Libderiv->deriv2_classes[1][1][32] = int_stack + 2412;
 Libderiv->deriv2_classes[1][2][32] = int_stack + 2421;
 Libderiv->deriv2_classes[2][1][32] = int_stack + 2439;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 2457;
 Libderiv->deriv2_classes[1][1][31] = int_stack + 2493;
 Libderiv->deriv2_classes[1][2][31] = int_stack + 2502;
 Libderiv->deriv2_classes[2][1][31] = int_stack + 2520;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 2538;
 Libderiv->deriv2_classes[1][1][30] = int_stack + 2574;
 Libderiv->deriv2_classes[1][2][30] = int_stack + 2583;
 Libderiv->deriv_classes[2][1][2] = int_stack + 2601;
 Libderiv->deriv2_classes[2][1][30] = int_stack + 2619;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 2637;
 Libderiv->deriv2_classes[1][1][26] = int_stack + 2673;
 Libderiv->deriv2_classes[1][2][26] = int_stack + 2682;
 Libderiv->deriv2_classes[2][1][26] = int_stack + 2700;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 2718;
 Libderiv->deriv2_classes[1][1][23] = int_stack + 2754;
 Libderiv->deriv2_classes[1][2][23] = int_stack + 2763;
 Libderiv->deriv2_classes[2][1][23] = int_stack + 2781;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 2799;
 Libderiv->deriv2_classes[1][1][22] = int_stack + 2835;
 Libderiv->deriv2_classes[1][2][22] = int_stack + 2844;
 Libderiv->deriv2_classes[2][1][22] = int_stack + 2862;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 2880;
 Libderiv->deriv2_classes[1][1][21] = int_stack + 2916;
 Libderiv->deriv2_classes[1][2][21] = int_stack + 2925;
 Libderiv->deriv2_classes[2][1][21] = int_stack + 2943;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 2961;
 Libderiv->deriv2_classes[1][1][20] = int_stack + 2997;
 Libderiv->deriv2_classes[1][2][20] = int_stack + 3006;
 Libderiv->deriv2_classes[2][1][20] = int_stack + 3024;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 3042;
 Libderiv->deriv2_classes[1][1][19] = int_stack + 3078;
 Libderiv->deriv2_classes[1][2][19] = int_stack + 3087;
 Libderiv->deriv2_classes[2][1][19] = int_stack + 3105;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 3123;
 Libderiv->deriv2_classes[1][1][18] = int_stack + 3159;
 Libderiv->deriv2_classes[1][2][18] = int_stack + 3168;
 Libderiv->deriv_classes[2][1][1] = int_stack + 3186;
 Libderiv->deriv2_classes[2][1][18] = int_stack + 3204;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 3222;
 Libderiv->deriv2_classes[1][1][14] = int_stack + 3258;
 Libderiv->deriv2_classes[1][2][14] = int_stack + 3267;
 Libderiv->deriv2_classes[2][1][14] = int_stack + 3285;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 3303;
 Libderiv->deriv2_classes[1][1][13] = int_stack + 3339;
 Libderiv->deriv2_classes[1][2][13] = int_stack + 3348;
 Libderiv->deriv2_classes[2][1][13] = int_stack + 3366;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 3384;
 Libderiv->deriv_classes[1][1][11] = int_stack + 3420;
 Libderiv->deriv_classes[1][2][11] = int_stack + 3429;
 Libderiv->deriv2_classes[1][1][11] = int_stack + 3447;
 Libderiv->deriv2_classes[1][2][11] = int_stack + 3456;
 Libderiv->deriv2_classes[2][1][11] = int_stack + 3474;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 3492;
 Libderiv->deriv_classes[1][1][10] = int_stack + 3528;
 Libderiv->deriv_classes[1][2][10] = int_stack + 3537;
 Libderiv->deriv2_classes[1][1][10] = int_stack + 3555;
 Libderiv->deriv2_classes[1][2][10] = int_stack + 3564;
 Libderiv->deriv2_classes[2][1][10] = int_stack + 3582;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 3600;
 Libderiv->deriv_classes[1][1][9] = int_stack + 3636;
 Libderiv->deriv_classes[1][2][9] = int_stack + 3645;
 Libderiv->deriv2_classes[1][1][9] = int_stack + 3663;
 Libderiv->deriv2_classes[1][2][9] = int_stack + 3672;
 Libderiv->deriv2_classes[2][1][9] = int_stack + 3690;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 3708;
 Libderiv->deriv_classes[1][1][8] = int_stack + 3744;
 Libderiv->deriv_classes[1][2][8] = int_stack + 3753;
 Libderiv->deriv2_classes[1][1][8] = int_stack + 3771;
 Libderiv->deriv2_classes[1][2][8] = int_stack + 3780;
 Libderiv->deriv2_classes[2][1][8] = int_stack + 3798;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 3816;
 Libderiv->deriv_classes[1][1][7] = int_stack + 3852;
 Libderiv->deriv_classes[1][2][7] = int_stack + 3861;
 Libderiv->deriv2_classes[1][1][7] = int_stack + 3879;
 Libderiv->deriv2_classes[1][2][7] = int_stack + 3888;
 Libderiv->deriv2_classes[2][1][7] = int_stack + 3906;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 3924;
 Libderiv->dvrr_classes[1][1] = int_stack + 3960;
 Libderiv->deriv_classes[1][1][6] = int_stack + 3969;
 Libderiv->deriv_classes[1][2][6] = int_stack + 3978;
 Libderiv->deriv2_classes[1][1][6] = int_stack + 3996;
 Libderiv->deriv2_classes[1][2][6] = int_stack + 4005;
 Libderiv->deriv_classes[2][1][0] = int_stack + 4023;
 Libderiv->deriv2_classes[2][1][6] = int_stack + 4041;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 4059;
 Libderiv->deriv_classes[1][1][2] = int_stack + 4095;
 Libderiv->deriv_classes[1][2][2] = int_stack + 4104;
 Libderiv->deriv2_classes[1][1][2] = int_stack + 4122;
 Libderiv->deriv2_classes[1][2][2] = int_stack + 4131;
 Libderiv->deriv2_classes[2][1][2] = int_stack + 4149;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 4167;
 Libderiv->deriv_classes[1][1][1] = int_stack + 4203;
 Libderiv->deriv_classes[1][2][1] = int_stack + 4212;
 Libderiv->deriv2_classes[1][1][1] = int_stack + 4230;
 Libderiv->deriv2_classes[1][2][1] = int_stack + 4239;
 Libderiv->deriv2_classes[2][1][1] = int_stack + 4257;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 4275;
 Libderiv->deriv_classes[1][1][0] = int_stack + 4311;
 Libderiv->deriv_classes[1][2][0] = int_stack + 4320;
 Libderiv->deriv2_classes[1][1][0] = int_stack + 4338;
 Libderiv->deriv2_classes[1][2][0] = int_stack + 4347;
 Libderiv->deriv2_classes[2][1][0] = int_stack + 4365;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 4383;
 memset(int_stack,0,35352);

 Libderiv->dvrr_stack = int_stack + 5337;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_pppp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+4419,int_stack+3429,int_stack+3420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3960,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+4446,int_stack+0,int_stack+1602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+0,int_stack+3537,int_stack+3528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3960, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+4500,int_stack+36,int_stack+1701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+27,int_stack+3645,int_stack+3636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+4554,int_stack+72,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+54,int_stack+3753,int_stack+3744, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+4608,int_stack+108,int_stack+1899, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+81,int_stack+3861,int_stack+3852, 0.0, zero_stack, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+4662,int_stack+144,int_stack+1998, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+108,int_stack+3978,int_stack+3969, 1.0, int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+4716,int_stack+198,int_stack+2097, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+135,int_stack+306,int_stack+3960,3);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+162,int_stack+4104,int_stack+4095,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4770,int_stack+234,int_stack+2601,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+189,int_stack+4212,int_stack+4203,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+216,int_stack+270,int_stack+3186,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+270,int_stack+4320,int_stack+4311,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4824,int_stack+324,int_stack+4023,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+297,int_stack+369,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3420,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+324,int_stack+405,int_stack+387, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1602,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+378,int_stack+450,int_stack+441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3420, 1.0, int_stack+3528,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+405,int_stack+486,int_stack+468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1602, 1.0, int_stack+1701,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+459,int_stack+531,int_stack+522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3528, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+486,int_stack+567,int_stack+549, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1701, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+540,int_stack+612,int_stack+603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3420, 0.0, zero_stack, 1.0, int_stack+3636,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+567,int_stack+648,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1602, 0.0, zero_stack, 1.0, int_stack+1800,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+621,int_stack+693,int_stack+684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3528, 1.0, int_stack+3636, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+648,int_stack+729,int_stack+711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1701, 1.0, int_stack+1800, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+702,int_stack+774,int_stack+765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3636, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+729,int_stack+810,int_stack+792, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1800, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+783,int_stack+855,int_stack+846, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3744,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+810,int_stack+891,int_stack+873, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1602, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1899,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+864,int_stack+936,int_stack+927, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3528, 0.0, zero_stack, 1.0, int_stack+3744, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+891,int_stack+972,int_stack+954, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1701, 0.0, zero_stack, 1.0, int_stack+1899, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+945,int_stack+1017,int_stack+1008, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3636, 1.0, int_stack+3744, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+972,int_stack+1053,int_stack+1035, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1800, 1.0, int_stack+1899, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1026,int_stack+1098,int_stack+1089, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3744, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1053,int_stack+1134,int_stack+1116, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1899, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1107,int_stack+1179,int_stack+1170, 0.0, zero_stack, 1.0, int_stack+3420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3852,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1134,int_stack+1215,int_stack+1197, 0.0, zero_stack, 1.0, int_stack+1602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1998,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1188,int_stack+1260,int_stack+1251, 0.0, zero_stack, 1.0, int_stack+3528, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3852, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1215,int_stack+1296,int_stack+1278, 0.0, zero_stack, 1.0, int_stack+1701, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1998, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1269,int_stack+1341,int_stack+1332, 0.0, zero_stack, 1.0, int_stack+3636, 0.0, zero_stack, 1.0, int_stack+3852, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1296,int_stack+1377,int_stack+1359, 0.0, zero_stack, 1.0, int_stack+1800, 0.0, zero_stack, 1.0, int_stack+1998, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1350,int_stack+1422,int_stack+1413, 0.0, zero_stack, 1.0, int_stack+3744, 1.0, int_stack+3852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1377,int_stack+1458,int_stack+1440, 0.0, zero_stack, 1.0, int_stack+1899, 1.0, int_stack+1998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1431,int_stack+1503,int_stack+1494, 0.0, zero_stack, 2.0, int_stack+3852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1458,int_stack+1539,int_stack+1521, 0.0, zero_stack, 2.0, int_stack+1998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1512,int_stack+1584,int_stack+1575, 1.0, int_stack+3420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3969,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1539,int_stack+1638,int_stack+1620, 1.0, int_stack+1602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2097,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3420,int_stack+1683,int_stack+1674, 1.0, int_stack+3528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3969, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1593,int_stack+1737,int_stack+1719, 1.0, int_stack+1701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2097, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3528,int_stack+1782,int_stack+1773, 1.0, int_stack+3636, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3969, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1647,int_stack+1836,int_stack+1818, 1.0, int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2097, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3636,int_stack+1881,int_stack+1872, 1.0, int_stack+3744, 0.0, zero_stack, 1.0, int_stack+3969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1701,int_stack+1935,int_stack+1917, 1.0, int_stack+1899, 0.0, zero_stack, 1.0, int_stack+2097, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3744,int_stack+1980,int_stack+1971, 1.0, int_stack+3852, 1.0, int_stack+3969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1755,int_stack+2034,int_stack+2016, 1.0, int_stack+1998, 1.0, int_stack+2097, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3852,int_stack+2079,int_stack+2070, 2.0, int_stack+3969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1809,int_stack+2133,int_stack+2115, 2.0, int_stack+2097, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1863,int_stack+2178,int_stack+2169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4095,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1890,int_stack+2214,int_stack+2196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2601,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1944,int_stack+2259,int_stack+2250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4095, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+1971,int_stack+2295,int_stack+2277, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2601, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2025,int_stack+2340,int_stack+2331, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4095, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2052,int_stack+2376,int_stack+2358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2601, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2106,int_stack+2421,int_stack+2412, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4095, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2133,int_stack+2457,int_stack+2439, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2187,int_stack+2502,int_stack+2493, 0.0, zero_stack, 1.0, int_stack+4095, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2214,int_stack+2538,int_stack+2520, 0.0, zero_stack, 1.0, int_stack+2601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2268,int_stack+2583,int_stack+2574, 1.0, int_stack+4095, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2295,int_stack+2637,int_stack+2619, 1.0, int_stack+2601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4095,int_stack+2682,int_stack+2673,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+2349,int_stack+2718,int_stack+2700,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2403,int_stack+2763,int_stack+2754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4203,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2430,int_stack+2799,int_stack+2781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3186,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2484,int_stack+2844,int_stack+2835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4203, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2511,int_stack+2880,int_stack+2862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3186, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2565,int_stack+2925,int_stack+2916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4203, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2592,int_stack+2961,int_stack+2943, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3186, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2646,int_stack+3006,int_stack+2997, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4203, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2673,int_stack+3042,int_stack+3024, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2727,int_stack+3087,int_stack+3078, 0.0, zero_stack, 1.0, int_stack+4203, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2754,int_stack+3123,int_stack+3105, 0.0, zero_stack, 1.0, int_stack+3186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2808,int_stack+3168,int_stack+3159, 1.0, int_stack+4203, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+2835,int_stack+3222,int_stack+3204, 1.0, int_stack+3186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4203,int_stack+3267,int_stack+3258,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+2889,int_stack+3303,int_stack+3285,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+2943,int_stack+3348,int_stack+3339,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+2970,int_stack+3384,int_stack+3366,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3024,int_stack+3456,int_stack+3447, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4311,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3051,int_stack+3492,int_stack+3474, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4023,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3447,int_stack+3564,int_stack+3555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4311, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3474,int_stack+3600,int_stack+3582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4023, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3555,int_stack+3672,int_stack+3663, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4311, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3582,int_stack+3708,int_stack+3690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4023, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3663,int_stack+3780,int_stack+3771, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3690,int_stack+3816,int_stack+3798, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4023, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3771,int_stack+3888,int_stack+3879, 0.0, zero_stack, 1.0, int_stack+4311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3798,int_stack+3924,int_stack+3906, 0.0, zero_stack, 1.0, int_stack+4023, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3879,int_stack+4005,int_stack+3996, 1.0, int_stack+4311, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+3906,int_stack+4059,int_stack+4041, 1.0, int_stack+4023, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4311,int_stack+4131,int_stack+4122,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+3960,int_stack+4167,int_stack+4149,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4122,int_stack+4239,int_stack+4230,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4149,int_stack+4275,int_stack+4257,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4230,int_stack+4347,int_stack+4338,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+4257,int_stack+4383,int_stack+4365,6);
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4338,int_stack+4446,int_stack+4419,9);
     Libderiv->ABCD[11] = int_stack + 4338;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4014,int_stack+4500,int_stack+0,9);
     Libderiv->ABCD[10] = int_stack + 4014;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4446,int_stack+4554,int_stack+27,9);
     Libderiv->ABCD[9] = int_stack + 4446;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4527,int_stack+4608,int_stack+54,9);
     Libderiv->ABCD[8] = int_stack + 4527;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3105,int_stack+4662,int_stack+81,9);
     Libderiv->ABCD[7] = int_stack + 3105;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4608,int_stack+4716,int_stack+108,9);
     Libderiv->ABCD[6] = int_stack + 4608;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4689,int_stack+4770,int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[2] = int_stack + 4689;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3186,int_stack+216,int_stack+189, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[1] = int_stack + 3186;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3267,int_stack+4824,int_stack+270, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[0] = int_stack + 3267;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4770,int_stack+324,int_stack+297,9);
     Libderiv->ABCD[155] = int_stack + 4770;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+297,int_stack+405,int_stack+378,9);
     Libderiv->ABCD[143] = int_stack + 297;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+378,int_stack+486,int_stack+459,9);
     Libderiv->ABCD[142] = int_stack + 378;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+459,int_stack+567,int_stack+540,9);
     Libderiv->ABCD[131] = int_stack + 459;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+540,int_stack+648,int_stack+621,9);
     Libderiv->ABCD[130] = int_stack + 540;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+621,int_stack+729,int_stack+702,9);
     Libderiv->ABCD[129] = int_stack + 621;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+702,int_stack+810,int_stack+783,9);
     Libderiv->ABCD[119] = int_stack + 702;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+783,int_stack+891,int_stack+864,9);
     Libderiv->ABCD[118] = int_stack + 783;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+864,int_stack+972,int_stack+945,9);
     Libderiv->ABCD[117] = int_stack + 864;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+945,int_stack+1053,int_stack+1026,9);
     Libderiv->ABCD[116] = int_stack + 945;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1026,int_stack+1134,int_stack+1107,9);
     Libderiv->ABCD[107] = int_stack + 1026;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1107,int_stack+1215,int_stack+1188,9);
     Libderiv->ABCD[106] = int_stack + 1107;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1188,int_stack+1296,int_stack+1269,9);
     Libderiv->ABCD[105] = int_stack + 1188;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1269,int_stack+1377,int_stack+1350,9);
     Libderiv->ABCD[104] = int_stack + 1269;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1350,int_stack+1458,int_stack+1431,9);
     Libderiv->ABCD[103] = int_stack + 1350;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1431,int_stack+1539,int_stack+1512,9);
     Libderiv->ABCD[95] = int_stack + 1431;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1512,int_stack+1593,int_stack+3420,9);
     Libderiv->ABCD[94] = int_stack + 1512;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4851,int_stack+1647,int_stack+3528,9);
     Libderiv->ABCD[93] = int_stack + 4851;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1593,int_stack+1701,int_stack+3636,9);
     Libderiv->ABCD[92] = int_stack + 1593;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1674,int_stack+1755,int_stack+3744,9);
     Libderiv->ABCD[91] = int_stack + 1674;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3348,int_stack+1809,int_stack+3852,9);
     Libderiv->ABCD[90] = int_stack + 3348;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1755,int_stack+1890,int_stack+1863, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4419, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[47] = int_stack + 1755;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1836,int_stack+1971,int_stack+1944, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[46] = int_stack + 1836;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1917,int_stack+2052,int_stack+2025, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[45] = int_stack + 1917;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1998,int_stack+2133,int_stack+2106, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[44] = int_stack + 1998;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2079,int_stack+2214,int_stack+2187, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[43] = int_stack + 2079;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2160,int_stack+2295,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[42] = int_stack + 2160;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2241,int_stack+2349,int_stack+4095, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[38] = int_stack + 2241;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2322,int_stack+2430,int_stack+2403, 0.0, zero_stack, 1.0, int_stack+4419, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[35] = int_stack + 2322;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2403,int_stack+2511,int_stack+2484, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[34] = int_stack + 2403;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2484,int_stack+2592,int_stack+2565, 0.0, zero_stack, 1.0, int_stack+27, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[33] = int_stack + 2484;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2565,int_stack+2673,int_stack+2646, 0.0, zero_stack, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[32] = int_stack + 2565;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2646,int_stack+2754,int_stack+2727, 0.0, zero_stack, 1.0, int_stack+81, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[31] = int_stack + 2646;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2727,int_stack+2835,int_stack+2808, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[30] = int_stack + 2727;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2808,int_stack+2889,int_stack+4203, 0.0, zero_stack, 1.0, int_stack+162, 1.0, int_stack+189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[26] = int_stack + 2808;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4932,int_stack+2970,int_stack+2943, 0.0, zero_stack, 2.0, int_stack+189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[25] = int_stack + 4932;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2889,int_stack+3051,int_stack+3024, 1.0, int_stack+4419, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[23] = int_stack + 2889;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2970,int_stack+3474,int_stack+3447, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[22] = int_stack + 2970;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5013,int_stack+3582,int_stack+3555, 1.0, int_stack+27, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[21] = int_stack + 5013;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5094,int_stack+3690,int_stack+3663, 1.0, int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[20] = int_stack + 5094;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+3798,int_stack+3771, 1.0, int_stack+81, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[19] = int_stack + 0;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5175,int_stack+3906,int_stack+3879, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[18] = int_stack + 5175;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+81,int_stack+3960,int_stack+4311, 1.0, int_stack+162, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[14] = int_stack + 81;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5256,int_stack+4149,int_stack+4122, 1.0, int_stack+189, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[13] = int_stack + 5256;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+162,int_stack+4257,int_stack+4230, 2.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[12] = int_stack + 162;

}
