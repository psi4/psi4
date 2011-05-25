#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|fd) integrals */

void d12hrr_order_p0fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][5][11] = int_stack + 0;
 Libderiv->deriv_classes[1][5][10] = int_stack + 63;
 Libderiv->deriv_classes[1][5][9] = int_stack + 126;
 Libderiv->deriv_classes[1][5][8] = int_stack + 189;
 Libderiv->deriv_classes[1][5][7] = int_stack + 252;
 Libderiv->dvrr_classes[1][4] = int_stack + 315;
 Libderiv->deriv_classes[1][5][6] = int_stack + 360;
 Libderiv->deriv_classes[1][5][2] = int_stack + 423;
 Libderiv->deriv_classes[1][5][1] = int_stack + 486;
 Libderiv->deriv_classes[1][5][0] = int_stack + 549;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 612;
 Libderiv->deriv2_classes[1][4][143] = int_stack + 642;
 Libderiv->deriv2_classes[1][5][143] = int_stack + 687;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 750;
 Libderiv->deriv2_classes[1][4][131] = int_stack + 780;
 Libderiv->deriv2_classes[1][5][131] = int_stack + 825;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 888;
 Libderiv->deriv2_classes[1][4][130] = int_stack + 918;
 Libderiv->deriv2_classes[1][5][130] = int_stack + 963;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 1026;
 Libderiv->deriv2_classes[1][4][119] = int_stack + 1056;
 Libderiv->deriv2_classes[1][5][119] = int_stack + 1101;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 1164;
 Libderiv->deriv2_classes[1][4][118] = int_stack + 1194;
 Libderiv->deriv2_classes[1][5][118] = int_stack + 1239;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 1302;
 Libderiv->deriv2_classes[1][4][117] = int_stack + 1332;
 Libderiv->deriv2_classes[1][5][117] = int_stack + 1377;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 1440;
 Libderiv->deriv2_classes[1][4][107] = int_stack + 1470;
 Libderiv->deriv2_classes[1][5][107] = int_stack + 1515;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 1578;
 Libderiv->deriv2_classes[1][4][106] = int_stack + 1608;
 Libderiv->deriv2_classes[1][5][106] = int_stack + 1653;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 1716;
 Libderiv->deriv2_classes[1][4][105] = int_stack + 1746;
 Libderiv->deriv2_classes[1][5][105] = int_stack + 1791;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 1854;
 Libderiv->deriv2_classes[1][4][104] = int_stack + 1884;
 Libderiv->deriv2_classes[1][5][104] = int_stack + 1929;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 1992;
 Libderiv->deriv2_classes[1][4][95] = int_stack + 2022;
 Libderiv->deriv2_classes[1][5][95] = int_stack + 2067;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 2130;
 Libderiv->deriv2_classes[1][4][94] = int_stack + 2160;
 Libderiv->deriv2_classes[1][5][94] = int_stack + 2205;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 2268;
 Libderiv->deriv2_classes[1][4][93] = int_stack + 2298;
 Libderiv->deriv2_classes[1][5][93] = int_stack + 2343;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 2406;
 Libderiv->deriv2_classes[1][4][92] = int_stack + 2436;
 Libderiv->deriv2_classes[1][5][92] = int_stack + 2481;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 2544;
 Libderiv->deriv2_classes[1][4][91] = int_stack + 2574;
 Libderiv->deriv2_classes[1][5][91] = int_stack + 2619;
 Libderiv->deriv_classes[1][3][11] = int_stack + 2682;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 2712;
 Libderiv->deriv_classes[1][4][11] = int_stack + 2742;
 Libderiv->deriv2_classes[1][4][83] = int_stack + 2787;
 Libderiv->deriv2_classes[1][5][83] = int_stack + 2832;
 Libderiv->deriv_classes[1][3][10] = int_stack + 2895;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 2925;
 Libderiv->deriv_classes[1][4][10] = int_stack + 2955;
 Libderiv->deriv2_classes[1][4][82] = int_stack + 3000;
 Libderiv->deriv2_classes[1][5][82] = int_stack + 3045;
 Libderiv->deriv_classes[1][3][9] = int_stack + 3108;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 3138;
 Libderiv->deriv_classes[1][4][9] = int_stack + 3168;
 Libderiv->deriv2_classes[1][4][81] = int_stack + 3213;
 Libderiv->deriv2_classes[1][5][81] = int_stack + 3258;
 Libderiv->deriv_classes[1][3][8] = int_stack + 3321;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 3351;
 Libderiv->deriv_classes[1][4][8] = int_stack + 3381;
 Libderiv->deriv2_classes[1][4][80] = int_stack + 3426;
 Libderiv->deriv2_classes[1][5][80] = int_stack + 3471;
 Libderiv->deriv_classes[1][3][7] = int_stack + 3534;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 3564;
 Libderiv->deriv_classes[1][4][7] = int_stack + 3594;
 Libderiv->deriv2_classes[1][4][79] = int_stack + 3639;
 Libderiv->deriv2_classes[1][5][79] = int_stack + 3684;
 Libderiv->dvrr_classes[1][3] = int_stack + 3747;
 Libderiv->deriv_classes[1][3][6] = int_stack + 3777;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 3807;
 Libderiv->deriv_classes[1][4][6] = int_stack + 3837;
 Libderiv->deriv2_classes[1][4][78] = int_stack + 3882;
 Libderiv->deriv2_classes[1][5][78] = int_stack + 3927;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 3990;
 Libderiv->deriv2_classes[1][4][35] = int_stack + 4020;
 Libderiv->deriv2_classes[1][5][35] = int_stack + 4065;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 4128;
 Libderiv->deriv2_classes[1][4][34] = int_stack + 4158;
 Libderiv->deriv2_classes[1][5][34] = int_stack + 4203;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 4266;
 Libderiv->deriv2_classes[1][4][33] = int_stack + 4296;
 Libderiv->deriv2_classes[1][5][33] = int_stack + 4341;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 4404;
 Libderiv->deriv2_classes[1][4][32] = int_stack + 4434;
 Libderiv->deriv2_classes[1][5][32] = int_stack + 4479;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 4542;
 Libderiv->deriv2_classes[1][4][31] = int_stack + 4572;
 Libderiv->deriv2_classes[1][5][31] = int_stack + 4617;
 Libderiv->deriv_classes[1][3][2] = int_stack + 4680;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 4710;
 Libderiv->deriv_classes[1][4][2] = int_stack + 4740;
 Libderiv->deriv2_classes[1][4][30] = int_stack + 4785;
 Libderiv->deriv2_classes[1][5][30] = int_stack + 4830;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 4893;
 Libderiv->deriv2_classes[1][4][26] = int_stack + 4923;
 Libderiv->deriv2_classes[1][5][26] = int_stack + 4968;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 5031;
 Libderiv->deriv2_classes[1][4][23] = int_stack + 5061;
 Libderiv->deriv2_classes[1][5][23] = int_stack + 5106;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 5169;
 Libderiv->deriv2_classes[1][4][22] = int_stack + 5199;
 Libderiv->deriv2_classes[1][5][22] = int_stack + 5244;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 5307;
 Libderiv->deriv2_classes[1][4][21] = int_stack + 5337;
 Libderiv->deriv2_classes[1][5][21] = int_stack + 5382;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 5445;
 Libderiv->deriv2_classes[1][4][20] = int_stack + 5475;
 Libderiv->deriv2_classes[1][5][20] = int_stack + 5520;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 5583;
 Libderiv->deriv2_classes[1][4][19] = int_stack + 5613;
 Libderiv->deriv2_classes[1][5][19] = int_stack + 5658;
 Libderiv->deriv_classes[1][3][1] = int_stack + 5721;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 5751;
 Libderiv->deriv_classes[1][4][1] = int_stack + 5781;
 Libderiv->deriv2_classes[1][4][18] = int_stack + 5826;
 Libderiv->deriv2_classes[1][5][18] = int_stack + 5871;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 5934;
 Libderiv->deriv2_classes[1][4][14] = int_stack + 5964;
 Libderiv->deriv2_classes[1][5][14] = int_stack + 6009;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 6072;
 Libderiv->deriv2_classes[1][4][13] = int_stack + 6102;
 Libderiv->deriv2_classes[1][5][13] = int_stack + 6147;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 6210;
 Libderiv->deriv2_classes[1][4][11] = int_stack + 6240;
 Libderiv->deriv2_classes[1][5][11] = int_stack + 6285;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 6348;
 Libderiv->deriv2_classes[1][4][10] = int_stack + 6378;
 Libderiv->deriv2_classes[1][5][10] = int_stack + 6423;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 6486;
 Libderiv->deriv2_classes[1][4][9] = int_stack + 6516;
 Libderiv->deriv2_classes[1][5][9] = int_stack + 6561;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 6624;
 Libderiv->deriv2_classes[1][4][8] = int_stack + 6654;
 Libderiv->deriv2_classes[1][5][8] = int_stack + 6699;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 6762;
 Libderiv->deriv2_classes[1][4][7] = int_stack + 6792;
 Libderiv->deriv2_classes[1][5][7] = int_stack + 6837;
 Libderiv->deriv_classes[1][3][0] = int_stack + 6900;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 6930;
 Libderiv->deriv_classes[1][4][0] = int_stack + 6960;
 Libderiv->deriv2_classes[1][4][6] = int_stack + 7005;
 Libderiv->deriv2_classes[1][5][6] = int_stack + 7050;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 7113;
 Libderiv->deriv2_classes[1][4][2] = int_stack + 7143;
 Libderiv->deriv2_classes[1][5][2] = int_stack + 7188;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 7251;
 Libderiv->deriv2_classes[1][4][1] = int_stack + 7281;
 Libderiv->deriv2_classes[1][5][1] = int_stack + 7326;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 7389;
 Libderiv->deriv2_classes[1][4][0] = int_stack + 7419;
 Libderiv->deriv2_classes[1][5][0] = int_stack + 7464;
 memset(int_stack,0,60216);

 Libderiv->dvrr_stack = int_stack + 13152;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7527,int_stack+315,int_stack+3747,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7617,int_stack+2742,int_stack+2682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3747,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7707,int_stack+0,int_stack+2742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7842,int_stack+2955,int_stack+2895, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3747, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7932,int_stack+63,int_stack+2955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+3168,int_stack+3108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3747, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8067,int_stack+126,int_stack+3168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+3381,int_stack+3321, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3747, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8202,int_stack+189,int_stack+3381, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8337,int_stack+3594,int_stack+3534, 0.0, zero_stack, 1.0, int_stack+3747, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8427,int_stack+252,int_stack+3594, 0.0, zero_stack, 1.0, int_stack+315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+3837,int_stack+3777, 1.0, int_stack+3747, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8562,int_stack+360,int_stack+3837, 1.0, int_stack+315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+270,int_stack+4740,int_stack+4680,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8697,int_stack+423,int_stack+4740,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+5781,int_stack+5721,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8832,int_stack+486,int_stack+5781,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+6960,int_stack+6900,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8967,int_stack+549,int_stack+6960,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9102,int_stack+642,int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2682,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9192,int_stack+687,int_stack+642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2742,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+780,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2682, 1.0, int_stack+2895,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+630,int_stack+825,int_stack+780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2742, 1.0, int_stack+2955,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+765,int_stack+918,int_stack+888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2895, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9327,int_stack+963,int_stack+918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+2955, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+855,int_stack+1056,int_stack+1026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2682, 0.0, zero_stack, 1.0, int_stack+3108,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9462,int_stack+1101,int_stack+1056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2742, 0.0, zero_stack, 1.0, int_stack+3168,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+945,int_stack+1194,int_stack+1164, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2895, 1.0, int_stack+3108, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1035,int_stack+1239,int_stack+1194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2955, 1.0, int_stack+3168, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1170,int_stack+1332,int_stack+1302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3108, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9597,int_stack+1377,int_stack+1332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3168, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1260,int_stack+1470,int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2682, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3321,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9732,int_stack+1515,int_stack+1470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3381,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1350,int_stack+1608,int_stack+1578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2895, 0.0, zero_stack, 1.0, int_stack+3321, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1440,int_stack+1653,int_stack+1608, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2955, 0.0, zero_stack, 1.0, int_stack+3381, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1575,int_stack+1746,int_stack+1716, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3108, 1.0, int_stack+3321, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9867,int_stack+1791,int_stack+1746, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3168, 1.0, int_stack+3381, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1665,int_stack+1884,int_stack+1854, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10002,int_stack+1929,int_stack+1884, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3381, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1755,int_stack+2022,int_stack+1992, 0.0, zero_stack, 1.0, int_stack+2682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3534,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1845,int_stack+2067,int_stack+2022, 0.0, zero_stack, 1.0, int_stack+2742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3594,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1980,int_stack+2160,int_stack+2130, 0.0, zero_stack, 1.0, int_stack+2895, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3534, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10137,int_stack+2205,int_stack+2160, 0.0, zero_stack, 1.0, int_stack+2955, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3594, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2070,int_stack+2298,int_stack+2268, 0.0, zero_stack, 1.0, int_stack+3108, 0.0, zero_stack, 1.0, int_stack+3534, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2160,int_stack+2343,int_stack+2298, 0.0, zero_stack, 1.0, int_stack+3168, 0.0, zero_stack, 1.0, int_stack+3594, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2295,int_stack+2436,int_stack+2406, 0.0, zero_stack, 1.0, int_stack+3321, 1.0, int_stack+3534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10272,int_stack+2481,int_stack+2436, 0.0, zero_stack, 1.0, int_stack+3381, 1.0, int_stack+3594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2385,int_stack+2574,int_stack+2544, 0.0, zero_stack, 2.0, int_stack+3534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10407,int_stack+2619,int_stack+2574, 0.0, zero_stack, 2.0, int_stack+3594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2475,int_stack+2787,int_stack+2712, 1.0, int_stack+2682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3777,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2565,int_stack+2832,int_stack+2787, 1.0, int_stack+2742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3837,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+3000,int_stack+2925, 1.0, int_stack+2895, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3777, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2790,int_stack+3045,int_stack+3000, 1.0, int_stack+2955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3837, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2925,int_stack+3213,int_stack+3138, 1.0, int_stack+3108, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3777, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3015,int_stack+3258,int_stack+3213, 1.0, int_stack+3168, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3837, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3150,int_stack+3426,int_stack+3351, 1.0, int_stack+3321, 0.0, zero_stack, 1.0, int_stack+3777, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3240,int_stack+3471,int_stack+3426, 1.0, int_stack+3381, 0.0, zero_stack, 1.0, int_stack+3837, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3375,int_stack+3639,int_stack+3564, 1.0, int_stack+3534, 1.0, int_stack+3777, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10542,int_stack+3684,int_stack+3639, 1.0, int_stack+3594, 1.0, int_stack+3837, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3465,int_stack+3882,int_stack+3807, 2.0, int_stack+3777, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3555,int_stack+3927,int_stack+3882, 2.0, int_stack+3837, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3690,int_stack+4020,int_stack+3990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4680,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3780,int_stack+4065,int_stack+4020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4740,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3915,int_stack+4158,int_stack+4128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4680, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4005,int_stack+4203,int_stack+4158, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4740, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4140,int_stack+4296,int_stack+4266, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4680, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10677,int_stack+4341,int_stack+4296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4740, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4230,int_stack+4434,int_stack+4404, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10812,int_stack+4479,int_stack+4434, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4320,int_stack+4572,int_stack+4542, 0.0, zero_stack, 1.0, int_stack+4680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4410,int_stack+4617,int_stack+4572, 0.0, zero_stack, 1.0, int_stack+4740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4545,int_stack+4785,int_stack+4710, 1.0, int_stack+4680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10947,int_stack+4830,int_stack+4785, 1.0, int_stack+4740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4635,int_stack+4923,int_stack+4893,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4725,int_stack+4968,int_stack+4923,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4860,int_stack+5061,int_stack+5031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5721,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11082,int_stack+5106,int_stack+5061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5781,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4950,int_stack+5199,int_stack+5169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5721, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5040,int_stack+5244,int_stack+5199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5781, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5175,int_stack+5337,int_stack+5307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5721, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11217,int_stack+5382,int_stack+5337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5781, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5265,int_stack+5475,int_stack+5445, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5721, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11352,int_stack+5520,int_stack+5475, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5355,int_stack+5613,int_stack+5583, 0.0, zero_stack, 1.0, int_stack+5721, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5445,int_stack+5658,int_stack+5613, 0.0, zero_stack, 1.0, int_stack+5781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5580,int_stack+5826,int_stack+5751, 1.0, int_stack+5721, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11487,int_stack+5871,int_stack+5826, 1.0, int_stack+5781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5670,int_stack+5964,int_stack+5934,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5760,int_stack+6009,int_stack+5964,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5895,int_stack+6102,int_stack+6072,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11622,int_stack+6147,int_stack+6102,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5985,int_stack+6240,int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6900,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6075,int_stack+6285,int_stack+6240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6960,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6210,int_stack+6378,int_stack+6348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6900, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11757,int_stack+6423,int_stack+6378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6960, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6300,int_stack+6516,int_stack+6486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6900, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11892,int_stack+6561,int_stack+6516, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6960, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6390,int_stack+6654,int_stack+6624, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6480,int_stack+6699,int_stack+6654, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6615,int_stack+6792,int_stack+6762, 0.0, zero_stack, 1.0, int_stack+6900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12027,int_stack+6837,int_stack+6792, 0.0, zero_stack, 1.0, int_stack+6960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6705,int_stack+7005,int_stack+6930, 1.0, int_stack+6900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6795,int_stack+7050,int_stack+7005, 1.0, int_stack+6960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6930,int_stack+7143,int_stack+7113,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12162,int_stack+7188,int_stack+7143,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7020,int_stack+7281,int_stack+7251,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7110,int_stack+7326,int_stack+7281,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7245,int_stack+7419,int_stack+7389,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12297,int_stack+7464,int_stack+7419,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7335,int_stack+7707,int_stack+7617, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7527,3);
     Libderiv->ABCD[11] = int_stack + 7335;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12432,int_stack+7932,int_stack+7842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7527, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 12432;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12612,int_stack+8067,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7527, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 12612;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7932,int_stack+8202,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7527, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 7932;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8112,int_stack+8427,int_stack+8337, 0.0, zero_stack, 1.0, int_stack+7527, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 8112;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12792,int_stack+8562,int_stack+180, 1.0, int_stack+7527, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 12792;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8427,int_stack+8697,int_stack+270,3);
     Libderiv->ABCD[2] = int_stack + 8427;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8607,int_stack+8832,int_stack+360,3);
     Libderiv->ABCD[1] = int_stack + 8607;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8787,int_stack+8967,int_stack+450,3);
     Libderiv->ABCD[0] = int_stack + 8787;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12972,int_stack+9192,int_stack+9102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7617,3);
     Libderiv->ABCD[155] = int_stack + 12972;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8967,int_stack+630,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7617, 1.0, int_stack+7842,3);
     Libderiv->ABCD[143] = int_stack + 8967;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9147,int_stack+9327,int_stack+765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7842, 0.0, zero_stack,3);
     Libderiv->ABCD[142] = int_stack + 9147;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+540,int_stack+9462,int_stack+855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7617, 0.0, zero_stack, 1.0, int_stack+0,3);
     Libderiv->ABCD[131] = int_stack + 540;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9327,int_stack+1035,int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7842, 1.0, int_stack+0, 0.0, zero_stack,3);
     Libderiv->ABCD[130] = int_stack + 9327;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+720,int_stack+9597,int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[129] = int_stack + 720;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9507,int_stack+9732,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7617, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90,3);
     Libderiv->ABCD[119] = int_stack + 9507;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9687,int_stack+1440,int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7842, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack,3);
     Libderiv->ABCD[118] = int_stack + 9687;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+900,int_stack+9867,int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[117] = int_stack + 900;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1080,int_stack+10002,int_stack+1665, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[116] = int_stack + 1080;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9867,int_stack+1845,int_stack+1755, 0.0, zero_stack, 1.0, int_stack+7617, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8337,3);
     Libderiv->ABCD[107] = int_stack + 9867;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1260,int_stack+10137,int_stack+1980, 0.0, zero_stack, 1.0, int_stack+7842, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8337, 0.0, zero_stack,3);
     Libderiv->ABCD[106] = int_stack + 1260;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10047,int_stack+2160,int_stack+2070, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+8337, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[105] = int_stack + 10047;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1440,int_stack+10272,int_stack+2295, 0.0, zero_stack, 1.0, int_stack+90, 1.0, int_stack+8337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[104] = int_stack + 1440;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10227,int_stack+10407,int_stack+2385, 0.0, zero_stack, 2.0, int_stack+8337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[103] = int_stack + 10227;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1620,int_stack+2565,int_stack+2475, 1.0, int_stack+7617, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180,3);
     Libderiv->ABCD[95] = int_stack + 1620;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1800,int_stack+2790,int_stack+2700, 1.0, int_stack+7842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack,3);
     Libderiv->ABCD[94] = int_stack + 1800;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1980,int_stack+3015,int_stack+2925, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[93] = int_stack + 1980;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2160,int_stack+3240,int_stack+3150, 1.0, int_stack+90, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[92] = int_stack + 2160;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+10542,int_stack+3375, 1.0, int_stack+8337, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[91] = int_stack + 0;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10407,int_stack+3555,int_stack+3465, 2.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[90] = int_stack + 10407;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2340,int_stack+3780,int_stack+3690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270,3);
     Libderiv->ABCD[47] = int_stack + 2340;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2520,int_stack+4005,int_stack+3915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack,3);
     Libderiv->ABCD[46] = int_stack + 2520;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2700,int_stack+10677,int_stack+4140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[45] = int_stack + 2700;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10587,int_stack+10812,int_stack+4230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[44] = int_stack + 10587;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10767,int_stack+4410,int_stack+4320, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[43] = int_stack + 10767;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2880,int_stack+10947,int_stack+4545, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[42] = int_stack + 2880;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+180,int_stack+4725,int_stack+4635,3);
     Libderiv->ABCD[38] = int_stack + 180;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3060,int_stack+11082,int_stack+4860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360,3);
     Libderiv->ABCD[35] = int_stack + 3060;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10947,int_stack+5040,int_stack+4950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack,3);
     Libderiv->ABCD[34] = int_stack + 10947;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3240,int_stack+11217,int_stack+5175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[33] = int_stack + 3240;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11127,int_stack+11352,int_stack+5265, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[32] = int_stack + 11127;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11307,int_stack+5445,int_stack+5355, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[31] = int_stack + 11307;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3420,int_stack+11487,int_stack+5580, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[30] = int_stack + 3420;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3600,int_stack+5760,int_stack+5670,3);
     Libderiv->ABCD[26] = int_stack + 3600;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3780,int_stack+11622,int_stack+5895,3);
     Libderiv->ABCD[25] = int_stack + 3780;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11487,int_stack+6075,int_stack+5985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,3);
     Libderiv->ABCD[23] = int_stack + 11487;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3960,int_stack+11757,int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,3);
     Libderiv->ABCD[22] = int_stack + 3960;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11667,int_stack+11892,int_stack+6300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[21] = int_stack + 11667;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11847,int_stack+6480,int_stack+6390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[20] = int_stack + 11847;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4140,int_stack+12027,int_stack+6615, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[19] = int_stack + 4140;
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4320,int_stack+6795,int_stack+6705, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[18] = int_stack + 4320;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+360,int_stack+12162,int_stack+6930,3);
     Libderiv->ABCD[14] = int_stack + 360;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12027,int_stack+7110,int_stack+7020,3);
     Libderiv->ABCD[13] = int_stack + 12027;
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4500,int_stack+12297,int_stack+7245,3);
     Libderiv->ABCD[12] = int_stack + 4500;

}
