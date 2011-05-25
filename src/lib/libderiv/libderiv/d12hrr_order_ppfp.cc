#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ppfp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|fp) integrals */

void d12hrr_order_ppfp(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->dvrr_classes[1][4] = int_stack + 780;
 Libderiv->deriv_classes[2][4][0] = int_stack + 825;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 915;
 Libderiv->deriv2_classes[1][4][143] = int_stack + 945;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 990;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 1050;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 1140;
 Libderiv->deriv2_classes[1][4][131] = int_stack + 1170;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1215;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 1275;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 1365;
 Libderiv->deriv2_classes[1][4][130] = int_stack + 1395;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 1440;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 1500;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 1590;
 Libderiv->deriv2_classes[1][4][119] = int_stack + 1620;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 1665;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 1725;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 1815;
 Libderiv->deriv2_classes[1][4][118] = int_stack + 1845;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 1890;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 1950;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 2040;
 Libderiv->deriv2_classes[1][4][117] = int_stack + 2070;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 2115;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 2175;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 2265;
 Libderiv->deriv2_classes[1][4][107] = int_stack + 2295;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 2340;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 2400;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 2490;
 Libderiv->deriv2_classes[1][4][106] = int_stack + 2520;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 2565;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 2625;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 2715;
 Libderiv->deriv2_classes[1][4][105] = int_stack + 2745;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 2790;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 2850;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 2940;
 Libderiv->deriv2_classes[1][4][104] = int_stack + 2970;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 3015;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 3075;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 3165;
 Libderiv->deriv2_classes[1][4][95] = int_stack + 3195;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 3240;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 3300;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 3390;
 Libderiv->deriv2_classes[1][4][94] = int_stack + 3420;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 3465;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 3525;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 3615;
 Libderiv->deriv2_classes[1][4][93] = int_stack + 3645;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 3690;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 3750;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 3840;
 Libderiv->deriv2_classes[1][4][92] = int_stack + 3870;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 3915;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 3975;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 4065;
 Libderiv->deriv2_classes[1][4][91] = int_stack + 4095;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 4140;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 4200;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 4290;
 Libderiv->deriv2_classes[1][4][83] = int_stack + 4320;
 Libderiv->deriv_classes[2][3][11] = int_stack + 4365;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 4425;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 4485;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 4575;
 Libderiv->deriv2_classes[1][4][82] = int_stack + 4605;
 Libderiv->deriv_classes[2][3][10] = int_stack + 4650;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 4710;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 4770;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 4860;
 Libderiv->deriv2_classes[1][4][81] = int_stack + 4890;
 Libderiv->deriv_classes[2][3][9] = int_stack + 4935;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 4995;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 5055;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 5145;
 Libderiv->deriv2_classes[1][4][80] = int_stack + 5175;
 Libderiv->deriv_classes[2][3][8] = int_stack + 5220;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 5280;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 5340;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 5430;
 Libderiv->deriv2_classes[1][4][79] = int_stack + 5460;
 Libderiv->deriv_classes[2][3][7] = int_stack + 5505;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 5565;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 5625;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 5715;
 Libderiv->deriv2_classes[1][4][78] = int_stack + 5745;
 Libderiv->deriv_classes[2][3][6] = int_stack + 5790;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 5850;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 5910;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 6000;
 Libderiv->deriv2_classes[1][4][35] = int_stack + 6030;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 6075;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 6135;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 6225;
 Libderiv->deriv2_classes[1][4][34] = int_stack + 6255;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 6300;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 6360;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 6450;
 Libderiv->deriv2_classes[1][4][33] = int_stack + 6480;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 6525;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 6585;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 6675;
 Libderiv->deriv2_classes[1][4][32] = int_stack + 6705;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 6750;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 6810;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 6900;
 Libderiv->deriv2_classes[1][4][31] = int_stack + 6930;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 6975;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 7035;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 7125;
 Libderiv->deriv2_classes[1][4][30] = int_stack + 7155;
 Libderiv->deriv_classes[2][3][2] = int_stack + 7200;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 7260;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 7320;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 7410;
 Libderiv->deriv2_classes[1][4][26] = int_stack + 7440;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 7485;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 7545;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 7635;
 Libderiv->deriv2_classes[1][4][23] = int_stack + 7665;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 7710;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 7770;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 7860;
 Libderiv->deriv2_classes[1][4][22] = int_stack + 7890;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 7935;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 7995;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 8085;
 Libderiv->deriv2_classes[1][4][21] = int_stack + 8115;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 8160;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 8220;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 8310;
 Libderiv->deriv2_classes[1][4][20] = int_stack + 8340;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 8385;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 8445;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 8535;
 Libderiv->deriv2_classes[1][4][19] = int_stack + 8565;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 8610;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 8670;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 8760;
 Libderiv->deriv2_classes[1][4][18] = int_stack + 8790;
 Libderiv->deriv_classes[2][3][1] = int_stack + 8835;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 8895;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 8955;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 9045;
 Libderiv->deriv2_classes[1][4][14] = int_stack + 9075;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 9120;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 9180;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 9270;
 Libderiv->deriv2_classes[1][4][13] = int_stack + 9300;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 9345;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 9405;
 Libderiv->deriv_classes[1][3][11] = int_stack + 9495;
 Libderiv->deriv_classes[1][4][11] = int_stack + 9525;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 9570;
 Libderiv->deriv2_classes[1][4][11] = int_stack + 9600;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 9645;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 9705;
 Libderiv->deriv_classes[1][3][10] = int_stack + 9795;
 Libderiv->deriv_classes[1][4][10] = int_stack + 9825;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 9870;
 Libderiv->deriv2_classes[1][4][10] = int_stack + 9900;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 9945;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 10005;
 Libderiv->deriv_classes[1][3][9] = int_stack + 10095;
 Libderiv->deriv_classes[1][4][9] = int_stack + 10125;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 10170;
 Libderiv->deriv2_classes[1][4][9] = int_stack + 10200;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 10245;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 10305;
 Libderiv->deriv_classes[1][3][8] = int_stack + 10395;
 Libderiv->deriv_classes[1][4][8] = int_stack + 10425;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 10470;
 Libderiv->deriv2_classes[1][4][8] = int_stack + 10500;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 10545;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 10605;
 Libderiv->deriv_classes[1][3][7] = int_stack + 10695;
 Libderiv->deriv_classes[1][4][7] = int_stack + 10725;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 10770;
 Libderiv->deriv2_classes[1][4][7] = int_stack + 10800;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 10845;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 10905;
 Libderiv->dvrr_classes[1][3] = int_stack + 10995;
 Libderiv->deriv_classes[1][3][6] = int_stack + 11025;
 Libderiv->deriv_classes[1][4][6] = int_stack + 11055;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 11100;
 Libderiv->deriv2_classes[1][4][6] = int_stack + 11130;
 Libderiv->deriv_classes[2][3][0] = int_stack + 11175;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 11235;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 11295;
 Libderiv->deriv_classes[1][3][2] = int_stack + 11385;
 Libderiv->deriv_classes[1][4][2] = int_stack + 11415;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 11460;
 Libderiv->deriv2_classes[1][4][2] = int_stack + 11490;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 11535;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 11595;
 Libderiv->deriv_classes[1][3][1] = int_stack + 11685;
 Libderiv->deriv_classes[1][4][1] = int_stack + 11715;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 11760;
 Libderiv->deriv2_classes[1][4][1] = int_stack + 11790;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 11835;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 11895;
 Libderiv->deriv_classes[1][3][0] = int_stack + 11985;
 Libderiv->deriv_classes[1][4][0] = int_stack + 12015;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 12060;
 Libderiv->deriv2_classes[1][4][0] = int_stack + 12090;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 12135;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 12195;
 memset(int_stack,0,98280);

 Libderiv->dvrr_stack = int_stack + 16155;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ppfp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12285,int_stack+9525,int_stack+9495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10995,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12375,int_stack+0,int_stack+4365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+9825,int_stack+9795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10995, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12555,int_stack+90,int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+10125,int_stack+10095, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10995, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12735,int_stack+180,int_stack+4935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+10425,int_stack+10395, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12915,int_stack+270,int_stack+5220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+270,int_stack+10725,int_stack+10695, 0.0, zero_stack, 1.0, int_stack+10995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13095,int_stack+360,int_stack+5505, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+11055,int_stack+11025, 1.0, int_stack+10995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13275,int_stack+510,int_stack+5790, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+780,int_stack+10995,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13455,int_stack+11415,int_stack+11385,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13545,int_stack+600,int_stack+7200,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+11715,int_stack+11685,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13725,int_stack+690,int_stack+8835,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+12015,int_stack+11985,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13905,int_stack+825,int_stack+11175,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+945,int_stack+915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9495,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+810,int_stack+1050,int_stack+990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4365,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+990,int_stack+1170,int_stack+1140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9495, 1.0, int_stack+9795,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14085,int_stack+1275,int_stack+1215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4365, 1.0, int_stack+4650,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+1395,int_stack+1365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9795, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1170,int_stack+1500,int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4650, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1350,int_stack+1620,int_stack+1590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9495, 0.0, zero_stack, 1.0, int_stack+10095,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1440,int_stack+1725,int_stack+1665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4365, 0.0, zero_stack, 1.0, int_stack+4935,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1620,int_stack+1845,int_stack+1815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9795, 1.0, int_stack+10095, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1710,int_stack+1950,int_stack+1890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4650, 1.0, int_stack+4935, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1890,int_stack+2070,int_stack+2040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10095, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14265,int_stack+2175,int_stack+2115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4935, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1980,int_stack+2295,int_stack+2265, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9495, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10395,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2070,int_stack+2400,int_stack+2340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4365, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5220,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2250,int_stack+2520,int_stack+2490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9795, 0.0, zero_stack, 1.0, int_stack+10395, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2340,int_stack+2625,int_stack+2565, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4650, 0.0, zero_stack, 1.0, int_stack+5220, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2520,int_stack+2745,int_stack+2715, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10095, 1.0, int_stack+10395, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2610,int_stack+2850,int_stack+2790, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4935, 1.0, int_stack+5220, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2790,int_stack+2970,int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14445,int_stack+3075,int_stack+3015, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2880,int_stack+3195,int_stack+3165, 0.0, zero_stack, 1.0, int_stack+9495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10695,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2970,int_stack+3300,int_stack+3240, 0.0, zero_stack, 1.0, int_stack+4365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5505,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3150,int_stack+3420,int_stack+3390, 0.0, zero_stack, 1.0, int_stack+9795, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10695, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3240,int_stack+3525,int_stack+3465, 0.0, zero_stack, 1.0, int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5505, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3420,int_stack+3645,int_stack+3615, 0.0, zero_stack, 1.0, int_stack+10095, 0.0, zero_stack, 1.0, int_stack+10695, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3510,int_stack+3750,int_stack+3690, 0.0, zero_stack, 1.0, int_stack+4935, 0.0, zero_stack, 1.0, int_stack+5505, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3690,int_stack+3870,int_stack+3840, 0.0, zero_stack, 1.0, int_stack+10395, 1.0, int_stack+10695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14625,int_stack+3975,int_stack+3915, 0.0, zero_stack, 1.0, int_stack+5220, 1.0, int_stack+5505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3780,int_stack+4095,int_stack+4065, 0.0, zero_stack, 2.0, int_stack+10695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3870,int_stack+4200,int_stack+4140, 0.0, zero_stack, 2.0, int_stack+5505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4050,int_stack+4320,int_stack+4290, 1.0, int_stack+9495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11025,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4140,int_stack+4485,int_stack+4425, 1.0, int_stack+4365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5790,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4320,int_stack+4605,int_stack+4575, 1.0, int_stack+9795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11025, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4410,int_stack+4770,int_stack+4710, 1.0, int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5790, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4590,int_stack+4890,int_stack+4860, 1.0, int_stack+10095, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11025, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4680,int_stack+5055,int_stack+4995, 1.0, int_stack+4935, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5790, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4860,int_stack+5175,int_stack+5145, 1.0, int_stack+10395, 0.0, zero_stack, 1.0, int_stack+11025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4950,int_stack+5340,int_stack+5280, 1.0, int_stack+5220, 0.0, zero_stack, 1.0, int_stack+5790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5130,int_stack+5460,int_stack+5430, 1.0, int_stack+10695, 1.0, int_stack+11025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5220,int_stack+5625,int_stack+5565, 1.0, int_stack+5505, 1.0, int_stack+5790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5400,int_stack+5745,int_stack+5715, 2.0, int_stack+11025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5490,int_stack+5910,int_stack+5850, 2.0, int_stack+5790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5670,int_stack+6030,int_stack+6000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11385,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5760,int_stack+6135,int_stack+6075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7200,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5940,int_stack+6255,int_stack+6225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11385, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6030,int_stack+6360,int_stack+6300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7200, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6210,int_stack+6480,int_stack+6450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11385, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6300,int_stack+6585,int_stack+6525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7200, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6480,int_stack+6705,int_stack+6675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6570,int_stack+6810,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6750,int_stack+6930,int_stack+6900, 0.0, zero_stack, 1.0, int_stack+11385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14805,int_stack+7035,int_stack+6975, 0.0, zero_stack, 1.0, int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6840,int_stack+7155,int_stack+7125, 1.0, int_stack+11385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6930,int_stack+7320,int_stack+7260, 1.0, int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7110,int_stack+7440,int_stack+7410,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+7545,int_stack+7485,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7380,int_stack+7665,int_stack+7635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11685,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7470,int_stack+7770,int_stack+7710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8835,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7650,int_stack+7890,int_stack+7860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11685, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7740,int_stack+7995,int_stack+7935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8835, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7920,int_stack+8115,int_stack+8085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11685, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14985,int_stack+8220,int_stack+8160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8835, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8010,int_stack+8340,int_stack+8310, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11685, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8100,int_stack+8445,int_stack+8385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8280,int_stack+8565,int_stack+8535, 0.0, zero_stack, 1.0, int_stack+11685, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8370,int_stack+8670,int_stack+8610, 0.0, zero_stack, 1.0, int_stack+8835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8550,int_stack+8790,int_stack+8760, 1.0, int_stack+11685, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8640,int_stack+8955,int_stack+8895, 1.0, int_stack+8835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8820,int_stack+9075,int_stack+9045,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8910,int_stack+9180,int_stack+9120,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+9090,int_stack+9300,int_stack+9270,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15165,int_stack+9405,int_stack+9345,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9180,int_stack+9600,int_stack+9570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11985,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9270,int_stack+9705,int_stack+9645, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11175,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9450,int_stack+9900,int_stack+9870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11985, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9540,int_stack+10005,int_stack+9945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11175, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9720,int_stack+10200,int_stack+10170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11985, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9810,int_stack+10305,int_stack+10245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11175, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9990,int_stack+10500,int_stack+10470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10080,int_stack+10605,int_stack+10545, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10260,int_stack+10800,int_stack+10770, 0.0, zero_stack, 1.0, int_stack+11985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10350,int_stack+10905,int_stack+10845, 0.0, zero_stack, 1.0, int_stack+11175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10530,int_stack+11130,int_stack+11100, 1.0, int_stack+11985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10620,int_stack+11295,int_stack+11235, 1.0, int_stack+11175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10800,int_stack+11490,int_stack+11460,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10890,int_stack+11595,int_stack+11535,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11070,int_stack+11790,int_stack+11760,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11160,int_stack+11895,int_stack+11835,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11340,int_stack+12090,int_stack+12060,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11430,int_stack+12195,int_stack+12135,6);
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+11610,int_stack+12375,int_stack+12285,30);
     Libderiv->ABCD[11] = int_stack + 11610;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+11880,int_stack+12555,int_stack+0,30);
     Libderiv->ABCD[10] = int_stack + 11880;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+12375,int_stack+12735,int_stack+90,30);
     Libderiv->ABCD[9] = int_stack + 12375;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+12645,int_stack+12915,int_stack+180,30);
     Libderiv->ABCD[8] = int_stack + 12645;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+15345,int_stack+13095,int_stack+270,30);
     Libderiv->ABCD[7] = int_stack + 15345;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+12915,int_stack+13275,int_stack+360,30);
     Libderiv->ABCD[6] = int_stack + 12915;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+13185,int_stack+13545,int_stack+13455, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[2] = int_stack + 13185;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+15615,int_stack+13725,int_stack+540, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[1] = int_stack + 15615;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+13545,int_stack+13905,int_stack+630, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[0] = int_stack + 13545;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+13815,int_stack+810,int_stack+720,30);
     Libderiv->ABCD[155] = int_stack + 13815;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+720,int_stack+14085,int_stack+990,30);
     Libderiv->ABCD[143] = int_stack + 720;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+15885,int_stack+1170,int_stack+1080,30);
     Libderiv->ABCD[142] = int_stack + 15885;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+990,int_stack+1440,int_stack+1350,30);
     Libderiv->ABCD[131] = int_stack + 990;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1260,int_stack+1710,int_stack+1620,30);
     Libderiv->ABCD[130] = int_stack + 1260;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1530,int_stack+14265,int_stack+1890,30);
     Libderiv->ABCD[129] = int_stack + 1530;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+14085,int_stack+2070,int_stack+1980,30);
     Libderiv->ABCD[119] = int_stack + 14085;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1800,int_stack+2340,int_stack+2250,30);
     Libderiv->ABCD[118] = int_stack + 1800;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2070,int_stack+2610,int_stack+2520,30);
     Libderiv->ABCD[117] = int_stack + 2070;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2340,int_stack+14445,int_stack+2790,30);
     Libderiv->ABCD[116] = int_stack + 2340;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2610,int_stack+2970,int_stack+2880,30);
     Libderiv->ABCD[107] = int_stack + 2610;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2880,int_stack+3240,int_stack+3150,30);
     Libderiv->ABCD[106] = int_stack + 2880;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3150,int_stack+3510,int_stack+3420,30);
     Libderiv->ABCD[105] = int_stack + 3150;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3420,int_stack+14625,int_stack+3690,30);
     Libderiv->ABCD[104] = int_stack + 3420;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+14355,int_stack+3870,int_stack+3780,30);
     Libderiv->ABCD[103] = int_stack + 14355;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3690,int_stack+4140,int_stack+4050,30);
     Libderiv->ABCD[95] = int_stack + 3690;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3960,int_stack+4410,int_stack+4320,30);
     Libderiv->ABCD[94] = int_stack + 3960;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4230,int_stack+4680,int_stack+4590,30);
     Libderiv->ABCD[93] = int_stack + 4230;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4500,int_stack+4950,int_stack+4860,30);
     Libderiv->ABCD[92] = int_stack + 4500;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4770,int_stack+5220,int_stack+5130,30);
     Libderiv->ABCD[91] = int_stack + 4770;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+5040,int_stack+5490,int_stack+5400,30);
     Libderiv->ABCD[90] = int_stack + 5040;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5310,int_stack+5760,int_stack+5670, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[47] = int_stack + 5310;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5580,int_stack+6030,int_stack+5940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[46] = int_stack + 5580;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5850,int_stack+6300,int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[45] = int_stack + 5850;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6120,int_stack+6570,int_stack+6480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[44] = int_stack + 6120;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6390,int_stack+14805,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[43] = int_stack + 6390;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+14625,int_stack+6930,int_stack+6840, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[42] = int_stack + 14625;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6660,int_stack+7200,int_stack+7110, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[38] = int_stack + 6660;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6930,int_stack+7470,int_stack+7380, 0.0, zero_stack, 1.0, int_stack+12285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[35] = int_stack + 6930;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7200,int_stack+7740,int_stack+7650, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[34] = int_stack + 7200;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7470,int_stack+14985,int_stack+7920, 0.0, zero_stack, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[33] = int_stack + 7470;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7740,int_stack+8100,int_stack+8010, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[32] = int_stack + 7740;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8010,int_stack+8370,int_stack+8280, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[31] = int_stack + 8010;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8280,int_stack+8640,int_stack+8550, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[30] = int_stack + 8280;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8550,int_stack+8910,int_stack+8820, 0.0, zero_stack, 1.0, int_stack+13455, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[26] = int_stack + 8550;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8820,int_stack+15165,int_stack+9090, 0.0, zero_stack, 2.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[25] = int_stack + 8820;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+14895,int_stack+9270,int_stack+9180, 1.0, int_stack+12285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[23] = int_stack + 14895;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9090,int_stack+9540,int_stack+9450, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[22] = int_stack + 9090;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9360,int_stack+9810,int_stack+9720, 1.0, int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[21] = int_stack + 9360;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9630,int_stack+10080,int_stack+9990, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[20] = int_stack + 9630;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+10350,int_stack+10260, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[19] = int_stack + 0;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9900,int_stack+10620,int_stack+10530, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[18] = int_stack + 9900;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+270,int_stack+10890,int_stack+10800, 1.0, int_stack+13455, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[14] = int_stack + 270;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+10170,int_stack+11160,int_stack+11070, 1.0, int_stack+540, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[13] = int_stack + 10170;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+10440,int_stack+11430,int_stack+11340, 2.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[12] = int_stack + 10440;

}
