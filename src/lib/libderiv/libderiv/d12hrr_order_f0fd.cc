#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_f0fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|fd) integrals */

void d12hrr_order_f0fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][5][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][10] = int_stack + 210;
 Libderiv->deriv_classes[3][5][9] = int_stack + 420;
 Libderiv->deriv_classes[3][5][8] = int_stack + 630;
 Libderiv->deriv_classes[3][5][7] = int_stack + 840;
 Libderiv->dvrr_classes[3][4] = int_stack + 1050;
 Libderiv->deriv_classes[3][5][6] = int_stack + 1200;
 Libderiv->deriv_classes[3][5][2] = int_stack + 1410;
 Libderiv->deriv_classes[3][5][1] = int_stack + 1620;
 Libderiv->deriv_classes[3][5][0] = int_stack + 1830;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 2040;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 2140;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 2290;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 2500;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 2600;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 2750;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 2960;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 3060;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 3210;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 3420;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 3520;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 3670;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 3880;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 3980;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 4130;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 4340;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 4440;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 4590;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 4800;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 4900;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 5050;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 5260;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 5360;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 5510;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 5720;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 5820;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 5970;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 6180;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 6280;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 6430;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 6640;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 6740;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 6890;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 7100;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 7200;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 7350;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 7560;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 7660;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 7810;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 8020;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 8120;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 8270;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 8480;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 8580;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 8730;
 Libderiv->deriv_classes[3][3][11] = int_stack + 8940;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 9040;
 Libderiv->deriv_classes[3][4][11] = int_stack + 9140;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 9290;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 9440;
 Libderiv->deriv_classes[3][3][10] = int_stack + 9650;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 9750;
 Libderiv->deriv_classes[3][4][10] = int_stack + 9850;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 10000;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 10150;
 Libderiv->deriv_classes[3][3][9] = int_stack + 10360;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 10460;
 Libderiv->deriv_classes[3][4][9] = int_stack + 10560;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 10710;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 10860;
 Libderiv->deriv_classes[3][3][8] = int_stack + 11070;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 11170;
 Libderiv->deriv_classes[3][4][8] = int_stack + 11270;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 11420;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 11570;
 Libderiv->deriv_classes[3][3][7] = int_stack + 11780;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 11880;
 Libderiv->deriv_classes[3][4][7] = int_stack + 11980;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 12130;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 12280;
 Libderiv->dvrr_classes[3][3] = int_stack + 12490;
 Libderiv->deriv_classes[3][3][6] = int_stack + 12590;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 12690;
 Libderiv->deriv_classes[3][4][6] = int_stack + 12790;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 12940;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 13090;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 13300;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 13400;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 13550;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 13760;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 13860;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 14010;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 14220;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 14320;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 14470;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 14680;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 14780;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 14930;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 15140;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 15240;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 15390;
 Libderiv->deriv_classes[3][3][2] = int_stack + 15600;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 15700;
 Libderiv->deriv_classes[3][4][2] = int_stack + 15800;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 15950;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 16100;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 16310;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 16410;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 16560;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 16770;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 16870;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 17020;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 17230;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 17330;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 17480;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 17690;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 17790;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 17940;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 18150;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 18250;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 18400;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 18610;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 18710;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 18860;
 Libderiv->deriv_classes[3][3][1] = int_stack + 19070;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 19170;
 Libderiv->deriv_classes[3][4][1] = int_stack + 19270;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 19420;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 19570;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 19780;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 19880;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 20030;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 20240;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 20340;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 20490;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 20700;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 20800;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 20950;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 21160;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 21260;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 21410;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 21620;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 21720;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 21870;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 22080;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 22180;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 22330;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 22540;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 22640;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 22790;
 Libderiv->deriv_classes[3][3][0] = int_stack + 23000;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 23100;
 Libderiv->deriv_classes[3][4][0] = int_stack + 23200;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 23350;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 23500;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 23710;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 23810;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 23960;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 24170;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 24270;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 24420;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 24630;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 24730;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 24880;
 memset(int_stack,0,200720);

 Libderiv->dvrr_stack = int_stack + 43840;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_f0fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+25090,int_stack+1050,int_stack+12490,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25390,int_stack+9140,int_stack+8940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12490,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25690,int_stack+0,int_stack+9140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26140,int_stack+9850,int_stack+9650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12490, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26440,int_stack+210,int_stack+9850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+10560,int_stack+10360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12490, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26890,int_stack+420,int_stack+10560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+11270,int_stack+11070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+27340,int_stack+630,int_stack+11270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27790,int_stack+11980,int_stack+11780, 0.0, zero_stack, 1.0, int_stack+12490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+28090,int_stack+840,int_stack+11980, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+12790,int_stack+12590, 1.0, int_stack+12490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+28540,int_stack+1200,int_stack+12790, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+15800,int_stack+15600,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28990,int_stack+1410,int_stack+15800,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+19270,int_stack+19070,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+29440,int_stack+1620,int_stack+19270,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1500,int_stack+23200,int_stack+23000,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+29890,int_stack+1830,int_stack+23200,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30340,int_stack+2140,int_stack+2040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8940,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+30640,int_stack+2290,int_stack+2140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9140,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2600,int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8940, 1.0, int_stack+9650,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2100,int_stack+2750,int_stack+2600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9140, 1.0, int_stack+9850,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+3060,int_stack+2960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9650, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31090,int_stack+3210,int_stack+3060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9850, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+3520,int_stack+3420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8940, 0.0, zero_stack, 1.0, int_stack+10360,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31540,int_stack+3670,int_stack+3520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9140, 0.0, zero_stack, 1.0, int_stack+10560,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3150,int_stack+3980,int_stack+3880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9650, 1.0, int_stack+10360, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3450,int_stack+4130,int_stack+3980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9850, 1.0, int_stack+10560, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3900,int_stack+4440,int_stack+4340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10360, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31990,int_stack+4590,int_stack+4440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10560, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4200,int_stack+4900,int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11070,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32440,int_stack+5050,int_stack+4900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11270,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+5360,int_stack+5260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9650, 0.0, zero_stack, 1.0, int_stack+11070, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4800,int_stack+5510,int_stack+5360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9850, 0.0, zero_stack, 1.0, int_stack+11270, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5250,int_stack+5820,int_stack+5720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10360, 1.0, int_stack+11070, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32890,int_stack+5970,int_stack+5820, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10560, 1.0, int_stack+11270, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5550,int_stack+6280,int_stack+6180, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33340,int_stack+6430,int_stack+6280, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5850,int_stack+6740,int_stack+6640, 0.0, zero_stack, 1.0, int_stack+8940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11780,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6150,int_stack+6890,int_stack+6740, 0.0, zero_stack, 1.0, int_stack+9140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11980,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6600,int_stack+7200,int_stack+7100, 0.0, zero_stack, 1.0, int_stack+9650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11780, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33790,int_stack+7350,int_stack+7200, 0.0, zero_stack, 1.0, int_stack+9850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11980, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6900,int_stack+7660,int_stack+7560, 0.0, zero_stack, 1.0, int_stack+10360, 0.0, zero_stack, 1.0, int_stack+11780, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7200,int_stack+7810,int_stack+7660, 0.0, zero_stack, 1.0, int_stack+10560, 0.0, zero_stack, 1.0, int_stack+11980, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7650,int_stack+8120,int_stack+8020, 0.0, zero_stack, 1.0, int_stack+11070, 1.0, int_stack+11780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34240,int_stack+8270,int_stack+8120, 0.0, zero_stack, 1.0, int_stack+11270, 1.0, int_stack+11980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7950,int_stack+8580,int_stack+8480, 0.0, zero_stack, 2.0, int_stack+11780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34690,int_stack+8730,int_stack+8580, 0.0, zero_stack, 2.0, int_stack+11980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8250,int_stack+9290,int_stack+9040, 1.0, int_stack+8940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12590,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8550,int_stack+9440,int_stack+9290, 1.0, int_stack+9140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12790,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9000,int_stack+10000,int_stack+9750, 1.0, int_stack+9650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12590, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9300,int_stack+10150,int_stack+10000, 1.0, int_stack+9850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12790, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9750,int_stack+10710,int_stack+10460, 1.0, int_stack+10360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12590, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10050,int_stack+10860,int_stack+10710, 1.0, int_stack+10560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12790, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10500,int_stack+11420,int_stack+11170, 1.0, int_stack+11070, 0.0, zero_stack, 1.0, int_stack+12590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10800,int_stack+11570,int_stack+11420, 1.0, int_stack+11270, 0.0, zero_stack, 1.0, int_stack+12790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11250,int_stack+12130,int_stack+11880, 1.0, int_stack+11780, 1.0, int_stack+12590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+35140,int_stack+12280,int_stack+12130, 1.0, int_stack+11980, 1.0, int_stack+12790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11550,int_stack+12940,int_stack+12690, 2.0, int_stack+12590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11850,int_stack+13090,int_stack+12940, 2.0, int_stack+12790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12300,int_stack+13400,int_stack+13300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15600,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+13550,int_stack+13400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15800,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13050,int_stack+13860,int_stack+13760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15600, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13350,int_stack+14010,int_stack+13860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15800, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13800,int_stack+14320,int_stack+14220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15600, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+35590,int_stack+14470,int_stack+14320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15800, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14100,int_stack+14780,int_stack+14680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36040,int_stack+14930,int_stack+14780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14400,int_stack+15240,int_stack+15140, 0.0, zero_stack, 1.0, int_stack+15600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14700,int_stack+15390,int_stack+15240, 0.0, zero_stack, 1.0, int_stack+15800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15150,int_stack+15950,int_stack+15700, 1.0, int_stack+15600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36490,int_stack+16100,int_stack+15950, 1.0, int_stack+15800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15450,int_stack+16410,int_stack+16310,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15750,int_stack+16560,int_stack+16410,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16200,int_stack+16870,int_stack+16770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19070,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36940,int_stack+17020,int_stack+16870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19270,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16500,int_stack+17330,int_stack+17230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19070, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16800,int_stack+17480,int_stack+17330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19270, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17250,int_stack+17790,int_stack+17690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19070, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37390,int_stack+17940,int_stack+17790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19270, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17550,int_stack+18250,int_stack+18150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37840,int_stack+18400,int_stack+18250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17850,int_stack+18710,int_stack+18610, 0.0, zero_stack, 1.0, int_stack+19070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18150,int_stack+18860,int_stack+18710, 0.0, zero_stack, 1.0, int_stack+19270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18600,int_stack+19420,int_stack+19170, 1.0, int_stack+19070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+38290,int_stack+19570,int_stack+19420, 1.0, int_stack+19270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+18900,int_stack+19880,int_stack+19780,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19200,int_stack+20030,int_stack+19880,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19650,int_stack+20340,int_stack+20240,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+38740,int_stack+20490,int_stack+20340,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+19950,int_stack+20800,int_stack+20700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23000,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20250,int_stack+20950,int_stack+20800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23200,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20700,int_stack+21260,int_stack+21160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23000, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+39190,int_stack+21410,int_stack+21260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23200, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21000,int_stack+21720,int_stack+21620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23000, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+39640,int_stack+21870,int_stack+21720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23200, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21300,int_stack+22180,int_stack+22080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21600,int_stack+22330,int_stack+22180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22050,int_stack+22640,int_stack+22540, 0.0, zero_stack, 1.0, int_stack+23000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+40090,int_stack+22790,int_stack+22640, 0.0, zero_stack, 1.0, int_stack+23200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22350,int_stack+23350,int_stack+23100, 1.0, int_stack+23000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22650,int_stack+23500,int_stack+23350, 1.0, int_stack+23200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+23100,int_stack+23810,int_stack+23710,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+40540,int_stack+23960,int_stack+23810,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+23400,int_stack+24270,int_stack+24170,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23700,int_stack+24420,int_stack+24270,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24150,int_stack+24730,int_stack+24630,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+40990,int_stack+24880,int_stack+24730,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24450,int_stack+25690,int_stack+25390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25090,10);
     Libderiv->ABCD[11] = int_stack + 24450;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41440,int_stack+26440,int_stack+26140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25090, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 41440;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42040,int_stack+26890,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25090, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 42040;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26440,int_stack+27340,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 26440;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+27040,int_stack+28090,int_stack+27790, 0.0, zero_stack, 1.0, int_stack+25090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 27040;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42640,int_stack+28540,int_stack+600, 1.0, int_stack+25090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 42640;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+28090,int_stack+28990,int_stack+900,10);
     Libderiv->ABCD[2] = int_stack + 28090;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+28690,int_stack+29440,int_stack+1200,10);
     Libderiv->ABCD[1] = int_stack + 28690;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+29290,int_stack+29890,int_stack+1500,10);
     Libderiv->ABCD[0] = int_stack + 29290;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+43240,int_stack+30640,int_stack+30340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25390,10);
     Libderiv->ABCD[155] = int_stack + 43240;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29890,int_stack+2100,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25390, 1.0, int_stack+26140,10);
     Libderiv->ABCD[143] = int_stack + 29890;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30490,int_stack+31090,int_stack+2550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+26140, 0.0, zero_stack,10);
     Libderiv->ABCD[142] = int_stack + 30490;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1800,int_stack+31540,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25390, 0.0, zero_stack, 1.0, int_stack+0,10);
     Libderiv->ABCD[131] = int_stack + 1800;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31090,int_stack+3450,int_stack+3150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26140, 1.0, int_stack+0, 0.0, zero_stack,10);
     Libderiv->ABCD[130] = int_stack + 31090;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2400,int_stack+31990,int_stack+3900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[129] = int_stack + 2400;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31690,int_stack+32440,int_stack+4200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300,10);
     Libderiv->ABCD[119] = int_stack + 31690;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32290,int_stack+4800,int_stack+4500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26140, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack,10);
     Libderiv->ABCD[118] = int_stack + 32290;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3000,int_stack+32890,int_stack+5250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[117] = int_stack + 3000;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3600,int_stack+33340,int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[116] = int_stack + 3600;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32890,int_stack+6150,int_stack+5850, 0.0, zero_stack, 1.0, int_stack+25390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27790,10);
     Libderiv->ABCD[107] = int_stack + 32890;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4200,int_stack+33790,int_stack+6600, 0.0, zero_stack, 1.0, int_stack+26140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27790, 0.0, zero_stack,10);
     Libderiv->ABCD[106] = int_stack + 4200;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33490,int_stack+7200,int_stack+6900, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+27790, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[105] = int_stack + 33490;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4800,int_stack+34240,int_stack+7650, 0.0, zero_stack, 1.0, int_stack+300, 1.0, int_stack+27790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[104] = int_stack + 4800;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+34090,int_stack+34690,int_stack+7950, 0.0, zero_stack, 2.0, int_stack+27790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[103] = int_stack + 34090;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5400,int_stack+8550,int_stack+8250, 1.0, int_stack+25390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+600,10);
     Libderiv->ABCD[95] = int_stack + 5400;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6000,int_stack+9300,int_stack+9000, 1.0, int_stack+26140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack,10);
     Libderiv->ABCD[94] = int_stack + 6000;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6600,int_stack+10050,int_stack+9750, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[93] = int_stack + 6600;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7200,int_stack+10800,int_stack+10500, 1.0, int_stack+300, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[92] = int_stack + 7200;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+35140,int_stack+11250, 1.0, int_stack+27790, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[91] = int_stack + 0;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+34690,int_stack+11850,int_stack+11550, 2.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[90] = int_stack + 34690;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7800,int_stack+12600,int_stack+12300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900,10);
     Libderiv->ABCD[47] = int_stack + 7800;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8400,int_stack+13350,int_stack+13050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack,10);
     Libderiv->ABCD[46] = int_stack + 8400;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9000,int_stack+35590,int_stack+13800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[45] = int_stack + 9000;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35290,int_stack+36040,int_stack+14100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[44] = int_stack + 35290;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35890,int_stack+14700,int_stack+14400, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[43] = int_stack + 35890;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9600,int_stack+36490,int_stack+15150, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[42] = int_stack + 9600;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+600,int_stack+15750,int_stack+15450,10);
     Libderiv->ABCD[38] = int_stack + 600;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10200,int_stack+36940,int_stack+16200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200,10);
     Libderiv->ABCD[35] = int_stack + 10200;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+36490,int_stack+16800,int_stack+16500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack,10);
     Libderiv->ABCD[34] = int_stack + 36490;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10800,int_stack+37390,int_stack+17250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[33] = int_stack + 10800;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+37090,int_stack+37840,int_stack+17550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[32] = int_stack + 37090;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+37690,int_stack+18150,int_stack+17850, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[31] = int_stack + 37690;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11400,int_stack+38290,int_stack+18600, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[30] = int_stack + 11400;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12000,int_stack+19200,int_stack+18900,10);
     Libderiv->ABCD[26] = int_stack + 12000;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+38740,int_stack+19650,10);
     Libderiv->ABCD[25] = int_stack + 12600;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+38290,int_stack+20250,int_stack+19950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500,10);
     Libderiv->ABCD[23] = int_stack + 38290;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13200,int_stack+39190,int_stack+20700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack,10);
     Libderiv->ABCD[22] = int_stack + 13200;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+38890,int_stack+39640,int_stack+21000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[21] = int_stack + 38890;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+39490,int_stack+21600,int_stack+21300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[20] = int_stack + 39490;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13800,int_stack+40090,int_stack+22050, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[19] = int_stack + 13800;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14400,int_stack+22650,int_stack+22350, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[18] = int_stack + 14400;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+40540,int_stack+23100,10);
     Libderiv->ABCD[14] = int_stack + 1200;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+40090,int_stack+23700,int_stack+23400,10);
     Libderiv->ABCD[13] = int_stack + 40090;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+15000,int_stack+40990,int_stack+24150,10);
     Libderiv->ABCD[12] = int_stack + 15000;

}
