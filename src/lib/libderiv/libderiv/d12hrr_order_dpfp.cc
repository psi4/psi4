#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_dpfp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|fp) integrals */

void d12hrr_order_dpfp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][10] = int_stack + 150;
 Libderiv->deriv_classes[3][4][9] = int_stack + 300;
 Libderiv->deriv_classes[3][4][8] = int_stack + 450;
 Libderiv->deriv_classes[3][4][7] = int_stack + 600;
 Libderiv->dvrr_classes[3][3] = int_stack + 750;
 Libderiv->deriv_classes[3][4][6] = int_stack + 850;
 Libderiv->deriv_classes[3][4][2] = int_stack + 1000;
 Libderiv->deriv_classes[3][4][1] = int_stack + 1150;
 Libderiv->dvrr_classes[2][4] = int_stack + 1300;
 Libderiv->deriv_classes[3][4][0] = int_stack + 1390;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 1540;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 1600;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 1690;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 1790;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1940;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 2000;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 2090;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 2190;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 2340;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 2400;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 2490;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 2590;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 2740;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 2800;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 2890;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 2990;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 3140;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 3200;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 3290;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 3390;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 3540;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 3600;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 3690;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 3790;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 3940;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 4000;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 4090;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 4190;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 4340;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 4400;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 4490;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 4590;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 4740;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 4800;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 4890;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 4990;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 5140;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 5200;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 5290;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 5390;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 5540;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 5600;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 5690;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 5790;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 5940;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 6000;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 6090;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 6190;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 6340;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 6400;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 6490;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 6590;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 6740;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 6800;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 6890;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 6990;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 7140;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 7200;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 7290;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 7390;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 7540;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 7600;
 Libderiv->deriv_classes[3][3][11] = int_stack + 7690;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 7790;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 7890;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 8040;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 8100;
 Libderiv->deriv_classes[3][3][10] = int_stack + 8190;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 8290;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 8390;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 8540;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 8600;
 Libderiv->deriv_classes[3][3][9] = int_stack + 8690;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 8790;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 8890;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 9040;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 9100;
 Libderiv->deriv_classes[3][3][8] = int_stack + 9190;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 9290;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 9390;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 9540;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 9600;
 Libderiv->deriv_classes[3][3][7] = int_stack + 9690;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 9790;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 9890;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 10040;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 10100;
 Libderiv->deriv_classes[3][3][6] = int_stack + 10190;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 10290;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 10390;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 10540;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 10600;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 10690;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 10790;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 10940;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 11000;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 11090;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 11190;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 11340;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 11400;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 11490;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 11590;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 11740;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 11800;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 11890;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 11990;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 12140;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 12200;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 12290;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 12390;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 12540;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 12600;
 Libderiv->deriv_classes[3][3][2] = int_stack + 12690;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 12790;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 12890;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 13040;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 13100;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 13190;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 13290;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 13440;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 13500;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 13590;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 13690;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 13840;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 13900;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 13990;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 14090;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 14240;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 14300;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 14390;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 14490;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 14640;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 14700;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 14790;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 14890;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 15040;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 15100;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 15190;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 15290;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 15440;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 15500;
 Libderiv->deriv_classes[3][3][1] = int_stack + 15590;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 15690;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 15790;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 15940;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 16000;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 16090;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 16190;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 16340;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 16400;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 16490;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 16590;
 Libderiv->deriv_classes[2][3][11] = int_stack + 16740;
 Libderiv->deriv_classes[2][4][11] = int_stack + 16800;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 16890;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 16950;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 17040;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 17140;
 Libderiv->deriv_classes[2][3][10] = int_stack + 17290;
 Libderiv->deriv_classes[2][4][10] = int_stack + 17350;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 17440;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 17500;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 17590;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 17690;
 Libderiv->deriv_classes[2][3][9] = int_stack + 17840;
 Libderiv->deriv_classes[2][4][9] = int_stack + 17900;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 17990;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 18050;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 18140;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 18240;
 Libderiv->deriv_classes[2][3][8] = int_stack + 18390;
 Libderiv->deriv_classes[2][4][8] = int_stack + 18450;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 18540;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 18600;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 18690;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 18790;
 Libderiv->deriv_classes[2][3][7] = int_stack + 18940;
 Libderiv->deriv_classes[2][4][7] = int_stack + 19000;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 19090;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 19150;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 19240;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 19340;
 Libderiv->dvrr_classes[2][3] = int_stack + 19490;
 Libderiv->deriv_classes[2][3][6] = int_stack + 19550;
 Libderiv->deriv_classes[2][4][6] = int_stack + 19610;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 19700;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 19760;
 Libderiv->deriv_classes[3][3][0] = int_stack + 19850;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 19950;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 20050;
 Libderiv->deriv_classes[2][3][2] = int_stack + 20200;
 Libderiv->deriv_classes[2][4][2] = int_stack + 20260;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 20350;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 20410;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 20500;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 20600;
 Libderiv->deriv_classes[2][3][1] = int_stack + 20750;
 Libderiv->deriv_classes[2][4][1] = int_stack + 20810;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 20900;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 20960;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 21050;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 21150;
 Libderiv->deriv_classes[2][3][0] = int_stack + 21300;
 Libderiv->deriv_classes[2][4][0] = int_stack + 21360;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 21450;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 21510;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 21600;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 21700;
 memset(int_stack,0,174800);

 Libderiv->dvrr_stack = int_stack + 33130;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_dpfp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21850,int_stack+16800,int_stack+16740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19490,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22030,int_stack+0,int_stack+7690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22330,int_stack+17350,int_stack+17290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19490, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22510,int_stack+150,int_stack+8190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+17900,int_stack+17840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19490, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22810,int_stack+300,int_stack+8690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+18450,int_stack+18390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23110,int_stack+450,int_stack+9190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+19000,int_stack+18940, 0.0, zero_stack, 1.0, int_stack+19490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23410,int_stack+600,int_stack+9690, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+19610,int_stack+19550, 1.0, int_stack+19490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23710,int_stack+850,int_stack+10190, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+1300,int_stack+19490,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24010,int_stack+20260,int_stack+20200,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24190,int_stack+1000,int_stack+12690,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+20810,int_stack+20750,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24490,int_stack+1150,int_stack+15590,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+21360,int_stack+21300,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24790,int_stack+1390,int_stack+19850,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1260,int_stack+1600,int_stack+1540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16740,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25090,int_stack+1790,int_stack+1690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7690,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1440,int_stack+2000,int_stack+1940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16740, 1.0, int_stack+17290,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1620,int_stack+2190,int_stack+2090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7690, 1.0, int_stack+8190,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1920,int_stack+2400,int_stack+2340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17290, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+2590,int_stack+2490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8190, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2400,int_stack+2800,int_stack+2740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16740, 0.0, zero_stack, 1.0, int_stack+17840,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2580,int_stack+2990,int_stack+2890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7690, 0.0, zero_stack, 1.0, int_stack+8690,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2880,int_stack+3200,int_stack+3140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17290, 1.0, int_stack+17840, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25390,int_stack+3390,int_stack+3290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8190, 1.0, int_stack+8690, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3060,int_stack+3600,int_stack+3540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17840, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3240,int_stack+3790,int_stack+3690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8690, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3540,int_stack+4000,int_stack+3940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18390,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3720,int_stack+4190,int_stack+4090, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7690, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9190,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4020,int_stack+4400,int_stack+4340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17290, 0.0, zero_stack, 1.0, int_stack+18390, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25690,int_stack+4590,int_stack+4490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8190, 0.0, zero_stack, 1.0, int_stack+9190, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4200,int_stack+4800,int_stack+4740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17840, 1.0, int_stack+18390, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4380,int_stack+4990,int_stack+4890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8690, 1.0, int_stack+9190, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4680,int_stack+5200,int_stack+5140, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+18390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4860,int_stack+5390,int_stack+5290, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5160,int_stack+5600,int_stack+5540, 0.0, zero_stack, 1.0, int_stack+16740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18940,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5340,int_stack+5790,int_stack+5690, 0.0, zero_stack, 1.0, int_stack+7690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9690,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5640,int_stack+6000,int_stack+5940, 0.0, zero_stack, 1.0, int_stack+17290, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18940, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25990,int_stack+6190,int_stack+6090, 0.0, zero_stack, 1.0, int_stack+8190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9690, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5820,int_stack+6400,int_stack+6340, 0.0, zero_stack, 1.0, int_stack+17840, 0.0, zero_stack, 1.0, int_stack+18940, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6000,int_stack+6590,int_stack+6490, 0.0, zero_stack, 1.0, int_stack+8690, 0.0, zero_stack, 1.0, int_stack+9690, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6300,int_stack+6800,int_stack+6740, 0.0, zero_stack, 1.0, int_stack+18390, 1.0, int_stack+18940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6480,int_stack+6990,int_stack+6890, 0.0, zero_stack, 1.0, int_stack+9190, 1.0, int_stack+9690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6780,int_stack+7200,int_stack+7140, 0.0, zero_stack, 2.0, int_stack+18940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6960,int_stack+7390,int_stack+7290, 0.0, zero_stack, 2.0, int_stack+9690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7260,int_stack+7600,int_stack+7540, 1.0, int_stack+16740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19550,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26290,int_stack+7890,int_stack+7790, 1.0, int_stack+7690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10190,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7440,int_stack+8100,int_stack+8040, 1.0, int_stack+17290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19550, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7620,int_stack+8390,int_stack+8290, 1.0, int_stack+8190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10190, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7920,int_stack+8600,int_stack+8540, 1.0, int_stack+17840, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19550, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8100,int_stack+8890,int_stack+8790, 1.0, int_stack+8690, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10190, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8400,int_stack+9100,int_stack+9040, 1.0, int_stack+18390, 0.0, zero_stack, 1.0, int_stack+19550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8580,int_stack+9390,int_stack+9290, 1.0, int_stack+9190, 0.0, zero_stack, 1.0, int_stack+10190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8880,int_stack+9600,int_stack+9540, 1.0, int_stack+18940, 1.0, int_stack+19550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9060,int_stack+9890,int_stack+9790, 1.0, int_stack+9690, 1.0, int_stack+10190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9360,int_stack+10100,int_stack+10040, 2.0, int_stack+19550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9540,int_stack+10390,int_stack+10290, 2.0, int_stack+10190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+19490,int_stack+10600,int_stack+10540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20200,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9840,int_stack+10790,int_stack+10690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12690,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10140,int_stack+11000,int_stack+10940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20200, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10320,int_stack+11190,int_stack+11090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12690, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10620,int_stack+11400,int_stack+11340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20200, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10800,int_stack+11590,int_stack+11490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12690, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11100,int_stack+11800,int_stack+11740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11280,int_stack+11990,int_stack+11890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11580,int_stack+12200,int_stack+12140, 0.0, zero_stack, 1.0, int_stack+20200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11760,int_stack+12390,int_stack+12290, 0.0, zero_stack, 1.0, int_stack+12690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12060,int_stack+12600,int_stack+12540, 1.0, int_stack+20200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12240,int_stack+12890,int_stack+12790, 1.0, int_stack+12690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12540,int_stack+13100,int_stack+13040,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12720,int_stack+13290,int_stack+13190,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13020,int_stack+13500,int_stack+13440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20750,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13200,int_stack+13690,int_stack+13590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15590,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13500,int_stack+13900,int_stack+13840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20750, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13680,int_stack+14090,int_stack+13990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15590, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13980,int_stack+14300,int_stack+14240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20750, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26590,int_stack+14490,int_stack+14390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15590, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14160,int_stack+14700,int_stack+14640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14340,int_stack+14890,int_stack+14790, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14640,int_stack+15100,int_stack+15040, 0.0, zero_stack, 1.0, int_stack+20750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14820,int_stack+15290,int_stack+15190, 0.0, zero_stack, 1.0, int_stack+15590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15120,int_stack+15500,int_stack+15440, 1.0, int_stack+20750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26890,int_stack+15790,int_stack+15690, 1.0, int_stack+15590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15300,int_stack+16000,int_stack+15940,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15480,int_stack+16190,int_stack+16090,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15780,int_stack+16400,int_stack+16340,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15960,int_stack+16590,int_stack+16490,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16260,int_stack+16950,int_stack+16890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21300,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16440,int_stack+17140,int_stack+17040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19850,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16740,int_stack+17500,int_stack+17440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21300, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16920,int_stack+17690,int_stack+17590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19850, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17220,int_stack+18050,int_stack+17990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21300, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17400,int_stack+18240,int_stack+18140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19850, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17700,int_stack+18600,int_stack+18540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17880,int_stack+18790,int_stack+18690, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18180,int_stack+19150,int_stack+19090, 0.0, zero_stack, 1.0, int_stack+21300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18360,int_stack+19340,int_stack+19240, 0.0, zero_stack, 1.0, int_stack+19850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18660,int_stack+19760,int_stack+19700, 1.0, int_stack+21300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18840,int_stack+20050,int_stack+19950, 1.0, int_stack+19850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19670,int_stack+20410,int_stack+20350,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19850,int_stack+20600,int_stack+20500,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+20150,int_stack+20960,int_stack+20900,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+20330,int_stack+21150,int_stack+21050,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+20630,int_stack+21510,int_stack+21450,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+20810,int_stack+21700,int_stack+21600,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+21110,int_stack+22030,int_stack+21850,30);
     Libderiv->ABCD[11] = int_stack + 21110;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+27190,int_stack+22510,int_stack+22330,30);
     Libderiv->ABCD[10] = int_stack + 27190;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+27730,int_stack+22810,int_stack+0,30);
     Libderiv->ABCD[9] = int_stack + 27730;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+22510,int_stack+23110,int_stack+180,30);
     Libderiv->ABCD[8] = int_stack + 22510;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+28270,int_stack+23410,int_stack+360,30);
     Libderiv->ABCD[7] = int_stack + 28270;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+23050,int_stack+23710,int_stack+540,30);
     Libderiv->ABCD[6] = int_stack + 23050;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+28810,int_stack+24190,int_stack+24010, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[2] = int_stack + 28810;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+29350,int_stack+24490,int_stack+900, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[1] = int_stack + 29350;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+24190,int_stack+24790,int_stack+1080, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[0] = int_stack + 24190;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+29890,int_stack+25090,int_stack+1260,30);
     Libderiv->ABCD[155] = int_stack + 29890;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+24730,int_stack+1620,int_stack+1440,30);
     Libderiv->ABCD[143] = int_stack + 24730;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1260,int_stack+2100,int_stack+1920,30);
     Libderiv->ABCD[142] = int_stack + 1260;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1800,int_stack+2580,int_stack+2400,30);
     Libderiv->ABCD[131] = int_stack + 1800;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2340,int_stack+25390,int_stack+2880,30);
     Libderiv->ABCD[130] = int_stack + 2340;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+30430,int_stack+3240,int_stack+3060,30);
     Libderiv->ABCD[129] = int_stack + 30430;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2880,int_stack+3720,int_stack+3540,30);
     Libderiv->ABCD[119] = int_stack + 2880;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3420,int_stack+25690,int_stack+4020,30);
     Libderiv->ABCD[118] = int_stack + 3420;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+25270,int_stack+4380,int_stack+4200,30);
     Libderiv->ABCD[117] = int_stack + 25270;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3960,int_stack+4860,int_stack+4680,30);
     Libderiv->ABCD[116] = int_stack + 3960;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4500,int_stack+5340,int_stack+5160,30);
     Libderiv->ABCD[107] = int_stack + 4500;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5040,int_stack+25990,int_stack+5640,30);
     Libderiv->ABCD[106] = int_stack + 5040;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+30970,int_stack+6000,int_stack+5820,30);
     Libderiv->ABCD[105] = int_stack + 30970;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5580,int_stack+6480,int_stack+6300,30);
     Libderiv->ABCD[104] = int_stack + 5580;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6120,int_stack+6960,int_stack+6780,30);
     Libderiv->ABCD[103] = int_stack + 6120;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6660,int_stack+26290,int_stack+7260,30);
     Libderiv->ABCD[95] = int_stack + 6660;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+25810,int_stack+7620,int_stack+7440,30);
     Libderiv->ABCD[94] = int_stack + 25810;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7200,int_stack+8100,int_stack+7920,30);
     Libderiv->ABCD[93] = int_stack + 7200;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7740,int_stack+8580,int_stack+8400,30);
     Libderiv->ABCD[92] = int_stack + 7740;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8280,int_stack+9060,int_stack+8880,30);
     Libderiv->ABCD[91] = int_stack + 8280;
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8820,int_stack+9540,int_stack+9360,30);
     Libderiv->ABCD[90] = int_stack + 8820;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+31510,int_stack+9840,int_stack+19490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[47] = int_stack + 31510;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9360,int_stack+10320,int_stack+10140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[46] = int_stack + 9360;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9900,int_stack+10800,int_stack+10620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[45] = int_stack + 9900;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+10440,int_stack+11280,int_stack+11100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[44] = int_stack + 10440;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+10980,int_stack+11760,int_stack+11580, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[43] = int_stack + 10980;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+11520,int_stack+12240,int_stack+12060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[42] = int_stack + 11520;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+32050,int_stack+12720,int_stack+12540, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[38] = int_stack + 32050;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+12060,int_stack+13200,int_stack+13020, 0.0, zero_stack, 1.0, int_stack+21850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[35] = int_stack + 12060;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+12600,int_stack+13680,int_stack+13500, 0.0, zero_stack, 1.0, int_stack+22330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[34] = int_stack + 12600;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+13140,int_stack+26590,int_stack+13980, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[33] = int_stack + 13140;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+26350,int_stack+14340,int_stack+14160, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[32] = int_stack + 26350;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+13680,int_stack+14820,int_stack+14640, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[31] = int_stack + 13680;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+14220,int_stack+26890,int_stack+15120, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[30] = int_stack + 14220;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+14760,int_stack+15480,int_stack+15300, 0.0, zero_stack, 1.0, int_stack+24010, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[26] = int_stack + 14760;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+32590,int_stack+15960,int_stack+15780, 0.0, zero_stack, 2.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[25] = int_stack + 32590;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+15300,int_stack+16440,int_stack+16260, 1.0, int_stack+21850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[23] = int_stack + 15300;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+15840,int_stack+16920,int_stack+16740, 1.0, int_stack+22330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[22] = int_stack + 15840;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+16380,int_stack+17400,int_stack+17220, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[21] = int_stack + 16380;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+16920,int_stack+17880,int_stack+17700, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[20] = int_stack + 16920;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+17460,int_stack+18360,int_stack+18180, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[19] = int_stack + 17460;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+18840,int_stack+18660, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[18] = int_stack + 0;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+18000,int_stack+19850,int_stack+19670, 1.0, int_stack+24010, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[14] = int_stack + 18000;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+18540,int_stack+20330,int_stack+20150, 1.0, int_stack+900, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[13] = int_stack + 18540;
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+540,int_stack+20810,int_stack+20630, 2.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[12] = int_stack + 540;

}
