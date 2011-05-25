#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gdff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gd|ff) integrals */

void d1hrr_order_gdff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][3][11] = int_stack + 0;
 Libderiv->deriv_classes[4][4][11] = int_stack + 150;
 Libderiv->deriv_classes[4][5][11] = int_stack + 375;
 Libderiv->deriv_classes[4][6][11] = int_stack + 690;
 Libderiv->deriv_classes[5][3][11] = int_stack + 1110;
 Libderiv->deriv_classes[5][4][11] = int_stack + 1320;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1635;
 Libderiv->deriv_classes[5][6][11] = int_stack + 2076;
 Libderiv->deriv_classes[6][3][11] = int_stack + 2664;
 Libderiv->deriv_classes[6][4][11] = int_stack + 2944;
 Libderiv->deriv_classes[6][5][11] = int_stack + 3364;
 Libderiv->deriv_classes[6][6][11] = int_stack + 3952;
 Libderiv->deriv_classes[4][3][10] = int_stack + 4736;
 Libderiv->deriv_classes[4][4][10] = int_stack + 4886;
 Libderiv->deriv_classes[4][5][10] = int_stack + 5111;
 Libderiv->deriv_classes[4][6][10] = int_stack + 5426;
 Libderiv->deriv_classes[5][3][10] = int_stack + 5846;
 Libderiv->deriv_classes[5][4][10] = int_stack + 6056;
 Libderiv->deriv_classes[5][5][10] = int_stack + 6371;
 Libderiv->deriv_classes[5][6][10] = int_stack + 6812;
 Libderiv->deriv_classes[6][3][10] = int_stack + 7400;
 Libderiv->deriv_classes[6][4][10] = int_stack + 7680;
 Libderiv->deriv_classes[6][5][10] = int_stack + 8100;
 Libderiv->deriv_classes[6][6][10] = int_stack + 8688;
 Libderiv->deriv_classes[4][3][9] = int_stack + 9472;
 Libderiv->deriv_classes[4][4][9] = int_stack + 9622;
 Libderiv->deriv_classes[4][5][9] = int_stack + 9847;
 Libderiv->deriv_classes[4][6][9] = int_stack + 10162;
 Libderiv->deriv_classes[5][3][9] = int_stack + 10582;
 Libderiv->deriv_classes[5][4][9] = int_stack + 10792;
 Libderiv->deriv_classes[5][5][9] = int_stack + 11107;
 Libderiv->deriv_classes[5][6][9] = int_stack + 11548;
 Libderiv->deriv_classes[6][3][9] = int_stack + 12136;
 Libderiv->deriv_classes[6][4][9] = int_stack + 12416;
 Libderiv->deriv_classes[6][5][9] = int_stack + 12836;
 Libderiv->deriv_classes[6][6][9] = int_stack + 13424;
 Libderiv->deriv_classes[4][3][8] = int_stack + 14208;
 Libderiv->deriv_classes[4][4][8] = int_stack + 14358;
 Libderiv->deriv_classes[4][5][8] = int_stack + 14583;
 Libderiv->deriv_classes[4][6][8] = int_stack + 14898;
 Libderiv->deriv_classes[5][3][8] = int_stack + 15318;
 Libderiv->deriv_classes[5][4][8] = int_stack + 15528;
 Libderiv->deriv_classes[5][5][8] = int_stack + 15843;
 Libderiv->deriv_classes[5][6][8] = int_stack + 16284;
 Libderiv->deriv_classes[6][3][8] = int_stack + 16872;
 Libderiv->deriv_classes[6][4][8] = int_stack + 17152;
 Libderiv->deriv_classes[6][5][8] = int_stack + 17572;
 Libderiv->deriv_classes[6][6][8] = int_stack + 18160;
 Libderiv->deriv_classes[4][3][7] = int_stack + 18944;
 Libderiv->deriv_classes[4][4][7] = int_stack + 19094;
 Libderiv->deriv_classes[4][5][7] = int_stack + 19319;
 Libderiv->deriv_classes[4][6][7] = int_stack + 19634;
 Libderiv->deriv_classes[5][3][7] = int_stack + 20054;
 Libderiv->deriv_classes[5][4][7] = int_stack + 20264;
 Libderiv->deriv_classes[5][5][7] = int_stack + 20579;
 Libderiv->deriv_classes[5][6][7] = int_stack + 21020;
 Libderiv->deriv_classes[6][3][7] = int_stack + 21608;
 Libderiv->deriv_classes[6][4][7] = int_stack + 21888;
 Libderiv->deriv_classes[6][5][7] = int_stack + 22308;
 Libderiv->deriv_classes[6][6][7] = int_stack + 22896;
 Libderiv->deriv_classes[4][3][6] = int_stack + 23680;
 Libderiv->deriv_classes[4][4][6] = int_stack + 23830;
 Libderiv->deriv_classes[4][5][6] = int_stack + 24055;
 Libderiv->deriv_classes[4][6][6] = int_stack + 24370;
 Libderiv->deriv_classes[5][3][6] = int_stack + 24790;
 Libderiv->deriv_classes[5][4][6] = int_stack + 25000;
 Libderiv->deriv_classes[5][5][6] = int_stack + 25315;
 Libderiv->deriv_classes[5][6][6] = int_stack + 25756;
 Libderiv->dvrr_classes[6][3] = int_stack + 26344;
 Libderiv->deriv_classes[6][3][6] = int_stack + 26624;
 Libderiv->dvrr_classes[6][4] = int_stack + 26904;
 Libderiv->deriv_classes[6][4][6] = int_stack + 27324;
 Libderiv->dvrr_classes[6][5] = int_stack + 27744;
 Libderiv->deriv_classes[6][5][6] = int_stack + 28332;
 Libderiv->deriv_classes[6][6][6] = int_stack + 28920;
 Libderiv->deriv_classes[4][3][2] = int_stack + 29704;
 Libderiv->deriv_classes[4][4][2] = int_stack + 29854;
 Libderiv->deriv_classes[4][5][2] = int_stack + 30079;
 Libderiv->deriv_classes[4][6][2] = int_stack + 30394;
 Libderiv->deriv_classes[5][3][2] = int_stack + 30814;
 Libderiv->deriv_classes[5][4][2] = int_stack + 31024;
 Libderiv->deriv_classes[5][5][2] = int_stack + 31339;
 Libderiv->deriv_classes[5][6][2] = int_stack + 31780;
 Libderiv->deriv_classes[6][3][2] = int_stack + 32368;
 Libderiv->deriv_classes[6][4][2] = int_stack + 32648;
 Libderiv->deriv_classes[6][5][2] = int_stack + 33068;
 Libderiv->deriv_classes[6][6][2] = int_stack + 33656;
 Libderiv->deriv_classes[4][3][1] = int_stack + 34440;
 Libderiv->deriv_classes[4][4][1] = int_stack + 34590;
 Libderiv->deriv_classes[4][5][1] = int_stack + 34815;
 Libderiv->deriv_classes[4][6][1] = int_stack + 35130;
 Libderiv->deriv_classes[5][3][1] = int_stack + 35550;
 Libderiv->deriv_classes[5][4][1] = int_stack + 35760;
 Libderiv->deriv_classes[5][5][1] = int_stack + 36075;
 Libderiv->deriv_classes[5][6][1] = int_stack + 36516;
 Libderiv->deriv_classes[6][3][1] = int_stack + 37104;
 Libderiv->deriv_classes[6][4][1] = int_stack + 37384;
 Libderiv->deriv_classes[6][5][1] = int_stack + 37804;
 Libderiv->deriv_classes[6][6][1] = int_stack + 38392;
 Libderiv->dvrr_classes[4][3] = int_stack + 39176;
 Libderiv->dvrr_classes[4][4] = int_stack + 39326;
 Libderiv->dvrr_classes[4][5] = int_stack + 39551;
 Libderiv->dvrr_classes[4][6] = int_stack + 39866;
 Libderiv->deriv_classes[4][3][0] = int_stack + 40286;
 Libderiv->deriv_classes[4][4][0] = int_stack + 40436;
 Libderiv->deriv_classes[4][5][0] = int_stack + 40661;
 Libderiv->deriv_classes[4][6][0] = int_stack + 40976;
 Libderiv->dvrr_classes[5][3] = int_stack + 41396;
 Libderiv->dvrr_classes[5][4] = int_stack + 41606;
 Libderiv->dvrr_classes[5][5] = int_stack + 41921;
 Libderiv->dvrr_classes[5][6] = int_stack + 42362;
 Libderiv->deriv_classes[5][3][0] = int_stack + 42950;
 Libderiv->deriv_classes[5][4][0] = int_stack + 43160;
 Libderiv->deriv_classes[5][5][0] = int_stack + 43475;
 Libderiv->deriv_classes[5][6][0] = int_stack + 43916;
 Libderiv->deriv_classes[6][3][0] = int_stack + 44504;
 Libderiv->deriv_classes[6][4][0] = int_stack + 44784;
 Libderiv->deriv_classes[6][5][0] = int_stack + 45204;
 Libderiv->deriv_classes[6][6][0] = int_stack + 45792;
 memset(int_stack,0,372608);

 Libderiv->dvrr_stack = int_stack + 124941;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gdff(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+46576,int_stack+39326,int_stack+39176,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+47026,int_stack+39551,int_stack+39326,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+47701,int_stack+47026,int_stack+46576,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+48601,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39176,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49051,int_stack+375,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39326,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+49726,int_stack+49051,int_stack+48601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+50626,int_stack+690,int_stack+375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39551,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+51571,int_stack+50626,int_stack+49051, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47026,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+52921,int_stack+51571,int_stack+49726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47701,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+48601,int_stack+41606,int_stack+41396,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+49231,int_stack+41921,int_stack+41606,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+50176,int_stack+49231,int_stack+48601,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51436,int_stack+1320,int_stack+1110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41396,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1635,int_stack+1320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41606,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+54421,int_stack+0,int_stack+51436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48601,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+51436,int_stack+2076,int_stack+1635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41921,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55681,int_stack+51436,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49231,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+55681,int_stack+54421, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50176,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+54421,int_stack+0,int_stack+52921,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51436,int_stack+26904,int_stack+26344,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52276,int_stack+27744,int_stack+26904,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+58921,int_stack+52276,int_stack+51436,28);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53536,int_stack+2944,int_stack+2664, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26344,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+60601,int_stack+3364,int_stack+2944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26904,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+61861,int_stack+60601,int_stack+53536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51436,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+63541,int_stack+3952,int_stack+3364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27744,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2100,int_stack+63541,int_stack+60601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52276,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+63541,int_stack+2100,int_stack+61861, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58921,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+66341,int_stack+63541,int_stack+0,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+4886,int_stack+4736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39176, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+5111,int_stack+4886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39326, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1125,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2025,int_stack+5426,int_stack+5111, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39551, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2970,int_stack+2025,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47026, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4320,int_stack+2970,int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47701, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+6056,int_stack+5846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41396, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+630,int_stack+6371,int_stack+6056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41606, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1575,int_stack+630,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48601, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2835,int_stack+6812,int_stack+6371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41921, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+60601,int_stack+2835,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49231, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+62491,int_stack+60601,int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50176, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+72641,int_stack+62491,int_stack+4320,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60601,int_stack+7680,int_stack+7400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26344, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+8100,int_stack+7680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26904, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1260,int_stack+0,int_stack+60601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51436, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60601,int_stack+8688,int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27744, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2940,int_stack+60601,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52276, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5460,int_stack+2940,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58921, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+77141,int_stack+5460,int_stack+62491,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+9622,int_stack+9472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39176, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+9847,int_stack+9622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39326, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1125,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2025,int_stack+10162,int_stack+9847, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39551, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2970,int_stack+2025,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47026, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4320,int_stack+2970,int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47701, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+10792,int_stack+10582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41396, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+630,int_stack+11107,int_stack+10792, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41606, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1575,int_stack+630,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48601, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2835,int_stack+11548,int_stack+11107, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41921, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5820,int_stack+2835,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49231, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7710,int_stack+5820,int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50176, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+60601,int_stack+7710,int_stack+4320,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+12416,int_stack+12136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26344, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+840,int_stack+12836,int_stack+12416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26904, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+840,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51436, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3780,int_stack+13424,int_stack+12836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27744, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9810,int_stack+3780,int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52276, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3780,int_stack+9810,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58921, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+83441,int_stack+3780,int_stack+7710,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+14358,int_stack+14208, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+14583,int_stack+14358, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1125,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2025,int_stack+14898,int_stack+14583, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2970,int_stack+2025,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4320,int_stack+2970,int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+15528,int_stack+15318, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+630,int_stack+15843,int_stack+15528, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1575,int_stack+630,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2835,int_stack+16284,int_stack+15843, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5820,int_stack+2835,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49231, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7710,int_stack+5820,int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+9810,int_stack+7710,int_stack+4320,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+17152,int_stack+16872, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+840,int_stack+17572,int_stack+17152, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+840,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3780,int_stack+18160,int_stack+17572, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27744, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14310,int_stack+3780,int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3780,int_stack+14310,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+89741,int_stack+3780,int_stack+7710,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14310,int_stack+19094,int_stack+18944, 0.0, zero_stack, 1.0, int_stack+39176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14760,int_stack+19319,int_stack+19094, 0.0, zero_stack, 1.0, int_stack+39326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15435,int_stack+14760,int_stack+14310, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16335,int_stack+19634,int_stack+19319, 0.0, zero_stack, 1.0, int_stack+39551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17280,int_stack+16335,int_stack+14760, 0.0, zero_stack, 1.0, int_stack+47026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+17280,int_stack+15435, 0.0, zero_stack, 1.0, int_stack+47701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14310,int_stack+20264,int_stack+20054, 0.0, zero_stack, 1.0, int_stack+41396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14940,int_stack+20579,int_stack+20264, 0.0, zero_stack, 1.0, int_stack+41606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15885,int_stack+14940,int_stack+14310, 0.0, zero_stack, 1.0, int_stack+48601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17145,int_stack+21020,int_stack+20579, 0.0, zero_stack, 1.0, int_stack+41921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18468,int_stack+17145,int_stack+14940, 0.0, zero_stack, 1.0, int_stack+49231, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1500,int_stack+18468,int_stack+15885, 0.0, zero_stack, 1.0, int_stack+50176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+14310,int_stack+1500,int_stack+0,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+21888,int_stack+21608, 0.0, zero_stack, 1.0, int_stack+26344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18810,int_stack+22308,int_stack+21888, 0.0, zero_stack, 1.0, int_stack+26904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20070,int_stack+18810,int_stack+0, 0.0, zero_stack, 1.0, int_stack+51436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3600,int_stack+22896,int_stack+22308, 0.0, zero_stack, 1.0, int_stack+27744, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5364,int_stack+3600,int_stack+18810, 0.0, zero_stack, 1.0, int_stack+52276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+96041,int_stack+5364,int_stack+20070, 0.0, zero_stack, 1.0, int_stack+58921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+98841,int_stack+96041,int_stack+1500,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+96041,int_stack+23830,int_stack+23680, 1.0, int_stack+39176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+96491,int_stack+24055,int_stack+23830, 1.0, int_stack+39326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97166,int_stack+96491,int_stack+96041, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18810,int_stack+24370,int_stack+24055, 1.0, int_stack+39551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19755,int_stack+18810,int_stack+96491, 1.0, int_stack+47026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+21105,int_stack+19755,int_stack+97166, 1.0, int_stack+47701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18810,int_stack+25000,int_stack+24790, 1.0, int_stack+41396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19440,int_stack+25315,int_stack+25000, 1.0, int_stack+41606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+96041,int_stack+19440,int_stack+18810, 1.0, int_stack+48601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+97301,int_stack+25756,int_stack+25315, 1.0, int_stack+41921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22605,int_stack+97301,int_stack+19440, 1.0, int_stack+49231, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18810,int_stack+22605,int_stack+96041, 1.0, int_stack+50176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+18810,int_stack+21105,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+96041,int_stack+27324,int_stack+26624, 1.0, int_stack+26344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+96881,int_stack+28332,int_stack+27324, 1.0, int_stack+26904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20910,int_stack+96881,int_stack+96041, 1.0, int_stack+51436, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22590,int_stack+28920,int_stack+28332, 1.0, int_stack+27744, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24354,int_stack+22590,int_stack+96881, 1.0, int_stack+52276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+96041,int_stack+24354,int_stack+20910, 1.0, int_stack+58921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+20910,int_stack+96041,int_stack+18810,100);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+18810,int_stack+39866,int_stack+39551,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+96041,int_stack+18810,int_stack+47026,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+18810,int_stack+96041,int_stack+47701,15);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+96041,int_stack+42362,int_stack+41921,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+27210,int_stack+96041,int_stack+49231,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+96041,int_stack+27210,int_stack+50176,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+46576,int_stack+96041,int_stack+18810,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27210,int_stack+29854,int_stack+29704,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27660,int_stack+30079,int_stack+29854,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+28335,int_stack+27660,int_stack+27210,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+58921,int_stack+30394,int_stack+30079,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+29235,int_stack+58921,int_stack+27660,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+58921,int_stack+29235,int_stack+28335,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27210,int_stack+31024,int_stack+30814,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27840,int_stack+31339,int_stack+31024,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+28785,int_stack+27840,int_stack+27210,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+41396,int_stack+31780,int_stack+31339,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+30045,int_stack+41396,int_stack+27840,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+51076,int_stack+30045,int_stack+28785,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+27210,int_stack+51076,int_stack+58921, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+58921,int_stack+32648,int_stack+32368,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+41396,int_stack+33068,int_stack+32648,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4500,int_stack+41396,int_stack+58921,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6180,int_stack+33656,int_stack+33068,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+31710,int_stack+6180,int_stack+41396,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+6180,int_stack+31710,int_stack+4500,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+105141,int_stack+6180,int_stack+51076, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96041, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51076,int_stack+34590,int_stack+34440,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51526,int_stack+34815,int_stack+34590,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+52201,int_stack+51526,int_stack+51076,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+53101,int_stack+35130,int_stack+34815,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4500,int_stack+53101,int_stack+51526,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+5850,int_stack+4500,int_stack+52201,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+35760,int_stack+35550,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51076,int_stack+36075,int_stack+35760,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+52021,int_stack+51076,int_stack+4500,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4500,int_stack+36516,int_stack+36075,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+7350,int_stack+4500,int_stack+51076,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+31710,int_stack+7350,int_stack+52021,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+111441,int_stack+31710,int_stack+5850, 0.0, zero_stack, 1.0, int_stack+18810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51076,int_stack+37384,int_stack+37104,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51916,int_stack+37804,int_stack+37384,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+58921,int_stack+51916,int_stack+51076,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4500,int_stack+38392,int_stack+37804,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6264,int_stack+4500,int_stack+51916,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+51076,int_stack+6264,int_stack+58921,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+33810,int_stack+51076,int_stack+31710, 0.0, zero_stack, 1.0, int_stack+96041, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+31710,int_stack+40436,int_stack+40286,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+32160,int_stack+40661,int_stack+40436,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+32835,int_stack+32160,int_stack+31710,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+51076,int_stack+40976,int_stack+40661,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+52021,int_stack+51076,int_stack+32160,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+58921,int_stack+52021,int_stack+32835,15);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51076,int_stack+43160,int_stack+42950,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51706,int_stack+43475,int_stack+43160,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+52651,int_stack+51706,int_stack+51076,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+31710,int_stack+43916,int_stack+43475,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4500,int_stack+31710,int_stack+51706,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+31710,int_stack+4500,int_stack+52651,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+4500,int_stack+31710,int_stack+58921, 1.0, int_stack+18810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+18810,int_stack+44784,int_stack+44504,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19650,int_stack+45204,int_stack+44784,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+58921,int_stack+19650,int_stack+18810,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+51076,int_stack+45792,int_stack+45204,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+40110,int_stack+51076,int_stack+19650,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+51076,int_stack+40110,int_stack+58921,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+40110,int_stack+51076,int_stack+31710, 1.0, int_stack+96041, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+115941,int_stack+66341,int_stack+54421,100);
     Libderiv->ABCD[11] = int_stack + 115941;
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+51076,int_stack+77141,int_stack+72641,100);
     Libderiv->ABCD[10] = int_stack + 51076;
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+65101,int_stack+83441,int_stack+60601,100);
     Libderiv->ABCD[9] = int_stack + 65101;
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+74101,int_stack+89741,int_stack+9810,100);
     Libderiv->ABCD[8] = int_stack + 74101;
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+83101,int_stack+98841,int_stack+14310,100);
     Libderiv->ABCD[7] = int_stack + 83101;
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+9000,int_stack+20910,int_stack+0,100);
     Libderiv->ABCD[6] = int_stack + 9000;
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+18000,int_stack+105141,int_stack+27210, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 18000;
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+92101,int_stack+33810,int_stack+111441, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 92101;
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+27000,int_stack+40110,int_stack+4500, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 27000;

}
