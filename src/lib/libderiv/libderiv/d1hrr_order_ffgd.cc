#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ffgd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (ff|gd) integrals */

void d1hrr_order_ffgd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][6][11] = int_stack + 360;
 Libderiv->deriv_classes[4][4][11] = int_stack + 640;
 Libderiv->deriv_classes[4][5][11] = int_stack + 865;
 Libderiv->deriv_classes[4][6][11] = int_stack + 1180;
 Libderiv->deriv_classes[5][4][11] = int_stack + 1600;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1915;
 Libderiv->deriv_classes[5][6][11] = int_stack + 2356;
 Libderiv->deriv_classes[6][4][11] = int_stack + 2944;
 Libderiv->deriv_classes[6][5][11] = int_stack + 3364;
 Libderiv->deriv_classes[6][6][11] = int_stack + 3952;
 Libderiv->deriv_classes[3][4][10] = int_stack + 4736;
 Libderiv->deriv_classes[3][5][10] = int_stack + 4886;
 Libderiv->deriv_classes[3][6][10] = int_stack + 5096;
 Libderiv->deriv_classes[4][4][10] = int_stack + 5376;
 Libderiv->deriv_classes[4][5][10] = int_stack + 5601;
 Libderiv->deriv_classes[4][6][10] = int_stack + 5916;
 Libderiv->deriv_classes[5][4][10] = int_stack + 6336;
 Libderiv->deriv_classes[5][5][10] = int_stack + 6651;
 Libderiv->deriv_classes[5][6][10] = int_stack + 7092;
 Libderiv->deriv_classes[6][4][10] = int_stack + 7680;
 Libderiv->deriv_classes[6][5][10] = int_stack + 8100;
 Libderiv->deriv_classes[6][6][10] = int_stack + 8688;
 Libderiv->deriv_classes[3][4][9] = int_stack + 9472;
 Libderiv->deriv_classes[3][5][9] = int_stack + 9622;
 Libderiv->deriv_classes[3][6][9] = int_stack + 9832;
 Libderiv->deriv_classes[4][4][9] = int_stack + 10112;
 Libderiv->deriv_classes[4][5][9] = int_stack + 10337;
 Libderiv->deriv_classes[4][6][9] = int_stack + 10652;
 Libderiv->deriv_classes[5][4][9] = int_stack + 11072;
 Libderiv->deriv_classes[5][5][9] = int_stack + 11387;
 Libderiv->deriv_classes[5][6][9] = int_stack + 11828;
 Libderiv->deriv_classes[6][4][9] = int_stack + 12416;
 Libderiv->deriv_classes[6][5][9] = int_stack + 12836;
 Libderiv->deriv_classes[6][6][9] = int_stack + 13424;
 Libderiv->deriv_classes[3][4][8] = int_stack + 14208;
 Libderiv->deriv_classes[3][5][8] = int_stack + 14358;
 Libderiv->deriv_classes[3][6][8] = int_stack + 14568;
 Libderiv->deriv_classes[4][4][8] = int_stack + 14848;
 Libderiv->deriv_classes[4][5][8] = int_stack + 15073;
 Libderiv->deriv_classes[4][6][8] = int_stack + 15388;
 Libderiv->deriv_classes[5][4][8] = int_stack + 15808;
 Libderiv->deriv_classes[5][5][8] = int_stack + 16123;
 Libderiv->deriv_classes[5][6][8] = int_stack + 16564;
 Libderiv->deriv_classes[6][4][8] = int_stack + 17152;
 Libderiv->deriv_classes[6][5][8] = int_stack + 17572;
 Libderiv->deriv_classes[6][6][8] = int_stack + 18160;
 Libderiv->deriv_classes[3][4][7] = int_stack + 18944;
 Libderiv->deriv_classes[3][5][7] = int_stack + 19094;
 Libderiv->deriv_classes[3][6][7] = int_stack + 19304;
 Libderiv->deriv_classes[4][4][7] = int_stack + 19584;
 Libderiv->deriv_classes[4][5][7] = int_stack + 19809;
 Libderiv->deriv_classes[4][6][7] = int_stack + 20124;
 Libderiv->deriv_classes[5][4][7] = int_stack + 20544;
 Libderiv->deriv_classes[5][5][7] = int_stack + 20859;
 Libderiv->deriv_classes[5][6][7] = int_stack + 21300;
 Libderiv->deriv_classes[6][4][7] = int_stack + 21888;
 Libderiv->deriv_classes[6][5][7] = int_stack + 22308;
 Libderiv->deriv_classes[6][6][7] = int_stack + 22896;
 Libderiv->deriv_classes[3][4][6] = int_stack + 23680;
 Libderiv->deriv_classes[3][5][6] = int_stack + 23830;
 Libderiv->deriv_classes[3][6][6] = int_stack + 24040;
 Libderiv->deriv_classes[4][4][6] = int_stack + 24320;
 Libderiv->deriv_classes[4][5][6] = int_stack + 24545;
 Libderiv->deriv_classes[4][6][6] = int_stack + 24860;
 Libderiv->deriv_classes[5][4][6] = int_stack + 25280;
 Libderiv->deriv_classes[5][5][6] = int_stack + 25595;
 Libderiv->deriv_classes[5][6][6] = int_stack + 26036;
 Libderiv->dvrr_classes[6][4] = int_stack + 26624;
 Libderiv->deriv_classes[6][4][6] = int_stack + 27044;
 Libderiv->dvrr_classes[6][5] = int_stack + 27464;
 Libderiv->deriv_classes[6][5][6] = int_stack + 28052;
 Libderiv->deriv_classes[6][6][6] = int_stack + 28640;
 Libderiv->deriv_classes[3][4][2] = int_stack + 29424;
 Libderiv->deriv_classes[3][5][2] = int_stack + 29574;
 Libderiv->deriv_classes[3][6][2] = int_stack + 29784;
 Libderiv->deriv_classes[4][4][2] = int_stack + 30064;
 Libderiv->deriv_classes[4][5][2] = int_stack + 30289;
 Libderiv->deriv_classes[4][6][2] = int_stack + 30604;
 Libderiv->deriv_classes[5][4][2] = int_stack + 31024;
 Libderiv->deriv_classes[5][5][2] = int_stack + 31339;
 Libderiv->deriv_classes[5][6][2] = int_stack + 31780;
 Libderiv->deriv_classes[6][4][2] = int_stack + 32368;
 Libderiv->deriv_classes[6][5][2] = int_stack + 32788;
 Libderiv->deriv_classes[6][6][2] = int_stack + 33376;
 Libderiv->deriv_classes[3][4][1] = int_stack + 34160;
 Libderiv->deriv_classes[3][5][1] = int_stack + 34310;
 Libderiv->deriv_classes[3][6][1] = int_stack + 34520;
 Libderiv->deriv_classes[4][4][1] = int_stack + 34800;
 Libderiv->deriv_classes[4][5][1] = int_stack + 35025;
 Libderiv->deriv_classes[4][6][1] = int_stack + 35340;
 Libderiv->deriv_classes[5][4][1] = int_stack + 35760;
 Libderiv->deriv_classes[5][5][1] = int_stack + 36075;
 Libderiv->deriv_classes[5][6][1] = int_stack + 36516;
 Libderiv->deriv_classes[6][4][1] = int_stack + 37104;
 Libderiv->deriv_classes[6][5][1] = int_stack + 37524;
 Libderiv->deriv_classes[6][6][1] = int_stack + 38112;
 Libderiv->dvrr_classes[3][4] = int_stack + 38896;
 Libderiv->dvrr_classes[3][5] = int_stack + 39046;
 Libderiv->dvrr_classes[3][6] = int_stack + 39256;
 Libderiv->deriv_classes[3][4][0] = int_stack + 39536;
 Libderiv->deriv_classes[3][5][0] = int_stack + 39686;
 Libderiv->deriv_classes[3][6][0] = int_stack + 39896;
 Libderiv->dvrr_classes[4][4] = int_stack + 40176;
 Libderiv->dvrr_classes[4][5] = int_stack + 40401;
 Libderiv->dvrr_classes[4][6] = int_stack + 40716;
 Libderiv->deriv_classes[4][4][0] = int_stack + 41136;
 Libderiv->deriv_classes[4][5][0] = int_stack + 41361;
 Libderiv->deriv_classes[4][6][0] = int_stack + 41676;
 Libderiv->dvrr_classes[5][4] = int_stack + 42096;
 Libderiv->dvrr_classes[5][5] = int_stack + 42411;
 Libderiv->dvrr_classes[5][6] = int_stack + 42852;
 Libderiv->deriv_classes[5][4][0] = int_stack + 43440;
 Libderiv->deriv_classes[5][5][0] = int_stack + 43755;
 Libderiv->deriv_classes[5][6][0] = int_stack + 44196;
 Libderiv->deriv_classes[6][4][0] = int_stack + 44784;
 Libderiv->deriv_classes[6][5][0] = int_stack + 45204;
 Libderiv->deriv_classes[6][6][0] = int_stack + 45792;
 memset(int_stack,0,372608);

 Libderiv->dvrr_stack = int_stack + 176284;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ffgd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+46576,int_stack+39046,int_stack+38896,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47026,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38896,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47476,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39046,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48106,int_stack+47476,int_stack+47026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+47026,int_stack+40401,int_stack+40176,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49006,int_stack+865,int_stack+640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40176,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49681,int_stack+1180,int_stack+865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40401,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+49681,int_stack+49006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47026,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+49006,int_stack+0,int_stack+48106,90);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+47701,int_stack+42411,int_stack+42096,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+51706,int_stack+1915,int_stack+1600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42096,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52651,int_stack+2356,int_stack+1915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42411,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+53974,int_stack+52651,int_stack+51706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47701,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+55864,int_stack+53974,int_stack+0,90);
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+59914,int_stack+55864,int_stack+49006,90);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+27464,int_stack+26624,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+3364,int_stack+2944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26624,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+48646,int_stack+3952,int_stack+3364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27464,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+50410,int_stack+48646,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+65314,int_stack+50410,int_stack+53974,90);
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+70984,int_stack+65314,int_stack+55864,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65314,int_stack+4886,int_stack+4736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38896, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+65764,int_stack+5096,int_stack+4886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39046, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+66394,int_stack+65764,int_stack+65314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65314,int_stack+5601,int_stack+5376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40176, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+67294,int_stack+5916,int_stack+5601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40401, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+68239,int_stack+67294,int_stack+65314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47026, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1260,int_stack+68239,int_stack+66394,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65314,int_stack+6651,int_stack+6336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42096, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66259,int_stack+7092,int_stack+6651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42411, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3960,int_stack+66259,int_stack+65314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47701, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+48646,int_stack+3960,int_stack+68239,90);
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+65314,int_stack+48646,int_stack+1260,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+8100,int_stack+7680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26624, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5850,int_stack+8688,int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27464, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+52696,int_stack+5850,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+79084,int_stack+52696,int_stack+3960,90);
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+1260,int_stack+79084,int_stack+48646,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48646,int_stack+9622,int_stack+9472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38896, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49096,int_stack+9832,int_stack+9622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39046, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+49726,int_stack+49096,int_stack+48646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48646,int_stack+10337,int_stack+10112, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40176, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+50626,int_stack+10652,int_stack+10337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40401, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+51571,int_stack+50626,int_stack+48646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47026, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+52921,int_stack+51571,int_stack+49726,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48646,int_stack+11387,int_stack+11072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42096, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49591,int_stack+11828,int_stack+11387, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42411, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55621,int_stack+49591,int_stack+48646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47701, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79084,int_stack+55621,int_stack+51571,90);
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+83134,int_stack+79084,int_stack+52921,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48646,int_stack+12836,int_stack+12416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26624, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49906,int_stack+13424,int_stack+12836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27464, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+51670,int_stack+49906,int_stack+48646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+88534,int_stack+51670,int_stack+55621,90);
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+48646,int_stack+88534,int_stack+79084,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+14358,int_stack+14208, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+79534,int_stack+14568,int_stack+14358, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+80164,int_stack+79534,int_stack+79084, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+15073,int_stack+14848, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81064,int_stack+15388,int_stack+15073, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+88534,int_stack+81064,int_stack+79084, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+89884,int_stack+88534,int_stack+80164,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+16123,int_stack+15808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80029,int_stack+16564,int_stack+16123, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42411, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+92584,int_stack+80029,int_stack+79084, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79084,int_stack+92584,int_stack+88534,90);
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+9360,int_stack+79084,int_stack+89884,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+88534,int_stack+17572,int_stack+17152, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+89794,int_stack+18160,int_stack+17572, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56746,int_stack+89794,int_stack+88534, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+94474,int_stack+56746,int_stack+92584,90);
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+100144,int_stack+94474,int_stack+79084,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+19094,int_stack+18944, 0.0, zero_stack, 1.0, int_stack+38896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+79534,int_stack+19304,int_stack+19094, 0.0, zero_stack, 1.0, int_stack+39046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+80164,int_stack+79534,int_stack+79084, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+19809,int_stack+19584, 0.0, zero_stack, 1.0, int_stack+40176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81064,int_stack+20124,int_stack+19809, 0.0, zero_stack, 1.0, int_stack+40401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56746,int_stack+81064,int_stack+79084, 0.0, zero_stack, 1.0, int_stack+47026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+88534,int_stack+56746,int_stack+80164,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+20859,int_stack+20544, 0.0, zero_stack, 1.0, int_stack+42096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80029,int_stack+21300,int_stack+20859, 0.0, zero_stack, 1.0, int_stack+42411, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+91234,int_stack+80029,int_stack+79084, 0.0, zero_stack, 1.0, int_stack+47701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79084,int_stack+91234,int_stack+56746,90);
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+93124,int_stack+79084,int_stack+88534,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+88534,int_stack+22308,int_stack+21888, 0.0, zero_stack, 1.0, int_stack+26624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+56746,int_stack+22896,int_stack+22308, 0.0, zero_stack, 1.0, int_stack+27464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14760,int_stack+56746,int_stack+88534, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+17280,int_stack+14760,int_stack+91234,90);
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+108244,int_stack+17280,int_stack+79084,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+23830,int_stack+23680, 1.0, int_stack+38896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+79534,int_stack+24040,int_stack+23830, 1.0, int_stack+39046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+80164,int_stack+79534,int_stack+79084, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+24545,int_stack+24320, 1.0, int_stack+40176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81064,int_stack+24860,int_stack+24545, 1.0, int_stack+40401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14760,int_stack+81064,int_stack+79084, 1.0, int_stack+47026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+16110,int_stack+14760,int_stack+80164,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79084,int_stack+25595,int_stack+25280, 1.0, int_stack+42096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80029,int_stack+26036,int_stack+25595, 1.0, int_stack+42411, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18810,int_stack+80029,int_stack+79084, 1.0, int_stack+47701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79084,int_stack+18810,int_stack+14760,90);
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+20700,int_stack+79084,int_stack+16110,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14760,int_stack+28052,int_stack+27044, 1.0, int_stack+26624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16020,int_stack+28640,int_stack+28052, 1.0, int_stack+27464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+26100,int_stack+16020,int_stack+14760, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+116344,int_stack+26100,int_stack+18810,90);
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+122014,int_stack+116344,int_stack+79084,90);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+79084,int_stack+39256,int_stack+39046,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+79714,int_stack+79084,int_stack+46576,10);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+80614,int_stack+40716,int_stack+40401,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81559,int_stack+80614,int_stack+47026,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+116344,int_stack+81559,int_stack+79714,90);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+119044,int_stack+42852,int_stack+42411,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+26100,int_stack+119044,int_stack+47701,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+14760,int_stack+26100,int_stack+81559,90);
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+130114,int_stack+14760,int_stack+116344,90);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+119044,int_stack+29574,int_stack+29424,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+79084,int_stack+29784,int_stack+29574,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+119494,int_stack+79084,int_stack+119044,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+120394,int_stack+30289,int_stack+30064,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+121069,int_stack+30604,int_stack+30289,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+46576,int_stack+121069,int_stack+120394,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+27990,int_stack+46576,int_stack+119494, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+80614,int_stack+31339,int_stack+31024,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+119044,int_stack+31780,int_stack+31339,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18810,int_stack+119044,int_stack+80614,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+88534,int_stack+18810,int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (fd|gd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+135514,int_stack+88534,int_stack+27990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+32788,int_stack+32368,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+27990,int_stack+33376,int_stack+32788,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+29754,int_stack+27990,int_stack+0,28);
 /*--- compute (hp|gd) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+140914,int_stack+29754,int_stack+18810, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (gd|gd) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+146584,int_stack+140914,int_stack+88534, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+88534,int_stack+34310,int_stack+34160,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+79084,int_stack+34520,int_stack+34310,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+88984,int_stack+79084,int_stack+88534,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+89884,int_stack+35025,int_stack+34800,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+80614,int_stack+35340,int_stack+35025,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+90559,int_stack+80614,int_stack+89884,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+140914,int_stack+90559,int_stack+88984, 0.0, zero_stack, 1.0, int_stack+79714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+80614,int_stack+36075,int_stack+35760,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+88534,int_stack+36516,int_stack+36075,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18810,int_stack+88534,int_stack+80614,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+27990,int_stack+18810,int_stack+90559, 0.0, zero_stack, 1.0, int_stack+81559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (fd|gd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+154684,int_stack+27990,int_stack+140914, 0.0, zero_stack, 1.0, int_stack+116344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+37524,int_stack+37104,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+140914,int_stack+38112,int_stack+37524,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+142678,int_stack+140914,int_stack+0,28);
 /*--- compute (hp|gd) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+32040,int_stack+142678,int_stack+18810, 0.0, zero_stack, 1.0, int_stack+26100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (gd|gd) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+160084,int_stack+32040,int_stack+27990, 0.0, zero_stack, 1.0, int_stack+14760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27990,int_stack+39686,int_stack+39536,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+79084,int_stack+39896,int_stack+39686,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+28440,int_stack+79084,int_stack+27990,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+29340,int_stack+41361,int_stack+41136,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+80614,int_stack+41676,int_stack+41361,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+30015,int_stack+80614,int_stack+29340,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+31365,int_stack+30015,int_stack+28440, 1.0, int_stack+79714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27990,int_stack+43755,int_stack+43440,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+79084,int_stack+44196,int_stack+43755,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18810,int_stack+79084,int_stack+27990,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+34065,int_stack+18810,int_stack+30015, 1.0, int_stack+81559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (fd|gd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+38115,int_stack+34065,int_stack+31365, 1.0, int_stack+116344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+45204,int_stack+44784,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+116344,int_stack+45792,int_stack+45204,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+118108,int_stack+116344,int_stack+0,28);
 /*--- compute (hp|gd) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+140914,int_stack+118108,int_stack+18810, 1.0, int_stack+26100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (gd|gd) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+168184,int_stack+140914,int_stack+34065, 1.0, int_stack+14760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (ff|gd) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+26100,int_stack+70984,int_stack+59914,90);
     Libderiv->ABCD[11] = int_stack + 26100;
 /*--- compute (ff|gd) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+70714,int_stack+1260,int_stack+65314,90);
     Libderiv->ABCD[10] = int_stack + 70714;
 /*--- compute (ff|gd) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+0,int_stack+48646,int_stack+83134,90);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (ff|gd) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+43515,int_stack+100144,int_stack+9360,90);
     Libderiv->ABCD[8] = int_stack + 43515;
 /*--- compute (ff|gd) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+9000,int_stack+108244,int_stack+93124,90);
     Libderiv->ABCD[7] = int_stack + 9000;
 /*--- compute (ff|gd) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+52515,int_stack+122014,int_stack+20700,90);
     Libderiv->ABCD[6] = int_stack + 52515;
 /*--- compute (ff|gd) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+61515,int_stack+146584,int_stack+135514, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[2] = int_stack + 61515;
 /*--- compute (ff|gd) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+135514,int_stack+160084,int_stack+154684, 0.0, zero_stack, 1.0, int_stack+130114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[1] = int_stack + 135514;
 /*--- compute (ff|gd) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+144514,int_stack+168184,int_stack+38115, 1.0, int_stack+130114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[0] = int_stack + 144514;

}
