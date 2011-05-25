#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ddgd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|gd) integrals */

void d1hrr_order_ddgd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[2][6][11] = int_stack + 216;
 Libderiv->deriv_classes[3][4][11] = int_stack + 384;
 Libderiv->deriv_classes[3][5][11] = int_stack + 534;
 Libderiv->deriv_classes[3][6][11] = int_stack + 744;
 Libderiv->deriv_classes[4][4][11] = int_stack + 1024;
 Libderiv->deriv_classes[4][5][11] = int_stack + 1249;
 Libderiv->deriv_classes[4][6][11] = int_stack + 1564;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1984;
 Libderiv->deriv_classes[2][5][10] = int_stack + 2074;
 Libderiv->deriv_classes[2][6][10] = int_stack + 2200;
 Libderiv->deriv_classes[3][4][10] = int_stack + 2368;
 Libderiv->deriv_classes[3][5][10] = int_stack + 2518;
 Libderiv->deriv_classes[3][6][10] = int_stack + 2728;
 Libderiv->deriv_classes[4][4][10] = int_stack + 3008;
 Libderiv->deriv_classes[4][5][10] = int_stack + 3233;
 Libderiv->deriv_classes[4][6][10] = int_stack + 3548;
 Libderiv->deriv_classes[2][4][9] = int_stack + 3968;
 Libderiv->deriv_classes[2][5][9] = int_stack + 4058;
 Libderiv->deriv_classes[2][6][9] = int_stack + 4184;
 Libderiv->deriv_classes[3][4][9] = int_stack + 4352;
 Libderiv->deriv_classes[3][5][9] = int_stack + 4502;
 Libderiv->deriv_classes[3][6][9] = int_stack + 4712;
 Libderiv->deriv_classes[4][4][9] = int_stack + 4992;
 Libderiv->deriv_classes[4][5][9] = int_stack + 5217;
 Libderiv->deriv_classes[4][6][9] = int_stack + 5532;
 Libderiv->deriv_classes[2][4][8] = int_stack + 5952;
 Libderiv->deriv_classes[2][5][8] = int_stack + 6042;
 Libderiv->deriv_classes[2][6][8] = int_stack + 6168;
 Libderiv->deriv_classes[3][4][8] = int_stack + 6336;
 Libderiv->deriv_classes[3][5][8] = int_stack + 6486;
 Libderiv->deriv_classes[3][6][8] = int_stack + 6696;
 Libderiv->deriv_classes[4][4][8] = int_stack + 6976;
 Libderiv->deriv_classes[4][5][8] = int_stack + 7201;
 Libderiv->deriv_classes[4][6][8] = int_stack + 7516;
 Libderiv->deriv_classes[2][4][7] = int_stack + 7936;
 Libderiv->deriv_classes[2][5][7] = int_stack + 8026;
 Libderiv->deriv_classes[2][6][7] = int_stack + 8152;
 Libderiv->deriv_classes[3][4][7] = int_stack + 8320;
 Libderiv->deriv_classes[3][5][7] = int_stack + 8470;
 Libderiv->deriv_classes[3][6][7] = int_stack + 8680;
 Libderiv->deriv_classes[4][4][7] = int_stack + 8960;
 Libderiv->deriv_classes[4][5][7] = int_stack + 9185;
 Libderiv->deriv_classes[4][6][7] = int_stack + 9500;
 Libderiv->deriv_classes[2][4][6] = int_stack + 9920;
 Libderiv->deriv_classes[2][5][6] = int_stack + 10010;
 Libderiv->deriv_classes[2][6][6] = int_stack + 10136;
 Libderiv->deriv_classes[3][4][6] = int_stack + 10304;
 Libderiv->deriv_classes[3][5][6] = int_stack + 10454;
 Libderiv->deriv_classes[3][6][6] = int_stack + 10664;
 Libderiv->dvrr_classes[4][4] = int_stack + 10944;
 Libderiv->deriv_classes[4][4][6] = int_stack + 11169;
 Libderiv->dvrr_classes[4][5] = int_stack + 11394;
 Libderiv->deriv_classes[4][5][6] = int_stack + 11709;
 Libderiv->deriv_classes[4][6][6] = int_stack + 12024;
 Libderiv->deriv_classes[2][4][2] = int_stack + 12444;
 Libderiv->deriv_classes[2][5][2] = int_stack + 12534;
 Libderiv->deriv_classes[2][6][2] = int_stack + 12660;
 Libderiv->deriv_classes[3][4][2] = int_stack + 12828;
 Libderiv->deriv_classes[3][5][2] = int_stack + 12978;
 Libderiv->deriv_classes[3][6][2] = int_stack + 13188;
 Libderiv->deriv_classes[4][4][2] = int_stack + 13468;
 Libderiv->deriv_classes[4][5][2] = int_stack + 13693;
 Libderiv->deriv_classes[4][6][2] = int_stack + 14008;
 Libderiv->deriv_classes[2][4][1] = int_stack + 14428;
 Libderiv->deriv_classes[2][5][1] = int_stack + 14518;
 Libderiv->deriv_classes[2][6][1] = int_stack + 14644;
 Libderiv->deriv_classes[3][4][1] = int_stack + 14812;
 Libderiv->deriv_classes[3][5][1] = int_stack + 14962;
 Libderiv->deriv_classes[3][6][1] = int_stack + 15172;
 Libderiv->deriv_classes[4][4][1] = int_stack + 15452;
 Libderiv->deriv_classes[4][5][1] = int_stack + 15677;
 Libderiv->deriv_classes[4][6][1] = int_stack + 15992;
 Libderiv->dvrr_classes[2][4] = int_stack + 16412;
 Libderiv->dvrr_classes[2][5] = int_stack + 16502;
 Libderiv->dvrr_classes[2][6] = int_stack + 16628;
 Libderiv->deriv_classes[2][4][0] = int_stack + 16796;
 Libderiv->deriv_classes[2][5][0] = int_stack + 16886;
 Libderiv->deriv_classes[2][6][0] = int_stack + 17012;
 Libderiv->dvrr_classes[3][4] = int_stack + 17180;
 Libderiv->dvrr_classes[3][5] = int_stack + 17330;
 Libderiv->dvrr_classes[3][6] = int_stack + 17540;
 Libderiv->deriv_classes[3][4][0] = int_stack + 17820;
 Libderiv->deriv_classes[3][5][0] = int_stack + 17970;
 Libderiv->deriv_classes[3][6][0] = int_stack + 18180;
 Libderiv->deriv_classes[4][4][0] = int_stack + 18460;
 Libderiv->deriv_classes[4][5][0] = int_stack + 18685;
 Libderiv->deriv_classes[4][6][0] = int_stack + 19000;
 memset(int_stack,0,155360);

 Libderiv->dvrr_stack = int_stack + 55168;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ddgd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19420,int_stack+16502,int_stack+16412,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19690,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16412,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+19960,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16502,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20338,int_stack+19960,int_stack+19690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19420,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19690,int_stack+17330,int_stack+17180,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20878,int_stack+534,int_stack+384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17180,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21328,int_stack+744,int_stack+534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17330,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+21328,int_stack+20878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19690,10);
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+20878,int_stack+0,int_stack+20338,90);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+20140,int_stack+11394,int_stack+10944,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22498,int_stack+1249,int_stack+1024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10944,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23173,int_stack+1564,int_stack+1249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11394,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24118,int_stack+23173,int_stack+22498, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20140,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+25468,int_stack+24118,int_stack+0,90);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+2074,int_stack+1984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16412, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+270,int_stack+2200,int_stack+2074, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16502, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+648,int_stack+270,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19420, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+2518,int_stack+2368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17180, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1188,int_stack+2728,int_stack+2518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17330, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1818,int_stack+1188,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19690, 0.0, zero_stack,10);
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+22498,int_stack+1818,int_stack+648,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3233,int_stack+3008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10944, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+675,int_stack+3548,int_stack+3233, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11394, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24118,int_stack+675,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20140, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28168,int_stack+24118,int_stack+1818,90);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+4058,int_stack+3968, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16412, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24388,int_stack+4184,int_stack+4058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16502, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24766,int_stack+24388,int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19420, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+4502,int_stack+4352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17180, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+4712,int_stack+4502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17330, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+630,int_stack+0,int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19690, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1530,int_stack+630,int_stack+24766,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+5217,int_stack+4992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10944, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3150,int_stack+5532,int_stack+5217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11394, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4095,int_stack+3150,int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20140, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+30868,int_stack+4095,int_stack+630,90);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+6042,int_stack+5952, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16412, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24388,int_stack+6168,int_stack+6042, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24766,int_stack+24388,int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+6486,int_stack+6336, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3150,int_stack+6696,int_stack+6486, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3780,int_stack+3150,int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4680,int_stack+3780,int_stack+24766,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+7201,int_stack+6976, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+7516,int_stack+7201, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6300,int_stack+0,int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+33568,int_stack+6300,int_stack+3780,90);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6300,int_stack+8026,int_stack+7936, 0.0, zero_stack, 1.0, int_stack+16412, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6570,int_stack+8152,int_stack+8026, 0.0, zero_stack, 1.0, int_stack+16502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6948,int_stack+6570,int_stack+6300, 0.0, zero_stack, 1.0, int_stack+19420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6300,int_stack+8470,int_stack+8320, 0.0, zero_stack, 1.0, int_stack+17180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7488,int_stack+8680,int_stack+8470, 0.0, zero_stack, 1.0, int_stack+17330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24118,int_stack+7488,int_stack+6300, 0.0, zero_stack, 1.0, int_stack+19690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+36268,int_stack+24118,int_stack+6948,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6300,int_stack+9185,int_stack+8960, 0.0, zero_stack, 1.0, int_stack+10944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6975,int_stack+9500,int_stack+9185, 0.0, zero_stack, 1.0, int_stack+11394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7920,int_stack+6975,int_stack+6300, 0.0, zero_stack, 1.0, int_stack+20140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+37888,int_stack+7920,int_stack+24118,90);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+10010,int_stack+9920, 1.0, int_stack+16412, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24388,int_stack+10136,int_stack+10010, 1.0, int_stack+16502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24766,int_stack+24388,int_stack+24118, 1.0, int_stack+19420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+10454,int_stack+10304, 1.0, int_stack+17180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6300,int_stack+10664,int_stack+10454, 1.0, int_stack+17330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6930,int_stack+6300,int_stack+24118, 1.0, int_stack+19690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7830,int_stack+6930,int_stack+24766,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+11709,int_stack+11169, 1.0, int_stack+10944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9450,int_stack+12024,int_stack+11709, 1.0, int_stack+11394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10395,int_stack+9450,int_stack+24118, 1.0, int_stack+20140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+40588,int_stack+10395,int_stack+6930,90);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+20140,int_stack+16628,int_stack+16502,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+24118,int_stack+20140,int_stack+19420,6);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+20140,int_stack+17540,int_stack+17330,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9450,int_stack+20140,int_stack+19690,10);
 /*--- compute (dp|gd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+10350,int_stack+9450,int_stack+24118,90);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19420,int_stack+12534,int_stack+12444,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+19690,int_stack+12660,int_stack+12534,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+20068,int_stack+19690,int_stack+19420,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19420,int_stack+12978,int_stack+12828,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24658,int_stack+13188,int_stack+12978,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11970,int_stack+24658,int_stack+19420,10);
 /*--- compute (dp|gd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+43288,int_stack+11970,int_stack+20068, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19420,int_stack+13693,int_stack+13468,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6300,int_stack+14008,int_stack+13693,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12870,int_stack+6300,int_stack+19420,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+44908,int_stack+12870,int_stack+11970, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11970,int_stack+14518,int_stack+14428,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12240,int_stack+14644,int_stack+14518,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12618,int_stack+12240,int_stack+11970,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11970,int_stack+14962,int_stack+14812,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+13158,int_stack+15172,int_stack+14962,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+13788,int_stack+13158,int_stack+11970,10);
 /*--- compute (dp|gd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+47608,int_stack+13788,int_stack+12618, 0.0, zero_stack, 1.0, int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11970,int_stack+15677,int_stack+15452,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12645,int_stack+15992,int_stack+15677,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+14688,int_stack+12645,int_stack+11970,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+49228,int_stack+14688,int_stack+13788, 0.0, zero_stack, 1.0, int_stack+9450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11970,int_stack+16886,int_stack+16796,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12240,int_stack+17012,int_stack+16886,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12618,int_stack+12240,int_stack+11970,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11970,int_stack+17970,int_stack+17820,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+13158,int_stack+18180,int_stack+17970,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+13788,int_stack+13158,int_stack+11970,10);
 /*--- compute (dp|gd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+14688,int_stack+13788,int_stack+12618, 1.0, int_stack+24118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24118,int_stack+18685,int_stack+18460,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11970,int_stack+19000,int_stack+18685,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16308,int_stack+11970,int_stack+24118,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+17658,int_stack+16308,int_stack+13788, 1.0, int_stack+9450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (dd|gd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+51928,int_stack+25468,int_stack+20878,90);
     Libderiv->ABCD[11] = int_stack + 51928;
 /*--- compute (dd|gd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+24118,int_stack+28168,int_stack+22498,90);
     Libderiv->ABCD[10] = int_stack + 24118;
 /*--- compute (dd|gd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+27358,int_stack+30868,int_stack+1530,90);
     Libderiv->ABCD[9] = int_stack + 27358;
 /*--- compute (dd|gd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+20358,int_stack+33568,int_stack+4680,90);
     Libderiv->ABCD[8] = int_stack + 20358;
 /*--- compute (dd|gd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+30598,int_stack+37888,int_stack+36268,90);
     Libderiv->ABCD[7] = int_stack + 30598;
 /*--- compute (dd|gd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+33838,int_stack+40588,int_stack+7830,90);
     Libderiv->ABCD[6] = int_stack + 33838;
 /*--- compute (dd|gd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+37078,int_stack+44908,int_stack+43288, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[2] = int_stack + 37078;
 /*--- compute (dd|gd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+40318,int_stack+49228,int_stack+47608, 0.0, zero_stack, 1.0, int_stack+10350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[1] = int_stack + 40318;
 /*--- compute (dd|gd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+43558,int_stack+17658,int_stack+14688, 1.0, int_stack+10350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[0] = int_stack + 43558;

}
