#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|gg) integrals */

void d1hrr_order_ppgg(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[1][5][11] = int_stack + 45;
 Libderiv->deriv_classes[1][6][11] = int_stack + 108;
 Libderiv->deriv_classes[1][7][11] = int_stack + 192;
 Libderiv->deriv_classes[1][8][11] = int_stack + 300;
 Libderiv->deriv_classes[2][4][11] = int_stack + 435;
 Libderiv->deriv_classes[2][5][11] = int_stack + 525;
 Libderiv->deriv_classes[2][6][11] = int_stack + 651;
 Libderiv->deriv_classes[2][7][11] = int_stack + 819;
 Libderiv->deriv_classes[2][8][11] = int_stack + 1035;
 Libderiv->deriv_classes[1][4][10] = int_stack + 1305;
 Libderiv->deriv_classes[1][5][10] = int_stack + 1350;
 Libderiv->deriv_classes[1][6][10] = int_stack + 1413;
 Libderiv->deriv_classes[1][7][10] = int_stack + 1497;
 Libderiv->deriv_classes[1][8][10] = int_stack + 1605;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1740;
 Libderiv->deriv_classes[2][5][10] = int_stack + 1830;
 Libderiv->deriv_classes[2][6][10] = int_stack + 1956;
 Libderiv->deriv_classes[2][7][10] = int_stack + 2124;
 Libderiv->deriv_classes[2][8][10] = int_stack + 2340;
 Libderiv->deriv_classes[1][4][9] = int_stack + 2610;
 Libderiv->deriv_classes[1][5][9] = int_stack + 2655;
 Libderiv->deriv_classes[1][6][9] = int_stack + 2718;
 Libderiv->deriv_classes[1][7][9] = int_stack + 2802;
 Libderiv->deriv_classes[1][8][9] = int_stack + 2910;
 Libderiv->deriv_classes[2][4][9] = int_stack + 3045;
 Libderiv->deriv_classes[2][5][9] = int_stack + 3135;
 Libderiv->deriv_classes[2][6][9] = int_stack + 3261;
 Libderiv->deriv_classes[2][7][9] = int_stack + 3429;
 Libderiv->deriv_classes[2][8][9] = int_stack + 3645;
 Libderiv->deriv_classes[1][4][8] = int_stack + 3915;
 Libderiv->deriv_classes[1][5][8] = int_stack + 3960;
 Libderiv->deriv_classes[1][6][8] = int_stack + 4023;
 Libderiv->deriv_classes[1][7][8] = int_stack + 4107;
 Libderiv->deriv_classes[1][8][8] = int_stack + 4215;
 Libderiv->deriv_classes[2][4][8] = int_stack + 4350;
 Libderiv->deriv_classes[2][5][8] = int_stack + 4440;
 Libderiv->deriv_classes[2][6][8] = int_stack + 4566;
 Libderiv->deriv_classes[2][7][8] = int_stack + 4734;
 Libderiv->deriv_classes[2][8][8] = int_stack + 4950;
 Libderiv->deriv_classes[1][4][7] = int_stack + 5220;
 Libderiv->deriv_classes[1][5][7] = int_stack + 5265;
 Libderiv->deriv_classes[1][6][7] = int_stack + 5328;
 Libderiv->deriv_classes[1][7][7] = int_stack + 5412;
 Libderiv->deriv_classes[1][8][7] = int_stack + 5520;
 Libderiv->deriv_classes[2][4][7] = int_stack + 5655;
 Libderiv->deriv_classes[2][5][7] = int_stack + 5745;
 Libderiv->deriv_classes[2][6][7] = int_stack + 5871;
 Libderiv->deriv_classes[2][7][7] = int_stack + 6039;
 Libderiv->deriv_classes[2][8][7] = int_stack + 6255;
 Libderiv->deriv_classes[1][4][6] = int_stack + 6525;
 Libderiv->deriv_classes[1][5][6] = int_stack + 6570;
 Libderiv->deriv_classes[1][6][6] = int_stack + 6633;
 Libderiv->deriv_classes[1][7][6] = int_stack + 6717;
 Libderiv->deriv_classes[1][8][6] = int_stack + 6825;
 Libderiv->dvrr_classes[2][4] = int_stack + 6960;
 Libderiv->deriv_classes[2][4][6] = int_stack + 7050;
 Libderiv->dvrr_classes[2][5] = int_stack + 7140;
 Libderiv->deriv_classes[2][5][6] = int_stack + 7266;
 Libderiv->dvrr_classes[2][6] = int_stack + 7392;
 Libderiv->deriv_classes[2][6][6] = int_stack + 7560;
 Libderiv->dvrr_classes[2][7] = int_stack + 7728;
 Libderiv->deriv_classes[2][7][6] = int_stack + 7944;
 Libderiv->deriv_classes[2][8][6] = int_stack + 8160;
 Libderiv->deriv_classes[1][4][2] = int_stack + 8430;
 Libderiv->deriv_classes[1][5][2] = int_stack + 8475;
 Libderiv->deriv_classes[1][6][2] = int_stack + 8538;
 Libderiv->deriv_classes[1][7][2] = int_stack + 8622;
 Libderiv->deriv_classes[1][8][2] = int_stack + 8730;
 Libderiv->deriv_classes[2][4][2] = int_stack + 8865;
 Libderiv->deriv_classes[2][5][2] = int_stack + 8955;
 Libderiv->deriv_classes[2][6][2] = int_stack + 9081;
 Libderiv->deriv_classes[2][7][2] = int_stack + 9249;
 Libderiv->deriv_classes[2][8][2] = int_stack + 9465;
 Libderiv->deriv_classes[1][4][1] = int_stack + 9735;
 Libderiv->deriv_classes[1][5][1] = int_stack + 9780;
 Libderiv->deriv_classes[1][6][1] = int_stack + 9843;
 Libderiv->deriv_classes[1][7][1] = int_stack + 9927;
 Libderiv->deriv_classes[1][8][1] = int_stack + 10035;
 Libderiv->deriv_classes[2][4][1] = int_stack + 10170;
 Libderiv->deriv_classes[2][5][1] = int_stack + 10260;
 Libderiv->deriv_classes[2][6][1] = int_stack + 10386;
 Libderiv->deriv_classes[2][7][1] = int_stack + 10554;
 Libderiv->deriv_classes[2][8][1] = int_stack + 10770;
 Libderiv->dvrr_classes[1][4] = int_stack + 11040;
 Libderiv->dvrr_classes[1][5] = int_stack + 11085;
 Libderiv->dvrr_classes[1][6] = int_stack + 11148;
 Libderiv->dvrr_classes[1][7] = int_stack + 11232;
 Libderiv->dvrr_classes[1][8] = int_stack + 11340;
 Libderiv->deriv_classes[1][4][0] = int_stack + 11475;
 Libderiv->deriv_classes[1][5][0] = int_stack + 11520;
 Libderiv->deriv_classes[1][6][0] = int_stack + 11583;
 Libderiv->deriv_classes[1][7][0] = int_stack + 11667;
 Libderiv->deriv_classes[1][8][0] = int_stack + 11775;
 Libderiv->deriv_classes[2][4][0] = int_stack + 11910;
 Libderiv->deriv_classes[2][5][0] = int_stack + 12000;
 Libderiv->deriv_classes[2][6][0] = int_stack + 12126;
 Libderiv->deriv_classes[2][7][0] = int_stack + 12294;
 Libderiv->deriv_classes[2][8][0] = int_stack + 12510;
 memset(int_stack,0,102240);

 Libderiv->dvrr_stack = int_stack + 28971;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12780,int_stack+11085,int_stack+11040,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12915,int_stack+11148,int_stack+11085,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+13104,int_stack+12915,int_stack+12780,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+13374,int_stack+11232,int_stack+11148,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+13626,int_stack+13374,int_stack+12915,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+14004,int_stack+13626,int_stack+13104,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14454,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11040,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14589,int_stack+108,int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11085,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14778,int_stack+14589,int_stack+14454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12780,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+15048,int_stack+192,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11148,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+15300,int_stack+15048,int_stack+14589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12915,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+15678,int_stack+15300,int_stack+14778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13104,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+14454,int_stack+300,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11232,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+16128,int_stack+14454,int_stack+15048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13374,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+14454,int_stack+16128,int_stack+15300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13626,3);
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+16128,int_stack+14454,int_stack+15678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14004,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14454,int_stack+7140,int_stack+6960,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+14724,int_stack+7392,int_stack+7140,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15102,int_stack+14724,int_stack+14454,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+16803,int_stack+7728,int_stack+7392,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+17307,int_stack+16803,int_stack+14724,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+18063,int_stack+17307,int_stack+15102,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15642,int_stack+525,int_stack+435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6960,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+651,int_stack+525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7140,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18963,int_stack+0,int_stack+15642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14454,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+19503,int_stack+819,int_stack+651, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7392,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+20007,int_stack+19503,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14724,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+20763,int_stack+20007,int_stack+18963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15102,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+1035,int_stack+819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7728,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+21663,int_stack+0,int_stack+19503, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16803,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+21663,int_stack+20007, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17307,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+18963,int_stack+0,int_stack+20763, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18063,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1350,int_stack+1305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11040, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+135,int_stack+1413,int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11085, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+324,int_stack+135,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12780, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+594,int_stack+1497,int_stack+1413, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11148, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+846,int_stack+594,int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12915, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+20313,int_stack+846,int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13104, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+1605,int_stack+1497, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11232, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1224,int_stack+0,int_stack+594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13374, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+1224,int_stack+846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13626, 0.0, zero_stack,3);
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+630,int_stack+0,int_stack+20313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14004, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20313,int_stack+1830,int_stack+1740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6960, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20583,int_stack+1956,int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7140, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20961,int_stack+20583,int_stack+20313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14454, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+21501,int_stack+2124,int_stack+1956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7392, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+22005,int_stack+21501,int_stack+20583, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14724, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+22761,int_stack+22005,int_stack+20961, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15102, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+20313,int_stack+2340,int_stack+2124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7728, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1305,int_stack+20313,int_stack+21501, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16803, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+20313,int_stack+1305,int_stack+22005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17307, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+23661,int_stack+20313,int_stack+22761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18063, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20313,int_stack+2655,int_stack+2610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11040, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20448,int_stack+2718,int_stack+2655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11085, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20637,int_stack+20448,int_stack+20313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12780, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20907,int_stack+2802,int_stack+2718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11148, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21159,int_stack+20907,int_stack+20448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12915, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+21537,int_stack+21159,int_stack+20637, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13104, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+20313,int_stack+2910,int_stack+2802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11232, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+21987,int_stack+20313,int_stack+20907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13374, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+21987,int_stack+21159, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13626, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+21987,int_stack+0,int_stack+21537, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14004, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3135,int_stack+3045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6960, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22662,int_stack+3261,int_stack+3135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7140, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23040,int_stack+22662,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14454, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+3429,int_stack+3261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7392, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+20313,int_stack+0,int_stack+22662, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14724, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+21069,int_stack+20313,int_stack+23040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15102, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+22662,int_stack+3645,int_stack+3429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7728, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1305,int_stack+22662,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16803, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+2313,int_stack+1305,int_stack+20313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17307, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+25011,int_stack+2313,int_stack+21069, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18063, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20313,int_stack+3960,int_stack+3915, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20448,int_stack+4023,int_stack+3960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20637,int_stack+20448,int_stack+20313, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20907,int_stack+4107,int_stack+4023, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21159,int_stack+20907,int_stack+20448, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+21537,int_stack+21159,int_stack+20637, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+20313,int_stack+4215,int_stack+4107, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1305,int_stack+20313,int_stack+20907, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+1305,int_stack+21159, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+1305,int_stack+0,int_stack+21537, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+4440,int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1980,int_stack+4566,int_stack+4440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2358,int_stack+1980,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+4734,int_stack+4566, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2898,int_stack+0,int_stack+1980, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14724, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3654,int_stack+2898,int_stack+2358, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+1980,int_stack+4950,int_stack+4734, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+20313,int_stack+1980,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16803, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+26361,int_stack+20313,int_stack+2898, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+20313,int_stack+26361,int_stack+3654, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18063, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26361,int_stack+5265,int_stack+5220, 0.0, zero_stack, 1.0, int_stack+11040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26496,int_stack+5328,int_stack+5265, 0.0, zero_stack, 1.0, int_stack+11085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+26685,int_stack+26496,int_stack+26361, 0.0, zero_stack, 1.0, int_stack+12780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+26955,int_stack+5412,int_stack+5328, 0.0, zero_stack, 1.0, int_stack+11148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+27207,int_stack+26955,int_stack+26496, 0.0, zero_stack, 1.0, int_stack+12915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+27207,int_stack+26685, 0.0, zero_stack, 1.0, int_stack+13104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+21663,int_stack+5520,int_stack+5412, 0.0, zero_stack, 1.0, int_stack+11232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+26361,int_stack+21663,int_stack+26955, 0.0, zero_stack, 1.0, int_stack+13374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+1980,int_stack+26361,int_stack+27207, 0.0, zero_stack, 1.0, int_stack+13626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+26361,int_stack+1980,int_stack+0, 0.0, zero_stack, 1.0, int_stack+14004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+5745,int_stack+5655, 0.0, zero_stack, 1.0, int_stack+6960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1980,int_stack+5871,int_stack+5745, 0.0, zero_stack, 1.0, int_stack+7140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2358,int_stack+1980,int_stack+0, 0.0, zero_stack, 1.0, int_stack+14454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+6039,int_stack+5871, 0.0, zero_stack, 1.0, int_stack+7392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2898,int_stack+0,int_stack+1980, 0.0, zero_stack, 1.0, int_stack+14724, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3654,int_stack+2898,int_stack+2358, 0.0, zero_stack, 1.0, int_stack+15102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+1980,int_stack+6255,int_stack+6039, 0.0, zero_stack, 1.0, int_stack+7728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+4554,int_stack+1980,int_stack+0, 0.0, zero_stack, 1.0, int_stack+16803, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+27036,int_stack+4554,int_stack+2898, 0.0, zero_stack, 1.0, int_stack+17307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+4554,int_stack+27036,int_stack+3654, 0.0, zero_stack, 1.0, int_stack+18063, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+27036,int_stack+6570,int_stack+6525, 1.0, int_stack+11040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+27171,int_stack+6633,int_stack+6570, 1.0, int_stack+11085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27360,int_stack+27171,int_stack+27036, 1.0, int_stack+12780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+27630,int_stack+6717,int_stack+6633, 1.0, int_stack+11148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+27882,int_stack+27630,int_stack+27171, 1.0, int_stack+12915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+28260,int_stack+27882,int_stack+27360, 1.0, int_stack+13104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+21663,int_stack+6825,int_stack+6717, 1.0, int_stack+11232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+12780,int_stack+21663,int_stack+27630, 1.0, int_stack+13374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+12780,int_stack+27882, 1.0, int_stack+13626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+27036,int_stack+0,int_stack+28260, 1.0, int_stack+14004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+7266,int_stack+7050, 1.0, int_stack+6960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12780,int_stack+7560,int_stack+7266, 1.0, int_stack+7140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27711,int_stack+12780,int_stack+0, 1.0, int_stack+14454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+7944,int_stack+7560, 1.0, int_stack+7392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+5904,int_stack+0,int_stack+12780, 1.0, int_stack+14724, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6660,int_stack+5904,int_stack+27711, 1.0, int_stack+15102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+27711,int_stack+8160,int_stack+7944, 1.0, int_stack+7728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+14454,int_stack+27711,int_stack+0, 1.0, int_stack+16803, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+27711,int_stack+14454,int_stack+5904, 1.0, int_stack+17307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+14454,int_stack+27711,int_stack+6660, 1.0, int_stack+18063, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+15804,int_stack+11340,int_stack+11232,3);
 /*--- compute (p0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+27711,int_stack+15804,int_stack+13374,3);
 /*--- compute (p0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+27711,int_stack+13626,3);
 /*--- compute (p0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+27711,int_stack+0,int_stack+14004,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+8475,int_stack+8430,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+135,int_stack+8538,int_stack+8475,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+324,int_stack+135,int_stack+0,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+28386,int_stack+8622,int_stack+8538,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+5904,int_stack+28386,int_stack+135,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+6282,int_stack+5904,int_stack+324,3);
 /*--- compute (p0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+15804,int_stack+8730,int_stack+8622,3);
 /*--- compute (p0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+0,int_stack+15804,int_stack+28386,3);
 /*--- compute (p0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+6732,int_stack+0,int_stack+5904,3);
 /*--- compute (p0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+7362,int_stack+6732,int_stack+6282,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5904,int_stack+8955,int_stack+8865,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6174,int_stack+9081,int_stack+8955,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6552,int_stack+6174,int_stack+5904,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+9249,int_stack+9081,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+8037,int_stack+0,int_stack+6174,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+16803,int_stack+8037,int_stack+6552,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+5904,int_stack+9465,int_stack+9249,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+17703,int_stack+5904,int_stack+0,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+5904,int_stack+17703,int_stack+8037,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+8037,int_stack+5904,int_stack+16803,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16803,int_stack+9780,int_stack+9735,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16938,int_stack+9843,int_stack+9780,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17127,int_stack+16938,int_stack+16803,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+17397,int_stack+9927,int_stack+9843,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+17649,int_stack+17397,int_stack+16938,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+18027,int_stack+17649,int_stack+17127,3);
 /*--- compute (p0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+15804,int_stack+10035,int_stack+9927,3);
 /*--- compute (p0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+16803,int_stack+15804,int_stack+17397,3);
 /*--- compute (p0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+16803,int_stack+17649,3);
 /*--- compute (p0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+16803,int_stack+0,int_stack+18027,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+10260,int_stack+10170,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+17478,int_stack+10386,int_stack+10260,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17856,int_stack+17478,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+10554,int_stack+10386,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+5904,int_stack+0,int_stack+17478,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+9387,int_stack+5904,int_stack+17856,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+17478,int_stack+10770,int_stack+10554,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+10287,int_stack+17478,int_stack+0,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+17478,int_stack+10287,int_stack+5904,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+5904,int_stack+17478,int_stack+9387,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9387,int_stack+11520,int_stack+11475,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+9522,int_stack+11583,int_stack+11520,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9711,int_stack+9522,int_stack+9387,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+9981,int_stack+11667,int_stack+11583,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+10233,int_stack+9981,int_stack+9522,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+10611,int_stack+10233,int_stack+9711,3);
 /*--- compute (p0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+15804,int_stack+11775,int_stack+11667,3);
 /*--- compute (p0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+9387,int_stack+15804,int_stack+9981,3);
 /*--- compute (p0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+9387,int_stack+10233,3);
 /*--- compute (p0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+9387,int_stack+0,int_stack+10611,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+12000,int_stack+11910,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10062,int_stack+12126,int_stack+12000,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10440,int_stack+10062,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+12294,int_stack+12126,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+10980,int_stack+0,int_stack+10062,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+17478,int_stack+10980,int_stack+10440,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+10062,int_stack+12510,int_stack+12294,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+11736,int_stack+10062,int_stack+0,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+12744,int_stack+11736,int_stack+10980,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+10062,int_stack+12744,int_stack+17478,6);
 /*--- compute (pp|gg) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+11412,int_stack+18963,int_stack+16128,225);
     Libderiv->ABCD[11] = int_stack + 11412;
 /*--- compute (pp|gg) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+17478,int_stack+23661,int_stack+630,225);
     Libderiv->ABCD[10] = int_stack + 17478;
 /*--- compute (pp|gg) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1980,int_stack+25011,int_stack+21987,225);
     Libderiv->ABCD[9] = int_stack + 1980;
 /*--- compute (pp|gg) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+21663,int_stack+20313,int_stack+1305,225);
     Libderiv->ABCD[8] = int_stack + 21663;
 /*--- compute (pp|gg) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+19503,int_stack+4554,int_stack+26361,225);
     Libderiv->ABCD[7] = int_stack + 19503;
 /*--- compute (pp|gg) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+23688,int_stack+14454,int_stack+27036,225);
     Libderiv->ABCD[6] = int_stack + 23688;
 /*--- compute (pp|gg) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+13437,int_stack+8037,int_stack+7362, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 13437;
 /*--- compute (pp|gg) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7254,int_stack+5904,int_stack+16803, 0.0, zero_stack, 1.0, int_stack+27711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 7254;
 /*--- compute (pp|gg) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4005,int_stack+10062,int_stack+9387, 1.0, int_stack+27711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 4005;

}
