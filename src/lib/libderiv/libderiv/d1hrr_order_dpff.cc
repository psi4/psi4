#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|ff) integrals */

void d1hrr_order_dpff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[2][5][11] = int_stack + 150;
 Libderiv->deriv_classes[2][6][11] = int_stack + 276;
 Libderiv->deriv_classes[3][3][11] = int_stack + 444;
 Libderiv->deriv_classes[3][4][11] = int_stack + 544;
 Libderiv->deriv_classes[3][5][11] = int_stack + 694;
 Libderiv->deriv_classes[3][6][11] = int_stack + 904;
 Libderiv->deriv_classes[2][3][10] = int_stack + 1184;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1244;
 Libderiv->deriv_classes[2][5][10] = int_stack + 1334;
 Libderiv->deriv_classes[2][6][10] = int_stack + 1460;
 Libderiv->deriv_classes[3][3][10] = int_stack + 1628;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1728;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1878;
 Libderiv->deriv_classes[3][6][10] = int_stack + 2088;
 Libderiv->deriv_classes[2][3][9] = int_stack + 2368;
 Libderiv->deriv_classes[2][4][9] = int_stack + 2428;
 Libderiv->deriv_classes[2][5][9] = int_stack + 2518;
 Libderiv->deriv_classes[2][6][9] = int_stack + 2644;
 Libderiv->deriv_classes[3][3][9] = int_stack + 2812;
 Libderiv->deriv_classes[3][4][9] = int_stack + 2912;
 Libderiv->deriv_classes[3][5][9] = int_stack + 3062;
 Libderiv->deriv_classes[3][6][9] = int_stack + 3272;
 Libderiv->deriv_classes[2][3][8] = int_stack + 3552;
 Libderiv->deriv_classes[2][4][8] = int_stack + 3612;
 Libderiv->deriv_classes[2][5][8] = int_stack + 3702;
 Libderiv->deriv_classes[2][6][8] = int_stack + 3828;
 Libderiv->deriv_classes[3][3][8] = int_stack + 3996;
 Libderiv->deriv_classes[3][4][8] = int_stack + 4096;
 Libderiv->deriv_classes[3][5][8] = int_stack + 4246;
 Libderiv->deriv_classes[3][6][8] = int_stack + 4456;
 Libderiv->deriv_classes[2][3][7] = int_stack + 4736;
 Libderiv->deriv_classes[2][4][7] = int_stack + 4796;
 Libderiv->deriv_classes[2][5][7] = int_stack + 4886;
 Libderiv->deriv_classes[2][6][7] = int_stack + 5012;
 Libderiv->deriv_classes[3][3][7] = int_stack + 5180;
 Libderiv->deriv_classes[3][4][7] = int_stack + 5280;
 Libderiv->deriv_classes[3][5][7] = int_stack + 5430;
 Libderiv->deriv_classes[3][6][7] = int_stack + 5640;
 Libderiv->deriv_classes[2][3][6] = int_stack + 5920;
 Libderiv->deriv_classes[2][4][6] = int_stack + 5980;
 Libderiv->deriv_classes[2][5][6] = int_stack + 6070;
 Libderiv->deriv_classes[2][6][6] = int_stack + 6196;
 Libderiv->dvrr_classes[3][3] = int_stack + 6364;
 Libderiv->deriv_classes[3][3][6] = int_stack + 6464;
 Libderiv->dvrr_classes[3][4] = int_stack + 6564;
 Libderiv->deriv_classes[3][4][6] = int_stack + 6714;
 Libderiv->dvrr_classes[3][5] = int_stack + 6864;
 Libderiv->deriv_classes[3][5][6] = int_stack + 7074;
 Libderiv->deriv_classes[3][6][6] = int_stack + 7284;
 Libderiv->deriv_classes[2][3][2] = int_stack + 7564;
 Libderiv->deriv_classes[2][4][2] = int_stack + 7624;
 Libderiv->deriv_classes[2][5][2] = int_stack + 7714;
 Libderiv->deriv_classes[2][6][2] = int_stack + 7840;
 Libderiv->deriv_classes[3][3][2] = int_stack + 8008;
 Libderiv->deriv_classes[3][4][2] = int_stack + 8108;
 Libderiv->deriv_classes[3][5][2] = int_stack + 8258;
 Libderiv->deriv_classes[3][6][2] = int_stack + 8468;
 Libderiv->deriv_classes[2][3][1] = int_stack + 8748;
 Libderiv->deriv_classes[2][4][1] = int_stack + 8808;
 Libderiv->deriv_classes[2][5][1] = int_stack + 8898;
 Libderiv->deriv_classes[2][6][1] = int_stack + 9024;
 Libderiv->deriv_classes[3][3][1] = int_stack + 9192;
 Libderiv->deriv_classes[3][4][1] = int_stack + 9292;
 Libderiv->deriv_classes[3][5][1] = int_stack + 9442;
 Libderiv->deriv_classes[3][6][1] = int_stack + 9652;
 Libderiv->dvrr_classes[2][3] = int_stack + 9932;
 Libderiv->dvrr_classes[2][4] = int_stack + 9992;
 Libderiv->dvrr_classes[2][5] = int_stack + 10082;
 Libderiv->dvrr_classes[2][6] = int_stack + 10208;
 Libderiv->deriv_classes[2][3][0] = int_stack + 10376;
 Libderiv->deriv_classes[2][4][0] = int_stack + 10436;
 Libderiv->deriv_classes[2][5][0] = int_stack + 10526;
 Libderiv->deriv_classes[2][6][0] = int_stack + 10652;
 Libderiv->deriv_classes[3][3][0] = int_stack + 10820;
 Libderiv->deriv_classes[3][4][0] = int_stack + 10920;
 Libderiv->deriv_classes[3][5][0] = int_stack + 11070;
 Libderiv->deriv_classes[3][6][0] = int_stack + 11280;
 memset(int_stack,0,92480);

 Libderiv->dvrr_stack = int_stack + 22548;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpff(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11560,int_stack+9992,int_stack+9932,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11740,int_stack+10082,int_stack+9992,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12010,int_stack+11740,int_stack+11560,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12370,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9932,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12550,int_stack+150,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9992,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12820,int_stack+12550,int_stack+12370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11560,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13180,int_stack+276,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10082,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13558,int_stack+13180,int_stack+12550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11740,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14098,int_stack+13558,int_stack+12820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12010,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12370,int_stack+6564,int_stack+6364,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12670,int_stack+6864,int_stack+6564,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+13120,int_stack+12670,int_stack+12370,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13720,int_stack+544,int_stack+444, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6364,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+694,int_stack+544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6564,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14698,int_stack+0,int_stack+13720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12370,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15298,int_stack+904,int_stack+694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6864,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15928,int_stack+15298,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12670,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+15928,int_stack+14698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13120,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14698,int_stack+1244,int_stack+1184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9932, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14878,int_stack+1334,int_stack+1244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9992, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15148,int_stack+14878,int_stack+14698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11560, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13720,int_stack+1460,int_stack+1334, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10082, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15508,int_stack+13720,int_stack+14878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11740, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16048,int_stack+15508,int_stack+15148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12010, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13720,int_stack+1728,int_stack+1628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6364, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14698,int_stack+1878,int_stack+1728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6564, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15148,int_stack+14698,int_stack+13720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12370, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1000,int_stack+2088,int_stack+1878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6864, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16648,int_stack+1000,int_stack+14698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12670, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1000,int_stack+16648,int_stack+15148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13120, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16648,int_stack+2428,int_stack+2368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9932, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16828,int_stack+2518,int_stack+2428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9992, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17098,int_stack+16828,int_stack+16648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11560, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13720,int_stack+2644,int_stack+2518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10082, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17458,int_stack+13720,int_stack+16828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11740, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14698,int_stack+17458,int_stack+17098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12010, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13720,int_stack+2912,int_stack+2812, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6364, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16648,int_stack+3062,int_stack+2912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6564, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17098,int_stack+16648,int_stack+13720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12370, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15298,int_stack+3272,int_stack+3062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6864, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2000,int_stack+15298,int_stack+16648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12670, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+17698,int_stack+2000,int_stack+17098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13120, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2000,int_stack+3612,int_stack+3552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2180,int_stack+3702,int_stack+3612, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2450,int_stack+2180,int_stack+2000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13720,int_stack+3828,int_stack+3702, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2810,int_stack+13720,int_stack+2180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3350,int_stack+2810,int_stack+2450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13720,int_stack+4096,int_stack+3996, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2000,int_stack+4246,int_stack+4096, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2450,int_stack+2000,int_stack+13720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16648,int_stack+4456,int_stack+4246, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18698,int_stack+16648,int_stack+2000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16648,int_stack+18698,int_stack+2450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18698,int_stack+4796,int_stack+4736, 0.0, zero_stack, 1.0, int_stack+9932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18878,int_stack+4886,int_stack+4796, 0.0, zero_stack, 1.0, int_stack+9992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19148,int_stack+18878,int_stack+18698, 0.0, zero_stack, 1.0, int_stack+11560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13720,int_stack+5012,int_stack+4886, 0.0, zero_stack, 1.0, int_stack+10082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2000,int_stack+13720,int_stack+18878, 0.0, zero_stack, 1.0, int_stack+11740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2540,int_stack+2000,int_stack+19148, 0.0, zero_stack, 1.0, int_stack+12010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2000,int_stack+5280,int_stack+5180, 0.0, zero_stack, 1.0, int_stack+6364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18698,int_stack+5430,int_stack+5280, 0.0, zero_stack, 1.0, int_stack+6564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19148,int_stack+18698,int_stack+2000, 0.0, zero_stack, 1.0, int_stack+12370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3950,int_stack+5640,int_stack+5430, 0.0, zero_stack, 1.0, int_stack+6864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4580,int_stack+3950,int_stack+18698, 0.0, zero_stack, 1.0, int_stack+12670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+19748,int_stack+4580,int_stack+19148, 0.0, zero_stack, 1.0, int_stack+13120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18698,int_stack+5980,int_stack+5920, 1.0, int_stack+9932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18878,int_stack+6070,int_stack+5980, 1.0, int_stack+9992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19148,int_stack+18878,int_stack+18698, 1.0, int_stack+11560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13720,int_stack+6196,int_stack+6070, 1.0, int_stack+10082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2000,int_stack+13720,int_stack+18878, 1.0, int_stack+11740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3950,int_stack+2000,int_stack+19148, 1.0, int_stack+12010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2000,int_stack+6714,int_stack+6464, 1.0, int_stack+6364, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18698,int_stack+7074,int_stack+6714, 1.0, int_stack+6564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19148,int_stack+18698,int_stack+2000, 1.0, int_stack+12370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4550,int_stack+7284,int_stack+7074, 1.0, int_stack+6864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5180,int_stack+4550,int_stack+18698, 1.0, int_stack+12670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6080,int_stack+5180,int_stack+19148, 1.0, int_stack+13120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+18698,int_stack+10208,int_stack+10082,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2000,int_stack+18698,int_stack+11740,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+18698,int_stack+2000,int_stack+12010,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2000,int_stack+7624,int_stack+7564,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2180,int_stack+7714,int_stack+7624,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+19298,int_stack+2180,int_stack+2000,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4550,int_stack+7840,int_stack+7714,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4928,int_stack+4550,int_stack+2180,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+5468,int_stack+4928,int_stack+19298,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19298,int_stack+8108,int_stack+8008,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4550,int_stack+8258,int_stack+8108,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11560,int_stack+4550,int_stack+19298,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12160,int_stack+8468,int_stack+8258,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12790,int_stack+12160,int_stack+4550,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+7080,int_stack+12790,int_stack+11560,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11560,int_stack+8808,int_stack+8748,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11740,int_stack+8898,int_stack+8808,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12010,int_stack+11740,int_stack+11560,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12370,int_stack+9024,int_stack+8898,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2000,int_stack+12370,int_stack+11740,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+12370,int_stack+2000,int_stack+12010,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2000,int_stack+9292,int_stack+9192,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19298,int_stack+9442,int_stack+9292,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12970,int_stack+19298,int_stack+2000,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11560,int_stack+9652,int_stack+9442,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4550,int_stack+11560,int_stack+19298,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+8080,int_stack+4550,int_stack+12970,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12970,int_stack+10436,int_stack+10376,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13150,int_stack+10526,int_stack+10436,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+13420,int_stack+13150,int_stack+12970,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4550,int_stack+10652,int_stack+10526,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4928,int_stack+4550,int_stack+13150,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+11560,int_stack+4928,int_stack+13420,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4550,int_stack+10920,int_stack+10820,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19298,int_stack+11070,int_stack+10920,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4850,int_stack+19298,int_stack+4550,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12970,int_stack+11280,int_stack+11070,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9080,int_stack+12970,int_stack+19298,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+12970,int_stack+9080,int_stack+4850,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9080,int_stack+0,int_stack+14098,100);
     Libderiv->ABCD[11] = int_stack + 9080;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+20748,int_stack+1000,int_stack+16048,100);
     Libderiv->ABCD[10] = int_stack + 20748;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+17698,int_stack+14698,100);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+13970,int_stack+16648,int_stack+3350,100);
     Libderiv->ABCD[8] = int_stack + 13970;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15770,int_stack+19748,int_stack+2540,100);
     Libderiv->ABCD[7] = int_stack + 15770;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1800,int_stack+6080,int_stack+3950,100);
     Libderiv->ABCD[6] = int_stack + 1800;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+3600,int_stack+7080,int_stack+5468, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 3600;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5400,int_stack+8080,int_stack+12370, 0.0, zero_stack, 1.0, int_stack+18698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 5400;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7200,int_stack+12970,int_stack+11560, 1.0, int_stack+18698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 7200;

}
