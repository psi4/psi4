#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|gf) integrals */

void d1hrr_order_ppgf(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[2][4][11] = int_stack + 300;
 Libderiv->deriv_classes[2][5][11] = int_stack + 390;
 Libderiv->deriv_classes[2][6][11] = int_stack + 516;
 Libderiv->deriv_classes[2][7][11] = int_stack + 684;
 Libderiv->deriv_classes[1][4][10] = int_stack + 900;
 Libderiv->deriv_classes[1][5][10] = int_stack + 945;
 Libderiv->deriv_classes[1][6][10] = int_stack + 1008;
 Libderiv->deriv_classes[1][7][10] = int_stack + 1092;
 Libderiv->deriv_classes[2][4][10] = int_stack + 1200;
 Libderiv->deriv_classes[2][5][10] = int_stack + 1290;
 Libderiv->deriv_classes[2][6][10] = int_stack + 1416;
 Libderiv->deriv_classes[2][7][10] = int_stack + 1584;
 Libderiv->deriv_classes[1][4][9] = int_stack + 1800;
 Libderiv->deriv_classes[1][5][9] = int_stack + 1845;
 Libderiv->deriv_classes[1][6][9] = int_stack + 1908;
 Libderiv->deriv_classes[1][7][9] = int_stack + 1992;
 Libderiv->deriv_classes[2][4][9] = int_stack + 2100;
 Libderiv->deriv_classes[2][5][9] = int_stack + 2190;
 Libderiv->deriv_classes[2][6][9] = int_stack + 2316;
 Libderiv->deriv_classes[2][7][9] = int_stack + 2484;
 Libderiv->deriv_classes[1][4][8] = int_stack + 2700;
 Libderiv->deriv_classes[1][5][8] = int_stack + 2745;
 Libderiv->deriv_classes[1][6][8] = int_stack + 2808;
 Libderiv->deriv_classes[1][7][8] = int_stack + 2892;
 Libderiv->deriv_classes[2][4][8] = int_stack + 3000;
 Libderiv->deriv_classes[2][5][8] = int_stack + 3090;
 Libderiv->deriv_classes[2][6][8] = int_stack + 3216;
 Libderiv->deriv_classes[2][7][8] = int_stack + 3384;
 Libderiv->deriv_classes[1][4][7] = int_stack + 3600;
 Libderiv->deriv_classes[1][5][7] = int_stack + 3645;
 Libderiv->deriv_classes[1][6][7] = int_stack + 3708;
 Libderiv->deriv_classes[1][7][7] = int_stack + 3792;
 Libderiv->deriv_classes[2][4][7] = int_stack + 3900;
 Libderiv->deriv_classes[2][5][7] = int_stack + 3990;
 Libderiv->deriv_classes[2][6][7] = int_stack + 4116;
 Libderiv->deriv_classes[2][7][7] = int_stack + 4284;
 Libderiv->deriv_classes[1][4][6] = int_stack + 4500;
 Libderiv->deriv_classes[1][5][6] = int_stack + 4545;
 Libderiv->deriv_classes[1][6][6] = int_stack + 4608;
 Libderiv->deriv_classes[1][7][6] = int_stack + 4692;
 Libderiv->dvrr_classes[2][4] = int_stack + 4800;
 Libderiv->deriv_classes[2][4][6] = int_stack + 4890;
 Libderiv->dvrr_classes[2][5] = int_stack + 4980;
 Libderiv->deriv_classes[2][5][6] = int_stack + 5106;
 Libderiv->dvrr_classes[2][6] = int_stack + 5232;
 Libderiv->deriv_classes[2][6][6] = int_stack + 5400;
 Libderiv->deriv_classes[2][7][6] = int_stack + 5568;
 Libderiv->deriv_classes[1][4][2] = int_stack + 5784;
 Libderiv->deriv_classes[1][5][2] = int_stack + 5829;
 Libderiv->deriv_classes[1][6][2] = int_stack + 5892;
 Libderiv->deriv_classes[1][7][2] = int_stack + 5976;
 Libderiv->deriv_classes[2][4][2] = int_stack + 6084;
 Libderiv->deriv_classes[2][5][2] = int_stack + 6174;
 Libderiv->deriv_classes[2][6][2] = int_stack + 6300;
 Libderiv->deriv_classes[2][7][2] = int_stack + 6468;
 Libderiv->deriv_classes[1][4][1] = int_stack + 6684;
 Libderiv->deriv_classes[1][5][1] = int_stack + 6729;
 Libderiv->deriv_classes[1][6][1] = int_stack + 6792;
 Libderiv->deriv_classes[1][7][1] = int_stack + 6876;
 Libderiv->deriv_classes[2][4][1] = int_stack + 6984;
 Libderiv->deriv_classes[2][5][1] = int_stack + 7074;
 Libderiv->deriv_classes[2][6][1] = int_stack + 7200;
 Libderiv->deriv_classes[2][7][1] = int_stack + 7368;
 Libderiv->dvrr_classes[1][4] = int_stack + 7584;
 Libderiv->dvrr_classes[1][5] = int_stack + 7629;
 Libderiv->dvrr_classes[1][6] = int_stack + 7692;
 Libderiv->dvrr_classes[1][7] = int_stack + 7776;
 Libderiv->deriv_classes[1][4][0] = int_stack + 7884;
 Libderiv->deriv_classes[1][5][0] = int_stack + 7929;
 Libderiv->deriv_classes[1][6][0] = int_stack + 7992;
 Libderiv->deriv_classes[1][7][0] = int_stack + 8076;
 Libderiv->deriv_classes[2][4][0] = int_stack + 8184;
 Libderiv->deriv_classes[2][5][0] = int_stack + 8274;
 Libderiv->deriv_classes[2][6][0] = int_stack + 8400;
 Libderiv->deriv_classes[2][7][0] = int_stack + 8568;
 memset(int_stack,0,70272);

 Libderiv->dvrr_stack = int_stack + 19602;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8784,int_stack+7629,int_stack+7584,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8919,int_stack+7692,int_stack+7629,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9108,int_stack+8919,int_stack+8784,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9378,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7584,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9513,int_stack+108,int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7629,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9702,int_stack+9513,int_stack+9378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8784,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9972,int_stack+192,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7692,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10224,int_stack+9972,int_stack+9513, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8919,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+10602,int_stack+10224,int_stack+9702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9108,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9378,int_stack+4980,int_stack+4800,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+9648,int_stack+5232,int_stack+4980,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10026,int_stack+9648,int_stack+9378,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+390,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11052,int_stack+516,int_stack+390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4980,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11430,int_stack+11052,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9378,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+0,int_stack+684,int_stack+516, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5232,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+11970,int_stack+0,int_stack+11052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9648,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+11970,int_stack+11430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10026,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11052,int_stack+945,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7584, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11187,int_stack+1008,int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7629, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11376,int_stack+11187,int_stack+11052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8784, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+11646,int_stack+1092,int_stack+1008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7692, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+11898,int_stack+11646,int_stack+11187, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8919, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+12276,int_stack+11898,int_stack+11376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9108, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11052,int_stack+1290,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11322,int_stack+1416,int_stack+1290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4980, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11700,int_stack+11322,int_stack+11052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9378, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+900,int_stack+1584,int_stack+1416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5232, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+12726,int_stack+900,int_stack+11322, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9648, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+900,int_stack+12726,int_stack+11700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10026, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12726,int_stack+1845,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7584, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12861,int_stack+1908,int_stack+1845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7629, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13050,int_stack+12861,int_stack+12726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8784, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+13320,int_stack+1992,int_stack+1908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7692, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13572,int_stack+13320,int_stack+12861, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8919, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+11052,int_stack+13572,int_stack+13050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9108, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12726,int_stack+2190,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12996,int_stack+2316,int_stack+2190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4980, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13374,int_stack+12996,int_stack+12726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9378, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1800,int_stack+2484,int_stack+2316, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5232, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+11502,int_stack+1800,int_stack+12996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9648, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1800,int_stack+11502,int_stack+13374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10026, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11502,int_stack+2745,int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11637,int_stack+2808,int_stack+2745, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7629, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11826,int_stack+11637,int_stack+11502, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+12726,int_stack+2892,int_stack+2808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+12978,int_stack+12726,int_stack+11637, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+13356,int_stack+12978,int_stack+11826, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12726,int_stack+3090,int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11502,int_stack+3216,int_stack+3090, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13806,int_stack+11502,int_stack+12726, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+12726,int_stack+3384,int_stack+3216, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2700,int_stack+12726,int_stack+11502, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+14346,int_stack+2700,int_stack+13806, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13806,int_stack+3645,int_stack+3600, 0.0, zero_stack, 1.0, int_stack+7584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13941,int_stack+3708,int_stack+3645, 0.0, zero_stack, 1.0, int_stack+7629, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2700,int_stack+13941,int_stack+13806, 0.0, zero_stack, 1.0, int_stack+8784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2970,int_stack+3792,int_stack+3708, 0.0, zero_stack, 1.0, int_stack+7692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3222,int_stack+2970,int_stack+13941, 0.0, zero_stack, 1.0, int_stack+8919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+13806,int_stack+3222,int_stack+2700, 0.0, zero_stack, 1.0, int_stack+9108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2700,int_stack+3990,int_stack+3900, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2970,int_stack+4116,int_stack+3990, 0.0, zero_stack, 1.0, int_stack+4980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3348,int_stack+2970,int_stack+2700, 0.0, zero_stack, 1.0, int_stack+9378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+11502,int_stack+4284,int_stack+4116, 0.0, zero_stack, 1.0, int_stack+5232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+15246,int_stack+11502,int_stack+2970, 0.0, zero_stack, 1.0, int_stack+9648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+16002,int_stack+15246,int_stack+3348, 0.0, zero_stack, 1.0, int_stack+10026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15246,int_stack+4545,int_stack+4500, 1.0, int_stack+7584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15381,int_stack+4608,int_stack+4545, 1.0, int_stack+7629, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15570,int_stack+15381,int_stack+15246, 1.0, int_stack+8784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+11502,int_stack+4692,int_stack+4608, 1.0, int_stack+7692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+11754,int_stack+11502,int_stack+15381, 1.0, int_stack+8919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2700,int_stack+11754,int_stack+15570, 1.0, int_stack+9108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11502,int_stack+5106,int_stack+4890, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11772,int_stack+5400,int_stack+5106, 1.0, int_stack+4980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15246,int_stack+11772,int_stack+11502, 1.0, int_stack+9378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3150,int_stack+5568,int_stack+5400, 1.0, int_stack+5232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3654,int_stack+3150,int_stack+11772, 1.0, int_stack+9648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4410,int_stack+3654,int_stack+15246, 1.0, int_stack+10026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+15246,int_stack+7776,int_stack+7692,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+15498,int_stack+15246,int_stack+8919,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+3150,int_stack+15498,int_stack+9108,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15246,int_stack+5829,int_stack+5784,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15381,int_stack+5892,int_stack+5829,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15570,int_stack+15381,int_stack+15246,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+3600,int_stack+5976,int_stack+5892,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+3852,int_stack+3600,int_stack+15381,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+11502,int_stack+3852,int_stack+15570,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3600,int_stack+6174,int_stack+6084,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3870,int_stack+6300,int_stack+6174,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15246,int_stack+3870,int_stack+3600,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8784,int_stack+6468,int_stack+6300,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+9288,int_stack+8784,int_stack+3870,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+5310,int_stack+9288,int_stack+15246,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15246,int_stack+6729,int_stack+6684,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15381,int_stack+6792,int_stack+6729,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15570,int_stack+15381,int_stack+15246,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8784,int_stack+6876,int_stack+6792,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+9036,int_stack+8784,int_stack+15381,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+9414,int_stack+9036,int_stack+15570,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8784,int_stack+7074,int_stack+6984,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15246,int_stack+7200,int_stack+7074,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9864,int_stack+15246,int_stack+8784,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8784,int_stack+7368,int_stack+7200,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+3600,int_stack+8784,int_stack+15246,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+6210,int_stack+3600,int_stack+9864,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9864,int_stack+7929,int_stack+7884,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+9999,int_stack+7992,int_stack+7929,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10188,int_stack+9999,int_stack+9864,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+3600,int_stack+8076,int_stack+7992,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+3852,int_stack+3600,int_stack+9999,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+15246,int_stack+3852,int_stack+10188,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3600,int_stack+8274,int_stack+8184,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3870,int_stack+8400,int_stack+8274,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9864,int_stack+3870,int_stack+3600,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8784,int_stack+8568,int_stack+8400,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+7110,int_stack+8784,int_stack+3870,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+7866,int_stack+7110,int_stack+9864,6);
 /*--- compute (pp|gf) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+16902,int_stack+0,int_stack+10602,150);
     Libderiv->ABCD[11] = int_stack + 16902;
 /*--- compute (pp|gf) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+18252,int_stack+900,int_stack+12276,150);
     Libderiv->ABCD[10] = int_stack + 18252;
 /*--- compute (pp|gf) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+1800,int_stack+11052,150);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (pp|gf) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1350,int_stack+14346,int_stack+13356,150);
     Libderiv->ABCD[8] = int_stack + 1350;
 /*--- compute (pp|gf) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+9864,int_stack+16002,int_stack+13806,150);
     Libderiv->ABCD[7] = int_stack + 9864;
 /*--- compute (pp|gf) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+11952,int_stack+4410,int_stack+2700,150);
     Libderiv->ABCD[6] = int_stack + 11952;
 /*--- compute (pp|gf) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3600,int_stack+5310,int_stack+11502, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 3600;
 /*--- compute (pp|gf) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+13302,int_stack+6210,int_stack+9414, 0.0, zero_stack, 1.0, int_stack+3150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 13302;
 /*--- compute (pp|gf) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+4950,int_stack+7866,int_stack+15246, 1.0, int_stack+3150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 4950;

}
