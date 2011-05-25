#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_g0gg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (g0|gg) integrals */

void d1hrr_order_g0gg(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 540;
 Libderiv->deriv_classes[4][7][11] = int_stack + 960;
 Libderiv->deriv_classes[4][8][11] = int_stack + 1500;
 Libderiv->deriv_classes[4][4][10] = int_stack + 2175;
 Libderiv->deriv_classes[4][5][10] = int_stack + 2400;
 Libderiv->deriv_classes[4][6][10] = int_stack + 2715;
 Libderiv->deriv_classes[4][7][10] = int_stack + 3135;
 Libderiv->deriv_classes[4][8][10] = int_stack + 3675;
 Libderiv->deriv_classes[4][4][9] = int_stack + 4350;
 Libderiv->deriv_classes[4][5][9] = int_stack + 4575;
 Libderiv->deriv_classes[4][6][9] = int_stack + 4890;
 Libderiv->deriv_classes[4][7][9] = int_stack + 5310;
 Libderiv->deriv_classes[4][8][9] = int_stack + 5850;
 Libderiv->deriv_classes[4][4][8] = int_stack + 6525;
 Libderiv->deriv_classes[4][5][8] = int_stack + 6750;
 Libderiv->deriv_classes[4][6][8] = int_stack + 7065;
 Libderiv->deriv_classes[4][7][8] = int_stack + 7485;
 Libderiv->deriv_classes[4][8][8] = int_stack + 8025;
 Libderiv->deriv_classes[4][4][7] = int_stack + 8700;
 Libderiv->deriv_classes[4][5][7] = int_stack + 8925;
 Libderiv->deriv_classes[4][6][7] = int_stack + 9240;
 Libderiv->deriv_classes[4][7][7] = int_stack + 9660;
 Libderiv->deriv_classes[4][8][7] = int_stack + 10200;
 Libderiv->dvrr_classes[4][4] = int_stack + 10875;
 Libderiv->deriv_classes[4][4][6] = int_stack + 11100;
 Libderiv->dvrr_classes[4][5] = int_stack + 11325;
 Libderiv->deriv_classes[4][5][6] = int_stack + 11640;
 Libderiv->dvrr_classes[4][6] = int_stack + 11955;
 Libderiv->deriv_classes[4][6][6] = int_stack + 12375;
 Libderiv->dvrr_classes[4][7] = int_stack + 12795;
 Libderiv->deriv_classes[4][7][6] = int_stack + 13335;
 Libderiv->deriv_classes[4][8][6] = int_stack + 13875;
 Libderiv->deriv_classes[4][4][2] = int_stack + 14550;
 Libderiv->deriv_classes[4][5][2] = int_stack + 14775;
 Libderiv->deriv_classes[4][6][2] = int_stack + 15090;
 Libderiv->deriv_classes[4][7][2] = int_stack + 15510;
 Libderiv->deriv_classes[4][8][2] = int_stack + 16050;
 Libderiv->deriv_classes[4][4][1] = int_stack + 16725;
 Libderiv->deriv_classes[4][5][1] = int_stack + 16950;
 Libderiv->deriv_classes[4][6][1] = int_stack + 17265;
 Libderiv->deriv_classes[4][7][1] = int_stack + 17685;
 Libderiv->deriv_classes[4][8][1] = int_stack + 18225;
 Libderiv->deriv_classes[4][4][0] = int_stack + 18900;
 Libderiv->deriv_classes[4][5][0] = int_stack + 19125;
 Libderiv->deriv_classes[4][6][0] = int_stack + 19440;
 Libderiv->deriv_classes[4][7][0] = int_stack + 19860;
 Libderiv->deriv_classes[4][8][0] = int_stack + 20400;
 memset(int_stack,0,168600);

 Libderiv->dvrr_stack = int_stack + 66165;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_g0gg(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+21075,int_stack+11325,int_stack+10875,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+21750,int_stack+11955,int_stack+11325,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+22695,int_stack+21750,int_stack+21075,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+24045,int_stack+12795,int_stack+11955,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+25305,int_stack+24045,int_stack+21750,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+27195,int_stack+25305,int_stack+22695,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29445,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10875,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30120,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11325,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31065,int_stack+30120,int_stack+29445, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21075,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+32415,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11955,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+33675,int_stack+32415,int_stack+30120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21750,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+35565,int_stack+33675,int_stack+31065, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22695,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+29445,int_stack+1500,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12795,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+37815,int_stack+29445,int_stack+32415, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24045,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+29445,int_stack+37815,int_stack+33675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25305,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37815,int_stack+2400,int_stack+2175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10875, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38490,int_stack+2715,int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11325, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+39435,int_stack+38490,int_stack+37815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21075, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+32595,int_stack+3135,int_stack+2715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11955, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+0,int_stack+32595,int_stack+38490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21750, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+40785,int_stack+0,int_stack+39435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22695, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+37815,int_stack+3675,int_stack+3135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12795, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+43035,int_stack+37815,int_stack+32595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24045, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+45555,int_stack+43035,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25305, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+4575,int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10875, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+675,int_stack+4890,int_stack+4575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11325, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1620,int_stack+675,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21075, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2970,int_stack+5310,int_stack+4890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11955, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+43035,int_stack+2970,int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21750, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+32595,int_stack+43035,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22695, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+5850,int_stack+5310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12795, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+37815,int_stack+0,int_stack+2970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24045, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+37815,int_stack+43035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25305, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43035,int_stack+6750,int_stack+6525, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43710,int_stack+7065,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+37815,int_stack+43710,int_stack+43035, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+39165,int_stack+7485,int_stack+7065, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3150,int_stack+39165,int_stack+43710, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+43035,int_stack+3150,int_stack+37815, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+5040,int_stack+8025,int_stack+7485, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+48705,int_stack+5040,int_stack+39165, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+5040,int_stack+48705,int_stack+3150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3150,int_stack+8925,int_stack+8700, 0.0, zero_stack, 1.0, int_stack+10875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3825,int_stack+9240,int_stack+8925, 0.0, zero_stack, 1.0, int_stack+11325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48705,int_stack+3825,int_stack+3150, 0.0, zero_stack, 1.0, int_stack+21075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+50055,int_stack+9660,int_stack+9240, 0.0, zero_stack, 1.0, int_stack+11955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+37815,int_stack+50055,int_stack+3825, 0.0, zero_stack, 1.0, int_stack+21750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+51315,int_stack+37815,int_stack+48705, 0.0, zero_stack, 1.0, int_stack+22695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+3150,int_stack+10200,int_stack+9660, 0.0, zero_stack, 1.0, int_stack+12795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+8190,int_stack+3150,int_stack+50055, 0.0, zero_stack, 1.0, int_stack+24045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+53565,int_stack+8190,int_stack+37815, 0.0, zero_stack, 1.0, int_stack+25305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37815,int_stack+11640,int_stack+11100, 1.0, int_stack+10875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38490,int_stack+12375,int_stack+11640, 1.0, int_stack+11325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+39435,int_stack+38490,int_stack+37815, 1.0, int_stack+21075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8190,int_stack+13335,int_stack+12375, 1.0, int_stack+11955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3150,int_stack+8190,int_stack+38490, 1.0, int_stack+21750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+9450,int_stack+3150,int_stack+39435, 1.0, int_stack+22695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+21075,int_stack+13875,int_stack+13335, 1.0, int_stack+12795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+37815,int_stack+21075,int_stack+8190, 1.0, int_stack+24045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+21075,int_stack+37815,int_stack+3150, 1.0, int_stack+25305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3150,int_stack+14775,int_stack+14550,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3825,int_stack+15090,int_stack+14775,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+37815,int_stack+3825,int_stack+3150,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8190,int_stack+15510,int_stack+15090,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+24225,int_stack+8190,int_stack+3825,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+11700,int_stack+24225,int_stack+37815,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+37815,int_stack+16050,int_stack+15510,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+13950,int_stack+37815,int_stack+8190,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+56715,int_stack+13950,int_stack+24225,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24225,int_stack+16950,int_stack+16725,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24900,int_stack+17265,int_stack+16950,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+25845,int_stack+24900,int_stack+24225,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8190,int_stack+17685,int_stack+17265,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+3150,int_stack+8190,int_stack+24900,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+13950,int_stack+3150,int_stack+25845,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+24225,int_stack+18225,int_stack+17685,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+16200,int_stack+24225,int_stack+8190,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+59865,int_stack+16200,int_stack+3150,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3150,int_stack+19125,int_stack+18900,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3825,int_stack+19440,int_stack+19125,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16200,int_stack+3825,int_stack+3150,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+8190,int_stack+19860,int_stack+19440,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+17550,int_stack+8190,int_stack+3825,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+24225,int_stack+17550,int_stack+16200,15);
 /*--- compute (g0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+3150,int_stack+20400,int_stack+19860,15);
 /*--- compute (g0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+37815,int_stack+3150,int_stack+8190,15);
 /*--- compute (g0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+63015,int_stack+37815,int_stack+17550,15);
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+16200,int_stack+29445,int_stack+35565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27195,15);
     Libderiv->ABCD[11] = int_stack + 16200;
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+34845,int_stack+45555,int_stack+40785, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27195, 0.0, zero_stack,15);
     Libderiv->ABCD[10] = int_stack + 34845;
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+45285,int_stack+0,int_stack+32595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27195, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[9] = int_stack + 45285;
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+0,int_stack+5040,int_stack+43035, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+3375,int_stack+53565,int_stack+51315, 0.0, zero_stack, 1.0, int_stack+27195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[7] = int_stack + 3375;
 /*--- compute (g0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+29445,int_stack+21075,int_stack+9450, 1.0, int_stack+27195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[6] = int_stack + 29445;
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+6750,int_stack+56715,int_stack+11700,15);
     Libderiv->ABCD[2] = int_stack + 6750;
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+10125,int_stack+59865,int_stack+13950,15);
     Libderiv->ABCD[1] = int_stack + 10125;
 /*--- compute (g0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+19575,int_stack+63015,int_stack+24225,15);
     Libderiv->ABCD[0] = int_stack + 19575;

}
