#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00gg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|gg) integrals */

void d1hrr_order_00gg(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][4][11] = int_stack + 0;
 Libderiv->deriv_classes[0][5][11] = int_stack + 15;
 Libderiv->deriv_classes[0][6][11] = int_stack + 36;
 Libderiv->deriv_classes[0][7][11] = int_stack + 64;
 Libderiv->deriv_classes[0][8][11] = int_stack + 100;
 Libderiv->deriv_classes[0][4][10] = int_stack + 145;
 Libderiv->deriv_classes[0][5][10] = int_stack + 160;
 Libderiv->deriv_classes[0][6][10] = int_stack + 181;
 Libderiv->deriv_classes[0][7][10] = int_stack + 209;
 Libderiv->deriv_classes[0][8][10] = int_stack + 245;
 Libderiv->deriv_classes[0][4][9] = int_stack + 290;
 Libderiv->deriv_classes[0][5][9] = int_stack + 305;
 Libderiv->deriv_classes[0][6][9] = int_stack + 326;
 Libderiv->deriv_classes[0][7][9] = int_stack + 354;
 Libderiv->deriv_classes[0][8][9] = int_stack + 390;
 Libderiv->deriv_classes[0][4][8] = int_stack + 435;
 Libderiv->deriv_classes[0][5][8] = int_stack + 450;
 Libderiv->deriv_classes[0][6][8] = int_stack + 471;
 Libderiv->deriv_classes[0][7][8] = int_stack + 499;
 Libderiv->deriv_classes[0][8][8] = int_stack + 535;
 Libderiv->deriv_classes[0][4][7] = int_stack + 580;
 Libderiv->deriv_classes[0][5][7] = int_stack + 595;
 Libderiv->deriv_classes[0][6][7] = int_stack + 616;
 Libderiv->deriv_classes[0][7][7] = int_stack + 644;
 Libderiv->deriv_classes[0][8][7] = int_stack + 680;
 Libderiv->dvrr_classes[0][4] = int_stack + 725;
 Libderiv->deriv_classes[0][4][6] = int_stack + 740;
 Libderiv->dvrr_classes[0][5] = int_stack + 755;
 Libderiv->deriv_classes[0][5][6] = int_stack + 776;
 Libderiv->dvrr_classes[0][6] = int_stack + 797;
 Libderiv->deriv_classes[0][6][6] = int_stack + 825;
 Libderiv->dvrr_classes[0][7] = int_stack + 853;
 Libderiv->deriv_classes[0][7][6] = int_stack + 889;
 Libderiv->deriv_classes[0][8][6] = int_stack + 925;
 Libderiv->deriv_classes[0][4][2] = int_stack + 970;
 Libderiv->deriv_classes[0][5][2] = int_stack + 985;
 Libderiv->deriv_classes[0][6][2] = int_stack + 1006;
 Libderiv->deriv_classes[0][7][2] = int_stack + 1034;
 Libderiv->deriv_classes[0][8][2] = int_stack + 1070;
 Libderiv->deriv_classes[0][4][1] = int_stack + 1115;
 Libderiv->deriv_classes[0][5][1] = int_stack + 1130;
 Libderiv->deriv_classes[0][6][1] = int_stack + 1151;
 Libderiv->deriv_classes[0][7][1] = int_stack + 1179;
 Libderiv->deriv_classes[0][8][1] = int_stack + 1215;
 Libderiv->deriv_classes[0][4][0] = int_stack + 1260;
 Libderiv->deriv_classes[0][5][0] = int_stack + 1275;
 Libderiv->deriv_classes[0][6][0] = int_stack + 1296;
 Libderiv->deriv_classes[0][7][0] = int_stack + 1324;
 Libderiv->deriv_classes[0][8][0] = int_stack + 1360;
 memset(int_stack,0,11240);

 Libderiv->dvrr_stack = int_stack + 4648;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00gg(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1405,int_stack+755,int_stack+725,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1450,int_stack+797,int_stack+755,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1513,int_stack+1450,int_stack+1405,1);
 /*--- compute (00|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+1603,int_stack+853,int_stack+797,1);
 /*--- compute (00|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+1687,int_stack+1603,int_stack+1450,1);
 /*--- compute (00|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+1813,int_stack+1687,int_stack+1513,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1963,int_stack+15,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+725,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2008,int_stack+36,int_stack+15, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+755,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2071,int_stack+2008,int_stack+1963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1405,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2161,int_stack+64,int_stack+36, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+797,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2245,int_stack+2161,int_stack+2008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1450,1);
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2371,int_stack+2245,int_stack+2071, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1513,1);
 /*--- compute (00|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+1963,int_stack+100,int_stack+64, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+853,1);
 /*--- compute (00|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+2521,int_stack+1963,int_stack+2161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1603,1);
 /*--- compute (00|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+1963,int_stack+2521,int_stack+2245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1687,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2521,int_stack+160,int_stack+145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+725, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2566,int_stack+181,int_stack+160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+755, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2629,int_stack+2566,int_stack+2521, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1405, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2719,int_stack+209,int_stack+181, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+797, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2803,int_stack+2719,int_stack+2566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1450, 0.0, zero_stack,1);
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2173,int_stack+2803,int_stack+2629, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1513, 0.0, zero_stack,1);
 /*--- compute (00|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+2521,int_stack+245,int_stack+209, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+853, 0.0, zero_stack,1);
 /*--- compute (00|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+0,int_stack+2521,int_stack+2719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1603, 0.0, zero_stack,1);
 /*--- compute (00|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+2521,int_stack+0,int_stack+2803, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1687, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+305,int_stack+290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+725, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45,int_stack+326,int_stack+305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+755, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+108,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1405, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+198,int_stack+354,int_stack+326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+797, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2731,int_stack+198,int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1450, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2857,int_stack+2731,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1513, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+390,int_stack+354, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+853, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3007,int_stack+0,int_stack+198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1603, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+3007,int_stack+2731, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1687, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2731,int_stack+450,int_stack+435, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2776,int_stack+471,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3007,int_stack+2776,int_stack+2731, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3097,int_stack+499,int_stack+471, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+797, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3181,int_stack+3097,int_stack+2776, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3307,int_stack+3181,int_stack+3007, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1513, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+2731,int_stack+535,int_stack+499, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+853, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3457,int_stack+2731,int_stack+3097, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+3625,int_stack+3457,int_stack+3181, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3457,int_stack+595,int_stack+580, 0.0, zero_stack, 1.0, int_stack+725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3502,int_stack+616,int_stack+595, 0.0, zero_stack, 1.0, int_stack+755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2731,int_stack+3502,int_stack+3457, 0.0, zero_stack, 1.0, int_stack+1405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3007,int_stack+644,int_stack+616, 0.0, zero_stack, 1.0, int_stack+797, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3091,int_stack+3007,int_stack+3502, 0.0, zero_stack, 1.0, int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3457,int_stack+3091,int_stack+2731, 0.0, zero_stack, 1.0, int_stack+1513, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+2731,int_stack+680,int_stack+644, 0.0, zero_stack, 1.0, int_stack+853, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+210,int_stack+2731,int_stack+3007, 0.0, zero_stack, 1.0, int_stack+1603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+378,int_stack+210,int_stack+3091, 0.0, zero_stack, 1.0, int_stack+1687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+210,int_stack+776,int_stack+740, 1.0, int_stack+725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+255,int_stack+825,int_stack+776, 1.0, int_stack+755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3007,int_stack+255,int_stack+210, 1.0, int_stack+1405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3097,int_stack+889,int_stack+825, 1.0, int_stack+797, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3181,int_stack+3097,int_stack+255, 1.0, int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+210,int_stack+3181,int_stack+3007, 1.0, int_stack+1513, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+1405,int_stack+925,int_stack+889, 1.0, int_stack+853, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+588,int_stack+1405,int_stack+3097, 1.0, int_stack+1603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+1405,int_stack+588,int_stack+3181, 1.0, int_stack+1687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+588,int_stack+985,int_stack+970,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+633,int_stack+1006,int_stack+985,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+696,int_stack+633,int_stack+588,1);
 /*--- compute (00|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+786,int_stack+1034,int_stack+1006,1);
 /*--- compute (00|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+2731,int_stack+786,int_stack+633,1);
 /*--- compute (00|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+870,int_stack+2731,int_stack+696,1);
 /*--- compute (00|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+588,int_stack+1070,int_stack+1034,1);
 /*--- compute (00|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+1615,int_stack+588,int_stack+786,1);
 /*--- compute (00|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+588,int_stack+1615,int_stack+2731,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2731,int_stack+1130,int_stack+1115,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2776,int_stack+1151,int_stack+1130,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1615,int_stack+2776,int_stack+2731,1);
 /*--- compute (00|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+1705,int_stack+1179,int_stack+1151,1);
 /*--- compute (00|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+1020,int_stack+1705,int_stack+2776,1);
 /*--- compute (00|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+3007,int_stack+1020,int_stack+1615,1);
 /*--- compute (00|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+2731,int_stack+1215,int_stack+1179,1);
 /*--- compute (00|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+3835,int_stack+2731,int_stack+1705,1);
 /*--- compute (00|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+4003,int_stack+3835,int_stack+1020,1);
 /*--- compute (00|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1020,int_stack+1275,int_stack+1260,1);
 /*--- compute (00|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1065,int_stack+1296,int_stack+1275,1);
 /*--- compute (00|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1128,int_stack+1065,int_stack+1020,1);
 /*--- compute (00|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+3835,int_stack+1324,int_stack+1296,1);
 /*--- compute (00|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+2731,int_stack+3835,int_stack+1065,1);
 /*--- compute (00|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+3157,int_stack+2731,int_stack+1128,1);
 /*--- compute (00|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+1020,int_stack+1360,int_stack+1324,1);
 /*--- compute (00|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+1128,int_stack+1020,int_stack+3835,1);
 /*--- compute (00|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+4213,int_stack+1128,int_stack+2731,1);
 /*--- compute (00|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+1020,int_stack+1963,int_stack+2371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1813,1);
     Libderiv->ABCD[11] = int_stack + 1020;
 /*--- compute (00|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+4423,int_stack+2521,int_stack+2173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1813, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 4423;
 /*--- compute (00|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+1963,int_stack+0,int_stack+2857, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1813, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 1963;
 /*--- compute (00|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+2188,int_stack+3625,int_stack+3307, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1813, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 2188;
 /*--- compute (00|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+2413,int_stack+378,int_stack+3457, 0.0, zero_stack, 1.0, int_stack+1813, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 2413;
 /*--- compute (00|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+3307,int_stack+1405,int_stack+210, 1.0, int_stack+1813, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 3307;
 /*--- compute (00|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+3532,int_stack+588,int_stack+870,1);
     Libderiv->ABCD[2] = int_stack + 3532;
 /*--- compute (00|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+3757,int_stack+4003,int_stack+3007,1);
     Libderiv->ABCD[1] = int_stack + 3757;
 /*--- compute (00|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+3982,int_stack+4213,int_stack+3157,1);
     Libderiv->ABCD[0] = int_stack + 3982;

}
