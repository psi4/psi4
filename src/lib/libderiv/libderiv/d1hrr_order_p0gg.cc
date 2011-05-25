#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0gg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|gg) integrals */

void d1hrr_order_p0gg(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[1][4][10] = int_stack + 435;
 Libderiv->deriv_classes[1][5][10] = int_stack + 480;
 Libderiv->deriv_classes[1][6][10] = int_stack + 543;
 Libderiv->deriv_classes[1][7][10] = int_stack + 627;
 Libderiv->deriv_classes[1][8][10] = int_stack + 735;
 Libderiv->deriv_classes[1][4][9] = int_stack + 870;
 Libderiv->deriv_classes[1][5][9] = int_stack + 915;
 Libderiv->deriv_classes[1][6][9] = int_stack + 978;
 Libderiv->deriv_classes[1][7][9] = int_stack + 1062;
 Libderiv->deriv_classes[1][8][9] = int_stack + 1170;
 Libderiv->deriv_classes[1][4][8] = int_stack + 1305;
 Libderiv->deriv_classes[1][5][8] = int_stack + 1350;
 Libderiv->deriv_classes[1][6][8] = int_stack + 1413;
 Libderiv->deriv_classes[1][7][8] = int_stack + 1497;
 Libderiv->deriv_classes[1][8][8] = int_stack + 1605;
 Libderiv->deriv_classes[1][4][7] = int_stack + 1740;
 Libderiv->deriv_classes[1][5][7] = int_stack + 1785;
 Libderiv->deriv_classes[1][6][7] = int_stack + 1848;
 Libderiv->deriv_classes[1][7][7] = int_stack + 1932;
 Libderiv->deriv_classes[1][8][7] = int_stack + 2040;
 Libderiv->dvrr_classes[1][4] = int_stack + 2175;
 Libderiv->deriv_classes[1][4][6] = int_stack + 2220;
 Libderiv->dvrr_classes[1][5] = int_stack + 2265;
 Libderiv->deriv_classes[1][5][6] = int_stack + 2328;
 Libderiv->dvrr_classes[1][6] = int_stack + 2391;
 Libderiv->deriv_classes[1][6][6] = int_stack + 2475;
 Libderiv->dvrr_classes[1][7] = int_stack + 2559;
 Libderiv->deriv_classes[1][7][6] = int_stack + 2667;
 Libderiv->deriv_classes[1][8][6] = int_stack + 2775;
 Libderiv->deriv_classes[1][4][2] = int_stack + 2910;
 Libderiv->deriv_classes[1][5][2] = int_stack + 2955;
 Libderiv->deriv_classes[1][6][2] = int_stack + 3018;
 Libderiv->deriv_classes[1][7][2] = int_stack + 3102;
 Libderiv->deriv_classes[1][8][2] = int_stack + 3210;
 Libderiv->deriv_classes[1][4][1] = int_stack + 3345;
 Libderiv->deriv_classes[1][5][1] = int_stack + 3390;
 Libderiv->deriv_classes[1][6][1] = int_stack + 3453;
 Libderiv->deriv_classes[1][7][1] = int_stack + 3537;
 Libderiv->deriv_classes[1][8][1] = int_stack + 3645;
 Libderiv->deriv_classes[1][4][0] = int_stack + 3780;
 Libderiv->deriv_classes[1][5][0] = int_stack + 3825;
 Libderiv->deriv_classes[1][6][0] = int_stack + 3888;
 Libderiv->deriv_classes[1][7][0] = int_stack + 3972;
 Libderiv->deriv_classes[1][8][0] = int_stack + 4080;
 memset(int_stack,0,33720);

 Libderiv->dvrr_stack = int_stack + 12936;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0gg(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4215,int_stack+2265,int_stack+2175,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4350,int_stack+2391,int_stack+2265,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4539,int_stack+4350,int_stack+4215,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+4809,int_stack+2559,int_stack+2391,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+5061,int_stack+4809,int_stack+4350,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+5439,int_stack+5061,int_stack+4539,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5889,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2175,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6024,int_stack+108,int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2265,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6213,int_stack+6024,int_stack+5889, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4215,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+6483,int_stack+192,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2391,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+6735,int_stack+6483,int_stack+6024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4350,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+7113,int_stack+6735,int_stack+6213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4539,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+5889,int_stack+300,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2559,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+7563,int_stack+5889,int_stack+6483, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4809,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+5889,int_stack+7563,int_stack+6735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5061,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7563,int_stack+480,int_stack+435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2175, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7698,int_stack+543,int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2265, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7887,int_stack+7698,int_stack+7563, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4215, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8157,int_stack+627,int_stack+543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2391, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8409,int_stack+8157,int_stack+7698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4350, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6519,int_stack+8409,int_stack+7887, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4539, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+7563,int_stack+735,int_stack+627, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2559, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+0,int_stack+7563,int_stack+8157, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4809, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+7563,int_stack+0,int_stack+8409, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5061, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+915,int_stack+870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2175, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+135,int_stack+978,int_stack+915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2265, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+324,int_stack+135,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4215, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+594,int_stack+1062,int_stack+978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2391, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8193,int_stack+594,int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4350, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8571,int_stack+8193,int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4539, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+1170,int_stack+1062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2559, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+9021,int_stack+0,int_stack+594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4809, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+9021,int_stack+8193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5061, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8193,int_stack+1350,int_stack+1305, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8328,int_stack+1413,int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9021,int_stack+8328,int_stack+8193, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9291,int_stack+1497,int_stack+1413, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+9543,int_stack+9291,int_stack+8328, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+630,int_stack+9543,int_stack+9021, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+8193,int_stack+1605,int_stack+1497, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1080,int_stack+8193,int_stack+9291, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4809, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+9921,int_stack+1080,int_stack+9543, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1080,int_stack+1785,int_stack+1740, 0.0, zero_stack, 1.0, int_stack+2175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1215,int_stack+1848,int_stack+1785, 0.0, zero_stack, 1.0, int_stack+2265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1404,int_stack+1215,int_stack+1080, 0.0, zero_stack, 1.0, int_stack+4215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8193,int_stack+1932,int_stack+1848, 0.0, zero_stack, 1.0, int_stack+2391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+9021,int_stack+8193,int_stack+1215, 0.0, zero_stack, 1.0, int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+10551,int_stack+9021,int_stack+1404, 0.0, zero_stack, 1.0, int_stack+4539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+1080,int_stack+2040,int_stack+1932, 0.0, zero_stack, 1.0, int_stack+2559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+1404,int_stack+1080,int_stack+8193, 0.0, zero_stack, 1.0, int_stack+4809, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+11001,int_stack+1404,int_stack+9021, 0.0, zero_stack, 1.0, int_stack+5061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9021,int_stack+2328,int_stack+2220, 1.0, int_stack+2175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9156,int_stack+2475,int_stack+2328, 1.0, int_stack+2265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9345,int_stack+9156,int_stack+9021, 1.0, int_stack+4215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9615,int_stack+2667,int_stack+2475, 1.0, int_stack+2391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8193,int_stack+9615,int_stack+9156, 1.0, int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1080,int_stack+8193,int_stack+9345, 1.0, int_stack+4539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+4215,int_stack+2775,int_stack+2667, 1.0, int_stack+2559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+9021,int_stack+4215,int_stack+9615, 1.0, int_stack+4809, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+4215,int_stack+9021,int_stack+8193, 1.0, int_stack+5061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8193,int_stack+2955,int_stack+2910,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8328,int_stack+3018,int_stack+2955,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9021,int_stack+8328,int_stack+8193,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+9291,int_stack+3102,int_stack+3018,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+9543,int_stack+9291,int_stack+8328,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+4845,int_stack+9543,int_stack+9021,3);
 /*--- compute (p0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+8193,int_stack+3210,int_stack+3102,3);
 /*--- compute (p0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+1530,int_stack+8193,int_stack+9291,3);
 /*--- compute (p0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+2034,int_stack+1530,int_stack+9543,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1530,int_stack+3390,int_stack+3345,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1665,int_stack+3453,int_stack+3390,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8193,int_stack+1665,int_stack+1530,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+9021,int_stack+3537,int_stack+3453,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+9273,int_stack+9021,int_stack+1665,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+1530,int_stack+9273,int_stack+8193,3);
 /*--- compute (p0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+8193,int_stack+3645,int_stack+3537,3);
 /*--- compute (p0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+2664,int_stack+8193,int_stack+9021,3);
 /*--- compute (p0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+11631,int_stack+2664,int_stack+9273,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2664,int_stack+3825,int_stack+3780,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2799,int_stack+3888,int_stack+3825,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2988,int_stack+2799,int_stack+2664,3);
 /*--- compute (p0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+3258,int_stack+3972,int_stack+3888,3);
 /*--- compute (p0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+8193,int_stack+3258,int_stack+2799,3);
 /*--- compute (p0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+3510,int_stack+8193,int_stack+2988,3);
 /*--- compute (p0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+2664,int_stack+4080,int_stack+3972,3);
 /*--- compute (p0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+9021,int_stack+2664,int_stack+3258,3);
 /*--- compute (p0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+2664,int_stack+9021,int_stack+8193,3);
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+9021,int_stack+5889,int_stack+7113, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5439,3);
     Libderiv->ABCD[11] = int_stack + 9021;
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+12261,int_stack+7563,int_stack+6519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5439, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 12261;
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+5889,int_stack+0,int_stack+8571, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5439, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 5889;
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+6564,int_stack+9921,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 6564;
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+0,int_stack+11001,int_stack+10551, 0.0, zero_stack, 1.0, int_stack+5439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (p0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+7239,int_stack+4215,int_stack+1080, 1.0, int_stack+5439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 7239;
 /*--- compute (p0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+675,int_stack+2034,int_stack+4845,3);
     Libderiv->ABCD[2] = int_stack + 675;
 /*--- compute (p0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+7914,int_stack+11631,int_stack+1530,3);
     Libderiv->ABCD[1] = int_stack + 7914;
 /*--- compute (p0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+1350,int_stack+2664,int_stack+3510,3);
     Libderiv->ABCD[0] = int_stack + 1350;

}
