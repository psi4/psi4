#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0gg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|gg) integrals */

void d1hrr_order_d0gg(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[2][6][11] = int_stack + 216;
 Libderiv->deriv_classes[2][7][11] = int_stack + 384;
 Libderiv->deriv_classes[2][8][11] = int_stack + 600;
 Libderiv->deriv_classes[2][4][10] = int_stack + 870;
 Libderiv->deriv_classes[2][5][10] = int_stack + 960;
 Libderiv->deriv_classes[2][6][10] = int_stack + 1086;
 Libderiv->deriv_classes[2][7][10] = int_stack + 1254;
 Libderiv->deriv_classes[2][8][10] = int_stack + 1470;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1740;
 Libderiv->deriv_classes[2][5][9] = int_stack + 1830;
 Libderiv->deriv_classes[2][6][9] = int_stack + 1956;
 Libderiv->deriv_classes[2][7][9] = int_stack + 2124;
 Libderiv->deriv_classes[2][8][9] = int_stack + 2340;
 Libderiv->deriv_classes[2][4][8] = int_stack + 2610;
 Libderiv->deriv_classes[2][5][8] = int_stack + 2700;
 Libderiv->deriv_classes[2][6][8] = int_stack + 2826;
 Libderiv->deriv_classes[2][7][8] = int_stack + 2994;
 Libderiv->deriv_classes[2][8][8] = int_stack + 3210;
 Libderiv->deriv_classes[2][4][7] = int_stack + 3480;
 Libderiv->deriv_classes[2][5][7] = int_stack + 3570;
 Libderiv->deriv_classes[2][6][7] = int_stack + 3696;
 Libderiv->deriv_classes[2][7][7] = int_stack + 3864;
 Libderiv->deriv_classes[2][8][7] = int_stack + 4080;
 Libderiv->dvrr_classes[2][4] = int_stack + 4350;
 Libderiv->deriv_classes[2][4][6] = int_stack + 4440;
 Libderiv->dvrr_classes[2][5] = int_stack + 4530;
 Libderiv->deriv_classes[2][5][6] = int_stack + 4656;
 Libderiv->dvrr_classes[2][6] = int_stack + 4782;
 Libderiv->deriv_classes[2][6][6] = int_stack + 4950;
 Libderiv->dvrr_classes[2][7] = int_stack + 5118;
 Libderiv->deriv_classes[2][7][6] = int_stack + 5334;
 Libderiv->deriv_classes[2][8][6] = int_stack + 5550;
 Libderiv->deriv_classes[2][4][2] = int_stack + 5820;
 Libderiv->deriv_classes[2][5][2] = int_stack + 5910;
 Libderiv->deriv_classes[2][6][2] = int_stack + 6036;
 Libderiv->deriv_classes[2][7][2] = int_stack + 6204;
 Libderiv->deriv_classes[2][8][2] = int_stack + 6420;
 Libderiv->deriv_classes[2][4][1] = int_stack + 6690;
 Libderiv->deriv_classes[2][5][1] = int_stack + 6780;
 Libderiv->deriv_classes[2][6][1] = int_stack + 6906;
 Libderiv->deriv_classes[2][7][1] = int_stack + 7074;
 Libderiv->deriv_classes[2][8][1] = int_stack + 7290;
 Libderiv->deriv_classes[2][4][0] = int_stack + 7560;
 Libderiv->deriv_classes[2][5][0] = int_stack + 7650;
 Libderiv->deriv_classes[2][6][0] = int_stack + 7776;
 Libderiv->deriv_classes[2][7][0] = int_stack + 7944;
 Libderiv->deriv_classes[2][8][0] = int_stack + 8160;
 memset(int_stack,0,67440);

 Libderiv->dvrr_stack = int_stack + 25296;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0gg(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8430,int_stack+4530,int_stack+4350,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8700,int_stack+4782,int_stack+4530,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9078,int_stack+8700,int_stack+8430,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+9618,int_stack+5118,int_stack+4782,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+10122,int_stack+9618,int_stack+8700,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+10878,int_stack+10122,int_stack+9078,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11778,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4350,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12048,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4530,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+12426,int_stack+12048,int_stack+11778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8430,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+12966,int_stack+384,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4782,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13470,int_stack+12966,int_stack+12048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8700,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+14226,int_stack+13470,int_stack+12426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9078,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+11778,int_stack+600,int_stack+384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5118,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+15126,int_stack+11778,int_stack+12966, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9618,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+11778,int_stack+15126,int_stack+13470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10122,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15126,int_stack+960,int_stack+870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4350, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15396,int_stack+1086,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4530, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15774,int_stack+15396,int_stack+15126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8430, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+16314,int_stack+1254,int_stack+1086, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4782, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13038,int_stack+16314,int_stack+15396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8700, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+13038,int_stack+15774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9078, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+15126,int_stack+1470,int_stack+1254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5118, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+16818,int_stack+15126,int_stack+16314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9618, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+15126,int_stack+16818,int_stack+13038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10122, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13038,int_stack+1830,int_stack+1740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4350, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13308,int_stack+1956,int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4530, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13686,int_stack+13308,int_stack+13038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8430, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+16386,int_stack+2124,int_stack+1956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4782, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+16890,int_stack+16386,int_stack+13308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8700, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+900,int_stack+16890,int_stack+13686, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9078, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+13038,int_stack+2340,int_stack+2124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5118, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+17646,int_stack+13038,int_stack+16386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9618, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+18654,int_stack+17646,int_stack+16890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10122, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16386,int_stack+2700,int_stack+2610, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16656,int_stack+2826,int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17034,int_stack+16656,int_stack+16386, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+17574,int_stack+2994,int_stack+2826, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13038,int_stack+17574,int_stack+16656, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1800,int_stack+13038,int_stack+17034, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+16386,int_stack+3210,int_stack+2994, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+19914,int_stack+16386,int_stack+17574, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9618, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+16386,int_stack+19914,int_stack+13038, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13038,int_stack+3570,int_stack+3480, 0.0, zero_stack, 1.0, int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13308,int_stack+3696,int_stack+3570, 0.0, zero_stack, 1.0, int_stack+4530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13686,int_stack+13308,int_stack+13038, 0.0, zero_stack, 1.0, int_stack+8430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+19914,int_stack+3864,int_stack+3696, 0.0, zero_stack, 1.0, int_stack+4782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+17646,int_stack+19914,int_stack+13308, 0.0, zero_stack, 1.0, int_stack+8700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2700,int_stack+17646,int_stack+13686, 0.0, zero_stack, 1.0, int_stack+9078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+13038,int_stack+4080,int_stack+3864, 0.0, zero_stack, 1.0, int_stack+5118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+20418,int_stack+13038,int_stack+19914, 0.0, zero_stack, 1.0, int_stack+9618, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+21426,int_stack+20418,int_stack+17646, 0.0, zero_stack, 1.0, int_stack+10122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17646,int_stack+4656,int_stack+4440, 1.0, int_stack+4350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17916,int_stack+4950,int_stack+4656, 1.0, int_stack+4530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19914,int_stack+17916,int_stack+17646, 1.0, int_stack+8430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20454,int_stack+5334,int_stack+4950, 1.0, int_stack+4782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13038,int_stack+20454,int_stack+17916, 1.0, int_stack+8700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+17646,int_stack+13038,int_stack+19914, 1.0, int_stack+9078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+8430,int_stack+5550,int_stack+5334, 1.0, int_stack+5118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3600,int_stack+8430,int_stack+20454, 1.0, int_stack+9618, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+8430,int_stack+3600,int_stack+13038, 1.0, int_stack+10122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13038,int_stack+5910,int_stack+5820,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+13308,int_stack+6036,int_stack+5910,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+13686,int_stack+13308,int_stack+13038,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+3600,int_stack+6204,int_stack+6036,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+4104,int_stack+3600,int_stack+13308,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+4860,int_stack+4104,int_stack+13686,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+13038,int_stack+6420,int_stack+6204,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+9690,int_stack+13038,int_stack+3600,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+19914,int_stack+9690,int_stack+4104,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9690,int_stack+6780,int_stack+6690,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+9960,int_stack+6906,int_stack+6780,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10338,int_stack+9960,int_stack+9690,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+3600,int_stack+7074,int_stack+6906,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+4104,int_stack+3600,int_stack+9960,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+13038,int_stack+4104,int_stack+10338,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+9690,int_stack+7290,int_stack+7074,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+5760,int_stack+9690,int_stack+3600,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+22686,int_stack+5760,int_stack+4104,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5760,int_stack+7650,int_stack+7560,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6030,int_stack+7776,int_stack+7650,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6408,int_stack+6030,int_stack+5760,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+6948,int_stack+7944,int_stack+7776,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+3600,int_stack+6948,int_stack+6030,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+9690,int_stack+3600,int_stack+6408,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+5760,int_stack+8160,int_stack+7944,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+23946,int_stack+5760,int_stack+6948,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+5760,int_stack+23946,int_stack+3600,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+7020,int_stack+11778,int_stack+14226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10878,6);
     Libderiv->ABCD[11] = int_stack + 7020;
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+23946,int_stack+15126,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10878, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 23946;
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+13938,int_stack+18654,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10878, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 13938;
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+0,int_stack+16386,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+1350,int_stack+21426,int_stack+2700, 0.0, zero_stack, 1.0, int_stack+10878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 1350;
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+2700,int_stack+8430,int_stack+17646, 1.0, int_stack+10878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 2700;
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+10590,int_stack+19914,int_stack+4860,6);
     Libderiv->ABCD[2] = int_stack + 10590;
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+4050,int_stack+22686,int_stack+13038,6);
     Libderiv->ABCD[1] = int_stack + 4050;
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+11940,int_stack+5760,int_stack+9690,6);
     Libderiv->ABCD[0] = int_stack + 11940;

}
