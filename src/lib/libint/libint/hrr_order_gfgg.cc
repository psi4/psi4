#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gfgg(Libint_t*, prim_data*);

  /* Computes quartets of (gf|gg) integrals */

REALTYPE *hrr_order_gfgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][4] = int_stack + 0;
 Libint->vrr_classes[4][5] = int_stack + 225;
 Libint->vrr_classes[4][6] = int_stack + 540;
 Libint->vrr_classes[4][7] = int_stack + 960;
 Libint->vrr_classes[4][8] = int_stack + 1500;
 Libint->vrr_classes[5][4] = int_stack + 2175;
 Libint->vrr_classes[5][5] = int_stack + 2490;
 Libint->vrr_classes[5][6] = int_stack + 2931;
 Libint->vrr_classes[5][7] = int_stack + 3519;
 Libint->vrr_classes[5][8] = int_stack + 4275;
 Libint->vrr_classes[6][4] = int_stack + 5220;
 Libint->vrr_classes[6][5] = int_stack + 5640;
 Libint->vrr_classes[6][6] = int_stack + 6228;
 Libint->vrr_classes[6][7] = int_stack + 7012;
 Libint->vrr_classes[6][8] = int_stack + 8020;
 Libint->vrr_classes[7][4] = int_stack + 9280;
 Libint->vrr_classes[7][5] = int_stack + 9820;
 Libint->vrr_classes[7][6] = int_stack + 10576;
 Libint->vrr_classes[7][7] = int_stack + 11584;
 Libint->vrr_classes[7][8] = int_stack + 12880;
 memset(int_stack,0,14500*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 14500;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gfgg(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+14500,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+15175,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+16120,int_stack+15175,int_stack+14500,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+17470,int_stack+960,int_stack+540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+18730,int_stack+17470,int_stack+15175,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+20620,int_stack+18730,int_stack+16120,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14500,int_stack+1500,int_stack+960,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+22870,int_stack+14500,int_stack+17470,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+14500,int_stack+22870,int_stack+18730,15);
 /*--- compute (g0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+22870,int_stack+14500,int_stack+20620,15);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+14500,int_stack+2490,int_stack+2175,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+15445,int_stack+2931,int_stack+2490,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+16768,int_stack+15445,int_stack+14500,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+18658,int_stack+3519,int_stack+2931,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+18658,int_stack+15445,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+26245,int_stack+0,int_stack+16768,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14500,int_stack+4275,int_stack+3519,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+29395,int_stack+14500,int_stack+18658,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+14500,int_stack+29395,int_stack+0,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+0,int_stack+14500,int_stack+26245,21);
 /*--- compute (gp|gg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+26245,int_stack+0,int_stack+22870,225);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+14500,int_stack+5640,int_stack+5220,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+15760,int_stack+6228,int_stack+5640,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+17524,int_stack+15760,int_stack+14500,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+20044,int_stack+7012,int_stack+6228,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+22396,int_stack+20044,int_stack+15760,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+36370,int_stack+22396,int_stack+17524,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14500,int_stack+8020,int_stack+7012,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+40570,int_stack+14500,int_stack+20044,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+14500,int_stack+40570,int_stack+22396,28);
 /*--- compute (i0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+40570,int_stack+14500,int_stack+36370,28);
 /*--- compute (hp|gg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+46870,int_stack+40570,int_stack+0,225);
 /*--- compute (gd|gg) ---*/
 hrr1_build_gd(Libint->AB,int_stack+61045,int_stack+46870,int_stack+26245,225);
 /*--- compute (k0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+0,int_stack+9820,int_stack+9280,36);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1620,int_stack+10576,int_stack+9820,36);
 /*--- compute (k0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3888,int_stack+1620,int_stack+0,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+7128,int_stack+11584,int_stack+10576,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+14500,int_stack+7128,int_stack+1620,36);
 /*--- compute (k0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+19036,int_stack+14500,int_stack+3888,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+12880,int_stack+11584,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+24436,int_stack+0,int_stack+7128,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+24436,int_stack+14500,36);
 /*--- compute (k0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+24436,int_stack+0,int_stack+19036,36);
 /*--- compute (ip|gg) ---*/
 hrr1_build_ip(Libint->AB,int_stack+0,int_stack+24436,int_stack+40570,225);
 /*--- compute (hd|gg) ---*/
 hrr1_build_hd(Libint->AB,int_stack+81295,int_stack+0,int_stack+46870,225);
 /*--- compute (gf|gg) ---*/
 hrr1_build_gf(Libint->AB,int_stack+0,int_stack+81295,int_stack+61045,225);
 return int_stack+0;}
