#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ffhh(Libint_t*, prim_data*);

  /* Computes quartets of (ff|hh) integrals */

REALTYPE *hrr_order_ffhh(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[3][7] = int_stack + 490;
 Libint->vrr_classes[3][8] = int_stack + 850;
 Libint->vrr_classes[3][9] = int_stack + 1300;
 Libint->vrr_classes[3][10] = int_stack + 1850;
 Libint->vrr_classes[4][5] = int_stack + 2510;
 Libint->vrr_classes[4][6] = int_stack + 2825;
 Libint->vrr_classes[4][7] = int_stack + 3245;
 Libint->vrr_classes[4][8] = int_stack + 3785;
 Libint->vrr_classes[4][9] = int_stack + 4460;
 Libint->vrr_classes[4][10] = int_stack + 5285;
 Libint->vrr_classes[5][5] = int_stack + 6275;
 Libint->vrr_classes[5][6] = int_stack + 6716;
 Libint->vrr_classes[5][7] = int_stack + 7304;
 Libint->vrr_classes[5][8] = int_stack + 8060;
 Libint->vrr_classes[5][9] = int_stack + 9005;
 Libint->vrr_classes[5][10] = int_stack + 10160;
 Libint->vrr_classes[6][5] = int_stack + 11546;
 Libint->vrr_classes[6][6] = int_stack + 12134;
 Libint->vrr_classes[6][7] = int_stack + 12918;
 Libint->vrr_classes[6][8] = int_stack + 13926;
 Libint->vrr_classes[6][9] = int_stack + 15186;
 Libint->vrr_classes[6][10] = int_stack + 16726;
 memset(int_stack,0,18574*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 18574;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ffhh(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+18574,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+19204,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+20044,int_stack+19204,int_stack+18574,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+21304,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+22384,int_stack+21304,int_stack+19204,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+24064,int_stack+22384,int_stack+20044,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+18574,int_stack+1300,int_stack+850,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+26164,int_stack+18574,int_stack+21304,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+28324,int_stack+26164,int_stack+22384,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+19924,int_stack+28324,int_stack+24064,10);
 /*--- compute (f0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+23074,int_stack+1850,int_stack+1300,10);
 /*--- compute (f0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+31124,int_stack+23074,int_stack+18574,10);
 /*--- compute (f0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+33824,int_stack+31124,int_stack+26164,10);
 /*--- compute (f0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+23074,int_stack+33824,int_stack+28324,10);
 /*--- compute (f0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+27274,int_stack+23074,int_stack+19924,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+18574,int_stack+2825,int_stack+2510,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+19519,int_stack+3245,int_stack+2825,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+20779,int_stack+19519,int_stack+18574,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+22669,int_stack+3785,int_stack+3245,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+24289,int_stack+22669,int_stack+19519,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+31684,int_stack+24289,int_stack+20779,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+18574,int_stack+4460,int_stack+3785,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+0,int_stack+18574,int_stack+22669,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+34834,int_stack+0,int_stack+24289,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+20599,int_stack+34834,int_stack+31684,15);
 /*--- compute (g0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+31684,int_stack+5285,int_stack+4460,15);
 /*--- compute (g0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+39034,int_stack+31684,int_stack+18574,15);
 /*--- compute (g0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+43084,int_stack+39034,int_stack+0,15);
 /*--- compute (g0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+48484,int_stack+43084,int_stack+34834,15);
 /*--- compute (g0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+31684,int_stack+48484,int_stack+20599,15);
 /*--- compute (fp|hh) ---*/
 hrr1_build_fp(Libint->AB,int_stack+38299,int_stack+31684,int_stack+27274,441);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+6716,int_stack+6275,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1323,int_stack+7304,int_stack+6716,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3087,int_stack+1323,int_stack+0,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18574,int_stack+8060,int_stack+7304,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+20842,int_stack+18574,int_stack+1323,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+24370,int_stack+20842,int_stack+3087,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+9005,int_stack+8060,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+2835,int_stack+0,int_stack+18574,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+51529,int_stack+2835,int_stack+20842,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+57409,int_stack+51529,int_stack+24370,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+18574,int_stack+10160,int_stack+9005,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+22039,int_stack+18574,int_stack+0,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+64024,int_stack+22039,int_stack+2835,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+0,int_stack+64024,int_stack+51529,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+18574,int_stack+0,int_stack+57409,21);
 /*--- compute (gp|hh) ---*/
 hrr1_build_gp(Libint->AB,int_stack+51529,int_stack+18574,int_stack+31684,441);
 /*--- compute (fd|hh) ---*/
 hrr1_build_fd(Libint->AB,int_stack+71374,int_stack+51529,int_stack+38299,441);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+12134,int_stack+11546,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+12918,int_stack+12134,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4116,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+7644,int_stack+13926,int_stack+12918,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+27835,int_stack+7644,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+32539,int_stack+27835,int_stack+4116,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+15186,int_stack+13926,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+38419,int_stack+0,int_stack+7644,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+3780,int_stack+38419,int_stack+27835,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+97834,int_stack+3780,int_stack+32539,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+27835,int_stack+16726,int_stack+15186,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+106654,int_stack+27835,int_stack+0,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+27835,int_stack+106654,int_stack+38419,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+37915,int_stack+27835,int_stack+3780,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+37915,int_stack+97834,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+97834,int_stack+0,int_stack+18574,441);
 /*--- compute (gd|hh) ---*/
 hrr1_build_gd(Libint->AB,int_stack+0,int_stack+97834,int_stack+51529,441);
 /*--- compute (ff|hh) ---*/
 hrr1_build_ff(Libint->AB,int_stack+97834,int_stack+0,int_stack+71374,441);
 return int_stack+97834;}
