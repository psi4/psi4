#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdhh(Libint_t*, prim_data*);

  /* Computes quartets of (fd|hh) integrals */

REALTYPE *hrr_order_fdhh(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,11546*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 11546;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdhh(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+11546,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12176,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13016,int_stack+12176,int_stack+11546,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14276,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+15356,int_stack+14276,int_stack+12176,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+17036,int_stack+15356,int_stack+13016,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+11546,int_stack+1300,int_stack+850,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+19136,int_stack+11546,int_stack+14276,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+21296,int_stack+19136,int_stack+15356,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+12896,int_stack+21296,int_stack+17036,10);
 /*--- compute (f0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+16046,int_stack+1850,int_stack+1300,10);
 /*--- compute (f0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+24096,int_stack+16046,int_stack+11546,10);
 /*--- compute (f0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+26796,int_stack+24096,int_stack+19136,10);
 /*--- compute (f0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+16046,int_stack+26796,int_stack+21296,10);
 /*--- compute (f0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+20246,int_stack+16046,int_stack+12896,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+11546,int_stack+2825,int_stack+2510,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12491,int_stack+3245,int_stack+2825,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13751,int_stack+12491,int_stack+11546,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+15641,int_stack+3785,int_stack+3245,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+17261,int_stack+15641,int_stack+12491,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+24656,int_stack+17261,int_stack+13751,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+11546,int_stack+4460,int_stack+3785,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+0,int_stack+11546,int_stack+15641,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+27806,int_stack+0,int_stack+17261,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+13571,int_stack+27806,int_stack+24656,15);
 /*--- compute (g0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+24656,int_stack+5285,int_stack+4460,15);
 /*--- compute (g0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+32006,int_stack+24656,int_stack+11546,15);
 /*--- compute (g0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+36056,int_stack+32006,int_stack+0,15);
 /*--- compute (g0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+41456,int_stack+36056,int_stack+27806,15);
 /*--- compute (g0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+24656,int_stack+41456,int_stack+13571,15);
 /*--- compute (fp|hh) ---*/
 hrr1_build_fp(Libint->AB,int_stack+31271,int_stack+24656,int_stack+20246,441);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+6716,int_stack+6275,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1323,int_stack+7304,int_stack+6716,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3087,int_stack+1323,int_stack+0,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+11546,int_stack+8060,int_stack+7304,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+13814,int_stack+11546,int_stack+1323,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+17342,int_stack+13814,int_stack+3087,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+9005,int_stack+8060,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+2835,int_stack+0,int_stack+11546,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+44501,int_stack+2835,int_stack+13814,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+50381,int_stack+44501,int_stack+17342,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+11546,int_stack+10160,int_stack+9005,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+15011,int_stack+11546,int_stack+0,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+7371,int_stack+15011,int_stack+2835,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+14931,int_stack+7371,int_stack+44501,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+14931,int_stack+50381,21);
 /*--- compute (gp|hh) ---*/
 hrr1_build_gp(Libint->AB,int_stack+44501,int_stack+0,int_stack+24656,441);
 /*--- compute (fd|hh) ---*/
 hrr1_build_fd(Libint->AB,int_stack+0,int_stack+44501,int_stack+31271,441);
 return int_stack+0;}
