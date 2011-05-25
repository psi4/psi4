#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ffhg(Libint_t*, prim_data*);

  /* Computes quartets of (ff|hg) integrals */

REALTYPE *hrr_order_ffhg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[3][7] = int_stack + 490;
 Libint->vrr_classes[3][8] = int_stack + 850;
 Libint->vrr_classes[3][9] = int_stack + 1300;
 Libint->vrr_classes[4][5] = int_stack + 1850;
 Libint->vrr_classes[4][6] = int_stack + 2165;
 Libint->vrr_classes[4][7] = int_stack + 2585;
 Libint->vrr_classes[4][8] = int_stack + 3125;
 Libint->vrr_classes[4][9] = int_stack + 3800;
 Libint->vrr_classes[5][5] = int_stack + 4625;
 Libint->vrr_classes[5][6] = int_stack + 5066;
 Libint->vrr_classes[5][7] = int_stack + 5654;
 Libint->vrr_classes[5][8] = int_stack + 6410;
 Libint->vrr_classes[5][9] = int_stack + 7355;
 Libint->vrr_classes[6][5] = int_stack + 8510;
 Libint->vrr_classes[6][6] = int_stack + 9098;
 Libint->vrr_classes[6][7] = int_stack + 9882;
 Libint->vrr_classes[6][8] = int_stack + 10890;
 Libint->vrr_classes[6][9] = int_stack + 12150;
 memset(int_stack,0,13690*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 13690;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ffhg(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+13690,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14320,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+15160,int_stack+14320,int_stack+13690,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+16420,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+17500,int_stack+16420,int_stack+14320,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+19180,int_stack+17500,int_stack+15160,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+13690,int_stack+1300,int_stack+850,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+21280,int_stack+13690,int_stack+16420,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+13690,int_stack+21280,int_stack+17500,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+21280,int_stack+13690,int_stack+19180,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+13690,int_stack+2165,int_stack+1850,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14635,int_stack+2585,int_stack+2165,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+15895,int_stack+14635,int_stack+13690,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+17785,int_stack+3125,int_stack+2585,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+17785,int_stack+14635,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+24430,int_stack+0,int_stack+15895,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+13690,int_stack+3800,int_stack+3125,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+27580,int_stack+13690,int_stack+17785,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+13690,int_stack+27580,int_stack+0,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+27580,int_stack+13690,int_stack+24430,15);
 /*--- compute (fp|hg) ---*/
 hrr1_build_fp(Libint->AB,int_stack+32305,int_stack+27580,int_stack+21280,315);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+13690,int_stack+5066,int_stack+4625,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+15013,int_stack+5654,int_stack+5066,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+16777,int_stack+15013,int_stack+13690,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+19423,int_stack+6410,int_stack+5654,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+21691,int_stack+19423,int_stack+15013,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+21691,int_stack+16777,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+13690,int_stack+7355,int_stack+6410,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+41755,int_stack+13690,int_stack+19423,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+13690,int_stack+41755,int_stack+21691,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+19570,int_stack+13690,int_stack+0,21);
 /*--- compute (gp|hg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+41755,int_stack+19570,int_stack+27580,315);
 /*--- compute (fd|hg) ---*/
 hrr1_build_fd(Libint->AB,int_stack+55930,int_stack+41755,int_stack+32305,315);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+9098,int_stack+8510,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+9882,int_stack+9098,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4116,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+13690,int_stack+10890,int_stack+9882,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+26185,int_stack+13690,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+30889,int_stack+26185,int_stack+4116,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+12150,int_stack+10890,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+3780,int_stack+0,int_stack+13690,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+9828,int_stack+3780,int_stack+26185,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+0,int_stack+9828,int_stack+30889,28);
 /*--- compute (hp|hg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+74830,int_stack+0,int_stack+19570,315);
 /*--- compute (gd|hg) ---*/
 hrr1_build_gd(Libint->AB,int_stack+0,int_stack+74830,int_stack+41755,315);
 /*--- compute (ff|hg) ---*/
 hrr1_build_ff(Libint->AB,int_stack+74830,int_stack+0,int_stack+55930,315);
 return int_stack+74830;}
