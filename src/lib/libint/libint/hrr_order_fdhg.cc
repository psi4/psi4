#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdhg(Libint_t*, prim_data*);

  /* Computes quartets of (fd|hg) integrals */

REALTYPE *hrr_order_fdhg(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,8510*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 8510;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdhg(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8510,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9140,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9980,int_stack+9140,int_stack+8510,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+11240,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+12320,int_stack+11240,int_stack+9140,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+14000,int_stack+12320,int_stack+9980,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+8510,int_stack+1300,int_stack+850,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+16100,int_stack+8510,int_stack+11240,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+8510,int_stack+16100,int_stack+12320,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+16100,int_stack+8510,int_stack+14000,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8510,int_stack+2165,int_stack+1850,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9455,int_stack+2585,int_stack+2165,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10715,int_stack+9455,int_stack+8510,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+12605,int_stack+3125,int_stack+2585,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+12605,int_stack+9455,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+19250,int_stack+0,int_stack+10715,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+8510,int_stack+3800,int_stack+3125,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+22400,int_stack+8510,int_stack+12605,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+8510,int_stack+22400,int_stack+0,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+22400,int_stack+8510,int_stack+19250,15);
 /*--- compute (fp|hg) ---*/
 hrr1_build_fp(Libint->AB,int_stack+27125,int_stack+22400,int_stack+16100,315);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8510,int_stack+5066,int_stack+4625,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9833,int_stack+5654,int_stack+5066,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+11597,int_stack+9833,int_stack+8510,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14243,int_stack+6410,int_stack+5654,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+16511,int_stack+14243,int_stack+9833,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+16511,int_stack+11597,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+8510,int_stack+7355,int_stack+6410,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+36575,int_stack+8510,int_stack+14243,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+4410,int_stack+36575,int_stack+16511,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+10290,int_stack+4410,int_stack+0,21);
 /*--- compute (gp|hg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+36575,int_stack+10290,int_stack+22400,315);
 /*--- compute (fd|hg) ---*/
 hrr1_build_fd(Libint->AB,int_stack+0,int_stack+36575,int_stack+27125,315);
 return int_stack+0;}
