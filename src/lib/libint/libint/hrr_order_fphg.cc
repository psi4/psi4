#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fphg(Libint_t*, prim_data*);

  /* Computes quartets of (fp|hg) integrals */

REALTYPE *hrr_order_fphg(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,4625*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4625;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fphg(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4625,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+5255,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6095,int_stack+5255,int_stack+4625,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+7355,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+8435,int_stack+7355,int_stack+5255,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+10115,int_stack+8435,int_stack+6095,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+4625,int_stack+1300,int_stack+850,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+12215,int_stack+4625,int_stack+7355,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+4625,int_stack+12215,int_stack+8435,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+12215,int_stack+4625,int_stack+10115,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4625,int_stack+2165,int_stack+1850,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+5570,int_stack+2585,int_stack+2165,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6830,int_stack+5570,int_stack+4625,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+8720,int_stack+3125,int_stack+2585,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+8720,int_stack+5570,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+15365,int_stack+0,int_stack+6830,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+4625,int_stack+3800,int_stack+3125,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+18515,int_stack+4625,int_stack+8720,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+2520,int_stack+18515,int_stack+0,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+6720,int_stack+2520,int_stack+15365,15);
 /*--- compute (fp|hg) ---*/
 hrr1_build_fp(Libint->AB,int_stack+15365,int_stack+6720,int_stack+12215,315);
 return int_stack+15365;}
