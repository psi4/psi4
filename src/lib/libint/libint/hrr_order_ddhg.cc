#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddhg(Libint_t*, prim_data*);

  /* Computes quartets of (dd|hg) integrals */

REALTYPE *hrr_order_ddhg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 Libint->vrr_classes[2][8] = int_stack + 510;
 Libint->vrr_classes[2][9] = int_stack + 780;
 Libint->vrr_classes[3][5] = int_stack + 1110;
 Libint->vrr_classes[3][6] = int_stack + 1320;
 Libint->vrr_classes[3][7] = int_stack + 1600;
 Libint->vrr_classes[3][8] = int_stack + 1960;
 Libint->vrr_classes[3][9] = int_stack + 2410;
 Libint->vrr_classes[4][5] = int_stack + 2960;
 Libint->vrr_classes[4][6] = int_stack + 3275;
 Libint->vrr_classes[4][7] = int_stack + 3695;
 Libint->vrr_classes[4][8] = int_stack + 4235;
 Libint->vrr_classes[4][9] = int_stack + 4910;
 memset(int_stack,0,5735*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 5735;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddhg(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5735,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6113,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6617,int_stack+6113,int_stack+5735,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+7373,int_stack+510,int_stack+294,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+8021,int_stack+7373,int_stack+6113,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+9029,int_stack+8021,int_stack+6617,6);
 /*--- compute (d0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+5735,int_stack+780,int_stack+510,6);
 /*--- compute (d0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+10289,int_stack+5735,int_stack+7373,6);
 /*--- compute (d0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+5735,int_stack+10289,int_stack+8021,6);
 /*--- compute (d0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+10289,int_stack+5735,int_stack+9029,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5735,int_stack+1320,int_stack+1110,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6365,int_stack+1600,int_stack+1320,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7205,int_stack+6365,int_stack+5735,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+8465,int_stack+1960,int_stack+1600,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+8465,int_stack+6365,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+12179,int_stack+0,int_stack+7205,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+5735,int_stack+2410,int_stack+1960,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+14279,int_stack+5735,int_stack+8465,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+5735,int_stack+14279,int_stack+0,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+14279,int_stack+5735,int_stack+12179,10);
 /*--- compute (dp|hg) ---*/
 hrr1_build_dp(Libint->AB,int_stack+17429,int_stack+14279,int_stack+10289,315);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5735,int_stack+3275,int_stack+2960,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6680,int_stack+3695,int_stack+3275,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7940,int_stack+6680,int_stack+5735,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9830,int_stack+4235,int_stack+3695,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+11450,int_stack+9830,int_stack+6680,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+11450,int_stack+7940,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+5735,int_stack+4910,int_stack+4235,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+23099,int_stack+5735,int_stack+9830,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+3150,int_stack+23099,int_stack+11450,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+7350,int_stack+3150,int_stack+0,15);
 /*--- compute (fp|hg) ---*/
 hrr1_build_fp(Libint->AB,int_stack+23099,int_stack+7350,int_stack+14279,315);
 /*--- compute (dd|hg) ---*/
 hrr1_build_dd(Libint->AB,int_stack+0,int_stack+23099,int_stack+17429,315);
 return int_stack+0;}
