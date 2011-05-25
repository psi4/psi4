#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddhh(Libint_t*, prim_data*);

  /* Computes quartets of (dd|hh) integrals */

REALTYPE *hrr_order_ddhh(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 Libint->vrr_classes[2][8] = int_stack + 510;
 Libint->vrr_classes[2][9] = int_stack + 780;
 Libint->vrr_classes[2][10] = int_stack + 1110;
 Libint->vrr_classes[3][5] = int_stack + 1506;
 Libint->vrr_classes[3][6] = int_stack + 1716;
 Libint->vrr_classes[3][7] = int_stack + 1996;
 Libint->vrr_classes[3][8] = int_stack + 2356;
 Libint->vrr_classes[3][9] = int_stack + 2806;
 Libint->vrr_classes[3][10] = int_stack + 3356;
 Libint->vrr_classes[4][5] = int_stack + 4016;
 Libint->vrr_classes[4][6] = int_stack + 4331;
 Libint->vrr_classes[4][7] = int_stack + 4751;
 Libint->vrr_classes[4][8] = int_stack + 5291;
 Libint->vrr_classes[4][9] = int_stack + 5966;
 Libint->vrr_classes[4][10] = int_stack + 6791;
 memset(int_stack,0,7781*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 7781;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddhh(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7781,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+8159,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+8663,int_stack+8159,int_stack+7781,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9419,int_stack+510,int_stack+294,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+10067,int_stack+9419,int_stack+8159,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+11075,int_stack+10067,int_stack+8663,6);
 /*--- compute (d0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+7781,int_stack+780,int_stack+510,6);
 /*--- compute (d0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+12335,int_stack+7781,int_stack+9419,6);
 /*--- compute (d0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+13631,int_stack+12335,int_stack+10067,6);
 /*--- compute (d0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+8591,int_stack+13631,int_stack+11075,6);
 /*--- compute (d0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+10481,int_stack+1110,int_stack+780,6);
 /*--- compute (d0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+15311,int_stack+10481,int_stack+7781,6);
 /*--- compute (d0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+16931,int_stack+15311,int_stack+12335,6);
 /*--- compute (d0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+10481,int_stack+16931,int_stack+13631,6);
 /*--- compute (d0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+13001,int_stack+10481,int_stack+8591,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7781,int_stack+1716,int_stack+1506,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+8411,int_stack+1996,int_stack+1716,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9251,int_stack+8411,int_stack+7781,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+10511,int_stack+2356,int_stack+1996,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+15647,int_stack+10511,int_stack+8411,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+17327,int_stack+15647,int_stack+9251,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+7781,int_stack+2806,int_stack+2356,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+0,int_stack+7781,int_stack+10511,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+9131,int_stack+0,int_stack+15647,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+19427,int_stack+9131,int_stack+17327,10);
 /*--- compute (f0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+15647,int_stack+3356,int_stack+2806,10);
 /*--- compute (f0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+22577,int_stack+15647,int_stack+7781,10);
 /*--- compute (f0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+15647,int_stack+22577,int_stack+0,10);
 /*--- compute (f0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+22577,int_stack+15647,int_stack+9131,10);
 /*--- compute (f0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+7781,int_stack+22577,int_stack+19427,10);
 /*--- compute (dp|hh) ---*/
 hrr1_build_dp(Libint->AB,int_stack+15647,int_stack+7781,int_stack+13001,441);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+23585,int_stack+4331,int_stack+4016,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+24530,int_stack+4751,int_stack+4331,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+24530,int_stack+23585,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+1890,int_stack+5291,int_stack+4751,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+12191,int_stack+1890,int_stack+24530,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+23585,int_stack+12191,int_stack+0,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+26735,int_stack+5966,int_stack+5291,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+28760,int_stack+26735,int_stack+1890,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+28760,int_stack+12191,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+32000,int_stack+0,int_stack+23585,15);
 /*--- compute (g0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+23585,int_stack+6791,int_stack+5966,15);
 /*--- compute (g0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+36725,int_stack+23585,int_stack+26735,15);
 /*--- compute (g0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+40775,int_stack+36725,int_stack+28760,15);
 /*--- compute (g0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+23585,int_stack+40775,int_stack+0,15);
 /*--- compute (g0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+23585,int_stack+32000,15);
 /*--- compute (fp|hh) ---*/
 hrr1_build_fp(Libint->AB,int_stack+23585,int_stack+0,int_stack+7781,441);
 /*--- compute (dd|hh) ---*/
 hrr1_build_dd(Libint->AB,int_stack+36815,int_stack+23585,int_stack+15647,441);
 return int_stack+36815;}
