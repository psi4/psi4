#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddhf(Libint_t*, prim_data*);

  /* Computes quartets of (dd|hf) integrals */

REALTYPE *hrr_order_ddhf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 Libint->vrr_classes[2][8] = int_stack + 510;
 Libint->vrr_classes[3][5] = int_stack + 780;
 Libint->vrr_classes[3][6] = int_stack + 990;
 Libint->vrr_classes[3][7] = int_stack + 1270;
 Libint->vrr_classes[3][8] = int_stack + 1630;
 Libint->vrr_classes[4][5] = int_stack + 2080;
 Libint->vrr_classes[4][6] = int_stack + 2395;
 Libint->vrr_classes[4][7] = int_stack + 2815;
 Libint->vrr_classes[4][8] = int_stack + 3355;
 memset(int_stack,0,4030*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4030;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddhf(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4030,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4408,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4912,int_stack+4408,int_stack+4030,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+5668,int_stack+510,int_stack+294,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+6316,int_stack+5668,int_stack+4408,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+7324,int_stack+6316,int_stack+4912,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4030,int_stack+990,int_stack+780,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4660,int_stack+1270,int_stack+990,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5500,int_stack+4660,int_stack+4030,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+1630,int_stack+1270,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+8584,int_stack+0,int_stack+4660,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+10264,int_stack+8584,int_stack+5500,10);
 /*--- compute (dp|hf) ---*/
 hrr1_build_dp(Libint->AB,int_stack+12364,int_stack+10264,int_stack+7324,210);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+2395,int_stack+2080,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+945,int_stack+2815,int_stack+2395,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4030,int_stack+945,int_stack+0,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+5920,int_stack+3355,int_stack+2815,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+7540,int_stack+5920,int_stack+945,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+7540,int_stack+4030,15);
 /*--- compute (fp|hf) ---*/
 hrr1_build_fp(Libint->AB,int_stack+3150,int_stack+0,int_stack+10264,210);
 /*--- compute (dd|hf) ---*/
 hrr1_build_dd(Libint->AB,int_stack+16144,int_stack+3150,int_stack+12364,210);
 return int_stack+16144;}
