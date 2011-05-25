#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_d0hh(Libint_t*, prim_data*);

  /* Computes quartets of (d0|hh) integrals */

REALTYPE *hrr_order_d0hh(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,1506*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1506;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_d0hh(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1506,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1884,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2388,int_stack+1884,int_stack+1506,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+3144,int_stack+510,int_stack+294,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+3792,int_stack+3144,int_stack+1884,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+4800,int_stack+3792,int_stack+2388,6);
 /*--- compute (d0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+1506,int_stack+780,int_stack+510,6);
 /*--- compute (d0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+6060,int_stack+1506,int_stack+3144,6);
 /*--- compute (d0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+7356,int_stack+6060,int_stack+3792,6);
 /*--- compute (d0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+2316,int_stack+7356,int_stack+4800,6);
 /*--- compute (d0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+4206,int_stack+1110,int_stack+780,6);
 /*--- compute (d0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+9036,int_stack+4206,int_stack+1506,6);
 /*--- compute (d0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+0,int_stack+9036,int_stack+6060,6);
 /*--- compute (d0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+4206,int_stack+0,int_stack+7356,6);
 /*--- compute (d0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+6726,int_stack+4206,int_stack+2316,6);
 return int_stack+6726;}
