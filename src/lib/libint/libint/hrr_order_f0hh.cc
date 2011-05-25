#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0hh(Libint_t*, prim_data*);

  /* Computes quartets of (f0|hh) integrals */

REALTYPE *hrr_order_f0hh(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,2510*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2510;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0hh(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2510,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3140,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3980,int_stack+3140,int_stack+2510,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+5240,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+6320,int_stack+5240,int_stack+3140,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+8000,int_stack+6320,int_stack+3980,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+2510,int_stack+1300,int_stack+850,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+10100,int_stack+2510,int_stack+5240,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+12260,int_stack+10100,int_stack+6320,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+3860,int_stack+12260,int_stack+8000,10);
 /*--- compute (f0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+7010,int_stack+1850,int_stack+1300,10);
 /*--- compute (f0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+15060,int_stack+7010,int_stack+2510,10);
 /*--- compute (f0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+0,int_stack+15060,int_stack+10100,10);
 /*--- compute (f0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+7010,int_stack+0,int_stack+12260,10);
 /*--- compute (f0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+11210,int_stack+7010,int_stack+3860,10);
 return int_stack+11210;}
