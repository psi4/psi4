#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hpgg(Libint_t*, prim_data*);

  /* Computes quartets of (hp|gg) integrals */

REALTYPE *hrr_order_hpgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 Libint->vrr_classes[5][6] = int_stack + 756;
 Libint->vrr_classes[5][7] = int_stack + 1344;
 Libint->vrr_classes[5][8] = int_stack + 2100;
 Libint->vrr_classes[6][4] = int_stack + 3045;
 Libint->vrr_classes[6][5] = int_stack + 3465;
 Libint->vrr_classes[6][6] = int_stack + 4053;
 Libint->vrr_classes[6][7] = int_stack + 4837;
 Libint->vrr_classes[6][8] = int_stack + 5845;
 memset(int_stack,0,7105*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 7105;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hpgg(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+7105,int_stack+315,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8050,int_stack+756,int_stack+315,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+9373,int_stack+8050,int_stack+7105,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+11263,int_stack+1344,int_stack+756,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13027,int_stack+11263,int_stack+8050,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+15673,int_stack+13027,int_stack+9373,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+7105,int_stack+2100,int_stack+1344,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+18823,int_stack+7105,int_stack+11263,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+7105,int_stack+18823,int_stack+13027,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+18823,int_stack+7105,int_stack+15673,21);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+7105,int_stack+3465,int_stack+3045,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8365,int_stack+4053,int_stack+3465,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+10129,int_stack+8365,int_stack+7105,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12649,int_stack+4837,int_stack+4053,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+15001,int_stack+12649,int_stack+8365,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+0,int_stack+15001,int_stack+10129,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+7105,int_stack+5845,int_stack+4837,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+23548,int_stack+7105,int_stack+12649,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+4200,int_stack+23548,int_stack+15001,28);
 /*--- compute (i0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+10080,int_stack+4200,int_stack+0,28);
 /*--- compute (hp|gg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+23548,int_stack+10080,int_stack+18823,225);
 return int_stack+23548;}
