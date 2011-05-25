#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hpgf(Libint_t*, prim_data*);

  /* Computes quartets of (hp|gf) integrals */

REALTYPE *hrr_order_hpgf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 Libint->vrr_classes[5][6] = int_stack + 756;
 Libint->vrr_classes[5][7] = int_stack + 1344;
 Libint->vrr_classes[6][4] = int_stack + 2100;
 Libint->vrr_classes[6][5] = int_stack + 2520;
 Libint->vrr_classes[6][6] = int_stack + 3108;
 Libint->vrr_classes[6][7] = int_stack + 3892;
 memset(int_stack,0,4900*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4900;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hpgf(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4900,int_stack+315,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5845,int_stack+756,int_stack+315,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7168,int_stack+5845,int_stack+4900,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9058,int_stack+1344,int_stack+756,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10822,int_stack+9058,int_stack+5845,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+13468,int_stack+10822,int_stack+7168,21);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4900,int_stack+2520,int_stack+2100,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6160,int_stack+3108,int_stack+2520,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7924,int_stack+6160,int_stack+4900,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10444,int_stack+3892,int_stack+3108,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+10444,int_stack+6160,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+3528,int_stack+0,int_stack+7924,28);
 /*--- compute (hp|gf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+16618,int_stack+3528,int_stack+13468,150);
 return int_stack+16618;}
