#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hpgd(Libint_t*, prim_data*);

  /* Computes quartets of (hp|gd) integrals */

REALTYPE *hrr_order_hpgd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 Libint->vrr_classes[5][6] = int_stack + 756;
 Libint->vrr_classes[6][4] = int_stack + 1344;
 Libint->vrr_classes[6][5] = int_stack + 1764;
 Libint->vrr_classes[6][6] = int_stack + 2352;
 memset(int_stack,0,3136*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3136;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hpgd(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3136,int_stack+315,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4081,int_stack+756,int_stack+315,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5404,int_stack+4081,int_stack+3136,21);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3136,int_stack+1764,int_stack+1344,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+2352,int_stack+1764,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7294,int_stack+0,int_stack+3136,28);
 /*--- compute (hp|gd) ---*/
 hrr1_build_hp(Libint->AB,int_stack+9814,int_stack+7294,int_stack+5404,90);
 return int_stack+9814;}
