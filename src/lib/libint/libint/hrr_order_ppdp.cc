#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppdp(Libint_t*, prim_data*);

  /* Computes quartets of (pp|dp) integrals */

REALTYPE *hrr_order_ppdp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][2] = int_stack + 0;
 Libint->vrr_classes[1][3] = int_stack + 18;
 Libint->vrr_classes[2][2] = int_stack + 48;
 Libint->vrr_classes[2][3] = int_stack + 84;
 memset(int_stack,0,144*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 144;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppdp(Libint, Data);
   Data++;
 }
 /*--- compute (p0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+144,int_stack+18,int_stack+0,3);
 /*--- compute (d0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+198,int_stack+84,int_stack+48,6);
 /*--- compute (pp|dp) ---*/
 hrr1_build_pp(Libint->AB,int_stack+306,int_stack+198,int_stack+144,18);
 return int_stack+306;}
