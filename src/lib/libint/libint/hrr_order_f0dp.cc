#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0dp(Libint_t*, prim_data*);

  /* Computes quartets of (f0|dp) integrals */

REALTYPE *hrr_order_f0dp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][2] = int_stack + 0;
 Libint->vrr_classes[3][3] = int_stack + 60;
 memset(int_stack,0,160*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 160;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0dp(Libint, Data);
   Data++;
 }
 /*--- compute (f0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+160,int_stack+60,int_stack+0,10);
 return int_stack+160;}
