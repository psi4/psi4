#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_d0dp(Libint_t*, prim_data*);

  /* Computes quartets of (d0|dp) integrals */

REALTYPE *hrr_order_d0dp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][2] = int_stack + 0;
 Libint->vrr_classes[2][3] = int_stack + 36;
 memset(int_stack,0,96*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 96;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_d0dp(Libint, Data);
   Data++;
 }
 /*--- compute (d0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+96,int_stack+36,int_stack+0,6);
 return int_stack+96;}
