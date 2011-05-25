#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_00dp(Libint_t*, prim_data*);

  /* Computes quartets of (00|dp) integrals */

REALTYPE *hrr_order_00dp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[0][2] = int_stack + 0;
 Libint->vrr_classes[0][3] = int_stack + 6;
 memset(int_stack,0,16*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 16;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_00dp(Libint, Data);
   Data++;
 }
 /*--- compute (00|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+16,int_stack+6,int_stack+0,1);
 return int_stack+16;}
