#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_d0dd(Libint_t*, prim_data*);

  /* Computes quartets of (d0|dd) integrals */

REALTYPE *hrr_order_d0dd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][2] = int_stack + 0;
 Libint->vrr_classes[2][3] = int_stack + 36;
 Libint->vrr_classes[2][4] = int_stack + 96;
 memset(int_stack,0,186*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 186;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_d0dd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+186,int_stack+36,int_stack+0,6);
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+294,int_stack+96,int_stack+36,6);
 /*--- compute (d0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+474,int_stack+294,int_stack+186,6);
 return int_stack+474;}
