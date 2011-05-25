#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_00dd(Libint_t*, prim_data*);

  /* Computes quartets of (00|dd) integrals */

REALTYPE *hrr_order_00dd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[0][2] = int_stack + 0;
 Libint->vrr_classes[0][3] = int_stack + 6;
 Libint->vrr_classes[0][4] = int_stack + 16;
 memset(int_stack,0,31*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 31;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_00dd(Libint, Data);
   Data++;
 }
 /*--- compute (00|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+31,int_stack+6,int_stack+0,1);
 /*--- compute (00|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+49,int_stack+16,int_stack+6,1);
 /*--- compute (00|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+79,int_stack+49,int_stack+31,1);
 return int_stack+79;}
