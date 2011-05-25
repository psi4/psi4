#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_00pp(Libint_t*, prim_data*);

  /* Computes quartets of (00|pp) integrals */

REALTYPE *hrr_order_00pp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[0][1] = int_stack + 0;
 Libint->vrr_classes[0][2] = int_stack + 3;
 memset(int_stack,0,9*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 9;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_00pp(Libint, Data);
   Data++;
 }
 /*--- compute (00|pp) ---*/
 hrr3_build_pp(Libint->CD,int_stack+9,int_stack+3,int_stack+0,1);
 return int_stack+9;}
