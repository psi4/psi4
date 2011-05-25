#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gph0(Libint_t*, prim_data*);

  /* Computes quartets of (gp|h0) integrals */

REALTYPE *hrr_order_gph0(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 memset(int_stack,0,756*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 756;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gph0(Libint, Data);
   Data++;
 }
 /*--- compute (gp|h0) ---*/
 hrr1_build_gp(Libint->AB,int_stack+756,int_stack+315,int_stack+0,21);
 return int_stack+756;}
