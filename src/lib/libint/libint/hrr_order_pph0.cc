#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_pph0(Libint_t*, prim_data*);

  /* Computes quartets of (pp|h0) integrals */

REALTYPE *hrr_order_pph0(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[2][5] = int_stack + 63;
 memset(int_stack,0,189*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 189;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_pph0(Libint, Data);
   Data++;
 }
 /*--- compute (pp|h0) ---*/
 hrr1_build_pp(Libint->AB,int_stack+189,int_stack+63,int_stack+0,21);
 return int_stack+189;}
