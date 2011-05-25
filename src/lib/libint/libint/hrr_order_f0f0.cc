#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0f0(Libint_t*, prim_data*);

  /* Computes quartets of (f0|f0) integrals */

REALTYPE *hrr_order_f0f0(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][3] = int_stack + 0;
 memset(int_stack,0,100*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 100;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0f0(Libint, Data);
   Data++;
 }
 /*--- compute (f0|f0) ---*/
 return int_stack+0;}
