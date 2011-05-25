#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppg0(Libint_t*, prim_data*);

  /* Computes quartets of (pp|g0) integrals */

REALTYPE *hrr_order_ppg0(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][4] = int_stack + 0;
 Libint->vrr_classes[2][4] = int_stack + 45;
 memset(int_stack,0,135*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 135;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppg0(Libint, Data);
   Data++;
 }
 /*--- compute (pp|g0) ---*/
 hrr1_build_pp(Libint->AB,int_stack+135,int_stack+45,int_stack+0,15);
 return int_stack+135;}
