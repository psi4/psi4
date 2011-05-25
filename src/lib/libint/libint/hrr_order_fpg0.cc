#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fpg0(Libint_t*, prim_data*);

  /* Computes quartets of (fp|g0) integrals */

REALTYPE *hrr_order_fpg0(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[4][4] = int_stack + 150;
 memset(int_stack,0,375*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 375;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fpg0(Libint, Data);
   Data++;
 }
 /*--- compute (fp|g0) ---*/
 hrr1_build_fp(Libint->AB,int_stack+375,int_stack+150,int_stack+0,15);
 return int_stack+375;}
