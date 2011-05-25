#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_h0gp(Libint_t*, prim_data*);

  /* Computes quartets of (h0|gp) integrals */

REALTYPE *hrr_order_h0gp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 memset(int_stack,0,756*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 756;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_h0gp(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+756,int_stack+315,int_stack+0,21);
 return int_stack+756;}
