#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppf0(Libint_t*, prim_data*);

  /* Computes quartets of (pp|f0) integrals */

REALTYPE *hrr_order_ppf0(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][3] = int_stack + 0;
 Libint->vrr_classes[2][3] = int_stack + 30;
 memset(int_stack,0,90*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 90;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppf0(Libint, Data);
   Data++;
 }
 /*--- compute (pp|f0) ---*/
 hrr1_build_pp(Libint->AB,int_stack+90,int_stack+30,int_stack+0,10);
 return int_stack+90;}
