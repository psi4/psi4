#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0pp(Libint_t*, prim_data*);

  /* Computes quartets of (p0|pp) integrals */

REALTYPE *hrr_order_p0pp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][1] = int_stack + 0;
 Libint->vrr_classes[1][2] = int_stack + 9;
 memset(int_stack,0,27*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 27;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0pp(Libint, Data);
   Data++;
 }
 /*--- compute (p0|pp) ---*/
 hrr3_build_pp(Libint->CD,int_stack+27,int_stack+9,int_stack+0,3);
 return int_stack+27;}
