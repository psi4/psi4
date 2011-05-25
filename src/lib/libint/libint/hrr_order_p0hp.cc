#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0hp(Libint_t*, prim_data*);

  /* Computes quartets of (p0|hp) integrals */

REALTYPE *hrr_order_p0hp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[1][6] = int_stack + 63;
 memset(int_stack,0,147*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 147;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0hp(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+147,int_stack+63,int_stack+0,3);
 return int_stack+147;}
