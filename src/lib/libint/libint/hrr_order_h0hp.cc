#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_h0hp(Libint_t*, prim_data*);

  /* Computes quartets of (h0|hp) integrals */

REALTYPE *hrr_order_h0hp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][5] = int_stack + 0;
 Libint->vrr_classes[5][6] = int_stack + 441;
 memset(int_stack,0,1029*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1029;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_h0hp(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1029,int_stack+441,int_stack+0,21);
 return int_stack+1029;}
