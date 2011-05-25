#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_g0hp(Libint_t*, prim_data*);

  /* Computes quartets of (g0|hp) integrals */

REALTYPE *hrr_order_g0hp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 memset(int_stack,0,735*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 735;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_g0hp(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+735,int_stack+315,int_stack+0,15);
 return int_stack+735;}
