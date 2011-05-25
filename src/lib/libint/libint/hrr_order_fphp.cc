#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fphp(Libint_t*, prim_data*);

  /* Computes quartets of (fp|hp) integrals */

REALTYPE *hrr_order_fphp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[4][5] = int_stack + 490;
 Libint->vrr_classes[4][6] = int_stack + 805;
 memset(int_stack,0,1225*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1225;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fphp(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1225,int_stack+210,int_stack+0,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1855,int_stack+805,int_stack+490,15);
 /*--- compute (fp|hp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+2800,int_stack+1855,int_stack+1225,63);
 return int_stack+2800;}
