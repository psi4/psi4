#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dphp(Libint_t*, prim_data*);

  /* Computes quartets of (dp|hp) integrals */

REALTYPE *hrr_order_dphp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[3][5] = int_stack + 294;
 Libint->vrr_classes[3][6] = int_stack + 504;
 memset(int_stack,0,784*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 784;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dphp(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+784,int_stack+126,int_stack+0,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1162,int_stack+504,int_stack+294,10);
 /*--- compute (dp|hp) ---*/
 hrr1_build_dp(Libint->AB,int_stack+1792,int_stack+1162,int_stack+784,63);
 return int_stack+1792;}
