#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_pphp(Libint_t*, prim_data*);

  /* Computes quartets of (pp|hp) integrals */

REALTYPE *hrr_order_pphp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[1][6] = int_stack + 63;
 Libint->vrr_classes[2][5] = int_stack + 147;
 Libint->vrr_classes[2][6] = int_stack + 273;
 memset(int_stack,0,441*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 441;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_pphp(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+441,int_stack+63,int_stack+0,3);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+630,int_stack+273,int_stack+147,6);
 /*--- compute (pp|hp) ---*/
 hrr1_build_pp(Libint->AB,int_stack+1008,int_stack+630,int_stack+441,63);
 return int_stack+1008;}
