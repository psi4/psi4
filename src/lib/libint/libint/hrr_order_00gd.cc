#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_00gd(Libint_t*, prim_data*);

  /* Computes quartets of (00|gd) integrals */

REALTYPE *hrr_order_00gd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[0][4] = int_stack + 0;
 Libint->vrr_classes[0][5] = int_stack + 15;
 Libint->vrr_classes[0][6] = int_stack + 36;
 memset(int_stack,0,64*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 64;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_00gd(Libint, Data);
   Data++;
 }
 /*--- compute (00|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+64,int_stack+15,int_stack+0,1);
 /*--- compute (00|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+109,int_stack+36,int_stack+15,1);
 /*--- compute (00|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+172,int_stack+109,int_stack+64,1);
 return int_stack+172;}
