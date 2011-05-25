#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0gd(Libint_t*, prim_data*);

  /* Computes quartets of (p0|gd) integrals */

REALTYPE *hrr_order_p0gd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][4] = int_stack + 0;
 Libint->vrr_classes[1][5] = int_stack + 45;
 Libint->vrr_classes[1][6] = int_stack + 108;
 memset(int_stack,0,192*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 192;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0gd(Libint, Data);
   Data++;
 }
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+192,int_stack+45,int_stack+0,3);
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+327,int_stack+108,int_stack+45,3);
 /*--- compute (p0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+516,int_stack+327,int_stack+192,3);
 return int_stack+516;}
