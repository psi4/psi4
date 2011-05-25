#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppgp(Libint_t*, prim_data*);

  /* Computes quartets of (pp|gp) integrals */

REALTYPE *hrr_order_ppgp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][4] = int_stack + 0;
 Libint->vrr_classes[1][5] = int_stack + 45;
 Libint->vrr_classes[2][4] = int_stack + 108;
 Libint->vrr_classes[2][5] = int_stack + 198;
 memset(int_stack,0,324*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 324;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppgp(Libint, Data);
   Data++;
 }
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+324,int_stack+45,int_stack+0,3);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+459,int_stack+198,int_stack+108,6);
 /*--- compute (pp|gp) ---*/
 hrr1_build_pp(Libint->AB,int_stack+729,int_stack+459,int_stack+324,45);
 return int_stack+729;}
