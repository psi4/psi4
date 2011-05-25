#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_pppp(Libint_t*, prim_data*);

  /* Computes quartets of (pp|pp) integrals */

REALTYPE *hrr_order_pppp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][1] = int_stack + 0;
 Libint->vrr_classes[1][2] = int_stack + 9;
 Libint->vrr_classes[2][1] = int_stack + 27;
 Libint->vrr_classes[2][2] = int_stack + 45;
 memset(int_stack,0,81*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 81;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_pppp(Libint, Data);
   Data++;
 }
 /*--- compute (p0|pp) ---*/
 hrr3_build_pp(Libint->CD,int_stack+81,int_stack+9,int_stack+0,3);
 /*--- compute (d0|pp) ---*/
 hrr3_build_pp(Libint->CD,int_stack+108,int_stack+45,int_stack+27,6);
 /*--- compute (pp|pp) ---*/
 hrr1_build_pp(Libint->AB,int_stack+0,int_stack+108,int_stack+81,9);
 return int_stack+0;}
