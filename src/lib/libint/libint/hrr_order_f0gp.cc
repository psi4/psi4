#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0gp(Libint_t*, prim_data*);

  /* Computes quartets of (f0|gp) integrals */

REALTYPE *hrr_order_f0gp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 memset(int_stack,0,360*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 360;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0gp(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+360,int_stack+150,int_stack+0,10);
 return int_stack+360;}
