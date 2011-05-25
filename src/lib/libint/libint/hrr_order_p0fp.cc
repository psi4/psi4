#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0fp(Libint_t*, prim_data*);

  /* Computes quartets of (p0|fp) integrals */

REALTYPE *hrr_order_p0fp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][3] = int_stack + 0;
 Libint->vrr_classes[1][4] = int_stack + 30;
 memset(int_stack,0,75*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 75;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0fp(Libint, Data);
   Data++;
 }
 /*--- compute (p0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+75,int_stack+30,int_stack+0,3);
 return int_stack+75;}
