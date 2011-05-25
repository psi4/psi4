#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppfp(Libint_t*, prim_data*);

  /* Computes quartets of (pp|fp) integrals */

REALTYPE *hrr_order_ppfp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][3] = int_stack + 0;
 Libint->vrr_classes[1][4] = int_stack + 30;
 Libint->vrr_classes[2][3] = int_stack + 75;
 Libint->vrr_classes[2][4] = int_stack + 135;
 memset(int_stack,0,225*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 225;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppfp(Libint, Data);
   Data++;
 }
 /*--- compute (p0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+225,int_stack+30,int_stack+0,3);
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+315,int_stack+135,int_stack+75,6);
 /*--- compute (pp|fp) ---*/
 hrr1_build_pp(Libint->AB,int_stack+495,int_stack+315,int_stack+225,30);
 return int_stack+495;}
