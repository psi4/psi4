#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_g0dd(Libint_t*, prim_data*);

  /* Computes quartets of (g0|dd) integrals */

REALTYPE *hrr_order_g0dd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][2] = int_stack + 0;
 Libint->vrr_classes[4][3] = int_stack + 90;
 Libint->vrr_classes[4][4] = int_stack + 240;
 memset(int_stack,0,465*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 465;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_g0dd(Libint, Data);
   Data++;
 }
 /*--- compute (g0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+465,int_stack+90,int_stack+0,15);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+735,int_stack+240,int_stack+90,15);
 /*--- compute (g0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+1185,int_stack+735,int_stack+465,15);
 return int_stack+1185;}
