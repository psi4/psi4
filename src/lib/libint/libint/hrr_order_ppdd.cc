#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppdd(Libint_t*, prim_data*);

  /* Computes quartets of (pp|dd) integrals */

REALTYPE *hrr_order_ppdd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][2] = int_stack + 0;
 Libint->vrr_classes[1][3] = int_stack + 18;
 Libint->vrr_classes[1][4] = int_stack + 48;
 Libint->vrr_classes[2][2] = int_stack + 93;
 Libint->vrr_classes[2][3] = int_stack + 129;
 Libint->vrr_classes[2][4] = int_stack + 189;
 memset(int_stack,0,279*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 279;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppdd(Libint, Data);
   Data++;
 }
 /*--- compute (p0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+279,int_stack+18,int_stack+0,3);
 /*--- compute (p0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+333,int_stack+48,int_stack+18,3);
 /*--- compute (p0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+423,int_stack+333,int_stack+279,3);
 /*--- compute (d0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+279,int_stack+129,int_stack+93,6);
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+531,int_stack+189,int_stack+129,6);
 /*--- compute (d0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+0,int_stack+531,int_stack+279,6);
 /*--- compute (pp|dd) ---*/
 hrr1_build_pp(Libint->AB,int_stack+531,int_stack+0,int_stack+423,36);
 return int_stack+531;}
