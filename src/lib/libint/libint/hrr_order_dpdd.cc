#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dpdd(Libint_t*, prim_data*);

  /* Computes quartets of (dp|dd) integrals */

REALTYPE *hrr_order_dpdd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][2] = int_stack + 0;
 Libint->vrr_classes[2][3] = int_stack + 36;
 Libint->vrr_classes[2][4] = int_stack + 96;
 Libint->vrr_classes[3][2] = int_stack + 186;
 Libint->vrr_classes[3][3] = int_stack + 246;
 Libint->vrr_classes[3][4] = int_stack + 346;
 memset(int_stack,0,496*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 496;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dpdd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+496,int_stack+36,int_stack+0,6);
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+604,int_stack+96,int_stack+36,6);
 /*--- compute (d0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+784,int_stack+604,int_stack+496,6);
 /*--- compute (f0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+496,int_stack+246,int_stack+186,10);
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1000,int_stack+346,int_stack+246,10);
 /*--- compute (f0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+0,int_stack+1000,int_stack+496,10);
 /*--- compute (dp|dd) ---*/
 hrr1_build_dp(Libint->AB,int_stack+1000,int_stack+0,int_stack+784,36);
 return int_stack+1000;}
