#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dddd(Libint_t*, prim_data*);

  /* Computes quartets of (dd|dd) integrals */

REALTYPE *hrr_order_dddd(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[4][2] = int_stack + 496;
 Libint->vrr_classes[4][3] = int_stack + 586;
 Libint->vrr_classes[4][4] = int_stack + 736;
 memset(int_stack,0,961*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 961;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dddd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+961,int_stack+36,int_stack+0,6);
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1069,int_stack+96,int_stack+36,6);
 /*--- compute (d0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+1249,int_stack+1069,int_stack+961,6);
 /*--- compute (f0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+961,int_stack+246,int_stack+186,10);
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1465,int_stack+346,int_stack+246,10);
 /*--- compute (f0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+0,int_stack+1465,int_stack+961,10);
 /*--- compute (dp|dd) ---*/
 hrr1_build_dp(Libint->AB,int_stack+1465,int_stack+0,int_stack+1249,36);
 /*--- compute (g0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+961,int_stack+586,int_stack+496,15);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+2113,int_stack+736,int_stack+586,15);
 /*--- compute (g0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+360,int_stack+2113,int_stack+961,15);
 /*--- compute (fp|dd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+2113,int_stack+360,int_stack+0,36);
 /*--- compute (dd|dd) ---*/
 hrr1_build_dd(Libint->AB,int_stack+0,int_stack+2113,int_stack+1465,36);
 return int_stack+0;}
