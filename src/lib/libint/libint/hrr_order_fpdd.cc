#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fpdd(Libint_t*, prim_data*);

  /* Computes quartets of (fp|dd) integrals */

REALTYPE *hrr_order_fpdd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][2] = int_stack + 0;
 Libint->vrr_classes[3][3] = int_stack + 60;
 Libint->vrr_classes[3][4] = int_stack + 160;
 Libint->vrr_classes[4][2] = int_stack + 310;
 Libint->vrr_classes[4][3] = int_stack + 400;
 Libint->vrr_classes[4][4] = int_stack + 550;
 memset(int_stack,0,775*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 775;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fpdd(Libint, Data);
   Data++;
 }
 /*--- compute (f0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+775,int_stack+60,int_stack+0,10);
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+955,int_stack+160,int_stack+60,10);
 /*--- compute (f0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+1255,int_stack+955,int_stack+775,10);
 /*--- compute (g0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+775,int_stack+400,int_stack+310,15);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1615,int_stack+550,int_stack+400,15);
 /*--- compute (g0|dd) ---*/
 hrr3_build_dd(Libint->CD,int_stack+0,int_stack+1615,int_stack+775,15);
 /*--- compute (fp|dd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+1615,int_stack+0,int_stack+1255,36);
 return int_stack+1615;}
