#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddfp(Libint_t*, prim_data*);

  /* Computes quartets of (dd|fp) integrals */

REALTYPE *hrr_order_ddfp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][3] = int_stack + 0;
 Libint->vrr_classes[2][4] = int_stack + 60;
 Libint->vrr_classes[3][3] = int_stack + 150;
 Libint->vrr_classes[3][4] = int_stack + 250;
 Libint->vrr_classes[4][3] = int_stack + 400;
 Libint->vrr_classes[4][4] = int_stack + 550;
 memset(int_stack,0,775*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 775;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddfp(Libint, Data);
   Data++;
 }
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+775,int_stack+60,int_stack+0,6);
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+955,int_stack+250,int_stack+150,10);
 /*--- compute (dp|fp) ---*/
 hrr1_build_dp(Libint->AB,int_stack+1255,int_stack+955,int_stack+775,30);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1795,int_stack+550,int_stack+400,15);
 /*--- compute (fp|fp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+0,int_stack+1795,int_stack+955,30);
 /*--- compute (dd|fp) ---*/
 hrr1_build_dd(Libint->AB,int_stack+1795,int_stack+0,int_stack+1255,30);
 return int_stack+1795;}
