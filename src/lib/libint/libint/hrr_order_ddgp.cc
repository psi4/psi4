#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddgp(Libint_t*, prim_data*);

  /* Computes quartets of (dd|gp) integrals */

REALTYPE *hrr_order_ddgp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][4] = int_stack + 0;
 Libint->vrr_classes[2][5] = int_stack + 90;
 Libint->vrr_classes[3][4] = int_stack + 216;
 Libint->vrr_classes[3][5] = int_stack + 366;
 Libint->vrr_classes[4][4] = int_stack + 576;
 Libint->vrr_classes[4][5] = int_stack + 801;
 memset(int_stack,0,1116*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1116;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddgp(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1116,int_stack+90,int_stack+0,6);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1386,int_stack+366,int_stack+216,10);
 /*--- compute (dp|gp) ---*/
 hrr1_build_dp(Libint->AB,int_stack+1836,int_stack+1386,int_stack+1116,45);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2646,int_stack+801,int_stack+576,15);
 /*--- compute (fp|gp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+0,int_stack+2646,int_stack+1386,45);
 /*--- compute (dd|gp) ---*/
 hrr1_build_dd(Libint->AB,int_stack+2646,int_stack+0,int_stack+1836,45);
 return int_stack+2646;}
