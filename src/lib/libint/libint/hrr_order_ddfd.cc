#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddfd(Libint_t*, prim_data*);

  /* Computes quartets of (dd|fd) integrals */

REALTYPE *hrr_order_ddfd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][3] = int_stack + 0;
 Libint->vrr_classes[2][4] = int_stack + 60;
 Libint->vrr_classes[2][5] = int_stack + 150;
 Libint->vrr_classes[3][3] = int_stack + 276;
 Libint->vrr_classes[3][4] = int_stack + 376;
 Libint->vrr_classes[3][5] = int_stack + 526;
 Libint->vrr_classes[4][3] = int_stack + 736;
 Libint->vrr_classes[4][4] = int_stack + 886;
 Libint->vrr_classes[4][5] = int_stack + 1111;
 memset(int_stack,0,1426*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1426;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddfd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1426,int_stack+60,int_stack+0,6);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1606,int_stack+150,int_stack+60,6);
 /*--- compute (d0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+1876,int_stack+1606,int_stack+1426,6);
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1426,int_stack+376,int_stack+276,10);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2236,int_stack+526,int_stack+376,10);
 /*--- compute (f0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+0,int_stack+2236,int_stack+1426,10);
 /*--- compute (dp|fd) ---*/
 hrr1_build_dp(Libint->AB,int_stack+2236,int_stack+0,int_stack+1876,60);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1426,int_stack+886,int_stack+736,15);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3316,int_stack+1111,int_stack+886,15);
 /*--- compute (g0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+3991,int_stack+3316,int_stack+1426,15);
 /*--- compute (fp|fd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+4891,int_stack+3991,int_stack+0,60);
 /*--- compute (dd|fd) ---*/
 hrr1_build_dd(Libint->AB,int_stack+0,int_stack+4891,int_stack+2236,60);
 return int_stack+0;}
