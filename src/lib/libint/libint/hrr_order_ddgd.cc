#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddgd(Libint_t*, prim_data*);

  /* Computes quartets of (dd|gd) integrals */

REALTYPE *hrr_order_ddgd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][4] = int_stack + 0;
 Libint->vrr_classes[2][5] = int_stack + 90;
 Libint->vrr_classes[2][6] = int_stack + 216;
 Libint->vrr_classes[3][4] = int_stack + 384;
 Libint->vrr_classes[3][5] = int_stack + 534;
 Libint->vrr_classes[3][6] = int_stack + 744;
 Libint->vrr_classes[4][4] = int_stack + 1024;
 Libint->vrr_classes[4][5] = int_stack + 1249;
 Libint->vrr_classes[4][6] = int_stack + 1564;
 memset(int_stack,0,1984*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1984;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddgd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1984,int_stack+90,int_stack+0,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2254,int_stack+216,int_stack+90,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2632,int_stack+2254,int_stack+1984,6);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1984,int_stack+534,int_stack+384,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3172,int_stack+744,int_stack+534,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+3172,int_stack+1984,10);
 /*--- compute (dp|gd) ---*/
 hrr1_build_dp(Libint->AB,int_stack+3172,int_stack+0,int_stack+2632,90);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1984,int_stack+1249,int_stack+1024,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4792,int_stack+1564,int_stack+1249,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5737,int_stack+4792,int_stack+1984,15);
 /*--- compute (fp|gd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+7087,int_stack+5737,int_stack+0,90);
 /*--- compute (dd|gd) ---*/
 hrr1_build_dd(Libint->AB,int_stack+9787,int_stack+7087,int_stack+3172,90);
 return int_stack+9787;}
