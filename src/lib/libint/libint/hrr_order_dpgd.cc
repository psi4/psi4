#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dpgd(Libint_t*, prim_data*);

  /* Computes quartets of (dp|gd) integrals */

REALTYPE *hrr_order_dpgd(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,1024*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1024;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dpgd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1024,int_stack+90,int_stack+0,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1294,int_stack+216,int_stack+90,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1672,int_stack+1294,int_stack+1024,6);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1024,int_stack+534,int_stack+384,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2212,int_stack+744,int_stack+534,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+2212,int_stack+1024,10);
 /*--- compute (dp|gd) ---*/
 hrr1_build_dp(Libint->AB,int_stack+2212,int_stack+0,int_stack+1672,90);
 return int_stack+2212;}
