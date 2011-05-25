#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_d0gf(Libint_t*, prim_data*);

  /* Computes quartets of (d0|gf) integrals */

REALTYPE *hrr_order_d0gf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][4] = int_stack + 0;
 Libint->vrr_classes[2][5] = int_stack + 90;
 Libint->vrr_classes[2][6] = int_stack + 216;
 Libint->vrr_classes[2][7] = int_stack + 384;
 memset(int_stack,0,600*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 600;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_d0gf(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+600,int_stack+90,int_stack+0,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+870,int_stack+216,int_stack+90,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1248,int_stack+870,int_stack+600,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1788,int_stack+384,int_stack+216,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+1788,int_stack+870,6);
 /*--- compute (d0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+1788,int_stack+0,int_stack+1248,6);
 return int_stack+1788;}
