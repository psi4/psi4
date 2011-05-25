#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppff(Libint_t*, prim_data*);

  /* Computes quartets of (pp|ff) integrals */

REALTYPE *hrr_order_ppff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][3] = int_stack + 0;
 Libint->vrr_classes[1][4] = int_stack + 30;
 Libint->vrr_classes[1][5] = int_stack + 75;
 Libint->vrr_classes[1][6] = int_stack + 138;
 Libint->vrr_classes[2][3] = int_stack + 222;
 Libint->vrr_classes[2][4] = int_stack + 282;
 Libint->vrr_classes[2][5] = int_stack + 372;
 Libint->vrr_classes[2][6] = int_stack + 498;
 memset(int_stack,0,666*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 666;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppff(Libint, Data);
   Data++;
 }
 /*--- compute (p0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+666,int_stack+30,int_stack+0,3);
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+756,int_stack+75,int_stack+30,3);
 /*--- compute (p0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+891,int_stack+756,int_stack+666,3);
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1071,int_stack+138,int_stack+75,3);
 /*--- compute (p0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1260,int_stack+1071,int_stack+756,3);
 /*--- compute (p0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+1530,int_stack+1260,int_stack+891,3);
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+666,int_stack+282,int_stack+222,6);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+846,int_stack+372,int_stack+282,6);
 /*--- compute (d0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+1116,int_stack+846,int_stack+666,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1830,int_stack+498,int_stack+372,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+1830,int_stack+846,6);
 /*--- compute (d0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+1830,int_stack+0,int_stack+1116,6);
 /*--- compute (pp|ff) ---*/
 hrr1_build_pp(Libint->AB,int_stack+0,int_stack+1830,int_stack+1530,100);
 return int_stack+0;}
