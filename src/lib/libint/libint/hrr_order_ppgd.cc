#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppgd(Libint_t*, prim_data*);

  /* Computes quartets of (pp|gd) integrals */

REALTYPE *hrr_order_ppgd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][4] = int_stack + 0;
 Libint->vrr_classes[1][5] = int_stack + 45;
 Libint->vrr_classes[1][6] = int_stack + 108;
 Libint->vrr_classes[2][4] = int_stack + 192;
 Libint->vrr_classes[2][5] = int_stack + 282;
 Libint->vrr_classes[2][6] = int_stack + 408;
 memset(int_stack,0,576*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 576;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppgd(Libint, Data);
   Data++;
 }
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+576,int_stack+45,int_stack+0,3);
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+711,int_stack+108,int_stack+45,3);
 /*--- compute (p0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+900,int_stack+711,int_stack+576,3);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+576,int_stack+282,int_stack+192,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1170,int_stack+408,int_stack+282,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+1170,int_stack+576,6);
 /*--- compute (pp|gd) ---*/
 hrr1_build_pp(Libint->AB,int_stack+1170,int_stack+0,int_stack+900,90);
 return int_stack+1170;}
