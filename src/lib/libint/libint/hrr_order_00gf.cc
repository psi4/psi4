#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_00gf(Libint_t*, prim_data*);

  /* Computes quartets of (00|gf) integrals */

REALTYPE *hrr_order_00gf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[0][4] = int_stack + 0;
 Libint->vrr_classes[0][5] = int_stack + 15;
 Libint->vrr_classes[0][6] = int_stack + 36;
 Libint->vrr_classes[0][7] = int_stack + 64;
 memset(int_stack,0,100*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 100;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_00gf(Libint, Data);
   Data++;
 }
 /*--- compute (00|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+100,int_stack+15,int_stack+0,1);
 /*--- compute (00|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+145,int_stack+36,int_stack+15,1);
 /*--- compute (00|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+208,int_stack+145,int_stack+100,1);
 /*--- compute (00|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+298,int_stack+64,int_stack+36,1);
 /*--- compute (00|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+298,int_stack+145,1);
 /*--- compute (00|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+298,int_stack+0,int_stack+208,1);
 return int_stack+298;}
