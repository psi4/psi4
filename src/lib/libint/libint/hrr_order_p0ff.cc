#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0ff(Libint_t*, prim_data*);

  /* Computes quartets of (p0|ff) integrals */

REALTYPE *hrr_order_p0ff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][3] = int_stack + 0;
 Libint->vrr_classes[1][4] = int_stack + 30;
 Libint->vrr_classes[1][5] = int_stack + 75;
 Libint->vrr_classes[1][6] = int_stack + 138;
 memset(int_stack,0,222*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 222;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0ff(Libint, Data);
   Data++;
 }
 /*--- compute (p0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+222,int_stack+30,int_stack+0,3);
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+312,int_stack+75,int_stack+30,3);
 /*--- compute (p0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+447,int_stack+312,int_stack+222,3);
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+627,int_stack+138,int_stack+75,3);
 /*--- compute (p0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+627,int_stack+312,3);
 /*--- compute (p0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+627,int_stack+0,int_stack+447,3);
 return int_stack+627;}
