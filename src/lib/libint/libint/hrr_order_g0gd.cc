#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_g0gd(Libint_t*, prim_data*);

  /* Computes quartets of (g0|gd) integrals */

REALTYPE *hrr_order_g0gd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][4] = int_stack + 0;
 Libint->vrr_classes[4][5] = int_stack + 225;
 Libint->vrr_classes[4][6] = int_stack + 540;
 memset(int_stack,0,960*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 960;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_g0gd(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+960,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1635,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2580,int_stack+1635,int_stack+960,15);
 return int_stack+2580;}
