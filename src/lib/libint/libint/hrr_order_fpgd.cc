#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fpgd(Libint_t*, prim_data*);

  /* Computes quartets of (fp|gd) integrals */

REALTYPE *hrr_order_fpgd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 Libint->vrr_classes[3][6] = int_stack + 360;
 Libint->vrr_classes[4][4] = int_stack + 640;
 Libint->vrr_classes[4][5] = int_stack + 865;
 Libint->vrr_classes[4][6] = int_stack + 1180;
 memset(int_stack,0,1600*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1600;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fpgd(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1600,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2050,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2680,int_stack+2050,int_stack+1600,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1600,int_stack+865,int_stack+640,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3580,int_stack+1180,int_stack+865,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+3580,int_stack+1600,15);
 /*--- compute (fp|gd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+3580,int_stack+0,int_stack+2680,90);
 return int_stack+3580;}
