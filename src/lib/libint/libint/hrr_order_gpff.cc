#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gpff(Libint_t*, prim_data*);

  /* Computes quartets of (gp|ff) integrals */

REALTYPE *hrr_order_gpff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][3] = int_stack + 0;
 Libint->vrr_classes[4][4] = int_stack + 150;
 Libint->vrr_classes[4][5] = int_stack + 375;
 Libint->vrr_classes[4][6] = int_stack + 690;
 Libint->vrr_classes[5][3] = int_stack + 1110;
 Libint->vrr_classes[5][4] = int_stack + 1320;
 Libint->vrr_classes[5][5] = int_stack + 1635;
 Libint->vrr_classes[5][6] = int_stack + 2076;
 memset(int_stack,0,2664*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2664;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gpff(Libint, Data);
   Data++;
 }
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+2664,int_stack+150,int_stack+0,15);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3114,int_stack+375,int_stack+150,15);
 /*--- compute (g0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+3789,int_stack+3114,int_stack+2664,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4689,int_stack+690,int_stack+375,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5634,int_stack+4689,int_stack+3114,15);
 /*--- compute (g0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+6984,int_stack+5634,int_stack+3789,15);
 /*--- compute (h0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+2664,int_stack+1320,int_stack+1110,21);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3294,int_stack+1635,int_stack+1320,21);
 /*--- compute (h0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+4239,int_stack+3294,int_stack+2664,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5499,int_stack+2076,int_stack+1635,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+5499,int_stack+3294,21);
 /*--- compute (h0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+1890,int_stack+0,int_stack+4239,21);
 /*--- compute (gp|ff) ---*/
 hrr1_build_gp(Libint->AB,int_stack+8484,int_stack+1890,int_stack+6984,100);
 return int_stack+8484;}
