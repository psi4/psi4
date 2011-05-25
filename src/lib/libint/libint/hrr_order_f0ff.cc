#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0ff(Libint_t*, prim_data*);

  /* Computes quartets of (f0|ff) integrals */

REALTYPE *hrr_order_f0ff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][3] = int_stack + 0;
 Libint->vrr_classes[3][4] = int_stack + 100;
 Libint->vrr_classes[3][5] = int_stack + 250;
 Libint->vrr_classes[3][6] = int_stack + 460;
 memset(int_stack,0,740*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 740;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0ff(Libint, Data);
   Data++;
 }
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+740,int_stack+100,int_stack+0,10);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1040,int_stack+250,int_stack+100,10);
 /*--- compute (f0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+1490,int_stack+1040,int_stack+740,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2090,int_stack+460,int_stack+250,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+2090,int_stack+1040,10);
 /*--- compute (f0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+2090,int_stack+0,int_stack+1490,10);
 return int_stack+2090;}
