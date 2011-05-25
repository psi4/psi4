#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_h0ff(Libint_t*, prim_data*);

  /* Computes quartets of (h0|ff) integrals */

REALTYPE *hrr_order_h0ff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][3] = int_stack + 0;
 Libint->vrr_classes[5][4] = int_stack + 210;
 Libint->vrr_classes[5][5] = int_stack + 525;
 Libint->vrr_classes[5][6] = int_stack + 966;
 memset(int_stack,0,1554*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1554;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_h0ff(Libint, Data);
   Data++;
 }
 /*--- compute (h0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1554,int_stack+210,int_stack+0,21);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2184,int_stack+525,int_stack+210,21);
 /*--- compute (h0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+3129,int_stack+2184,int_stack+1554,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4389,int_stack+966,int_stack+525,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+4389,int_stack+2184,21);
 /*--- compute (h0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+4389,int_stack+0,int_stack+3129,21);
 return int_stack+4389;}
