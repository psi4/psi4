#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_h0gf(Libint_t*, prim_data*);

  /* Computes quartets of (h0|gf) integrals */

REALTYPE *hrr_order_h0gf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 Libint->vrr_classes[5][6] = int_stack + 756;
 Libint->vrr_classes[5][7] = int_stack + 1344;
 memset(int_stack,0,2100*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2100;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_h0gf(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2100,int_stack+315,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3045,int_stack+756,int_stack+315,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+4368,int_stack+3045,int_stack+2100,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6258,int_stack+1344,int_stack+756,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+6258,int_stack+3045,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+6258,int_stack+0,int_stack+4368,21);
 return int_stack+6258;}
