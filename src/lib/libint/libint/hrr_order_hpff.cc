#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hpff(Libint_t*, prim_data*);

  /* Computes quartets of (hp|ff) integrals */

REALTYPE *hrr_order_hpff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][3] = int_stack + 0;
 Libint->vrr_classes[5][4] = int_stack + 210;
 Libint->vrr_classes[5][5] = int_stack + 525;
 Libint->vrr_classes[5][6] = int_stack + 966;
 Libint->vrr_classes[6][3] = int_stack + 1554;
 Libint->vrr_classes[6][4] = int_stack + 1834;
 Libint->vrr_classes[6][5] = int_stack + 2254;
 Libint->vrr_classes[6][6] = int_stack + 2842;
 memset(int_stack,0,3626*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3626;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hpff(Libint, Data);
   Data++;
 }
 /*--- compute (h0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+3626,int_stack+210,int_stack+0,21);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4256,int_stack+525,int_stack+210,21);
 /*--- compute (h0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+5201,int_stack+4256,int_stack+3626,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6461,int_stack+966,int_stack+525,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7784,int_stack+6461,int_stack+4256,21);
 /*--- compute (h0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+9674,int_stack+7784,int_stack+5201,21);
 /*--- compute (i0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+3626,int_stack+1834,int_stack+1554,28);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4466,int_stack+2254,int_stack+1834,28);
 /*--- compute (i0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+5726,int_stack+4466,int_stack+3626,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7406,int_stack+2842,int_stack+2254,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+7406,int_stack+4466,28);
 /*--- compute (i0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+2520,int_stack+0,int_stack+5726,28);
 /*--- compute (hp|ff) ---*/
 hrr1_build_hp(Libint->AB,int_stack+11774,int_stack+2520,int_stack+9674,100);
 return int_stack+11774;}
