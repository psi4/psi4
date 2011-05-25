#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ffhp(Libint_t*, prim_data*);

  /* Computes quartets of (ff|hp) integrals */

REALTYPE *hrr_order_ffhp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[4][5] = int_stack + 490;
 Libint->vrr_classes[4][6] = int_stack + 805;
 Libint->vrr_classes[5][5] = int_stack + 1225;
 Libint->vrr_classes[5][6] = int_stack + 1666;
 Libint->vrr_classes[6][5] = int_stack + 2254;
 Libint->vrr_classes[6][6] = int_stack + 2842;
 memset(int_stack,0,3626*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3626;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ffhp(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3626,int_stack+210,int_stack+0,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4256,int_stack+805,int_stack+490,15);
 /*--- compute (fp|hp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+5201,int_stack+4256,int_stack+3626,63);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7091,int_stack+1666,int_stack+1225,21);
 /*--- compute (gp|hp) ---*/
 hrr1_build_gp(Libint->AB,int_stack+8414,int_stack+7091,int_stack+4256,63);
 /*--- compute (fd|hp) ---*/
 hrr1_build_fd(Libint->AB,int_stack+11249,int_stack+8414,int_stack+5201,63);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3626,int_stack+2842,int_stack+2254,28);
 /*--- compute (hp|hp) ---*/
 hrr1_build_hp(Libint->AB,int_stack+15029,int_stack+3626,int_stack+7091,63);
 /*--- compute (gd|hp) ---*/
 hrr1_build_gd(Libint->AB,int_stack+0,int_stack+15029,int_stack+8414,63);
 /*--- compute (ff|hp) ---*/
 hrr1_build_ff(Libint->AB,int_stack+15029,int_stack+0,int_stack+11249,63);
 return int_stack+15029;}
