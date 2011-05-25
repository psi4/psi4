#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdhp(Libint_t*, prim_data*);

  /* Computes quartets of (fd|hp) integrals */

REALTYPE *hrr_order_fdhp(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,2254*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2254;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdhp(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2254,int_stack+210,int_stack+0,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2884,int_stack+805,int_stack+490,15);
 /*--- compute (fp|hp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+3829,int_stack+2884,int_stack+2254,63);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5719,int_stack+1666,int_stack+1225,21);
 /*--- compute (gp|hp) ---*/
 hrr1_build_gp(Libint->AB,int_stack+0,int_stack+5719,int_stack+2884,63);
 /*--- compute (fd|hp) ---*/
 hrr1_build_fd(Libint->AB,int_stack+5719,int_stack+0,int_stack+3829,63);
 return int_stack+5719;}
