#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gphp(Libint_t*, prim_data*);

  /* Computes quartets of (gp|hp) integrals */

REALTYPE *hrr_order_gphp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[5][5] = int_stack + 735;
 Libint->vrr_classes[5][6] = int_stack + 1176;
 memset(int_stack,0,1764*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1764;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gphp(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1764,int_stack+315,int_stack+0,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2709,int_stack+1176,int_stack+735,21);
 /*--- compute (gp|hp) ---*/
 hrr1_build_gp(Libint->AB,int_stack+4032,int_stack+2709,int_stack+1764,63);
 return int_stack+4032;}
