#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hphp(Libint_t*, prim_data*);

  /* Computes quartets of (hp|hp) integrals */

REALTYPE *hrr_order_hphp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][5] = int_stack + 0;
 Libint->vrr_classes[5][6] = int_stack + 441;
 Libint->vrr_classes[6][5] = int_stack + 1029;
 Libint->vrr_classes[6][6] = int_stack + 1617;
 memset(int_stack,0,2401*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2401;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hphp(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2401,int_stack+441,int_stack+0,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3724,int_stack+1617,int_stack+1029,28);
 /*--- compute (hp|hp) ---*/
 hrr1_build_hp(Libint->AB,int_stack+5488,int_stack+3724,int_stack+2401,63);
 return int_stack+5488;}
