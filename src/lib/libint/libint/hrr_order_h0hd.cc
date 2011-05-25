#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_h0hd(Libint_t*, prim_data*);

  /* Computes quartets of (h0|hd) integrals */

REALTYPE *hrr_order_h0hd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][5] = int_stack + 0;
 Libint->vrr_classes[5][6] = int_stack + 441;
 Libint->vrr_classes[5][7] = int_stack + 1029;
 memset(int_stack,0,1785*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1785;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_h0hd(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1785,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3108,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4872,int_stack+3108,int_stack+1785,21);
 return int_stack+4872;}
