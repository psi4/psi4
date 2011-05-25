#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0hd(Libint_t*, prim_data*);

  /* Computes quartets of (f0|hd) integrals */

REALTYPE *hrr_order_f0hd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[3][7] = int_stack + 490;
 memset(int_stack,0,850*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 850;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0hd(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+850,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1480,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2320,int_stack+1480,int_stack+850,10);
 return int_stack+2320;}
