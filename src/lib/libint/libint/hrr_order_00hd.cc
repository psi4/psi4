#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_00hd(Libint_t*, prim_data*);

  /* Computes quartets of (00|hd) integrals */

REALTYPE *hrr_order_00hd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[0][5] = int_stack + 0;
 Libint->vrr_classes[0][6] = int_stack + 21;
 Libint->vrr_classes[0][7] = int_stack + 49;
 memset(int_stack,0,85*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 85;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_00hd(Libint, Data);
   Data++;
 }
 /*--- compute (00|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+85,int_stack+21,int_stack+0,1);
 /*--- compute (00|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+148,int_stack+49,int_stack+21,1);
 /*--- compute (00|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+232,int_stack+148,int_stack+85,1);
 return int_stack+232;}
