#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_d0fd(Libint_t*, prim_data*);

  /* Computes quartets of (d0|fd) integrals */

REALTYPE *hrr_order_d0fd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][3] = int_stack + 0;
 Libint->vrr_classes[2][4] = int_stack + 60;
 Libint->vrr_classes[2][5] = int_stack + 150;
 memset(int_stack,0,276*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 276;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_d0fd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+276,int_stack+60,int_stack+0,6);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+456,int_stack+150,int_stack+60,6);
 /*--- compute (d0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+726,int_stack+456,int_stack+276,6);
 return int_stack+726;}
