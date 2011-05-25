#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0fd(Libint_t*, prim_data*);

  /* Computes quartets of (f0|fd) integrals */

REALTYPE *hrr_order_f0fd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][3] = int_stack + 0;
 Libint->vrr_classes[3][4] = int_stack + 100;
 Libint->vrr_classes[3][5] = int_stack + 250;
 memset(int_stack,0,460*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 460;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0fd(Libint, Data);
   Data++;
 }
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+460,int_stack+100,int_stack+0,10);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+760,int_stack+250,int_stack+100,10);
 /*--- compute (f0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+1210,int_stack+760,int_stack+460,10);
 return int_stack+1210;}
