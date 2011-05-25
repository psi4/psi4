#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdh0(Libint_t*, prim_data*);

  /* Computes quartets of (fd|h0) integrals */

REALTYPE *hrr_order_fdh0(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[4][5] = int_stack + 210;
 Libint->vrr_classes[5][5] = int_stack + 525;
 memset(int_stack,0,966*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 966;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdh0(Libint, Data);
   Data++;
 }
 /*--- compute (fp|h0) ---*/
 hrr1_build_fp(Libint->AB,int_stack+966,int_stack+210,int_stack+0,21);
 /*--- compute (gp|h0) ---*/
 hrr1_build_gp(Libint->AB,int_stack+1596,int_stack+525,int_stack+210,21);
 /*--- compute (fd|h0) ---*/
 hrr1_build_fd(Libint->AB,int_stack+2541,int_stack+1596,int_stack+966,21);
 return int_stack+2541;}
