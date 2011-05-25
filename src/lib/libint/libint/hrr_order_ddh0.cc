#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddh0(Libint_t*, prim_data*);

  /* Computes quartets of (dd|h0) integrals */

REALTYPE *hrr_order_ddh0(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 126;
 Libint->vrr_classes[4][5] = int_stack + 336;
 memset(int_stack,0,651*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 651;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddh0(Libint, Data);
   Data++;
 }
 /*--- compute (dp|h0) ---*/
 hrr1_build_dp(Libint->AB,int_stack+651,int_stack+126,int_stack+0,21);
 /*--- compute (fp|h0) ---*/
 hrr1_build_fp(Libint->AB,int_stack+1029,int_stack+336,int_stack+126,21);
 /*--- compute (dd|h0) ---*/
 hrr1_build_dd(Libint->AB,int_stack+1659,int_stack+1029,int_stack+651,21);
 return int_stack+1659;}
