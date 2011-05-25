#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dpdp(Libint_t*, prim_data*);

  /* Computes quartets of (dp|dp) integrals */

REALTYPE *hrr_order_dpdp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][2] = int_stack + 0;
 Libint->vrr_classes[2][3] = int_stack + 36;
 Libint->vrr_classes[3][2] = int_stack + 96;
 Libint->vrr_classes[3][3] = int_stack + 156;
 memset(int_stack,0,256*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 256;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dpdp(Libint, Data);
   Data++;
 }
 /*--- compute (d0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+256,int_stack+36,int_stack+0,6);
 /*--- compute (f0|dp) ---*/
 hrr3_build_dp(Libint->CD,int_stack+364,int_stack+156,int_stack+96,10);
 /*--- compute (dp|dp) ---*/
 hrr1_build_dp(Libint->AB,int_stack+544,int_stack+364,int_stack+256,18);
 return int_stack+544;}
