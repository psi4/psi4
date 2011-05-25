#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddhp(Libint_t*, prim_data*);

  /* Computes quartets of (dd|hp) integrals */

REALTYPE *hrr_order_ddhp(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[3][5] = int_stack + 294;
 Libint->vrr_classes[3][6] = int_stack + 504;
 Libint->vrr_classes[4][5] = int_stack + 784;
 Libint->vrr_classes[4][6] = int_stack + 1099;
 memset(int_stack,0,1519*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1519;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddhp(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1519,int_stack+126,int_stack+0,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1897,int_stack+504,int_stack+294,10);
 /*--- compute (dp|hp) ---*/
 hrr1_build_dp(Libint->AB,int_stack+2527,int_stack+1897,int_stack+1519,63);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3661,int_stack+1099,int_stack+784,15);
 /*--- compute (fp|hp) ---*/
 hrr1_build_fp(Libint->AB,int_stack+0,int_stack+3661,int_stack+1897,63);
 /*--- compute (dd|hp) ---*/
 hrr1_build_dd(Libint->AB,int_stack+3661,int_stack+0,int_stack+2527,63);
 return int_stack+3661;}
