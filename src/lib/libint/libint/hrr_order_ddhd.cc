#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddhd(Libint_t*, prim_data*);

  /* Computes quartets of (dd|hd) integrals */

REALTYPE *hrr_order_ddhd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 Libint->vrr_classes[3][5] = int_stack + 510;
 Libint->vrr_classes[3][6] = int_stack + 720;
 Libint->vrr_classes[3][7] = int_stack + 1000;
 Libint->vrr_classes[4][5] = int_stack + 1360;
 Libint->vrr_classes[4][6] = int_stack + 1675;
 Libint->vrr_classes[4][7] = int_stack + 2095;
 memset(int_stack,0,2635*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2635;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddhd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2635,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3013,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3517,int_stack+3013,int_stack+2635,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2635,int_stack+720,int_stack+510,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4273,int_stack+1000,int_stack+720,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+4273,int_stack+2635,10);
 /*--- compute (dp|hd) ---*/
 hrr1_build_dp(Libint->AB,int_stack+4273,int_stack+0,int_stack+3517,126);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2635,int_stack+1675,int_stack+1360,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6541,int_stack+2095,int_stack+1675,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7801,int_stack+6541,int_stack+2635,15);
 /*--- compute (fp|hd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+9691,int_stack+7801,int_stack+0,126);
 /*--- compute (dd|hd) ---*/
 hrr1_build_dd(Libint->AB,int_stack+13471,int_stack+9691,int_stack+4273,126);
 return int_stack+13471;}
