#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dphf(Libint_t*, prim_data*);

  /* Computes quartets of (dp|hf) integrals */

REALTYPE *hrr_order_dphf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 Libint->vrr_classes[2][8] = int_stack + 510;
 Libint->vrr_classes[3][5] = int_stack + 780;
 Libint->vrr_classes[3][6] = int_stack + 990;
 Libint->vrr_classes[3][7] = int_stack + 1270;
 Libint->vrr_classes[3][8] = int_stack + 1630;
 memset(int_stack,0,2080*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2080;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dphf(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2080,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2458,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2962,int_stack+2458,int_stack+2080,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+3718,int_stack+510,int_stack+294,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+4366,int_stack+3718,int_stack+2458,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+5374,int_stack+4366,int_stack+2962,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2080,int_stack+990,int_stack+780,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2710,int_stack+1270,int_stack+990,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3550,int_stack+2710,int_stack+2080,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+1630,int_stack+1270,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+6634,int_stack+0,int_stack+2710,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+6634,int_stack+3550,10);
 /*--- compute (dp|hf) ---*/
 hrr1_build_dp(Libint->AB,int_stack+6634,int_stack+0,int_stack+5374,210);
 return int_stack+6634;}
