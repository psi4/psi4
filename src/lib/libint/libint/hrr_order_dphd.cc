#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dphd(Libint_t*, prim_data*);

  /* Computes quartets of (dp|hd) integrals */

REALTYPE *hrr_order_dphd(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,1360*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1360;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dphd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1360,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1738,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2242,int_stack+1738,int_stack+1360,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1360,int_stack+720,int_stack+510,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2998,int_stack+1000,int_stack+720,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+2998,int_stack+1360,10);
 /*--- compute (dp|hd) ---*/
 hrr1_build_dp(Libint->AB,int_stack+2998,int_stack+0,int_stack+2242,126);
 return int_stack+2998;}
