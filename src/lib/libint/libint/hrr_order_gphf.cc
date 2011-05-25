#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gphf(Libint_t*, prim_data*);

  /* Computes quartets of (gp|hf) integrals */

REALTYPE *hrr_order_gphf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[4][8] = int_stack + 1275;
 Libint->vrr_classes[5][5] = int_stack + 1950;
 Libint->vrr_classes[5][6] = int_stack + 2391;
 Libint->vrr_classes[5][7] = int_stack + 2979;
 Libint->vrr_classes[5][8] = int_stack + 3735;
 memset(int_stack,0,4680*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4680;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gphf(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4680,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+5625,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6885,int_stack+5625,int_stack+4680,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+8775,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+10395,int_stack+8775,int_stack+5625,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+12915,int_stack+10395,int_stack+6885,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4680,int_stack+2391,int_stack+1950,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6003,int_stack+2979,int_stack+2391,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7767,int_stack+6003,int_stack+4680,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+10413,int_stack+3735,int_stack+2979,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+10413,int_stack+6003,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+16065,int_stack+0,int_stack+7767,21);
 /*--- compute (gp|hf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+0,int_stack+16065,int_stack+12915,210);
 return int_stack+0;}
