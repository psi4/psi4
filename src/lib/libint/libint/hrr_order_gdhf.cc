#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gdhf(Libint_t*, prim_data*);

  /* Computes quartets of (gd|hf) integrals */

REALTYPE *hrr_order_gdhf(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[6][5] = int_stack + 4680;
 Libint->vrr_classes[6][6] = int_stack + 5268;
 Libint->vrr_classes[6][7] = int_stack + 6052;
 Libint->vrr_classes[6][8] = int_stack + 7060;
 memset(int_stack,0,8320*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 8320;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gdhf(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8320,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9265,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10525,int_stack+9265,int_stack+8320,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+12415,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+14035,int_stack+12415,int_stack+9265,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+16555,int_stack+14035,int_stack+10525,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8320,int_stack+2391,int_stack+1950,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9643,int_stack+2979,int_stack+2391,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+11407,int_stack+9643,int_stack+8320,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14053,int_stack+3735,int_stack+2979,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+14053,int_stack+9643,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+19705,int_stack+0,int_stack+11407,21);
 /*--- compute (gp|hf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+24115,int_stack+19705,int_stack+16555,210);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+5268,int_stack+4680,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+6052,int_stack+5268,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+8320,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+11848,int_stack+7060,int_stack+6052,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+14872,int_stack+11848,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+14872,int_stack+8320,28);
 /*--- compute (hp|hf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+5880,int_stack+0,int_stack+19705,210);
 /*--- compute (gd|hf) ---*/
 hrr1_build_gd(Libint->AB,int_stack+33565,int_stack+5880,int_stack+24115,210);
 return int_stack+33565;}
