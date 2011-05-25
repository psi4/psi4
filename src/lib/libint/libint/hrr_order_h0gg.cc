#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_h0gg(Libint_t*, prim_data*);

  /* Computes quartets of (h0|gg) integrals */

REALTYPE *hrr_order_h0gg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 Libint->vrr_classes[5][6] = int_stack + 756;
 Libint->vrr_classes[5][7] = int_stack + 1344;
 Libint->vrr_classes[5][8] = int_stack + 2100;
 memset(int_stack,0,3045*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3045;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_h0gg(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3045,int_stack+315,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3990,int_stack+756,int_stack+315,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5313,int_stack+3990,int_stack+3045,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+7203,int_stack+1344,int_stack+756,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+8967,int_stack+7203,int_stack+3990,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+11613,int_stack+8967,int_stack+5313,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+3045,int_stack+2100,int_stack+1344,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+14763,int_stack+3045,int_stack+7203,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+14763,int_stack+8967,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+4410,int_stack+0,int_stack+11613,21);
 return int_stack+4410;}
