#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_d0gg(Libint_t*, prim_data*);

  /* Computes quartets of (d0|gg) integrals */

REALTYPE *hrr_order_d0gg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][4] = int_stack + 0;
 Libint->vrr_classes[2][5] = int_stack + 90;
 Libint->vrr_classes[2][6] = int_stack + 216;
 Libint->vrr_classes[2][7] = int_stack + 384;
 Libint->vrr_classes[2][8] = int_stack + 600;
 memset(int_stack,0,870*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 870;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_d0gg(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+870,int_stack+90,int_stack+0,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1140,int_stack+216,int_stack+90,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1518,int_stack+1140,int_stack+870,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2058,int_stack+384,int_stack+216,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2562,int_stack+2058,int_stack+1140,6);
 /*--- compute (d0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+3318,int_stack+2562,int_stack+1518,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+870,int_stack+600,int_stack+384,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+4218,int_stack+870,int_stack+2058,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+4218,int_stack+2562,6);
 /*--- compute (d0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+4218,int_stack+0,int_stack+3318,6);
 return int_stack+4218;}
