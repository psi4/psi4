#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_pphf(Libint_t*, prim_data*);

  /* Computes quartets of (pp|hf) integrals */

REALTYPE *hrr_order_pphf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[1][6] = int_stack + 63;
 Libint->vrr_classes[1][7] = int_stack + 147;
 Libint->vrr_classes[1][8] = int_stack + 255;
 Libint->vrr_classes[2][5] = int_stack + 390;
 Libint->vrr_classes[2][6] = int_stack + 516;
 Libint->vrr_classes[2][7] = int_stack + 684;
 Libint->vrr_classes[2][8] = int_stack + 900;
 memset(int_stack,0,1170*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1170;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_pphf(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1170,int_stack+63,int_stack+0,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1359,int_stack+147,int_stack+63,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+1611,int_stack+1359,int_stack+1170,3);
 /*--- compute (p0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+1989,int_stack+255,int_stack+147,3);
 /*--- compute (p0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+2313,int_stack+1989,int_stack+1359,3);
 /*--- compute (p0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+2817,int_stack+2313,int_stack+1611,3);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1170,int_stack+516,int_stack+390,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1548,int_stack+684,int_stack+516,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2052,int_stack+1548,int_stack+1170,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+900,int_stack+684,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+3447,int_stack+0,int_stack+1548,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+3447,int_stack+2052,6);
 /*--- compute (pp|hf) ---*/
 hrr1_build_pp(Libint->AB,int_stack+3447,int_stack+0,int_stack+2817,210);
 return int_stack+3447;}
