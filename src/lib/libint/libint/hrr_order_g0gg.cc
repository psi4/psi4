#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_g0gg(Libint_t*, prim_data*);

  /* Computes quartets of (g0|gg) integrals */

REALTYPE *hrr_order_g0gg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][4] = int_stack + 0;
 Libint->vrr_classes[4][5] = int_stack + 225;
 Libint->vrr_classes[4][6] = int_stack + 540;
 Libint->vrr_classes[4][7] = int_stack + 960;
 Libint->vrr_classes[4][8] = int_stack + 1500;
 memset(int_stack,0,2175*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2175;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_g0gg(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2175,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2850,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3795,int_stack+2850,int_stack+2175,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+5145,int_stack+960,int_stack+540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6405,int_stack+5145,int_stack+2850,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+8295,int_stack+6405,int_stack+3795,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+2175,int_stack+1500,int_stack+960,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+10545,int_stack+2175,int_stack+5145,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+10545,int_stack+6405,15);
 /*--- compute (g0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+10545,int_stack+0,int_stack+8295,15);
 return int_stack+10545;}
