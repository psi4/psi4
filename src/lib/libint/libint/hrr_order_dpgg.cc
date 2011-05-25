#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dpgg(Libint_t*, prim_data*);

  /* Computes quartets of (dp|gg) integrals */

REALTYPE *hrr_order_dpgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][4] = int_stack + 0;
 Libint->vrr_classes[2][5] = int_stack + 90;
 Libint->vrr_classes[2][6] = int_stack + 216;
 Libint->vrr_classes[2][7] = int_stack + 384;
 Libint->vrr_classes[2][8] = int_stack + 600;
 Libint->vrr_classes[3][4] = int_stack + 870;
 Libint->vrr_classes[3][5] = int_stack + 1020;
 Libint->vrr_classes[3][6] = int_stack + 1230;
 Libint->vrr_classes[3][7] = int_stack + 1510;
 Libint->vrr_classes[3][8] = int_stack + 1870;
 memset(int_stack,0,2320*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2320;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dpgg(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2320,int_stack+90,int_stack+0,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2590,int_stack+216,int_stack+90,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2968,int_stack+2590,int_stack+2320,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3508,int_stack+384,int_stack+216,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4012,int_stack+3508,int_stack+2590,6);
 /*--- compute (d0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+4768,int_stack+4012,int_stack+2968,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+2320,int_stack+600,int_stack+384,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+5668,int_stack+2320,int_stack+3508,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+2320,int_stack+5668,int_stack+4012,6);
 /*--- compute (d0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+5668,int_stack+2320,int_stack+4768,6);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2320,int_stack+1020,int_stack+870,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2770,int_stack+1230,int_stack+1020,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3400,int_stack+2770,int_stack+2320,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4300,int_stack+1510,int_stack+1230,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+4300,int_stack+2770,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+7018,int_stack+0,int_stack+3400,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+2320,int_stack+1870,int_stack+1510,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+8518,int_stack+2320,int_stack+4300,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+1260,int_stack+8518,int_stack+0,10);
 /*--- compute (f0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+8518,int_stack+1260,int_stack+7018,10);
 /*--- compute (dp|gg) ---*/
 hrr1_build_dp(Libint->AB,int_stack+0,int_stack+8518,int_stack+5668,225);
 return int_stack+0;}
