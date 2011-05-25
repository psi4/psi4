#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0gg(Libint_t*, prim_data*);

  /* Computes quartets of (f0|gg) integrals */

REALTYPE *hrr_order_f0gg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 Libint->vrr_classes[3][6] = int_stack + 360;
 Libint->vrr_classes[3][7] = int_stack + 640;
 Libint->vrr_classes[3][8] = int_stack + 1000;
 memset(int_stack,0,1450*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1450;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0gg(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1450,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1900,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2530,int_stack+1900,int_stack+1450,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3430,int_stack+640,int_stack+360,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4270,int_stack+3430,int_stack+1900,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+5530,int_stack+4270,int_stack+2530,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+1450,int_stack+1000,int_stack+640,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+7030,int_stack+1450,int_stack+3430,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+7030,int_stack+4270,10);
 /*--- compute (f0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+2100,int_stack+0,int_stack+5530,10);
 return int_stack+2100;}
