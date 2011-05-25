#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fphf(Libint_t*, prim_data*);

  /* Computes quartets of (fp|hf) integrals */

REALTYPE *hrr_order_fphf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[3][7] = int_stack + 490;
 Libint->vrr_classes[3][8] = int_stack + 850;
 Libint->vrr_classes[4][5] = int_stack + 1300;
 Libint->vrr_classes[4][6] = int_stack + 1615;
 Libint->vrr_classes[4][7] = int_stack + 2035;
 Libint->vrr_classes[4][8] = int_stack + 2575;
 memset(int_stack,0,3250*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3250;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fphf(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3250,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3880,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4720,int_stack+3880,int_stack+3250,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+5980,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+7060,int_stack+5980,int_stack+3880,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+8740,int_stack+7060,int_stack+4720,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3250,int_stack+1615,int_stack+1300,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4195,int_stack+2035,int_stack+1615,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5455,int_stack+4195,int_stack+3250,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+2575,int_stack+2035,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+1620,int_stack+0,int_stack+4195,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+10840,int_stack+1620,int_stack+5455,15);
 /*--- compute (fp|hf) ---*/
 hrr1_build_fp(Libint->AB,int_stack+0,int_stack+10840,int_stack+8740,210);
 return int_stack+0;}
