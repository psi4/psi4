#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdhf(Libint_t*, prim_data*);

  /* Computes quartets of (fd|hf) integrals */

REALTYPE *hrr_order_fdhf(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[5][5] = int_stack + 3250;
 Libint->vrr_classes[5][6] = int_stack + 3691;
 Libint->vrr_classes[5][7] = int_stack + 4279;
 Libint->vrr_classes[5][8] = int_stack + 5035;
 memset(int_stack,0,5980*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 5980;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdhf(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5980,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6610,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7450,int_stack+6610,int_stack+5980,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+8710,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+9790,int_stack+8710,int_stack+6610,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+11470,int_stack+9790,int_stack+7450,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5980,int_stack+1615,int_stack+1300,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6925,int_stack+2035,int_stack+1615,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+8185,int_stack+6925,int_stack+5980,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+2575,int_stack+2035,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+13570,int_stack+0,int_stack+6925,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+13570,int_stack+8185,15);
 /*--- compute (fp|hf) ---*/
 hrr1_build_fp(Libint->AB,int_stack+13570,int_stack+0,int_stack+11470,210);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5980,int_stack+3691,int_stack+3250,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+7303,int_stack+4279,int_stack+3691,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9067,int_stack+7303,int_stack+5980,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+19870,int_stack+5035,int_stack+4279,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+3150,int_stack+19870,int_stack+7303,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+19870,int_stack+3150,int_stack+9067,21);
 /*--- compute (gp|hf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+3150,int_stack+19870,int_stack+0,210);
 /*--- compute (fd|hf) ---*/
 hrr1_build_fd(Libint->AB,int_stack+19870,int_stack+3150,int_stack+13570,210);
 return int_stack+19870;}
