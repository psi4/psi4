#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ffhf(Libint_t*, prim_data*);

  /* Computes quartets of (ff|hf) integrals */

REALTYPE *hrr_order_ffhf(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[6][5] = int_stack + 5980;
 Libint->vrr_classes[6][6] = int_stack + 6568;
 Libint->vrr_classes[6][7] = int_stack + 7352;
 Libint->vrr_classes[6][8] = int_stack + 8360;
 memset(int_stack,0,9620*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 9620;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ffhf(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9620,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10250,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+11090,int_stack+10250,int_stack+9620,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+12350,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+13430,int_stack+12350,int_stack+10250,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+15110,int_stack+13430,int_stack+11090,10);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9620,int_stack+1615,int_stack+1300,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10565,int_stack+2035,int_stack+1615,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+11825,int_stack+10565,int_stack+9620,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+2575,int_stack+2035,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+17210,int_stack+0,int_stack+10565,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+17210,int_stack+11825,15);
 /*--- compute (fp|hf) ---*/
 hrr1_build_fp(Libint->AB,int_stack+17210,int_stack+0,int_stack+15110,210);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9620,int_stack+3691,int_stack+3250,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10943,int_stack+4279,int_stack+3691,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+12707,int_stack+10943,int_stack+9620,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+23510,int_stack+5035,int_stack+4279,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+25778,int_stack+23510,int_stack+10943,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+29306,int_stack+25778,int_stack+12707,21);
 /*--- compute (gp|hf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+33716,int_stack+29306,int_stack+0,210);
 /*--- compute (fd|hf) ---*/
 hrr1_build_fd(Libint->AB,int_stack+43166,int_stack+33716,int_stack+17210,210);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+6568,int_stack+5980,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+7352,int_stack+6568,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9620,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+4116,int_stack+8360,int_stack+7352,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+13148,int_stack+4116,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+13148,int_stack+9620,28);
 /*--- compute (hp|hf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+5880,int_stack+0,int_stack+29306,210);
 /*--- compute (gd|hf) ---*/
 hrr1_build_gd(Libint->AB,int_stack+55766,int_stack+5880,int_stack+33716,210);
 /*--- compute (ff|hf) ---*/
 hrr1_build_ff(Libint->AB,int_stack+0,int_stack+55766,int_stack+43166,210);
 return int_stack+0;}
