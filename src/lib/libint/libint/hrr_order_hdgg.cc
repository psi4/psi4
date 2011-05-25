#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hdgg(Libint_t*, prim_data*);

  /* Computes quartets of (hd|gg) integrals */

REALTYPE *hrr_order_hdgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 Libint->vrr_classes[5][6] = int_stack + 756;
 Libint->vrr_classes[5][7] = int_stack + 1344;
 Libint->vrr_classes[5][8] = int_stack + 2100;
 Libint->vrr_classes[6][4] = int_stack + 3045;
 Libint->vrr_classes[6][5] = int_stack + 3465;
 Libint->vrr_classes[6][6] = int_stack + 4053;
 Libint->vrr_classes[6][7] = int_stack + 4837;
 Libint->vrr_classes[6][8] = int_stack + 5845;
 Libint->vrr_classes[7][4] = int_stack + 7105;
 Libint->vrr_classes[7][5] = int_stack + 7645;
 Libint->vrr_classes[7][6] = int_stack + 8401;
 Libint->vrr_classes[7][7] = int_stack + 9409;
 Libint->vrr_classes[7][8] = int_stack + 10705;
 memset(int_stack,0,12325*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 12325;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hdgg(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+12325,int_stack+315,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+13270,int_stack+756,int_stack+315,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+14593,int_stack+13270,int_stack+12325,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+16483,int_stack+1344,int_stack+756,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+18247,int_stack+16483,int_stack+13270,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+20893,int_stack+18247,int_stack+14593,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+12325,int_stack+2100,int_stack+1344,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+24043,int_stack+12325,int_stack+16483,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+12325,int_stack+24043,int_stack+18247,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+24043,int_stack+12325,int_stack+20893,21);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+12325,int_stack+3465,int_stack+3045,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+13585,int_stack+4053,int_stack+3465,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+15349,int_stack+13585,int_stack+12325,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+17869,int_stack+4837,int_stack+4053,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+20221,int_stack+17869,int_stack+13585,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+0,int_stack+20221,int_stack+15349,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+12325,int_stack+5845,int_stack+4837,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+28768,int_stack+12325,int_stack+17869,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+12325,int_stack+28768,int_stack+20221,28);
 /*--- compute (i0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+28768,int_stack+12325,int_stack+0,28);
 /*--- compute (hp|gg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+35068,int_stack+28768,int_stack+24043,225);
 /*--- compute (k0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+0,int_stack+7645,int_stack+7105,36);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1620,int_stack+8401,int_stack+7645,36);
 /*--- compute (k0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3888,int_stack+1620,int_stack+0,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12325,int_stack+9409,int_stack+8401,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+15349,int_stack+12325,int_stack+1620,36);
 /*--- compute (k0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+19885,int_stack+15349,int_stack+3888,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+10705,int_stack+9409,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+3888,int_stack+0,int_stack+12325,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+49243,int_stack+3888,int_stack+15349,36);
 /*--- compute (k0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+0,int_stack+49243,int_stack+19885,36);
 /*--- compute (ip|gg) ---*/
 hrr1_build_ip(Libint->AB,int_stack+8100,int_stack+0,int_stack+28768,225);
 /*--- compute (hd|gg) ---*/
 hrr1_build_hd(Libint->AB,int_stack+49243,int_stack+8100,int_stack+35068,225);
 return int_stack+49243;}
