#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dphh(Libint_t*, prim_data*);

  /* Computes quartets of (dp|hh) integrals */

REALTYPE *hrr_order_dphh(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 Libint->vrr_classes[2][8] = int_stack + 510;
 Libint->vrr_classes[2][9] = int_stack + 780;
 Libint->vrr_classes[2][10] = int_stack + 1110;
 Libint->vrr_classes[3][5] = int_stack + 1506;
 Libint->vrr_classes[3][6] = int_stack + 1716;
 Libint->vrr_classes[3][7] = int_stack + 1996;
 Libint->vrr_classes[3][8] = int_stack + 2356;
 Libint->vrr_classes[3][9] = int_stack + 2806;
 Libint->vrr_classes[3][10] = int_stack + 3356;
 memset(int_stack,0,4016*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4016;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dphh(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4016,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4394,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4898,int_stack+4394,int_stack+4016,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+5654,int_stack+510,int_stack+294,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+6302,int_stack+5654,int_stack+4394,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+7310,int_stack+6302,int_stack+4898,6);
 /*--- compute (d0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+4016,int_stack+780,int_stack+510,6);
 /*--- compute (d0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+8570,int_stack+4016,int_stack+5654,6);
 /*--- compute (d0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+9866,int_stack+8570,int_stack+6302,6);
 /*--- compute (d0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+4826,int_stack+9866,int_stack+7310,6);
 /*--- compute (d0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+6716,int_stack+1110,int_stack+780,6);
 /*--- compute (d0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+11546,int_stack+6716,int_stack+4016,6);
 /*--- compute (d0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+13166,int_stack+11546,int_stack+8570,6);
 /*--- compute (d0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+6716,int_stack+13166,int_stack+9866,6);
 /*--- compute (d0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+9236,int_stack+6716,int_stack+4826,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4016,int_stack+1716,int_stack+1506,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4646,int_stack+1996,int_stack+1716,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5486,int_stack+4646,int_stack+4016,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+6746,int_stack+2356,int_stack+1996,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+11882,int_stack+6746,int_stack+4646,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+13562,int_stack+11882,int_stack+5486,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+4016,int_stack+2806,int_stack+2356,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+0,int_stack+4016,int_stack+6746,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+5366,int_stack+0,int_stack+11882,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+15662,int_stack+5366,int_stack+13562,10);
 /*--- compute (f0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+11882,int_stack+3356,int_stack+2806,10);
 /*--- compute (f0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+18812,int_stack+11882,int_stack+4016,10);
 /*--- compute (f0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+11882,int_stack+18812,int_stack+0,10);
 /*--- compute (f0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+0,int_stack+11882,int_stack+5366,10);
 /*--- compute (f0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+4200,int_stack+0,int_stack+15662,10);
 /*--- compute (dp|hh) ---*/
 hrr1_build_dp(Libint->AB,int_stack+11882,int_stack+4200,int_stack+9236,441);
 return int_stack+11882;}
