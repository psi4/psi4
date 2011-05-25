#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_pphh(Libint_t*, prim_data*);

  /* Computes quartets of (pp|hh) integrals */

REALTYPE *hrr_order_pphh(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[1][6] = int_stack + 63;
 Libint->vrr_classes[1][7] = int_stack + 147;
 Libint->vrr_classes[1][8] = int_stack + 255;
 Libint->vrr_classes[1][9] = int_stack + 390;
 Libint->vrr_classes[1][10] = int_stack + 555;
 Libint->vrr_classes[2][5] = int_stack + 753;
 Libint->vrr_classes[2][6] = int_stack + 879;
 Libint->vrr_classes[2][7] = int_stack + 1047;
 Libint->vrr_classes[2][8] = int_stack + 1263;
 Libint->vrr_classes[2][9] = int_stack + 1533;
 Libint->vrr_classes[2][10] = int_stack + 1863;
 memset(int_stack,0,2259*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2259;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_pphh(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2259,int_stack+63,int_stack+0,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2448,int_stack+147,int_stack+63,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2700,int_stack+2448,int_stack+2259,3);
 /*--- compute (p0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+3078,int_stack+255,int_stack+147,3);
 /*--- compute (p0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+3402,int_stack+3078,int_stack+2448,3);
 /*--- compute (p0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+3906,int_stack+3402,int_stack+2700,3);
 /*--- compute (p0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+2259,int_stack+390,int_stack+255,3);
 /*--- compute (p0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+4536,int_stack+2259,int_stack+3078,3);
 /*--- compute (p0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+5184,int_stack+4536,int_stack+3402,3);
 /*--- compute (p0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+2664,int_stack+5184,int_stack+3906,3);
 /*--- compute (p0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+3609,int_stack+555,int_stack+390,3);
 /*--- compute (p0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+6024,int_stack+3609,int_stack+2259,3);
 /*--- compute (p0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+6834,int_stack+6024,int_stack+4536,3);
 /*--- compute (p0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+3609,int_stack+6834,int_stack+5184,3);
 /*--- compute (p0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+4869,int_stack+3609,int_stack+2664,3);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2259,int_stack+879,int_stack+753,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2637,int_stack+1047,int_stack+879,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3141,int_stack+2637,int_stack+2259,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+3897,int_stack+1263,int_stack+1047,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+6192,int_stack+3897,int_stack+2637,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+6192,int_stack+3141,6);
 /*--- compute (d0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+2259,int_stack+1533,int_stack+1263,6);
 /*--- compute (d0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+7200,int_stack+2259,int_stack+3897,6);
 /*--- compute (d0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+3069,int_stack+7200,int_stack+6192,6);
 /*--- compute (d0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+8496,int_stack+3069,int_stack+0,6);
 /*--- compute (d0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+0,int_stack+1863,int_stack+1533,6);
 /*--- compute (d0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+10386,int_stack+0,int_stack+2259,6);
 /*--- compute (d0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+0,int_stack+10386,int_stack+7200,6);
 /*--- compute (d0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+10386,int_stack+0,int_stack+3069,6);
 /*--- compute (d0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+10386,int_stack+8496,6);
 /*--- compute (pp|hh) ---*/
 hrr1_build_pp(Libint->AB,int_stack+6192,int_stack+0,int_stack+4869,441);
 return int_stack+6192;}
