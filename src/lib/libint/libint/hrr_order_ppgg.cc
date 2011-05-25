#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppgg(Libint_t*, prim_data*);

  /* Computes quartets of (pp|gg) integrals */

REALTYPE *hrr_order_ppgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][4] = int_stack + 0;
 Libint->vrr_classes[1][5] = int_stack + 45;
 Libint->vrr_classes[1][6] = int_stack + 108;
 Libint->vrr_classes[1][7] = int_stack + 192;
 Libint->vrr_classes[1][8] = int_stack + 300;
 Libint->vrr_classes[2][4] = int_stack + 435;
 Libint->vrr_classes[2][5] = int_stack + 525;
 Libint->vrr_classes[2][6] = int_stack + 651;
 Libint->vrr_classes[2][7] = int_stack + 819;
 Libint->vrr_classes[2][8] = int_stack + 1035;
 memset(int_stack,0,1305*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1305;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppgg(Libint, Data);
   Data++;
 }
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1305,int_stack+45,int_stack+0,3);
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1440,int_stack+108,int_stack+45,3);
 /*--- compute (p0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1629,int_stack+1440,int_stack+1305,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1899,int_stack+192,int_stack+108,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2151,int_stack+1899,int_stack+1440,3);
 /*--- compute (p0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+2529,int_stack+2151,int_stack+1629,3);
 /*--- compute (p0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+1305,int_stack+300,int_stack+192,3);
 /*--- compute (p0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+2979,int_stack+1305,int_stack+1899,3);
 /*--- compute (p0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+1305,int_stack+2979,int_stack+2151,3);
 /*--- compute (p0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+2979,int_stack+1305,int_stack+2529,3);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1305,int_stack+525,int_stack+435,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1575,int_stack+651,int_stack+525,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1953,int_stack+1575,int_stack+1305,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+0,int_stack+819,int_stack+651,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3654,int_stack+0,int_stack+1575,6);
 /*--- compute (d0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+4410,int_stack+3654,int_stack+1953,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+1305,int_stack+1035,int_stack+819,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+1953,int_stack+1305,int_stack+0,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+1953,int_stack+3654,6);
 /*--- compute (d0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+1260,int_stack+0,int_stack+4410,6);
 /*--- compute (pp|gg) ---*/
 hrr1_build_pp(Libint->AB,int_stack+3654,int_stack+1260,int_stack+2979,225);
 return int_stack+3654;}
