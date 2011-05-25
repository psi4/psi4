#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0gg(Libint_t*, prim_data*);

  /* Computes quartets of (p0|gg) integrals */

REALTYPE *hrr_order_p0gg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][4] = int_stack + 0;
 Libint->vrr_classes[1][5] = int_stack + 45;
 Libint->vrr_classes[1][6] = int_stack + 108;
 Libint->vrr_classes[1][7] = int_stack + 192;
 Libint->vrr_classes[1][8] = int_stack + 300;
 memset(int_stack,0,435*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 435;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0gg(Libint, Data);
   Data++;
 }
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+435,int_stack+45,int_stack+0,3);
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+570,int_stack+108,int_stack+45,3);
 /*--- compute (p0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+759,int_stack+570,int_stack+435,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1029,int_stack+192,int_stack+108,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+1281,int_stack+1029,int_stack+570,3);
 /*--- compute (p0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+1659,int_stack+1281,int_stack+759,3);
 /*--- compute (p0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+435,int_stack+300,int_stack+192,3);
 /*--- compute (p0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+2109,int_stack+435,int_stack+1029,3);
 /*--- compute (p0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+2109,int_stack+1281,3);
 /*--- compute (p0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+2109,int_stack+0,int_stack+1659,3);
 return int_stack+2109;}
