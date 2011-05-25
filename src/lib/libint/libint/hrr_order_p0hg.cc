#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0hg(Libint_t*, prim_data*);

  /* Computes quartets of (p0|hg) integrals */

REALTYPE *hrr_order_p0hg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[1][6] = int_stack + 63;
 Libint->vrr_classes[1][7] = int_stack + 147;
 Libint->vrr_classes[1][8] = int_stack + 255;
 Libint->vrr_classes[1][9] = int_stack + 390;
 memset(int_stack,0,555*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 555;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0hg(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+555,int_stack+63,int_stack+0,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+744,int_stack+147,int_stack+63,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+996,int_stack+744,int_stack+555,3);
 /*--- compute (p0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+1374,int_stack+255,int_stack+147,3);
 /*--- compute (p0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+1698,int_stack+1374,int_stack+744,3);
 /*--- compute (p0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+2202,int_stack+1698,int_stack+996,3);
 /*--- compute (p0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+555,int_stack+390,int_stack+255,3);
 /*--- compute (p0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+2832,int_stack+555,int_stack+1374,3);
 /*--- compute (p0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+2832,int_stack+1698,3);
 /*--- compute (p0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+2832,int_stack+0,int_stack+2202,3);
 return int_stack+2832;}
