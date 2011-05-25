#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0hh(Libint_t*, prim_data*);

  /* Computes quartets of (p0|hh) integrals */

REALTYPE *hrr_order_p0hh(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,753*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 753;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0hh(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+753,int_stack+63,int_stack+0,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+942,int_stack+147,int_stack+63,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+1194,int_stack+942,int_stack+753,3);
 /*--- compute (p0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+1572,int_stack+255,int_stack+147,3);
 /*--- compute (p0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+1896,int_stack+1572,int_stack+942,3);
 /*--- compute (p0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+2400,int_stack+1896,int_stack+1194,3);
 /*--- compute (p0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+753,int_stack+390,int_stack+255,3);
 /*--- compute (p0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+3030,int_stack+753,int_stack+1572,3);
 /*--- compute (p0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+3678,int_stack+3030,int_stack+1896,3);
 /*--- compute (p0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+1158,int_stack+3678,int_stack+2400,3);
 /*--- compute (p0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+2103,int_stack+555,int_stack+390,3);
 /*--- compute (p0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+4518,int_stack+2103,int_stack+753,3);
 /*--- compute (p0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+0,int_stack+4518,int_stack+3030,3);
 /*--- compute (p0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+4518,int_stack+0,int_stack+3678,3);
 /*--- compute (p0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+2103,int_stack+4518,int_stack+1158,3);
 return int_stack+2103;}
