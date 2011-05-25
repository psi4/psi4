#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_00hg(Libint_t*, prim_data*);

  /* Computes quartets of (00|hg) integrals */

REALTYPE *hrr_order_00hg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[0][5] = int_stack + 0;
 Libint->vrr_classes[0][6] = int_stack + 21;
 Libint->vrr_classes[0][7] = int_stack + 49;
 Libint->vrr_classes[0][8] = int_stack + 85;
 Libint->vrr_classes[0][9] = int_stack + 130;
 memset(int_stack,0,185*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 185;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_00hg(Libint, Data);
   Data++;
 }
 /*--- compute (00|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+185,int_stack+21,int_stack+0,1);
 /*--- compute (00|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+248,int_stack+49,int_stack+21,1);
 /*--- compute (00|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+332,int_stack+248,int_stack+185,1);
 /*--- compute (00|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+458,int_stack+85,int_stack+49,1);
 /*--- compute (00|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+566,int_stack+458,int_stack+248,1);
 /*--- compute (00|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+734,int_stack+566,int_stack+332,1);
 /*--- compute (00|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+185,int_stack+130,int_stack+85,1);
 /*--- compute (00|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+944,int_stack+185,int_stack+458,1);
 /*--- compute (00|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+944,int_stack+566,1);
 /*--- compute (00|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+944,int_stack+0,int_stack+734,1);
 return int_stack+944;}
