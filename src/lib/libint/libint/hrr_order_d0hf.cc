#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_d0hf(Libint_t*, prim_data*);

  /* Computes quartets of (d0|hf) integrals */

REALTYPE *hrr_order_d0hf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 Libint->vrr_classes[2][8] = int_stack + 510;
 memset(int_stack,0,780*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 780;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_d0hf(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+780,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1158,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+1662,int_stack+1158,int_stack+780,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+2418,int_stack+510,int_stack+294,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+2418,int_stack+1158,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+2418,int_stack+0,int_stack+1662,6);
 return int_stack+2418;}
