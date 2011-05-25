#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hphf(Libint_t*, prim_data*);

  /* Computes quartets of (hp|hf) integrals */

REALTYPE *hrr_order_hphf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][5] = int_stack + 0;
 Libint->vrr_classes[5][6] = int_stack + 441;
 Libint->vrr_classes[5][7] = int_stack + 1029;
 Libint->vrr_classes[5][8] = int_stack + 1785;
 Libint->vrr_classes[6][5] = int_stack + 2730;
 Libint->vrr_classes[6][6] = int_stack + 3318;
 Libint->vrr_classes[6][7] = int_stack + 4102;
 Libint->vrr_classes[6][8] = int_stack + 5110;
 memset(int_stack,0,6370*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 6370;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hphf(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6370,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+7693,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9457,int_stack+7693,int_stack+6370,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+12103,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+14371,int_stack+12103,int_stack+7693,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+17899,int_stack+14371,int_stack+9457,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6370,int_stack+3318,int_stack+2730,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+8134,int_stack+4102,int_stack+3318,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10486,int_stack+8134,int_stack+6370,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14014,int_stack+5110,int_stack+4102,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+14014,int_stack+8134,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+22309,int_stack+0,int_stack+10486,28);
 /*--- compute (hp|hf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+0,int_stack+22309,int_stack+17899,210);
 return int_stack+0;}
