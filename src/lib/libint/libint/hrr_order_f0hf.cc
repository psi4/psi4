#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_f0hf(Libint_t*, prim_data*);

  /* Computes quartets of (f0|hf) integrals */

REALTYPE *hrr_order_f0hf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][5] = int_stack + 0;
 Libint->vrr_classes[3][6] = int_stack + 210;
 Libint->vrr_classes[3][7] = int_stack + 490;
 Libint->vrr_classes[3][8] = int_stack + 850;
 memset(int_stack,0,1300*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1300;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_f0hf(Libint, Data);
   Data++;
 }
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1300,int_stack+210,int_stack+0,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1930,int_stack+490,int_stack+210,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2770,int_stack+1930,int_stack+1300,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+4030,int_stack+850,int_stack+490,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+4030,int_stack+1930,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+4030,int_stack+0,int_stack+2770,10);
 return int_stack+4030;}
