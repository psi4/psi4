#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_00hf(Libint_t*, prim_data*);

  /* Computes quartets of (00|hf) integrals */

REALTYPE *hrr_order_00hf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[0][5] = int_stack + 0;
 Libint->vrr_classes[0][6] = int_stack + 21;
 Libint->vrr_classes[0][7] = int_stack + 49;
 Libint->vrr_classes[0][8] = int_stack + 85;
 memset(int_stack,0,130*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 130;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_00hf(Libint, Data);
   Data++;
 }
 /*--- compute (00|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+130,int_stack+21,int_stack+0,1);
 /*--- compute (00|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+193,int_stack+49,int_stack+21,1);
 /*--- compute (00|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+277,int_stack+193,int_stack+130,1);
 /*--- compute (00|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+403,int_stack+85,int_stack+49,1);
 /*--- compute (00|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+403,int_stack+193,1);
 /*--- compute (00|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+403,int_stack+0,int_stack+277,1);
 return int_stack+403;}
