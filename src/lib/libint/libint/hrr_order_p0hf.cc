#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_p0hf(Libint_t*, prim_data*);

  /* Computes quartets of (p0|hf) integrals */

REALTYPE *hrr_order_p0hf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[1][6] = int_stack + 63;
 Libint->vrr_classes[1][7] = int_stack + 147;
 Libint->vrr_classes[1][8] = int_stack + 255;
 memset(int_stack,0,390*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 390;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_p0hf(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+390,int_stack+63,int_stack+0,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+579,int_stack+147,int_stack+63,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+831,int_stack+579,int_stack+390,3);
 /*--- compute (p0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+1209,int_stack+255,int_stack+147,3);
 /*--- compute (p0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+1209,int_stack+579,3);
 /*--- compute (p0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+1209,int_stack+0,int_stack+831,3);
 return int_stack+1209;}
