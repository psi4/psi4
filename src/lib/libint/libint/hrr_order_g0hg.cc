#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_g0hg(Libint_t*, prim_data*);

  /* Computes quartets of (g0|hg) integrals */

REALTYPE *hrr_order_g0hg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[4][8] = int_stack + 1275;
 Libint->vrr_classes[4][9] = int_stack + 1950;
 memset(int_stack,0,2775*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2775;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_g0hg(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2775,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3720,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4980,int_stack+3720,int_stack+2775,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+6870,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+8490,int_stack+6870,int_stack+3720,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+11010,int_stack+8490,int_stack+4980,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+2775,int_stack+1950,int_stack+1275,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+14160,int_stack+2775,int_stack+6870,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+14160,int_stack+8490,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+4200,int_stack+0,int_stack+11010,15);
 return int_stack+4200;}
