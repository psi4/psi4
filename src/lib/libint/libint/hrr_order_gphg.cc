#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gphg(Libint_t*, prim_data*);

  /* Computes quartets of (gp|hg) integrals */

REALTYPE *hrr_order_gphg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[4][8] = int_stack + 1275;
 Libint->vrr_classes[4][9] = int_stack + 1950;
 Libint->vrr_classes[5][5] = int_stack + 2775;
 Libint->vrr_classes[5][6] = int_stack + 3216;
 Libint->vrr_classes[5][7] = int_stack + 3804;
 Libint->vrr_classes[5][8] = int_stack + 4560;
 Libint->vrr_classes[5][9] = int_stack + 5505;
 memset(int_stack,0,6660*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 6660;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gphg(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6660,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+7605,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+8865,int_stack+7605,int_stack+6660,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+10755,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+12375,int_stack+10755,int_stack+7605,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+14895,int_stack+12375,int_stack+8865,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+6660,int_stack+1950,int_stack+1275,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+18045,int_stack+6660,int_stack+10755,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+6660,int_stack+18045,int_stack+12375,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+18045,int_stack+6660,int_stack+14895,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6660,int_stack+3216,int_stack+2775,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+7983,int_stack+3804,int_stack+3216,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9747,int_stack+7983,int_stack+6660,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+12393,int_stack+4560,int_stack+3804,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+12393,int_stack+7983,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+22770,int_stack+0,int_stack+9747,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+6660,int_stack+5505,int_stack+4560,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+27180,int_stack+6660,int_stack+12393,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+3528,int_stack+27180,int_stack+0,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+9408,int_stack+3528,int_stack+22770,21);
 /*--- compute (gp|hg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+22770,int_stack+9408,int_stack+18045,315);
 return int_stack+22770;}
