#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gphh(Libint_t*, prim_data*);

  /* Computes quartets of (gp|hh) integrals */

REALTYPE *hrr_order_gphh(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[4][8] = int_stack + 1275;
 Libint->vrr_classes[4][9] = int_stack + 1950;
 Libint->vrr_classes[4][10] = int_stack + 2775;
 Libint->vrr_classes[5][5] = int_stack + 3765;
 Libint->vrr_classes[5][6] = int_stack + 4206;
 Libint->vrr_classes[5][7] = int_stack + 4794;
 Libint->vrr_classes[5][8] = int_stack + 5550;
 Libint->vrr_classes[5][9] = int_stack + 6495;
 Libint->vrr_classes[5][10] = int_stack + 7650;
 memset(int_stack,0,9036*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 9036;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gphh(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9036,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9981,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+11241,int_stack+9981,int_stack+9036,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+13131,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+14751,int_stack+13131,int_stack+9981,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+17271,int_stack+14751,int_stack+11241,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+9036,int_stack+1950,int_stack+1275,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+20421,int_stack+9036,int_stack+13131,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+23661,int_stack+20421,int_stack+14751,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+11061,int_stack+23661,int_stack+17271,15);
 /*--- compute (g0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+15786,int_stack+2775,int_stack+1950,15);
 /*--- compute (g0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+27861,int_stack+15786,int_stack+9036,15);
 /*--- compute (g0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+31911,int_stack+27861,int_stack+20421,15);
 /*--- compute (g0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+15786,int_stack+31911,int_stack+23661,15);
 /*--- compute (g0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+22086,int_stack+15786,int_stack+11061,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9036,int_stack+4206,int_stack+3765,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10359,int_stack+4794,int_stack+4206,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+12123,int_stack+10359,int_stack+9036,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14769,int_stack+5550,int_stack+4794,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+17037,int_stack+14769,int_stack+10359,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+28701,int_stack+17037,int_stack+12123,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+9036,int_stack+6495,int_stack+5550,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+33111,int_stack+9036,int_stack+14769,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+33111,int_stack+17037,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+11871,int_stack+0,int_stack+28701,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+28701,int_stack+7650,int_stack+6495,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+37647,int_stack+28701,int_stack+9036,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+43317,int_stack+37647,int_stack+33111,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+28701,int_stack+43317,int_stack+0,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+28701,int_stack+11871,21);
 /*--- compute (gp|hh) ---*/
 hrr1_build_gp(Libint->AB,int_stack+28701,int_stack+0,int_stack+22086,441);
 return int_stack+28701;}
