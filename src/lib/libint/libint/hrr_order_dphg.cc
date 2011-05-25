#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dphg(Libint_t*, prim_data*);

  /* Computes quartets of (dp|hg) integrals */

REALTYPE *hrr_order_dphg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 Libint->vrr_classes[2][8] = int_stack + 510;
 Libint->vrr_classes[2][9] = int_stack + 780;
 Libint->vrr_classes[3][5] = int_stack + 1110;
 Libint->vrr_classes[3][6] = int_stack + 1320;
 Libint->vrr_classes[3][7] = int_stack + 1600;
 Libint->vrr_classes[3][8] = int_stack + 1960;
 Libint->vrr_classes[3][9] = int_stack + 2410;
 memset(int_stack,0,2960*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2960;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dphg(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2960,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3338,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3842,int_stack+3338,int_stack+2960,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+4598,int_stack+510,int_stack+294,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+5246,int_stack+4598,int_stack+3338,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+6254,int_stack+5246,int_stack+3842,6);
 /*--- compute (d0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+2960,int_stack+780,int_stack+510,6);
 /*--- compute (d0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+7514,int_stack+2960,int_stack+4598,6);
 /*--- compute (d0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+2960,int_stack+7514,int_stack+5246,6);
 /*--- compute (d0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+7514,int_stack+2960,int_stack+6254,6);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2960,int_stack+1320,int_stack+1110,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3590,int_stack+1600,int_stack+1320,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4430,int_stack+3590,int_stack+2960,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+5690,int_stack+1960,int_stack+1600,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+5690,int_stack+3590,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+9404,int_stack+0,int_stack+4430,10);
 /*--- compute (f0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+2960,int_stack+2410,int_stack+1960,10);
 /*--- compute (f0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+11504,int_stack+2960,int_stack+5690,10);
 /*--- compute (f0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+1680,int_stack+11504,int_stack+0,10);
 /*--- compute (f0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+11504,int_stack+1680,int_stack+9404,10);
 /*--- compute (dp|hg) ---*/
 hrr1_build_dp(Libint->AB,int_stack+0,int_stack+11504,int_stack+7514,315);
 return int_stack+0;}
