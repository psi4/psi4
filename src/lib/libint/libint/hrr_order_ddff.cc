#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddff(Libint_t*, prim_data*);

  /* Computes quartets of (dd|ff) integrals */

REALTYPE *hrr_order_ddff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][3] = int_stack + 0;
 Libint->vrr_classes[2][4] = int_stack + 60;
 Libint->vrr_classes[2][5] = int_stack + 150;
 Libint->vrr_classes[2][6] = int_stack + 276;
 Libint->vrr_classes[3][3] = int_stack + 444;
 Libint->vrr_classes[3][4] = int_stack + 544;
 Libint->vrr_classes[3][5] = int_stack + 694;
 Libint->vrr_classes[3][6] = int_stack + 904;
 Libint->vrr_classes[4][3] = int_stack + 1184;
 Libint->vrr_classes[4][4] = int_stack + 1334;
 Libint->vrr_classes[4][5] = int_stack + 1559;
 Libint->vrr_classes[4][6] = int_stack + 1874;
 memset(int_stack,0,2294*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2294;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddff(Libint, Data);
   Data++;
 }
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+2294,int_stack+60,int_stack+0,6);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2474,int_stack+150,int_stack+60,6);
 /*--- compute (d0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+2744,int_stack+2474,int_stack+2294,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3104,int_stack+276,int_stack+150,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3482,int_stack+3104,int_stack+2474,6);
 /*--- compute (d0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+4022,int_stack+3482,int_stack+2744,6);
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+2294,int_stack+544,int_stack+444,10);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2594,int_stack+694,int_stack+544,10);
 /*--- compute (f0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+3044,int_stack+2594,int_stack+2294,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+904,int_stack+694,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+4622,int_stack+0,int_stack+2594,10);
 /*--- compute (f0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+0,int_stack+4622,int_stack+3044,10);
 /*--- compute (dp|ff) ---*/
 hrr1_build_dp(Libint->AB,int_stack+4622,int_stack+0,int_stack+4022,100);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+6422,int_stack+1334,int_stack+1184,15);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2294,int_stack+1559,int_stack+1334,15);
 /*--- compute (g0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+2969,int_stack+2294,int_stack+6422,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6422,int_stack+1874,int_stack+1559,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7367,int_stack+6422,int_stack+2294,15);
 /*--- compute (g0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+1000,int_stack+7367,int_stack+2969,15);
 /*--- compute (fp|ff) ---*/
 hrr1_build_fp(Libint->AB,int_stack+6422,int_stack+1000,int_stack+0,100);
 /*--- compute (dd|ff) ---*/
 hrr1_build_dd(Libint->AB,int_stack+0,int_stack+6422,int_stack+4622,100);
 return int_stack+0;}
