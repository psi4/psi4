#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dpgf(Libint_t*, prim_data*);

  /* Computes quartets of (dp|gf) integrals */

REALTYPE *hrr_order_dpgf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][4] = int_stack + 0;
 Libint->vrr_classes[2][5] = int_stack + 90;
 Libint->vrr_classes[2][6] = int_stack + 216;
 Libint->vrr_classes[2][7] = int_stack + 384;
 Libint->vrr_classes[3][4] = int_stack + 600;
 Libint->vrr_classes[3][5] = int_stack + 750;
 Libint->vrr_classes[3][6] = int_stack + 960;
 Libint->vrr_classes[3][7] = int_stack + 1240;
 memset(int_stack,0,1600*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1600;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dpgf(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1600,int_stack+90,int_stack+0,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1870,int_stack+216,int_stack+90,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2248,int_stack+1870,int_stack+1600,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2788,int_stack+384,int_stack+216,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3292,int_stack+2788,int_stack+1870,6);
 /*--- compute (d0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+4048,int_stack+3292,int_stack+2248,6);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1600,int_stack+750,int_stack+600,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+2050,int_stack+960,int_stack+750,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2680,int_stack+2050,int_stack+1600,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+0,int_stack+1240,int_stack+960,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4948,int_stack+0,int_stack+2050,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+0,int_stack+4948,int_stack+2680,10);
 /*--- compute (dp|gf) ---*/
 hrr1_build_dp(Libint->AB,int_stack+4948,int_stack+0,int_stack+4048,150);
 return int_stack+4948;}
