#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddgf(Libint_t*, prim_data*);

  /* Computes quartets of (dd|gf) integrals */

REALTYPE *hrr_order_ddgf(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[4][4] = int_stack + 1600;
 Libint->vrr_classes[4][5] = int_stack + 1825;
 Libint->vrr_classes[4][6] = int_stack + 2140;
 Libint->vrr_classes[4][7] = int_stack + 2560;
 memset(int_stack,0,3100*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3100;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddgf(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3100,int_stack+90,int_stack+0,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3370,int_stack+216,int_stack+90,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3748,int_stack+3370,int_stack+3100,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4288,int_stack+384,int_stack+216,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4792,int_stack+4288,int_stack+3370,6);
 /*--- compute (d0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+5548,int_stack+4792,int_stack+3748,6);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3100,int_stack+750,int_stack+600,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3550,int_stack+960,int_stack+750,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+4180,int_stack+3550,int_stack+3100,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+0,int_stack+1240,int_stack+960,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6448,int_stack+0,int_stack+3550,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+0,int_stack+6448,int_stack+4180,10);
 /*--- compute (dp|gf) ---*/
 hrr1_build_dp(Libint->AB,int_stack+6448,int_stack+0,int_stack+5548,150);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+9148,int_stack+1825,int_stack+1600,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3100,int_stack+2140,int_stack+1825,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+4045,int_stack+3100,int_stack+9148,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9148,int_stack+2560,int_stack+2140,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10408,int_stack+9148,int_stack+3100,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+1500,int_stack+10408,int_stack+4045,15);
 /*--- compute (fp|gf) ---*/
 hrr1_build_fp(Libint->AB,int_stack+9148,int_stack+1500,int_stack+0,150);
 /*--- compute (dd|gf) ---*/
 hrr1_build_dd(Libint->AB,int_stack+0,int_stack+9148,int_stack+6448,150);
 return int_stack+0;}
