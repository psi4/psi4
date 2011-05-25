#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ddgg(Libint_t*, prim_data*);

  /* Computes quartets of (dd|gg) integrals */

REALTYPE *hrr_order_ddgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][4] = int_stack + 0;
 Libint->vrr_classes[2][5] = int_stack + 90;
 Libint->vrr_classes[2][6] = int_stack + 216;
 Libint->vrr_classes[2][7] = int_stack + 384;
 Libint->vrr_classes[2][8] = int_stack + 600;
 Libint->vrr_classes[3][4] = int_stack + 870;
 Libint->vrr_classes[3][5] = int_stack + 1020;
 Libint->vrr_classes[3][6] = int_stack + 1230;
 Libint->vrr_classes[3][7] = int_stack + 1510;
 Libint->vrr_classes[3][8] = int_stack + 1870;
 Libint->vrr_classes[4][4] = int_stack + 2320;
 Libint->vrr_classes[4][5] = int_stack + 2545;
 Libint->vrr_classes[4][6] = int_stack + 2860;
 Libint->vrr_classes[4][7] = int_stack + 3280;
 Libint->vrr_classes[4][8] = int_stack + 3820;
 memset(int_stack,0,4495*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4495;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ddgg(Libint, Data);
   Data++;
 }
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4495,int_stack+90,int_stack+0,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4765,int_stack+216,int_stack+90,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5143,int_stack+4765,int_stack+4495,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+5683,int_stack+384,int_stack+216,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6187,int_stack+5683,int_stack+4765,6);
 /*--- compute (d0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+6943,int_stack+6187,int_stack+5143,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+4495,int_stack+600,int_stack+384,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+7843,int_stack+4495,int_stack+5683,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+4495,int_stack+7843,int_stack+6187,6);
 /*--- compute (d0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+7843,int_stack+4495,int_stack+6943,6);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4495,int_stack+1020,int_stack+870,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4945,int_stack+1230,int_stack+1020,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5575,int_stack+4945,int_stack+4495,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6475,int_stack+1510,int_stack+1230,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+6475,int_stack+4945,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+9193,int_stack+0,int_stack+5575,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+4495,int_stack+1870,int_stack+1510,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+10693,int_stack+4495,int_stack+6475,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+4495,int_stack+10693,int_stack+0,10);
 /*--- compute (f0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+0,int_stack+4495,int_stack+9193,10);
 /*--- compute (dp|gg) ---*/
 hrr1_build_dp(Libint->AB,int_stack+9193,int_stack+0,int_stack+7843,225);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+13243,int_stack+2545,int_stack+2320,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4495,int_stack+2860,int_stack+2545,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5440,int_stack+4495,int_stack+13243,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6790,int_stack+3280,int_stack+2860,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13243,int_stack+6790,int_stack+4495,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+15133,int_stack+13243,int_stack+5440,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+4495,int_stack+3820,int_stack+3280,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+17383,int_stack+4495,int_stack+6790,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+2250,int_stack+17383,int_stack+13243,15);
 /*--- compute (g0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+5400,int_stack+2250,int_stack+15133,15);
 /*--- compute (fp|gg) ---*/
 hrr1_build_fp(Libint->AB,int_stack+13243,int_stack+5400,int_stack+0,225);
 /*--- compute (dd|gg) ---*/
 hrr1_build_dd(Libint->AB,int_stack+0,int_stack+13243,int_stack+9193,225);
 return int_stack+0;}
