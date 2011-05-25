#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gpgg(Libint_t*, prim_data*);

  /* Computes quartets of (gp|gg) integrals */

REALTYPE *hrr_order_gpgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][4] = int_stack + 0;
 Libint->vrr_classes[4][5] = int_stack + 225;
 Libint->vrr_classes[4][6] = int_stack + 540;
 Libint->vrr_classes[4][7] = int_stack + 960;
 Libint->vrr_classes[4][8] = int_stack + 1500;
 Libint->vrr_classes[5][4] = int_stack + 2175;
 Libint->vrr_classes[5][5] = int_stack + 2490;
 Libint->vrr_classes[5][6] = int_stack + 2931;
 Libint->vrr_classes[5][7] = int_stack + 3519;
 Libint->vrr_classes[5][8] = int_stack + 4275;
 memset(int_stack,0,5220*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 5220;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gpgg(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+5220,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5895,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+6840,int_stack+5895,int_stack+5220,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+8190,int_stack+960,int_stack+540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9450,int_stack+8190,int_stack+5895,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+11340,int_stack+9450,int_stack+6840,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+5220,int_stack+1500,int_stack+960,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+13590,int_stack+5220,int_stack+8190,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+5220,int_stack+13590,int_stack+9450,15);
 /*--- compute (g0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+13590,int_stack+5220,int_stack+11340,15);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+5220,int_stack+2490,int_stack+2175,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6165,int_stack+2931,int_stack+2490,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7488,int_stack+6165,int_stack+5220,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9378,int_stack+3519,int_stack+2931,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+9378,int_stack+6165,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+16965,int_stack+0,int_stack+7488,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+5220,int_stack+4275,int_stack+3519,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+20115,int_stack+5220,int_stack+9378,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+2646,int_stack+20115,int_stack+0,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+7056,int_stack+2646,int_stack+16965,21);
 /*--- compute (gp|gg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+16965,int_stack+7056,int_stack+13590,225);
 return int_stack+16965;}
