#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gdgg(Libint_t*, prim_data*);

  /* Computes quartets of (gd|gg) integrals */

REALTYPE *hrr_order_gdgg(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[6][4] = int_stack + 5220;
 Libint->vrr_classes[6][5] = int_stack + 5640;
 Libint->vrr_classes[6][6] = int_stack + 6228;
 Libint->vrr_classes[6][7] = int_stack + 7012;
 Libint->vrr_classes[6][8] = int_stack + 8020;
 memset(int_stack,0,9280*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 9280;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gdgg(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+9280,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9955,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+10900,int_stack+9955,int_stack+9280,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12250,int_stack+960,int_stack+540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13510,int_stack+12250,int_stack+9955,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+15400,int_stack+13510,int_stack+10900,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9280,int_stack+1500,int_stack+960,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+17650,int_stack+9280,int_stack+12250,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+9280,int_stack+17650,int_stack+13510,15);
 /*--- compute (g0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+17650,int_stack+9280,int_stack+15400,15);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+9280,int_stack+2490,int_stack+2175,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+10225,int_stack+2931,int_stack+2490,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+11548,int_stack+10225,int_stack+9280,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+13438,int_stack+3519,int_stack+2931,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+13438,int_stack+10225,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+21025,int_stack+0,int_stack+11548,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9280,int_stack+4275,int_stack+3519,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+24175,int_stack+9280,int_stack+13438,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+9280,int_stack+24175,int_stack+0,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+0,int_stack+9280,int_stack+21025,21);
 /*--- compute (gp|gg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+21025,int_stack+0,int_stack+17650,225);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+9280,int_stack+5640,int_stack+5220,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+10540,int_stack+6228,int_stack+5640,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+12304,int_stack+10540,int_stack+9280,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14824,int_stack+7012,int_stack+6228,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+17176,int_stack+14824,int_stack+10540,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+31150,int_stack+17176,int_stack+12304,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9280,int_stack+8020,int_stack+7012,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+35350,int_stack+9280,int_stack+14824,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+4725,int_stack+35350,int_stack+17176,28);
 /*--- compute (i0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+10605,int_stack+4725,int_stack+31150,28);
 /*--- compute (hp|gg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+31150,int_stack+10605,int_stack+0,225);
 /*--- compute (gd|gg) ---*/
 hrr1_build_gd(Libint->AB,int_stack+0,int_stack+31150,int_stack+21025,225);
 return int_stack+0;}
