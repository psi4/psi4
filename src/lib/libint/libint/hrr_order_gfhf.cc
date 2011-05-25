#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gfhf(Libint_t*, prim_data*);

  /* Computes quartets of (gf|hf) integrals */

REALTYPE *hrr_order_gfhf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[4][8] = int_stack + 1275;
 Libint->vrr_classes[5][5] = int_stack + 1950;
 Libint->vrr_classes[5][6] = int_stack + 2391;
 Libint->vrr_classes[5][7] = int_stack + 2979;
 Libint->vrr_classes[5][8] = int_stack + 3735;
 Libint->vrr_classes[6][5] = int_stack + 4680;
 Libint->vrr_classes[6][6] = int_stack + 5268;
 Libint->vrr_classes[6][7] = int_stack + 6052;
 Libint->vrr_classes[6][8] = int_stack + 7060;
 Libint->vrr_classes[7][5] = int_stack + 8320;
 Libint->vrr_classes[7][6] = int_stack + 9076;
 Libint->vrr_classes[7][7] = int_stack + 10084;
 Libint->vrr_classes[7][8] = int_stack + 11380;
 memset(int_stack,0,13000*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 13000;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gfhf(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+13000,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+13945,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+15205,int_stack+13945,int_stack+13000,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+17095,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+18715,int_stack+17095,int_stack+13945,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+21235,int_stack+18715,int_stack+15205,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+13000,int_stack+2391,int_stack+1950,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14323,int_stack+2979,int_stack+2391,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+16087,int_stack+14323,int_stack+13000,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18733,int_stack+3735,int_stack+2979,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+18733,int_stack+14323,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+24385,int_stack+0,int_stack+16087,21);
 /*--- compute (gp|hf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+28795,int_stack+24385,int_stack+21235,210);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+5268,int_stack+4680,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+6052,int_stack+5268,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13000,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+16528,int_stack+7060,int_stack+6052,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+19552,int_stack+16528,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+19552,int_stack+13000,28);
 /*--- compute (hp|hf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+38245,int_stack+0,int_stack+24385,210);
 /*--- compute (gd|hf) ---*/
 hrr1_build_gd(Libint->AB,int_stack+51475,int_stack+38245,int_stack+28795,210);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+13000,int_stack+9076,int_stack+8320,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+15268,int_stack+10084,int_stack+9076,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+18292,int_stack+15268,int_stack+13000,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+22828,int_stack+11380,int_stack+10084,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+26716,int_stack+22828,int_stack+15268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+5880,int_stack+26716,int_stack+18292,36);
 /*--- compute (ip|hf) ---*/
 hrr1_build_ip(Libint->AB,int_stack+13440,int_stack+5880,int_stack+0,210);
 /*--- compute (hd|hf) ---*/
 hrr1_build_hd(Libint->AB,int_stack+70375,int_stack+13440,int_stack+38245,210);
 /*--- compute (gf|hf) ---*/
 hrr1_build_gf(Libint->AB,int_stack+0,int_stack+70375,int_stack+51475,210);
 return int_stack+0;}
