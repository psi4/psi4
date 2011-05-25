#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fpgg(Libint_t*, prim_data*);

  /* Computes quartets of (fp|gg) integrals */

REALTYPE *hrr_order_fpgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 Libint->vrr_classes[3][6] = int_stack + 360;
 Libint->vrr_classes[3][7] = int_stack + 640;
 Libint->vrr_classes[3][8] = int_stack + 1000;
 Libint->vrr_classes[4][4] = int_stack + 1450;
 Libint->vrr_classes[4][5] = int_stack + 1675;
 Libint->vrr_classes[4][6] = int_stack + 1990;
 Libint->vrr_classes[4][7] = int_stack + 2410;
 Libint->vrr_classes[4][8] = int_stack + 2950;
 memset(int_stack,0,3625*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3625;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fpgg(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3625,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4075,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+4705,int_stack+4075,int_stack+3625,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+5605,int_stack+640,int_stack+360,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6445,int_stack+5605,int_stack+4075,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+7705,int_stack+6445,int_stack+4705,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+3625,int_stack+1000,int_stack+640,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+9205,int_stack+3625,int_stack+5605,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+3625,int_stack+9205,int_stack+6445,10);
 /*--- compute (f0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+9205,int_stack+3625,int_stack+7705,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3625,int_stack+1675,int_stack+1450,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4300,int_stack+1990,int_stack+1675,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5245,int_stack+4300,int_stack+3625,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6595,int_stack+2410,int_stack+1990,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+6595,int_stack+4300,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+11455,int_stack+0,int_stack+5245,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+3625,int_stack+2950,int_stack+2410,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+13705,int_stack+3625,int_stack+6595,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+1890,int_stack+13705,int_stack+0,15);
 /*--- compute (g0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+5040,int_stack+1890,int_stack+11455,15);
 /*--- compute (fp|gg) ---*/
 hrr1_build_fp(Libint->AB,int_stack+11455,int_stack+5040,int_stack+9205,225);
 return int_stack+11455;}
