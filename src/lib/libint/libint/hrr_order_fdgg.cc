#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdgg(Libint_t*, prim_data*);

  /* Computes quartets of (fd|gg) integrals */

REALTYPE *hrr_order_fdgg(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[5][4] = int_stack + 3625;
 Libint->vrr_classes[5][5] = int_stack + 3940;
 Libint->vrr_classes[5][6] = int_stack + 4381;
 Libint->vrr_classes[5][7] = int_stack + 4969;
 Libint->vrr_classes[5][8] = int_stack + 5725;
 memset(int_stack,0,6670*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 6670;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdgg(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+6670,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7120,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7750,int_stack+7120,int_stack+6670,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+8650,int_stack+640,int_stack+360,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9490,int_stack+8650,int_stack+7120,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+10750,int_stack+9490,int_stack+7750,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+6670,int_stack+1000,int_stack+640,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+12250,int_stack+6670,int_stack+8650,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+6670,int_stack+12250,int_stack+9490,10);
 /*--- compute (f0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+12250,int_stack+6670,int_stack+10750,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+6670,int_stack+1675,int_stack+1450,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7345,int_stack+1990,int_stack+1675,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+8290,int_stack+7345,int_stack+6670,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9640,int_stack+2410,int_stack+1990,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+9640,int_stack+7345,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+14500,int_stack+0,int_stack+8290,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+6670,int_stack+2950,int_stack+2410,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+16750,int_stack+6670,int_stack+9640,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+6670,int_stack+16750,int_stack+0,15);
 /*--- compute (g0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+0,int_stack+6670,int_stack+14500,15);
 /*--- compute (fp|gg) ---*/
 hrr1_build_fp(Libint->AB,int_stack+14500,int_stack+0,int_stack+12250,225);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+6670,int_stack+3940,int_stack+3625,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7615,int_stack+4381,int_stack+3940,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+8938,int_stack+7615,int_stack+6670,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10828,int_stack+4969,int_stack+4381,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+21250,int_stack+10828,int_stack+7615,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+23896,int_stack+21250,int_stack+8938,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+6670,int_stack+5725,int_stack+4969,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+27046,int_stack+6670,int_stack+10828,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+3375,int_stack+27046,int_stack+21250,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+7785,int_stack+3375,int_stack+23896,21);
 /*--- compute (gp|gg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+21250,int_stack+7785,int_stack+0,225);
 /*--- compute (fd|gg) ---*/
 hrr1_build_fd(Libint->AB,int_stack+0,int_stack+21250,int_stack+14500,225);
 return int_stack+0;}
