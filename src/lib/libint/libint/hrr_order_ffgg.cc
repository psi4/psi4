#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ffgg(Libint_t*, prim_data*);

  /* Computes quartets of (ff|gg) integrals */

REALTYPE *hrr_order_ffgg(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[6][4] = int_stack + 6670;
 Libint->vrr_classes[6][5] = int_stack + 7090;
 Libint->vrr_classes[6][6] = int_stack + 7678;
 Libint->vrr_classes[6][7] = int_stack + 8462;
 Libint->vrr_classes[6][8] = int_stack + 9470;
 memset(int_stack,0,10730*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 10730;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ffgg(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+10730,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+11180,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+11810,int_stack+11180,int_stack+10730,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12710,int_stack+640,int_stack+360,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13550,int_stack+12710,int_stack+11180,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+14810,int_stack+13550,int_stack+11810,10);
 /*--- compute (f0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+10730,int_stack+1000,int_stack+640,10);
 /*--- compute (f0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+16310,int_stack+10730,int_stack+12710,10);
 /*--- compute (f0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+10730,int_stack+16310,int_stack+13550,10);
 /*--- compute (f0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+16310,int_stack+10730,int_stack+14810,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+10730,int_stack+1675,int_stack+1450,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+11405,int_stack+1990,int_stack+1675,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+12350,int_stack+11405,int_stack+10730,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+13700,int_stack+2410,int_stack+1990,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+13700,int_stack+11405,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+18560,int_stack+0,int_stack+12350,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+10730,int_stack+2950,int_stack+2410,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+20810,int_stack+10730,int_stack+13700,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+10730,int_stack+20810,int_stack+0,15);
 /*--- compute (g0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+0,int_stack+10730,int_stack+18560,15);
 /*--- compute (fp|gg) ---*/
 hrr1_build_fp(Libint->AB,int_stack+18560,int_stack+0,int_stack+16310,225);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+10730,int_stack+3940,int_stack+3625,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+11675,int_stack+4381,int_stack+3940,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+12998,int_stack+11675,int_stack+10730,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14888,int_stack+4969,int_stack+4381,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+25310,int_stack+14888,int_stack+11675,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+27956,int_stack+25310,int_stack+12998,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+10730,int_stack+5725,int_stack+4969,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+31106,int_stack+10730,int_stack+14888,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+10730,int_stack+31106,int_stack+25310,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+31106,int_stack+10730,int_stack+27956,21);
 /*--- compute (gp|gg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+35831,int_stack+31106,int_stack+0,225);
 /*--- compute (fd|gg) ---*/
 hrr1_build_fd(Libint->AB,int_stack+45956,int_stack+35831,int_stack+18560,225);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+0,int_stack+7090,int_stack+6670,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1260,int_stack+7678,int_stack+7090,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3024,int_stack+1260,int_stack+0,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10730,int_stack+8462,int_stack+7678,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13082,int_stack+10730,int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+16610,int_stack+13082,int_stack+3024,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+9470,int_stack+8462,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+3024,int_stack+0,int_stack+10730,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+20810,int_stack+3024,int_stack+13082,28);
 /*--- compute (i0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+0,int_stack+20810,int_stack+16610,28);
 /*--- compute (hp|gg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+6300,int_stack+0,int_stack+31106,225);
 /*--- compute (gd|gg) ---*/
 hrr1_build_gd(Libint->AB,int_stack+59456,int_stack+6300,int_stack+35831,225);
 /*--- compute (ff|gg) ---*/
 hrr1_build_ff(Libint->AB,int_stack+0,int_stack+59456,int_stack+45956,225);
 return int_stack+0;}
