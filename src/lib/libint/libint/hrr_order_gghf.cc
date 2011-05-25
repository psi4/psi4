#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gghf(Libint_t*, prim_data*);

  /* Computes quartets of (gg|hf) integrals */

REALTYPE *hrr_order_gghf(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[8][5] = int_stack + 13000;
 Libint->vrr_classes[8][6] = int_stack + 13945;
 Libint->vrr_classes[8][7] = int_stack + 15205;
 Libint->vrr_classes[8][8] = int_stack + 16825;
 memset(int_stack,0,18850*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 18850;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gghf(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+18850,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+19795,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+21055,int_stack+19795,int_stack+18850,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+22945,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+24565,int_stack+22945,int_stack+19795,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+27085,int_stack+24565,int_stack+21055,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+18850,int_stack+2391,int_stack+1950,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+20173,int_stack+2979,int_stack+2391,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+21937,int_stack+20173,int_stack+18850,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+24583,int_stack+3735,int_stack+2979,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+24583,int_stack+20173,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+30235,int_stack+0,int_stack+21937,21);
 /*--- compute (gp|hf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+34645,int_stack+30235,int_stack+27085,210);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+5268,int_stack+4680,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+6052,int_stack+5268,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+18850,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+22378,int_stack+7060,int_stack+6052,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+25402,int_stack+22378,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+25402,int_stack+18850,28);
 /*--- compute (hp|hf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+44095,int_stack+0,int_stack+30235,210);
 /*--- compute (gd|hf) ---*/
 hrr1_build_gd(Libint->AB,int_stack+57325,int_stack+44095,int_stack+34645,210);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+18850,int_stack+9076,int_stack+8320,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+21118,int_stack+10084,int_stack+9076,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+24142,int_stack+21118,int_stack+18850,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+28678,int_stack+11380,int_stack+10084,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+32566,int_stack+28678,int_stack+21118,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+76225,int_stack+32566,int_stack+24142,36);
 /*--- compute (ip|hf) ---*/
 hrr1_build_ip(Libint->AB,int_stack+18850,int_stack+76225,int_stack+0,210);
 /*--- compute (hd|hf) ---*/
 hrr1_build_hd(Libint->AB,int_stack+83785,int_stack+18850,int_stack+44095,210);
 /*--- compute (gf|hf) ---*/
 hrr1_build_gf(Libint->AB,int_stack+110245,int_stack+83785,int_stack+57325,210);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+13945,int_stack+13000,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2835,int_stack+15205,int_stack+13945,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+6615,int_stack+2835,int_stack+0,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+36490,int_stack+16825,int_stack+15205,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+41350,int_stack+36490,int_stack+2835,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+48910,int_stack+41350,int_stack+6615,45);
 /*--- compute (kp|hf) ---*/
 hrr1_build_kp(Libint->AB,int_stack+141745,int_stack+48910,int_stack+76225,210);
 /*--- compute (id|hf) ---*/
 hrr1_build_id(Libint->AB,int_stack+36490,int_stack+141745,int_stack+18850,210);
 /*--- compute (hf|hf) ---*/
 hrr1_build_hf(Libint->AB,int_stack+141745,int_stack+36490,int_stack+83785,210);
 /*--- compute (gg|hf) ---*/
 hrr1_build_gg(Libint->AB,int_stack+0,int_stack+141745,int_stack+110245,210);
 return int_stack+0;}
