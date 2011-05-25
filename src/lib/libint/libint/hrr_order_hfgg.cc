#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hfgg(Libint_t*, prim_data*);

  /* Computes quartets of (hf|gg) integrals */

REALTYPE *hrr_order_hfgg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 Libint->vrr_classes[5][6] = int_stack + 756;
 Libint->vrr_classes[5][7] = int_stack + 1344;
 Libint->vrr_classes[5][8] = int_stack + 2100;
 Libint->vrr_classes[6][4] = int_stack + 3045;
 Libint->vrr_classes[6][5] = int_stack + 3465;
 Libint->vrr_classes[6][6] = int_stack + 4053;
 Libint->vrr_classes[6][7] = int_stack + 4837;
 Libint->vrr_classes[6][8] = int_stack + 5845;
 Libint->vrr_classes[7][4] = int_stack + 7105;
 Libint->vrr_classes[7][5] = int_stack + 7645;
 Libint->vrr_classes[7][6] = int_stack + 8401;
 Libint->vrr_classes[7][7] = int_stack + 9409;
 Libint->vrr_classes[7][8] = int_stack + 10705;
 Libint->vrr_classes[8][4] = int_stack + 12325;
 Libint->vrr_classes[8][5] = int_stack + 13000;
 Libint->vrr_classes[8][6] = int_stack + 13945;
 Libint->vrr_classes[8][7] = int_stack + 15205;
 Libint->vrr_classes[8][8] = int_stack + 16825;
 memset(int_stack,0,18850*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 18850;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hfgg(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+18850,int_stack+315,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+19795,int_stack+756,int_stack+315,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+21118,int_stack+19795,int_stack+18850,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+23008,int_stack+1344,int_stack+756,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+24772,int_stack+23008,int_stack+19795,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+27418,int_stack+24772,int_stack+21118,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18850,int_stack+2100,int_stack+1344,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+30568,int_stack+18850,int_stack+23008,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+18850,int_stack+30568,int_stack+24772,21);
 /*--- compute (h0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+30568,int_stack+18850,int_stack+27418,21);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+18850,int_stack+3465,int_stack+3045,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+20110,int_stack+4053,int_stack+3465,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+21874,int_stack+20110,int_stack+18850,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+24394,int_stack+4837,int_stack+4053,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+26746,int_stack+24394,int_stack+20110,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+0,int_stack+26746,int_stack+21874,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18850,int_stack+5845,int_stack+4837,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+35293,int_stack+18850,int_stack+24394,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+18850,int_stack+35293,int_stack+26746,28);
 /*--- compute (i0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+35293,int_stack+18850,int_stack+0,28);
 /*--- compute (hp|gg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+41593,int_stack+35293,int_stack+30568,225);
 /*--- compute (k0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+0,int_stack+7645,int_stack+7105,36);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1620,int_stack+8401,int_stack+7645,36);
 /*--- compute (k0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3888,int_stack+1620,int_stack+0,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+18850,int_stack+9409,int_stack+8401,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+21874,int_stack+18850,int_stack+1620,36);
 /*--- compute (k0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+26410,int_stack+21874,int_stack+3888,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+10705,int_stack+9409,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+3888,int_stack+0,int_stack+18850,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+55768,int_stack+3888,int_stack+21874,36);
 /*--- compute (k0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+0,int_stack+55768,int_stack+26410,36);
 /*--- compute (ip|gg) ---*/
 hrr1_build_ip(Libint->AB,int_stack+55768,int_stack+0,int_stack+35293,225);
 /*--- compute (hd|gg) ---*/
 hrr1_build_hd(Libint->AB,int_stack+74668,int_stack+55768,int_stack+41593,225);
 /*--- compute (l0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+18850,int_stack+13000,int_stack+12325,45);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+20875,int_stack+13945,int_stack+13000,45);
 /*--- compute (l0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+23710,int_stack+20875,int_stack+18850,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+27760,int_stack+15205,int_stack+13945,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+31540,int_stack+27760,int_stack+20875,45);
 /*--- compute (l0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+37210,int_stack+31540,int_stack+23710,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18850,int_stack+16825,int_stack+15205,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+43960,int_stack+18850,int_stack+27760,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+8100,int_stack+43960,int_stack+31540,45);
 /*--- compute (l0|gg) ---*/
 hrr3_build_gg(Libint->CD,int_stack+43960,int_stack+8100,int_stack+37210,45);
 /*--- compute (kp|gg) ---*/
 hrr1_build_kp(Libint->AB,int_stack+8100,int_stack+43960,int_stack+0,225);
 /*--- compute (id|gg) ---*/
 hrr1_build_id(Libint->AB,int_stack+103018,int_stack+8100,int_stack+55768,225);
 /*--- compute (hf|gg) ---*/
 hrr1_build_hf(Libint->AB,int_stack+0,int_stack+103018,int_stack+74668,225);
 return int_stack+0;}
