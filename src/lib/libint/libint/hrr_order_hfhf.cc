#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hfhf(Libint_t*, prim_data*);

  /* Computes quartets of (hf|hf) integrals */

REALTYPE *hrr_order_hfhf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][5] = int_stack + 0;
 Libint->vrr_classes[5][6] = int_stack + 441;
 Libint->vrr_classes[5][7] = int_stack + 1029;
 Libint->vrr_classes[5][8] = int_stack + 1785;
 Libint->vrr_classes[6][5] = int_stack + 2730;
 Libint->vrr_classes[6][6] = int_stack + 3318;
 Libint->vrr_classes[6][7] = int_stack + 4102;
 Libint->vrr_classes[6][8] = int_stack + 5110;
 Libint->vrr_classes[7][5] = int_stack + 6370;
 Libint->vrr_classes[7][6] = int_stack + 7126;
 Libint->vrr_classes[7][7] = int_stack + 8134;
 Libint->vrr_classes[7][8] = int_stack + 9430;
 Libint->vrr_classes[8][5] = int_stack + 11050;
 Libint->vrr_classes[8][6] = int_stack + 11995;
 Libint->vrr_classes[8][7] = int_stack + 13255;
 Libint->vrr_classes[8][8] = int_stack + 14875;
 memset(int_stack,0,16900*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 16900;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hfhf(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+16900,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+18223,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+19987,int_stack+18223,int_stack+16900,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+22633,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+24901,int_stack+22633,int_stack+18223,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+28429,int_stack+24901,int_stack+19987,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+16900,int_stack+3318,int_stack+2730,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+18664,int_stack+4102,int_stack+3318,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+21016,int_stack+18664,int_stack+16900,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+24544,int_stack+5110,int_stack+4102,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+24544,int_stack+18664,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+32839,int_stack+0,int_stack+21016,28);
 /*--- compute (hp|hf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+38719,int_stack+32839,int_stack+28429,210);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+7126,int_stack+6370,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+8134,int_stack+7126,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+16900,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+21436,int_stack+9430,int_stack+8134,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+25324,int_stack+21436,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+25324,int_stack+16900,36);
 /*--- compute (ip|hf) ---*/
 hrr1_build_ip(Libint->AB,int_stack+51949,int_stack+0,int_stack+32839,210);
 /*--- compute (hd|hf) ---*/
 hrr1_build_hd(Libint->AB,int_stack+69589,int_stack+51949,int_stack+38719,210);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+16900,int_stack+11995,int_stack+11050,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+19735,int_stack+13255,int_stack+11995,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+23515,int_stack+19735,int_stack+16900,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+29185,int_stack+14875,int_stack+13255,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+34045,int_stack+29185,int_stack+19735,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+41605,int_stack+34045,int_stack+23515,45);
 /*--- compute (kp|hf) ---*/
 hrr1_build_kp(Libint->AB,int_stack+7560,int_stack+41605,int_stack+0,210);
 /*--- compute (id|hf) ---*/
 hrr1_build_id(Libint->AB,int_stack+96049,int_stack+7560,int_stack+51949,210);
 /*--- compute (hf|hf) ---*/
 hrr1_build_hf(Libint->AB,int_stack+0,int_stack+96049,int_stack+69589,210);
 return int_stack+0;}
