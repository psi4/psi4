#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hdhf(Libint_t*, prim_data*);

  /* Computes quartets of (hd|hf) integrals */

REALTYPE *hrr_order_hdhf(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,11050*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 11050;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hdhf(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+11050,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12373,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+14137,int_stack+12373,int_stack+11050,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+16783,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+19051,int_stack+16783,int_stack+12373,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+22579,int_stack+19051,int_stack+14137,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+11050,int_stack+3318,int_stack+2730,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12814,int_stack+4102,int_stack+3318,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+15166,int_stack+12814,int_stack+11050,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18694,int_stack+5110,int_stack+4102,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+18694,int_stack+12814,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+26989,int_stack+0,int_stack+15166,28);
 /*--- compute (hp|hf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+32869,int_stack+26989,int_stack+22579,210);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+7126,int_stack+6370,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+8134,int_stack+7126,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+11050,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+15586,int_stack+9430,int_stack+8134,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+19474,int_stack+15586,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+19474,int_stack+11050,36);
 /*--- compute (ip|hf) ---*/
 hrr1_build_ip(Libint->AB,int_stack+7560,int_stack+0,int_stack+26989,210);
 /*--- compute (hd|hf) ---*/
 hrr1_build_hd(Libint->AB,int_stack+46099,int_stack+7560,int_stack+32869,210);
 return int_stack+46099;}
