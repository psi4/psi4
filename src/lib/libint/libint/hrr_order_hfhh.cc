#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hfhh(Libint_t*, prim_data*);

  /* Computes quartets of (hf|hh) integrals */

REALTYPE *hrr_order_hfhh(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][5] = int_stack + 0;
 Libint->vrr_classes[5][6] = int_stack + 441;
 Libint->vrr_classes[5][7] = int_stack + 1029;
 Libint->vrr_classes[5][8] = int_stack + 1785;
 Libint->vrr_classes[5][9] = int_stack + 2730;
 Libint->vrr_classes[5][10] = int_stack + 3885;
 Libint->vrr_classes[6][5] = int_stack + 5271;
 Libint->vrr_classes[6][6] = int_stack + 5859;
 Libint->vrr_classes[6][7] = int_stack + 6643;
 Libint->vrr_classes[6][8] = int_stack + 7651;
 Libint->vrr_classes[6][9] = int_stack + 8911;
 Libint->vrr_classes[6][10] = int_stack + 10451;
 Libint->vrr_classes[7][5] = int_stack + 12299;
 Libint->vrr_classes[7][6] = int_stack + 13055;
 Libint->vrr_classes[7][7] = int_stack + 14063;
 Libint->vrr_classes[7][8] = int_stack + 15359;
 Libint->vrr_classes[7][9] = int_stack + 16979;
 Libint->vrr_classes[7][10] = int_stack + 18959;
 Libint->vrr_classes[8][5] = int_stack + 21335;
 Libint->vrr_classes[8][6] = int_stack + 22280;
 Libint->vrr_classes[8][7] = int_stack + 23540;
 Libint->vrr_classes[8][8] = int_stack + 25160;
 Libint->vrr_classes[8][9] = int_stack + 27185;
 Libint->vrr_classes[8][10] = int_stack + 29660;
 memset(int_stack,0,32630*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 32630;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hfhh(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+32630,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+33953,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+35717,int_stack+33953,int_stack+32630,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+38363,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+40631,int_stack+38363,int_stack+33953,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+44159,int_stack+40631,int_stack+35717,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+32630,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+48569,int_stack+32630,int_stack+38363,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+53105,int_stack+48569,int_stack+40631,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+35465,int_stack+53105,int_stack+44159,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+42080,int_stack+3885,int_stack+2730,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+58985,int_stack+42080,int_stack+32630,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+64655,int_stack+58985,int_stack+48569,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+42080,int_stack+64655,int_stack+53105,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+50900,int_stack+42080,int_stack+35465,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+32630,int_stack+5859,int_stack+5271,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+34394,int_stack+6643,int_stack+5859,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+36746,int_stack+34394,int_stack+32630,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+40274,int_stack+7651,int_stack+6643,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+43298,int_stack+40274,int_stack+34394,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+60161,int_stack+43298,int_stack+36746,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+32630,int_stack+8911,int_stack+7651,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+66041,int_stack+32630,int_stack+40274,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+66041,int_stack+43298,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+36410,int_stack+0,int_stack+60161,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+60161,int_stack+10451,int_stack+8911,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+72089,int_stack+60161,int_stack+32630,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+79649,int_stack+72089,int_stack+66041,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+60161,int_stack+79649,int_stack+0,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+71921,int_stack+60161,int_stack+36410,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+84269,int_stack+71921,int_stack+50900,441);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+13055,int_stack+12299,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+14063,int_stack+13055,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9828,int_stack+15359,int_stack+14063,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+32630,int_stack+9828,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+38678,int_stack+32630,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+16979,int_stack+15359,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+46238,int_stack+0,int_stack+9828,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+4860,int_stack+46238,int_stack+32630,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+54014,int_stack+4860,int_stack+38678,36);
 /*--- compute (k0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+32630,int_stack+18959,int_stack+16979,36);
 /*--- compute (k0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+112052,int_stack+32630,int_stack+0,36);
 /*--- compute (k0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+32630,int_stack+112052,int_stack+46238,36);
 /*--- compute (k0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+112052,int_stack+32630,int_stack+4860,36);
 /*--- compute (k0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+32630,int_stack+112052,int_stack+54014,36);
 /*--- compute (ip|hh) ---*/
 hrr1_build_ip(Libint->AB,int_stack+112052,int_stack+32630,int_stack+71921,441);
 /*--- compute (hd|hh) ---*/
 hrr1_build_hd(Libint->AB,int_stack+149096,int_stack+112052,int_stack+84269,441);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+48506,int_stack+22280,int_stack+21335,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+51341,int_stack+23540,int_stack+22280,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+55121,int_stack+51341,int_stack+48506,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+60791,int_stack+25160,int_stack+23540,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+65651,int_stack+60791,int_stack+51341,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+73211,int_stack+65651,int_stack+55121,45);
 /*--- compute (l0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+48506,int_stack+27185,int_stack+25160,45);
 /*--- compute (l0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+82661,int_stack+48506,int_stack+60791,45);
 /*--- compute (l0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+92381,int_stack+82661,int_stack+65651,45);
 /*--- compute (l0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+54581,int_stack+92381,int_stack+73211,45);
 /*--- compute (l0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+68756,int_stack+29660,int_stack+27185,45);
 /*--- compute (l0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+0,int_stack+68756,int_stack+48506,45);
 /*--- compute (l0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+12150,int_stack+0,int_stack+82661,45);
 /*--- compute (l0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+68756,int_stack+12150,int_stack+92381,45);
 /*--- compute (l0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+68756,int_stack+54581,45);
 /*--- compute (kp|hh) ---*/
 hrr1_build_kp(Libint->AB,int_stack+48506,int_stack+0,int_stack+32630,441);
 /*--- compute (id|hh) ---*/
 hrr1_build_id(Libint->AB,int_stack+204662,int_stack+48506,int_stack+112052,441);
 /*--- compute (hf|hh) ---*/
 hrr1_build_hf(Libint->AB,int_stack+0,int_stack+204662,int_stack+149096,441);
 return int_stack+0;}
