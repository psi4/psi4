#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hhhh(Libint_t*, prim_data*);

  /* Computes quartets of (hh|hh) integrals */

REALTYPE *hrr_order_hhhh(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[9][5] = int_stack + 32630;
 Libint->vrr_classes[9][6] = int_stack + 33785;
 Libint->vrr_classes[9][7] = int_stack + 35325;
 Libint->vrr_classes[9][8] = int_stack + 37305;
 Libint->vrr_classes[9][9] = int_stack + 39780;
 Libint->vrr_classes[9][10] = int_stack + 42805;
 Libint->vrr_classes[10][5] = int_stack + 46435;
 Libint->vrr_classes[10][6] = int_stack + 47821;
 Libint->vrr_classes[10][7] = int_stack + 49669;
 Libint->vrr_classes[10][8] = int_stack + 52045;
 Libint->vrr_classes[10][9] = int_stack + 55015;
 Libint->vrr_classes[10][10] = int_stack + 58645;
 memset(int_stack,0,63001*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 63001;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hhhh(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+63001,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+64324,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+66088,int_stack+64324,int_stack+63001,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+68734,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+71002,int_stack+68734,int_stack+64324,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+74530,int_stack+71002,int_stack+66088,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+63001,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+78940,int_stack+63001,int_stack+68734,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+83476,int_stack+78940,int_stack+71002,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+65836,int_stack+83476,int_stack+74530,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+72451,int_stack+3885,int_stack+2730,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+89356,int_stack+72451,int_stack+63001,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+95026,int_stack+89356,int_stack+78940,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+72451,int_stack+95026,int_stack+83476,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+81271,int_stack+72451,int_stack+65836,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+63001,int_stack+5859,int_stack+5271,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+64765,int_stack+6643,int_stack+5859,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+67117,int_stack+64765,int_stack+63001,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+70645,int_stack+7651,int_stack+6643,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+73669,int_stack+70645,int_stack+64765,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+90532,int_stack+73669,int_stack+67117,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+63001,int_stack+8911,int_stack+7651,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+96412,int_stack+63001,int_stack+70645,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+96412,int_stack+73669,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+66781,int_stack+0,int_stack+90532,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+90532,int_stack+10451,int_stack+8911,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+102460,int_stack+90532,int_stack+63001,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+110020,int_stack+102460,int_stack+96412,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+90532,int_stack+110020,int_stack+0,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+102292,int_stack+90532,int_stack+66781,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+114640,int_stack+102292,int_stack+81271,441);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+13055,int_stack+12299,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+14063,int_stack+13055,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9828,int_stack+15359,int_stack+14063,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+63001,int_stack+9828,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+69049,int_stack+63001,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+16979,int_stack+15359,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+76609,int_stack+0,int_stack+9828,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+4860,int_stack+76609,int_stack+63001,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+84385,int_stack+4860,int_stack+69049,36);
 /*--- compute (k0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+63001,int_stack+18959,int_stack+16979,36);
 /*--- compute (k0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+142423,int_stack+63001,int_stack+0,36);
 /*--- compute (k0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+63001,int_stack+142423,int_stack+76609,36);
 /*--- compute (k0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+142423,int_stack+63001,int_stack+4860,36);
 /*--- compute (k0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+63001,int_stack+142423,int_stack+84385,36);
 /*--- compute (ip|hh) ---*/
 hrr1_build_ip(Libint->AB,int_stack+142423,int_stack+63001,int_stack+102292,441);
 /*--- compute (hd|hh) ---*/
 hrr1_build_hd(Libint->AB,int_stack+179467,int_stack+142423,int_stack+114640,441);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+78877,int_stack+22280,int_stack+21335,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+81712,int_stack+23540,int_stack+22280,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+85492,int_stack+81712,int_stack+78877,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+91162,int_stack+25160,int_stack+23540,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+96022,int_stack+91162,int_stack+81712,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+103582,int_stack+96022,int_stack+85492,45);
 /*--- compute (l0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+78877,int_stack+27185,int_stack+25160,45);
 /*--- compute (l0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+113032,int_stack+78877,int_stack+91162,45);
 /*--- compute (l0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+122752,int_stack+113032,int_stack+96022,45);
 /*--- compute (l0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+84952,int_stack+122752,int_stack+103582,45);
 /*--- compute (l0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+99127,int_stack+29660,int_stack+27185,45);
 /*--- compute (l0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+0,int_stack+99127,int_stack+78877,45);
 /*--- compute (l0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+12150,int_stack+0,int_stack+113032,45);
 /*--- compute (l0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+99127,int_stack+12150,int_stack+122752,45);
 /*--- compute (l0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+99127,int_stack+84952,45);
 /*--- compute (kp|hh) ---*/
 hrr1_build_kp(Libint->AB,int_stack+78877,int_stack+0,int_stack+63001,441);
 /*--- compute (id|hh) ---*/
 hrr1_build_id(Libint->AB,int_stack+235033,int_stack+78877,int_stack+142423,441);
 /*--- compute (hf|hh) ---*/
 hrr1_build_hf(Libint->AB,int_stack+309121,int_stack+235033,int_stack+179467,441);
 /*--- compute (m0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+63001,int_stack+33785,int_stack+32630,55);
 /*--- compute (m0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+66466,int_stack+35325,int_stack+33785,55);
 /*--- compute (m0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+71086,int_stack+66466,int_stack+63001,55);
 /*--- compute (m0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+19845,int_stack+37305,int_stack+35325,55);
 /*--- compute (m0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+25785,int_stack+19845,int_stack+66466,55);
 /*--- compute (m0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+126505,int_stack+25785,int_stack+71086,55);
 /*--- compute (m0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+63001,int_stack+39780,int_stack+37305,55);
 /*--- compute (m0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+138055,int_stack+63001,int_stack+19845,55);
 /*--- compute (m0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+149935,int_stack+138055,int_stack+25785,55);
 /*--- compute (m0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+19845,int_stack+149935,int_stack+126505,55);
 /*--- compute (m0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+126505,int_stack+42805,int_stack+39780,55);
 /*--- compute (m0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+165335,int_stack+126505,int_stack+63001,55);
 /*--- compute (m0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+180185,int_stack+165335,int_stack+138055,55);
 /*--- compute (m0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+126505,int_stack+180185,int_stack+149935,55);
 /*--- compute (m0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+149605,int_stack+126505,int_stack+19845,55);
 /*--- compute (lp|hh) ---*/
 hrr1_build_lp(Libint->AB,int_stack+173860,int_stack+149605,int_stack+0,441);
 /*--- compute (kd|hh) ---*/
 hrr1_build_kd(Libint->AB,int_stack+401731,int_stack+173860,int_stack+78877,441);
 /*--- compute (if|hh) ---*/
 hrr1_build_if(Libint->AB,int_stack+496987,int_stack+401731,int_stack+235033,441);
 /*--- compute (hg|hh) ---*/
 hrr1_build_hg(Libint->AB,int_stack+620467,int_stack+496987,int_stack+309121,441);
 /*--- compute (n0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+47821,int_stack+46435,66);
 /*--- compute (n0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+4158,int_stack+49669,int_stack+47821,66);
 /*--- compute (n0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9702,int_stack+4158,int_stack+0,66);
 /*--- compute (n0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18018,int_stack+52045,int_stack+49669,66);
 /*--- compute (n0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+25146,int_stack+18018,int_stack+4158,66);
 /*--- compute (n0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+36234,int_stack+25146,int_stack+9702,66);
 /*--- compute (n0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+55015,int_stack+52045,66);
 /*--- compute (n0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+63001,int_stack+0,int_stack+18018,66);
 /*--- compute (n0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+77257,int_stack+63001,int_stack+25146,66);
 /*--- compute (n0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+8910,int_stack+77257,int_stack+36234,66);
 /*--- compute (n0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+29700,int_stack+58645,int_stack+55015,66);
 /*--- compute (n0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+40590,int_stack+29700,int_stack+0,66);
 /*--- compute (n0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+95737,int_stack+40590,int_stack+63001,66);
 /*--- compute (n0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+29700,int_stack+95737,int_stack+77257,66);
 /*--- compute (n0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+57420,int_stack+29700,int_stack+8910,66);
 /*--- compute (mp|hh) ---*/
 hrr1_build_mp(Libint->AB,int_stack+233395,int_stack+57420,int_stack+149605,441);
 /*--- compute (ld|hh) ---*/
 hrr1_build_ld(Libint->AB,int_stack+0,int_stack+233395,int_stack+173860,441);
 /*--- compute (kf|hh) ---*/
 hrr1_build_kf(Libint->AB,int_stack+119070,int_stack+0,int_stack+401731,441);
 /*--- compute (ig|hh) ---*/
 hrr1_build_ig(Libint->AB,int_stack+277830,int_stack+119070,int_stack+496987,441);
 /*--- compute (hh|hh) ---*/
 hrr1_build_hh(Libint->AB,int_stack+0,int_stack+277830,int_stack+620467,441);
 return int_stack+0;}
