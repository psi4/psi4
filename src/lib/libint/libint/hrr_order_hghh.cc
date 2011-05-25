#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hghh(Libint_t*, prim_data*);

  /* Computes quartets of (hg|hh) integrals */

REALTYPE *hrr_order_hghh(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,46435*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 46435;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hghh(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+46435,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+47758,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+49522,int_stack+47758,int_stack+46435,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+52168,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+54436,int_stack+52168,int_stack+47758,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+57964,int_stack+54436,int_stack+49522,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+46435,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+62374,int_stack+46435,int_stack+52168,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+66910,int_stack+62374,int_stack+54436,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+49270,int_stack+66910,int_stack+57964,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+55885,int_stack+3885,int_stack+2730,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+72790,int_stack+55885,int_stack+46435,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+78460,int_stack+72790,int_stack+62374,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+55885,int_stack+78460,int_stack+66910,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+64705,int_stack+55885,int_stack+49270,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+46435,int_stack+5859,int_stack+5271,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+48199,int_stack+6643,int_stack+5859,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+50551,int_stack+48199,int_stack+46435,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+54079,int_stack+7651,int_stack+6643,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+57103,int_stack+54079,int_stack+48199,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+73966,int_stack+57103,int_stack+50551,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+46435,int_stack+8911,int_stack+7651,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+79846,int_stack+46435,int_stack+54079,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+79846,int_stack+57103,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+50215,int_stack+0,int_stack+73966,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+73966,int_stack+10451,int_stack+8911,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+85894,int_stack+73966,int_stack+46435,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+93454,int_stack+85894,int_stack+79846,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+73966,int_stack+93454,int_stack+0,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+85726,int_stack+73966,int_stack+50215,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+98074,int_stack+85726,int_stack+64705,441);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+13055,int_stack+12299,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+14063,int_stack+13055,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9828,int_stack+15359,int_stack+14063,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+46435,int_stack+9828,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+52483,int_stack+46435,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+16979,int_stack+15359,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+60043,int_stack+0,int_stack+9828,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+4860,int_stack+60043,int_stack+46435,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+67819,int_stack+4860,int_stack+52483,36);
 /*--- compute (k0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+46435,int_stack+18959,int_stack+16979,36);
 /*--- compute (k0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+125857,int_stack+46435,int_stack+0,36);
 /*--- compute (k0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+46435,int_stack+125857,int_stack+60043,36);
 /*--- compute (k0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+125857,int_stack+46435,int_stack+4860,36);
 /*--- compute (k0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+46435,int_stack+125857,int_stack+67819,36);
 /*--- compute (ip|hh) ---*/
 hrr1_build_ip(Libint->AB,int_stack+125857,int_stack+46435,int_stack+85726,441);
 /*--- compute (hd|hh) ---*/
 hrr1_build_hd(Libint->AB,int_stack+162901,int_stack+125857,int_stack+98074,441);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+62311,int_stack+22280,int_stack+21335,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+65146,int_stack+23540,int_stack+22280,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+68926,int_stack+65146,int_stack+62311,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+74596,int_stack+25160,int_stack+23540,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+79456,int_stack+74596,int_stack+65146,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+87016,int_stack+79456,int_stack+68926,45);
 /*--- compute (l0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+62311,int_stack+27185,int_stack+25160,45);
 /*--- compute (l0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+96466,int_stack+62311,int_stack+74596,45);
 /*--- compute (l0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+106186,int_stack+96466,int_stack+79456,45);
 /*--- compute (l0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+68386,int_stack+106186,int_stack+87016,45);
 /*--- compute (l0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+82561,int_stack+29660,int_stack+27185,45);
 /*--- compute (l0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+0,int_stack+82561,int_stack+62311,45);
 /*--- compute (l0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+12150,int_stack+0,int_stack+96466,45);
 /*--- compute (l0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+82561,int_stack+12150,int_stack+106186,45);
 /*--- compute (l0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+82561,int_stack+68386,45);
 /*--- compute (kp|hh) ---*/
 hrr1_build_kp(Libint->AB,int_stack+62311,int_stack+0,int_stack+46435,441);
 /*--- compute (id|hh) ---*/
 hrr1_build_id(Libint->AB,int_stack+218467,int_stack+62311,int_stack+125857,441);
 /*--- compute (hf|hh) ---*/
 hrr1_build_hf(Libint->AB,int_stack+292555,int_stack+218467,int_stack+162901,441);
 /*--- compute (m0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+46435,int_stack+33785,int_stack+32630,55);
 /*--- compute (m0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+49900,int_stack+35325,int_stack+33785,55);
 /*--- compute (m0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+54520,int_stack+49900,int_stack+46435,55);
 /*--- compute (m0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+19845,int_stack+37305,int_stack+35325,55);
 /*--- compute (m0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+25785,int_stack+19845,int_stack+49900,55);
 /*--- compute (m0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+109939,int_stack+25785,int_stack+54520,55);
 /*--- compute (m0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+46435,int_stack+39780,int_stack+37305,55);
 /*--- compute (m0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+121489,int_stack+46435,int_stack+19845,55);
 /*--- compute (m0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+133369,int_stack+121489,int_stack+25785,55);
 /*--- compute (m0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+19845,int_stack+133369,int_stack+109939,55);
 /*--- compute (m0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+109939,int_stack+42805,int_stack+39780,55);
 /*--- compute (m0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+148769,int_stack+109939,int_stack+46435,55);
 /*--- compute (m0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+37170,int_stack+148769,int_stack+121489,55);
 /*--- compute (m0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+148769,int_stack+37170,int_stack+133369,55);
 /*--- compute (m0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+37170,int_stack+148769,int_stack+19845,55);
 /*--- compute (lp|hh) ---*/
 hrr1_build_lp(Libint->AB,int_stack+109939,int_stack+37170,int_stack+0,441);
 /*--- compute (kd|hh) ---*/
 hrr1_build_kd(Libint->AB,int_stack+385165,int_stack+109939,int_stack+62311,441);
 /*--- compute (if|hh) ---*/
 hrr1_build_if(Libint->AB,int_stack+0,int_stack+385165,int_stack+218467,441);
 /*--- compute (hg|hh) ---*/
 hrr1_build_hg(Libint->AB,int_stack+123480,int_stack+0,int_stack+292555,441);
 return int_stack+123480;}
