#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gghh(Libint_t*, prim_data*);

  /* Computes quartets of (gg|hh) integrals */

REALTYPE *hrr_order_gghh(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[4][8] = int_stack + 1275;
 Libint->vrr_classes[4][9] = int_stack + 1950;
 Libint->vrr_classes[4][10] = int_stack + 2775;
 Libint->vrr_classes[5][5] = int_stack + 3765;
 Libint->vrr_classes[5][6] = int_stack + 4206;
 Libint->vrr_classes[5][7] = int_stack + 4794;
 Libint->vrr_classes[5][8] = int_stack + 5550;
 Libint->vrr_classes[5][9] = int_stack + 6495;
 Libint->vrr_classes[5][10] = int_stack + 7650;
 Libint->vrr_classes[6][5] = int_stack + 9036;
 Libint->vrr_classes[6][6] = int_stack + 9624;
 Libint->vrr_classes[6][7] = int_stack + 10408;
 Libint->vrr_classes[6][8] = int_stack + 11416;
 Libint->vrr_classes[6][9] = int_stack + 12676;
 Libint->vrr_classes[6][10] = int_stack + 14216;
 Libint->vrr_classes[7][5] = int_stack + 16064;
 Libint->vrr_classes[7][6] = int_stack + 16820;
 Libint->vrr_classes[7][7] = int_stack + 17828;
 Libint->vrr_classes[7][8] = int_stack + 19124;
 Libint->vrr_classes[7][9] = int_stack + 20744;
 Libint->vrr_classes[7][10] = int_stack + 22724;
 Libint->vrr_classes[8][5] = int_stack + 25100;
 Libint->vrr_classes[8][6] = int_stack + 26045;
 Libint->vrr_classes[8][7] = int_stack + 27305;
 Libint->vrr_classes[8][8] = int_stack + 28925;
 Libint->vrr_classes[8][9] = int_stack + 30950;
 Libint->vrr_classes[8][10] = int_stack + 33425;
 memset(int_stack,0,36395*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 36395;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gghh(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+36395,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+37340,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+38600,int_stack+37340,int_stack+36395,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+40490,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+42110,int_stack+40490,int_stack+37340,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+44630,int_stack+42110,int_stack+38600,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+36395,int_stack+1950,int_stack+1275,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+47780,int_stack+36395,int_stack+40490,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+51020,int_stack+47780,int_stack+42110,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+38420,int_stack+51020,int_stack+44630,15);
 /*--- compute (g0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+43145,int_stack+2775,int_stack+1950,15);
 /*--- compute (g0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+55220,int_stack+43145,int_stack+36395,15);
 /*--- compute (g0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+59270,int_stack+55220,int_stack+47780,15);
 /*--- compute (g0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+43145,int_stack+59270,int_stack+51020,15);
 /*--- compute (g0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+49445,int_stack+43145,int_stack+38420,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+36395,int_stack+4206,int_stack+3765,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+37718,int_stack+4794,int_stack+4206,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+39482,int_stack+37718,int_stack+36395,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+42128,int_stack+5550,int_stack+4794,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+44396,int_stack+42128,int_stack+37718,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+56060,int_stack+44396,int_stack+39482,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+36395,int_stack+6495,int_stack+5550,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+0,int_stack+36395,int_stack+42128,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+60470,int_stack+0,int_stack+44396,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+39230,int_stack+60470,int_stack+56060,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+56060,int_stack+7650,int_stack+6495,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+66350,int_stack+56060,int_stack+36395,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+72020,int_stack+66350,int_stack+0,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+0,int_stack+72020,int_stack+60470,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+56060,int_stack+0,int_stack+39230,21);
 /*--- compute (gp|hh) ---*/
 hrr1_build_gp(Libint->AB,int_stack+65321,int_stack+56060,int_stack+49445,441);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+9624,int_stack+9036,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+10408,int_stack+9624,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4116,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+36395,int_stack+11416,int_stack+10408,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+39419,int_stack+36395,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+44123,int_stack+39419,int_stack+4116,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+12676,int_stack+11416,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+3780,int_stack+0,int_stack+36395,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+85166,int_stack+3780,int_stack+39419,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+93006,int_stack+85166,int_stack+44123,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+36395,int_stack+14216,int_stack+12676,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+41015,int_stack+36395,int_stack+0,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+101826,int_stack+41015,int_stack+3780,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+0,int_stack+101826,int_stack+85166,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+36395,int_stack+0,int_stack+93006,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+85166,int_stack+36395,int_stack+56060,441);
 /*--- compute (gd|hh) ---*/
 hrr1_build_gd(Libint->AB,int_stack+112949,int_stack+85166,int_stack+65321,441);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+16820,int_stack+16064,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+17828,int_stack+16820,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9828,int_stack+19124,int_stack+17828,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+48743,int_stack+9828,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+54791,int_stack+48743,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+20744,int_stack+19124,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+62351,int_stack+0,int_stack+9828,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+4860,int_stack+62351,int_stack+48743,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+70127,int_stack+4860,int_stack+54791,36);
 /*--- compute (k0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+48743,int_stack+22724,int_stack+20744,36);
 /*--- compute (k0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+14940,int_stack+48743,int_stack+0,36);
 /*--- compute (k0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+48743,int_stack+14940,int_stack+62351,36);
 /*--- compute (k0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+152639,int_stack+48743,int_stack+4860,36);
 /*--- compute (k0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+48743,int_stack+152639,int_stack+70127,36);
 /*--- compute (ip|hh) ---*/
 hrr1_build_ip(Libint->AB,int_stack+152639,int_stack+48743,int_stack+36395,441);
 /*--- compute (hd|hh) ---*/
 hrr1_build_hd(Libint->AB,int_stack+189683,int_stack+152639,int_stack+85166,441);
 /*--- compute (gf|hh) ---*/
 hrr1_build_gf(Libint->AB,int_stack+245249,int_stack+189683,int_stack+112949,441);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+36395,int_stack+26045,int_stack+25100,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+39230,int_stack+27305,int_stack+26045,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+43010,int_stack+39230,int_stack+36395,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+64619,int_stack+28925,int_stack+27305,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+69479,int_stack+64619,int_stack+39230,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+77039,int_stack+69479,int_stack+43010,45);
 /*--- compute (l0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+36395,int_stack+30950,int_stack+28925,45);
 /*--- compute (l0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+86489,int_stack+36395,int_stack+64619,45);
 /*--- compute (l0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+96209,int_stack+86489,int_stack+69479,45);
 /*--- compute (l0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+108809,int_stack+96209,int_stack+77039,45);
 /*--- compute (l0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+64619,int_stack+33425,int_stack+30950,45);
 /*--- compute (l0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+72044,int_stack+64619,int_stack+36395,45);
 /*--- compute (l0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+122984,int_stack+72044,int_stack+86489,45);
 /*--- compute (l0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+64619,int_stack+122984,int_stack+96209,45);
 /*--- compute (l0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+122984,int_stack+64619,int_stack+108809,45);
 /*--- compute (kp|hh) ---*/
 hrr1_build_kp(Libint->AB,int_stack+64619,int_stack+122984,int_stack+48743,441);
 /*--- compute (id|hh) ---*/
 hrr1_build_id(Libint->AB,int_stack+311399,int_stack+64619,int_stack+152639,441);
 /*--- compute (hf|hh) ---*/
 hrr1_build_hf(Libint->AB,int_stack+0,int_stack+311399,int_stack+189683,441);
 /*--- compute (gg|hh) ---*/
 hrr1_build_gg(Libint->AB,int_stack+92610,int_stack+0,int_stack+245249,441);
 return int_stack+92610;}
