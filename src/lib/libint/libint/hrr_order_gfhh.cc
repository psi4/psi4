#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gfhh(Libint_t*, prim_data*);

  /* Computes quartets of (gf|hh) integrals */

REALTYPE *hrr_order_gfhh(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,25100*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 25100;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gfhh(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+25100,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+26045,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+27305,int_stack+26045,int_stack+25100,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+29195,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+30815,int_stack+29195,int_stack+26045,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+33335,int_stack+30815,int_stack+27305,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+25100,int_stack+1950,int_stack+1275,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+36485,int_stack+25100,int_stack+29195,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+39725,int_stack+36485,int_stack+30815,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+27125,int_stack+39725,int_stack+33335,15);
 /*--- compute (g0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+31850,int_stack+2775,int_stack+1950,15);
 /*--- compute (g0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+43925,int_stack+31850,int_stack+25100,15);
 /*--- compute (g0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+47975,int_stack+43925,int_stack+36485,15);
 /*--- compute (g0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+31850,int_stack+47975,int_stack+39725,15);
 /*--- compute (g0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+38150,int_stack+31850,int_stack+27125,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+25100,int_stack+4206,int_stack+3765,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+26423,int_stack+4794,int_stack+4206,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+28187,int_stack+26423,int_stack+25100,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+30833,int_stack+5550,int_stack+4794,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+33101,int_stack+30833,int_stack+26423,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+44765,int_stack+33101,int_stack+28187,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+25100,int_stack+6495,int_stack+5550,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+49175,int_stack+25100,int_stack+30833,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+49175,int_stack+33101,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+27935,int_stack+0,int_stack+44765,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+44765,int_stack+7650,int_stack+6495,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+53711,int_stack+44765,int_stack+25100,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+59381,int_stack+53711,int_stack+49175,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+44765,int_stack+59381,int_stack+0,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+53585,int_stack+44765,int_stack+27935,21);
 /*--- compute (gp|hh) ---*/
 hrr1_build_gp(Libint->AB,int_stack+62846,int_stack+53585,int_stack+38150,441);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+9624,int_stack+9036,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1764,int_stack+10408,int_stack+9624,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+4116,int_stack+1764,int_stack+0,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+25100,int_stack+11416,int_stack+10408,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+28124,int_stack+25100,int_stack+1764,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+32828,int_stack+28124,int_stack+4116,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+12676,int_stack+11416,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+3780,int_stack+0,int_stack+25100,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+38708,int_stack+3780,int_stack+28124,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+82691,int_stack+38708,int_stack+32828,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+25100,int_stack+14216,int_stack+12676,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+29720,int_stack+25100,int_stack+0,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+91511,int_stack+29720,int_stack+3780,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+0,int_stack+91511,int_stack+38708,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+25100,int_stack+0,int_stack+82691,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+82691,int_stack+25100,int_stack+53585,441);
 /*--- compute (gd|hh) ---*/
 hrr1_build_gd(Libint->AB,int_stack+110474,int_stack+82691,int_stack+62846,441);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+16820,int_stack+16064,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+17828,int_stack+16820,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9828,int_stack+19124,int_stack+17828,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+37448,int_stack+9828,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+43496,int_stack+37448,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+20744,int_stack+19124,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+51056,int_stack+0,int_stack+9828,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+4860,int_stack+51056,int_stack+37448,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+58832,int_stack+4860,int_stack+43496,36);
 /*--- compute (k0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+37448,int_stack+22724,int_stack+20744,36);
 /*--- compute (k0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+14940,int_stack+37448,int_stack+0,36);
 /*--- compute (k0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+37448,int_stack+14940,int_stack+51056,36);
 /*--- compute (k0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+150164,int_stack+37448,int_stack+4860,36);
 /*--- compute (k0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+37448,int_stack+150164,int_stack+58832,36);
 /*--- compute (ip|hh) ---*/
 hrr1_build_ip(Libint->AB,int_stack+150164,int_stack+37448,int_stack+25100,441);
 /*--- compute (hd|hh) ---*/
 hrr1_build_hd(Libint->AB,int_stack+0,int_stack+150164,int_stack+82691,441);
 /*--- compute (gf|hh) ---*/
 hrr1_build_gf(Libint->AB,int_stack+150164,int_stack+0,int_stack+110474,441);
 return int_stack+150164;}
