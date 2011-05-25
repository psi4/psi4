#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hdhh(Libint_t*, prim_data*);

  /* Computes quartets of (hd|hh) integrals */

REALTYPE *hrr_order_hdhh(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,21335*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 21335;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hdhh(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+21335,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+22658,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+24422,int_stack+22658,int_stack+21335,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+27068,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+29336,int_stack+27068,int_stack+22658,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+32864,int_stack+29336,int_stack+24422,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+21335,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+37274,int_stack+21335,int_stack+27068,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+41810,int_stack+37274,int_stack+29336,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+24170,int_stack+41810,int_stack+32864,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+30785,int_stack+3885,int_stack+2730,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+47690,int_stack+30785,int_stack+21335,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+53360,int_stack+47690,int_stack+37274,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+30785,int_stack+53360,int_stack+41810,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+39605,int_stack+30785,int_stack+24170,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+21335,int_stack+5859,int_stack+5271,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+23099,int_stack+6643,int_stack+5859,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+25451,int_stack+23099,int_stack+21335,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+28979,int_stack+7651,int_stack+6643,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+32003,int_stack+28979,int_stack+23099,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+48866,int_stack+32003,int_stack+25451,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+21335,int_stack+8911,int_stack+7651,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+54746,int_stack+21335,int_stack+28979,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+54746,int_stack+32003,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+25115,int_stack+0,int_stack+48866,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+48866,int_stack+10451,int_stack+8911,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+60794,int_stack+48866,int_stack+21335,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+68354,int_stack+60794,int_stack+54746,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+48866,int_stack+68354,int_stack+0,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+60626,int_stack+48866,int_stack+25115,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+72974,int_stack+60626,int_stack+39605,441);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+13055,int_stack+12299,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+14063,int_stack+13055,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+9828,int_stack+15359,int_stack+14063,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+21335,int_stack+9828,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+27383,int_stack+21335,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+16979,int_stack+15359,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+34943,int_stack+0,int_stack+9828,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+4860,int_stack+34943,int_stack+21335,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+42719,int_stack+4860,int_stack+27383,36);
 /*--- compute (k0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+21335,int_stack+18959,int_stack+16979,36);
 /*--- compute (k0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+100757,int_stack+21335,int_stack+0,36);
 /*--- compute (k0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+14940,int_stack+100757,int_stack+34943,36);
 /*--- compute (k0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+100757,int_stack+14940,int_stack+4860,36);
 /*--- compute (k0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+100757,int_stack+42719,36);
 /*--- compute (ip|hh) ---*/
 hrr1_build_ip(Libint->AB,int_stack+15876,int_stack+0,int_stack+60626,441);
 /*--- compute (hd|hh) ---*/
 hrr1_build_hd(Libint->AB,int_stack+100757,int_stack+15876,int_stack+72974,441);
 return int_stack+100757;}
