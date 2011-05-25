#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hfhg(Libint_t*, prim_data*);

  /* Computes quartets of (hf|hg) integrals */

REALTYPE *hrr_order_hfhg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][5] = int_stack + 0;
 Libint->vrr_classes[5][6] = int_stack + 441;
 Libint->vrr_classes[5][7] = int_stack + 1029;
 Libint->vrr_classes[5][8] = int_stack + 1785;
 Libint->vrr_classes[5][9] = int_stack + 2730;
 Libint->vrr_classes[6][5] = int_stack + 3885;
 Libint->vrr_classes[6][6] = int_stack + 4473;
 Libint->vrr_classes[6][7] = int_stack + 5257;
 Libint->vrr_classes[6][8] = int_stack + 6265;
 Libint->vrr_classes[6][9] = int_stack + 7525;
 Libint->vrr_classes[7][5] = int_stack + 9065;
 Libint->vrr_classes[7][6] = int_stack + 9821;
 Libint->vrr_classes[7][7] = int_stack + 10829;
 Libint->vrr_classes[7][8] = int_stack + 12125;
 Libint->vrr_classes[7][9] = int_stack + 13745;
 Libint->vrr_classes[8][5] = int_stack + 15725;
 Libint->vrr_classes[8][6] = int_stack + 16670;
 Libint->vrr_classes[8][7] = int_stack + 17930;
 Libint->vrr_classes[8][8] = int_stack + 19550;
 Libint->vrr_classes[8][9] = int_stack + 21575;
 memset(int_stack,0,24050*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 24050;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hfhg(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+24050,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+25373,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+27137,int_stack+25373,int_stack+24050,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+29783,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+32051,int_stack+29783,int_stack+25373,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+35579,int_stack+32051,int_stack+27137,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+24050,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+39989,int_stack+24050,int_stack+29783,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+24050,int_stack+39989,int_stack+32051,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+39989,int_stack+24050,int_stack+35579,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+24050,int_stack+4473,int_stack+3885,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+25814,int_stack+5257,int_stack+4473,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+28166,int_stack+25814,int_stack+24050,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+31694,int_stack+6265,int_stack+5257,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+34718,int_stack+31694,int_stack+25814,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+34718,int_stack+28166,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+24050,int_stack+7525,int_stack+6265,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+46604,int_stack+24050,int_stack+31694,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+24050,int_stack+46604,int_stack+34718,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+46604,int_stack+24050,int_stack+0,28);
 /*--- compute (hp|hg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+55424,int_stack+46604,int_stack+39989,315);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+9821,int_stack+9065,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+10829,int_stack+9821,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+24050,int_stack+12125,int_stack+10829,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+27938,int_stack+24050,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+33986,int_stack+27938,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+13745,int_stack+12125,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+4860,int_stack+0,int_stack+24050,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+75269,int_stack+4860,int_stack+27938,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+0,int_stack+75269,int_stack+33986,36);
 /*--- compute (ip|hg) ---*/
 hrr1_build_ip(Libint->AB,int_stack+75269,int_stack+0,int_stack+46604,315);
 /*--- compute (hd|hg) ---*/
 hrr1_build_hd(Libint->AB,int_stack+101729,int_stack+75269,int_stack+55424,315);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+24050,int_stack+16670,int_stack+15725,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+26885,int_stack+17930,int_stack+16670,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+30665,int_stack+26885,int_stack+24050,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+36335,int_stack+19550,int_stack+17930,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+41195,int_stack+36335,int_stack+26885,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+48755,int_stack+41195,int_stack+30665,45);
 /*--- compute (l0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+24050,int_stack+21575,int_stack+19550,45);
 /*--- compute (l0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+58205,int_stack+24050,int_stack+36335,45);
 /*--- compute (l0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+11340,int_stack+58205,int_stack+41195,45);
 /*--- compute (l0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+58205,int_stack+11340,int_stack+48755,45);
 /*--- compute (kp|hg) ---*/
 hrr1_build_kp(Libint->AB,int_stack+11340,int_stack+58205,int_stack+0,315);
 /*--- compute (id|hg) ---*/
 hrr1_build_id(Libint->AB,int_stack+141419,int_stack+11340,int_stack+75269,315);
 /*--- compute (hf|hg) ---*/
 hrr1_build_hf(Libint->AB,int_stack+0,int_stack+141419,int_stack+101729,315);
 return int_stack+0;}
