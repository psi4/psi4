#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hghg(Libint_t*, prim_data*);

  /* Computes quartets of (hg|hg) integrals */

REALTYPE *hrr_order_hghg(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[9][5] = int_stack + 24050;
 Libint->vrr_classes[9][6] = int_stack + 25205;
 Libint->vrr_classes[9][7] = int_stack + 26745;
 Libint->vrr_classes[9][8] = int_stack + 28725;
 Libint->vrr_classes[9][9] = int_stack + 31200;
 memset(int_stack,0,34225*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 34225;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hghg(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+34225,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+35548,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+37312,int_stack+35548,int_stack+34225,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+39958,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+42226,int_stack+39958,int_stack+35548,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+45754,int_stack+42226,int_stack+37312,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+34225,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+50164,int_stack+34225,int_stack+39958,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+34225,int_stack+50164,int_stack+42226,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+50164,int_stack+34225,int_stack+45754,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+34225,int_stack+4473,int_stack+3885,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+35989,int_stack+5257,int_stack+4473,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+38341,int_stack+35989,int_stack+34225,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+41869,int_stack+6265,int_stack+5257,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+44893,int_stack+41869,int_stack+35989,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+44893,int_stack+38341,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+34225,int_stack+7525,int_stack+6265,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+56779,int_stack+34225,int_stack+41869,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+34225,int_stack+56779,int_stack+44893,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+56779,int_stack+34225,int_stack+0,28);
 /*--- compute (hp|hg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+65599,int_stack+56779,int_stack+50164,315);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+9821,int_stack+9065,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+10829,int_stack+9821,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+34225,int_stack+12125,int_stack+10829,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+38113,int_stack+34225,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+44161,int_stack+38113,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+13745,int_stack+12125,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+4860,int_stack+0,int_stack+34225,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+85444,int_stack+4860,int_stack+38113,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+0,int_stack+85444,int_stack+44161,36);
 /*--- compute (ip|hg) ---*/
 hrr1_build_ip(Libint->AB,int_stack+85444,int_stack+0,int_stack+56779,315);
 /*--- compute (hd|hg) ---*/
 hrr1_build_hd(Libint->AB,int_stack+111904,int_stack+85444,int_stack+65599,315);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+34225,int_stack+16670,int_stack+15725,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+37060,int_stack+17930,int_stack+16670,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+40840,int_stack+37060,int_stack+34225,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+46510,int_stack+19550,int_stack+17930,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+51370,int_stack+46510,int_stack+37060,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+58930,int_stack+51370,int_stack+40840,45);
 /*--- compute (l0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+34225,int_stack+21575,int_stack+19550,45);
 /*--- compute (l0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+68380,int_stack+34225,int_stack+46510,45);
 /*--- compute (l0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+34225,int_stack+68380,int_stack+51370,45);
 /*--- compute (l0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+68380,int_stack+34225,int_stack+58930,45);
 /*--- compute (kp|hg) ---*/
 hrr1_build_kp(Libint->AB,int_stack+34225,int_stack+68380,int_stack+0,315);
 /*--- compute (id|hg) ---*/
 hrr1_build_id(Libint->AB,int_stack+151594,int_stack+34225,int_stack+85444,315);
 /*--- compute (hf|hg) ---*/
 hrr1_build_hf(Libint->AB,int_stack+204514,int_stack+151594,int_stack+111904,315);
 /*--- compute (m0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+25205,int_stack+24050,55);
 /*--- compute (m0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+3465,int_stack+26745,int_stack+25205,55);
 /*--- compute (m0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+8085,int_stack+3465,int_stack+0,55);
 /*--- compute (m0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+15015,int_stack+28725,int_stack+26745,55);
 /*--- compute (m0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+82555,int_stack+15015,int_stack+3465,55);
 /*--- compute (m0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+91795,int_stack+82555,int_stack+8085,55);
 /*--- compute (m0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+31200,int_stack+28725,55);
 /*--- compute (m0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+20955,int_stack+0,int_stack+15015,55);
 /*--- compute (m0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+20955,int_stack+82555,55);
 /*--- compute (m0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+15400,int_stack+0,int_stack+91795,55);
 /*--- compute (lp|hg) ---*/
 hrr1_build_lp(Libint->AB,int_stack+82555,int_stack+15400,int_stack+68380,315);
 /*--- compute (kd|hg) ---*/
 hrr1_build_kd(Libint->AB,int_stack+270664,int_stack+82555,int_stack+34225,315);
 /*--- compute (if|hg) ---*/
 hrr1_build_if(Libint->AB,int_stack+0,int_stack+270664,int_stack+151594,315);
 /*--- compute (hg|hg) ---*/
 hrr1_build_hg(Libint->AB,int_stack+88200,int_stack+0,int_stack+204514,315);
 return int_stack+88200;}
