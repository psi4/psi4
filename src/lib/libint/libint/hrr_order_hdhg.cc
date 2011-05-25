#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hdhg(Libint_t*, prim_data*);

  /* Computes quartets of (hd|hg) integrals */

REALTYPE *hrr_order_hdhg(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,15725*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 15725;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hdhg(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+15725,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+17048,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+18812,int_stack+17048,int_stack+15725,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+21458,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+23726,int_stack+21458,int_stack+17048,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+27254,int_stack+23726,int_stack+18812,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+15725,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+31664,int_stack+15725,int_stack+21458,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+15725,int_stack+31664,int_stack+23726,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+31664,int_stack+15725,int_stack+27254,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+15725,int_stack+4473,int_stack+3885,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+17489,int_stack+5257,int_stack+4473,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+19841,int_stack+17489,int_stack+15725,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+23369,int_stack+6265,int_stack+5257,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+26393,int_stack+23369,int_stack+17489,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+26393,int_stack+19841,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+15725,int_stack+7525,int_stack+6265,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+38279,int_stack+15725,int_stack+23369,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+15725,int_stack+38279,int_stack+26393,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+38279,int_stack+15725,int_stack+0,28);
 /*--- compute (hp|hg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+47099,int_stack+38279,int_stack+31664,315);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+9821,int_stack+9065,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+10829,int_stack+9821,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+15725,int_stack+12125,int_stack+10829,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+19613,int_stack+15725,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+25661,int_stack+19613,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+13745,int_stack+12125,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+4860,int_stack+0,int_stack+15725,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+66944,int_stack+4860,int_stack+19613,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+0,int_stack+66944,int_stack+25661,36);
 /*--- compute (ip|hg) ---*/
 hrr1_build_ip(Libint->AB,int_stack+11340,int_stack+0,int_stack+38279,315);
 /*--- compute (hd|hg) ---*/
 hrr1_build_hd(Libint->AB,int_stack+66944,int_stack+11340,int_stack+47099,315);
 return int_stack+66944;}
