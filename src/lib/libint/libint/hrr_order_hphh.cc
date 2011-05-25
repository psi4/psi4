#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hphh(Libint_t*, prim_data*);

  /* Computes quartets of (hp|hh) integrals */

REALTYPE *hrr_order_hphh(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,12299*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 12299;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hphh(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+12299,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+13622,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+15386,int_stack+13622,int_stack+12299,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18032,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+20300,int_stack+18032,int_stack+13622,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+23828,int_stack+20300,int_stack+15386,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+12299,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+28238,int_stack+12299,int_stack+18032,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+32774,int_stack+28238,int_stack+20300,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+15134,int_stack+32774,int_stack+23828,21);
 /*--- compute (h0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+21749,int_stack+3885,int_stack+2730,21);
 /*--- compute (h0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+38654,int_stack+21749,int_stack+12299,21);
 /*--- compute (h0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+44324,int_stack+38654,int_stack+28238,21);
 /*--- compute (h0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+21749,int_stack+44324,int_stack+32774,21);
 /*--- compute (h0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+30569,int_stack+21749,int_stack+15134,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+12299,int_stack+5859,int_stack+5271,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14063,int_stack+6643,int_stack+5859,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+16415,int_stack+14063,int_stack+12299,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+19943,int_stack+7651,int_stack+6643,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+22967,int_stack+19943,int_stack+14063,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+39830,int_stack+22967,int_stack+16415,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+12299,int_stack+8911,int_stack+7651,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+45710,int_stack+12299,int_stack+19943,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+45710,int_stack+22967,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+16079,int_stack+0,int_stack+39830,28);
 /*--- compute (i0|mp) ---*/
 hrr3_build_mp(Libint->CD,int_stack+39830,int_stack+10451,int_stack+8911,28);
 /*--- compute (i0|ld) ---*/
 hrr3_build_ld(Libint->CD,int_stack+51758,int_stack+39830,int_stack+12299,28);
 /*--- compute (i0|kf) ---*/
 hrr3_build_kf(Libint->CD,int_stack+59318,int_stack+51758,int_stack+45710,28);
 /*--- compute (i0|ig) ---*/
 hrr3_build_ig(Libint->CD,int_stack+39830,int_stack+59318,int_stack+0,28);
 /*--- compute (i0|hh) ---*/
 hrr3_build_hh(Libint->CD,int_stack+0,int_stack+39830,int_stack+16079,28);
 /*--- compute (hp|hh) ---*/
 hrr1_build_hp(Libint->AB,int_stack+39830,int_stack+0,int_stack+30569,441);
 return int_stack+39830;}
