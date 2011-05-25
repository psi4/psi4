#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hphg(Libint_t*, prim_data*);

  /* Computes quartets of (hp|hg) integrals */

REALTYPE *hrr_order_hphg(Libint_t *Libint, int num_prim_comb)
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
 memset(int_stack,0,9065*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 9065;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hphg(Libint, Data);
   Data++;
 }
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9065,int_stack+441,int_stack+0,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10388,int_stack+1029,int_stack+441,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+12152,int_stack+10388,int_stack+9065,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+14798,int_stack+1785,int_stack+1029,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+17066,int_stack+14798,int_stack+10388,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+20594,int_stack+17066,int_stack+12152,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+9065,int_stack+2730,int_stack+1785,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+25004,int_stack+9065,int_stack+14798,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+9065,int_stack+25004,int_stack+17066,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+25004,int_stack+9065,int_stack+20594,21);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9065,int_stack+4473,int_stack+3885,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10829,int_stack+5257,int_stack+4473,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+13181,int_stack+10829,int_stack+9065,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+16709,int_stack+6265,int_stack+5257,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+19733,int_stack+16709,int_stack+10829,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+0,int_stack+19733,int_stack+13181,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+9065,int_stack+7525,int_stack+6265,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+31619,int_stack+9065,int_stack+16709,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+5880,int_stack+31619,int_stack+19733,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+13720,int_stack+5880,int_stack+0,28);
 /*--- compute (hp|hg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+31619,int_stack+13720,int_stack+25004,315);
 return int_stack+31619;}
