#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_pphg(Libint_t*, prim_data*);

  /* Computes quartets of (pp|hg) integrals */

REALTYPE *hrr_order_pphg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][5] = int_stack + 0;
 Libint->vrr_classes[1][6] = int_stack + 63;
 Libint->vrr_classes[1][7] = int_stack + 147;
 Libint->vrr_classes[1][8] = int_stack + 255;
 Libint->vrr_classes[1][9] = int_stack + 390;
 Libint->vrr_classes[2][5] = int_stack + 555;
 Libint->vrr_classes[2][6] = int_stack + 681;
 Libint->vrr_classes[2][7] = int_stack + 849;
 Libint->vrr_classes[2][8] = int_stack + 1065;
 Libint->vrr_classes[2][9] = int_stack + 1335;
 memset(int_stack,0,1665*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1665;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_pphg(Libint, Data);
   Data++;
 }
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1665,int_stack+63,int_stack+0,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1854,int_stack+147,int_stack+63,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2106,int_stack+1854,int_stack+1665,3);
 /*--- compute (p0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+2484,int_stack+255,int_stack+147,3);
 /*--- compute (p0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+2808,int_stack+2484,int_stack+1854,3);
 /*--- compute (p0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+3312,int_stack+2808,int_stack+2106,3);
 /*--- compute (p0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+1665,int_stack+390,int_stack+255,3);
 /*--- compute (p0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+3942,int_stack+1665,int_stack+2484,3);
 /*--- compute (p0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+1665,int_stack+3942,int_stack+2808,3);
 /*--- compute (p0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+3942,int_stack+1665,int_stack+3312,3);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1665,int_stack+681,int_stack+555,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2043,int_stack+849,int_stack+681,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2547,int_stack+2043,int_stack+1665,6);
 /*--- compute (d0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+0,int_stack+1065,int_stack+849,6);
 /*--- compute (d0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+4887,int_stack+0,int_stack+2043,6);
 /*--- compute (d0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+5895,int_stack+4887,int_stack+2547,6);
 /*--- compute (d0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+1665,int_stack+1335,int_stack+1065,6);
 /*--- compute (d0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+2475,int_stack+1665,int_stack+0,6);
 /*--- compute (d0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+0,int_stack+2475,int_stack+4887,6);
 /*--- compute (d0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+1680,int_stack+0,int_stack+5895,6);
 /*--- compute (pp|hg) ---*/
 hrr1_build_pp(Libint->AB,int_stack+4887,int_stack+1680,int_stack+3942,315);
 return int_stack+4887;}
