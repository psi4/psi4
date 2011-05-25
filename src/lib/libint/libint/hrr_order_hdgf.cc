#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_hdgf(Libint_t*, prim_data*);

  /* Computes quartets of (hd|gf) integrals */

REALTYPE *hrr_order_hdgf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[5][4] = int_stack + 0;
 Libint->vrr_classes[5][5] = int_stack + 315;
 Libint->vrr_classes[5][6] = int_stack + 756;
 Libint->vrr_classes[5][7] = int_stack + 1344;
 Libint->vrr_classes[6][4] = int_stack + 2100;
 Libint->vrr_classes[6][5] = int_stack + 2520;
 Libint->vrr_classes[6][6] = int_stack + 3108;
 Libint->vrr_classes[6][7] = int_stack + 3892;
 Libint->vrr_classes[7][4] = int_stack + 4900;
 Libint->vrr_classes[7][5] = int_stack + 5440;
 Libint->vrr_classes[7][6] = int_stack + 6196;
 Libint->vrr_classes[7][7] = int_stack + 7204;
 memset(int_stack,0,8500*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 8500;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_hdgf(Libint, Data);
   Data++;
 }
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+8500,int_stack+315,int_stack+0,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9445,int_stack+756,int_stack+315,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+10768,int_stack+9445,int_stack+8500,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12658,int_stack+1344,int_stack+756,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+14422,int_stack+12658,int_stack+9445,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+17068,int_stack+14422,int_stack+10768,21);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+8500,int_stack+2520,int_stack+2100,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+9760,int_stack+3108,int_stack+2520,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+11524,int_stack+9760,int_stack+8500,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14044,int_stack+3892,int_stack+3108,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+14044,int_stack+9760,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+20218,int_stack+0,int_stack+11524,28);
 /*--- compute (hp|gf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+24418,int_stack+20218,int_stack+17068,150);
 /*--- compute (k0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+0,int_stack+5440,int_stack+4900,36);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1620,int_stack+6196,int_stack+5440,36);
 /*--- compute (k0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+8500,int_stack+1620,int_stack+0,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+11740,int_stack+7204,int_stack+6196,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+3888,int_stack+11740,int_stack+1620,36);
 /*--- compute (k0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+11740,int_stack+3888,int_stack+8500,36);
 /*--- compute (ip|gf) ---*/
 hrr1_build_ip(Libint->AB,int_stack+33868,int_stack+11740,int_stack+20218,150);
 /*--- compute (hd|gf) ---*/
 hrr1_build_hd(Libint->AB,int_stack+0,int_stack+33868,int_stack+24418,150);
 return int_stack+0;}
