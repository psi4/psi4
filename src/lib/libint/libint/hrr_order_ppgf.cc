#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ppgf(Libint_t*, prim_data*);

  /* Computes quartets of (pp|gf) integrals */

REALTYPE *hrr_order_ppgf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[1][4] = int_stack + 0;
 Libint->vrr_classes[1][5] = int_stack + 45;
 Libint->vrr_classes[1][6] = int_stack + 108;
 Libint->vrr_classes[1][7] = int_stack + 192;
 Libint->vrr_classes[2][4] = int_stack + 300;
 Libint->vrr_classes[2][5] = int_stack + 390;
 Libint->vrr_classes[2][6] = int_stack + 516;
 Libint->vrr_classes[2][7] = int_stack + 684;
 memset(int_stack,0,900*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 900;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ppgf(Libint, Data);
   Data++;
 }
 /*--- compute (p0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+900,int_stack+45,int_stack+0,3);
 /*--- compute (p0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1035,int_stack+108,int_stack+45,3);
 /*--- compute (p0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1224,int_stack+1035,int_stack+900,3);
 /*--- compute (p0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+1494,int_stack+192,int_stack+108,3);
 /*--- compute (p0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+1746,int_stack+1494,int_stack+1035,3);
 /*--- compute (p0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+2124,int_stack+1746,int_stack+1224,3);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+900,int_stack+390,int_stack+300,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1170,int_stack+516,int_stack+390,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1548,int_stack+1170,int_stack+900,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+0,int_stack+684,int_stack+516,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2574,int_stack+0,int_stack+1170,6);
 /*--- compute (d0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+0,int_stack+2574,int_stack+1548,6);
 /*--- compute (pp|gf) ---*/
 hrr1_build_pp(Libint->AB,int_stack+2574,int_stack+0,int_stack+2124,150);
 return int_stack+2574;}
