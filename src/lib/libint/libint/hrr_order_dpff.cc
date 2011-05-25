#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_dpff(Libint_t*, prim_data*);

  /* Computes quartets of (dp|ff) integrals */

REALTYPE *hrr_order_dpff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][3] = int_stack + 0;
 Libint->vrr_classes[2][4] = int_stack + 60;
 Libint->vrr_classes[2][5] = int_stack + 150;
 Libint->vrr_classes[2][6] = int_stack + 276;
 Libint->vrr_classes[3][3] = int_stack + 444;
 Libint->vrr_classes[3][4] = int_stack + 544;
 Libint->vrr_classes[3][5] = int_stack + 694;
 Libint->vrr_classes[3][6] = int_stack + 904;
 memset(int_stack,0,1184*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 1184;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_dpff(Libint, Data);
   Data++;
 }
 /*--- compute (d0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1184,int_stack+60,int_stack+0,6);
 /*--- compute (d0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1364,int_stack+150,int_stack+60,6);
 /*--- compute (d0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+1634,int_stack+1364,int_stack+1184,6);
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1994,int_stack+276,int_stack+150,6);
 /*--- compute (d0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2372,int_stack+1994,int_stack+1364,6);
 /*--- compute (d0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+2912,int_stack+2372,int_stack+1634,6);
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+1184,int_stack+544,int_stack+444,10);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+1484,int_stack+694,int_stack+544,10);
 /*--- compute (f0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+1934,int_stack+1484,int_stack+1184,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+904,int_stack+694,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3512,int_stack+0,int_stack+1484,10);
 /*--- compute (f0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+0,int_stack+3512,int_stack+1934,10);
 /*--- compute (dp|ff) ---*/
 hrr1_build_dp(Libint->AB,int_stack+1000,int_stack+0,int_stack+2912,100);
 return int_stack+1000;}
