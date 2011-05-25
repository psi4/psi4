#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gdff(Libint_t*, prim_data*);

  /* Computes quartets of (gd|ff) integrals */

REALTYPE *hrr_order_gdff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][3] = int_stack + 0;
 Libint->vrr_classes[4][4] = int_stack + 150;
 Libint->vrr_classes[4][5] = int_stack + 375;
 Libint->vrr_classes[4][6] = int_stack + 690;
 Libint->vrr_classes[5][3] = int_stack + 1110;
 Libint->vrr_classes[5][4] = int_stack + 1320;
 Libint->vrr_classes[5][5] = int_stack + 1635;
 Libint->vrr_classes[5][6] = int_stack + 2076;
 Libint->vrr_classes[6][3] = int_stack + 2664;
 Libint->vrr_classes[6][4] = int_stack + 2944;
 Libint->vrr_classes[6][5] = int_stack + 3364;
 Libint->vrr_classes[6][6] = int_stack + 3952;
 memset(int_stack,0,4736*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4736;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gdff(Libint, Data);
   Data++;
 }
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+4736,int_stack+150,int_stack+0,15);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+5186,int_stack+375,int_stack+150,15);
 /*--- compute (g0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+5861,int_stack+5186,int_stack+4736,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6761,int_stack+690,int_stack+375,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7706,int_stack+6761,int_stack+5186,15);
 /*--- compute (g0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+9056,int_stack+7706,int_stack+5861,15);
 /*--- compute (h0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+4736,int_stack+1320,int_stack+1110,21);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+5366,int_stack+1635,int_stack+1320,21);
 /*--- compute (h0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+6311,int_stack+5366,int_stack+4736,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7571,int_stack+2076,int_stack+1635,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+7571,int_stack+5366,21);
 /*--- compute (h0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+10556,int_stack+0,int_stack+6311,21);
 /*--- compute (gp|ff) ---*/
 hrr1_build_gp(Libint->AB,int_stack+12656,int_stack+10556,int_stack+9056,100);
 /*--- compute (i0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+0,int_stack+2944,int_stack+2664,28);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+840,int_stack+3364,int_stack+2944,28);
 /*--- compute (i0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+4736,int_stack+840,int_stack+0,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+6416,int_stack+3952,int_stack+3364,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+2100,int_stack+6416,int_stack+840,28);
 /*--- compute (i0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+6416,int_stack+2100,int_stack+4736,28);
 /*--- compute (hp|ff) ---*/
 hrr1_build_hp(Libint->AB,int_stack+0,int_stack+6416,int_stack+10556,100);
 /*--- compute (gd|ff) ---*/
 hrr1_build_gd(Libint->AB,int_stack+17156,int_stack+0,int_stack+12656,100);
 return int_stack+17156;}
