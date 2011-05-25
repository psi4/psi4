#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdff(Libint_t*, prim_data*);

  /* Computes quartets of (fd|ff) integrals */

REALTYPE *hrr_order_fdff(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][3] = int_stack + 0;
 Libint->vrr_classes[3][4] = int_stack + 100;
 Libint->vrr_classes[3][5] = int_stack + 250;
 Libint->vrr_classes[3][6] = int_stack + 460;
 Libint->vrr_classes[4][3] = int_stack + 740;
 Libint->vrr_classes[4][4] = int_stack + 890;
 Libint->vrr_classes[4][5] = int_stack + 1115;
 Libint->vrr_classes[4][6] = int_stack + 1430;
 Libint->vrr_classes[5][3] = int_stack + 1850;
 Libint->vrr_classes[5][4] = int_stack + 2060;
 Libint->vrr_classes[5][5] = int_stack + 2375;
 Libint->vrr_classes[5][6] = int_stack + 2816;
 memset(int_stack,0,3404*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3404;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdff(Libint, Data);
   Data++;
 }
 /*--- compute (f0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+3404,int_stack+100,int_stack+0,10);
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3704,int_stack+250,int_stack+100,10);
 /*--- compute (f0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+4154,int_stack+3704,int_stack+3404,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4754,int_stack+460,int_stack+250,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5384,int_stack+4754,int_stack+3704,10);
 /*--- compute (f0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+6284,int_stack+5384,int_stack+4154,10);
 /*--- compute (g0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+3404,int_stack+890,int_stack+740,15);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3854,int_stack+1115,int_stack+890,15);
 /*--- compute (g0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+4529,int_stack+3854,int_stack+3404,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+1430,int_stack+1115,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+7284,int_stack+0,int_stack+3854,15);
 /*--- compute (g0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+0,int_stack+7284,int_stack+4529,15);
 /*--- compute (fp|ff) ---*/
 hrr1_build_fp(Libint->AB,int_stack+7284,int_stack+0,int_stack+6284,100);
 /*--- compute (h0|fp) ---*/
 hrr3_build_fp(Libint->CD,int_stack+10284,int_stack+2060,int_stack+1850,21);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3404,int_stack+2375,int_stack+2060,21);
 /*--- compute (h0|fd) ---*/
 hrr3_build_fd(Libint->CD,int_stack+4349,int_stack+3404,int_stack+10284,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5609,int_stack+2816,int_stack+2375,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+1500,int_stack+5609,int_stack+3404,21);
 /*--- compute (h0|ff) ---*/
 hrr3_build_ff(Libint->CD,int_stack+10284,int_stack+1500,int_stack+4349,21);
 /*--- compute (gp|ff) ---*/
 hrr1_build_gp(Libint->AB,int_stack+1500,int_stack+10284,int_stack+0,100);
 /*--- compute (fd|ff) ---*/
 hrr1_build_fd(Libint->AB,int_stack+10284,int_stack+1500,int_stack+7284,100);
 return int_stack+10284;}
