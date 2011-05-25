#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdgd(Libint_t*, prim_data*);

  /* Computes quartets of (fd|gd) integrals */

REALTYPE *hrr_order_fdgd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 Libint->vrr_classes[3][6] = int_stack + 360;
 Libint->vrr_classes[4][4] = int_stack + 640;
 Libint->vrr_classes[4][5] = int_stack + 865;
 Libint->vrr_classes[4][6] = int_stack + 1180;
 Libint->vrr_classes[5][4] = int_stack + 1600;
 Libint->vrr_classes[5][5] = int_stack + 1915;
 Libint->vrr_classes[5][6] = int_stack + 2356;
 memset(int_stack,0,2944*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 2944;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdgd(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2944,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+3394,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+4024,int_stack+3394,int_stack+2944,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2944,int_stack+865,int_stack+640,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4924,int_stack+1180,int_stack+865,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+0,int_stack+4924,int_stack+2944,15);
 /*--- compute (fp|gd) ---*/
 hrr1_build_fp(Libint->AB,int_stack+4924,int_stack+0,int_stack+4024,90);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+2944,int_stack+1915,int_stack+1600,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7624,int_stack+2356,int_stack+1915,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+8947,int_stack+7624,int_stack+2944,21);
 /*--- compute (gp|gd) ---*/
 hrr1_build_gp(Libint->AB,int_stack+10837,int_stack+8947,int_stack+0,90);
 /*--- compute (fd|gd) ---*/
 hrr1_build_fd(Libint->AB,int_stack+14887,int_stack+10837,int_stack+4924,90);
 return int_stack+14887;}
