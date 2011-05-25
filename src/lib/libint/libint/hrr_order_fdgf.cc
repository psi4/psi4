#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_fdgf(Libint_t*, prim_data*);

  /* Computes quartets of (fd|gf) integrals */

REALTYPE *hrr_order_fdgf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[3][4] = int_stack + 0;
 Libint->vrr_classes[3][5] = int_stack + 150;
 Libint->vrr_classes[3][6] = int_stack + 360;
 Libint->vrr_classes[3][7] = int_stack + 640;
 Libint->vrr_classes[4][4] = int_stack + 1000;
 Libint->vrr_classes[4][5] = int_stack + 1225;
 Libint->vrr_classes[4][6] = int_stack + 1540;
 Libint->vrr_classes[4][7] = int_stack + 1960;
 Libint->vrr_classes[5][4] = int_stack + 2500;
 Libint->vrr_classes[5][5] = int_stack + 2815;
 Libint->vrr_classes[5][6] = int_stack + 3256;
 Libint->vrr_classes[5][7] = int_stack + 3844;
 memset(int_stack,0,4600*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 4600;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_fdgf(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4600,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5050,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5680,int_stack+5050,int_stack+4600,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6580,int_stack+640,int_stack+360,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7420,int_stack+6580,int_stack+5050,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+8680,int_stack+7420,int_stack+5680,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4600,int_stack+1225,int_stack+1000,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5275,int_stack+1540,int_stack+1225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+6220,int_stack+5275,int_stack+4600,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+0,int_stack+1960,int_stack+1540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10180,int_stack+0,int_stack+5275,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+0,int_stack+10180,int_stack+6220,15);
 /*--- compute (fp|gf) ---*/
 hrr1_build_fp(Libint->AB,int_stack+10180,int_stack+0,int_stack+8680,150);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+4600,int_stack+2815,int_stack+2500,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+5545,int_stack+3256,int_stack+2815,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+6868,int_stack+5545,int_stack+4600,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14680,int_stack+3844,int_stack+3256,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+2250,int_stack+14680,int_stack+5545,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+14680,int_stack+2250,int_stack+6868,21);
 /*--- compute (gp|gf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+2250,int_stack+14680,int_stack+0,150);
 /*--- compute (fd|gf) ---*/
 hrr1_build_fd(Libint->AB,int_stack+14680,int_stack+2250,int_stack+10180,150);
 return int_stack+14680;}
