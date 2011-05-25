#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gpgf(Libint_t*, prim_data*);

  /* Computes quartets of (gp|gf) integrals */

REALTYPE *hrr_order_gpgf(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][4] = int_stack + 0;
 Libint->vrr_classes[4][5] = int_stack + 225;
 Libint->vrr_classes[4][6] = int_stack + 540;
 Libint->vrr_classes[4][7] = int_stack + 960;
 Libint->vrr_classes[5][4] = int_stack + 1500;
 Libint->vrr_classes[5][5] = int_stack + 1815;
 Libint->vrr_classes[5][6] = int_stack + 2256;
 Libint->vrr_classes[5][7] = int_stack + 2844;
 memset(int_stack,0,3600*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 3600;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gpgf(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3600,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4275,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5220,int_stack+4275,int_stack+3600,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+6570,int_stack+960,int_stack+540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+7830,int_stack+6570,int_stack+4275,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+9720,int_stack+7830,int_stack+5220,15);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+3600,int_stack+1815,int_stack+1500,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+4545,int_stack+2256,int_stack+1815,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+5868,int_stack+4545,int_stack+3600,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+7758,int_stack+2844,int_stack+2256,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+7758,int_stack+4545,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+2646,int_stack+0,int_stack+5868,21);
 /*--- compute (gp|gf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+11970,int_stack+2646,int_stack+9720,150);
 return int_stack+11970;}
