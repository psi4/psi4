#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gdgf(Libint_t*, prim_data*);

  /* Computes quartets of (gd|gf) integrals */

REALTYPE *hrr_order_gdgf(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[6][4] = int_stack + 3600;
 Libint->vrr_classes[6][5] = int_stack + 4020;
 Libint->vrr_classes[6][6] = int_stack + 4608;
 Libint->vrr_classes[6][7] = int_stack + 5392;
 memset(int_stack,0,6400*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 6400;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gdgf(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+6400,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7075,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+8020,int_stack+7075,int_stack+6400,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9370,int_stack+960,int_stack+540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10630,int_stack+9370,int_stack+7075,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+12520,int_stack+10630,int_stack+8020,15);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+6400,int_stack+1815,int_stack+1500,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7345,int_stack+2256,int_stack+1815,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+8668,int_stack+7345,int_stack+6400,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+10558,int_stack+2844,int_stack+2256,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+10558,int_stack+7345,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+14770,int_stack+0,int_stack+8668,21);
 /*--- compute (gp|gf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+17920,int_stack+14770,int_stack+12520,150);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+0,int_stack+4020,int_stack+3600,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1260,int_stack+4608,int_stack+4020,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+6400,int_stack+1260,int_stack+0,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+8920,int_stack+5392,int_stack+4608,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+24670,int_stack+8920,int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+8920,int_stack+24670,int_stack+6400,28);
 /*--- compute (hp|gf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+24670,int_stack+8920,int_stack+14770,150);
 /*--- compute (gd|gf) ---*/
 hrr1_build_gd(Libint->AB,int_stack+0,int_stack+24670,int_stack+17920,150);
 return int_stack+0;}
