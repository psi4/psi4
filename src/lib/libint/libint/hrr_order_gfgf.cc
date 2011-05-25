#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gfgf(Libint_t*, prim_data*);

  /* Computes quartets of (gf|gf) integrals */

REALTYPE *hrr_order_gfgf(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[7][4] = int_stack + 6400;
 Libint->vrr_classes[7][5] = int_stack + 6940;
 Libint->vrr_classes[7][6] = int_stack + 7696;
 Libint->vrr_classes[7][7] = int_stack + 8704;
 memset(int_stack,0,10000*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 10000;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gfgf(Libint, Data);
   Data++;
 }
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+10000,int_stack+225,int_stack+0,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+10675,int_stack+540,int_stack+225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+11620,int_stack+10675,int_stack+10000,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12970,int_stack+960,int_stack+540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+14230,int_stack+12970,int_stack+10675,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+16120,int_stack+14230,int_stack+11620,15);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+10000,int_stack+1815,int_stack+1500,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+10945,int_stack+2256,int_stack+1815,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+12268,int_stack+10945,int_stack+10000,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+14158,int_stack+2844,int_stack+2256,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+14158,int_stack+10945,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+18370,int_stack+0,int_stack+12268,21);
 /*--- compute (gp|gf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+21520,int_stack+18370,int_stack+16120,150);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+0,int_stack+4020,int_stack+3600,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1260,int_stack+4608,int_stack+4020,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+10000,int_stack+1260,int_stack+0,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+12520,int_stack+5392,int_stack+4608,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+28270,int_stack+12520,int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+12520,int_stack+28270,int_stack+10000,28);
 /*--- compute (hp|gf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+28270,int_stack+12520,int_stack+18370,150);
 /*--- compute (gd|gf) ---*/
 hrr1_build_gd(Libint->AB,int_stack+37720,int_stack+28270,int_stack+21520,150);
 /*--- compute (k0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+10000,int_stack+6940,int_stack+6400,36);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+16720,int_stack+7696,int_stack+6940,36);
 /*--- compute (k0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+18988,int_stack+16720,int_stack+10000,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+22228,int_stack+8704,int_stack+7696,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+0,int_stack+22228,int_stack+16720,36);
 /*--- compute (k0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+22228,int_stack+0,int_stack+18988,36);
 /*--- compute (ip|gf) ---*/
 hrr1_build_ip(Libint->AB,int_stack+51220,int_stack+22228,int_stack+12520,150);
 /*--- compute (hd|gf) ---*/
 hrr1_build_hd(Libint->AB,int_stack+0,int_stack+51220,int_stack+28270,150);
 /*--- compute (gf|gf) ---*/
 hrr1_build_gf(Libint->AB,int_stack+51220,int_stack+0,int_stack+37720,150);
 return int_stack+51220;}
