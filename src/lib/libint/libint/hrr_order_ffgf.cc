#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_ffgf(Libint_t*, prim_data*);

  /* Computes quartets of (ff|gf) integrals */

REALTYPE *hrr_order_ffgf(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[6][4] = int_stack + 4600;
 Libint->vrr_classes[6][5] = int_stack + 5020;
 Libint->vrr_classes[6][6] = int_stack + 5608;
 Libint->vrr_classes[6][7] = int_stack + 6392;
 memset(int_stack,0,7400*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 7400;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_ffgf(Libint, Data);
   Data++;
 }
 /*--- compute (f0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+7400,int_stack+150,int_stack+0,10);
 /*--- compute (f0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+7850,int_stack+360,int_stack+150,10);
 /*--- compute (f0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+8480,int_stack+7850,int_stack+7400,10);
 /*--- compute (f0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+9380,int_stack+640,int_stack+360,10);
 /*--- compute (f0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+10220,int_stack+9380,int_stack+7850,10);
 /*--- compute (f0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+11480,int_stack+10220,int_stack+8480,10);
 /*--- compute (g0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+7400,int_stack+1225,int_stack+1000,15);
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8075,int_stack+1540,int_stack+1225,15);
 /*--- compute (g0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+9020,int_stack+8075,int_stack+7400,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+0,int_stack+1960,int_stack+1540,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+12980,int_stack+0,int_stack+8075,15);
 /*--- compute (g0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+0,int_stack+12980,int_stack+9020,15);
 /*--- compute (fp|gf) ---*/
 hrr1_build_fp(Libint->AB,int_stack+12980,int_stack+0,int_stack+11480,150);
 /*--- compute (h0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+7400,int_stack+2815,int_stack+2500,21);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+8345,int_stack+3256,int_stack+2815,21);
 /*--- compute (h0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+9668,int_stack+8345,int_stack+7400,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+17480,int_stack+3844,int_stack+3256,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+19244,int_stack+17480,int_stack+8345,21);
 /*--- compute (h0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+21890,int_stack+19244,int_stack+9668,21);
 /*--- compute (gp|gf) ---*/
 hrr1_build_gp(Libint->AB,int_stack+25040,int_stack+21890,int_stack+0,150);
 /*--- compute (fd|gf) ---*/
 hrr1_build_fd(Libint->AB,int_stack+31790,int_stack+25040,int_stack+12980,150);
 /*--- compute (i0|gp) ---*/
 hrr3_build_gp(Libint->CD,int_stack+0,int_stack+5020,int_stack+4600,28);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+1260,int_stack+5608,int_stack+5020,28);
 /*--- compute (i0|gd) ---*/
 hrr3_build_gd(Libint->CD,int_stack+3024,int_stack+1260,int_stack+0,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+7400,int_stack+6392,int_stack+5608,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+9752,int_stack+7400,int_stack+1260,28);
 /*--- compute (i0|gf) ---*/
 hrr3_build_gf(Libint->CD,int_stack+5544,int_stack+9752,int_stack+3024,28);
 /*--- compute (hp|gf) ---*/
 hrr1_build_hp(Libint->AB,int_stack+9744,int_stack+5544,int_stack+21890,150);
 /*--- compute (gd|gf) ---*/
 hrr1_build_gd(Libint->AB,int_stack+40790,int_stack+9744,int_stack+25040,150);
 /*--- compute (ff|gf) ---*/
 hrr1_build_ff(Libint->AB,int_stack+0,int_stack+40790,int_stack+31790,150);
 return int_stack+0;}
