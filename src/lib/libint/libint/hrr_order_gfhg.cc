#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gfhg(Libint_t*, prim_data*);

  /* Computes quartets of (gf|hg) integrals */

REALTYPE *hrr_order_gfhg(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[4][5] = int_stack + 0;
 Libint->vrr_classes[4][6] = int_stack + 315;
 Libint->vrr_classes[4][7] = int_stack + 735;
 Libint->vrr_classes[4][8] = int_stack + 1275;
 Libint->vrr_classes[4][9] = int_stack + 1950;
 Libint->vrr_classes[5][5] = int_stack + 2775;
 Libint->vrr_classes[5][6] = int_stack + 3216;
 Libint->vrr_classes[5][7] = int_stack + 3804;
 Libint->vrr_classes[5][8] = int_stack + 4560;
 Libint->vrr_classes[5][9] = int_stack + 5505;
 Libint->vrr_classes[6][5] = int_stack + 6660;
 Libint->vrr_classes[6][6] = int_stack + 7248;
 Libint->vrr_classes[6][7] = int_stack + 8032;
 Libint->vrr_classes[6][8] = int_stack + 9040;
 Libint->vrr_classes[6][9] = int_stack + 10300;
 Libint->vrr_classes[7][5] = int_stack + 11840;
 Libint->vrr_classes[7][6] = int_stack + 12596;
 Libint->vrr_classes[7][7] = int_stack + 13604;
 Libint->vrr_classes[7][8] = int_stack + 14900;
 Libint->vrr_classes[7][9] = int_stack + 16520;
 memset(int_stack,0,18500*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 18500;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gfhg(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+18500,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+19445,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+20705,int_stack+19445,int_stack+18500,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+22595,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+24215,int_stack+22595,int_stack+19445,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+26735,int_stack+24215,int_stack+20705,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+18500,int_stack+1950,int_stack+1275,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+29885,int_stack+18500,int_stack+22595,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+18500,int_stack+29885,int_stack+24215,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+29885,int_stack+18500,int_stack+26735,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+18500,int_stack+3216,int_stack+2775,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+19823,int_stack+3804,int_stack+3216,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+21587,int_stack+19823,int_stack+18500,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+24233,int_stack+4560,int_stack+3804,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+24233,int_stack+19823,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+34610,int_stack+0,int_stack+21587,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+18500,int_stack+5505,int_stack+4560,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+39020,int_stack+18500,int_stack+24233,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+18500,int_stack+39020,int_stack+0,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+0,int_stack+18500,int_stack+34610,21);
 /*--- compute (gp|hg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+34610,int_stack+0,int_stack+29885,315);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+18500,int_stack+7248,int_stack+6660,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+20264,int_stack+8032,int_stack+7248,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+22616,int_stack+20264,int_stack+18500,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+26144,int_stack+9040,int_stack+8032,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+29168,int_stack+26144,int_stack+20264,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+48785,int_stack+29168,int_stack+22616,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+18500,int_stack+10300,int_stack+9040,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+54665,int_stack+18500,int_stack+26144,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+18500,int_stack+54665,int_stack+29168,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+54665,int_stack+18500,int_stack+48785,28);
 /*--- compute (hp|hg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+63485,int_stack+54665,int_stack+0,315);
 /*--- compute (gd|hg) ---*/
 hrr1_build_gd(Libint->AB,int_stack+83330,int_stack+63485,int_stack+34610,315);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+12596,int_stack+11840,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+13604,int_stack+12596,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+18500,int_stack+14900,int_stack+13604,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+22388,int_stack+18500,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+28436,int_stack+22388,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+16520,int_stack+14900,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+4860,int_stack+0,int_stack+18500,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+35996,int_stack+4860,int_stack+22388,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+0,int_stack+35996,int_stack+28436,36);
 /*--- compute (ip|hg) ---*/
 hrr1_build_ip(Libint->AB,int_stack+11340,int_stack+0,int_stack+54665,315);
 /*--- compute (hd|hg) ---*/
 hrr1_build_hd(Libint->AB,int_stack+111680,int_stack+11340,int_stack+63485,315);
 /*--- compute (gf|hg) ---*/
 hrr1_build_gf(Libint->AB,int_stack+0,int_stack+111680,int_stack+83330,315);
 return int_stack+0;}
