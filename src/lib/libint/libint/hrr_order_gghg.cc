#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_gghg(Libint_t*, prim_data*);

  /* Computes quartets of (gg|hg) integrals */

REALTYPE *hrr_order_gghg(Libint_t *Libint, int num_prim_comb)
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
 Libint->vrr_classes[8][5] = int_stack + 18500;
 Libint->vrr_classes[8][6] = int_stack + 19445;
 Libint->vrr_classes[8][7] = int_stack + 20705;
 Libint->vrr_classes[8][8] = int_stack + 22325;
 Libint->vrr_classes[8][9] = int_stack + 24350;
 memset(int_stack,0,26825*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 26825;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_gghg(Libint, Data);
   Data++;
 }
 /*--- compute (g0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+26825,int_stack+315,int_stack+0,15);
 /*--- compute (g0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+27770,int_stack+735,int_stack+315,15);
 /*--- compute (g0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+29030,int_stack+27770,int_stack+26825,15);
 /*--- compute (g0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+30920,int_stack+1275,int_stack+735,15);
 /*--- compute (g0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+32540,int_stack+30920,int_stack+27770,15);
 /*--- compute (g0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+35060,int_stack+32540,int_stack+29030,15);
 /*--- compute (g0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+26825,int_stack+1950,int_stack+1275,15);
 /*--- compute (g0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+38210,int_stack+26825,int_stack+30920,15);
 /*--- compute (g0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+26825,int_stack+38210,int_stack+32540,15);
 /*--- compute (g0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+38210,int_stack+26825,int_stack+35060,15);
 /*--- compute (h0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+26825,int_stack+3216,int_stack+2775,21);
 /*--- compute (h0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+28148,int_stack+3804,int_stack+3216,21);
 /*--- compute (h0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+29912,int_stack+28148,int_stack+26825,21);
 /*--- compute (h0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+32558,int_stack+4560,int_stack+3804,21);
 /*--- compute (h0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+0,int_stack+32558,int_stack+28148,21);
 /*--- compute (h0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+42935,int_stack+0,int_stack+29912,21);
 /*--- compute (h0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+26825,int_stack+5505,int_stack+4560,21);
 /*--- compute (h0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+47345,int_stack+26825,int_stack+32558,21);
 /*--- compute (h0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+26825,int_stack+47345,int_stack+0,21);
 /*--- compute (h0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+0,int_stack+26825,int_stack+42935,21);
 /*--- compute (gp|hg) ---*/
 hrr1_build_gp(Libint->AB,int_stack+42935,int_stack+0,int_stack+38210,315);
 /*--- compute (i0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+26825,int_stack+7248,int_stack+6660,28);
 /*--- compute (i0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+28589,int_stack+8032,int_stack+7248,28);
 /*--- compute (i0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+30941,int_stack+28589,int_stack+26825,28);
 /*--- compute (i0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+34469,int_stack+9040,int_stack+8032,28);
 /*--- compute (i0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+37493,int_stack+34469,int_stack+28589,28);
 /*--- compute (i0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+57110,int_stack+37493,int_stack+30941,28);
 /*--- compute (i0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+26825,int_stack+10300,int_stack+9040,28);
 /*--- compute (i0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+62990,int_stack+26825,int_stack+34469,28);
 /*--- compute (i0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+26825,int_stack+62990,int_stack+37493,28);
 /*--- compute (i0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+62990,int_stack+26825,int_stack+57110,28);
 /*--- compute (hp|hg) ---*/
 hrr1_build_hp(Libint->AB,int_stack+71810,int_stack+62990,int_stack+0,315);
 /*--- compute (gd|hg) ---*/
 hrr1_build_gd(Libint->AB,int_stack+91655,int_stack+71810,int_stack+42935,315);
 /*--- compute (k0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+0,int_stack+12596,int_stack+11840,36);
 /*--- compute (k0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+2268,int_stack+13604,int_stack+12596,36);
 /*--- compute (k0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+5292,int_stack+2268,int_stack+0,36);
 /*--- compute (k0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+26825,int_stack+14900,int_stack+13604,36);
 /*--- compute (k0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+30713,int_stack+26825,int_stack+2268,36);
 /*--- compute (k0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+36761,int_stack+30713,int_stack+5292,36);
 /*--- compute (k0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+0,int_stack+16520,int_stack+14900,36);
 /*--- compute (k0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+4860,int_stack+0,int_stack+26825,36);
 /*--- compute (k0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+44321,int_stack+4860,int_stack+30713,36);
 /*--- compute (k0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+0,int_stack+44321,int_stack+36761,36);
 /*--- compute (ip|hg) ---*/
 hrr1_build_ip(Libint->AB,int_stack+26825,int_stack+0,int_stack+62990,315);
 /*--- compute (hd|hg) ---*/
 hrr1_build_hd(Libint->AB,int_stack+120005,int_stack+26825,int_stack+71810,315);
 /*--- compute (gf|hg) ---*/
 hrr1_build_gf(Libint->AB,int_stack+159695,int_stack+120005,int_stack+91655,315);
 /*--- compute (l0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+53285,int_stack+19445,int_stack+18500,45);
 /*--- compute (l0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+56120,int_stack+20705,int_stack+19445,45);
 /*--- compute (l0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+59900,int_stack+56120,int_stack+53285,45);
 /*--- compute (l0|kp) ---*/
 hrr3_build_kp(Libint->CD,int_stack+65570,int_stack+22325,int_stack+20705,45);
 /*--- compute (l0|id) ---*/
 hrr3_build_id(Libint->CD,int_stack+70430,int_stack+65570,int_stack+56120,45);
 /*--- compute (l0|hf) ---*/
 hrr3_build_hf(Libint->CD,int_stack+77990,int_stack+70430,int_stack+59900,45);
 /*--- compute (l0|lp) ---*/
 hrr3_build_lp(Libint->CD,int_stack+53285,int_stack+24350,int_stack+22325,45);
 /*--- compute (l0|kd) ---*/
 hrr3_build_kd(Libint->CD,int_stack+87440,int_stack+53285,int_stack+65570,45);
 /*--- compute (l0|if) ---*/
 hrr3_build_if(Libint->CD,int_stack+53285,int_stack+87440,int_stack+70430,45);
 /*--- compute (l0|hg) ---*/
 hrr3_build_hg(Libint->CD,int_stack+87440,int_stack+53285,int_stack+77990,45);
 /*--- compute (kp|hg) ---*/
 hrr1_build_kp(Libint->AB,int_stack+53285,int_stack+87440,int_stack+0,315);
 /*--- compute (id|hg) ---*/
 hrr1_build_id(Libint->AB,int_stack+206945,int_stack+53285,int_stack+26825,315);
 /*--- compute (hf|hg) ---*/
 hrr1_build_hf(Libint->AB,int_stack+0,int_stack+206945,int_stack+120005,315);
 /*--- compute (gg|hg) ---*/
 hrr1_build_gg(Libint->AB,int_stack+66150,int_stack+0,int_stack+159695,315);
 return int_stack+66150;}
