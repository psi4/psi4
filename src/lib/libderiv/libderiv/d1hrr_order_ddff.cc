#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ddff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|ff) integrals */

void d1hrr_order_ddff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[2][5][11] = int_stack + 150;
 Libderiv->deriv_classes[2][6][11] = int_stack + 276;
 Libderiv->deriv_classes[3][3][11] = int_stack + 444;
 Libderiv->deriv_classes[3][4][11] = int_stack + 544;
 Libderiv->deriv_classes[3][5][11] = int_stack + 694;
 Libderiv->deriv_classes[3][6][11] = int_stack + 904;
 Libderiv->deriv_classes[4][3][11] = int_stack + 1184;
 Libderiv->deriv_classes[4][4][11] = int_stack + 1334;
 Libderiv->deriv_classes[4][5][11] = int_stack + 1559;
 Libderiv->deriv_classes[4][6][11] = int_stack + 1874;
 Libderiv->deriv_classes[2][3][10] = int_stack + 2294;
 Libderiv->deriv_classes[2][4][10] = int_stack + 2354;
 Libderiv->deriv_classes[2][5][10] = int_stack + 2444;
 Libderiv->deriv_classes[2][6][10] = int_stack + 2570;
 Libderiv->deriv_classes[3][3][10] = int_stack + 2738;
 Libderiv->deriv_classes[3][4][10] = int_stack + 2838;
 Libderiv->deriv_classes[3][5][10] = int_stack + 2988;
 Libderiv->deriv_classes[3][6][10] = int_stack + 3198;
 Libderiv->deriv_classes[4][3][10] = int_stack + 3478;
 Libderiv->deriv_classes[4][4][10] = int_stack + 3628;
 Libderiv->deriv_classes[4][5][10] = int_stack + 3853;
 Libderiv->deriv_classes[4][6][10] = int_stack + 4168;
 Libderiv->deriv_classes[2][3][9] = int_stack + 4588;
 Libderiv->deriv_classes[2][4][9] = int_stack + 4648;
 Libderiv->deriv_classes[2][5][9] = int_stack + 4738;
 Libderiv->deriv_classes[2][6][9] = int_stack + 4864;
 Libderiv->deriv_classes[3][3][9] = int_stack + 5032;
 Libderiv->deriv_classes[3][4][9] = int_stack + 5132;
 Libderiv->deriv_classes[3][5][9] = int_stack + 5282;
 Libderiv->deriv_classes[3][6][9] = int_stack + 5492;
 Libderiv->deriv_classes[4][3][9] = int_stack + 5772;
 Libderiv->deriv_classes[4][4][9] = int_stack + 5922;
 Libderiv->deriv_classes[4][5][9] = int_stack + 6147;
 Libderiv->deriv_classes[4][6][9] = int_stack + 6462;
 Libderiv->deriv_classes[2][3][8] = int_stack + 6882;
 Libderiv->deriv_classes[2][4][8] = int_stack + 6942;
 Libderiv->deriv_classes[2][5][8] = int_stack + 7032;
 Libderiv->deriv_classes[2][6][8] = int_stack + 7158;
 Libderiv->deriv_classes[3][3][8] = int_stack + 7326;
 Libderiv->deriv_classes[3][4][8] = int_stack + 7426;
 Libderiv->deriv_classes[3][5][8] = int_stack + 7576;
 Libderiv->deriv_classes[3][6][8] = int_stack + 7786;
 Libderiv->deriv_classes[4][3][8] = int_stack + 8066;
 Libderiv->deriv_classes[4][4][8] = int_stack + 8216;
 Libderiv->deriv_classes[4][5][8] = int_stack + 8441;
 Libderiv->deriv_classes[4][6][8] = int_stack + 8756;
 Libderiv->deriv_classes[2][3][7] = int_stack + 9176;
 Libderiv->deriv_classes[2][4][7] = int_stack + 9236;
 Libderiv->deriv_classes[2][5][7] = int_stack + 9326;
 Libderiv->deriv_classes[2][6][7] = int_stack + 9452;
 Libderiv->deriv_classes[3][3][7] = int_stack + 9620;
 Libderiv->deriv_classes[3][4][7] = int_stack + 9720;
 Libderiv->deriv_classes[3][5][7] = int_stack + 9870;
 Libderiv->deriv_classes[3][6][7] = int_stack + 10080;
 Libderiv->deriv_classes[4][3][7] = int_stack + 10360;
 Libderiv->deriv_classes[4][4][7] = int_stack + 10510;
 Libderiv->deriv_classes[4][5][7] = int_stack + 10735;
 Libderiv->deriv_classes[4][6][7] = int_stack + 11050;
 Libderiv->deriv_classes[2][3][6] = int_stack + 11470;
 Libderiv->deriv_classes[2][4][6] = int_stack + 11530;
 Libderiv->deriv_classes[2][5][6] = int_stack + 11620;
 Libderiv->deriv_classes[2][6][6] = int_stack + 11746;
 Libderiv->deriv_classes[3][3][6] = int_stack + 11914;
 Libderiv->deriv_classes[3][4][6] = int_stack + 12014;
 Libderiv->deriv_classes[3][5][6] = int_stack + 12164;
 Libderiv->deriv_classes[3][6][6] = int_stack + 12374;
 Libderiv->dvrr_classes[4][3] = int_stack + 12654;
 Libderiv->deriv_classes[4][3][6] = int_stack + 12804;
 Libderiv->dvrr_classes[4][4] = int_stack + 12954;
 Libderiv->deriv_classes[4][4][6] = int_stack + 13179;
 Libderiv->dvrr_classes[4][5] = int_stack + 13404;
 Libderiv->deriv_classes[4][5][6] = int_stack + 13719;
 Libderiv->deriv_classes[4][6][6] = int_stack + 14034;
 Libderiv->deriv_classes[2][3][2] = int_stack + 14454;
 Libderiv->deriv_classes[2][4][2] = int_stack + 14514;
 Libderiv->deriv_classes[2][5][2] = int_stack + 14604;
 Libderiv->deriv_classes[2][6][2] = int_stack + 14730;
 Libderiv->deriv_classes[3][3][2] = int_stack + 14898;
 Libderiv->deriv_classes[3][4][2] = int_stack + 14998;
 Libderiv->deriv_classes[3][5][2] = int_stack + 15148;
 Libderiv->deriv_classes[3][6][2] = int_stack + 15358;
 Libderiv->deriv_classes[4][3][2] = int_stack + 15638;
 Libderiv->deriv_classes[4][4][2] = int_stack + 15788;
 Libderiv->deriv_classes[4][5][2] = int_stack + 16013;
 Libderiv->deriv_classes[4][6][2] = int_stack + 16328;
 Libderiv->deriv_classes[2][3][1] = int_stack + 16748;
 Libderiv->deriv_classes[2][4][1] = int_stack + 16808;
 Libderiv->deriv_classes[2][5][1] = int_stack + 16898;
 Libderiv->deriv_classes[2][6][1] = int_stack + 17024;
 Libderiv->deriv_classes[3][3][1] = int_stack + 17192;
 Libderiv->deriv_classes[3][4][1] = int_stack + 17292;
 Libderiv->deriv_classes[3][5][1] = int_stack + 17442;
 Libderiv->deriv_classes[3][6][1] = int_stack + 17652;
 Libderiv->deriv_classes[4][3][1] = int_stack + 17932;
 Libderiv->deriv_classes[4][4][1] = int_stack + 18082;
 Libderiv->deriv_classes[4][5][1] = int_stack + 18307;
 Libderiv->deriv_classes[4][6][1] = int_stack + 18622;
 Libderiv->dvrr_classes[2][3] = int_stack + 19042;
 Libderiv->dvrr_classes[2][4] = int_stack + 19102;
 Libderiv->dvrr_classes[2][5] = int_stack + 19192;
 Libderiv->dvrr_classes[2][6] = int_stack + 19318;
 Libderiv->deriv_classes[2][3][0] = int_stack + 19486;
 Libderiv->deriv_classes[2][4][0] = int_stack + 19546;
 Libderiv->deriv_classes[2][5][0] = int_stack + 19636;
 Libderiv->deriv_classes[2][6][0] = int_stack + 19762;
 Libderiv->dvrr_classes[3][3] = int_stack + 19930;
 Libderiv->dvrr_classes[3][4] = int_stack + 20030;
 Libderiv->dvrr_classes[3][5] = int_stack + 20180;
 Libderiv->dvrr_classes[3][6] = int_stack + 20390;
 Libderiv->deriv_classes[3][3][0] = int_stack + 20670;
 Libderiv->deriv_classes[3][4][0] = int_stack + 20770;
 Libderiv->deriv_classes[3][5][0] = int_stack + 20920;
 Libderiv->deriv_classes[3][6][0] = int_stack + 21130;
 Libderiv->deriv_classes[4][3][0] = int_stack + 21410;
 Libderiv->deriv_classes[4][4][0] = int_stack + 21560;
 Libderiv->deriv_classes[4][5][0] = int_stack + 21785;
 Libderiv->deriv_classes[4][6][0] = int_stack + 22100;
 memset(int_stack,0,180160);

 Libderiv->dvrr_stack = int_stack + 59183;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ddff(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+22520,int_stack+19102,int_stack+19042,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22700,int_stack+19192,int_stack+19102,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+22970,int_stack+22700,int_stack+22520,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23330,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19042,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23510,int_stack+150,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19102,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23780,int_stack+23510,int_stack+23330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22520,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24140,int_stack+276,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19192,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24518,int_stack+24140,int_stack+23510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22700,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+25058,int_stack+24518,int_stack+23780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22970,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+23330,int_stack+20030,int_stack+19930,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23630,int_stack+20180,int_stack+20030,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+24080,int_stack+23630,int_stack+23330,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24680,int_stack+544,int_stack+444, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19930,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+694,int_stack+544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20030,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25658,int_stack+0,int_stack+24680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23330,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26258,int_stack+904,int_stack+694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20180,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+26888,int_stack+26258,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23630,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+26888,int_stack+25658, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24080,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+25658,int_stack+0,int_stack+25058,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27458,int_stack+12954,int_stack+12654,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24680,int_stack+13404,int_stack+12954,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+27908,int_stack+24680,int_stack+27458,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28808,int_stack+1334,int_stack+1184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12654,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29258,int_stack+1559,int_stack+1334, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12954,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29933,int_stack+29258,int_stack+28808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27458,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30833,int_stack+1874,int_stack+1559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13404,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31778,int_stack+30833,int_stack+29258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24680,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+33128,int_stack+31778,int_stack+29933, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27908,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28808,int_stack+33128,int_stack+0,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+2354,int_stack+2294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19042, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180,int_stack+2444,int_stack+2354, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19102, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+450,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22520, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+810,int_stack+2570,int_stack+2444, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19192, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1188,int_stack+810,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22700, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1728,int_stack+1188,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22970, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+2838,int_stack+2738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19930, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+300,int_stack+2988,int_stack+2838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20030, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+750,int_stack+300,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23330, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2328,int_stack+3198,int_stack+2988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20180, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31808,int_stack+2328,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23630, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2328,int_stack+31808,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24080, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+31808,int_stack+2328,int_stack+1728,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33608,int_stack+3628,int_stack+3478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12654, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34058,int_stack+3853,int_stack+3628, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12954, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+34058,int_stack+33608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27458, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+900,int_stack+4168,int_stack+3853, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13404, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+34733,int_stack+900,int_stack+34058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24680, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+36083,int_stack+34733,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27908, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+37583,int_stack+36083,int_stack+2328,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+4648,int_stack+4588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19042, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180,int_stack+4738,int_stack+4648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19102, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+450,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22520, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+810,int_stack+4864,int_stack+4738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19192, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1188,int_stack+810,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22700, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1728,int_stack+1188,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22970, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+5132,int_stack+5032, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19930, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+300,int_stack+5282,int_stack+5132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20030, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+750,int_stack+300,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23330, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2328,int_stack+5492,int_stack+5282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20180, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2958,int_stack+2328,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23630, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3858,int_stack+2958,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24080, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+33608,int_stack+3858,int_stack+1728,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+5922,int_stack+5772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12654, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+6147,int_stack+5922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12954, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1125,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27458, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2025,int_stack+6462,int_stack+6147, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13404, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4858,int_stack+2025,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24680, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2025,int_stack+4858,int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27908, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+40583,int_stack+2025,int_stack+3858,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+6942,int_stack+6882, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180,int_stack+7032,int_stack+6942, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+450,int_stack+180,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+810,int_stack+7158,int_stack+7032, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1188,int_stack+810,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1728,int_stack+1188,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+7426,int_stack+7326, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+300,int_stack+7576,int_stack+7426, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+750,int_stack+300,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2328,int_stack+7786,int_stack+7576, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2958,int_stack+2328,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3858,int_stack+2958,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4858,int_stack+3858,int_stack+1728,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+8216,int_stack+8066, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+8441,int_stack+8216, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12954, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1125,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2025,int_stack+8756,int_stack+8441, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6658,int_stack+2025,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2025,int_stack+6658,int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+43583,int_stack+2025,int_stack+3858,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6658,int_stack+9236,int_stack+9176, 0.0, zero_stack, 1.0, int_stack+19042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6838,int_stack+9326,int_stack+9236, 0.0, zero_stack, 1.0, int_stack+19102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7108,int_stack+6838,int_stack+6658, 0.0, zero_stack, 1.0, int_stack+22520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7468,int_stack+9452,int_stack+9326, 0.0, zero_stack, 1.0, int_stack+19192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7846,int_stack+7468,int_stack+6838, 0.0, zero_stack, 1.0, int_stack+22700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8386,int_stack+7846,int_stack+7108, 0.0, zero_stack, 1.0, int_stack+22970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6658,int_stack+9720,int_stack+9620, 0.0, zero_stack, 1.0, int_stack+19930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6958,int_stack+9870,int_stack+9720, 0.0, zero_stack, 1.0, int_stack+20030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7408,int_stack+6958,int_stack+6658, 0.0, zero_stack, 1.0, int_stack+23330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8986,int_stack+10080,int_stack+9870, 0.0, zero_stack, 1.0, int_stack+20180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+8986,int_stack+6958, 0.0, zero_stack, 1.0, int_stack+23630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8986,int_stack+0,int_stack+7408, 0.0, zero_stack, 1.0, int_stack+24080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+8986,int_stack+8386,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+10510,int_stack+10360, 0.0, zero_stack, 1.0, int_stack+12654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2250,int_stack+10735,int_stack+10510, 0.0, zero_stack, 1.0, int_stack+12954, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2925,int_stack+2250,int_stack+1800, 0.0, zero_stack, 1.0, int_stack+27458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3825,int_stack+11050,int_stack+10735, 0.0, zero_stack, 1.0, int_stack+13404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9986,int_stack+3825,int_stack+2250, 0.0, zero_stack, 1.0, int_stack+24680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6658,int_stack+9986,int_stack+2925, 0.0, zero_stack, 1.0, int_stack+27908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1800,int_stack+6658,int_stack+8986,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6658,int_stack+11530,int_stack+11470, 1.0, int_stack+19042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6838,int_stack+11620,int_stack+11530, 1.0, int_stack+19102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7108,int_stack+6838,int_stack+6658, 1.0, int_stack+22520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7468,int_stack+11746,int_stack+11620, 1.0, int_stack+19192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7846,int_stack+7468,int_stack+6838, 1.0, int_stack+22700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8386,int_stack+7846,int_stack+7108, 1.0, int_stack+22970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6658,int_stack+12014,int_stack+11914, 1.0, int_stack+19930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6958,int_stack+12164,int_stack+12014, 1.0, int_stack+20030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7408,int_stack+6958,int_stack+6658, 1.0, int_stack+23330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8986,int_stack+12374,int_stack+12164, 1.0, int_stack+20180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9616,int_stack+8986,int_stack+6958, 1.0, int_stack+23630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10516,int_stack+9616,int_stack+7408, 1.0, int_stack+24080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+35408,int_stack+10516,int_stack+8386,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6658,int_stack+13179,int_stack+12804, 1.0, int_stack+12654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7108,int_stack+13719,int_stack+13179, 1.0, int_stack+12954, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7783,int_stack+7108,int_stack+6658, 1.0, int_stack+27458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8683,int_stack+14034,int_stack+13719, 1.0, int_stack+13404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11516,int_stack+8683,int_stack+7108, 1.0, int_stack+24680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8683,int_stack+11516,int_stack+7783, 1.0, int_stack+27908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+46583,int_stack+8683,int_stack+10516,100);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24680,int_stack+19318,int_stack+19192,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+25058,int_stack+24680,int_stack+22700,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+27458,int_stack+25058,int_stack+22970,6);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24680,int_stack+20390,int_stack+20180,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6658,int_stack+24680,int_stack+23630,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+7558,int_stack+6658,int_stack+24080,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8558,int_stack+7558,int_stack+27458,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6658,int_stack+14514,int_stack+14454,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6838,int_stack+14604,int_stack+14514,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7108,int_stack+6838,int_stack+6658,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+28058,int_stack+14730,int_stack+14604,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10358,int_stack+28058,int_stack+6838,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+28058,int_stack+10358,int_stack+7108,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10358,int_stack+14998,int_stack+14898,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10658,int_stack+15148,int_stack+14998,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11108,int_stack+10658,int_stack+10358,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11708,int_stack+15358,int_stack+15148,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6658,int_stack+11708,int_stack+10658,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+11708,int_stack+6658,int_stack+11108,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+12708,int_stack+11708,int_stack+28058, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+28058,int_stack+15788,int_stack+15638,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6658,int_stack+16013,int_stack+15788,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14508,int_stack+6658,int_stack+28058,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10358,int_stack+16328,int_stack+16013,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+22520,int_stack+10358,int_stack+6658,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+23870,int_stack+22520,int_stack+14508,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+49583,int_stack+23870,int_stack+11708, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14508,int_stack+16808,int_stack+16748,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14688,int_stack+16898,int_stack+16808,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14958,int_stack+14688,int_stack+14508,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15318,int_stack+17024,int_stack+16898,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15696,int_stack+15318,int_stack+14688,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+16236,int_stack+15696,int_stack+14958,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14508,int_stack+17292,int_stack+17192,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14808,int_stack+17442,int_stack+17292,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+15258,int_stack+14808,int_stack+14508,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+22520,int_stack+17652,int_stack+17442,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6658,int_stack+22520,int_stack+14808,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+22520,int_stack+6658,int_stack+15258,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+23520,int_stack+22520,int_stack+16236, 0.0, zero_stack, 1.0, int_stack+27458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6658,int_stack+18082,int_stack+17932,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14508,int_stack+18307,int_stack+18082,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+15183,int_stack+14508,int_stack+6658,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+16083,int_stack+18622,int_stack+18307,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17028,int_stack+16083,int_stack+14508,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+10358,int_stack+17028,int_stack+15183,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+14508,int_stack+10358,int_stack+22520, 0.0, zero_stack, 1.0, int_stack+7558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+22520,int_stack+19546,int_stack+19486,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22700,int_stack+19636,int_stack+19546,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+22970,int_stack+22700,int_stack+22520,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10358,int_stack+19762,int_stack+19636,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10736,int_stack+10358,int_stack+22700,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+11276,int_stack+10736,int_stack+22970,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10358,int_stack+20770,int_stack+20670,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10658,int_stack+20920,int_stack+20770,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+22520,int_stack+10658,int_stack+10358,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11876,int_stack+21130,int_stack+20920,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6658,int_stack+11876,int_stack+10658,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+17508,int_stack+6658,int_stack+22520,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+18508,int_stack+17508,int_stack+11276, 1.0, int_stack+27458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27458,int_stack+21560,int_stack+21410,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27908,int_stack+21785,int_stack+21560,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+6658,int_stack+27908,int_stack+27458,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+22520,int_stack+22100,int_stack+21785,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10358,int_stack+22520,int_stack+27908,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+20308,int_stack+10358,int_stack+6658,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+52583,int_stack+20308,int_stack+17508, 1.0, int_stack+7558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+55583,int_stack+28808,int_stack+25658,100);
     Libderiv->ABCD[11] = int_stack + 55583;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+25320,int_stack+37583,int_stack+31808,100);
     Libderiv->ABCD[10] = int_stack + 25320;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+28920,int_stack+40583,int_stack+33608,100);
     Libderiv->ABCD[9] = int_stack + 28920;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+37208,int_stack+43583,int_stack+4858,100);
     Libderiv->ABCD[8] = int_stack + 37208;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+4800,int_stack+1800,int_stack+0,100);
     Libderiv->ABCD[7] = int_stack + 4800;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+46583,int_stack+35408,100);
     Libderiv->ABCD[6] = int_stack + 0;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+32520,int_stack+49583,int_stack+12708, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 32520;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+10358,int_stack+14508,int_stack+23520, 0.0, zero_stack, 1.0, int_stack+8558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 10358;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+20308,int_stack+52583,int_stack+18508, 1.0, int_stack+8558, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 20308;

}
