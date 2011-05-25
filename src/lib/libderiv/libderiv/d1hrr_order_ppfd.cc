#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|fd) integrals */

void d1hrr_order_ppfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[1][4][11] = int_stack + 30;
 Libderiv->deriv_classes[1][5][11] = int_stack + 75;
 Libderiv->deriv_classes[2][3][11] = int_stack + 138;
 Libderiv->deriv_classes[2][4][11] = int_stack + 198;
 Libderiv->deriv_classes[2][5][11] = int_stack + 288;
 Libderiv->deriv_classes[1][3][10] = int_stack + 414;
 Libderiv->deriv_classes[1][4][10] = int_stack + 444;
 Libderiv->deriv_classes[1][5][10] = int_stack + 489;
 Libderiv->deriv_classes[2][3][10] = int_stack + 552;
 Libderiv->deriv_classes[2][4][10] = int_stack + 612;
 Libderiv->deriv_classes[2][5][10] = int_stack + 702;
 Libderiv->deriv_classes[1][3][9] = int_stack + 828;
 Libderiv->deriv_classes[1][4][9] = int_stack + 858;
 Libderiv->deriv_classes[1][5][9] = int_stack + 903;
 Libderiv->deriv_classes[2][3][9] = int_stack + 966;
 Libderiv->deriv_classes[2][4][9] = int_stack + 1026;
 Libderiv->deriv_classes[2][5][9] = int_stack + 1116;
 Libderiv->deriv_classes[1][3][8] = int_stack + 1242;
 Libderiv->deriv_classes[1][4][8] = int_stack + 1272;
 Libderiv->deriv_classes[1][5][8] = int_stack + 1317;
 Libderiv->deriv_classes[2][3][8] = int_stack + 1380;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1440;
 Libderiv->deriv_classes[2][5][8] = int_stack + 1530;
 Libderiv->deriv_classes[1][3][7] = int_stack + 1656;
 Libderiv->deriv_classes[1][4][7] = int_stack + 1686;
 Libderiv->deriv_classes[1][5][7] = int_stack + 1731;
 Libderiv->deriv_classes[2][3][7] = int_stack + 1794;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1854;
 Libderiv->deriv_classes[2][5][7] = int_stack + 1944;
 Libderiv->deriv_classes[1][3][6] = int_stack + 2070;
 Libderiv->deriv_classes[1][4][6] = int_stack + 2100;
 Libderiv->deriv_classes[1][5][6] = int_stack + 2145;
 Libderiv->dvrr_classes[2][3] = int_stack + 2208;
 Libderiv->deriv_classes[2][3][6] = int_stack + 2268;
 Libderiv->dvrr_classes[2][4] = int_stack + 2328;
 Libderiv->deriv_classes[2][4][6] = int_stack + 2418;
 Libderiv->deriv_classes[2][5][6] = int_stack + 2508;
 Libderiv->deriv_classes[1][3][2] = int_stack + 2634;
 Libderiv->deriv_classes[1][4][2] = int_stack + 2664;
 Libderiv->deriv_classes[1][5][2] = int_stack + 2709;
 Libderiv->deriv_classes[2][3][2] = int_stack + 2772;
 Libderiv->deriv_classes[2][4][2] = int_stack + 2832;
 Libderiv->deriv_classes[2][5][2] = int_stack + 2922;
 Libderiv->deriv_classes[1][3][1] = int_stack + 3048;
 Libderiv->deriv_classes[1][4][1] = int_stack + 3078;
 Libderiv->deriv_classes[1][5][1] = int_stack + 3123;
 Libderiv->deriv_classes[2][3][1] = int_stack + 3186;
 Libderiv->deriv_classes[2][4][1] = int_stack + 3246;
 Libderiv->deriv_classes[2][5][1] = int_stack + 3336;
 Libderiv->dvrr_classes[1][3] = int_stack + 3462;
 Libderiv->dvrr_classes[1][4] = int_stack + 3492;
 Libderiv->dvrr_classes[1][5] = int_stack + 3537;
 Libderiv->deriv_classes[1][3][0] = int_stack + 3600;
 Libderiv->deriv_classes[1][4][0] = int_stack + 3630;
 Libderiv->deriv_classes[1][5][0] = int_stack + 3675;
 Libderiv->deriv_classes[2][3][0] = int_stack + 3738;
 Libderiv->deriv_classes[2][4][0] = int_stack + 3798;
 Libderiv->deriv_classes[2][5][0] = int_stack + 3888;
 memset(int_stack,0,32112);

 Libderiv->dvrr_stack = int_stack + 7524;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4014,int_stack+3492,int_stack+3462,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4104,int_stack+30,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3462,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4194,int_stack+75,int_stack+30, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3492,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4329,int_stack+4194,int_stack+4104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4014,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4104,int_stack+2328,int_stack+2208,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4509,int_stack+198,int_stack+138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2208,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4689,int_stack+288,int_stack+198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2328,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+4689,int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4104,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4509,int_stack+444,int_stack+414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3462, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4599,int_stack+489,int_stack+444, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3492, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4734,int_stack+4599,int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4014, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4509,int_stack+612,int_stack+552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2208, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4914,int_stack+702,int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2328, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5184,int_stack+4914,int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4104, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4509,int_stack+858,int_stack+828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3462, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4599,int_stack+903,int_stack+858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3492, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4914,int_stack+4599,int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4014, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4509,int_stack+1026,int_stack+966, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2208, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5544,int_stack+1116,int_stack+1026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2328, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+360,int_stack+5544,int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4104, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5094,int_stack+1272,int_stack+1242, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4509,int_stack+1317,int_stack+1272, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5544,int_stack+4509,int_stack+5094, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4509,int_stack+1440,int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5724,int_stack+1530,int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+720,int_stack+5724,int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5094,int_stack+1686,int_stack+1656, 0.0, zero_stack, 1.0, int_stack+3462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4509,int_stack+1731,int_stack+1686, 0.0, zero_stack, 1.0, int_stack+3492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5724,int_stack+4509,int_stack+5094, 0.0, zero_stack, 1.0, int_stack+4014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4509,int_stack+1854,int_stack+1794, 0.0, zero_stack, 1.0, int_stack+2208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1080,int_stack+1944,int_stack+1854, 0.0, zero_stack, 1.0, int_stack+2328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1350,int_stack+1080,int_stack+4509, 0.0, zero_stack, 1.0, int_stack+4104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5094,int_stack+2100,int_stack+2070, 1.0, int_stack+3462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4509,int_stack+2145,int_stack+2100, 1.0, int_stack+3492, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1080,int_stack+4509,int_stack+5094, 1.0, int_stack+4014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4509,int_stack+2418,int_stack+2268, 1.0, int_stack+2208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1710,int_stack+2508,int_stack+2418, 1.0, int_stack+2328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1980,int_stack+1710,int_stack+4509, 1.0, int_stack+4104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4104,int_stack+3537,int_stack+3492,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4509,int_stack+4104,int_stack+4014,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5094,int_stack+2664,int_stack+2634,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4014,int_stack+2709,int_stack+2664,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+4149,int_stack+4014,int_stack+5094,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1710,int_stack+2832,int_stack+2772,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2340,int_stack+2922,int_stack+2832,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2610,int_stack+2340,int_stack+1710,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5094,int_stack+3078,int_stack+3048,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4014,int_stack+3123,int_stack+3078,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1710,int_stack+4014,int_stack+5094,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2340,int_stack+3246,int_stack+3186,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2970,int_stack+3336,int_stack+3246,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3240,int_stack+2970,int_stack+2340,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5094,int_stack+3630,int_stack+3600,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4014,int_stack+3675,int_stack+3630,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2340,int_stack+4014,int_stack+5094,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2970,int_stack+3798,int_stack+3738,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5904,int_stack+3888,int_stack+3798,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3600,int_stack+5904,int_stack+2970,6);
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+5904,int_stack+0,int_stack+4329,60);
     Libderiv->ABCD[11] = int_stack + 5904;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+6444,int_stack+5184,int_stack+4734,60);
     Libderiv->ABCD[10] = int_stack + 6444;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+6984,int_stack+360,int_stack+4914,60);
     Libderiv->ABCD[9] = int_stack + 6984;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+720,int_stack+5544,60);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+540,int_stack+1350,int_stack+5724,60);
     Libderiv->ABCD[7] = int_stack + 540;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4689,int_stack+1980,int_stack+1080,60);
     Libderiv->ABCD[6] = int_stack + 4689;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1080,int_stack+2610,int_stack+4149, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 1080;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3960,int_stack+3240,int_stack+1710, 0.0, zero_stack, 1.0, int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 3960;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1620,int_stack+3600,int_stack+2340, 1.0, int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 1620;

}
