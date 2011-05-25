#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fdgd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fd|gd) integrals */

void d1hrr_order_fdgd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][6][11] = int_stack + 360;
 Libderiv->deriv_classes[4][4][11] = int_stack + 640;
 Libderiv->deriv_classes[4][5][11] = int_stack + 865;
 Libderiv->deriv_classes[4][6][11] = int_stack + 1180;
 Libderiv->deriv_classes[5][4][11] = int_stack + 1600;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1915;
 Libderiv->deriv_classes[5][6][11] = int_stack + 2356;
 Libderiv->deriv_classes[3][4][10] = int_stack + 2944;
 Libderiv->deriv_classes[3][5][10] = int_stack + 3094;
 Libderiv->deriv_classes[3][6][10] = int_stack + 3304;
 Libderiv->deriv_classes[4][4][10] = int_stack + 3584;
 Libderiv->deriv_classes[4][5][10] = int_stack + 3809;
 Libderiv->deriv_classes[4][6][10] = int_stack + 4124;
 Libderiv->deriv_classes[5][4][10] = int_stack + 4544;
 Libderiv->deriv_classes[5][5][10] = int_stack + 4859;
 Libderiv->deriv_classes[5][6][10] = int_stack + 5300;
 Libderiv->deriv_classes[3][4][9] = int_stack + 5888;
 Libderiv->deriv_classes[3][5][9] = int_stack + 6038;
 Libderiv->deriv_classes[3][6][9] = int_stack + 6248;
 Libderiv->deriv_classes[4][4][9] = int_stack + 6528;
 Libderiv->deriv_classes[4][5][9] = int_stack + 6753;
 Libderiv->deriv_classes[4][6][9] = int_stack + 7068;
 Libderiv->deriv_classes[5][4][9] = int_stack + 7488;
 Libderiv->deriv_classes[5][5][9] = int_stack + 7803;
 Libderiv->deriv_classes[5][6][9] = int_stack + 8244;
 Libderiv->deriv_classes[3][4][8] = int_stack + 8832;
 Libderiv->deriv_classes[3][5][8] = int_stack + 8982;
 Libderiv->deriv_classes[3][6][8] = int_stack + 9192;
 Libderiv->deriv_classes[4][4][8] = int_stack + 9472;
 Libderiv->deriv_classes[4][5][8] = int_stack + 9697;
 Libderiv->deriv_classes[4][6][8] = int_stack + 10012;
 Libderiv->deriv_classes[5][4][8] = int_stack + 10432;
 Libderiv->deriv_classes[5][5][8] = int_stack + 10747;
 Libderiv->deriv_classes[5][6][8] = int_stack + 11188;
 Libderiv->deriv_classes[3][4][7] = int_stack + 11776;
 Libderiv->deriv_classes[3][5][7] = int_stack + 11926;
 Libderiv->deriv_classes[3][6][7] = int_stack + 12136;
 Libderiv->deriv_classes[4][4][7] = int_stack + 12416;
 Libderiv->deriv_classes[4][5][7] = int_stack + 12641;
 Libderiv->deriv_classes[4][6][7] = int_stack + 12956;
 Libderiv->deriv_classes[5][4][7] = int_stack + 13376;
 Libderiv->deriv_classes[5][5][7] = int_stack + 13691;
 Libderiv->deriv_classes[5][6][7] = int_stack + 14132;
 Libderiv->deriv_classes[3][4][6] = int_stack + 14720;
 Libderiv->deriv_classes[3][5][6] = int_stack + 14870;
 Libderiv->deriv_classes[3][6][6] = int_stack + 15080;
 Libderiv->deriv_classes[4][4][6] = int_stack + 15360;
 Libderiv->deriv_classes[4][5][6] = int_stack + 15585;
 Libderiv->deriv_classes[4][6][6] = int_stack + 15900;
 Libderiv->dvrr_classes[5][4] = int_stack + 16320;
 Libderiv->deriv_classes[5][4][6] = int_stack + 16635;
 Libderiv->dvrr_classes[5][5] = int_stack + 16950;
 Libderiv->deriv_classes[5][5][6] = int_stack + 17391;
 Libderiv->deriv_classes[5][6][6] = int_stack + 17832;
 Libderiv->deriv_classes[3][4][2] = int_stack + 18420;
 Libderiv->deriv_classes[3][5][2] = int_stack + 18570;
 Libderiv->deriv_classes[3][6][2] = int_stack + 18780;
 Libderiv->deriv_classes[4][4][2] = int_stack + 19060;
 Libderiv->deriv_classes[4][5][2] = int_stack + 19285;
 Libderiv->deriv_classes[4][6][2] = int_stack + 19600;
 Libderiv->deriv_classes[5][4][2] = int_stack + 20020;
 Libderiv->deriv_classes[5][5][2] = int_stack + 20335;
 Libderiv->deriv_classes[5][6][2] = int_stack + 20776;
 Libderiv->deriv_classes[3][4][1] = int_stack + 21364;
 Libderiv->deriv_classes[3][5][1] = int_stack + 21514;
 Libderiv->deriv_classes[3][6][1] = int_stack + 21724;
 Libderiv->deriv_classes[4][4][1] = int_stack + 22004;
 Libderiv->deriv_classes[4][5][1] = int_stack + 22229;
 Libderiv->deriv_classes[4][6][1] = int_stack + 22544;
 Libderiv->deriv_classes[5][4][1] = int_stack + 22964;
 Libderiv->deriv_classes[5][5][1] = int_stack + 23279;
 Libderiv->deriv_classes[5][6][1] = int_stack + 23720;
 Libderiv->dvrr_classes[3][4] = int_stack + 24308;
 Libderiv->dvrr_classes[3][5] = int_stack + 24458;
 Libderiv->dvrr_classes[3][6] = int_stack + 24668;
 Libderiv->deriv_classes[3][4][0] = int_stack + 24948;
 Libderiv->deriv_classes[3][5][0] = int_stack + 25098;
 Libderiv->deriv_classes[3][6][0] = int_stack + 25308;
 Libderiv->dvrr_classes[4][4] = int_stack + 25588;
 Libderiv->dvrr_classes[4][5] = int_stack + 25813;
 Libderiv->dvrr_classes[4][6] = int_stack + 26128;
 Libderiv->deriv_classes[4][4][0] = int_stack + 26548;
 Libderiv->deriv_classes[4][5][0] = int_stack + 26773;
 Libderiv->deriv_classes[4][6][0] = int_stack + 27088;
 Libderiv->deriv_classes[5][4][0] = int_stack + 27508;
 Libderiv->deriv_classes[5][5][0] = int_stack + 27823;
 Libderiv->deriv_classes[5][6][0] = int_stack + 28264;
 memset(int_stack,0,230816);

 Libderiv->dvrr_stack = int_stack + 82420;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fdgd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28852,int_stack+24458,int_stack+24308,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29302,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24308,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29752,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24458,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+30382,int_stack+29752,int_stack+29302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28852,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+29302,int_stack+25813,int_stack+25588,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31282,int_stack+865,int_stack+640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25588,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+31957,int_stack+1180,int_stack+865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25813,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+31957,int_stack+31282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29302,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+31282,int_stack+0,int_stack+30382,90);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+29977,int_stack+16950,int_stack+16320,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33982,int_stack+1915,int_stack+1600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16320,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+34927,int_stack+2356,int_stack+1915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16950,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+36250,int_stack+34927,int_stack+33982, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29977,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+38140,int_stack+36250,int_stack+0,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3094,int_stack+2944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24308, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+450,int_stack+3304,int_stack+3094, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24458, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28852, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+3809,int_stack+3584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25588, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1980,int_stack+4124,int_stack+3809, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25813, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2925,int_stack+1980,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29302, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+33982,int_stack+2925,int_stack+1080,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+4859,int_stack+4544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16320, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+945,int_stack+5300,int_stack+4859, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16950, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42190,int_stack+945,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29977, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+44080,int_stack+42190,int_stack+2925,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42190,int_stack+6038,int_stack+5888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24308, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42640,int_stack+6248,int_stack+6038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24458, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+42640,int_stack+42190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28852, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42190,int_stack+6753,int_stack+6528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25588, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42865,int_stack+7068,int_stack+6753, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25813, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+900,int_stack+42865,int_stack+42190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29302, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+2250,int_stack+900,int_stack+0,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42190,int_stack+7803,int_stack+7488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16320, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4950,int_stack+8244,int_stack+7803, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16950, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6273,int_stack+4950,int_stack+42190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29977, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+48130,int_stack+6273,int_stack+900,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42190,int_stack+8982,int_stack+8832, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42640,int_stack+9192,int_stack+8982, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4950,int_stack+42640,int_stack+42190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42190,int_stack+9697,int_stack+9472, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42865,int_stack+10012,int_stack+9697, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25813, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5850,int_stack+42865,int_stack+42190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+7200,int_stack+5850,int_stack+4950,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42190,int_stack+10747,int_stack+10432, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+11188,int_stack+10747, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+52180,int_stack+0,int_stack+42190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29977, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+54070,int_stack+52180,int_stack+5850,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52180,int_stack+11926,int_stack+11776, 0.0, zero_stack, 1.0, int_stack+24308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52630,int_stack+12136,int_stack+11926, 0.0, zero_stack, 1.0, int_stack+24458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42190,int_stack+52630,int_stack+52180, 0.0, zero_stack, 1.0, int_stack+28852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52180,int_stack+12641,int_stack+12416, 0.0, zero_stack, 1.0, int_stack+25588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52855,int_stack+12956,int_stack+12641, 0.0, zero_stack, 1.0, int_stack+25813, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+52855,int_stack+52180, 0.0, zero_stack, 1.0, int_stack+29302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+9900,int_stack+0,int_stack+42190,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42190,int_stack+13691,int_stack+13376, 0.0, zero_stack, 1.0, int_stack+16320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52180,int_stack+14132,int_stack+13691, 0.0, zero_stack, 1.0, int_stack+16950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4950,int_stack+52180,int_stack+42190, 0.0, zero_stack, 1.0, int_stack+29977, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+58120,int_stack+4950,int_stack+0,90);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+14870,int_stack+14720, 1.0, int_stack+24308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+450,int_stack+15080,int_stack+14870, 1.0, int_stack+24458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+450,int_stack+0, 1.0, int_stack+28852, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+15585,int_stack+15360, 1.0, int_stack+25588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4950,int_stack+15900,int_stack+15585, 1.0, int_stack+25813, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42190,int_stack+4950,int_stack+0, 1.0, int_stack+29302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12600,int_stack+42190,int_stack+1080,90);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+17391,int_stack+16635, 1.0, int_stack+16320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4950,int_stack+17832,int_stack+17391, 1.0, int_stack+16950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+52180,int_stack+4950,int_stack+0, 1.0, int_stack+29977, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+62170,int_stack+52180,int_stack+42190,90);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+42190,int_stack+24668,int_stack+24458,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+42820,int_stack+42190,int_stack+28852,10);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52180,int_stack+26128,int_stack+25813,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+52180,int_stack+29302,15);
 /*--- compute (fp|gd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15300,int_stack+0,int_stack+42820,90);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52180,int_stack+18570,int_stack+18420,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+42190,int_stack+18780,int_stack+18570,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1350,int_stack+42190,int_stack+52180,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52180,int_stack+19285,int_stack+19060,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52855,int_stack+19600,int_stack+19285,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+28852,int_stack+52855,int_stack+52180,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+66220,int_stack+28852,int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52180,int_stack+20335,int_stack+20020,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4950,int_stack+20776,int_stack+20335,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18000,int_stack+4950,int_stack+52180,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+68920,int_stack+18000,int_stack+28852, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28852,int_stack+21514,int_stack+21364,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+42190,int_stack+21724,int_stack+21514,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1350,int_stack+42190,int_stack+28852,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28852,int_stack+22229,int_stack+22004,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+29527,int_stack+22544,int_stack+22229,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18000,int_stack+29527,int_stack+28852,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+19350,int_stack+18000,int_stack+1350, 0.0, zero_stack, 1.0, int_stack+42820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28852,int_stack+23279,int_stack+22964,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+29797,int_stack+23720,int_stack+23279,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+52180,int_stack+29797,int_stack+28852,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+72970,int_stack+52180,int_stack+18000, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18000,int_stack+25098,int_stack+24948,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+42190,int_stack+25308,int_stack+25098,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18450,int_stack+42190,int_stack+18000,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52180,int_stack+26773,int_stack+26548,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52855,int_stack+27088,int_stack+26773,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+28852,int_stack+52855,int_stack+52180,15);
 /*--- compute (fp|gd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+22050,int_stack+28852,int_stack+18450, 1.0, int_stack+42820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52180,int_stack+27823,int_stack+27508,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+18000,int_stack+28264,int_stack+27823,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+42190,int_stack+18000,int_stack+52180,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+24750,int_stack+42190,int_stack+28852, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+77020,int_stack+38140,int_stack+31282,90);
     Libderiv->ABCD[11] = int_stack + 77020;
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+36682,int_stack+44080,int_stack+33982,90);
     Libderiv->ABCD[10] = int_stack + 36682;
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+28800,int_stack+48130,int_stack+2250,90);
     Libderiv->ABCD[9] = int_stack + 28800;
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+0,int_stack+54070,int_stack+7200,90);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+42082,int_stack+58120,int_stack+9900,90);
     Libderiv->ABCD[7] = int_stack + 42082;
 /*--- compute (fd|gd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+5400,int_stack+62170,int_stack+12600,90);
     Libderiv->ABCD[6] = int_stack + 5400;
 /*--- compute (fd|gd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+47482,int_stack+68920,int_stack+66220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[2] = int_stack + 47482;
 /*--- compute (fd|gd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+52882,int_stack+72970,int_stack+19350, 0.0, zero_stack, 1.0, int_stack+15300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[1] = int_stack + 52882;
 /*--- compute (fd|gd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+58282,int_stack+24750,int_stack+22050, 1.0, int_stack+15300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[0] = int_stack + 58282;

}
