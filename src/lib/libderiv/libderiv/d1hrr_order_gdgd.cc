#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gdgd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gd|gd) integrals */

void d1hrr_order_gdgd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 540;
 Libderiv->deriv_classes[5][4][11] = int_stack + 960;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1275;
 Libderiv->deriv_classes[5][6][11] = int_stack + 1716;
 Libderiv->deriv_classes[6][4][11] = int_stack + 2304;
 Libderiv->deriv_classes[6][5][11] = int_stack + 2724;
 Libderiv->deriv_classes[6][6][11] = int_stack + 3312;
 Libderiv->deriv_classes[4][4][10] = int_stack + 4096;
 Libderiv->deriv_classes[4][5][10] = int_stack + 4321;
 Libderiv->deriv_classes[4][6][10] = int_stack + 4636;
 Libderiv->deriv_classes[5][4][10] = int_stack + 5056;
 Libderiv->deriv_classes[5][5][10] = int_stack + 5371;
 Libderiv->deriv_classes[5][6][10] = int_stack + 5812;
 Libderiv->deriv_classes[6][4][10] = int_stack + 6400;
 Libderiv->deriv_classes[6][5][10] = int_stack + 6820;
 Libderiv->deriv_classes[6][6][10] = int_stack + 7408;
 Libderiv->deriv_classes[4][4][9] = int_stack + 8192;
 Libderiv->deriv_classes[4][5][9] = int_stack + 8417;
 Libderiv->deriv_classes[4][6][9] = int_stack + 8732;
 Libderiv->deriv_classes[5][4][9] = int_stack + 9152;
 Libderiv->deriv_classes[5][5][9] = int_stack + 9467;
 Libderiv->deriv_classes[5][6][9] = int_stack + 9908;
 Libderiv->deriv_classes[6][4][9] = int_stack + 10496;
 Libderiv->deriv_classes[6][5][9] = int_stack + 10916;
 Libderiv->deriv_classes[6][6][9] = int_stack + 11504;
 Libderiv->deriv_classes[4][4][8] = int_stack + 12288;
 Libderiv->deriv_classes[4][5][8] = int_stack + 12513;
 Libderiv->deriv_classes[4][6][8] = int_stack + 12828;
 Libderiv->deriv_classes[5][4][8] = int_stack + 13248;
 Libderiv->deriv_classes[5][5][8] = int_stack + 13563;
 Libderiv->deriv_classes[5][6][8] = int_stack + 14004;
 Libderiv->deriv_classes[6][4][8] = int_stack + 14592;
 Libderiv->deriv_classes[6][5][8] = int_stack + 15012;
 Libderiv->deriv_classes[6][6][8] = int_stack + 15600;
 Libderiv->deriv_classes[4][4][7] = int_stack + 16384;
 Libderiv->deriv_classes[4][5][7] = int_stack + 16609;
 Libderiv->deriv_classes[4][6][7] = int_stack + 16924;
 Libderiv->deriv_classes[5][4][7] = int_stack + 17344;
 Libderiv->deriv_classes[5][5][7] = int_stack + 17659;
 Libderiv->deriv_classes[5][6][7] = int_stack + 18100;
 Libderiv->deriv_classes[6][4][7] = int_stack + 18688;
 Libderiv->deriv_classes[6][5][7] = int_stack + 19108;
 Libderiv->deriv_classes[6][6][7] = int_stack + 19696;
 Libderiv->deriv_classes[4][4][6] = int_stack + 20480;
 Libderiv->deriv_classes[4][5][6] = int_stack + 20705;
 Libderiv->deriv_classes[4][6][6] = int_stack + 21020;
 Libderiv->deriv_classes[5][4][6] = int_stack + 21440;
 Libderiv->deriv_classes[5][5][6] = int_stack + 21755;
 Libderiv->deriv_classes[5][6][6] = int_stack + 22196;
 Libderiv->dvrr_classes[6][4] = int_stack + 22784;
 Libderiv->deriv_classes[6][4][6] = int_stack + 23204;
 Libderiv->dvrr_classes[6][5] = int_stack + 23624;
 Libderiv->deriv_classes[6][5][6] = int_stack + 24212;
 Libderiv->deriv_classes[6][6][6] = int_stack + 24800;
 Libderiv->deriv_classes[4][4][2] = int_stack + 25584;
 Libderiv->deriv_classes[4][5][2] = int_stack + 25809;
 Libderiv->deriv_classes[4][6][2] = int_stack + 26124;
 Libderiv->deriv_classes[5][4][2] = int_stack + 26544;
 Libderiv->deriv_classes[5][5][2] = int_stack + 26859;
 Libderiv->deriv_classes[5][6][2] = int_stack + 27300;
 Libderiv->deriv_classes[6][4][2] = int_stack + 27888;
 Libderiv->deriv_classes[6][5][2] = int_stack + 28308;
 Libderiv->deriv_classes[6][6][2] = int_stack + 28896;
 Libderiv->deriv_classes[4][4][1] = int_stack + 29680;
 Libderiv->deriv_classes[4][5][1] = int_stack + 29905;
 Libderiv->deriv_classes[4][6][1] = int_stack + 30220;
 Libderiv->deriv_classes[5][4][1] = int_stack + 30640;
 Libderiv->deriv_classes[5][5][1] = int_stack + 30955;
 Libderiv->deriv_classes[5][6][1] = int_stack + 31396;
 Libderiv->deriv_classes[6][4][1] = int_stack + 31984;
 Libderiv->deriv_classes[6][5][1] = int_stack + 32404;
 Libderiv->deriv_classes[6][6][1] = int_stack + 32992;
 Libderiv->dvrr_classes[4][4] = int_stack + 33776;
 Libderiv->dvrr_classes[4][5] = int_stack + 34001;
 Libderiv->dvrr_classes[4][6] = int_stack + 34316;
 Libderiv->deriv_classes[4][4][0] = int_stack + 34736;
 Libderiv->deriv_classes[4][5][0] = int_stack + 34961;
 Libderiv->deriv_classes[4][6][0] = int_stack + 35276;
 Libderiv->dvrr_classes[5][4] = int_stack + 35696;
 Libderiv->dvrr_classes[5][5] = int_stack + 36011;
 Libderiv->dvrr_classes[5][6] = int_stack + 36452;
 Libderiv->deriv_classes[5][4][0] = int_stack + 37040;
 Libderiv->deriv_classes[5][5][0] = int_stack + 37355;
 Libderiv->deriv_classes[5][6][0] = int_stack + 37796;
 Libderiv->deriv_classes[6][4][0] = int_stack + 38384;
 Libderiv->deriv_classes[6][5][0] = int_stack + 38804;
 Libderiv->deriv_classes[6][6][0] = int_stack + 39392;
 memset(int_stack,0,321408);

 Libderiv->dvrr_stack = int_stack + 122598;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gdgd(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+40176,int_stack+34001,int_stack+33776,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+40851,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33776,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41526,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34001,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42471,int_stack+41526,int_stack+40851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40176,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+40851,int_stack+36011,int_stack+35696,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1275,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35696,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43821,int_stack+1716,int_stack+1275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36011,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45144,int_stack+43821,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40851,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+47034,int_stack+45144,int_stack+42471,90);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+23624,int_stack+22784,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41796,int_stack+2724,int_stack+2304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22784,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43056,int_stack+3312,int_stack+2724, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23624,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1260,int_stack+43056,int_stack+41796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+51084,int_stack+1260,int_stack+45144,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+4321,int_stack+4096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33776, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1935,int_stack+4636,int_stack+4321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34001, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2880,int_stack+1935,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40176, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+5371,int_stack+5056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35696, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41796,int_stack+5812,int_stack+5371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36011, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4230,int_stack+41796,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40851, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+41796,int_stack+4230,int_stack+2880,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+6820,int_stack+6400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22784, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+56754,int_stack+7408,int_stack+6820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23624, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+58518,int_stack+56754,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+61038,int_stack+58518,int_stack+4230,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+8417,int_stack+8192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33776, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1935,int_stack+8732,int_stack+8417, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34001, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2880,int_stack+1935,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40176, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+9467,int_stack+9152, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35696, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4230,int_stack+9908,int_stack+9467, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36011, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5553,int_stack+4230,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40851, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+56754,int_stack+5553,int_stack+2880,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+10916,int_stack+10496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22784, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2520,int_stack+11504,int_stack+10916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23624, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7443,int_stack+2520,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+66708,int_stack+7443,int_stack+5553,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+12513,int_stack+12288, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1935,int_stack+12828,int_stack+12513, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34001, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2880,int_stack+1935,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+13563,int_stack+13248, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4230,int_stack+14004,int_stack+13563, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5553,int_stack+4230,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+7443,int_stack+5553,int_stack+2880,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+15012,int_stack+14592, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2520,int_stack+15600,int_stack+15012, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11493,int_stack+2520,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+72378,int_stack+11493,int_stack+5553,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11493,int_stack+16609,int_stack+16384, 0.0, zero_stack, 1.0, int_stack+33776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12168,int_stack+16924,int_stack+16609, 0.0, zero_stack, 1.0, int_stack+34001, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13113,int_stack+12168,int_stack+11493, 0.0, zero_stack, 1.0, int_stack+40176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11493,int_stack+17659,int_stack+17344, 0.0, zero_stack, 1.0, int_stack+35696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14463,int_stack+18100,int_stack+17659, 0.0, zero_stack, 1.0, int_stack+36011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15786,int_stack+14463,int_stack+11493, 0.0, zero_stack, 1.0, int_stack+40851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+1260,int_stack+15786,int_stack+13113,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11493,int_stack+19108,int_stack+18688, 0.0, zero_stack, 1.0, int_stack+22784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12753,int_stack+19696,int_stack+19108, 0.0, zero_stack, 1.0, int_stack+23624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17676,int_stack+12753,int_stack+11493, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+78048,int_stack+17676,int_stack+15786,90);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11493,int_stack+20705,int_stack+20480, 1.0, int_stack+33776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12168,int_stack+21020,int_stack+20705, 1.0, int_stack+34001, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13113,int_stack+12168,int_stack+11493, 1.0, int_stack+40176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11493,int_stack+21755,int_stack+21440, 1.0, int_stack+35696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14463,int_stack+22196,int_stack+21755, 1.0, int_stack+36011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15786,int_stack+14463,int_stack+11493, 1.0, int_stack+40851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+17676,int_stack+15786,int_stack+13113,90);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11493,int_stack+24212,int_stack+23204, 1.0, int_stack+22784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12753,int_stack+24800,int_stack+24212, 1.0, int_stack+23624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21726,int_stack+12753,int_stack+11493, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gd) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+83718,int_stack+21726,int_stack+15786,90);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+21726,int_stack+34316,int_stack+34001,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+22671,int_stack+21726,int_stack+40176,15);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24021,int_stack+36452,int_stack+36011,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11493,int_stack+24021,int_stack+40851,21);
 /*--- compute (gp|gd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+13383,int_stack+11493,int_stack+22671,90);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24021,int_stack+25809,int_stack+25584,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+21726,int_stack+26124,int_stack+25809,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+24696,int_stack+21726,int_stack+24021,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+21726,int_stack+26859,int_stack+26544,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+40176,int_stack+27300,int_stack+26859,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5310,int_stack+40176,int_stack+21726,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+89388,int_stack+5310,int_stack+24696, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22671, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+28308,int_stack+27888,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24021,int_stack+28896,int_stack+28308,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+25785,int_stack+24021,int_stack+0,28);
 /*--- compute (hp|gd) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+93438,int_stack+25785,int_stack+5310, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5310,int_stack+29905,int_stack+29680,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+21726,int_stack+30220,int_stack+29905,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5985,int_stack+21726,int_stack+5310,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+21726,int_stack+30955,int_stack+30640,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24021,int_stack+31396,int_stack+30955,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+25344,int_stack+24021,int_stack+21726,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+27234,int_stack+25344,int_stack+5985, 0.0, zero_stack, 1.0, int_stack+22671, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+32404,int_stack+31984,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+5310,int_stack+32992,int_stack+32404,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+31284,int_stack+5310,int_stack+0,28);
 /*--- compute (hp|gd) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+99108,int_stack+31284,int_stack+25344, 0.0, zero_stack, 1.0, int_stack+11493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+31284,int_stack+34961,int_stack+34736,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+21726,int_stack+35276,int_stack+34961,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+31959,int_stack+21726,int_stack+31284,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+21726,int_stack+37355,int_stack+37040,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+33309,int_stack+37796,int_stack+37355,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+34632,int_stack+33309,int_stack+21726,21);
 /*--- compute (gp|gd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+104778,int_stack+34632,int_stack+31959, 1.0, int_stack+22671, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+38804,int_stack+38384,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+21726,int_stack+39392,int_stack+38804,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23490,int_stack+21726,int_stack+0,28);
 /*--- compute (hp|gd) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+108828,int_stack+23490,int_stack+34632, 1.0, int_stack+11493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+31284,int_stack+51084,int_stack+47034,90);
     Libderiv->ABCD[11] = int_stack + 31284;
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+45846,int_stack+61038,int_stack+41796,90);
     Libderiv->ABCD[10] = int_stack + 45846;
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+114498,int_stack+66708,int_stack+56754,90);
     Libderiv->ABCD[9] = int_stack + 114498;
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+53946,int_stack+72378,int_stack+7443,90);
     Libderiv->ABCD[8] = int_stack + 53946;
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+62046,int_stack+78048,int_stack+1260,90);
     Libderiv->ABCD[7] = int_stack + 62046;
 /*--- compute (gd|gd) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+0,int_stack+83718,int_stack+17676,90);
     Libderiv->ABCD[6] = int_stack + 0;
 /*--- compute (gd|gd) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+17433,int_stack+93438,int_stack+89388, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[2] = int_stack + 17433;
 /*--- compute (gd|gd) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+70146,int_stack+99108,int_stack+27234, 0.0, zero_stack, 1.0, int_stack+13383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[1] = int_stack + 70146;
 /*--- compute (gd|gd) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+78246,int_stack+108828,int_stack+104778, 1.0, int_stack+13383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,90);
     Libderiv->ABCD[0] = int_stack + 78246;

}
