#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fdfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fd|fd) integrals */

void d1hrr_order_fdfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 100;
 Libderiv->deriv_classes[3][5][11] = int_stack + 250;
 Libderiv->deriv_classes[4][3][11] = int_stack + 460;
 Libderiv->deriv_classes[4][4][11] = int_stack + 610;
 Libderiv->deriv_classes[4][5][11] = int_stack + 835;
 Libderiv->deriv_classes[5][3][11] = int_stack + 1150;
 Libderiv->deriv_classes[5][4][11] = int_stack + 1360;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1675;
 Libderiv->deriv_classes[3][3][10] = int_stack + 2116;
 Libderiv->deriv_classes[3][4][10] = int_stack + 2216;
 Libderiv->deriv_classes[3][5][10] = int_stack + 2366;
 Libderiv->deriv_classes[4][3][10] = int_stack + 2576;
 Libderiv->deriv_classes[4][4][10] = int_stack + 2726;
 Libderiv->deriv_classes[4][5][10] = int_stack + 2951;
 Libderiv->deriv_classes[5][3][10] = int_stack + 3266;
 Libderiv->deriv_classes[5][4][10] = int_stack + 3476;
 Libderiv->deriv_classes[5][5][10] = int_stack + 3791;
 Libderiv->deriv_classes[3][3][9] = int_stack + 4232;
 Libderiv->deriv_classes[3][4][9] = int_stack + 4332;
 Libderiv->deriv_classes[3][5][9] = int_stack + 4482;
 Libderiv->deriv_classes[4][3][9] = int_stack + 4692;
 Libderiv->deriv_classes[4][4][9] = int_stack + 4842;
 Libderiv->deriv_classes[4][5][9] = int_stack + 5067;
 Libderiv->deriv_classes[5][3][9] = int_stack + 5382;
 Libderiv->deriv_classes[5][4][9] = int_stack + 5592;
 Libderiv->deriv_classes[5][5][9] = int_stack + 5907;
 Libderiv->deriv_classes[3][3][8] = int_stack + 6348;
 Libderiv->deriv_classes[3][4][8] = int_stack + 6448;
 Libderiv->deriv_classes[3][5][8] = int_stack + 6598;
 Libderiv->deriv_classes[4][3][8] = int_stack + 6808;
 Libderiv->deriv_classes[4][4][8] = int_stack + 6958;
 Libderiv->deriv_classes[4][5][8] = int_stack + 7183;
 Libderiv->deriv_classes[5][3][8] = int_stack + 7498;
 Libderiv->deriv_classes[5][4][8] = int_stack + 7708;
 Libderiv->deriv_classes[5][5][8] = int_stack + 8023;
 Libderiv->deriv_classes[3][3][7] = int_stack + 8464;
 Libderiv->deriv_classes[3][4][7] = int_stack + 8564;
 Libderiv->deriv_classes[3][5][7] = int_stack + 8714;
 Libderiv->deriv_classes[4][3][7] = int_stack + 8924;
 Libderiv->deriv_classes[4][4][7] = int_stack + 9074;
 Libderiv->deriv_classes[4][5][7] = int_stack + 9299;
 Libderiv->deriv_classes[5][3][7] = int_stack + 9614;
 Libderiv->deriv_classes[5][4][7] = int_stack + 9824;
 Libderiv->deriv_classes[5][5][7] = int_stack + 10139;
 Libderiv->deriv_classes[3][3][6] = int_stack + 10580;
 Libderiv->deriv_classes[3][4][6] = int_stack + 10680;
 Libderiv->deriv_classes[3][5][6] = int_stack + 10830;
 Libderiv->deriv_classes[4][3][6] = int_stack + 11040;
 Libderiv->deriv_classes[4][4][6] = int_stack + 11190;
 Libderiv->deriv_classes[4][5][6] = int_stack + 11415;
 Libderiv->dvrr_classes[5][3] = int_stack + 11730;
 Libderiv->deriv_classes[5][3][6] = int_stack + 11940;
 Libderiv->dvrr_classes[5][4] = int_stack + 12150;
 Libderiv->deriv_classes[5][4][6] = int_stack + 12465;
 Libderiv->deriv_classes[5][5][6] = int_stack + 12780;
 Libderiv->deriv_classes[3][3][2] = int_stack + 13221;
 Libderiv->deriv_classes[3][4][2] = int_stack + 13321;
 Libderiv->deriv_classes[3][5][2] = int_stack + 13471;
 Libderiv->deriv_classes[4][3][2] = int_stack + 13681;
 Libderiv->deriv_classes[4][4][2] = int_stack + 13831;
 Libderiv->deriv_classes[4][5][2] = int_stack + 14056;
 Libderiv->deriv_classes[5][3][2] = int_stack + 14371;
 Libderiv->deriv_classes[5][4][2] = int_stack + 14581;
 Libderiv->deriv_classes[5][5][2] = int_stack + 14896;
 Libderiv->deriv_classes[3][3][1] = int_stack + 15337;
 Libderiv->deriv_classes[3][4][1] = int_stack + 15437;
 Libderiv->deriv_classes[3][5][1] = int_stack + 15587;
 Libderiv->deriv_classes[4][3][1] = int_stack + 15797;
 Libderiv->deriv_classes[4][4][1] = int_stack + 15947;
 Libderiv->deriv_classes[4][5][1] = int_stack + 16172;
 Libderiv->deriv_classes[5][3][1] = int_stack + 16487;
 Libderiv->deriv_classes[5][4][1] = int_stack + 16697;
 Libderiv->deriv_classes[5][5][1] = int_stack + 17012;
 Libderiv->dvrr_classes[3][3] = int_stack + 17453;
 Libderiv->dvrr_classes[3][4] = int_stack + 17553;
 Libderiv->dvrr_classes[3][5] = int_stack + 17703;
 Libderiv->deriv_classes[3][3][0] = int_stack + 17913;
 Libderiv->deriv_classes[3][4][0] = int_stack + 18013;
 Libderiv->deriv_classes[3][5][0] = int_stack + 18163;
 Libderiv->dvrr_classes[4][3] = int_stack + 18373;
 Libderiv->dvrr_classes[4][4] = int_stack + 18523;
 Libderiv->dvrr_classes[4][5] = int_stack + 18748;
 Libderiv->deriv_classes[4][3][0] = int_stack + 19063;
 Libderiv->deriv_classes[4][4][0] = int_stack + 19213;
 Libderiv->deriv_classes[4][5][0] = int_stack + 19438;
 Libderiv->deriv_classes[5][3][0] = int_stack + 19753;
 Libderiv->deriv_classes[5][4][0] = int_stack + 19963;
 Libderiv->deriv_classes[5][5][0] = int_stack + 20278;
 memset(int_stack,0,165752);

 Libderiv->dvrr_stack = int_stack + 53104;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fdfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+20719,int_stack+17553,int_stack+17453,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21019,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17453,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21319,int_stack+250,int_stack+100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17553,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21769,int_stack+21319,int_stack+21019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20719,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+21019,int_stack+18523,int_stack+18373,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+610,int_stack+460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18373,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22369,int_stack+835,int_stack+610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18523,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23044,int_stack+22369,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21019,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+23944,int_stack+23044,int_stack+21769,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+12150,int_stack+11730,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21469,int_stack+1360,int_stack+1150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11730,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22099,int_stack+1675,int_stack+1360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12150,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+630,int_stack+22099,int_stack+21469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+25744,int_stack+630,int_stack+23044,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+2216,int_stack+2116, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17453, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+930,int_stack+2366,int_stack+2216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17553, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1380,int_stack+930,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20719, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+2726,int_stack+2576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18373, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1980,int_stack+2951,int_stack+2726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18523, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21469,int_stack+1980,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21019, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28444,int_stack+21469,int_stack+1380,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+3476,int_stack+3266, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11730, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1260,int_stack+3791,int_stack+3476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12150, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2205,int_stack+1260,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+30244,int_stack+2205,int_stack+21469,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21469,int_stack+4332,int_stack+4232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17453, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21769,int_stack+4482,int_stack+4332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17553, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22219,int_stack+21769,int_stack+21469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20719, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21469,int_stack+4842,int_stack+4692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18373, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22819,int_stack+5067,int_stack+4842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18523, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+630,int_stack+22819,int_stack+21469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21019, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1530,int_stack+630,int_stack+22219,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21469,int_stack+5592,int_stack+5382, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11730, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22099,int_stack+5907,int_stack+5592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12150, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3330,int_stack+22099,int_stack+21469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+32944,int_stack+3330,int_stack+630,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+6448,int_stack+6348, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17453, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+930,int_stack+6598,int_stack+6448, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3330,int_stack+930,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+6958,int_stack+6808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18373, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3930,int_stack+7183,int_stack+6958, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4605,int_stack+3930,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5505,int_stack+4605,int_stack+3330,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3330,int_stack+7708,int_stack+7498, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21469,int_stack+8023,int_stack+7708, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22414,int_stack+21469,int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+35644,int_stack+22414,int_stack+4605,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3330,int_stack+8564,int_stack+8464, 0.0, zero_stack, 1.0, int_stack+17453, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3630,int_stack+8714,int_stack+8564, 0.0, zero_stack, 1.0, int_stack+17553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4080,int_stack+3630,int_stack+3330, 0.0, zero_stack, 1.0, int_stack+20719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3330,int_stack+9074,int_stack+8924, 0.0, zero_stack, 1.0, int_stack+18373, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4680,int_stack+9299,int_stack+9074, 0.0, zero_stack, 1.0, int_stack+18523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+630,int_stack+4680,int_stack+3330, 0.0, zero_stack, 1.0, int_stack+21019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+21469,int_stack+630,int_stack+4080,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3330,int_stack+9824,int_stack+9614, 0.0, zero_stack, 1.0, int_stack+11730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3960,int_stack+10139,int_stack+9824, 0.0, zero_stack, 1.0, int_stack+12150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7305,int_stack+3960,int_stack+3330, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+38344,int_stack+7305,int_stack+630,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+10680,int_stack+10580, 1.0, int_stack+17453, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+930,int_stack+10830,int_stack+10680, 1.0, int_stack+17553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7305,int_stack+930,int_stack+630, 1.0, int_stack+20719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+11190,int_stack+11040, 1.0, int_stack+18373, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23269,int_stack+11415,int_stack+11190, 1.0, int_stack+18523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7905,int_stack+23269,int_stack+630, 1.0, int_stack+21019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8805,int_stack+7905,int_stack+7305,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+12465,int_stack+11940, 1.0, int_stack+11730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10605,int_stack+12780,int_stack+12465, 1.0, int_stack+12150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11550,int_stack+10605,int_stack+630, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+41044,int_stack+11550,int_stack+7905,60);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+17703,int_stack+17553,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+450,int_stack+0,int_stack+20719,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23269,int_stack+18748,int_stack+18523,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10605,int_stack+23269,int_stack+21019,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3330,int_stack+10605,int_stack+450,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+23269,int_stack+13321,int_stack+13221,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+13471,int_stack+13321,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+20719,int_stack+0,int_stack+23269,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+13831,int_stack+13681,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23269,int_stack+14056,int_stack+13831,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11505,int_stack+23269,int_stack+0,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+12405,int_stack+11505,int_stack+20719, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+20719,int_stack+14581,int_stack+14371,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7305,int_stack+14896,int_stack+14581,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+43744,int_stack+7305,int_stack+20719,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+45004,int_stack+43744,int_stack+11505, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11505,int_stack+15437,int_stack+15337,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+15587,int_stack+15437,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11805,int_stack+0,int_stack+11505,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+15947,int_stack+15797,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23269,int_stack+16172,int_stack+15947,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+43744,int_stack+23269,int_stack+0,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+14205,int_stack+43744,int_stack+11805, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+23269,int_stack+16697,int_stack+16487,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7305,int_stack+17012,int_stack+16697,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+16005,int_stack+7305,int_stack+23269,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+47704,int_stack+16005,int_stack+43744, 0.0, zero_stack, 1.0, int_stack+10605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+43744,int_stack+18013,int_stack+17913,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+18163,int_stack+18013,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+44044,int_stack+0,int_stack+43744,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+19213,int_stack+19063,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23269,int_stack+19438,int_stack+19213,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11505,int_stack+23269,int_stack+0,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+16005,int_stack+11505,int_stack+44044, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+19963,int_stack+19753,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+43744,int_stack+20278,int_stack+19963,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17805,int_stack+43744,int_stack+0,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+50404,int_stack+17805,int_stack+11505, 1.0, int_stack+10605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+17805,int_stack+25744,int_stack+23944,60);
     Libderiv->ABCD[11] = int_stack + 17805;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+23269,int_stack+30244,int_stack+28444,60);
     Libderiv->ABCD[10] = int_stack + 23269;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+26869,int_stack+32944,int_stack+1530,60);
     Libderiv->ABCD[9] = int_stack + 26869;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+30469,int_stack+35644,int_stack+5505,60);
     Libderiv->ABCD[8] = int_stack + 30469;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+34069,int_stack+38344,int_stack+21469,60);
     Libderiv->ABCD[7] = int_stack + 34069;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+5130,int_stack+41044,int_stack+8805,60);
     Libderiv->ABCD[6] = int_stack + 5130;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+37669,int_stack+45004,int_stack+12405, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 37669;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+41269,int_stack+47704,int_stack+14205, 0.0, zero_stack, 1.0, int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 41269;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+44869,int_stack+50404,int_stack+16005, 1.0, int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 44869;

}
