#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpgg(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|gg) integrals */

void d1hrr_order_dpgg(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[2][6][11] = int_stack + 216;
 Libderiv->deriv_classes[2][7][11] = int_stack + 384;
 Libderiv->deriv_classes[2][8][11] = int_stack + 600;
 Libderiv->deriv_classes[3][4][11] = int_stack + 870;
 Libderiv->deriv_classes[3][5][11] = int_stack + 1020;
 Libderiv->deriv_classes[3][6][11] = int_stack + 1230;
 Libderiv->deriv_classes[3][7][11] = int_stack + 1510;
 Libderiv->deriv_classes[3][8][11] = int_stack + 1870;
 Libderiv->deriv_classes[2][4][10] = int_stack + 2320;
 Libderiv->deriv_classes[2][5][10] = int_stack + 2410;
 Libderiv->deriv_classes[2][6][10] = int_stack + 2536;
 Libderiv->deriv_classes[2][7][10] = int_stack + 2704;
 Libderiv->deriv_classes[2][8][10] = int_stack + 2920;
 Libderiv->deriv_classes[3][4][10] = int_stack + 3190;
 Libderiv->deriv_classes[3][5][10] = int_stack + 3340;
 Libderiv->deriv_classes[3][6][10] = int_stack + 3550;
 Libderiv->deriv_classes[3][7][10] = int_stack + 3830;
 Libderiv->deriv_classes[3][8][10] = int_stack + 4190;
 Libderiv->deriv_classes[2][4][9] = int_stack + 4640;
 Libderiv->deriv_classes[2][5][9] = int_stack + 4730;
 Libderiv->deriv_classes[2][6][9] = int_stack + 4856;
 Libderiv->deriv_classes[2][7][9] = int_stack + 5024;
 Libderiv->deriv_classes[2][8][9] = int_stack + 5240;
 Libderiv->deriv_classes[3][4][9] = int_stack + 5510;
 Libderiv->deriv_classes[3][5][9] = int_stack + 5660;
 Libderiv->deriv_classes[3][6][9] = int_stack + 5870;
 Libderiv->deriv_classes[3][7][9] = int_stack + 6150;
 Libderiv->deriv_classes[3][8][9] = int_stack + 6510;
 Libderiv->deriv_classes[2][4][8] = int_stack + 6960;
 Libderiv->deriv_classes[2][5][8] = int_stack + 7050;
 Libderiv->deriv_classes[2][6][8] = int_stack + 7176;
 Libderiv->deriv_classes[2][7][8] = int_stack + 7344;
 Libderiv->deriv_classes[2][8][8] = int_stack + 7560;
 Libderiv->deriv_classes[3][4][8] = int_stack + 7830;
 Libderiv->deriv_classes[3][5][8] = int_stack + 7980;
 Libderiv->deriv_classes[3][6][8] = int_stack + 8190;
 Libderiv->deriv_classes[3][7][8] = int_stack + 8470;
 Libderiv->deriv_classes[3][8][8] = int_stack + 8830;
 Libderiv->deriv_classes[2][4][7] = int_stack + 9280;
 Libderiv->deriv_classes[2][5][7] = int_stack + 9370;
 Libderiv->deriv_classes[2][6][7] = int_stack + 9496;
 Libderiv->deriv_classes[2][7][7] = int_stack + 9664;
 Libderiv->deriv_classes[2][8][7] = int_stack + 9880;
 Libderiv->deriv_classes[3][4][7] = int_stack + 10150;
 Libderiv->deriv_classes[3][5][7] = int_stack + 10300;
 Libderiv->deriv_classes[3][6][7] = int_stack + 10510;
 Libderiv->deriv_classes[3][7][7] = int_stack + 10790;
 Libderiv->deriv_classes[3][8][7] = int_stack + 11150;
 Libderiv->deriv_classes[2][4][6] = int_stack + 11600;
 Libderiv->deriv_classes[2][5][6] = int_stack + 11690;
 Libderiv->deriv_classes[2][6][6] = int_stack + 11816;
 Libderiv->deriv_classes[2][7][6] = int_stack + 11984;
 Libderiv->deriv_classes[2][8][6] = int_stack + 12200;
 Libderiv->dvrr_classes[3][4] = int_stack + 12470;
 Libderiv->deriv_classes[3][4][6] = int_stack + 12620;
 Libderiv->dvrr_classes[3][5] = int_stack + 12770;
 Libderiv->deriv_classes[3][5][6] = int_stack + 12980;
 Libderiv->dvrr_classes[3][6] = int_stack + 13190;
 Libderiv->deriv_classes[3][6][6] = int_stack + 13470;
 Libderiv->dvrr_classes[3][7] = int_stack + 13750;
 Libderiv->deriv_classes[3][7][6] = int_stack + 14110;
 Libderiv->deriv_classes[3][8][6] = int_stack + 14470;
 Libderiv->deriv_classes[2][4][2] = int_stack + 14920;
 Libderiv->deriv_classes[2][5][2] = int_stack + 15010;
 Libderiv->deriv_classes[2][6][2] = int_stack + 15136;
 Libderiv->deriv_classes[2][7][2] = int_stack + 15304;
 Libderiv->deriv_classes[2][8][2] = int_stack + 15520;
 Libderiv->deriv_classes[3][4][2] = int_stack + 15790;
 Libderiv->deriv_classes[3][5][2] = int_stack + 15940;
 Libderiv->deriv_classes[3][6][2] = int_stack + 16150;
 Libderiv->deriv_classes[3][7][2] = int_stack + 16430;
 Libderiv->deriv_classes[3][8][2] = int_stack + 16790;
 Libderiv->deriv_classes[2][4][1] = int_stack + 17240;
 Libderiv->deriv_classes[2][5][1] = int_stack + 17330;
 Libderiv->deriv_classes[2][6][1] = int_stack + 17456;
 Libderiv->deriv_classes[2][7][1] = int_stack + 17624;
 Libderiv->deriv_classes[2][8][1] = int_stack + 17840;
 Libderiv->deriv_classes[3][4][1] = int_stack + 18110;
 Libderiv->deriv_classes[3][5][1] = int_stack + 18260;
 Libderiv->deriv_classes[3][6][1] = int_stack + 18470;
 Libderiv->deriv_classes[3][7][1] = int_stack + 18750;
 Libderiv->deriv_classes[3][8][1] = int_stack + 19110;
 Libderiv->dvrr_classes[2][4] = int_stack + 19560;
 Libderiv->dvrr_classes[2][5] = int_stack + 19650;
 Libderiv->dvrr_classes[2][6] = int_stack + 19776;
 Libderiv->dvrr_classes[2][7] = int_stack + 19944;
 Libderiv->dvrr_classes[2][8] = int_stack + 20160;
 Libderiv->deriv_classes[2][4][0] = int_stack + 20430;
 Libderiv->deriv_classes[2][5][0] = int_stack + 20520;
 Libderiv->deriv_classes[2][6][0] = int_stack + 20646;
 Libderiv->deriv_classes[2][7][0] = int_stack + 20814;
 Libderiv->deriv_classes[2][8][0] = int_stack + 21030;
 Libderiv->deriv_classes[3][4][0] = int_stack + 21300;
 Libderiv->deriv_classes[3][5][0] = int_stack + 21450;
 Libderiv->deriv_classes[3][6][0] = int_stack + 21660;
 Libderiv->deriv_classes[3][7][0] = int_stack + 21940;
 Libderiv->deriv_classes[3][8][0] = int_stack + 22300;
 memset(int_stack,0,182000);

 Libderiv->dvrr_stack = int_stack + 61666;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpgg(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22750,int_stack+19650,int_stack+19560,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+23020,int_stack+19776,int_stack+19650,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23398,int_stack+23020,int_stack+22750,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+23938,int_stack+19944,int_stack+19776,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+24442,int_stack+23938,int_stack+23020,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+25198,int_stack+24442,int_stack+23398,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26098,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19560,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26368,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19650,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+26746,int_stack+26368,int_stack+26098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22750,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+27286,int_stack+384,int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19776,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+27790,int_stack+27286,int_stack+26368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23020,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+28546,int_stack+27790,int_stack+26746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23398,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+26098,int_stack+600,int_stack+384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19944,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+29446,int_stack+26098,int_stack+27286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23938,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+26098,int_stack+29446,int_stack+27790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24442,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+29446,int_stack+26098,int_stack+28546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25198,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+26098,int_stack+12770,int_stack+12470,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+26548,int_stack+13190,int_stack+12770,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+27178,int_stack+26548,int_stack+26098,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+28078,int_stack+13750,int_stack+13190,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+30796,int_stack+28078,int_stack+26548,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+32056,int_stack+30796,int_stack+27178,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+28918,int_stack+1020,int_stack+870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12470,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+1230,int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12770,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33556,int_stack+0,int_stack+28918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26098,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+34456,int_stack+1510,int_stack+1230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13190,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+35296,int_stack+34456,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26548,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+35296,int_stack+33556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27178,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+36556,int_stack+1870,int_stack+1510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13750,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+37636,int_stack+36556,int_stack+34456, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28078,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+39316,int_stack+37636,int_stack+35296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30796,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+33556,int_stack+39316,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32056,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+2410,int_stack+2320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19560, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+270,int_stack+2536,int_stack+2410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19650, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+648,int_stack+270,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22750, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+1188,int_stack+2704,int_stack+2536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19776, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1692,int_stack+1188,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23020, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+35806,int_stack+1692,int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23398, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+2920,int_stack+2704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19944, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+36706,int_stack+0,int_stack+1188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23938, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+0,int_stack+36706,int_stack+1692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24442, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+36706,int_stack+0,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25198, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+35806,int_stack+3340,int_stack+3190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12470, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+3550,int_stack+3340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12770, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+630,int_stack+0,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26098, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+35806,int_stack+3830,int_stack+3550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13190, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+1530,int_stack+35806,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26548, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+38056,int_stack+1530,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27178, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+0,int_stack+4190,int_stack+3830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13750, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+2790,int_stack+0,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28078, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+39556,int_stack+2790,int_stack+1530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30796, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+0,int_stack+39556,int_stack+38056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32056, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+38056,int_stack+4730,int_stack+4640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19560, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38326,int_stack+4856,int_stack+4730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19650, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+38704,int_stack+38326,int_stack+38056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22750, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+39244,int_stack+5024,int_stack+4856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19776, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+39748,int_stack+39244,int_stack+38326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23020, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+35806,int_stack+39748,int_stack+38704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23398, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+38056,int_stack+5240,int_stack+5024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19944, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+40504,int_stack+38056,int_stack+39244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23938, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+38056,int_stack+40504,int_stack+39748, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24442, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+39316,int_stack+38056,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25198, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+35806,int_stack+5660,int_stack+5510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12470, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38056,int_stack+5870,int_stack+5660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12770, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+40666,int_stack+38056,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26098, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+35806,int_stack+6150,int_stack+5870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13190, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2250,int_stack+35806,int_stack+38056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26548, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3510,int_stack+2250,int_stack+40666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27178, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+40666,int_stack+6510,int_stack+6150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13750, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+5010,int_stack+40666,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28078, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+40666,int_stack+5010,int_stack+2250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30796, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+42766,int_stack+40666,int_stack+3510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32056, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+40666,int_stack+7050,int_stack+6960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+40936,int_stack+7176,int_stack+7050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41314,int_stack+40936,int_stack+40666, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+41854,int_stack+7344,int_stack+7176, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2250,int_stack+41854,int_stack+40936, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+35806,int_stack+2250,int_stack+41314, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+40666,int_stack+7560,int_stack+7344, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3006,int_stack+40666,int_stack+41854, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+38056,int_stack+3006,int_stack+2250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+2250,int_stack+38056,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+35806,int_stack+7980,int_stack+7830, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38056,int_stack+8190,int_stack+7980, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3600,int_stack+38056,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+35806,int_stack+8470,int_stack+8190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+4500,int_stack+35806,int_stack+38056, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+5760,int_stack+4500,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+38056,int_stack+8830,int_stack+8470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+7260,int_stack+38056,int_stack+35806, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+40666,int_stack+7260,int_stack+4500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+45016,int_stack+40666,int_stack+5760, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+40666,int_stack+9370,int_stack+9280, 0.0, zero_stack, 1.0, int_stack+19560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+40936,int_stack+9496,int_stack+9370, 0.0, zero_stack, 1.0, int_stack+19650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41314,int_stack+40936,int_stack+40666, 0.0, zero_stack, 1.0, int_stack+22750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+41854,int_stack+9664,int_stack+9496, 0.0, zero_stack, 1.0, int_stack+19776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+35806,int_stack+41854,int_stack+40936, 0.0, zero_stack, 1.0, int_stack+23020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+38056,int_stack+35806,int_stack+41314, 0.0, zero_stack, 1.0, int_stack+23398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+40666,int_stack+9880,int_stack+9664, 0.0, zero_stack, 1.0, int_stack+19944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+3600,int_stack+40666,int_stack+41854, 0.0, zero_stack, 1.0, int_stack+23938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+40666,int_stack+3600,int_stack+35806, 0.0, zero_stack, 1.0, int_stack+24442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+3600,int_stack+40666,int_stack+38056, 0.0, zero_stack, 1.0, int_stack+25198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+38056,int_stack+10300,int_stack+10150, 0.0, zero_stack, 1.0, int_stack+12470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38506,int_stack+10510,int_stack+10300, 0.0, zero_stack, 1.0, int_stack+12770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+35806,int_stack+38506,int_stack+38056, 0.0, zero_stack, 1.0, int_stack+26098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+40666,int_stack+10790,int_stack+10510, 0.0, zero_stack, 1.0, int_stack+13190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+41506,int_stack+40666,int_stack+38506, 0.0, zero_stack, 1.0, int_stack+26548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+4950,int_stack+41506,int_stack+35806, 0.0, zero_stack, 1.0, int_stack+27178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+38056,int_stack+11150,int_stack+10790, 0.0, zero_stack, 1.0, int_stack+13750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+6450,int_stack+38056,int_stack+40666, 0.0, zero_stack, 1.0, int_stack+28078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+8130,int_stack+6450,int_stack+41506, 0.0, zero_stack, 1.0, int_stack+30796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+47266,int_stack+8130,int_stack+4950, 0.0, zero_stack, 1.0, int_stack+32056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4950,int_stack+11690,int_stack+11600, 1.0, int_stack+19560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5220,int_stack+11816,int_stack+11690, 1.0, int_stack+19650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5598,int_stack+5220,int_stack+4950, 1.0, int_stack+22750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+6138,int_stack+11984,int_stack+11816, 1.0, int_stack+19776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+6642,int_stack+6138,int_stack+5220, 1.0, int_stack+23020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+35806,int_stack+6642,int_stack+5598, 1.0, int_stack+23398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+22750,int_stack+12200,int_stack+11984, 1.0, int_stack+19944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+4950,int_stack+22750,int_stack+6138, 1.0, int_stack+23938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+38056,int_stack+4950,int_stack+6642, 1.0, int_stack+24442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+4950,int_stack+38056,int_stack+35806, 1.0, int_stack+25198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+35806,int_stack+12980,int_stack+12620, 1.0, int_stack+12470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38056,int_stack+13470,int_stack+12980, 1.0, int_stack+12770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6300,int_stack+38056,int_stack+35806, 1.0, int_stack+26098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+35806,int_stack+14110,int_stack+13470, 1.0, int_stack+13190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7200,int_stack+35806,int_stack+38056, 1.0, int_stack+26548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8460,int_stack+7200,int_stack+6300, 1.0, int_stack+27178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|kp) ---*/
   d1hrr3_build_kp(Libderiv->CD,int_stack+38056,int_stack+14470,int_stack+14110, 1.0, int_stack+13750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|id) ---*/
   d1hrr3_build_id(Libderiv->CD,int_stack+26098,int_stack+38056,int_stack+35806, 1.0, int_stack+28078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hf) ---*/
   d1hrr3_build_hf(Libderiv->CD,int_stack+40666,int_stack+26098,int_stack+7200, 1.0, int_stack+30796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gg) ---*/
   d1hrr3_build_gg(Libderiv->CD,int_stack+26098,int_stack+40666,int_stack+8460, 1.0, int_stack+32056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+40666,int_stack+20160,int_stack+19944,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+41314,int_stack+40666,int_stack+23938,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+38056,int_stack+41314,int_stack+24442,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+40666,int_stack+38056,int_stack+25198,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+38056,int_stack+15010,int_stack+14920,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+38326,int_stack+15136,int_stack+15010,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+38704,int_stack+38326,int_stack+38056,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+42016,int_stack+15304,int_stack+15136,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+30796,int_stack+42016,int_stack+38326,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+35806,int_stack+30796,int_stack+38704,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+38056,int_stack+15520,int_stack+15304,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+31552,int_stack+38056,int_stack+42016,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+38056,int_stack+31552,int_stack+30796,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+30796,int_stack+38056,int_stack+35806,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35806,int_stack+15940,int_stack+15790,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+38056,int_stack+16150,int_stack+15940,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+32146,int_stack+38056,int_stack+35806,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+35806,int_stack+16430,int_stack+16150,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+6300,int_stack+35806,int_stack+38056,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+7560,int_stack+6300,int_stack+32146,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+32146,int_stack+16790,int_stack+16430,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+9060,int_stack+32146,int_stack+35806,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+10740,int_stack+9060,int_stack+6300,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+12840,int_stack+10740,int_stack+7560,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6300,int_stack+17330,int_stack+17240,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6570,int_stack+17456,int_stack+17330,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+6948,int_stack+6570,int_stack+6300,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+7488,int_stack+17624,int_stack+17456,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+7992,int_stack+7488,int_stack+6570,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+35806,int_stack+7992,int_stack+6948,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+6300,int_stack+17840,int_stack+17624,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+8748,int_stack+6300,int_stack+7488,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+38056,int_stack+8748,int_stack+7992,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+6300,int_stack+38056,int_stack+35806,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35806,int_stack+18260,int_stack+18110,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+38056,int_stack+18470,int_stack+18260,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+7650,int_stack+38056,int_stack+35806,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+35806,int_stack+18750,int_stack+18470,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+8550,int_stack+35806,int_stack+38056,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+9810,int_stack+8550,int_stack+7650,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+38056,int_stack+19110,int_stack+18750,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+15090,int_stack+38056,int_stack+35806,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+16770,int_stack+15090,int_stack+8550,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+22750,int_stack+16770,int_stack+9810,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+15090,int_stack+20520,int_stack+20430,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15360,int_stack+20646,int_stack+20520,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15738,int_stack+15360,int_stack+15090,6);
 /*--- compute (d0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+16278,int_stack+20814,int_stack+20646,6);
 /*--- compute (d0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+16782,int_stack+16278,int_stack+15360,6);
 /*--- compute (d0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+35806,int_stack+16782,int_stack+15738,6);
 /*--- compute (d0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+15090,int_stack+21030,int_stack+20814,6);
 /*--- compute (d0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+17538,int_stack+15090,int_stack+16278,6);
 /*--- compute (d0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+38056,int_stack+17538,int_stack+16782,6);
 /*--- compute (d0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+15090,int_stack+38056,int_stack+35806,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35806,int_stack+21450,int_stack+21300,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+38056,int_stack+21660,int_stack+21450,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16440,int_stack+38056,int_stack+35806,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+35806,int_stack+21940,int_stack+21660,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+17340,int_stack+35806,int_stack+38056,10);
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+18600,int_stack+17340,int_stack+16440,10);
 /*--- compute (f0|kp) ---*/
   hrr3_build_kp(Libderiv->CD,int_stack+38056,int_stack+22300,int_stack+21940,10);
 /*--- compute (f0|id) ---*/
   hrr3_build_id(Libderiv->CD,int_stack+20100,int_stack+38056,int_stack+35806,10);
 /*--- compute (f0|hf) ---*/
   hrr3_build_hf(Libderiv->CD,int_stack+7650,int_stack+20100,int_stack+17340,10);
 /*--- compute (f0|gg) ---*/
   hrr3_build_gg(Libderiv->CD,int_stack+20100,int_stack+7650,int_stack+18600,10);
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7650,int_stack+33556,int_stack+29446,225);
     Libderiv->ABCD[11] = int_stack + 7650;
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+32146,int_stack+0,int_stack+36706,225);
     Libderiv->ABCD[10] = int_stack + 32146;
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+49516,int_stack+42766,int_stack+39316,225);
     Libderiv->ABCD[9] = int_stack + 49516;
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+36196,int_stack+45016,int_stack+2250,225);
     Libderiv->ABCD[8] = int_stack + 36196;
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+42016,int_stack+47266,int_stack+3600,225);
     Libderiv->ABCD[7] = int_stack + 42016;
 /*--- compute (dp|gg) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+26098,int_stack+4950,225);
     Libderiv->ABCD[6] = int_stack + 0;
 /*--- compute (dp|gg) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+25000,int_stack+12840,int_stack+30796, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[2] = int_stack + 25000;
 /*--- compute (dp|gg) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+53566,int_stack+22750,int_stack+6300, 0.0, zero_stack, 1.0, int_stack+40666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[1] = int_stack + 53566;
 /*--- compute (dp|gg) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+57616,int_stack+20100,int_stack+15090, 1.0, int_stack+40666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,225);
     Libderiv->ABCD[0] = int_stack + 57616;

}
