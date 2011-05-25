#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|ff) integrals */

void d1hrr_order_fpff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 100;
 Libderiv->deriv_classes[3][5][11] = int_stack + 250;
 Libderiv->deriv_classes[3][6][11] = int_stack + 460;
 Libderiv->deriv_classes[4][3][11] = int_stack + 740;
 Libderiv->deriv_classes[4][4][11] = int_stack + 890;
 Libderiv->deriv_classes[4][5][11] = int_stack + 1115;
 Libderiv->deriv_classes[4][6][11] = int_stack + 1430;
 Libderiv->deriv_classes[3][3][10] = int_stack + 1850;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1950;
 Libderiv->deriv_classes[3][5][10] = int_stack + 2100;
 Libderiv->deriv_classes[3][6][10] = int_stack + 2310;
 Libderiv->deriv_classes[4][3][10] = int_stack + 2590;
 Libderiv->deriv_classes[4][4][10] = int_stack + 2740;
 Libderiv->deriv_classes[4][5][10] = int_stack + 2965;
 Libderiv->deriv_classes[4][6][10] = int_stack + 3280;
 Libderiv->deriv_classes[3][3][9] = int_stack + 3700;
 Libderiv->deriv_classes[3][4][9] = int_stack + 3800;
 Libderiv->deriv_classes[3][5][9] = int_stack + 3950;
 Libderiv->deriv_classes[3][6][9] = int_stack + 4160;
 Libderiv->deriv_classes[4][3][9] = int_stack + 4440;
 Libderiv->deriv_classes[4][4][9] = int_stack + 4590;
 Libderiv->deriv_classes[4][5][9] = int_stack + 4815;
 Libderiv->deriv_classes[4][6][9] = int_stack + 5130;
 Libderiv->deriv_classes[3][3][8] = int_stack + 5550;
 Libderiv->deriv_classes[3][4][8] = int_stack + 5650;
 Libderiv->deriv_classes[3][5][8] = int_stack + 5800;
 Libderiv->deriv_classes[3][6][8] = int_stack + 6010;
 Libderiv->deriv_classes[4][3][8] = int_stack + 6290;
 Libderiv->deriv_classes[4][4][8] = int_stack + 6440;
 Libderiv->deriv_classes[4][5][8] = int_stack + 6665;
 Libderiv->deriv_classes[4][6][8] = int_stack + 6980;
 Libderiv->deriv_classes[3][3][7] = int_stack + 7400;
 Libderiv->deriv_classes[3][4][7] = int_stack + 7500;
 Libderiv->deriv_classes[3][5][7] = int_stack + 7650;
 Libderiv->deriv_classes[3][6][7] = int_stack + 7860;
 Libderiv->deriv_classes[4][3][7] = int_stack + 8140;
 Libderiv->deriv_classes[4][4][7] = int_stack + 8290;
 Libderiv->deriv_classes[4][5][7] = int_stack + 8515;
 Libderiv->deriv_classes[4][6][7] = int_stack + 8830;
 Libderiv->deriv_classes[3][3][6] = int_stack + 9250;
 Libderiv->deriv_classes[3][4][6] = int_stack + 9350;
 Libderiv->deriv_classes[3][5][6] = int_stack + 9500;
 Libderiv->deriv_classes[3][6][6] = int_stack + 9710;
 Libderiv->dvrr_classes[4][3] = int_stack + 9990;
 Libderiv->deriv_classes[4][3][6] = int_stack + 10140;
 Libderiv->dvrr_classes[4][4] = int_stack + 10290;
 Libderiv->deriv_classes[4][4][6] = int_stack + 10515;
 Libderiv->dvrr_classes[4][5] = int_stack + 10740;
 Libderiv->deriv_classes[4][5][6] = int_stack + 11055;
 Libderiv->deriv_classes[4][6][6] = int_stack + 11370;
 Libderiv->deriv_classes[3][3][2] = int_stack + 11790;
 Libderiv->deriv_classes[3][4][2] = int_stack + 11890;
 Libderiv->deriv_classes[3][5][2] = int_stack + 12040;
 Libderiv->deriv_classes[3][6][2] = int_stack + 12250;
 Libderiv->deriv_classes[4][3][2] = int_stack + 12530;
 Libderiv->deriv_classes[4][4][2] = int_stack + 12680;
 Libderiv->deriv_classes[4][5][2] = int_stack + 12905;
 Libderiv->deriv_classes[4][6][2] = int_stack + 13220;
 Libderiv->deriv_classes[3][3][1] = int_stack + 13640;
 Libderiv->deriv_classes[3][4][1] = int_stack + 13740;
 Libderiv->deriv_classes[3][5][1] = int_stack + 13890;
 Libderiv->deriv_classes[3][6][1] = int_stack + 14100;
 Libderiv->deriv_classes[4][3][1] = int_stack + 14380;
 Libderiv->deriv_classes[4][4][1] = int_stack + 14530;
 Libderiv->deriv_classes[4][5][1] = int_stack + 14755;
 Libderiv->deriv_classes[4][6][1] = int_stack + 15070;
 Libderiv->dvrr_classes[3][3] = int_stack + 15490;
 Libderiv->dvrr_classes[3][4] = int_stack + 15590;
 Libderiv->dvrr_classes[3][5] = int_stack + 15740;
 Libderiv->dvrr_classes[3][6] = int_stack + 15950;
 Libderiv->deriv_classes[3][3][0] = int_stack + 16230;
 Libderiv->deriv_classes[3][4][0] = int_stack + 16330;
 Libderiv->deriv_classes[3][5][0] = int_stack + 16480;
 Libderiv->deriv_classes[3][6][0] = int_stack + 16690;
 Libderiv->deriv_classes[4][3][0] = int_stack + 16970;
 Libderiv->deriv_classes[4][4][0] = int_stack + 17120;
 Libderiv->deriv_classes[4][5][0] = int_stack + 17345;
 Libderiv->deriv_classes[4][6][0] = int_stack + 17660;
 memset(int_stack,0,144640);

 Libderiv->dvrr_stack = int_stack + 41130;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpff(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+18080,int_stack+15590,int_stack+15490,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18380,int_stack+15740,int_stack+15590,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+18830,int_stack+18380,int_stack+18080,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+19430,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15490,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19730,int_stack+250,int_stack+100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15590,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20180,int_stack+19730,int_stack+19430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18080,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20780,int_stack+460,int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15740,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21410,int_stack+20780,int_stack+19730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18380,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+22310,int_stack+21410,int_stack+20180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18830,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19430,int_stack+10290,int_stack+9990,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+19880,int_stack+10740,int_stack+10290,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+20555,int_stack+19880,int_stack+19430,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21455,int_stack+890,int_stack+740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9990,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1115,int_stack+890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10290,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23310,int_stack+0,int_stack+21455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19430,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24210,int_stack+1430,int_stack+1115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10740,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+25155,int_stack+24210,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19880,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+25155,int_stack+23310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20555,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23310,int_stack+1950,int_stack+1850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15490, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1500,int_stack+2100,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15590, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23610,int_stack+1500,int_stack+23310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18080, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24210,int_stack+2310,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15740, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24840,int_stack+24210,int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18380, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1500,int_stack+24840,int_stack+23610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18830, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23310,int_stack+2740,int_stack+2590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9990, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23760,int_stack+2965,int_stack+2740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10290, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24435,int_stack+23760,int_stack+23310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19430, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+25335,int_stack+3280,int_stack+2965, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10740, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+26280,int_stack+25335,int_stack+23760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19880, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+27630,int_stack+26280,int_stack+24435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20555, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23310,int_stack+3800,int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15490, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23610,int_stack+3950,int_stack+3800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15590, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24060,int_stack+23610,int_stack+23310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18080, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24660,int_stack+4160,int_stack+3950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15740, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+25290,int_stack+24660,int_stack+23610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18380, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+26190,int_stack+25290,int_stack+24060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18830, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23310,int_stack+4590,int_stack+4440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9990, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23760,int_stack+4815,int_stack+4590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10290, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24435,int_stack+23760,int_stack+23310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19430, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2500,int_stack+5130,int_stack+4815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10740, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3445,int_stack+2500,int_stack+23760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19880, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+29130,int_stack+3445,int_stack+24435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20555, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2500,int_stack+5650,int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2800,int_stack+5800,int_stack+5650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3250,int_stack+2800,int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3850,int_stack+6010,int_stack+5800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4480,int_stack+3850,int_stack+2800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+23310,int_stack+4480,int_stack+3250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2500,int_stack+6440,int_stack+6290, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2950,int_stack+6665,int_stack+6440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3625,int_stack+2950,int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4525,int_stack+6980,int_stack+6665, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5470,int_stack+4525,int_stack+2950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+24310,int_stack+5470,int_stack+3625, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2500,int_stack+7500,int_stack+7400, 0.0, zero_stack, 1.0, int_stack+15490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2800,int_stack+7650,int_stack+7500, 0.0, zero_stack, 1.0, int_stack+15590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3250,int_stack+2800,int_stack+2500, 0.0, zero_stack, 1.0, int_stack+18080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3850,int_stack+7860,int_stack+7650, 0.0, zero_stack, 1.0, int_stack+15740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4480,int_stack+3850,int_stack+2800, 0.0, zero_stack, 1.0, int_stack+18380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5380,int_stack+4480,int_stack+3250, 0.0, zero_stack, 1.0, int_stack+18830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2500,int_stack+8290,int_stack+8140, 0.0, zero_stack, 1.0, int_stack+9990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2950,int_stack+8515,int_stack+8290, 0.0, zero_stack, 1.0, int_stack+10290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3625,int_stack+2950,int_stack+2500, 0.0, zero_stack, 1.0, int_stack+19430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6380,int_stack+8830,int_stack+8515, 0.0, zero_stack, 1.0, int_stack+10740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7325,int_stack+6380,int_stack+2950, 0.0, zero_stack, 1.0, int_stack+19880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+30630,int_stack+7325,int_stack+3625, 0.0, zero_stack, 1.0, int_stack+20555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6380,int_stack+9350,int_stack+9250, 1.0, int_stack+15490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6680,int_stack+9500,int_stack+9350, 1.0, int_stack+15590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7130,int_stack+6680,int_stack+6380, 1.0, int_stack+18080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7730,int_stack+9710,int_stack+9500, 1.0, int_stack+15740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8360,int_stack+7730,int_stack+6680, 1.0, int_stack+18380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2500,int_stack+8360,int_stack+7130, 1.0, int_stack+18830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6380,int_stack+10515,int_stack+10140, 1.0, int_stack+9990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6830,int_stack+11055,int_stack+10515, 1.0, int_stack+10290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7505,int_stack+6830,int_stack+6380, 1.0, int_stack+19430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8405,int_stack+11370,int_stack+11055, 1.0, int_stack+10740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9350,int_stack+8405,int_stack+6830, 1.0, int_stack+19880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3500,int_stack+9350,int_stack+7505, 1.0, int_stack+20555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+19430,int_stack+15950,int_stack+15740,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+20060,int_stack+19430,int_stack+18380,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+20960,int_stack+20060,int_stack+18830,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+21960,int_stack+11890,int_stack+11790,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6380,int_stack+12040,int_stack+11890,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+6830,int_stack+6380,int_stack+21960,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7430,int_stack+12250,int_stack+12040,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8060,int_stack+7430,int_stack+6380,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+8960,int_stack+8060,int_stack+6830,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6380,int_stack+12680,int_stack+12530,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6830,int_stack+12905,int_stack+12680,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+7505,int_stack+6830,int_stack+6380,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+9960,int_stack+13220,int_stack+12905,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10905,int_stack+9960,int_stack+6830,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+18080,int_stack+10905,int_stack+7505,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+9960,int_stack+13740,int_stack+13640,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10260,int_stack+13890,int_stack+13740,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+10710,int_stack+10260,int_stack+9960,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11310,int_stack+14100,int_stack+13890,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11940,int_stack+11310,int_stack+10260,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+12840,int_stack+11940,int_stack+10710,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+9960,int_stack+14530,int_stack+14380,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10410,int_stack+14755,int_stack+14530,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11085,int_stack+10410,int_stack+9960,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6380,int_stack+15070,int_stack+14755,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+13840,int_stack+6380,int_stack+10410,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+6380,int_stack+13840,int_stack+11085,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13840,int_stack+16330,int_stack+16230,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14140,int_stack+16480,int_stack+16330,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14590,int_stack+14140,int_stack+13840,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15190,int_stack+16690,int_stack+16480,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+15820,int_stack+15190,int_stack+14140,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+7880,int_stack+15820,int_stack+14590,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13840,int_stack+17120,int_stack+16970,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14290,int_stack+17345,int_stack+17120,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+14965,int_stack+14290,int_stack+13840,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15865,int_stack+17660,int_stack+17345,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9960,int_stack+15865,int_stack+14290,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+15865,int_stack+9960,int_stack+14965,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+32130,int_stack+0,int_stack+22310,100);
     Libderiv->ABCD[11] = int_stack + 32130;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+35130,int_stack+27630,int_stack+1500,100);
     Libderiv->ABCD[10] = int_stack + 35130;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+38130,int_stack+29130,int_stack+26190,100);
     Libderiv->ABCD[9] = int_stack + 38130;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+25810,int_stack+24310,int_stack+23310,100);
     Libderiv->ABCD[8] = int_stack + 25810;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+21960,int_stack+30630,int_stack+5380,100);
     Libderiv->ABCD[7] = int_stack + 21960;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28810,int_stack+3500,int_stack+2500,100);
     Libderiv->ABCD[6] = int_stack + 28810;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+18080,int_stack+8960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 0;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+3000,int_stack+6380,int_stack+12840, 0.0, zero_stack, 1.0, int_stack+20960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 3000;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+17365,int_stack+15865,int_stack+7880, 1.0, int_stack+20960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 17365;

}
