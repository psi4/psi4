#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fdff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fd|ff) integrals */

void d1hrr_order_fdff(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[5][3][11] = int_stack + 1850;
 Libderiv->deriv_classes[5][4][11] = int_stack + 2060;
 Libderiv->deriv_classes[5][5][11] = int_stack + 2375;
 Libderiv->deriv_classes[5][6][11] = int_stack + 2816;
 Libderiv->deriv_classes[3][3][10] = int_stack + 3404;
 Libderiv->deriv_classes[3][4][10] = int_stack + 3504;
 Libderiv->deriv_classes[3][5][10] = int_stack + 3654;
 Libderiv->deriv_classes[3][6][10] = int_stack + 3864;
 Libderiv->deriv_classes[4][3][10] = int_stack + 4144;
 Libderiv->deriv_classes[4][4][10] = int_stack + 4294;
 Libderiv->deriv_classes[4][5][10] = int_stack + 4519;
 Libderiv->deriv_classes[4][6][10] = int_stack + 4834;
 Libderiv->deriv_classes[5][3][10] = int_stack + 5254;
 Libderiv->deriv_classes[5][4][10] = int_stack + 5464;
 Libderiv->deriv_classes[5][5][10] = int_stack + 5779;
 Libderiv->deriv_classes[5][6][10] = int_stack + 6220;
 Libderiv->deriv_classes[3][3][9] = int_stack + 6808;
 Libderiv->deriv_classes[3][4][9] = int_stack + 6908;
 Libderiv->deriv_classes[3][5][9] = int_stack + 7058;
 Libderiv->deriv_classes[3][6][9] = int_stack + 7268;
 Libderiv->deriv_classes[4][3][9] = int_stack + 7548;
 Libderiv->deriv_classes[4][4][9] = int_stack + 7698;
 Libderiv->deriv_classes[4][5][9] = int_stack + 7923;
 Libderiv->deriv_classes[4][6][9] = int_stack + 8238;
 Libderiv->deriv_classes[5][3][9] = int_stack + 8658;
 Libderiv->deriv_classes[5][4][9] = int_stack + 8868;
 Libderiv->deriv_classes[5][5][9] = int_stack + 9183;
 Libderiv->deriv_classes[5][6][9] = int_stack + 9624;
 Libderiv->deriv_classes[3][3][8] = int_stack + 10212;
 Libderiv->deriv_classes[3][4][8] = int_stack + 10312;
 Libderiv->deriv_classes[3][5][8] = int_stack + 10462;
 Libderiv->deriv_classes[3][6][8] = int_stack + 10672;
 Libderiv->deriv_classes[4][3][8] = int_stack + 10952;
 Libderiv->deriv_classes[4][4][8] = int_stack + 11102;
 Libderiv->deriv_classes[4][5][8] = int_stack + 11327;
 Libderiv->deriv_classes[4][6][8] = int_stack + 11642;
 Libderiv->deriv_classes[5][3][8] = int_stack + 12062;
 Libderiv->deriv_classes[5][4][8] = int_stack + 12272;
 Libderiv->deriv_classes[5][5][8] = int_stack + 12587;
 Libderiv->deriv_classes[5][6][8] = int_stack + 13028;
 Libderiv->deriv_classes[3][3][7] = int_stack + 13616;
 Libderiv->deriv_classes[3][4][7] = int_stack + 13716;
 Libderiv->deriv_classes[3][5][7] = int_stack + 13866;
 Libderiv->deriv_classes[3][6][7] = int_stack + 14076;
 Libderiv->deriv_classes[4][3][7] = int_stack + 14356;
 Libderiv->deriv_classes[4][4][7] = int_stack + 14506;
 Libderiv->deriv_classes[4][5][7] = int_stack + 14731;
 Libderiv->deriv_classes[4][6][7] = int_stack + 15046;
 Libderiv->deriv_classes[5][3][7] = int_stack + 15466;
 Libderiv->deriv_classes[5][4][7] = int_stack + 15676;
 Libderiv->deriv_classes[5][5][7] = int_stack + 15991;
 Libderiv->deriv_classes[5][6][7] = int_stack + 16432;
 Libderiv->deriv_classes[3][3][6] = int_stack + 17020;
 Libderiv->deriv_classes[3][4][6] = int_stack + 17120;
 Libderiv->deriv_classes[3][5][6] = int_stack + 17270;
 Libderiv->deriv_classes[3][6][6] = int_stack + 17480;
 Libderiv->deriv_classes[4][3][6] = int_stack + 17760;
 Libderiv->deriv_classes[4][4][6] = int_stack + 17910;
 Libderiv->deriv_classes[4][5][6] = int_stack + 18135;
 Libderiv->deriv_classes[4][6][6] = int_stack + 18450;
 Libderiv->dvrr_classes[5][3] = int_stack + 18870;
 Libderiv->deriv_classes[5][3][6] = int_stack + 19080;
 Libderiv->dvrr_classes[5][4] = int_stack + 19290;
 Libderiv->deriv_classes[5][4][6] = int_stack + 19605;
 Libderiv->dvrr_classes[5][5] = int_stack + 19920;
 Libderiv->deriv_classes[5][5][6] = int_stack + 20361;
 Libderiv->deriv_classes[5][6][6] = int_stack + 20802;
 Libderiv->deriv_classes[3][3][2] = int_stack + 21390;
 Libderiv->deriv_classes[3][4][2] = int_stack + 21490;
 Libderiv->deriv_classes[3][5][2] = int_stack + 21640;
 Libderiv->deriv_classes[3][6][2] = int_stack + 21850;
 Libderiv->deriv_classes[4][3][2] = int_stack + 22130;
 Libderiv->deriv_classes[4][4][2] = int_stack + 22280;
 Libderiv->deriv_classes[4][5][2] = int_stack + 22505;
 Libderiv->deriv_classes[4][6][2] = int_stack + 22820;
 Libderiv->deriv_classes[5][3][2] = int_stack + 23240;
 Libderiv->deriv_classes[5][4][2] = int_stack + 23450;
 Libderiv->deriv_classes[5][5][2] = int_stack + 23765;
 Libderiv->deriv_classes[5][6][2] = int_stack + 24206;
 Libderiv->deriv_classes[3][3][1] = int_stack + 24794;
 Libderiv->deriv_classes[3][4][1] = int_stack + 24894;
 Libderiv->deriv_classes[3][5][1] = int_stack + 25044;
 Libderiv->deriv_classes[3][6][1] = int_stack + 25254;
 Libderiv->deriv_classes[4][3][1] = int_stack + 25534;
 Libderiv->deriv_classes[4][4][1] = int_stack + 25684;
 Libderiv->deriv_classes[4][5][1] = int_stack + 25909;
 Libderiv->deriv_classes[4][6][1] = int_stack + 26224;
 Libderiv->deriv_classes[5][3][1] = int_stack + 26644;
 Libderiv->deriv_classes[5][4][1] = int_stack + 26854;
 Libderiv->deriv_classes[5][5][1] = int_stack + 27169;
 Libderiv->deriv_classes[5][6][1] = int_stack + 27610;
 Libderiv->dvrr_classes[3][3] = int_stack + 28198;
 Libderiv->dvrr_classes[3][4] = int_stack + 28298;
 Libderiv->dvrr_classes[3][5] = int_stack + 28448;
 Libderiv->dvrr_classes[3][6] = int_stack + 28658;
 Libderiv->deriv_classes[3][3][0] = int_stack + 28938;
 Libderiv->deriv_classes[3][4][0] = int_stack + 29038;
 Libderiv->deriv_classes[3][5][0] = int_stack + 29188;
 Libderiv->deriv_classes[3][6][0] = int_stack + 29398;
 Libderiv->dvrr_classes[4][3] = int_stack + 29678;
 Libderiv->dvrr_classes[4][4] = int_stack + 29828;
 Libderiv->dvrr_classes[4][5] = int_stack + 30053;
 Libderiv->dvrr_classes[4][6] = int_stack + 30368;
 Libderiv->deriv_classes[4][3][0] = int_stack + 30788;
 Libderiv->deriv_classes[4][4][0] = int_stack + 30938;
 Libderiv->deriv_classes[4][5][0] = int_stack + 31163;
 Libderiv->deriv_classes[4][6][0] = int_stack + 31478;
 Libderiv->deriv_classes[5][3][0] = int_stack + 31898;
 Libderiv->deriv_classes[5][4][0] = int_stack + 32108;
 Libderiv->deriv_classes[5][5][0] = int_stack + 32423;
 Libderiv->deriv_classes[5][6][0] = int_stack + 32864;
 memset(int_stack,0,267616);

 Libderiv->dvrr_stack = int_stack + 87477;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fdff(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+33452,int_stack+28298,int_stack+28198,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+33752,int_stack+28448,int_stack+28298,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+34202,int_stack+33752,int_stack+33452,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+34802,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28198,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+35102,int_stack+250,int_stack+100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28298,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35552,int_stack+35102,int_stack+34802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33452,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36152,int_stack+460,int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28448,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+36782,int_stack+36152,int_stack+35102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33752,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+37682,int_stack+36782,int_stack+35552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34202,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+34802,int_stack+29828,int_stack+29678,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35252,int_stack+30053,int_stack+29828,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+35927,int_stack+35252,int_stack+34802,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36827,int_stack+890,int_stack+740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29678,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1115,int_stack+890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29828,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+38682,int_stack+0,int_stack+36827, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34802,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+39582,int_stack+1430,int_stack+1115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30053,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+40527,int_stack+39582,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35252,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+40527,int_stack+38682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35927,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+38682,int_stack+0,int_stack+37682,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+36827,int_stack+19290,int_stack+18870,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+37457,int_stack+19920,int_stack+19290,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+41682,int_stack+37457,int_stack+36827,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42942,int_stack+2060,int_stack+1850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18870,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43572,int_stack+2375,int_stack+2060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19290,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+44517,int_stack+43572,int_stack+42942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36827,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45777,int_stack+2816,int_stack+2375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19920,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1500,int_stack+45777,int_stack+43572, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37457,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+45777,int_stack+1500,int_stack+44517, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41682,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+47877,int_stack+45777,int_stack+0,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+3504,int_stack+3404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28198, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+300,int_stack+3654,int_stack+3504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28298, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+750,int_stack+300,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33452, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1350,int_stack+3864,int_stack+3654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28448, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1980,int_stack+1350,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33752, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2880,int_stack+1980,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34202, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+4294,int_stack+4144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29678, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+4519,int_stack+4294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29828, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1125,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34802, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42942,int_stack+4834,int_stack+4519, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30053, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3880,int_stack+42942,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35252, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+42942,int_stack+3880,int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35927, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+44442,int_stack+42942,int_stack+2880,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+5464,int_stack+5254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18870, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+630,int_stack+5779,int_stack+5464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19290, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1575,int_stack+630,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36827, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2835,int_stack+6220,int_stack+5779, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19920, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4158,int_stack+2835,int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37457, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+52377,int_stack+4158,int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41682, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+52377,int_stack+42942,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42942,int_stack+6908,int_stack+6808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28198, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43242,int_stack+7058,int_stack+6908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28298, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+43692,int_stack+43242,int_stack+42942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33452, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52377,int_stack+7268,int_stack+7058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28448, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+53007,int_stack+52377,int_stack+43242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33752, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+53907,int_stack+53007,int_stack+43692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34202, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+52377,int_stack+7698,int_stack+7548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29678, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52827,int_stack+7923,int_stack+7698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29828, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42942,int_stack+52827,int_stack+52377, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34802, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4500,int_stack+8238,int_stack+7923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30053, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5445,int_stack+4500,int_stack+52827, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35252, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+52377,int_stack+5445,int_stack+42942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35927, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+4500,int_stack+52377,int_stack+53907,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42942,int_stack+8868,int_stack+8658, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18870, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7500,int_stack+9183,int_stack+8868, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19290, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+53877,int_stack+7500,int_stack+42942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36827, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42942,int_stack+9624,int_stack+9183, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19920, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55137,int_stack+42942,int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37457, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7500,int_stack+55137,int_stack+53877, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41682, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+53877,int_stack+7500,int_stack+52377,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+52377,int_stack+10312,int_stack+10212, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52677,int_stack+10462,int_stack+10312, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+53127,int_stack+52677,int_stack+52377, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7500,int_stack+10672,int_stack+10462, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8130,int_stack+7500,int_stack+52677, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33752, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9030,int_stack+8130,int_stack+53127, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34202, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7500,int_stack+11102,int_stack+10952, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7950,int_stack+11327,int_stack+11102, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+52377,int_stack+7950,int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10030,int_stack+11642,int_stack+11327, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42942,int_stack+10030,int_stack+7950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10030,int_stack+42942,int_stack+52377, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+58377,int_stack+10030,int_stack+9030,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+52377,int_stack+12272,int_stack+12062, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42942,int_stack+12587,int_stack+12272, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7500,int_stack+42942,int_stack+52377, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36827, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52377,int_stack+13028,int_stack+12587, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11530,int_stack+52377,int_stack+42942, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37457, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+61377,int_stack+11530,int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+63477,int_stack+61377,int_stack+10030,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61377,int_stack+13716,int_stack+13616, 0.0, zero_stack, 1.0, int_stack+28198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61677,int_stack+13866,int_stack+13716, 0.0, zero_stack, 1.0, int_stack+28298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+62127,int_stack+61677,int_stack+61377, 0.0, zero_stack, 1.0, int_stack+33452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+62727,int_stack+14076,int_stack+13866, 0.0, zero_stack, 1.0, int_stack+28448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7500,int_stack+62727,int_stack+61677, 0.0, zero_stack, 1.0, int_stack+33752, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8400,int_stack+7500,int_stack+62127, 0.0, zero_stack, 1.0, int_stack+34202, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7500,int_stack+14506,int_stack+14356, 0.0, zero_stack, 1.0, int_stack+29678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61377,int_stack+14731,int_stack+14506, 0.0, zero_stack, 1.0, int_stack+29828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+62052,int_stack+61377,int_stack+7500, 0.0, zero_stack, 1.0, int_stack+34802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9400,int_stack+15046,int_stack+14731, 0.0, zero_stack, 1.0, int_stack+30053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10345,int_stack+9400,int_stack+61377, 0.0, zero_stack, 1.0, int_stack+35252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+42942,int_stack+10345,int_stack+62052, 0.0, zero_stack, 1.0, int_stack+35927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+9400,int_stack+42942,int_stack+8400,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61377,int_stack+15676,int_stack+15466, 0.0, zero_stack, 1.0, int_stack+18870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+62007,int_stack+15991,int_stack+15676, 0.0, zero_stack, 1.0, int_stack+19290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12400,int_stack+62007,int_stack+61377, 0.0, zero_stack, 1.0, int_stack+36827, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13660,int_stack+16432,int_stack+15991, 0.0, zero_stack, 1.0, int_stack+19920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14983,int_stack+13660,int_stack+62007, 0.0, zero_stack, 1.0, int_stack+37457, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+61377,int_stack+14983,int_stack+12400, 0.0, zero_stack, 1.0, int_stack+41682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+12400,int_stack+61377,int_stack+42942,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42942,int_stack+17120,int_stack+17020, 1.0, int_stack+28198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43242,int_stack+17270,int_stack+17120, 1.0, int_stack+28298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+43692,int_stack+43242,int_stack+42942, 1.0, int_stack+33452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+61377,int_stack+17480,int_stack+17270, 1.0, int_stack+28448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+62007,int_stack+61377,int_stack+43242, 1.0, int_stack+33752, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7500,int_stack+62007,int_stack+43692, 1.0, int_stack+34202, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+61377,int_stack+17910,int_stack+17760, 1.0, int_stack+29678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+61827,int_stack+18135,int_stack+17910, 1.0, int_stack+29828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8500,int_stack+61827,int_stack+61377, 1.0, int_stack+34802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+62502,int_stack+18450,int_stack+18135, 1.0, int_stack+30053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42942,int_stack+62502,int_stack+61827, 1.0, int_stack+35252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+52377,int_stack+42942,int_stack+8500, 1.0, int_stack+35927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+67977,int_stack+52377,int_stack+7500,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7500,int_stack+19605,int_stack+19080, 1.0, int_stack+18870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8130,int_stack+20361,int_stack+19605, 1.0, int_stack+19290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42942,int_stack+8130,int_stack+7500, 1.0, int_stack+36827, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+61377,int_stack+20802,int_stack+20361, 1.0, int_stack+19920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16900,int_stack+61377,int_stack+8130, 1.0, int_stack+37457, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+61377,int_stack+16900,int_stack+42942, 1.0, int_stack+41682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+70977,int_stack+61377,int_stack+52377,100);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52377,int_stack+28658,int_stack+28448,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+61377,int_stack+52377,int_stack+33752,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+52377,int_stack+61377,int_stack+34202,10);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+61377,int_stack+30368,int_stack+30053,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+41682,int_stack+61377,int_stack+35252,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+61377,int_stack+41682,int_stack+35927,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+16900,int_stack+61377,int_stack+52377,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+41682,int_stack+21490,int_stack+21390,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+41982,int_stack+21640,int_stack+21490,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+62877,int_stack+41982,int_stack+41682,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+42432,int_stack+21850,int_stack+21640,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+43062,int_stack+42432,int_stack+41982,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+41682,int_stack+43062,int_stack+62877,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+62877,int_stack+22280,int_stack+22130,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+42682,int_stack+22505,int_stack+22280,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+43357,int_stack+42682,int_stack+62877,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+19900,int_stack+22820,int_stack+22505,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+20845,int_stack+19900,int_stack+42682,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+7500,int_stack+20845,int_stack+43357,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+19900,int_stack+7500,int_stack+41682, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52377, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+41682,int_stack+23450,int_stack+23240,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+42312,int_stack+23765,int_stack+23450,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+33452,int_stack+42312,int_stack+41682,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+34712,int_stack+24206,int_stack+23765,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+22900,int_stack+34712,int_stack+42312,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+34712,int_stack+22900,int_stack+33452,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+75477,int_stack+34712,int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61377, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7500,int_stack+24894,int_stack+24794,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7800,int_stack+25044,int_stack+24894,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+62877,int_stack+7800,int_stack+7500,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+8250,int_stack+25254,int_stack+25044,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+33452,int_stack+8250,int_stack+7800,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+7500,int_stack+33452,int_stack+62877,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+62877,int_stack+25684,int_stack+25534,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+33452,int_stack+25909,int_stack+25684,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+8500,int_stack+33452,int_stack+62877,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+34127,int_stack+26224,int_stack+25909,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+35072,int_stack+34127,int_stack+33452,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+33452,int_stack+35072,int_stack+8500,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+34952,int_stack+33452,int_stack+7500, 0.0, zero_stack, 1.0, int_stack+52377, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+7500,int_stack+26854,int_stack+26644,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8130,int_stack+27169,int_stack+26854,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+22900,int_stack+8130,int_stack+7500,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24160,int_stack+27610,int_stack+27169,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+25483,int_stack+24160,int_stack+8130,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+41682,int_stack+25483,int_stack+22900,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+22900,int_stack+41682,int_stack+33452, 0.0, zero_stack, 1.0, int_stack+61377, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+33452,int_stack+29038,int_stack+28938,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+33752,int_stack+29188,int_stack+29038,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+62877,int_stack+33752,int_stack+33452,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+34202,int_stack+29398,int_stack+29188,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+41682,int_stack+34202,int_stack+33752,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+33452,int_stack+41682,int_stack+62877,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+62877,int_stack+30938,int_stack+30788,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+41682,int_stack+31163,int_stack+30938,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+42357,int_stack+41682,int_stack+62877,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+43257,int_stack+31478,int_stack+31163,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+27400,int_stack+43257,int_stack+41682,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+28750,int_stack+27400,int_stack+42357,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+79977,int_stack+28750,int_stack+33452, 1.0, int_stack+52377, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52377,int_stack+32108,int_stack+31898,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+33452,int_stack+32423,int_stack+32108,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+27400,int_stack+33452,int_stack+52377,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52377,int_stack+32864,int_stack+32423,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+41682,int_stack+52377,int_stack+33452,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+30250,int_stack+41682,int_stack+27400,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+82977,int_stack+30250,int_stack+28750, 1.0, int_stack+61377, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+27400,int_stack+47877,int_stack+38682,100);
     Libderiv->ABCD[11] = int_stack + 27400;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+37952,int_stack+0,int_stack+44442,100);
     Libderiv->ABCD[10] = int_stack + 37952;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+43952,int_stack+53877,int_stack+4500,100);
     Libderiv->ABCD[9] = int_stack + 43952;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+0,int_stack+63477,int_stack+58377,100);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+49952,int_stack+12400,int_stack+9400,100);
     Libderiv->ABCD[7] = int_stack + 49952;
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+6000,int_stack+70977,int_stack+67977,100);
     Libderiv->ABCD[6] = int_stack + 6000;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+55952,int_stack+75477,int_stack+19900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 55952;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+61952,int_stack+22900,int_stack+34952, 0.0, zero_stack, 1.0, int_stack+16900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 61952;
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+19900,int_stack+82977,int_stack+79977, 1.0, int_stack+16900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 19900;

}
