#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gpgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gp|gf) integrals */

void d1hrr_order_gpgf(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 540;
 Libderiv->deriv_classes[4][7][11] = int_stack + 960;
 Libderiv->deriv_classes[5][4][11] = int_stack + 1500;
 Libderiv->deriv_classes[5][5][11] = int_stack + 1815;
 Libderiv->deriv_classes[5][6][11] = int_stack + 2256;
 Libderiv->deriv_classes[5][7][11] = int_stack + 2844;
 Libderiv->deriv_classes[4][4][10] = int_stack + 3600;
 Libderiv->deriv_classes[4][5][10] = int_stack + 3825;
 Libderiv->deriv_classes[4][6][10] = int_stack + 4140;
 Libderiv->deriv_classes[4][7][10] = int_stack + 4560;
 Libderiv->deriv_classes[5][4][10] = int_stack + 5100;
 Libderiv->deriv_classes[5][5][10] = int_stack + 5415;
 Libderiv->deriv_classes[5][6][10] = int_stack + 5856;
 Libderiv->deriv_classes[5][7][10] = int_stack + 6444;
 Libderiv->deriv_classes[4][4][9] = int_stack + 7200;
 Libderiv->deriv_classes[4][5][9] = int_stack + 7425;
 Libderiv->deriv_classes[4][6][9] = int_stack + 7740;
 Libderiv->deriv_classes[4][7][9] = int_stack + 8160;
 Libderiv->deriv_classes[5][4][9] = int_stack + 8700;
 Libderiv->deriv_classes[5][5][9] = int_stack + 9015;
 Libderiv->deriv_classes[5][6][9] = int_stack + 9456;
 Libderiv->deriv_classes[5][7][9] = int_stack + 10044;
 Libderiv->deriv_classes[4][4][8] = int_stack + 10800;
 Libderiv->deriv_classes[4][5][8] = int_stack + 11025;
 Libderiv->deriv_classes[4][6][8] = int_stack + 11340;
 Libderiv->deriv_classes[4][7][8] = int_stack + 11760;
 Libderiv->deriv_classes[5][4][8] = int_stack + 12300;
 Libderiv->deriv_classes[5][5][8] = int_stack + 12615;
 Libderiv->deriv_classes[5][6][8] = int_stack + 13056;
 Libderiv->deriv_classes[5][7][8] = int_stack + 13644;
 Libderiv->deriv_classes[4][4][7] = int_stack + 14400;
 Libderiv->deriv_classes[4][5][7] = int_stack + 14625;
 Libderiv->deriv_classes[4][6][7] = int_stack + 14940;
 Libderiv->deriv_classes[4][7][7] = int_stack + 15360;
 Libderiv->deriv_classes[5][4][7] = int_stack + 15900;
 Libderiv->deriv_classes[5][5][7] = int_stack + 16215;
 Libderiv->deriv_classes[5][6][7] = int_stack + 16656;
 Libderiv->deriv_classes[5][7][7] = int_stack + 17244;
 Libderiv->deriv_classes[4][4][6] = int_stack + 18000;
 Libderiv->deriv_classes[4][5][6] = int_stack + 18225;
 Libderiv->deriv_classes[4][6][6] = int_stack + 18540;
 Libderiv->deriv_classes[4][7][6] = int_stack + 18960;
 Libderiv->dvrr_classes[5][4] = int_stack + 19500;
 Libderiv->deriv_classes[5][4][6] = int_stack + 19815;
 Libderiv->dvrr_classes[5][5] = int_stack + 20130;
 Libderiv->deriv_classes[5][5][6] = int_stack + 20571;
 Libderiv->dvrr_classes[5][6] = int_stack + 21012;
 Libderiv->deriv_classes[5][6][6] = int_stack + 21600;
 Libderiv->deriv_classes[5][7][6] = int_stack + 22188;
 Libderiv->deriv_classes[4][4][2] = int_stack + 22944;
 Libderiv->deriv_classes[4][5][2] = int_stack + 23169;
 Libderiv->deriv_classes[4][6][2] = int_stack + 23484;
 Libderiv->deriv_classes[4][7][2] = int_stack + 23904;
 Libderiv->deriv_classes[5][4][2] = int_stack + 24444;
 Libderiv->deriv_classes[5][5][2] = int_stack + 24759;
 Libderiv->deriv_classes[5][6][2] = int_stack + 25200;
 Libderiv->deriv_classes[5][7][2] = int_stack + 25788;
 Libderiv->deriv_classes[4][4][1] = int_stack + 26544;
 Libderiv->deriv_classes[4][5][1] = int_stack + 26769;
 Libderiv->deriv_classes[4][6][1] = int_stack + 27084;
 Libderiv->deriv_classes[4][7][1] = int_stack + 27504;
 Libderiv->deriv_classes[5][4][1] = int_stack + 28044;
 Libderiv->deriv_classes[5][5][1] = int_stack + 28359;
 Libderiv->deriv_classes[5][6][1] = int_stack + 28800;
 Libderiv->deriv_classes[5][7][1] = int_stack + 29388;
 Libderiv->dvrr_classes[4][4] = int_stack + 30144;
 Libderiv->dvrr_classes[4][5] = int_stack + 30369;
 Libderiv->dvrr_classes[4][6] = int_stack + 30684;
 Libderiv->dvrr_classes[4][7] = int_stack + 31104;
 Libderiv->deriv_classes[4][4][0] = int_stack + 31644;
 Libderiv->deriv_classes[4][5][0] = int_stack + 31869;
 Libderiv->deriv_classes[4][6][0] = int_stack + 32184;
 Libderiv->deriv_classes[4][7][0] = int_stack + 32604;
 Libderiv->deriv_classes[5][4][0] = int_stack + 33144;
 Libderiv->deriv_classes[5][5][0] = int_stack + 33459;
 Libderiv->deriv_classes[5][6][0] = int_stack + 33900;
 Libderiv->deriv_classes[5][7][0] = int_stack + 34488;
 memset(int_stack,0,281952);

 Libderiv->dvrr_stack = int_stack + 86418;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gpgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35244,int_stack+30369,int_stack+30144,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+35919,int_stack+30684,int_stack+30369,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+36864,int_stack+35919,int_stack+35244,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+38214,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30144,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+38889,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30369,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+39834,int_stack+38889,int_stack+38214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35244,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+41184,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30684,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+42444,int_stack+41184,int_stack+38889, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35919,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+44334,int_stack+42444,int_stack+39834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36864,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+38214,int_stack+20130,int_stack+19500,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+39159,int_stack+21012,int_stack+20130,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+40482,int_stack+39159,int_stack+38214,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42372,int_stack+1815,int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19500,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+2256,int_stack+1815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20130,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+46584,int_stack+0,int_stack+42372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38214,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+42372,int_stack+2844,int_stack+2256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21012,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+48474,int_stack+42372,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39159,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+48474,int_stack+46584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40482,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+46584,int_stack+3825,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30144, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47259,int_stack+4140,int_stack+3825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30369, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48204,int_stack+47259,int_stack+46584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35244, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+49554,int_stack+4560,int_stack+4140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30684, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+3150,int_stack+49554,int_stack+47259, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35919, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+49554,int_stack+3150,int_stack+48204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36864, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3150,int_stack+5415,int_stack+5100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19500, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+46584,int_stack+5856,int_stack+5415, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20130, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42372,int_stack+46584,int_stack+3150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38214, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3150,int_stack+6444,int_stack+5856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21012, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+51804,int_stack+3150,int_stack+46584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39159, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3150,int_stack+51804,int_stack+42372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40482, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42372,int_stack+7425,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30144, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43047,int_stack+7740,int_stack+7425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30369, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+51804,int_stack+43047,int_stack+42372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35244, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+53154,int_stack+8160,int_stack+7740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30684, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+46584,int_stack+53154,int_stack+43047, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35919, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6300,int_stack+46584,int_stack+51804, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36864, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+51804,int_stack+9015,int_stack+8700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19500, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+52749,int_stack+9456,int_stack+9015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20130, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+46584,int_stack+52749,int_stack+51804, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38214, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+42372,int_stack+10044,int_stack+9456, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21012, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+54072,int_stack+42372,int_stack+52749, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39159, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+56718,int_stack+54072,int_stack+46584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40482, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+46584,int_stack+11025,int_stack+10800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47259,int_stack+11340,int_stack+11025, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48204,int_stack+47259,int_stack+46584, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+42372,int_stack+11760,int_stack+11340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+51804,int_stack+42372,int_stack+47259, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+53694,int_stack+51804,int_stack+48204, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+51804,int_stack+12615,int_stack+12300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+42372,int_stack+13056,int_stack+12615, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+46584,int_stack+42372,int_stack+51804, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+51804,int_stack+13644,int_stack+13056, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8550,int_stack+51804,int_stack+42372, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39159, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+11196,int_stack+8550,int_stack+46584, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+46584,int_stack+14625,int_stack+14400, 0.0, zero_stack, 1.0, int_stack+30144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47259,int_stack+14940,int_stack+14625, 0.0, zero_stack, 1.0, int_stack+30369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+48204,int_stack+47259,int_stack+46584, 0.0, zero_stack, 1.0, int_stack+35244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8550,int_stack+15360,int_stack+14940, 0.0, zero_stack, 1.0, int_stack+30684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+51804,int_stack+8550,int_stack+47259, 0.0, zero_stack, 1.0, int_stack+35919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+8550,int_stack+51804,int_stack+48204, 0.0, zero_stack, 1.0, int_stack+36864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+51804,int_stack+16215,int_stack+15900, 0.0, zero_stack, 1.0, int_stack+19500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+46584,int_stack+16656,int_stack+16215, 0.0, zero_stack, 1.0, int_stack+20130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42372,int_stack+46584,int_stack+51804, 0.0, zero_stack, 1.0, int_stack+38214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+51804,int_stack+17244,int_stack+16656, 0.0, zero_stack, 1.0, int_stack+21012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+14346,int_stack+51804,int_stack+46584, 0.0, zero_stack, 1.0, int_stack+39159, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+59868,int_stack+14346,int_stack+42372, 0.0, zero_stack, 1.0, int_stack+40482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42372,int_stack+18225,int_stack+18000, 1.0, int_stack+30144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43047,int_stack+18540,int_stack+18225, 1.0, int_stack+30369, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14346,int_stack+43047,int_stack+42372, 1.0, int_stack+35244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+15696,int_stack+18960,int_stack+18540, 1.0, int_stack+30684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+51804,int_stack+15696,int_stack+43047, 1.0, int_stack+35919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+15696,int_stack+51804,int_stack+14346, 1.0, int_stack+36864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14346,int_stack+20571,int_stack+19815, 1.0, int_stack+19500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+51804,int_stack+21600,int_stack+20571, 1.0, int_stack+20130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17946,int_stack+51804,int_stack+14346, 1.0, int_stack+38214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+42372,int_stack+22188,int_stack+21600, 1.0, int_stack+21012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+19836,int_stack+42372,int_stack+51804, 1.0, int_stack+39159, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+63018,int_stack+19836,int_stack+17946, 1.0, int_stack+40482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+17946,int_stack+31104,int_stack+30684,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+51804,int_stack+17946,int_stack+35919,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+17946,int_stack+51804,int_stack+36864,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51804,int_stack+23169,int_stack+22944,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52479,int_stack+23484,int_stack+23169,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+14346,int_stack+52479,int_stack+51804,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+20196,int_stack+23904,int_stack+23484,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+21456,int_stack+20196,int_stack+52479,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+35244,int_stack+21456,int_stack+14346,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14346,int_stack+24759,int_stack+24444,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+20196,int_stack+25200,int_stack+24759,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+51804,int_stack+20196,int_stack+14346,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+21519,int_stack+25788,int_stack+25200,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+23283,int_stack+21519,int_stack+20196,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+37494,int_stack+23283,int_stack+51804,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51804,int_stack+26769,int_stack+26544,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52479,int_stack+27084,int_stack+26769,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+14346,int_stack+52479,int_stack+51804,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+20196,int_stack+27504,int_stack+27084,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+21456,int_stack+20196,int_stack+52479,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+23346,int_stack+21456,int_stack+14346,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14346,int_stack+28359,int_stack+28044,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+20196,int_stack+28800,int_stack+28359,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+51804,int_stack+20196,int_stack+14346,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+21519,int_stack+29388,int_stack+28800,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+25596,int_stack+21519,int_stack+20196,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+20196,int_stack+25596,int_stack+51804,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51804,int_stack+31869,int_stack+31644,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52479,int_stack+32184,int_stack+31869,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+14346,int_stack+52479,int_stack+51804,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+25596,int_stack+32604,int_stack+32184,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+26856,int_stack+25596,int_stack+52479,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+28746,int_stack+26856,int_stack+14346,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14346,int_stack+33459,int_stack+33144,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+25596,int_stack+33900,int_stack+33459,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+51804,int_stack+25596,int_stack+14346,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+26919,int_stack+34488,int_stack+33900,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+30996,int_stack+26919,int_stack+25596,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+25596,int_stack+30996,int_stack+51804,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+66168,int_stack+0,int_stack+44334,150);
     Libderiv->ABCD[11] = int_stack + 66168;
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+40644,int_stack+3150,int_stack+49554,150);
     Libderiv->ABCD[10] = int_stack + 40644;
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+72918,int_stack+56718,int_stack+6300,150);
     Libderiv->ABCD[9] = int_stack + 72918;
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+11196,int_stack+53694,150);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+47394,int_stack+59868,int_stack+8550,150);
     Libderiv->ABCD[7] = int_stack + 47394;
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+6750,int_stack+63018,int_stack+15696,150);
     Libderiv->ABCD[6] = int_stack + 6750;
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+54144,int_stack+37494,int_stack+35244, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 54144;
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+30996,int_stack+20196,int_stack+23346, 0.0, zero_stack, 1.0, int_stack+17946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 30996;
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+79668,int_stack+25596,int_stack+28746, 1.0, int_stack+17946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 79668;

}
