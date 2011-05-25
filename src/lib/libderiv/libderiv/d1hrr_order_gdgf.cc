#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gdgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gd|gf) integrals */

void d1hrr_order_gdgf(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[6][4][11] = int_stack + 3600;
 Libderiv->deriv_classes[6][5][11] = int_stack + 4020;
 Libderiv->deriv_classes[6][6][11] = int_stack + 4608;
 Libderiv->deriv_classes[6][7][11] = int_stack + 5392;
 Libderiv->deriv_classes[4][4][10] = int_stack + 6400;
 Libderiv->deriv_classes[4][5][10] = int_stack + 6625;
 Libderiv->deriv_classes[4][6][10] = int_stack + 6940;
 Libderiv->deriv_classes[4][7][10] = int_stack + 7360;
 Libderiv->deriv_classes[5][4][10] = int_stack + 7900;
 Libderiv->deriv_classes[5][5][10] = int_stack + 8215;
 Libderiv->deriv_classes[5][6][10] = int_stack + 8656;
 Libderiv->deriv_classes[5][7][10] = int_stack + 9244;
 Libderiv->deriv_classes[6][4][10] = int_stack + 10000;
 Libderiv->deriv_classes[6][5][10] = int_stack + 10420;
 Libderiv->deriv_classes[6][6][10] = int_stack + 11008;
 Libderiv->deriv_classes[6][7][10] = int_stack + 11792;
 Libderiv->deriv_classes[4][4][9] = int_stack + 12800;
 Libderiv->deriv_classes[4][5][9] = int_stack + 13025;
 Libderiv->deriv_classes[4][6][9] = int_stack + 13340;
 Libderiv->deriv_classes[4][7][9] = int_stack + 13760;
 Libderiv->deriv_classes[5][4][9] = int_stack + 14300;
 Libderiv->deriv_classes[5][5][9] = int_stack + 14615;
 Libderiv->deriv_classes[5][6][9] = int_stack + 15056;
 Libderiv->deriv_classes[5][7][9] = int_stack + 15644;
 Libderiv->deriv_classes[6][4][9] = int_stack + 16400;
 Libderiv->deriv_classes[6][5][9] = int_stack + 16820;
 Libderiv->deriv_classes[6][6][9] = int_stack + 17408;
 Libderiv->deriv_classes[6][7][9] = int_stack + 18192;
 Libderiv->deriv_classes[4][4][8] = int_stack + 19200;
 Libderiv->deriv_classes[4][5][8] = int_stack + 19425;
 Libderiv->deriv_classes[4][6][8] = int_stack + 19740;
 Libderiv->deriv_classes[4][7][8] = int_stack + 20160;
 Libderiv->deriv_classes[5][4][8] = int_stack + 20700;
 Libderiv->deriv_classes[5][5][8] = int_stack + 21015;
 Libderiv->deriv_classes[5][6][8] = int_stack + 21456;
 Libderiv->deriv_classes[5][7][8] = int_stack + 22044;
 Libderiv->deriv_classes[6][4][8] = int_stack + 22800;
 Libderiv->deriv_classes[6][5][8] = int_stack + 23220;
 Libderiv->deriv_classes[6][6][8] = int_stack + 23808;
 Libderiv->deriv_classes[6][7][8] = int_stack + 24592;
 Libderiv->deriv_classes[4][4][7] = int_stack + 25600;
 Libderiv->deriv_classes[4][5][7] = int_stack + 25825;
 Libderiv->deriv_classes[4][6][7] = int_stack + 26140;
 Libderiv->deriv_classes[4][7][7] = int_stack + 26560;
 Libderiv->deriv_classes[5][4][7] = int_stack + 27100;
 Libderiv->deriv_classes[5][5][7] = int_stack + 27415;
 Libderiv->deriv_classes[5][6][7] = int_stack + 27856;
 Libderiv->deriv_classes[5][7][7] = int_stack + 28444;
 Libderiv->deriv_classes[6][4][7] = int_stack + 29200;
 Libderiv->deriv_classes[6][5][7] = int_stack + 29620;
 Libderiv->deriv_classes[6][6][7] = int_stack + 30208;
 Libderiv->deriv_classes[6][7][7] = int_stack + 30992;
 Libderiv->deriv_classes[4][4][6] = int_stack + 32000;
 Libderiv->deriv_classes[4][5][6] = int_stack + 32225;
 Libderiv->deriv_classes[4][6][6] = int_stack + 32540;
 Libderiv->deriv_classes[4][7][6] = int_stack + 32960;
 Libderiv->deriv_classes[5][4][6] = int_stack + 33500;
 Libderiv->deriv_classes[5][5][6] = int_stack + 33815;
 Libderiv->deriv_classes[5][6][6] = int_stack + 34256;
 Libderiv->deriv_classes[5][7][6] = int_stack + 34844;
 Libderiv->dvrr_classes[6][4] = int_stack + 35600;
 Libderiv->deriv_classes[6][4][6] = int_stack + 36020;
 Libderiv->dvrr_classes[6][5] = int_stack + 36440;
 Libderiv->deriv_classes[6][5][6] = int_stack + 37028;
 Libderiv->dvrr_classes[6][6] = int_stack + 37616;
 Libderiv->deriv_classes[6][6][6] = int_stack + 38400;
 Libderiv->deriv_classes[6][7][6] = int_stack + 39184;
 Libderiv->deriv_classes[4][4][2] = int_stack + 40192;
 Libderiv->deriv_classes[4][5][2] = int_stack + 40417;
 Libderiv->deriv_classes[4][6][2] = int_stack + 40732;
 Libderiv->deriv_classes[4][7][2] = int_stack + 41152;
 Libderiv->deriv_classes[5][4][2] = int_stack + 41692;
 Libderiv->deriv_classes[5][5][2] = int_stack + 42007;
 Libderiv->deriv_classes[5][6][2] = int_stack + 42448;
 Libderiv->deriv_classes[5][7][2] = int_stack + 43036;
 Libderiv->deriv_classes[6][4][2] = int_stack + 43792;
 Libderiv->deriv_classes[6][5][2] = int_stack + 44212;
 Libderiv->deriv_classes[6][6][2] = int_stack + 44800;
 Libderiv->deriv_classes[6][7][2] = int_stack + 45584;
 Libderiv->deriv_classes[4][4][1] = int_stack + 46592;
 Libderiv->deriv_classes[4][5][1] = int_stack + 46817;
 Libderiv->deriv_classes[4][6][1] = int_stack + 47132;
 Libderiv->deriv_classes[4][7][1] = int_stack + 47552;
 Libderiv->deriv_classes[5][4][1] = int_stack + 48092;
 Libderiv->deriv_classes[5][5][1] = int_stack + 48407;
 Libderiv->deriv_classes[5][6][1] = int_stack + 48848;
 Libderiv->deriv_classes[5][7][1] = int_stack + 49436;
 Libderiv->deriv_classes[6][4][1] = int_stack + 50192;
 Libderiv->deriv_classes[6][5][1] = int_stack + 50612;
 Libderiv->deriv_classes[6][6][1] = int_stack + 51200;
 Libderiv->deriv_classes[6][7][1] = int_stack + 51984;
 Libderiv->dvrr_classes[4][4] = int_stack + 52992;
 Libderiv->dvrr_classes[4][5] = int_stack + 53217;
 Libderiv->dvrr_classes[4][6] = int_stack + 53532;
 Libderiv->dvrr_classes[4][7] = int_stack + 53952;
 Libderiv->deriv_classes[4][4][0] = int_stack + 54492;
 Libderiv->deriv_classes[4][5][0] = int_stack + 54717;
 Libderiv->deriv_classes[4][6][0] = int_stack + 55032;
 Libderiv->deriv_classes[4][7][0] = int_stack + 55452;
 Libderiv->dvrr_classes[5][4] = int_stack + 55992;
 Libderiv->dvrr_classes[5][5] = int_stack + 56307;
 Libderiv->dvrr_classes[5][6] = int_stack + 56748;
 Libderiv->dvrr_classes[5][7] = int_stack + 57336;
 Libderiv->deriv_classes[5][4][0] = int_stack + 58092;
 Libderiv->deriv_classes[5][5][0] = int_stack + 58407;
 Libderiv->deriv_classes[5][6][0] = int_stack + 58848;
 Libderiv->deriv_classes[5][7][0] = int_stack + 59436;
 Libderiv->deriv_classes[6][4][0] = int_stack + 60192;
 Libderiv->deriv_classes[6][5][0] = int_stack + 60612;
 Libderiv->deriv_classes[6][6][0] = int_stack + 61200;
 Libderiv->deriv_classes[6][7][0] = int_stack + 61984;
 memset(int_stack,0,503936);

 Libderiv->dvrr_stack = int_stack + 202372;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gdgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+62992,int_stack+53217,int_stack+52992,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+63667,int_stack+53532,int_stack+53217,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+64612,int_stack+63667,int_stack+62992,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65962,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52992,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66637,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53217,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67582,int_stack+66637,int_stack+65962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+68932,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53532,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+70192,int_stack+68932,int_stack+66637, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63667,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+72082,int_stack+70192,int_stack+67582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64612,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+65962,int_stack+56307,int_stack+55992,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+66907,int_stack+56748,int_stack+56307,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+68230,int_stack+66907,int_stack+65962,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+70120,int_stack+1815,int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55992,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+2256,int_stack+1815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56307,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+74332,int_stack+0,int_stack+70120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65962,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+70120,int_stack+2844,int_stack+2256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56748,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+76222,int_stack+70120,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66907,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+76222,int_stack+74332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68230,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+74332,int_stack+0,int_stack+72082,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+70120,int_stack+36440,int_stack+35600,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+71380,int_stack+37616,int_stack+36440,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81082,int_stack+71380,int_stack+70120,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83602,int_stack+4020,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35600,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84862,int_stack+4608,int_stack+4020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36440,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+86626,int_stack+84862,int_stack+83602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70120,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+89146,int_stack+5392,int_stack+4608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37616,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+91498,int_stack+89146,int_stack+84862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71380,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+95026,int_stack+91498,int_stack+86626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81082,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+83602,int_stack+95026,int_stack+0,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+6625,int_stack+6400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52992, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+675,int_stack+6940,int_stack+6625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53217, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1620,int_stack+675,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+2970,int_stack+7360,int_stack+6940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53532, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+4230,int_stack+2970,int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63667, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+93052,int_stack+4230,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64612, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+8215,int_stack+7900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55992, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+945,int_stack+8656,int_stack+8215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56307, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2268,int_stack+945,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65962, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+4158,int_stack+9244,int_stack+8656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56748, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+5922,int_stack+4158,int_stack+945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66907, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+95302,int_stack+5922,int_stack+2268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68230, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+95302,int_stack+93052,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+93052,int_stack+10420,int_stack+10000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35600, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6750,int_stack+11008,int_stack+10420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36440, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+98452,int_stack+6750,int_stack+93052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70120, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8514,int_stack+11792,int_stack+11008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37616, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+100972,int_stack+8514,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71380, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+6750,int_stack+100972,int_stack+98452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81082, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+98452,int_stack+6750,int_stack+95302,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6750,int_stack+13025,int_stack+12800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52992, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7425,int_stack+13340,int_stack+13025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53217, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8370,int_stack+7425,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9720,int_stack+13760,int_stack+13340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53532, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+10980,int_stack+9720,int_stack+7425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63667, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+93052,int_stack+10980,int_stack+8370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64612, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6750,int_stack+14615,int_stack+14300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55992, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7695,int_stack+15056,int_stack+14615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56307, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9018,int_stack+7695,int_stack+6750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65962, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+10908,int_stack+15644,int_stack+15056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56748, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+12672,int_stack+10908,int_stack+7695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66907, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+95302,int_stack+12672,int_stack+9018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68230, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+6750,int_stack+95302,int_stack+93052,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+93052,int_stack+16820,int_stack+16400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35600, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13500,int_stack+17408,int_stack+16820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36440, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+107902,int_stack+13500,int_stack+93052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70120, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+110422,int_stack+18192,int_stack+17408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37616, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+15264,int_stack+110422,int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71380, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+110422,int_stack+15264,int_stack+107902, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81082, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+114622,int_stack+110422,int_stack+95302,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107902,int_stack+19425,int_stack+19200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108577,int_stack+19740,int_stack+19425, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109522,int_stack+108577,int_stack+107902, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+110872,int_stack+20160,int_stack+19740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+112132,int_stack+110872,int_stack+108577, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63667, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+13500,int_stack+112132,int_stack+109522, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107902,int_stack+21015,int_stack+20700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108847,int_stack+21456,int_stack+21015, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+110170,int_stack+108847,int_stack+107902, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+112060,int_stack+22044,int_stack+21456, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56748, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+15750,int_stack+112060,int_stack+108847, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+18396,int_stack+15750,int_stack+110170, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+124072,int_stack+18396,int_stack+13500,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13500,int_stack+23220,int_stack+22800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14760,int_stack+23808,int_stack+23220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+107902,int_stack+14760,int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+110422,int_stack+24592,int_stack+23808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21546,int_stack+110422,int_stack+14760, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+110422,int_stack+21546,int_stack+107902, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+130822,int_stack+110422,int_stack+18396,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107902,int_stack+25825,int_stack+25600, 0.0, zero_stack, 1.0, int_stack+52992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108577,int_stack+26140,int_stack+25825, 0.0, zero_stack, 1.0, int_stack+53217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109522,int_stack+108577,int_stack+107902, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+110872,int_stack+26560,int_stack+26140, 0.0, zero_stack, 1.0, int_stack+53532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+112132,int_stack+110872,int_stack+108577, 0.0, zero_stack, 1.0, int_stack+63667, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+13500,int_stack+112132,int_stack+109522, 0.0, zero_stack, 1.0, int_stack+64612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107902,int_stack+27415,int_stack+27100, 0.0, zero_stack, 1.0, int_stack+55992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108847,int_stack+27856,int_stack+27415, 0.0, zero_stack, 1.0, int_stack+56307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+110170,int_stack+108847,int_stack+107902, 0.0, zero_stack, 1.0, int_stack+65962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+112060,int_stack+28444,int_stack+27856, 0.0, zero_stack, 1.0, int_stack+56748, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+15750,int_stack+112060,int_stack+108847, 0.0, zero_stack, 1.0, int_stack+66907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+18396,int_stack+15750,int_stack+110170, 0.0, zero_stack, 1.0, int_stack+68230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+21546,int_stack+18396,int_stack+13500,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13500,int_stack+29620,int_stack+29200, 0.0, zero_stack, 1.0, int_stack+35600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14760,int_stack+30208,int_stack+29620, 0.0, zero_stack, 1.0, int_stack+36440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+107902,int_stack+14760,int_stack+13500, 0.0, zero_stack, 1.0, int_stack+70120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+110422,int_stack+30992,int_stack+30208, 0.0, zero_stack, 1.0, int_stack+37616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+28296,int_stack+110422,int_stack+14760, 0.0, zero_stack, 1.0, int_stack+71380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+110422,int_stack+28296,int_stack+107902, 0.0, zero_stack, 1.0, int_stack+81082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+140272,int_stack+110422,int_stack+18396,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107902,int_stack+32225,int_stack+32000, 1.0, int_stack+52992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108577,int_stack+32540,int_stack+32225, 1.0, int_stack+53217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+109522,int_stack+108577,int_stack+107902, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+110872,int_stack+32960,int_stack+32540, 1.0, int_stack+53532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+112132,int_stack+110872,int_stack+108577, 1.0, int_stack+63667, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+28296,int_stack+112132,int_stack+109522, 1.0, int_stack+64612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107902,int_stack+33815,int_stack+33500, 1.0, int_stack+55992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+108847,int_stack+34256,int_stack+33815, 1.0, int_stack+56307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+110170,int_stack+108847,int_stack+107902, 1.0, int_stack+65962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+112060,int_stack+34844,int_stack+34256, 1.0, int_stack+56748, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+30546,int_stack+112060,int_stack+108847, 1.0, int_stack+66907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+13500,int_stack+30546,int_stack+110170, 1.0, int_stack+68230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+149722,int_stack+13500,int_stack+28296,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+28296,int_stack+37028,int_stack+36020, 1.0, int_stack+35600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29556,int_stack+38400,int_stack+37028, 1.0, int_stack+36440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31320,int_stack+29556,int_stack+28296, 1.0, int_stack+70120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+33840,int_stack+39184,int_stack+38400, 1.0, int_stack+37616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+36192,int_stack+33840,int_stack+29556, 1.0, int_stack+71380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+70120,int_stack+36192,int_stack+31320, 1.0, int_stack+81082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+28296,int_stack+70120,int_stack+13500,150);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+13500,int_stack+53952,int_stack+53532,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+14760,int_stack+13500,int_stack+63667,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+16650,int_stack+14760,int_stack+64612,15);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+13500,int_stack+57336,int_stack+56748,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+18900,int_stack+13500,int_stack+66907,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+13500,int_stack+18900,int_stack+68230,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+62992,int_stack+13500,int_stack+16650,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18900,int_stack+40417,int_stack+40192,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+19575,int_stack+40732,int_stack+40417,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81082,int_stack+19575,int_stack+18900,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+37746,int_stack+41152,int_stack+40732,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+39006,int_stack+37746,int_stack+19575,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+18900,int_stack+39006,int_stack+81082,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+81082,int_stack+42007,int_stack+41692,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+82027,int_stack+42448,int_stack+42007,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+37746,int_stack+82027,int_stack+81082,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+39636,int_stack+43036,int_stack+42448,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+107902,int_stack+39636,int_stack+82027,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+39636,int_stack+107902,int_stack+37746,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+156472,int_stack+39636,int_stack+18900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+18900,int_stack+44212,int_stack+43792,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+37746,int_stack+44800,int_stack+44212,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81082,int_stack+37746,int_stack+18900,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+18900,int_stack+45584,int_stack+44800,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+107902,int_stack+18900,int_stack+37746,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+69742,int_stack+107902,int_stack+81082,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+163222,int_stack+69742,int_stack+39636, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+69742,int_stack+46817,int_stack+46592,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+70417,int_stack+47132,int_stack+46817,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+71362,int_stack+70417,int_stack+69742,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+72712,int_stack+47552,int_stack+47132,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+81082,int_stack+72712,int_stack+70417,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+107902,int_stack+81082,int_stack+71362,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+81082,int_stack+48407,int_stack+48092,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+82027,int_stack+48848,int_stack+48407,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+69742,int_stack+82027,int_stack+81082,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+71632,int_stack+49436,int_stack+48848,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+18900,int_stack+71632,int_stack+82027,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+110152,int_stack+18900,int_stack+69742,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+37746,int_stack+110152,int_stack+107902, 0.0, zero_stack, 1.0, int_stack+16650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+107902,int_stack+50612,int_stack+50192,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+69742,int_stack+51200,int_stack+50612,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81082,int_stack+69742,int_stack+107902,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+71506,int_stack+51984,int_stack+51200,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+44496,int_stack+71506,int_stack+69742,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+69742,int_stack+44496,int_stack+81082,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+44496,int_stack+69742,int_stack+110152, 0.0, zero_stack, 1.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+69742,int_stack+54717,int_stack+54492,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+70417,int_stack+55032,int_stack+54717,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+71362,int_stack+70417,int_stack+69742,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+72712,int_stack+55452,int_stack+55032,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+81082,int_stack+72712,int_stack+70417,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+53946,int_stack+81082,int_stack+71362,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+81082,int_stack+58407,int_stack+58092,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+82027,int_stack+58848,int_stack+58407,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+69742,int_stack+82027,int_stack+81082,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+71632,int_stack+59436,int_stack+58848,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+18900,int_stack+71632,int_stack+82027,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+56196,int_stack+18900,int_stack+69742,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+172672,int_stack+56196,int_stack+53946, 1.0, int_stack+16650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16650,int_stack+60612,int_stack+60192,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+17910,int_stack+61200,int_stack+60612,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+81082,int_stack+17910,int_stack+16650,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+69742,int_stack+61984,int_stack+61200,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+59346,int_stack+69742,int_stack+17910,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+69742,int_stack+59346,int_stack+81082,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+179422,int_stack+69742,int_stack+56196, 1.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+188872,int_stack+83602,int_stack+74332,150);
     Libderiv->ABCD[11] = int_stack + 188872;
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+69742,int_stack+98452,int_stack+0,150);
     Libderiv->ABCD[10] = int_stack + 69742;
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+83242,int_stack+114622,int_stack+6750,150);
     Libderiv->ABCD[9] = int_stack + 83242;
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+0,int_stack+130822,int_stack+124072,150);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+96742,int_stack+140272,int_stack+21546,150);
     Libderiv->ABCD[7] = int_stack + 96742;
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+13500,int_stack+28296,int_stack+149722,150);
     Libderiv->ABCD[6] = int_stack + 13500;
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+110242,int_stack+163222,int_stack+156472, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 110242;
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+123742,int_stack+44496,int_stack+37746, 0.0, zero_stack, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 123742;
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+27000,int_stack+179422,int_stack+172672, 1.0, int_stack+62992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 27000;

}
