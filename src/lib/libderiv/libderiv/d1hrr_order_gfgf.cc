#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gfgf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gf|gf) integrals */

void d1hrr_order_gfgf(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[7][4][11] = int_stack + 6400;
 Libderiv->deriv_classes[7][5][11] = int_stack + 6940;
 Libderiv->deriv_classes[7][6][11] = int_stack + 7696;
 Libderiv->deriv_classes[7][7][11] = int_stack + 8704;
 Libderiv->deriv_classes[4][4][10] = int_stack + 10000;
 Libderiv->deriv_classes[4][5][10] = int_stack + 10225;
 Libderiv->deriv_classes[4][6][10] = int_stack + 10540;
 Libderiv->deriv_classes[4][7][10] = int_stack + 10960;
 Libderiv->deriv_classes[5][4][10] = int_stack + 11500;
 Libderiv->deriv_classes[5][5][10] = int_stack + 11815;
 Libderiv->deriv_classes[5][6][10] = int_stack + 12256;
 Libderiv->deriv_classes[5][7][10] = int_stack + 12844;
 Libderiv->deriv_classes[6][4][10] = int_stack + 13600;
 Libderiv->deriv_classes[6][5][10] = int_stack + 14020;
 Libderiv->deriv_classes[6][6][10] = int_stack + 14608;
 Libderiv->deriv_classes[6][7][10] = int_stack + 15392;
 Libderiv->deriv_classes[7][4][10] = int_stack + 16400;
 Libderiv->deriv_classes[7][5][10] = int_stack + 16940;
 Libderiv->deriv_classes[7][6][10] = int_stack + 17696;
 Libderiv->deriv_classes[7][7][10] = int_stack + 18704;
 Libderiv->deriv_classes[4][4][9] = int_stack + 20000;
 Libderiv->deriv_classes[4][5][9] = int_stack + 20225;
 Libderiv->deriv_classes[4][6][9] = int_stack + 20540;
 Libderiv->deriv_classes[4][7][9] = int_stack + 20960;
 Libderiv->deriv_classes[5][4][9] = int_stack + 21500;
 Libderiv->deriv_classes[5][5][9] = int_stack + 21815;
 Libderiv->deriv_classes[5][6][9] = int_stack + 22256;
 Libderiv->deriv_classes[5][7][9] = int_stack + 22844;
 Libderiv->deriv_classes[6][4][9] = int_stack + 23600;
 Libderiv->deriv_classes[6][5][9] = int_stack + 24020;
 Libderiv->deriv_classes[6][6][9] = int_stack + 24608;
 Libderiv->deriv_classes[6][7][9] = int_stack + 25392;
 Libderiv->deriv_classes[7][4][9] = int_stack + 26400;
 Libderiv->deriv_classes[7][5][9] = int_stack + 26940;
 Libderiv->deriv_classes[7][6][9] = int_stack + 27696;
 Libderiv->deriv_classes[7][7][9] = int_stack + 28704;
 Libderiv->deriv_classes[4][4][8] = int_stack + 30000;
 Libderiv->deriv_classes[4][5][8] = int_stack + 30225;
 Libderiv->deriv_classes[4][6][8] = int_stack + 30540;
 Libderiv->deriv_classes[4][7][8] = int_stack + 30960;
 Libderiv->deriv_classes[5][4][8] = int_stack + 31500;
 Libderiv->deriv_classes[5][5][8] = int_stack + 31815;
 Libderiv->deriv_classes[5][6][8] = int_stack + 32256;
 Libderiv->deriv_classes[5][7][8] = int_stack + 32844;
 Libderiv->deriv_classes[6][4][8] = int_stack + 33600;
 Libderiv->deriv_classes[6][5][8] = int_stack + 34020;
 Libderiv->deriv_classes[6][6][8] = int_stack + 34608;
 Libderiv->deriv_classes[6][7][8] = int_stack + 35392;
 Libderiv->deriv_classes[7][4][8] = int_stack + 36400;
 Libderiv->deriv_classes[7][5][8] = int_stack + 36940;
 Libderiv->deriv_classes[7][6][8] = int_stack + 37696;
 Libderiv->deriv_classes[7][7][8] = int_stack + 38704;
 Libderiv->deriv_classes[4][4][7] = int_stack + 40000;
 Libderiv->deriv_classes[4][5][7] = int_stack + 40225;
 Libderiv->deriv_classes[4][6][7] = int_stack + 40540;
 Libderiv->deriv_classes[4][7][7] = int_stack + 40960;
 Libderiv->deriv_classes[5][4][7] = int_stack + 41500;
 Libderiv->deriv_classes[5][5][7] = int_stack + 41815;
 Libderiv->deriv_classes[5][6][7] = int_stack + 42256;
 Libderiv->deriv_classes[5][7][7] = int_stack + 42844;
 Libderiv->deriv_classes[6][4][7] = int_stack + 43600;
 Libderiv->deriv_classes[6][5][7] = int_stack + 44020;
 Libderiv->deriv_classes[6][6][7] = int_stack + 44608;
 Libderiv->deriv_classes[6][7][7] = int_stack + 45392;
 Libderiv->deriv_classes[7][4][7] = int_stack + 46400;
 Libderiv->deriv_classes[7][5][7] = int_stack + 46940;
 Libderiv->deriv_classes[7][6][7] = int_stack + 47696;
 Libderiv->deriv_classes[7][7][7] = int_stack + 48704;
 Libderiv->deriv_classes[4][4][6] = int_stack + 50000;
 Libderiv->deriv_classes[4][5][6] = int_stack + 50225;
 Libderiv->deriv_classes[4][6][6] = int_stack + 50540;
 Libderiv->deriv_classes[4][7][6] = int_stack + 50960;
 Libderiv->deriv_classes[5][4][6] = int_stack + 51500;
 Libderiv->deriv_classes[5][5][6] = int_stack + 51815;
 Libderiv->deriv_classes[5][6][6] = int_stack + 52256;
 Libderiv->deriv_classes[5][7][6] = int_stack + 52844;
 Libderiv->deriv_classes[6][4][6] = int_stack + 53600;
 Libderiv->deriv_classes[6][5][6] = int_stack + 54020;
 Libderiv->deriv_classes[6][6][6] = int_stack + 54608;
 Libderiv->deriv_classes[6][7][6] = int_stack + 55392;
 Libderiv->dvrr_classes[7][4] = int_stack + 56400;
 Libderiv->deriv_classes[7][4][6] = int_stack + 56940;
 Libderiv->dvrr_classes[7][5] = int_stack + 57480;
 Libderiv->deriv_classes[7][5][6] = int_stack + 58236;
 Libderiv->dvrr_classes[7][6] = int_stack + 58992;
 Libderiv->deriv_classes[7][6][6] = int_stack + 60000;
 Libderiv->deriv_classes[7][7][6] = int_stack + 61008;
 Libderiv->deriv_classes[4][4][2] = int_stack + 62304;
 Libderiv->deriv_classes[4][5][2] = int_stack + 62529;
 Libderiv->deriv_classes[4][6][2] = int_stack + 62844;
 Libderiv->deriv_classes[4][7][2] = int_stack + 63264;
 Libderiv->deriv_classes[5][4][2] = int_stack + 63804;
 Libderiv->deriv_classes[5][5][2] = int_stack + 64119;
 Libderiv->deriv_classes[5][6][2] = int_stack + 64560;
 Libderiv->deriv_classes[5][7][2] = int_stack + 65148;
 Libderiv->deriv_classes[6][4][2] = int_stack + 65904;
 Libderiv->deriv_classes[6][5][2] = int_stack + 66324;
 Libderiv->deriv_classes[6][6][2] = int_stack + 66912;
 Libderiv->deriv_classes[6][7][2] = int_stack + 67696;
 Libderiv->deriv_classes[7][4][2] = int_stack + 68704;
 Libderiv->deriv_classes[7][5][2] = int_stack + 69244;
 Libderiv->deriv_classes[7][6][2] = int_stack + 70000;
 Libderiv->deriv_classes[7][7][2] = int_stack + 71008;
 Libderiv->deriv_classes[4][4][1] = int_stack + 72304;
 Libderiv->deriv_classes[4][5][1] = int_stack + 72529;
 Libderiv->deriv_classes[4][6][1] = int_stack + 72844;
 Libderiv->deriv_classes[4][7][1] = int_stack + 73264;
 Libderiv->deriv_classes[5][4][1] = int_stack + 73804;
 Libderiv->deriv_classes[5][5][1] = int_stack + 74119;
 Libderiv->deriv_classes[5][6][1] = int_stack + 74560;
 Libderiv->deriv_classes[5][7][1] = int_stack + 75148;
 Libderiv->deriv_classes[6][4][1] = int_stack + 75904;
 Libderiv->deriv_classes[6][5][1] = int_stack + 76324;
 Libderiv->deriv_classes[6][6][1] = int_stack + 76912;
 Libderiv->deriv_classes[6][7][1] = int_stack + 77696;
 Libderiv->deriv_classes[7][4][1] = int_stack + 78704;
 Libderiv->deriv_classes[7][5][1] = int_stack + 79244;
 Libderiv->deriv_classes[7][6][1] = int_stack + 80000;
 Libderiv->deriv_classes[7][7][1] = int_stack + 81008;
 Libderiv->dvrr_classes[4][4] = int_stack + 82304;
 Libderiv->dvrr_classes[4][5] = int_stack + 82529;
 Libderiv->dvrr_classes[4][6] = int_stack + 82844;
 Libderiv->dvrr_classes[4][7] = int_stack + 83264;
 Libderiv->deriv_classes[4][4][0] = int_stack + 83804;
 Libderiv->deriv_classes[4][5][0] = int_stack + 84029;
 Libderiv->deriv_classes[4][6][0] = int_stack + 84344;
 Libderiv->deriv_classes[4][7][0] = int_stack + 84764;
 Libderiv->dvrr_classes[5][4] = int_stack + 85304;
 Libderiv->dvrr_classes[5][5] = int_stack + 85619;
 Libderiv->dvrr_classes[5][6] = int_stack + 86060;
 Libderiv->dvrr_classes[5][7] = int_stack + 86648;
 Libderiv->deriv_classes[5][4][0] = int_stack + 87404;
 Libderiv->deriv_classes[5][5][0] = int_stack + 87719;
 Libderiv->deriv_classes[5][6][0] = int_stack + 88160;
 Libderiv->deriv_classes[5][7][0] = int_stack + 88748;
 Libderiv->dvrr_classes[6][4] = int_stack + 89504;
 Libderiv->dvrr_classes[6][5] = int_stack + 89924;
 Libderiv->dvrr_classes[6][6] = int_stack + 90512;
 Libderiv->dvrr_classes[6][7] = int_stack + 91296;
 Libderiv->deriv_classes[6][4][0] = int_stack + 92304;
 Libderiv->deriv_classes[6][5][0] = int_stack + 92724;
 Libderiv->deriv_classes[6][6][0] = int_stack + 93312;
 Libderiv->deriv_classes[6][7][0] = int_stack + 94096;
 Libderiv->deriv_classes[7][4][0] = int_stack + 95104;
 Libderiv->deriv_classes[7][5][0] = int_stack + 95644;
 Libderiv->deriv_classes[7][6][0] = int_stack + 96400;
 Libderiv->deriv_classes[7][7][0] = int_stack + 97408;
 memset(int_stack,0,789632);

 Libderiv->dvrr_stack = int_stack + 413938;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gfgf(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+98704,int_stack+82529,int_stack+82304,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+99379,int_stack+82844,int_stack+82529,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+100324,int_stack+99379,int_stack+98704,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101674,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82304,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+102349,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82529,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+103294,int_stack+102349,int_stack+101674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98704,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+104644,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82844,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+105904,int_stack+104644,int_stack+102349, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99379,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+107794,int_stack+105904,int_stack+103294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100324,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+101674,int_stack+85619,int_stack+85304,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+102619,int_stack+86060,int_stack+85619,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+103942,int_stack+102619,int_stack+101674,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+105832,int_stack+1815,int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85304,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+2256,int_stack+1815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85619,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+110044,int_stack+0,int_stack+105832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101674,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+105832,int_stack+2844,int_stack+2256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86060,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+111934,int_stack+105832,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102619,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+111934,int_stack+110044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103942,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+110044,int_stack+0,int_stack+107794,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+105832,int_stack+89924,int_stack+89504,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+107092,int_stack+90512,int_stack+89924,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+116794,int_stack+107092,int_stack+105832,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+119314,int_stack+4020,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89504,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+120574,int_stack+4608,int_stack+4020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89924,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+122338,int_stack+120574,int_stack+119314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105832,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+124858,int_stack+5392,int_stack+4608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90512,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+127210,int_stack+124858,int_stack+120574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107092,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+130738,int_stack+127210,int_stack+122338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116794,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+119314,int_stack+130738,int_stack+0,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+134938,int_stack+119314,int_stack+110044,150);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+57480,int_stack+56400,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1620,int_stack+58992,int_stack+57480,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+108856,int_stack+1620,int_stack+0,36);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3888,int_stack+6940,int_stack+6400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56400,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+112096,int_stack+7696,int_stack+6940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57480,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+148438,int_stack+112096,int_stack+3888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+3888,int_stack+8704,int_stack+7696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58992,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+151678,int_stack+3888,int_stack+112096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1620,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3888,int_stack+151678,int_stack+148438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108856,36);
 /*--- compute (ip|gf) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+148438,int_stack+3888,int_stack+130738,150);
 /*--- compute (hd|gf) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+161038,int_stack+148438,int_stack+119314,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+119314,int_stack+10225,int_stack+10000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82304, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+119989,int_stack+10540,int_stack+10225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82529, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+120934,int_stack+119989,int_stack+119314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98704, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+122284,int_stack+10960,int_stack+10540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82844, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+123544,int_stack+122284,int_stack+119989, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99379, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+125434,int_stack+123544,int_stack+120934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100324, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+119314,int_stack+11815,int_stack+11500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85304, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+120259,int_stack+12256,int_stack+11815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85619, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+121582,int_stack+120259,int_stack+119314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101674, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+123472,int_stack+12844,int_stack+12256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86060, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+127684,int_stack+123472,int_stack+120259, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102619, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+130330,int_stack+127684,int_stack+121582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103942, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+148438,int_stack+130330,int_stack+125434,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+119314,int_stack+14020,int_stack+13600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89504, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+120574,int_stack+14608,int_stack+14020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89924, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+122338,int_stack+120574,int_stack+119314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105832, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+124858,int_stack+15392,int_stack+14608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90512, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+155188,int_stack+124858,int_stack+120574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107092, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+124858,int_stack+155188,int_stack+122338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116794, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+3888,int_stack+124858,int_stack+130330,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+179938,int_stack+3888,int_stack+148438,150);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+16940,int_stack+16400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56400, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+150058,int_stack+17696,int_stack+16940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57480, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+152326,int_stack+150058,int_stack+148438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+155566,int_stack+18704,int_stack+17696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58992, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+129058,int_stack+155566,int_stack+150058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1620, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+155566,int_stack+129058,int_stack+152326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108856, 0.0, zero_stack,36);
 /*--- compute (ip|gf) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+193438,int_stack+155566,int_stack+124858,150);
 /*--- compute (hd|gf) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+206038,int_stack+193438,int_stack+3888,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3888,int_stack+20225,int_stack+20000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82304, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4563,int_stack+20540,int_stack+20225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82529, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5508,int_stack+4563,int_stack+3888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98704, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+6858,int_stack+20960,int_stack+20540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82844, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+8118,int_stack+6858,int_stack+4563, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99379, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+10008,int_stack+8118,int_stack+5508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100324, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3888,int_stack+21815,int_stack+21500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85304, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4833,int_stack+22256,int_stack+21815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85619, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6156,int_stack+4833,int_stack+3888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101674, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+8046,int_stack+22844,int_stack+22256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86060, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+12258,int_stack+8046,int_stack+4833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102619, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+14904,int_stack+12258,int_stack+6156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103942, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+193438,int_stack+14904,int_stack+10008,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3888,int_stack+24020,int_stack+23600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89504, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5148,int_stack+24608,int_stack+24020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89924, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6912,int_stack+5148,int_stack+3888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105832, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+9432,int_stack+25392,int_stack+24608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90512, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+18054,int_stack+9432,int_stack+5148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107092, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+9432,int_stack+18054,int_stack+6912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116794, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+148438,int_stack+9432,int_stack+14904,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+119314,int_stack+148438,int_stack+193438,150);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+193438,int_stack+26940,int_stack+26400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56400, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+195058,int_stack+27696,int_stack+26940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57480, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+197326,int_stack+195058,int_stack+193438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+200566,int_stack+28704,int_stack+27696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58992, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+13632,int_stack+200566,int_stack+195058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1620, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+200566,int_stack+13632,int_stack+197326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108856, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gf) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+13632,int_stack+200566,int_stack+9432,150);
 /*--- compute (hd|gf) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+224938,int_stack+13632,int_stack+148438,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+30225,int_stack+30000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149113,int_stack+30540,int_stack+30225, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+150058,int_stack+149113,int_stack+148438, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+151408,int_stack+30960,int_stack+30540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+152668,int_stack+151408,int_stack+149113, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+154558,int_stack+152668,int_stack+150058, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+31815,int_stack+31500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149383,int_stack+32256,int_stack+31815, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+150706,int_stack+149383,int_stack+148438, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+152596,int_stack+32844,int_stack+32256, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+156808,int_stack+152596,int_stack+149383, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+193438,int_stack+156808,int_stack+150706, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+196588,int_stack+193438,int_stack+154558,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+34020,int_stack+33600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149698,int_stack+34608,int_stack+34020, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89924, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+151462,int_stack+149698,int_stack+148438, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+153982,int_stack+35392,int_stack+34608, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+156334,int_stack+153982,int_stack+149698, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3888,int_stack+156334,int_stack+151462, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+148438,int_stack+3888,int_stack+193438,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+8088,int_stack+148438,int_stack+196588,150);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+193438,int_stack+36940,int_stack+36400, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+195058,int_stack+37696,int_stack+36940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+197326,int_stack+195058,int_stack+193438, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+200566,int_stack+38704,int_stack+37696, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21588,int_stack+200566,int_stack+195058, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+200566,int_stack+21588,int_stack+197326, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gf) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+21588,int_stack+200566,int_stack+3888,150);
 /*--- compute (hd|gf) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+243838,int_stack+21588,int_stack+148438,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+40225,int_stack+40000, 0.0, zero_stack, 1.0, int_stack+82304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149113,int_stack+40540,int_stack+40225, 0.0, zero_stack, 1.0, int_stack+82529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+150058,int_stack+149113,int_stack+148438, 0.0, zero_stack, 1.0, int_stack+98704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+151408,int_stack+40960,int_stack+40540, 0.0, zero_stack, 1.0, int_stack+82844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+152668,int_stack+151408,int_stack+149113, 0.0, zero_stack, 1.0, int_stack+99379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+154558,int_stack+152668,int_stack+150058, 0.0, zero_stack, 1.0, int_stack+100324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+41815,int_stack+41500, 0.0, zero_stack, 1.0, int_stack+85304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149383,int_stack+42256,int_stack+41815, 0.0, zero_stack, 1.0, int_stack+85619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+150706,int_stack+149383,int_stack+148438, 0.0, zero_stack, 1.0, int_stack+101674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+152596,int_stack+42844,int_stack+42256, 0.0, zero_stack, 1.0, int_stack+86060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+156808,int_stack+152596,int_stack+149383, 0.0, zero_stack, 1.0, int_stack+102619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+21588,int_stack+156808,int_stack+150706, 0.0, zero_stack, 1.0, int_stack+103942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+24738,int_stack+21588,int_stack+154558,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+44020,int_stack+43600, 0.0, zero_stack, 1.0, int_stack+89504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149698,int_stack+44608,int_stack+44020, 0.0, zero_stack, 1.0, int_stack+89924, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+151462,int_stack+149698,int_stack+148438, 0.0, zero_stack, 1.0, int_stack+105832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+153982,int_stack+45392,int_stack+44608, 0.0, zero_stack, 1.0, int_stack+90512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+156334,int_stack+153982,int_stack+149698, 0.0, zero_stack, 1.0, int_stack+107092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3888,int_stack+156334,int_stack+151462, 0.0, zero_stack, 1.0, int_stack+116794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+148438,int_stack+3888,int_stack+21588,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+31488,int_stack+148438,int_stack+24738,150);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21588,int_stack+46940,int_stack+46400, 0.0, zero_stack, 1.0, int_stack+56400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23208,int_stack+47696,int_stack+46940, 0.0, zero_stack, 1.0, int_stack+57480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+25476,int_stack+23208,int_stack+21588, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+157888,int_stack+48704,int_stack+47696, 0.0, zero_stack, 1.0, int_stack+58992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+44988,int_stack+157888,int_stack+23208, 0.0, zero_stack, 1.0, int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+193438,int_stack+44988,int_stack+25476, 0.0, zero_stack, 1.0, int_stack+108856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gf) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+262738,int_stack+193438,int_stack+3888,150);
 /*--- compute (hd|gf) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+275338,int_stack+262738,int_stack+148438,150);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+50225,int_stack+50000, 1.0, int_stack+82304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149113,int_stack+50540,int_stack+50225, 1.0, int_stack+82529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+150058,int_stack+149113,int_stack+148438, 1.0, int_stack+98704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+151408,int_stack+50960,int_stack+50540, 1.0, int_stack+82844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+152668,int_stack+151408,int_stack+149113, 1.0, int_stack+99379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+154558,int_stack+152668,int_stack+150058, 1.0, int_stack+100324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+51815,int_stack+51500, 1.0, int_stack+85304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149383,int_stack+52256,int_stack+51815, 1.0, int_stack+85619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+150706,int_stack+149383,int_stack+148438, 1.0, int_stack+101674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+152596,int_stack+52844,int_stack+52256, 1.0, int_stack+86060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+156808,int_stack+152596,int_stack+149383, 1.0, int_stack+102619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+262738,int_stack+156808,int_stack+150706, 1.0, int_stack+103942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+265888,int_stack+262738,int_stack+154558,150);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+54020,int_stack+53600, 1.0, int_stack+89504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149698,int_stack+54608,int_stack+54020, 1.0, int_stack+89924, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+151462,int_stack+149698,int_stack+148438, 1.0, int_stack+105832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+153982,int_stack+55392,int_stack+54608, 1.0, int_stack+90512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+156334,int_stack+153982,int_stack+149698, 1.0, int_stack+107092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+3888,int_stack+156334,int_stack+151462, 1.0, int_stack+116794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+148438,int_stack+3888,int_stack+262738,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+294238,int_stack+148438,int_stack+265888,150);
 /*--- compute (k0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+262738,int_stack+58236,int_stack+56940, 1.0, int_stack+56400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+264358,int_stack+60000,int_stack+58236, 1.0, int_stack+57480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+266626,int_stack+264358,int_stack+262738, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+269866,int_stack+61008,int_stack+60000, 1.0, int_stack+58992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+193438,int_stack+269866,int_stack+264358, 1.0, int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (k0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+269866,int_stack+193438,int_stack+266626, 1.0, int_stack+108856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (ip|gf) ---*/
   hrr1_build_ip(Libderiv->AB,int_stack+193438,int_stack+269866,int_stack+3888,150);
 /*--- compute (hd|gf) ---*/
   hrr1_build_hd(Libderiv->AB,int_stack+307738,int_stack+193438,int_stack+148438,150);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+105832,int_stack+83264,int_stack+82844,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+148438,int_stack+105832,int_stack+99379,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+150328,int_stack+148438,int_stack+100324,15);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+148438,int_stack+86648,int_stack+86060,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+152578,int_stack+148438,int_stack+102619,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+155224,int_stack+152578,int_stack+103942,21);
 /*--- compute (gp|gf) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+193438,int_stack+155224,int_stack+150328,150);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+152578,int_stack+91296,int_stack+90512,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+200188,int_stack+152578,int_stack+107092,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+200188,int_stack+116794,28);
 /*--- compute (hp|gf) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+262738,int_stack+0,int_stack+155224,150);
 /*--- compute (gd|gf) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+98704,int_stack+262738,int_stack+193438,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+200188,int_stack+62529,int_stack+62304,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+200863,int_stack+62844,int_stack+62529,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+201808,int_stack+200863,int_stack+200188,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+203158,int_stack+63264,int_stack+62844,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+148438,int_stack+203158,int_stack+200863,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+203158,int_stack+148438,int_stack+201808,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+64119,int_stack+63804,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+200188,int_stack+64560,int_stack+64119,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+152578,int_stack+200188,int_stack+148438,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+148438,int_stack+65148,int_stack+64560,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+158374,int_stack+148438,int_stack+200188,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+272188,int_stack+158374,int_stack+152578,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+112204,int_stack+272188,int_stack+203158, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+152578,int_stack+66324,int_stack+65904,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+158374,int_stack+66912,int_stack+66324,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+200188,int_stack+158374,int_stack+152578,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+152578,int_stack+67696,int_stack+66912,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+4200,int_stack+152578,int_stack+158374,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+44988,int_stack+4200,int_stack+200188,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+49188,int_stack+44988,int_stack+272188, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+326638,int_stack+49188,int_stack+112204, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+193438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+112204,int_stack+69244,int_stack+68704,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+113824,int_stack+70000,int_stack+69244,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+200188,int_stack+113824,int_stack+112204,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+116092,int_stack+71008,int_stack+70000,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+58638,int_stack+116092,int_stack+113824,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+112204,int_stack+58638,int_stack+200188,36);
 /*--- compute (ip|gf) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+58638,int_stack+112204,int_stack+44988, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (hd|gf) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+340138,int_stack+58638,int_stack+49188, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+262738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+44988,int_stack+72529,int_stack+72304,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+45663,int_stack+72844,int_stack+72529,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+46608,int_stack+45663,int_stack+44988,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+47958,int_stack+73264,int_stack+72844,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+148438,int_stack+47958,int_stack+45663,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+47958,int_stack+148438,int_stack+46608,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+74119,int_stack+73804,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+50208,int_stack+74560,int_stack+74119,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+51531,int_stack+50208,int_stack+148438,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+148438,int_stack+75148,int_stack+74560,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+152578,int_stack+148438,int_stack+50208,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+272188,int_stack+152578,int_stack+51531,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+50208,int_stack+272188,int_stack+47958, 0.0, zero_stack, 1.0, int_stack+150328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+152578,int_stack+76324,int_stack+75904,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+56958,int_stack+76912,int_stack+76324,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+58722,int_stack+56958,int_stack+152578,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+152578,int_stack+77696,int_stack+76912,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+61242,int_stack+152578,int_stack+56958,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+64770,int_stack+61242,int_stack+58722,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+68970,int_stack+64770,int_stack+272188, 0.0, zero_stack, 1.0, int_stack+155224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+359038,int_stack+68970,int_stack+50208, 0.0, zero_stack, 1.0, int_stack+193438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+272188,int_stack+79244,int_stack+78704,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+152578,int_stack+80000,int_stack+79244,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+44988,int_stack+152578,int_stack+272188,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+272188,int_stack+81008,int_stack+80000,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+78420,int_stack+272188,int_stack+152578,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+48228,int_stack+78420,int_stack+44988,36);
 /*--- compute (ip|gf) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+372538,int_stack+48228,int_stack+64770, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (hd|gf) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+44988,int_stack+372538,int_stack+68970, 0.0, zero_stack, 1.0, int_stack+262738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+372538,int_stack+84029,int_stack+83804,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+373213,int_stack+84344,int_stack+84029,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+374158,int_stack+373213,int_stack+372538,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+375508,int_stack+84764,int_stack+84344,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+148438,int_stack+375508,int_stack+373213,15);
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+375508,int_stack+148438,int_stack+374158,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+148438,int_stack+87719,int_stack+87404,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+377758,int_stack+88160,int_stack+87719,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+379081,int_stack+377758,int_stack+148438,21);
 /*--- compute (h0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+148438,int_stack+88748,int_stack+88160,21);
 /*--- compute (h0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+152578,int_stack+148438,int_stack+377758,21);
 /*--- compute (h0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+272188,int_stack+152578,int_stack+379081,21);
 /*--- compute (gp|gf) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+377758,int_stack+272188,int_stack+375508, 1.0, int_stack+150328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+384508,int_stack+92724,int_stack+92304,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+148438,int_stack+93312,int_stack+92724,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+150202,int_stack+148438,int_stack+384508,28);
 /*--- compute (i0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+152722,int_stack+94096,int_stack+93312,28);
 /*--- compute (i0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+372538,int_stack+152722,int_stack+148438,28);
 /*--- compute (i0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+63888,int_stack+372538,int_stack+150202,28);
 /*--- compute (hp|gf) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+68088,int_stack+63888,int_stack+272188, 1.0, int_stack+155224, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (gd|gf) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+77538,int_stack+68088,int_stack+377758, 1.0, int_stack+193438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (k0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+193438,int_stack+95644,int_stack+95104,36);
 /*--- compute (k0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+195058,int_stack+96400,int_stack+95644,36);
 /*--- compute (k0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+197326,int_stack+195058,int_stack+193438,36);
 /*--- compute (k0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+200566,int_stack+97408,int_stack+96400,36);
 /*--- compute (k0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+372538,int_stack+200566,int_stack+195058,36);
 /*--- compute (k0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+200566,int_stack+372538,int_stack+197326,36);
 /*--- compute (ip|gf) ---*/
   d1hrr1_build_ip(Libderiv->AB,int_stack+148438,int_stack+200566,int_stack+63888, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (hd|gf) ---*/
   d1hrr1_build_hd(Libderiv->AB,int_stack+372538,int_stack+148438,int_stack+68088, 1.0, int_stack+262738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
 /*--- compute (gf|gf) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+391438,int_stack+161038,int_stack+134938,150);
     Libderiv->ABCD[11] = int_stack + 391438;
 /*--- compute (gf|gf) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+132814,int_stack+206038,int_stack+179938,150);
     Libderiv->ABCD[10] = int_stack + 132814;
 /*--- compute (gf|gf) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+155314,int_stack+224938,int_stack+119314,150);
     Libderiv->ABCD[9] = int_stack + 155314;
 /*--- compute (gf|gf) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+177814,int_stack+243838,int_stack+8088,150);
     Libderiv->ABCD[8] = int_stack + 177814;
 /*--- compute (gf|gf) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+0,int_stack+275338,int_stack+31488,150);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (gf|gf) ---*/
   hrr1_build_gf(Libderiv->AB,int_stack+200314,int_stack+307738,int_stack+294238,150);
     Libderiv->ABCD[6] = int_stack + 200314;
 /*--- compute (gf|gf) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+222814,int_stack+340138,int_stack+326638, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[2] = int_stack + 222814;
 /*--- compute (gf|gf) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+245314,int_stack+44988,int_stack+359038, 0.0, zero_stack, 1.0, int_stack+98704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[1] = int_stack + 245314;
 /*--- compute (gf|gf) ---*/
   d1hrr1_build_gf(Libderiv->AB,int_stack+22500,int_stack+372538,int_stack+77538, 1.0, int_stack+98704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,150);
     Libderiv->ABCD[0] = int_stack + 22500;

}
