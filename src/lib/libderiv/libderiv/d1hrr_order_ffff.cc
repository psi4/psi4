#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ffff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (ff|ff) integrals */

void d1hrr_order_ffff(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[6][3][11] = int_stack + 3404;
 Libderiv->deriv_classes[6][4][11] = int_stack + 3684;
 Libderiv->deriv_classes[6][5][11] = int_stack + 4104;
 Libderiv->deriv_classes[6][6][11] = int_stack + 4692;
 Libderiv->deriv_classes[3][3][10] = int_stack + 5476;
 Libderiv->deriv_classes[3][4][10] = int_stack + 5576;
 Libderiv->deriv_classes[3][5][10] = int_stack + 5726;
 Libderiv->deriv_classes[3][6][10] = int_stack + 5936;
 Libderiv->deriv_classes[4][3][10] = int_stack + 6216;
 Libderiv->deriv_classes[4][4][10] = int_stack + 6366;
 Libderiv->deriv_classes[4][5][10] = int_stack + 6591;
 Libderiv->deriv_classes[4][6][10] = int_stack + 6906;
 Libderiv->deriv_classes[5][3][10] = int_stack + 7326;
 Libderiv->deriv_classes[5][4][10] = int_stack + 7536;
 Libderiv->deriv_classes[5][5][10] = int_stack + 7851;
 Libderiv->deriv_classes[5][6][10] = int_stack + 8292;
 Libderiv->deriv_classes[6][3][10] = int_stack + 8880;
 Libderiv->deriv_classes[6][4][10] = int_stack + 9160;
 Libderiv->deriv_classes[6][5][10] = int_stack + 9580;
 Libderiv->deriv_classes[6][6][10] = int_stack + 10168;
 Libderiv->deriv_classes[3][3][9] = int_stack + 10952;
 Libderiv->deriv_classes[3][4][9] = int_stack + 11052;
 Libderiv->deriv_classes[3][5][9] = int_stack + 11202;
 Libderiv->deriv_classes[3][6][9] = int_stack + 11412;
 Libderiv->deriv_classes[4][3][9] = int_stack + 11692;
 Libderiv->deriv_classes[4][4][9] = int_stack + 11842;
 Libderiv->deriv_classes[4][5][9] = int_stack + 12067;
 Libderiv->deriv_classes[4][6][9] = int_stack + 12382;
 Libderiv->deriv_classes[5][3][9] = int_stack + 12802;
 Libderiv->deriv_classes[5][4][9] = int_stack + 13012;
 Libderiv->deriv_classes[5][5][9] = int_stack + 13327;
 Libderiv->deriv_classes[5][6][9] = int_stack + 13768;
 Libderiv->deriv_classes[6][3][9] = int_stack + 14356;
 Libderiv->deriv_classes[6][4][9] = int_stack + 14636;
 Libderiv->deriv_classes[6][5][9] = int_stack + 15056;
 Libderiv->deriv_classes[6][6][9] = int_stack + 15644;
 Libderiv->deriv_classes[3][3][8] = int_stack + 16428;
 Libderiv->deriv_classes[3][4][8] = int_stack + 16528;
 Libderiv->deriv_classes[3][5][8] = int_stack + 16678;
 Libderiv->deriv_classes[3][6][8] = int_stack + 16888;
 Libderiv->deriv_classes[4][3][8] = int_stack + 17168;
 Libderiv->deriv_classes[4][4][8] = int_stack + 17318;
 Libderiv->deriv_classes[4][5][8] = int_stack + 17543;
 Libderiv->deriv_classes[4][6][8] = int_stack + 17858;
 Libderiv->deriv_classes[5][3][8] = int_stack + 18278;
 Libderiv->deriv_classes[5][4][8] = int_stack + 18488;
 Libderiv->deriv_classes[5][5][8] = int_stack + 18803;
 Libderiv->deriv_classes[5][6][8] = int_stack + 19244;
 Libderiv->deriv_classes[6][3][8] = int_stack + 19832;
 Libderiv->deriv_classes[6][4][8] = int_stack + 20112;
 Libderiv->deriv_classes[6][5][8] = int_stack + 20532;
 Libderiv->deriv_classes[6][6][8] = int_stack + 21120;
 Libderiv->deriv_classes[3][3][7] = int_stack + 21904;
 Libderiv->deriv_classes[3][4][7] = int_stack + 22004;
 Libderiv->deriv_classes[3][5][7] = int_stack + 22154;
 Libderiv->deriv_classes[3][6][7] = int_stack + 22364;
 Libderiv->deriv_classes[4][3][7] = int_stack + 22644;
 Libderiv->deriv_classes[4][4][7] = int_stack + 22794;
 Libderiv->deriv_classes[4][5][7] = int_stack + 23019;
 Libderiv->deriv_classes[4][6][7] = int_stack + 23334;
 Libderiv->deriv_classes[5][3][7] = int_stack + 23754;
 Libderiv->deriv_classes[5][4][7] = int_stack + 23964;
 Libderiv->deriv_classes[5][5][7] = int_stack + 24279;
 Libderiv->deriv_classes[5][6][7] = int_stack + 24720;
 Libderiv->deriv_classes[6][3][7] = int_stack + 25308;
 Libderiv->deriv_classes[6][4][7] = int_stack + 25588;
 Libderiv->deriv_classes[6][5][7] = int_stack + 26008;
 Libderiv->deriv_classes[6][6][7] = int_stack + 26596;
 Libderiv->deriv_classes[3][3][6] = int_stack + 27380;
 Libderiv->deriv_classes[3][4][6] = int_stack + 27480;
 Libderiv->deriv_classes[3][5][6] = int_stack + 27630;
 Libderiv->deriv_classes[3][6][6] = int_stack + 27840;
 Libderiv->deriv_classes[4][3][6] = int_stack + 28120;
 Libderiv->deriv_classes[4][4][6] = int_stack + 28270;
 Libderiv->deriv_classes[4][5][6] = int_stack + 28495;
 Libderiv->deriv_classes[4][6][6] = int_stack + 28810;
 Libderiv->deriv_classes[5][3][6] = int_stack + 29230;
 Libderiv->deriv_classes[5][4][6] = int_stack + 29440;
 Libderiv->deriv_classes[5][5][6] = int_stack + 29755;
 Libderiv->deriv_classes[5][6][6] = int_stack + 30196;
 Libderiv->dvrr_classes[6][3] = int_stack + 30784;
 Libderiv->deriv_classes[6][3][6] = int_stack + 31064;
 Libderiv->dvrr_classes[6][4] = int_stack + 31344;
 Libderiv->deriv_classes[6][4][6] = int_stack + 31764;
 Libderiv->dvrr_classes[6][5] = int_stack + 32184;
 Libderiv->deriv_classes[6][5][6] = int_stack + 32772;
 Libderiv->deriv_classes[6][6][6] = int_stack + 33360;
 Libderiv->deriv_classes[3][3][2] = int_stack + 34144;
 Libderiv->deriv_classes[3][4][2] = int_stack + 34244;
 Libderiv->deriv_classes[3][5][2] = int_stack + 34394;
 Libderiv->deriv_classes[3][6][2] = int_stack + 34604;
 Libderiv->deriv_classes[4][3][2] = int_stack + 34884;
 Libderiv->deriv_classes[4][4][2] = int_stack + 35034;
 Libderiv->deriv_classes[4][5][2] = int_stack + 35259;
 Libderiv->deriv_classes[4][6][2] = int_stack + 35574;
 Libderiv->deriv_classes[5][3][2] = int_stack + 35994;
 Libderiv->deriv_classes[5][4][2] = int_stack + 36204;
 Libderiv->deriv_classes[5][5][2] = int_stack + 36519;
 Libderiv->deriv_classes[5][6][2] = int_stack + 36960;
 Libderiv->deriv_classes[6][3][2] = int_stack + 37548;
 Libderiv->deriv_classes[6][4][2] = int_stack + 37828;
 Libderiv->deriv_classes[6][5][2] = int_stack + 38248;
 Libderiv->deriv_classes[6][6][2] = int_stack + 38836;
 Libderiv->deriv_classes[3][3][1] = int_stack + 39620;
 Libderiv->deriv_classes[3][4][1] = int_stack + 39720;
 Libderiv->deriv_classes[3][5][1] = int_stack + 39870;
 Libderiv->deriv_classes[3][6][1] = int_stack + 40080;
 Libderiv->deriv_classes[4][3][1] = int_stack + 40360;
 Libderiv->deriv_classes[4][4][1] = int_stack + 40510;
 Libderiv->deriv_classes[4][5][1] = int_stack + 40735;
 Libderiv->deriv_classes[4][6][1] = int_stack + 41050;
 Libderiv->deriv_classes[5][3][1] = int_stack + 41470;
 Libderiv->deriv_classes[5][4][1] = int_stack + 41680;
 Libderiv->deriv_classes[5][5][1] = int_stack + 41995;
 Libderiv->deriv_classes[5][6][1] = int_stack + 42436;
 Libderiv->deriv_classes[6][3][1] = int_stack + 43024;
 Libderiv->deriv_classes[6][4][1] = int_stack + 43304;
 Libderiv->deriv_classes[6][5][1] = int_stack + 43724;
 Libderiv->deriv_classes[6][6][1] = int_stack + 44312;
 Libderiv->dvrr_classes[3][3] = int_stack + 45096;
 Libderiv->dvrr_classes[3][4] = int_stack + 45196;
 Libderiv->dvrr_classes[3][5] = int_stack + 45346;
 Libderiv->dvrr_classes[3][6] = int_stack + 45556;
 Libderiv->deriv_classes[3][3][0] = int_stack + 45836;
 Libderiv->deriv_classes[3][4][0] = int_stack + 45936;
 Libderiv->deriv_classes[3][5][0] = int_stack + 46086;
 Libderiv->deriv_classes[3][6][0] = int_stack + 46296;
 Libderiv->dvrr_classes[4][3] = int_stack + 46576;
 Libderiv->dvrr_classes[4][4] = int_stack + 46726;
 Libderiv->dvrr_classes[4][5] = int_stack + 46951;
 Libderiv->dvrr_classes[4][6] = int_stack + 47266;
 Libderiv->deriv_classes[4][3][0] = int_stack + 47686;
 Libderiv->deriv_classes[4][4][0] = int_stack + 47836;
 Libderiv->deriv_classes[4][5][0] = int_stack + 48061;
 Libderiv->deriv_classes[4][6][0] = int_stack + 48376;
 Libderiv->dvrr_classes[5][3] = int_stack + 48796;
 Libderiv->dvrr_classes[5][4] = int_stack + 49006;
 Libderiv->dvrr_classes[5][5] = int_stack + 49321;
 Libderiv->dvrr_classes[5][6] = int_stack + 49762;
 Libderiv->deriv_classes[5][3][0] = int_stack + 50350;
 Libderiv->deriv_classes[5][4][0] = int_stack + 50560;
 Libderiv->deriv_classes[5][5][0] = int_stack + 50875;
 Libderiv->deriv_classes[5][6][0] = int_stack + 51316;
 Libderiv->deriv_classes[6][3][0] = int_stack + 51904;
 Libderiv->deriv_classes[6][4][0] = int_stack + 52184;
 Libderiv->deriv_classes[6][5][0] = int_stack + 52604;
 Libderiv->deriv_classes[6][6][0] = int_stack + 53192;
 memset(int_stack,0,431808);

 Libderiv->dvrr_stack = int_stack + 188931;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ffff(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53976,int_stack+45196,int_stack+45096,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+54276,int_stack+45346,int_stack+45196,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+54726,int_stack+54276,int_stack+53976,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+55326,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45096,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+55626,int_stack+250,int_stack+100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45196,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56076,int_stack+55626,int_stack+55326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53976,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+56676,int_stack+460,int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45346,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+57306,int_stack+56676,int_stack+55626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54276,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+58206,int_stack+57306,int_stack+56076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54726,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+55326,int_stack+46726,int_stack+46576,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+55776,int_stack+46951,int_stack+46726,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56451,int_stack+55776,int_stack+55326,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+57351,int_stack+890,int_stack+740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1115,int_stack+890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46726,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+59206,int_stack+0,int_stack+57351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55326,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60106,int_stack+1430,int_stack+1115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46951,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+61051,int_stack+60106,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55776,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+61051,int_stack+59206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56451,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+59206,int_stack+0,int_stack+58206,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+62206,int_stack+49006,int_stack+48796,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+57351,int_stack+49321,int_stack+49006,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+62836,int_stack+57351,int_stack+62206,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+64096,int_stack+2060,int_stack+1850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48796,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+64726,int_stack+2375,int_stack+2060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49006,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65671,int_stack+64726,int_stack+64096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62206,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66931,int_stack+2816,int_stack+2375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49321,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1500,int_stack+66931,int_stack+64726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57351,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+66931,int_stack+1500,int_stack+65671, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62836,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+69031,int_stack+66931,int_stack+0,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+73531,int_stack+69031,int_stack+59206,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+31344,int_stack+30784,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+840,int_stack+32184,int_stack+31344,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+64096,int_stack+840,int_stack+0,28);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+3684,int_stack+3404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30784,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+58296,int_stack+4104,int_stack+3684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31344,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+59556,int_stack+58296,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2100,int_stack+4692,int_stack+4104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32184,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+79531,int_stack+2100,int_stack+58296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2100,int_stack+79531,int_stack+59556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64096,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+79531,int_stack+2100,int_stack+66931,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+85831,int_stack+79531,int_stack+69031,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+5576,int_stack+5476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45096, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79831,int_stack+5726,int_stack+5576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45196, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80281,int_stack+79831,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53976, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80881,int_stack+5936,int_stack+5726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45346, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81511,int_stack+80881,int_stack+79831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54276, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+82411,int_stack+81511,int_stack+80281, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54726, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+6366,int_stack+6216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79981,int_stack+6591,int_stack+6366, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46726, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80656,int_stack+79981,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55326, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83411,int_stack+6906,int_stack+6591, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46951, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84356,int_stack+83411,int_stack+79981, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55776, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2100,int_stack+84356,int_stack+80656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56451, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3600,int_stack+2100,int_stack+82411,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+7536,int_stack+7326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48796, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80161,int_stack+7851,int_stack+7536, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49006, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81106,int_stack+80161,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62206, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+82366,int_stack+8292,int_stack+7851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49321, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83689,int_stack+82366,int_stack+80161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57351, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6600,int_stack+83689,int_stack+81106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62836, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79531,int_stack+6600,int_stack+2100,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+65776,int_stack+79531,int_stack+3600,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+9160,int_stack+8880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30784, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2940,int_stack+9580,int_stack+9160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31344, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4200,int_stack+2940,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84031,int_stack+10168,int_stack+9580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32184, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+58296,int_stack+84031,int_stack+2940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+94831,int_stack+58296,int_stack+4200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64096, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+97631,int_stack+94831,int_stack+6600,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+103931,int_stack+97631,int_stack+79531,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+11052,int_stack+10952, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45096, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79831,int_stack+11202,int_stack+11052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45196, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80281,int_stack+79831,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53976, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80881,int_stack+11412,int_stack+11202, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45346, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81511,int_stack+80881,int_stack+79831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54276, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+82411,int_stack+81511,int_stack+80281, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54726, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+11842,int_stack+11692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79981,int_stack+12067,int_stack+11842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46726, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80656,int_stack+79981,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55326, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83411,int_stack+12382,int_stack+12067, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46951, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84356,int_stack+83411,int_stack+79981, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55776, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+94831,int_stack+84356,int_stack+80656, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56451, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+96331,int_stack+94831,int_stack+82411,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+13012,int_stack+12802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48796, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80161,int_stack+13327,int_stack+13012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49006, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81106,int_stack+80161,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62206, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+82366,int_stack+13768,int_stack+13327, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49321, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83689,int_stack+82366,int_stack+80161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57351, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+99331,int_stack+83689,int_stack+81106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62836, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79531,int_stack+99331,int_stack+94831,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+2100,int_stack+79531,int_stack+96331,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+94831,int_stack+14636,int_stack+14356, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30784, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+95671,int_stack+15056,int_stack+14636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31344, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+96931,int_stack+95671,int_stack+94831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84031,int_stack+15644,int_stack+15056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32184, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+58296,int_stack+84031,int_stack+95671, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8100,int_stack+58296,int_stack+96931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64096, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+112931,int_stack+8100,int_stack+99331,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+94831,int_stack+112931,int_stack+79531,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+16528,int_stack+16428, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79831,int_stack+16678,int_stack+16528, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80281,int_stack+79831,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80881,int_stack+16888,int_stack+16678, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81511,int_stack+80881,int_stack+79831, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+82411,int_stack+81511,int_stack+80281, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+17318,int_stack+17168, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79981,int_stack+17543,int_stack+17318, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80656,int_stack+79981,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83411,int_stack+17858,int_stack+17543, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46951, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84356,int_stack+83411,int_stack+79981, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+112931,int_stack+84356,int_stack+80656, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+114431,int_stack+112931,int_stack+82411,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+18488,int_stack+18278, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80161,int_stack+18803,int_stack+18488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81106,int_stack+80161,int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+82366,int_stack+19244,int_stack+18803, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83689,int_stack+82366,int_stack+80161, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+117431,int_stack+83689,int_stack+81106, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79531,int_stack+117431,int_stack+112931,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+8100,int_stack+79531,int_stack+114431,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+112931,int_stack+20112,int_stack+19832, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113771,int_stack+20532,int_stack+20112, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+115031,int_stack+113771,int_stack+112931, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84031,int_stack+21120,int_stack+20532, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14100,int_stack+84031,int_stack+113771, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16620,int_stack+14100,int_stack+115031, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+119531,int_stack+16620,int_stack+117431,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+125831,int_stack+119531,int_stack+79531,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+22004,int_stack+21904, 0.0, zero_stack, 1.0, int_stack+45096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79831,int_stack+22154,int_stack+22004, 0.0, zero_stack, 1.0, int_stack+45196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80281,int_stack+79831,int_stack+79531, 0.0, zero_stack, 1.0, int_stack+53976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80881,int_stack+22364,int_stack+22154, 0.0, zero_stack, 1.0, int_stack+45346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81511,int_stack+80881,int_stack+79831, 0.0, zero_stack, 1.0, int_stack+54276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+82411,int_stack+81511,int_stack+80281, 0.0, zero_stack, 1.0, int_stack+54726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+22794,int_stack+22644, 0.0, zero_stack, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79981,int_stack+23019,int_stack+22794, 0.0, zero_stack, 1.0, int_stack+46726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80656,int_stack+79981,int_stack+79531, 0.0, zero_stack, 1.0, int_stack+55326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83411,int_stack+23334,int_stack+23019, 0.0, zero_stack, 1.0, int_stack+46951, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84356,int_stack+83411,int_stack+79981, 0.0, zero_stack, 1.0, int_stack+55776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14100,int_stack+84356,int_stack+80656, 0.0, zero_stack, 1.0, int_stack+56451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15600,int_stack+14100,int_stack+82411,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+23964,int_stack+23754, 0.0, zero_stack, 1.0, int_stack+48796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80161,int_stack+24279,int_stack+23964, 0.0, zero_stack, 1.0, int_stack+49006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81106,int_stack+80161,int_stack+79531, 0.0, zero_stack, 1.0, int_stack+62206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+82366,int_stack+24720,int_stack+24279, 0.0, zero_stack, 1.0, int_stack+49321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83689,int_stack+82366,int_stack+80161, 0.0, zero_stack, 1.0, int_stack+57351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18600,int_stack+83689,int_stack+81106, 0.0, zero_stack, 1.0, int_stack+62836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79531,int_stack+18600,int_stack+14100,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+112931,int_stack+79531,int_stack+15600,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14100,int_stack+25588,int_stack+25308, 0.0, zero_stack, 1.0, int_stack+30784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14940,int_stack+26008,int_stack+25588, 0.0, zero_stack, 1.0, int_stack+31344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16200,int_stack+14940,int_stack+14100, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84031,int_stack+26596,int_stack+26008, 0.0, zero_stack, 1.0, int_stack+32184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20700,int_stack+84031,int_stack+14940, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+23220,int_stack+20700,int_stack+16200, 0.0, zero_stack, 1.0, int_stack+64096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+118931,int_stack+23220,int_stack+18600,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+14100,int_stack+118931,int_stack+79531,100);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+27480,int_stack+27380, 1.0, int_stack+45096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79831,int_stack+27630,int_stack+27480, 1.0, int_stack+45196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80281,int_stack+79831,int_stack+79531, 1.0, int_stack+53976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80881,int_stack+27840,int_stack+27630, 1.0, int_stack+45346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81511,int_stack+80881,int_stack+79831, 1.0, int_stack+54276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+82411,int_stack+81511,int_stack+80281, 1.0, int_stack+54726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+28270,int_stack+28120, 1.0, int_stack+46576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79981,int_stack+28495,int_stack+28270, 1.0, int_stack+46726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80656,int_stack+79981,int_stack+79531, 1.0, int_stack+55326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83411,int_stack+28810,int_stack+28495, 1.0, int_stack+46951, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+84356,int_stack+83411,int_stack+79981, 1.0, int_stack+55776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+118931,int_stack+84356,int_stack+80656, 1.0, int_stack+56451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+120431,int_stack+118931,int_stack+82411,100);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+29440,int_stack+29230, 1.0, int_stack+48796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80161,int_stack+29755,int_stack+29440, 1.0, int_stack+49006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81106,int_stack+80161,int_stack+79531, 1.0, int_stack+62206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+82366,int_stack+30196,int_stack+29755, 1.0, int_stack+49321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83689,int_stack+82366,int_stack+80161, 1.0, int_stack+57351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+123431,int_stack+83689,int_stack+81106, 1.0, int_stack+62836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+79531,int_stack+123431,int_stack+118931,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+23100,int_stack+79531,int_stack+120431,100);
 /*--- compute (i0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+118931,int_stack+31764,int_stack+31064, 1.0, int_stack+30784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+119771,int_stack+32772,int_stack+31764, 1.0, int_stack+31344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121031,int_stack+119771,int_stack+118931, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84031,int_stack+33360,int_stack+32772, 1.0, int_stack+32184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29100,int_stack+84031,int_stack+119771, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (i0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+58296,int_stack+29100,int_stack+121031, 1.0, int_stack+64096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,28);
 /*--- compute (hp|ff) ---*/
   hrr1_build_hp(Libderiv->AB,int_stack+134831,int_stack+58296,int_stack+123431,100);
 /*--- compute (gd|ff) ---*/
   hrr1_build_gd(Libderiv->AB,int_stack+141131,int_stack+134831,int_stack+79531,100);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+79531,int_stack+45556,int_stack+45346,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+80161,int_stack+79531,int_stack+54276,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+81061,int_stack+80161,int_stack+54726,10);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+79531,int_stack+47266,int_stack+46951,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+82061,int_stack+79531,int_stack+55776,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+79531,int_stack+82061,int_stack+56451,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+82061,int_stack+79531,int_stack+81061,100);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+134831,int_stack+49762,int_stack+49321,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+136154,int_stack+134831,int_stack+57351,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+136154,int_stack+62836,21);
 /*--- compute (gp|ff) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+134831,int_stack+0,int_stack+79531,100);
 /*--- compute (fd|ff) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+118931,int_stack+134831,int_stack+82061,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+139331,int_stack+34244,int_stack+34144,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+139631,int_stack+34394,int_stack+34244,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+140081,int_stack+139631,int_stack+139331,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+85061,int_stack+34604,int_stack+34394,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+124931,int_stack+85061,int_stack+139631,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+29100,int_stack+124931,int_stack+140081,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+124931,int_stack+35034,int_stack+34884,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+85061,int_stack+35259,int_stack+35034,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+139331,int_stack+85061,int_stack+124931,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+30100,int_stack+35574,int_stack+35259,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+31045,int_stack+30100,int_stack+85061,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+32395,int_stack+31045,int_stack+139331,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+53976,int_stack+32395,int_stack+29100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+29100,int_stack+36204,int_stack+35994,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+29730,int_stack+36519,int_stack+36204,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+30675,int_stack+29730,int_stack+29100,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+139331,int_stack+36960,int_stack+36519,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+33895,int_stack+139331,int_stack+29730,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+56976,int_stack+33895,int_stack+30675,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+59076,int_stack+56976,int_stack+32395, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+29100,int_stack+59076,int_stack+53976, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53976,int_stack+37828,int_stack+37548,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+54816,int_stack+38248,int_stack+37828,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+139331,int_stack+54816,int_stack+53976,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+35100,int_stack+38836,int_stack+38248,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+36864,int_stack+35100,int_stack+54816,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+53976,int_stack+36864,int_stack+139331,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+150131,int_stack+53976,int_stack+56976, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+156431,int_stack+150131,int_stack+59076, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+134831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+150131,int_stack+39720,int_stack+39620,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+150431,int_stack+39870,int_stack+39720,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+150881,int_stack+150431,int_stack+150131,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+151481,int_stack+40080,int_stack+39870,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+124931,int_stack+151481,int_stack+150431,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+151481,int_stack+124931,int_stack+150881,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+124931,int_stack+40510,int_stack+40360,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+152481,int_stack+40735,int_stack+40510,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+153156,int_stack+152481,int_stack+124931,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+154056,int_stack+41050,int_stack+40735,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+150131,int_stack+154056,int_stack+152481,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+154056,int_stack+150131,int_stack+153156,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+53976,int_stack+154056,int_stack+151481, 0.0, zero_stack, 1.0, int_stack+81061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+150131,int_stack+41680,int_stack+41470,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+150761,int_stack+41995,int_stack+41680,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+151706,int_stack+150761,int_stack+150131,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+56976,int_stack+42436,int_stack+41995,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+58299,int_stack+56976,int_stack+150761,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+60189,int_stack+58299,int_stack+151706,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+35100,int_stack+60189,int_stack+154056, 0.0, zero_stack, 1.0, int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+150131,int_stack+35100,int_stack+53976, 0.0, zero_stack, 1.0, int_stack+82061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53976,int_stack+43304,int_stack+43024,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+54816,int_stack+43724,int_stack+43304,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56076,int_stack+54816,int_stack+53976,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+57756,int_stack+44312,int_stack+43724,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+62289,int_stack+57756,int_stack+54816,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+39600,int_stack+62289,int_stack+56076,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+165431,int_stack+39600,int_stack+60189, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+53976,int_stack+165431,int_stack+35100, 0.0, zero_stack, 1.0, int_stack+134831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+156131,int_stack+45936,int_stack+45836,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35100,int_stack+46086,int_stack+45936,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+35550,int_stack+35100,int_stack+156131,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+36150,int_stack+46296,int_stack+46086,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+124931,int_stack+36150,int_stack+35100,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+36150,int_stack+124931,int_stack+35550,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+124931,int_stack+47836,int_stack+47686,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35100,int_stack+48061,int_stack+47836,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+37150,int_stack+35100,int_stack+124931,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+38050,int_stack+48376,int_stack+48061,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+38995,int_stack+38050,int_stack+35100,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+40345,int_stack+38995,int_stack+37150,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+37150,int_stack+40345,int_stack+36150, 1.0, int_stack+81061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+35100,int_stack+50560,int_stack+50350,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+35730,int_stack+50875,int_stack+50560,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+41845,int_stack+35730,int_stack+35100,21);
 /*--- compute (h0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+43105,int_stack+51316,int_stack+50875,21);
 /*--- compute (h0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+44428,int_stack+43105,int_stack+35730,21);
 /*--- compute (h0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+46318,int_stack+44428,int_stack+41845,21);
 /*--- compute (gp|ff) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+165431,int_stack+46318,int_stack+40345, 1.0, int_stack+79531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (fd|ff) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+40150,int_stack+165431,int_stack+37150, 1.0, int_stack+82061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (i0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+79531,int_stack+52184,int_stack+51904,28);
 /*--- compute (i0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+80371,int_stack+52604,int_stack+52184,28);
 /*--- compute (i0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+81631,int_stack+80371,int_stack+79531,28);
 /*--- compute (i0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+83311,int_stack+53192,int_stack+52604,28);
 /*--- compute (i0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+35100,int_stack+83311,int_stack+80371,28);
 /*--- compute (i0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+62976,int_stack+35100,int_stack+81631,28);
 /*--- compute (hp|ff) ---*/
   d1hrr1_build_hp(Libderiv->AB,int_stack+79531,int_stack+62976,int_stack+46318, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (gd|ff) ---*/
   d1hrr1_build_gd(Libderiv->AB,int_stack+169931,int_stack+79531,int_stack+165431, 1.0, int_stack+134831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+178931,int_stack+85831,int_stack+73531,100);
     Libderiv->ABCD[11] = int_stack + 178931;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+71776,int_stack+103931,int_stack+65776,100);
     Libderiv->ABCD[10] = int_stack + 71776;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+81776,int_stack+94831,int_stack+2100,100);
     Libderiv->ABCD[9] = int_stack + 81776;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+91776,int_stack+125831,int_stack+8100,100);
     Libderiv->ABCD[8] = int_stack + 91776;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+0,int_stack+14100,int_stack+112931,100);
     Libderiv->ABCD[7] = int_stack + 0;
 /*--- compute (ff|ff) ---*/
   hrr1_build_ff(Libderiv->AB,int_stack+10000,int_stack+141131,int_stack+23100,100);
     Libderiv->ABCD[6] = int_stack + 10000;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+124931,int_stack+156431,int_stack+29100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 124931;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+20000,int_stack+53976,int_stack+150131, 0.0, zero_stack, 1.0, int_stack+118931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 20000;
 /*--- compute (ff|ff) ---*/
   d1hrr1_build_ff(Libderiv->AB,int_stack+30000,int_stack+169931,int_stack+40150, 1.0, int_stack+118931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 30000;

}
