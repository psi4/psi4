#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_dpfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|fd) integrals */

void d12hrr_order_dpfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][5][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][10] = int_stack + 210;
 Libderiv->deriv_classes[3][5][9] = int_stack + 420;
 Libderiv->deriv_classes[3][5][8] = int_stack + 630;
 Libderiv->deriv_classes[3][5][7] = int_stack + 840;
 Libderiv->dvrr_classes[3][4] = int_stack + 1050;
 Libderiv->deriv_classes[3][5][6] = int_stack + 1200;
 Libderiv->deriv_classes[3][5][2] = int_stack + 1410;
 Libderiv->deriv_classes[3][5][1] = int_stack + 1620;
 Libderiv->dvrr_classes[2][5] = int_stack + 1830;
 Libderiv->deriv_classes[3][5][0] = int_stack + 1956;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 2166;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 2226;
 Libderiv->deriv2_classes[2][5][143] = int_stack + 2316;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 2442;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 2542;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 2692;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 2902;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 2962;
 Libderiv->deriv2_classes[2][5][131] = int_stack + 3052;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 3178;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 3278;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 3428;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 3638;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 3698;
 Libderiv->deriv2_classes[2][5][130] = int_stack + 3788;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 3914;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 4014;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 4164;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 4374;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 4434;
 Libderiv->deriv2_classes[2][5][119] = int_stack + 4524;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 4650;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 4750;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 4900;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 5110;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 5170;
 Libderiv->deriv2_classes[2][5][118] = int_stack + 5260;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 5386;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 5486;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 5636;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 5846;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 5906;
 Libderiv->deriv2_classes[2][5][117] = int_stack + 5996;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 6122;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 6222;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 6372;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 6582;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 6642;
 Libderiv->deriv2_classes[2][5][107] = int_stack + 6732;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 6858;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 6958;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 7108;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 7318;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 7378;
 Libderiv->deriv2_classes[2][5][106] = int_stack + 7468;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 7594;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 7694;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 7844;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 8054;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 8114;
 Libderiv->deriv2_classes[2][5][105] = int_stack + 8204;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 8330;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 8430;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 8580;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 8790;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 8850;
 Libderiv->deriv2_classes[2][5][104] = int_stack + 8940;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 9066;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 9166;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 9316;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 9526;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 9586;
 Libderiv->deriv2_classes[2][5][95] = int_stack + 9676;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 9802;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 9902;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 10052;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 10262;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 10322;
 Libderiv->deriv2_classes[2][5][94] = int_stack + 10412;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 10538;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 10638;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 10788;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 10998;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 11058;
 Libderiv->deriv2_classes[2][5][93] = int_stack + 11148;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 11274;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 11374;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 11524;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 11734;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 11794;
 Libderiv->deriv2_classes[2][5][92] = int_stack + 11884;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 12010;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 12110;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 12260;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 12470;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 12530;
 Libderiv->deriv2_classes[2][5][91] = int_stack + 12620;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 12746;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 12846;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 12996;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 13206;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 13266;
 Libderiv->deriv2_classes[2][5][83] = int_stack + 13356;
 Libderiv->deriv_classes[3][3][11] = int_stack + 13482;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 13582;
 Libderiv->deriv_classes[3][4][11] = int_stack + 13682;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 13832;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 13982;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 14192;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 14252;
 Libderiv->deriv2_classes[2][5][82] = int_stack + 14342;
 Libderiv->deriv_classes[3][3][10] = int_stack + 14468;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 14568;
 Libderiv->deriv_classes[3][4][10] = int_stack + 14668;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 14818;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 14968;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 15178;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 15238;
 Libderiv->deriv2_classes[2][5][81] = int_stack + 15328;
 Libderiv->deriv_classes[3][3][9] = int_stack + 15454;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 15554;
 Libderiv->deriv_classes[3][4][9] = int_stack + 15654;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 15804;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 15954;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 16164;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 16224;
 Libderiv->deriv2_classes[2][5][80] = int_stack + 16314;
 Libderiv->deriv_classes[3][3][8] = int_stack + 16440;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 16540;
 Libderiv->deriv_classes[3][4][8] = int_stack + 16640;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 16790;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 16940;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 17150;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 17210;
 Libderiv->deriv2_classes[2][5][79] = int_stack + 17300;
 Libderiv->deriv_classes[3][3][7] = int_stack + 17426;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 17526;
 Libderiv->deriv_classes[3][4][7] = int_stack + 17626;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 17776;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 17926;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 18136;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 18196;
 Libderiv->deriv2_classes[2][5][78] = int_stack + 18286;
 Libderiv->dvrr_classes[3][3] = int_stack + 18412;
 Libderiv->deriv_classes[3][3][6] = int_stack + 18512;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 18612;
 Libderiv->deriv_classes[3][4][6] = int_stack + 18712;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 18862;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 19012;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 19222;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 19282;
 Libderiv->deriv2_classes[2][5][35] = int_stack + 19372;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 19498;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 19598;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 19748;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 19958;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 20018;
 Libderiv->deriv2_classes[2][5][34] = int_stack + 20108;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 20234;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 20334;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 20484;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 20694;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 20754;
 Libderiv->deriv2_classes[2][5][33] = int_stack + 20844;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 20970;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 21070;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 21220;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 21430;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 21490;
 Libderiv->deriv2_classes[2][5][32] = int_stack + 21580;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 21706;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 21806;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 21956;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 22166;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 22226;
 Libderiv->deriv2_classes[2][5][31] = int_stack + 22316;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 22442;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 22542;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 22692;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 22902;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 22962;
 Libderiv->deriv2_classes[2][5][30] = int_stack + 23052;
 Libderiv->deriv_classes[3][3][2] = int_stack + 23178;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 23278;
 Libderiv->deriv_classes[3][4][2] = int_stack + 23378;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 23528;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 23678;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 23888;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 23948;
 Libderiv->deriv2_classes[2][5][26] = int_stack + 24038;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 24164;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 24264;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 24414;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 24624;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 24684;
 Libderiv->deriv2_classes[2][5][23] = int_stack + 24774;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 24900;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 25000;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 25150;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 25360;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 25420;
 Libderiv->deriv2_classes[2][5][22] = int_stack + 25510;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 25636;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 25736;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 25886;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 26096;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 26156;
 Libderiv->deriv2_classes[2][5][21] = int_stack + 26246;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 26372;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 26472;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 26622;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 26832;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 26892;
 Libderiv->deriv2_classes[2][5][20] = int_stack + 26982;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 27108;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 27208;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 27358;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 27568;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 27628;
 Libderiv->deriv2_classes[2][5][19] = int_stack + 27718;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 27844;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 27944;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 28094;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 28304;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 28364;
 Libderiv->deriv2_classes[2][5][18] = int_stack + 28454;
 Libderiv->deriv_classes[3][3][1] = int_stack + 28580;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 28680;
 Libderiv->deriv_classes[3][4][1] = int_stack + 28780;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 28930;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 29080;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 29290;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 29350;
 Libderiv->deriv2_classes[2][5][14] = int_stack + 29440;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 29566;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 29666;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 29816;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 30026;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 30086;
 Libderiv->deriv2_classes[2][5][13] = int_stack + 30176;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 30302;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 30402;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 30552;
 Libderiv->deriv_classes[2][3][11] = int_stack + 30762;
 Libderiv->deriv_classes[2][4][11] = int_stack + 30822;
 Libderiv->deriv_classes[2][5][11] = int_stack + 30912;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 31038;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 31098;
 Libderiv->deriv2_classes[2][5][11] = int_stack + 31188;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 31314;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 31414;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 31564;
 Libderiv->deriv_classes[2][3][10] = int_stack + 31774;
 Libderiv->deriv_classes[2][4][10] = int_stack + 31834;
 Libderiv->deriv_classes[2][5][10] = int_stack + 31924;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 32050;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 32110;
 Libderiv->deriv2_classes[2][5][10] = int_stack + 32200;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 32326;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 32426;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 32576;
 Libderiv->deriv_classes[2][3][9] = int_stack + 32786;
 Libderiv->deriv_classes[2][4][9] = int_stack + 32846;
 Libderiv->deriv_classes[2][5][9] = int_stack + 32936;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 33062;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 33122;
 Libderiv->deriv2_classes[2][5][9] = int_stack + 33212;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 33338;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 33438;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 33588;
 Libderiv->deriv_classes[2][3][8] = int_stack + 33798;
 Libderiv->deriv_classes[2][4][8] = int_stack + 33858;
 Libderiv->deriv_classes[2][5][8] = int_stack + 33948;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 34074;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 34134;
 Libderiv->deriv2_classes[2][5][8] = int_stack + 34224;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 34350;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 34450;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 34600;
 Libderiv->deriv_classes[2][3][7] = int_stack + 34810;
 Libderiv->deriv_classes[2][4][7] = int_stack + 34870;
 Libderiv->deriv_classes[2][5][7] = int_stack + 34960;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 35086;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 35146;
 Libderiv->deriv2_classes[2][5][7] = int_stack + 35236;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 35362;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 35462;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 35612;
 Libderiv->dvrr_classes[2][3] = int_stack + 35822;
 Libderiv->deriv_classes[2][3][6] = int_stack + 35882;
 Libderiv->dvrr_classes[2][4] = int_stack + 35942;
 Libderiv->deriv_classes[2][4][6] = int_stack + 36032;
 Libderiv->deriv_classes[2][5][6] = int_stack + 36122;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 36248;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 36308;
 Libderiv->deriv2_classes[2][5][6] = int_stack + 36398;
 Libderiv->deriv_classes[3][3][0] = int_stack + 36524;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 36624;
 Libderiv->deriv_classes[3][4][0] = int_stack + 36724;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 36874;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 37024;
 Libderiv->deriv_classes[2][3][2] = int_stack + 37234;
 Libderiv->deriv_classes[2][4][2] = int_stack + 37294;
 Libderiv->deriv_classes[2][5][2] = int_stack + 37384;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 37510;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 37570;
 Libderiv->deriv2_classes[2][5][2] = int_stack + 37660;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 37786;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 37886;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 38036;
 Libderiv->deriv_classes[2][3][1] = int_stack + 38246;
 Libderiv->deriv_classes[2][4][1] = int_stack + 38306;
 Libderiv->deriv_classes[2][5][1] = int_stack + 38396;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 38522;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 38582;
 Libderiv->deriv2_classes[2][5][1] = int_stack + 38672;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 38798;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 38898;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 39048;
 Libderiv->deriv_classes[2][3][0] = int_stack + 39258;
 Libderiv->deriv_classes[2][4][0] = int_stack + 39318;
 Libderiv->deriv_classes[2][5][0] = int_stack + 39408;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 39534;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 39594;
 Libderiv->deriv2_classes[2][5][0] = int_stack + 39684;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 39810;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 39910;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 40060;
 memset(int_stack,0,322160);

 Libderiv->dvrr_stack = int_stack + 70390;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_dpfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+35942,int_stack+35822,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40450,int_stack+30822,int_stack+30762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35822,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+40630,int_stack+30912,int_stack+30822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35942,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+40900,int_stack+40630,int_stack+40450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+1050,int_stack+18412,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41560,int_stack+13682,int_stack+13482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18412,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41860,int_stack+0,int_stack+13682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42310,int_stack+41860,int_stack+41560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41260,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41860,int_stack+31834,int_stack+31774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35822, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42040,int_stack+31924,int_stack+31834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35942, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42910,int_stack+42040,int_stack+41860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40270, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43270,int_stack+14668,int_stack+14468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18412, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43570,int_stack+210,int_stack+14668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+44020,int_stack+43570,int_stack+43270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41260, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43570,int_stack+32846,int_stack+32786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35822, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43750,int_stack+32936,int_stack+32846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35942, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+43750,int_stack+43570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40270, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44620,int_stack+15654,int_stack+15454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18412, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+44920,int_stack+420,int_stack+15654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45370,int_stack+44920,int_stack+44620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41260, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44920,int_stack+33858,int_stack+33798, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45100,int_stack+33948,int_stack+33858, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45970,int_stack+45100,int_stack+44920, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+46330,int_stack+16640,int_stack+16440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18412, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+46630,int_stack+630,int_stack+16640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+47080,int_stack+46630,int_stack+46330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+46630,int_stack+34870,int_stack+34810, 0.0, zero_stack, 1.0, int_stack+35822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+46810,int_stack+34960,int_stack+34870, 0.0, zero_stack, 1.0, int_stack+35942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+360,int_stack+46810,int_stack+46630, 0.0, zero_stack, 1.0, int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47680,int_stack+17626,int_stack+17426, 0.0, zero_stack, 1.0, int_stack+18412, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47980,int_stack+840,int_stack+17626, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+48430,int_stack+47980,int_stack+47680, 0.0, zero_stack, 1.0, int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47980,int_stack+36032,int_stack+35882, 1.0, int_stack+35822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48160,int_stack+36122,int_stack+36032, 1.0, int_stack+35942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+49030,int_stack+48160,int_stack+47980, 1.0, int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+18712,int_stack+18512, 1.0, int_stack+18412, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49390,int_stack+1200,int_stack+18712, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+49840,int_stack+49390,int_stack+720, 1.0, int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+48160,int_stack+1830,int_stack+35942,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+49390,int_stack+48160,int_stack+40270,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+37294,int_stack+37234,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+48160,int_stack+37384,int_stack+37294,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1020,int_stack+48160,int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+23378,int_stack+23178,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+50440,int_stack+1410,int_stack+23378,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+50890,int_stack+50440,int_stack+41260,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+50440,int_stack+38306,int_stack+38246,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+50620,int_stack+38396,int_stack+38306,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+51490,int_stack+50620,int_stack+50440,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51850,int_stack+28780,int_stack+28580,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52150,int_stack+1620,int_stack+28780,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+52600,int_stack+52150,int_stack+51850,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52150,int_stack+39318,int_stack+39258,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+39408,int_stack+39318,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1380,int_stack+52330,int_stack+52150,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53200,int_stack+36724,int_stack+36524,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+53500,int_stack+1956,int_stack+36724,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+53950,int_stack+53500,int_stack+53200,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+2226,int_stack+2166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+30762,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+2316,int_stack+2226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+30822,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1740,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40450,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+2542,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13482,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+54550,int_stack+2692,int_stack+2542, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13682,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+54550,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41560,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+2962,int_stack+2902, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30762, 1.0, int_stack+31774,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+3052,int_stack+2962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30822, 1.0, int_stack+31834,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+54550,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40450, 1.0, int_stack+41860,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+3278,int_stack+3178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13482, 1.0, int_stack+14468,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2700,int_stack+3428,int_stack+3278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13682, 1.0, int_stack+14668,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+54910,int_stack+2700,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41560, 1.0, int_stack+43270,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+3698,int_stack+3638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+31774, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+3788,int_stack+3698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+31834, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2700,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41860, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+4014,int_stack+3914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14468, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3060,int_stack+4164,int_stack+4014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14668, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3510,int_stack+3060,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43270, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+4434,int_stack+4374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30762, 0.0, zero_stack, 1.0, int_stack+32786,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+4524,int_stack+4434, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30822, 0.0, zero_stack, 1.0, int_stack+32846,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3060,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40450, 0.0, zero_stack, 1.0, int_stack+43570,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+4750,int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13482, 0.0, zero_stack, 1.0, int_stack+15454,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4110,int_stack+4900,int_stack+4750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13682, 0.0, zero_stack, 1.0, int_stack+15654,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+55510,int_stack+4110,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41560, 0.0, zero_stack, 1.0, int_stack+44620,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+5170,int_stack+5110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31774, 1.0, int_stack+32786, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+5260,int_stack+5170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31834, 1.0, int_stack+32846, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4110,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41860, 1.0, int_stack+43570, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+5486,int_stack+5386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14468, 1.0, int_stack+15454, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4470,int_stack+5636,int_stack+5486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14668, 1.0, int_stack+15654, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4920,int_stack+4470,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43270, 1.0, int_stack+44620, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+5906,int_stack+5846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+32786, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+5996,int_stack+5906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+32846, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4470,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43570, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+6222,int_stack+6122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15454, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5520,int_stack+6372,int_stack+6222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15654, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5970,int_stack+5520,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+44620, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+6642,int_stack+6582, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33798,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+6732,int_stack+6642, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30822, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33858,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5520,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44920,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+6958,int_stack+6858, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13482, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16440,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56110,int_stack+7108,int_stack+6958, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13682, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16640,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6570,int_stack+56110,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46330,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+7378,int_stack+7318, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31774, 0.0, zero_stack, 1.0, int_stack+33798, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+7468,int_stack+7378, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31834, 0.0, zero_stack, 1.0, int_stack+33858, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56110,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41860, 0.0, zero_stack, 1.0, int_stack+44920, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+7694,int_stack+7594, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14468, 0.0, zero_stack, 1.0, int_stack+16440, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56470,int_stack+7844,int_stack+7694, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14668, 0.0, zero_stack, 1.0, int_stack+16640, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7170,int_stack+56470,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43270, 0.0, zero_stack, 1.0, int_stack+46330, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+8114,int_stack+8054, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32786, 1.0, int_stack+33798, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+8204,int_stack+8114, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32846, 1.0, int_stack+33858, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56470,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43570, 1.0, int_stack+44920, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+8430,int_stack+8330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15454, 1.0, int_stack+16440, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7770,int_stack+8580,int_stack+8430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15654, 1.0, int_stack+16640, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56830,int_stack+7770,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44620, 1.0, int_stack+46330, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+8850,int_stack+8790, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+33798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+8940,int_stack+8850, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+33858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7770,int_stack+53680,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+44920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+9166,int_stack+9066, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8130,int_stack+9316,int_stack+9166, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8580,int_stack+8130,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+46330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+9586,int_stack+9526, 0.0, zero_stack, 1.0, int_stack+30762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34810,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+9676,int_stack+9586, 0.0, zero_stack, 1.0, int_stack+30822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34870,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8130,int_stack+53680,int_stack+53500, 0.0, zero_stack, 1.0, int_stack+40450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46630,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+9902,int_stack+9802, 0.0, zero_stack, 1.0, int_stack+13482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17426,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9180,int_stack+10052,int_stack+9902, 0.0, zero_stack, 1.0, int_stack+13682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17626,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9630,int_stack+9180,int_stack+53500, 0.0, zero_stack, 1.0, int_stack+41560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47680,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+10322,int_stack+10262, 0.0, zero_stack, 1.0, int_stack+31774, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34810, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+10412,int_stack+10322, 0.0, zero_stack, 1.0, int_stack+31834, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34870, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9180,int_stack+53680,int_stack+53500, 0.0, zero_stack, 1.0, int_stack+41860, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46630, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+10638,int_stack+10538, 0.0, zero_stack, 1.0, int_stack+14468, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17426, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+57430,int_stack+10788,int_stack+10638, 0.0, zero_stack, 1.0, int_stack+14668, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17626, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10230,int_stack+57430,int_stack+53500, 0.0, zero_stack, 1.0, int_stack+43270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47680, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+11058,int_stack+10998, 0.0, zero_stack, 1.0, int_stack+32786, 0.0, zero_stack, 1.0, int_stack+34810, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+11148,int_stack+11058, 0.0, zero_stack, 1.0, int_stack+32846, 0.0, zero_stack, 1.0, int_stack+34870, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57430,int_stack+53680,int_stack+53500, 0.0, zero_stack, 1.0, int_stack+43570, 0.0, zero_stack, 1.0, int_stack+46630, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+11374,int_stack+11274, 0.0, zero_stack, 1.0, int_stack+15454, 0.0, zero_stack, 1.0, int_stack+17426, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10830,int_stack+11524,int_stack+11374, 0.0, zero_stack, 1.0, int_stack+15654, 0.0, zero_stack, 1.0, int_stack+17626, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57790,int_stack+10830,int_stack+53500, 0.0, zero_stack, 1.0, int_stack+44620, 0.0, zero_stack, 1.0, int_stack+47680, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+11794,int_stack+11734, 0.0, zero_stack, 1.0, int_stack+33798, 1.0, int_stack+34810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+11884,int_stack+11794, 0.0, zero_stack, 1.0, int_stack+33858, 1.0, int_stack+34870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10830,int_stack+53680,int_stack+53500, 0.0, zero_stack, 1.0, int_stack+44920, 1.0, int_stack+46630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+12110,int_stack+12010, 0.0, zero_stack, 1.0, int_stack+16440, 1.0, int_stack+17426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11190,int_stack+12260,int_stack+12110, 0.0, zero_stack, 1.0, int_stack+16640, 1.0, int_stack+17626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11640,int_stack+11190,int_stack+53500, 0.0, zero_stack, 1.0, int_stack+46330, 1.0, int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+12530,int_stack+12470, 0.0, zero_stack, 2.0, int_stack+34810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+12620,int_stack+12530, 0.0, zero_stack, 2.0, int_stack+34870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11190,int_stack+53680,int_stack+53500, 0.0, zero_stack, 2.0, int_stack+46630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+12846,int_stack+12746, 0.0, zero_stack, 2.0, int_stack+17426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12240,int_stack+12996,int_stack+12846, 0.0, zero_stack, 2.0, int_stack+17626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+58390,int_stack+12240,int_stack+53500, 0.0, zero_stack, 2.0, int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53500,int_stack+13266,int_stack+13206, 1.0, int_stack+30762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35882,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53680,int_stack+13356,int_stack+13266, 1.0, int_stack+30822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36032,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12240,int_stack+53680,int_stack+53500, 1.0, int_stack+40450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47980,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40450,int_stack+13832,int_stack+13582, 1.0, int_stack+13482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18512,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53500,int_stack+13982,int_stack+13832, 1.0, int_stack+13682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18712,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+53500,int_stack+40450, 1.0, int_stack+41560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41560,int_stack+14252,int_stack+14192, 1.0, int_stack+31774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35882, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+14342,int_stack+14252, 1.0, int_stack+31834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36032, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+40450,int_stack+52330,int_stack+41560, 1.0, int_stack+41860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47980, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41560,int_stack+14818,int_stack+14568, 1.0, int_stack+14468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18512, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41860,int_stack+14968,int_stack+14818, 1.0, int_stack+14668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18712, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13200,int_stack+41860,int_stack+41560, 1.0, int_stack+43270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43270,int_stack+15238,int_stack+15178, 1.0, int_stack+32786, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35882, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+15328,int_stack+15238, 1.0, int_stack+32846, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36032, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41560,int_stack+52330,int_stack+43270, 1.0, int_stack+43570, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47980, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43270,int_stack+15804,int_stack+15554, 1.0, int_stack+15454, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18512, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43570,int_stack+15954,int_stack+15804, 1.0, int_stack+15654, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18712, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13800,int_stack+43570,int_stack+43270, 1.0, int_stack+44620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44620,int_stack+16224,int_stack+16164, 1.0, int_stack+33798, 0.0, zero_stack, 1.0, int_stack+35882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+16314,int_stack+16224, 1.0, int_stack+33858, 0.0, zero_stack, 1.0, int_stack+36032, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+43270,int_stack+52330,int_stack+44620, 1.0, int_stack+44920, 0.0, zero_stack, 1.0, int_stack+47980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44620,int_stack+16790,int_stack+16540, 1.0, int_stack+16440, 0.0, zero_stack, 1.0, int_stack+18512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+44920,int_stack+16940,int_stack+16790, 1.0, int_stack+16640, 0.0, zero_stack, 1.0, int_stack+18712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14400,int_stack+44920,int_stack+44620, 1.0, int_stack+46330, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+46330,int_stack+17210,int_stack+17150, 1.0, int_stack+34810, 1.0, int_stack+35882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+17300,int_stack+17210, 1.0, int_stack+34870, 1.0, int_stack+36032, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+44620,int_stack+52330,int_stack+46330, 1.0, int_stack+46630, 1.0, int_stack+47980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+46330,int_stack+17776,int_stack+17526, 1.0, int_stack+17426, 1.0, int_stack+18512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+46630,int_stack+17926,int_stack+17776, 1.0, int_stack+17626, 1.0, int_stack+18712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15000,int_stack+46630,int_stack+46330, 1.0, int_stack+47680, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47680,int_stack+18196,int_stack+18136, 2.0, int_stack+35882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+18286,int_stack+18196, 2.0, int_stack+36032, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+46330,int_stack+52330,int_stack+47680, 2.0, int_stack+47980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47680,int_stack+18862,int_stack+18612, 2.0, int_stack+18512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47980,int_stack+19012,int_stack+18862, 2.0, int_stack+18712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15600,int_stack+47980,int_stack+47680, 2.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+19282,int_stack+19222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37234,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+19372,int_stack+19282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37294,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+47680,int_stack+52330,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+19598,int_stack+19498, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23178,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53500,int_stack+19748,int_stack+19598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23378,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16200,int_stack+53500,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41260,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+20018,int_stack+19958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37234, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+20108,int_stack+20018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37294, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+53500,int_stack+52330,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40270, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+20334,int_stack+20234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23178, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16800,int_stack+20484,int_stack+20334, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23378, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17250,int_stack+16800,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41260, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+20754,int_stack+20694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37234, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+20844,int_stack+20754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37294, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16800,int_stack+52330,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40270, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+21070,int_stack+20970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23178, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17850,int_stack+21220,int_stack+21070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23378, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18300,int_stack+17850,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41260, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+21490,int_stack+21430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+21580,int_stack+21490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17850,int_stack+52330,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+21806,int_stack+21706, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18900,int_stack+21956,int_stack+21806, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19350,int_stack+18900,int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+22226,int_stack+22166, 0.0, zero_stack, 1.0, int_stack+37234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+22316,int_stack+22226, 0.0, zero_stack, 1.0, int_stack+37294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18900,int_stack+52330,int_stack+720, 0.0, zero_stack, 1.0, int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+22542,int_stack+22442, 0.0, zero_stack, 1.0, int_stack+23178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19950,int_stack+22692,int_stack+22542, 0.0, zero_stack, 1.0, int_stack+23378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20400,int_stack+19950,int_stack+720, 0.0, zero_stack, 1.0, int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+22962,int_stack+22902, 1.0, int_stack+37234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+23052,int_stack+22962, 1.0, int_stack+37294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19950,int_stack+52330,int_stack+720, 1.0, int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+23528,int_stack+23278, 1.0, int_stack+23178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21000,int_stack+23678,int_stack+23528, 1.0, int_stack+23378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21450,int_stack+21000,int_stack+720, 1.0, int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+23948,int_stack+23888,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+24038,int_stack+23948,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+21000,int_stack+52330,int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+24264,int_stack+24164,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22050,int_stack+24414,int_stack+24264,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+22500,int_stack+22050,int_stack+41260,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+24684,int_stack+24624, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38246,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+24774,int_stack+24684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38306,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22050,int_stack+52330,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50440,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+25000,int_stack+24900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28580,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23100,int_stack+25150,int_stack+25000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28780,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23550,int_stack+23100,int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51850,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+25420,int_stack+25360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38246, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+25510,int_stack+25420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38306, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23100,int_stack+52330,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50440, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+25736,int_stack+25636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28580, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24150,int_stack+25886,int_stack+25736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28780, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24600,int_stack+24150,int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51850, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+26156,int_stack+26096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38246, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+26246,int_stack+26156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38306, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24150,int_stack+52330,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50440, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+26472,int_stack+26372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28580, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25200,int_stack+26622,int_stack+26472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28780, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25650,int_stack+25200,int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51850, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+26892,int_stack+26832, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+26982,int_stack+26892, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38306, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25200,int_stack+52330,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+27208,int_stack+27108, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26250,int_stack+27358,int_stack+27208, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26700,int_stack+26250,int_stack+41260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+27628,int_stack+27568, 0.0, zero_stack, 1.0, int_stack+38246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+27718,int_stack+27628, 0.0, zero_stack, 1.0, int_stack+38306, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26250,int_stack+52330,int_stack+40270, 0.0, zero_stack, 1.0, int_stack+50440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+27944,int_stack+27844, 0.0, zero_stack, 1.0, int_stack+28580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+27300,int_stack+28094,int_stack+27944, 0.0, zero_stack, 1.0, int_stack+28780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+58990,int_stack+27300,int_stack+41260, 0.0, zero_stack, 1.0, int_stack+51850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+28364,int_stack+28304, 1.0, int_stack+38246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+28454,int_stack+28364, 1.0, int_stack+38306, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+27300,int_stack+52330,int_stack+40270, 1.0, int_stack+50440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+28930,int_stack+28680, 1.0, int_stack+28580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+50440,int_stack+29080,int_stack+28930, 1.0, int_stack+28780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+27660,int_stack+50440,int_stack+41260, 1.0, int_stack+51850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+29350,int_stack+29290,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+29440,int_stack+29350,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+50440,int_stack+52330,int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51850,int_stack+29666,int_stack+29566,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28260,int_stack+29816,int_stack+29666,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+28710,int_stack+28260,int_stack+51850,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+30086,int_stack+30026,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+30176,int_stack+30086,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+28260,int_stack+52330,int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51850,int_stack+30402,int_stack+30302,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+29310,int_stack+30552,int_stack+30402,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+29760,int_stack+29310,int_stack+51850,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+31098,int_stack+31038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39258,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+31188,int_stack+31098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39318,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29310,int_stack+52330,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52150,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51850,int_stack+31414,int_stack+31314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36524,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+30360,int_stack+31564,int_stack+31414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36724,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30810,int_stack+30360,int_stack+51850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53200,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+32110,int_stack+32050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39258, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+32200,int_stack+32110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39318, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30360,int_stack+52330,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52150, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51850,int_stack+32426,int_stack+32326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36524, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31410,int_stack+32576,int_stack+32426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36724, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31860,int_stack+31410,int_stack+51850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53200, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+33122,int_stack+33062, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39258, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+33212,int_stack+33122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39318, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31410,int_stack+52330,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52150, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51850,int_stack+33438,int_stack+33338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36524, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32460,int_stack+33588,int_stack+33438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36724, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32910,int_stack+32460,int_stack+51850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53200, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+34134,int_stack+34074, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+34224,int_stack+34134, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32460,int_stack+52330,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51850,int_stack+34450,int_stack+34350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33510,int_stack+34600,int_stack+34450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36724, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33960,int_stack+33510,int_stack+51850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+35146,int_stack+35086, 0.0, zero_stack, 1.0, int_stack+39258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+35236,int_stack+35146, 0.0, zero_stack, 1.0, int_stack+39318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33510,int_stack+52330,int_stack+40270, 0.0, zero_stack, 1.0, int_stack+52150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51850,int_stack+35462,int_stack+35362, 0.0, zero_stack, 1.0, int_stack+36524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34560,int_stack+35612,int_stack+35462, 0.0, zero_stack, 1.0, int_stack+36724, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35010,int_stack+34560,int_stack+51850, 0.0, zero_stack, 1.0, int_stack+53200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+36308,int_stack+36248, 1.0, int_stack+39258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+52330,int_stack+36398,int_stack+36308, 1.0, int_stack+39318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+34560,int_stack+52330,int_stack+40270, 1.0, int_stack+52150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41260,int_stack+36874,int_stack+36624, 1.0, int_stack+36524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+51850,int_stack+37024,int_stack+36874, 1.0, int_stack+36724, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35610,int_stack+51850,int_stack+41260, 1.0, int_stack+53200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+37570,int_stack+37510,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+53200,int_stack+37660,int_stack+37570,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+51850,int_stack+53200,int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53200,int_stack+37886,int_stack+37786,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+36210,int_stack+38036,int_stack+37886,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+36660,int_stack+36210,int_stack+53200,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+38582,int_stack+38522,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+53200,int_stack+38672,int_stack+38582,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+36210,int_stack+53200,int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53200,int_stack+38898,int_stack+38798,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+37260,int_stack+39048,int_stack+38898,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+37710,int_stack+37260,int_stack+53200,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+40270,int_stack+39594,int_stack+39534,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+53200,int_stack+39684,int_stack+39594,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+37260,int_stack+53200,int_stack+40270,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53200,int_stack+39910,int_stack+39810,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+38310,int_stack+40060,int_stack+39910,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+38760,int_stack+38310,int_stack+53200,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+39360,int_stack+42310,int_stack+40900,60);
     Libderiv->ABCD[11] = int_stack + 39360;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+59590,int_stack+44020,int_stack+42910,60);
     Libderiv->ABCD[10] = int_stack + 59590;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+60670,int_stack+45370,int_stack+0,60);
     Libderiv->ABCD[9] = int_stack + 60670;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+61750,int_stack+47080,int_stack+45970,60);
     Libderiv->ABCD[8] = int_stack + 61750;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+62830,int_stack+48430,int_stack+360,60);
     Libderiv->ABCD[7] = int_stack + 62830;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+63910,int_stack+49840,int_stack+49030,60);
     Libderiv->ABCD[6] = int_stack + 63910;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+64990,int_stack+50890,int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 64990;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+66070,int_stack+52600,int_stack+51490, 0.0, zero_stack, 1.0, int_stack+49390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 66070;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+52210,int_stack+53950,int_stack+1380, 1.0, int_stack+49390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 52210;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+67150,int_stack+2100,int_stack+1740,60);
     Libderiv->ABCD[155] = int_stack + 67150;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+68230,int_stack+54910,int_stack+54550,60);
     Libderiv->ABCD[143] = int_stack + 68230;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+53860,int_stack+3510,int_stack+2700,60);
     Libderiv->ABCD[142] = int_stack + 53860;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1740,int_stack+55510,int_stack+3060,60);
     Libderiv->ABCD[131] = int_stack + 1740;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2820,int_stack+4920,int_stack+4110,60);
     Libderiv->ABCD[130] = int_stack + 2820;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+54940,int_stack+5970,int_stack+4470,60);
     Libderiv->ABCD[129] = int_stack + 54940;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3900,int_stack+6570,int_stack+5520,60);
     Libderiv->ABCD[119] = int_stack + 3900;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4980,int_stack+7170,int_stack+56110,60);
     Libderiv->ABCD[118] = int_stack + 4980;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6060,int_stack+56830,int_stack+56470,60);
     Libderiv->ABCD[117] = int_stack + 6060;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+56020,int_stack+8580,int_stack+7770,60);
     Libderiv->ABCD[116] = int_stack + 56020;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+69310,int_stack+9630,int_stack+8130,60);
     Libderiv->ABCD[107] = int_stack + 69310;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7140,int_stack+10230,int_stack+9180,60);
     Libderiv->ABCD[106] = int_stack + 7140;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8220,int_stack+57790,int_stack+57430,60);
     Libderiv->ABCD[105] = int_stack + 8220;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9300,int_stack+11640,int_stack+10830,60);
     Libderiv->ABCD[104] = int_stack + 9300;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+57100,int_stack+58390,int_stack+11190,60);
     Libderiv->ABCD[103] = int_stack + 57100;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+10380,int_stack+12600,int_stack+12240,60);
     Libderiv->ABCD[95] = int_stack + 10380;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+11460,int_stack+13200,int_stack+40450,60);
     Libderiv->ABCD[94] = int_stack + 11460;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+12540,int_stack+13800,int_stack+41560,60);
     Libderiv->ABCD[93] = int_stack + 12540;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+41260,int_stack+14400,int_stack+43270,60);
     Libderiv->ABCD[92] = int_stack + 41260;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+43270,int_stack+15000,int_stack+44620,60);
     Libderiv->ABCD[91] = int_stack + 43270;
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+44350,int_stack+15600,int_stack+46330,60);
     Libderiv->ABCD[90] = int_stack + 44350;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+46330,int_stack+16200,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[47] = int_stack + 46330;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+47410,int_stack+17250,int_stack+53500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[46] = int_stack + 47410;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+13620,int_stack+18300,int_stack+16800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[45] = int_stack + 13620;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+14700,int_stack+19350,int_stack+17850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[44] = int_stack + 14700;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+15780,int_stack+20400,int_stack+18900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[43] = int_stack + 15780;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+16860,int_stack+21450,int_stack+19950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[42] = int_stack + 16860;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+17940,int_stack+22500,int_stack+21000, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[38] = int_stack + 17940;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+19020,int_stack+23550,int_stack+22050, 0.0, zero_stack, 1.0, int_stack+40900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[35] = int_stack + 19020;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+20100,int_stack+24600,int_stack+23100, 0.0, zero_stack, 1.0, int_stack+42910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[34] = int_stack + 20100;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+21180,int_stack+25650,int_stack+24150, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[33] = int_stack + 21180;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+22260,int_stack+26700,int_stack+25200, 0.0, zero_stack, 1.0, int_stack+45970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[32] = int_stack + 22260;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+23340,int_stack+58990,int_stack+26250, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[31] = int_stack + 23340;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+24420,int_stack+27660,int_stack+27300, 0.0, zero_stack, 1.0, int_stack+49030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[30] = int_stack + 24420;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+25500,int_stack+28710,int_stack+50440, 0.0, zero_stack, 1.0, int_stack+1020, 1.0, int_stack+51490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[26] = int_stack + 25500;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+26580,int_stack+29760,int_stack+28260, 0.0, zero_stack, 2.0, int_stack+51490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[25] = int_stack + 26580;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+27660,int_stack+30810,int_stack+29310, 1.0, int_stack+40900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[23] = int_stack + 27660;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+28740,int_stack+31860,int_stack+30360, 1.0, int_stack+42910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[22] = int_stack + 28740;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+29820,int_stack+32910,int_stack+31410, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[21] = int_stack + 29820;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+30900,int_stack+33960,int_stack+32460, 1.0, int_stack+45970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[20] = int_stack + 30900;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+31980,int_stack+35010,int_stack+33510, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[19] = int_stack + 31980;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+33060,int_stack+35610,int_stack+34560, 1.0, int_stack+49030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[18] = int_stack + 33060;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+48490,int_stack+36660,int_stack+51850, 1.0, int_stack+1020, 0.0, zero_stack, 1.0, int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[14] = int_stack + 48490;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+37710,int_stack+36210, 1.0, int_stack+51490, 1.0, int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[13] = int_stack + 0;
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+49570,int_stack+38760,int_stack+37260, 2.0, int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[12] = int_stack + 49570;

}
