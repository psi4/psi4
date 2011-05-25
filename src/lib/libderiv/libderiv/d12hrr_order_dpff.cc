#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_dpff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|ff) integrals */

void d12hrr_order_dpff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][6][11] = int_stack + 0;
 Libderiv->deriv_classes[3][6][10] = int_stack + 280;
 Libderiv->deriv_classes[3][6][9] = int_stack + 560;
 Libderiv->deriv_classes[3][6][8] = int_stack + 840;
 Libderiv->deriv_classes[3][6][7] = int_stack + 1120;
 Libderiv->dvrr_classes[3][5] = int_stack + 1400;
 Libderiv->deriv_classes[3][6][6] = int_stack + 1610;
 Libderiv->deriv_classes[3][6][2] = int_stack + 1890;
 Libderiv->deriv_classes[3][6][1] = int_stack + 2170;
 Libderiv->dvrr_classes[2][6] = int_stack + 2450;
 Libderiv->deriv_classes[3][6][0] = int_stack + 2618;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 2898;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 2958;
 Libderiv->deriv2_classes[2][5][143] = int_stack + 3048;
 Libderiv->deriv2_classes[2][6][143] = int_stack + 3174;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 3342;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 3442;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 3592;
 Libderiv->deriv2_classes[3][6][143] = int_stack + 3802;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 4082;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 4142;
 Libderiv->deriv2_classes[2][5][131] = int_stack + 4232;
 Libderiv->deriv2_classes[2][6][131] = int_stack + 4358;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 4526;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 4626;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 4776;
 Libderiv->deriv2_classes[3][6][131] = int_stack + 4986;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 5266;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 5326;
 Libderiv->deriv2_classes[2][5][130] = int_stack + 5416;
 Libderiv->deriv2_classes[2][6][130] = int_stack + 5542;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 5710;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 5810;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 5960;
 Libderiv->deriv2_classes[3][6][130] = int_stack + 6170;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 6450;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 6510;
 Libderiv->deriv2_classes[2][5][119] = int_stack + 6600;
 Libderiv->deriv2_classes[2][6][119] = int_stack + 6726;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 6894;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 6994;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 7144;
 Libderiv->deriv2_classes[3][6][119] = int_stack + 7354;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 7634;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 7694;
 Libderiv->deriv2_classes[2][5][118] = int_stack + 7784;
 Libderiv->deriv2_classes[2][6][118] = int_stack + 7910;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 8078;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 8178;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 8328;
 Libderiv->deriv2_classes[3][6][118] = int_stack + 8538;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 8818;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 8878;
 Libderiv->deriv2_classes[2][5][117] = int_stack + 8968;
 Libderiv->deriv2_classes[2][6][117] = int_stack + 9094;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 9262;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 9362;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 9512;
 Libderiv->deriv2_classes[3][6][117] = int_stack + 9722;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 10002;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 10062;
 Libderiv->deriv2_classes[2][5][107] = int_stack + 10152;
 Libderiv->deriv2_classes[2][6][107] = int_stack + 10278;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 10446;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 10546;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 10696;
 Libderiv->deriv2_classes[3][6][107] = int_stack + 10906;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 11186;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 11246;
 Libderiv->deriv2_classes[2][5][106] = int_stack + 11336;
 Libderiv->deriv2_classes[2][6][106] = int_stack + 11462;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 11630;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 11730;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 11880;
 Libderiv->deriv2_classes[3][6][106] = int_stack + 12090;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 12370;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 12430;
 Libderiv->deriv2_classes[2][5][105] = int_stack + 12520;
 Libderiv->deriv2_classes[2][6][105] = int_stack + 12646;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 12814;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 12914;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 13064;
 Libderiv->deriv2_classes[3][6][105] = int_stack + 13274;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 13554;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 13614;
 Libderiv->deriv2_classes[2][5][104] = int_stack + 13704;
 Libderiv->deriv2_classes[2][6][104] = int_stack + 13830;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 13998;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 14098;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 14248;
 Libderiv->deriv2_classes[3][6][104] = int_stack + 14458;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 14738;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 14798;
 Libderiv->deriv2_classes[2][5][95] = int_stack + 14888;
 Libderiv->deriv2_classes[2][6][95] = int_stack + 15014;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 15182;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 15282;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 15432;
 Libderiv->deriv2_classes[3][6][95] = int_stack + 15642;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 15922;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 15982;
 Libderiv->deriv2_classes[2][5][94] = int_stack + 16072;
 Libderiv->deriv2_classes[2][6][94] = int_stack + 16198;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 16366;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 16466;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 16616;
 Libderiv->deriv2_classes[3][6][94] = int_stack + 16826;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 17106;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 17166;
 Libderiv->deriv2_classes[2][5][93] = int_stack + 17256;
 Libderiv->deriv2_classes[2][6][93] = int_stack + 17382;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 17550;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 17650;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 17800;
 Libderiv->deriv2_classes[3][6][93] = int_stack + 18010;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 18290;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 18350;
 Libderiv->deriv2_classes[2][5][92] = int_stack + 18440;
 Libderiv->deriv2_classes[2][6][92] = int_stack + 18566;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 18734;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 18834;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 18984;
 Libderiv->deriv2_classes[3][6][92] = int_stack + 19194;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 19474;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 19534;
 Libderiv->deriv2_classes[2][5][91] = int_stack + 19624;
 Libderiv->deriv2_classes[2][6][91] = int_stack + 19750;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 19918;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 20018;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 20168;
 Libderiv->deriv2_classes[3][6][91] = int_stack + 20378;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 20658;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 20718;
 Libderiv->deriv2_classes[2][5][83] = int_stack + 20808;
 Libderiv->deriv2_classes[2][6][83] = int_stack + 20934;
 Libderiv->deriv_classes[3][3][11] = int_stack + 21102;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 21202;
 Libderiv->deriv_classes[3][4][11] = int_stack + 21302;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 21452;
 Libderiv->deriv_classes[3][5][11] = int_stack + 21602;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 21812;
 Libderiv->deriv2_classes[3][6][83] = int_stack + 22022;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 22302;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 22362;
 Libderiv->deriv2_classes[2][5][82] = int_stack + 22452;
 Libderiv->deriv2_classes[2][6][82] = int_stack + 22578;
 Libderiv->deriv_classes[3][3][10] = int_stack + 22746;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 22846;
 Libderiv->deriv_classes[3][4][10] = int_stack + 22946;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 23096;
 Libderiv->deriv_classes[3][5][10] = int_stack + 23246;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 23456;
 Libderiv->deriv2_classes[3][6][82] = int_stack + 23666;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 23946;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 24006;
 Libderiv->deriv2_classes[2][5][81] = int_stack + 24096;
 Libderiv->deriv2_classes[2][6][81] = int_stack + 24222;
 Libderiv->deriv_classes[3][3][9] = int_stack + 24390;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 24490;
 Libderiv->deriv_classes[3][4][9] = int_stack + 24590;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 24740;
 Libderiv->deriv_classes[3][5][9] = int_stack + 24890;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 25100;
 Libderiv->deriv2_classes[3][6][81] = int_stack + 25310;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 25590;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 25650;
 Libderiv->deriv2_classes[2][5][80] = int_stack + 25740;
 Libderiv->deriv2_classes[2][6][80] = int_stack + 25866;
 Libderiv->deriv_classes[3][3][8] = int_stack + 26034;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 26134;
 Libderiv->deriv_classes[3][4][8] = int_stack + 26234;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 26384;
 Libderiv->deriv_classes[3][5][8] = int_stack + 26534;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 26744;
 Libderiv->deriv2_classes[3][6][80] = int_stack + 26954;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 27234;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 27294;
 Libderiv->deriv2_classes[2][5][79] = int_stack + 27384;
 Libderiv->deriv2_classes[2][6][79] = int_stack + 27510;
 Libderiv->deriv_classes[3][3][7] = int_stack + 27678;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 27778;
 Libderiv->deriv_classes[3][4][7] = int_stack + 27878;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 28028;
 Libderiv->deriv_classes[3][5][7] = int_stack + 28178;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 28388;
 Libderiv->deriv2_classes[3][6][79] = int_stack + 28598;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 28878;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 28938;
 Libderiv->deriv2_classes[2][5][78] = int_stack + 29028;
 Libderiv->deriv2_classes[2][6][78] = int_stack + 29154;
 Libderiv->dvrr_classes[3][3] = int_stack + 29322;
 Libderiv->deriv_classes[3][3][6] = int_stack + 29422;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 29522;
 Libderiv->dvrr_classes[3][4] = int_stack + 29622;
 Libderiv->deriv_classes[3][4][6] = int_stack + 29772;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 29922;
 Libderiv->deriv_classes[3][5][6] = int_stack + 30072;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 30282;
 Libderiv->deriv2_classes[3][6][78] = int_stack + 30492;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 30772;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 30832;
 Libderiv->deriv2_classes[2][5][35] = int_stack + 30922;
 Libderiv->deriv2_classes[2][6][35] = int_stack + 31048;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 31216;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 31316;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 31466;
 Libderiv->deriv2_classes[3][6][35] = int_stack + 31676;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 31956;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 32016;
 Libderiv->deriv2_classes[2][5][34] = int_stack + 32106;
 Libderiv->deriv2_classes[2][6][34] = int_stack + 32232;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 32400;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 32500;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 32650;
 Libderiv->deriv2_classes[3][6][34] = int_stack + 32860;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 33140;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 33200;
 Libderiv->deriv2_classes[2][5][33] = int_stack + 33290;
 Libderiv->deriv2_classes[2][6][33] = int_stack + 33416;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 33584;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 33684;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 33834;
 Libderiv->deriv2_classes[3][6][33] = int_stack + 34044;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 34324;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 34384;
 Libderiv->deriv2_classes[2][5][32] = int_stack + 34474;
 Libderiv->deriv2_classes[2][6][32] = int_stack + 34600;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 34768;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 34868;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 35018;
 Libderiv->deriv2_classes[3][6][32] = int_stack + 35228;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 35508;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 35568;
 Libderiv->deriv2_classes[2][5][31] = int_stack + 35658;
 Libderiv->deriv2_classes[2][6][31] = int_stack + 35784;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 35952;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 36052;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 36202;
 Libderiv->deriv2_classes[3][6][31] = int_stack + 36412;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 36692;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 36752;
 Libderiv->deriv2_classes[2][5][30] = int_stack + 36842;
 Libderiv->deriv2_classes[2][6][30] = int_stack + 36968;
 Libderiv->deriv_classes[3][3][2] = int_stack + 37136;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 37236;
 Libderiv->deriv_classes[3][4][2] = int_stack + 37336;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 37486;
 Libderiv->deriv_classes[3][5][2] = int_stack + 37636;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 37846;
 Libderiv->deriv2_classes[3][6][30] = int_stack + 38056;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 38336;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 38396;
 Libderiv->deriv2_classes[2][5][26] = int_stack + 38486;
 Libderiv->deriv2_classes[2][6][26] = int_stack + 38612;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 38780;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 38880;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 39030;
 Libderiv->deriv2_classes[3][6][26] = int_stack + 39240;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 39520;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 39580;
 Libderiv->deriv2_classes[2][5][23] = int_stack + 39670;
 Libderiv->deriv2_classes[2][6][23] = int_stack + 39796;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 39964;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 40064;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 40214;
 Libderiv->deriv2_classes[3][6][23] = int_stack + 40424;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 40704;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 40764;
 Libderiv->deriv2_classes[2][5][22] = int_stack + 40854;
 Libderiv->deriv2_classes[2][6][22] = int_stack + 40980;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 41148;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 41248;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 41398;
 Libderiv->deriv2_classes[3][6][22] = int_stack + 41608;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 41888;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 41948;
 Libderiv->deriv2_classes[2][5][21] = int_stack + 42038;
 Libderiv->deriv2_classes[2][6][21] = int_stack + 42164;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 42332;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 42432;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 42582;
 Libderiv->deriv2_classes[3][6][21] = int_stack + 42792;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 43072;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 43132;
 Libderiv->deriv2_classes[2][5][20] = int_stack + 43222;
 Libderiv->deriv2_classes[2][6][20] = int_stack + 43348;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 43516;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 43616;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 43766;
 Libderiv->deriv2_classes[3][6][20] = int_stack + 43976;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 44256;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 44316;
 Libderiv->deriv2_classes[2][5][19] = int_stack + 44406;
 Libderiv->deriv2_classes[2][6][19] = int_stack + 44532;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 44700;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 44800;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 44950;
 Libderiv->deriv2_classes[3][6][19] = int_stack + 45160;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 45440;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 45500;
 Libderiv->deriv2_classes[2][5][18] = int_stack + 45590;
 Libderiv->deriv2_classes[2][6][18] = int_stack + 45716;
 Libderiv->deriv_classes[3][3][1] = int_stack + 45884;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 45984;
 Libderiv->deriv_classes[3][4][1] = int_stack + 46084;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 46234;
 Libderiv->deriv_classes[3][5][1] = int_stack + 46384;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 46594;
 Libderiv->deriv2_classes[3][6][18] = int_stack + 46804;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 47084;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 47144;
 Libderiv->deriv2_classes[2][5][14] = int_stack + 47234;
 Libderiv->deriv2_classes[2][6][14] = int_stack + 47360;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 47528;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 47628;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 47778;
 Libderiv->deriv2_classes[3][6][14] = int_stack + 47988;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 48268;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 48328;
 Libderiv->deriv2_classes[2][5][13] = int_stack + 48418;
 Libderiv->deriv2_classes[2][6][13] = int_stack + 48544;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 48712;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 48812;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 48962;
 Libderiv->deriv2_classes[3][6][13] = int_stack + 49172;
 Libderiv->deriv_classes[2][3][11] = int_stack + 49452;
 Libderiv->deriv_classes[2][4][11] = int_stack + 49512;
 Libderiv->deriv_classes[2][5][11] = int_stack + 49602;
 Libderiv->deriv_classes[2][6][11] = int_stack + 49728;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 49896;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 49956;
 Libderiv->deriv2_classes[2][5][11] = int_stack + 50046;
 Libderiv->deriv2_classes[2][6][11] = int_stack + 50172;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 50340;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 50440;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 50590;
 Libderiv->deriv2_classes[3][6][11] = int_stack + 50800;
 Libderiv->deriv_classes[2][3][10] = int_stack + 51080;
 Libderiv->deriv_classes[2][4][10] = int_stack + 51140;
 Libderiv->deriv_classes[2][5][10] = int_stack + 51230;
 Libderiv->deriv_classes[2][6][10] = int_stack + 51356;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 51524;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 51584;
 Libderiv->deriv2_classes[2][5][10] = int_stack + 51674;
 Libderiv->deriv2_classes[2][6][10] = int_stack + 51800;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 51968;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 52068;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 52218;
 Libderiv->deriv2_classes[3][6][10] = int_stack + 52428;
 Libderiv->deriv_classes[2][3][9] = int_stack + 52708;
 Libderiv->deriv_classes[2][4][9] = int_stack + 52768;
 Libderiv->deriv_classes[2][5][9] = int_stack + 52858;
 Libderiv->deriv_classes[2][6][9] = int_stack + 52984;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 53152;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 53212;
 Libderiv->deriv2_classes[2][5][9] = int_stack + 53302;
 Libderiv->deriv2_classes[2][6][9] = int_stack + 53428;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 53596;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 53696;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 53846;
 Libderiv->deriv2_classes[3][6][9] = int_stack + 54056;
 Libderiv->deriv_classes[2][3][8] = int_stack + 54336;
 Libderiv->deriv_classes[2][4][8] = int_stack + 54396;
 Libderiv->deriv_classes[2][5][8] = int_stack + 54486;
 Libderiv->deriv_classes[2][6][8] = int_stack + 54612;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 54780;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 54840;
 Libderiv->deriv2_classes[2][5][8] = int_stack + 54930;
 Libderiv->deriv2_classes[2][6][8] = int_stack + 55056;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 55224;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 55324;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 55474;
 Libderiv->deriv2_classes[3][6][8] = int_stack + 55684;
 Libderiv->deriv_classes[2][3][7] = int_stack + 55964;
 Libderiv->deriv_classes[2][4][7] = int_stack + 56024;
 Libderiv->deriv_classes[2][5][7] = int_stack + 56114;
 Libderiv->deriv_classes[2][6][7] = int_stack + 56240;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 56408;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 56468;
 Libderiv->deriv2_classes[2][5][7] = int_stack + 56558;
 Libderiv->deriv2_classes[2][6][7] = int_stack + 56684;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 56852;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 56952;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 57102;
 Libderiv->deriv2_classes[3][6][7] = int_stack + 57312;
 Libderiv->dvrr_classes[2][3] = int_stack + 57592;
 Libderiv->deriv_classes[2][3][6] = int_stack + 57652;
 Libderiv->dvrr_classes[2][4] = int_stack + 57712;
 Libderiv->deriv_classes[2][4][6] = int_stack + 57802;
 Libderiv->dvrr_classes[2][5] = int_stack + 57892;
 Libderiv->deriv_classes[2][5][6] = int_stack + 58018;
 Libderiv->deriv_classes[2][6][6] = int_stack + 58144;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 58312;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 58372;
 Libderiv->deriv2_classes[2][5][6] = int_stack + 58462;
 Libderiv->deriv2_classes[2][6][6] = int_stack + 58588;
 Libderiv->deriv_classes[3][3][0] = int_stack + 58756;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 58856;
 Libderiv->deriv_classes[3][4][0] = int_stack + 58956;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 59106;
 Libderiv->deriv_classes[3][5][0] = int_stack + 59256;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 59466;
 Libderiv->deriv2_classes[3][6][6] = int_stack + 59676;
 Libderiv->deriv_classes[2][3][2] = int_stack + 59956;
 Libderiv->deriv_classes[2][4][2] = int_stack + 60016;
 Libderiv->deriv_classes[2][5][2] = int_stack + 60106;
 Libderiv->deriv_classes[2][6][2] = int_stack + 60232;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 60400;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 60460;
 Libderiv->deriv2_classes[2][5][2] = int_stack + 60550;
 Libderiv->deriv2_classes[2][6][2] = int_stack + 60676;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 60844;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 60944;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 61094;
 Libderiv->deriv2_classes[3][6][2] = int_stack + 61304;
 Libderiv->deriv_classes[2][3][1] = int_stack + 61584;
 Libderiv->deriv_classes[2][4][1] = int_stack + 61644;
 Libderiv->deriv_classes[2][5][1] = int_stack + 61734;
 Libderiv->deriv_classes[2][6][1] = int_stack + 61860;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 62028;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 62088;
 Libderiv->deriv2_classes[2][5][1] = int_stack + 62178;
 Libderiv->deriv2_classes[2][6][1] = int_stack + 62304;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 62472;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 62572;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 62722;
 Libderiv->deriv2_classes[3][6][1] = int_stack + 62932;
 Libderiv->deriv_classes[2][3][0] = int_stack + 63212;
 Libderiv->deriv_classes[2][4][0] = int_stack + 63272;
 Libderiv->deriv_classes[2][5][0] = int_stack + 63362;
 Libderiv->deriv_classes[2][6][0] = int_stack + 63488;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 63656;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 63716;
 Libderiv->deriv2_classes[2][5][0] = int_stack + 63806;
 Libderiv->deriv2_classes[2][6][0] = int_stack + 63932;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 64100;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 64200;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 64350;
 Libderiv->deriv2_classes[3][6][0] = int_stack + 64560;
 memset(int_stack,0,518720);

 Libderiv->dvrr_stack = int_stack + 127224;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_dpff(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+64840,int_stack+57712,int_stack+57592,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+65020,int_stack+57892,int_stack+57712,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+65290,int_stack+65020,int_stack+64840,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+65650,int_stack+49512,int_stack+49452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57592,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65830,int_stack+49602,int_stack+49512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57712,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+66100,int_stack+65830,int_stack+65650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64840,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66460,int_stack+49728,int_stack+49602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57892,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+66838,int_stack+66460,int_stack+65830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65020,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+67378,int_stack+66838,int_stack+66100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65290,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+66460,int_stack+29622,int_stack+29322,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+66760,int_stack+1400,int_stack+29622,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+67978,int_stack+66760,int_stack+66460,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+68578,int_stack+21302,int_stack+21102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29322,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+68878,int_stack+21602,int_stack+21302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29622,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+69328,int_stack+68878,int_stack+68578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66460,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+69928,int_stack+0,int_stack+21602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1400,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+70558,int_stack+69928,int_stack+68878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66760,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+71458,int_stack+70558,int_stack+69328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67978,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+69928,int_stack+51140,int_stack+51080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57592, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+70108,int_stack+51230,int_stack+51140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57712, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+70378,int_stack+70108,int_stack+69928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64840, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+70738,int_stack+51356,int_stack+51230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57892, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+70738,int_stack+70108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65020, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+70738,int_stack+72458,int_stack+70378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65290, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72458,int_stack+22946,int_stack+22746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29322, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+72758,int_stack+23246,int_stack+22946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29622, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+73208,int_stack+72758,int_stack+72458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66460, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+73808,int_stack+280,int_stack+23246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1400, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+74438,int_stack+73808,int_stack+72758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66760, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+75338,int_stack+74438,int_stack+73208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67978, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73808,int_stack+52768,int_stack+52708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57592, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73988,int_stack+52858,int_stack+52768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57712, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+74258,int_stack+73988,int_stack+73808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64840, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+74618,int_stack+52984,int_stack+52858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57892, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+76338,int_stack+74618,int_stack+73988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65020, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+74618,int_stack+76338,int_stack+74258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65290, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+76338,int_stack+24590,int_stack+24390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29322, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+24890,int_stack+24590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29622, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+76638,int_stack+0,int_stack+76338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66460, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+77238,int_stack+560,int_stack+24890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1400, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+77868,int_stack+77238,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66760, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+78768,int_stack+77868,int_stack+76638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67978, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+77238,int_stack+54396,int_stack+54336, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+77418,int_stack+54486,int_stack+54396, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+77688,int_stack+77418,int_stack+77238, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+78048,int_stack+54612,int_stack+54486, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+79768,int_stack+78048,int_stack+77418, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+78048,int_stack+79768,int_stack+77688, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79768,int_stack+26234,int_stack+26034, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29322, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80068,int_stack+26534,int_stack+26234, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80518,int_stack+80068,int_stack+79768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81118,int_stack+840,int_stack+26534, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+81748,int_stack+81118,int_stack+80068, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+82648,int_stack+81748,int_stack+80518, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+81118,int_stack+56024,int_stack+55964, 0.0, zero_stack, 1.0, int_stack+57592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+81298,int_stack+56114,int_stack+56024, 0.0, zero_stack, 1.0, int_stack+57712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+81298,int_stack+81118, 0.0, zero_stack, 1.0, int_stack+64840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+81928,int_stack+56240,int_stack+56114, 0.0, zero_stack, 1.0, int_stack+57892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+450,int_stack+81928,int_stack+81298, 0.0, zero_stack, 1.0, int_stack+65020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+81928,int_stack+450,int_stack+81568, 0.0, zero_stack, 1.0, int_stack+65290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+27878,int_stack+27678, 0.0, zero_stack, 1.0, int_stack+29322, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83648,int_stack+28178,int_stack+27878, 0.0, zero_stack, 1.0, int_stack+29622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+84098,int_stack+83648,int_stack+450, 0.0, zero_stack, 1.0, int_stack+66460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84698,int_stack+1120,int_stack+28178, 0.0, zero_stack, 1.0, int_stack+1400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+85328,int_stack+84698,int_stack+83648, 0.0, zero_stack, 1.0, int_stack+66760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+86228,int_stack+85328,int_stack+84098, 0.0, zero_stack, 1.0, int_stack+67978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+84698,int_stack+57802,int_stack+57652, 1.0, int_stack+57592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+84878,int_stack+58018,int_stack+57802, 1.0, int_stack+57712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+85148,int_stack+84878,int_stack+84698, 1.0, int_stack+64840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+85508,int_stack+58144,int_stack+58018, 1.0, int_stack+57892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+750,int_stack+85508,int_stack+84878, 1.0, int_stack+65020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+85508,int_stack+750,int_stack+85148, 1.0, int_stack+65290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+750,int_stack+29772,int_stack+29422, 1.0, int_stack+29322, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87228,int_stack+30072,int_stack+29772, 1.0, int_stack+29622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+87678,int_stack+87228,int_stack+750, 1.0, int_stack+66460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88278,int_stack+1610,int_stack+30072, 1.0, int_stack+1400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+88908,int_stack+88278,int_stack+87228, 1.0, int_stack+66760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+89808,int_stack+88908,int_stack+87678, 1.0, int_stack+67978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+67978,int_stack+2450,int_stack+57892,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+88278,int_stack+67978,int_stack+65020,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+67978,int_stack+88278,int_stack+65290,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+88278,int_stack+60016,int_stack+59956,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+88458,int_stack+60106,int_stack+60016,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+88728,int_stack+88458,int_stack+88278,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+89088,int_stack+60232,int_stack+60106,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+66460,int_stack+89088,int_stack+88458,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+89088,int_stack+66460,int_stack+88728,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+66460,int_stack+37336,int_stack+37136,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+66760,int_stack+37636,int_stack+37336,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1050,int_stack+66760,int_stack+66460,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+64840,int_stack+1890,int_stack+37636,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+90808,int_stack+64840,int_stack+66760,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+91708,int_stack+90808,int_stack+1050,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+90808,int_stack+61644,int_stack+61584,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+90988,int_stack+61734,int_stack+61644,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+91258,int_stack+90988,int_stack+90808,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+64840,int_stack+61860,int_stack+61734,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+92708,int_stack+64840,int_stack+90988,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+64840,int_stack+92708,int_stack+91258,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+92708,int_stack+46084,int_stack+45884,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+93008,int_stack+46384,int_stack+46084,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+93458,int_stack+93008,int_stack+92708,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+94058,int_stack+2170,int_stack+46384,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1650,int_stack+94058,int_stack+93008,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+94058,int_stack+1650,int_stack+93458,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1650,int_stack+63272,int_stack+63212,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1830,int_stack+63362,int_stack+63272,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+1830,int_stack+1650,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+95058,int_stack+63488,int_stack+63362,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+95436,int_stack+95058,int_stack+1830,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+95976,int_stack+95436,int_stack+2100,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+95058,int_stack+58956,int_stack+58756,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+95358,int_stack+59256,int_stack+58956,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+96576,int_stack+95358,int_stack+95058,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+97176,int_stack+2618,int_stack+59256,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+97806,int_stack+97176,int_stack+95358,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+98706,int_stack+97806,int_stack+96576,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97176,int_stack+2958,int_stack+2898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49452,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97356,int_stack+3048,int_stack+2958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49512,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97626,int_stack+97356,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+65650,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+97986,int_stack+3174,int_stack+3048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49602,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2460,int_stack+97986,int_stack+97356, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+65830,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+97986,int_stack+2460,int_stack+97626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+66100,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2460,int_stack+3442,int_stack+3342, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+21102,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2760,int_stack+3592,int_stack+3442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+21302,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97176,int_stack+2760,int_stack+2460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+68578,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+99706,int_stack+3802,int_stack+3592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+21602,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+100336,int_stack+99706,int_stack+2760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+68878,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2460,int_stack+100336,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+69328,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97176,int_stack+4142,int_stack+4082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49452, 1.0, int_stack+51080,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97356,int_stack+4232,int_stack+4142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49512, 1.0, int_stack+51140,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97626,int_stack+97356,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65650, 1.0, int_stack+69928,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+99706,int_stack+4358,int_stack+4232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49602, 1.0, int_stack+51230,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+100084,int_stack+99706,int_stack+97356, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65830, 1.0, int_stack+70108,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+100624,int_stack+100084,int_stack+97626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66100, 1.0, int_stack+70378,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+4626,int_stack+4526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21102, 1.0, int_stack+22746,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+100006,int_stack+4776,int_stack+4626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21302, 1.0, int_stack+22946,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97176,int_stack+100006,int_stack+99706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68578, 1.0, int_stack+72458,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+101224,int_stack+4986,int_stack+4776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21602, 1.0, int_stack+23246,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3460,int_stack+101224,int_stack+100006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68878, 1.0, int_stack+72758,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+101224,int_stack+3460,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69328, 1.0, int_stack+73208,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97176,int_stack+5326,int_stack+5266, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+51080, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97356,int_stack+5416,int_stack+5326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+51140, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97626,int_stack+97356,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+69928, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3460,int_stack+5542,int_stack+5416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+51230, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3838,int_stack+3460,int_stack+97356, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+70108, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4378,int_stack+3838,int_stack+97626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+70378, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3460,int_stack+5810,int_stack+5710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22746, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3760,int_stack+5960,int_stack+5810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22946, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97176,int_stack+3760,int_stack+3460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+72458, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4978,int_stack+6170,int_stack+5960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+23246, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+99706,int_stack+4978,int_stack+3760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+72758, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4978,int_stack+99706,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+73208, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97176,int_stack+6510,int_stack+6450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49452, 0.0, zero_stack, 1.0, int_stack+52708,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97356,int_stack+6600,int_stack+6510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49512, 0.0, zero_stack, 1.0, int_stack+52768,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97626,int_stack+97356,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65650, 0.0, zero_stack, 1.0, int_stack+73808,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+99706,int_stack+6726,int_stack+6600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49602, 0.0, zero_stack, 1.0, int_stack+52858,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+100084,int_stack+99706,int_stack+97356, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65830, 0.0, zero_stack, 1.0, int_stack+73988,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5978,int_stack+100084,int_stack+97626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66100, 0.0, zero_stack, 1.0, int_stack+74258,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+6994,int_stack+6894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21102, 0.0, zero_stack, 1.0, int_stack+24390,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+100006,int_stack+7144,int_stack+6994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21302, 0.0, zero_stack, 1.0, int_stack+24590,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97176,int_stack+100006,int_stack+99706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68578, 0.0, zero_stack, 1.0, int_stack+76338,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3460,int_stack+7354,int_stack+7144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21602, 0.0, zero_stack, 1.0, int_stack+24890,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6578,int_stack+3460,int_stack+100006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68878, 0.0, zero_stack, 1.0, int_stack+0,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+102224,int_stack+6578,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69328, 0.0, zero_stack, 1.0, int_stack+76638,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97176,int_stack+7694,int_stack+7634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51080, 1.0, int_stack+52708, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97356,int_stack+7784,int_stack+7694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51140, 1.0, int_stack+52768, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97626,int_stack+97356,int_stack+97176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69928, 1.0, int_stack+73808, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6578,int_stack+7910,int_stack+7784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51230, 1.0, int_stack+52858, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6956,int_stack+6578,int_stack+97356, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70108, 1.0, int_stack+73988, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3460,int_stack+6956,int_stack+97626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70378, 1.0, int_stack+74258, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6578,int_stack+8178,int_stack+8078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22746, 1.0, int_stack+24390, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6878,int_stack+8328,int_stack+8178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22946, 1.0, int_stack+24590, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7328,int_stack+6878,int_stack+6578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72458, 1.0, int_stack+76338, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+97176,int_stack+8538,int_stack+8328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23246, 1.0, int_stack+24890, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+99706,int_stack+97176,int_stack+6878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72758, 1.0, int_stack+0, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+103224,int_stack+99706,int_stack+7328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73208, 1.0, int_stack+76638, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+8878,int_stack+8818, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+52708, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+99886,int_stack+8968,int_stack+8878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+52768, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+100156,int_stack+99886,int_stack+99706, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+73808, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+97176,int_stack+9094,int_stack+8968, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+52858, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6578,int_stack+97176,int_stack+99886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+73988, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+97176,int_stack+6578,int_stack+100156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+74258, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6578,int_stack+9362,int_stack+9262, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24390, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6878,int_stack+9512,int_stack+9362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24590, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7328,int_stack+6878,int_stack+6578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+76338, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7928,int_stack+9722,int_stack+9512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24890, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8558,int_stack+7928,int_stack+6878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104224,int_stack+8558,int_stack+7328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+76638, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6578,int_stack+10062,int_stack+10002, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49452, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54336,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6758,int_stack+10152,int_stack+10062, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49512, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54396,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7028,int_stack+6758,int_stack+6578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77238,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7388,int_stack+10278,int_stack+10152, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49602, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54486,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7766,int_stack+7388,int_stack+6758, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65830, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77418,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8306,int_stack+7766,int_stack+7028, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77688,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6578,int_stack+10546,int_stack+10446, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21102, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26034,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6878,int_stack+10696,int_stack+10546, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21302, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26234,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7328,int_stack+6878,int_stack+6578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79768,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8906,int_stack+10906,int_stack+10696, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21602, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26534,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9536,int_stack+8906,int_stack+6878, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68878, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80068,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+105224,int_stack+9536,int_stack+7328, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69328, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80518,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8906,int_stack+11246,int_stack+11186, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51080, 0.0, zero_stack, 1.0, int_stack+54336, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9086,int_stack+11336,int_stack+11246, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51140, 0.0, zero_stack, 1.0, int_stack+54396, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9356,int_stack+9086,int_stack+8906, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69928, 0.0, zero_stack, 1.0, int_stack+77238, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9716,int_stack+11462,int_stack+11336, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51230, 0.0, zero_stack, 1.0, int_stack+54486, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10094,int_stack+9716,int_stack+9086, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70108, 0.0, zero_stack, 1.0, int_stack+77418, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10634,int_stack+10094,int_stack+9356, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70378, 0.0, zero_stack, 1.0, int_stack+77688, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8906,int_stack+11730,int_stack+11630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22746, 0.0, zero_stack, 1.0, int_stack+26034, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9206,int_stack+11880,int_stack+11730, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22946, 0.0, zero_stack, 1.0, int_stack+26234, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9656,int_stack+9206,int_stack+8906, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72458, 0.0, zero_stack, 1.0, int_stack+79768, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11234,int_stack+12090,int_stack+11880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23246, 0.0, zero_stack, 1.0, int_stack+26534, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6578,int_stack+11234,int_stack+9206, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72758, 0.0, zero_stack, 1.0, int_stack+80068, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11234,int_stack+6578,int_stack+9656, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73208, 0.0, zero_stack, 1.0, int_stack+80518, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6578,int_stack+12430,int_stack+12370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52708, 1.0, int_stack+54336, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6758,int_stack+12520,int_stack+12430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52768, 1.0, int_stack+54396, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7028,int_stack+6758,int_stack+6578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73808, 1.0, int_stack+77238, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7388,int_stack+12646,int_stack+12520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52858, 1.0, int_stack+54486, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7766,int_stack+7388,int_stack+6758, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73988, 1.0, int_stack+77418, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8906,int_stack+7766,int_stack+7028, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74258, 1.0, int_stack+77688, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6578,int_stack+12914,int_stack+12814, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24390, 1.0, int_stack+26034, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6878,int_stack+13064,int_stack+12914, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24590, 1.0, int_stack+26234, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7328,int_stack+6878,int_stack+6578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76338, 1.0, int_stack+79768, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12234,int_stack+13274,int_stack+13064, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24890, 1.0, int_stack+26534, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9506,int_stack+12234,int_stack+6878, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+80068, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12234,int_stack+9506,int_stack+7328, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76638, 1.0, int_stack+80518, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9506,int_stack+13614,int_stack+13554, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+54336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9686,int_stack+13704,int_stack+13614, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+54396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9956,int_stack+9686,int_stack+9506, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+77238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13234,int_stack+13830,int_stack+13704, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+54486, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6578,int_stack+13234,int_stack+9686, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+77418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13234,int_stack+6578,int_stack+9956, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+77688, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6578,int_stack+14098,int_stack+13998, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+26034, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6878,int_stack+14248,int_stack+14098, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+26234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7328,int_stack+6878,int_stack+6578, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+79768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9506,int_stack+14458,int_stack+14248, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+26534, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13834,int_stack+9506,int_stack+6878, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+80068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9506,int_stack+13834,int_stack+7328, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+80518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13834,int_stack+14798,int_stack+14738, 0.0, zero_stack, 1.0, int_stack+49452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55964,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14014,int_stack+14888,int_stack+14798, 0.0, zero_stack, 1.0, int_stack+49512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56024,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14284,int_stack+14014,int_stack+13834, 0.0, zero_stack, 1.0, int_stack+65650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81118,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6578,int_stack+15014,int_stack+14888, 0.0, zero_stack, 1.0, int_stack+49602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56114,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6956,int_stack+6578,int_stack+14014, 0.0, zero_stack, 1.0, int_stack+65830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81298,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7496,int_stack+6956,int_stack+14284, 0.0, zero_stack, 1.0, int_stack+66100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81568,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6578,int_stack+15282,int_stack+15182, 0.0, zero_stack, 1.0, int_stack+21102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27678,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6878,int_stack+15432,int_stack+15282, 0.0, zero_stack, 1.0, int_stack+21302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27878,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13834,int_stack+6878,int_stack+6578, 0.0, zero_stack, 1.0, int_stack+68578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14434,int_stack+15642,int_stack+15432, 0.0, zero_stack, 1.0, int_stack+21602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28178,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+99706,int_stack+14434,int_stack+6878, 0.0, zero_stack, 1.0, int_stack+68878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+83648,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14434,int_stack+99706,int_stack+13834, 0.0, zero_stack, 1.0, int_stack+69328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84098,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13834,int_stack+15982,int_stack+15922, 0.0, zero_stack, 1.0, int_stack+51080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55964, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14014,int_stack+16072,int_stack+15982, 0.0, zero_stack, 1.0, int_stack+51140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56024, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+99706,int_stack+14014,int_stack+13834, 0.0, zero_stack, 1.0, int_stack+69928, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81118, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+100066,int_stack+16198,int_stack+16072, 0.0, zero_stack, 1.0, int_stack+51230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56114, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15434,int_stack+100066,int_stack+14014, 0.0, zero_stack, 1.0, int_stack+70108, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81298, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13834,int_stack+15434,int_stack+99706, 0.0, zero_stack, 1.0, int_stack+70378, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81568, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+16466,int_stack+16366, 0.0, zero_stack, 1.0, int_stack+22746, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27678, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+100006,int_stack+16616,int_stack+16466, 0.0, zero_stack, 1.0, int_stack+22946, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27878, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15434,int_stack+100006,int_stack+99706, 0.0, zero_stack, 1.0, int_stack+72458, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6578,int_stack+16826,int_stack+16616, 0.0, zero_stack, 1.0, int_stack+23246, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28178, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16034,int_stack+6578,int_stack+100006, 0.0, zero_stack, 1.0, int_stack+72758, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+83648, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+106224,int_stack+16034,int_stack+15434, 0.0, zero_stack, 1.0, int_stack+73208, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84098, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15434,int_stack+17166,int_stack+17106, 0.0, zero_stack, 1.0, int_stack+52708, 0.0, zero_stack, 1.0, int_stack+55964, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15614,int_stack+17256,int_stack+17166, 0.0, zero_stack, 1.0, int_stack+52768, 0.0, zero_stack, 1.0, int_stack+56024, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15884,int_stack+15614,int_stack+15434, 0.0, zero_stack, 1.0, int_stack+73808, 0.0, zero_stack, 1.0, int_stack+81118, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16244,int_stack+17382,int_stack+17256, 0.0, zero_stack, 1.0, int_stack+52858, 0.0, zero_stack, 1.0, int_stack+56114, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16622,int_stack+16244,int_stack+15614, 0.0, zero_stack, 1.0, int_stack+73988, 0.0, zero_stack, 1.0, int_stack+81298, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6578,int_stack+16622,int_stack+15884, 0.0, zero_stack, 1.0, int_stack+74258, 0.0, zero_stack, 1.0, int_stack+81568, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15434,int_stack+17650,int_stack+17550, 0.0, zero_stack, 1.0, int_stack+24390, 0.0, zero_stack, 1.0, int_stack+27678, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15734,int_stack+17800,int_stack+17650, 0.0, zero_stack, 1.0, int_stack+24590, 0.0, zero_stack, 1.0, int_stack+27878, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16184,int_stack+15734,int_stack+15434, 0.0, zero_stack, 1.0, int_stack+76338, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16784,int_stack+18010,int_stack+17800, 0.0, zero_stack, 1.0, int_stack+24890, 0.0, zero_stack, 1.0, int_stack+28178, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+99706,int_stack+16784,int_stack+15734, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+83648, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16784,int_stack+99706,int_stack+16184, 0.0, zero_stack, 1.0, int_stack+76638, 0.0, zero_stack, 1.0, int_stack+84098, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+18350,int_stack+18290, 0.0, zero_stack, 1.0, int_stack+54336, 1.0, int_stack+55964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+99886,int_stack+18440,int_stack+18350, 0.0, zero_stack, 1.0, int_stack+54396, 1.0, int_stack+56024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+100156,int_stack+99886,int_stack+99706, 0.0, zero_stack, 1.0, int_stack+77238, 1.0, int_stack+81118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17784,int_stack+18566,int_stack+18440, 0.0, zero_stack, 1.0, int_stack+54486, 1.0, int_stack+56114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18162,int_stack+17784,int_stack+99886, 0.0, zero_stack, 1.0, int_stack+77418, 1.0, int_stack+81298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+15434,int_stack+18162,int_stack+100156, 0.0, zero_stack, 1.0, int_stack+77688, 1.0, int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17784,int_stack+18834,int_stack+18734, 0.0, zero_stack, 1.0, int_stack+26034, 1.0, int_stack+27678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18084,int_stack+18984,int_stack+18834, 0.0, zero_stack, 1.0, int_stack+26234, 1.0, int_stack+27878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+99706,int_stack+18084,int_stack+17784, 0.0, zero_stack, 1.0, int_stack+79768, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16034,int_stack+19194,int_stack+18984, 0.0, zero_stack, 1.0, int_stack+26534, 1.0, int_stack+28178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18534,int_stack+16034,int_stack+18084, 0.0, zero_stack, 1.0, int_stack+80068, 1.0, int_stack+83648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+107224,int_stack+18534,int_stack+99706, 0.0, zero_stack, 1.0, int_stack+80518, 1.0, int_stack+84098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+19534,int_stack+19474, 0.0, zero_stack, 2.0, int_stack+55964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+99886,int_stack+19624,int_stack+19534, 0.0, zero_stack, 2.0, int_stack+56024, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+100156,int_stack+99886,int_stack+99706, 0.0, zero_stack, 2.0, int_stack+81118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16034,int_stack+19750,int_stack+19624, 0.0, zero_stack, 2.0, int_stack+56114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17784,int_stack+16034,int_stack+99886, 0.0, zero_stack, 2.0, int_stack+81298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16034,int_stack+17784,int_stack+100156, 0.0, zero_stack, 2.0, int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17784,int_stack+20018,int_stack+19918, 0.0, zero_stack, 2.0, int_stack+27678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18084,int_stack+20168,int_stack+20018, 0.0, zero_stack, 2.0, int_stack+27878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18534,int_stack+18084,int_stack+17784, 0.0, zero_stack, 2.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+19134,int_stack+20378,int_stack+20168, 0.0, zero_stack, 2.0, int_stack+28178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+99706,int_stack+19134,int_stack+18084, 0.0, zero_stack, 2.0, int_stack+83648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+19134,int_stack+99706,int_stack+18534, 0.0, zero_stack, 2.0, int_stack+84098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+20718,int_stack+20658, 1.0, int_stack+49452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57652,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+99886,int_stack+20808,int_stack+20718, 1.0, int_stack+49512, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57802,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+100156,int_stack+99886,int_stack+99706, 1.0, int_stack+65650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84698,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20134,int_stack+20934,int_stack+20808, 1.0, int_stack+49602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58018,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20512,int_stack+20134,int_stack+99886, 1.0, int_stack+65830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84878,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+17784,int_stack+20512,int_stack+100156, 1.0, int_stack+66100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85148,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20134,int_stack+21452,int_stack+21202, 1.0, int_stack+21102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29422,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20434,int_stack+21812,int_stack+21452, 1.0, int_stack+21302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29772,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20884,int_stack+20434,int_stack+20134, 1.0, int_stack+68578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+99706,int_stack+22022,int_stack+21812, 1.0, int_stack+21602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30072,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+65440,int_stack+99706,int_stack+20434, 1.0, int_stack+68878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87228,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+108224,int_stack+65440,int_stack+20884, 1.0, int_stack+69328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87678,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+65440,int_stack+22362,int_stack+22302, 1.0, int_stack+51080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57652, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65620,int_stack+22452,int_stack+22362, 1.0, int_stack+51140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57802, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65890,int_stack+65620,int_stack+65440, 1.0, int_stack+69928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84698, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+99706,int_stack+22578,int_stack+22452, 1.0, int_stack+51230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58018, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+100084,int_stack+99706,int_stack+65620, 1.0, int_stack+70108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84878, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+68578,int_stack+100084,int_stack+65890, 1.0, int_stack+70378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85148, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+23096,int_stack+22846, 1.0, int_stack+22746, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29422, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+100006,int_stack+23456,int_stack+23096, 1.0, int_stack+22946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29772, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65440,int_stack+100006,int_stack+99706, 1.0, int_stack+72458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+69178,int_stack+23666,int_stack+23456, 1.0, int_stack+23246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30072, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+69808,int_stack+69178,int_stack+100006, 1.0, int_stack+72758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87228, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+20134,int_stack+69808,int_stack+65440, 1.0, int_stack+73208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87678, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+65440,int_stack+24006,int_stack+23946, 1.0, int_stack+52708, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57652, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65620,int_stack+24096,int_stack+24006, 1.0, int_stack+52768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57802, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65890,int_stack+65620,int_stack+65440, 1.0, int_stack+73808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84698, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+69178,int_stack+24222,int_stack+24096, 1.0, int_stack+52858, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58018, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+69556,int_stack+69178,int_stack+65620, 1.0, int_stack+73988, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84878, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+70096,int_stack+69556,int_stack+65890, 1.0, int_stack+74258, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85148, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+69178,int_stack+24740,int_stack+24490, 1.0, int_stack+24390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29422, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+69478,int_stack+25100,int_stack+24740, 1.0, int_stack+24590, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29772, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65440,int_stack+69478,int_stack+69178, 1.0, int_stack+76338, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+72458,int_stack+25310,int_stack+25100, 1.0, int_stack+24890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30072, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+73088,int_stack+72458,int_stack+69478, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87228, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+21134,int_stack+73088,int_stack+65440, 1.0, int_stack+76638, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87678, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+65440,int_stack+25650,int_stack+25590, 1.0, int_stack+54336, 0.0, zero_stack, 1.0, int_stack+57652, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65620,int_stack+25740,int_stack+25650, 1.0, int_stack+54396, 0.0, zero_stack, 1.0, int_stack+57802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65890,int_stack+65620,int_stack+65440, 1.0, int_stack+77238, 0.0, zero_stack, 1.0, int_stack+84698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+25866,int_stack+25740, 1.0, int_stack+54486, 0.0, zero_stack, 1.0, int_stack+58018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+0,int_stack+65620, 1.0, int_stack+77418, 0.0, zero_stack, 1.0, int_stack+84878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+72998,int_stack+72458,int_stack+65890, 1.0, int_stack+77688, 0.0, zero_stack, 1.0, int_stack+85148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72458,int_stack+26384,int_stack+26134, 1.0, int_stack+26034, 0.0, zero_stack, 1.0, int_stack+29422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+26744,int_stack+26384, 1.0, int_stack+26234, 0.0, zero_stack, 1.0, int_stack+29772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65440,int_stack+0,int_stack+72458, 1.0, int_stack+79768, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+73598,int_stack+26954,int_stack+26744, 1.0, int_stack+26534, 0.0, zero_stack, 1.0, int_stack+30072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+76338,int_stack+73598,int_stack+0, 1.0, int_stack+80068, 0.0, zero_stack, 1.0, int_stack+87228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+73598,int_stack+76338,int_stack+65440, 1.0, int_stack+80518, 0.0, zero_stack, 1.0, int_stack+87678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+65440,int_stack+27294,int_stack+27234, 1.0, int_stack+55964, 1.0, int_stack+57652, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65620,int_stack+27384,int_stack+27294, 1.0, int_stack+56024, 1.0, int_stack+57802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65890,int_stack+65620,int_stack+65440, 1.0, int_stack+81118, 1.0, int_stack+84698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+76338,int_stack+27510,int_stack+27384, 1.0, int_stack+56114, 1.0, int_stack+58018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+76338,int_stack+65620, 1.0, int_stack+81298, 1.0, int_stack+84878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+76338,int_stack+72458,int_stack+65890, 1.0, int_stack+81568, 1.0, int_stack+85148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72458,int_stack+28028,int_stack+27778, 1.0, int_stack+27678, 1.0, int_stack+29422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+28388,int_stack+28028, 1.0, int_stack+27878, 1.0, int_stack+29772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+76938,int_stack+0,int_stack+72458, 1.0, int_stack+450, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+65440,int_stack+28598,int_stack+28388, 1.0, int_stack+28178, 1.0, int_stack+30072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+79768,int_stack+65440,int_stack+0, 1.0, int_stack+83648, 1.0, int_stack+87228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+65440,int_stack+79768,int_stack+76938, 1.0, int_stack+84098, 1.0, int_stack+87678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+76938,int_stack+28938,int_stack+28878, 2.0, int_stack+57652, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+77118,int_stack+29028,int_stack+28938, 2.0, int_stack+57802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+77388,int_stack+77118,int_stack+76938, 2.0, int_stack+84698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+79768,int_stack+29154,int_stack+29028, 2.0, int_stack+58018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+79768,int_stack+77118, 2.0, int_stack+84878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+79768,int_stack+72458,int_stack+77388, 2.0, int_stack+85148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72458,int_stack+29922,int_stack+29522, 2.0, int_stack+29422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80368,int_stack+30282,int_stack+29922, 2.0, int_stack+29772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80818,int_stack+80368,int_stack+72458, 2.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+76938,int_stack+30492,int_stack+30282, 2.0, int_stack+30072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83648,int_stack+76938,int_stack+80368, 2.0, int_stack+87228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+76938,int_stack+83648,int_stack+80818, 2.0, int_stack+87678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+83648,int_stack+30832,int_stack+30772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59956,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83828,int_stack+30922,int_stack+30832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60016,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+84098,int_stack+83828,int_stack+83648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88278,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+84458,int_stack+31048,int_stack+30922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60106,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+84458,int_stack+83828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88458,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+84458,int_stack+72458,int_stack+84098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88728,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72458,int_stack+31316,int_stack+31216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37136,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+31466,int_stack+31316, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37336,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+83648,int_stack+85058,int_stack+72458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66460,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+87228,int_stack+31676,int_stack+31466, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37636,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+80368,int_stack+87228,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66760,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+87228,int_stack+80368,int_stack+83648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+83648,int_stack+32016,int_stack+31956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59956, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83828,int_stack+32106,int_stack+32016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60016, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+84098,int_stack+83828,int_stack+83648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88278, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80368,int_stack+32232,int_stack+32106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60106, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+80368,int_stack+83828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88458, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+80368,int_stack+72458,int_stack+84098, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88728, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72458,int_stack+32500,int_stack+32400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37136, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+32650,int_stack+32500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37336, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80968,int_stack+85058,int_stack+72458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66460, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83648,int_stack+32860,int_stack+32650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37636, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+83648,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66760, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+22134,int_stack+0,int_stack+80968, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+80968,int_stack+33200,int_stack+33140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59956, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+81148,int_stack+33290,int_stack+33200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60016, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81418,int_stack+81148,int_stack+80968, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88278, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+0,int_stack+33416,int_stack+33290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60106, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+0,int_stack+81148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88458, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+72458,int_stack+81418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88728, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72458,int_stack+33684,int_stack+33584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37136, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+600,int_stack+33834,int_stack+33684, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37336, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80968,int_stack+600,int_stack+72458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66460, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83648,int_stack+34044,int_stack+33834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37636, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+69178,int_stack+83648,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66760, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+23134,int_stack+69178,int_stack+80968, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+80968,int_stack+34384,int_stack+34324, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+81148,int_stack+34474,int_stack+34384, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81418,int_stack+81148,int_stack+80968, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+69178,int_stack+34600,int_stack+34474, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+69556,int_stack+69178,int_stack+81148, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+83648,int_stack+69556,int_stack+81418, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+69178,int_stack+34868,int_stack+34768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+600,int_stack+35018,int_stack+34868, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+69478,int_stack+600,int_stack+69178, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80968,int_stack+35228,int_stack+35018, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+99706,int_stack+80968,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+24134,int_stack+99706,int_stack+69478, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+35568,int_stack+35508, 0.0, zero_stack, 1.0, int_stack+59956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+99886,int_stack+35658,int_stack+35568, 0.0, zero_stack, 1.0, int_stack+60016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+100156,int_stack+99886,int_stack+99706, 0.0, zero_stack, 1.0, int_stack+88278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+600,int_stack+35784,int_stack+35658, 0.0, zero_stack, 1.0, int_stack+60106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+600,int_stack+99886, 0.0, zero_stack, 1.0, int_stack+88458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+80968,int_stack+72458,int_stack+100156, 0.0, zero_stack, 1.0, int_stack+88728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72458,int_stack+36052,int_stack+35952, 0.0, zero_stack, 1.0, int_stack+37136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+600,int_stack+36202,int_stack+36052, 0.0, zero_stack, 1.0, int_stack+37336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+99706,int_stack+600,int_stack+72458, 0.0, zero_stack, 1.0, int_stack+66460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+69178,int_stack+36412,int_stack+36202, 0.0, zero_stack, 1.0, int_stack+37636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+25134,int_stack+69178,int_stack+600, 0.0, zero_stack, 1.0, int_stack+66760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+26034,int_stack+25134,int_stack+99706, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+99706,int_stack+36752,int_stack+36692, 1.0, int_stack+59956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+99886,int_stack+36842,int_stack+36752, 1.0, int_stack+60016, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+99886,int_stack+99706, 1.0, int_stack+88278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+100156,int_stack+36968,int_stack+36842, 1.0, int_stack+60106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+100156,int_stack+99886, 1.0, int_stack+88458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+99706,int_stack+72458,int_stack+81568, 1.0, int_stack+88728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+81568,int_stack+37486,int_stack+37236, 1.0, int_stack+37136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+600,int_stack+37846,int_stack+37486, 1.0, int_stack+37336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25134,int_stack+600,int_stack+81568, 1.0, int_stack+66460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+69178,int_stack+38056,int_stack+37846, 1.0, int_stack+37636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27034,int_stack+69178,int_stack+600, 1.0, int_stack+66760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+27934,int_stack+27034,int_stack+25134, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+25134,int_stack+38396,int_stack+38336,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+25314,int_stack+38486,int_stack+38396,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+25314,int_stack+25134,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+25584,int_stack+38612,int_stack+38486,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+25584,int_stack+25314,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+25134,int_stack+72458,int_stack+81568,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+25734,int_stack+38880,int_stack+38780,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+39030,int_stack+38880,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+27034,int_stack+85058,int_stack+25734,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+600,int_stack+39240,int_stack+39030,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+69178,int_stack+600,int_stack+85058,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+600,int_stack+69178,int_stack+27034,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27034,int_stack+39580,int_stack+39520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61584,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+27214,int_stack+39670,int_stack+39580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61644,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+27214,int_stack+27034, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90808,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+27484,int_stack+39796,int_stack+39670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61734,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+27484,int_stack+27214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90988,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+27034,int_stack+72458,int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91258,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27634,int_stack+40064,int_stack+39964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45884,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+40214,int_stack+40064, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46084,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+69178,int_stack+85058,int_stack+27634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92708,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88228,int_stack+40424,int_stack+40214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46384,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+66440,int_stack+88228,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93008,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+28934,int_stack+66440,int_stack+69178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93458,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+69178,int_stack+40764,int_stack+40704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61584, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+69358,int_stack+40854,int_stack+40764, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61644, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+69358,int_stack+69178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90808, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+69628,int_stack+40980,int_stack+40854, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61734, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+69628,int_stack+69358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90988, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+69178,int_stack+72458,int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91258, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27634,int_stack+41248,int_stack+41148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45884, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+41398,int_stack+41248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46084, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+66440,int_stack+85058,int_stack+27634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92708, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88228,int_stack+41608,int_stack+41398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46384, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29934,int_stack+88228,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93008, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+30834,int_stack+29934,int_stack+66440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93458, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+66440,int_stack+41948,int_stack+41888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61584, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66620,int_stack+42038,int_stack+41948, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61644, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+66620,int_stack+66440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90808, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66890,int_stack+42164,int_stack+42038, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61734, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+66890,int_stack+66620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90988, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+66440,int_stack+72458,int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91258, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27634,int_stack+42432,int_stack+42332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45884, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+42582,int_stack+42432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46084, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29934,int_stack+85058,int_stack+27634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92708, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88228,int_stack+42792,int_stack+42582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46384, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31834,int_stack+88228,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93008, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+32734,int_stack+31834,int_stack+29934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93458, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29934,int_stack+43132,int_stack+43072, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+30114,int_stack+43222,int_stack+43132, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+30114,int_stack+29934, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30384,int_stack+43348,int_stack+43222, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+30384,int_stack+30114, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+29934,int_stack+72458,int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30534,int_stack+43616,int_stack+43516, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45884, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+43766,int_stack+43616, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31834,int_stack+85058,int_stack+30534, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88228,int_stack+43976,int_stack+43766, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33734,int_stack+88228,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+34634,int_stack+33734,int_stack+31834, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+31834,int_stack+44316,int_stack+44256, 0.0, zero_stack, 1.0, int_stack+61584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32014,int_stack+44406,int_stack+44316, 0.0, zero_stack, 1.0, int_stack+61644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+32014,int_stack+31834, 0.0, zero_stack, 1.0, int_stack+90808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+32284,int_stack+44532,int_stack+44406, 0.0, zero_stack, 1.0, int_stack+61734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+32284,int_stack+32014, 0.0, zero_stack, 1.0, int_stack+90988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+31834,int_stack+72458,int_stack+81568, 0.0, zero_stack, 1.0, int_stack+91258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+32434,int_stack+44800,int_stack+44700, 0.0, zero_stack, 1.0, int_stack+45884, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+44950,int_stack+44800, 0.0, zero_stack, 1.0, int_stack+46084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33734,int_stack+85058,int_stack+32434, 0.0, zero_stack, 1.0, int_stack+92708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+88228,int_stack+45160,int_stack+44950, 0.0, zero_stack, 1.0, int_stack+46384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+35634,int_stack+88228,int_stack+85058, 0.0, zero_stack, 1.0, int_stack+93008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+36534,int_stack+35634,int_stack+33734, 0.0, zero_stack, 1.0, int_stack+93458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33734,int_stack+45500,int_stack+45440, 1.0, int_stack+61584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33914,int_stack+45590,int_stack+45500, 1.0, int_stack+61644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+33914,int_stack+33734, 1.0, int_stack+90808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+34184,int_stack+45716,int_stack+45590, 1.0, int_stack+61734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+34184,int_stack+33914, 1.0, int_stack+90988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+33734,int_stack+72458,int_stack+81568, 1.0, int_stack+91258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+34334,int_stack+46234,int_stack+45984, 1.0, int_stack+45884, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+46594,int_stack+46234, 1.0, int_stack+46084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+90808,int_stack+85058,int_stack+34334, 1.0, int_stack+92708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+46804,int_stack+46594, 1.0, int_stack+46384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+37534,int_stack+35634,int_stack+85058, 1.0, int_stack+93008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+38434,int_stack+37534,int_stack+90808, 1.0, int_stack+93458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+90808,int_stack+47144,int_stack+47084,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+90988,int_stack+47234,int_stack+47144,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+90988,int_stack+90808,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+91258,int_stack+47360,int_stack+47234,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+91258,int_stack+90988,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+90808,int_stack+72458,int_stack+81568,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+91408,int_stack+47628,int_stack+47528,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+47778,int_stack+47628,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+37534,int_stack+85058,int_stack+91408,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+47988,int_stack+47778,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+92708,int_stack+35634,int_stack+85058,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+39434,int_stack+92708,int_stack+37534,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+37534,int_stack+48328,int_stack+48268,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+37714,int_stack+48418,int_stack+48328,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+37714,int_stack+37534,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+37984,int_stack+48544,int_stack+48418,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+37984,int_stack+37714,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+37534,int_stack+72458,int_stack+81568,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+38134,int_stack+48812,int_stack+48712,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+48962,int_stack+48812,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+92708,int_stack+85058,int_stack+38134,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+93308,int_stack+49172,int_stack+48962,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+35634,int_stack+93308,int_stack+85058,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+40434,int_stack+35634,int_stack+92708,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+92708,int_stack+49956,int_stack+49896, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63212,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+92888,int_stack+50046,int_stack+49956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63272,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+92888,int_stack+92708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1650,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+93158,int_stack+50172,int_stack+50046, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63362,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+93158,int_stack+92888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1830,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+92708,int_stack+72458,int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+38134,int_stack+50440,int_stack+50340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58756,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+50590,int_stack+50440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58956,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+93308,int_stack+85058,int_stack+38134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95058,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+50800,int_stack+50590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59256,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41434,int_stack+35634,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95358,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+42334,int_stack+41434,int_stack+93308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96576,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+93308,int_stack+51584,int_stack+51524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63212, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+93488,int_stack+51674,int_stack+51584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63272, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+93488,int_stack+93308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1650, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41434,int_stack+51800,int_stack+51674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63362, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+41434,int_stack+93488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1830, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+41434,int_stack+72458,int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42034,int_stack+52068,int_stack+51968, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58756, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+52218,int_stack+52068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58956, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+93308,int_stack+85058,int_stack+42034, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95058, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+52428,int_stack+52218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59256, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+43334,int_stack+35634,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95358, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+44234,int_stack+43334,int_stack+93308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96576, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+93308,int_stack+53212,int_stack+53152, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63212, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+93488,int_stack+53302,int_stack+53212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63272, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+93488,int_stack+93308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1650, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43334,int_stack+53428,int_stack+53302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63362, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+43334,int_stack+93488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1830, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+43334,int_stack+72458,int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43934,int_stack+53696,int_stack+53596, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58756, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+53846,int_stack+53696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58956, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+93308,int_stack+85058,int_stack+43934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95058, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+54056,int_stack+53846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59256, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45234,int_stack+35634,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95358, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+46134,int_stack+45234,int_stack+93308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96576, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+93308,int_stack+54840,int_stack+54780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+93488,int_stack+54930,int_stack+54840, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+93488,int_stack+93308, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45234,int_stack+55056,int_stack+54930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+45234,int_stack+93488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+45234,int_stack+72458,int_stack+81568, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45834,int_stack+55324,int_stack+55224, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+55474,int_stack+55324, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+93308,int_stack+85058,int_stack+45834, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+55684,int_stack+55474, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+47134,int_stack+35634,int_stack+85058, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+48034,int_stack+47134,int_stack+93308, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+93308,int_stack+56468,int_stack+56408, 0.0, zero_stack, 1.0, int_stack+63212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+93488,int_stack+56558,int_stack+56468, 0.0, zero_stack, 1.0, int_stack+63272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+93488,int_stack+93308, 0.0, zero_stack, 1.0, int_stack+1650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47134,int_stack+56684,int_stack+56558, 0.0, zero_stack, 1.0, int_stack+63362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+47134,int_stack+93488, 0.0, zero_stack, 1.0, int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+47134,int_stack+72458,int_stack+81568, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47734,int_stack+56952,int_stack+56852, 0.0, zero_stack, 1.0, int_stack+58756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+57102,int_stack+56952, 0.0, zero_stack, 1.0, int_stack+58956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+93308,int_stack+85058,int_stack+47734, 0.0, zero_stack, 1.0, int_stack+95058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+57312,int_stack+57102, 0.0, zero_stack, 1.0, int_stack+59256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+49034,int_stack+35634,int_stack+85058, 0.0, zero_stack, 1.0, int_stack+95358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+49934,int_stack+49034,int_stack+93308, 0.0, zero_stack, 1.0, int_stack+96576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+93308,int_stack+58372,int_stack+58312, 1.0, int_stack+63212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+93488,int_stack+58462,int_stack+58372, 1.0, int_stack+63272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+93488,int_stack+93308, 1.0, int_stack+1650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49034,int_stack+58588,int_stack+58462, 1.0, int_stack+63362, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+49034,int_stack+93488, 1.0, int_stack+1830, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+49034,int_stack+72458,int_stack+81568, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49634,int_stack+59106,int_stack+58856, 1.0, int_stack+58756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+59466,int_stack+59106, 1.0, int_stack+58956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+93308,int_stack+85058,int_stack+49634, 1.0, int_stack+95058, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+59676,int_stack+59466, 1.0, int_stack+59256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+50934,int_stack+35634,int_stack+85058, 1.0, int_stack+95358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+51834,int_stack+50934,int_stack+93308, 1.0, int_stack+96576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+96576,int_stack+60460,int_stack+60400,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+96756,int_stack+60550,int_stack+60460,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+96756,int_stack+96576,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+93308,int_stack+60676,int_stack+60550,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+93308,int_stack+96756,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+96576,int_stack+72458,int_stack+81568,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+49634,int_stack+60944,int_stack+60844,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+61094,int_stack+60944,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+93308,int_stack+85058,int_stack+49634,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+50934,int_stack+61304,int_stack+61094,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+35634,int_stack+50934,int_stack+85058,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+52834,int_stack+35634,int_stack+93308,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+93308,int_stack+62088,int_stack+62028,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+93488,int_stack+62178,int_stack+62088,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+93488,int_stack+93308,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+35634,int_stack+62304,int_stack+62178,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+72458,int_stack+35634,int_stack+93488,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+35634,int_stack+72458,int_stack+81568,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+36234,int_stack+62572,int_stack+62472,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+62722,int_stack+62572,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+93308,int_stack+85058,int_stack+36234,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+50934,int_stack+62932,int_stack+62722,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+95058,int_stack+50934,int_stack+85058,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+53834,int_stack+95058,int_stack+93308,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+93308,int_stack+63716,int_stack+63656,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+93488,int_stack+63806,int_stack+63716,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+81568,int_stack+93488,int_stack+93308,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+95058,int_stack+63932,int_stack+63806,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+95436,int_stack+95058,int_stack+93488,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+93308,int_stack+95436,int_stack+81568,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+36234,int_stack+64200,int_stack+64100,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+85058,int_stack+64350,int_stack+64200,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+95058,int_stack+85058,int_stack+36234,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+50934,int_stack+64560,int_stack+64350,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+54834,int_stack+50934,int_stack+85058,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+55734,int_stack+54834,int_stack+95058,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+56734,int_stack+71458,int_stack+67378,100);
     Libderiv->ABCD[11] = int_stack + 56734;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+58534,int_stack+75338,int_stack+70738,100);
     Libderiv->ABCD[10] = int_stack + 58534;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+60334,int_stack+78768,int_stack+74618,100);
     Libderiv->ABCD[9] = int_stack + 60334;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+62134,int_stack+82648,int_stack+78048,100);
     Libderiv->ABCD[8] = int_stack + 62134;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+109224,int_stack+86228,int_stack+81928,100);
     Libderiv->ABCD[7] = int_stack + 109224;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+111024,int_stack+89808,int_stack+85508,100);
     Libderiv->ABCD[6] = int_stack + 111024;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+112824,int_stack+91708,int_stack+89088, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 112824;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+114624,int_stack+94058,int_stack+64840, 0.0, zero_stack, 1.0, int_stack+67978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 114624;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+93908,int_stack+98706,int_stack+95976, 1.0, int_stack+67978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 93908;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+116424,int_stack+2460,int_stack+97986,100);
     Libderiv->ABCD[155] = int_stack + 116424;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1600,int_stack+101224,int_stack+100624,100);
     Libderiv->ABCD[143] = int_stack + 1600;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+100306,int_stack+4978,int_stack+4378,100);
     Libderiv->ABCD[142] = int_stack + 100306;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+97776,int_stack+102224,int_stack+5978,100);
     Libderiv->ABCD[131] = int_stack + 97776;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4060,int_stack+103224,int_stack+3460,100);
     Libderiv->ABCD[130] = int_stack + 4060;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+102106,int_stack+104224,int_stack+97176,100);
     Libderiv->ABCD[129] = int_stack + 102106;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+118224,int_stack+105224,int_stack+8306,100);
     Libderiv->ABCD[119] = int_stack + 118224;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+103906,int_stack+11234,int_stack+10634,100);
     Libderiv->ABCD[118] = int_stack + 103906;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+120024,int_stack+12234,int_stack+8906,100);
     Libderiv->ABCD[117] = int_stack + 120024;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+10506,int_stack+9506,int_stack+13234,100);
     Libderiv->ABCD[116] = int_stack + 10506;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8096,int_stack+14434,int_stack+7496,100);
     Libderiv->ABCD[107] = int_stack + 8096;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+121824,int_stack+106224,int_stack+13834,100);
     Libderiv->ABCD[106] = int_stack + 121824;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+12306,int_stack+16784,int_stack+6578,100);
     Libderiv->ABCD[105] = int_stack + 12306;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5860,int_stack+107224,int_stack+15434,100);
     Libderiv->ABCD[104] = int_stack + 5860;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+105706,int_stack+19134,int_stack+16034,100);
     Libderiv->ABCD[103] = int_stack + 105706;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14106,int_stack+108224,int_stack+17784,100);
     Libderiv->ABCD[95] = int_stack + 14106;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15906,int_stack+20134,int_stack+68578,100);
     Libderiv->ABCD[94] = int_stack + 15906;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+17706,int_stack+21134,int_stack+70096,100);
     Libderiv->ABCD[93] = int_stack + 17706;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+19506,int_stack+73598,int_stack+72998,100);
     Libderiv->ABCD[92] = int_stack + 19506;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+71338,int_stack+65440,int_stack+76338,100);
     Libderiv->ABCD[91] = int_stack + 71338;
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+123624,int_stack+76938,int_stack+79768,100);
     Libderiv->ABCD[90] = int_stack + 123624;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+75218,int_stack+87228,int_stack+84458, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[47] = int_stack + 75218;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+86108,int_stack+22134,int_stack+80368, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[46] = int_stack + 86108;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+21306,int_stack+23134,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74618, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[45] = int_stack + 21306;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+78648,int_stack+24134,int_stack+83648, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[44] = int_stack + 78648;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+23106,int_stack+26034,int_stack+80968, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[43] = int_stack + 23106;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+82528,int_stack+27934,int_stack+99706, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[42] = int_stack + 82528;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+125424,int_stack+600,int_stack+25134, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+89088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[38] = int_stack + 125424;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+24906,int_stack+28934,int_stack+27034, 0.0, zero_stack, 1.0, int_stack+67378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[35] = int_stack + 24906;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+26706,int_stack+30834,int_stack+69178, 0.0, zero_stack, 1.0, int_stack+70738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[34] = int_stack + 26706;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+67978,int_stack+32734,int_stack+66440, 0.0, zero_stack, 1.0, int_stack+74618, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[33] = int_stack + 67978;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+65440,int_stack+34634,int_stack+29934, 0.0, zero_stack, 1.0, int_stack+78048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[32] = int_stack + 65440;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+28506,int_stack+36534,int_stack+31834, 0.0, zero_stack, 1.0, int_stack+81928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[31] = int_stack + 28506;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+30306,int_stack+38434,int_stack+33734, 0.0, zero_stack, 1.0, int_stack+85508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[30] = int_stack + 30306;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+32106,int_stack+39434,int_stack+90808, 0.0, zero_stack, 1.0, int_stack+89088, 1.0, int_stack+64840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[26] = int_stack + 32106;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+38134,int_stack+40434,int_stack+37534, 0.0, zero_stack, 2.0, int_stack+64840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[25] = int_stack + 38134;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+36234,int_stack+42334,int_stack+92708, 1.0, int_stack+67378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[23] = int_stack + 36234;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+89688,int_stack+44234,int_stack+41434, 1.0, int_stack+70738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[22] = int_stack + 89688;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+39934,int_stack+46134,int_stack+43334, 1.0, int_stack+74618, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[21] = int_stack + 39934;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+41734,int_stack+48034,int_stack+45234, 1.0, int_stack+78048, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[20] = int_stack + 41734;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+43534,int_stack+49934,int_stack+47134, 1.0, int_stack+81928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[19] = int_stack + 43534;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+49634,int_stack+51834,int_stack+49034, 1.0, int_stack+85508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[18] = int_stack + 49634;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+45334,int_stack+52834,int_stack+96576, 1.0, int_stack+89088, 0.0, zero_stack, 1.0, int_stack+95976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[14] = int_stack + 45334;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+51434,int_stack+53834,int_stack+35634, 1.0, int_stack+64840, 1.0, int_stack+95976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[13] = int_stack + 51434;
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+53234,int_stack+55734,int_stack+93308, 2.0, int_stack+95976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[12] = int_stack + 53234;

}
