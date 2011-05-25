#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_dpdd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|dd) integrals */

void d12hrr_order_dpdd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][10] = int_stack + 150;
 Libderiv->deriv_classes[3][4][9] = int_stack + 300;
 Libderiv->deriv_classes[3][4][8] = int_stack + 450;
 Libderiv->deriv_classes[3][4][7] = int_stack + 600;
 Libderiv->dvrr_classes[3][3] = int_stack + 750;
 Libderiv->deriv_classes[3][4][6] = int_stack + 850;
 Libderiv->deriv_classes[3][4][2] = int_stack + 1000;
 Libderiv->deriv_classes[3][4][1] = int_stack + 1150;
 Libderiv->dvrr_classes[2][4] = int_stack + 1300;
 Libderiv->deriv_classes[3][4][0] = int_stack + 1390;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 1540;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 1576;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 1636;
 Libderiv->deriv2_classes[3][2][143] = int_stack + 1726;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 1786;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 1886;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 2036;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 2072;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 2132;
 Libderiv->deriv2_classes[3][2][131] = int_stack + 2222;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 2282;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 2382;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 2532;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 2568;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 2628;
 Libderiv->deriv2_classes[3][2][130] = int_stack + 2718;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 2778;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 2878;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 3028;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 3064;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 3124;
 Libderiv->deriv2_classes[3][2][119] = int_stack + 3214;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 3274;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 3374;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 3524;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 3560;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 3620;
 Libderiv->deriv2_classes[3][2][118] = int_stack + 3710;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 3770;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 3870;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 4020;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 4056;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 4116;
 Libderiv->deriv2_classes[3][2][117] = int_stack + 4206;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 4266;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 4366;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 4516;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 4552;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 4612;
 Libderiv->deriv2_classes[3][2][107] = int_stack + 4702;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 4762;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 4862;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 5012;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 5048;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 5108;
 Libderiv->deriv2_classes[3][2][106] = int_stack + 5198;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 5258;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 5358;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 5508;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 5544;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 5604;
 Libderiv->deriv2_classes[3][2][105] = int_stack + 5694;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 5754;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 5854;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 6004;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 6040;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 6100;
 Libderiv->deriv2_classes[3][2][104] = int_stack + 6190;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 6250;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 6350;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 6500;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 6536;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 6596;
 Libderiv->deriv2_classes[3][2][95] = int_stack + 6686;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 6746;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 6846;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 6996;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 7032;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 7092;
 Libderiv->deriv2_classes[3][2][94] = int_stack + 7182;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 7242;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 7342;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 7492;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 7528;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 7588;
 Libderiv->deriv2_classes[3][2][93] = int_stack + 7678;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 7738;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 7838;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 7988;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 8024;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 8084;
 Libderiv->deriv2_classes[3][2][92] = int_stack + 8174;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 8234;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 8334;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 8484;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 8520;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 8580;
 Libderiv->deriv2_classes[3][2][91] = int_stack + 8670;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 8730;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 8830;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 8980;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 9016;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 9076;
 Libderiv->deriv_classes[3][2][11] = int_stack + 9166;
 Libderiv->deriv2_classes[3][2][83] = int_stack + 9226;
 Libderiv->deriv_classes[3][3][11] = int_stack + 9286;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 9386;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 9486;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 9636;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 9672;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 9732;
 Libderiv->deriv_classes[3][2][10] = int_stack + 9822;
 Libderiv->deriv2_classes[3][2][82] = int_stack + 9882;
 Libderiv->deriv_classes[3][3][10] = int_stack + 9942;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 10042;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 10142;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 10292;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 10328;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 10388;
 Libderiv->deriv_classes[3][2][9] = int_stack + 10478;
 Libderiv->deriv2_classes[3][2][81] = int_stack + 10538;
 Libderiv->deriv_classes[3][3][9] = int_stack + 10598;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 10698;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 10798;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 10948;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 10984;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 11044;
 Libderiv->deriv_classes[3][2][8] = int_stack + 11134;
 Libderiv->deriv2_classes[3][2][80] = int_stack + 11194;
 Libderiv->deriv_classes[3][3][8] = int_stack + 11254;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 11354;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 11454;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 11604;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 11640;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 11700;
 Libderiv->deriv_classes[3][2][7] = int_stack + 11790;
 Libderiv->deriv2_classes[3][2][79] = int_stack + 11850;
 Libderiv->deriv_classes[3][3][7] = int_stack + 11910;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 12010;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 12110;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 12260;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 12296;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 12356;
 Libderiv->dvrr_classes[3][2] = int_stack + 12446;
 Libderiv->deriv_classes[3][2][6] = int_stack + 12506;
 Libderiv->deriv2_classes[3][2][78] = int_stack + 12566;
 Libderiv->deriv_classes[3][3][6] = int_stack + 12626;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 12726;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 12826;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 12976;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 13012;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 13072;
 Libderiv->deriv2_classes[3][2][35] = int_stack + 13162;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 13222;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 13322;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 13472;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 13508;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 13568;
 Libderiv->deriv2_classes[3][2][34] = int_stack + 13658;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 13718;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 13818;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 13968;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 14004;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 14064;
 Libderiv->deriv2_classes[3][2][33] = int_stack + 14154;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 14214;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 14314;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 14464;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 14500;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 14560;
 Libderiv->deriv2_classes[3][2][32] = int_stack + 14650;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 14710;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 14810;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 14960;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 14996;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 15056;
 Libderiv->deriv2_classes[3][2][31] = int_stack + 15146;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 15206;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 15306;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 15456;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 15492;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 15552;
 Libderiv->deriv_classes[3][2][2] = int_stack + 15642;
 Libderiv->deriv2_classes[3][2][30] = int_stack + 15702;
 Libderiv->deriv_classes[3][3][2] = int_stack + 15762;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 15862;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 15962;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 16112;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 16148;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 16208;
 Libderiv->deriv2_classes[3][2][26] = int_stack + 16298;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 16358;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 16458;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 16608;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 16644;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 16704;
 Libderiv->deriv2_classes[3][2][23] = int_stack + 16794;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 16854;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 16954;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 17104;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 17140;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 17200;
 Libderiv->deriv2_classes[3][2][22] = int_stack + 17290;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 17350;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 17450;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 17600;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 17636;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 17696;
 Libderiv->deriv2_classes[3][2][21] = int_stack + 17786;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 17846;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 17946;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 18096;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 18132;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 18192;
 Libderiv->deriv2_classes[3][2][20] = int_stack + 18282;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 18342;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 18442;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 18592;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 18628;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 18688;
 Libderiv->deriv2_classes[3][2][19] = int_stack + 18778;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 18838;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 18938;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 19088;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 19124;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 19184;
 Libderiv->deriv_classes[3][2][1] = int_stack + 19274;
 Libderiv->deriv2_classes[3][2][18] = int_stack + 19334;
 Libderiv->deriv_classes[3][3][1] = int_stack + 19394;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 19494;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 19594;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 19744;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 19780;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 19840;
 Libderiv->deriv2_classes[3][2][14] = int_stack + 19930;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 19990;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 20090;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 20240;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 20276;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 20336;
 Libderiv->deriv2_classes[3][2][13] = int_stack + 20426;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 20486;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 20586;
 Libderiv->deriv_classes[2][2][11] = int_stack + 20736;
 Libderiv->deriv_classes[2][3][11] = int_stack + 20772;
 Libderiv->deriv_classes[2][4][11] = int_stack + 20832;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 20922;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 20958;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 21018;
 Libderiv->deriv2_classes[3][2][11] = int_stack + 21108;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 21168;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 21268;
 Libderiv->deriv_classes[2][2][10] = int_stack + 21418;
 Libderiv->deriv_classes[2][3][10] = int_stack + 21454;
 Libderiv->deriv_classes[2][4][10] = int_stack + 21514;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 21604;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 21640;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 21700;
 Libderiv->deriv2_classes[3][2][10] = int_stack + 21790;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 21850;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 21950;
 Libderiv->deriv_classes[2][2][9] = int_stack + 22100;
 Libderiv->deriv_classes[2][3][9] = int_stack + 22136;
 Libderiv->deriv_classes[2][4][9] = int_stack + 22196;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 22286;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 22322;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 22382;
 Libderiv->deriv2_classes[3][2][9] = int_stack + 22472;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 22532;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 22632;
 Libderiv->deriv_classes[2][2][8] = int_stack + 22782;
 Libderiv->deriv_classes[2][3][8] = int_stack + 22818;
 Libderiv->deriv_classes[2][4][8] = int_stack + 22878;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 22968;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 23004;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 23064;
 Libderiv->deriv2_classes[3][2][8] = int_stack + 23154;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 23214;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 23314;
 Libderiv->deriv_classes[2][2][7] = int_stack + 23464;
 Libderiv->deriv_classes[2][3][7] = int_stack + 23500;
 Libderiv->deriv_classes[2][4][7] = int_stack + 23560;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 23650;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 23686;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 23746;
 Libderiv->deriv2_classes[3][2][7] = int_stack + 23836;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 23896;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 23996;
 Libderiv->dvrr_classes[2][2] = int_stack + 24146;
 Libderiv->deriv_classes[2][2][6] = int_stack + 24182;
 Libderiv->dvrr_classes[2][3] = int_stack + 24218;
 Libderiv->deriv_classes[2][3][6] = int_stack + 24278;
 Libderiv->deriv_classes[2][4][6] = int_stack + 24338;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 24428;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 24464;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 24524;
 Libderiv->deriv_classes[3][2][0] = int_stack + 24614;
 Libderiv->deriv2_classes[3][2][6] = int_stack + 24674;
 Libderiv->deriv_classes[3][3][0] = int_stack + 24734;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 24834;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 24934;
 Libderiv->deriv_classes[2][2][2] = int_stack + 25084;
 Libderiv->deriv_classes[2][3][2] = int_stack + 25120;
 Libderiv->deriv_classes[2][4][2] = int_stack + 25180;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 25270;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 25306;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 25366;
 Libderiv->deriv2_classes[3][2][2] = int_stack + 25456;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 25516;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 25616;
 Libderiv->deriv_classes[2][2][1] = int_stack + 25766;
 Libderiv->deriv_classes[2][3][1] = int_stack + 25802;
 Libderiv->deriv_classes[2][4][1] = int_stack + 25862;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 25952;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 25988;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 26048;
 Libderiv->deriv2_classes[3][2][1] = int_stack + 26138;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 26198;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 26298;
 Libderiv->deriv_classes[2][2][0] = int_stack + 26448;
 Libderiv->deriv_classes[2][3][0] = int_stack + 26484;
 Libderiv->deriv_classes[2][4][0] = int_stack + 26544;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 26634;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 26670;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 26730;
 Libderiv->deriv2_classes[3][2][0] = int_stack + 26820;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 26880;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 26980;
 memset(int_stack,0,217040);

 Libderiv->dvrr_stack = int_stack + 42718;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_dpdd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27130,int_stack+24218,int_stack+24146,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+20772,int_stack+20736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24146,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+20832,int_stack+20772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24218,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+27526,int_stack+27346,int_stack+27238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27130,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27346,int_stack+750,int_stack+12446,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+9286,int_stack+9166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12446,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27922,int_stack+0,int_stack+9286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+28222,int_stack+27922,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27346,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27922,int_stack+21454,int_stack+21418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24146, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28030,int_stack+21514,int_stack+21454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24218, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+28582,int_stack+28030,int_stack+27922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27130, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+28030,int_stack+9942,int_stack+9822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12446, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28798,int_stack+150,int_stack+9942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+29098,int_stack+28798,int_stack+28030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27346, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+28798,int_stack+22136,int_stack+22100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24146, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28906,int_stack+22196,int_stack+22136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24218, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+28906,int_stack+28798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27130, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+28906,int_stack+10598,int_stack+10478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12446, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29458,int_stack+300,int_stack+10598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+29758,int_stack+29458,int_stack+28906, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27346, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+29458,int_stack+22818,int_stack+22782, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29566,int_stack+22878,int_stack+22818, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+216,int_stack+29566,int_stack+29458, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+29566,int_stack+11254,int_stack+11134, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30118,int_stack+450,int_stack+11254, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+30418,int_stack+30118,int_stack+29566, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+30118,int_stack+23500,int_stack+23464, 0.0, zero_stack, 1.0, int_stack+24146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30226,int_stack+23560,int_stack+23500, 0.0, zero_stack, 1.0, int_stack+24218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+30778,int_stack+30226,int_stack+30118, 0.0, zero_stack, 1.0, int_stack+27130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+30226,int_stack+11910,int_stack+11790, 0.0, zero_stack, 1.0, int_stack+12446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30994,int_stack+600,int_stack+11910, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+31294,int_stack+30994,int_stack+30226, 0.0, zero_stack, 1.0, int_stack+27346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+30994,int_stack+24278,int_stack+24182, 1.0, int_stack+24146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+31102,int_stack+24338,int_stack+24278, 1.0, int_stack+24218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+432,int_stack+31102,int_stack+30994, 1.0, int_stack+27130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+31102,int_stack+12626,int_stack+12506, 1.0, int_stack+12446, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+31654,int_stack+850,int_stack+12626, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+31954,int_stack+31654,int_stack+31102, 1.0, int_stack+27346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+1300,int_stack+24218,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+31654,int_stack+27346,int_stack+27130,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27130,int_stack+25120,int_stack+25084,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+25180,int_stack+25120,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+648,int_stack+27346,int_stack+27130,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27346,int_stack+15762,int_stack+15642,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+32314,int_stack+1000,int_stack+15762,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+32614,int_stack+32314,int_stack+27346,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+32314,int_stack+25802,int_stack+25766,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+32422,int_stack+25862,int_stack+25802,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+864,int_stack+32422,int_stack+32314,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+32422,int_stack+19394,int_stack+19274,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+32974,int_stack+1150,int_stack+19394,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+33274,int_stack+32974,int_stack+32422,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+32974,int_stack+26484,int_stack+26448,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+33082,int_stack+26544,int_stack+26484,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1080,int_stack+33082,int_stack+32974,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+33082,int_stack+24734,int_stack+24614,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+33634,int_stack+1390,int_stack+24734,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+33934,int_stack+33634,int_stack+33082,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+1576,int_stack+1540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+20736,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+1636,int_stack+1576, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+20772,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1296,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27238,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+1786,int_stack+1726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9166,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+34294,int_stack+1886,int_stack+1786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9286,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1512,int_stack+34294,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27742,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+2072,int_stack+2036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20736, 1.0, int_stack+21418,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+2132,int_stack+2072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20772, 1.0, int_stack+21454,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+34294,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27238, 1.0, int_stack+27922,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+2282,int_stack+2222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9166, 1.0, int_stack+9822,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+34510,int_stack+2382,int_stack+2282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9286, 1.0, int_stack+9942,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1872,int_stack+34510,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27742, 1.0, int_stack+28030,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+2568,int_stack+2532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+21418, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+2628,int_stack+2568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+21454, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+34510,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27922, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+2778,int_stack+2718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9822, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2232,int_stack+2878,int_stack+2778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9942, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2532,int_stack+2232,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28030, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+3064,int_stack+3028, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20736, 0.0, zero_stack, 1.0, int_stack+22100,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+3124,int_stack+3064, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20772, 0.0, zero_stack, 1.0, int_stack+22136,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2232,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27238, 0.0, zero_stack, 1.0, int_stack+28798,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+3274,int_stack+3214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9166, 0.0, zero_stack, 1.0, int_stack+10478,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2892,int_stack+3374,int_stack+3274, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9286, 0.0, zero_stack, 1.0, int_stack+10598,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+34726,int_stack+2892,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27742, 0.0, zero_stack, 1.0, int_stack+28906,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+3560,int_stack+3524, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21418, 1.0, int_stack+22100, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+3620,int_stack+3560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21454, 1.0, int_stack+22136, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2892,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27922, 1.0, int_stack+28798, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+3770,int_stack+3710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9822, 1.0, int_stack+10478, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3108,int_stack+3870,int_stack+3770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9942, 1.0, int_stack+10598, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3408,int_stack+3108,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28030, 1.0, int_stack+28906, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+4056,int_stack+4020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22100, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+4116,int_stack+4056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22136, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3108,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28798, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+4266,int_stack+4206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10478, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3768,int_stack+4366,int_stack+4266, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+10598, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4068,int_stack+3768,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28906, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+4552,int_stack+4516, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20736, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22782,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+4612,int_stack+4552, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20772, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22818,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3768,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27238, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29458,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+4762,int_stack+4702, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9166, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11134,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4428,int_stack+4862,int_stack+4762, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9286, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11254,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+35086,int_stack+4428,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29566,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+5048,int_stack+5012, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21418, 0.0, zero_stack, 1.0, int_stack+22782, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+5108,int_stack+5048, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21454, 0.0, zero_stack, 1.0, int_stack+22818, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4428,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27922, 0.0, zero_stack, 1.0, int_stack+29458, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+5258,int_stack+5198, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9822, 0.0, zero_stack, 1.0, int_stack+11134, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4644,int_stack+5358,int_stack+5258, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9942, 0.0, zero_stack, 1.0, int_stack+11254, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4944,int_stack+4644,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28030, 0.0, zero_stack, 1.0, int_stack+29566, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+5544,int_stack+5508, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22100, 1.0, int_stack+22782, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+5604,int_stack+5544, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22136, 1.0, int_stack+22818, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4644,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28798, 1.0, int_stack+29458, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+5754,int_stack+5694, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10478, 1.0, int_stack+11134, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5304,int_stack+5854,int_stack+5754, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10598, 1.0, int_stack+11254, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5604,int_stack+5304,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28906, 1.0, int_stack+29566, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+6040,int_stack+6004, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+6100,int_stack+6040, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22818, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5304,int_stack+33742,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+6250,int_stack+6190, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+35446,int_stack+6350,int_stack+6250, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5964,int_stack+35446,int_stack+33634, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+6536,int_stack+6500, 0.0, zero_stack, 1.0, int_stack+20736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23464,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+6596,int_stack+6536, 0.0, zero_stack, 1.0, int_stack+20772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23500,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+35446,int_stack+33742,int_stack+33634, 0.0, zero_stack, 1.0, int_stack+27238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30118,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+6746,int_stack+6686, 0.0, zero_stack, 1.0, int_stack+9166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11790,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+35662,int_stack+6846,int_stack+6746, 0.0, zero_stack, 1.0, int_stack+9286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11910,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6324,int_stack+35662,int_stack+33634, 0.0, zero_stack, 1.0, int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30226,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+7032,int_stack+6996, 0.0, zero_stack, 1.0, int_stack+21418, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23464, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+7092,int_stack+7032, 0.0, zero_stack, 1.0, int_stack+21454, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23500, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+35662,int_stack+33742,int_stack+33634, 0.0, zero_stack, 1.0, int_stack+27922, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30118, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+7242,int_stack+7182, 0.0, zero_stack, 1.0, int_stack+9822, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11790, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6684,int_stack+7342,int_stack+7242, 0.0, zero_stack, 1.0, int_stack+9942, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11910, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6984,int_stack+6684,int_stack+33634, 0.0, zero_stack, 1.0, int_stack+28030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30226, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+7528,int_stack+7492, 0.0, zero_stack, 1.0, int_stack+22100, 0.0, zero_stack, 1.0, int_stack+23464, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+7588,int_stack+7528, 0.0, zero_stack, 1.0, int_stack+22136, 0.0, zero_stack, 1.0, int_stack+23500, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6684,int_stack+33742,int_stack+33634, 0.0, zero_stack, 1.0, int_stack+28798, 0.0, zero_stack, 1.0, int_stack+30118, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+7738,int_stack+7678, 0.0, zero_stack, 1.0, int_stack+10478, 0.0, zero_stack, 1.0, int_stack+11790, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7344,int_stack+7838,int_stack+7738, 0.0, zero_stack, 1.0, int_stack+10598, 0.0, zero_stack, 1.0, int_stack+11910, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+35878,int_stack+7344,int_stack+33634, 0.0, zero_stack, 1.0, int_stack+28906, 0.0, zero_stack, 1.0, int_stack+30226, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+8024,int_stack+7988, 0.0, zero_stack, 1.0, int_stack+22782, 1.0, int_stack+23464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+8084,int_stack+8024, 0.0, zero_stack, 1.0, int_stack+22818, 1.0, int_stack+23500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7344,int_stack+33742,int_stack+33634, 0.0, zero_stack, 1.0, int_stack+29458, 1.0, int_stack+30118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+8234,int_stack+8174, 0.0, zero_stack, 1.0, int_stack+11134, 1.0, int_stack+11790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7560,int_stack+8334,int_stack+8234, 0.0, zero_stack, 1.0, int_stack+11254, 1.0, int_stack+11910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7860,int_stack+7560,int_stack+33634, 0.0, zero_stack, 1.0, int_stack+29566, 1.0, int_stack+30226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+8520,int_stack+8484, 0.0, zero_stack, 2.0, int_stack+23464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+8580,int_stack+8520, 0.0, zero_stack, 2.0, int_stack+23500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7560,int_stack+33742,int_stack+33634, 0.0, zero_stack, 2.0, int_stack+30118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+8730,int_stack+8670, 0.0, zero_stack, 2.0, int_stack+11790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8220,int_stack+8830,int_stack+8730, 0.0, zero_stack, 2.0, int_stack+11910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8520,int_stack+8220,int_stack+33634, 0.0, zero_stack, 2.0, int_stack+30226, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+9016,int_stack+8980, 1.0, int_stack+20736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24182,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33742,int_stack+9076,int_stack+9016, 1.0, int_stack+20772, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24278,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8220,int_stack+33742,int_stack+33634, 1.0, int_stack+27238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30994,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+33634,int_stack+9386,int_stack+9226, 1.0, int_stack+9166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12506,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8880,int_stack+9486,int_stack+9386, 1.0, int_stack+9286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12626,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9180,int_stack+8880,int_stack+33634, 1.0, int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31102,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+9672,int_stack+9636, 1.0, int_stack+21418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24182, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+9732,int_stack+9672, 1.0, int_stack+21454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24278, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+33634,int_stack+27742,int_stack+27238, 1.0, int_stack+27922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30994, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+10042,int_stack+9882, 1.0, int_stack+9822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12506, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8880,int_stack+10142,int_stack+10042, 1.0, int_stack+9942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12626, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9540,int_stack+8880,int_stack+27742, 1.0, int_stack+28030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31102, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+10328,int_stack+10292, 1.0, int_stack+22100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24182, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+10388,int_stack+10328, 1.0, int_stack+22136, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24278, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+27922,int_stack+27742,int_stack+27238, 1.0, int_stack+28798, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30994, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+10698,int_stack+10538, 1.0, int_stack+10478, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12506, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8880,int_stack+10798,int_stack+10698, 1.0, int_stack+10598, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12626, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9900,int_stack+8880,int_stack+27742, 1.0, int_stack+28906, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31102, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+10984,int_stack+10948, 1.0, int_stack+22782, 0.0, zero_stack, 1.0, int_stack+24182, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+11044,int_stack+10984, 1.0, int_stack+22818, 0.0, zero_stack, 1.0, int_stack+24278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8880,int_stack+27742,int_stack+27238, 1.0, int_stack+29458, 0.0, zero_stack, 1.0, int_stack+30994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+11354,int_stack+11194, 1.0, int_stack+11134, 0.0, zero_stack, 1.0, int_stack+12506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28798,int_stack+11454,int_stack+11354, 1.0, int_stack+11254, 0.0, zero_stack, 1.0, int_stack+12626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10260,int_stack+28798,int_stack+27742, 1.0, int_stack+29566, 0.0, zero_stack, 1.0, int_stack+31102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+11640,int_stack+11604, 1.0, int_stack+23464, 1.0, int_stack+24182, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+11700,int_stack+11640, 1.0, int_stack+23500, 1.0, int_stack+24278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+28798,int_stack+27742,int_stack+27238, 1.0, int_stack+30118, 1.0, int_stack+30994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+12010,int_stack+11850, 1.0, int_stack+11790, 1.0, int_stack+12506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29458,int_stack+12110,int_stack+12010, 1.0, int_stack+11910, 1.0, int_stack+12626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10620,int_stack+29458,int_stack+27742, 1.0, int_stack+30226, 1.0, int_stack+31102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+12296,int_stack+12260, 2.0, int_stack+24182, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+12356,int_stack+12296, 2.0, int_stack+24278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+29458,int_stack+27742,int_stack+27238, 2.0, int_stack+30994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+12726,int_stack+12566, 2.0, int_stack+12506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30118,int_stack+12826,int_stack+12726, 2.0, int_stack+12626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10980,int_stack+30118,int_stack+27742, 2.0, int_stack+31102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+13012,int_stack+12976, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25084,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+13072,int_stack+13012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25120,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+30118,int_stack+27742,int_stack+27238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27130,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+13222,int_stack+13162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15642,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30994,int_stack+13322,int_stack+13222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15762,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11340,int_stack+30994,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27346,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+13508,int_stack+13472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25084, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+13568,int_stack+13508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25120, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+30994,int_stack+27742,int_stack+27238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27130, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+13718,int_stack+13658, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15642, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11700,int_stack+13818,int_stack+13718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15762, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12000,int_stack+11700,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27346, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+14004,int_stack+13968, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25084, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+14064,int_stack+14004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25120, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11700,int_stack+27742,int_stack+27238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27130, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+14214,int_stack+14154, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15642, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12360,int_stack+14314,int_stack+14214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15762, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12660,int_stack+12360,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27346, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+14500,int_stack+14464, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+14560,int_stack+14500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12360,int_stack+27742,int_stack+27238, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+14710,int_stack+14650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13020,int_stack+14810,int_stack+14710, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+13320,int_stack+13020,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+14996,int_stack+14960, 0.0, zero_stack, 1.0, int_stack+25084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+15056,int_stack+14996, 0.0, zero_stack, 1.0, int_stack+25120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+13020,int_stack+27742,int_stack+27238, 0.0, zero_stack, 1.0, int_stack+27130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+15206,int_stack+15146, 0.0, zero_stack, 1.0, int_stack+15642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13680,int_stack+15306,int_stack+15206, 0.0, zero_stack, 1.0, int_stack+15762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+13980,int_stack+13680,int_stack+27742, 0.0, zero_stack, 1.0, int_stack+27346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27238,int_stack+15492,int_stack+15456, 1.0, int_stack+25084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27742,int_stack+15552,int_stack+15492, 1.0, int_stack+25120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+13680,int_stack+27742,int_stack+27238, 1.0, int_stack+27130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+15862,int_stack+15702, 1.0, int_stack+15642, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14340,int_stack+15962,int_stack+15862, 1.0, int_stack+15762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+14640,int_stack+14340,int_stack+27742, 1.0, int_stack+27346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+16148,int_stack+16112,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14340,int_stack+16208,int_stack+16148,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+27130,int_stack+14340,int_stack+27742,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+16358,int_stack+16298,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14340,int_stack+16458,int_stack+16358,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+15000,int_stack+14340,int_stack+27742,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+16644,int_stack+16608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25766,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+16704,int_stack+16644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25802,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+14340,int_stack+27346,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32314,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+16854,int_stack+16794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19274,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15360,int_stack+16954,int_stack+16854, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19394,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15660,int_stack+15360,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32422,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+17140,int_stack+17104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25766, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+17200,int_stack+17140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25802, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15360,int_stack+27346,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32314, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+17350,int_stack+17290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19274, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16020,int_stack+17450,int_stack+17350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19394, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16320,int_stack+16020,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32422, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+17636,int_stack+17600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25766, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+17696,int_stack+17636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25802, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16020,int_stack+27346,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32314, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+17846,int_stack+17786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19274, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16680,int_stack+17946,int_stack+17846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19394, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16980,int_stack+16680,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32422, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+18132,int_stack+18096, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+18192,int_stack+18132, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16680,int_stack+27346,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+18342,int_stack+18282, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19274, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17340,int_stack+18442,int_stack+18342, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17640,int_stack+17340,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+18628,int_stack+18592, 0.0, zero_stack, 1.0, int_stack+25766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+18688,int_stack+18628, 0.0, zero_stack, 1.0, int_stack+25802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17340,int_stack+27346,int_stack+27742, 0.0, zero_stack, 1.0, int_stack+32314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+18838,int_stack+18778, 0.0, zero_stack, 1.0, int_stack+19274, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18000,int_stack+18938,int_stack+18838, 0.0, zero_stack, 1.0, int_stack+19394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+18300,int_stack+18000,int_stack+27742, 0.0, zero_stack, 1.0, int_stack+32422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+19124,int_stack+19088, 1.0, int_stack+25766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+19184,int_stack+19124, 1.0, int_stack+25802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+18000,int_stack+27346,int_stack+27742, 1.0, int_stack+32314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+19494,int_stack+19334, 1.0, int_stack+19274, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18660,int_stack+19594,int_stack+19494, 1.0, int_stack+19394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+18960,int_stack+18660,int_stack+27742, 1.0, int_stack+32422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+19780,int_stack+19744,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+19840,int_stack+19780,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+18660,int_stack+27346,int_stack+27742,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+19990,int_stack+19930,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+32314,int_stack+20090,int_stack+19990,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+19320,int_stack+32314,int_stack+27742,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+20276,int_stack+20240,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+20336,int_stack+20276,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+32314,int_stack+27346,int_stack+27742,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+20486,int_stack+20426,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19680,int_stack+20586,int_stack+20486,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+19980,int_stack+19680,int_stack+27742,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+20958,int_stack+20922, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26448,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+21018,int_stack+20958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26484,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19680,int_stack+27346,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32974,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+21168,int_stack+21108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24614,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20340,int_stack+21268,int_stack+21168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24734,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20640,int_stack+20340,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33082,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+21640,int_stack+21604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26448, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+21700,int_stack+21640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26484, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20340,int_stack+27346,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32974, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+21850,int_stack+21790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24614, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21000,int_stack+21950,int_stack+21850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24734, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21300,int_stack+21000,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33082, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+22322,int_stack+22286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26448, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+22382,int_stack+22322, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26484, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21000,int_stack+27346,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32974, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+22532,int_stack+22472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24614, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21660,int_stack+22632,int_stack+22532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24734, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21960,int_stack+21660,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33082, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+23004,int_stack+22968, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+23064,int_stack+23004, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21660,int_stack+27346,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32974, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+23214,int_stack+23154, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22320,int_stack+23314,int_stack+23214, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+22620,int_stack+22320,int_stack+27742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+23686,int_stack+23650, 0.0, zero_stack, 1.0, int_stack+26448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+23746,int_stack+23686, 0.0, zero_stack, 1.0, int_stack+26484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+22320,int_stack+27346,int_stack+27742, 0.0, zero_stack, 1.0, int_stack+32974, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+23896,int_stack+23836, 0.0, zero_stack, 1.0, int_stack+24614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22980,int_stack+23996,int_stack+23896, 0.0, zero_stack, 1.0, int_stack+24734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+23280,int_stack+22980,int_stack+27742, 0.0, zero_stack, 1.0, int_stack+33082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+24464,int_stack+24428, 1.0, int_stack+26448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+24524,int_stack+24464, 1.0, int_stack+26484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+22980,int_stack+27346,int_stack+27742, 1.0, int_stack+32974, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+24834,int_stack+24674, 1.0, int_stack+24614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23640,int_stack+24934,int_stack+24834, 1.0, int_stack+24734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+23940,int_stack+23640,int_stack+27742, 1.0, int_stack+33082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+25306,int_stack+25270,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+25366,int_stack+25306,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+23640,int_stack+27346,int_stack+27742,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+25516,int_stack+25456,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+32974,int_stack+25616,int_stack+25516,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+24300,int_stack+32974,int_stack+27742,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+25988,int_stack+25952,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+26048,int_stack+25988,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+32974,int_stack+27346,int_stack+27742,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+26198,int_stack+26138,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24660,int_stack+26298,int_stack+26198,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+24960,int_stack+24660,int_stack+27742,10);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+26670,int_stack+26634,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27346,int_stack+26730,int_stack+26670,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+24660,int_stack+27346,int_stack+27742,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+27742,int_stack+26880,int_stack+26820,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+25320,int_stack+26980,int_stack+26880,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+25620,int_stack+25320,int_stack+27742,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+25980,int_stack+28222,int_stack+27526,36);
     Libderiv->ABCD[11] = int_stack + 25980;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+36238,int_stack+29098,int_stack+28582,36);
     Libderiv->ABCD[10] = int_stack + 36238;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+36886,int_stack+29758,int_stack+0,36);
     Libderiv->ABCD[9] = int_stack + 36886;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+37534,int_stack+30418,int_stack+216,36);
     Libderiv->ABCD[8] = int_stack + 37534;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+38182,int_stack+31294,int_stack+30778,36);
     Libderiv->ABCD[7] = int_stack + 38182;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+38830,int_stack+31954,int_stack+432,36);
     Libderiv->ABCD[6] = int_stack + 38830;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+39478,int_stack+32614,int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[2] = int_stack + 39478;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+40126,int_stack+33274,int_stack+864, 0.0, zero_stack, 1.0, int_stack+31654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[1] = int_stack + 40126;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+40774,int_stack+33934,int_stack+1080, 1.0, int_stack+31654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[0] = int_stack + 40774;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+31210,int_stack+1512,int_stack+1296,36);
     Libderiv->ABCD[155] = int_stack + 31210;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+41422,int_stack+1872,int_stack+34294,36);
     Libderiv->ABCD[143] = int_stack + 41422;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1296,int_stack+2532,int_stack+34510,36);
     Libderiv->ABCD[142] = int_stack + 1296;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+33850,int_stack+34726,int_stack+2232,36);
     Libderiv->ABCD[131] = int_stack + 33850;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1944,int_stack+3408,int_stack+2892,36);
     Libderiv->ABCD[130] = int_stack + 1944;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+42070,int_stack+4068,int_stack+3108,36);
     Libderiv->ABCD[129] = int_stack + 42070;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2592,int_stack+35086,int_stack+3768,36);
     Libderiv->ABCD[119] = int_stack + 2592;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3240,int_stack+4944,int_stack+4428,36);
     Libderiv->ABCD[118] = int_stack + 3240;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3888,int_stack+5604,int_stack+4644,36);
     Libderiv->ABCD[117] = int_stack + 3888;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4536,int_stack+5964,int_stack+5304,36);
     Libderiv->ABCD[116] = int_stack + 4536;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5184,int_stack+6324,int_stack+35446,36);
     Libderiv->ABCD[107] = int_stack + 5184;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5832,int_stack+6984,int_stack+35662,36);
     Libderiv->ABCD[106] = int_stack + 5832;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+34498,int_stack+35878,int_stack+6684,36);
     Libderiv->ABCD[105] = int_stack + 34498;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6480,int_stack+7860,int_stack+7344,36);
     Libderiv->ABCD[104] = int_stack + 6480;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+35146,int_stack+8520,int_stack+7560,36);
     Libderiv->ABCD[103] = int_stack + 35146;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7128,int_stack+9180,int_stack+8220,36);
     Libderiv->ABCD[95] = int_stack + 7128;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7776,int_stack+9540,int_stack+33634,36);
     Libderiv->ABCD[94] = int_stack + 7776;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+33190,int_stack+9900,int_stack+27922,36);
     Libderiv->ABCD[93] = int_stack + 33190;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+27742,int_stack+10260,int_stack+8880,36);
     Libderiv->ABCD[92] = int_stack + 27742;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8424,int_stack+10620,int_stack+28798,36);
     Libderiv->ABCD[91] = int_stack + 8424;
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+28798,int_stack+10980,int_stack+29458,36);
     Libderiv->ABCD[90] = int_stack + 28798;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+29446,int_stack+11340,int_stack+30118, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[47] = int_stack + 29446;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+30094,int_stack+12000,int_stack+30994, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[46] = int_stack + 30094;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9072,int_stack+12660,int_stack+11700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[45] = int_stack + 9072;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9720,int_stack+13320,int_stack+12360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[44] = int_stack + 9720;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+10368,int_stack+13980,int_stack+13020, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[43] = int_stack + 10368;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+11016,int_stack+14640,int_stack+13680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[42] = int_stack + 11016;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+11664,int_stack+15000,int_stack+27130, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[38] = int_stack + 11664;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+12312,int_stack+15660,int_stack+14340, 0.0, zero_stack, 1.0, int_stack+27526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[35] = int_stack + 12312;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+12960,int_stack+16320,int_stack+15360, 0.0, zero_stack, 1.0, int_stack+28582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[34] = int_stack + 12960;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+13608,int_stack+16980,int_stack+16020, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[33] = int_stack + 13608;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+14256,int_stack+17640,int_stack+16680, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[32] = int_stack + 14256;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+14904,int_stack+18300,int_stack+17340, 0.0, zero_stack, 1.0, int_stack+30778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[31] = int_stack + 14904;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+15552,int_stack+18960,int_stack+18000, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[30] = int_stack + 15552;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+16200,int_stack+19320,int_stack+18660, 0.0, zero_stack, 1.0, int_stack+648, 1.0, int_stack+864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[26] = int_stack + 16200;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+16848,int_stack+19980,int_stack+32314, 0.0, zero_stack, 2.0, int_stack+864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[25] = int_stack + 16848;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+17496,int_stack+20640,int_stack+19680, 1.0, int_stack+27526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[23] = int_stack + 17496;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+18144,int_stack+21300,int_stack+20340, 1.0, int_stack+28582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[22] = int_stack + 18144;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+18792,int_stack+21960,int_stack+21000, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[21] = int_stack + 18792;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+19440,int_stack+22620,int_stack+21660, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[20] = int_stack + 19440;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+20088,int_stack+23280,int_stack+22320, 1.0, int_stack+30778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[19] = int_stack + 20088;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+20736,int_stack+23940,int_stack+22980, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[18] = int_stack + 20736;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+24300,int_stack+23640, 1.0, int_stack+648, 0.0, zero_stack, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[14] = int_stack + 0;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+21384,int_stack+24960,int_stack+32974, 1.0, int_stack+864, 1.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[13] = int_stack + 21384;
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+22032,int_stack+25620,int_stack+24660, 2.0, int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[12] = int_stack + 22032;

}
