#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ppdd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|dd) integrals */

void d12hrr_order_ppdd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][10] = int_stack + 90;
 Libderiv->deriv_classes[2][4][9] = int_stack + 180;
 Libderiv->deriv_classes[2][4][8] = int_stack + 270;
 Libderiv->deriv_classes[2][4][7] = int_stack + 360;
 Libderiv->dvrr_classes[2][3] = int_stack + 450;
 Libderiv->deriv_classes[2][4][6] = int_stack + 510;
 Libderiv->deriv_classes[2][4][2] = int_stack + 600;
 Libderiv->deriv_classes[2][4][1] = int_stack + 690;
 Libderiv->dvrr_classes[1][4] = int_stack + 780;
 Libderiv->deriv_classes[2][4][0] = int_stack + 825;
 Libderiv->deriv2_classes[1][2][143] = int_stack + 915;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 933;
 Libderiv->deriv2_classes[1][4][143] = int_stack + 963;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 1008;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 1044;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 1104;
 Libderiv->deriv2_classes[1][2][131] = int_stack + 1194;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 1212;
 Libderiv->deriv2_classes[1][4][131] = int_stack + 1242;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 1287;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1323;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 1383;
 Libderiv->deriv2_classes[1][2][130] = int_stack + 1473;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 1491;
 Libderiv->deriv2_classes[1][4][130] = int_stack + 1521;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 1566;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 1602;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 1662;
 Libderiv->deriv2_classes[1][2][119] = int_stack + 1752;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 1770;
 Libderiv->deriv2_classes[1][4][119] = int_stack + 1800;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 1845;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 1881;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 1941;
 Libderiv->deriv2_classes[1][2][118] = int_stack + 2031;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 2049;
 Libderiv->deriv2_classes[1][4][118] = int_stack + 2079;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 2124;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 2160;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 2220;
 Libderiv->deriv2_classes[1][2][117] = int_stack + 2310;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 2328;
 Libderiv->deriv2_classes[1][4][117] = int_stack + 2358;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 2403;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 2439;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 2499;
 Libderiv->deriv2_classes[1][2][107] = int_stack + 2589;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 2607;
 Libderiv->deriv2_classes[1][4][107] = int_stack + 2637;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 2682;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 2718;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 2778;
 Libderiv->deriv2_classes[1][2][106] = int_stack + 2868;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 2886;
 Libderiv->deriv2_classes[1][4][106] = int_stack + 2916;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 2961;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 2997;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 3057;
 Libderiv->deriv2_classes[1][2][105] = int_stack + 3147;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 3165;
 Libderiv->deriv2_classes[1][4][105] = int_stack + 3195;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 3240;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 3276;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 3336;
 Libderiv->deriv2_classes[1][2][104] = int_stack + 3426;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 3444;
 Libderiv->deriv2_classes[1][4][104] = int_stack + 3474;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 3519;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 3555;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 3615;
 Libderiv->deriv2_classes[1][2][95] = int_stack + 3705;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 3723;
 Libderiv->deriv2_classes[1][4][95] = int_stack + 3753;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 3798;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 3834;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 3894;
 Libderiv->deriv2_classes[1][2][94] = int_stack + 3984;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 4002;
 Libderiv->deriv2_classes[1][4][94] = int_stack + 4032;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 4077;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 4113;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 4173;
 Libderiv->deriv2_classes[1][2][93] = int_stack + 4263;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 4281;
 Libderiv->deriv2_classes[1][4][93] = int_stack + 4311;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 4356;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 4392;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 4452;
 Libderiv->deriv2_classes[1][2][92] = int_stack + 4542;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 4560;
 Libderiv->deriv2_classes[1][4][92] = int_stack + 4590;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 4635;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 4671;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 4731;
 Libderiv->deriv2_classes[1][2][91] = int_stack + 4821;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 4839;
 Libderiv->deriv2_classes[1][4][91] = int_stack + 4869;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 4914;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 4950;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 5010;
 Libderiv->deriv2_classes[1][2][83] = int_stack + 5100;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 5118;
 Libderiv->deriv2_classes[1][4][83] = int_stack + 5148;
 Libderiv->deriv_classes[2][2][11] = int_stack + 5193;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 5229;
 Libderiv->deriv_classes[2][3][11] = int_stack + 5265;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 5325;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 5385;
 Libderiv->deriv2_classes[1][2][82] = int_stack + 5475;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 5493;
 Libderiv->deriv2_classes[1][4][82] = int_stack + 5523;
 Libderiv->deriv_classes[2][2][10] = int_stack + 5568;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 5604;
 Libderiv->deriv_classes[2][3][10] = int_stack + 5640;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 5700;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 5760;
 Libderiv->deriv2_classes[1][2][81] = int_stack + 5850;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 5868;
 Libderiv->deriv2_classes[1][4][81] = int_stack + 5898;
 Libderiv->deriv_classes[2][2][9] = int_stack + 5943;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 5979;
 Libderiv->deriv_classes[2][3][9] = int_stack + 6015;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 6075;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 6135;
 Libderiv->deriv2_classes[1][2][80] = int_stack + 6225;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 6243;
 Libderiv->deriv2_classes[1][4][80] = int_stack + 6273;
 Libderiv->deriv_classes[2][2][8] = int_stack + 6318;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 6354;
 Libderiv->deriv_classes[2][3][8] = int_stack + 6390;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 6450;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 6510;
 Libderiv->deriv2_classes[1][2][79] = int_stack + 6600;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 6618;
 Libderiv->deriv2_classes[1][4][79] = int_stack + 6648;
 Libderiv->deriv_classes[2][2][7] = int_stack + 6693;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 6729;
 Libderiv->deriv_classes[2][3][7] = int_stack + 6765;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 6825;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 6885;
 Libderiv->deriv2_classes[1][2][78] = int_stack + 6975;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 6993;
 Libderiv->deriv2_classes[1][4][78] = int_stack + 7023;
 Libderiv->dvrr_classes[2][2] = int_stack + 7068;
 Libderiv->deriv_classes[2][2][6] = int_stack + 7104;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 7140;
 Libderiv->deriv_classes[2][3][6] = int_stack + 7176;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 7236;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 7296;
 Libderiv->deriv2_classes[1][2][35] = int_stack + 7386;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 7404;
 Libderiv->deriv2_classes[1][4][35] = int_stack + 7434;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 7479;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 7515;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 7575;
 Libderiv->deriv2_classes[1][2][34] = int_stack + 7665;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 7683;
 Libderiv->deriv2_classes[1][4][34] = int_stack + 7713;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 7758;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 7794;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 7854;
 Libderiv->deriv2_classes[1][2][33] = int_stack + 7944;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 7962;
 Libderiv->deriv2_classes[1][4][33] = int_stack + 7992;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 8037;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 8073;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 8133;
 Libderiv->deriv2_classes[1][2][32] = int_stack + 8223;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 8241;
 Libderiv->deriv2_classes[1][4][32] = int_stack + 8271;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 8316;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 8352;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 8412;
 Libderiv->deriv2_classes[1][2][31] = int_stack + 8502;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 8520;
 Libderiv->deriv2_classes[1][4][31] = int_stack + 8550;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 8595;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 8631;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 8691;
 Libderiv->deriv2_classes[1][2][30] = int_stack + 8781;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 8799;
 Libderiv->deriv2_classes[1][4][30] = int_stack + 8829;
 Libderiv->deriv_classes[2][2][2] = int_stack + 8874;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 8910;
 Libderiv->deriv_classes[2][3][2] = int_stack + 8946;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 9006;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 9066;
 Libderiv->deriv2_classes[1][2][26] = int_stack + 9156;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 9174;
 Libderiv->deriv2_classes[1][4][26] = int_stack + 9204;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 9249;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 9285;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 9345;
 Libderiv->deriv2_classes[1][2][23] = int_stack + 9435;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 9453;
 Libderiv->deriv2_classes[1][4][23] = int_stack + 9483;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 9528;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 9564;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 9624;
 Libderiv->deriv2_classes[1][2][22] = int_stack + 9714;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 9732;
 Libderiv->deriv2_classes[1][4][22] = int_stack + 9762;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 9807;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 9843;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 9903;
 Libderiv->deriv2_classes[1][2][21] = int_stack + 9993;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 10011;
 Libderiv->deriv2_classes[1][4][21] = int_stack + 10041;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 10086;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 10122;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 10182;
 Libderiv->deriv2_classes[1][2][20] = int_stack + 10272;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 10290;
 Libderiv->deriv2_classes[1][4][20] = int_stack + 10320;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 10365;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 10401;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 10461;
 Libderiv->deriv2_classes[1][2][19] = int_stack + 10551;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 10569;
 Libderiv->deriv2_classes[1][4][19] = int_stack + 10599;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 10644;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 10680;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 10740;
 Libderiv->deriv2_classes[1][2][18] = int_stack + 10830;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 10848;
 Libderiv->deriv2_classes[1][4][18] = int_stack + 10878;
 Libderiv->deriv_classes[2][2][1] = int_stack + 10923;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 10959;
 Libderiv->deriv_classes[2][3][1] = int_stack + 10995;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 11055;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 11115;
 Libderiv->deriv2_classes[1][2][14] = int_stack + 11205;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 11223;
 Libderiv->deriv2_classes[1][4][14] = int_stack + 11253;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 11298;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 11334;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 11394;
 Libderiv->deriv2_classes[1][2][13] = int_stack + 11484;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 11502;
 Libderiv->deriv2_classes[1][4][13] = int_stack + 11532;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 11577;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 11613;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 11673;
 Libderiv->deriv_classes[1][2][11] = int_stack + 11763;
 Libderiv->deriv_classes[1][3][11] = int_stack + 11781;
 Libderiv->deriv_classes[1][4][11] = int_stack + 11811;
 Libderiv->deriv2_classes[1][2][11] = int_stack + 11856;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 11874;
 Libderiv->deriv2_classes[1][4][11] = int_stack + 11904;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 11949;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 11985;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 12045;
 Libderiv->deriv_classes[1][2][10] = int_stack + 12135;
 Libderiv->deriv_classes[1][3][10] = int_stack + 12153;
 Libderiv->deriv_classes[1][4][10] = int_stack + 12183;
 Libderiv->deriv2_classes[1][2][10] = int_stack + 12228;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 12246;
 Libderiv->deriv2_classes[1][4][10] = int_stack + 12276;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 12321;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 12357;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 12417;
 Libderiv->deriv_classes[1][2][9] = int_stack + 12507;
 Libderiv->deriv_classes[1][3][9] = int_stack + 12525;
 Libderiv->deriv_classes[1][4][9] = int_stack + 12555;
 Libderiv->deriv2_classes[1][2][9] = int_stack + 12600;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 12618;
 Libderiv->deriv2_classes[1][4][9] = int_stack + 12648;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 12693;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 12729;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 12789;
 Libderiv->deriv_classes[1][2][8] = int_stack + 12879;
 Libderiv->deriv_classes[1][3][8] = int_stack + 12897;
 Libderiv->deriv_classes[1][4][8] = int_stack + 12927;
 Libderiv->deriv2_classes[1][2][8] = int_stack + 12972;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 12990;
 Libderiv->deriv2_classes[1][4][8] = int_stack + 13020;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 13065;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 13101;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 13161;
 Libderiv->deriv_classes[1][2][7] = int_stack + 13251;
 Libderiv->deriv_classes[1][3][7] = int_stack + 13269;
 Libderiv->deriv_classes[1][4][7] = int_stack + 13299;
 Libderiv->deriv2_classes[1][2][7] = int_stack + 13344;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 13362;
 Libderiv->deriv2_classes[1][4][7] = int_stack + 13392;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 13437;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 13473;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 13533;
 Libderiv->dvrr_classes[1][2] = int_stack + 13623;
 Libderiv->deriv_classes[1][2][6] = int_stack + 13641;
 Libderiv->dvrr_classes[1][3] = int_stack + 13659;
 Libderiv->deriv_classes[1][3][6] = int_stack + 13689;
 Libderiv->deriv_classes[1][4][6] = int_stack + 13719;
 Libderiv->deriv2_classes[1][2][6] = int_stack + 13764;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 13782;
 Libderiv->deriv2_classes[1][4][6] = int_stack + 13812;
 Libderiv->deriv_classes[2][2][0] = int_stack + 13857;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 13893;
 Libderiv->deriv_classes[2][3][0] = int_stack + 13929;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 13989;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 14049;
 Libderiv->deriv_classes[1][2][2] = int_stack + 14139;
 Libderiv->deriv_classes[1][3][2] = int_stack + 14157;
 Libderiv->deriv_classes[1][4][2] = int_stack + 14187;
 Libderiv->deriv2_classes[1][2][2] = int_stack + 14232;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 14250;
 Libderiv->deriv2_classes[1][4][2] = int_stack + 14280;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 14325;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 14361;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 14421;
 Libderiv->deriv_classes[1][2][1] = int_stack + 14511;
 Libderiv->deriv_classes[1][3][1] = int_stack + 14529;
 Libderiv->deriv_classes[1][4][1] = int_stack + 14559;
 Libderiv->deriv2_classes[1][2][1] = int_stack + 14604;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 14622;
 Libderiv->deriv2_classes[1][4][1] = int_stack + 14652;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 14697;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 14733;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 14793;
 Libderiv->deriv_classes[1][2][0] = int_stack + 14883;
 Libderiv->deriv_classes[1][3][0] = int_stack + 14901;
 Libderiv->deriv_classes[1][4][0] = int_stack + 14931;
 Libderiv->deriv2_classes[1][2][0] = int_stack + 14976;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 14994;
 Libderiv->deriv2_classes[1][4][0] = int_stack + 15024;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 15069;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 15105;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 15165;
 memset(int_stack,0,122040);

 Libderiv->dvrr_stack = int_stack + 24813;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ppdd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+13659,int_stack+13623,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15309,int_stack+11781,int_stack+11763, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13623,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+11811,int_stack+11781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13659,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15453,int_stack+15363,int_stack+15309, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15255,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15561,int_stack+450,int_stack+7068,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15669,int_stack+5265,int_stack+5193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7068,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15777,int_stack+0,int_stack+5265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15957,int_stack+15777,int_stack+15669, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15561,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15777,int_stack+12153,int_stack+12135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13623, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+12183,int_stack+12153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13659, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15831,int_stack+0,int_stack+15777, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15255, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+16173,int_stack+5640,int_stack+5568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7068, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16281,int_stack+90,int_stack+5640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16461,int_stack+16281,int_stack+16173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15561, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+16281,int_stack+12525,int_stack+12507, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13623, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+12555,int_stack+12525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13659, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16335,int_stack+15363,int_stack+16281, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15255, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+6015,int_stack+5943, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7068, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16677,int_stack+180,int_stack+6015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16857,int_stack+16677,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15561, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+16677,int_stack+12897,int_stack+12879, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+12927,int_stack+12897, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16731,int_stack+15363,int_stack+16677, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+6390,int_stack+6318, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17073,int_stack+270,int_stack+6390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17253,int_stack+17073,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15561, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17073,int_stack+13269,int_stack+13251, 0.0, zero_stack, 1.0, int_stack+13623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+13299,int_stack+13269, 0.0, zero_stack, 1.0, int_stack+13659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17127,int_stack+15363,int_stack+17073, 0.0, zero_stack, 1.0, int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+6765,int_stack+6693, 0.0, zero_stack, 1.0, int_stack+7068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17469,int_stack+360,int_stack+6765, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17649,int_stack+17469,int_stack+216, 0.0, zero_stack, 1.0, int_stack+15561, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17469,int_stack+13689,int_stack+13641, 1.0, int_stack+13623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+13719,int_stack+13689, 1.0, int_stack+13659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17523,int_stack+15363,int_stack+17469, 1.0, int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+7176,int_stack+7104, 1.0, int_stack+7068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17865,int_stack+510,int_stack+7176, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+18045,int_stack+17865,int_stack+324, 1.0, int_stack+15561, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+780,int_stack+13659,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+15561,int_stack+15363,int_stack+15255,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+14157,int_stack+14139,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+14187,int_stack+14157,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+17865,int_stack+15363,int_stack+15255,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+8946,int_stack+8874,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+18261,int_stack+600,int_stack+8946,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+18441,int_stack+18261,int_stack+432,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+14529,int_stack+14511,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+14559,int_stack+14529,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+18315,int_stack+15363,int_stack+18261,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+10995,int_stack+10923,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+18657,int_stack+690,int_stack+10995,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+18837,int_stack+18657,int_stack+540,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+18657,int_stack+14901,int_stack+14883,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+14931,int_stack+14901,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+18711,int_stack+15363,int_stack+18657,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+13929,int_stack+13857,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19053,int_stack+825,int_stack+13929,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+19233,int_stack+19053,int_stack+648,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+933,int_stack+915, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11763,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+963,int_stack+933, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11781,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19107,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15309,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+756,int_stack+1044,int_stack+1008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5193,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+864,int_stack+1104,int_stack+1044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5265,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19449,int_stack+864,int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15669,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+1212,int_stack+1194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11763, 1.0, int_stack+12135,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+1242,int_stack+1212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11781, 1.0, int_stack+12153,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+756,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15309, 1.0, int_stack+15777,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+864,int_stack+1323,int_stack+1287, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5193, 1.0, int_stack+5568,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+972,int_stack+1383,int_stack+1323, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5265, 1.0, int_stack+5640,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1152,int_stack+972,int_stack+864, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15669, 1.0, int_stack+16173,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+1491,int_stack+1473, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12135, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+1521,int_stack+1491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12153, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+864,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15777, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+972,int_stack+1602,int_stack+1566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5568, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1368,int_stack+1662,int_stack+1602, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5640, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19665,int_stack+1368,int_stack+972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16173, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+1770,int_stack+1752, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11763, 0.0, zero_stack, 1.0, int_stack+12507,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+1800,int_stack+1770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11781, 0.0, zero_stack, 1.0, int_stack+12525,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+972,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15309, 0.0, zero_stack, 1.0, int_stack+16281,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1368,int_stack+1881,int_stack+1845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5193, 0.0, zero_stack, 1.0, int_stack+5943,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1476,int_stack+1941,int_stack+1881, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5265, 0.0, zero_stack, 1.0, int_stack+6015,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1656,int_stack+1476,int_stack+1368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15669, 0.0, zero_stack, 1.0, int_stack+0,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+2049,int_stack+2031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12135, 1.0, int_stack+12507, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+2079,int_stack+2049, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12153, 1.0, int_stack+12525, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1368,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15777, 1.0, int_stack+16281, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1476,int_stack+2160,int_stack+2124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5568, 1.0, int_stack+5943, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1872,int_stack+2220,int_stack+2160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5640, 1.0, int_stack+6015, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2052,int_stack+1872,int_stack+1476, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16173, 1.0, int_stack+0, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+2328,int_stack+2310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12507, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+2358,int_stack+2328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12525, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1476,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16281, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1872,int_stack+2439,int_stack+2403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5943, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+19881,int_stack+2499,int_stack+2439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6015, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2268,int_stack+19881,int_stack+1872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+2607,int_stack+2589, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11763, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12879,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+2637,int_stack+2607, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11781, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12897,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1872,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15309, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16677,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19881,int_stack+2718,int_stack+2682, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5193, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6318,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+19989,int_stack+2778,int_stack+2718, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5265, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6390,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20169,int_stack+19989,int_stack+19881, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15669, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+2886,int_stack+2868, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12135, 0.0, zero_stack, 1.0, int_stack+12879, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+2916,int_stack+2886, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12153, 0.0, zero_stack, 1.0, int_stack+12897, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19881,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15777, 0.0, zero_stack, 1.0, int_stack+16677, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19989,int_stack+2997,int_stack+2961, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5568, 0.0, zero_stack, 1.0, int_stack+6318, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20385,int_stack+3057,int_stack+2997, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5640, 0.0, zero_stack, 1.0, int_stack+6390, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20565,int_stack+20385,int_stack+19989, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16173, 0.0, zero_stack, 1.0, int_stack+108, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+3165,int_stack+3147, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12507, 1.0, int_stack+12879, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+3195,int_stack+3165, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12525, 1.0, int_stack+12897, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19989,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16281, 1.0, int_stack+16677, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+20385,int_stack+3276,int_stack+3240, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5943, 1.0, int_stack+6318, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20781,int_stack+3336,int_stack+3276, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6015, 1.0, int_stack+6390, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2484,int_stack+20781,int_stack+20385, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+3444,int_stack+3426, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12879, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+3474,int_stack+3444, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20385,int_stack+15363,int_stack+19053, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16677, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+20781,int_stack+3555,int_stack+3519, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+3615,int_stack+3555, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2880,int_stack+2700,int_stack+20781, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+3723,int_stack+3705, 0.0, zero_stack, 1.0, int_stack+11763, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13251,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+3753,int_stack+3723, 0.0, zero_stack, 1.0, int_stack+11781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13269,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20781,int_stack+15363,int_stack+19053, 0.0, zero_stack, 1.0, int_stack+15309, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17073,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+20889,int_stack+3834,int_stack+3798, 0.0, zero_stack, 1.0, int_stack+5193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6693,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+3894,int_stack+3834, 0.0, zero_stack, 1.0, int_stack+5265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6765,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3096,int_stack+2700,int_stack+20889, 0.0, zero_stack, 1.0, int_stack+15669, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+4002,int_stack+3984, 0.0, zero_stack, 1.0, int_stack+12135, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13251, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+4032,int_stack+4002, 0.0, zero_stack, 1.0, int_stack+12153, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13269, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20889,int_stack+15363,int_stack+19053, 0.0, zero_stack, 1.0, int_stack+15777, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17073, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+4113,int_stack+4077, 0.0, zero_stack, 1.0, int_stack+5568, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6693, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3312,int_stack+4173,int_stack+4113, 0.0, zero_stack, 1.0, int_stack+5640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6765, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3492,int_stack+3312,int_stack+2700, 0.0, zero_stack, 1.0, int_stack+16173, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+4281,int_stack+4263, 0.0, zero_stack, 1.0, int_stack+12507, 0.0, zero_stack, 1.0, int_stack+13251, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+4311,int_stack+4281, 0.0, zero_stack, 1.0, int_stack+12525, 0.0, zero_stack, 1.0, int_stack+13269, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2700,int_stack+15363,int_stack+19053, 0.0, zero_stack, 1.0, int_stack+16281, 0.0, zero_stack, 1.0, int_stack+17073, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3312,int_stack+4392,int_stack+4356, 0.0, zero_stack, 1.0, int_stack+5943, 0.0, zero_stack, 1.0, int_stack+6693, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3708,int_stack+4452,int_stack+4392, 0.0, zero_stack, 1.0, int_stack+6015, 0.0, zero_stack, 1.0, int_stack+6765, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3888,int_stack+3708,int_stack+3312, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+4560,int_stack+4542, 0.0, zero_stack, 1.0, int_stack+12879, 1.0, int_stack+13251, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+4590,int_stack+4560, 0.0, zero_stack, 1.0, int_stack+12897, 1.0, int_stack+13269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3312,int_stack+15363,int_stack+19053, 0.0, zero_stack, 1.0, int_stack+16677, 1.0, int_stack+17073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3708,int_stack+4671,int_stack+4635, 0.0, zero_stack, 1.0, int_stack+6318, 1.0, int_stack+6693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4104,int_stack+4731,int_stack+4671, 0.0, zero_stack, 1.0, int_stack+6390, 1.0, int_stack+6765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4284,int_stack+4104,int_stack+3708, 0.0, zero_stack, 1.0, int_stack+108, 1.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+4839,int_stack+4821, 0.0, zero_stack, 2.0, int_stack+13251, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+4869,int_stack+4839, 0.0, zero_stack, 2.0, int_stack+13269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3708,int_stack+15363,int_stack+19053, 0.0, zero_stack, 2.0, int_stack+17073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4104,int_stack+4950,int_stack+4914, 0.0, zero_stack, 2.0, int_stack+6693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+5010,int_stack+4950, 0.0, zero_stack, 2.0, int_stack+6765, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4680,int_stack+4500,int_stack+4104, 0.0, zero_stack, 2.0, int_stack+216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+5118,int_stack+5100, 1.0, int_stack+11763, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13641,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15363,int_stack+5148,int_stack+5118, 1.0, int_stack+11781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13689,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4104,int_stack+15363,int_stack+19053, 1.0, int_stack+15309, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17469,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15309,int_stack+5325,int_stack+5229, 1.0, int_stack+5193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7104,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+5385,int_stack+5325, 1.0, int_stack+5265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7176,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4896,int_stack+4500,int_stack+15309, 1.0, int_stack+15669, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+5493,int_stack+5475, 1.0, int_stack+12135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13641, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15669,int_stack+5523,int_stack+5493, 1.0, int_stack+12153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13689, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15309,int_stack+15669,int_stack+19053, 1.0, int_stack+15777, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17469, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15669,int_stack+5700,int_stack+5604, 1.0, int_stack+5568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7104, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+5760,int_stack+5700, 1.0, int_stack+5640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7176, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5112,int_stack+4500,int_stack+15669, 1.0, int_stack+16173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+5868,int_stack+5850, 1.0, int_stack+12507, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13641, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16173,int_stack+5898,int_stack+5868, 1.0, int_stack+12525, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13689, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15669,int_stack+16173,int_stack+19053, 1.0, int_stack+16281, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17469, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+16173,int_stack+6075,int_stack+5979, 1.0, int_stack+5943, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7104, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+6135,int_stack+6075, 1.0, int_stack+6015, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7176, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5328,int_stack+4500,int_stack+16173, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+19053,int_stack+6243,int_stack+6225, 1.0, int_stack+12879, 0.0, zero_stack, 1.0, int_stack+13641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+6273,int_stack+6243, 1.0, int_stack+12897, 0.0, zero_stack, 1.0, int_stack+13689, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16173,int_stack+0,int_stack+19053, 1.0, int_stack+16677, 0.0, zero_stack, 1.0, int_stack+17469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+6450,int_stack+6354, 1.0, int_stack+6318, 0.0, zero_stack, 1.0, int_stack+7104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+6510,int_stack+6450, 1.0, int_stack+6390, 0.0, zero_stack, 1.0, int_stack+7176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5544,int_stack+4500,int_stack+0, 1.0, int_stack+108, 0.0, zero_stack, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+16677,int_stack+6618,int_stack+6600, 1.0, int_stack+13251, 1.0, int_stack+13641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+6648,int_stack+6618, 1.0, int_stack+13269, 1.0, int_stack+13689, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+90,int_stack+0,int_stack+16677, 1.0, int_stack+17073, 1.0, int_stack+17469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4500,int_stack+6825,int_stack+6729, 1.0, int_stack+6693, 1.0, int_stack+7104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5760,int_stack+6885,int_stack+6825, 1.0, int_stack+6765, 1.0, int_stack+7176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5940,int_stack+5760,int_stack+4500, 1.0, int_stack+216, 1.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17073,int_stack+6993,int_stack+6975, 2.0, int_stack+13641, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+7023,int_stack+6993, 2.0, int_stack+13689, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4500,int_stack+0,int_stack+17073, 2.0, int_stack+17469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5760,int_stack+7236,int_stack+7140, 2.0, int_stack+7104, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6156,int_stack+7296,int_stack+7236, 2.0, int_stack+7176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6336,int_stack+6156,int_stack+5760, 2.0, int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17469,int_stack+7404,int_stack+7386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14139,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+7434,int_stack+7404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14157,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5760,int_stack+0,int_stack+17469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15255,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6156,int_stack+7515,int_stack+7479, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8874,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+198,int_stack+7575,int_stack+7515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8946,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6552,int_stack+198,int_stack+6156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17469,int_stack+7683,int_stack+7665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14139, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+7713,int_stack+7683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14157, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6156,int_stack+0,int_stack+17469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15255, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+198,int_stack+7794,int_stack+7758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8874, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6768,int_stack+7854,int_stack+7794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8946, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6948,int_stack+6768,int_stack+198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17469,int_stack+7962,int_stack+7944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14139, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+7992,int_stack+7962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14157, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+198,int_stack+0,int_stack+17469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15255, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+306,int_stack+8073,int_stack+8037, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8874, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6768,int_stack+8133,int_stack+8073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8946, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7164,int_stack+6768,int_stack+306, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17469,int_stack+8241,int_stack+8223, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+8271,int_stack+8241, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14157, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+306,int_stack+0,int_stack+17469, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6768,int_stack+8352,int_stack+8316, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8874, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7380,int_stack+8412,int_stack+8352, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7560,int_stack+7380,int_stack+6768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17469,int_stack+8520,int_stack+8502, 0.0, zero_stack, 1.0, int_stack+14139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+8550,int_stack+8520, 0.0, zero_stack, 1.0, int_stack+14157, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6768,int_stack+0,int_stack+17469, 0.0, zero_stack, 1.0, int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7380,int_stack+8631,int_stack+8595, 0.0, zero_stack, 1.0, int_stack+8874, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7776,int_stack+8691,int_stack+8631, 0.0, zero_stack, 1.0, int_stack+8946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7956,int_stack+7776,int_stack+7380, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+17469,int_stack+8799,int_stack+8781, 1.0, int_stack+14139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+8829,int_stack+8799, 1.0, int_stack+14157, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7380,int_stack+0,int_stack+17469, 1.0, int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7776,int_stack+9006,int_stack+8910, 1.0, int_stack+8874, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8172,int_stack+9066,int_stack+9006, 1.0, int_stack+8946, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8352,int_stack+8172,int_stack+7776, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+9174,int_stack+9156,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+9204,int_stack+9174,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+7776,int_stack+0,int_stack+15255,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+8172,int_stack+9285,int_stack+9249,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8568,int_stack+9345,int_stack+9285,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+8748,int_stack+8568,int_stack+8172,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+9453,int_stack+9435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14511,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+9483,int_stack+9453, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14529,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8172,int_stack+0,int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18261,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8568,int_stack+9564,int_stack+9528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10923,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8964,int_stack+9624,int_stack+9564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10995,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9144,int_stack+8964,int_stack+8568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+9732,int_stack+9714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14511, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+9762,int_stack+9732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14529, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8568,int_stack+0,int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18261, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8964,int_stack+9843,int_stack+9807, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10923, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9360,int_stack+9903,int_stack+9843, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10995, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9540,int_stack+9360,int_stack+8964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+10011,int_stack+9993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14511, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+10041,int_stack+10011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14529, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8964,int_stack+0,int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18261, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9360,int_stack+10122,int_stack+10086, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10923, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9756,int_stack+10182,int_stack+10122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10995, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9936,int_stack+9756,int_stack+9360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+10290,int_stack+10272, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+10320,int_stack+10290, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9360,int_stack+0,int_stack+15255, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9756,int_stack+10401,int_stack+10365, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10152,int_stack+10461,int_stack+10401, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10332,int_stack+10152,int_stack+9756, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+10569,int_stack+10551, 0.0, zero_stack, 1.0, int_stack+14511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+10599,int_stack+10569, 0.0, zero_stack, 1.0, int_stack+14529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9756,int_stack+0,int_stack+15255, 0.0, zero_stack, 1.0, int_stack+18261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+10152,int_stack+10680,int_stack+10644, 0.0, zero_stack, 1.0, int_stack+10923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20997,int_stack+10740,int_stack+10680, 0.0, zero_stack, 1.0, int_stack+10995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21177,int_stack+20997,int_stack+10152, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+15255,int_stack+10848,int_stack+10830, 1.0, int_stack+14511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+10878,int_stack+10848, 1.0, int_stack+14529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10152,int_stack+0,int_stack+15255, 1.0, int_stack+18261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+20997,int_stack+11055,int_stack+10959, 1.0, int_stack+10923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21393,int_stack+11115,int_stack+11055, 1.0, int_stack+10995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21573,int_stack+21393,int_stack+20997, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+11223,int_stack+11205,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+11253,int_stack+11223,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+20997,int_stack+0,int_stack+18261,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+21393,int_stack+11334,int_stack+11298,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+414,int_stack+11394,int_stack+11334,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+10548,int_stack+414,int_stack+21393,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+11502,int_stack+11484,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+11532,int_stack+11502,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+21393,int_stack+0,int_stack+18261,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+414,int_stack+11613,int_stack+11577,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+21789,int_stack+11673,int_stack+11613,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+10764,int_stack+21789,int_stack+414,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+11874,int_stack+11856, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14883,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+11904,int_stack+11874, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14901,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+414,int_stack+0,int_stack+18261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18657,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+522,int_stack+11985,int_stack+11949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13857,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21789,int_stack+12045,int_stack+11985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13929,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10980,int_stack+21789,int_stack+522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+12246,int_stack+12228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14883, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+12276,int_stack+12246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14901, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+522,int_stack+0,int_stack+18261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18657, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+21789,int_stack+12357,int_stack+12321, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13857, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11196,int_stack+12417,int_stack+12357, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13929, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11376,int_stack+11196,int_stack+21789, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+12618,int_stack+12600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14883, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+12648,int_stack+12618, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14901, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21789,int_stack+0,int_stack+18261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18657, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11196,int_stack+12729,int_stack+12693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13857, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11592,int_stack+12789,int_stack+12729, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13929, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11772,int_stack+11592,int_stack+11196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+12990,int_stack+12972, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+13020,int_stack+12990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11196,int_stack+0,int_stack+18261, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18657, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11592,int_stack+13101,int_stack+13065, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13857, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11988,int_stack+13161,int_stack+13101, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12168,int_stack+11988,int_stack+11592, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+13362,int_stack+13344, 0.0, zero_stack, 1.0, int_stack+14883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+13392,int_stack+13362, 0.0, zero_stack, 1.0, int_stack+14901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11592,int_stack+0,int_stack+18261, 0.0, zero_stack, 1.0, int_stack+18657, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+11988,int_stack+13473,int_stack+13437, 0.0, zero_stack, 1.0, int_stack+13857, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12384,int_stack+13533,int_stack+13473, 0.0, zero_stack, 1.0, int_stack+13929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12564,int_stack+12384,int_stack+11988, 0.0, zero_stack, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18261,int_stack+13782,int_stack+13764, 1.0, int_stack+14883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+13812,int_stack+13782, 1.0, int_stack+14901, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11988,int_stack+0,int_stack+18261, 1.0, int_stack+18657, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+12384,int_stack+13989,int_stack+13893, 1.0, int_stack+13857, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12780,int_stack+14049,int_stack+13989, 1.0, int_stack+13929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12960,int_stack+12780,int_stack+12384, 1.0, int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+18657,int_stack+14250,int_stack+14232,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+14280,int_stack+14250,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+12384,int_stack+0,int_stack+18657,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+12780,int_stack+14361,int_stack+14325,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13176,int_stack+14421,int_stack+14361,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+13356,int_stack+13176,int_stack+12780,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+18657,int_stack+14622,int_stack+14604,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+14652,int_stack+14622,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+12780,int_stack+0,int_stack+18657,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+13176,int_stack+14733,int_stack+14697,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13572,int_stack+14793,int_stack+14733,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+13752,int_stack+13572,int_stack+13176,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+18657,int_stack+14994,int_stack+14976,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+15024,int_stack+14994,3);
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+13176,int_stack+0,int_stack+18657,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+13572,int_stack+15105,int_stack+15069,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+13968,int_stack+15165,int_stack+15105,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+14148,int_stack+13968,int_stack+13572,6);
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+14364,int_stack+15957,int_stack+15453,36);
     Libderiv->ABCD[11] = int_stack + 14364;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+14688,int_stack+16461,int_stack+15831,36);
     Libderiv->ABCD[10] = int_stack + 14688;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+21897,int_stack+16857,int_stack+16335,36);
     Libderiv->ABCD[9] = int_stack + 21897;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+22221,int_stack+17253,int_stack+16731,36);
     Libderiv->ABCD[8] = int_stack + 22221;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+22545,int_stack+17649,int_stack+17127,36);
     Libderiv->ABCD[7] = int_stack + 22545;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+22869,int_stack+18045,int_stack+17523,36);
     Libderiv->ABCD[6] = int_stack + 22869;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+23193,int_stack+18441,int_stack+17865, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15561, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[2] = int_stack + 23193;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+23517,int_stack+18837,int_stack+18315, 0.0, zero_stack, 1.0, int_stack+15561, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[1] = int_stack + 23517;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+17973,int_stack+19233,int_stack+18711, 1.0, int_stack+15561, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[0] = int_stack + 17973;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+23841,int_stack+19449,int_stack+19107,36);
     Libderiv->ABCD[155] = int_stack + 23841;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+24165,int_stack+1152,int_stack+756,36);
     Libderiv->ABCD[143] = int_stack + 24165;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+24489,int_stack+19665,int_stack+864,36);
     Libderiv->ABCD[142] = int_stack + 24489;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+630,int_stack+1656,int_stack+972,36);
     Libderiv->ABCD[131] = int_stack + 630;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+954,int_stack+2052,int_stack+1368,36);
     Libderiv->ABCD[130] = int_stack + 954;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+18819,int_stack+2268,int_stack+1476,36);
     Libderiv->ABCD[129] = int_stack + 18819;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1278,int_stack+20169,int_stack+1872,36);
     Libderiv->ABCD[119] = int_stack + 1278;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1602,int_stack+20565,int_stack+19881,36);
     Libderiv->ABCD[118] = int_stack + 1602;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1926,int_stack+2484,int_stack+19989,36);
     Libderiv->ABCD[117] = int_stack + 1926;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2250,int_stack+2880,int_stack+20385,36);
     Libderiv->ABCD[116] = int_stack + 2250;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+19143,int_stack+3096,int_stack+20781,36);
     Libderiv->ABCD[107] = int_stack + 19143;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2808,int_stack+3492,int_stack+20889,36);
     Libderiv->ABCD[106] = int_stack + 2808;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+19467,int_stack+3888,int_stack+2700,36);
     Libderiv->ABCD[105] = int_stack + 19467;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+19791,int_stack+4284,int_stack+3312,36);
     Libderiv->ABCD[104] = int_stack + 19791;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3132,int_stack+4680,int_stack+3708,36);
     Libderiv->ABCD[103] = int_stack + 3132;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3456,int_stack+4896,int_stack+4104,36);
     Libderiv->ABCD[95] = int_stack + 3456;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4608,int_stack+5112,int_stack+15309,36);
     Libderiv->ABCD[94] = int_stack + 4608;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4932,int_stack+5328,int_stack+15669,36);
     Libderiv->ABCD[93] = int_stack + 4932;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3780,int_stack+5544,int_stack+16173,36);
     Libderiv->ABCD[92] = int_stack + 3780;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+5256,int_stack+5940,int_stack+90,36);
     Libderiv->ABCD[91] = int_stack + 5256;
 /*--- compute (pp|dd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4104,int_stack+6336,int_stack+4500,36);
     Libderiv->ABCD[90] = int_stack + 4104;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+15012,int_stack+6552,int_stack+5760, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15453, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[47] = int_stack + 15012;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6264,int_stack+6948,int_stack+6156, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[46] = int_stack + 6264;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5580,int_stack+7164,int_stack+198, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[45] = int_stack + 5580;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+6876,int_stack+7560,int_stack+306, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16731, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[44] = int_stack + 6876;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+7956,int_stack+6768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[43] = int_stack + 0;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+5904,int_stack+8352,int_stack+7380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[42] = int_stack + 5904;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7200,int_stack+8748,int_stack+7776, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[38] = int_stack + 7200;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7524,int_stack+9144,int_stack+8172, 0.0, zero_stack, 1.0, int_stack+15453, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[35] = int_stack + 7524;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7848,int_stack+9540,int_stack+8568, 0.0, zero_stack, 1.0, int_stack+15831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[34] = int_stack + 7848;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8172,int_stack+9936,int_stack+8964, 0.0, zero_stack, 1.0, int_stack+16335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[33] = int_stack + 8172;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8496,int_stack+10332,int_stack+9360, 0.0, zero_stack, 1.0, int_stack+16731, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[32] = int_stack + 8496;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8820,int_stack+21177,int_stack+9756, 0.0, zero_stack, 1.0, int_stack+17127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[31] = int_stack + 8820;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9144,int_stack+21573,int_stack+10152, 0.0, zero_stack, 1.0, int_stack+17523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[30] = int_stack + 9144;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9468,int_stack+10548,int_stack+20997, 0.0, zero_stack, 1.0, int_stack+17865, 1.0, int_stack+18315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[26] = int_stack + 9468;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9792,int_stack+10764,int_stack+21393, 0.0, zero_stack, 2.0, int_stack+18315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[25] = int_stack + 9792;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+10116,int_stack+10980,int_stack+414, 1.0, int_stack+15453, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[23] = int_stack + 10116;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+10440,int_stack+11376,int_stack+522, 1.0, int_stack+15831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[22] = int_stack + 10440;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+10764,int_stack+11772,int_stack+21789, 1.0, int_stack+16335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[21] = int_stack + 10764;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+15336,int_stack+12168,int_stack+11196, 1.0, int_stack+16731, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[20] = int_stack + 15336;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+11088,int_stack+12564,int_stack+11592, 1.0, int_stack+17127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[19] = int_stack + 11088;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+11412,int_stack+12960,int_stack+11988, 1.0, int_stack+17523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[18] = int_stack + 11412;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+11736,int_stack+13356,int_stack+12384, 1.0, int_stack+17865, 0.0, zero_stack, 1.0, int_stack+18711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[14] = int_stack + 11736;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+13284,int_stack+13752,int_stack+12780, 1.0, int_stack+18315, 1.0, int_stack+18711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[13] = int_stack + 13284;
 /*--- compute (pp|dd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+13608,int_stack+14148,int_stack+13176, 2.0, int_stack+18711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[12] = int_stack + 13608;

}
