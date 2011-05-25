#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ppff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|ff) integrals */

void d12hrr_order_ppff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][6][11] = int_stack + 0;
 Libderiv->deriv_classes[2][6][10] = int_stack + 168;
 Libderiv->deriv_classes[2][6][9] = int_stack + 336;
 Libderiv->deriv_classes[2][6][8] = int_stack + 504;
 Libderiv->deriv_classes[2][6][7] = int_stack + 672;
 Libderiv->dvrr_classes[2][5] = int_stack + 840;
 Libderiv->deriv_classes[2][6][6] = int_stack + 966;
 Libderiv->deriv_classes[2][6][2] = int_stack + 1134;
 Libderiv->deriv_classes[2][6][1] = int_stack + 1302;
 Libderiv->dvrr_classes[1][6] = int_stack + 1470;
 Libderiv->deriv_classes[2][6][0] = int_stack + 1554;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 1722;
 Libderiv->deriv2_classes[1][4][143] = int_stack + 1752;
 Libderiv->deriv2_classes[1][5][143] = int_stack + 1797;
 Libderiv->deriv2_classes[1][6][143] = int_stack + 1860;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 1944;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 2004;
 Libderiv->deriv2_classes[2][5][143] = int_stack + 2094;
 Libderiv->deriv2_classes[2][6][143] = int_stack + 2220;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 2388;
 Libderiv->deriv2_classes[1][4][131] = int_stack + 2418;
 Libderiv->deriv2_classes[1][5][131] = int_stack + 2463;
 Libderiv->deriv2_classes[1][6][131] = int_stack + 2526;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 2610;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 2670;
 Libderiv->deriv2_classes[2][5][131] = int_stack + 2760;
 Libderiv->deriv2_classes[2][6][131] = int_stack + 2886;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 3054;
 Libderiv->deriv2_classes[1][4][130] = int_stack + 3084;
 Libderiv->deriv2_classes[1][5][130] = int_stack + 3129;
 Libderiv->deriv2_classes[1][6][130] = int_stack + 3192;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 3276;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 3336;
 Libderiv->deriv2_classes[2][5][130] = int_stack + 3426;
 Libderiv->deriv2_classes[2][6][130] = int_stack + 3552;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 3720;
 Libderiv->deriv2_classes[1][4][119] = int_stack + 3750;
 Libderiv->deriv2_classes[1][5][119] = int_stack + 3795;
 Libderiv->deriv2_classes[1][6][119] = int_stack + 3858;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 3942;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 4002;
 Libderiv->deriv2_classes[2][5][119] = int_stack + 4092;
 Libderiv->deriv2_classes[2][6][119] = int_stack + 4218;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 4386;
 Libderiv->deriv2_classes[1][4][118] = int_stack + 4416;
 Libderiv->deriv2_classes[1][5][118] = int_stack + 4461;
 Libderiv->deriv2_classes[1][6][118] = int_stack + 4524;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 4608;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 4668;
 Libderiv->deriv2_classes[2][5][118] = int_stack + 4758;
 Libderiv->deriv2_classes[2][6][118] = int_stack + 4884;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 5052;
 Libderiv->deriv2_classes[1][4][117] = int_stack + 5082;
 Libderiv->deriv2_classes[1][5][117] = int_stack + 5127;
 Libderiv->deriv2_classes[1][6][117] = int_stack + 5190;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 5274;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 5334;
 Libderiv->deriv2_classes[2][5][117] = int_stack + 5424;
 Libderiv->deriv2_classes[2][6][117] = int_stack + 5550;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 5718;
 Libderiv->deriv2_classes[1][4][107] = int_stack + 5748;
 Libderiv->deriv2_classes[1][5][107] = int_stack + 5793;
 Libderiv->deriv2_classes[1][6][107] = int_stack + 5856;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 5940;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 6000;
 Libderiv->deriv2_classes[2][5][107] = int_stack + 6090;
 Libderiv->deriv2_classes[2][6][107] = int_stack + 6216;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 6384;
 Libderiv->deriv2_classes[1][4][106] = int_stack + 6414;
 Libderiv->deriv2_classes[1][5][106] = int_stack + 6459;
 Libderiv->deriv2_classes[1][6][106] = int_stack + 6522;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 6606;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 6666;
 Libderiv->deriv2_classes[2][5][106] = int_stack + 6756;
 Libderiv->deriv2_classes[2][6][106] = int_stack + 6882;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 7050;
 Libderiv->deriv2_classes[1][4][105] = int_stack + 7080;
 Libderiv->deriv2_classes[1][5][105] = int_stack + 7125;
 Libderiv->deriv2_classes[1][6][105] = int_stack + 7188;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 7272;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 7332;
 Libderiv->deriv2_classes[2][5][105] = int_stack + 7422;
 Libderiv->deriv2_classes[2][6][105] = int_stack + 7548;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 7716;
 Libderiv->deriv2_classes[1][4][104] = int_stack + 7746;
 Libderiv->deriv2_classes[1][5][104] = int_stack + 7791;
 Libderiv->deriv2_classes[1][6][104] = int_stack + 7854;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 7938;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 7998;
 Libderiv->deriv2_classes[2][5][104] = int_stack + 8088;
 Libderiv->deriv2_classes[2][6][104] = int_stack + 8214;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 8382;
 Libderiv->deriv2_classes[1][4][95] = int_stack + 8412;
 Libderiv->deriv2_classes[1][5][95] = int_stack + 8457;
 Libderiv->deriv2_classes[1][6][95] = int_stack + 8520;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 8604;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 8664;
 Libderiv->deriv2_classes[2][5][95] = int_stack + 8754;
 Libderiv->deriv2_classes[2][6][95] = int_stack + 8880;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 9048;
 Libderiv->deriv2_classes[1][4][94] = int_stack + 9078;
 Libderiv->deriv2_classes[1][5][94] = int_stack + 9123;
 Libderiv->deriv2_classes[1][6][94] = int_stack + 9186;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 9270;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 9330;
 Libderiv->deriv2_classes[2][5][94] = int_stack + 9420;
 Libderiv->deriv2_classes[2][6][94] = int_stack + 9546;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 9714;
 Libderiv->deriv2_classes[1][4][93] = int_stack + 9744;
 Libderiv->deriv2_classes[1][5][93] = int_stack + 9789;
 Libderiv->deriv2_classes[1][6][93] = int_stack + 9852;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 9936;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 9996;
 Libderiv->deriv2_classes[2][5][93] = int_stack + 10086;
 Libderiv->deriv2_classes[2][6][93] = int_stack + 10212;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 10380;
 Libderiv->deriv2_classes[1][4][92] = int_stack + 10410;
 Libderiv->deriv2_classes[1][5][92] = int_stack + 10455;
 Libderiv->deriv2_classes[1][6][92] = int_stack + 10518;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 10602;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 10662;
 Libderiv->deriv2_classes[2][5][92] = int_stack + 10752;
 Libderiv->deriv2_classes[2][6][92] = int_stack + 10878;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 11046;
 Libderiv->deriv2_classes[1][4][91] = int_stack + 11076;
 Libderiv->deriv2_classes[1][5][91] = int_stack + 11121;
 Libderiv->deriv2_classes[1][6][91] = int_stack + 11184;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 11268;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 11328;
 Libderiv->deriv2_classes[2][5][91] = int_stack + 11418;
 Libderiv->deriv2_classes[2][6][91] = int_stack + 11544;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 11712;
 Libderiv->deriv2_classes[1][4][83] = int_stack + 11742;
 Libderiv->deriv2_classes[1][5][83] = int_stack + 11787;
 Libderiv->deriv2_classes[1][6][83] = int_stack + 11850;
 Libderiv->deriv_classes[2][3][11] = int_stack + 11934;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 11994;
 Libderiv->deriv_classes[2][4][11] = int_stack + 12054;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 12144;
 Libderiv->deriv_classes[2][5][11] = int_stack + 12234;
 Libderiv->deriv2_classes[2][5][83] = int_stack + 12360;
 Libderiv->deriv2_classes[2][6][83] = int_stack + 12486;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 12654;
 Libderiv->deriv2_classes[1][4][82] = int_stack + 12684;
 Libderiv->deriv2_classes[1][5][82] = int_stack + 12729;
 Libderiv->deriv2_classes[1][6][82] = int_stack + 12792;
 Libderiv->deriv_classes[2][3][10] = int_stack + 12876;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 12936;
 Libderiv->deriv_classes[2][4][10] = int_stack + 12996;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 13086;
 Libderiv->deriv_classes[2][5][10] = int_stack + 13176;
 Libderiv->deriv2_classes[2][5][82] = int_stack + 13302;
 Libderiv->deriv2_classes[2][6][82] = int_stack + 13428;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 13596;
 Libderiv->deriv2_classes[1][4][81] = int_stack + 13626;
 Libderiv->deriv2_classes[1][5][81] = int_stack + 13671;
 Libderiv->deriv2_classes[1][6][81] = int_stack + 13734;
 Libderiv->deriv_classes[2][3][9] = int_stack + 13818;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 13878;
 Libderiv->deriv_classes[2][4][9] = int_stack + 13938;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 14028;
 Libderiv->deriv_classes[2][5][9] = int_stack + 14118;
 Libderiv->deriv2_classes[2][5][81] = int_stack + 14244;
 Libderiv->deriv2_classes[2][6][81] = int_stack + 14370;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 14538;
 Libderiv->deriv2_classes[1][4][80] = int_stack + 14568;
 Libderiv->deriv2_classes[1][5][80] = int_stack + 14613;
 Libderiv->deriv2_classes[1][6][80] = int_stack + 14676;
 Libderiv->deriv_classes[2][3][8] = int_stack + 14760;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 14820;
 Libderiv->deriv_classes[2][4][8] = int_stack + 14880;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 14970;
 Libderiv->deriv_classes[2][5][8] = int_stack + 15060;
 Libderiv->deriv2_classes[2][5][80] = int_stack + 15186;
 Libderiv->deriv2_classes[2][6][80] = int_stack + 15312;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 15480;
 Libderiv->deriv2_classes[1][4][79] = int_stack + 15510;
 Libderiv->deriv2_classes[1][5][79] = int_stack + 15555;
 Libderiv->deriv2_classes[1][6][79] = int_stack + 15618;
 Libderiv->deriv_classes[2][3][7] = int_stack + 15702;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 15762;
 Libderiv->deriv_classes[2][4][7] = int_stack + 15822;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 15912;
 Libderiv->deriv_classes[2][5][7] = int_stack + 16002;
 Libderiv->deriv2_classes[2][5][79] = int_stack + 16128;
 Libderiv->deriv2_classes[2][6][79] = int_stack + 16254;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 16422;
 Libderiv->deriv2_classes[1][4][78] = int_stack + 16452;
 Libderiv->deriv2_classes[1][5][78] = int_stack + 16497;
 Libderiv->deriv2_classes[1][6][78] = int_stack + 16560;
 Libderiv->dvrr_classes[2][3] = int_stack + 16644;
 Libderiv->deriv_classes[2][3][6] = int_stack + 16704;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 16764;
 Libderiv->dvrr_classes[2][4] = int_stack + 16824;
 Libderiv->deriv_classes[2][4][6] = int_stack + 16914;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 17004;
 Libderiv->deriv_classes[2][5][6] = int_stack + 17094;
 Libderiv->deriv2_classes[2][5][78] = int_stack + 17220;
 Libderiv->deriv2_classes[2][6][78] = int_stack + 17346;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 17514;
 Libderiv->deriv2_classes[1][4][35] = int_stack + 17544;
 Libderiv->deriv2_classes[1][5][35] = int_stack + 17589;
 Libderiv->deriv2_classes[1][6][35] = int_stack + 17652;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 17736;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 17796;
 Libderiv->deriv2_classes[2][5][35] = int_stack + 17886;
 Libderiv->deriv2_classes[2][6][35] = int_stack + 18012;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 18180;
 Libderiv->deriv2_classes[1][4][34] = int_stack + 18210;
 Libderiv->deriv2_classes[1][5][34] = int_stack + 18255;
 Libderiv->deriv2_classes[1][6][34] = int_stack + 18318;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 18402;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 18462;
 Libderiv->deriv2_classes[2][5][34] = int_stack + 18552;
 Libderiv->deriv2_classes[2][6][34] = int_stack + 18678;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 18846;
 Libderiv->deriv2_classes[1][4][33] = int_stack + 18876;
 Libderiv->deriv2_classes[1][5][33] = int_stack + 18921;
 Libderiv->deriv2_classes[1][6][33] = int_stack + 18984;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 19068;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 19128;
 Libderiv->deriv2_classes[2][5][33] = int_stack + 19218;
 Libderiv->deriv2_classes[2][6][33] = int_stack + 19344;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 19512;
 Libderiv->deriv2_classes[1][4][32] = int_stack + 19542;
 Libderiv->deriv2_classes[1][5][32] = int_stack + 19587;
 Libderiv->deriv2_classes[1][6][32] = int_stack + 19650;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 19734;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 19794;
 Libderiv->deriv2_classes[2][5][32] = int_stack + 19884;
 Libderiv->deriv2_classes[2][6][32] = int_stack + 20010;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 20178;
 Libderiv->deriv2_classes[1][4][31] = int_stack + 20208;
 Libderiv->deriv2_classes[1][5][31] = int_stack + 20253;
 Libderiv->deriv2_classes[1][6][31] = int_stack + 20316;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 20400;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 20460;
 Libderiv->deriv2_classes[2][5][31] = int_stack + 20550;
 Libderiv->deriv2_classes[2][6][31] = int_stack + 20676;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 20844;
 Libderiv->deriv2_classes[1][4][30] = int_stack + 20874;
 Libderiv->deriv2_classes[1][5][30] = int_stack + 20919;
 Libderiv->deriv2_classes[1][6][30] = int_stack + 20982;
 Libderiv->deriv_classes[2][3][2] = int_stack + 21066;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 21126;
 Libderiv->deriv_classes[2][4][2] = int_stack + 21186;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 21276;
 Libderiv->deriv_classes[2][5][2] = int_stack + 21366;
 Libderiv->deriv2_classes[2][5][30] = int_stack + 21492;
 Libderiv->deriv2_classes[2][6][30] = int_stack + 21618;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 21786;
 Libderiv->deriv2_classes[1][4][26] = int_stack + 21816;
 Libderiv->deriv2_classes[1][5][26] = int_stack + 21861;
 Libderiv->deriv2_classes[1][6][26] = int_stack + 21924;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 22008;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 22068;
 Libderiv->deriv2_classes[2][5][26] = int_stack + 22158;
 Libderiv->deriv2_classes[2][6][26] = int_stack + 22284;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 22452;
 Libderiv->deriv2_classes[1][4][23] = int_stack + 22482;
 Libderiv->deriv2_classes[1][5][23] = int_stack + 22527;
 Libderiv->deriv2_classes[1][6][23] = int_stack + 22590;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 22674;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 22734;
 Libderiv->deriv2_classes[2][5][23] = int_stack + 22824;
 Libderiv->deriv2_classes[2][6][23] = int_stack + 22950;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 23118;
 Libderiv->deriv2_classes[1][4][22] = int_stack + 23148;
 Libderiv->deriv2_classes[1][5][22] = int_stack + 23193;
 Libderiv->deriv2_classes[1][6][22] = int_stack + 23256;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 23340;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 23400;
 Libderiv->deriv2_classes[2][5][22] = int_stack + 23490;
 Libderiv->deriv2_classes[2][6][22] = int_stack + 23616;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 23784;
 Libderiv->deriv2_classes[1][4][21] = int_stack + 23814;
 Libderiv->deriv2_classes[1][5][21] = int_stack + 23859;
 Libderiv->deriv2_classes[1][6][21] = int_stack + 23922;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 24006;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 24066;
 Libderiv->deriv2_classes[2][5][21] = int_stack + 24156;
 Libderiv->deriv2_classes[2][6][21] = int_stack + 24282;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 24450;
 Libderiv->deriv2_classes[1][4][20] = int_stack + 24480;
 Libderiv->deriv2_classes[1][5][20] = int_stack + 24525;
 Libderiv->deriv2_classes[1][6][20] = int_stack + 24588;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 24672;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 24732;
 Libderiv->deriv2_classes[2][5][20] = int_stack + 24822;
 Libderiv->deriv2_classes[2][6][20] = int_stack + 24948;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 25116;
 Libderiv->deriv2_classes[1][4][19] = int_stack + 25146;
 Libderiv->deriv2_classes[1][5][19] = int_stack + 25191;
 Libderiv->deriv2_classes[1][6][19] = int_stack + 25254;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 25338;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 25398;
 Libderiv->deriv2_classes[2][5][19] = int_stack + 25488;
 Libderiv->deriv2_classes[2][6][19] = int_stack + 25614;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 25782;
 Libderiv->deriv2_classes[1][4][18] = int_stack + 25812;
 Libderiv->deriv2_classes[1][5][18] = int_stack + 25857;
 Libderiv->deriv2_classes[1][6][18] = int_stack + 25920;
 Libderiv->deriv_classes[2][3][1] = int_stack + 26004;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 26064;
 Libderiv->deriv_classes[2][4][1] = int_stack + 26124;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 26214;
 Libderiv->deriv_classes[2][5][1] = int_stack + 26304;
 Libderiv->deriv2_classes[2][5][18] = int_stack + 26430;
 Libderiv->deriv2_classes[2][6][18] = int_stack + 26556;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 26724;
 Libderiv->deriv2_classes[1][4][14] = int_stack + 26754;
 Libderiv->deriv2_classes[1][5][14] = int_stack + 26799;
 Libderiv->deriv2_classes[1][6][14] = int_stack + 26862;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 26946;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 27006;
 Libderiv->deriv2_classes[2][5][14] = int_stack + 27096;
 Libderiv->deriv2_classes[2][6][14] = int_stack + 27222;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 27390;
 Libderiv->deriv2_classes[1][4][13] = int_stack + 27420;
 Libderiv->deriv2_classes[1][5][13] = int_stack + 27465;
 Libderiv->deriv2_classes[1][6][13] = int_stack + 27528;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 27612;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 27672;
 Libderiv->deriv2_classes[2][5][13] = int_stack + 27762;
 Libderiv->deriv2_classes[2][6][13] = int_stack + 27888;
 Libderiv->deriv_classes[1][3][11] = int_stack + 28056;
 Libderiv->deriv_classes[1][4][11] = int_stack + 28086;
 Libderiv->deriv_classes[1][5][11] = int_stack + 28131;
 Libderiv->deriv_classes[1][6][11] = int_stack + 28194;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 28278;
 Libderiv->deriv2_classes[1][4][11] = int_stack + 28308;
 Libderiv->deriv2_classes[1][5][11] = int_stack + 28353;
 Libderiv->deriv2_classes[1][6][11] = int_stack + 28416;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 28500;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 28560;
 Libderiv->deriv2_classes[2][5][11] = int_stack + 28650;
 Libderiv->deriv2_classes[2][6][11] = int_stack + 28776;
 Libderiv->deriv_classes[1][3][10] = int_stack + 28944;
 Libderiv->deriv_classes[1][4][10] = int_stack + 28974;
 Libderiv->deriv_classes[1][5][10] = int_stack + 29019;
 Libderiv->deriv_classes[1][6][10] = int_stack + 29082;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 29166;
 Libderiv->deriv2_classes[1][4][10] = int_stack + 29196;
 Libderiv->deriv2_classes[1][5][10] = int_stack + 29241;
 Libderiv->deriv2_classes[1][6][10] = int_stack + 29304;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 29388;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 29448;
 Libderiv->deriv2_classes[2][5][10] = int_stack + 29538;
 Libderiv->deriv2_classes[2][6][10] = int_stack + 29664;
 Libderiv->deriv_classes[1][3][9] = int_stack + 29832;
 Libderiv->deriv_classes[1][4][9] = int_stack + 29862;
 Libderiv->deriv_classes[1][5][9] = int_stack + 29907;
 Libderiv->deriv_classes[1][6][9] = int_stack + 29970;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 30054;
 Libderiv->deriv2_classes[1][4][9] = int_stack + 30084;
 Libderiv->deriv2_classes[1][5][9] = int_stack + 30129;
 Libderiv->deriv2_classes[1][6][9] = int_stack + 30192;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 30276;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 30336;
 Libderiv->deriv2_classes[2][5][9] = int_stack + 30426;
 Libderiv->deriv2_classes[2][6][9] = int_stack + 30552;
 Libderiv->deriv_classes[1][3][8] = int_stack + 30720;
 Libderiv->deriv_classes[1][4][8] = int_stack + 30750;
 Libderiv->deriv_classes[1][5][8] = int_stack + 30795;
 Libderiv->deriv_classes[1][6][8] = int_stack + 30858;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 30942;
 Libderiv->deriv2_classes[1][4][8] = int_stack + 30972;
 Libderiv->deriv2_classes[1][5][8] = int_stack + 31017;
 Libderiv->deriv2_classes[1][6][8] = int_stack + 31080;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 31164;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 31224;
 Libderiv->deriv2_classes[2][5][8] = int_stack + 31314;
 Libderiv->deriv2_classes[2][6][8] = int_stack + 31440;
 Libderiv->deriv_classes[1][3][7] = int_stack + 31608;
 Libderiv->deriv_classes[1][4][7] = int_stack + 31638;
 Libderiv->deriv_classes[1][5][7] = int_stack + 31683;
 Libderiv->deriv_classes[1][6][7] = int_stack + 31746;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 31830;
 Libderiv->deriv2_classes[1][4][7] = int_stack + 31860;
 Libderiv->deriv2_classes[1][5][7] = int_stack + 31905;
 Libderiv->deriv2_classes[1][6][7] = int_stack + 31968;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 32052;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 32112;
 Libderiv->deriv2_classes[2][5][7] = int_stack + 32202;
 Libderiv->deriv2_classes[2][6][7] = int_stack + 32328;
 Libderiv->dvrr_classes[1][3] = int_stack + 32496;
 Libderiv->deriv_classes[1][3][6] = int_stack + 32526;
 Libderiv->dvrr_classes[1][4] = int_stack + 32556;
 Libderiv->deriv_classes[1][4][6] = int_stack + 32601;
 Libderiv->dvrr_classes[1][5] = int_stack + 32646;
 Libderiv->deriv_classes[1][5][6] = int_stack + 32709;
 Libderiv->deriv_classes[1][6][6] = int_stack + 32772;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 32856;
 Libderiv->deriv2_classes[1][4][6] = int_stack + 32886;
 Libderiv->deriv2_classes[1][5][6] = int_stack + 32931;
 Libderiv->deriv2_classes[1][6][6] = int_stack + 32994;
 Libderiv->deriv_classes[2][3][0] = int_stack + 33078;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 33138;
 Libderiv->deriv_classes[2][4][0] = int_stack + 33198;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 33288;
 Libderiv->deriv_classes[2][5][0] = int_stack + 33378;
 Libderiv->deriv2_classes[2][5][6] = int_stack + 33504;
 Libderiv->deriv2_classes[2][6][6] = int_stack + 33630;
 Libderiv->deriv_classes[1][3][2] = int_stack + 33798;
 Libderiv->deriv_classes[1][4][2] = int_stack + 33828;
 Libderiv->deriv_classes[1][5][2] = int_stack + 33873;
 Libderiv->deriv_classes[1][6][2] = int_stack + 33936;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 34020;
 Libderiv->deriv2_classes[1][4][2] = int_stack + 34050;
 Libderiv->deriv2_classes[1][5][2] = int_stack + 34095;
 Libderiv->deriv2_classes[1][6][2] = int_stack + 34158;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 34242;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 34302;
 Libderiv->deriv2_classes[2][5][2] = int_stack + 34392;
 Libderiv->deriv2_classes[2][6][2] = int_stack + 34518;
 Libderiv->deriv_classes[1][3][1] = int_stack + 34686;
 Libderiv->deriv_classes[1][4][1] = int_stack + 34716;
 Libderiv->deriv_classes[1][5][1] = int_stack + 34761;
 Libderiv->deriv_classes[1][6][1] = int_stack + 34824;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 34908;
 Libderiv->deriv2_classes[1][4][1] = int_stack + 34938;
 Libderiv->deriv2_classes[1][5][1] = int_stack + 34983;
 Libderiv->deriv2_classes[1][6][1] = int_stack + 35046;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 35130;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 35190;
 Libderiv->deriv2_classes[2][5][1] = int_stack + 35280;
 Libderiv->deriv2_classes[2][6][1] = int_stack + 35406;
 Libderiv->deriv_classes[1][3][0] = int_stack + 35574;
 Libderiv->deriv_classes[1][4][0] = int_stack + 35604;
 Libderiv->deriv_classes[1][5][0] = int_stack + 35649;
 Libderiv->deriv_classes[1][6][0] = int_stack + 35712;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 35796;
 Libderiv->deriv2_classes[1][4][0] = int_stack + 35826;
 Libderiv->deriv2_classes[1][5][0] = int_stack + 35871;
 Libderiv->deriv2_classes[1][6][0] = int_stack + 35934;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 36018;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 36078;
 Libderiv->deriv2_classes[2][5][0] = int_stack + 36168;
 Libderiv->deriv2_classes[2][6][0] = int_stack + 36294;
 memset(int_stack,0,291696);

 Libderiv->dvrr_stack = int_stack + 64428;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ppff(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+36462,int_stack+32556,int_stack+32496,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+36552,int_stack+32646,int_stack+32556,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+36687,int_stack+36552,int_stack+36462,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36867,int_stack+28086,int_stack+28056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32496,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36957,int_stack+28131,int_stack+28086, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32556,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+37092,int_stack+36957,int_stack+36867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36462,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+37272,int_stack+28194,int_stack+28131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32646,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+37461,int_stack+37272,int_stack+36957, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36552,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+37731,int_stack+37461,int_stack+37092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36687,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+37272,int_stack+16824,int_stack+16644,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+37452,int_stack+840,int_stack+16824,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+38031,int_stack+37452,int_stack+37272,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+38391,int_stack+12054,int_stack+11934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16644,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+38571,int_stack+12234,int_stack+12054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16824,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+38841,int_stack+38571,int_stack+38391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37272,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+39201,int_stack+0,int_stack+12234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+39579,int_stack+39201,int_stack+38571, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37452,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+40119,int_stack+39579,int_stack+38841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38031,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+39201,int_stack+28974,int_stack+28944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32496, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+39291,int_stack+29019,int_stack+28974, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32556, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+39426,int_stack+39291,int_stack+39201, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36462, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+39606,int_stack+29082,int_stack+29019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32646, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+39795,int_stack+39606,int_stack+39291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36552, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+40719,int_stack+39795,int_stack+39426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36687, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+39606,int_stack+12996,int_stack+12876, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16644, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+39786,int_stack+13176,int_stack+12996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16824, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41019,int_stack+39786,int_stack+39606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37272, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41379,int_stack+168,int_stack+13176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41757,int_stack+41379,int_stack+39786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37452, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+42297,int_stack+41757,int_stack+41019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38031, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41379,int_stack+29862,int_stack+29832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32496, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41469,int_stack+29907,int_stack+29862, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32556, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41604,int_stack+41469,int_stack+41379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36462, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41784,int_stack+29970,int_stack+29907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32646, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41973,int_stack+41784,int_stack+41469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36552, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+41973,int_stack+41604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36687, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41784,int_stack+13938,int_stack+13818, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16644, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41964,int_stack+14118,int_stack+13938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16824, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42897,int_stack+41964,int_stack+41784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37272, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43257,int_stack+336,int_stack+14118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+43635,int_stack+43257,int_stack+41964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37452, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+44175,int_stack+43635,int_stack+42897, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38031, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43257,int_stack+30750,int_stack+30720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43347,int_stack+30795,int_stack+30750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+43482,int_stack+43347,int_stack+43257, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43662,int_stack+30858,int_stack+30795, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+43851,int_stack+43662,int_stack+43347, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+44775,int_stack+43851,int_stack+43482, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43662,int_stack+14880,int_stack+14760, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+43842,int_stack+15060,int_stack+14880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45075,int_stack+43842,int_stack+43662, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45435,int_stack+504,int_stack+15060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+45813,int_stack+45435,int_stack+43842, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+46353,int_stack+45813,int_stack+45075, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45435,int_stack+31638,int_stack+31608, 0.0, zero_stack, 1.0, int_stack+32496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45525,int_stack+31683,int_stack+31638, 0.0, zero_stack, 1.0, int_stack+32556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45660,int_stack+45525,int_stack+45435, 0.0, zero_stack, 1.0, int_stack+36462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45840,int_stack+31746,int_stack+31683, 0.0, zero_stack, 1.0, int_stack+32646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+46029,int_stack+45840,int_stack+45525, 0.0, zero_stack, 1.0, int_stack+36552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+300,int_stack+46029,int_stack+45660, 0.0, zero_stack, 1.0, int_stack+36687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45840,int_stack+15822,int_stack+15702, 0.0, zero_stack, 1.0, int_stack+16644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+46020,int_stack+16002,int_stack+15822, 0.0, zero_stack, 1.0, int_stack+16824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+46953,int_stack+46020,int_stack+45840, 0.0, zero_stack, 1.0, int_stack+37272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47313,int_stack+672,int_stack+16002, 0.0, zero_stack, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+47691,int_stack+47313,int_stack+46020, 0.0, zero_stack, 1.0, int_stack+37452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+48231,int_stack+47691,int_stack+46953, 0.0, zero_stack, 1.0, int_stack+38031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47313,int_stack+32601,int_stack+32526, 1.0, int_stack+32496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47403,int_stack+32709,int_stack+32601, 1.0, int_stack+32556, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+47538,int_stack+47403,int_stack+47313, 1.0, int_stack+36462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47718,int_stack+32772,int_stack+32709, 1.0, int_stack+32646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+47907,int_stack+47718,int_stack+47403, 1.0, int_stack+36552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+48831,int_stack+47907,int_stack+47538, 1.0, int_stack+36687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47718,int_stack+16914,int_stack+16704, 1.0, int_stack+16644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47898,int_stack+17094,int_stack+16914, 1.0, int_stack+16824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+49131,int_stack+47898,int_stack+47718, 1.0, int_stack+37272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49491,int_stack+966,int_stack+17094, 1.0, int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+49869,int_stack+49491,int_stack+47898, 1.0, int_stack+37452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+50409,int_stack+49869,int_stack+49131, 1.0, int_stack+38031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+38031,int_stack+1470,int_stack+32646,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+49491,int_stack+38031,int_stack+36552,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+38031,int_stack+49491,int_stack+36687,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16824,int_stack+33828,int_stack+33798,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+49491,int_stack+33873,int_stack+33828,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+49626,int_stack+49491,int_stack+16824,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+49806,int_stack+33936,int_stack+33873,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+49995,int_stack+49806,int_stack+49491,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+37272,int_stack+49995,int_stack+49626,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+49806,int_stack+21186,int_stack+21066,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+49986,int_stack+21366,int_stack+21186,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+36462,int_stack+49986,int_stack+49806,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+600,int_stack+1134,int_stack+21366,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+51009,int_stack+600,int_stack+49986,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+600,int_stack+51009,int_stack+36462,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51009,int_stack+34716,int_stack+34686,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51099,int_stack+34761,int_stack+34716,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+51234,int_stack+51099,int_stack+51009,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+51414,int_stack+34824,int_stack+34761,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+51603,int_stack+51414,int_stack+51099,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+51873,int_stack+51603,int_stack+51234,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51414,int_stack+26124,int_stack+26004,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+51594,int_stack+26304,int_stack+26124,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+52173,int_stack+51594,int_stack+51414,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52533,int_stack+1302,int_stack+26304,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+52911,int_stack+52533,int_stack+51594,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+53451,int_stack+52911,int_stack+52173,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52533,int_stack+35604,int_stack+35574,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52623,int_stack+35649,int_stack+35604,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+52758,int_stack+52623,int_stack+52533,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+52938,int_stack+35712,int_stack+35649,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+53127,int_stack+52938,int_stack+52623,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+54051,int_stack+53127,int_stack+52758,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52938,int_stack+33198,int_stack+33078,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+53118,int_stack+33378,int_stack+33198,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+54351,int_stack+53118,int_stack+52938,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+54711,int_stack+1554,int_stack+33378,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+55089,int_stack+54711,int_stack+53118,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+55629,int_stack+55089,int_stack+54351,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+54711,int_stack+1752,int_stack+1722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28056,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+54801,int_stack+1797,int_stack+1752, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28086,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+54936,int_stack+54801,int_stack+54711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+36867,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55116,int_stack+1860,int_stack+1797, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28131,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55305,int_stack+55116,int_stack+54801, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+36957,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+56229,int_stack+55305,int_stack+54936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+37092,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+54711,int_stack+2004,int_stack+1944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11934,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+54891,int_stack+2094,int_stack+2004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12054,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+55161,int_stack+54891,int_stack+54711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+38391,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+56529,int_stack+2220,int_stack+2094, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12234,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1200,int_stack+56529,int_stack+54891, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+38571,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1740,int_stack+1200,int_stack+55161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+38841,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+2418,int_stack+2388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28056, 1.0, int_stack+28944,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1290,int_stack+2463,int_stack+2418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28086, 1.0, int_stack+28974,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1425,int_stack+1290,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36867, 1.0, int_stack+39201,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+56529,int_stack+2526,int_stack+2463, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28131, 1.0, int_stack+29019,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2340,int_stack+56529,int_stack+1290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36957, 1.0, int_stack+39291,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+56529,int_stack+2340,int_stack+1425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37092, 1.0, int_stack+39426,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2340,int_stack+2670,int_stack+2610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11934, 1.0, int_stack+12876,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1200,int_stack+2760,int_stack+2670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12054, 1.0, int_stack+12996,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+54711,int_stack+1200,int_stack+2340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38391, 1.0, int_stack+39606,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2340,int_stack+2886,int_stack+2760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12234, 1.0, int_stack+13176,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55071,int_stack+2340,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38571, 1.0, int_stack+39786,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2340,int_stack+55071,int_stack+54711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38841, 1.0, int_stack+41019,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+54711,int_stack+3084,int_stack+3054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28944, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+54801,int_stack+3129,int_stack+3084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28974, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+54936,int_stack+54801,int_stack+54711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+39201, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2940,int_stack+3192,int_stack+3129, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29019, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+55116,int_stack+2940,int_stack+54801, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+39291, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2940,int_stack+55116,int_stack+54936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+39426, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+54711,int_stack+3336,int_stack+3276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12876, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+54891,int_stack+3426,int_stack+3336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12996, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+55161,int_stack+54891,int_stack+54711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+39606, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1200,int_stack+3552,int_stack+3426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13176, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56829,int_stack+1200,int_stack+54891, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+39786, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+57369,int_stack+56829,int_stack+55161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41019, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+3750,int_stack+3720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28056, 0.0, zero_stack, 1.0, int_stack+29832,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+3795,int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28086, 0.0, zero_stack, 1.0, int_stack+29862,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36867, 0.0, zero_stack, 1.0, int_stack+41379,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1200,int_stack+3858,int_stack+3795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28131, 0.0, zero_stack, 1.0, int_stack+29907,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1389,int_stack+1200,int_stack+56919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36957, 0.0, zero_stack, 1.0, int_stack+41469,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+54711,int_stack+1389,int_stack+57054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37092, 0.0, zero_stack, 1.0, int_stack+41604,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+4002,int_stack+3942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11934, 0.0, zero_stack, 1.0, int_stack+13818,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1380,int_stack+4092,int_stack+4002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12054, 0.0, zero_stack, 1.0, int_stack+13938,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+1380,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38391, 0.0, zero_stack, 1.0, int_stack+41784,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55011,int_stack+4218,int_stack+4092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12234, 0.0, zero_stack, 1.0, int_stack+14118,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3240,int_stack+55011,int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38571, 0.0, zero_stack, 1.0, int_stack+41964,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+55011,int_stack+3240,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38841, 0.0, zero_stack, 1.0, int_stack+42897,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+4416,int_stack+4386, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28944, 1.0, int_stack+29832, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+4461,int_stack+4416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28974, 1.0, int_stack+29862, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39201, 1.0, int_stack+41379, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3240,int_stack+4524,int_stack+4461, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29019, 1.0, int_stack+29907, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3429,int_stack+3240,int_stack+56919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39291, 1.0, int_stack+41469, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3699,int_stack+3429,int_stack+57054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39426, 1.0, int_stack+41604, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3240,int_stack+4668,int_stack+4608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12876, 1.0, int_stack+13818, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3420,int_stack+4758,int_stack+4668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12996, 1.0, int_stack+13938, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+3420,int_stack+3240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39606, 1.0, int_stack+41784, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+3999,int_stack+4884,int_stack+4758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13176, 1.0, int_stack+14118, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1200,int_stack+3999,int_stack+3420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39786, 1.0, int_stack+41964, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3999,int_stack+1200,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41019, 1.0, int_stack+42897, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+5082,int_stack+5052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29832, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+5127,int_stack+5082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29862, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41379, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1200,int_stack+5190,int_stack+5127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29907, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1389,int_stack+1200,int_stack+56919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41469, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4599,int_stack+1389,int_stack+57054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41604, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+5334,int_stack+5274, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13818, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1380,int_stack+5424,int_stack+5334, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13938, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+1380,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41784, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4899,int_stack+5550,int_stack+5424, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14118, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+57969,int_stack+4899,int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41964, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4899,int_stack+57969,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+42897, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+5748,int_stack+5718, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28056, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30720,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+5793,int_stack+5748, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28086, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30750,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36867, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43257,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+57969,int_stack+5856,int_stack+5793, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28131, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30795,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+58158,int_stack+57969,int_stack+56919, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36957, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43347,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+58428,int_stack+58158,int_stack+57054, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37092, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43482,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+57969,int_stack+6000,int_stack+5940, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11934, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14760,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+58149,int_stack+6090,int_stack+6000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12054, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14880,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+58149,int_stack+57969, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38391, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43662,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5499,int_stack+6216,int_stack+6090, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12234, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15060,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1200,int_stack+5499,int_stack+58149, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38571, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43842,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5499,int_stack+1200,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38841, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45075,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+6414,int_stack+6384, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28944, 0.0, zero_stack, 1.0, int_stack+30720, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+6459,int_stack+6414, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28974, 0.0, zero_stack, 1.0, int_stack+30750, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39201, 0.0, zero_stack, 1.0, int_stack+43257, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1200,int_stack+6522,int_stack+6459, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29019, 0.0, zero_stack, 1.0, int_stack+30795, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1389,int_stack+1200,int_stack+56919, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39291, 0.0, zero_stack, 1.0, int_stack+43347, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6099,int_stack+1389,int_stack+57054, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39426, 0.0, zero_stack, 1.0, int_stack+43482, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+6666,int_stack+6606, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12876, 0.0, zero_stack, 1.0, int_stack+14760, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1380,int_stack+6756,int_stack+6666, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12996, 0.0, zero_stack, 1.0, int_stack+14880, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+1380,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39606, 0.0, zero_stack, 1.0, int_stack+43662, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+57969,int_stack+6882,int_stack+6756, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13176, 0.0, zero_stack, 1.0, int_stack+15060, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6399,int_stack+57969,int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39786, 0.0, zero_stack, 1.0, int_stack+43842, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+58728,int_stack+6399,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41019, 0.0, zero_stack, 1.0, int_stack+45075, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+7080,int_stack+7050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29832, 1.0, int_stack+30720, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+7125,int_stack+7080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29862, 1.0, int_stack+30750, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41379, 1.0, int_stack+43257, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6399,int_stack+7188,int_stack+7125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29907, 1.0, int_stack+30795, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6588,int_stack+6399,int_stack+56919, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41469, 1.0, int_stack+43347, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6858,int_stack+6588,int_stack+57054, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41604, 1.0, int_stack+43482, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6399,int_stack+7332,int_stack+7272, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13818, 1.0, int_stack+14760, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6579,int_stack+7422,int_stack+7332, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13938, 1.0, int_stack+14880, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+6579,int_stack+6399, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41784, 1.0, int_stack+43662, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+57969,int_stack+7548,int_stack+7422, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14118, 1.0, int_stack+15060, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1200,int_stack+57969,int_stack+6579, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41964, 1.0, int_stack+43842, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+59328,int_stack+1200,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42897, 1.0, int_stack+45075, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+7746,int_stack+7716, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+30720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+7791,int_stack+7746, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+30750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43257, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1200,int_stack+7854,int_stack+7791, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+30795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1389,int_stack+1200,int_stack+56919, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43347, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+57969,int_stack+1389,int_stack+57054, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+7998,int_stack+7938, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1380,int_stack+8088,int_stack+7998, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+1380,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43662, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6399,int_stack+8214,int_stack+8088, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+6399,int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7698,int_stack+7158,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+45075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+8412,int_stack+8382, 0.0, zero_stack, 1.0, int_stack+28056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31608,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+8457,int_stack+8412, 0.0, zero_stack, 1.0, int_stack+28086, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31638,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+36867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45435,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7158,int_stack+8520,int_stack+8457, 0.0, zero_stack, 1.0, int_stack+28131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31683,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7347,int_stack+7158,int_stack+56919, 0.0, zero_stack, 1.0, int_stack+36957, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45525,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6399,int_stack+7347,int_stack+57054, 0.0, zero_stack, 1.0, int_stack+37092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45660,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7158,int_stack+8664,int_stack+8604, 0.0, zero_stack, 1.0, int_stack+11934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15702,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7338,int_stack+8754,int_stack+8664, 0.0, zero_stack, 1.0, int_stack+12054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15822,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+7338,int_stack+7158, 0.0, zero_stack, 1.0, int_stack+38391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45840,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1200,int_stack+8880,int_stack+8754, 0.0, zero_stack, 1.0, int_stack+12234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16002,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8298,int_stack+1200,int_stack+7338, 0.0, zero_stack, 1.0, int_stack+38571, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46020,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+59928,int_stack+8298,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+38841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46953,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+9078,int_stack+9048, 0.0, zero_stack, 1.0, int_stack+28944, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31608, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+9123,int_stack+9078, 0.0, zero_stack, 1.0, int_stack+28974, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31638, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+39201, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45435, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8298,int_stack+9186,int_stack+9123, 0.0, zero_stack, 1.0, int_stack+29019, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31683, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8487,int_stack+8298,int_stack+56919, 0.0, zero_stack, 1.0, int_stack+39291, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45525, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8757,int_stack+8487,int_stack+57054, 0.0, zero_stack, 1.0, int_stack+39426, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45660, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8298,int_stack+9330,int_stack+9270, 0.0, zero_stack, 1.0, int_stack+12876, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15702, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8478,int_stack+9420,int_stack+9330, 0.0, zero_stack, 1.0, int_stack+12996, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15822, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+8478,int_stack+8298, 0.0, zero_stack, 1.0, int_stack+39606, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45840, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1200,int_stack+9546,int_stack+9420, 0.0, zero_stack, 1.0, int_stack+13176, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16002, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+1200,int_stack+8478, 0.0, zero_stack, 1.0, int_stack+39786, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46020, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9057,int_stack+7158,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+41019, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46953, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+9744,int_stack+9714, 0.0, zero_stack, 1.0, int_stack+29832, 0.0, zero_stack, 1.0, int_stack+31608, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+9789,int_stack+9744, 0.0, zero_stack, 1.0, int_stack+29862, 0.0, zero_stack, 1.0, int_stack+31638, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+41379, 0.0, zero_stack, 1.0, int_stack+45435, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7158,int_stack+9852,int_stack+9789, 0.0, zero_stack, 1.0, int_stack+29907, 0.0, zero_stack, 1.0, int_stack+31683, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7347,int_stack+7158,int_stack+56919, 0.0, zero_stack, 1.0, int_stack+41469, 0.0, zero_stack, 1.0, int_stack+45525, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1200,int_stack+7347,int_stack+57054, 0.0, zero_stack, 1.0, int_stack+41604, 0.0, zero_stack, 1.0, int_stack+45660, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7158,int_stack+9996,int_stack+9936, 0.0, zero_stack, 1.0, int_stack+13818, 0.0, zero_stack, 1.0, int_stack+15702, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7338,int_stack+10086,int_stack+9996, 0.0, zero_stack, 1.0, int_stack+13938, 0.0, zero_stack, 1.0, int_stack+15822, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+7338,int_stack+7158, 0.0, zero_stack, 1.0, int_stack+41784, 0.0, zero_stack, 1.0, int_stack+45840, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8298,int_stack+10212,int_stack+10086, 0.0, zero_stack, 1.0, int_stack+14118, 0.0, zero_stack, 1.0, int_stack+16002, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9657,int_stack+8298,int_stack+7338, 0.0, zero_stack, 1.0, int_stack+41964, 0.0, zero_stack, 1.0, int_stack+46020, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+60528,int_stack+9657,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+42897, 0.0, zero_stack, 1.0, int_stack+46953, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+10410,int_stack+10380, 0.0, zero_stack, 1.0, int_stack+30720, 1.0, int_stack+31608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+10455,int_stack+10410, 0.0, zero_stack, 1.0, int_stack+30750, 1.0, int_stack+31638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+43257, 1.0, int_stack+45435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9657,int_stack+10518,int_stack+10455, 0.0, zero_stack, 1.0, int_stack+30795, 1.0, int_stack+31683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9846,int_stack+9657,int_stack+56919, 0.0, zero_stack, 1.0, int_stack+43347, 1.0, int_stack+45525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10116,int_stack+9846,int_stack+57054, 0.0, zero_stack, 1.0, int_stack+43482, 1.0, int_stack+45660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9657,int_stack+10662,int_stack+10602, 0.0, zero_stack, 1.0, int_stack+14760, 1.0, int_stack+15702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9837,int_stack+10752,int_stack+10662, 0.0, zero_stack, 1.0, int_stack+14880, 1.0, int_stack+15822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+9837,int_stack+9657, 0.0, zero_stack, 1.0, int_stack+43662, 1.0, int_stack+45840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8298,int_stack+10878,int_stack+10752, 0.0, zero_stack, 1.0, int_stack+15060, 1.0, int_stack+16002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+8298,int_stack+9837, 0.0, zero_stack, 1.0, int_stack+43842, 1.0, int_stack+46020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+10416,int_stack+7158,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+45075, 1.0, int_stack+46953, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+11076,int_stack+11046, 0.0, zero_stack, 2.0, int_stack+31608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+11121,int_stack+11076, 0.0, zero_stack, 2.0, int_stack+31638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 0.0, zero_stack, 2.0, int_stack+45435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7158,int_stack+11184,int_stack+11121, 0.0, zero_stack, 2.0, int_stack+31683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7347,int_stack+7158,int_stack+56919, 0.0, zero_stack, 2.0, int_stack+45525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8298,int_stack+7347,int_stack+57054, 0.0, zero_stack, 2.0, int_stack+45660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7158,int_stack+11328,int_stack+11268, 0.0, zero_stack, 2.0, int_stack+15702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7338,int_stack+11418,int_stack+11328, 0.0, zero_stack, 2.0, int_stack+15822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+7338,int_stack+7158, 0.0, zero_stack, 2.0, int_stack+45840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9657,int_stack+11544,int_stack+11418, 0.0, zero_stack, 2.0, int_stack+16002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11016,int_stack+9657,int_stack+7338, 0.0, zero_stack, 2.0, int_stack+46020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+61128,int_stack+11016,int_stack+56829, 0.0, zero_stack, 2.0, int_stack+46953, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+11742,int_stack+11712, 1.0, int_stack+28056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32526,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+11787,int_stack+11742, 1.0, int_stack+28086, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32601,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 1.0, int_stack+36867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47313,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11016,int_stack+11850,int_stack+11787, 1.0, int_stack+28131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32709,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11205,int_stack+11016,int_stack+56919, 1.0, int_stack+36957, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47403,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11475,int_stack+11205,int_stack+57054, 1.0, int_stack+37092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47538,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11016,int_stack+12144,int_stack+11994, 1.0, int_stack+11934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16704,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11196,int_stack+12360,int_stack+12144, 1.0, int_stack+12054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16914,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+11196,int_stack+11016, 1.0, int_stack+38391, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47718,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11775,int_stack+12486,int_stack+12360, 1.0, int_stack+12234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17094,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+11775,int_stack+11196, 1.0, int_stack+38571, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47898,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11775,int_stack+7158,int_stack+56829, 1.0, int_stack+38841, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49131,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+12684,int_stack+12654, 1.0, int_stack+28944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32526, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56919,int_stack+12729,int_stack+12684, 1.0, int_stack+28974, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32601, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+57054,int_stack+56919,int_stack+56829, 1.0, int_stack+39201, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47313, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7158,int_stack+12792,int_stack+12729, 1.0, int_stack+29019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32709, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7347,int_stack+7158,int_stack+56919, 1.0, int_stack+39291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47403, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12375,int_stack+7347,int_stack+57054, 1.0, int_stack+39426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47538, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7158,int_stack+13086,int_stack+12936, 1.0, int_stack+12876, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16704, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7338,int_stack+13302,int_stack+13086, 1.0, int_stack+12996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16914, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+7338,int_stack+7158, 1.0, int_stack+39606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47718, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12675,int_stack+13428,int_stack+13302, 1.0, int_stack+13176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17094, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13053,int_stack+12675,int_stack+7338, 1.0, int_stack+39786, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47898, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+38331,int_stack+13053,int_stack+56829, 1.0, int_stack+41019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49131, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41019,int_stack+13626,int_stack+13596, 1.0, int_stack+29832, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32526, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41109,int_stack+13671,int_stack+13626, 1.0, int_stack+29862, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32601, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+41109,int_stack+41019, 1.0, int_stack+41379, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47313, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41244,int_stack+13734,int_stack+13671, 1.0, int_stack+29907, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32709, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+57009,int_stack+41244,int_stack+41109, 1.0, int_stack+41469, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47403, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+41019,int_stack+57009,int_stack+56829, 1.0, int_stack+41604, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47538, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+14028,int_stack+13878, 1.0, int_stack+13818, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16704, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+57009,int_stack+14244,int_stack+14028, 1.0, int_stack+13938, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16914, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41319,int_stack+57009,int_stack+56829, 1.0, int_stack+41784, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47718, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12675,int_stack+14370,int_stack+14244, 1.0, int_stack+14118, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17094, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+12675,int_stack+57009, 1.0, int_stack+41964, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47898, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12675,int_stack+7158,int_stack+41319, 1.0, int_stack+42897, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49131, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42897,int_stack+14568,int_stack+14538, 1.0, int_stack+30720, 0.0, zero_stack, 1.0, int_stack+32526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+42987,int_stack+14613,int_stack+14568, 1.0, int_stack+30750, 0.0, zero_stack, 1.0, int_stack+32601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41319,int_stack+42987,int_stack+42897, 1.0, int_stack+43257, 0.0, zero_stack, 1.0, int_stack+47313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43122,int_stack+14676,int_stack+14613, 1.0, int_stack+30795, 0.0, zero_stack, 1.0, int_stack+32709, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41499,int_stack+43122,int_stack+42987, 1.0, int_stack+43347, 0.0, zero_stack, 1.0, int_stack+47403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+42897,int_stack+41499,int_stack+41319, 1.0, int_stack+43482, 0.0, zero_stack, 1.0, int_stack+47538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41319,int_stack+14970,int_stack+14820, 1.0, int_stack+14760, 0.0, zero_stack, 1.0, int_stack+16704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41499,int_stack+15186,int_stack+14970, 1.0, int_stack+14880, 0.0, zero_stack, 1.0, int_stack+16914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41769,int_stack+41499,int_stack+41319, 1.0, int_stack+43662, 0.0, zero_stack, 1.0, int_stack+47718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43197,int_stack+15312,int_stack+15186, 1.0, int_stack+15060, 0.0, zero_stack, 1.0, int_stack+17094, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+43197,int_stack+41499, 1.0, int_stack+43842, 0.0, zero_stack, 1.0, int_stack+47898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+43197,int_stack+7158,int_stack+41769, 1.0, int_stack+45075, 0.0, zero_stack, 1.0, int_stack+49131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45075,int_stack+15510,int_stack+15480, 1.0, int_stack+31608, 1.0, int_stack+32526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45165,int_stack+15555,int_stack+15510, 1.0, int_stack+31638, 1.0, int_stack+32601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7158,int_stack+45165,int_stack+45075, 1.0, int_stack+45435, 1.0, int_stack+47313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+45300,int_stack+15618,int_stack+15555, 1.0, int_stack+31683, 1.0, int_stack+32709, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7338,int_stack+45300,int_stack+45165, 1.0, int_stack+45525, 1.0, int_stack+47403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+45075,int_stack+7338,int_stack+7158, 1.0, int_stack+45660, 1.0, int_stack+47538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7158,int_stack+15912,int_stack+15762, 1.0, int_stack+15702, 1.0, int_stack+16704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7338,int_stack+16128,int_stack+15912, 1.0, int_stack+15822, 1.0, int_stack+16914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45375,int_stack+7338,int_stack+7158, 1.0, int_stack+45840, 1.0, int_stack+47718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+16254,int_stack+16128, 1.0, int_stack+16002, 1.0, int_stack+17094, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56829,int_stack+43797,int_stack+7338, 1.0, int_stack+46020, 1.0, int_stack+47898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+45735,int_stack+56829,int_stack+45375, 1.0, int_stack+46953, 1.0, int_stack+49131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+46953,int_stack+16452,int_stack+16422, 2.0, int_stack+32526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+47043,int_stack+16497,int_stack+16452, 2.0, int_stack+32601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45375,int_stack+47043,int_stack+46953, 2.0, int_stack+47313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+47178,int_stack+16560,int_stack+16497, 2.0, int_stack+32709, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56829,int_stack+47178,int_stack+47043, 2.0, int_stack+47403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+46953,int_stack+56829,int_stack+45375, 2.0, int_stack+47538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45375,int_stack+17004,int_stack+16764, 2.0, int_stack+16704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+56829,int_stack+17220,int_stack+17004, 2.0, int_stack+16914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32496,int_stack+56829,int_stack+45375, 2.0, int_stack+47718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+17346,int_stack+17220, 2.0, int_stack+17094, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+43797,int_stack+56829, 2.0, int_stack+47898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16914,int_stack+7158,int_stack+32496, 2.0, int_stack+49131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49131,int_stack+17544,int_stack+17514, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33798,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49221,int_stack+17589,int_stack+17544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33828,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32496,int_stack+49221,int_stack+49131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16824,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7158,int_stack+17652,int_stack+17589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33873,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7347,int_stack+7158,int_stack+49221, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49491,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+49131,int_stack+7347,int_stack+32496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49626,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+32496,int_stack+17796,int_stack+17736, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21066,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7158,int_stack+17886,int_stack+17796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21186,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45375,int_stack+7158,int_stack+32496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49806,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+18012,int_stack+17886, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21366,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56829,int_stack+43797,int_stack+7158, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49986,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+17514,int_stack+56829,int_stack+45375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36462,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45375,int_stack+18210,int_stack+18180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33798, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45465,int_stack+18255,int_stack+18210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33828, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+45465,int_stack+45375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16824, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+57009,int_stack+18318,int_stack+18255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33873, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+57009,int_stack+45465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49491, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+57009,int_stack+7158,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49626, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+18462,int_stack+18402, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21066, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7158,int_stack+18552,int_stack+18462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21186, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45375,int_stack+7158,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49806, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+18678,int_stack+18552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21366, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18114,int_stack+43797,int_stack+7158, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49986, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+47253,int_stack+18114,int_stack+45375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36462, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45375,int_stack+18876,int_stack+18846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33798, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+45465,int_stack+18921,int_stack+18876, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33828, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+45465,int_stack+45375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16824, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18114,int_stack+18984,int_stack+18921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33873, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18303,int_stack+18114,int_stack+45465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49491, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+45375,int_stack+18303,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49626, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+19128,int_stack+19068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21066, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18114,int_stack+19218,int_stack+19128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21186, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32496,int_stack+18114,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49806, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+19344,int_stack+19218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21366, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+43797,int_stack+18114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49986, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18114,int_stack+7158,int_stack+32496, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36462, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+32496,int_stack+19542,int_stack+19512, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32586,int_stack+19587,int_stack+19542, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+32586,int_stack+32496, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7158,int_stack+19650,int_stack+19587, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33873, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7347,int_stack+7158,int_stack+32586, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+32496,int_stack+7347,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+19794,int_stack+19734, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21066, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7158,int_stack+19884,int_stack+19794, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18714,int_stack+7158,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+20010,int_stack+19884, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21366, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19074,int_stack+43797,int_stack+7158, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+41319,int_stack+19074,int_stack+18714, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18714,int_stack+20208,int_stack+20178, 0.0, zero_stack, 1.0, int_stack+33798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18804,int_stack+20253,int_stack+20208, 0.0, zero_stack, 1.0, int_stack+33828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+18804,int_stack+18714, 0.0, zero_stack, 1.0, int_stack+16824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18939,int_stack+20316,int_stack+20253, 0.0, zero_stack, 1.0, int_stack+33873, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19128,int_stack+18939,int_stack+18804, 0.0, zero_stack, 1.0, int_stack+49491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18714,int_stack+19128,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+49626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+20460,int_stack+20400, 0.0, zero_stack, 1.0, int_stack+21066, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19014,int_stack+20550,int_stack+20460, 0.0, zero_stack, 1.0, int_stack+21186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19284,int_stack+19014,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+49806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+20676,int_stack+20550, 0.0, zero_stack, 1.0, int_stack+21366, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7158,int_stack+43797,int_stack+19014, 0.0, zero_stack, 1.0, int_stack+49986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+19644,int_stack+7158,int_stack+19284, 0.0, zero_stack, 1.0, int_stack+36462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7158,int_stack+20874,int_stack+20844, 1.0, int_stack+33798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7248,int_stack+20919,int_stack+20874, 1.0, int_stack+33828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+7248,int_stack+7158, 1.0, int_stack+16824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7383,int_stack+20982,int_stack+20919, 1.0, int_stack+33873, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19014,int_stack+7383,int_stack+7248, 1.0, int_stack+49491, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7158,int_stack+19014,int_stack+56829, 1.0, int_stack+49626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+21276,int_stack+21126, 1.0, int_stack+21066, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19014,int_stack+21492,int_stack+21276, 1.0, int_stack+21186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19284,int_stack+19014,int_stack+56829, 1.0, int_stack+49806, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+21618,int_stack+21492, 1.0, int_stack+21366, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20244,int_stack+43797,int_stack+19014, 1.0, int_stack+49986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+20784,int_stack+20244,int_stack+19284, 1.0, int_stack+36462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+36462,int_stack+21816,int_stack+21786,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+36552,int_stack+21861,int_stack+21816,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+36552,int_stack+36462,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+36687,int_stack+21924,int_stack+21861,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+36876,int_stack+36687,int_stack+36552,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+36462,int_stack+36876,int_stack+56829,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+22068,int_stack+22008,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+36762,int_stack+22158,int_stack+22068,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+20244,int_stack+36762,int_stack+56829,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+22284,int_stack+22158,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+19014,int_stack+43797,int_stack+36762,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+21384,int_stack+19014,int_stack+20244,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20244,int_stack+22482,int_stack+22452, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34686,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20334,int_stack+22527,int_stack+22482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34716,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+20334,int_stack+20244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51009,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20469,int_stack+22590,int_stack+22527, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34761,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19014,int_stack+20469,int_stack+20334, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51099,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+20244,int_stack+19014,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51234,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+22734,int_stack+22674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26004,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19014,int_stack+22824,int_stack+22734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26124,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19284,int_stack+19014,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51414,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+22950,int_stack+22824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26304,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21984,int_stack+43797,int_stack+19014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51594,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+49431,int_stack+21984,int_stack+19284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52173,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21984,int_stack+23148,int_stack+23118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34686, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22074,int_stack+23193,int_stack+23148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34716, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+22074,int_stack+21984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51009, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22209,int_stack+23256,int_stack+23193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34761, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22398,int_stack+22209,int_stack+22074, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51099, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+21984,int_stack+22398,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51234, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+23400,int_stack+23340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26004, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22284,int_stack+23490,int_stack+23400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26124, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22554,int_stack+22284,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51414, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+23616,int_stack+23490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26304, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22914,int_stack+43797,int_stack+22284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51594, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+19014,int_stack+22914,int_stack+22554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52173, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22284,int_stack+23814,int_stack+23784, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34686, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22374,int_stack+23859,int_stack+23814, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34716, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+22374,int_stack+22284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51009, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22509,int_stack+23922,int_stack+23859, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34761, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22698,int_stack+22509,int_stack+22374, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51099, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+22284,int_stack+22698,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51234, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+24066,int_stack+24006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26004, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22584,int_stack+24156,int_stack+24066, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26124, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22854,int_stack+22584,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51414, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+24282,int_stack+24156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26304, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23214,int_stack+43797,int_stack+22584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51594, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+23754,int_stack+23214,int_stack+22854, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52173, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22584,int_stack+24480,int_stack+24450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34686, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22674,int_stack+24525,int_stack+24480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+22674,int_stack+22584, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51009, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22809,int_stack+24588,int_stack+24525, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22998,int_stack+22809,int_stack+22674, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51099, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+22584,int_stack+22998,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+24732,int_stack+24672, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22884,int_stack+24822,int_stack+24732, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23154,int_stack+22884,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+24948,int_stack+24822, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24354,int_stack+43797,int_stack+22884, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13275,int_stack+24354,int_stack+23154, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24354,int_stack+25146,int_stack+25116, 0.0, zero_stack, 1.0, int_stack+34686, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24444,int_stack+25191,int_stack+25146, 0.0, zero_stack, 1.0, int_stack+34716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+24444,int_stack+24354, 0.0, zero_stack, 1.0, int_stack+51009, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24579,int_stack+25254,int_stack+25191, 0.0, zero_stack, 1.0, int_stack+34761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24768,int_stack+24579,int_stack+24444, 0.0, zero_stack, 1.0, int_stack+51099, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+25038,int_stack+24768,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+51234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+25398,int_stack+25338, 0.0, zero_stack, 1.0, int_stack+26004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24354,int_stack+25488,int_stack+25398, 0.0, zero_stack, 1.0, int_stack+26124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24624,int_stack+24354,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+51414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+25614,int_stack+25488, 0.0, zero_stack, 1.0, int_stack+26304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22884,int_stack+43797,int_stack+24354, 0.0, zero_stack, 1.0, int_stack+51594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13875,int_stack+22884,int_stack+24624, 0.0, zero_stack, 1.0, int_stack+52173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22884,int_stack+25812,int_stack+25782, 1.0, int_stack+34686, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22974,int_stack+25857,int_stack+25812, 1.0, int_stack+34716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+22974,int_stack+22884, 1.0, int_stack+51009, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23109,int_stack+25920,int_stack+25857, 1.0, int_stack+34761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23298,int_stack+23109,int_stack+22974, 1.0, int_stack+51099, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+22884,int_stack+23298,int_stack+56829, 1.0, int_stack+51234, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+26214,int_stack+26064, 1.0, int_stack+26004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+51009,int_stack+26430,int_stack+26214, 1.0, int_stack+26124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23184,int_stack+51009,int_stack+56829, 1.0, int_stack+51414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+26556,int_stack+26430, 1.0, int_stack+26304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24354,int_stack+43797,int_stack+51009, 1.0, int_stack+51594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+51009,int_stack+24354,int_stack+23184, 1.0, int_stack+52173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52173,int_stack+26754,int_stack+26724,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+52263,int_stack+26799,int_stack+26754,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+52263,int_stack+52173,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+23184,int_stack+26862,int_stack+26799,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23373,int_stack+23184,int_stack+52263,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+52173,int_stack+23373,int_stack+56829,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+27006,int_stack+26946,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+27096,int_stack+27006,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+24354,int_stack+23184,int_stack+56829,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+27222,int_stack+27096,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+25338,int_stack+43797,int_stack+23184,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+25878,int_stack+25338,int_stack+24354,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24354,int_stack+27420,int_stack+27390,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24444,int_stack+27465,int_stack+27420,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+24444,int_stack+24354,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+24579,int_stack+27528,int_stack+27465,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+24768,int_stack+24579,int_stack+24444,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+24354,int_stack+24768,int_stack+56829,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+27672,int_stack+27612,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+24654,int_stack+27762,int_stack+27672,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+25338,int_stack+24654,int_stack+56829,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+27888,int_stack+27762,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+43797,int_stack+24654,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+26478,int_stack+23184,int_stack+25338,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25338,int_stack+28308,int_stack+28278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35574,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25428,int_stack+28353,int_stack+28308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35604,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+25428,int_stack+25338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52533,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+25563,int_stack+28416,int_stack+28353, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35649,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+25563,int_stack+25428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52623,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+23454,int_stack+23184,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52758,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+28560,int_stack+28500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33078,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+28650,int_stack+28560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33198,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25338,int_stack+23184,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52938,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+28776,int_stack+28650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33378,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27078,int_stack+43797,int_stack+23184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53118,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+27618,int_stack+27078,int_stack+25338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54351,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25338,int_stack+29196,int_stack+29166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35574, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25428,int_stack+29241,int_stack+29196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35604, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+25428,int_stack+25338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52533, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+25563,int_stack+29304,int_stack+29241, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35649, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+25563,int_stack+25428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52623, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+25338,int_stack+23184,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52758, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+29448,int_stack+29388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33078, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+29538,int_stack+29448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33198, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+27078,int_stack+23184,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52938, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+29664,int_stack+29538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33378, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+28218,int_stack+43797,int_stack+23184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53118, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+28758,int_stack+28218,int_stack+27078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54351, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27078,int_stack+30084,int_stack+30054, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35574, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+27168,int_stack+30129,int_stack+30084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35604, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+27168,int_stack+27078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52533, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+27303,int_stack+30192,int_stack+30129, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35649, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+27303,int_stack+27168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52623, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+27078,int_stack+23184,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52758, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+30336,int_stack+30276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33078, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+30426,int_stack+30336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33198, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+28218,int_stack+23184,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52938, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+30552,int_stack+30426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33378, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29358,int_stack+43797,int_stack+23184, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53118, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+29898,int_stack+29358,int_stack+28218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54351, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28218,int_stack+30972,int_stack+30942, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+28308,int_stack+31017,int_stack+30972, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+28308,int_stack+28218, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52533, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+28443,int_stack+31080,int_stack+31017, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35649, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+28443,int_stack+28308, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+28218,int_stack+23184,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+31224,int_stack+31164, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+31314,int_stack+31224, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29358,int_stack+23184,int_stack+56829, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+31440,int_stack+31314, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+30498,int_stack+43797,int_stack+23184, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+31038,int_stack+30498,int_stack+29358, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29358,int_stack+31860,int_stack+31830, 0.0, zero_stack, 1.0, int_stack+35574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29448,int_stack+31905,int_stack+31860, 0.0, zero_stack, 1.0, int_stack+35604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+29448,int_stack+29358, 0.0, zero_stack, 1.0, int_stack+52533, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29583,int_stack+31968,int_stack+31905, 0.0, zero_stack, 1.0, int_stack+35649, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+29583,int_stack+29448, 0.0, zero_stack, 1.0, int_stack+52623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+29358,int_stack+23184,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+52758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+32112,int_stack+32052, 0.0, zero_stack, 1.0, int_stack+33078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+32202,int_stack+32112, 0.0, zero_stack, 1.0, int_stack+33198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30498,int_stack+23184,int_stack+56829, 0.0, zero_stack, 1.0, int_stack+52938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+32328,int_stack+32202, 0.0, zero_stack, 1.0, int_stack+33378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+31638,int_stack+43797,int_stack+23184, 0.0, zero_stack, 1.0, int_stack+53118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14475,int_stack+31638,int_stack+30498, 0.0, zero_stack, 1.0, int_stack+54351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+30498,int_stack+32886,int_stack+32856, 1.0, int_stack+35574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+30588,int_stack+32931,int_stack+32886, 1.0, int_stack+35604, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+30588,int_stack+30498, 1.0, int_stack+52533, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+30723,int_stack+32994,int_stack+32931, 1.0, int_stack+35649, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+30723,int_stack+30588, 1.0, int_stack+52623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+30498,int_stack+23184,int_stack+56829, 1.0, int_stack+52758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+33288,int_stack+33138, 1.0, int_stack+33078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+33504,int_stack+33288, 1.0, int_stack+33198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31638,int_stack+23184,int_stack+56829, 1.0, int_stack+52938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+33630,int_stack+33504, 1.0, int_stack+33378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+52473,int_stack+43797,int_stack+23184, 1.0, int_stack+53118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+32796,int_stack+52473,int_stack+31638, 1.0, int_stack+54351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+54351,int_stack+34050,int_stack+34020,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+54441,int_stack+34095,int_stack+34050,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+54441,int_stack+54351,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+31638,int_stack+34158,int_stack+34095,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+31638,int_stack+54441,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+31638,int_stack+23184,int_stack+56829,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+34302,int_stack+34242,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+34392,int_stack+34302,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+54351,int_stack+23184,int_stack+56829,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+34518,int_stack+34392,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+31938,int_stack+43797,int_stack+23184,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+52473,int_stack+31938,int_stack+54351,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+54351,int_stack+34938,int_stack+34908,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+54441,int_stack+34983,int_stack+34938,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+54441,int_stack+54351,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+31938,int_stack+35046,int_stack+34983,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+31938,int_stack+54441,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+31938,int_stack+23184,int_stack+56829,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+35190,int_stack+35130,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+35280,int_stack+35190,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+54351,int_stack+23184,int_stack+56829,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+35406,int_stack+35280,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+33396,int_stack+43797,int_stack+23184,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+33936,int_stack+33396,int_stack+54351,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+54351,int_stack+35826,int_stack+35796,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+54441,int_stack+35871,int_stack+35826,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56829,int_stack+54441,int_stack+54351,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+33396,int_stack+35934,int_stack+35871,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+23184,int_stack+33396,int_stack+54441,3);
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+33396,int_stack+23184,int_stack+56829,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+56829,int_stack+36078,int_stack+36018,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+23184,int_stack+36168,int_stack+36078,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+54351,int_stack+23184,int_stack+56829,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+43797,int_stack+36294,int_stack+36168,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+34536,int_stack+43797,int_stack+23184,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+35076,int_stack+34536,int_stack+54351,6);
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+15075,int_stack+40119,int_stack+37731,100);
     Libderiv->ABCD[11] = int_stack + 15075;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+15975,int_stack+42297,int_stack+40719,100);
     Libderiv->ABCD[10] = int_stack + 15975;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+41919,int_stack+44175,int_stack+0,100);
     Libderiv->ABCD[9] = int_stack + 41919;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+43797,int_stack+46353,int_stack+44775,100);
     Libderiv->ABCD[8] = int_stack + 43797;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+38931,int_stack+48231,int_stack+300,100);
     Libderiv->ABCD[7] = int_stack + 38931;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+47853,int_stack+50409,int_stack+48831,100);
     Libderiv->ABCD[6] = int_stack + 47853;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+50031,int_stack+600,int_stack+37272, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 50031;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+61728,int_stack+53451,int_stack+51873, 0.0, zero_stack, 1.0, int_stack+38031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 61728;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+53073,int_stack+55629,int_stack+54051, 1.0, int_stack+38031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 53073;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+62628,int_stack+1740,int_stack+56229,100);
     Libderiv->ABCD[155] = int_stack + 62628;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+55611,int_stack+2340,int_stack+56529,100);
     Libderiv->ABCD[143] = int_stack + 55611;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1500,int_stack+57369,int_stack+2940,100);
     Libderiv->ABCD[142] = int_stack + 1500;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2400,int_stack+55011,int_stack+54711,100);
     Libderiv->ABCD[131] = int_stack + 2400;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+54351,int_stack+3999,int_stack+3699,100);
     Libderiv->ABCD[130] = int_stack + 54351;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3300,int_stack+4899,int_stack+4599,100);
     Libderiv->ABCD[129] = int_stack + 3300;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4200,int_stack+5499,int_stack+58428,100);
     Libderiv->ABCD[119] = int_stack + 4200;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+5100,int_stack+58728,int_stack+6099,100);
     Libderiv->ABCD[118] = int_stack + 5100;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+58269,int_stack+59328,int_stack+6858,100);
     Libderiv->ABCD[117] = int_stack + 58269;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+63528,int_stack+7698,int_stack+57969,100);
     Libderiv->ABCD[116] = int_stack + 63528;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+57309,int_stack+59928,int_stack+6399,100);
     Libderiv->ABCD[107] = int_stack + 57309;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+6000,int_stack+9057,int_stack+8757,100);
     Libderiv->ABCD[106] = int_stack + 6000;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+8598,int_stack+60528,int_stack+1200,100);
     Libderiv->ABCD[105] = int_stack + 8598;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+600,int_stack+10416,int_stack+10116,100);
     Libderiv->ABCD[104] = int_stack + 600;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+9498,int_stack+61128,int_stack+8298,100);
     Libderiv->ABCD[103] = int_stack + 9498;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+7458,int_stack+11775,int_stack+11475,100);
     Libderiv->ABCD[95] = int_stack + 7458;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+10398,int_stack+38331,int_stack+12375,100);
     Libderiv->ABCD[94] = int_stack + 10398;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+38031,int_stack+12675,int_stack+41019,100);
     Libderiv->ABCD[93] = int_stack + 38031;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+11298,int_stack+43197,int_stack+42897,100);
     Libderiv->ABCD[92] = int_stack + 11298;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+42819,int_stack+45735,int_stack+45075,100);
     Libderiv->ABCD[91] = int_stack + 42819;
 /*--- compute (pp|ff) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+45675,int_stack+16914,int_stack+46953,100);
     Libderiv->ABCD[90] = int_stack + 45675;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+12198,int_stack+17514,int_stack+49131, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37731, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[47] = int_stack + 12198;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+16875,int_stack+47253,int_stack+57009, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[46] = int_stack + 16875;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+46575,int_stack+18114,int_stack+45375, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[45] = int_stack + 46575;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+17775,int_stack+41319,int_stack+32496, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[44] = int_stack + 17775;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+41019,int_stack+19644,int_stack+18714, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[43] = int_stack + 41019;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+59169,int_stack+20784,int_stack+7158, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[42] = int_stack + 59169;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+60069,int_stack+21384,int_stack+36462, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+37272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[38] = int_stack + 60069;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+20544,int_stack+49431,int_stack+20244, 0.0, zero_stack, 1.0, int_stack+37731, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[35] = int_stack + 20544;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+49131,int_stack+19014,int_stack+21984, 0.0, zero_stack, 1.0, int_stack+40719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[34] = int_stack + 49131;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+35676,int_stack+23754,int_stack+22284, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[33] = int_stack + 35676;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+21444,int_stack+13275,int_stack+22584, 0.0, zero_stack, 1.0, int_stack+44775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[32] = int_stack + 21444;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+18675,int_stack+13875,int_stack+25038, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[31] = int_stack + 18675;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+19575,int_stack+51009,int_stack+22884, 0.0, zero_stack, 1.0, int_stack+48831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[30] = int_stack + 19575;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+22344,int_stack+25878,int_stack+52173, 0.0, zero_stack, 1.0, int_stack+37272, 1.0, int_stack+51873, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[26] = int_stack + 22344;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+50931,int_stack+26478,int_stack+24354, 0.0, zero_stack, 2.0, int_stack+51873, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[25] = int_stack + 50931;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+23754,int_stack+27618,int_stack+23454, 1.0, int_stack+37731, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[23] = int_stack + 23754;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+25638,int_stack+28758,int_stack+25338, 1.0, int_stack+40719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[22] = int_stack + 25638;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+24654,int_stack+29898,int_stack+27078, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[21] = int_stack + 24654;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+26538,int_stack+31038,int_stack+28218, 1.0, int_stack+44775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[20] = int_stack + 26538;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+44697,int_stack+14475,int_stack+29358, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[19] = int_stack + 44697;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+27438,int_stack+32796,int_stack+30498, 1.0, int_stack+48831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[18] = int_stack + 27438;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+32238,int_stack+52473,int_stack+31638, 1.0, int_stack+37272, 0.0, zero_stack, 1.0, int_stack+54051, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[14] = int_stack + 32238;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+52173,int_stack+33936,int_stack+31938, 1.0, int_stack+51873, 1.0, int_stack+54051, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[13] = int_stack + 52173;
 /*--- compute (pp|ff) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+33696,int_stack+35076,int_stack+33396, 2.0, int_stack+54051, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[12] = int_stack + 33696;

}
