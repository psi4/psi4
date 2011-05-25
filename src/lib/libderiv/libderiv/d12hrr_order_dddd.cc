#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_dddd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|dd) integrals */

void d12hrr_order_dddd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][4][10] = int_stack + 225;
 Libderiv->deriv_classes[4][4][9] = int_stack + 450;
 Libderiv->deriv_classes[4][4][8] = int_stack + 675;
 Libderiv->deriv_classes[4][4][7] = int_stack + 900;
 Libderiv->dvrr_classes[4][3] = int_stack + 1125;
 Libderiv->deriv_classes[4][4][6] = int_stack + 1275;
 Libderiv->deriv_classes[4][4][2] = int_stack + 1500;
 Libderiv->deriv_classes[4][4][1] = int_stack + 1725;
 Libderiv->dvrr_classes[3][4] = int_stack + 1950;
 Libderiv->deriv_classes[4][4][0] = int_stack + 2100;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 2325;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 2361;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 2421;
 Libderiv->deriv2_classes[3][2][143] = int_stack + 2511;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 2571;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 2671;
 Libderiv->deriv2_classes[4][2][143] = int_stack + 2821;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 2911;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 3061;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 3286;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 3322;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 3382;
 Libderiv->deriv2_classes[3][2][131] = int_stack + 3472;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 3532;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 3632;
 Libderiv->deriv2_classes[4][2][131] = int_stack + 3782;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 3872;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 4022;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 4247;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 4283;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 4343;
 Libderiv->deriv2_classes[3][2][130] = int_stack + 4433;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 4493;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 4593;
 Libderiv->deriv2_classes[4][2][130] = int_stack + 4743;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 4833;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 4983;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 5208;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 5244;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 5304;
 Libderiv->deriv2_classes[3][2][119] = int_stack + 5394;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 5454;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 5554;
 Libderiv->deriv2_classes[4][2][119] = int_stack + 5704;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 5794;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 5944;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 6169;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 6205;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 6265;
 Libderiv->deriv2_classes[3][2][118] = int_stack + 6355;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 6415;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 6515;
 Libderiv->deriv2_classes[4][2][118] = int_stack + 6665;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 6755;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 6905;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 7130;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 7166;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 7226;
 Libderiv->deriv2_classes[3][2][117] = int_stack + 7316;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 7376;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 7476;
 Libderiv->deriv2_classes[4][2][117] = int_stack + 7626;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 7716;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 7866;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 8091;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 8127;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 8187;
 Libderiv->deriv2_classes[3][2][107] = int_stack + 8277;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 8337;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 8437;
 Libderiv->deriv2_classes[4][2][107] = int_stack + 8587;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 8677;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 8827;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 9052;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 9088;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 9148;
 Libderiv->deriv2_classes[3][2][106] = int_stack + 9238;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 9298;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 9398;
 Libderiv->deriv2_classes[4][2][106] = int_stack + 9548;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 9638;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 9788;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 10013;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 10049;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 10109;
 Libderiv->deriv2_classes[3][2][105] = int_stack + 10199;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 10259;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 10359;
 Libderiv->deriv2_classes[4][2][105] = int_stack + 10509;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 10599;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 10749;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 10974;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 11010;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 11070;
 Libderiv->deriv2_classes[3][2][104] = int_stack + 11160;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 11220;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 11320;
 Libderiv->deriv2_classes[4][2][104] = int_stack + 11470;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 11560;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 11710;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 11935;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 11971;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 12031;
 Libderiv->deriv2_classes[3][2][95] = int_stack + 12121;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 12181;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 12281;
 Libderiv->deriv2_classes[4][2][95] = int_stack + 12431;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 12521;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 12671;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 12896;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 12932;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 12992;
 Libderiv->deriv2_classes[3][2][94] = int_stack + 13082;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 13142;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 13242;
 Libderiv->deriv2_classes[4][2][94] = int_stack + 13392;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 13482;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 13632;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 13857;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 13893;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 13953;
 Libderiv->deriv2_classes[3][2][93] = int_stack + 14043;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 14103;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 14203;
 Libderiv->deriv2_classes[4][2][93] = int_stack + 14353;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 14443;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 14593;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 14818;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 14854;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 14914;
 Libderiv->deriv2_classes[3][2][92] = int_stack + 15004;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 15064;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 15164;
 Libderiv->deriv2_classes[4][2][92] = int_stack + 15314;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 15404;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 15554;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 15779;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 15815;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 15875;
 Libderiv->deriv2_classes[3][2][91] = int_stack + 15965;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 16025;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 16125;
 Libderiv->deriv2_classes[4][2][91] = int_stack + 16275;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 16365;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 16515;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 16740;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 16776;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 16836;
 Libderiv->deriv2_classes[3][2][83] = int_stack + 16926;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 16986;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 17086;
 Libderiv->deriv_classes[4][2][11] = int_stack + 17236;
 Libderiv->deriv2_classes[4][2][83] = int_stack + 17326;
 Libderiv->deriv_classes[4][3][11] = int_stack + 17416;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 17566;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 17716;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 17941;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 17977;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 18037;
 Libderiv->deriv2_classes[3][2][82] = int_stack + 18127;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 18187;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 18287;
 Libderiv->deriv_classes[4][2][10] = int_stack + 18437;
 Libderiv->deriv2_classes[4][2][82] = int_stack + 18527;
 Libderiv->deriv_classes[4][3][10] = int_stack + 18617;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 18767;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 18917;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 19142;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 19178;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 19238;
 Libderiv->deriv2_classes[3][2][81] = int_stack + 19328;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 19388;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 19488;
 Libderiv->deriv_classes[4][2][9] = int_stack + 19638;
 Libderiv->deriv2_classes[4][2][81] = int_stack + 19728;
 Libderiv->deriv_classes[4][3][9] = int_stack + 19818;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 19968;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 20118;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 20343;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 20379;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 20439;
 Libderiv->deriv2_classes[3][2][80] = int_stack + 20529;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 20589;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 20689;
 Libderiv->deriv_classes[4][2][8] = int_stack + 20839;
 Libderiv->deriv2_classes[4][2][80] = int_stack + 20929;
 Libderiv->deriv_classes[4][3][8] = int_stack + 21019;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 21169;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 21319;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 21544;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 21580;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 21640;
 Libderiv->deriv2_classes[3][2][79] = int_stack + 21730;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 21790;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 21890;
 Libderiv->deriv_classes[4][2][7] = int_stack + 22040;
 Libderiv->deriv2_classes[4][2][79] = int_stack + 22130;
 Libderiv->deriv_classes[4][3][7] = int_stack + 22220;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 22370;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 22520;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 22745;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 22781;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 22841;
 Libderiv->deriv2_classes[3][2][78] = int_stack + 22931;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 22991;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 23091;
 Libderiv->dvrr_classes[4][2] = int_stack + 23241;
 Libderiv->deriv_classes[4][2][6] = int_stack + 23331;
 Libderiv->deriv2_classes[4][2][78] = int_stack + 23421;
 Libderiv->deriv_classes[4][3][6] = int_stack + 23511;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 23661;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 23811;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 24036;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 24072;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 24132;
 Libderiv->deriv2_classes[3][2][35] = int_stack + 24222;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 24282;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 24382;
 Libderiv->deriv2_classes[4][2][35] = int_stack + 24532;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 24622;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 24772;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 24997;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 25033;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 25093;
 Libderiv->deriv2_classes[3][2][34] = int_stack + 25183;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 25243;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 25343;
 Libderiv->deriv2_classes[4][2][34] = int_stack + 25493;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 25583;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 25733;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 25958;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 25994;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 26054;
 Libderiv->deriv2_classes[3][2][33] = int_stack + 26144;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 26204;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 26304;
 Libderiv->deriv2_classes[4][2][33] = int_stack + 26454;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 26544;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 26694;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 26919;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 26955;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 27015;
 Libderiv->deriv2_classes[3][2][32] = int_stack + 27105;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 27165;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 27265;
 Libderiv->deriv2_classes[4][2][32] = int_stack + 27415;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 27505;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 27655;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 27880;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 27916;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 27976;
 Libderiv->deriv2_classes[3][2][31] = int_stack + 28066;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 28126;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 28226;
 Libderiv->deriv2_classes[4][2][31] = int_stack + 28376;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 28466;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 28616;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 28841;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 28877;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 28937;
 Libderiv->deriv2_classes[3][2][30] = int_stack + 29027;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 29087;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 29187;
 Libderiv->deriv_classes[4][2][2] = int_stack + 29337;
 Libderiv->deriv2_classes[4][2][30] = int_stack + 29427;
 Libderiv->deriv_classes[4][3][2] = int_stack + 29517;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 29667;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 29817;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 30042;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 30078;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 30138;
 Libderiv->deriv2_classes[3][2][26] = int_stack + 30228;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 30288;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 30388;
 Libderiv->deriv2_classes[4][2][26] = int_stack + 30538;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 30628;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 30778;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 31003;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 31039;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 31099;
 Libderiv->deriv2_classes[3][2][23] = int_stack + 31189;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 31249;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 31349;
 Libderiv->deriv2_classes[4][2][23] = int_stack + 31499;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 31589;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 31739;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 31964;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 32000;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 32060;
 Libderiv->deriv2_classes[3][2][22] = int_stack + 32150;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 32210;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 32310;
 Libderiv->deriv2_classes[4][2][22] = int_stack + 32460;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 32550;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 32700;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 32925;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 32961;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 33021;
 Libderiv->deriv2_classes[3][2][21] = int_stack + 33111;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 33171;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 33271;
 Libderiv->deriv2_classes[4][2][21] = int_stack + 33421;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 33511;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 33661;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 33886;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 33922;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 33982;
 Libderiv->deriv2_classes[3][2][20] = int_stack + 34072;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 34132;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 34232;
 Libderiv->deriv2_classes[4][2][20] = int_stack + 34382;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 34472;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 34622;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 34847;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 34883;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 34943;
 Libderiv->deriv2_classes[3][2][19] = int_stack + 35033;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 35093;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 35193;
 Libderiv->deriv2_classes[4][2][19] = int_stack + 35343;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 35433;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 35583;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 35808;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 35844;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 35904;
 Libderiv->deriv2_classes[3][2][18] = int_stack + 35994;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 36054;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 36154;
 Libderiv->deriv_classes[4][2][1] = int_stack + 36304;
 Libderiv->deriv2_classes[4][2][18] = int_stack + 36394;
 Libderiv->deriv_classes[4][3][1] = int_stack + 36484;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 36634;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 36784;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 37009;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 37045;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 37105;
 Libderiv->deriv2_classes[3][2][14] = int_stack + 37195;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 37255;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 37355;
 Libderiv->deriv2_classes[4][2][14] = int_stack + 37505;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 37595;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 37745;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 37970;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 38006;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 38066;
 Libderiv->deriv2_classes[3][2][13] = int_stack + 38156;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 38216;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 38316;
 Libderiv->deriv2_classes[4][2][13] = int_stack + 38466;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 38556;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 38706;
 Libderiv->deriv_classes[2][2][11] = int_stack + 38931;
 Libderiv->deriv_classes[2][3][11] = int_stack + 38967;
 Libderiv->deriv_classes[2][4][11] = int_stack + 39027;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 39117;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 39153;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 39213;
 Libderiv->deriv_classes[3][2][11] = int_stack + 39303;
 Libderiv->deriv_classes[3][3][11] = int_stack + 39363;
 Libderiv->deriv_classes[3][4][11] = int_stack + 39463;
 Libderiv->deriv2_classes[3][2][11] = int_stack + 39613;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 39673;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 39773;
 Libderiv->deriv2_classes[4][2][11] = int_stack + 39923;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 40013;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 40163;
 Libderiv->deriv_classes[2][2][10] = int_stack + 40388;
 Libderiv->deriv_classes[2][3][10] = int_stack + 40424;
 Libderiv->deriv_classes[2][4][10] = int_stack + 40484;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 40574;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 40610;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 40670;
 Libderiv->deriv_classes[3][2][10] = int_stack + 40760;
 Libderiv->deriv_classes[3][3][10] = int_stack + 40820;
 Libderiv->deriv_classes[3][4][10] = int_stack + 40920;
 Libderiv->deriv2_classes[3][2][10] = int_stack + 41070;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 41130;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 41230;
 Libderiv->deriv2_classes[4][2][10] = int_stack + 41380;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 41470;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 41620;
 Libderiv->deriv_classes[2][2][9] = int_stack + 41845;
 Libderiv->deriv_classes[2][3][9] = int_stack + 41881;
 Libderiv->deriv_classes[2][4][9] = int_stack + 41941;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 42031;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 42067;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 42127;
 Libderiv->deriv_classes[3][2][9] = int_stack + 42217;
 Libderiv->deriv_classes[3][3][9] = int_stack + 42277;
 Libderiv->deriv_classes[3][4][9] = int_stack + 42377;
 Libderiv->deriv2_classes[3][2][9] = int_stack + 42527;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 42587;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 42687;
 Libderiv->deriv2_classes[4][2][9] = int_stack + 42837;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 42927;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 43077;
 Libderiv->deriv_classes[2][2][8] = int_stack + 43302;
 Libderiv->deriv_classes[2][3][8] = int_stack + 43338;
 Libderiv->deriv_classes[2][4][8] = int_stack + 43398;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 43488;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 43524;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 43584;
 Libderiv->deriv_classes[3][2][8] = int_stack + 43674;
 Libderiv->deriv_classes[3][3][8] = int_stack + 43734;
 Libderiv->deriv_classes[3][4][8] = int_stack + 43834;
 Libderiv->deriv2_classes[3][2][8] = int_stack + 43984;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 44044;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 44144;
 Libderiv->deriv2_classes[4][2][8] = int_stack + 44294;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 44384;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 44534;
 Libderiv->deriv_classes[2][2][7] = int_stack + 44759;
 Libderiv->deriv_classes[2][3][7] = int_stack + 44795;
 Libderiv->deriv_classes[2][4][7] = int_stack + 44855;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 44945;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 44981;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 45041;
 Libderiv->deriv_classes[3][2][7] = int_stack + 45131;
 Libderiv->deriv_classes[3][3][7] = int_stack + 45191;
 Libderiv->deriv_classes[3][4][7] = int_stack + 45291;
 Libderiv->deriv2_classes[3][2][7] = int_stack + 45441;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 45501;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 45601;
 Libderiv->deriv2_classes[4][2][7] = int_stack + 45751;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 45841;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 45991;
 Libderiv->deriv_classes[2][2][6] = int_stack + 46216;
 Libderiv->deriv_classes[2][3][6] = int_stack + 46252;
 Libderiv->deriv_classes[2][4][6] = int_stack + 46312;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 46402;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 46438;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 46498;
 Libderiv->dvrr_classes[3][2] = int_stack + 46588;
 Libderiv->deriv_classes[3][2][6] = int_stack + 46648;
 Libderiv->dvrr_classes[3][3] = int_stack + 46708;
 Libderiv->deriv_classes[3][3][6] = int_stack + 46808;
 Libderiv->deriv_classes[3][4][6] = int_stack + 46908;
 Libderiv->deriv2_classes[3][2][6] = int_stack + 47058;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 47118;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 47218;
 Libderiv->deriv_classes[4][2][0] = int_stack + 47368;
 Libderiv->deriv2_classes[4][2][6] = int_stack + 47458;
 Libderiv->deriv_classes[4][3][0] = int_stack + 47548;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 47698;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 47848;
 Libderiv->deriv_classes[2][2][2] = int_stack + 48073;
 Libderiv->deriv_classes[2][3][2] = int_stack + 48109;
 Libderiv->deriv_classes[2][4][2] = int_stack + 48169;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 48259;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 48295;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 48355;
 Libderiv->deriv_classes[3][2][2] = int_stack + 48445;
 Libderiv->deriv_classes[3][3][2] = int_stack + 48505;
 Libderiv->deriv_classes[3][4][2] = int_stack + 48605;
 Libderiv->deriv2_classes[3][2][2] = int_stack + 48755;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 48815;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 48915;
 Libderiv->deriv2_classes[4][2][2] = int_stack + 49065;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 49155;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 49305;
 Libderiv->deriv_classes[2][2][1] = int_stack + 49530;
 Libderiv->deriv_classes[2][3][1] = int_stack + 49566;
 Libderiv->deriv_classes[2][4][1] = int_stack + 49626;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 49716;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 49752;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 49812;
 Libderiv->deriv_classes[3][2][1] = int_stack + 49902;
 Libderiv->deriv_classes[3][3][1] = int_stack + 49962;
 Libderiv->deriv_classes[3][4][1] = int_stack + 50062;
 Libderiv->deriv2_classes[3][2][1] = int_stack + 50212;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 50272;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 50372;
 Libderiv->deriv2_classes[4][2][1] = int_stack + 50522;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 50612;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 50762;
 Libderiv->dvrr_classes[2][2] = int_stack + 50987;
 Libderiv->dvrr_classes[2][3] = int_stack + 51023;
 Libderiv->dvrr_classes[2][4] = int_stack + 51083;
 Libderiv->deriv_classes[2][2][0] = int_stack + 51173;
 Libderiv->deriv_classes[2][3][0] = int_stack + 51209;
 Libderiv->deriv_classes[2][4][0] = int_stack + 51269;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 51359;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 51395;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 51455;
 Libderiv->deriv_classes[3][2][0] = int_stack + 51545;
 Libderiv->deriv_classes[3][3][0] = int_stack + 51605;
 Libderiv->deriv_classes[3][4][0] = int_stack + 51705;
 Libderiv->deriv2_classes[3][2][0] = int_stack + 51855;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 51915;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 52015;
 Libderiv->deriv2_classes[4][2][0] = int_stack + 52165;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 52255;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 52405;
 memset(int_stack,0,421040);

 Libderiv->dvrr_stack = int_stack + 106750;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_dddd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+52630,int_stack+51023,int_stack+50987,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+38967,int_stack+38931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50987,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+52846,int_stack+39027,int_stack+38967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51023,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+53026,int_stack+52846,int_stack+52738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+52846,int_stack+46708,int_stack+46588,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+39363,int_stack+39303, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46588,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53422,int_stack+39463,int_stack+39363, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46708,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+53722,int_stack+53422,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52846,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+54082,int_stack+53722,int_stack+53026,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+53422,int_stack+1125,int_stack+23241,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54730,int_stack+17416,int_stack+17236, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23241,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+55000,int_stack+0,int_stack+17416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+55450,int_stack+55000,int_stack+54730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53422,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+55990,int_stack+55450,int_stack+53722,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+55000,int_stack+40424,int_stack+40388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50987, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+55108,int_stack+40484,int_stack+40424, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51023, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+55288,int_stack+55108,int_stack+55000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+55108,int_stack+40820,int_stack+40760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46588, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+55504,int_stack+40920,int_stack+40820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46708, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+57070,int_stack+55504,int_stack+55108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52846, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+57430,int_stack+57070,int_stack+55288,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+55504,int_stack+18617,int_stack+18437, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23241, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+58078,int_stack+225,int_stack+18617, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+58528,int_stack+58078,int_stack+55504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53422, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+59068,int_stack+58528,int_stack+57070,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+58078,int_stack+41881,int_stack+41845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50987, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+58186,int_stack+41941,int_stack+41881, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51023, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+55774,int_stack+58186,int_stack+58078, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+58186,int_stack+42277,int_stack+42217, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46588, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+58366,int_stack+42377,int_stack+42277, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46708, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+58666,int_stack+58366,int_stack+58186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52846, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+60148,int_stack+58666,int_stack+55774,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+58366,int_stack+19818,int_stack+19638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23241, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+450,int_stack+19818, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+60796,int_stack+0,int_stack+58366, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53422, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+61336,int_stack+60796,int_stack+58666,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+60796,int_stack+43338,int_stack+43302, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+60904,int_stack+43398,int_stack+43338, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51023, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+61084,int_stack+60904,int_stack+60796, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+60904,int_stack+43734,int_stack+43674, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+43834,int_stack+43734, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+300,int_stack+0,int_stack+60904, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+62416,int_stack+300,int_stack+61084,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+21019,int_stack+20839, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23241, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63064,int_stack+675,int_stack+21019, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+63514,int_stack+63064,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+64054,int_stack+63514,int_stack+300,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63064,int_stack+44795,int_stack+44759, 0.0, zero_stack, 1.0, int_stack+50987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63172,int_stack+44855,int_stack+44795, 0.0, zero_stack, 1.0, int_stack+51023, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+63352,int_stack+63172,int_stack+63064, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63172,int_stack+45191,int_stack+45131, 0.0, zero_stack, 1.0, int_stack+46588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63568,int_stack+45291,int_stack+45191, 0.0, zero_stack, 1.0, int_stack+46708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+65134,int_stack+63568,int_stack+63172, 0.0, zero_stack, 1.0, int_stack+52846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+65494,int_stack+65134,int_stack+63352,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+22220,int_stack+22040, 0.0, zero_stack, 1.0, int_stack+23241, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+66142,int_stack+900,int_stack+22220, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+66592,int_stack+66142,int_stack+63568, 0.0, zero_stack, 1.0, int_stack+53422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+67132,int_stack+66592,int_stack+65134,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+66142,int_stack+46252,int_stack+46216, 1.0, int_stack+50987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+66250,int_stack+46312,int_stack+46252, 1.0, int_stack+51023, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+63838,int_stack+66250,int_stack+66142, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+66250,int_stack+46808,int_stack+46648, 1.0, int_stack+46588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+66430,int_stack+46908,int_stack+46808, 1.0, int_stack+46708, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+66730,int_stack+66430,int_stack+66250, 1.0, int_stack+52846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+68212,int_stack+66730,int_stack+63838,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+66430,int_stack+23511,int_stack+23331, 1.0, int_stack+23241, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+660,int_stack+1275,int_stack+23511, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+68860,int_stack+660,int_stack+66430, 1.0, int_stack+53422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+69400,int_stack+68860,int_stack+66730,36);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+51083,int_stack+51023,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+69040,int_stack+68860,int_stack+52630,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53422,int_stack+1950,int_stack+46708,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+53422,int_stack+52846,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+70480,int_stack+660,int_stack+69040,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+52630,int_stack+48109,int_stack+48073,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52846,int_stack+48169,int_stack+48109,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+53422,int_stack+52846,int_stack+52630,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+52846,int_stack+48505,int_stack+48445,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1020,int_stack+48605,int_stack+48505,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+71128,int_stack+1020,int_stack+52846,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+71488,int_stack+71128,int_stack+53422, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1020,int_stack+29517,int_stack+29337,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+72136,int_stack+1500,int_stack+29517,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+72586,int_stack+72136,int_stack+1020,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+73126,int_stack+72586,int_stack+71128, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+49566,int_stack+49530,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+49626,int_stack+49566,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+72244,int_stack+68860,int_stack+72136,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+49962,int_stack+49902,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+50062,int_stack+49962,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+72760,int_stack+72460,int_stack+68860,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+74206,int_stack+72760,int_stack+72244, 0.0, zero_stack, 1.0, int_stack+69040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+72460,int_stack+36484,int_stack+36304,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+74854,int_stack+1725,int_stack+36484,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1290,int_stack+74854,int_stack+72460,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+74854,int_stack+1290,int_stack+72760, 0.0, zero_stack, 1.0, int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1290,int_stack+51209,int_stack+51173,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1398,int_stack+51269,int_stack+51209,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1578,int_stack+1398,int_stack+1290,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1398,int_stack+51605,int_stack+51545,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1794,int_stack+51705,int_stack+51605,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+75934,int_stack+1794,int_stack+1398,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+76294,int_stack+75934,int_stack+1578, 1.0, int_stack+69040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+69040,int_stack+47548,int_stack+47368,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+76942,int_stack+2100,int_stack+47548,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+77392,int_stack+76942,int_stack+69040,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+77932,int_stack+77392,int_stack+75934, 1.0, int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+2361,int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+38931,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+768,int_stack+2421,int_stack+2361, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+38967,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+76942,int_stack+768,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+52738,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+2571,int_stack+2511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+39303,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+77158,int_stack+2671,int_stack+2571, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+39363,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+77458,int_stack+77158,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+53242,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1794,int_stack+77458,int_stack+76942,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+76942,int_stack+2911,int_stack+2821, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17236,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2442,int_stack+3061,int_stack+2911, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17416,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+2442,int_stack+76942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+54730,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+79552,int_stack+79012,int_stack+77458,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+79012,int_stack+3322,int_stack+3286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38931, 1.0, int_stack+40388,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79120,int_stack+3382,int_stack+3322, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38967, 1.0, int_stack+40424,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79300,int_stack+79120,int_stack+79012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52738, 1.0, int_stack+55000,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+79012,int_stack+3532,int_stack+3472, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39303, 1.0, int_stack+40760,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+76942,int_stack+3632,int_stack+3532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39363, 1.0, int_stack+40820,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+76942,int_stack+79012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53242, 1.0, int_stack+55108,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+76942,int_stack+660,int_stack+79300,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+79012,int_stack+3872,int_stack+3782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17236, 1.0, int_stack+18437,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2442,int_stack+4022,int_stack+3872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17416, 1.0, int_stack+18617,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2892,int_stack+2442,int_stack+79012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54730, 1.0, int_stack+55504,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+80632,int_stack+2892,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+4283,int_stack+4247, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40388, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+768,int_stack+4343,int_stack+4283, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40424, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+768,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+55000, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+4493,int_stack+4433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40760, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79228,int_stack+4593,int_stack+4493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40820, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+79228,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+55108, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2802,int_stack+2442,int_stack+79012,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+79012,int_stack+4833,int_stack+4743, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+18437, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3450,int_stack+4983,int_stack+4833, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+18617, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3900,int_stack+3450,int_stack+79012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+55504, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+81712,int_stack+3900,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+5244,int_stack+5208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38931, 0.0, zero_stack, 1.0, int_stack+41845,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+5304,int_stack+5244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38967, 0.0, zero_stack, 1.0, int_stack+41881,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+2550,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52738, 0.0, zero_stack, 1.0, int_stack+58078,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+5454,int_stack+5394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39303, 0.0, zero_stack, 1.0, int_stack+42217,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79228,int_stack+5554,int_stack+5454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39363, 0.0, zero_stack, 1.0, int_stack+42277,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+79228,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53242, 0.0, zero_stack, 1.0, int_stack+58186,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3450,int_stack+660,int_stack+79012,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+79012,int_stack+5794,int_stack+5704, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17236, 0.0, zero_stack, 1.0, int_stack+19638,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4098,int_stack+5944,int_stack+5794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17416, 0.0, zero_stack, 1.0, int_stack+19818,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4548,int_stack+4098,int_stack+79012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54730, 0.0, zero_stack, 1.0, int_stack+58366,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5088,int_stack+4548,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+6205,int_stack+6169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40388, 1.0, int_stack+41845, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+768,int_stack+6265,int_stack+6205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40424, 1.0, int_stack+41881, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+768,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55000, 1.0, int_stack+58078, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+6415,int_stack+6355, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40760, 1.0, int_stack+42217, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79228,int_stack+6515,int_stack+6415, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40820, 1.0, int_stack+42277, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+79228,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55108, 1.0, int_stack+58186, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+4098,int_stack+2442,int_stack+79012,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+79012,int_stack+6755,int_stack+6665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18437, 1.0, int_stack+19638, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6168,int_stack+6905,int_stack+6755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18617, 1.0, int_stack+19818, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+6168,int_stack+79012, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55504, 1.0, int_stack+58366, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+83332,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+7166,int_stack+7130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41845, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+7226,int_stack+7166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41881, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+2550,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+58078, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+7376,int_stack+7316, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+42217, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+83008,int_stack+7476,int_stack+7376, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+42277, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+83008,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+58186, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6168,int_stack+660,int_stack+82792,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+82792,int_stack+7716,int_stack+7626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+19638, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+7866,int_stack+7716, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+19818, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+79012,int_stack+82792, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+58366, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+84412,int_stack+6816,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+8127,int_stack+8091, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38931, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43302,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+768,int_stack+8187,int_stack+8127, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38967, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43338,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+768,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52738, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60796,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+8337,int_stack+8277, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39303, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43674,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7032,int_stack+8437,int_stack+8337, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39363, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43734,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+7032,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60904,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7032,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+8677,int_stack+8587, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17236, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20839,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7680,int_stack+8827,int_stack+8677, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17416, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21019,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+7680,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54730, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+7680,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+9088,int_stack+9052, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40388, 0.0, zero_stack, 1.0, int_stack+43302, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+9148,int_stack+9088, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40424, 0.0, zero_stack, 1.0, int_stack+43338, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+2550,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55000, 0.0, zero_stack, 1.0, int_stack+60796, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+9298,int_stack+9238, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40760, 0.0, zero_stack, 1.0, int_stack+43674, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+9398,int_stack+9298, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40820, 0.0, zero_stack, 1.0, int_stack+43734, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+82792,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55108, 0.0, zero_stack, 1.0, int_stack+60904, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8760,int_stack+660,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+9638,int_stack+9548, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18437, 0.0, zero_stack, 1.0, int_stack+20839, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+9788,int_stack+9638, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18617, 0.0, zero_stack, 1.0, int_stack+21019, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55504, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+85492,int_stack+79012,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+10049,int_stack+10013, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41845, 1.0, int_stack+43302, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+768,int_stack+10109,int_stack+10049, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41881, 1.0, int_stack+43338, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+768,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58078, 1.0, int_stack+60796, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+10259,int_stack+10199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42217, 1.0, int_stack+43674, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+10359,int_stack+10259, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42277, 1.0, int_stack+43734, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+79012,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58186, 1.0, int_stack+60904, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9408,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+10599,int_stack+10509, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19638, 1.0, int_stack+20839, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+10749,int_stack+10599, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19818, 1.0, int_stack+21019, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58366, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+86572,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+11010,int_stack+10974, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+11070,int_stack+11010, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+2550,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+60796, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+11220,int_stack+11160, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43674, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+11320,int_stack+11220, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+82792,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+60904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+10056,int_stack+660,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+11560,int_stack+11470, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+20839, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+11710,int_stack+11560, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+21019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+10704,int_stack+79012,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+11971,int_stack+11935, 0.0, zero_stack, 1.0, int_stack+38931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44759,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+768,int_stack+12031,int_stack+11971, 0.0, zero_stack, 1.0, int_stack+38967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44795,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+768,int_stack+660, 0.0, zero_stack, 1.0, int_stack+52738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63064,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+12181,int_stack+12121, 0.0, zero_stack, 1.0, int_stack+39303, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45131,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+12281,int_stack+12181, 0.0, zero_stack, 1.0, int_stack+39363, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45191,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+79012,int_stack+660, 0.0, zero_stack, 1.0, int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63172,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+87652,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+12521,int_stack+12431, 0.0, zero_stack, 1.0, int_stack+17236, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22040,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+12671,int_stack+12521, 0.0, zero_stack, 1.0, int_stack+17416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22220,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+660, 0.0, zero_stack, 1.0, int_stack+54730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63568,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11784,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+12932,int_stack+12896, 0.0, zero_stack, 1.0, int_stack+40388, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44759, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+12992,int_stack+12932, 0.0, zero_stack, 1.0, int_stack+40424, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44795, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+2550,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+55000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63064, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+13142,int_stack+13082, 0.0, zero_stack, 1.0, int_stack+40760, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45131, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+13242,int_stack+13142, 0.0, zero_stack, 1.0, int_stack+40820, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45191, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+82792,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+55108, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63172, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+88300,int_stack+660,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+13482,int_stack+13392, 0.0, zero_stack, 1.0, int_stack+18437, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22040, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+13632,int_stack+13482, 0.0, zero_stack, 1.0, int_stack+18617, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22220, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+55504, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63568, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+88948,int_stack+79012,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+13893,int_stack+13857, 0.0, zero_stack, 1.0, int_stack+41845, 0.0, zero_stack, 1.0, int_stack+44759, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+768,int_stack+13953,int_stack+13893, 0.0, zero_stack, 1.0, int_stack+41881, 0.0, zero_stack, 1.0, int_stack+44795, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+768,int_stack+660, 0.0, zero_stack, 1.0, int_stack+58078, 0.0, zero_stack, 1.0, int_stack+63064, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+14103,int_stack+14043, 0.0, zero_stack, 1.0, int_stack+42217, 0.0, zero_stack, 1.0, int_stack+45131, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+14203,int_stack+14103, 0.0, zero_stack, 1.0, int_stack+42277, 0.0, zero_stack, 1.0, int_stack+45191, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+79012,int_stack+660, 0.0, zero_stack, 1.0, int_stack+58186, 0.0, zero_stack, 1.0, int_stack+63172, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+12864,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+14443,int_stack+14353, 0.0, zero_stack, 1.0, int_stack+19638, 0.0, zero_stack, 1.0, int_stack+22040, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+14593,int_stack+14443, 0.0, zero_stack, 1.0, int_stack+19818, 0.0, zero_stack, 1.0, int_stack+22220, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+660, 0.0, zero_stack, 1.0, int_stack+58366, 0.0, zero_stack, 1.0, int_stack+63568, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+13512,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+14854,int_stack+14818, 0.0, zero_stack, 1.0, int_stack+43302, 1.0, int_stack+44759, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+14914,int_stack+14854, 0.0, zero_stack, 1.0, int_stack+43338, 1.0, int_stack+44795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+2550,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+60796, 1.0, int_stack+63064, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+15064,int_stack+15004, 0.0, zero_stack, 1.0, int_stack+43674, 1.0, int_stack+45131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+15164,int_stack+15064, 0.0, zero_stack, 1.0, int_stack+43734, 1.0, int_stack+45191, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+82792,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+60904, 1.0, int_stack+63172, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14592,int_stack+660,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+15404,int_stack+15314, 0.0, zero_stack, 1.0, int_stack+20839, 1.0, int_stack+22040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+15554,int_stack+15404, 0.0, zero_stack, 1.0, int_stack+21019, 1.0, int_stack+22220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+90028,int_stack+79012,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+15815,int_stack+15779, 0.0, zero_stack, 2.0, int_stack+44759, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+768,int_stack+15875,int_stack+15815, 0.0, zero_stack, 2.0, int_stack+44795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+768,int_stack+660, 0.0, zero_stack, 2.0, int_stack+63064, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+16025,int_stack+15965, 0.0, zero_stack, 2.0, int_stack+45131, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+16125,int_stack+16025, 0.0, zero_stack, 2.0, int_stack+45191, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+79012,int_stack+660, 0.0, zero_stack, 2.0, int_stack+63172, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15240,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+16365,int_stack+16275, 0.0, zero_stack, 2.0, int_stack+22040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+16515,int_stack+16365, 0.0, zero_stack, 2.0, int_stack+22220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+660, 0.0, zero_stack, 2.0, int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+91108,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+16776,int_stack+16740, 1.0, int_stack+38931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46216,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+16836,int_stack+16776, 1.0, int_stack+38967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46252,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+2550,int_stack+2442, 1.0, int_stack+52738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66142,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+16986,int_stack+16926, 1.0, int_stack+39303, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46648,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+17086,int_stack+16986, 1.0, int_stack+39363, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46808,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+82792,int_stack+2442, 1.0, int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66250,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15888,int_stack+660,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+17566,int_stack+17326, 1.0, int_stack+17236, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23331,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+17716,int_stack+17566, 1.0, int_stack+17416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23511,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+2442, 1.0, int_stack+54730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66430,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+16536,int_stack+79012,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+17977,int_stack+17941, 1.0, int_stack+40388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46216, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+18037,int_stack+17977, 1.0, int_stack+40424, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46252, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 1.0, int_stack+55000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66142, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+18187,int_stack+18127, 1.0, int_stack+40760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46648, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+660,int_stack+18287,int_stack+18187, 1.0, int_stack+40820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46808, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+660,int_stack+53242, 1.0, int_stack+55108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66250, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+17616,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+18767,int_stack+18527, 1.0, int_stack+18437, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23331, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+18917,int_stack+18767, 1.0, int_stack+18617, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23511, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+660, 1.0, int_stack+55504, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66430, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+92188,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+19178,int_stack+19142, 1.0, int_stack+41845, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46216, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+19238,int_stack+19178, 1.0, int_stack+41881, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46252, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 1.0, int_stack+58078, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66142, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+19388,int_stack+19328, 1.0, int_stack+42217, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46648, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2442,int_stack+19488,int_stack+19388, 1.0, int_stack+42277, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46808, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+660,int_stack+2442,int_stack+53242, 1.0, int_stack+58186, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66250, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+18264,int_stack+660,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+55504,int_stack+19968,int_stack+19728, 1.0, int_stack+19638, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23331, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+20118,int_stack+19968, 1.0, int_stack+19818, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23511, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+55504, 1.0, int_stack+58366, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66430, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+18912,int_stack+79012,int_stack+660,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+20379,int_stack+20343, 1.0, int_stack+43302, 0.0, zero_stack, 1.0, int_stack+46216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+20439,int_stack+20379, 1.0, int_stack+43338, 0.0, zero_stack, 1.0, int_stack+46252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 1.0, int_stack+60796, 0.0, zero_stack, 1.0, int_stack+66142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+20589,int_stack+20529, 1.0, int_stack+43674, 0.0, zero_stack, 1.0, int_stack+46648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+660,int_stack+20689,int_stack+20589, 1.0, int_stack+43734, 0.0, zero_stack, 1.0, int_stack+46808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+660,int_stack+53242, 1.0, int_stack+60904, 0.0, zero_stack, 1.0, int_stack+66250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+19992,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+55504,int_stack+21169,int_stack+20929, 1.0, int_stack+20839, 0.0, zero_stack, 1.0, int_stack+23331, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+21319,int_stack+21169, 1.0, int_stack+21019, 0.0, zero_stack, 1.0, int_stack+23511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+55504, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+66430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+93268,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+21580,int_stack+21544, 1.0, int_stack+44759, 1.0, int_stack+46216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+21640,int_stack+21580, 1.0, int_stack+44795, 1.0, int_stack+46252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 1.0, int_stack+63064, 1.0, int_stack+66142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+21790,int_stack+21730, 1.0, int_stack+45131, 1.0, int_stack+46648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+21890,int_stack+21790, 1.0, int_stack+45191, 1.0, int_stack+46808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 1.0, int_stack+63172, 1.0, int_stack+66250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+20640,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+55504,int_stack+22370,int_stack+22130, 1.0, int_stack+22040, 1.0, int_stack+23331, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+22520,int_stack+22370, 1.0, int_stack+22220, 1.0, int_stack+23511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+55504, 1.0, int_stack+63568, 1.0, int_stack+66430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+21288,int_stack+79012,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+22781,int_stack+22745, 2.0, int_stack+46216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+22841,int_stack+22781, 2.0, int_stack+46252, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 2.0, int_stack+66142, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+22991,int_stack+22931, 2.0, int_stack+46648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+23091,int_stack+22991, 2.0, int_stack+46808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 2.0, int_stack+66250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+22368,int_stack+2442,int_stack+6816,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+23661,int_stack+23421, 2.0, int_stack+23331, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+23811,int_stack+23661, 2.0, int_stack+23511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 2.0, int_stack+66430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+94348,int_stack+82792,int_stack+2442,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+24072,int_stack+24036, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48073,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+24132,int_stack+24072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48109,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+24282,int_stack+24222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48445,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+24382,int_stack+24282, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48505,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52846,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+23016,int_stack+2442,int_stack+6816, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+24622,int_stack+24532, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29337,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+24772,int_stack+24622, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29517,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1020,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+23664,int_stack+79012,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+25033,int_stack+24997, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48073, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+25093,int_stack+25033, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48109, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+25243,int_stack+25183, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48445, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+25343,int_stack+25243, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48505, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52846, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+24744,int_stack+2442,int_stack+6816, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+25583,int_stack+25493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29337, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+25733,int_stack+25583, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29517, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1020, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+95428,int_stack+82792,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+25994,int_stack+25958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48073, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+26054,int_stack+25994, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48109, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+26204,int_stack+26144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48445, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+26304,int_stack+26204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48505, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52846, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+25392,int_stack+2442,int_stack+6816, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+26544,int_stack+26454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29337, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+26694,int_stack+26544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29517, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1020, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+96508,int_stack+79012,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+26955,int_stack+26919, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+27015,int_stack+26955, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48109, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+27165,int_stack+27105, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48445, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+27265,int_stack+27165, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+26040,int_stack+2442,int_stack+6816, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+27505,int_stack+27415, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+27655,int_stack+27505, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29517, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+26688,int_stack+82792,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+27916,int_stack+27880, 0.0, zero_stack, 1.0, int_stack+48073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+27976,int_stack+27916, 0.0, zero_stack, 1.0, int_stack+48109, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+28126,int_stack+28066, 0.0, zero_stack, 1.0, int_stack+48445, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+28226,int_stack+28126, 0.0, zero_stack, 1.0, int_stack+48505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 1.0, int_stack+52846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+97588,int_stack+2442,int_stack+6816, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+28466,int_stack+28376, 0.0, zero_stack, 1.0, int_stack+29337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+28616,int_stack+28466, 0.0, zero_stack, 1.0, int_stack+29517, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568, 0.0, zero_stack, 1.0, int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+98236,int_stack+79012,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52738,int_stack+28877,int_stack+28841, 1.0, int_stack+48073, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+28937,int_stack+28877, 1.0, int_stack+48109, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+52738, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+29087,int_stack+29027, 1.0, int_stack+48445, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+29187,int_stack+29087, 1.0, int_stack+48505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 1.0, int_stack+52846, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+27768,int_stack+2442,int_stack+6816, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+29667,int_stack+29427, 1.0, int_stack+29337, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+29817,int_stack+29667, 1.0, int_stack+29517, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 1.0, int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+28416,int_stack+82792,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+30078,int_stack+30042,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+30138,int_stack+30078,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+2442,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+30288,int_stack+30228,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+30388,int_stack+30288,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+29496,int_stack+2442,int_stack+6816, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+53422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+30628,int_stack+30538,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+30778,int_stack+30628,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+99316,int_stack+79012,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+71128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+31039,int_stack+31003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49530,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+31099,int_stack+31039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49566,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72136,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+31249,int_stack+31189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49902,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+31349,int_stack+31249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49962,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68860,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+30144,int_stack+2442,int_stack+6816, 0.0, zero_stack, 1.0, int_stack+53026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+31589,int_stack+31499, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36304,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+31739,int_stack+31589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36484,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72460,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+30792,int_stack+82792,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+53722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+32000,int_stack+31964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49530, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+32060,int_stack+32000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49566, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72136, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+32210,int_stack+32150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49902, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+32310,int_stack+32210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49962, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68860, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+100396,int_stack+2442,int_stack+6816, 0.0, zero_stack, 1.0, int_stack+55288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+32550,int_stack+32460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36304, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+32700,int_stack+32550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36484, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72460, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+101044,int_stack+79012,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+57070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+32961,int_stack+32925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49530, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+33021,int_stack+32961, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49566, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72136, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+33171,int_stack+33111, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49902, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+33271,int_stack+33171, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49962, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68860, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+31872,int_stack+2442,int_stack+6816, 0.0, zero_stack, 1.0, int_stack+55774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+33511,int_stack+33421, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36304, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+33661,int_stack+33511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36484, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72460, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+32520,int_stack+82792,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+58666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+33922,int_stack+33886, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+33982,int_stack+33922, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+2442, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+34132,int_stack+34072, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49902, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+34232,int_stack+34132, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+33600,int_stack+2442,int_stack+6816, 0.0, zero_stack, 1.0, int_stack+61084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+34472,int_stack+34382, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+34622,int_stack+34472, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+102124,int_stack+79012,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+34883,int_stack+34847, 0.0, zero_stack, 1.0, int_stack+49530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+34943,int_stack+34883, 0.0, zero_stack, 1.0, int_stack+49566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+72136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+35093,int_stack+35033, 0.0, zero_stack, 1.0, int_stack+49902, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+35193,int_stack+35093, 0.0, zero_stack, 1.0, int_stack+49962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 0.0, zero_stack, 1.0, int_stack+68860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+34248,int_stack+2442,int_stack+6816, 0.0, zero_stack, 1.0, int_stack+63352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+35433,int_stack+35343, 0.0, zero_stack, 1.0, int_stack+36304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+35583,int_stack+35433, 0.0, zero_stack, 1.0, int_stack+36484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 0.0, zero_stack, 1.0, int_stack+72460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+103204,int_stack+82792,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+65134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2442,int_stack+35844,int_stack+35808, 1.0, int_stack+49530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53242,int_stack+35904,int_stack+35844, 1.0, int_stack+49566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+53242,int_stack+2442, 1.0, int_stack+72136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+53242,int_stack+36054,int_stack+35994, 1.0, int_stack+49902, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+36154,int_stack+36054, 1.0, int_stack+49962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+0,int_stack+53242, 1.0, int_stack+68860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+34896,int_stack+2442,int_stack+6816, 0.0, zero_stack, 1.0, int_stack+63838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+36634,int_stack+36394, 1.0, int_stack+36304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+36784,int_stack+36634, 1.0, int_stack+36484, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568, 1.0, int_stack+72460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+35544,int_stack+79012,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+66730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+37045,int_stack+37009,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+37105,int_stack+37045,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+68860,int_stack+72136,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+37255,int_stack+37195,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+37355,int_stack+37255,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+72460,int_stack+68860,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+36624,int_stack+2442,int_stack+6816, 0.0, zero_stack, 1.0, int_stack+53422, 1.0, int_stack+72244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+37595,int_stack+37505,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+37745,int_stack+37595,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+104284,int_stack+82792,int_stack+2442, 0.0, zero_stack, 1.0, int_stack+71128, 1.0, int_stack+72760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+38006,int_stack+37970,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+38066,int_stack+38006,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+68860,int_stack+72136,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+38216,int_stack+38156,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+38316,int_stack+38216,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+72460,int_stack+68860,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+37272,int_stack+2442,int_stack+6816, 0.0, zero_stack, 2.0, int_stack+72244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+38556,int_stack+38466,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+38706,int_stack+38556,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+37920,int_stack+79012,int_stack+2442, 0.0, zero_stack, 2.0, int_stack+72760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+39153,int_stack+39117, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51173,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+39213,int_stack+39153, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51209,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+68860,int_stack+72136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1290,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+39673,int_stack+39613, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51545,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+39773,int_stack+39673, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51605,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+72460,int_stack+68860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1398,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+39000,int_stack+2442,int_stack+6816, 1.0, int_stack+53026, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+40013,int_stack+39923, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47368,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+40163,int_stack+40013, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47548,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69040,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+105364,int_stack+82792,int_stack+2442, 1.0, int_stack+53722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+40610,int_stack+40574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51173, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+40670,int_stack+40610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51209, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+68860,int_stack+72136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1290, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+41130,int_stack+41070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51545, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+41230,int_stack+41130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51605, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2442,int_stack+72460,int_stack+68860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1398, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+52630,int_stack+2442,int_stack+6816, 1.0, int_stack+55288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+41470,int_stack+41380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47368, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+41620,int_stack+41470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47548, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69040, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+39648,int_stack+79012,int_stack+2442, 1.0, int_stack+57070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+42067,int_stack+42031, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51173, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+42127,int_stack+42067, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51209, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+68860,int_stack+72136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1290, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+42587,int_stack+42527, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51545, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+42687,int_stack+42587, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51605, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+57070,int_stack+72460,int_stack+68860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1398, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+54730,int_stack+57070,int_stack+6816, 1.0, int_stack+55774, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+42927,int_stack+42837, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47368, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+43077,int_stack+42927, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47548, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69040, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+40728,int_stack+82792,int_stack+57070, 1.0, int_stack+58666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+43524,int_stack+43488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+43584,int_stack+43524, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51209, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+68860,int_stack+72136, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+44044,int_stack+43984, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+44144,int_stack+44044, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+57070,int_stack+72460,int_stack+68860, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+58078,int_stack+57070,int_stack+6816, 1.0, int_stack+61084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+63568,int_stack+44384,int_stack+44294, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+44534,int_stack+44384, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+63568, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+41808,int_stack+79012,int_stack+57070, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+44981,int_stack+44945, 0.0, zero_stack, 1.0, int_stack+51173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+45041,int_stack+44981, 0.0, zero_stack, 1.0, int_stack+51209, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+68860,int_stack+72136, 0.0, zero_stack, 1.0, int_stack+1290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+45501,int_stack+45441, 0.0, zero_stack, 1.0, int_stack+51545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+45601,int_stack+45501, 0.0, zero_stack, 1.0, int_stack+51605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+57070,int_stack+72460,int_stack+68860, 0.0, zero_stack, 1.0, int_stack+1398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+57070,int_stack+6816, 1.0, int_stack+63352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72460,int_stack+45841,int_stack+45751, 0.0, zero_stack, 1.0, int_stack+47368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+45991,int_stack+45841, 0.0, zero_stack, 1.0, int_stack+47548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+72460, 0.0, zero_stack, 1.0, int_stack+69040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+42888,int_stack+82792,int_stack+57070, 1.0, int_stack+65134, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+46438,int_stack+46402, 1.0, int_stack+51173, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+68860,int_stack+46498,int_stack+46438, 1.0, int_stack+51209, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+68860,int_stack+72136, 1.0, int_stack+1290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+68860,int_stack+47118,int_stack+47058, 1.0, int_stack+51545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+47218,int_stack+47118, 1.0, int_stack+51605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+65134,int_stack+72460,int_stack+68860, 1.0, int_stack+1398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+648,int_stack+65134,int_stack+6816, 1.0, int_stack+63838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+72460,int_stack+47698,int_stack+47458, 1.0, int_stack+47368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+47848,int_stack+47698, 1.0, int_stack+47548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+72460, 1.0, int_stack+69040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+43968,int_stack+79012,int_stack+65134, 1.0, int_stack+66730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+48295,int_stack+48259,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+65134,int_stack+48355,int_stack+48295,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+65134,int_stack+72136,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+65134,int_stack+48815,int_stack+48755,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+48915,int_stack+48815,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+57070,int_stack+72460,int_stack+65134,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+66142,int_stack+57070,int_stack+6816, 1.0, int_stack+53422, 0.0, zero_stack, 1.0, int_stack+1578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+65134,int_stack+49155,int_stack+49065,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+79012,int_stack+49305,int_stack+49155,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+79012,int_stack+65134,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+45048,int_stack+82792,int_stack+57070, 1.0, int_stack+71128, 0.0, zero_stack, 1.0, int_stack+75934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+72136,int_stack+49752,int_stack+49716,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+71128,int_stack+49812,int_stack+49752,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+71128,int_stack+72136,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+71128,int_stack+50272,int_stack+50212,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+72460,int_stack+50372,int_stack+50272,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+57070,int_stack+72460,int_stack+71128,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+53278,int_stack+57070,int_stack+6816, 1.0, int_stack+72244, 1.0, int_stack+1578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+71128,int_stack+50612,int_stack+50522,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+72136,int_stack+50762,int_stack+50612,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+82792,int_stack+72136,int_stack+71128,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+46128,int_stack+82792,int_stack+57070, 1.0, int_stack+72760, 1.0, int_stack+75934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+57070,int_stack+51395,int_stack+51359,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+57178,int_stack+51455,int_stack+51395,6);
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+6816,int_stack+57178,int_stack+57070,6);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+57070,int_stack+51915,int_stack+51855,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+52015,int_stack+51915,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+71128,int_stack+82792,int_stack+57070,10);
 /*--- compute (dp|dd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+72136,int_stack+71128,int_stack+6816, 2.0, int_stack+1578, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+57070,int_stack+52255,int_stack+52165,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+82792,int_stack+52405,int_stack+52255,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+79012,int_stack+82792,int_stack+57070,15);
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+47208,int_stack+79012,int_stack+71128, 2.0, int_stack+75934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+48288,int_stack+55990,int_stack+54082,36);
     Libderiv->ABCD[11] = int_stack + 48288;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+55378,int_stack+59068,int_stack+57430,36);
     Libderiv->ABCD[10] = int_stack + 55378;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+58726,int_stack+61336,int_stack+60148,36);
     Libderiv->ABCD[9] = int_stack + 58726;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+60796,int_stack+64054,int_stack+62416,36);
     Libderiv->ABCD[8] = int_stack + 60796;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+63064,int_stack+67132,int_stack+65494,36);
     Libderiv->ABCD[7] = int_stack + 63064;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+66790,int_stack+69400,int_stack+68212,36);
     Libderiv->ABCD[6] = int_stack + 66790;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+68860,int_stack+73126,int_stack+71488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[2] = int_stack + 68860;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+72784,int_stack+74854,int_stack+74206, 0.0, zero_stack, 1.0, int_stack+70480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[1] = int_stack + 72784;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+74854,int_stack+77932,int_stack+76294, 1.0, int_stack+70480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[0] = int_stack + 74854;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+70156,int_stack+79552,int_stack+1794,36);
     Libderiv->ABCD[155] = int_stack + 70156;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+1296,int_stack+80632,int_stack+76942,36);
     Libderiv->ABCD[143] = int_stack + 1296;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+76942,int_stack+81712,int_stack+2802,36);
     Libderiv->ABCD[142] = int_stack + 76942;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+78238,int_stack+5088,int_stack+3450,36);
     Libderiv->ABCD[131] = int_stack + 78238;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+79534,int_stack+83332,int_stack+4098,36);
     Libderiv->ABCD[130] = int_stack + 79534;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+80830,int_stack+84412,int_stack+6168,36);
     Libderiv->ABCD[129] = int_stack + 80830;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+82126,int_stack+7680,int_stack+7032,36);
     Libderiv->ABCD[119] = int_stack + 82126;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+83422,int_stack+85492,int_stack+8760,36);
     Libderiv->ABCD[118] = int_stack + 83422;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+84718,int_stack+86572,int_stack+9408,36);
     Libderiv->ABCD[117] = int_stack + 84718;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+86014,int_stack+10704,int_stack+10056,36);
     Libderiv->ABCD[116] = int_stack + 86014;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+2592,int_stack+11784,int_stack+87652,36);
     Libderiv->ABCD[107] = int_stack + 2592;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+3888,int_stack+88948,int_stack+88300,36);
     Libderiv->ABCD[106] = int_stack + 3888;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+87310,int_stack+13512,int_stack+12864,36);
     Libderiv->ABCD[105] = int_stack + 87310;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+88606,int_stack+90028,int_stack+14592,36);
     Libderiv->ABCD[104] = int_stack + 88606;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+5184,int_stack+91108,int_stack+15240,36);
     Libderiv->ABCD[103] = int_stack + 5184;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+89902,int_stack+16536,int_stack+15888,36);
     Libderiv->ABCD[95] = int_stack + 89902;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+6480,int_stack+92188,int_stack+17616,36);
     Libderiv->ABCD[94] = int_stack + 6480;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+91198,int_stack+18912,int_stack+18264,36);
     Libderiv->ABCD[93] = int_stack + 91198;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+7776,int_stack+93268,int_stack+19992,36);
     Libderiv->ABCD[92] = int_stack + 7776;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+92494,int_stack+21288,int_stack+20640,36);
     Libderiv->ABCD[91] = int_stack + 92494;
 /*--- compute (dd|dd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+9072,int_stack+94348,int_stack+22368,36);
     Libderiv->ABCD[90] = int_stack + 9072;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+93790,int_stack+23664,int_stack+23016, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[47] = int_stack + 93790;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+10368,int_stack+95428,int_stack+24744, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[46] = int_stack + 10368;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+95086,int_stack+96508,int_stack+25392, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[45] = int_stack + 95086;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+11664,int_stack+26688,int_stack+26040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[44] = int_stack + 11664;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+12960,int_stack+98236,int_stack+97588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[43] = int_stack + 12960;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+96382,int_stack+28416,int_stack+27768, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[42] = int_stack + 96382;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+97678,int_stack+99316,int_stack+29496, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+71488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[38] = int_stack + 97678;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+98974,int_stack+30792,int_stack+30144, 0.0, zero_stack, 1.0, int_stack+54082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[35] = int_stack + 98974;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+14256,int_stack+101044,int_stack+100396, 0.0, zero_stack, 1.0, int_stack+57430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[34] = int_stack + 14256;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+100270,int_stack+32520,int_stack+31872, 0.0, zero_stack, 1.0, int_stack+60148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[33] = int_stack + 100270;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+15552,int_stack+102124,int_stack+33600, 0.0, zero_stack, 1.0, int_stack+62416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[32] = int_stack + 15552;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+101566,int_stack+103204,int_stack+34248, 0.0, zero_stack, 1.0, int_stack+65494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[31] = int_stack + 101566;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+102862,int_stack+35544,int_stack+34896, 0.0, zero_stack, 1.0, int_stack+68212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[30] = int_stack + 102862;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+16848,int_stack+104284,int_stack+36624, 0.0, zero_stack, 1.0, int_stack+71488, 1.0, int_stack+74206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[26] = int_stack + 16848;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+18144,int_stack+37920,int_stack+37272, 0.0, zero_stack, 2.0, int_stack+74206, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[25] = int_stack + 18144;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+19440,int_stack+105364,int_stack+39000, 1.0, int_stack+54082, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[23] = int_stack + 19440;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+104158,int_stack+39648,int_stack+52630, 1.0, int_stack+57430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[22] = int_stack + 104158;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+105454,int_stack+40728,int_stack+54730, 1.0, int_stack+60148, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[21] = int_stack + 105454;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+20736,int_stack+41808,int_stack+58078, 1.0, int_stack+62416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[20] = int_stack + 20736;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+22032,int_stack+42888,int_stack+0, 1.0, int_stack+65494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[19] = int_stack + 22032;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+23328,int_stack+43968,int_stack+648, 1.0, int_stack+68212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[18] = int_stack + 23328;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+45048,int_stack+66142, 1.0, int_stack+71488, 0.0, zero_stack, 1.0, int_stack+76294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[14] = int_stack + 0;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+24624,int_stack+46128,int_stack+53278, 1.0, int_stack+74206, 1.0, int_stack+76294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[13] = int_stack + 24624;
 /*--- compute (dd|dd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+25920,int_stack+47208,int_stack+72136, 2.0, int_stack+76294, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[12] = int_stack + 25920;

}
