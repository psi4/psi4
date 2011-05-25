#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ddfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|fd) integrals */

void d12hrr_order_ddfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][5][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][10] = int_stack + 315;
 Libderiv->deriv_classes[4][5][9] = int_stack + 630;
 Libderiv->deriv_classes[4][5][8] = int_stack + 945;
 Libderiv->deriv_classes[4][5][7] = int_stack + 1260;
 Libderiv->dvrr_classes[4][4] = int_stack + 1575;
 Libderiv->deriv_classes[4][5][6] = int_stack + 1800;
 Libderiv->deriv_classes[4][5][2] = int_stack + 2115;
 Libderiv->deriv_classes[4][5][1] = int_stack + 2430;
 Libderiv->dvrr_classes[3][5] = int_stack + 2745;
 Libderiv->deriv_classes[4][5][0] = int_stack + 2955;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 3270;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 3330;
 Libderiv->deriv2_classes[2][5][143] = int_stack + 3420;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 3546;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 3646;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 3796;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 4006;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 4156;
 Libderiv->deriv2_classes[4][5][143] = int_stack + 4381;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 4696;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 4756;
 Libderiv->deriv2_classes[2][5][131] = int_stack + 4846;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 4972;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 5072;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 5222;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 5432;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 5582;
 Libderiv->deriv2_classes[4][5][131] = int_stack + 5807;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 6122;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 6182;
 Libderiv->deriv2_classes[2][5][130] = int_stack + 6272;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 6398;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 6498;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 6648;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 6858;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 7008;
 Libderiv->deriv2_classes[4][5][130] = int_stack + 7233;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 7548;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 7608;
 Libderiv->deriv2_classes[2][5][119] = int_stack + 7698;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 7824;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 7924;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 8074;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 8284;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 8434;
 Libderiv->deriv2_classes[4][5][119] = int_stack + 8659;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 8974;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 9034;
 Libderiv->deriv2_classes[2][5][118] = int_stack + 9124;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 9250;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 9350;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 9500;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 9710;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 9860;
 Libderiv->deriv2_classes[4][5][118] = int_stack + 10085;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 10400;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 10460;
 Libderiv->deriv2_classes[2][5][117] = int_stack + 10550;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 10676;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 10776;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 10926;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 11136;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 11286;
 Libderiv->deriv2_classes[4][5][117] = int_stack + 11511;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 11826;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 11886;
 Libderiv->deriv2_classes[2][5][107] = int_stack + 11976;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 12102;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 12202;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 12352;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 12562;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 12712;
 Libderiv->deriv2_classes[4][5][107] = int_stack + 12937;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 13252;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 13312;
 Libderiv->deriv2_classes[2][5][106] = int_stack + 13402;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 13528;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 13628;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 13778;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 13988;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 14138;
 Libderiv->deriv2_classes[4][5][106] = int_stack + 14363;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 14678;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 14738;
 Libderiv->deriv2_classes[2][5][105] = int_stack + 14828;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 14954;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 15054;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 15204;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 15414;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 15564;
 Libderiv->deriv2_classes[4][5][105] = int_stack + 15789;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 16104;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 16164;
 Libderiv->deriv2_classes[2][5][104] = int_stack + 16254;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 16380;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 16480;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 16630;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 16840;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 16990;
 Libderiv->deriv2_classes[4][5][104] = int_stack + 17215;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 17530;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 17590;
 Libderiv->deriv2_classes[2][5][95] = int_stack + 17680;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 17806;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 17906;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 18056;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 18266;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 18416;
 Libderiv->deriv2_classes[4][5][95] = int_stack + 18641;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 18956;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 19016;
 Libderiv->deriv2_classes[2][5][94] = int_stack + 19106;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 19232;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 19332;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 19482;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 19692;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 19842;
 Libderiv->deriv2_classes[4][5][94] = int_stack + 20067;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 20382;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 20442;
 Libderiv->deriv2_classes[2][5][93] = int_stack + 20532;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 20658;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 20758;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 20908;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 21118;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 21268;
 Libderiv->deriv2_classes[4][5][93] = int_stack + 21493;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 21808;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 21868;
 Libderiv->deriv2_classes[2][5][92] = int_stack + 21958;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 22084;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 22184;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 22334;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 22544;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 22694;
 Libderiv->deriv2_classes[4][5][92] = int_stack + 22919;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 23234;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 23294;
 Libderiv->deriv2_classes[2][5][91] = int_stack + 23384;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 23510;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 23610;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 23760;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 23970;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 24120;
 Libderiv->deriv2_classes[4][5][91] = int_stack + 24345;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 24660;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 24720;
 Libderiv->deriv2_classes[2][5][83] = int_stack + 24810;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 24936;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 25036;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 25186;
 Libderiv->deriv_classes[4][3][11] = int_stack + 25396;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 25546;
 Libderiv->deriv_classes[4][4][11] = int_stack + 25696;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 25921;
 Libderiv->deriv2_classes[4][5][83] = int_stack + 26146;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 26461;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 26521;
 Libderiv->deriv2_classes[2][5][82] = int_stack + 26611;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 26737;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 26837;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 26987;
 Libderiv->deriv_classes[4][3][10] = int_stack + 27197;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 27347;
 Libderiv->deriv_classes[4][4][10] = int_stack + 27497;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 27722;
 Libderiv->deriv2_classes[4][5][82] = int_stack + 27947;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 28262;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 28322;
 Libderiv->deriv2_classes[2][5][81] = int_stack + 28412;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 28538;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 28638;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 28788;
 Libderiv->deriv_classes[4][3][9] = int_stack + 28998;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 29148;
 Libderiv->deriv_classes[4][4][9] = int_stack + 29298;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 29523;
 Libderiv->deriv2_classes[4][5][81] = int_stack + 29748;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 30063;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 30123;
 Libderiv->deriv2_classes[2][5][80] = int_stack + 30213;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 30339;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 30439;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 30589;
 Libderiv->deriv_classes[4][3][8] = int_stack + 30799;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 30949;
 Libderiv->deriv_classes[4][4][8] = int_stack + 31099;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 31324;
 Libderiv->deriv2_classes[4][5][80] = int_stack + 31549;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 31864;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 31924;
 Libderiv->deriv2_classes[2][5][79] = int_stack + 32014;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 32140;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 32240;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 32390;
 Libderiv->deriv_classes[4][3][7] = int_stack + 32600;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 32750;
 Libderiv->deriv_classes[4][4][7] = int_stack + 32900;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 33125;
 Libderiv->deriv2_classes[4][5][79] = int_stack + 33350;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 33665;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 33725;
 Libderiv->deriv2_classes[2][5][78] = int_stack + 33815;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 33941;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 34041;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 34191;
 Libderiv->dvrr_classes[4][3] = int_stack + 34401;
 Libderiv->deriv_classes[4][3][6] = int_stack + 34551;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 34701;
 Libderiv->deriv_classes[4][4][6] = int_stack + 34851;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 35076;
 Libderiv->deriv2_classes[4][5][78] = int_stack + 35301;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 35616;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 35676;
 Libderiv->deriv2_classes[2][5][35] = int_stack + 35766;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 35892;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 35992;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 36142;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 36352;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 36502;
 Libderiv->deriv2_classes[4][5][35] = int_stack + 36727;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 37042;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 37102;
 Libderiv->deriv2_classes[2][5][34] = int_stack + 37192;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 37318;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 37418;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 37568;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 37778;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 37928;
 Libderiv->deriv2_classes[4][5][34] = int_stack + 38153;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 38468;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 38528;
 Libderiv->deriv2_classes[2][5][33] = int_stack + 38618;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 38744;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 38844;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 38994;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 39204;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 39354;
 Libderiv->deriv2_classes[4][5][33] = int_stack + 39579;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 39894;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 39954;
 Libderiv->deriv2_classes[2][5][32] = int_stack + 40044;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 40170;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 40270;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 40420;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 40630;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 40780;
 Libderiv->deriv2_classes[4][5][32] = int_stack + 41005;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 41320;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 41380;
 Libderiv->deriv2_classes[2][5][31] = int_stack + 41470;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 41596;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 41696;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 41846;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 42056;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 42206;
 Libderiv->deriv2_classes[4][5][31] = int_stack + 42431;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 42746;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 42806;
 Libderiv->deriv2_classes[2][5][30] = int_stack + 42896;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 43022;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 43122;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 43272;
 Libderiv->deriv_classes[4][3][2] = int_stack + 43482;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 43632;
 Libderiv->deriv_classes[4][4][2] = int_stack + 43782;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 44007;
 Libderiv->deriv2_classes[4][5][30] = int_stack + 44232;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 44547;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 44607;
 Libderiv->deriv2_classes[2][5][26] = int_stack + 44697;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 44823;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 44923;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 45073;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 45283;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 45433;
 Libderiv->deriv2_classes[4][5][26] = int_stack + 45658;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 45973;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 46033;
 Libderiv->deriv2_classes[2][5][23] = int_stack + 46123;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 46249;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 46349;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 46499;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 46709;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 46859;
 Libderiv->deriv2_classes[4][5][23] = int_stack + 47084;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 47399;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 47459;
 Libderiv->deriv2_classes[2][5][22] = int_stack + 47549;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 47675;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 47775;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 47925;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 48135;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 48285;
 Libderiv->deriv2_classes[4][5][22] = int_stack + 48510;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 48825;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 48885;
 Libderiv->deriv2_classes[2][5][21] = int_stack + 48975;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 49101;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 49201;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 49351;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 49561;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 49711;
 Libderiv->deriv2_classes[4][5][21] = int_stack + 49936;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 50251;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 50311;
 Libderiv->deriv2_classes[2][5][20] = int_stack + 50401;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 50527;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 50627;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 50777;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 50987;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 51137;
 Libderiv->deriv2_classes[4][5][20] = int_stack + 51362;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 51677;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 51737;
 Libderiv->deriv2_classes[2][5][19] = int_stack + 51827;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 51953;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 52053;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 52203;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 52413;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 52563;
 Libderiv->deriv2_classes[4][5][19] = int_stack + 52788;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 53103;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 53163;
 Libderiv->deriv2_classes[2][5][18] = int_stack + 53253;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 53379;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 53479;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 53629;
 Libderiv->deriv_classes[4][3][1] = int_stack + 53839;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 53989;
 Libderiv->deriv_classes[4][4][1] = int_stack + 54139;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 54364;
 Libderiv->deriv2_classes[4][5][18] = int_stack + 54589;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 54904;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 54964;
 Libderiv->deriv2_classes[2][5][14] = int_stack + 55054;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 55180;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 55280;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 55430;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 55640;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 55790;
 Libderiv->deriv2_classes[4][5][14] = int_stack + 56015;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 56330;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 56390;
 Libderiv->deriv2_classes[2][5][13] = int_stack + 56480;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 56606;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 56706;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 56856;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 57066;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 57216;
 Libderiv->deriv2_classes[4][5][13] = int_stack + 57441;
 Libderiv->deriv_classes[2][3][11] = int_stack + 57756;
 Libderiv->deriv_classes[2][4][11] = int_stack + 57816;
 Libderiv->deriv_classes[2][5][11] = int_stack + 57906;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 58032;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 58092;
 Libderiv->deriv2_classes[2][5][11] = int_stack + 58182;
 Libderiv->deriv_classes[3][3][11] = int_stack + 58308;
 Libderiv->deriv_classes[3][4][11] = int_stack + 58408;
 Libderiv->deriv_classes[3][5][11] = int_stack + 58558;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 58768;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 58868;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 59018;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 59228;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 59378;
 Libderiv->deriv2_classes[4][5][11] = int_stack + 59603;
 Libderiv->deriv_classes[2][3][10] = int_stack + 59918;
 Libderiv->deriv_classes[2][4][10] = int_stack + 59978;
 Libderiv->deriv_classes[2][5][10] = int_stack + 60068;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 60194;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 60254;
 Libderiv->deriv2_classes[2][5][10] = int_stack + 60344;
 Libderiv->deriv_classes[3][3][10] = int_stack + 60470;
 Libderiv->deriv_classes[3][4][10] = int_stack + 60570;
 Libderiv->deriv_classes[3][5][10] = int_stack + 60720;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 60930;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 61030;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 61180;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 61390;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 61540;
 Libderiv->deriv2_classes[4][5][10] = int_stack + 61765;
 Libderiv->deriv_classes[2][3][9] = int_stack + 62080;
 Libderiv->deriv_classes[2][4][9] = int_stack + 62140;
 Libderiv->deriv_classes[2][5][9] = int_stack + 62230;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 62356;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 62416;
 Libderiv->deriv2_classes[2][5][9] = int_stack + 62506;
 Libderiv->deriv_classes[3][3][9] = int_stack + 62632;
 Libderiv->deriv_classes[3][4][9] = int_stack + 62732;
 Libderiv->deriv_classes[3][5][9] = int_stack + 62882;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 63092;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 63192;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 63342;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 63552;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 63702;
 Libderiv->deriv2_classes[4][5][9] = int_stack + 63927;
 Libderiv->deriv_classes[2][3][8] = int_stack + 64242;
 Libderiv->deriv_classes[2][4][8] = int_stack + 64302;
 Libderiv->deriv_classes[2][5][8] = int_stack + 64392;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 64518;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 64578;
 Libderiv->deriv2_classes[2][5][8] = int_stack + 64668;
 Libderiv->deriv_classes[3][3][8] = int_stack + 64794;
 Libderiv->deriv_classes[3][4][8] = int_stack + 64894;
 Libderiv->deriv_classes[3][5][8] = int_stack + 65044;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 65254;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 65354;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 65504;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 65714;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 65864;
 Libderiv->deriv2_classes[4][5][8] = int_stack + 66089;
 Libderiv->deriv_classes[2][3][7] = int_stack + 66404;
 Libderiv->deriv_classes[2][4][7] = int_stack + 66464;
 Libderiv->deriv_classes[2][5][7] = int_stack + 66554;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 66680;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 66740;
 Libderiv->deriv2_classes[2][5][7] = int_stack + 66830;
 Libderiv->deriv_classes[3][3][7] = int_stack + 66956;
 Libderiv->deriv_classes[3][4][7] = int_stack + 67056;
 Libderiv->deriv_classes[3][5][7] = int_stack + 67206;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 67416;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 67516;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 67666;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 67876;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 68026;
 Libderiv->deriv2_classes[4][5][7] = int_stack + 68251;
 Libderiv->deriv_classes[2][3][6] = int_stack + 68566;
 Libderiv->deriv_classes[2][4][6] = int_stack + 68626;
 Libderiv->deriv_classes[2][5][6] = int_stack + 68716;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 68842;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 68902;
 Libderiv->deriv2_classes[2][5][6] = int_stack + 68992;
 Libderiv->dvrr_classes[3][3] = int_stack + 69118;
 Libderiv->deriv_classes[3][3][6] = int_stack + 69218;
 Libderiv->dvrr_classes[3][4] = int_stack + 69318;
 Libderiv->deriv_classes[3][4][6] = int_stack + 69468;
 Libderiv->deriv_classes[3][5][6] = int_stack + 69618;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 69828;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 69928;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 70078;
 Libderiv->deriv_classes[4][3][0] = int_stack + 70288;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 70438;
 Libderiv->deriv_classes[4][4][0] = int_stack + 70588;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 70813;
 Libderiv->deriv2_classes[4][5][6] = int_stack + 71038;
 Libderiv->deriv_classes[2][3][2] = int_stack + 71353;
 Libderiv->deriv_classes[2][4][2] = int_stack + 71413;
 Libderiv->deriv_classes[2][5][2] = int_stack + 71503;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 71629;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 71689;
 Libderiv->deriv2_classes[2][5][2] = int_stack + 71779;
 Libderiv->deriv_classes[3][3][2] = int_stack + 71905;
 Libderiv->deriv_classes[3][4][2] = int_stack + 72005;
 Libderiv->deriv_classes[3][5][2] = int_stack + 72155;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 72365;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 72465;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 72615;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 72825;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 72975;
 Libderiv->deriv2_classes[4][5][2] = int_stack + 73200;
 Libderiv->deriv_classes[2][3][1] = int_stack + 73515;
 Libderiv->deriv_classes[2][4][1] = int_stack + 73575;
 Libderiv->deriv_classes[2][5][1] = int_stack + 73665;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 73791;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 73851;
 Libderiv->deriv2_classes[2][5][1] = int_stack + 73941;
 Libderiv->deriv_classes[3][3][1] = int_stack + 74067;
 Libderiv->deriv_classes[3][4][1] = int_stack + 74167;
 Libderiv->deriv_classes[3][5][1] = int_stack + 74317;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 74527;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 74627;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 74777;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 74987;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 75137;
 Libderiv->deriv2_classes[4][5][1] = int_stack + 75362;
 Libderiv->dvrr_classes[2][3] = int_stack + 75677;
 Libderiv->dvrr_classes[2][4] = int_stack + 75737;
 Libderiv->dvrr_classes[2][5] = int_stack + 75827;
 Libderiv->deriv_classes[2][3][0] = int_stack + 75953;
 Libderiv->deriv_classes[2][4][0] = int_stack + 76013;
 Libderiv->deriv_classes[2][5][0] = int_stack + 76103;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 76229;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 76289;
 Libderiv->deriv2_classes[2][5][0] = int_stack + 76379;
 Libderiv->deriv_classes[3][3][0] = int_stack + 76505;
 Libderiv->deriv_classes[3][4][0] = int_stack + 76605;
 Libderiv->deriv_classes[3][5][0] = int_stack + 76755;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 76965;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 77065;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 77215;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 77425;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 77575;
 Libderiv->deriv2_classes[4][5][0] = int_stack + 77800;
 memset(int_stack,0,624920);

 Libderiv->dvrr_stack = int_stack + 177475;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ddfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+78115,int_stack+75737,int_stack+75677,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+57816,int_stack+57756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75677,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+57906,int_stack+57816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75737,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+78745,int_stack+78475,int_stack+78295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78115,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+79105,int_stack+69318,int_stack+69118,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+58408,int_stack+58308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69118,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79705,int_stack+58558,int_stack+58408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69318,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+80155,int_stack+79705,int_stack+79405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79105,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+80755,int_stack+80155,int_stack+78745,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+79705,int_stack+1575,int_stack+34401,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+81835,int_stack+25696,int_stack+25396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34401,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82285,int_stack+0,int_stack+25696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1575,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+82960,int_stack+82285,int_stack+81835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79705,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+83860,int_stack+82960,int_stack+80155,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82285,int_stack+59978,int_stack+59918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75677, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+60068,int_stack+59978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75737, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+82465,int_stack+78475,int_stack+82285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78115, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82825,int_stack+60570,int_stack+60470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69118, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+83125,int_stack+60720,int_stack+60570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69318, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+85660,int_stack+83125,int_stack+82825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79105, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+86260,int_stack+85660,int_stack+82465,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+83125,int_stack+27497,int_stack+27197, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34401, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+87340,int_stack+315,int_stack+27497, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1575, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+88015,int_stack+87340,int_stack+83125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79705, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+88915,int_stack+88015,int_stack+85660,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+87340,int_stack+62140,int_stack+62080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75677, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+62230,int_stack+62140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75737, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+87520,int_stack+78475,int_stack+87340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78115, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+87880,int_stack+62732,int_stack+62632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69118, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+88180,int_stack+62882,int_stack+62732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69318, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+88180,int_stack+87880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79105, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+90715,int_stack+0,int_stack+87520,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+88180,int_stack+29298,int_stack+28998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34401, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+91795,int_stack+630,int_stack+29298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+92470,int_stack+91795,int_stack+88180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79705, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+93370,int_stack+92470,int_stack+0,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+91795,int_stack+64302,int_stack+64242, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75677, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+64392,int_stack+64302, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75737, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+91975,int_stack+78475,int_stack+91795, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+92335,int_stack+64894,int_stack+64794, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+92635,int_stack+65044,int_stack+64894, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+95170,int_stack+92635,int_stack+92335, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+95770,int_stack+95170,int_stack+91975,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+92635,int_stack+31099,int_stack+30799, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+96850,int_stack+945,int_stack+31099, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97525,int_stack+96850,int_stack+92635, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+98425,int_stack+97525,int_stack+95170,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+96850,int_stack+66464,int_stack+66404, 0.0, zero_stack, 1.0, int_stack+75677, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+66554,int_stack+66464, 0.0, zero_stack, 1.0, int_stack+75737, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97030,int_stack+78475,int_stack+96850, 0.0, zero_stack, 1.0, int_stack+78115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97390,int_stack+67056,int_stack+66956, 0.0, zero_stack, 1.0, int_stack+69118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97690,int_stack+67206,int_stack+67056, 0.0, zero_stack, 1.0, int_stack+69318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+600,int_stack+97690,int_stack+97390, 0.0, zero_stack, 1.0, int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+100225,int_stack+600,int_stack+97030,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97690,int_stack+32900,int_stack+32600, 0.0, zero_stack, 1.0, int_stack+34401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101305,int_stack+1260,int_stack+32900, 0.0, zero_stack, 1.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101980,int_stack+101305,int_stack+97690, 0.0, zero_stack, 1.0, int_stack+79705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+102880,int_stack+101980,int_stack+600,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+68626,int_stack+68566, 1.0, int_stack+75677, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+68716,int_stack+68626, 1.0, int_stack+75737, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101485,int_stack+78475,int_stack+101305, 1.0, int_stack+78115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101845,int_stack+69468,int_stack+69218, 1.0, int_stack+69118, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102145,int_stack+69618,int_stack+69468, 1.0, int_stack+69318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+104680,int_stack+102145,int_stack+101845, 1.0, int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+105280,int_stack+104680,int_stack+101485,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102145,int_stack+34851,int_stack+34551, 1.0, int_stack+34401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+106360,int_stack+1800,int_stack+34851, 1.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+106360,int_stack+102145, 1.0, int_stack+79705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+106360,int_stack+1200,int_stack+104680,60);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+75827,int_stack+75737,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+78115,6);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+79705,int_stack+2745,int_stack+69318,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+79705,int_stack+79105,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+108760,int_stack+108160,int_stack+1200,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+78115,int_stack+71413,int_stack+71353,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+71503,int_stack+71413,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+79705,int_stack+78475,int_stack+78115,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+79105,int_stack+72005,int_stack+71905,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1560,int_stack+72155,int_stack+72005,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+109840,int_stack+1560,int_stack+79105,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+110440,int_stack+109840,int_stack+79705, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1560,int_stack+43782,int_stack+43482,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+111520,int_stack+2115,int_stack+43782,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+112195,int_stack+111520,int_stack+1560,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+113095,int_stack+112195,int_stack+109840, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+73575,int_stack+73515,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+73665,int_stack+73575,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+111700,int_stack+78475,int_stack+111520,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+112060,int_stack+74167,int_stack+74067,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+112360,int_stack+74317,int_stack+74167,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+114895,int_stack+112360,int_stack+112060,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+115495,int_stack+114895,int_stack+111700, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+112360,int_stack+54139,int_stack+53839,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+116575,int_stack+2430,int_stack+54139,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2010,int_stack+116575,int_stack+112360,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+116575,int_stack+2010,int_stack+114895, 0.0, zero_stack, 1.0, int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2010,int_stack+76013,int_stack+75953,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+76103,int_stack+76013,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2190,int_stack+78475,int_stack+2010,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+76605,int_stack+76505,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+118375,int_stack+76755,int_stack+76605,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+118825,int_stack+118375,int_stack+2550,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+119425,int_stack+118825,int_stack+2190, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+118375,int_stack+70588,int_stack+70288,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+120505,int_stack+2955,int_stack+70588,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+121180,int_stack+120505,int_stack+118375,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+122080,int_stack+121180,int_stack+118825, 1.0, int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+3330,int_stack+3270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+57756,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+3420,int_stack+3330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+57816,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+78295,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+3646,int_stack+3546, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+58308,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+120505,int_stack+3796,int_stack+3646, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+58408,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+120955,int_stack+120505,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+79405,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2850,int_stack+120955,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+120505,int_stack+4156,int_stack+4006, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25396,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+123880,int_stack+4381,int_stack+4156, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25696,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+124555,int_stack+123880,int_stack+120505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+81835,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+125455,int_stack+124555,int_stack+120955,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+120505,int_stack+4756,int_stack+4696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57756, 1.0, int_stack+59918,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+4846,int_stack+4756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57816, 1.0, int_stack+59978,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+120505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78295, 1.0, int_stack+82285,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+120505,int_stack+5072,int_stack+4972, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58308, 1.0, int_stack+60470,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+120805,int_stack+5222,int_stack+5072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58408, 1.0, int_stack+60570,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+120805,int_stack+120505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79405, 1.0, int_stack+82825,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+120505,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+121585,int_stack+5582,int_stack+5432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25396, 1.0, int_stack+27197,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+123880,int_stack+5807,int_stack+5582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25696, 1.0, int_stack+27497,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+124555,int_stack+123880,int_stack+121585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81835, 1.0, int_stack+83125,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3930,int_stack+124555,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+6182,int_stack+6122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+59918, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+6272,int_stack+6182, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+59978, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+82285, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+6498,int_stack+6398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+60470, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+121585,int_stack+6648,int_stack+6498, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+60570, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+123880,int_stack+121585,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+82825, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5730,int_stack+123880,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+7008,int_stack+6858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27197, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+124480,int_stack+7233,int_stack+7008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27497, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+124480,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+83125, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+128155,int_stack+127255,int_stack+123880,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+123880,int_stack+7608,int_stack+7548, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57756, 0.0, zero_stack, 1.0, int_stack+62080,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+7698,int_stack+7608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57816, 0.0, zero_stack, 1.0, int_stack+62140,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+123880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78295, 0.0, zero_stack, 1.0, int_stack+87340,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+123880,int_stack+7924,int_stack+7824, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58308, 0.0, zero_stack, 1.0, int_stack+62632,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+124180,int_stack+8074,int_stack+7924, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58408, 0.0, zero_stack, 1.0, int_stack+62732,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+124180,int_stack+123880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79405, 0.0, zero_stack, 1.0, int_stack+87880,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+123880,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+124960,int_stack+8434,int_stack+8284, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25396, 0.0, zero_stack, 1.0, int_stack+28998,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+8659,int_stack+8434, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25696, 0.0, zero_stack, 1.0, int_stack+29298,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6810,int_stack+127255,int_stack+124960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81835, 0.0, zero_stack, 1.0, int_stack+88180,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+129955,int_stack+6810,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+9034,int_stack+8974, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59918, 1.0, int_stack+62080, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+9124,int_stack+9034, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59978, 1.0, int_stack+62140, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82285, 1.0, int_stack+87340, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+9350,int_stack+9250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60470, 1.0, int_stack+62632, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6810,int_stack+9500,int_stack+9350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60570, 1.0, int_stack+62732, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7260,int_stack+6810,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82825, 1.0, int_stack+87880, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7860,int_stack+7260,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6810,int_stack+9860,int_stack+9710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27197, 1.0, int_stack+28998, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8940,int_stack+10085,int_stack+9860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27497, 1.0, int_stack+29298, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+8940,int_stack+6810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+83125, 1.0, int_stack+88180, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+131755,int_stack+127255,int_stack+7260,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+127255,int_stack+10460,int_stack+10400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+62080, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+10550,int_stack+10460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+62140, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+127255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+87340, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+127255,int_stack+10776,int_stack+10676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+62632, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127555,int_stack+10926,int_stack+10776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+62732, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+127555,int_stack+127255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+87880, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8940,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+127255,int_stack+11286,int_stack+11136, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28998, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6810,int_stack+11511,int_stack+11286, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+29298, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10020,int_stack+6810,int_stack+127255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+88180, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+133555,int_stack+10020,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+11886,int_stack+11826, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57756, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64242,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+11976,int_stack+11886, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57816, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64302,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78295, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91795,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+12202,int_stack+12102, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58308, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64794,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10020,int_stack+12352,int_stack+12202, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58408, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64894,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10470,int_stack+10020,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79405, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92335,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+11070,int_stack+10470,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10020,int_stack+12712,int_stack+12562, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25396, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30799,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+12937,int_stack+12712, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25696, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31099,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12150,int_stack+127255,int_stack+10020, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+81835, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92635,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+135355,int_stack+12150,int_stack+10470,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12150,int_stack+13312,int_stack+13252, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59918, 0.0, zero_stack, 1.0, int_stack+64242, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+13402,int_stack+13312, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59978, 0.0, zero_stack, 1.0, int_stack+64302, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+12150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82285, 0.0, zero_stack, 1.0, int_stack+91795, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12150,int_stack+13628,int_stack+13528, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60470, 0.0, zero_stack, 1.0, int_stack+64794, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12450,int_stack+13778,int_stack+13628, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+60570, 0.0, zero_stack, 1.0, int_stack+64894, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+12450,int_stack+12150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82825, 0.0, zero_stack, 1.0, int_stack+92335, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+12150,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+14138,int_stack+13988, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27197, 0.0, zero_stack, 1.0, int_stack+30799, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10020,int_stack+14363,int_stack+14138, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27497, 0.0, zero_stack, 1.0, int_stack+31099, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+10020,int_stack+13230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+83125, 0.0, zero_stack, 1.0, int_stack+92635, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+137155,int_stack+127255,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+14738,int_stack+14678, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62080, 1.0, int_stack+64242, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+14828,int_stack+14738, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62140, 1.0, int_stack+64302, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87340, 1.0, int_stack+91795, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+15054,int_stack+14954, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62632, 1.0, int_stack+64794, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+15204,int_stack+15054, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62732, 1.0, int_stack+64894, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+127255,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87880, 1.0, int_stack+92335, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+13830,int_stack+13230,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+15564,int_stack+15414, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28998, 1.0, int_stack+30799, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+15789,int_stack+15564, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29298, 1.0, int_stack+31099, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14910,int_stack+127255,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88180, 1.0, int_stack+92635, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+138955,int_stack+14910,int_stack+13230,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+16164,int_stack+16104, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+64242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+16254,int_stack+16164, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+64302, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+13230, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+91795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+16480,int_stack+16380, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+64794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14910,int_stack+16630,int_stack+16480, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+64894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+14910,int_stack+13230, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+92335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14910,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+16990,int_stack+16840, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+30799, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15990,int_stack+17215,int_stack+16990, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+31099, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+15990,int_stack+13230, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+92635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+140755,int_stack+127255,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+17590,int_stack+17530, 0.0, zero_stack, 1.0, int_stack+57756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66404,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+17680,int_stack+17590, 0.0, zero_stack, 1.0, int_stack+57816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66464,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+78295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96850,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+17906,int_stack+17806, 0.0, zero_stack, 1.0, int_stack+58308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66956,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+18056,int_stack+17906, 0.0, zero_stack, 1.0, int_stack+58408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67056,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+127255,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+79405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97390,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15990,int_stack+13230,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+18416,int_stack+18266, 0.0, zero_stack, 1.0, int_stack+25396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32600,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+18641,int_stack+18416, 0.0, zero_stack, 1.0, int_stack+25696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32900,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17070,int_stack+127255,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+81835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97690,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+142555,int_stack+17070,int_stack+13230,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+19016,int_stack+18956, 0.0, zero_stack, 1.0, int_stack+59918, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66404, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+19106,int_stack+19016, 0.0, zero_stack, 1.0, int_stack+59978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66464, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+13230, 0.0, zero_stack, 1.0, int_stack+82285, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96850, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+19332,int_stack+19232, 0.0, zero_stack, 1.0, int_stack+60470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66956, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17070,int_stack+19482,int_stack+19332, 0.0, zero_stack, 1.0, int_stack+60570, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67056, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+17070,int_stack+13230, 0.0, zero_stack, 1.0, int_stack+82825, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97390, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+17070,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+19842,int_stack+19692, 0.0, zero_stack, 1.0, int_stack+27197, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32600, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18150,int_stack+20067,int_stack+19842, 0.0, zero_stack, 1.0, int_stack+27497, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32900, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+18150,int_stack+13230, 0.0, zero_stack, 1.0, int_stack+83125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97690, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+18150,int_stack+127255,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+20442,int_stack+20382, 0.0, zero_stack, 1.0, int_stack+62080, 0.0, zero_stack, 1.0, int_stack+66404, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+20532,int_stack+20442, 0.0, zero_stack, 1.0, int_stack+62140, 0.0, zero_stack, 1.0, int_stack+66464, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+87340, 0.0, zero_stack, 1.0, int_stack+96850, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+20758,int_stack+20658, 0.0, zero_stack, 1.0, int_stack+62632, 0.0, zero_stack, 1.0, int_stack+66956, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+20908,int_stack+20758, 0.0, zero_stack, 1.0, int_stack+62732, 0.0, zero_stack, 1.0, int_stack+67056, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+127255,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+87880, 0.0, zero_stack, 1.0, int_stack+97390, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+19950,int_stack+13230,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+21268,int_stack+21118, 0.0, zero_stack, 1.0, int_stack+28998, 0.0, zero_stack, 1.0, int_stack+32600, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+21493,int_stack+21268, 0.0, zero_stack, 1.0, int_stack+29298, 0.0, zero_stack, 1.0, int_stack+32900, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10020,int_stack+127255,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+88180, 0.0, zero_stack, 1.0, int_stack+97690, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+144355,int_stack+10020,int_stack+13230,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+21868,int_stack+21808, 0.0, zero_stack, 1.0, int_stack+64242, 1.0, int_stack+66404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+21958,int_stack+21868, 0.0, zero_stack, 1.0, int_stack+64302, 1.0, int_stack+66464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+13230, 0.0, zero_stack, 1.0, int_stack+91795, 1.0, int_stack+96850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+22184,int_stack+22084, 0.0, zero_stack, 1.0, int_stack+64794, 1.0, int_stack+66956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10020,int_stack+22334,int_stack+22184, 0.0, zero_stack, 1.0, int_stack+64894, 1.0, int_stack+67056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10470,int_stack+10020,int_stack+13230, 0.0, zero_stack, 1.0, int_stack+92335, 1.0, int_stack+97390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+21030,int_stack+10470,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10020,int_stack+22694,int_stack+22544, 0.0, zero_stack, 1.0, int_stack+30799, 1.0, int_stack+32600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+22919,int_stack+22694, 0.0, zero_stack, 1.0, int_stack+31099, 1.0, int_stack+32900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22110,int_stack+127255,int_stack+10020, 0.0, zero_stack, 1.0, int_stack+92635, 1.0, int_stack+97690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+146155,int_stack+22110,int_stack+10470,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22110,int_stack+23294,int_stack+23234, 0.0, zero_stack, 2.0, int_stack+66404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+23384,int_stack+23294, 0.0, zero_stack, 2.0, int_stack+66464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+22110, 0.0, zero_stack, 2.0, int_stack+96850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22110,int_stack+23610,int_stack+23510, 0.0, zero_stack, 2.0, int_stack+66956, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22410,int_stack+23760,int_stack+23610, 0.0, zero_stack, 2.0, int_stack+67056, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+22410,int_stack+22110, 0.0, zero_stack, 2.0, int_stack+97390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+22110,int_stack+13230,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23190,int_stack+24120,int_stack+23970, 0.0, zero_stack, 2.0, int_stack+32600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10020,int_stack+24345,int_stack+24120, 0.0, zero_stack, 2.0, int_stack+32900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+10020,int_stack+23190, 0.0, zero_stack, 2.0, int_stack+97690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+147955,int_stack+127255,int_stack+13230,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+24720,int_stack+24660, 1.0, int_stack+57756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68566,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78475,int_stack+24810,int_stack+24720, 1.0, int_stack+57816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68626,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+78475,int_stack+13230, 1.0, int_stack+78295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101305,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+25036,int_stack+24936, 1.0, int_stack+58308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69218,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13230,int_stack+25186,int_stack+25036, 1.0, int_stack+58408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69468,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+13230,int_stack+78295, 1.0, int_stack+79405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101845,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+23190,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+25921,int_stack+25546, 1.0, int_stack+25396, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34551,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+26146,int_stack+25921, 1.0, int_stack+25696, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34851,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24270,int_stack+127255,int_stack+78295, 1.0, int_stack+81835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102145,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+149755,int_stack+24270,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+26521,int_stack+26461, 1.0, int_stack+59918, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68566, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108340,int_stack+26611,int_stack+26521, 1.0, int_stack+59978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68626, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108340,int_stack+108160, 1.0, int_stack+82285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101305, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+26837,int_stack+26737, 1.0, int_stack+60470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69218, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+26987,int_stack+26837, 1.0, int_stack+60570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69468, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 1.0, int_stack+82825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101845, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+24270,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+27722,int_stack+27347, 1.0, int_stack+27197, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34551, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25350,int_stack+27947,int_stack+27722, 1.0, int_stack+27497, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34851, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+25350,int_stack+78295, 1.0, int_stack+83125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102145, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+25350,int_stack+127255,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+28322,int_stack+28262, 1.0, int_stack+62080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68566, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108340,int_stack+28412,int_stack+28322, 1.0, int_stack+62140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+68626, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108340,int_stack+108160, 1.0, int_stack+87340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101305, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+28638,int_stack+28538, 1.0, int_stack+62632, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69218, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+28788,int_stack+28638, 1.0, int_stack+62732, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69468, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 1.0, int_stack+87880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101845, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+27150,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+29523,int_stack+29148, 1.0, int_stack+28998, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34551, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+29748,int_stack+29523, 1.0, int_stack+29298, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34851, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+28230,int_stack+127255,int_stack+78295, 1.0, int_stack+88180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102145, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+151555,int_stack+28230,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+87340,int_stack+30123,int_stack+30063, 1.0, int_stack+64242, 0.0, zero_stack, 1.0, int_stack+68566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+30213,int_stack+30123, 1.0, int_stack+64302, 0.0, zero_stack, 1.0, int_stack+68626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+87340, 1.0, int_stack+91795, 0.0, zero_stack, 1.0, int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+30439,int_stack+30339, 1.0, int_stack+64794, 0.0, zero_stack, 1.0, int_stack+69218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+30589,int_stack+30439, 1.0, int_stack+64894, 0.0, zero_stack, 1.0, int_stack+69468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 1.0, int_stack+92335, 0.0, zero_stack, 1.0, int_stack+101845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+28230,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+31324,int_stack+30949, 1.0, int_stack+30799, 0.0, zero_stack, 1.0, int_stack+34551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29310,int_stack+31549,int_stack+31324, 1.0, int_stack+31099, 0.0, zero_stack, 1.0, int_stack+34851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+29310,int_stack+78295, 1.0, int_stack+92635, 0.0, zero_stack, 1.0, int_stack+102145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+29310,int_stack+127255,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+91795,int_stack+31924,int_stack+31864, 1.0, int_stack+66404, 1.0, int_stack+68566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+32014,int_stack+31924, 1.0, int_stack+66464, 1.0, int_stack+68626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+91795, 1.0, int_stack+96850, 1.0, int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+32240,int_stack+32140, 1.0, int_stack+66956, 1.0, int_stack+69218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+32390,int_stack+32240, 1.0, int_stack+67056, 1.0, int_stack+69468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 1.0, int_stack+97390, 1.0, int_stack+101845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+31110,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+33125,int_stack+32750, 1.0, int_stack+32600, 1.0, int_stack+34551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+33350,int_stack+33125, 1.0, int_stack+32900, 1.0, int_stack+34851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32190,int_stack+127255,int_stack+78295, 1.0, int_stack+97690, 1.0, int_stack+102145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+153355,int_stack+32190,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+96850,int_stack+33725,int_stack+33665, 2.0, int_stack+68566, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+33815,int_stack+33725, 2.0, int_stack+68626, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+96850, 2.0, int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+34041,int_stack+33941, 2.0, int_stack+69218, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+34191,int_stack+34041, 2.0, int_stack+69468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 2.0, int_stack+101845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+32190,int_stack+108160,int_stack+1200,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+35076,int_stack+34701, 2.0, int_stack+34551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33270,int_stack+35301,int_stack+35076, 2.0, int_stack+34851, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+33270,int_stack+78295, 2.0, int_stack+102145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+33270,int_stack+127255,int_stack+108160,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+35676,int_stack+35616, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71353,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+35766,int_stack+35676, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71413,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78115,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+35992,int_stack+35892, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71905,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+36142,int_stack+35992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72005,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79105,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+35070,int_stack+108160,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78745, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+36502,int_stack+36352, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43482,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+36727,int_stack+36502, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43782,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+127255,int_stack+78295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+155155,int_stack+101845,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+37102,int_stack+37042, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71353, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+37192,int_stack+37102, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71413, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78115, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+37418,int_stack+37318, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71905, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+37568,int_stack+37418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72005, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79105, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+36150,int_stack+108160,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+37928,int_stack+37778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43482, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+38153,int_stack+37928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43782, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+101845,int_stack+78295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+156955,int_stack+127255,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+38528,int_stack+38468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71353, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+38618,int_stack+38528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71413, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78115, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+38844,int_stack+38744, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71905, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+38994,int_stack+38844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72005, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79105, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+37230,int_stack+108160,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+39354,int_stack+39204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43482, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+39579,int_stack+39354, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43782, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+127255,int_stack+78295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+158755,int_stack+101845,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+39954,int_stack+39894, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71353, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+40044,int_stack+39954, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71413, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+40270,int_stack+40170, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+40420,int_stack+40270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+38310,int_stack+108160,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+40780,int_stack+40630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+41005,int_stack+40780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+101845,int_stack+78295, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+39390,int_stack+127255,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+41380,int_stack+41320, 0.0, zero_stack, 1.0, int_stack+71353, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+41470,int_stack+41380, 0.0, zero_stack, 1.0, int_stack+71413, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+101305, 0.0, zero_stack, 1.0, int_stack+78115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+41696,int_stack+41596, 0.0, zero_stack, 1.0, int_stack+71905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78295,int_stack+41846,int_stack+41696, 0.0, zero_stack, 1.0, int_stack+72005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78295,int_stack+79405, 0.0, zero_stack, 1.0, int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+160555,int_stack+108160,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+78295,int_stack+42206,int_stack+42056, 0.0, zero_stack, 1.0, int_stack+43482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+42431,int_stack+42206, 0.0, zero_stack, 1.0, int_stack+43782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+127255,int_stack+78295, 0.0, zero_stack, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+161635,int_stack+101845,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+42806,int_stack+42746, 1.0, int_stack+71353, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+42896,int_stack+42806, 1.0, int_stack+71413, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+108160,int_stack+101305, 1.0, int_stack+78115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79405,int_stack+43122,int_stack+43022, 1.0, int_stack+71905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78115,int_stack+43272,int_stack+43122, 1.0, int_stack+72005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+78115,int_stack+79405, 1.0, int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+41190,int_stack+108160,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79105,int_stack+44007,int_stack+43632, 1.0, int_stack+43482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+44232,int_stack+44007, 1.0, int_stack+43782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+101845,int_stack+79105, 1.0, int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+42270,int_stack+127255,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+44607,int_stack+44547,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+44697,int_stack+44607,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+108160,int_stack+101305,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+44923,int_stack+44823,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+127615,int_stack+45073,int_stack+44923,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+79105,int_stack+127615,int_stack+108160,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+44070,int_stack+79105,int_stack+127255, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+79705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+127255,int_stack+45433,int_stack+45283,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+45658,int_stack+45433,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+97390,int_stack+101845,int_stack+127255,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+163435,int_stack+97390,int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+109840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+46033,int_stack+45973, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73515,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79105,int_stack+46123,int_stack+46033, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73575,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97390,int_stack+79105,int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111520,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79105,int_stack+46349,int_stack+46249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74067,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97750,int_stack+46499,int_stack+46349, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74167,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+97750,int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112060,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+45150,int_stack+108160,int_stack+97390, 0.0, zero_stack, 1.0, int_stack+78745, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97390,int_stack+46859,int_stack+46709, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53839,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+47084,int_stack+46859, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54139,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+127255,int_stack+97390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112360,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+165235,int_stack+101845,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+80155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+47459,int_stack+47399, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73515, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+47549,int_stack+47459, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73575, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+108160,int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111520, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+47775,int_stack+47675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74067, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102205,int_stack+47925,int_stack+47775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74167, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+79105,int_stack+102205,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112060, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+46230,int_stack+79105,int_stack+101845, 0.0, zero_stack, 1.0, int_stack+82465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101845,int_stack+48285,int_stack+48135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53839, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97390,int_stack+48510,int_stack+48285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54139, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+97390,int_stack+101845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112360, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+167035,int_stack+127255,int_stack+79105, 0.0, zero_stack, 1.0, int_stack+85660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+48885,int_stack+48825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73515, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79105,int_stack+48975,int_stack+48885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73575, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+79105,int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111520, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79105,int_stack+49201,int_stack+49101, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74067, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127615,int_stack+49351,int_stack+49201, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74167, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+127615,int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112060, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+47310,int_stack+108160,int_stack+127255, 0.0, zero_stack, 1.0, int_stack+87520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+127255,int_stack+49711,int_stack+49561, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53839, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+49936,int_stack+49711, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54139, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97390,int_stack+101845,int_stack+127255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112360, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+48390,int_stack+97390,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+50311,int_stack+50251, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+50401,int_stack+50311, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97390,int_stack+108160,int_stack+101305, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+50627,int_stack+50527, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74067, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97750,int_stack+50777,int_stack+50627, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74167, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+79105,int_stack+97750,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+168835,int_stack+79105,int_stack+97390, 0.0, zero_stack, 1.0, int_stack+91975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97390,int_stack+51137,int_stack+50987, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53839, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+51362,int_stack+51137, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+127255,int_stack+97390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+169915,int_stack+101845,int_stack+79105, 0.0, zero_stack, 1.0, int_stack+95170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+51737,int_stack+51677, 0.0, zero_stack, 1.0, int_stack+73515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79105,int_stack+51827,int_stack+51737, 0.0, zero_stack, 1.0, int_stack+73575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+79105,int_stack+101305, 0.0, zero_stack, 1.0, int_stack+111520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79105,int_stack+52053,int_stack+51953, 0.0, zero_stack, 1.0, int_stack+74067, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102205,int_stack+52203,int_stack+52053, 0.0, zero_stack, 1.0, int_stack+74167, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+102205,int_stack+79105, 0.0, zero_stack, 1.0, int_stack+112060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+50190,int_stack+108160,int_stack+101845, 0.0, zero_stack, 1.0, int_stack+97030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101845,int_stack+52563,int_stack+52413, 0.0, zero_stack, 1.0, int_stack+53839, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97390,int_stack+52788,int_stack+52563, 0.0, zero_stack, 1.0, int_stack+54139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+97390,int_stack+101845, 0.0, zero_stack, 1.0, int_stack+112360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+51270,int_stack+127255,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+101305,int_stack+53163,int_stack+53103, 1.0, int_stack+73515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+53253,int_stack+53163, 1.0, int_stack+73575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+108160,int_stack+101305, 1.0, int_stack+111520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+53479,int_stack+53379, 1.0, int_stack+74067, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127615,int_stack+53629,int_stack+53479, 1.0, int_stack+74167, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+79105,int_stack+127615,int_stack+108160, 1.0, int_stack+112060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+171715,int_stack+79105,int_stack+127255, 0.0, zero_stack, 1.0, int_stack+101485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+127255,int_stack+54364,int_stack+53989, 1.0, int_stack+53839, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+54589,int_stack+54364, 1.0, int_stack+54139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97390,int_stack+101845,int_stack+127255, 1.0, int_stack+112360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+53070,int_stack+97390,int_stack+79105, 0.0, zero_stack, 1.0, int_stack+104680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+54964,int_stack+54904,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+79105,int_stack+55054,int_stack+54964,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+97390,int_stack+79105,int_stack+111520,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+79105,int_stack+55280,int_stack+55180,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+97750,int_stack+55430,int_stack+55280,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+97750,int_stack+79105,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+172795,int_stack+108160,int_stack+97390, 0.0, zero_stack, 1.0, int_stack+79705, 1.0, int_stack+111700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+97390,int_stack+55790,int_stack+55640,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+56015,int_stack+55790,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+127255,int_stack+97390,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+173875,int_stack+101845,int_stack+108160, 0.0, zero_stack, 1.0, int_stack+109840, 1.0, int_stack+114895, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+56390,int_stack+56330,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+56480,int_stack+56390,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+108160,int_stack+111520,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+56706,int_stack+56606,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+102205,int_stack+56856,int_stack+56706,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+79105,int_stack+102205,int_stack+108160,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+54870,int_stack+79105,int_stack+101845, 0.0, zero_stack, 2.0, int_stack+111700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+101845,int_stack+57216,int_stack+57066,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+97390,int_stack+57441,int_stack+57216,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+97390,int_stack+101845,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+55950,int_stack+127255,int_stack+79105, 0.0, zero_stack, 2.0, int_stack+114895, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+58092,int_stack+58032, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75953,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+79105,int_stack+58182,int_stack+58092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76013,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+79105,int_stack+111520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2010,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+79105,int_stack+58868,int_stack+58768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76505,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127615,int_stack+59018,int_stack+58868, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76605,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+127615,int_stack+79105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2550,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+57750,int_stack+108160,int_stack+127255, 1.0, int_stack+78745, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+127255,int_stack+59378,int_stack+59228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70288,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+59603,int_stack+59378, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70588,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97390,int_stack+101845,int_stack+127255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118375,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+175675,int_stack+97390,int_stack+108160, 1.0, int_stack+80155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+60254,int_stack+60194, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75953, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+60344,int_stack+60254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76013, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+97390,int_stack+108160,int_stack+111520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2010, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+61030,int_stack+60930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76505, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+97750,int_stack+61180,int_stack+61030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76605, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+97750,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2550, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+78115,int_stack+13230,int_stack+97390, 1.0, int_stack+82465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+97390,int_stack+61540,int_stack+61390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70288, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+61765,int_stack+61540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70588, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+127255,int_stack+97390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118375, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+58830,int_stack+101845,int_stack+13230, 1.0, int_stack+85660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+62416,int_stack+62356, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75953, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+85660,int_stack+62506,int_stack+62416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76013, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+85660,int_stack+111520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2010, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+85660,int_stack+63192,int_stack+63092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76505, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+63342,int_stack+63192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76605, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+101845,int_stack+85660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2550, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+60630,int_stack+108160,int_stack+13230, 1.0, int_stack+87520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+63702,int_stack+63552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70288, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+63927,int_stack+63702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70588, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+101845,int_stack+13230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118375, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+61710,int_stack+127255,int_stack+108160, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+64578,int_stack+64518, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+75953, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+64668,int_stack+64578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76013, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+0,int_stack+111520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+65354,int_stack+65254, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+65504,int_stack+65354, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+127255,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+63510,int_stack+13230,int_stack+108160, 1.0, int_stack+91975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+65864,int_stack+65714, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+66089,int_stack+65864, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+101845,int_stack+127255,int_stack+108160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+64590,int_stack+101845,int_stack+13230, 1.0, int_stack+95170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+66740,int_stack+66680, 0.0, zero_stack, 1.0, int_stack+75953, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+95170,int_stack+66830,int_stack+66740, 0.0, zero_stack, 1.0, int_stack+76013, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+95170,int_stack+111520, 0.0, zero_stack, 1.0, int_stack+2010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+95170,int_stack+67516,int_stack+67416, 0.0, zero_stack, 1.0, int_stack+76505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+67666,int_stack+67516, 0.0, zero_stack, 1.0, int_stack+76605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108160,int_stack+101845,int_stack+95170, 0.0, zero_stack, 1.0, int_stack+2550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+66390,int_stack+108160,int_stack+13230, 1.0, int_stack+97030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+68026,int_stack+67876, 0.0, zero_stack, 1.0, int_stack+70288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+101845,int_stack+68251,int_stack+68026, 0.0, zero_stack, 1.0, int_stack+70588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+101845,int_stack+13230, 0.0, zero_stack, 1.0, int_stack+118375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+81835,int_stack+127255,int_stack+108160, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+68902,int_stack+68842, 1.0, int_stack+75953, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+108160,int_stack+68992,int_stack+68902, 1.0, int_stack+76013, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+108160,int_stack+111520, 1.0, int_stack+2010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108160,int_stack+69928,int_stack+69828, 1.0, int_stack+76505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+127615,int_stack+70078,int_stack+69928, 1.0, int_stack+76605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+127615,int_stack+108160, 1.0, int_stack+2550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+13230,int_stack+127255, 1.0, int_stack+101485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+127255,int_stack+70813,int_stack+70438, 1.0, int_stack+70288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1080,int_stack+71038,int_stack+70813, 1.0, int_stack+70588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+112060,int_stack+1080,int_stack+127255, 1.0, int_stack+118375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+67470,int_stack+112060,int_stack+13230, 1.0, int_stack+104680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+71689,int_stack+71629,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+104680,int_stack+71779,int_stack+71689,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+104680,int_stack+111520,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+72465,int_stack+72365,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+118375,int_stack+72615,int_stack+72465,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+104680,int_stack+118375,int_stack+2550,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1080,int_stack+104680,int_stack+13230, 1.0, int_stack+79705, 0.0, zero_stack, 1.0, int_stack+2190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+118375,int_stack+72975,int_stack+72825,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+112060,int_stack+73200,int_stack+72975,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+127255,int_stack+112060,int_stack+118375,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+69270,int_stack+127255,int_stack+104680, 1.0, int_stack+109840, 0.0, zero_stack, 1.0, int_stack+118825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+111520,int_stack+73851,int_stack+73791,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+109840,int_stack+73941,int_stack+73851,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+104680,int_stack+109840,int_stack+111520,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+74627,int_stack+74527,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+118375,int_stack+74777,int_stack+74627,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+109840,int_stack+118375,int_stack+2550,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+101305,int_stack+109840,int_stack+104680, 1.0, int_stack+111700, 1.0, int_stack+2190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+118375,int_stack+75137,int_stack+74987,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+111520,int_stack+75362,int_stack+75137,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+112195,int_stack+111520,int_stack+118375,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+71070,int_stack+112195,int_stack+109840, 1.0, int_stack+114895, 1.0, int_stack+118825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+114895,int_stack+76289,int_stack+76229,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+115075,int_stack+76379,int_stack+76289,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+109840,int_stack+115075,int_stack+114895,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+77065,int_stack+76965,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+118375,int_stack+77215,int_stack+77065,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+114895,int_stack+118375,int_stack+2550,10);
 /*--- compute (dp|fd) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+111520,int_stack+114895,int_stack+109840, 2.0, int_stack+2190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+118375,int_stack+77575,int_stack+77425,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+127255,int_stack+77800,int_stack+77575,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+79195,int_stack+127255,int_stack+118375,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+72870,int_stack+79195,int_stack+114895, 2.0, int_stack+118825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+74670,int_stack+83860,int_stack+80755,60);
     Libderiv->ABCD[11] = int_stack + 74670;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+83635,int_stack+88915,int_stack+86260,60);
     Libderiv->ABCD[10] = int_stack + 83635;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+87340,int_stack+93370,int_stack+90715,60);
     Libderiv->ABCD[9] = int_stack + 87340;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+91795,int_stack+98425,int_stack+95770,60);
     Libderiv->ABCD[8] = int_stack + 91795;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+96850,int_stack+102880,int_stack+100225,60);
     Libderiv->ABCD[7] = int_stack + 96850;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+102385,int_stack+106360,int_stack+105280,60);
     Libderiv->ABCD[6] = int_stack + 102385;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+106360,int_stack+113095,int_stack+110440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 106360;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+112600,int_stack+116575,int_stack+115495, 0.0, zero_stack, 1.0, int_stack+108760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 112600;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+116575,int_stack+122080,int_stack+119425, 1.0, int_stack+108760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 116575;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+121585,int_stack+125455,int_stack+2850,60);
     Libderiv->ABCD[155] = int_stack + 121585;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+124960,int_stack+3930,int_stack+120505,60);
     Libderiv->ABCD[143] = int_stack + 124960;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+2160,int_stack+128155,int_stack+5730,60);
     Libderiv->ABCD[142] = int_stack + 2160;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+4320,int_stack+129955,int_stack+123880,60);
     Libderiv->ABCD[131] = int_stack + 4320;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+127120,int_stack+131755,int_stack+7860,60);
     Libderiv->ABCD[130] = int_stack + 127120;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+6480,int_stack+133555,int_stack+8940,60);
     Libderiv->ABCD[129] = int_stack + 6480;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+8640,int_stack+135355,int_stack+11070,60);
     Libderiv->ABCD[119] = int_stack + 8640;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+129280,int_stack+137155,int_stack+12150,60);
     Libderiv->ABCD[118] = int_stack + 129280;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+10800,int_stack+138955,int_stack+13830,60);
     Libderiv->ABCD[117] = int_stack + 10800;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+131440,int_stack+140755,int_stack+14910,60);
     Libderiv->ABCD[116] = int_stack + 131440;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+12960,int_stack+142555,int_stack+15990,60);
     Libderiv->ABCD[107] = int_stack + 12960;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+133600,int_stack+18150,int_stack+17070,60);
     Libderiv->ABCD[106] = int_stack + 133600;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+15120,int_stack+144355,int_stack+19950,60);
     Libderiv->ABCD[105] = int_stack + 15120;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+17280,int_stack+146155,int_stack+21030,60);
     Libderiv->ABCD[104] = int_stack + 17280;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+19440,int_stack+147955,int_stack+22110,60);
     Libderiv->ABCD[103] = int_stack + 19440;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+135760,int_stack+149755,int_stack+23190,60);
     Libderiv->ABCD[95] = int_stack + 135760;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+21600,int_stack+25350,int_stack+24270,60);
     Libderiv->ABCD[94] = int_stack + 21600;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+23760,int_stack+151555,int_stack+27150,60);
     Libderiv->ABCD[93] = int_stack + 23760;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+25920,int_stack+29310,int_stack+28230,60);
     Libderiv->ABCD[92] = int_stack + 25920;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+28080,int_stack+153355,int_stack+31110,60);
     Libderiv->ABCD[91] = int_stack + 28080;
 /*--- compute (dd|fd) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+137920,int_stack+33270,int_stack+32190,60);
     Libderiv->ABCD[90] = int_stack + 137920;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+30240,int_stack+155155,int_stack+35070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[47] = int_stack + 30240;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+32400,int_stack+156955,int_stack+36150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[46] = int_stack + 32400;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+34560,int_stack+158755,int_stack+37230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[45] = int_stack + 34560;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+140080,int_stack+39390,int_stack+38310, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[44] = int_stack + 140080;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+36720,int_stack+161635,int_stack+160555, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[43] = int_stack + 36720;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+38880,int_stack+42270,int_stack+41190, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[42] = int_stack + 38880;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+41040,int_stack+163435,int_stack+44070, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+110440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[38] = int_stack + 41040;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+142240,int_stack+165235,int_stack+45150, 0.0, zero_stack, 1.0, int_stack+80755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[35] = int_stack + 142240;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+43200,int_stack+167035,int_stack+46230, 0.0, zero_stack, 1.0, int_stack+86260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[34] = int_stack + 43200;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+144400,int_stack+48390,int_stack+47310, 0.0, zero_stack, 1.0, int_stack+90715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[33] = int_stack + 144400;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+45360,int_stack+169915,int_stack+168835, 0.0, zero_stack, 1.0, int_stack+95770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[32] = int_stack + 45360;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+47520,int_stack+51270,int_stack+50190, 0.0, zero_stack, 1.0, int_stack+100225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[31] = int_stack + 47520;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+49680,int_stack+53070,int_stack+171715, 0.0, zero_stack, 1.0, int_stack+105280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[30] = int_stack + 49680;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+51840,int_stack+173875,int_stack+172795, 0.0, zero_stack, 1.0, int_stack+110440, 1.0, int_stack+115495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[26] = int_stack + 51840;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+146560,int_stack+55950,int_stack+54870, 0.0, zero_stack, 2.0, int_stack+115495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[25] = int_stack + 146560;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+54000,int_stack+175675,int_stack+57750, 1.0, int_stack+80755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[23] = int_stack + 54000;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+79195,int_stack+58830,int_stack+78115, 1.0, int_stack+86260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[22] = int_stack + 79195;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+56160,int_stack+61710,int_stack+60630, 1.0, int_stack+90715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[21] = int_stack + 56160;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+58320,int_stack+64590,int_stack+63510, 1.0, int_stack+95770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[20] = int_stack + 58320;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+60480,int_stack+81835,int_stack+66390, 1.0, int_stack+100225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[19] = int_stack + 60480;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+81355,int_stack+67470,int_stack+0, 1.0, int_stack+105280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[18] = int_stack + 81355;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+62640,int_stack+69270,int_stack+1080, 1.0, int_stack+110440, 0.0, zero_stack, 1.0, int_stack+119425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[14] = int_stack + 62640;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+71070,int_stack+101305, 1.0, int_stack+115495, 1.0, int_stack+119425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[13] = int_stack + 0;
 /*--- compute (dd|fd) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+108520,int_stack+72870,int_stack+111520, 2.0, int_stack+119425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[12] = int_stack + 108520;

}
