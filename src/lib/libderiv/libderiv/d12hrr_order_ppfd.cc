#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ppfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|fd) integrals */

void d12hrr_order_ppfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][5][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][10] = int_stack + 126;
 Libderiv->deriv_classes[2][5][9] = int_stack + 252;
 Libderiv->deriv_classes[2][5][8] = int_stack + 378;
 Libderiv->deriv_classes[2][5][7] = int_stack + 504;
 Libderiv->dvrr_classes[2][4] = int_stack + 630;
 Libderiv->deriv_classes[2][5][6] = int_stack + 720;
 Libderiv->deriv_classes[2][5][2] = int_stack + 846;
 Libderiv->deriv_classes[2][5][1] = int_stack + 972;
 Libderiv->dvrr_classes[1][5] = int_stack + 1098;
 Libderiv->deriv_classes[2][5][0] = int_stack + 1161;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 1287;
 Libderiv->deriv2_classes[1][4][143] = int_stack + 1317;
 Libderiv->deriv2_classes[1][5][143] = int_stack + 1362;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 1425;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 1485;
 Libderiv->deriv2_classes[2][5][143] = int_stack + 1575;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 1701;
 Libderiv->deriv2_classes[1][4][131] = int_stack + 1731;
 Libderiv->deriv2_classes[1][5][131] = int_stack + 1776;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1839;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 1899;
 Libderiv->deriv2_classes[2][5][131] = int_stack + 1989;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 2115;
 Libderiv->deriv2_classes[1][4][130] = int_stack + 2145;
 Libderiv->deriv2_classes[1][5][130] = int_stack + 2190;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 2253;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 2313;
 Libderiv->deriv2_classes[2][5][130] = int_stack + 2403;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 2529;
 Libderiv->deriv2_classes[1][4][119] = int_stack + 2559;
 Libderiv->deriv2_classes[1][5][119] = int_stack + 2604;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 2667;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 2727;
 Libderiv->deriv2_classes[2][5][119] = int_stack + 2817;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 2943;
 Libderiv->deriv2_classes[1][4][118] = int_stack + 2973;
 Libderiv->deriv2_classes[1][5][118] = int_stack + 3018;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 3081;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 3141;
 Libderiv->deriv2_classes[2][5][118] = int_stack + 3231;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 3357;
 Libderiv->deriv2_classes[1][4][117] = int_stack + 3387;
 Libderiv->deriv2_classes[1][5][117] = int_stack + 3432;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 3495;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 3555;
 Libderiv->deriv2_classes[2][5][117] = int_stack + 3645;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 3771;
 Libderiv->deriv2_classes[1][4][107] = int_stack + 3801;
 Libderiv->deriv2_classes[1][5][107] = int_stack + 3846;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 3909;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 3969;
 Libderiv->deriv2_classes[2][5][107] = int_stack + 4059;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 4185;
 Libderiv->deriv2_classes[1][4][106] = int_stack + 4215;
 Libderiv->deriv2_classes[1][5][106] = int_stack + 4260;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 4323;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 4383;
 Libderiv->deriv2_classes[2][5][106] = int_stack + 4473;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 4599;
 Libderiv->deriv2_classes[1][4][105] = int_stack + 4629;
 Libderiv->deriv2_classes[1][5][105] = int_stack + 4674;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 4737;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 4797;
 Libderiv->deriv2_classes[2][5][105] = int_stack + 4887;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 5013;
 Libderiv->deriv2_classes[1][4][104] = int_stack + 5043;
 Libderiv->deriv2_classes[1][5][104] = int_stack + 5088;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 5151;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 5211;
 Libderiv->deriv2_classes[2][5][104] = int_stack + 5301;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 5427;
 Libderiv->deriv2_classes[1][4][95] = int_stack + 5457;
 Libderiv->deriv2_classes[1][5][95] = int_stack + 5502;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 5565;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 5625;
 Libderiv->deriv2_classes[2][5][95] = int_stack + 5715;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 5841;
 Libderiv->deriv2_classes[1][4][94] = int_stack + 5871;
 Libderiv->deriv2_classes[1][5][94] = int_stack + 5916;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 5979;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 6039;
 Libderiv->deriv2_classes[2][5][94] = int_stack + 6129;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 6255;
 Libderiv->deriv2_classes[1][4][93] = int_stack + 6285;
 Libderiv->deriv2_classes[1][5][93] = int_stack + 6330;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 6393;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 6453;
 Libderiv->deriv2_classes[2][5][93] = int_stack + 6543;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 6669;
 Libderiv->deriv2_classes[1][4][92] = int_stack + 6699;
 Libderiv->deriv2_classes[1][5][92] = int_stack + 6744;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 6807;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 6867;
 Libderiv->deriv2_classes[2][5][92] = int_stack + 6957;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 7083;
 Libderiv->deriv2_classes[1][4][91] = int_stack + 7113;
 Libderiv->deriv2_classes[1][5][91] = int_stack + 7158;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 7221;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 7281;
 Libderiv->deriv2_classes[2][5][91] = int_stack + 7371;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 7497;
 Libderiv->deriv2_classes[1][4][83] = int_stack + 7527;
 Libderiv->deriv2_classes[1][5][83] = int_stack + 7572;
 Libderiv->deriv_classes[2][3][11] = int_stack + 7635;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 7695;
 Libderiv->deriv_classes[2][4][11] = int_stack + 7755;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 7845;
 Libderiv->deriv2_classes[2][5][83] = int_stack + 7935;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 8061;
 Libderiv->deriv2_classes[1][4][82] = int_stack + 8091;
 Libderiv->deriv2_classes[1][5][82] = int_stack + 8136;
 Libderiv->deriv_classes[2][3][10] = int_stack + 8199;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 8259;
 Libderiv->deriv_classes[2][4][10] = int_stack + 8319;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 8409;
 Libderiv->deriv2_classes[2][5][82] = int_stack + 8499;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 8625;
 Libderiv->deriv2_classes[1][4][81] = int_stack + 8655;
 Libderiv->deriv2_classes[1][5][81] = int_stack + 8700;
 Libderiv->deriv_classes[2][3][9] = int_stack + 8763;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 8823;
 Libderiv->deriv_classes[2][4][9] = int_stack + 8883;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 8973;
 Libderiv->deriv2_classes[2][5][81] = int_stack + 9063;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 9189;
 Libderiv->deriv2_classes[1][4][80] = int_stack + 9219;
 Libderiv->deriv2_classes[1][5][80] = int_stack + 9264;
 Libderiv->deriv_classes[2][3][8] = int_stack + 9327;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 9387;
 Libderiv->deriv_classes[2][4][8] = int_stack + 9447;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 9537;
 Libderiv->deriv2_classes[2][5][80] = int_stack + 9627;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 9753;
 Libderiv->deriv2_classes[1][4][79] = int_stack + 9783;
 Libderiv->deriv2_classes[1][5][79] = int_stack + 9828;
 Libderiv->deriv_classes[2][3][7] = int_stack + 9891;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 9951;
 Libderiv->deriv_classes[2][4][7] = int_stack + 10011;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 10101;
 Libderiv->deriv2_classes[2][5][79] = int_stack + 10191;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 10317;
 Libderiv->deriv2_classes[1][4][78] = int_stack + 10347;
 Libderiv->deriv2_classes[1][5][78] = int_stack + 10392;
 Libderiv->dvrr_classes[2][3] = int_stack + 10455;
 Libderiv->deriv_classes[2][3][6] = int_stack + 10515;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 10575;
 Libderiv->deriv_classes[2][4][6] = int_stack + 10635;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 10725;
 Libderiv->deriv2_classes[2][5][78] = int_stack + 10815;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 10941;
 Libderiv->deriv2_classes[1][4][35] = int_stack + 10971;
 Libderiv->deriv2_classes[1][5][35] = int_stack + 11016;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 11079;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 11139;
 Libderiv->deriv2_classes[2][5][35] = int_stack + 11229;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 11355;
 Libderiv->deriv2_classes[1][4][34] = int_stack + 11385;
 Libderiv->deriv2_classes[1][5][34] = int_stack + 11430;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 11493;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 11553;
 Libderiv->deriv2_classes[2][5][34] = int_stack + 11643;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 11769;
 Libderiv->deriv2_classes[1][4][33] = int_stack + 11799;
 Libderiv->deriv2_classes[1][5][33] = int_stack + 11844;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 11907;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 11967;
 Libderiv->deriv2_classes[2][5][33] = int_stack + 12057;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 12183;
 Libderiv->deriv2_classes[1][4][32] = int_stack + 12213;
 Libderiv->deriv2_classes[1][5][32] = int_stack + 12258;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 12321;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 12381;
 Libderiv->deriv2_classes[2][5][32] = int_stack + 12471;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 12597;
 Libderiv->deriv2_classes[1][4][31] = int_stack + 12627;
 Libderiv->deriv2_classes[1][5][31] = int_stack + 12672;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 12735;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 12795;
 Libderiv->deriv2_classes[2][5][31] = int_stack + 12885;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 13011;
 Libderiv->deriv2_classes[1][4][30] = int_stack + 13041;
 Libderiv->deriv2_classes[1][5][30] = int_stack + 13086;
 Libderiv->deriv_classes[2][3][2] = int_stack + 13149;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 13209;
 Libderiv->deriv_classes[2][4][2] = int_stack + 13269;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 13359;
 Libderiv->deriv2_classes[2][5][30] = int_stack + 13449;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 13575;
 Libderiv->deriv2_classes[1][4][26] = int_stack + 13605;
 Libderiv->deriv2_classes[1][5][26] = int_stack + 13650;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 13713;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 13773;
 Libderiv->deriv2_classes[2][5][26] = int_stack + 13863;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 13989;
 Libderiv->deriv2_classes[1][4][23] = int_stack + 14019;
 Libderiv->deriv2_classes[1][5][23] = int_stack + 14064;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 14127;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 14187;
 Libderiv->deriv2_classes[2][5][23] = int_stack + 14277;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 14403;
 Libderiv->deriv2_classes[1][4][22] = int_stack + 14433;
 Libderiv->deriv2_classes[1][5][22] = int_stack + 14478;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 14541;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 14601;
 Libderiv->deriv2_classes[2][5][22] = int_stack + 14691;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 14817;
 Libderiv->deriv2_classes[1][4][21] = int_stack + 14847;
 Libderiv->deriv2_classes[1][5][21] = int_stack + 14892;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 14955;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 15015;
 Libderiv->deriv2_classes[2][5][21] = int_stack + 15105;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 15231;
 Libderiv->deriv2_classes[1][4][20] = int_stack + 15261;
 Libderiv->deriv2_classes[1][5][20] = int_stack + 15306;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 15369;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 15429;
 Libderiv->deriv2_classes[2][5][20] = int_stack + 15519;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 15645;
 Libderiv->deriv2_classes[1][4][19] = int_stack + 15675;
 Libderiv->deriv2_classes[1][5][19] = int_stack + 15720;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 15783;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 15843;
 Libderiv->deriv2_classes[2][5][19] = int_stack + 15933;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 16059;
 Libderiv->deriv2_classes[1][4][18] = int_stack + 16089;
 Libderiv->deriv2_classes[1][5][18] = int_stack + 16134;
 Libderiv->deriv_classes[2][3][1] = int_stack + 16197;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 16257;
 Libderiv->deriv_classes[2][4][1] = int_stack + 16317;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 16407;
 Libderiv->deriv2_classes[2][5][18] = int_stack + 16497;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 16623;
 Libderiv->deriv2_classes[1][4][14] = int_stack + 16653;
 Libderiv->deriv2_classes[1][5][14] = int_stack + 16698;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 16761;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 16821;
 Libderiv->deriv2_classes[2][5][14] = int_stack + 16911;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 17037;
 Libderiv->deriv2_classes[1][4][13] = int_stack + 17067;
 Libderiv->deriv2_classes[1][5][13] = int_stack + 17112;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 17175;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 17235;
 Libderiv->deriv2_classes[2][5][13] = int_stack + 17325;
 Libderiv->deriv_classes[1][3][11] = int_stack + 17451;
 Libderiv->deriv_classes[1][4][11] = int_stack + 17481;
 Libderiv->deriv_classes[1][5][11] = int_stack + 17526;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 17589;
 Libderiv->deriv2_classes[1][4][11] = int_stack + 17619;
 Libderiv->deriv2_classes[1][5][11] = int_stack + 17664;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 17727;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 17787;
 Libderiv->deriv2_classes[2][5][11] = int_stack + 17877;
 Libderiv->deriv_classes[1][3][10] = int_stack + 18003;
 Libderiv->deriv_classes[1][4][10] = int_stack + 18033;
 Libderiv->deriv_classes[1][5][10] = int_stack + 18078;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 18141;
 Libderiv->deriv2_classes[1][4][10] = int_stack + 18171;
 Libderiv->deriv2_classes[1][5][10] = int_stack + 18216;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 18279;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 18339;
 Libderiv->deriv2_classes[2][5][10] = int_stack + 18429;
 Libderiv->deriv_classes[1][3][9] = int_stack + 18555;
 Libderiv->deriv_classes[1][4][9] = int_stack + 18585;
 Libderiv->deriv_classes[1][5][9] = int_stack + 18630;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 18693;
 Libderiv->deriv2_classes[1][4][9] = int_stack + 18723;
 Libderiv->deriv2_classes[1][5][9] = int_stack + 18768;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 18831;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 18891;
 Libderiv->deriv2_classes[2][5][9] = int_stack + 18981;
 Libderiv->deriv_classes[1][3][8] = int_stack + 19107;
 Libderiv->deriv_classes[1][4][8] = int_stack + 19137;
 Libderiv->deriv_classes[1][5][8] = int_stack + 19182;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 19245;
 Libderiv->deriv2_classes[1][4][8] = int_stack + 19275;
 Libderiv->deriv2_classes[1][5][8] = int_stack + 19320;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 19383;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 19443;
 Libderiv->deriv2_classes[2][5][8] = int_stack + 19533;
 Libderiv->deriv_classes[1][3][7] = int_stack + 19659;
 Libderiv->deriv_classes[1][4][7] = int_stack + 19689;
 Libderiv->deriv_classes[1][5][7] = int_stack + 19734;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 19797;
 Libderiv->deriv2_classes[1][4][7] = int_stack + 19827;
 Libderiv->deriv2_classes[1][5][7] = int_stack + 19872;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 19935;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 19995;
 Libderiv->deriv2_classes[2][5][7] = int_stack + 20085;
 Libderiv->dvrr_classes[1][3] = int_stack + 20211;
 Libderiv->deriv_classes[1][3][6] = int_stack + 20241;
 Libderiv->dvrr_classes[1][4] = int_stack + 20271;
 Libderiv->deriv_classes[1][4][6] = int_stack + 20316;
 Libderiv->deriv_classes[1][5][6] = int_stack + 20361;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 20424;
 Libderiv->deriv2_classes[1][4][6] = int_stack + 20454;
 Libderiv->deriv2_classes[1][5][6] = int_stack + 20499;
 Libderiv->deriv_classes[2][3][0] = int_stack + 20562;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 20622;
 Libderiv->deriv_classes[2][4][0] = int_stack + 20682;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 20772;
 Libderiv->deriv2_classes[2][5][6] = int_stack + 20862;
 Libderiv->deriv_classes[1][3][2] = int_stack + 20988;
 Libderiv->deriv_classes[1][4][2] = int_stack + 21018;
 Libderiv->deriv_classes[1][5][2] = int_stack + 21063;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 21126;
 Libderiv->deriv2_classes[1][4][2] = int_stack + 21156;
 Libderiv->deriv2_classes[1][5][2] = int_stack + 21201;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 21264;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 21324;
 Libderiv->deriv2_classes[2][5][2] = int_stack + 21414;
 Libderiv->deriv_classes[1][3][1] = int_stack + 21540;
 Libderiv->deriv_classes[1][4][1] = int_stack + 21570;
 Libderiv->deriv_classes[1][5][1] = int_stack + 21615;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 21678;
 Libderiv->deriv2_classes[1][4][1] = int_stack + 21708;
 Libderiv->deriv2_classes[1][5][1] = int_stack + 21753;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 21816;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 21876;
 Libderiv->deriv2_classes[2][5][1] = int_stack + 21966;
 Libderiv->deriv_classes[1][3][0] = int_stack + 22092;
 Libderiv->deriv_classes[1][4][0] = int_stack + 22122;
 Libderiv->deriv_classes[1][5][0] = int_stack + 22167;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 22230;
 Libderiv->deriv2_classes[1][4][0] = int_stack + 22260;
 Libderiv->deriv2_classes[1][5][0] = int_stack + 22305;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 22368;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 22428;
 Libderiv->deriv2_classes[2][5][0] = int_stack + 22518;
 memset(int_stack,0,181152);

 Libderiv->dvrr_stack = int_stack + 38799;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ppfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+20271,int_stack+20211,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22734,int_stack+17481,int_stack+17451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20211,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+17526,int_stack+17481, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20271,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22959,int_stack+22824,int_stack+22734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22644,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+23139,int_stack+630,int_stack+10455,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23319,int_stack+7755,int_stack+7635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10455,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23499,int_stack+0,int_stack+7755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+630,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23769,int_stack+23499,int_stack+23319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23139,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23499,int_stack+18033,int_stack+18003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20211, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+18078,int_stack+18033, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20271, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23589,int_stack+22824,int_stack+23499, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22644, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24129,int_stack+8319,int_stack+8199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10455, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24309,int_stack+126,int_stack+8319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24579,int_stack+24309,int_stack+24129, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23139, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24309,int_stack+18585,int_stack+18555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20211, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+18630,int_stack+18585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20271, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24399,int_stack+22824,int_stack+24309, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22644, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+8883,int_stack+8763, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10455, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24939,int_stack+252,int_stack+8883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25209,int_stack+24939,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23139, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24939,int_stack+19137,int_stack+19107, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+19182,int_stack+19137, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20271, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25029,int_stack+22824,int_stack+24939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+9447,int_stack+9327, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+25569,int_stack+378,int_stack+9447, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25839,int_stack+25569,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25569,int_stack+19689,int_stack+19659, 0.0, zero_stack, 1.0, int_stack+20211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+19734,int_stack+19689, 0.0, zero_stack, 1.0, int_stack+20271, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25659,int_stack+22824,int_stack+25569, 0.0, zero_stack, 1.0, int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26199,int_stack+10011,int_stack+9891, 0.0, zero_stack, 1.0, int_stack+10455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+26379,int_stack+504,int_stack+10011, 0.0, zero_stack, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26649,int_stack+26379,int_stack+26199, 0.0, zero_stack, 1.0, int_stack+23139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26379,int_stack+20316,int_stack+20241, 1.0, int_stack+20211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+20361,int_stack+20316, 1.0, int_stack+20271, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26469,int_stack+22824,int_stack+26379, 1.0, int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+10635,int_stack+10515, 1.0, int_stack+10455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+27009,int_stack+720,int_stack+10635, 1.0, int_stack+630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+27279,int_stack+27009,int_stack+360, 1.0, int_stack+23139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+1098,int_stack+20271,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+23139,int_stack+22824,int_stack+22644,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+21018,int_stack+20988,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+21063,int_stack+21018,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+27009,int_stack+22824,int_stack+22644,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+540,int_stack+13269,int_stack+13149,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+27639,int_stack+846,int_stack+13269,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+27909,int_stack+27639,int_stack+540,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+21570,int_stack+21540,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+21615,int_stack+21570,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+27639,int_stack+22824,int_stack+27189,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+16317,int_stack+16197,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28269,int_stack+972,int_stack+16317,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+28539,int_stack+28269,int_stack+720,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27819,int_stack+22122,int_stack+22092,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+22167,int_stack+22122,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+28269,int_stack+22824,int_stack+27819,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+20682,int_stack+20562,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+28899,int_stack+1161,int_stack+20682,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+29169,int_stack+28899,int_stack+900,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+1317,int_stack+1287, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17451,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+1362,int_stack+1317, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17481,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+28899,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22734,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1080,int_stack+1485,int_stack+1425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7635,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29529,int_stack+1575,int_stack+1485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+7755,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1260,int_stack+29529,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+23319,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+1731,int_stack+1701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17451, 1.0, int_stack+18003,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+1776,int_stack+1731, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17481, 1.0, int_stack+18033,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1080,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22734, 1.0, int_stack+23499,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29529,int_stack+1899,int_stack+1839, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7635, 1.0, int_stack+8199,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29709,int_stack+1989,int_stack+1899, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7755, 1.0, int_stack+8319,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1620,int_stack+29709,int_stack+29529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23319, 1.0, int_stack+24129,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+2145,int_stack+2115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+18003, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+2190,int_stack+2145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+18033, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29529,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+23499, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29709,int_stack+2313,int_stack+2253, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8199, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1980,int_stack+2403,int_stack+2313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8319, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29889,int_stack+1980,int_stack+29709, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24129, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+2559,int_stack+2529, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17451, 0.0, zero_stack, 1.0, int_stack+18555,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+2604,int_stack+2559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17481, 0.0, zero_stack, 1.0, int_stack+18585,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29709,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22734, 0.0, zero_stack, 1.0, int_stack+24309,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1980,int_stack+2727,int_stack+2667, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7635, 0.0, zero_stack, 1.0, int_stack+8763,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2160,int_stack+2817,int_stack+2727, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7755, 0.0, zero_stack, 1.0, int_stack+8883,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2430,int_stack+2160,int_stack+1980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23319, 0.0, zero_stack, 1.0, int_stack+0,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+2973,int_stack+2943, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18003, 1.0, int_stack+18555, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+3018,int_stack+2973, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18033, 1.0, int_stack+18585, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1980,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23499, 1.0, int_stack+24309, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2160,int_stack+3141,int_stack+3081, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8199, 1.0, int_stack+8763, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2790,int_stack+3231,int_stack+3141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8319, 1.0, int_stack+8883, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30249,int_stack+2790,int_stack+2160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24129, 1.0, int_stack+0, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+3387,int_stack+3357, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+18555, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+3432,int_stack+3387, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+18585, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2160,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24309, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2790,int_stack+3555,int_stack+3495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8763, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2970,int_stack+3645,int_stack+3555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8883, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3240,int_stack+2970,int_stack+2790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+3801,int_stack+3771, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17451, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19107,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+3846,int_stack+3801, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17481, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19137,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2790,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22734, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24939,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2970,int_stack+3969,int_stack+3909, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7635, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9327,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3600,int_stack+4059,int_stack+3969, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9447,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30609,int_stack+3600,int_stack+2970, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23319, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+4215,int_stack+4185, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18003, 0.0, zero_stack, 1.0, int_stack+19107, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+4260,int_stack+4215, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18033, 0.0, zero_stack, 1.0, int_stack+19137, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2970,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23499, 0.0, zero_stack, 1.0, int_stack+24939, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4383,int_stack+4323, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8199, 0.0, zero_stack, 1.0, int_stack+9327, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3780,int_stack+4473,int_stack+4383, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8319, 0.0, zero_stack, 1.0, int_stack+9447, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4050,int_stack+3780,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24129, 0.0, zero_stack, 1.0, int_stack+180, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+4629,int_stack+4599, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18555, 1.0, int_stack+19107, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+4674,int_stack+4629, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18585, 1.0, int_stack+19137, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3600,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24309, 1.0, int_stack+24939, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3780,int_stack+4797,int_stack+4737, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8763, 1.0, int_stack+9327, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4410,int_stack+4887,int_stack+4797, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8883, 1.0, int_stack+9447, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30969,int_stack+4410,int_stack+3780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+5043,int_stack+5013, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+19107, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+5088,int_stack+5043, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+19137, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3780,int_stack+22824,int_stack+28449, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4410,int_stack+5211,int_stack+5151, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9327, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4590,int_stack+5301,int_stack+5211, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+9447, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4860,int_stack+4590,int_stack+4410, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+5457,int_stack+5427, 0.0, zero_stack, 1.0, int_stack+17451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19659,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+5502,int_stack+5457, 0.0, zero_stack, 1.0, int_stack+17481, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19689,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4410,int_stack+22824,int_stack+28449, 0.0, zero_stack, 1.0, int_stack+22734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25569,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4590,int_stack+5625,int_stack+5565, 0.0, zero_stack, 1.0, int_stack+7635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9891,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5220,int_stack+5715,int_stack+5625, 0.0, zero_stack, 1.0, int_stack+7755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10011,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31329,int_stack+5220,int_stack+4590, 0.0, zero_stack, 1.0, int_stack+23319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26199,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+5871,int_stack+5841, 0.0, zero_stack, 1.0, int_stack+18003, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19659, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+5916,int_stack+5871, 0.0, zero_stack, 1.0, int_stack+18033, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19689, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4590,int_stack+22824,int_stack+28449, 0.0, zero_stack, 1.0, int_stack+23499, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25569, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5220,int_stack+6039,int_stack+5979, 0.0, zero_stack, 1.0, int_stack+8199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9891, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5400,int_stack+6129,int_stack+6039, 0.0, zero_stack, 1.0, int_stack+8319, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10011, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5670,int_stack+5400,int_stack+5220, 0.0, zero_stack, 1.0, int_stack+24129, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26199, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+6285,int_stack+6255, 0.0, zero_stack, 1.0, int_stack+18555, 0.0, zero_stack, 1.0, int_stack+19659, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+6330,int_stack+6285, 0.0, zero_stack, 1.0, int_stack+18585, 0.0, zero_stack, 1.0, int_stack+19689, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5220,int_stack+22824,int_stack+28449, 0.0, zero_stack, 1.0, int_stack+24309, 0.0, zero_stack, 1.0, int_stack+25569, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5400,int_stack+6453,int_stack+6393, 0.0, zero_stack, 1.0, int_stack+8763, 0.0, zero_stack, 1.0, int_stack+9891, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6030,int_stack+6543,int_stack+6453, 0.0, zero_stack, 1.0, int_stack+8883, 0.0, zero_stack, 1.0, int_stack+10011, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6300,int_stack+6030,int_stack+5400, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+26199, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+6699,int_stack+6669, 0.0, zero_stack, 1.0, int_stack+19107, 1.0, int_stack+19659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+6744,int_stack+6699, 0.0, zero_stack, 1.0, int_stack+19137, 1.0, int_stack+19689, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5400,int_stack+22824,int_stack+28449, 0.0, zero_stack, 1.0, int_stack+24939, 1.0, int_stack+25569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6030,int_stack+6867,int_stack+6807, 0.0, zero_stack, 1.0, int_stack+9327, 1.0, int_stack+9891, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31689,int_stack+6957,int_stack+6867, 0.0, zero_stack, 1.0, int_stack+9447, 1.0, int_stack+10011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6660,int_stack+31689,int_stack+6030, 0.0, zero_stack, 1.0, int_stack+180, 1.0, int_stack+26199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+7113,int_stack+7083, 0.0, zero_stack, 2.0, int_stack+19659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+7158,int_stack+7113, 0.0, zero_stack, 2.0, int_stack+19689, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6030,int_stack+22824,int_stack+28449, 0.0, zero_stack, 2.0, int_stack+25569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+31689,int_stack+7281,int_stack+7221, 0.0, zero_stack, 2.0, int_stack+9891, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31869,int_stack+7371,int_stack+7281, 0.0, zero_stack, 2.0, int_stack+10011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32139,int_stack+31869,int_stack+31689, 0.0, zero_stack, 2.0, int_stack+26199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+7527,int_stack+7497, 1.0, int_stack+17451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20241,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22824,int_stack+7572,int_stack+7527, 1.0, int_stack+17481, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20316,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31689,int_stack+22824,int_stack+28449, 1.0, int_stack+22734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26379,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22734,int_stack+7845,int_stack+7695, 1.0, int_stack+7635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10515,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31869,int_stack+7935,int_stack+7845, 1.0, int_stack+7755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10635,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32499,int_stack+31869,int_stack+22734, 1.0, int_stack+23319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+8091,int_stack+8061, 1.0, int_stack+18003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20241, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23319,int_stack+8136,int_stack+8091, 1.0, int_stack+18033, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20316, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22734,int_stack+23319,int_stack+28449, 1.0, int_stack+23499, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26379, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23319,int_stack+8409,int_stack+8259, 1.0, int_stack+8199, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10515, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31869,int_stack+8499,int_stack+8409, 1.0, int_stack+8319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10635, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7020,int_stack+31869,int_stack+23319, 1.0, int_stack+24129, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+8655,int_stack+8625, 1.0, int_stack+18555, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20241, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24129,int_stack+8700,int_stack+8655, 1.0, int_stack+18585, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20316, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23319,int_stack+24129,int_stack+28449, 1.0, int_stack+24309, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26379, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24129,int_stack+8973,int_stack+8823, 1.0, int_stack+8763, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10515, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31869,int_stack+9063,int_stack+8973, 1.0, int_stack+8883, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10635, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7380,int_stack+31869,int_stack+24129, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28449,int_stack+9219,int_stack+9189, 1.0, int_stack+19107, 0.0, zero_stack, 1.0, int_stack+20241, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+9264,int_stack+9219, 1.0, int_stack+19137, 0.0, zero_stack, 1.0, int_stack+20316, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24129,int_stack+0,int_stack+28449, 1.0, int_stack+24939, 0.0, zero_stack, 1.0, int_stack+26379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+9537,int_stack+9387, 1.0, int_stack+9327, 0.0, zero_stack, 1.0, int_stack+10515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+31869,int_stack+9627,int_stack+9537, 1.0, int_stack+9447, 0.0, zero_stack, 1.0, int_stack+10635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7740,int_stack+31869,int_stack+0, 1.0, int_stack+180, 0.0, zero_stack, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24939,int_stack+9783,int_stack+9753, 1.0, int_stack+19659, 1.0, int_stack+20241, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+9828,int_stack+9783, 1.0, int_stack+19689, 1.0, int_stack+20316, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+135,int_stack+0,int_stack+24939, 1.0, int_stack+25569, 1.0, int_stack+26379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+31869,int_stack+10101,int_stack+9951, 1.0, int_stack+9891, 1.0, int_stack+10515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+10191,int_stack+10101, 1.0, int_stack+10011, 1.0, int_stack+10635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8370,int_stack+8100,int_stack+31869, 1.0, int_stack+26199, 1.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25569,int_stack+10347,int_stack+10317, 2.0, int_stack+20241, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+10392,int_stack+10347, 2.0, int_stack+20316, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26199,int_stack+0,int_stack+25569, 2.0, int_stack+26379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+31869,int_stack+10725,int_stack+10575, 2.0, int_stack+10515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8100,int_stack+10815,int_stack+10725, 2.0, int_stack+10635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8730,int_stack+8100,int_stack+31869, 2.0, int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26379,int_stack+10971,int_stack+10941, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20988,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+11016,int_stack+10971, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21018,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31869,int_stack+0,int_stack+26379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22644,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8100,int_stack+11139,int_stack+11079, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13149,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9090,int_stack+11229,int_stack+11139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13269,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9360,int_stack+9090,int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26379,int_stack+11385,int_stack+11355, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20988, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+11430,int_stack+11385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21018, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8100,int_stack+0,int_stack+26379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22644, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9090,int_stack+11553,int_stack+11493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13149, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9720,int_stack+11643,int_stack+11553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13269, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9990,int_stack+9720,int_stack+9090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26379,int_stack+11799,int_stack+11769, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20988, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+11844,int_stack+11799, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21018, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9090,int_stack+0,int_stack+26379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22644, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9720,int_stack+11967,int_stack+11907, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13149, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10350,int_stack+12057,int_stack+11967, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13269, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10620,int_stack+10350,int_stack+9720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26379,int_stack+12213,int_stack+12183, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+12258,int_stack+12213, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9720,int_stack+0,int_stack+26379, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10350,int_stack+12381,int_stack+12321, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13149, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10980,int_stack+12471,int_stack+12381, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11250,int_stack+10980,int_stack+10350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26379,int_stack+12627,int_stack+12597, 0.0, zero_stack, 1.0, int_stack+20988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+12672,int_stack+12627, 0.0, zero_stack, 1.0, int_stack+21018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10350,int_stack+0,int_stack+26379, 0.0, zero_stack, 1.0, int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10980,int_stack+12795,int_stack+12735, 0.0, zero_stack, 1.0, int_stack+13149, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11610,int_stack+12885,int_stack+12795, 0.0, zero_stack, 1.0, int_stack+13269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11880,int_stack+11610,int_stack+10980, 0.0, zero_stack, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26379,int_stack+13041,int_stack+13011, 1.0, int_stack+20988, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+13086,int_stack+13041, 1.0, int_stack+21018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10980,int_stack+0,int_stack+26379, 1.0, int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11610,int_stack+13359,int_stack+13209, 1.0, int_stack+13149, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12240,int_stack+13449,int_stack+13359, 1.0, int_stack+13269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12510,int_stack+12240,int_stack+11610, 1.0, int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+13605,int_stack+13575,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+13650,int_stack+13605,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11610,int_stack+0,int_stack+22644,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12240,int_stack+13773,int_stack+13713,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+315,int_stack+13863,int_stack+13773,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12870,int_stack+315,int_stack+12240,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+14019,int_stack+13989, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21540,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+14064,int_stack+14019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21570,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12240,int_stack+0,int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27189,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+315,int_stack+14187,int_stack+14127, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16197,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13230,int_stack+14277,int_stack+14187, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16317,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13500,int_stack+13230,int_stack+315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+14433,int_stack+14403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21540, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+14478,int_stack+14433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21570, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+315,int_stack+0,int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27189, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+495,int_stack+14601,int_stack+14541, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16197, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13230,int_stack+14691,int_stack+14601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16317, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13860,int_stack+13230,int_stack+495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+14847,int_stack+14817, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21540, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+14892,int_stack+14847, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21570, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+495,int_stack+0,int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27189, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13230,int_stack+15015,int_stack+14955, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16197, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14220,int_stack+15105,int_stack+15015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16317, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14490,int_stack+14220,int_stack+13230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+15261,int_stack+15231, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+15306,int_stack+15261, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13230,int_stack+0,int_stack+22644, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14220,int_stack+15429,int_stack+15369, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16197, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14850,int_stack+15519,int_stack+15429, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16317, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15120,int_stack+14850,int_stack+14220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+15675,int_stack+15645, 0.0, zero_stack, 1.0, int_stack+21540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+15720,int_stack+15675, 0.0, zero_stack, 1.0, int_stack+21570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14220,int_stack+0,int_stack+22644, 0.0, zero_stack, 1.0, int_stack+27189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14850,int_stack+15843,int_stack+15783, 0.0, zero_stack, 1.0, int_stack+16197, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15480,int_stack+15933,int_stack+15843, 0.0, zero_stack, 1.0, int_stack+16317, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+32859,int_stack+15480,int_stack+14850, 0.0, zero_stack, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22644,int_stack+16089,int_stack+16059, 1.0, int_stack+21540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+16134,int_stack+16089, 1.0, int_stack+21570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14850,int_stack+0,int_stack+22644, 1.0, int_stack+27189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15480,int_stack+16407,int_stack+16257, 1.0, int_stack+16197, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15660,int_stack+16497,int_stack+16407, 1.0, int_stack+16317, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15930,int_stack+15660,int_stack+15480, 1.0, int_stack+720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+16653,int_stack+16623,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+16698,int_stack+16653,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+15480,int_stack+0,int_stack+27189,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+15660,int_stack+16821,int_stack+16761,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16290,int_stack+16911,int_stack+16821,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+16560,int_stack+16290,int_stack+15660,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+17067,int_stack+17037,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+17112,int_stack+17067,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+15660,int_stack+0,int_stack+27189,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16290,int_stack+17235,int_stack+17175,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16920,int_stack+17325,int_stack+17235,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17190,int_stack+16920,int_stack+16290,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+17619,int_stack+17589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22092,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+17664,int_stack+17619, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22122,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16290,int_stack+0,int_stack+27189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27819,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16920,int_stack+17787,int_stack+17727, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20562,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33219,int_stack+17877,int_stack+17787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20682,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17550,int_stack+33219,int_stack+16920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+18171,int_stack+18141, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22092, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+18216,int_stack+18171, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22122, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16920,int_stack+0,int_stack+27189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27819, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33219,int_stack+18339,int_stack+18279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20562, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33399,int_stack+18429,int_stack+18339, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20682, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17910,int_stack+33399,int_stack+33219, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+18723,int_stack+18693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22092, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+18768,int_stack+18723, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22122, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33219,int_stack+0,int_stack+27189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27819, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33399,int_stack+18891,int_stack+18831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20562, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+33579,int_stack+18981,int_stack+18891, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20682, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18270,int_stack+33579,int_stack+33399, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+19275,int_stack+19245, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+19320,int_stack+19275, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33399,int_stack+0,int_stack+27189, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33579,int_stack+19443,int_stack+19383, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18630,int_stack+19533,int_stack+19443, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18900,int_stack+18630,int_stack+33579, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+19827,int_stack+19797, 0.0, zero_stack, 1.0, int_stack+22092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+19872,int_stack+19827, 0.0, zero_stack, 1.0, int_stack+22122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33579,int_stack+0,int_stack+27189, 0.0, zero_stack, 1.0, int_stack+27819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33759,int_stack+19995,int_stack+19935, 0.0, zero_stack, 1.0, int_stack+20562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18630,int_stack+20085,int_stack+19995, 0.0, zero_stack, 1.0, int_stack+20682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19260,int_stack+18630,int_stack+33759, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27189,int_stack+20454,int_stack+20424, 1.0, int_stack+22092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+20499,int_stack+20454, 1.0, int_stack+22122, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33759,int_stack+0,int_stack+27189, 1.0, int_stack+27819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18630,int_stack+20772,int_stack+20622, 1.0, int_stack+20562, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19620,int_stack+20862,int_stack+20772, 1.0, int_stack+20682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19890,int_stack+19620,int_stack+18630, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27819,int_stack+21156,int_stack+21126,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+21201,int_stack+21156,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+18630,int_stack+0,int_stack+27819,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19620,int_stack+21324,int_stack+21264,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+20250,int_stack+21414,int_stack+21324,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+20520,int_stack+20250,int_stack+19620,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27819,int_stack+21708,int_stack+21678,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+21753,int_stack+21708,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+19620,int_stack+0,int_stack+27819,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+20250,int_stack+21876,int_stack+21816,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+20880,int_stack+21966,int_stack+21876,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+21150,int_stack+20880,int_stack+20250,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+27819,int_stack+22260,int_stack+22230,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+22305,int_stack+22260,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+20250,int_stack+0,int_stack+27819,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+20880,int_stack+22428,int_stack+22368,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+21510,int_stack+22518,int_stack+22428,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+21780,int_stack+21510,int_stack+20880,6);
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+22140,int_stack+23769,int_stack+22959,60);
     Libderiv->ABCD[11] = int_stack + 22140;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+33939,int_stack+24579,int_stack+23589,60);
     Libderiv->ABCD[10] = int_stack + 33939;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+34479,int_stack+25209,int_stack+24399,60);
     Libderiv->ABCD[9] = int_stack + 34479;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+35019,int_stack+25839,int_stack+25029,60);
     Libderiv->ABCD[8] = int_stack + 35019;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+35559,int_stack+26649,int_stack+25659,60);
     Libderiv->ABCD[7] = int_stack + 35559;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+36099,int_stack+27279,int_stack+26469,60);
     Libderiv->ABCD[6] = int_stack + 36099;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+36639,int_stack+27909,int_stack+27009, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 36639;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+37179,int_stack+28539,int_stack+27639, 0.0, zero_stack, 1.0, int_stack+23139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 37179;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+37719,int_stack+29169,int_stack+28269, 1.0, int_stack+23139, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 37719;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+38259,int_stack+1260,int_stack+28899,60);
     Libderiv->ABCD[155] = int_stack + 38259;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+28449,int_stack+1620,int_stack+1080,60);
     Libderiv->ABCD[143] = int_stack + 28449;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+28989,int_stack+29889,int_stack+29529,60);
     Libderiv->ABCD[142] = int_stack + 28989;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+675,int_stack+2430,int_stack+29709,60);
     Libderiv->ABCD[131] = int_stack + 675;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+29529,int_stack+30249,int_stack+1980,60);
     Libderiv->ABCD[130] = int_stack + 29529;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+30069,int_stack+3240,int_stack+2160,60);
     Libderiv->ABCD[129] = int_stack + 30069;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1215,int_stack+30609,int_stack+2790,60);
     Libderiv->ABCD[119] = int_stack + 1215;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1755,int_stack+4050,int_stack+2970,60);
     Libderiv->ABCD[118] = int_stack + 1755;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2295,int_stack+30969,int_stack+3600,60);
     Libderiv->ABCD[117] = int_stack + 2295;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+30609,int_stack+4860,int_stack+3780,60);
     Libderiv->ABCD[116] = int_stack + 30609;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2835,int_stack+31329,int_stack+4410,60);
     Libderiv->ABCD[107] = int_stack + 2835;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+31149,int_stack+5670,int_stack+4590,60);
     Libderiv->ABCD[106] = int_stack + 31149;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3375,int_stack+6300,int_stack+5220,60);
     Libderiv->ABCD[105] = int_stack + 3375;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3915,int_stack+6660,int_stack+5400,60);
     Libderiv->ABCD[104] = int_stack + 3915;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4455,int_stack+32139,int_stack+6030,60);
     Libderiv->ABCD[103] = int_stack + 4455;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4995,int_stack+32499,int_stack+31689,60);
     Libderiv->ABCD[95] = int_stack + 4995;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+5535,int_stack+7020,int_stack+22734,60);
     Libderiv->ABCD[94] = int_stack + 5535;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+6075,int_stack+7380,int_stack+23319,60);
     Libderiv->ABCD[93] = int_stack + 6075;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+6615,int_stack+7740,int_stack+24129,60);
     Libderiv->ABCD[92] = int_stack + 6615;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+23769,int_stack+8370,int_stack+135,60);
     Libderiv->ABCD[91] = int_stack + 23769;
 /*--- compute (pp|fd) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+7155,int_stack+8730,int_stack+26199,60);
     Libderiv->ABCD[90] = int_stack + 7155;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+25839,int_stack+9360,int_stack+31869, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22959, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[47] = int_stack + 25839;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+31689,int_stack+9990,int_stack+8100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[46] = int_stack + 31689;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+32229,int_stack+10620,int_stack+9090, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24399, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[45] = int_stack + 32229;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+7695,int_stack+11250,int_stack+9720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25029, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[44] = int_stack + 7695;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8235,int_stack+11880,int_stack+10350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[43] = int_stack + 8235;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+8775,int_stack+12510,int_stack+10980, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[42] = int_stack + 8775;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9315,int_stack+12870,int_stack+11610, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27009, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[38] = int_stack + 9315;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+9855,int_stack+13500,int_stack+12240, 0.0, zero_stack, 1.0, int_stack+22959, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[35] = int_stack + 9855;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+10395,int_stack+13860,int_stack+315, 0.0, zero_stack, 1.0, int_stack+23589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[34] = int_stack + 10395;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+13410,int_stack+14490,int_stack+495, 0.0, zero_stack, 1.0, int_stack+24399, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[33] = int_stack + 13410;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+15120,int_stack+13230, 0.0, zero_stack, 1.0, int_stack+25029, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[32] = int_stack + 0;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+10935,int_stack+32859,int_stack+14220, 0.0, zero_stack, 1.0, int_stack+25659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[31] = int_stack + 10935;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+13950,int_stack+15930,int_stack+14850, 0.0, zero_stack, 1.0, int_stack+26469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[30] = int_stack + 13950;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+14490,int_stack+16560,int_stack+15480, 0.0, zero_stack, 1.0, int_stack+27009, 1.0, int_stack+27639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[26] = int_stack + 14490;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+15030,int_stack+17190,int_stack+15660, 0.0, zero_stack, 2.0, int_stack+27639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[25] = int_stack + 15030;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+15570,int_stack+17550,int_stack+16290, 1.0, int_stack+22959, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[23] = int_stack + 15570;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+17100,int_stack+17910,int_stack+16920, 1.0, int_stack+23589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[22] = int_stack + 17100;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+17640,int_stack+18270,int_stack+33219, 1.0, int_stack+24399, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[21] = int_stack + 17640;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+32769,int_stack+18900,int_stack+33399, 1.0, int_stack+25029, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[20] = int_stack + 32769;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+24309,int_stack+19260,int_stack+33579, 1.0, int_stack+25659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[19] = int_stack + 24309;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+24849,int_stack+19890,int_stack+33759, 1.0, int_stack+26469, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[18] = int_stack + 24849;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+33309,int_stack+20520,int_stack+18630, 1.0, int_stack+27009, 0.0, zero_stack, 1.0, int_stack+28269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[14] = int_stack + 33309;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+26379,int_stack+21150,int_stack+19620, 1.0, int_stack+27639, 1.0, int_stack+28269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[13] = int_stack + 26379;
 /*--- compute (pp|fd) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+26919,int_stack+21780,int_stack+20250, 2.0, int_stack+28269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[12] = int_stack + 26919;

}
