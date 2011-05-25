#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_fdfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fd|fd) integrals */

void d12hrr_order_fdfd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[5][5][11] = int_stack + 0;
 Libderiv->deriv_classes[5][5][10] = int_stack + 441;
 Libderiv->deriv_classes[5][5][9] = int_stack + 882;
 Libderiv->deriv_classes[5][5][8] = int_stack + 1323;
 Libderiv->deriv_classes[5][5][7] = int_stack + 1764;
 Libderiv->dvrr_classes[5][4] = int_stack + 2205;
 Libderiv->deriv_classes[5][5][6] = int_stack + 2520;
 Libderiv->deriv_classes[5][5][2] = int_stack + 2961;
 Libderiv->deriv_classes[5][5][1] = int_stack + 3402;
 Libderiv->dvrr_classes[4][5] = int_stack + 3843;
 Libderiv->deriv_classes[5][5][0] = int_stack + 4158;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 4599;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 4699;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 4849;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 5059;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 5209;
 Libderiv->deriv2_classes[4][5][143] = int_stack + 5434;
 Libderiv->deriv2_classes[5][3][143] = int_stack + 5749;
 Libderiv->deriv2_classes[5][4][143] = int_stack + 5959;
 Libderiv->deriv2_classes[5][5][143] = int_stack + 6274;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 6715;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 6815;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 6965;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 7175;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 7325;
 Libderiv->deriv2_classes[4][5][131] = int_stack + 7550;
 Libderiv->deriv2_classes[5][3][131] = int_stack + 7865;
 Libderiv->deriv2_classes[5][4][131] = int_stack + 8075;
 Libderiv->deriv2_classes[5][5][131] = int_stack + 8390;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 8831;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 8931;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 9081;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 9291;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 9441;
 Libderiv->deriv2_classes[4][5][130] = int_stack + 9666;
 Libderiv->deriv2_classes[5][3][130] = int_stack + 9981;
 Libderiv->deriv2_classes[5][4][130] = int_stack + 10191;
 Libderiv->deriv2_classes[5][5][130] = int_stack + 10506;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 10947;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 11047;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 11197;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 11407;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 11557;
 Libderiv->deriv2_classes[4][5][119] = int_stack + 11782;
 Libderiv->deriv2_classes[5][3][119] = int_stack + 12097;
 Libderiv->deriv2_classes[5][4][119] = int_stack + 12307;
 Libderiv->deriv2_classes[5][5][119] = int_stack + 12622;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 13063;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 13163;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 13313;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 13523;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 13673;
 Libderiv->deriv2_classes[4][5][118] = int_stack + 13898;
 Libderiv->deriv2_classes[5][3][118] = int_stack + 14213;
 Libderiv->deriv2_classes[5][4][118] = int_stack + 14423;
 Libderiv->deriv2_classes[5][5][118] = int_stack + 14738;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 15179;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 15279;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 15429;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 15639;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 15789;
 Libderiv->deriv2_classes[4][5][117] = int_stack + 16014;
 Libderiv->deriv2_classes[5][3][117] = int_stack + 16329;
 Libderiv->deriv2_classes[5][4][117] = int_stack + 16539;
 Libderiv->deriv2_classes[5][5][117] = int_stack + 16854;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 17295;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 17395;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 17545;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 17755;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 17905;
 Libderiv->deriv2_classes[4][5][107] = int_stack + 18130;
 Libderiv->deriv2_classes[5][3][107] = int_stack + 18445;
 Libderiv->deriv2_classes[5][4][107] = int_stack + 18655;
 Libderiv->deriv2_classes[5][5][107] = int_stack + 18970;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 19411;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 19511;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 19661;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 19871;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 20021;
 Libderiv->deriv2_classes[4][5][106] = int_stack + 20246;
 Libderiv->deriv2_classes[5][3][106] = int_stack + 20561;
 Libderiv->deriv2_classes[5][4][106] = int_stack + 20771;
 Libderiv->deriv2_classes[5][5][106] = int_stack + 21086;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 21527;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 21627;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 21777;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 21987;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 22137;
 Libderiv->deriv2_classes[4][5][105] = int_stack + 22362;
 Libderiv->deriv2_classes[5][3][105] = int_stack + 22677;
 Libderiv->deriv2_classes[5][4][105] = int_stack + 22887;
 Libderiv->deriv2_classes[5][5][105] = int_stack + 23202;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 23643;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 23743;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 23893;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 24103;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 24253;
 Libderiv->deriv2_classes[4][5][104] = int_stack + 24478;
 Libderiv->deriv2_classes[5][3][104] = int_stack + 24793;
 Libderiv->deriv2_classes[5][4][104] = int_stack + 25003;
 Libderiv->deriv2_classes[5][5][104] = int_stack + 25318;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 25759;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 25859;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 26009;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 26219;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 26369;
 Libderiv->deriv2_classes[4][5][95] = int_stack + 26594;
 Libderiv->deriv2_classes[5][3][95] = int_stack + 26909;
 Libderiv->deriv2_classes[5][4][95] = int_stack + 27119;
 Libderiv->deriv2_classes[5][5][95] = int_stack + 27434;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 27875;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 27975;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 28125;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 28335;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 28485;
 Libderiv->deriv2_classes[4][5][94] = int_stack + 28710;
 Libderiv->deriv2_classes[5][3][94] = int_stack + 29025;
 Libderiv->deriv2_classes[5][4][94] = int_stack + 29235;
 Libderiv->deriv2_classes[5][5][94] = int_stack + 29550;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 29991;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 30091;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 30241;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 30451;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 30601;
 Libderiv->deriv2_classes[4][5][93] = int_stack + 30826;
 Libderiv->deriv2_classes[5][3][93] = int_stack + 31141;
 Libderiv->deriv2_classes[5][4][93] = int_stack + 31351;
 Libderiv->deriv2_classes[5][5][93] = int_stack + 31666;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 32107;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 32207;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 32357;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 32567;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 32717;
 Libderiv->deriv2_classes[4][5][92] = int_stack + 32942;
 Libderiv->deriv2_classes[5][3][92] = int_stack + 33257;
 Libderiv->deriv2_classes[5][4][92] = int_stack + 33467;
 Libderiv->deriv2_classes[5][5][92] = int_stack + 33782;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 34223;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 34323;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 34473;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 34683;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 34833;
 Libderiv->deriv2_classes[4][5][91] = int_stack + 35058;
 Libderiv->deriv2_classes[5][3][91] = int_stack + 35373;
 Libderiv->deriv2_classes[5][4][91] = int_stack + 35583;
 Libderiv->deriv2_classes[5][5][91] = int_stack + 35898;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 36339;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 36439;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 36589;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 36799;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 36949;
 Libderiv->deriv2_classes[4][5][83] = int_stack + 37174;
 Libderiv->deriv_classes[5][3][11] = int_stack + 37489;
 Libderiv->deriv2_classes[5][3][83] = int_stack + 37699;
 Libderiv->deriv_classes[5][4][11] = int_stack + 37909;
 Libderiv->deriv2_classes[5][4][83] = int_stack + 38224;
 Libderiv->deriv2_classes[5][5][83] = int_stack + 38539;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 38980;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 39080;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 39230;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 39440;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 39590;
 Libderiv->deriv2_classes[4][5][82] = int_stack + 39815;
 Libderiv->deriv_classes[5][3][10] = int_stack + 40130;
 Libderiv->deriv2_classes[5][3][82] = int_stack + 40340;
 Libderiv->deriv_classes[5][4][10] = int_stack + 40550;
 Libderiv->deriv2_classes[5][4][82] = int_stack + 40865;
 Libderiv->deriv2_classes[5][5][82] = int_stack + 41180;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 41621;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 41721;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 41871;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 42081;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 42231;
 Libderiv->deriv2_classes[4][5][81] = int_stack + 42456;
 Libderiv->deriv_classes[5][3][9] = int_stack + 42771;
 Libderiv->deriv2_classes[5][3][81] = int_stack + 42981;
 Libderiv->deriv_classes[5][4][9] = int_stack + 43191;
 Libderiv->deriv2_classes[5][4][81] = int_stack + 43506;
 Libderiv->deriv2_classes[5][5][81] = int_stack + 43821;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 44262;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 44362;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 44512;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 44722;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 44872;
 Libderiv->deriv2_classes[4][5][80] = int_stack + 45097;
 Libderiv->deriv_classes[5][3][8] = int_stack + 45412;
 Libderiv->deriv2_classes[5][3][80] = int_stack + 45622;
 Libderiv->deriv_classes[5][4][8] = int_stack + 45832;
 Libderiv->deriv2_classes[5][4][80] = int_stack + 46147;
 Libderiv->deriv2_classes[5][5][80] = int_stack + 46462;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 46903;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 47003;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 47153;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 47363;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 47513;
 Libderiv->deriv2_classes[4][5][79] = int_stack + 47738;
 Libderiv->deriv_classes[5][3][7] = int_stack + 48053;
 Libderiv->deriv2_classes[5][3][79] = int_stack + 48263;
 Libderiv->deriv_classes[5][4][7] = int_stack + 48473;
 Libderiv->deriv2_classes[5][4][79] = int_stack + 48788;
 Libderiv->deriv2_classes[5][5][79] = int_stack + 49103;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 49544;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 49644;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 49794;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 50004;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 50154;
 Libderiv->deriv2_classes[4][5][78] = int_stack + 50379;
 Libderiv->dvrr_classes[5][3] = int_stack + 50694;
 Libderiv->deriv_classes[5][3][6] = int_stack + 50904;
 Libderiv->deriv2_classes[5][3][78] = int_stack + 51114;
 Libderiv->deriv_classes[5][4][6] = int_stack + 51324;
 Libderiv->deriv2_classes[5][4][78] = int_stack + 51639;
 Libderiv->deriv2_classes[5][5][78] = int_stack + 51954;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 52395;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 52495;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 52645;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 52855;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 53005;
 Libderiv->deriv2_classes[4][5][35] = int_stack + 53230;
 Libderiv->deriv2_classes[5][3][35] = int_stack + 53545;
 Libderiv->deriv2_classes[5][4][35] = int_stack + 53755;
 Libderiv->deriv2_classes[5][5][35] = int_stack + 54070;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 54511;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 54611;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 54761;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 54971;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 55121;
 Libderiv->deriv2_classes[4][5][34] = int_stack + 55346;
 Libderiv->deriv2_classes[5][3][34] = int_stack + 55661;
 Libderiv->deriv2_classes[5][4][34] = int_stack + 55871;
 Libderiv->deriv2_classes[5][5][34] = int_stack + 56186;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 56627;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 56727;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 56877;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 57087;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 57237;
 Libderiv->deriv2_classes[4][5][33] = int_stack + 57462;
 Libderiv->deriv2_classes[5][3][33] = int_stack + 57777;
 Libderiv->deriv2_classes[5][4][33] = int_stack + 57987;
 Libderiv->deriv2_classes[5][5][33] = int_stack + 58302;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 58743;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 58843;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 58993;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 59203;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 59353;
 Libderiv->deriv2_classes[4][5][32] = int_stack + 59578;
 Libderiv->deriv2_classes[5][3][32] = int_stack + 59893;
 Libderiv->deriv2_classes[5][4][32] = int_stack + 60103;
 Libderiv->deriv2_classes[5][5][32] = int_stack + 60418;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 60859;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 60959;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 61109;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 61319;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 61469;
 Libderiv->deriv2_classes[4][5][31] = int_stack + 61694;
 Libderiv->deriv2_classes[5][3][31] = int_stack + 62009;
 Libderiv->deriv2_classes[5][4][31] = int_stack + 62219;
 Libderiv->deriv2_classes[5][5][31] = int_stack + 62534;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 62975;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 63075;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 63225;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 63435;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 63585;
 Libderiv->deriv2_classes[4][5][30] = int_stack + 63810;
 Libderiv->deriv_classes[5][3][2] = int_stack + 64125;
 Libderiv->deriv2_classes[5][3][30] = int_stack + 64335;
 Libderiv->deriv_classes[5][4][2] = int_stack + 64545;
 Libderiv->deriv2_classes[5][4][30] = int_stack + 64860;
 Libderiv->deriv2_classes[5][5][30] = int_stack + 65175;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 65616;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 65716;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 65866;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 66076;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 66226;
 Libderiv->deriv2_classes[4][5][26] = int_stack + 66451;
 Libderiv->deriv2_classes[5][3][26] = int_stack + 66766;
 Libderiv->deriv2_classes[5][4][26] = int_stack + 66976;
 Libderiv->deriv2_classes[5][5][26] = int_stack + 67291;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 67732;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 67832;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 67982;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 68192;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 68342;
 Libderiv->deriv2_classes[4][5][23] = int_stack + 68567;
 Libderiv->deriv2_classes[5][3][23] = int_stack + 68882;
 Libderiv->deriv2_classes[5][4][23] = int_stack + 69092;
 Libderiv->deriv2_classes[5][5][23] = int_stack + 69407;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 69848;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 69948;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 70098;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 70308;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 70458;
 Libderiv->deriv2_classes[4][5][22] = int_stack + 70683;
 Libderiv->deriv2_classes[5][3][22] = int_stack + 70998;
 Libderiv->deriv2_classes[5][4][22] = int_stack + 71208;
 Libderiv->deriv2_classes[5][5][22] = int_stack + 71523;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 71964;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 72064;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 72214;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 72424;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 72574;
 Libderiv->deriv2_classes[4][5][21] = int_stack + 72799;
 Libderiv->deriv2_classes[5][3][21] = int_stack + 73114;
 Libderiv->deriv2_classes[5][4][21] = int_stack + 73324;
 Libderiv->deriv2_classes[5][5][21] = int_stack + 73639;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 74080;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 74180;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 74330;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 74540;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 74690;
 Libderiv->deriv2_classes[4][5][20] = int_stack + 74915;
 Libderiv->deriv2_classes[5][3][20] = int_stack + 75230;
 Libderiv->deriv2_classes[5][4][20] = int_stack + 75440;
 Libderiv->deriv2_classes[5][5][20] = int_stack + 75755;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 76196;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 76296;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 76446;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 76656;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 76806;
 Libderiv->deriv2_classes[4][5][19] = int_stack + 77031;
 Libderiv->deriv2_classes[5][3][19] = int_stack + 77346;
 Libderiv->deriv2_classes[5][4][19] = int_stack + 77556;
 Libderiv->deriv2_classes[5][5][19] = int_stack + 77871;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 78312;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 78412;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 78562;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 78772;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 78922;
 Libderiv->deriv2_classes[4][5][18] = int_stack + 79147;
 Libderiv->deriv_classes[5][3][1] = int_stack + 79462;
 Libderiv->deriv2_classes[5][3][18] = int_stack + 79672;
 Libderiv->deriv_classes[5][4][1] = int_stack + 79882;
 Libderiv->deriv2_classes[5][4][18] = int_stack + 80197;
 Libderiv->deriv2_classes[5][5][18] = int_stack + 80512;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 80953;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 81053;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 81203;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 81413;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 81563;
 Libderiv->deriv2_classes[4][5][14] = int_stack + 81788;
 Libderiv->deriv2_classes[5][3][14] = int_stack + 82103;
 Libderiv->deriv2_classes[5][4][14] = int_stack + 82313;
 Libderiv->deriv2_classes[5][5][14] = int_stack + 82628;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 83069;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 83169;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 83319;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 83529;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 83679;
 Libderiv->deriv2_classes[4][5][13] = int_stack + 83904;
 Libderiv->deriv2_classes[5][3][13] = int_stack + 84219;
 Libderiv->deriv2_classes[5][4][13] = int_stack + 84429;
 Libderiv->deriv2_classes[5][5][13] = int_stack + 84744;
 Libderiv->deriv_classes[3][3][11] = int_stack + 85185;
 Libderiv->deriv_classes[3][4][11] = int_stack + 85285;
 Libderiv->deriv_classes[3][5][11] = int_stack + 85435;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 85645;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 85745;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 85895;
 Libderiv->deriv_classes[4][3][11] = int_stack + 86105;
 Libderiv->deriv_classes[4][4][11] = int_stack + 86255;
 Libderiv->deriv_classes[4][5][11] = int_stack + 86480;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 86795;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 86945;
 Libderiv->deriv2_classes[4][5][11] = int_stack + 87170;
 Libderiv->deriv2_classes[5][3][11] = int_stack + 87485;
 Libderiv->deriv2_classes[5][4][11] = int_stack + 87695;
 Libderiv->deriv2_classes[5][5][11] = int_stack + 88010;
 Libderiv->deriv_classes[3][3][10] = int_stack + 88451;
 Libderiv->deriv_classes[3][4][10] = int_stack + 88551;
 Libderiv->deriv_classes[3][5][10] = int_stack + 88701;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 88911;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 89011;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 89161;
 Libderiv->deriv_classes[4][3][10] = int_stack + 89371;
 Libderiv->deriv_classes[4][4][10] = int_stack + 89521;
 Libderiv->deriv_classes[4][5][10] = int_stack + 89746;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 90061;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 90211;
 Libderiv->deriv2_classes[4][5][10] = int_stack + 90436;
 Libderiv->deriv2_classes[5][3][10] = int_stack + 90751;
 Libderiv->deriv2_classes[5][4][10] = int_stack + 90961;
 Libderiv->deriv2_classes[5][5][10] = int_stack + 91276;
 Libderiv->deriv_classes[3][3][9] = int_stack + 91717;
 Libderiv->deriv_classes[3][4][9] = int_stack + 91817;
 Libderiv->deriv_classes[3][5][9] = int_stack + 91967;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 92177;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 92277;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 92427;
 Libderiv->deriv_classes[4][3][9] = int_stack + 92637;
 Libderiv->deriv_classes[4][4][9] = int_stack + 92787;
 Libderiv->deriv_classes[4][5][9] = int_stack + 93012;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 93327;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 93477;
 Libderiv->deriv2_classes[4][5][9] = int_stack + 93702;
 Libderiv->deriv2_classes[5][3][9] = int_stack + 94017;
 Libderiv->deriv2_classes[5][4][9] = int_stack + 94227;
 Libderiv->deriv2_classes[5][5][9] = int_stack + 94542;
 Libderiv->deriv_classes[3][3][8] = int_stack + 94983;
 Libderiv->deriv_classes[3][4][8] = int_stack + 95083;
 Libderiv->deriv_classes[3][5][8] = int_stack + 95233;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 95443;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 95543;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 95693;
 Libderiv->deriv_classes[4][3][8] = int_stack + 95903;
 Libderiv->deriv_classes[4][4][8] = int_stack + 96053;
 Libderiv->deriv_classes[4][5][8] = int_stack + 96278;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 96593;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 96743;
 Libderiv->deriv2_classes[4][5][8] = int_stack + 96968;
 Libderiv->deriv2_classes[5][3][8] = int_stack + 97283;
 Libderiv->deriv2_classes[5][4][8] = int_stack + 97493;
 Libderiv->deriv2_classes[5][5][8] = int_stack + 97808;
 Libderiv->deriv_classes[3][3][7] = int_stack + 98249;
 Libderiv->deriv_classes[3][4][7] = int_stack + 98349;
 Libderiv->deriv_classes[3][5][7] = int_stack + 98499;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 98709;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 98809;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 98959;
 Libderiv->deriv_classes[4][3][7] = int_stack + 99169;
 Libderiv->deriv_classes[4][4][7] = int_stack + 99319;
 Libderiv->deriv_classes[4][5][7] = int_stack + 99544;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 99859;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 100009;
 Libderiv->deriv2_classes[4][5][7] = int_stack + 100234;
 Libderiv->deriv2_classes[5][3][7] = int_stack + 100549;
 Libderiv->deriv2_classes[5][4][7] = int_stack + 100759;
 Libderiv->deriv2_classes[5][5][7] = int_stack + 101074;
 Libderiv->deriv_classes[3][3][6] = int_stack + 101515;
 Libderiv->deriv_classes[3][4][6] = int_stack + 101615;
 Libderiv->deriv_classes[3][5][6] = int_stack + 101765;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 101975;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 102075;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 102225;
 Libderiv->dvrr_classes[4][3] = int_stack + 102435;
 Libderiv->deriv_classes[4][3][6] = int_stack + 102585;
 Libderiv->dvrr_classes[4][4] = int_stack + 102735;
 Libderiv->deriv_classes[4][4][6] = int_stack + 102960;
 Libderiv->deriv_classes[4][5][6] = int_stack + 103185;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 103500;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 103650;
 Libderiv->deriv2_classes[4][5][6] = int_stack + 103875;
 Libderiv->deriv_classes[5][3][0] = int_stack + 104190;
 Libderiv->deriv2_classes[5][3][6] = int_stack + 104400;
 Libderiv->deriv_classes[5][4][0] = int_stack + 104610;
 Libderiv->deriv2_classes[5][4][6] = int_stack + 104925;
 Libderiv->deriv2_classes[5][5][6] = int_stack + 105240;
 Libderiv->deriv_classes[3][3][2] = int_stack + 105681;
 Libderiv->deriv_classes[3][4][2] = int_stack + 105781;
 Libderiv->deriv_classes[3][5][2] = int_stack + 105931;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 106141;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 106241;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 106391;
 Libderiv->deriv_classes[4][3][2] = int_stack + 106601;
 Libderiv->deriv_classes[4][4][2] = int_stack + 106751;
 Libderiv->deriv_classes[4][5][2] = int_stack + 106976;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 107291;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 107441;
 Libderiv->deriv2_classes[4][5][2] = int_stack + 107666;
 Libderiv->deriv2_classes[5][3][2] = int_stack + 107981;
 Libderiv->deriv2_classes[5][4][2] = int_stack + 108191;
 Libderiv->deriv2_classes[5][5][2] = int_stack + 108506;
 Libderiv->deriv_classes[3][3][1] = int_stack + 108947;
 Libderiv->deriv_classes[3][4][1] = int_stack + 109047;
 Libderiv->deriv_classes[3][5][1] = int_stack + 109197;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 109407;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 109507;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 109657;
 Libderiv->deriv_classes[4][3][1] = int_stack + 109867;
 Libderiv->deriv_classes[4][4][1] = int_stack + 110017;
 Libderiv->deriv_classes[4][5][1] = int_stack + 110242;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 110557;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 110707;
 Libderiv->deriv2_classes[4][5][1] = int_stack + 110932;
 Libderiv->deriv2_classes[5][3][1] = int_stack + 111247;
 Libderiv->deriv2_classes[5][4][1] = int_stack + 111457;
 Libderiv->deriv2_classes[5][5][1] = int_stack + 111772;
 Libderiv->dvrr_classes[3][3] = int_stack + 112213;
 Libderiv->dvrr_classes[3][4] = int_stack + 112313;
 Libderiv->dvrr_classes[3][5] = int_stack + 112463;
 Libderiv->deriv_classes[3][3][0] = int_stack + 112673;
 Libderiv->deriv_classes[3][4][0] = int_stack + 112773;
 Libderiv->deriv_classes[3][5][0] = int_stack + 112923;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 113133;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 113233;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 113383;
 Libderiv->deriv_classes[4][3][0] = int_stack + 113593;
 Libderiv->deriv_classes[4][4][0] = int_stack + 113743;
 Libderiv->deriv_classes[4][5][0] = int_stack + 113968;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 114283;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 114433;
 Libderiv->deriv2_classes[4][5][0] = int_stack + 114658;
 Libderiv->deriv2_classes[5][3][0] = int_stack + 114973;
 Libderiv->deriv2_classes[5][4][0] = int_stack + 115183;
 Libderiv->deriv2_classes[5][5][0] = int_stack + 115498;
 memset(int_stack,0,927512);

 Libderiv->dvrr_stack = int_stack + 280324;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_fdfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+115939,int_stack+112313,int_stack+112213,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+85285,int_stack+85185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112213,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+116539,int_stack+85435,int_stack+85285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112313,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+116989,int_stack+116539,int_stack+116239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115939,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+116539,int_stack+102735,int_stack+102435,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+86255,int_stack+86105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102435,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+118039,int_stack+86480,int_stack+86255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102735,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+118714,int_stack+118039,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116539,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+119614,int_stack+118714,int_stack+116989,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+118039,int_stack+2205,int_stack+50694,21);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+121414,int_stack+37909,int_stack+37489, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50694,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+122044,int_stack+0,int_stack+37909, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2205,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+122989,int_stack+122044,int_stack+121414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118039,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+124249,int_stack+122989,int_stack+118714,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+122044,int_stack+88551,int_stack+88451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112213, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+122344,int_stack+88701,int_stack+88551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112313, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+122794,int_stack+122344,int_stack+122044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115939, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+122344,int_stack+89521,int_stack+89371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102435, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+123394,int_stack+89746,int_stack+89521, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102735, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+126949,int_stack+123394,int_stack+122344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116539, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+127849,int_stack+126949,int_stack+122794,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+123394,int_stack+40550,int_stack+40130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50694, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+129649,int_stack+441,int_stack+40550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2205, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130594,int_stack+129649,int_stack+123394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118039, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+131854,int_stack+130594,int_stack+126949,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+129649,int_stack+91817,int_stack+91717, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112213, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+129949,int_stack+91967,int_stack+91817, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112313, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130399,int_stack+129949,int_stack+129649, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115939, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+129949,int_stack+92787,int_stack+92637, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102435, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130999,int_stack+93012,int_stack+92787, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102735, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+134554,int_stack+130999,int_stack+129949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116539, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+135454,int_stack+134554,int_stack+130399,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+130999,int_stack+43191,int_stack+42771, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50694, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+137254,int_stack+882,int_stack+43191, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2205, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+137254,int_stack+130999, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118039, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+137254,int_stack+0,int_stack+134554,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+95083,int_stack+94983, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+300,int_stack+95233,int_stack+95083, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+139954,int_stack+300,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+96053,int_stack+95903, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+96278,int_stack+96053, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+141229,int_stack+140554,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+142129,int_stack+141229,int_stack+139954,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+45832,int_stack+45412, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+143929,int_stack+1323,int_stack+45832, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+144874,int_stack+143929,int_stack+140554, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+146134,int_stack+144874,int_stack+141229,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+143929,int_stack+98349,int_stack+98249, 0.0, zero_stack, 1.0, int_stack+112213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+144229,int_stack+98499,int_stack+98349, 0.0, zero_stack, 1.0, int_stack+112313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+144679,int_stack+144229,int_stack+143929, 0.0, zero_stack, 1.0, int_stack+115939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+144229,int_stack+99319,int_stack+99169, 0.0, zero_stack, 1.0, int_stack+102435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+145279,int_stack+99544,int_stack+99319, 0.0, zero_stack, 1.0, int_stack+102735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+750,int_stack+145279,int_stack+144229, 0.0, zero_stack, 1.0, int_stack+116539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+148834,int_stack+750,int_stack+144679,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+145279,int_stack+48473,int_stack+48053, 0.0, zero_stack, 1.0, int_stack+50694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+150634,int_stack+1764,int_stack+48473, 0.0, zero_stack, 1.0, int_stack+2205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+151579,int_stack+150634,int_stack+145279, 0.0, zero_stack, 1.0, int_stack+118039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+152839,int_stack+151579,int_stack+750,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+150634,int_stack+101615,int_stack+101515, 1.0, int_stack+112213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+150934,int_stack+101765,int_stack+101615, 1.0, int_stack+112313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+151384,int_stack+150934,int_stack+150634, 1.0, int_stack+115939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+150934,int_stack+102960,int_stack+102585, 1.0, int_stack+102435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+151984,int_stack+103185,int_stack+102960, 1.0, int_stack+102735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+155539,int_stack+151984,int_stack+150934, 1.0, int_stack+116539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+156439,int_stack+155539,int_stack+151384,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+151984,int_stack+51324,int_stack+50904, 1.0, int_stack+50694, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+158239,int_stack+2520,int_stack+51324, 1.0, int_stack+2205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1650,int_stack+158239,int_stack+151984, 1.0, int_stack+118039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+158239,int_stack+1650,int_stack+155539,60);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1650,int_stack+112463,int_stack+112313,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+1650,int_stack+115939,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+118039,int_stack+3843,int_stack+102735,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+118039,int_stack+116539,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+161839,int_stack+160939,int_stack+2100,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+115939,int_stack+105781,int_stack+105681,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+116539,int_stack+105931,int_stack+105781,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+118039,int_stack+116539,int_stack+115939,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+116539,int_stack+106751,int_stack+106601,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+106976,int_stack+106751,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+164314,int_stack+163639,int_stack+116539,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+165214,int_stack+164314,int_stack+118039, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+163639,int_stack+64545,int_stack+64125,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+167014,int_stack+2961,int_stack+64545,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+167959,int_stack+167014,int_stack+163639,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+169219,int_stack+167959,int_stack+164314, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+109047,int_stack+108947,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1650,int_stack+109197,int_stack+109047,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+167314,int_stack+1650,int_stack+167014,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1650,int_stack+110017,int_stack+109867,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+167914,int_stack+110242,int_stack+110017,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+171919,int_stack+167914,int_stack+1650,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+172819,int_stack+171919,int_stack+167314, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+79882,int_stack+79462,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+174619,int_stack+3402,int_stack+79882,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2700,int_stack+174619,int_stack+167914,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+174619,int_stack+2700,int_stack+171919, 0.0, zero_stack, 1.0, int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+112773,int_stack+112673,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3000,int_stack+112923,int_stack+112773,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3450,int_stack+3000,int_stack+2700,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3000,int_stack+113743,int_stack+113593,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+168544,int_stack+113968,int_stack+113743,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+177319,int_stack+168544,int_stack+3000,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+178219,int_stack+177319,int_stack+3450, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+168544,int_stack+104610,int_stack+104190,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+180019,int_stack+4158,int_stack+104610,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+180964,int_stack+180019,int_stack+168544,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+182224,int_stack+180964,int_stack+177319, 1.0, int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+4699,int_stack+4599, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+85185,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+4849,int_stack+4699, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+85285,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+116239,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+5209,int_stack+5059, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+86105,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180019,int_stack+5434,int_stack+5209, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+86255,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+180694,int_stack+180019,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+117589,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+184924,int_stack+180694,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+181594,int_stack+5959,int_stack+5749, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+37489,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4050,int_stack+6274,int_stack+5959, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+37909,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4995,int_stack+4050,int_stack+181594, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+121414,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+186724,int_stack+4995,int_stack+180694,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4050,int_stack+6815,int_stack+6715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85185, 1.0, int_stack+88451,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4350,int_stack+6965,int_stack+6815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85285, 1.0, int_stack+88551,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+4350,int_stack+4050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116239, 1.0, int_stack+122044,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4050,int_stack+7325,int_stack+7175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86105, 1.0, int_stack+89371,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4500,int_stack+7550,int_stack+7325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86255, 1.0, int_stack+89521,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+4500,int_stack+4050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+117589, 1.0, int_stack+122344,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+4050,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5850,int_stack+8075,int_stack+7865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37489, 1.0, int_stack+40130,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6480,int_stack+8390,int_stack+8075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37909, 1.0, int_stack+40550,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7425,int_stack+6480,int_stack+5850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121414, 1.0, int_stack+123394,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+189424,int_stack+7425,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+8931,int_stack+8831, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+88451, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+9081,int_stack+8931, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+88551, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+122044, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+9441,int_stack+9291, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+89371, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5850,int_stack+9666,int_stack+9441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+89521, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6525,int_stack+5850,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+122344, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+7425,int_stack+6525,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+10191,int_stack+9981, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40130, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+10506,int_stack+10191, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40550, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+180019,int_stack+9225,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+123394, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+192124,int_stack+180019,int_stack+6525,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180019,int_stack+11047,int_stack+10947, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85185, 0.0, zero_stack, 1.0, int_stack+91717,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180319,int_stack+11197,int_stack+11047, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85285, 0.0, zero_stack, 1.0, int_stack+91817,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+180319,int_stack+180019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116239, 0.0, zero_stack, 1.0, int_stack+129649,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180019,int_stack+11557,int_stack+11407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86105, 0.0, zero_stack, 1.0, int_stack+92637,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180469,int_stack+11782,int_stack+11557, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86255, 0.0, zero_stack, 1.0, int_stack+92787,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+180469,int_stack+180019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+117589, 0.0, zero_stack, 1.0, int_stack+129949,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+180019,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+12307,int_stack+12097, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37489, 0.0, zero_stack, 1.0, int_stack+42771,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9855,int_stack+12622,int_stack+12307, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37909, 0.0, zero_stack, 1.0, int_stack+43191,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10800,int_stack+9855,int_stack+9225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121414, 0.0, zero_stack, 1.0, int_stack+130999,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+194824,int_stack+10800,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+13163,int_stack+13063, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88451, 1.0, int_stack+91717, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+13313,int_stack+13163, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88551, 1.0, int_stack+91817, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122044, 1.0, int_stack+129649, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+13673,int_stack+13523, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89371, 1.0, int_stack+92637, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+13898,int_stack+13673, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89521, 1.0, int_stack+92787, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9900,int_stack+9225,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122344, 1.0, int_stack+129949, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+10800,int_stack+9900,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+14423,int_stack+14213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40130, 1.0, int_stack+42771, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+14738,int_stack+14423, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40550, 1.0, int_stack+43191, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13545,int_stack+12600,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123394, 1.0, int_stack+130999, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+197524,int_stack+13545,int_stack+9900,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+15279,int_stack+15179, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+91717, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+15429,int_stack+15279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+91817, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+129649, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+15789,int_stack+15639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+92637, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+16014,int_stack+15789, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+92787, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13275,int_stack+12600,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+129949, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+14175,int_stack+13275,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+16539,int_stack+16329, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+42771, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+16854,int_stack+16539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43191, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15975,int_stack+9225,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+130999, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+200224,int_stack+15975,int_stack+13275,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15975,int_stack+17395,int_stack+17295, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85185, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+94983,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16275,int_stack+17545,int_stack+17395, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+85285, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95083,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+16275,int_stack+15975, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116239, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15975,int_stack+17905,int_stack+17755, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86105, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+95903,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16425,int_stack+18130,int_stack+17905, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86255, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96053,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+16425,int_stack+15975, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15975,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17775,int_stack+18655,int_stack+18445, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37489, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45412,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+18970,int_stack+18655, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37909, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45832,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+9225,int_stack+17775, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121414, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140554,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+202924,int_stack+12600,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+19511,int_stack+19411, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88451, 0.0, zero_stack, 1.0, int_stack+94983, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+19661,int_stack+19511, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+88551, 0.0, zero_stack, 1.0, int_stack+95083, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122044, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+20021,int_stack+19871, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89371, 0.0, zero_stack, 1.0, int_stack+95903, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+20246,int_stack+20021, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89521, 0.0, zero_stack, 1.0, int_stack+96053, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13275,int_stack+12600,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122344, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+17775,int_stack+13275,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+20771,int_stack+20561, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40130, 0.0, zero_stack, 1.0, int_stack+45412, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19575,int_stack+21086,int_stack+20771, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40550, 0.0, zero_stack, 1.0, int_stack+45832, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+19575,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123394, 0.0, zero_stack, 1.0, int_stack+140554, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+205624,int_stack+9225,int_stack+13275,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+21627,int_stack+21527, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91717, 1.0, int_stack+94983, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9525,int_stack+21777,int_stack+21627, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91817, 1.0, int_stack+95083, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+9525,int_stack+9225, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129649, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+22137,int_stack+21987, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92637, 1.0, int_stack+95903, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9675,int_stack+22362,int_stack+22137, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92787, 1.0, int_stack+96053, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+9675,int_stack+9225, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129949, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+19575,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+22887,int_stack+22677, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42771, 1.0, int_stack+45412, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9855,int_stack+23202,int_stack+22887, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43191, 1.0, int_stack+45832, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21375,int_stack+9855,int_stack+9225, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130999, 1.0, int_stack+140554, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+208324,int_stack+21375,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+23743,int_stack+23643, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+94983, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+23893,int_stack+23743, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+95083, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+24253,int_stack+24103, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+95903, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21375,int_stack+24478,int_stack+24253, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+96053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22050,int_stack+21375,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+22950,int_stack+22050,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+25003,int_stack+24793, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+45412, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+25318,int_stack+25003, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+45832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+9225,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+140554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+211024,int_stack+12600,int_stack+22050,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12600,int_stack+25859,int_stack+25759, 0.0, zero_stack, 1.0, int_stack+85185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98249,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12900,int_stack+26009,int_stack+25859, 0.0, zero_stack, 1.0, int_stack+85285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98349,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+12900,int_stack+12600, 0.0, zero_stack, 1.0, int_stack+116239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143929,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12600,int_stack+26369,int_stack+26219, 0.0, zero_stack, 1.0, int_stack+86105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99169,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13050,int_stack+26594,int_stack+26369, 0.0, zero_stack, 1.0, int_stack+86255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99319,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+13050,int_stack+12600, 0.0, zero_stack, 1.0, int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144229,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+24750,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12600,int_stack+27119,int_stack+26909, 0.0, zero_stack, 1.0, int_stack+37489, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48053,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13230,int_stack+27434,int_stack+27119, 0.0, zero_stack, 1.0, int_stack+37909, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48473,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+13230,int_stack+12600, 0.0, zero_stack, 1.0, int_stack+121414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145279,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+213724,int_stack+9225,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+27975,int_stack+27875, 0.0, zero_stack, 1.0, int_stack+88451, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98249, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+28125,int_stack+27975, 0.0, zero_stack, 1.0, int_stack+88551, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98349, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+122044, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+143929, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+28485,int_stack+28335, 0.0, zero_stack, 1.0, int_stack+89371, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99169, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+28710,int_stack+28485, 0.0, zero_stack, 1.0, int_stack+89521, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99319, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9900,int_stack+9225,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+122344, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144229, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+26550,int_stack+9900,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+29235,int_stack+29025, 0.0, zero_stack, 1.0, int_stack+40130, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48053, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+29550,int_stack+29235, 0.0, zero_stack, 1.0, int_stack+40550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48473, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21375,int_stack+12600,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+123394, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145279, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+216424,int_stack+21375,int_stack+9900,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21375,int_stack+30091,int_stack+29991, 0.0, zero_stack, 1.0, int_stack+91717, 0.0, zero_stack, 1.0, int_stack+98249, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21675,int_stack+30241,int_stack+30091, 0.0, zero_stack, 1.0, int_stack+91817, 0.0, zero_stack, 1.0, int_stack+98349, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+21675,int_stack+21375, 0.0, zero_stack, 1.0, int_stack+129649, 0.0, zero_stack, 1.0, int_stack+143929, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21375,int_stack+30601,int_stack+30451, 0.0, zero_stack, 1.0, int_stack+92637, 0.0, zero_stack, 1.0, int_stack+99169, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21825,int_stack+30826,int_stack+30601, 0.0, zero_stack, 1.0, int_stack+92787, 0.0, zero_stack, 1.0, int_stack+99319, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+21825,int_stack+21375, 0.0, zero_stack, 1.0, int_stack+129949, 0.0, zero_stack, 1.0, int_stack+144229, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+28350,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21375,int_stack+31351,int_stack+31141, 0.0, zero_stack, 1.0, int_stack+42771, 0.0, zero_stack, 1.0, int_stack+48053, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+22005,int_stack+31666,int_stack+31351, 0.0, zero_stack, 1.0, int_stack+43191, 0.0, zero_stack, 1.0, int_stack+48473, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+22005,int_stack+21375, 0.0, zero_stack, 1.0, int_stack+130999, 0.0, zero_stack, 1.0, int_stack+145279, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+219124,int_stack+12600,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+32207,int_stack+32107, 0.0, zero_stack, 1.0, int_stack+94983, 1.0, int_stack+98249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+32357,int_stack+32207, 0.0, zero_stack, 1.0, int_stack+95083, 1.0, int_stack+98349, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+143929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+32717,int_stack+32567, 0.0, zero_stack, 1.0, int_stack+95903, 1.0, int_stack+99169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+32942,int_stack+32717, 0.0, zero_stack, 1.0, int_stack+96053, 1.0, int_stack+99319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13275,int_stack+12600,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+300, 1.0, int_stack+144229, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+30150,int_stack+13275,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+33467,int_stack+33257, 0.0, zero_stack, 1.0, int_stack+45412, 1.0, int_stack+48053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21375,int_stack+33782,int_stack+33467, 0.0, zero_stack, 1.0, int_stack+45832, 1.0, int_stack+48473, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+21375,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+140554, 1.0, int_stack+145279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+221824,int_stack+9225,int_stack+13275,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+34323,int_stack+34223, 0.0, zero_stack, 2.0, int_stack+98249, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9525,int_stack+34473,int_stack+34323, 0.0, zero_stack, 2.0, int_stack+98349, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+9525,int_stack+9225, 0.0, zero_stack, 2.0, int_stack+143929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+34833,int_stack+34683, 0.0, zero_stack, 2.0, int_stack+99169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9675,int_stack+35058,int_stack+34833, 0.0, zero_stack, 2.0, int_stack+99319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+9675,int_stack+9225, 0.0, zero_stack, 2.0, int_stack+144229, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+31950,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+35583,int_stack+35373, 0.0, zero_stack, 2.0, int_stack+48053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9855,int_stack+35898,int_stack+35583, 0.0, zero_stack, 2.0, int_stack+48473, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21375,int_stack+9855,int_stack+9225, 0.0, zero_stack, 2.0, int_stack+145279, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+224524,int_stack+21375,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+36439,int_stack+36339, 1.0, int_stack+85185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101515,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+161239,int_stack+36589,int_stack+36439, 1.0, int_stack+85285, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101615,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+161239,int_stack+160939, 1.0, int_stack+116239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150634,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+36949,int_stack+36799, 1.0, int_stack+86105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102585,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+21375,int_stack+37174,int_stack+36949, 1.0, int_stack+86255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102960,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22050,int_stack+21375,int_stack+160939, 1.0, int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150934,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+33750,int_stack+22050,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+38224,int_stack+37699, 1.0, int_stack+37489, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50904,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+38539,int_stack+38224, 1.0, int_stack+37909, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51324,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+9225,int_stack+160939, 1.0, int_stack+121414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151984,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+35550,int_stack+12600,int_stack+22050,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+39080,int_stack+38980, 1.0, int_stack+88451, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101515, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+39230,int_stack+39080, 1.0, int_stack+88551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101615, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 1.0, int_stack+122044, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150634, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+39590,int_stack+39440, 1.0, int_stack+89371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102585, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+39815,int_stack+39590, 1.0, int_stack+89521, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102960, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13275,int_stack+12600,int_stack+117589, 1.0, int_stack+122344, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150934, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+38250,int_stack+13275,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12600,int_stack+40865,int_stack+40340, 1.0, int_stack+40130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50904, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+41180,int_stack+40865, 1.0, int_stack+40550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51324, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+121414,int_stack+12600, 1.0, int_stack+123394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151984, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+227224,int_stack+9225,int_stack+13275,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+41721,int_stack+41621, 1.0, int_stack+91717, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101515, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+41871,int_stack+41721, 1.0, int_stack+91817, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101615, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 1.0, int_stack+129649, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150634, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+42231,int_stack+42081, 1.0, int_stack+92637, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102585, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+42456,int_stack+42231, 1.0, int_stack+92787, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102960, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9900,int_stack+9225,int_stack+117589, 1.0, int_stack+129949, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150934, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+40050,int_stack+9900,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+43506,int_stack+42981, 1.0, int_stack+42771, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50904, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+43821,int_stack+43506, 1.0, int_stack+43191, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51324, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+12600,int_stack+9225, 1.0, int_stack+130999, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151984, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+229924,int_stack+121414,int_stack+9900,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+44362,int_stack+44262, 1.0, int_stack+94983, 0.0, zero_stack, 1.0, int_stack+101515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+44512,int_stack+44362, 1.0, int_stack+95083, 0.0, zero_stack, 1.0, int_stack+101615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+150634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+44872,int_stack+44722, 1.0, int_stack+95903, 0.0, zero_stack, 1.0, int_stack+102585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+45097,int_stack+44872, 1.0, int_stack+96053, 0.0, zero_stack, 1.0, int_stack+102960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+121414,int_stack+117589, 1.0, int_stack+300, 0.0, zero_stack, 1.0, int_stack+150934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+41850,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+121414,int_stack+46147,int_stack+45622, 1.0, int_stack+45412, 0.0, zero_stack, 1.0, int_stack+50904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+46462,int_stack+46147, 1.0, int_stack+45832, 0.0, zero_stack, 1.0, int_stack+51324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+9225,int_stack+121414, 1.0, int_stack+140554, 0.0, zero_stack, 1.0, int_stack+151984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+43650,int_stack+12600,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+47003,int_stack+46903, 1.0, int_stack+98249, 1.0, int_stack+101515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+47153,int_stack+47003, 1.0, int_stack+98349, 1.0, int_stack+101615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 1.0, int_stack+143929, 1.0, int_stack+150634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+47513,int_stack+47363, 1.0, int_stack+99169, 1.0, int_stack+102585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+47738,int_stack+47513, 1.0, int_stack+99319, 1.0, int_stack+102960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+140554,int_stack+117589, 1.0, int_stack+144229, 1.0, int_stack+150934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+232624,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+48788,int_stack+48263, 1.0, int_stack+48053, 1.0, int_stack+50904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+49103,int_stack+48788, 1.0, int_stack+48473, 1.0, int_stack+51324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+12600,int_stack+140554, 1.0, int_stack+145279, 1.0, int_stack+151984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+46350,int_stack+121414,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+49644,int_stack+49544, 2.0, int_stack+101515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+49794,int_stack+49644, 2.0, int_stack+101615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 2.0, int_stack+150634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+50154,int_stack+50004, 2.0, int_stack+102585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+50379,int_stack+50154, 2.0, int_stack+102960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+140554,int_stack+117589, 2.0, int_stack+150934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+49050,int_stack+160939,int_stack+2100,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+51639,int_stack+51114, 2.0, int_stack+50904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+51954,int_stack+51639, 2.0, int_stack+51324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+121414,int_stack+140554, 2.0, int_stack+151984, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+234424,int_stack+12600,int_stack+160939,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+52495,int_stack+52395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105681,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+52645,int_stack+52495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105781,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115939,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+53005,int_stack+52855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106601,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+53230,int_stack+53005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106751,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+140554,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116539,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+50850,int_stack+160939,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116989, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+53755,int_stack+53545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64125,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+54070,int_stack+53755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64545,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+12600,int_stack+140554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163639,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+237124,int_stack+121414,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+54611,int_stack+54511, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105681, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+54761,int_stack+54611, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105781, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115939, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+55121,int_stack+54971, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106601, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+55346,int_stack+55121, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106751, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+140554,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116539, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+52650,int_stack+160939,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+55871,int_stack+55661, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64125, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+56186,int_stack+55871, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64545, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+121414,int_stack+140554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163639, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+239824,int_stack+12600,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+56727,int_stack+56627, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105681, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+56877,int_stack+56727, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105781, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115939, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+57237,int_stack+57087, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106601, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+57462,int_stack+57237, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106751, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+140554,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116539, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+54450,int_stack+160939,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130399, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+57987,int_stack+57777, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64125, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+58302,int_stack+57987, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64545, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+12600,int_stack+140554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163639, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+242524,int_stack+121414,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+134554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+58843,int_stack+58743, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+58993,int_stack+58843, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+59353,int_stack+59203, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+59578,int_stack+59353, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106751, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+140554,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+56250,int_stack+160939,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139954, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+60103,int_stack+59893, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+60418,int_stack+60103, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+121414,int_stack+140554, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+163639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+58050,int_stack+12600,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141229, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+60959,int_stack+60859, 0.0, zero_stack, 1.0, int_stack+105681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+61109,int_stack+60959, 0.0, zero_stack, 1.0, int_stack+105781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 0.0, zero_stack, 1.0, int_stack+115939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+61469,int_stack+61319, 0.0, zero_stack, 1.0, int_stack+106601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+61694,int_stack+61469, 0.0, zero_stack, 1.0, int_stack+106751, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+140554,int_stack+117589, 0.0, zero_stack, 1.0, int_stack+116539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+245224,int_stack+160939,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144679, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+62219,int_stack+62009, 0.0, zero_stack, 1.0, int_stack+64125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+62534,int_stack+62219, 0.0, zero_stack, 1.0, int_stack+64545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+12600,int_stack+140554, 0.0, zero_stack, 1.0, int_stack+163639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+247024,int_stack+121414,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116239,int_stack+63075,int_stack+62975, 1.0, int_stack+105681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+63225,int_stack+63075, 1.0, int_stack+105781, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+116239, 1.0, int_stack+115939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+63585,int_stack+63435, 1.0, int_stack+106601, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140554,int_stack+63810,int_stack+63585, 1.0, int_stack+106751, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+140554,int_stack+117589, 1.0, int_stack+116539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+60750,int_stack+160939,int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+140554,int_stack+64860,int_stack+64335, 1.0, int_stack+64125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+115939,int_stack+65175,int_stack+64860, 1.0, int_stack+64545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+115939,int_stack+140554, 1.0, int_stack+163639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+62550,int_stack+121414,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+65716,int_stack+65616,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+65866,int_stack+65716,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+161239,int_stack+117589,int_stack+160939,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+66226,int_stack+66076,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+66451,int_stack+66226,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+163639,int_stack+117589,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+249724,int_stack+121414,int_stack+161239, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+118039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+163639,int_stack+66976,int_stack+66766,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+115939,int_stack+67291,int_stack+66976,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+115939,int_stack+163639,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+251524,int_stack+12600,int_stack+121414, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+164314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+121414,int_stack+67832,int_stack+67732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108947,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+67982,int_stack+67832, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109047,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+121414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167014,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+68342,int_stack+68192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109867,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+68567,int_stack+68342, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110017,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1650,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+65250,int_stack+160939,int_stack+2100, 0.0, zero_stack, 1.0, int_stack+116989, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+163639,int_stack+69092,int_stack+68882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79462,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+69407,int_stack+69092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79882,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+121414,int_stack+163639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167914,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+67050,int_stack+12600,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+118714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+69948,int_stack+69848, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108947, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+70098,int_stack+69948, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109047, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+161239,int_stack+117589,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167014, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+70458,int_stack+70308, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109867, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+70683,int_stack+70458, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110017, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+163639,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1650, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+254224,int_stack+12600,int_stack+161239, 0.0, zero_stack, 1.0, int_stack+122794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+163639,int_stack+71208,int_stack+70998, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79462, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+71523,int_stack+71208, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79882, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+121414,int_stack+163639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167914, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+256024,int_stack+9225,int_stack+12600, 0.0, zero_stack, 1.0, int_stack+126949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12600,int_stack+72064,int_stack+71964, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108947, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+72214,int_stack+72064, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109047, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+12600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167014, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+72574,int_stack+72424, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109867, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+72799,int_stack+72574, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110017, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1650, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+69750,int_stack+160939,int_stack+2100, 0.0, zero_stack, 1.0, int_stack+130399, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+163639,int_stack+73324,int_stack+73114, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79462, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+73639,int_stack+73324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79882, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+12600,int_stack+163639, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167914, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+258724,int_stack+9225,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+134554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+74180,int_stack+74080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108947, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+74330,int_stack+74180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109047, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+161239,int_stack+117589,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+74690,int_stack+74540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+74915,int_stack+74690, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110017, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+163639,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+71550,int_stack+9225,int_stack+161239, 0.0, zero_stack, 1.0, int_stack+139954, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+163639,int_stack+75440,int_stack+75230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+75755,int_stack+75440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+12600,int_stack+163639, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+73350,int_stack+121414,int_stack+9225, 0.0, zero_stack, 1.0, int_stack+141229, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9225,int_stack+76296,int_stack+76196, 0.0, zero_stack, 1.0, int_stack+108947, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+76446,int_stack+76296, 0.0, zero_stack, 1.0, int_stack+109047, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2100,int_stack+117589,int_stack+9225, 0.0, zero_stack, 1.0, int_stack+167014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+76806,int_stack+76656, 0.0, zero_stack, 1.0, int_stack+109867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+77031,int_stack+76806, 0.0, zero_stack, 1.0, int_stack+110017, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+117589, 0.0, zero_stack, 1.0, int_stack+1650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+261424,int_stack+160939,int_stack+2100, 0.0, zero_stack, 1.0, int_stack+144679, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+163639,int_stack+77556,int_stack+77346, 0.0, zero_stack, 1.0, int_stack+79462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+77871,int_stack+77556, 0.0, zero_stack, 1.0, int_stack+79882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+9225,int_stack+163639, 0.0, zero_stack, 1.0, int_stack+167914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+263224,int_stack+121414,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+78412,int_stack+78312, 1.0, int_stack+108947, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+78562,int_stack+78412, 1.0, int_stack+109047, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+161239,int_stack+117589,int_stack+160939, 1.0, int_stack+167014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+78922,int_stack+78772, 1.0, int_stack+109867, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+79147,int_stack+78922, 1.0, int_stack+110017, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+163639,int_stack+117589, 1.0, int_stack+1650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+76050,int_stack+121414,int_stack+161239, 0.0, zero_stack, 1.0, int_stack+151384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1650,int_stack+80197,int_stack+79672, 1.0, int_stack+79462, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+80512,int_stack+80197, 1.0, int_stack+79882, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+9225,int_stack+1650, 1.0, int_stack+167914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+77850,int_stack+12600,int_stack+121414, 0.0, zero_stack, 1.0, int_stack+155539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+81053,int_stack+80953,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+81203,int_stack+81053,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+117589,int_stack+167014,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+81563,int_stack+81413,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+81788,int_stack+81563,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+117589,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+265924,int_stack+160939,int_stack+121414, 0.0, zero_stack, 1.0, int_stack+118039, 1.0, int_stack+167314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+82313,int_stack+82103,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+82628,int_stack+82313,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+121414,int_stack+167914,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+267724,int_stack+12600,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+164314, 1.0, int_stack+171919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+83169,int_stack+83069,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+83319,int_stack+83169,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+117589,int_stack+167014,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+83679,int_stack+83529,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+83904,int_stack+83679,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+163639,int_stack+117589,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+80550,int_stack+12600,int_stack+160939, 0.0, zero_stack, 2.0, int_stack+167314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+84429,int_stack+84219,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+121414,int_stack+84744,int_stack+84429,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+121414,int_stack+167914,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+82350,int_stack+9225,int_stack+12600, 0.0, zero_stack, 2.0, int_stack+171919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+85745,int_stack+85645, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112673,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+117589,int_stack+85895,int_stack+85745, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112773,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+117589,int_stack+167014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+117589,int_stack+86945,int_stack+86795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113593,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+87170,int_stack+86945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113743,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+117589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+85050,int_stack+160939,int_stack+12600, 1.0, int_stack+116989, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+87695,int_stack+87485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104190,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+88010,int_stack+87695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104610,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+12600,int_stack+167914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168544,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+270424,int_stack+9225,int_stack+160939, 1.0, int_stack+118714, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+89011,int_stack+88911, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112673, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+160939,int_stack+89161,int_stack+89011, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112773, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+160939,int_stack+167014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+90211,int_stack+90061, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113593, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+90436,int_stack+90211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113743, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9825,int_stack+163639,int_stack+160939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+115939,int_stack+9825,int_stack+9225, 1.0, int_stack+122794, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+90961,int_stack+90751, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104190, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+91276,int_stack+90961, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104610, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+12600,int_stack+167914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168544, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+86850,int_stack+121414,int_stack+9825, 1.0, int_stack+126949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+92277,int_stack+92177, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112673, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+126949,int_stack+92427,int_stack+92277, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112773, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+126949,int_stack+167014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+126949,int_stack+93477,int_stack+93327, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113593, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+93702,int_stack+93477, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113743, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+126949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+122014,int_stack+160939,int_stack+121414, 1.0, int_stack+130399, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+94227,int_stack+94017, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104190, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+94542,int_stack+94227, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104610, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+12600,int_stack+167914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168544, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+89550,int_stack+9225,int_stack+160939, 1.0, int_stack+134554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+95543,int_stack+95443, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112673, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+134554,int_stack+95693,int_stack+95543, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+112773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+134554,int_stack+167014, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+134554,int_stack+96743,int_stack+96593, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+96968,int_stack+96743, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113743, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+134554, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+129649,int_stack+160939,int_stack+121414, 1.0, int_stack+139954, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+97493,int_stack+97283, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+139954,int_stack+97808,int_stack+97493, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+139954,int_stack+167914, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+92250,int_stack+9225,int_stack+160939, 1.0, int_stack+141229, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+98809,int_stack+98709, 0.0, zero_stack, 1.0, int_stack+112673, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+160939,int_stack+98959,int_stack+98809, 0.0, zero_stack, 1.0, int_stack+112773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+160939,int_stack+167014, 0.0, zero_stack, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+100009,int_stack+99859, 0.0, zero_stack, 1.0, int_stack+113593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+100234,int_stack+100009, 0.0, zero_stack, 1.0, int_stack+113743, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+134554,int_stack+163639,int_stack+160939, 0.0, zero_stack, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+139954,int_stack+134554,int_stack+121414, 1.0, int_stack+144679, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+100759,int_stack+100549, 0.0, zero_stack, 1.0, int_stack+104190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+101074,int_stack+100759, 0.0, zero_stack, 1.0, int_stack+104610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12600,int_stack+9225,int_stack+167914, 0.0, zero_stack, 1.0, int_stack+168544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+94950,int_stack+12600,int_stack+134554, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+102075,int_stack+101975, 1.0, int_stack+112673, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+134554,int_stack+102225,int_stack+102075, 1.0, int_stack+112773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+134554,int_stack+167014, 1.0, int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+134554,int_stack+103650,int_stack+103500, 1.0, int_stack+113593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+103875,int_stack+103650, 1.0, int_stack+113743, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+134554, 1.0, int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+143929,int_stack+160939,int_stack+121414, 1.0, int_stack+151384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+167914,int_stack+104925,int_stack+104400, 1.0, int_stack+104190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12600,int_stack+105240,int_stack+104925, 1.0, int_stack+104610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (h0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9225,int_stack+12600,int_stack+167914, 1.0, int_stack+168544, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+9225,int_stack+160939, 1.0, int_stack+155539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+106241,int_stack+106141,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+155539,int_stack+106391,int_stack+106241,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+155539,int_stack+167014,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+155539,int_stack+107441,int_stack+107291,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+163639,int_stack+107666,int_stack+107441,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+163639,int_stack+155539,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+150634,int_stack+160939,int_stack+121414, 1.0, int_stack+118039, 0.0, zero_stack, 1.0, int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+155539,int_stack+108191,int_stack+107981,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+108506,int_stack+108191,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+167914,int_stack+9225,int_stack+155539,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+97650,int_stack+167914,int_stack+160939, 1.0, int_stack+164314, 0.0, zero_stack, 1.0, int_stack+177319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+167014,int_stack+109507,int_stack+109407,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+160939,int_stack+109657,int_stack+109507,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+160939,int_stack+167014,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+110707,int_stack+110557,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+167914,int_stack+110932,int_stack+110707,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+155539,int_stack+167914,int_stack+160939,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+117739,int_stack+155539,int_stack+121414, 1.0, int_stack+167314, 1.0, int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+160939,int_stack+111457,int_stack+111247,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+167014,int_stack+111772,int_stack+111457,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+167959,int_stack+167014,int_stack+160939,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+100350,int_stack+167959,int_stack+155539, 1.0, int_stack+171919, 1.0, int_stack+177319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+171919,int_stack+113233,int_stack+113133,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+172219,int_stack+113383,int_stack+113233,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+121414,int_stack+172219,int_stack+171919,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+171919,int_stack+114433,int_stack+114283,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+155539,int_stack+114658,int_stack+114433,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+160939,int_stack+155539,int_stack+171919,15);
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+167014,int_stack+160939,int_stack+121414, 2.0, int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (h0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+171919,int_stack+115183,int_stack+114973,21);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9225,int_stack+115498,int_stack+115183,21);
 /*--- compute (h0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+163639,int_stack+9225,int_stack+171919,21);
 /*--- compute (gp|fd) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+103050,int_stack+163639,int_stack+160939, 2.0, int_stack+177319, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+105750,int_stack+124249,int_stack+119614,60);
     Libderiv->ABCD[11] = int_stack + 105750;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+123814,int_stack+131854,int_stack+127849,60);
     Libderiv->ABCD[10] = int_stack + 123814;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+131449,int_stack+137254,int_stack+135454,60);
     Libderiv->ABCD[9] = int_stack + 131449;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+109350,int_stack+146134,int_stack+142129,60);
     Libderiv->ABCD[8] = int_stack + 109350;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+273124,int_stack+152839,int_stack+148834,60);
     Libderiv->ABCD[7] = int_stack + 273124;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+152434,int_stack+158239,int_stack+156439,60);
     Libderiv->ABCD[6] = int_stack + 152434;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+158239,int_stack+169219,int_stack+165214, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+161839, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 158239;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+168814,int_stack+174619,int_stack+172819, 0.0, zero_stack, 1.0, int_stack+161839, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 168814;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+174619,int_stack+182224,int_stack+178219, 1.0, int_stack+161839, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 174619;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+276724,int_stack+186724,int_stack+184924,60);
     Libderiv->ABCD[155] = int_stack + 276724;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+181819,int_stack+189424,int_stack+4050,60);
     Libderiv->ABCD[143] = int_stack + 181819;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+2700,int_stack+192124,int_stack+7425,60);
     Libderiv->ABCD[142] = int_stack + 2700;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+6300,int_stack+194824,int_stack+180019,60);
     Libderiv->ABCD[131] = int_stack + 6300;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+185419,int_stack+197524,int_stack+10800,60);
     Libderiv->ABCD[130] = int_stack + 185419;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+9900,int_stack+200224,int_stack+14175,60);
     Libderiv->ABCD[129] = int_stack + 9900;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+189019,int_stack+202924,int_stack+15975,60);
     Libderiv->ABCD[119] = int_stack + 189019;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+13500,int_stack+205624,int_stack+17775,60);
     Libderiv->ABCD[118] = int_stack + 13500;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+192619,int_stack+208324,int_stack+19575,60);
     Libderiv->ABCD[117] = int_stack + 192619;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+17100,int_stack+211024,int_stack+22950,60);
     Libderiv->ABCD[116] = int_stack + 17100;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+20700,int_stack+213724,int_stack+24750,60);
     Libderiv->ABCD[107] = int_stack + 20700;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+196219,int_stack+216424,int_stack+26550,60);
     Libderiv->ABCD[106] = int_stack + 196219;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+24300,int_stack+219124,int_stack+28350,60);
     Libderiv->ABCD[105] = int_stack + 24300;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+199819,int_stack+221824,int_stack+30150,60);
     Libderiv->ABCD[104] = int_stack + 199819;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+27900,int_stack+224524,int_stack+31950,60);
     Libderiv->ABCD[103] = int_stack + 27900;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+203419,int_stack+35550,int_stack+33750,60);
     Libderiv->ABCD[95] = int_stack + 203419;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+31500,int_stack+227224,int_stack+38250,60);
     Libderiv->ABCD[94] = int_stack + 31500;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+35100,int_stack+229924,int_stack+40050,60);
     Libderiv->ABCD[93] = int_stack + 35100;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+207019,int_stack+43650,int_stack+41850,60);
     Libderiv->ABCD[92] = int_stack + 207019;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+38700,int_stack+46350,int_stack+232624,60);
     Libderiv->ABCD[91] = int_stack + 38700;
 /*--- compute (fd|fd) ---*/
   hrr1_build_fd(Libderiv->AB,int_stack+42300,int_stack+234424,int_stack+49050,60);
     Libderiv->ABCD[90] = int_stack + 42300;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+45900,int_stack+237124,int_stack+50850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[47] = int_stack + 45900;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+210619,int_stack+239824,int_stack+52650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127849, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[46] = int_stack + 210619;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+49500,int_stack+242524,int_stack+54450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[45] = int_stack + 49500;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+214219,int_stack+58050,int_stack+56250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+142129, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[44] = int_stack + 214219;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+53100,int_stack+247024,int_stack+245224, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+148834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[43] = int_stack + 53100;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+56700,int_stack+62550,int_stack+60750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+156439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[42] = int_stack + 56700;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+60300,int_stack+251524,int_stack+249724, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+165214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[38] = int_stack + 60300;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+217819,int_stack+67050,int_stack+65250, 0.0, zero_stack, 1.0, int_stack+119614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[35] = int_stack + 217819;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+63900,int_stack+256024,int_stack+254224, 0.0, zero_stack, 1.0, int_stack+127849, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[34] = int_stack + 63900;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+221419,int_stack+258724,int_stack+69750, 0.0, zero_stack, 1.0, int_stack+135454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[33] = int_stack + 221419;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+67500,int_stack+73350,int_stack+71550, 0.0, zero_stack, 1.0, int_stack+142129, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[32] = int_stack + 67500;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+71100,int_stack+263224,int_stack+261424, 0.0, zero_stack, 1.0, int_stack+148834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[31] = int_stack + 71100;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+225019,int_stack+77850,int_stack+76050, 0.0, zero_stack, 1.0, int_stack+156439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[30] = int_stack + 225019;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+74700,int_stack+267724,int_stack+265924, 0.0, zero_stack, 1.0, int_stack+165214, 1.0, int_stack+172819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[26] = int_stack + 74700;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+228619,int_stack+82350,int_stack+80550, 0.0, zero_stack, 2.0, int_stack+172819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[25] = int_stack + 228619;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+78300,int_stack+270424,int_stack+85050, 1.0, int_stack+119614, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[23] = int_stack + 78300;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+81900,int_stack+86850,int_stack+115939, 1.0, int_stack+127849, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[22] = int_stack + 81900;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+85500,int_stack+89550,int_stack+122014, 1.0, int_stack+135454, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[21] = int_stack + 85500;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+119539,int_stack+92250,int_stack+129649, 1.0, int_stack+142129, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[20] = int_stack + 119539;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+127414,int_stack+94950,int_stack+139954, 1.0, int_stack+148834, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[19] = int_stack + 127414;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+145729,int_stack+0,int_stack+143929, 1.0, int_stack+156439, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[18] = int_stack + 145729;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+89100,int_stack+97650,int_stack+150634, 1.0, int_stack+165214, 0.0, zero_stack, 1.0, int_stack+178219, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[14] = int_stack + 89100;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+161839,int_stack+100350,int_stack+117739, 1.0, int_stack+172819, 1.0, int_stack+178219, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[13] = int_stack + 161839;
 /*--- compute (fd|fd) ---*/
   d1hrr1_build_fd(Libderiv->AB,int_stack+92700,int_stack+103050,int_stack+167014, 2.0, int_stack+178219, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[12] = int_stack + 92700;

}
