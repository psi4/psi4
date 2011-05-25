#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ddff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|ff) integrals */

void d12hrr_order_ddff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][6][11] = int_stack + 0;
 Libderiv->deriv_classes[4][6][10] = int_stack + 420;
 Libderiv->deriv_classes[4][6][9] = int_stack + 840;
 Libderiv->deriv_classes[4][6][8] = int_stack + 1260;
 Libderiv->deriv_classes[4][6][7] = int_stack + 1680;
 Libderiv->dvrr_classes[4][5] = int_stack + 2100;
 Libderiv->deriv_classes[4][6][6] = int_stack + 2415;
 Libderiv->deriv_classes[4][6][2] = int_stack + 2835;
 Libderiv->deriv_classes[4][6][1] = int_stack + 3255;
 Libderiv->dvrr_classes[3][6] = int_stack + 3675;
 Libderiv->deriv_classes[4][6][0] = int_stack + 3955;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 4375;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 4435;
 Libderiv->deriv2_classes[2][5][143] = int_stack + 4525;
 Libderiv->deriv2_classes[2][6][143] = int_stack + 4651;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 4819;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 4919;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 5069;
 Libderiv->deriv2_classes[3][6][143] = int_stack + 5279;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 5559;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 5709;
 Libderiv->deriv2_classes[4][5][143] = int_stack + 5934;
 Libderiv->deriv2_classes[4][6][143] = int_stack + 6249;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 6669;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 6729;
 Libderiv->deriv2_classes[2][5][131] = int_stack + 6819;
 Libderiv->deriv2_classes[2][6][131] = int_stack + 6945;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 7113;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 7213;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 7363;
 Libderiv->deriv2_classes[3][6][131] = int_stack + 7573;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 7853;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 8003;
 Libderiv->deriv2_classes[4][5][131] = int_stack + 8228;
 Libderiv->deriv2_classes[4][6][131] = int_stack + 8543;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 8963;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 9023;
 Libderiv->deriv2_classes[2][5][130] = int_stack + 9113;
 Libderiv->deriv2_classes[2][6][130] = int_stack + 9239;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 9407;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 9507;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 9657;
 Libderiv->deriv2_classes[3][6][130] = int_stack + 9867;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 10147;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 10297;
 Libderiv->deriv2_classes[4][5][130] = int_stack + 10522;
 Libderiv->deriv2_classes[4][6][130] = int_stack + 10837;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 11257;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 11317;
 Libderiv->deriv2_classes[2][5][119] = int_stack + 11407;
 Libderiv->deriv2_classes[2][6][119] = int_stack + 11533;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 11701;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 11801;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 11951;
 Libderiv->deriv2_classes[3][6][119] = int_stack + 12161;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 12441;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 12591;
 Libderiv->deriv2_classes[4][5][119] = int_stack + 12816;
 Libderiv->deriv2_classes[4][6][119] = int_stack + 13131;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 13551;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 13611;
 Libderiv->deriv2_classes[2][5][118] = int_stack + 13701;
 Libderiv->deriv2_classes[2][6][118] = int_stack + 13827;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 13995;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 14095;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 14245;
 Libderiv->deriv2_classes[3][6][118] = int_stack + 14455;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 14735;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 14885;
 Libderiv->deriv2_classes[4][5][118] = int_stack + 15110;
 Libderiv->deriv2_classes[4][6][118] = int_stack + 15425;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 15845;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 15905;
 Libderiv->deriv2_classes[2][5][117] = int_stack + 15995;
 Libderiv->deriv2_classes[2][6][117] = int_stack + 16121;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 16289;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 16389;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 16539;
 Libderiv->deriv2_classes[3][6][117] = int_stack + 16749;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 17029;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 17179;
 Libderiv->deriv2_classes[4][5][117] = int_stack + 17404;
 Libderiv->deriv2_classes[4][6][117] = int_stack + 17719;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 18139;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 18199;
 Libderiv->deriv2_classes[2][5][107] = int_stack + 18289;
 Libderiv->deriv2_classes[2][6][107] = int_stack + 18415;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 18583;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 18683;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 18833;
 Libderiv->deriv2_classes[3][6][107] = int_stack + 19043;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 19323;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 19473;
 Libderiv->deriv2_classes[4][5][107] = int_stack + 19698;
 Libderiv->deriv2_classes[4][6][107] = int_stack + 20013;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 20433;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 20493;
 Libderiv->deriv2_classes[2][5][106] = int_stack + 20583;
 Libderiv->deriv2_classes[2][6][106] = int_stack + 20709;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 20877;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 20977;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 21127;
 Libderiv->deriv2_classes[3][6][106] = int_stack + 21337;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 21617;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 21767;
 Libderiv->deriv2_classes[4][5][106] = int_stack + 21992;
 Libderiv->deriv2_classes[4][6][106] = int_stack + 22307;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 22727;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 22787;
 Libderiv->deriv2_classes[2][5][105] = int_stack + 22877;
 Libderiv->deriv2_classes[2][6][105] = int_stack + 23003;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 23171;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 23271;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 23421;
 Libderiv->deriv2_classes[3][6][105] = int_stack + 23631;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 23911;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 24061;
 Libderiv->deriv2_classes[4][5][105] = int_stack + 24286;
 Libderiv->deriv2_classes[4][6][105] = int_stack + 24601;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 25021;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 25081;
 Libderiv->deriv2_classes[2][5][104] = int_stack + 25171;
 Libderiv->deriv2_classes[2][6][104] = int_stack + 25297;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 25465;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 25565;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 25715;
 Libderiv->deriv2_classes[3][6][104] = int_stack + 25925;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 26205;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 26355;
 Libderiv->deriv2_classes[4][5][104] = int_stack + 26580;
 Libderiv->deriv2_classes[4][6][104] = int_stack + 26895;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 27315;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 27375;
 Libderiv->deriv2_classes[2][5][95] = int_stack + 27465;
 Libderiv->deriv2_classes[2][6][95] = int_stack + 27591;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 27759;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 27859;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 28009;
 Libderiv->deriv2_classes[3][6][95] = int_stack + 28219;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 28499;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 28649;
 Libderiv->deriv2_classes[4][5][95] = int_stack + 28874;
 Libderiv->deriv2_classes[4][6][95] = int_stack + 29189;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 29609;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 29669;
 Libderiv->deriv2_classes[2][5][94] = int_stack + 29759;
 Libderiv->deriv2_classes[2][6][94] = int_stack + 29885;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 30053;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 30153;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 30303;
 Libderiv->deriv2_classes[3][6][94] = int_stack + 30513;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 30793;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 30943;
 Libderiv->deriv2_classes[4][5][94] = int_stack + 31168;
 Libderiv->deriv2_classes[4][6][94] = int_stack + 31483;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 31903;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 31963;
 Libderiv->deriv2_classes[2][5][93] = int_stack + 32053;
 Libderiv->deriv2_classes[2][6][93] = int_stack + 32179;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 32347;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 32447;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 32597;
 Libderiv->deriv2_classes[3][6][93] = int_stack + 32807;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 33087;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 33237;
 Libderiv->deriv2_classes[4][5][93] = int_stack + 33462;
 Libderiv->deriv2_classes[4][6][93] = int_stack + 33777;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 34197;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 34257;
 Libderiv->deriv2_classes[2][5][92] = int_stack + 34347;
 Libderiv->deriv2_classes[2][6][92] = int_stack + 34473;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 34641;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 34741;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 34891;
 Libderiv->deriv2_classes[3][6][92] = int_stack + 35101;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 35381;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 35531;
 Libderiv->deriv2_classes[4][5][92] = int_stack + 35756;
 Libderiv->deriv2_classes[4][6][92] = int_stack + 36071;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 36491;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 36551;
 Libderiv->deriv2_classes[2][5][91] = int_stack + 36641;
 Libderiv->deriv2_classes[2][6][91] = int_stack + 36767;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 36935;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 37035;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 37185;
 Libderiv->deriv2_classes[3][6][91] = int_stack + 37395;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 37675;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 37825;
 Libderiv->deriv2_classes[4][5][91] = int_stack + 38050;
 Libderiv->deriv2_classes[4][6][91] = int_stack + 38365;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 38785;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 38845;
 Libderiv->deriv2_classes[2][5][83] = int_stack + 38935;
 Libderiv->deriv2_classes[2][6][83] = int_stack + 39061;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 39229;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 39329;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 39479;
 Libderiv->deriv2_classes[3][6][83] = int_stack + 39689;
 Libderiv->deriv_classes[4][3][11] = int_stack + 39969;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 40119;
 Libderiv->deriv_classes[4][4][11] = int_stack + 40269;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 40494;
 Libderiv->deriv_classes[4][5][11] = int_stack + 40719;
 Libderiv->deriv2_classes[4][5][83] = int_stack + 41034;
 Libderiv->deriv2_classes[4][6][83] = int_stack + 41349;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 41769;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 41829;
 Libderiv->deriv2_classes[2][5][82] = int_stack + 41919;
 Libderiv->deriv2_classes[2][6][82] = int_stack + 42045;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 42213;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 42313;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 42463;
 Libderiv->deriv2_classes[3][6][82] = int_stack + 42673;
 Libderiv->deriv_classes[4][3][10] = int_stack + 42953;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 43103;
 Libderiv->deriv_classes[4][4][10] = int_stack + 43253;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 43478;
 Libderiv->deriv_classes[4][5][10] = int_stack + 43703;
 Libderiv->deriv2_classes[4][5][82] = int_stack + 44018;
 Libderiv->deriv2_classes[4][6][82] = int_stack + 44333;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 44753;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 44813;
 Libderiv->deriv2_classes[2][5][81] = int_stack + 44903;
 Libderiv->deriv2_classes[2][6][81] = int_stack + 45029;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 45197;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 45297;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 45447;
 Libderiv->deriv2_classes[3][6][81] = int_stack + 45657;
 Libderiv->deriv_classes[4][3][9] = int_stack + 45937;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 46087;
 Libderiv->deriv_classes[4][4][9] = int_stack + 46237;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 46462;
 Libderiv->deriv_classes[4][5][9] = int_stack + 46687;
 Libderiv->deriv2_classes[4][5][81] = int_stack + 47002;
 Libderiv->deriv2_classes[4][6][81] = int_stack + 47317;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 47737;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 47797;
 Libderiv->deriv2_classes[2][5][80] = int_stack + 47887;
 Libderiv->deriv2_classes[2][6][80] = int_stack + 48013;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 48181;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 48281;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 48431;
 Libderiv->deriv2_classes[3][6][80] = int_stack + 48641;
 Libderiv->deriv_classes[4][3][8] = int_stack + 48921;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 49071;
 Libderiv->deriv_classes[4][4][8] = int_stack + 49221;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 49446;
 Libderiv->deriv_classes[4][5][8] = int_stack + 49671;
 Libderiv->deriv2_classes[4][5][80] = int_stack + 49986;
 Libderiv->deriv2_classes[4][6][80] = int_stack + 50301;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 50721;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 50781;
 Libderiv->deriv2_classes[2][5][79] = int_stack + 50871;
 Libderiv->deriv2_classes[2][6][79] = int_stack + 50997;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 51165;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 51265;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 51415;
 Libderiv->deriv2_classes[3][6][79] = int_stack + 51625;
 Libderiv->deriv_classes[4][3][7] = int_stack + 51905;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 52055;
 Libderiv->deriv_classes[4][4][7] = int_stack + 52205;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 52430;
 Libderiv->deriv_classes[4][5][7] = int_stack + 52655;
 Libderiv->deriv2_classes[4][5][79] = int_stack + 52970;
 Libderiv->deriv2_classes[4][6][79] = int_stack + 53285;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 53705;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 53765;
 Libderiv->deriv2_classes[2][5][78] = int_stack + 53855;
 Libderiv->deriv2_classes[2][6][78] = int_stack + 53981;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 54149;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 54249;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 54399;
 Libderiv->deriv2_classes[3][6][78] = int_stack + 54609;
 Libderiv->dvrr_classes[4][3] = int_stack + 54889;
 Libderiv->deriv_classes[4][3][6] = int_stack + 55039;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 55189;
 Libderiv->dvrr_classes[4][4] = int_stack + 55339;
 Libderiv->deriv_classes[4][4][6] = int_stack + 55564;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 55789;
 Libderiv->deriv_classes[4][5][6] = int_stack + 56014;
 Libderiv->deriv2_classes[4][5][78] = int_stack + 56329;
 Libderiv->deriv2_classes[4][6][78] = int_stack + 56644;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 57064;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 57124;
 Libderiv->deriv2_classes[2][5][35] = int_stack + 57214;
 Libderiv->deriv2_classes[2][6][35] = int_stack + 57340;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 57508;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 57608;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 57758;
 Libderiv->deriv2_classes[3][6][35] = int_stack + 57968;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 58248;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 58398;
 Libderiv->deriv2_classes[4][5][35] = int_stack + 58623;
 Libderiv->deriv2_classes[4][6][35] = int_stack + 58938;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 59358;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 59418;
 Libderiv->deriv2_classes[2][5][34] = int_stack + 59508;
 Libderiv->deriv2_classes[2][6][34] = int_stack + 59634;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 59802;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 59902;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 60052;
 Libderiv->deriv2_classes[3][6][34] = int_stack + 60262;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 60542;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 60692;
 Libderiv->deriv2_classes[4][5][34] = int_stack + 60917;
 Libderiv->deriv2_classes[4][6][34] = int_stack + 61232;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 61652;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 61712;
 Libderiv->deriv2_classes[2][5][33] = int_stack + 61802;
 Libderiv->deriv2_classes[2][6][33] = int_stack + 61928;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 62096;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 62196;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 62346;
 Libderiv->deriv2_classes[3][6][33] = int_stack + 62556;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 62836;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 62986;
 Libderiv->deriv2_classes[4][5][33] = int_stack + 63211;
 Libderiv->deriv2_classes[4][6][33] = int_stack + 63526;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 63946;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 64006;
 Libderiv->deriv2_classes[2][5][32] = int_stack + 64096;
 Libderiv->deriv2_classes[2][6][32] = int_stack + 64222;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 64390;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 64490;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 64640;
 Libderiv->deriv2_classes[3][6][32] = int_stack + 64850;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 65130;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 65280;
 Libderiv->deriv2_classes[4][5][32] = int_stack + 65505;
 Libderiv->deriv2_classes[4][6][32] = int_stack + 65820;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 66240;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 66300;
 Libderiv->deriv2_classes[2][5][31] = int_stack + 66390;
 Libderiv->deriv2_classes[2][6][31] = int_stack + 66516;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 66684;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 66784;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 66934;
 Libderiv->deriv2_classes[3][6][31] = int_stack + 67144;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 67424;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 67574;
 Libderiv->deriv2_classes[4][5][31] = int_stack + 67799;
 Libderiv->deriv2_classes[4][6][31] = int_stack + 68114;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 68534;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 68594;
 Libderiv->deriv2_classes[2][5][30] = int_stack + 68684;
 Libderiv->deriv2_classes[2][6][30] = int_stack + 68810;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 68978;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 69078;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 69228;
 Libderiv->deriv2_classes[3][6][30] = int_stack + 69438;
 Libderiv->deriv_classes[4][3][2] = int_stack + 69718;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 69868;
 Libderiv->deriv_classes[4][4][2] = int_stack + 70018;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 70243;
 Libderiv->deriv_classes[4][5][2] = int_stack + 70468;
 Libderiv->deriv2_classes[4][5][30] = int_stack + 70783;
 Libderiv->deriv2_classes[4][6][30] = int_stack + 71098;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 71518;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 71578;
 Libderiv->deriv2_classes[2][5][26] = int_stack + 71668;
 Libderiv->deriv2_classes[2][6][26] = int_stack + 71794;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 71962;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 72062;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 72212;
 Libderiv->deriv2_classes[3][6][26] = int_stack + 72422;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 72702;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 72852;
 Libderiv->deriv2_classes[4][5][26] = int_stack + 73077;
 Libderiv->deriv2_classes[4][6][26] = int_stack + 73392;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 73812;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 73872;
 Libderiv->deriv2_classes[2][5][23] = int_stack + 73962;
 Libderiv->deriv2_classes[2][6][23] = int_stack + 74088;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 74256;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 74356;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 74506;
 Libderiv->deriv2_classes[3][6][23] = int_stack + 74716;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 74996;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 75146;
 Libderiv->deriv2_classes[4][5][23] = int_stack + 75371;
 Libderiv->deriv2_classes[4][6][23] = int_stack + 75686;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 76106;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 76166;
 Libderiv->deriv2_classes[2][5][22] = int_stack + 76256;
 Libderiv->deriv2_classes[2][6][22] = int_stack + 76382;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 76550;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 76650;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 76800;
 Libderiv->deriv2_classes[3][6][22] = int_stack + 77010;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 77290;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 77440;
 Libderiv->deriv2_classes[4][5][22] = int_stack + 77665;
 Libderiv->deriv2_classes[4][6][22] = int_stack + 77980;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 78400;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 78460;
 Libderiv->deriv2_classes[2][5][21] = int_stack + 78550;
 Libderiv->deriv2_classes[2][6][21] = int_stack + 78676;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 78844;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 78944;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 79094;
 Libderiv->deriv2_classes[3][6][21] = int_stack + 79304;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 79584;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 79734;
 Libderiv->deriv2_classes[4][5][21] = int_stack + 79959;
 Libderiv->deriv2_classes[4][6][21] = int_stack + 80274;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 80694;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 80754;
 Libderiv->deriv2_classes[2][5][20] = int_stack + 80844;
 Libderiv->deriv2_classes[2][6][20] = int_stack + 80970;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 81138;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 81238;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 81388;
 Libderiv->deriv2_classes[3][6][20] = int_stack + 81598;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 81878;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 82028;
 Libderiv->deriv2_classes[4][5][20] = int_stack + 82253;
 Libderiv->deriv2_classes[4][6][20] = int_stack + 82568;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 82988;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 83048;
 Libderiv->deriv2_classes[2][5][19] = int_stack + 83138;
 Libderiv->deriv2_classes[2][6][19] = int_stack + 83264;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 83432;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 83532;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 83682;
 Libderiv->deriv2_classes[3][6][19] = int_stack + 83892;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 84172;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 84322;
 Libderiv->deriv2_classes[4][5][19] = int_stack + 84547;
 Libderiv->deriv2_classes[4][6][19] = int_stack + 84862;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 85282;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 85342;
 Libderiv->deriv2_classes[2][5][18] = int_stack + 85432;
 Libderiv->deriv2_classes[2][6][18] = int_stack + 85558;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 85726;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 85826;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 85976;
 Libderiv->deriv2_classes[3][6][18] = int_stack + 86186;
 Libderiv->deriv_classes[4][3][1] = int_stack + 86466;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 86616;
 Libderiv->deriv_classes[4][4][1] = int_stack + 86766;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 86991;
 Libderiv->deriv_classes[4][5][1] = int_stack + 87216;
 Libderiv->deriv2_classes[4][5][18] = int_stack + 87531;
 Libderiv->deriv2_classes[4][6][18] = int_stack + 87846;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 88266;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 88326;
 Libderiv->deriv2_classes[2][5][14] = int_stack + 88416;
 Libderiv->deriv2_classes[2][6][14] = int_stack + 88542;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 88710;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 88810;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 88960;
 Libderiv->deriv2_classes[3][6][14] = int_stack + 89170;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 89450;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 89600;
 Libderiv->deriv2_classes[4][5][14] = int_stack + 89825;
 Libderiv->deriv2_classes[4][6][14] = int_stack + 90140;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 90560;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 90620;
 Libderiv->deriv2_classes[2][5][13] = int_stack + 90710;
 Libderiv->deriv2_classes[2][6][13] = int_stack + 90836;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 91004;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 91104;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 91254;
 Libderiv->deriv2_classes[3][6][13] = int_stack + 91464;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 91744;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 91894;
 Libderiv->deriv2_classes[4][5][13] = int_stack + 92119;
 Libderiv->deriv2_classes[4][6][13] = int_stack + 92434;
 Libderiv->deriv_classes[2][3][11] = int_stack + 92854;
 Libderiv->deriv_classes[2][4][11] = int_stack + 92914;
 Libderiv->deriv_classes[2][5][11] = int_stack + 93004;
 Libderiv->deriv_classes[2][6][11] = int_stack + 93130;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 93298;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 93358;
 Libderiv->deriv2_classes[2][5][11] = int_stack + 93448;
 Libderiv->deriv2_classes[2][6][11] = int_stack + 93574;
 Libderiv->deriv_classes[3][3][11] = int_stack + 93742;
 Libderiv->deriv_classes[3][4][11] = int_stack + 93842;
 Libderiv->deriv_classes[3][5][11] = int_stack + 93992;
 Libderiv->deriv_classes[3][6][11] = int_stack + 94202;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 94482;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 94582;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 94732;
 Libderiv->deriv2_classes[3][6][11] = int_stack + 94942;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 95222;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 95372;
 Libderiv->deriv2_classes[4][5][11] = int_stack + 95597;
 Libderiv->deriv2_classes[4][6][11] = int_stack + 95912;
 Libderiv->deriv_classes[2][3][10] = int_stack + 96332;
 Libderiv->deriv_classes[2][4][10] = int_stack + 96392;
 Libderiv->deriv_classes[2][5][10] = int_stack + 96482;
 Libderiv->deriv_classes[2][6][10] = int_stack + 96608;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 96776;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 96836;
 Libderiv->deriv2_classes[2][5][10] = int_stack + 96926;
 Libderiv->deriv2_classes[2][6][10] = int_stack + 97052;
 Libderiv->deriv_classes[3][3][10] = int_stack + 97220;
 Libderiv->deriv_classes[3][4][10] = int_stack + 97320;
 Libderiv->deriv_classes[3][5][10] = int_stack + 97470;
 Libderiv->deriv_classes[3][6][10] = int_stack + 97680;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 97960;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 98060;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 98210;
 Libderiv->deriv2_classes[3][6][10] = int_stack + 98420;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 98700;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 98850;
 Libderiv->deriv2_classes[4][5][10] = int_stack + 99075;
 Libderiv->deriv2_classes[4][6][10] = int_stack + 99390;
 Libderiv->deriv_classes[2][3][9] = int_stack + 99810;
 Libderiv->deriv_classes[2][4][9] = int_stack + 99870;
 Libderiv->deriv_classes[2][5][9] = int_stack + 99960;
 Libderiv->deriv_classes[2][6][9] = int_stack + 100086;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 100254;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 100314;
 Libderiv->deriv2_classes[2][5][9] = int_stack + 100404;
 Libderiv->deriv2_classes[2][6][9] = int_stack + 100530;
 Libderiv->deriv_classes[3][3][9] = int_stack + 100698;
 Libderiv->deriv_classes[3][4][9] = int_stack + 100798;
 Libderiv->deriv_classes[3][5][9] = int_stack + 100948;
 Libderiv->deriv_classes[3][6][9] = int_stack + 101158;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 101438;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 101538;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 101688;
 Libderiv->deriv2_classes[3][6][9] = int_stack + 101898;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 102178;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 102328;
 Libderiv->deriv2_classes[4][5][9] = int_stack + 102553;
 Libderiv->deriv2_classes[4][6][9] = int_stack + 102868;
 Libderiv->deriv_classes[2][3][8] = int_stack + 103288;
 Libderiv->deriv_classes[2][4][8] = int_stack + 103348;
 Libderiv->deriv_classes[2][5][8] = int_stack + 103438;
 Libderiv->deriv_classes[2][6][8] = int_stack + 103564;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 103732;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 103792;
 Libderiv->deriv2_classes[2][5][8] = int_stack + 103882;
 Libderiv->deriv2_classes[2][6][8] = int_stack + 104008;
 Libderiv->deriv_classes[3][3][8] = int_stack + 104176;
 Libderiv->deriv_classes[3][4][8] = int_stack + 104276;
 Libderiv->deriv_classes[3][5][8] = int_stack + 104426;
 Libderiv->deriv_classes[3][6][8] = int_stack + 104636;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 104916;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 105016;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 105166;
 Libderiv->deriv2_classes[3][6][8] = int_stack + 105376;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 105656;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 105806;
 Libderiv->deriv2_classes[4][5][8] = int_stack + 106031;
 Libderiv->deriv2_classes[4][6][8] = int_stack + 106346;
 Libderiv->deriv_classes[2][3][7] = int_stack + 106766;
 Libderiv->deriv_classes[2][4][7] = int_stack + 106826;
 Libderiv->deriv_classes[2][5][7] = int_stack + 106916;
 Libderiv->deriv_classes[2][6][7] = int_stack + 107042;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 107210;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 107270;
 Libderiv->deriv2_classes[2][5][7] = int_stack + 107360;
 Libderiv->deriv2_classes[2][6][7] = int_stack + 107486;
 Libderiv->deriv_classes[3][3][7] = int_stack + 107654;
 Libderiv->deriv_classes[3][4][7] = int_stack + 107754;
 Libderiv->deriv_classes[3][5][7] = int_stack + 107904;
 Libderiv->deriv_classes[3][6][7] = int_stack + 108114;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 108394;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 108494;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 108644;
 Libderiv->deriv2_classes[3][6][7] = int_stack + 108854;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 109134;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 109284;
 Libderiv->deriv2_classes[4][5][7] = int_stack + 109509;
 Libderiv->deriv2_classes[4][6][7] = int_stack + 109824;
 Libderiv->deriv_classes[2][3][6] = int_stack + 110244;
 Libderiv->deriv_classes[2][4][6] = int_stack + 110304;
 Libderiv->deriv_classes[2][5][6] = int_stack + 110394;
 Libderiv->deriv_classes[2][6][6] = int_stack + 110520;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 110688;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 110748;
 Libderiv->deriv2_classes[2][5][6] = int_stack + 110838;
 Libderiv->deriv2_classes[2][6][6] = int_stack + 110964;
 Libderiv->dvrr_classes[3][3] = int_stack + 111132;
 Libderiv->deriv_classes[3][3][6] = int_stack + 111232;
 Libderiv->dvrr_classes[3][4] = int_stack + 111332;
 Libderiv->deriv_classes[3][4][6] = int_stack + 111482;
 Libderiv->dvrr_classes[3][5] = int_stack + 111632;
 Libderiv->deriv_classes[3][5][6] = int_stack + 111842;
 Libderiv->deriv_classes[3][6][6] = int_stack + 112052;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 112332;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 112432;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 112582;
 Libderiv->deriv2_classes[3][6][6] = int_stack + 112792;
 Libderiv->deriv_classes[4][3][0] = int_stack + 113072;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 113222;
 Libderiv->deriv_classes[4][4][0] = int_stack + 113372;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 113597;
 Libderiv->deriv_classes[4][5][0] = int_stack + 113822;
 Libderiv->deriv2_classes[4][5][6] = int_stack + 114137;
 Libderiv->deriv2_classes[4][6][6] = int_stack + 114452;
 Libderiv->deriv_classes[2][3][2] = int_stack + 114872;
 Libderiv->deriv_classes[2][4][2] = int_stack + 114932;
 Libderiv->deriv_classes[2][5][2] = int_stack + 115022;
 Libderiv->deriv_classes[2][6][2] = int_stack + 115148;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 115316;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 115376;
 Libderiv->deriv2_classes[2][5][2] = int_stack + 115466;
 Libderiv->deriv2_classes[2][6][2] = int_stack + 115592;
 Libderiv->deriv_classes[3][3][2] = int_stack + 115760;
 Libderiv->deriv_classes[3][4][2] = int_stack + 115860;
 Libderiv->deriv_classes[3][5][2] = int_stack + 116010;
 Libderiv->deriv_classes[3][6][2] = int_stack + 116220;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 116500;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 116600;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 116750;
 Libderiv->deriv2_classes[3][6][2] = int_stack + 116960;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 117240;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 117390;
 Libderiv->deriv2_classes[4][5][2] = int_stack + 117615;
 Libderiv->deriv2_classes[4][6][2] = int_stack + 117930;
 Libderiv->deriv_classes[2][3][1] = int_stack + 118350;
 Libderiv->deriv_classes[2][4][1] = int_stack + 118410;
 Libderiv->deriv_classes[2][5][1] = int_stack + 118500;
 Libderiv->deriv_classes[2][6][1] = int_stack + 118626;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 118794;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 118854;
 Libderiv->deriv2_classes[2][5][1] = int_stack + 118944;
 Libderiv->deriv2_classes[2][6][1] = int_stack + 119070;
 Libderiv->deriv_classes[3][3][1] = int_stack + 119238;
 Libderiv->deriv_classes[3][4][1] = int_stack + 119338;
 Libderiv->deriv_classes[3][5][1] = int_stack + 119488;
 Libderiv->deriv_classes[3][6][1] = int_stack + 119698;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 119978;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 120078;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 120228;
 Libderiv->deriv2_classes[3][6][1] = int_stack + 120438;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 120718;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 120868;
 Libderiv->deriv2_classes[4][5][1] = int_stack + 121093;
 Libderiv->deriv2_classes[4][6][1] = int_stack + 121408;
 Libderiv->dvrr_classes[2][3] = int_stack + 121828;
 Libderiv->dvrr_classes[2][4] = int_stack + 121888;
 Libderiv->dvrr_classes[2][5] = int_stack + 121978;
 Libderiv->dvrr_classes[2][6] = int_stack + 122104;
 Libderiv->deriv_classes[2][3][0] = int_stack + 122272;
 Libderiv->deriv_classes[2][4][0] = int_stack + 122332;
 Libderiv->deriv_classes[2][5][0] = int_stack + 122422;
 Libderiv->deriv_classes[2][6][0] = int_stack + 122548;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 122716;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 122776;
 Libderiv->deriv2_classes[2][5][0] = int_stack + 122866;
 Libderiv->deriv2_classes[2][6][0] = int_stack + 122992;
 Libderiv->deriv_classes[3][3][0] = int_stack + 123160;
 Libderiv->deriv_classes[3][4][0] = int_stack + 123260;
 Libderiv->deriv_classes[3][5][0] = int_stack + 123410;
 Libderiv->deriv_classes[3][6][0] = int_stack + 123620;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 123900;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 124000;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 124150;
 Libderiv->deriv2_classes[3][6][0] = int_stack + 124360;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 124640;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 124790;
 Libderiv->deriv2_classes[4][5][0] = int_stack + 125015;
 Libderiv->deriv2_classes[4][6][0] = int_stack + 125330;
 memset(int_stack,0,1006000);

 Libderiv->dvrr_stack = int_stack + 301578;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ddff(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+125750,int_stack+121888,int_stack+121828,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+125930,int_stack+121978,int_stack+121888,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+126200,int_stack+125930,int_stack+125750,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+126560,int_stack+92914,int_stack+92854, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121828,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+126740,int_stack+93004,int_stack+92914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121888,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127010,int_stack+126740,int_stack+126560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125750,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+127370,int_stack+93130,int_stack+93004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121978,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+127748,int_stack+127370,int_stack+126740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125930,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+128288,int_stack+127748,int_stack+127010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126200,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+127370,int_stack+111332,int_stack+111132,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+127670,int_stack+111632,int_stack+111332,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+128888,int_stack+127670,int_stack+127370,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+129488,int_stack+93842,int_stack+93742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111132,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+129788,int_stack+93992,int_stack+93842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111332,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130238,int_stack+129788,int_stack+129488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127370,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+130838,int_stack+94202,int_stack+93992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111632,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+131468,int_stack+130838,int_stack+129788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127670,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+132368,int_stack+131468,int_stack+130238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128888,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+133368,int_stack+132368,int_stack+128288,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+130838,int_stack+55339,int_stack+54889,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+131288,int_stack+2100,int_stack+55339,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+131288,int_stack+130838,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+136068,int_stack+40269,int_stack+39969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54889,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+136518,int_stack+40719,int_stack+40269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55339,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+137193,int_stack+136518,int_stack+136068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130838,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+138093,int_stack+0,int_stack+40719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+139038,int_stack+138093,int_stack+136518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131288,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+140388,int_stack+139038,int_stack+137193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135168,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+141888,int_stack+140388,int_stack+132368,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+138093,int_stack+96392,int_stack+96332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121828, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+138273,int_stack+96482,int_stack+96392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121888, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+138543,int_stack+138273,int_stack+138093, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125750, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+138903,int_stack+96608,int_stack+96482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121978, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+139281,int_stack+138903,int_stack+138273, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125930, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+139821,int_stack+139281,int_stack+138543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126200, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+138903,int_stack+97320,int_stack+97220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111132, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+139203,int_stack+97470,int_stack+97320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111332, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+140421,int_stack+139203,int_stack+138903, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127370, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+141021,int_stack+97680,int_stack+97470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111632, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+144888,int_stack+141021,int_stack+139203, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127670, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+145788,int_stack+144888,int_stack+140421, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128888, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+146788,int_stack+145788,int_stack+139821,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+144888,int_stack+43253,int_stack+42953, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54889, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+141021,int_stack+43703,int_stack+43253, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55339, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+148588,int_stack+141021,int_stack+144888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130838, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+149488,int_stack+420,int_stack+43703, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+150433,int_stack+149488,int_stack+141021, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131288, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+151783,int_stack+150433,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135168, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+153283,int_stack+151783,int_stack+145788,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+149488,int_stack+99870,int_stack+99810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121828, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+149668,int_stack+99960,int_stack+99870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121888, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+149938,int_stack+149668,int_stack+149488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125750, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+150298,int_stack+100086,int_stack+99960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121978, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+150676,int_stack+150298,int_stack+149668, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125930, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+151216,int_stack+150676,int_stack+149938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126200, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+150298,int_stack+100798,int_stack+100698, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111132, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+145338,int_stack+100948,int_stack+100798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111332, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+150598,int_stack+145338,int_stack+150298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127370, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+151816,int_stack+101158,int_stack+100948, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111632, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+156283,int_stack+151816,int_stack+145338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127670, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+151816,int_stack+156283,int_stack+150598, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128888, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+156283,int_stack+151816,int_stack+151216,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+158083,int_stack+46237,int_stack+45937, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54889, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+46687,int_stack+46237, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55339, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+158533,int_stack+0,int_stack+158083, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130838, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+159433,int_stack+840,int_stack+46687, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+160378,int_stack+159433,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131288, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+161728,int_stack+160378,int_stack+158533, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135168, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+163228,int_stack+161728,int_stack+151816,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+159433,int_stack+103348,int_stack+103288, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+159613,int_stack+103438,int_stack+103348, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+159883,int_stack+159613,int_stack+159433, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+160243,int_stack+103564,int_stack+103438, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+160621,int_stack+160243,int_stack+159613, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+161161,int_stack+160621,int_stack+159883, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+160243,int_stack+104276,int_stack+104176, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+160543,int_stack+104426,int_stack+104276, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+161761,int_stack+160543,int_stack+160243, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+162361,int_stack+104636,int_stack+104426, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+166228,int_stack+162361,int_stack+160543, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+167128,int_stack+166228,int_stack+161761, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+168128,int_stack+167128,int_stack+161161,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+166228,int_stack+49221,int_stack+48921, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54889, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+162361,int_stack+49671,int_stack+49221, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55339, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169928,int_stack+162361,int_stack+166228, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+170828,int_stack+1260,int_stack+49671, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+171773,int_stack+170828,int_stack+162361, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+173123,int_stack+171773,int_stack+169928, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+174623,int_stack+173123,int_stack+167128,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+170828,int_stack+106826,int_stack+106766, 0.0, zero_stack, 1.0, int_stack+121828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+171008,int_stack+106916,int_stack+106826, 0.0, zero_stack, 1.0, int_stack+121888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+171278,int_stack+171008,int_stack+170828, 0.0, zero_stack, 1.0, int_stack+125750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+177623,int_stack+107042,int_stack+106916, 0.0, zero_stack, 1.0, int_stack+121978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+171638,int_stack+177623,int_stack+171008, 0.0, zero_stack, 1.0, int_stack+125930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+177623,int_stack+171638,int_stack+171278, 0.0, zero_stack, 1.0, int_stack+126200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+171638,int_stack+107754,int_stack+107654, 0.0, zero_stack, 1.0, int_stack+111132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+166678,int_stack+107904,int_stack+107754, 0.0, zero_stack, 1.0, int_stack+111332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+171938,int_stack+166678,int_stack+171638, 0.0, zero_stack, 1.0, int_stack+127370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+172538,int_stack+108114,int_stack+107904, 0.0, zero_stack, 1.0, int_stack+111632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+173168,int_stack+172538,int_stack+166678, 0.0, zero_stack, 1.0, int_stack+127670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+675,int_stack+173168,int_stack+171938, 0.0, zero_stack, 1.0, int_stack+128888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+172538,int_stack+675,int_stack+177623,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+152816,int_stack+52205,int_stack+51905, 0.0, zero_stack, 1.0, int_stack+54889, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+178223,int_stack+52655,int_stack+52205, 0.0, zero_stack, 1.0, int_stack+55339, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+178898,int_stack+178223,int_stack+152816, 0.0, zero_stack, 1.0, int_stack+130838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+179798,int_stack+1680,int_stack+52655, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+180743,int_stack+179798,int_stack+178223, 0.0, zero_stack, 1.0, int_stack+131288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+182093,int_stack+180743,int_stack+178898, 0.0, zero_stack, 1.0, int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+183593,int_stack+182093,int_stack+675,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+179798,int_stack+110304,int_stack+110244, 1.0, int_stack+121828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+179978,int_stack+110394,int_stack+110304, 1.0, int_stack+121888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+180248,int_stack+179978,int_stack+179798, 1.0, int_stack+125750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+180608,int_stack+110520,int_stack+110394, 1.0, int_stack+121978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+180986,int_stack+180608,int_stack+179978, 1.0, int_stack+125930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+181526,int_stack+180986,int_stack+180248, 1.0, int_stack+126200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180608,int_stack+111482,int_stack+111232, 1.0, int_stack+111132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+180908,int_stack+111842,int_stack+111482, 1.0, int_stack+111332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+182126,int_stack+180908,int_stack+180608, 1.0, int_stack+127370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+182726,int_stack+112052,int_stack+111842, 1.0, int_stack+111632, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+186593,int_stack+182726,int_stack+180908, 1.0, int_stack+127670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+187493,int_stack+186593,int_stack+182126, 1.0, int_stack+128888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+188493,int_stack+187493,int_stack+181526,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+186593,int_stack+55564,int_stack+55039, 1.0, int_stack+54889, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+182726,int_stack+56014,int_stack+55564, 1.0, int_stack+55339, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+190293,int_stack+182726,int_stack+186593, 1.0, int_stack+130838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+191193,int_stack+2415,int_stack+56014, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+192138,int_stack+191193,int_stack+182726, 1.0, int_stack+131288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+130838,int_stack+192138,int_stack+190293, 1.0, int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+191193,int_stack+130838,int_stack+187493,100);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130838,int_stack+122104,int_stack+121978,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+131216,int_stack+130838,int_stack+125930,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+131756,int_stack+131216,int_stack+126200,6);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130838,int_stack+3675,int_stack+111632,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+135168,int_stack+130838,int_stack+127670,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+135168,int_stack+128888,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+194193,int_stack+1675,int_stack+131756,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+128888,int_stack+114932,int_stack+114872,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+129068,int_stack+115022,int_stack+114932,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+129068,int_stack+128888,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+135528,int_stack+115148,int_stack+115022,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+130838,int_stack+135528,int_stack+129068,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+127370,int_stack+130838,int_stack+135168,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+130838,int_stack+115860,int_stack+115760,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+187043,int_stack+116010,int_stack+115860,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+131138,int_stack+187043,int_stack+130838,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+125750,int_stack+116220,int_stack+116010,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+195993,int_stack+125750,int_stack+187043,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+196893,int_stack+195993,int_stack+131138,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+197893,int_stack+196893,int_stack+127370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+195993,int_stack+70018,int_stack+69718,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+125750,int_stack+70468,int_stack+70018,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+199693,int_stack+125750,int_stack+195993,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+200593,int_stack+2835,int_stack+70468,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+201538,int_stack+200593,int_stack+125750,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+202888,int_stack+201538,int_stack+199693,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+204388,int_stack+202888,int_stack+196893, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+200593,int_stack+118410,int_stack+118350,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+200773,int_stack+118500,int_stack+118410,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+201043,int_stack+200773,int_stack+200593,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+201403,int_stack+118626,int_stack+118500,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+135528,int_stack+201403,int_stack+200773,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+201403,int_stack+135528,int_stack+201043,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+135528,int_stack+119338,int_stack+119238,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+196443,int_stack+119488,int_stack+119338,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+202003,int_stack+196443,int_stack+135528,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+202603,int_stack+119698,int_stack+119488,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+203233,int_stack+202603,int_stack+196443,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+207388,int_stack+203233,int_stack+202003,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+208388,int_stack+207388,int_stack+201403, 0.0, zero_stack, 1.0, int_stack+131756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+202603,int_stack+86766,int_stack+86466,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+203053,int_stack+87216,int_stack+86766,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+210188,int_stack+203053,int_stack+202603,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+211088,int_stack+3255,int_stack+87216,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+212033,int_stack+211088,int_stack+203053,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+213383,int_stack+212033,int_stack+210188,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+214883,int_stack+213383,int_stack+207388, 0.0, zero_stack, 1.0, int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+211088,int_stack+122332,int_stack+122272,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+211268,int_stack+122422,int_stack+122332,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+211538,int_stack+211268,int_stack+211088,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+211898,int_stack+122548,int_stack+122422,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+212276,int_stack+211898,int_stack+211268,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+212816,int_stack+212276,int_stack+211538,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+211898,int_stack+123260,int_stack+123160,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+212198,int_stack+123410,int_stack+123260,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+213416,int_stack+212198,int_stack+211898,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+214016,int_stack+123620,int_stack+123410,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2675,int_stack+214016,int_stack+212198,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+217883,int_stack+2675,int_stack+213416,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+218883,int_stack+217883,int_stack+212816, 1.0, int_stack+131756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2675,int_stack+113372,int_stack+113072,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3125,int_stack+113822,int_stack+113372,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+220683,int_stack+3125,int_stack+2675,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+221583,int_stack+3955,int_stack+113822,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+222528,int_stack+221583,int_stack+3125,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+223878,int_stack+222528,int_stack+220683,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+225378,int_stack+223878,int_stack+217883, 1.0, int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+4435,int_stack+4375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+92854,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+4525,int_stack+4435, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+92914,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+126560,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+221583,int_stack+4651,int_stack+4525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+93004,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+221961,int_stack+221583,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+126740,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+222501,int_stack+221961,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+127010,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+221583,int_stack+4919,int_stack+4819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+93742,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+221883,int_stack+5069,int_stack+4919, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+93842,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1675,int_stack+221883,int_stack+221583, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+129488,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+5279,int_stack+5069, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+93992,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+223101,int_stack+131738,int_stack+221883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+129788,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+224001,int_stack+223101,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+130238,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+228378,int_stack+224001,int_stack+222501,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+5709,int_stack+5559, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+39969,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+221583,int_stack+5934,int_stack+5709, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40269,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+222258,int_stack+221583,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+136068,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+6249,int_stack+5934, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40719,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3800,int_stack+1675,int_stack+221583, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+136518,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5150,int_stack+3800,int_stack+222258, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+137193,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+230178,int_stack+5150,int_stack+224001,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3800,int_stack+6729,int_stack+6669, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92854, 1.0, int_stack+96332,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3980,int_stack+6819,int_stack+6729, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92914, 1.0, int_stack+96392,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4250,int_stack+3980,int_stack+3800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126560, 1.0, int_stack+138093,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4610,int_stack+6945,int_stack+6819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93004, 1.0, int_stack+96482,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4988,int_stack+4610,int_stack+3980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126740, 1.0, int_stack+138273,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5528,int_stack+4988,int_stack+4250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127010, 1.0, int_stack+138543,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3800,int_stack+7213,int_stack+7113, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93742, 1.0, int_stack+97220,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4100,int_stack+7363,int_stack+7213, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93842, 1.0, int_stack+97320,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4550,int_stack+4100,int_stack+3800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129488, 1.0, int_stack+138903,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+7573,int_stack+7363, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93992, 1.0, int_stack+97470,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6128,int_stack+131738,int_stack+4100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129788, 1.0, int_stack+139203,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+6128,int_stack+4550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130238, 1.0, int_stack+140421,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+221583,int_stack+1675,int_stack+5528,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+131738,int_stack+8003,int_stack+7853, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39969, 1.0, int_stack+42953,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3800,int_stack+8228,int_stack+8003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40269, 1.0, int_stack+43253,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4475,int_stack+3800,int_stack+131738, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136068, 1.0, int_stack+144888,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5375,int_stack+8543,int_stack+8228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40719, 1.0, int_stack+43703,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6320,int_stack+5375,int_stack+3800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136518, 1.0, int_stack+141021,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+223383,int_stack+6320,int_stack+4475, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137193, 1.0, int_stack+148588,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3800,int_stack+223383,int_stack+1675,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+9023,int_stack+8963, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+96332, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+9113,int_stack+9023, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+96392, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+138093, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+223383,int_stack+9239,int_stack+9113, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+96482, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+223761,int_stack+223383,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+138273, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+224301,int_stack+223761,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+138543, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+223383,int_stack+9507,int_stack+9407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+97220, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+223683,int_stack+9657,int_stack+9507, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+97320, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1675,int_stack+223683,int_stack+223383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+138903, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+9867,int_stack+9657, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+97470, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6800,int_stack+131738,int_stack+223683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+139203, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7700,int_stack+6800,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+140421, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+233178,int_stack+7700,int_stack+224301,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+10297,int_stack+10147, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+42953, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6800,int_stack+10522,int_stack+10297, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43253, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+223383,int_stack+6800,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+144888, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+10837,int_stack+10522, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43703, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8700,int_stack+1675,int_stack+6800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+141021, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+8700,int_stack+223383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+148588, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+236478,int_stack+234978,int_stack+7700,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+11317,int_stack+11257, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92854, 0.0, zero_stack, 1.0, int_stack+99810,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235158,int_stack+11407,int_stack+11317, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92914, 0.0, zero_stack, 1.0, int_stack+99870,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235428,int_stack+235158,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126560, 0.0, zero_stack, 1.0, int_stack+149488,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+235788,int_stack+11533,int_stack+11407, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93004, 0.0, zero_stack, 1.0, int_stack+99960,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+223383,int_stack+235788,int_stack+235158, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126740, 0.0, zero_stack, 1.0, int_stack+149668,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+235788,int_stack+223383,int_stack+235428, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127010, 0.0, zero_stack, 1.0, int_stack+149938,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+223383,int_stack+11801,int_stack+11701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93742, 0.0, zero_stack, 1.0, int_stack+100698,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+223683,int_stack+11951,int_stack+11801, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93842, 0.0, zero_stack, 1.0, int_stack+100798,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+224133,int_stack+223683,int_stack+223383, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129488, 0.0, zero_stack, 1.0, int_stack+150298,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+12161,int_stack+11951, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93992, 0.0, zero_stack, 1.0, int_stack+100948,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6800,int_stack+131738,int_stack+223683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129788, 0.0, zero_stack, 1.0, int_stack+145338,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+6800,int_stack+224133, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130238, 0.0, zero_stack, 1.0, int_stack+150598,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6800,int_stack+1675,int_stack+235788,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8600,int_stack+12591,int_stack+12441, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39969, 0.0, zero_stack, 1.0, int_stack+45937,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9050,int_stack+12816,int_stack+12591, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40269, 0.0, zero_stack, 1.0, int_stack+46237,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9725,int_stack+9050,int_stack+8600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136068, 0.0, zero_stack, 1.0, int_stack+158083,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10625,int_stack+13131,int_stack+12816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40719, 0.0, zero_stack, 1.0, int_stack+46687,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11570,int_stack+10625,int_stack+9050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136518, 0.0, zero_stack, 1.0, int_stack+0,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+11570,int_stack+9725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137193, 0.0, zero_stack, 1.0, int_stack+158533,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8600,int_stack+234978,int_stack+1675,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+13611,int_stack+13551, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96332, 1.0, int_stack+99810, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+13701,int_stack+13611, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96392, 1.0, int_stack+99870, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138093, 1.0, int_stack+149488, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+13827,int_stack+13701, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96482, 1.0, int_stack+99960, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138273, 1.0, int_stack+149668, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+235356,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138543, 1.0, int_stack+149938, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+14095,int_stack+13995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97220, 1.0, int_stack+100698, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+14245,int_stack+14095, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97320, 1.0, int_stack+100798, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138903, 1.0, int_stack+150298, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+14455,int_stack+14245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97470, 1.0, int_stack+100948, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1675,int_stack+131738,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139203, 1.0, int_stack+145338, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12200,int_stack+1675,int_stack+235728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140421, 1.0, int_stack+150598, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+223383,int_stack+12200,int_stack+11600,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+14885,int_stack+14735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42953, 1.0, int_stack+45937, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1675,int_stack+15110,int_stack+14885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43253, 1.0, int_stack+46237, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+234978,int_stack+1675,int_stack+11600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144888, 1.0, int_stack+158083, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13200,int_stack+15425,int_stack+15110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43703, 1.0, int_stack+46687, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14145,int_stack+13200,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141021, 1.0, int_stack+0, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+14145,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+148588, 1.0, int_stack+158533, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+240978,int_stack+239478,int_stack+12200,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+15905,int_stack+15845, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+99810, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239658,int_stack+15995,int_stack+15905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+99870, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239928,int_stack+239658,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+149488, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+240288,int_stack+16121,int_stack+15995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+99960, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+240288,int_stack+239658, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+149668, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+240288,int_stack+234978,int_stack+239928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+149938, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+16389,int_stack+16289, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+100698, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+16539,int_stack+16389, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+100798, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+150298, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+16749,int_stack+16539, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+100948, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1675,int_stack+131738,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+145338, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+1675,int_stack+235728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+150598, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+12600,int_stack+11600,int_stack+240288,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+17179,int_stack+17029, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+45937, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+234978,int_stack+17404,int_stack+17179, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+46237, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+158083, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+17719,int_stack+17404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+46687, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14400,int_stack+1675,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+14400,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+158533, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+14400,int_stack+234978,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+18199,int_stack+18139, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92854, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103288,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+18289,int_stack+18199, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92914, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103348,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159433,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+18415,int_stack+18289, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93004, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103438,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+11780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159613,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+235356,int_stack+12050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127010, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+159883,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+18683,int_stack+18583, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104176,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+18833,int_stack+18683, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93842, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104276,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+160243,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+19043,int_stack+18833, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93992, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104426,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+240078,int_stack+131738,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129788, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+160543,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+240078,int_stack+235728, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130238, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+161761,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+17400,int_stack+11600,int_stack+239478,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+19473,int_stack+19323, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39969, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48921,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239928,int_stack+19698,int_stack+19473, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40269, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49221,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+234978,int_stack+239928,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136068, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+166228,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+20013,int_stack+19698, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40719, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49671,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+243978,int_stack+1675,int_stack+239928, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+136518, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+162361,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+243978,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137193, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+169928,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+243978,int_stack+239478,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+20493,int_stack+20433, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96332, 0.0, zero_stack, 1.0, int_stack+103288, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+20583,int_stack+20493, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96392, 0.0, zero_stack, 1.0, int_stack+103348, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138093, 0.0, zero_stack, 1.0, int_stack+159433, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+20709,int_stack+20583, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96482, 0.0, zero_stack, 1.0, int_stack+103438, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+11780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138273, 0.0, zero_stack, 1.0, int_stack+159613, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239856,int_stack+12050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138543, 0.0, zero_stack, 1.0, int_stack+159883, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+20977,int_stack+20877, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97220, 0.0, zero_stack, 1.0, int_stack+104176, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239778,int_stack+21127,int_stack+20977, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97320, 0.0, zero_stack, 1.0, int_stack+104276, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+240228,int_stack+239778,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138903, 0.0, zero_stack, 1.0, int_stack+160243, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+21337,int_stack+21127, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+97470, 0.0, zero_stack, 1.0, int_stack+104426, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235578,int_stack+131738,int_stack+239778, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139203, 0.0, zero_stack, 1.0, int_stack+160543, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+235578,int_stack+240228, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140421, 0.0, zero_stack, 1.0, int_stack+161761, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+19200,int_stack+11600,int_stack+234978,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+21767,int_stack+21617, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42953, 0.0, zero_stack, 1.0, int_stack+48921, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235428,int_stack+21992,int_stack+21767, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43253, 0.0, zero_stack, 1.0, int_stack+49221, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+235428,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+144888, 0.0, zero_stack, 1.0, int_stack+166228, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+22307,int_stack+21992, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43703, 0.0, zero_stack, 1.0, int_stack+49671, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21000,int_stack+1675,int_stack+235428, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141021, 0.0, zero_stack, 1.0, int_stack+162361, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+21000,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+148588, 0.0, zero_stack, 1.0, int_stack+169928, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+246978,int_stack+234978,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+22787,int_stack+22727, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99810, 1.0, int_stack+103288, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+22877,int_stack+22787, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99870, 1.0, int_stack+103348, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149488, 1.0, int_stack+159433, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+23003,int_stack+22877, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+99960, 1.0, int_stack+103438, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+11780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149668, 1.0, int_stack+159613, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+235356,int_stack+12050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149938, 1.0, int_stack+159883, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+23271,int_stack+23171, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100698, 1.0, int_stack+104176, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+23421,int_stack+23271, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100798, 1.0, int_stack+104276, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150298, 1.0, int_stack+160243, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+23631,int_stack+23421, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+100948, 1.0, int_stack+104426, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+240078,int_stack+131738,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145338, 1.0, int_stack+160543, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+240078,int_stack+235728, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+150598, 1.0, int_stack+161761, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+21000,int_stack+11600,int_stack+239478,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+24061,int_stack+23911, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45937, 1.0, int_stack+48921, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239928,int_stack+24286,int_stack+24061, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46237, 1.0, int_stack+49221, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+234978,int_stack+239928,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+158083, 1.0, int_stack+166228, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22800,int_stack+24601,int_stack+24286, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46687, 1.0, int_stack+49671, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+249978,int_stack+22800,int_stack+239928, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+162361, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+249978,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+158533, 1.0, int_stack+169928, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+249978,int_stack+239478,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+25081,int_stack+25021, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+103288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+25171,int_stack+25081, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+103348, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+159433, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+25297,int_stack+25171, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+103438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+11780, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+159613, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239856,int_stack+12050, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+159883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+25565,int_stack+25465, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+104176, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239778,int_stack+25715,int_stack+25565, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+104276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+240228,int_stack+239778,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+160243, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+25925,int_stack+25715, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+104426, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235578,int_stack+131738,int_stack+239778, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+160543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+235578,int_stack+240228, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+161761, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+22800,int_stack+11600,int_stack+234978,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+26355,int_stack+26205, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+48921, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235428,int_stack+26580,int_stack+26355, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49221, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+235428,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+166228, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24600,int_stack+26895,int_stack+26580, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49671, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+25545,int_stack+24600,int_stack+235428, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+162361, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+25545,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+169928, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+252978,int_stack+234978,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+27375,int_stack+27315, 0.0, zero_stack, 1.0, int_stack+92854, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106766,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+27465,int_stack+27375, 0.0, zero_stack, 1.0, int_stack+92914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106826,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 1.0, int_stack+126560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170828,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+27591,int_stack+27465, 0.0, zero_stack, 1.0, int_stack+93004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106916,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+11780, 0.0, zero_stack, 1.0, int_stack+126740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171008,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+235356,int_stack+12050, 0.0, zero_stack, 1.0, int_stack+127010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171278,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+27859,int_stack+27759, 0.0, zero_stack, 1.0, int_stack+93742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107654,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+28009,int_stack+27859, 0.0, zero_stack, 1.0, int_stack+93842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107754,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 1.0, int_stack+129488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171638,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+28219,int_stack+28009, 0.0, zero_stack, 1.0, int_stack+93992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107904,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+240078,int_stack+131738,int_stack+235278, 0.0, zero_stack, 1.0, int_stack+129788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+166678,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+240078,int_stack+235728, 0.0, zero_stack, 1.0, int_stack+130238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171938,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+24600,int_stack+11600,int_stack+239478,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+28649,int_stack+28499, 0.0, zero_stack, 1.0, int_stack+39969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51905,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239928,int_stack+28874,int_stack+28649, 0.0, zero_stack, 1.0, int_stack+40269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52205,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+234978,int_stack+239928,int_stack+239478, 0.0, zero_stack, 1.0, int_stack+136068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152816,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26400,int_stack+29189,int_stack+28874, 0.0, zero_stack, 1.0, int_stack+40719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52655,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27345,int_stack+26400,int_stack+239928, 0.0, zero_stack, 1.0, int_stack+136518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178223,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+27345,int_stack+234978, 0.0, zero_stack, 1.0, int_stack+137193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178898,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+26400,int_stack+239478,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+29669,int_stack+29609, 0.0, zero_stack, 1.0, int_stack+96332, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106766, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+29759,int_stack+29669, 0.0, zero_stack, 1.0, int_stack+96392, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106826, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 1.0, int_stack+138093, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+170828, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+29885,int_stack+29759, 0.0, zero_stack, 1.0, int_stack+96482, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106916, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+11780, 0.0, zero_stack, 1.0, int_stack+138273, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171008, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239856,int_stack+12050, 0.0, zero_stack, 1.0, int_stack+138543, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171278, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+30153,int_stack+30053, 0.0, zero_stack, 1.0, int_stack+97220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107654, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239778,int_stack+30303,int_stack+30153, 0.0, zero_stack, 1.0, int_stack+97320, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107754, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+240228,int_stack+239778,int_stack+239478, 0.0, zero_stack, 1.0, int_stack+138903, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171638, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+30513,int_stack+30303, 0.0, zero_stack, 1.0, int_stack+97470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107904, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235578,int_stack+131738,int_stack+239778, 0.0, zero_stack, 1.0, int_stack+139203, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+166678, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+235578,int_stack+240228, 0.0, zero_stack, 1.0, int_stack+140421, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+171938, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+255978,int_stack+11600,int_stack+234978,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+30943,int_stack+30793, 0.0, zero_stack, 1.0, int_stack+42953, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51905, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235428,int_stack+31168,int_stack+30943, 0.0, zero_stack, 1.0, int_stack+43253, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52205, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+235428,int_stack+234978, 0.0, zero_stack, 1.0, int_stack+144888, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+152816, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29400,int_stack+31483,int_stack+31168, 0.0, zero_stack, 1.0, int_stack+43703, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52655, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+30345,int_stack+29400,int_stack+235428, 0.0, zero_stack, 1.0, int_stack+141021, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178223, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+30345,int_stack+239478, 0.0, zero_stack, 1.0, int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+178898, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+257778,int_stack+234978,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+31963,int_stack+31903, 0.0, zero_stack, 1.0, int_stack+99810, 0.0, zero_stack, 1.0, int_stack+106766, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+32053,int_stack+31963, 0.0, zero_stack, 1.0, int_stack+99870, 0.0, zero_stack, 1.0, int_stack+106826, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 1.0, int_stack+149488, 0.0, zero_stack, 1.0, int_stack+170828, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+32179,int_stack+32053, 0.0, zero_stack, 1.0, int_stack+99960, 0.0, zero_stack, 1.0, int_stack+106916, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+11780, 0.0, zero_stack, 1.0, int_stack+149668, 0.0, zero_stack, 1.0, int_stack+171008, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+235356,int_stack+12050, 0.0, zero_stack, 1.0, int_stack+149938, 0.0, zero_stack, 1.0, int_stack+171278, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+32447,int_stack+32347, 0.0, zero_stack, 1.0, int_stack+100698, 0.0, zero_stack, 1.0, int_stack+107654, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+32597,int_stack+32447, 0.0, zero_stack, 1.0, int_stack+100798, 0.0, zero_stack, 1.0, int_stack+107754, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 1.0, int_stack+150298, 0.0, zero_stack, 1.0, int_stack+171638, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+32807,int_stack+32597, 0.0, zero_stack, 1.0, int_stack+100948, 0.0, zero_stack, 1.0, int_stack+107904, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+240078,int_stack+131738,int_stack+235278, 0.0, zero_stack, 1.0, int_stack+145338, 0.0, zero_stack, 1.0, int_stack+166678, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+240078,int_stack+235728, 0.0, zero_stack, 1.0, int_stack+150598, 0.0, zero_stack, 1.0, int_stack+171938, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+29400,int_stack+11600,int_stack+239478,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+33237,int_stack+33087, 0.0, zero_stack, 1.0, int_stack+45937, 0.0, zero_stack, 1.0, int_stack+51905, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239928,int_stack+33462,int_stack+33237, 0.0, zero_stack, 1.0, int_stack+46237, 0.0, zero_stack, 1.0, int_stack+52205, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+234978,int_stack+239928,int_stack+239478, 0.0, zero_stack, 1.0, int_stack+158083, 0.0, zero_stack, 1.0, int_stack+152816, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+31200,int_stack+33777,int_stack+33462, 0.0, zero_stack, 1.0, int_stack+46687, 0.0, zero_stack, 1.0, int_stack+52655, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+32145,int_stack+31200,int_stack+239928, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+178223, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+32145,int_stack+234978, 0.0, zero_stack, 1.0, int_stack+158533, 0.0, zero_stack, 1.0, int_stack+178898, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+260778,int_stack+239478,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+34257,int_stack+34197, 0.0, zero_stack, 1.0, int_stack+103288, 1.0, int_stack+106766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+34347,int_stack+34257, 0.0, zero_stack, 1.0, int_stack+103348, 1.0, int_stack+106826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 1.0, int_stack+159433, 1.0, int_stack+170828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+34473,int_stack+34347, 0.0, zero_stack, 1.0, int_stack+103438, 1.0, int_stack+106916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+11780, 0.0, zero_stack, 1.0, int_stack+159613, 1.0, int_stack+171008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239856,int_stack+12050, 0.0, zero_stack, 1.0, int_stack+159883, 1.0, int_stack+171278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+34741,int_stack+34641, 0.0, zero_stack, 1.0, int_stack+104176, 1.0, int_stack+107654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239778,int_stack+34891,int_stack+34741, 0.0, zero_stack, 1.0, int_stack+104276, 1.0, int_stack+107754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+240228,int_stack+239778,int_stack+239478, 0.0, zero_stack, 1.0, int_stack+160243, 1.0, int_stack+171638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+35101,int_stack+34891, 0.0, zero_stack, 1.0, int_stack+104426, 1.0, int_stack+107904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235578,int_stack+131738,int_stack+239778, 0.0, zero_stack, 1.0, int_stack+160543, 1.0, int_stack+166678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+235578,int_stack+240228, 0.0, zero_stack, 1.0, int_stack+161761, 1.0, int_stack+171938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+31200,int_stack+11600,int_stack+234978,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+35531,int_stack+35381, 0.0, zero_stack, 1.0, int_stack+48921, 1.0, int_stack+51905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235428,int_stack+35756,int_stack+35531, 0.0, zero_stack, 1.0, int_stack+49221, 1.0, int_stack+52205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+235428,int_stack+234978, 0.0, zero_stack, 1.0, int_stack+166228, 1.0, int_stack+152816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+33000,int_stack+36071,int_stack+35756, 0.0, zero_stack, 1.0, int_stack+49671, 1.0, int_stack+52655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33945,int_stack+33000,int_stack+235428, 0.0, zero_stack, 1.0, int_stack+162361, 1.0, int_stack+178223, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+33945,int_stack+239478, 0.0, zero_stack, 1.0, int_stack+169928, 1.0, int_stack+178898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+33000,int_stack+234978,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+36551,int_stack+36491, 0.0, zero_stack, 2.0, int_stack+106766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+36641,int_stack+36551, 0.0, zero_stack, 2.0, int_stack+106826, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 0.0, zero_stack, 2.0, int_stack+170828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+36767,int_stack+36641, 0.0, zero_stack, 2.0, int_stack+106916, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+11780, 0.0, zero_stack, 2.0, int_stack+171008, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+235356,int_stack+12050, 0.0, zero_stack, 2.0, int_stack+171278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+37035,int_stack+36935, 0.0, zero_stack, 2.0, int_stack+107654, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+37185,int_stack+37035, 0.0, zero_stack, 2.0, int_stack+107754, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 2.0, int_stack+171638, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+37395,int_stack+37185, 0.0, zero_stack, 2.0, int_stack+107904, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+240078,int_stack+131738,int_stack+235278, 0.0, zero_stack, 2.0, int_stack+166678, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+240078,int_stack+235728, 0.0, zero_stack, 2.0, int_stack+171938, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+263778,int_stack+11600,int_stack+239478,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+37825,int_stack+37675, 0.0, zero_stack, 2.0, int_stack+51905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239928,int_stack+38050,int_stack+37825, 0.0, zero_stack, 2.0, int_stack+52205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+234978,int_stack+239928,int_stack+239478, 0.0, zero_stack, 2.0, int_stack+152816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36000,int_stack+38365,int_stack+38050, 0.0, zero_stack, 2.0, int_stack+52655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+36945,int_stack+36000,int_stack+239928, 0.0, zero_stack, 2.0, int_stack+178223, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+36945,int_stack+234978, 0.0, zero_stack, 2.0, int_stack+178898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+265578,int_stack+239478,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+38845,int_stack+38785, 1.0, int_stack+92854, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110244,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+38935,int_stack+38845, 1.0, int_stack+92914, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110304,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 1.0, int_stack+126560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+179798,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+39061,int_stack+38935, 1.0, int_stack+93004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110394,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+11780, 1.0, int_stack+126740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+179978,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239856,int_stack+12050, 1.0, int_stack+127010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180248,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+39329,int_stack+39229, 1.0, int_stack+93742, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111232,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239778,int_stack+39479,int_stack+39329, 1.0, int_stack+93842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111482,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+240228,int_stack+239778,int_stack+239478, 1.0, int_stack+129488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180608,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+39689,int_stack+39479, 1.0, int_stack+93992, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111842,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235578,int_stack+131738,int_stack+239778, 1.0, int_stack+129788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180908,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+235578,int_stack+240228, 1.0, int_stack+130238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182126,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+36000,int_stack+11600,int_stack+234978,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+40494,int_stack+40119, 1.0, int_stack+39969, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55039,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235428,int_stack+41034,int_stack+40494, 1.0, int_stack+40269, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55564,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+235428,int_stack+234978, 1.0, int_stack+136068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186593,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+41349,int_stack+41034, 1.0, int_stack+40719, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56014,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+37800,int_stack+126425,int_stack+235428, 1.0, int_stack+136518, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182726,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+37800,int_stack+239478, 1.0, int_stack+137193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+190293,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+37800,int_stack+234978,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+41829,int_stack+41769, 1.0, int_stack+96332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110244, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+41919,int_stack+41829, 1.0, int_stack+96392, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110304, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 1.0, int_stack+138093, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+179798, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+42045,int_stack+41919, 1.0, int_stack+96482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110394, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+11780, 1.0, int_stack+138273, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+179978, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+235356,int_stack+12050, 1.0, int_stack+138543, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180248, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+42313,int_stack+42213, 1.0, int_stack+97220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111232, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+42463,int_stack+42313, 1.0, int_stack+97320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111482, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 1.0, int_stack+138903, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180608, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+42673,int_stack+42463, 1.0, int_stack+97470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111842, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+240078,int_stack+131738,int_stack+235278, 1.0, int_stack+139203, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180908, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+240078,int_stack+235728, 1.0, int_stack+140421, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182126, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+40800,int_stack+11600,int_stack+239478,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+43478,int_stack+43103, 1.0, int_stack+42953, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55039, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+239928,int_stack+44018,int_stack+43478, 1.0, int_stack+43253, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55564, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+234978,int_stack+239928,int_stack+239478, 1.0, int_stack+144888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186593, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+44333,int_stack+44018, 1.0, int_stack+43703, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56014, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+42600,int_stack+126425,int_stack+239928, 1.0, int_stack+141021, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182726, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+42600,int_stack+234978, 1.0, int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+190293, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+135828,int_stack+239478,int_stack+11600,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+44813,int_stack+44753, 1.0, int_stack+99810, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110244, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+44903,int_stack+44813, 1.0, int_stack+99870, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110304, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600, 1.0, int_stack+149488, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+179798, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+45029,int_stack+44903, 1.0, int_stack+99960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110394, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+11780, 1.0, int_stack+149668, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+179978, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239856,int_stack+12050, 1.0, int_stack+149938, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180248, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+45297,int_stack+45197, 1.0, int_stack+100698, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111232, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+144888,int_stack+45447,int_stack+45297, 1.0, int_stack+100798, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111482, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239778,int_stack+144888,int_stack+239478, 1.0, int_stack+150298, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180608, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+45657,int_stack+45447, 1.0, int_stack+100948, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+111842, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11600,int_stack+131738,int_stack+144888, 1.0, int_stack+145338, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+180908, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+11600,int_stack+239778, 1.0, int_stack+150598, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182126, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+149188,int_stack+1675,int_stack+148588,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+46462,int_stack+46087, 1.0, int_stack+45937, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55039, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11600,int_stack+47002,int_stack+46462, 1.0, int_stack+46237, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55564, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+144888,int_stack+11600,int_stack+148588, 1.0, int_stack+158083, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+186593, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+47317,int_stack+47002, 1.0, int_stack+46687, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56014, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+126425,int_stack+11600, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+182726, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+144888, 1.0, int_stack+158533, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+190293, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+42600,int_stack+234978,int_stack+1675,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+47797,int_stack+47737, 1.0, int_stack+103288, 0.0, zero_stack, 1.0, int_stack+110244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+47887,int_stack+47797, 1.0, int_stack+103348, 0.0, zero_stack, 1.0, int_stack+110304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 1.0, int_stack+159433, 0.0, zero_stack, 1.0, int_stack+179798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+48013,int_stack+47887, 1.0, int_stack+103438, 0.0, zero_stack, 1.0, int_stack+110394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 1.0, int_stack+159613, 0.0, zero_stack, 1.0, int_stack+179978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 1.0, int_stack+159883, 0.0, zero_stack, 1.0, int_stack+180248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+48281,int_stack+48181, 1.0, int_stack+104176, 0.0, zero_stack, 1.0, int_stack+111232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+48431,int_stack+48281, 1.0, int_stack+104276, 0.0, zero_stack, 1.0, int_stack+111482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 1.0, int_stack+160243, 0.0, zero_stack, 1.0, int_stack+180608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+48641,int_stack+48431, 1.0, int_stack+104426, 0.0, zero_stack, 1.0, int_stack+111842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+144888,int_stack+131738,int_stack+235278, 1.0, int_stack+160543, 0.0, zero_stack, 1.0, int_stack+180908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+144888,int_stack+235728, 1.0, int_stack+161761, 0.0, zero_stack, 1.0, int_stack+182126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+158083,int_stack+1675,int_stack+148588,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+49446,int_stack+49071, 1.0, int_stack+48921, 0.0, zero_stack, 1.0, int_stack+55039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+49986,int_stack+49446, 1.0, int_stack+49221, 0.0, zero_stack, 1.0, int_stack+55564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+144888,int_stack+0,int_stack+148588, 1.0, int_stack+166228, 0.0, zero_stack, 1.0, int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+50301,int_stack+49986, 1.0, int_stack+49671, 0.0, zero_stack, 1.0, int_stack+56014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+126425,int_stack+0, 1.0, int_stack+162361, 0.0, zero_stack, 1.0, int_stack+182726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+144888, 1.0, int_stack+169928, 0.0, zero_stack, 1.0, int_stack+190293, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+45600,int_stack+239478,int_stack+1675,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+50781,int_stack+50721, 1.0, int_stack+106766, 1.0, int_stack+110244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+50871,int_stack+50781, 1.0, int_stack+106826, 1.0, int_stack+110304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 1.0, int_stack+170828, 1.0, int_stack+179798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+50997,int_stack+50871, 1.0, int_stack+106916, 1.0, int_stack+110394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+1855, 1.0, int_stack+171008, 1.0, int_stack+179978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239856,int_stack+2125, 1.0, int_stack+171278, 1.0, int_stack+180248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+51265,int_stack+51165, 1.0, int_stack+107654, 1.0, int_stack+111232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+166228,int_stack+51415,int_stack+51265, 1.0, int_stack+107754, 1.0, int_stack+111482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239778,int_stack+166228,int_stack+239478, 1.0, int_stack+171638, 1.0, int_stack+180608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+51625,int_stack+51415, 1.0, int_stack+107904, 1.0, int_stack+111842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+144888,int_stack+131738,int_stack+166228, 1.0, int_stack+166678, 1.0, int_stack+180908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+144888,int_stack+239778, 1.0, int_stack+171938, 1.0, int_stack+182126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+169928,int_stack+1675,int_stack+148588,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+52430,int_stack+52055, 1.0, int_stack+51905, 1.0, int_stack+55039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+52970,int_stack+52430, 1.0, int_stack+52205, 1.0, int_stack+55564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+144888,int_stack+0,int_stack+148588, 1.0, int_stack+152816, 1.0, int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+53285,int_stack+52970, 1.0, int_stack+52655, 1.0, int_stack+56014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+126425,int_stack+0, 1.0, int_stack+178223, 1.0, int_stack+182726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+144888, 1.0, int_stack+178898, 1.0, int_stack+190293, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+48600,int_stack+234978,int_stack+1675,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+53765,int_stack+53705, 2.0, int_stack+110244, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+53855,int_stack+53765, 2.0, int_stack+110304, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 2.0, int_stack+179798, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+53981,int_stack+53855, 2.0, int_stack+110394, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 2.0, int_stack+179978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 2.0, int_stack+180248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+54249,int_stack+54149, 2.0, int_stack+111232, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+54399,int_stack+54249, 2.0, int_stack+111482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 2.0, int_stack+180608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+54609,int_stack+54399, 2.0, int_stack+111842, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+144888,int_stack+131738,int_stack+235278, 2.0, int_stack+180908, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+144888,int_stack+235728, 2.0, int_stack+182126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+178223,int_stack+1675,int_stack+148588,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+55789,int_stack+55189, 2.0, int_stack+55039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+56329,int_stack+55789, 2.0, int_stack+55564, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+144888,int_stack+0,int_stack+148588, 2.0, int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+56644,int_stack+56329, 2.0, int_stack+56014, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+126425,int_stack+0, 2.0, int_stack+182726, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+144888, 2.0, int_stack+190293, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+51600,int_stack+239478,int_stack+1675,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+57124,int_stack+57064, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114872,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+57214,int_stack+57124, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114932,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128888,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+57340,int_stack+57214, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115022,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129068,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239856,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135168,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+57608,int_stack+57508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115760,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+186593,int_stack+57758,int_stack+57608, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115860,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239778,int_stack+186593,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130838,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+57968,int_stack+57758, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116010,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+190293,int_stack+131738,int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187043,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+190293,int_stack+239778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131138,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+54600,int_stack+1675,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+186593,int_stack+58398,int_stack+58248, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69718,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+58623,int_stack+58398, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70018,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+190293,int_stack+0,int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195993,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+58938,int_stack+58623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70468,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+126425,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125750,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+190293, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+199693,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+268578,int_stack+234978,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+132368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+59418,int_stack+59358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114872, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+59508,int_stack+59418, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114932, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128888, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+59634,int_stack+59508, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115022, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129068, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135168, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+59902,int_stack+59802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115760, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+186593,int_stack+60052,int_stack+59902, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115860, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235278,int_stack+186593,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130838, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+60262,int_stack+60052, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116010, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+190293,int_stack+131738,int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187043, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+190293,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131138, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+56400,int_stack+1675,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139821, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+186593,int_stack+60692,int_stack+60542, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69718, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+60917,int_stack+60692, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70018, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+190293,int_stack+0,int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195993, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+61232,int_stack+60917, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70468, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+126425,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125750, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+190293, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+199693, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+58200,int_stack+239478,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+61712,int_stack+61652, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114872, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+61802,int_stack+61712, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114932, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128888, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+61928,int_stack+61802, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115022, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129068, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239856,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135168, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+62196,int_stack+62096, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115760, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+186593,int_stack+62346,int_stack+62196, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115860, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239778,int_stack+186593,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130838, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+62556,int_stack+62346, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116010, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+190293,int_stack+131738,int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187043, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+190293,int_stack+239778, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131138, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+271578,int_stack+1675,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+186593,int_stack+62986,int_stack+62836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69718, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+63211,int_stack+62986, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70018, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+190293,int_stack+0,int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195993, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+63526,int_stack+63211, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70468, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+126425,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125750, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+190293, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+199693, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+273378,int_stack+234978,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+64006,int_stack+63946, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+64096,int_stack+64006, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+64222,int_stack+64096, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+129068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+64490,int_stack+64390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+186593,int_stack+64640,int_stack+64490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235278,int_stack+186593,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+64850,int_stack+64640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+190293,int_stack+131738,int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187043, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+190293,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+61200,int_stack+1675,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+161161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+186593,int_stack+65280,int_stack+65130, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+65505,int_stack+65280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+190293,int_stack+0,int_stack+186593, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+65820,int_stack+65505, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+70468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+126425,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+190293, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+199693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+63000,int_stack+239478,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+66300,int_stack+66240, 0.0, zero_stack, 1.0, int_stack+114872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+66390,int_stack+66300, 0.0, zero_stack, 1.0, int_stack+114932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+128888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+66516,int_stack+66390, 0.0, zero_stack, 1.0, int_stack+115022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239856,int_stack+239478,int_stack+1855, 0.0, zero_stack, 1.0, int_stack+129068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239856,int_stack+2125, 0.0, zero_stack, 1.0, int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+239478,int_stack+66784,int_stack+66684, 0.0, zero_stack, 1.0, int_stack+115760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+186593,int_stack+66934,int_stack+66784, 0.0, zero_stack, 1.0, int_stack+115860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239778,int_stack+186593,int_stack+239478, 0.0, zero_stack, 1.0, int_stack+130838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+67144,int_stack+66934, 0.0, zero_stack, 1.0, int_stack+116010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+190293,int_stack+131738,int_stack+186593, 0.0, zero_stack, 1.0, int_stack+187043, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+190293,int_stack+239778, 0.0, zero_stack, 1.0, int_stack+131138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+276378,int_stack+1675,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+177623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+186593,int_stack+67574,int_stack+67424, 0.0, zero_stack, 1.0, int_stack+69718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+67799,int_stack+67574, 0.0, zero_stack, 1.0, int_stack+70018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+190293,int_stack+0,int_stack+186593, 0.0, zero_stack, 1.0, int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+68114,int_stack+67799, 0.0, zero_stack, 1.0, int_stack+70468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+126425,int_stack+0, 0.0, zero_stack, 1.0, int_stack+125750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+190293, 0.0, zero_stack, 1.0, int_stack+199693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+278178,int_stack+234978,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+68594,int_stack+68534, 1.0, int_stack+114872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+68684,int_stack+68594, 1.0, int_stack+114932, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 1.0, int_stack+128888, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+68810,int_stack+68684, 1.0, int_stack+115022, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 1.0, int_stack+129068, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 1.0, int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+135168,int_stack+69078,int_stack+68978, 1.0, int_stack+115760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+186593,int_stack+69228,int_stack+69078, 1.0, int_stack+115860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+234978,int_stack+186593,int_stack+135168, 1.0, int_stack+130838, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+131738,int_stack+69438,int_stack+69228, 1.0, int_stack+116010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235578,int_stack+131738,int_stack+186593, 1.0, int_stack+187043, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+235578,int_stack+234978, 1.0, int_stack+131138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+128888,int_stack+1675,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+181526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+70243,int_stack+69868, 1.0, int_stack+69718, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+70783,int_stack+70243, 1.0, int_stack+70018, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+186593,int_stack+0,int_stack+148588, 1.0, int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126425,int_stack+71098,int_stack+70783, 1.0, int_stack+70468, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+126425,int_stack+0, 1.0, int_stack+125750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+186593, 1.0, int_stack+199693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+66000,int_stack+239478,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+187493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+71578,int_stack+71518,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+71668,int_stack+71578,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+1855,int_stack+1675,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2125,int_stack+71794,int_stack+71668,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+2125,int_stack+1855,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239478,int_stack+135168,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+135168,int_stack+72062,int_stack+71962,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+195993,int_stack+72212,int_stack+72062,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+195993,int_stack+135168,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+240078,int_stack+72422,int_stack+72212,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+199693,int_stack+240078,int_stack+195993,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+199693,int_stack+239478,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+69000,int_stack+1675,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+127370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+195993,int_stack+72852,int_stack+72702,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+73077,int_stack+72852,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+199693,int_stack+0,int_stack+195993,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+73392,int_stack+73077,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+0,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+199693,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+70800,int_stack+239478,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+196893, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+73872,int_stack+73812, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118350,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+73962,int_stack+73872, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118410,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+200593,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2125,int_stack+74088,int_stack+73962, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118500,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+2125,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+200773,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239478,int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201043,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+135168,int_stack+74356,int_stack+74256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119238,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195993,int_stack+74506,int_stack+74356, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119338,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+195993,int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135528,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+240078,int_stack+74716,int_stack+74506, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119488,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+199693,int_stack+240078,int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+196443,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+199693,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202003,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+281178,int_stack+1675,int_stack+148588, 0.0, zero_stack, 1.0, int_stack+128288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+195993,int_stack+75146,int_stack+74996, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86466,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+75371,int_stack+75146, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86766,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+199693,int_stack+0,int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202603,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+75686,int_stack+75371, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87216,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203053,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+199693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+210188,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+282978,int_stack+239478,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+132368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+76166,int_stack+76106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118350, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+76256,int_stack+76166, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118410, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+200593, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2125,int_stack+76382,int_stack+76256, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118500, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+2125,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+200773, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239478,int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201043, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+135168,int_stack+76650,int_stack+76550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119238, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195993,int_stack+76800,int_stack+76650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119338, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+195993,int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135528, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+240078,int_stack+77010,int_stack+76800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119488, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+199693,int_stack+240078,int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+196443, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+199693,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202003, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+73800,int_stack+1675,int_stack+148588, 0.0, zero_stack, 1.0, int_stack+139821, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+195993,int_stack+77440,int_stack+77290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86466, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+77665,int_stack+77440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86766, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+199693,int_stack+0,int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202603, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+77980,int_stack+77665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87216, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203053, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+199693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+210188, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+285978,int_stack+239478,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+145788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+78460,int_stack+78400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118350, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+78550,int_stack+78460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118410, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+200593, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2125,int_stack+78676,int_stack+78550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118500, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+2125,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+200773, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239478,int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201043, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+135168,int_stack+78944,int_stack+78844, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119238, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195993,int_stack+79094,int_stack+78944, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119338, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+195993,int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135528, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+240078,int_stack+79304,int_stack+79094, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119488, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+199693,int_stack+240078,int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+196443, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+199693,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202003, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+75600,int_stack+1675,int_stack+148588, 0.0, zero_stack, 1.0, int_stack+151216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+195993,int_stack+79734,int_stack+79584, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86466, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+79959,int_stack+79734, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86766, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+199693,int_stack+0,int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202603, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+80274,int_stack+79959, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87216, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203053, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+199693, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+210188, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+77400,int_stack+239478,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+151816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+80754,int_stack+80694, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+80844,int_stack+80754, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+200593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2125,int_stack+80970,int_stack+80844, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+118500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+2125,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+200773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239478,int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+201043, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+135168,int_stack+81238,int_stack+81138, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195993,int_stack+81388,int_stack+81238, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+195993,int_stack+135168, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+240078,int_stack+81598,int_stack+81388, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+199693,int_stack+240078,int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+196443, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+199693,int_stack+239478, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+288978,int_stack+1675,int_stack+148588, 0.0, zero_stack, 1.0, int_stack+161161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+195993,int_stack+82028,int_stack+81878, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86466, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+82253,int_stack+82028, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+86766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+199693,int_stack+0,int_stack+195993, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+202603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+82568,int_stack+82253, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+203053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+199693, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+210188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+290778,int_stack+239478,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+167128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+83048,int_stack+82988, 0.0, zero_stack, 1.0, int_stack+118350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+83138,int_stack+83048, 0.0, zero_stack, 1.0, int_stack+118410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+1855,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+200593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2125,int_stack+83264,int_stack+83138, 0.0, zero_stack, 1.0, int_stack+118500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+2125,int_stack+1855, 0.0, zero_stack, 1.0, int_stack+200773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239478,int_stack+135168, 0.0, zero_stack, 1.0, int_stack+201043, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+135168,int_stack+83532,int_stack+83432, 0.0, zero_stack, 1.0, int_stack+119238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195993,int_stack+83682,int_stack+83532, 0.0, zero_stack, 1.0, int_stack+119338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+195993,int_stack+135168, 0.0, zero_stack, 1.0, int_stack+135528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+240078,int_stack+83892,int_stack+83682, 0.0, zero_stack, 1.0, int_stack+119488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+199693,int_stack+240078,int_stack+195993, 0.0, zero_stack, 1.0, int_stack+196443, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+199693,int_stack+239478, 0.0, zero_stack, 1.0, int_stack+202003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+80400,int_stack+1675,int_stack+148588, 0.0, zero_stack, 1.0, int_stack+177623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+195993,int_stack+84322,int_stack+84172, 0.0, zero_stack, 1.0, int_stack+86466, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+84547,int_stack+84322, 0.0, zero_stack, 1.0, int_stack+86766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+199693,int_stack+0,int_stack+195993, 0.0, zero_stack, 1.0, int_stack+202603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+239478,int_stack+84862,int_stack+84547, 0.0, zero_stack, 1.0, int_stack+87216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+0, 0.0, zero_stack, 1.0, int_stack+203053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+199693, 0.0, zero_stack, 1.0, int_stack+210188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+82200,int_stack+239478,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+85342,int_stack+85282, 1.0, int_stack+118350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+85432,int_stack+85342, 1.0, int_stack+118410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+135168,int_stack+1855,int_stack+1675, 1.0, int_stack+200593, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2125,int_stack+85558,int_stack+85432, 1.0, int_stack+118500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+2125,int_stack+1855, 1.0, int_stack+200773, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+239478,int_stack+135168, 1.0, int_stack+201043, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+135168,int_stack+85826,int_stack+85726, 1.0, int_stack+119238, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+195993,int_stack+85976,int_stack+85826, 1.0, int_stack+119338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+239478,int_stack+195993,int_stack+135168, 1.0, int_stack+135528, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+135168,int_stack+86186,int_stack+85976, 1.0, int_stack+119488, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+240078,int_stack+135168,int_stack+195993, 1.0, int_stack+196443, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+240078,int_stack+239478, 1.0, int_stack+202003, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+293778,int_stack+1675,int_stack+148588, 0.0, zero_stack, 1.0, int_stack+181526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+86991,int_stack+86616, 1.0, int_stack+86466, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+87531,int_stack+86991, 1.0, int_stack+86766, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+195993,int_stack+0,int_stack+148588, 1.0, int_stack+202603, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+202003,int_stack+87846,int_stack+87531, 1.0, int_stack+87216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+202003,int_stack+0, 1.0, int_stack+203053, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+195993, 1.0, int_stack+210188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+85200,int_stack+234978,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+187493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+88326,int_stack+88266,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+88416,int_stack+88326,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+88542,int_stack+88416,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+88810,int_stack+88710,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+88960,int_stack+88810,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+89170,int_stack+88960,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+210188,int_stack+1675,int_stack+235278,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+210188,int_stack+235728,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+202003,int_stack+1675,int_stack+148588, 0.0, zero_stack, 1.0, int_stack+127370, 1.0, int_stack+201403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+89600,int_stack+89450,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+89825,int_stack+89600,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+210188,int_stack+0,int_stack+148588,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+90140,int_stack+89825,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+0,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+210188,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+295578,int_stack+234978,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+196893, 1.0, int_stack+207388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+90620,int_stack+90560,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+90710,int_stack+90620,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+90836,int_stack+90710,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+91104,int_stack+91004,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+91254,int_stack+91104,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+91464,int_stack+91254,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+210188,int_stack+1675,int_stack+235278,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+210188,int_stack+235728,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+88200,int_stack+1675,int_stack+148588, 0.0, zero_stack, 2.0, int_stack+201403, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+91894,int_stack+91744,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+92119,int_stack+91894,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+210188,int_stack+0,int_stack+148588,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+92434,int_stack+92119,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+0,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+210188,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+90000,int_stack+234978,int_stack+1675, 0.0, zero_stack, 2.0, int_stack+207388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+93358,int_stack+93298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122272,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+93448,int_stack+93358, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122332,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211088,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+93574,int_stack+93448, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122422,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211268,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211538,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+94582,int_stack+94482, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123160,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+94732,int_stack+94582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123260,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211898,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+94942,int_stack+94732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123410,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+210188,int_stack+1675,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+212198,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+210188,int_stack+235728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+213416,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+93000,int_stack+1675,int_stack+148588, 1.0, int_stack+128288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+95372,int_stack+95222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113072,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+95597,int_stack+95372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113372,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+210188,int_stack+0,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2675,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+95912,int_stack+95597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113822,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3125,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+210188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220683,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+298578,int_stack+234978,int_stack+1675, 1.0, int_stack+132368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+96836,int_stack+96776, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122272, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+96926,int_stack+96836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122332, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211088, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+97052,int_stack+96926, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122422, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211268, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211538, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+98060,int_stack+97960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123160, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+98210,int_stack+98060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123260, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211898, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+98420,int_stack+98210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123410, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+210188,int_stack+1675,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+212198, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+210188,int_stack+235728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+213416, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+130688,int_stack+1675,int_stack+148588, 1.0, int_stack+139821, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+98850,int_stack+98700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113072, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+99075,int_stack+98850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113372, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+210188,int_stack+0,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2675, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+99390,int_stack+99075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113822, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3125, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+210188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220683, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+94800,int_stack+234978,int_stack+1675, 1.0, int_stack+145788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+100314,int_stack+100254, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122272, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+100404,int_stack+100314, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122332, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211088, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+100530,int_stack+100404, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122422, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211268, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211538, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+101538,int_stack+101438, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123160, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+101688,int_stack+101538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123260, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211898, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+101898,int_stack+101688, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123410, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+210188,int_stack+1675,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+212198, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+210188,int_stack+235728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+213416, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+144888,int_stack+1675,int_stack+148588, 1.0, int_stack+151216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+102328,int_stack+102178, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113072, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+102553,int_stack+102328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113372, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+210188,int_stack+0,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2675, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+102868,int_stack+102553, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113822, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3125, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+210188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220683, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+97800,int_stack+234978,int_stack+1675, 1.0, int_stack+151816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+103792,int_stack+103732, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+103882,int_stack+103792, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+104008,int_stack+103882, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+105016,int_stack+104916, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+105166,int_stack+105016, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+211898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+105376,int_stack+105166, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+210188,int_stack+1675,int_stack+235278, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+212198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+210188,int_stack+235728, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+213416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+150988,int_stack+1675,int_stack+148588, 1.0, int_stack+161161, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+148588,int_stack+105806,int_stack+105656, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+106031,int_stack+105806, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+210188,int_stack+0,int_stack+148588, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+106346,int_stack+106031, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+210188, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+220683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+159883,int_stack+234978,int_stack+1675, 1.0, int_stack+167128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1675,int_stack+107270,int_stack+107210, 0.0, zero_stack, 1.0, int_stack+122272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1855,int_stack+107360,int_stack+107270, 0.0, zero_stack, 1.0, int_stack+122332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2125,int_stack+1855,int_stack+1675, 0.0, zero_stack, 1.0, int_stack+211088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+107486,int_stack+107360, 0.0, zero_stack, 1.0, int_stack+122422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+1855, 0.0, zero_stack, 1.0, int_stack+211268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148588,int_stack+235356,int_stack+2125, 0.0, zero_stack, 1.0, int_stack+211538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+108494,int_stack+108394, 0.0, zero_stack, 1.0, int_stack+123160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+108644,int_stack+108494, 0.0, zero_stack, 1.0, int_stack+123260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235728,int_stack+235278,int_stack+234978, 0.0, zero_stack, 1.0, int_stack+211898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1675,int_stack+108854,int_stack+108644, 0.0, zero_stack, 1.0, int_stack+123410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+210188,int_stack+1675,int_stack+235278, 0.0, zero_stack, 1.0, int_stack+212198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1675,int_stack+210188,int_stack+235728, 0.0, zero_stack, 1.0, int_stack+213416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+166228,int_stack+1675,int_stack+148588, 1.0, int_stack+177623, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+177623,int_stack+109284,int_stack+109134, 0.0, zero_stack, 1.0, int_stack+113072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+109509,int_stack+109284, 0.0, zero_stack, 1.0, int_stack+113372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+210188,int_stack+0,int_stack+177623, 0.0, zero_stack, 1.0, int_stack+2675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+109824,int_stack+109509, 0.0, zero_stack, 1.0, int_stack+113822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+0, 0.0, zero_stack, 1.0, int_stack+3125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+210188, 0.0, zero_stack, 1.0, int_stack+220683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+100800,int_stack+234978,int_stack+1675, 1.0, int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+110748,int_stack+110688, 1.0, int_stack+122272, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+235158,int_stack+110838,int_stack+110748, 1.0, int_stack+122332, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+235428,int_stack+235158,int_stack+234978, 1.0, int_stack+211088, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+235788,int_stack+110964,int_stack+110838, 1.0, int_stack+122422, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+210188,int_stack+235788,int_stack+235158, 1.0, int_stack+211268, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+177623,int_stack+210188,int_stack+235428, 1.0, int_stack+211538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+210188,int_stack+112432,int_stack+112332, 1.0, int_stack+123160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+210488,int_stack+112582,int_stack+112432, 1.0, int_stack+123260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+148588,int_stack+210488,int_stack+210188, 1.0, int_stack+211898, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+210938,int_stack+112792,int_stack+112582, 1.0, int_stack+123410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+195993,int_stack+210938,int_stack+210488, 1.0, int_stack+212198, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+195993,int_stack+148588, 1.0, int_stack+213416, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+210188,int_stack+11600,int_stack+177623, 1.0, int_stack+181526, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+177623,int_stack+113597,int_stack+113222, 1.0, int_stack+113072, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+213416,int_stack+114137,int_stack+113597, 1.0, int_stack+113372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+195993,int_stack+213416,int_stack+177623, 1.0, int_stack+2675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+114452,int_stack+114137, 1.0, int_stack+113822, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+213416, 1.0, int_stack+3125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+195993, 1.0, int_stack+220683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+234978,int_stack+11600, 1.0, int_stack+187493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+115376,int_stack+115316,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+115466,int_stack+115376,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+115592,int_stack+115466,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+11780,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+177623,int_stack+235356,int_stack+12050,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+116600,int_stack+116500,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+116750,int_stack+116600,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+148588,int_stack+235278,int_stack+234978,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+235728,int_stack+116960,int_stack+116750,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+220683,int_stack+235728,int_stack+235278,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+220683,int_stack+148588,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+186593,int_stack+11600,int_stack+177623, 1.0, int_stack+127370, 0.0, zero_stack, 1.0, int_stack+212816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+177623,int_stack+117390,int_stack+117240,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+220683,int_stack+117615,int_stack+117390,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+195993,int_stack+220683,int_stack+177623,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+117930,int_stack+117615,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+220683,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+195993,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+125750,int_stack+234978,int_stack+11600, 1.0, int_stack+196893, 0.0, zero_stack, 1.0, int_stack+217883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+118854,int_stack+118794,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11780,int_stack+118944,int_stack+118854,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12050,int_stack+11780,int_stack+11600,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+119070,int_stack+118944,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+235356,int_stack+234978,int_stack+11780,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+177623,int_stack+235356,int_stack+12050,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+234978,int_stack+120078,int_stack+119978,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+235278,int_stack+120228,int_stack+120078,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+148588,int_stack+235278,int_stack+234978,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+235728,int_stack+120438,int_stack+120228,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+220683,int_stack+235728,int_stack+235278,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+11600,int_stack+220683,int_stack+148588,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+195993,int_stack+11600,int_stack+177623, 1.0, int_stack+201403, 1.0, int_stack+212816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+177623,int_stack+120868,int_stack+120718,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+220683,int_stack+121093,int_stack+120868,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+190293,int_stack+220683,int_stack+177623,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+234978,int_stack+121408,int_stack+121093,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+220683,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+234978,int_stack+239478,int_stack+190293,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+180023,int_stack+234978,int_stack+11600, 1.0, int_stack+207388, 1.0, int_stack+217883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+207388,int_stack+122776,int_stack+122716,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+207568,int_stack+122866,int_stack+122776,6);
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+207838,int_stack+207568,int_stack+207388,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11600,int_stack+122992,int_stack+122866,6);
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11978,int_stack+11600,int_stack+207568,6);
 /*--- compute (d0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+177623,int_stack+11978,int_stack+207838,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+11600,int_stack+124000,int_stack+123900,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11900,int_stack+124150,int_stack+124000,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+148588,int_stack+11900,int_stack+11600,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+207388,int_stack+124360,int_stack+124150,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+190293,int_stack+207388,int_stack+11900,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+207388,int_stack+190293,int_stack+148588,10);
 /*--- compute (dp|ff) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+199693,int_stack+207388,int_stack+177623, 2.0, int_stack+212816, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+177623,int_stack+124790,int_stack+124640,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+190293,int_stack+125015,int_stack+124790,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+220683,int_stack+190293,int_stack+177623,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11600,int_stack+125330,int_stack+125015,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+234978,int_stack+11600,int_stack+190293,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+239478,int_stack+234978,int_stack+220683,15);
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+103800,int_stack+239478,int_stack+207388, 2.0, int_stack+217883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+106800,int_stack+141888,int_stack+133368,100);
     Libderiv->ABCD[11] = int_stack + 106800;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+110400,int_stack+153283,int_stack+146788,100);
     Libderiv->ABCD[10] = int_stack + 110400;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+114000,int_stack+163228,int_stack+156283,100);
     Libderiv->ABCD[9] = int_stack + 114000;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+117600,int_stack+174623,int_stack+168128,100);
     Libderiv->ABCD[8] = int_stack + 117600;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+121200,int_stack+183593,int_stack+172538,100);
     Libderiv->ABCD[7] = int_stack + 121200;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+138828,int_stack+191193,int_stack+188493,100);
     Libderiv->ABCD[6] = int_stack + 138828;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+190293,int_stack+204388,int_stack+197893, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+194193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 190293;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+203803,int_stack+214883,int_stack+208388, 0.0, zero_stack, 1.0, int_stack+194193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 203803;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+211988,int_stack+225378,int_stack+218883, 1.0, int_stack+194193, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 211988;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+174338,int_stack+230178,int_stack+228378,100);
     Libderiv->ABCD[155] = int_stack + 174338;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+225183,int_stack+3800,int_stack+221583,100);
     Libderiv->ABCD[143] = int_stack + 225183;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+3000,int_stack+236478,int_stack+233178,100);
     Libderiv->ABCD[142] = int_stack + 3000;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+228783,int_stack+8600,int_stack+6800,100);
     Libderiv->ABCD[131] = int_stack + 228783;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+6600,int_stack+240978,int_stack+223383,100);
     Libderiv->ABCD[130] = int_stack + 6600;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+220683,int_stack+14400,int_stack+12600,100);
     Libderiv->ABCD[129] = int_stack + 220683;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+10200,int_stack+243978,int_stack+17400,100);
     Libderiv->ABCD[119] = int_stack + 10200;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+13800,int_stack+246978,int_stack+19200,100);
     Libderiv->ABCD[118] = int_stack + 13800;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+17400,int_stack+249978,int_stack+21000,100);
     Libderiv->ABCD[117] = int_stack + 17400;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+232383,int_stack+252978,int_stack+22800,100);
     Libderiv->ABCD[116] = int_stack + 232383;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+21000,int_stack+26400,int_stack+24600,100);
     Libderiv->ABCD[107] = int_stack + 21000;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+24600,int_stack+257778,int_stack+255978,100);
     Libderiv->ABCD[106] = int_stack + 24600;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+235983,int_stack+260778,int_stack+29400,100);
     Libderiv->ABCD[105] = int_stack + 235983;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+239583,int_stack+33000,int_stack+31200,100);
     Libderiv->ABCD[104] = int_stack + 239583;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+28200,int_stack+265578,int_stack+263778,100);
     Libderiv->ABCD[103] = int_stack + 28200;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+31800,int_stack+37800,int_stack+36000,100);
     Libderiv->ABCD[95] = int_stack + 31800;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+35400,int_stack+135828,int_stack+40800,100);
     Libderiv->ABCD[94] = int_stack + 35400;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+39000,int_stack+42600,int_stack+149188,100);
     Libderiv->ABCD[93] = int_stack + 39000;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+135168,int_stack+45600,int_stack+158083,100);
     Libderiv->ABCD[92] = int_stack + 135168;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+42600,int_stack+48600,int_stack+169928,100);
     Libderiv->ABCD[91] = int_stack + 42600;
 /*--- compute (dd|ff) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+46200,int_stack+51600,int_stack+178223,100);
     Libderiv->ABCD[90] = int_stack + 46200;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+49800,int_stack+268578,int_stack+54600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[47] = int_stack + 49800;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+243183,int_stack+58200,int_stack+56400, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[46] = int_stack + 243183;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+53400,int_stack+273378,int_stack+271578, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+156283, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[45] = int_stack + 53400;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+57000,int_stack+63000,int_stack+61200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+168128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[44] = int_stack + 57000;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+60600,int_stack+278178,int_stack+276378, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+172538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[43] = int_stack + 60600;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+246783,int_stack+66000,int_stack+128888, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+188493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[42] = int_stack + 246783;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+64200,int_stack+70800,int_stack+69000, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+197893, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[38] = int_stack + 64200;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+67800,int_stack+282978,int_stack+281178, 0.0, zero_stack, 1.0, int_stack+133368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[35] = int_stack + 67800;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+250383,int_stack+285978,int_stack+73800, 0.0, zero_stack, 1.0, int_stack+146788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[34] = int_stack + 250383;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+71400,int_stack+77400,int_stack+75600, 0.0, zero_stack, 1.0, int_stack+156283, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[33] = int_stack + 71400;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+75000,int_stack+290778,int_stack+288978, 0.0, zero_stack, 1.0, int_stack+168128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[32] = int_stack + 75000;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+253983,int_stack+82200,int_stack+80400, 0.0, zero_stack, 1.0, int_stack+172538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[31] = int_stack + 253983;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+78600,int_stack+85200,int_stack+293778, 0.0, zero_stack, 1.0, int_stack+188493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[30] = int_stack + 78600;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+82200,int_stack+295578,int_stack+202003, 0.0, zero_stack, 1.0, int_stack+197893, 1.0, int_stack+208388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[26] = int_stack + 82200;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+257583,int_stack+90000,int_stack+88200, 0.0, zero_stack, 2.0, int_stack+208388, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[25] = int_stack + 257583;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+85800,int_stack+298578,int_stack+93000, 1.0, int_stack+133368, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[23] = int_stack + 85800;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+89400,int_stack+94800,int_stack+130688, 1.0, int_stack+146788, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[22] = int_stack + 89400;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+93000,int_stack+97800,int_stack+144888, 1.0, int_stack+156283, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[21] = int_stack + 93000;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+96600,int_stack+159883,int_stack+150988, 1.0, int_stack+168128, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[20] = int_stack + 96600;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+128750,int_stack+100800,int_stack+166228, 1.0, int_stack+172538, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[19] = int_stack + 128750;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+100200,int_stack+0,int_stack+210188, 1.0, int_stack+188493, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[18] = int_stack + 100200;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+261183,int_stack+125750,int_stack+186593, 1.0, int_stack+197893, 0.0, zero_stack, 1.0, int_stack+218883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[14] = int_stack + 261183;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+183023,int_stack+180023,int_stack+195993, 1.0, int_stack+208388, 1.0, int_stack+218883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[13] = int_stack + 183023;
 /*--- compute (dd|ff) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+193893,int_stack+103800,int_stack+199693, 2.0, int_stack+218883, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[12] = int_stack + 193893;

}
