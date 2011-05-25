#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_fpff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|ff) integrals */

void d12hrr_order_fpff(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv2_classes[3][3][143] = int_stack + 4375;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 4475;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 4625;
 Libderiv->deriv2_classes[3][6][143] = int_stack + 4835;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 5115;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 5265;
 Libderiv->deriv2_classes[4][5][143] = int_stack + 5490;
 Libderiv->deriv2_classes[4][6][143] = int_stack + 5805;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 6225;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 6325;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 6475;
 Libderiv->deriv2_classes[3][6][131] = int_stack + 6685;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 6965;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 7115;
 Libderiv->deriv2_classes[4][5][131] = int_stack + 7340;
 Libderiv->deriv2_classes[4][6][131] = int_stack + 7655;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 8075;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 8175;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 8325;
 Libderiv->deriv2_classes[3][6][130] = int_stack + 8535;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 8815;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 8965;
 Libderiv->deriv2_classes[4][5][130] = int_stack + 9190;
 Libderiv->deriv2_classes[4][6][130] = int_stack + 9505;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 9925;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 10025;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 10175;
 Libderiv->deriv2_classes[3][6][119] = int_stack + 10385;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 10665;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 10815;
 Libderiv->deriv2_classes[4][5][119] = int_stack + 11040;
 Libderiv->deriv2_classes[4][6][119] = int_stack + 11355;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 11775;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 11875;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 12025;
 Libderiv->deriv2_classes[3][6][118] = int_stack + 12235;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 12515;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 12665;
 Libderiv->deriv2_classes[4][5][118] = int_stack + 12890;
 Libderiv->deriv2_classes[4][6][118] = int_stack + 13205;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 13625;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 13725;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 13875;
 Libderiv->deriv2_classes[3][6][117] = int_stack + 14085;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 14365;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 14515;
 Libderiv->deriv2_classes[4][5][117] = int_stack + 14740;
 Libderiv->deriv2_classes[4][6][117] = int_stack + 15055;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 15475;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 15575;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 15725;
 Libderiv->deriv2_classes[3][6][107] = int_stack + 15935;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 16215;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 16365;
 Libderiv->deriv2_classes[4][5][107] = int_stack + 16590;
 Libderiv->deriv2_classes[4][6][107] = int_stack + 16905;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 17325;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 17425;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 17575;
 Libderiv->deriv2_classes[3][6][106] = int_stack + 17785;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 18065;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 18215;
 Libderiv->deriv2_classes[4][5][106] = int_stack + 18440;
 Libderiv->deriv2_classes[4][6][106] = int_stack + 18755;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 19175;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 19275;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 19425;
 Libderiv->deriv2_classes[3][6][105] = int_stack + 19635;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 19915;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 20065;
 Libderiv->deriv2_classes[4][5][105] = int_stack + 20290;
 Libderiv->deriv2_classes[4][6][105] = int_stack + 20605;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 21025;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 21125;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 21275;
 Libderiv->deriv2_classes[3][6][104] = int_stack + 21485;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 21765;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 21915;
 Libderiv->deriv2_classes[4][5][104] = int_stack + 22140;
 Libderiv->deriv2_classes[4][6][104] = int_stack + 22455;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 22875;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 22975;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 23125;
 Libderiv->deriv2_classes[3][6][95] = int_stack + 23335;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 23615;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 23765;
 Libderiv->deriv2_classes[4][5][95] = int_stack + 23990;
 Libderiv->deriv2_classes[4][6][95] = int_stack + 24305;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 24725;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 24825;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 24975;
 Libderiv->deriv2_classes[3][6][94] = int_stack + 25185;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 25465;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 25615;
 Libderiv->deriv2_classes[4][5][94] = int_stack + 25840;
 Libderiv->deriv2_classes[4][6][94] = int_stack + 26155;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 26575;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 26675;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 26825;
 Libderiv->deriv2_classes[3][6][93] = int_stack + 27035;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 27315;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 27465;
 Libderiv->deriv2_classes[4][5][93] = int_stack + 27690;
 Libderiv->deriv2_classes[4][6][93] = int_stack + 28005;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 28425;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 28525;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 28675;
 Libderiv->deriv2_classes[3][6][92] = int_stack + 28885;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 29165;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 29315;
 Libderiv->deriv2_classes[4][5][92] = int_stack + 29540;
 Libderiv->deriv2_classes[4][6][92] = int_stack + 29855;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 30275;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 30375;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 30525;
 Libderiv->deriv2_classes[3][6][91] = int_stack + 30735;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 31015;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 31165;
 Libderiv->deriv2_classes[4][5][91] = int_stack + 31390;
 Libderiv->deriv2_classes[4][6][91] = int_stack + 31705;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 32125;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 32225;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 32375;
 Libderiv->deriv2_classes[3][6][83] = int_stack + 32585;
 Libderiv->deriv_classes[4][3][11] = int_stack + 32865;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 33015;
 Libderiv->deriv_classes[4][4][11] = int_stack + 33165;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 33390;
 Libderiv->deriv_classes[4][5][11] = int_stack + 33615;
 Libderiv->deriv2_classes[4][5][83] = int_stack + 33930;
 Libderiv->deriv2_classes[4][6][83] = int_stack + 34245;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 34665;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 34765;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 34915;
 Libderiv->deriv2_classes[3][6][82] = int_stack + 35125;
 Libderiv->deriv_classes[4][3][10] = int_stack + 35405;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 35555;
 Libderiv->deriv_classes[4][4][10] = int_stack + 35705;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 35930;
 Libderiv->deriv_classes[4][5][10] = int_stack + 36155;
 Libderiv->deriv2_classes[4][5][82] = int_stack + 36470;
 Libderiv->deriv2_classes[4][6][82] = int_stack + 36785;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 37205;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 37305;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 37455;
 Libderiv->deriv2_classes[3][6][81] = int_stack + 37665;
 Libderiv->deriv_classes[4][3][9] = int_stack + 37945;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 38095;
 Libderiv->deriv_classes[4][4][9] = int_stack + 38245;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 38470;
 Libderiv->deriv_classes[4][5][9] = int_stack + 38695;
 Libderiv->deriv2_classes[4][5][81] = int_stack + 39010;
 Libderiv->deriv2_classes[4][6][81] = int_stack + 39325;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 39745;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 39845;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 39995;
 Libderiv->deriv2_classes[3][6][80] = int_stack + 40205;
 Libderiv->deriv_classes[4][3][8] = int_stack + 40485;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 40635;
 Libderiv->deriv_classes[4][4][8] = int_stack + 40785;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 41010;
 Libderiv->deriv_classes[4][5][8] = int_stack + 41235;
 Libderiv->deriv2_classes[4][5][80] = int_stack + 41550;
 Libderiv->deriv2_classes[4][6][80] = int_stack + 41865;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 42285;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 42385;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 42535;
 Libderiv->deriv2_classes[3][6][79] = int_stack + 42745;
 Libderiv->deriv_classes[4][3][7] = int_stack + 43025;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 43175;
 Libderiv->deriv_classes[4][4][7] = int_stack + 43325;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 43550;
 Libderiv->deriv_classes[4][5][7] = int_stack + 43775;
 Libderiv->deriv2_classes[4][5][79] = int_stack + 44090;
 Libderiv->deriv2_classes[4][6][79] = int_stack + 44405;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 44825;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 44925;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 45075;
 Libderiv->deriv2_classes[3][6][78] = int_stack + 45285;
 Libderiv->dvrr_classes[4][3] = int_stack + 45565;
 Libderiv->deriv_classes[4][3][6] = int_stack + 45715;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 45865;
 Libderiv->dvrr_classes[4][4] = int_stack + 46015;
 Libderiv->deriv_classes[4][4][6] = int_stack + 46240;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 46465;
 Libderiv->deriv_classes[4][5][6] = int_stack + 46690;
 Libderiv->deriv2_classes[4][5][78] = int_stack + 47005;
 Libderiv->deriv2_classes[4][6][78] = int_stack + 47320;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 47740;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 47840;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 47990;
 Libderiv->deriv2_classes[3][6][35] = int_stack + 48200;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 48480;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 48630;
 Libderiv->deriv2_classes[4][5][35] = int_stack + 48855;
 Libderiv->deriv2_classes[4][6][35] = int_stack + 49170;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 49590;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 49690;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 49840;
 Libderiv->deriv2_classes[3][6][34] = int_stack + 50050;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 50330;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 50480;
 Libderiv->deriv2_classes[4][5][34] = int_stack + 50705;
 Libderiv->deriv2_classes[4][6][34] = int_stack + 51020;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 51440;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 51540;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 51690;
 Libderiv->deriv2_classes[3][6][33] = int_stack + 51900;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 52180;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 52330;
 Libderiv->deriv2_classes[4][5][33] = int_stack + 52555;
 Libderiv->deriv2_classes[4][6][33] = int_stack + 52870;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 53290;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 53390;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 53540;
 Libderiv->deriv2_classes[3][6][32] = int_stack + 53750;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 54030;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 54180;
 Libderiv->deriv2_classes[4][5][32] = int_stack + 54405;
 Libderiv->deriv2_classes[4][6][32] = int_stack + 54720;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 55140;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 55240;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 55390;
 Libderiv->deriv2_classes[3][6][31] = int_stack + 55600;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 55880;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 56030;
 Libderiv->deriv2_classes[4][5][31] = int_stack + 56255;
 Libderiv->deriv2_classes[4][6][31] = int_stack + 56570;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 56990;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 57090;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 57240;
 Libderiv->deriv2_classes[3][6][30] = int_stack + 57450;
 Libderiv->deriv_classes[4][3][2] = int_stack + 57730;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 57880;
 Libderiv->deriv_classes[4][4][2] = int_stack + 58030;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 58255;
 Libderiv->deriv_classes[4][5][2] = int_stack + 58480;
 Libderiv->deriv2_classes[4][5][30] = int_stack + 58795;
 Libderiv->deriv2_classes[4][6][30] = int_stack + 59110;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 59530;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 59630;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 59780;
 Libderiv->deriv2_classes[3][6][26] = int_stack + 59990;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 60270;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 60420;
 Libderiv->deriv2_classes[4][5][26] = int_stack + 60645;
 Libderiv->deriv2_classes[4][6][26] = int_stack + 60960;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 61380;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 61480;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 61630;
 Libderiv->deriv2_classes[3][6][23] = int_stack + 61840;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 62120;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 62270;
 Libderiv->deriv2_classes[4][5][23] = int_stack + 62495;
 Libderiv->deriv2_classes[4][6][23] = int_stack + 62810;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 63230;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 63330;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 63480;
 Libderiv->deriv2_classes[3][6][22] = int_stack + 63690;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 63970;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 64120;
 Libderiv->deriv2_classes[4][5][22] = int_stack + 64345;
 Libderiv->deriv2_classes[4][6][22] = int_stack + 64660;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 65080;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 65180;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 65330;
 Libderiv->deriv2_classes[3][6][21] = int_stack + 65540;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 65820;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 65970;
 Libderiv->deriv2_classes[4][5][21] = int_stack + 66195;
 Libderiv->deriv2_classes[4][6][21] = int_stack + 66510;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 66930;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 67030;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 67180;
 Libderiv->deriv2_classes[3][6][20] = int_stack + 67390;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 67670;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 67820;
 Libderiv->deriv2_classes[4][5][20] = int_stack + 68045;
 Libderiv->deriv2_classes[4][6][20] = int_stack + 68360;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 68780;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 68880;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 69030;
 Libderiv->deriv2_classes[3][6][19] = int_stack + 69240;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 69520;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 69670;
 Libderiv->deriv2_classes[4][5][19] = int_stack + 69895;
 Libderiv->deriv2_classes[4][6][19] = int_stack + 70210;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 70630;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 70730;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 70880;
 Libderiv->deriv2_classes[3][6][18] = int_stack + 71090;
 Libderiv->deriv_classes[4][3][1] = int_stack + 71370;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 71520;
 Libderiv->deriv_classes[4][4][1] = int_stack + 71670;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 71895;
 Libderiv->deriv_classes[4][5][1] = int_stack + 72120;
 Libderiv->deriv2_classes[4][5][18] = int_stack + 72435;
 Libderiv->deriv2_classes[4][6][18] = int_stack + 72750;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 73170;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 73270;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 73420;
 Libderiv->deriv2_classes[3][6][14] = int_stack + 73630;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 73910;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 74060;
 Libderiv->deriv2_classes[4][5][14] = int_stack + 74285;
 Libderiv->deriv2_classes[4][6][14] = int_stack + 74600;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 75020;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 75120;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 75270;
 Libderiv->deriv2_classes[3][6][13] = int_stack + 75480;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 75760;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 75910;
 Libderiv->deriv2_classes[4][5][13] = int_stack + 76135;
 Libderiv->deriv2_classes[4][6][13] = int_stack + 76450;
 Libderiv->deriv_classes[3][3][11] = int_stack + 76870;
 Libderiv->deriv_classes[3][4][11] = int_stack + 76970;
 Libderiv->deriv_classes[3][5][11] = int_stack + 77120;
 Libderiv->deriv_classes[3][6][11] = int_stack + 77330;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 77610;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 77710;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 77860;
 Libderiv->deriv2_classes[3][6][11] = int_stack + 78070;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 78350;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 78500;
 Libderiv->deriv2_classes[4][5][11] = int_stack + 78725;
 Libderiv->deriv2_classes[4][6][11] = int_stack + 79040;
 Libderiv->deriv_classes[3][3][10] = int_stack + 79460;
 Libderiv->deriv_classes[3][4][10] = int_stack + 79560;
 Libderiv->deriv_classes[3][5][10] = int_stack + 79710;
 Libderiv->deriv_classes[3][6][10] = int_stack + 79920;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 80200;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 80300;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 80450;
 Libderiv->deriv2_classes[3][6][10] = int_stack + 80660;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 80940;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 81090;
 Libderiv->deriv2_classes[4][5][10] = int_stack + 81315;
 Libderiv->deriv2_classes[4][6][10] = int_stack + 81630;
 Libderiv->deriv_classes[3][3][9] = int_stack + 82050;
 Libderiv->deriv_classes[3][4][9] = int_stack + 82150;
 Libderiv->deriv_classes[3][5][9] = int_stack + 82300;
 Libderiv->deriv_classes[3][6][9] = int_stack + 82510;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 82790;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 82890;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 83040;
 Libderiv->deriv2_classes[3][6][9] = int_stack + 83250;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 83530;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 83680;
 Libderiv->deriv2_classes[4][5][9] = int_stack + 83905;
 Libderiv->deriv2_classes[4][6][9] = int_stack + 84220;
 Libderiv->deriv_classes[3][3][8] = int_stack + 84640;
 Libderiv->deriv_classes[3][4][8] = int_stack + 84740;
 Libderiv->deriv_classes[3][5][8] = int_stack + 84890;
 Libderiv->deriv_classes[3][6][8] = int_stack + 85100;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 85380;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 85480;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 85630;
 Libderiv->deriv2_classes[3][6][8] = int_stack + 85840;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 86120;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 86270;
 Libderiv->deriv2_classes[4][5][8] = int_stack + 86495;
 Libderiv->deriv2_classes[4][6][8] = int_stack + 86810;
 Libderiv->deriv_classes[3][3][7] = int_stack + 87230;
 Libderiv->deriv_classes[3][4][7] = int_stack + 87330;
 Libderiv->deriv_classes[3][5][7] = int_stack + 87480;
 Libderiv->deriv_classes[3][6][7] = int_stack + 87690;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 87970;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 88070;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 88220;
 Libderiv->deriv2_classes[3][6][7] = int_stack + 88430;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 88710;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 88860;
 Libderiv->deriv2_classes[4][5][7] = int_stack + 89085;
 Libderiv->deriv2_classes[4][6][7] = int_stack + 89400;
 Libderiv->dvrr_classes[3][3] = int_stack + 89820;
 Libderiv->deriv_classes[3][3][6] = int_stack + 89920;
 Libderiv->dvrr_classes[3][4] = int_stack + 90020;
 Libderiv->deriv_classes[3][4][6] = int_stack + 90170;
 Libderiv->dvrr_classes[3][5] = int_stack + 90320;
 Libderiv->deriv_classes[3][5][6] = int_stack + 90530;
 Libderiv->deriv_classes[3][6][6] = int_stack + 90740;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 91020;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 91120;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 91270;
 Libderiv->deriv2_classes[3][6][6] = int_stack + 91480;
 Libderiv->deriv_classes[4][3][0] = int_stack + 91760;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 91910;
 Libderiv->deriv_classes[4][4][0] = int_stack + 92060;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 92285;
 Libderiv->deriv_classes[4][5][0] = int_stack + 92510;
 Libderiv->deriv2_classes[4][5][6] = int_stack + 92825;
 Libderiv->deriv2_classes[4][6][6] = int_stack + 93140;
 Libderiv->deriv_classes[3][3][2] = int_stack + 93560;
 Libderiv->deriv_classes[3][4][2] = int_stack + 93660;
 Libderiv->deriv_classes[3][5][2] = int_stack + 93810;
 Libderiv->deriv_classes[3][6][2] = int_stack + 94020;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 94300;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 94400;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 94550;
 Libderiv->deriv2_classes[3][6][2] = int_stack + 94760;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 95040;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 95190;
 Libderiv->deriv2_classes[4][5][2] = int_stack + 95415;
 Libderiv->deriv2_classes[4][6][2] = int_stack + 95730;
 Libderiv->deriv_classes[3][3][1] = int_stack + 96150;
 Libderiv->deriv_classes[3][4][1] = int_stack + 96250;
 Libderiv->deriv_classes[3][5][1] = int_stack + 96400;
 Libderiv->deriv_classes[3][6][1] = int_stack + 96610;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 96890;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 96990;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 97140;
 Libderiv->deriv2_classes[3][6][1] = int_stack + 97350;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 97630;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 97780;
 Libderiv->deriv2_classes[4][5][1] = int_stack + 98005;
 Libderiv->deriv2_classes[4][6][1] = int_stack + 98320;
 Libderiv->deriv_classes[3][3][0] = int_stack + 98740;
 Libderiv->deriv_classes[3][4][0] = int_stack + 98840;
 Libderiv->deriv_classes[3][5][0] = int_stack + 98990;
 Libderiv->deriv_classes[3][6][0] = int_stack + 99200;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 99480;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 99580;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 99730;
 Libderiv->deriv2_classes[3][6][0] = int_stack + 99940;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 100220;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 100370;
 Libderiv->deriv2_classes[4][5][0] = int_stack + 100595;
 Libderiv->deriv2_classes[4][6][0] = int_stack + 100910;
 memset(int_stack,0,810640);

 Libderiv->dvrr_stack = int_stack + 206855;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_fpff(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+101330,int_stack+90020,int_stack+89820,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+101630,int_stack+90320,int_stack+90020,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+102080,int_stack+101630,int_stack+101330,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+102680,int_stack+76970,int_stack+76870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89820,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+102980,int_stack+77120,int_stack+76970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90020,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+103430,int_stack+102980,int_stack+102680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101330,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+104030,int_stack+77330,int_stack+77120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90320,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+104660,int_stack+104030,int_stack+102980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101630,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+105560,int_stack+104660,int_stack+103430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102080,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+104030,int_stack+46015,int_stack+45565,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+104480,int_stack+2100,int_stack+46015,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+104480,int_stack+104030,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+107460,int_stack+33165,int_stack+32865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45565,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+107910,int_stack+33615,int_stack+33165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46015,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+108585,int_stack+107910,int_stack+107460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104030,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+109485,int_stack+0,int_stack+33615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+110430,int_stack+109485,int_stack+107910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104480,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+111780,int_stack+110430,int_stack+108585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106560,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+109485,int_stack+79560,int_stack+79460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89820, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+109785,int_stack+79710,int_stack+79560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90020, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+110235,int_stack+109785,int_stack+109485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101330, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+110835,int_stack+79920,int_stack+79710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90320, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+113280,int_stack+110835,int_stack+109785, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101630, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+114180,int_stack+113280,int_stack+110235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102080, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+113280,int_stack+35705,int_stack+35405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45565, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+110835,int_stack+36155,int_stack+35705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46015, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+115180,int_stack+110835,int_stack+113280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104030, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+116080,int_stack+420,int_stack+36155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+117025,int_stack+116080,int_stack+110835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104480, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+118375,int_stack+117025,int_stack+115180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106560, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116080,int_stack+82150,int_stack+82050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89820, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113730,int_stack+82300,int_stack+82150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90020, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+116380,int_stack+113730,int_stack+116080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101330, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+116980,int_stack+82510,int_stack+82300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90320, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+119875,int_stack+116980,int_stack+113730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101630, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+116980,int_stack+119875,int_stack+116380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102080, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+119875,int_stack+38245,int_stack+37945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45565, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+120325,int_stack+38695,int_stack+38245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46015, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+121000,int_stack+120325,int_stack+119875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104030, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+121900,int_stack+840,int_stack+38695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+122845,int_stack+121900,int_stack+120325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104480, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+124195,int_stack+122845,int_stack+121000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106560, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+121900,int_stack+84740,int_stack+84640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+122200,int_stack+84890,int_stack+84740, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+122650,int_stack+122200,int_stack+121900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+123250,int_stack+85100,int_stack+84890, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+123250,int_stack+122200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+125695,int_stack+0,int_stack+122650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+40785,int_stack+40485, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+41235,int_stack+40785, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+123250,int_stack+450,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126695,int_stack+1260,int_stack+41235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+127640,int_stack+126695,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+104480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+128990,int_stack+127640,int_stack+123250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+126695,int_stack+87330,int_stack+87230, 0.0, zero_stack, 1.0, int_stack+89820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+126995,int_stack+87480,int_stack+87330, 0.0, zero_stack, 1.0, int_stack+90020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+127445,int_stack+126995,int_stack+126695, 0.0, zero_stack, 1.0, int_stack+101330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+128045,int_stack+87690,int_stack+87480, 0.0, zero_stack, 1.0, int_stack+90320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+130490,int_stack+128045,int_stack+126995, 0.0, zero_stack, 1.0, int_stack+101630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+131390,int_stack+130490,int_stack+127445, 0.0, zero_stack, 1.0, int_stack+102080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+130490,int_stack+43325,int_stack+43025, 0.0, zero_stack, 1.0, int_stack+45565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+128045,int_stack+43775,int_stack+43325, 0.0, zero_stack, 1.0, int_stack+46015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+132390,int_stack+128045,int_stack+130490, 0.0, zero_stack, 1.0, int_stack+104030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+133290,int_stack+1680,int_stack+43775, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+134235,int_stack+133290,int_stack+128045, 0.0, zero_stack, 1.0, int_stack+104480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+135585,int_stack+134235,int_stack+132390, 0.0, zero_stack, 1.0, int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+133290,int_stack+90170,int_stack+89920, 1.0, int_stack+89820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130940,int_stack+90530,int_stack+90170, 1.0, int_stack+90020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+133590,int_stack+130940,int_stack+133290, 1.0, int_stack+101330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+134190,int_stack+90740,int_stack+90530, 1.0, int_stack+90320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1125,int_stack+134190,int_stack+130940, 1.0, int_stack+101630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+134190,int_stack+1125,int_stack+133590, 1.0, int_stack+102080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1125,int_stack+46240,int_stack+45715, 1.0, int_stack+45565, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+137085,int_stack+46690,int_stack+46240, 1.0, int_stack+46015, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+137760,int_stack+137085,int_stack+1125, 1.0, int_stack+104030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+138660,int_stack+2415,int_stack+46690, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+139605,int_stack+138660,int_stack+137085, 1.0, int_stack+104480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+104030,int_stack+139605,int_stack+137760, 1.0, int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+106560,int_stack+3675,int_stack+90320,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+138660,int_stack+106560,int_stack+101630,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+139560,int_stack+138660,int_stack+102080,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+138660,int_stack+93660,int_stack+93560,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+138960,int_stack+93810,int_stack+93660,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+138960,int_stack+138660,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1575,int_stack+94020,int_stack+93810,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+101330,int_stack+1575,int_stack+138960,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+1575,int_stack+101330,int_stack+106560,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+101330,int_stack+58030,int_stack+57730,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+101780,int_stack+58480,int_stack+58030,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+140560,int_stack+101780,int_stack+101330,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+141460,int_stack+2835,int_stack+58480,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+142405,int_stack+141460,int_stack+101780,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+143755,int_stack+142405,int_stack+140560,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+107160,int_stack+96250,int_stack+96150,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+141460,int_stack+96400,int_stack+96250,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+141910,int_stack+141460,int_stack+107160,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+142510,int_stack+96610,int_stack+96400,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+145255,int_stack+142510,int_stack+141460,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+142510,int_stack+145255,int_stack+141910,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+145255,int_stack+71670,int_stack+71370,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+145705,int_stack+72120,int_stack+71670,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+146380,int_stack+145705,int_stack+145255,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+147280,int_stack+3255,int_stack+72120,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2575,int_stack+147280,int_stack+145705,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+147280,int_stack+2575,int_stack+146380,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2575,int_stack+98840,int_stack+98740,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2875,int_stack+98990,int_stack+98840,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3325,int_stack+2875,int_stack+2575,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+148780,int_stack+99200,int_stack+98990,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+149410,int_stack+148780,int_stack+2875,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+150310,int_stack+149410,int_stack+3325,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+148780,int_stack+92060,int_stack+91760,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+149230,int_stack+92510,int_stack+92060,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+151310,int_stack+149230,int_stack+148780,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+152210,int_stack+3955,int_stack+92510,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+153155,int_stack+152210,int_stack+149230,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+154505,int_stack+153155,int_stack+151310,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+152210,int_stack+4475,int_stack+4375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+76870,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152510,int_stack+4625,int_stack+4475, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+76970,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152960,int_stack+152510,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+102680,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+153560,int_stack+4835,int_stack+4625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+77120,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3925,int_stack+153560,int_stack+152510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+102980,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+156005,int_stack+3925,int_stack+152960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+103430,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3925,int_stack+5265,int_stack+5115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+32865,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4375,int_stack+5490,int_stack+5265, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+33165,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152210,int_stack+4375,int_stack+3925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+107460,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+153110,int_stack+5805,int_stack+5490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+33615,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+153110,int_stack+4375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+107910,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3925,int_stack+157005,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+108585,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+152210,int_stack+6325,int_stack+6225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76870, 1.0, int_stack+79460,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152510,int_stack+6475,int_stack+6325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76970, 1.0, int_stack+79560,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152960,int_stack+152510,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102680, 1.0, int_stack+109485,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+153560,int_stack+6685,int_stack+6475, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77120, 1.0, int_stack+79710,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+153560,int_stack+152510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102980, 1.0, int_stack+109785,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+157905,int_stack+157005,int_stack+152960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103430, 1.0, int_stack+110235,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+7115,int_stack+6965, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32865, 1.0, int_stack+35405,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152210,int_stack+7340,int_stack+7115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33165, 1.0, int_stack+35705,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152885,int_stack+152210,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107460, 1.0, int_stack+113280,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5425,int_stack+7655,int_stack+7340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33615, 1.0, int_stack+36155,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6370,int_stack+5425,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107910, 1.0, int_stack+110835,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+158905,int_stack+6370,int_stack+152885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108585, 1.0, int_stack+115180,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+152210,int_stack+8175,int_stack+8075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+79460, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152510,int_stack+8325,int_stack+8175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+79560, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152960,int_stack+152510,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+109485, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+153560,int_stack+8535,int_stack+8325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+79710, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+153560,int_stack+152510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+109785, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5425,int_stack+157005,int_stack+152960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+110235, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+8965,int_stack+8815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+35405, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152210,int_stack+9190,int_stack+8965, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+35705, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152885,int_stack+152210,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+113280, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6425,int_stack+9505,int_stack+9190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+36155, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7370,int_stack+6425,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+110835, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+160405,int_stack+7370,int_stack+152885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+115180, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+152210,int_stack+10025,int_stack+9925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76870, 0.0, zero_stack, 1.0, int_stack+82050,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152510,int_stack+10175,int_stack+10025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76970, 0.0, zero_stack, 1.0, int_stack+82150,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152960,int_stack+152510,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102680, 0.0, zero_stack, 1.0, int_stack+116080,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+153560,int_stack+10385,int_stack+10175, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77120, 0.0, zero_stack, 1.0, int_stack+82300,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+153560,int_stack+152510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102980, 0.0, zero_stack, 1.0, int_stack+113730,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6425,int_stack+157005,int_stack+152960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103430, 0.0, zero_stack, 1.0, int_stack+116380,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+10815,int_stack+10665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32865, 0.0, zero_stack, 1.0, int_stack+37945,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152210,int_stack+11040,int_stack+10815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33165, 0.0, zero_stack, 1.0, int_stack+38245,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152885,int_stack+152210,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107460, 0.0, zero_stack, 1.0, int_stack+119875,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7425,int_stack+11355,int_stack+11040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33615, 0.0, zero_stack, 1.0, int_stack+38695,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8370,int_stack+7425,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107910, 0.0, zero_stack, 1.0, int_stack+120325,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9720,int_stack+8370,int_stack+152885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108585, 0.0, zero_stack, 1.0, int_stack+121000,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+152210,int_stack+11875,int_stack+11775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79460, 1.0, int_stack+82050, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152510,int_stack+12025,int_stack+11875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79560, 1.0, int_stack+82150, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152960,int_stack+152510,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109485, 1.0, int_stack+116080, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+153560,int_stack+12235,int_stack+12025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79710, 1.0, int_stack+82300, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+153560,int_stack+152510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109785, 1.0, int_stack+113730, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7425,int_stack+157005,int_stack+152960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110235, 1.0, int_stack+116380, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+12665,int_stack+12515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35405, 1.0, int_stack+37945, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152210,int_stack+12890,int_stack+12665, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35705, 1.0, int_stack+38245, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152885,int_stack+152210,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113280, 1.0, int_stack+119875, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8425,int_stack+13205,int_stack+12890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36155, 1.0, int_stack+38695, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11220,int_stack+8425,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110835, 1.0, int_stack+120325, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+161905,int_stack+11220,int_stack+152885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115180, 1.0, int_stack+121000, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11220,int_stack+13725,int_stack+13625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+82050, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11520,int_stack+13875,int_stack+13725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+82150, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11970,int_stack+11520,int_stack+11220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+116080, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12570,int_stack+14085,int_stack+13875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+82300, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+12570,int_stack+11520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+113730, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12570,int_stack+157005,int_stack+11970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+116380, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+14515,int_stack+14365, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+37945, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13570,int_stack+14740,int_stack+14515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+38245, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11220,int_stack+13570,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+119875, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+152210,int_stack+15055,int_stack+14740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+38695, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+153155,int_stack+152210,int_stack+13570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+120325, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13570,int_stack+153155,int_stack+11220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+121000, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11220,int_stack+15575,int_stack+15475, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76870, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84640,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11520,int_stack+15725,int_stack+15575, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76970, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84740,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11970,int_stack+11520,int_stack+11220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121900,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15070,int_stack+15935,int_stack+15725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+77120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+84890,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+15070,int_stack+11520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+102980, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122200,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+15070,int_stack+157005,int_stack+11970, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+103430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+122650,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+16365,int_stack+16215, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32865, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40485,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11220,int_stack+16590,int_stack+16365, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33165, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40785,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152210,int_stack+11220,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+153110,int_stack+16905,int_stack+16590, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33615, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41235,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+163405,int_stack+153110,int_stack+11220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107910, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+164755,int_stack+163405,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+108585, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+123250,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+152210,int_stack+17425,int_stack+17325, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79460, 0.0, zero_stack, 1.0, int_stack+84640, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152510,int_stack+17575,int_stack+17425, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79560, 0.0, zero_stack, 1.0, int_stack+84740, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152960,int_stack+152510,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109485, 0.0, zero_stack, 1.0, int_stack+121900, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+153560,int_stack+17785,int_stack+17575, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+79710, 0.0, zero_stack, 1.0, int_stack+84890, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+153560,int_stack+152510, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+109785, 0.0, zero_stack, 1.0, int_stack+122200, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+163405,int_stack+157005,int_stack+152960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110235, 0.0, zero_stack, 1.0, int_stack+122650, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+18215,int_stack+18065, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35405, 0.0, zero_stack, 1.0, int_stack+40485, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+152210,int_stack+18440,int_stack+18215, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35705, 0.0, zero_stack, 1.0, int_stack+40785, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+152885,int_stack+152210,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113280, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11220,int_stack+18755,int_stack+18440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36155, 0.0, zero_stack, 1.0, int_stack+41235, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16070,int_stack+11220,int_stack+152210, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+110835, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+17420,int_stack+16070,int_stack+152885, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+115180, 0.0, zero_stack, 1.0, int_stack+123250, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16070,int_stack+19275,int_stack+19175, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82050, 1.0, int_stack+84640, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16370,int_stack+19425,int_stack+19275, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82150, 1.0, int_stack+84740, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16820,int_stack+16370,int_stack+16070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116080, 1.0, int_stack+121900, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+152210,int_stack+19635,int_stack+19425, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+82300, 1.0, int_stack+84890, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+152210,int_stack+16370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+113730, 1.0, int_stack+122200, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+152210,int_stack+157005,int_stack+16820, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116380, 1.0, int_stack+122650, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+20065,int_stack+19915, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37945, 1.0, int_stack+40485, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+153210,int_stack+20290,int_stack+20065, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38245, 1.0, int_stack+40785, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16070,int_stack+153210,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+119875, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11220,int_stack+20605,int_stack+20290, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38695, 1.0, int_stack+41235, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18920,int_stack+11220,int_stack+153210, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+120325, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+166255,int_stack+18920,int_stack+16070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+121000, 1.0, int_stack+123250, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16070,int_stack+21125,int_stack+21025, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+84640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16370,int_stack+21275,int_stack+21125, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+84740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16820,int_stack+16370,int_stack+16070, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+121900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18920,int_stack+21485,int_stack+21275, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+84890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+18920,int_stack+16370, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+122200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18920,int_stack+157005,int_stack+16820, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+122650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+21915,int_stack+21765, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+19920,int_stack+22140,int_stack+21915, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+40785, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20595,int_stack+19920,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16070,int_stack+22455,int_stack+22140, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+41235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11220,int_stack+16070,int_stack+19920, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+167755,int_stack+11220,int_stack+20595, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+123250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11220,int_stack+22975,int_stack+22875, 0.0, zero_stack, 1.0, int_stack+76870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87230,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11520,int_stack+23125,int_stack+22975, 0.0, zero_stack, 1.0, int_stack+76970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87330,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11970,int_stack+11520,int_stack+11220, 0.0, zero_stack, 1.0, int_stack+102680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126695,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+19920,int_stack+23335,int_stack+23125, 0.0, zero_stack, 1.0, int_stack+77120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87480,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+19920,int_stack+11520, 0.0, zero_stack, 1.0, int_stack+102980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126995,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+19920,int_stack+157005,int_stack+11970, 0.0, zero_stack, 1.0, int_stack+103430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127445,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+23765,int_stack+23615, 0.0, zero_stack, 1.0, int_stack+32865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43025,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20920,int_stack+23990,int_stack+23765, 0.0, zero_stack, 1.0, int_stack+33165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43325,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21595,int_stack+20920,int_stack+157005, 0.0, zero_stack, 1.0, int_stack+107460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130490,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22495,int_stack+24305,int_stack+23990, 0.0, zero_stack, 1.0, int_stack+33615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43775,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11220,int_stack+22495,int_stack+20920, 0.0, zero_stack, 1.0, int_stack+107910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128045,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+22495,int_stack+11220,int_stack+21595, 0.0, zero_stack, 1.0, int_stack+108585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+132390,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11220,int_stack+24825,int_stack+24725, 0.0, zero_stack, 1.0, int_stack+79460, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87230, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11520,int_stack+24975,int_stack+24825, 0.0, zero_stack, 1.0, int_stack+79560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87330, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11970,int_stack+11520,int_stack+11220, 0.0, zero_stack, 1.0, int_stack+109485, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126695, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20920,int_stack+25185,int_stack+24975, 0.0, zero_stack, 1.0, int_stack+79710, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+87480, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+20920,int_stack+11520, 0.0, zero_stack, 1.0, int_stack+109785, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+126995, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+20920,int_stack+157005,int_stack+11970, 0.0, zero_stack, 1.0, int_stack+110235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+127445, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+25615,int_stack+25465, 0.0, zero_stack, 1.0, int_stack+35405, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43025, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11220,int_stack+25840,int_stack+25615, 0.0, zero_stack, 1.0, int_stack+35705, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43325, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23995,int_stack+11220,int_stack+157005, 0.0, zero_stack, 1.0, int_stack+113280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130490, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+24895,int_stack+26155,int_stack+25840, 0.0, zero_stack, 1.0, int_stack+36155, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43775, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16070,int_stack+24895,int_stack+11220, 0.0, zero_stack, 1.0, int_stack+110835, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+128045, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+24895,int_stack+16070,int_stack+23995, 0.0, zero_stack, 1.0, int_stack+115180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+132390, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23995,int_stack+26675,int_stack+26575, 0.0, zero_stack, 1.0, int_stack+82050, 0.0, zero_stack, 1.0, int_stack+87230, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24295,int_stack+26825,int_stack+26675, 0.0, zero_stack, 1.0, int_stack+82150, 0.0, zero_stack, 1.0, int_stack+87330, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16070,int_stack+24295,int_stack+23995, 0.0, zero_stack, 1.0, int_stack+116080, 0.0, zero_stack, 1.0, int_stack+126695, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16670,int_stack+27035,int_stack+26825, 0.0, zero_stack, 1.0, int_stack+82300, 0.0, zero_stack, 1.0, int_stack+87480, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+16670,int_stack+24295, 0.0, zero_stack, 1.0, int_stack+113730, 0.0, zero_stack, 1.0, int_stack+126995, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11220,int_stack+157005,int_stack+16070, 0.0, zero_stack, 1.0, int_stack+116380, 0.0, zero_stack, 1.0, int_stack+127445, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16070,int_stack+27465,int_stack+27315, 0.0, zero_stack, 1.0, int_stack+37945, 0.0, zero_stack, 1.0, int_stack+43025, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+16520,int_stack+27690,int_stack+27465, 0.0, zero_stack, 1.0, int_stack+38245, 0.0, zero_stack, 1.0, int_stack+43325, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+157005,int_stack+16520,int_stack+16070, 0.0, zero_stack, 1.0, int_stack+119875, 0.0, zero_stack, 1.0, int_stack+130490, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+26395,int_stack+28005,int_stack+27690, 0.0, zero_stack, 1.0, int_stack+38695, 0.0, zero_stack, 1.0, int_stack+43775, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+169255,int_stack+26395,int_stack+16520, 0.0, zero_stack, 1.0, int_stack+120325, 0.0, zero_stack, 1.0, int_stack+128045, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+26395,int_stack+169255,int_stack+157005, 0.0, zero_stack, 1.0, int_stack+121000, 0.0, zero_stack, 1.0, int_stack+132390, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+28525,int_stack+28425, 0.0, zero_stack, 1.0, int_stack+84640, 1.0, int_stack+87230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157305,int_stack+28675,int_stack+28525, 0.0, zero_stack, 1.0, int_stack+84740, 1.0, int_stack+87330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+157305,int_stack+157005, 0.0, zero_stack, 1.0, int_stack+121900, 1.0, int_stack+126695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+169855,int_stack+28885,int_stack+28675, 0.0, zero_stack, 1.0, int_stack+84890, 1.0, int_stack+87480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23995,int_stack+169855,int_stack+157305, 0.0, zero_stack, 1.0, int_stack+122200, 1.0, int_stack+126995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+169855,int_stack+23995,int_stack+169255, 0.0, zero_stack, 1.0, int_stack+122650, 1.0, int_stack+127445, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+29315,int_stack+29165, 0.0, zero_stack, 1.0, int_stack+40485, 1.0, int_stack+43025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23995,int_stack+29540,int_stack+29315, 0.0, zero_stack, 1.0, int_stack+40785, 1.0, int_stack+43325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+157005,int_stack+23995,int_stack+169255, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+27895,int_stack+29855,int_stack+29540, 0.0, zero_stack, 1.0, int_stack+41235, 1.0, int_stack+43775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16070,int_stack+27895,int_stack+23995, 0.0, zero_stack, 1.0, int_stack+450, 1.0, int_stack+128045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+27895,int_stack+16070,int_stack+157005, 0.0, zero_stack, 1.0, int_stack+123250, 1.0, int_stack+132390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+30375,int_stack+30275, 0.0, zero_stack, 2.0, int_stack+87230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157305,int_stack+30525,int_stack+30375, 0.0, zero_stack, 2.0, int_stack+87330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+157305,int_stack+157005, 0.0, zero_stack, 2.0, int_stack+126695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16070,int_stack+30735,int_stack+30525, 0.0, zero_stack, 2.0, int_stack+87480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23995,int_stack+16070,int_stack+157305, 0.0, zero_stack, 2.0, int_stack+126995, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+16070,int_stack+23995,int_stack+169255, 0.0, zero_stack, 2.0, int_stack+127445, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+31165,int_stack+31015, 0.0, zero_stack, 2.0, int_stack+43025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+23995,int_stack+31390,int_stack+31165, 0.0, zero_stack, 2.0, int_stack+43325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+157005,int_stack+23995,int_stack+169255, 0.0, zero_stack, 2.0, int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29395,int_stack+31705,int_stack+31390, 0.0, zero_stack, 2.0, int_stack+43775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+30340,int_stack+29395,int_stack+23995, 0.0, zero_stack, 2.0, int_stack+128045, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+170855,int_stack+30340,int_stack+157005, 0.0, zero_stack, 2.0, int_stack+132390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+32225,int_stack+32125, 1.0, int_stack+76870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89920,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157305,int_stack+32375,int_stack+32225, 1.0, int_stack+76970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90170,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+157305,int_stack+157005, 1.0, int_stack+102680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133290,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23995,int_stack+32585,int_stack+32375, 1.0, int_stack+77120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90530,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29395,int_stack+23995,int_stack+157305, 1.0, int_stack+102980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130940,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+30295,int_stack+29395,int_stack+169255, 1.0, int_stack+103430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133590,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+33390,int_stack+33015, 1.0, int_stack+32865, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45715,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29395,int_stack+33930,int_stack+33390, 1.0, int_stack+33165, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46240,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23995,int_stack+29395,int_stack+169255, 1.0, int_stack+107460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+31295,int_stack+34245,int_stack+33930, 1.0, int_stack+33615, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46690,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+32240,int_stack+31295,int_stack+29395, 1.0, int_stack+107910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137085,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+102455,int_stack+32240,int_stack+23995, 1.0, int_stack+108585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137760,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23995,int_stack+34765,int_stack+34665, 1.0, int_stack+79460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89920, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+24295,int_stack+34915,int_stack+34765, 1.0, int_stack+79560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90170, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+24295,int_stack+23995, 1.0, int_stack+109485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133290, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29395,int_stack+35125,int_stack+34915, 1.0, int_stack+79710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90530, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+29395,int_stack+24295, 1.0, int_stack+109785, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130940, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+31295,int_stack+157005,int_stack+169255, 1.0, int_stack+110235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133590, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+35930,int_stack+35555, 1.0, int_stack+35405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45715, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+36470,int_stack+35930, 1.0, int_stack+35705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46240, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29395,int_stack+157005,int_stack+169255, 1.0, int_stack+113280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+32295,int_stack+36785,int_stack+36470, 1.0, int_stack+36155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46690, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+33240,int_stack+32295,int_stack+157005, 1.0, int_stack+110835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137085, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+34590,int_stack+33240,int_stack+29395, 1.0, int_stack+115180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137760, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+115180,int_stack+37305,int_stack+37205, 1.0, int_stack+82050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+89920, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113280,int_stack+37455,int_stack+37305, 1.0, int_stack+82150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90170, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+115480,int_stack+113280,int_stack+115180, 1.0, int_stack+116080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133290, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29395,int_stack+37665,int_stack+37455, 1.0, int_stack+82300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+90530, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+29395,int_stack+113280, 1.0, int_stack+113730, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+130940, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+32295,int_stack+157005,int_stack+115480, 1.0, int_stack+116380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+133590, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+157005,int_stack+38470,int_stack+38095, 1.0, int_stack+37945, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45715, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113280,int_stack+39010,int_stack+38470, 1.0, int_stack+38245, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46240, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29395,int_stack+113280,int_stack+157005, 1.0, int_stack+119875, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+115180,int_stack+39325,int_stack+39010, 1.0, int_stack+38695, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46690, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+36090,int_stack+115180,int_stack+113280, 1.0, int_stack+120325, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137085, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+115180,int_stack+36090,int_stack+29395, 1.0, int_stack+121000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+137760, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+39845,int_stack+39745, 1.0, int_stack+84640, 0.0, zero_stack, 1.0, int_stack+89920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29395,int_stack+39995,int_stack+39845, 1.0, int_stack+84740, 0.0, zero_stack, 1.0, int_stack+90170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+29395,int_stack+116680, 1.0, int_stack+121900, 0.0, zero_stack, 1.0, int_stack+133290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+36090,int_stack+40205,int_stack+39995, 1.0, int_stack+84890, 0.0, zero_stack, 1.0, int_stack+90530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+113280,int_stack+36090,int_stack+29395, 1.0, int_stack+122200, 0.0, zero_stack, 1.0, int_stack+130940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+36090,int_stack+113280,int_stack+169255, 1.0, int_stack+122650, 0.0, zero_stack, 1.0, int_stack+133590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+41010,int_stack+40635, 1.0, int_stack+40485, 0.0, zero_stack, 1.0, int_stack+45715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+113280,int_stack+41550,int_stack+41010, 1.0, int_stack+40785, 0.0, zero_stack, 1.0, int_stack+46240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29395,int_stack+113280,int_stack+169255, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+37090,int_stack+41865,int_stack+41550, 1.0, int_stack+41235, 0.0, zero_stack, 1.0, int_stack+46690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+38035,int_stack+37090,int_stack+113280, 1.0, int_stack+450, 0.0, zero_stack, 1.0, int_stack+137085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+39385,int_stack+38035,int_stack+29395, 1.0, int_stack+123250, 0.0, zero_stack, 1.0, int_stack+137760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+42385,int_stack+42285, 1.0, int_stack+87230, 1.0, int_stack+89920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29395,int_stack+42535,int_stack+42385, 1.0, int_stack+87330, 1.0, int_stack+90170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+29395,int_stack+116680, 1.0, int_stack+126695, 1.0, int_stack+133290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+113280,int_stack+42745,int_stack+42535, 1.0, int_stack+87480, 1.0, int_stack+90530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+113280,int_stack+29395, 1.0, int_stack+126995, 1.0, int_stack+130940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+37090,int_stack+157005,int_stack+169255, 1.0, int_stack+127445, 1.0, int_stack+133590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+43550,int_stack+43175, 1.0, int_stack+43025, 1.0, int_stack+45715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+44090,int_stack+43550, 1.0, int_stack+43325, 1.0, int_stack+46240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29395,int_stack+157005,int_stack+169255, 1.0, int_stack+130490, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+126695,int_stack+44405,int_stack+44090, 1.0, int_stack+43775, 1.0, int_stack+46690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+40885,int_stack+126695,int_stack+157005, 1.0, int_stack+128045, 1.0, int_stack+137085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+126695,int_stack+40885,int_stack+29395, 1.0, int_stack+132390, 1.0, int_stack+137760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+44925,int_stack+44825, 2.0, int_stack+89920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+45075,int_stack+44925, 2.0, int_stack+90170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 2.0, int_stack+133290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+132390,int_stack+45285,int_stack+45075, 2.0, int_stack+90530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29395,int_stack+132390,int_stack+130490, 2.0, int_stack+130940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+132390,int_stack+29395,int_stack+169255, 2.0, int_stack+133590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+46465,int_stack+45865, 2.0, int_stack+45715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29395,int_stack+47005,int_stack+46465, 2.0, int_stack+46240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+29395,int_stack+169255, 2.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+40885,int_stack+47320,int_stack+47005, 2.0, int_stack+46690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41830,int_stack+40885,int_stack+29395, 2.0, int_stack+137085, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+43180,int_stack+41830,int_stack+130490, 2.0, int_stack+137760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+47840,int_stack+47740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93560,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+47990,int_stack+47840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93660,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138660,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+137085,int_stack+48200,int_stack+47990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93810,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+29395,int_stack+137085,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138960,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+137085,int_stack+29395,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106560,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+48630,int_stack+48480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57730,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29395,int_stack+48855,int_stack+48630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58030,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+29395,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101330,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+40885,int_stack+49170,int_stack+48855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58480,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+41830,int_stack+40885,int_stack+29395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101780,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+44680,int_stack+41830,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140560,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+49690,int_stack+49590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93560, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+49840,int_stack+49690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93660, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138660, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+29395,int_stack+50050,int_stack+49840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93810, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+29395,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138960, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+40885,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106560, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+50480,int_stack+50330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57730, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+50705,int_stack+50480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58030, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101330, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41885,int_stack+51020,int_stack+50705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58480, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+46180,int_stack+41885,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101780, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+47530,int_stack+46180,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140560, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+51540,int_stack+51440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93560, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+51690,int_stack+51540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93660, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138660, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+46180,int_stack+51900,int_stack+51690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93810, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+46180,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138960, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+46180,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106560, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+52330,int_stack+52180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57730, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+52555,int_stack+52330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58030, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101330, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41885,int_stack+52870,int_stack+52555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58480, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+49030,int_stack+41885,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101780, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+50380,int_stack+49030,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140560, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+53390,int_stack+53290, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+53540,int_stack+53390, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+49030,int_stack+53750,int_stack+53540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+93810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+49030,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+138960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+49030,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+54180,int_stack+54030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+54405,int_stack+54180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41885,int_stack+54720,int_stack+54405, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+51880,int_stack+41885,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+101780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+53230,int_stack+51880,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+55240,int_stack+55140, 0.0, zero_stack, 1.0, int_stack+93560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+55390,int_stack+55240, 0.0, zero_stack, 1.0, int_stack+93660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 1.0, int_stack+138660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+51880,int_stack+55600,int_stack+55390, 0.0, zero_stack, 1.0, int_stack+93810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+51880,int_stack+130490, 0.0, zero_stack, 1.0, int_stack+138960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+51880,int_stack+157005,int_stack+169255, 0.0, zero_stack, 1.0, int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+56030,int_stack+55880, 0.0, zero_stack, 1.0, int_stack+57730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+56255,int_stack+56030, 0.0, zero_stack, 1.0, int_stack+58030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 1.0, int_stack+101330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+41885,int_stack+56570,int_stack+56255, 0.0, zero_stack, 1.0, int_stack+58480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+54730,int_stack+41885,int_stack+157005, 0.0, zero_stack, 1.0, int_stack+101780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+54730,int_stack+130490, 0.0, zero_stack, 1.0, int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+57090,int_stack+56990, 1.0, int_stack+93560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+57240,int_stack+57090, 1.0, int_stack+93660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 1.0, int_stack+138660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+54730,int_stack+57450,int_stack+57240, 1.0, int_stack+93810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+54730,int_stack+130490, 1.0, int_stack+138960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+54730,int_stack+157005,int_stack+169255, 1.0, int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+106560,int_stack+58255,int_stack+57880, 1.0, int_stack+57730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+58795,int_stack+58255, 1.0, int_stack+58030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+106560, 1.0, int_stack+101330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+55730,int_stack+59110,int_stack+58795, 1.0, int_stack+58480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+56675,int_stack+55730,int_stack+157005, 1.0, int_stack+101780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+58025,int_stack+56675,int_stack+130490, 1.0, int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+59630,int_stack+59530,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+140560,int_stack+59780,int_stack+59630,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+140560,int_stack+116680,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130490,int_stack+59990,int_stack+59780,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+130490,int_stack+140560,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+55730,int_stack+157005,int_stack+106560,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+106560,int_stack+60420,int_stack+60270,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+60645,int_stack+60420,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+140560,int_stack+157005,int_stack+106560,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+56730,int_stack+60960,int_stack+60645,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+59525,int_stack+56730,int_stack+157005,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+119875,int_stack+59525,int_stack+140560,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+61480,int_stack+61380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96150,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140560,int_stack+61630,int_stack+61480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96250,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+140560,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107160,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+59525,int_stack+61840,int_stack+61630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96400,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+59525,int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141460,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+59525,int_stack+157005,int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141910,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+106560,int_stack+62270,int_stack+62120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+62495,int_stack+62270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71670,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+140560,int_stack+157005,int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145255,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+60525,int_stack+62810,int_stack+62495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72120,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+61470,int_stack+60525,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145705,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+121375,int_stack+61470,int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146380,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+63330,int_stack+63230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96150, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140560,int_stack+63480,int_stack+63330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96250, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+140560,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107160, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+157005,int_stack+63690,int_stack+63480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96400, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141460, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+60525,int_stack+130490,int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141910, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+106560,int_stack+64120,int_stack+63970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+64345,int_stack+64120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71670, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+140560,int_stack+130490,int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145255, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+61525,int_stack+64660,int_stack+64345, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72120, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+62470,int_stack+61525,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145705, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+107460,int_stack+62470,int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146380, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+65180,int_stack+65080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96150, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140560,int_stack+65330,int_stack+65180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96250, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+140560,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107160, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+130490,int_stack+65540,int_stack+65330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96400, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+130490,int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141460, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+61525,int_stack+157005,int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141910, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+106560,int_stack+65970,int_stack+65820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+66195,int_stack+65970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71670, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+140560,int_stack+157005,int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145255, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+62525,int_stack+66510,int_stack+66195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72120, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+63470,int_stack+62525,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145705, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+64820,int_stack+63470,int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146380, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+67030,int_stack+66930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140560,int_stack+67180,int_stack+67030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+140560,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+107160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+157005,int_stack+67390,int_stack+67180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+96400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+62525,int_stack+130490,int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+141910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+106560,int_stack+67820,int_stack+67670, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+68045,int_stack+67820, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+140560,int_stack+130490,int_stack+106560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+63525,int_stack+68360,int_stack+68045, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+66320,int_stack+63525,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+145705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+108960,int_stack+66320,int_stack+140560, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+146380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+68880,int_stack+68780, 0.0, zero_stack, 1.0, int_stack+96150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140560,int_stack+69030,int_stack+68880, 0.0, zero_stack, 1.0, int_stack+96250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+140560,int_stack+116680, 0.0, zero_stack, 1.0, int_stack+107160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+66320,int_stack+69240,int_stack+69030, 0.0, zero_stack, 1.0, int_stack+96400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+130490,int_stack+66320,int_stack+140560, 0.0, zero_stack, 1.0, int_stack+141460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+66320,int_stack+130490,int_stack+106560, 0.0, zero_stack, 1.0, int_stack+141910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+106560,int_stack+69670,int_stack+69520, 0.0, zero_stack, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+69895,int_stack+69670, 0.0, zero_stack, 1.0, int_stack+71670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+140560,int_stack+130490,int_stack+106560, 0.0, zero_stack, 1.0, int_stack+145255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+67320,int_stack+70210,int_stack+69895, 0.0, zero_stack, 1.0, int_stack+72120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+68265,int_stack+67320,int_stack+130490, 0.0, zero_stack, 1.0, int_stack+145705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+172355,int_stack+68265,int_stack+140560, 0.0, zero_stack, 1.0, int_stack+146380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+70730,int_stack+70630, 1.0, int_stack+96150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+140560,int_stack+70880,int_stack+70730, 1.0, int_stack+96250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+106560,int_stack+140560,int_stack+116680, 1.0, int_stack+107160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+130490,int_stack+71090,int_stack+70880, 1.0, int_stack+96400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+130490,int_stack+140560, 1.0, int_stack+141460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+140560,int_stack+157005,int_stack+106560, 1.0, int_stack+141910, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+106560,int_stack+71895,int_stack+71520, 1.0, int_stack+71370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+72435,int_stack+71895, 1.0, int_stack+71670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+106560, 1.0, int_stack+145255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+141560,int_stack+72750,int_stack+72435, 1.0, int_stack+72120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+67320,int_stack+141560,int_stack+157005, 1.0, int_stack+145705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+68670,int_stack+67320,int_stack+130490, 1.0, int_stack+146380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+73270,int_stack+73170,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+73420,int_stack+73270,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+67320,int_stack+73630,int_stack+73420,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+67320,int_stack+130490,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+67320,int_stack+157005,int_stack+169255,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+74060,int_stack+73910,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+74285,int_stack+74060,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+141560,int_stack+74600,int_stack+74285,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+145255,int_stack+141560,int_stack+157005,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+70170,int_stack+145255,int_stack+130490,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+75120,int_stack+75020,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+75270,int_stack+75120,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+145255,int_stack+75480,int_stack+75270,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+145255,int_stack+130490,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+145255,int_stack+157005,int_stack+169255,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+75910,int_stack+75760,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+76135,int_stack+75910,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+146255,int_stack+76450,int_stack+76135,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+71670,int_stack+146255,int_stack+157005,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+73020,int_stack+71670,int_stack+130490,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+77710,int_stack+77610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98740,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+77860,int_stack+77710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98840,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2575,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+71670,int_stack+78070,int_stack+77860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98990,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+71670,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2875,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+71670,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3325,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+78500,int_stack+78350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91760,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+78725,int_stack+78500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92060,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+148780,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+146255,int_stack+79040,int_stack+78725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92510,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+74520,int_stack+146255,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149230,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+75870,int_stack+74520,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151310,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+80300,int_stack+80200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98740, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+80450,int_stack+80300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98840, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2575, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+74520,int_stack+80660,int_stack+80450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98990, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+74520,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2875, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+74520,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3325, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+81090,int_stack+80940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91760, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+81315,int_stack+81090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92060, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+148780, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+146255,int_stack+81630,int_stack+81315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92510, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+77370,int_stack+146255,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149230, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+78720,int_stack+77370,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151310, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+82890,int_stack+82790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98740, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+83040,int_stack+82890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98840, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2575, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+77370,int_stack+83250,int_stack+83040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98990, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+77370,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2875, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+77370,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3325, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+83680,int_stack+83530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91760, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+83905,int_stack+83680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92060, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+148780, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+146255,int_stack+84220,int_stack+83905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92510, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+80220,int_stack+146255,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149230, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+81570,int_stack+80220,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151310, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+85480,int_stack+85380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+85630,int_stack+85480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+80220,int_stack+85840,int_stack+85630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+98990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+80220,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+80220,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+86270,int_stack+86120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+91760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+86495,int_stack+86270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+148780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+146255,int_stack+86810,int_stack+86495, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+92510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+83070,int_stack+146255,int_stack+157005, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+149230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+84420,int_stack+83070,int_stack+130490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+151310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+88070,int_stack+87970, 0.0, zero_stack, 1.0, int_stack+98740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+88220,int_stack+88070, 0.0, zero_stack, 1.0, int_stack+98840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 0.0, zero_stack, 1.0, int_stack+2575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+83070,int_stack+88430,int_stack+88220, 0.0, zero_stack, 1.0, int_stack+98990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+83070,int_stack+130490, 0.0, zero_stack, 1.0, int_stack+2875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+83070,int_stack+157005,int_stack+169255, 0.0, zero_stack, 1.0, int_stack+3325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+88860,int_stack+88710, 0.0, zero_stack, 1.0, int_stack+91760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+89085,int_stack+88860, 0.0, zero_stack, 1.0, int_stack+92060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 0.0, zero_stack, 1.0, int_stack+148780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+146255,int_stack+89400,int_stack+89085, 0.0, zero_stack, 1.0, int_stack+92510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+85920,int_stack+146255,int_stack+157005, 0.0, zero_stack, 1.0, int_stack+149230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+87270,int_stack+85920,int_stack+130490, 0.0, zero_stack, 1.0, int_stack+151310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+91120,int_stack+91020, 1.0, int_stack+98740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+91270,int_stack+91120, 1.0, int_stack+98840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+130490,int_stack+116680, 1.0, int_stack+2575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+85920,int_stack+91480,int_stack+91270, 1.0, int_stack+98990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+85920,int_stack+130490, 1.0, int_stack+2875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+85920,int_stack+157005,int_stack+169255, 1.0, int_stack+3325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+92285,int_stack+91910, 1.0, int_stack+91760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+92825,int_stack+92285, 1.0, int_stack+92060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+169255, 1.0, int_stack+148780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2575,int_stack+93140,int_stack+92825, 1.0, int_stack+92510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+88770,int_stack+2575,int_stack+157005, 1.0, int_stack+149230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+148780,int_stack+88770,int_stack+130490, 1.0, int_stack+151310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+94400,int_stack+94300,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+151310,int_stack+94550,int_stack+94400,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+151310,int_stack+116680,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+130490,int_stack+94760,int_stack+94550,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+130490,int_stack+151310,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+88770,int_stack+157005,int_stack+169255,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+95190,int_stack+95040,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+95415,int_stack+95190,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+151310,int_stack+157005,int_stack+169255,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+89770,int_stack+95730,int_stack+95415,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2575,int_stack+89770,int_stack+157005,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+89770,int_stack+2575,int_stack+151310,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+96990,int_stack+96890,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+151310,int_stack+97140,int_stack+96990,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+151310,int_stack+116680,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2575,int_stack+97350,int_stack+97140,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+157005,int_stack+2575,int_stack+151310,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+2575,int_stack+157005,int_stack+169255,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+97780,int_stack+97630,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+157005,int_stack+98005,int_stack+97780,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+151310,int_stack+157005,int_stack+169255,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+91270,int_stack+98320,int_stack+98005,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+92215,int_stack+91270,int_stack+157005,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+93565,int_stack+92215,int_stack+151310,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+116680,int_stack+99580,int_stack+99480,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+151310,int_stack+99730,int_stack+99580,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+169255,int_stack+151310,int_stack+116680,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+157005,int_stack+99940,int_stack+99730,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+130490,int_stack+157005,int_stack+151310,10);
 /*--- compute (f0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+91270,int_stack+130490,int_stack+169255,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+169255,int_stack+100370,int_stack+100220,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+130490,int_stack+100595,int_stack+100370,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+151310,int_stack+130490,int_stack+169255,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+92270,int_stack+100910,int_stack+100595,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+95065,int_stack+92270,int_stack+130490,15);
 /*--- compute (g0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+96415,int_stack+95065,int_stack+151310,15);
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+97915,int_stack+111780,int_stack+105560,100);
     Libderiv->ABCD[11] = int_stack + 97915;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+110460,int_stack+118375,int_stack+114180,100);
     Libderiv->ABCD[10] = int_stack + 110460;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+173855,int_stack+124195,int_stack+116980,100);
     Libderiv->ABCD[9] = int_stack + 173855;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+176855,int_stack+128990,int_stack+125695,100);
     Libderiv->ABCD[8] = int_stack + 176855;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+128195,int_stack+135585,int_stack+131390,100);
     Libderiv->ABCD[7] = int_stack + 128195;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+179855,int_stack+104030,int_stack+134190,100);
     Libderiv->ABCD[6] = int_stack + 179855;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+182855,int_stack+143755,int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+139560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[2] = int_stack + 182855;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+185855,int_stack+147280,int_stack+142510, 0.0, zero_stack, 1.0, int_stack+139560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[1] = int_stack + 185855;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+188855,int_stack+154505,int_stack+150310, 1.0, int_stack+139560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[0] = int_stack + 188855;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+191855,int_stack+3925,int_stack+156005,100);
     Libderiv->ABCD[155] = int_stack + 191855;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+153210,int_stack+158905,int_stack+157905,100);
     Libderiv->ABCD[143] = int_stack + 153210;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+156210,int_stack+160405,int_stack+5425,100);
     Libderiv->ABCD[142] = int_stack + 156210;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+194855,int_stack+9720,int_stack+6425,100);
     Libderiv->ABCD[131] = int_stack + 194855;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3575,int_stack+161905,int_stack+7425,100);
     Libderiv->ABCD[130] = int_stack + 3575;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+6575,int_stack+13570,int_stack+12570,100);
     Libderiv->ABCD[129] = int_stack + 6575;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+159210,int_stack+164755,int_stack+15070,100);
     Libderiv->ABCD[119] = int_stack + 159210;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12220,int_stack+17420,int_stack+163405,100);
     Libderiv->ABCD[118] = int_stack + 12220;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+162210,int_stack+166255,int_stack+152210,100);
     Libderiv->ABCD[117] = int_stack + 162210;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+197855,int_stack+167755,int_stack+18920,100);
     Libderiv->ABCD[116] = int_stack + 197855;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+165210,int_stack+22495,int_stack+19920,100);
     Libderiv->ABCD[107] = int_stack + 165210;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+17070,int_stack+24895,int_stack+20920,100);
     Libderiv->ABCD[106] = int_stack + 17070;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+20070,int_stack+26395,int_stack+11220,100);
     Libderiv->ABCD[105] = int_stack + 20070;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+23070,int_stack+27895,int_stack+169855,100);
     Libderiv->ABCD[104] = int_stack + 23070;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+26070,int_stack+170855,int_stack+16070,100);
     Libderiv->ABCD[103] = int_stack + 26070;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+168210,int_stack+102455,int_stack+30295,100);
     Libderiv->ABCD[95] = int_stack + 168210;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+100915,int_stack+34590,int_stack+31295,100);
     Libderiv->ABCD[94] = int_stack + 100915;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+29070,int_stack+115180,int_stack+32295,100);
     Libderiv->ABCD[93] = int_stack + 29070;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+32070,int_stack+39385,int_stack+36090,100);
     Libderiv->ABCD[92] = int_stack + 32070;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+200855,int_stack+126695,int_stack+37090,100);
     Libderiv->ABCD[91] = int_stack + 200855;
 /*--- compute (fp|ff) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+35070,int_stack+43180,int_stack+132390,100);
     Libderiv->ABCD[90] = int_stack + 35070;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+203855,int_stack+44680,int_stack+137085, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+105560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[47] = int_stack + 203855;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+41885,int_stack+47530,int_stack+40885, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+114180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[46] = int_stack + 41885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+38070,int_stack+50380,int_stack+46180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+116980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[45] = int_stack + 38070;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+44885,int_stack+53230,int_stack+49030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+125695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[44] = int_stack + 44885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+47885,int_stack+0,int_stack+51880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+131390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[43] = int_stack + 47885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+50885,int_stack+58025,int_stack+54730, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+134190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[42] = int_stack + 50885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+135190,int_stack+119875,int_stack+55730, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[38] = int_stack + 135190;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+53885,int_stack+121375,int_stack+59525, 0.0, zero_stack, 1.0, int_stack+105560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[35] = int_stack + 53885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+56885,int_stack+107460,int_stack+60525, 0.0, zero_stack, 1.0, int_stack+114180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[34] = int_stack + 56885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+117980,int_stack+64820,int_stack+61525, 0.0, zero_stack, 1.0, int_stack+116980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[33] = int_stack + 117980;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+120980,int_stack+108960,int_stack+62525, 0.0, zero_stack, 1.0, int_stack+125695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[32] = int_stack + 120980;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+106560,int_stack+172355,int_stack+66320, 0.0, zero_stack, 1.0, int_stack+131390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[31] = int_stack + 106560;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+59885,int_stack+68670,int_stack+140560, 0.0, zero_stack, 1.0, int_stack+134190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[30] = int_stack + 59885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+62885,int_stack+70170,int_stack+67320, 0.0, zero_stack, 1.0, int_stack+1575, 1.0, int_stack+142510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[26] = int_stack + 62885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+65885,int_stack+73020,int_stack+145255, 0.0, zero_stack, 2.0, int_stack+142510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[25] = int_stack + 65885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+143510,int_stack+75870,int_stack+71670, 1.0, int_stack+105560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[23] = int_stack + 143510;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+68885,int_stack+78720,int_stack+74520, 1.0, int_stack+114180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[22] = int_stack + 68885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+71885,int_stack+81570,int_stack+77370, 1.0, int_stack+116980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[21] = int_stack + 71885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+74885,int_stack+84420,int_stack+80220, 1.0, int_stack+125695, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[20] = int_stack + 74885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+77885,int_stack+87270,int_stack+83070, 1.0, int_stack+131390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[19] = int_stack + 77885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+80885,int_stack+148780,int_stack+85920, 1.0, int_stack+134190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[18] = int_stack + 80885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+83885,int_stack+89770,int_stack+88770, 1.0, int_stack+1575, 0.0, zero_stack, 1.0, int_stack+150310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[14] = int_stack + 83885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+86885,int_stack+93565,int_stack+2575, 1.0, int_stack+142510, 1.0, int_stack+150310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[13] = int_stack + 86885;
 /*--- compute (fp|ff) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+96415,int_stack+91270, 2.0, int_stack+150310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,100);
     Libderiv->ABCD[12] = int_stack + 0;

}
