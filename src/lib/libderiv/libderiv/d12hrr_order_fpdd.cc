#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_fpdd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|dd) integrals */

void d12hrr_order_fpdd(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv2_classes[3][2][143] = int_stack + 2325;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 2385;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 2485;
 Libderiv->deriv2_classes[4][2][143] = int_stack + 2635;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 2725;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 2875;
 Libderiv->deriv2_classes[3][2][131] = int_stack + 3100;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 3160;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 3260;
 Libderiv->deriv2_classes[4][2][131] = int_stack + 3410;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 3500;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 3650;
 Libderiv->deriv2_classes[3][2][130] = int_stack + 3875;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 3935;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 4035;
 Libderiv->deriv2_classes[4][2][130] = int_stack + 4185;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 4275;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 4425;
 Libderiv->deriv2_classes[3][2][119] = int_stack + 4650;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 4710;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 4810;
 Libderiv->deriv2_classes[4][2][119] = int_stack + 4960;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 5050;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 5200;
 Libderiv->deriv2_classes[3][2][118] = int_stack + 5425;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 5485;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 5585;
 Libderiv->deriv2_classes[4][2][118] = int_stack + 5735;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 5825;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 5975;
 Libderiv->deriv2_classes[3][2][117] = int_stack + 6200;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 6260;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 6360;
 Libderiv->deriv2_classes[4][2][117] = int_stack + 6510;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 6600;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 6750;
 Libderiv->deriv2_classes[3][2][107] = int_stack + 6975;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 7035;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 7135;
 Libderiv->deriv2_classes[4][2][107] = int_stack + 7285;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 7375;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 7525;
 Libderiv->deriv2_classes[3][2][106] = int_stack + 7750;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 7810;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 7910;
 Libderiv->deriv2_classes[4][2][106] = int_stack + 8060;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 8150;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 8300;
 Libderiv->deriv2_classes[3][2][105] = int_stack + 8525;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 8585;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 8685;
 Libderiv->deriv2_classes[4][2][105] = int_stack + 8835;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 8925;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 9075;
 Libderiv->deriv2_classes[3][2][104] = int_stack + 9300;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 9360;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 9460;
 Libderiv->deriv2_classes[4][2][104] = int_stack + 9610;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 9700;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 9850;
 Libderiv->deriv2_classes[3][2][95] = int_stack + 10075;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 10135;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 10235;
 Libderiv->deriv2_classes[4][2][95] = int_stack + 10385;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 10475;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 10625;
 Libderiv->deriv2_classes[3][2][94] = int_stack + 10850;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 10910;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 11010;
 Libderiv->deriv2_classes[4][2][94] = int_stack + 11160;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 11250;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 11400;
 Libderiv->deriv2_classes[3][2][93] = int_stack + 11625;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 11685;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 11785;
 Libderiv->deriv2_classes[4][2][93] = int_stack + 11935;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 12025;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 12175;
 Libderiv->deriv2_classes[3][2][92] = int_stack + 12400;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 12460;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 12560;
 Libderiv->deriv2_classes[4][2][92] = int_stack + 12710;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 12800;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 12950;
 Libderiv->deriv2_classes[3][2][91] = int_stack + 13175;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 13235;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 13335;
 Libderiv->deriv2_classes[4][2][91] = int_stack + 13485;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 13575;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 13725;
 Libderiv->deriv2_classes[3][2][83] = int_stack + 13950;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 14010;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 14110;
 Libderiv->deriv_classes[4][2][11] = int_stack + 14260;
 Libderiv->deriv2_classes[4][2][83] = int_stack + 14350;
 Libderiv->deriv_classes[4][3][11] = int_stack + 14440;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 14590;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 14740;
 Libderiv->deriv2_classes[3][2][82] = int_stack + 14965;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 15025;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 15125;
 Libderiv->deriv_classes[4][2][10] = int_stack + 15275;
 Libderiv->deriv2_classes[4][2][82] = int_stack + 15365;
 Libderiv->deriv_classes[4][3][10] = int_stack + 15455;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 15605;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 15755;
 Libderiv->deriv2_classes[3][2][81] = int_stack + 15980;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 16040;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 16140;
 Libderiv->deriv_classes[4][2][9] = int_stack + 16290;
 Libderiv->deriv2_classes[4][2][81] = int_stack + 16380;
 Libderiv->deriv_classes[4][3][9] = int_stack + 16470;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 16620;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 16770;
 Libderiv->deriv2_classes[3][2][80] = int_stack + 16995;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 17055;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 17155;
 Libderiv->deriv_classes[4][2][8] = int_stack + 17305;
 Libderiv->deriv2_classes[4][2][80] = int_stack + 17395;
 Libderiv->deriv_classes[4][3][8] = int_stack + 17485;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 17635;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 17785;
 Libderiv->deriv2_classes[3][2][79] = int_stack + 18010;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 18070;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 18170;
 Libderiv->deriv_classes[4][2][7] = int_stack + 18320;
 Libderiv->deriv2_classes[4][2][79] = int_stack + 18410;
 Libderiv->deriv_classes[4][3][7] = int_stack + 18500;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 18650;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 18800;
 Libderiv->deriv2_classes[3][2][78] = int_stack + 19025;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 19085;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 19185;
 Libderiv->dvrr_classes[4][2] = int_stack + 19335;
 Libderiv->deriv_classes[4][2][6] = int_stack + 19425;
 Libderiv->deriv2_classes[4][2][78] = int_stack + 19515;
 Libderiv->deriv_classes[4][3][6] = int_stack + 19605;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 19755;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 19905;
 Libderiv->deriv2_classes[3][2][35] = int_stack + 20130;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 20190;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 20290;
 Libderiv->deriv2_classes[4][2][35] = int_stack + 20440;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 20530;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 20680;
 Libderiv->deriv2_classes[3][2][34] = int_stack + 20905;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 20965;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 21065;
 Libderiv->deriv2_classes[4][2][34] = int_stack + 21215;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 21305;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 21455;
 Libderiv->deriv2_classes[3][2][33] = int_stack + 21680;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 21740;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 21840;
 Libderiv->deriv2_classes[4][2][33] = int_stack + 21990;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 22080;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 22230;
 Libderiv->deriv2_classes[3][2][32] = int_stack + 22455;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 22515;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 22615;
 Libderiv->deriv2_classes[4][2][32] = int_stack + 22765;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 22855;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 23005;
 Libderiv->deriv2_classes[3][2][31] = int_stack + 23230;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 23290;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 23390;
 Libderiv->deriv2_classes[4][2][31] = int_stack + 23540;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 23630;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 23780;
 Libderiv->deriv2_classes[3][2][30] = int_stack + 24005;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 24065;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 24165;
 Libderiv->deriv_classes[4][2][2] = int_stack + 24315;
 Libderiv->deriv2_classes[4][2][30] = int_stack + 24405;
 Libderiv->deriv_classes[4][3][2] = int_stack + 24495;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 24645;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 24795;
 Libderiv->deriv2_classes[3][2][26] = int_stack + 25020;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 25080;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 25180;
 Libderiv->deriv2_classes[4][2][26] = int_stack + 25330;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 25420;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 25570;
 Libderiv->deriv2_classes[3][2][23] = int_stack + 25795;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 25855;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 25955;
 Libderiv->deriv2_classes[4][2][23] = int_stack + 26105;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 26195;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 26345;
 Libderiv->deriv2_classes[3][2][22] = int_stack + 26570;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 26630;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 26730;
 Libderiv->deriv2_classes[4][2][22] = int_stack + 26880;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 26970;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 27120;
 Libderiv->deriv2_classes[3][2][21] = int_stack + 27345;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 27405;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 27505;
 Libderiv->deriv2_classes[4][2][21] = int_stack + 27655;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 27745;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 27895;
 Libderiv->deriv2_classes[3][2][20] = int_stack + 28120;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 28180;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 28280;
 Libderiv->deriv2_classes[4][2][20] = int_stack + 28430;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 28520;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 28670;
 Libderiv->deriv2_classes[3][2][19] = int_stack + 28895;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 28955;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 29055;
 Libderiv->deriv2_classes[4][2][19] = int_stack + 29205;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 29295;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 29445;
 Libderiv->deriv2_classes[3][2][18] = int_stack + 29670;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 29730;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 29830;
 Libderiv->deriv_classes[4][2][1] = int_stack + 29980;
 Libderiv->deriv2_classes[4][2][18] = int_stack + 30070;
 Libderiv->deriv_classes[4][3][1] = int_stack + 30160;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 30310;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 30460;
 Libderiv->deriv2_classes[3][2][14] = int_stack + 30685;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 30745;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 30845;
 Libderiv->deriv2_classes[4][2][14] = int_stack + 30995;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 31085;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 31235;
 Libderiv->deriv2_classes[3][2][13] = int_stack + 31460;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 31520;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 31620;
 Libderiv->deriv2_classes[4][2][13] = int_stack + 31770;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 31860;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 32010;
 Libderiv->deriv_classes[3][2][11] = int_stack + 32235;
 Libderiv->deriv_classes[3][3][11] = int_stack + 32295;
 Libderiv->deriv_classes[3][4][11] = int_stack + 32395;
 Libderiv->deriv2_classes[3][2][11] = int_stack + 32545;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 32605;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 32705;
 Libderiv->deriv2_classes[4][2][11] = int_stack + 32855;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 32945;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 33095;
 Libderiv->deriv_classes[3][2][10] = int_stack + 33320;
 Libderiv->deriv_classes[3][3][10] = int_stack + 33380;
 Libderiv->deriv_classes[3][4][10] = int_stack + 33480;
 Libderiv->deriv2_classes[3][2][10] = int_stack + 33630;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 33690;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 33790;
 Libderiv->deriv2_classes[4][2][10] = int_stack + 33940;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 34030;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 34180;
 Libderiv->deriv_classes[3][2][9] = int_stack + 34405;
 Libderiv->deriv_classes[3][3][9] = int_stack + 34465;
 Libderiv->deriv_classes[3][4][9] = int_stack + 34565;
 Libderiv->deriv2_classes[3][2][9] = int_stack + 34715;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 34775;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 34875;
 Libderiv->deriv2_classes[4][2][9] = int_stack + 35025;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 35115;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 35265;
 Libderiv->deriv_classes[3][2][8] = int_stack + 35490;
 Libderiv->deriv_classes[3][3][8] = int_stack + 35550;
 Libderiv->deriv_classes[3][4][8] = int_stack + 35650;
 Libderiv->deriv2_classes[3][2][8] = int_stack + 35800;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 35860;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 35960;
 Libderiv->deriv2_classes[4][2][8] = int_stack + 36110;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 36200;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 36350;
 Libderiv->deriv_classes[3][2][7] = int_stack + 36575;
 Libderiv->deriv_classes[3][3][7] = int_stack + 36635;
 Libderiv->deriv_classes[3][4][7] = int_stack + 36735;
 Libderiv->deriv2_classes[3][2][7] = int_stack + 36885;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 36945;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 37045;
 Libderiv->deriv2_classes[4][2][7] = int_stack + 37195;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 37285;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 37435;
 Libderiv->dvrr_classes[3][2] = int_stack + 37660;
 Libderiv->deriv_classes[3][2][6] = int_stack + 37720;
 Libderiv->dvrr_classes[3][3] = int_stack + 37780;
 Libderiv->deriv_classes[3][3][6] = int_stack + 37880;
 Libderiv->deriv_classes[3][4][6] = int_stack + 37980;
 Libderiv->deriv2_classes[3][2][6] = int_stack + 38130;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 38190;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 38290;
 Libderiv->deriv_classes[4][2][0] = int_stack + 38440;
 Libderiv->deriv2_classes[4][2][6] = int_stack + 38530;
 Libderiv->deriv_classes[4][3][0] = int_stack + 38620;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 38770;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 38920;
 Libderiv->deriv_classes[3][2][2] = int_stack + 39145;
 Libderiv->deriv_classes[3][3][2] = int_stack + 39205;
 Libderiv->deriv_classes[3][4][2] = int_stack + 39305;
 Libderiv->deriv2_classes[3][2][2] = int_stack + 39455;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 39515;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 39615;
 Libderiv->deriv2_classes[4][2][2] = int_stack + 39765;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 39855;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 40005;
 Libderiv->deriv_classes[3][2][1] = int_stack + 40230;
 Libderiv->deriv_classes[3][3][1] = int_stack + 40290;
 Libderiv->deriv_classes[3][4][1] = int_stack + 40390;
 Libderiv->deriv2_classes[3][2][1] = int_stack + 40540;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 40600;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 40700;
 Libderiv->deriv2_classes[4][2][1] = int_stack + 40850;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 40940;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 41090;
 Libderiv->deriv_classes[3][2][0] = int_stack + 41315;
 Libderiv->deriv_classes[3][3][0] = int_stack + 41375;
 Libderiv->deriv_classes[3][4][0] = int_stack + 41475;
 Libderiv->deriv2_classes[3][2][0] = int_stack + 41625;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 41685;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 41785;
 Libderiv->deriv2_classes[4][2][0] = int_stack + 41935;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 42025;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 42175;
 memset(int_stack,0,339200);

 Libderiv->dvrr_stack = int_stack + 68500;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_fpdd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+42400,int_stack+37780,int_stack+37660,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+42580,int_stack+32295,int_stack+32235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37660,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42760,int_stack+32395,int_stack+32295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37780,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+43060,int_stack+42760,int_stack+42580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+42760,int_stack+1125,int_stack+19335,15);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+43420,int_stack+14440,int_stack+14260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19335,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43690,int_stack+0,int_stack+14440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+44140,int_stack+43690,int_stack+43420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42760,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+43690,int_stack+33380,int_stack+33320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37660, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+33480,int_stack+33380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37780, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+44980,int_stack+44680,int_stack+43690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+43870,int_stack+15455,int_stack+15275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19335, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45340,int_stack+225,int_stack+15455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+45790,int_stack+45340,int_stack+43870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42760, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+45340,int_stack+34465,int_stack+34405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37660, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+34565,int_stack+34465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37780, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+44680,int_stack+45340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+45520,int_stack+16470,int_stack+16290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19335, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+46330,int_stack+450,int_stack+16470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+46780,int_stack+46330,int_stack+45520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42760, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+46330,int_stack+35550,int_stack+35490, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+35650,int_stack+35550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+47320,int_stack+44680,int_stack+46330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+46510,int_stack+17485,int_stack+17305, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47680,int_stack+675,int_stack+17485, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+360,int_stack+47680,int_stack+46510, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+36635,int_stack+36575, 0.0, zero_stack, 1.0, int_stack+37660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+36735,int_stack+36635, 0.0, zero_stack, 1.0, int_stack+37780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+47860,int_stack+44680,int_stack+47680, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+18500,int_stack+18320, 0.0, zero_stack, 1.0, int_stack+19335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+48220,int_stack+900,int_stack+18500, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+48670,int_stack+48220,int_stack+44680, 0.0, zero_stack, 1.0, int_stack+42760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+48220,int_stack+37880,int_stack+37720, 1.0, int_stack+37660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+37980,int_stack+37880, 1.0, int_stack+37780, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+49510,int_stack+49210,int_stack+48220, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+48400,int_stack+19605,int_stack+19425, 1.0, int_stack+19335, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49870,int_stack+1275,int_stack+19605, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+900,int_stack+49870,int_stack+48400, 1.0, int_stack+42760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+42760,int_stack+1950,int_stack+37780,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+49870,int_stack+42760,int_stack+42400,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+42400,int_stack+39205,int_stack+39145,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+42760,int_stack+39305,int_stack+39205,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+50230,int_stack+42760,int_stack+42400,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+42760,int_stack+24495,int_stack+24315,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+50590,int_stack+1500,int_stack+24495,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+51040,int_stack+50590,int_stack+42760,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+50590,int_stack+40290,int_stack+40230,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+40390,int_stack+40290,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+51580,int_stack+49210,int_stack+50590,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+50770,int_stack+30160,int_stack+29980,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51940,int_stack+1725,int_stack+30160,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+52390,int_stack+51940,int_stack+50770,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+51940,int_stack+41375,int_stack+41315,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+41475,int_stack+41375,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1440,int_stack+49210,int_stack+51940,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+52120,int_stack+38620,int_stack+38440,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52930,int_stack+2100,int_stack+38620,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+53380,int_stack+52930,int_stack+52120,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+2385,int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+32235,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+2485,int_stack+2385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+32295,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1800,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+42580,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+2725,int_stack+2635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14260,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2160,int_stack+2875,int_stack+2725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14440,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+53920,int_stack+2160,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43420,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+3160,int_stack+3100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32235, 1.0, int_stack+33320,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+3260,int_stack+3160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32295, 1.0, int_stack+33380,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2160,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42580, 1.0, int_stack+43690,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+3500,int_stack+3410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14260, 1.0, int_stack+15275,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2520,int_stack+3650,int_stack+3500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14440, 1.0, int_stack+15455,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2970,int_stack+2520,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43420, 1.0, int_stack+43870,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+3935,int_stack+3875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+33320, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+4035,int_stack+3935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+33380, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2520,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43690, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+4275,int_stack+4185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15275, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3510,int_stack+4425,int_stack+4275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15455, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3960,int_stack+3510,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+43870, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+4710,int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32235, 0.0, zero_stack, 1.0, int_stack+34405,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+4810,int_stack+4710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32295, 0.0, zero_stack, 1.0, int_stack+34465,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3510,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42580, 0.0, zero_stack, 1.0, int_stack+45340,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+5050,int_stack+4960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14260, 0.0, zero_stack, 1.0, int_stack+16290,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+5200,int_stack+5050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14440, 0.0, zero_stack, 1.0, int_stack+16470,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+54460,int_stack+4500,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43420, 0.0, zero_stack, 1.0, int_stack+45520,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+5485,int_stack+5425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33320, 1.0, int_stack+34405, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+5585,int_stack+5485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33380, 1.0, int_stack+34465, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4500,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43690, 1.0, int_stack+45340, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+5825,int_stack+5735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15275, 1.0, int_stack+16290, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4860,int_stack+5975,int_stack+5825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15455, 1.0, int_stack+16470, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5310,int_stack+4860,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43870, 1.0, int_stack+45520, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+6260,int_stack+6200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+34405, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+6360,int_stack+6260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+34465, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4860,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+45340, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+6600,int_stack+6510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16290, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5850,int_stack+6750,int_stack+6600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16470, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6300,int_stack+5850,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+45520, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+7035,int_stack+6975, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35490,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+7135,int_stack+7035, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32295, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35550,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+5850,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42580, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46330,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+7375,int_stack+7285, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17305,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6840,int_stack+7525,int_stack+7375, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17485,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+55000,int_stack+6840,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46510,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+7810,int_stack+7750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33320, 0.0, zero_stack, 1.0, int_stack+35490, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+7910,int_stack+7810, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33380, 0.0, zero_stack, 1.0, int_stack+35550, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+6840,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43690, 0.0, zero_stack, 1.0, int_stack+46330, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+8150,int_stack+8060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15275, 0.0, zero_stack, 1.0, int_stack+17305, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+8300,int_stack+8150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15455, 0.0, zero_stack, 1.0, int_stack+17485, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7650,int_stack+7200,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43870, 0.0, zero_stack, 1.0, int_stack+46510, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+8585,int_stack+8525, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34405, 1.0, int_stack+35490, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+8685,int_stack+8585, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34465, 1.0, int_stack+35550, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+7200,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45340, 1.0, int_stack+46330, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+8925,int_stack+8835, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16290, 1.0, int_stack+17305, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8190,int_stack+9075,int_stack+8925, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16470, 1.0, int_stack+17485, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8640,int_stack+8190,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45520, 1.0, int_stack+46510, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+9360,int_stack+9300, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+35490, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+9460,int_stack+9360, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+35550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+8190,int_stack+49210,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+46330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+9700,int_stack+9610, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9180,int_stack+9850,int_stack+9700, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+55540,int_stack+9180,int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+46510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+10135,int_stack+10075, 0.0, zero_stack, 1.0, int_stack+32235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36575,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+10235,int_stack+10135, 0.0, zero_stack, 1.0, int_stack+32295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36635,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9180,int_stack+49210,int_stack+52930, 0.0, zero_stack, 1.0, int_stack+42580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47680,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+10475,int_stack+10385, 0.0, zero_stack, 1.0, int_stack+14260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18320,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9540,int_stack+10625,int_stack+10475, 0.0, zero_stack, 1.0, int_stack+14440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18500,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9990,int_stack+9540,int_stack+52930, 0.0, zero_stack, 1.0, int_stack+43420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44680,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+10910,int_stack+10850, 0.0, zero_stack, 1.0, int_stack+33320, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36575, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+11010,int_stack+10910, 0.0, zero_stack, 1.0, int_stack+33380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36635, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+9540,int_stack+49210,int_stack+52930, 0.0, zero_stack, 1.0, int_stack+43690, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47680, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+11250,int_stack+11160, 0.0, zero_stack, 1.0, int_stack+15275, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18320, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10530,int_stack+11400,int_stack+11250, 0.0, zero_stack, 1.0, int_stack+15455, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18500, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10980,int_stack+10530,int_stack+52930, 0.0, zero_stack, 1.0, int_stack+43870, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44680, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+11685,int_stack+11625, 0.0, zero_stack, 1.0, int_stack+34405, 0.0, zero_stack, 1.0, int_stack+36575, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+11785,int_stack+11685, 0.0, zero_stack, 1.0, int_stack+34465, 0.0, zero_stack, 1.0, int_stack+36635, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+10530,int_stack+49210,int_stack+52930, 0.0, zero_stack, 1.0, int_stack+45340, 0.0, zero_stack, 1.0, int_stack+47680, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+12025,int_stack+11935, 0.0, zero_stack, 1.0, int_stack+16290, 0.0, zero_stack, 1.0, int_stack+18320, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11520,int_stack+12175,int_stack+12025, 0.0, zero_stack, 1.0, int_stack+16470, 0.0, zero_stack, 1.0, int_stack+18500, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+56080,int_stack+11520,int_stack+52930, 0.0, zero_stack, 1.0, int_stack+45520, 0.0, zero_stack, 1.0, int_stack+44680, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+12460,int_stack+12400, 0.0, zero_stack, 1.0, int_stack+35490, 1.0, int_stack+36575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+12560,int_stack+12460, 0.0, zero_stack, 1.0, int_stack+35550, 1.0, int_stack+36635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11520,int_stack+49210,int_stack+52930, 0.0, zero_stack, 1.0, int_stack+46330, 1.0, int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+12800,int_stack+12710, 0.0, zero_stack, 1.0, int_stack+17305, 1.0, int_stack+18320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11880,int_stack+12950,int_stack+12800, 0.0, zero_stack, 1.0, int_stack+17485, 1.0, int_stack+18500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12330,int_stack+11880,int_stack+52930, 0.0, zero_stack, 1.0, int_stack+46510, 1.0, int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+13235,int_stack+13175, 0.0, zero_stack, 2.0, int_stack+36575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+13335,int_stack+13235, 0.0, zero_stack, 2.0, int_stack+36635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+11880,int_stack+49210,int_stack+52930, 0.0, zero_stack, 2.0, int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+13575,int_stack+13485, 0.0, zero_stack, 2.0, int_stack+18320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12870,int_stack+13725,int_stack+13575, 0.0, zero_stack, 2.0, int_stack+18500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+13320,int_stack+12870,int_stack+52930, 0.0, zero_stack, 2.0, int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+14010,int_stack+13950, 1.0, int_stack+32235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37720,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+14110,int_stack+14010, 1.0, int_stack+32295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37880,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+12870,int_stack+49210,int_stack+52930, 1.0, int_stack+42580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48220,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+52930,int_stack+14590,int_stack+14350, 1.0, int_stack+14260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19425,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13860,int_stack+14740,int_stack+14590, 1.0, int_stack+14440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19605,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+14310,int_stack+13860,int_stack+52930, 1.0, int_stack+43420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48400,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+42580,int_stack+15025,int_stack+14965, 1.0, int_stack+33320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37720, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+15125,int_stack+15025, 1.0, int_stack+33380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37880, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+52930,int_stack+49210,int_stack+42580, 1.0, int_stack+43690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48220, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+49210,int_stack+15605,int_stack+15365, 1.0, int_stack+15275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19425, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43420,int_stack+15755,int_stack+15605, 1.0, int_stack+15455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19605, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+14850,int_stack+43420,int_stack+49210, 1.0, int_stack+43870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48400, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+42580,int_stack+16040,int_stack+15980, 1.0, int_stack+34405, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37720, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+16140,int_stack+16040, 1.0, int_stack+34465, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37880, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+43420,int_stack+49210,int_stack+42580, 1.0, int_stack+45340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48220, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+49210,int_stack+16620,int_stack+16380, 1.0, int_stack+16290, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19425, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13860,int_stack+16770,int_stack+16620, 1.0, int_stack+16470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19605, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15390,int_stack+13860,int_stack+49210, 1.0, int_stack+45520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48400, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+42580,int_stack+17055,int_stack+16995, 1.0, int_stack+35490, 0.0, zero_stack, 1.0, int_stack+37720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+17155,int_stack+17055, 1.0, int_stack+35550, 0.0, zero_stack, 1.0, int_stack+37880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+43780,int_stack+49210,int_stack+42580, 1.0, int_stack+46330, 0.0, zero_stack, 1.0, int_stack+48220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+49210,int_stack+17635,int_stack+17395, 1.0, int_stack+17305, 0.0, zero_stack, 1.0, int_stack+19425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13860,int_stack+17785,int_stack+17635, 1.0, int_stack+17485, 0.0, zero_stack, 1.0, int_stack+19605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+15930,int_stack+13860,int_stack+49210, 1.0, int_stack+46510, 0.0, zero_stack, 1.0, int_stack+48400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+42580,int_stack+18070,int_stack+18010, 1.0, int_stack+36575, 1.0, int_stack+37720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49210,int_stack+18170,int_stack+18070, 1.0, int_stack+36635, 1.0, int_stack+37880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+13860,int_stack+49210,int_stack+42580, 1.0, int_stack+47680, 1.0, int_stack+48220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+49210,int_stack+18650,int_stack+18410, 1.0, int_stack+18320, 1.0, int_stack+19425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+46330,int_stack+18800,int_stack+18650, 1.0, int_stack+18500, 1.0, int_stack+19605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+16470,int_stack+46330,int_stack+49210, 1.0, int_stack+44680, 1.0, int_stack+48400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+19085,int_stack+19025, 2.0, int_stack+37720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+19185,int_stack+19085, 2.0, int_stack+37880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+46330,int_stack+44680,int_stack+47680, 2.0, int_stack+48220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+19755,int_stack+19515, 2.0, int_stack+19425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45340,int_stack+19905,int_stack+19755, 2.0, int_stack+19605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17010,int_stack+45340,int_stack+44680, 2.0, int_stack+48400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+20190,int_stack+20130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39145,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+20290,int_stack+20190, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39205,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+45340,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+20530,int_stack+20440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24315,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+48220,int_stack+20680,int_stack+20530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24495,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+17550,int_stack+48220,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42760,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+20965,int_stack+20905, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39145, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+21065,int_stack+20965, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39205, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+48220,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+21305,int_stack+21215, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24315, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18090,int_stack+21455,int_stack+21305, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24495, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+18540,int_stack+18090,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42760, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+21740,int_stack+21680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39145, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+21840,int_stack+21740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39205, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+18090,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+22080,int_stack+21990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24315, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+19080,int_stack+22230,int_stack+22080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24495, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19530,int_stack+19080,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42760, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+22515,int_stack+22455, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+22615,int_stack+22515, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+19080,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+22855,int_stack+22765, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20070,int_stack+23005,int_stack+22855, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20520,int_stack+20070,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+23290,int_stack+23230, 0.0, zero_stack, 1.0, int_stack+39145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+23390,int_stack+23290, 0.0, zero_stack, 1.0, int_stack+39205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+20070,int_stack+44680,int_stack+47680, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+23630,int_stack+23540, 0.0, zero_stack, 1.0, int_stack+24315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21060,int_stack+23780,int_stack+23630, 0.0, zero_stack, 1.0, int_stack+24495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21510,int_stack+21060,int_stack+44680, 0.0, zero_stack, 1.0, int_stack+42760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+24065,int_stack+24005, 1.0, int_stack+39145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+24165,int_stack+24065, 1.0, int_stack+39205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+21060,int_stack+44680,int_stack+47680, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+42400,int_stack+24645,int_stack+24405, 1.0, int_stack+24315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22050,int_stack+24795,int_stack+24645, 1.0, int_stack+24495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+22500,int_stack+22050,int_stack+42400, 1.0, int_stack+42760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+25080,int_stack+25020,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+25180,int_stack+25080,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+42400,int_stack+44680,int_stack+47680,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+25420,int_stack+25330,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+22050,int_stack+25570,int_stack+25420,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+23040,int_stack+22050,int_stack+44680,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+25855,int_stack+25795, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40230,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+25955,int_stack+25855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40290,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+22050,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50590,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+26195,int_stack+26105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29980,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23580,int_stack+26345,int_stack+26195, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30160,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+24030,int_stack+23580,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50770,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+26630,int_stack+26570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40230, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+26730,int_stack+26630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40290, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+23580,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50590, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+26970,int_stack+26880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29980, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+24570,int_stack+27120,int_stack+26970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30160, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+25020,int_stack+24570,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50770, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+27405,int_stack+27345, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40230, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+27505,int_stack+27405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40290, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+24570,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50590, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+27745,int_stack+27655, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29980, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25560,int_stack+27895,int_stack+27745, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30160, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+26010,int_stack+25560,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50770, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+28180,int_stack+28120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+28280,int_stack+28180, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+25560,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+28520,int_stack+28430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26550,int_stack+28670,int_stack+28520, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+27000,int_stack+26550,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+28955,int_stack+28895, 0.0, zero_stack, 1.0, int_stack+40230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+29055,int_stack+28955, 0.0, zero_stack, 1.0, int_stack+40290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+26550,int_stack+44680,int_stack+47680, 0.0, zero_stack, 1.0, int_stack+50590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+29295,int_stack+29205, 0.0, zero_stack, 1.0, int_stack+29980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27540,int_stack+29445,int_stack+29295, 0.0, zero_stack, 1.0, int_stack+30160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+27990,int_stack+27540,int_stack+44680, 0.0, zero_stack, 1.0, int_stack+50770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+29730,int_stack+29670, 1.0, int_stack+40230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+29830,int_stack+29730, 1.0, int_stack+40290, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+27540,int_stack+44680,int_stack+47680, 1.0, int_stack+50590, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+30310,int_stack+30070, 1.0, int_stack+29980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28530,int_stack+30460,int_stack+30310, 1.0, int_stack+30160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+28980,int_stack+28530,int_stack+44680, 1.0, int_stack+50770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+30745,int_stack+30685,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+30845,int_stack+30745,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+28530,int_stack+44680,int_stack+47680,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+31085,int_stack+30995,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+50590,int_stack+31235,int_stack+31085,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+29520,int_stack+50590,int_stack+44680,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+31520,int_stack+31460,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+31620,int_stack+31520,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+50590,int_stack+44680,int_stack+47680,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+31860,int_stack+31770,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+30060,int_stack+32010,int_stack+31860,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+30510,int_stack+30060,int_stack+44680,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+32605,int_stack+32545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41315,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+32705,int_stack+32605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41375,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+30060,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51940,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+32945,int_stack+32855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38440,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+31050,int_stack+33095,int_stack+32945, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38620,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+31500,int_stack+31050,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52120,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+33690,int_stack+33630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41315, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+33790,int_stack+33690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41375, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+31050,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51940, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+34030,int_stack+33940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38440, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+32040,int_stack+34180,int_stack+34030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38620, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+32490,int_stack+32040,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52120, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+34775,int_stack+34715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41315, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+34875,int_stack+34775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41375, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+32040,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51940, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+35115,int_stack+35025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38440, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+33030,int_stack+35265,int_stack+35115, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38620, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+33480,int_stack+33030,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52120, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+35860,int_stack+35800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+35960,int_stack+35860, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+33030,int_stack+44680,int_stack+47680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+36200,int_stack+36110, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+34020,int_stack+36350,int_stack+36200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+34470,int_stack+34020,int_stack+44680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+36945,int_stack+36885, 0.0, zero_stack, 1.0, int_stack+41315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+37045,int_stack+36945, 0.0, zero_stack, 1.0, int_stack+41375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+34020,int_stack+44680,int_stack+47680, 0.0, zero_stack, 1.0, int_stack+51940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+37285,int_stack+37195, 0.0, zero_stack, 1.0, int_stack+38440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+35010,int_stack+37435,int_stack+37285, 0.0, zero_stack, 1.0, int_stack+38620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+35460,int_stack+35010,int_stack+44680, 0.0, zero_stack, 1.0, int_stack+52120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+38190,int_stack+38130, 1.0, int_stack+41315, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+38290,int_stack+38190, 1.0, int_stack+41375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+35010,int_stack+44680,int_stack+47680, 1.0, int_stack+51940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+38770,int_stack+38530, 1.0, int_stack+38440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36000,int_stack+38920,int_stack+38770, 1.0, int_stack+38620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+36450,int_stack+36000,int_stack+44680, 1.0, int_stack+52120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+39515,int_stack+39455,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+39615,int_stack+39515,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+36000,int_stack+44680,int_stack+47680,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+39855,int_stack+39765,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+51940,int_stack+40005,int_stack+39855,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+36990,int_stack+51940,int_stack+44680,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+40600,int_stack+40540,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+40700,int_stack+40600,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+51940,int_stack+44680,int_stack+47680,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+40940,int_stack+40850,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+37530,int_stack+41090,int_stack+40940,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+37980,int_stack+37530,int_stack+44680,15);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+47680,int_stack+41685,int_stack+41625,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+44680,int_stack+41785,int_stack+41685,10);
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+37530,int_stack+44680,int_stack+47680,10);
 /*--- compute (g0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+44680,int_stack+42025,int_stack+41935,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+38520,int_stack+42175,int_stack+42025,15);
 /*--- compute (g0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+38970,int_stack+38520,int_stack+44680,15);
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+39510,int_stack+44140,int_stack+43060,36);
     Libderiv->ABCD[11] = int_stack + 39510;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+40590,int_stack+45790,int_stack+44980,36);
     Libderiv->ABCD[10] = int_stack + 40590;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+56620,int_stack+46780,int_stack+0,36);
     Libderiv->ABCD[9] = int_stack + 56620;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+57700,int_stack+360,int_stack+47320,36);
     Libderiv->ABCD[8] = int_stack + 57700;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+58780,int_stack+48670,int_stack+47860,36);
     Libderiv->ABCD[7] = int_stack + 58780;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+59860,int_stack+900,int_stack+49510,36);
     Libderiv->ABCD[6] = int_stack + 59860;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+360,int_stack+51040,int_stack+50230, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[2] = int_stack + 360;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+60940,int_stack+52390,int_stack+51580, 0.0, zero_stack, 1.0, int_stack+49870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[1] = int_stack + 60940;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+62020,int_stack+53380,int_stack+1440, 1.0, int_stack+49870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[0] = int_stack + 62020;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+63100,int_stack+53920,int_stack+1800,36);
     Libderiv->ABCD[155] = int_stack + 63100;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+53290,int_stack+2970,int_stack+2160,36);
     Libderiv->ABCD[143] = int_stack + 53290;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+64180,int_stack+3960,int_stack+2520,36);
     Libderiv->ABCD[142] = int_stack + 64180;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1800,int_stack+54460,int_stack+3510,36);
     Libderiv->ABCD[131] = int_stack + 1800;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+2880,int_stack+5310,int_stack+4500,36);
     Libderiv->ABCD[130] = int_stack + 2880;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+65260,int_stack+6300,int_stack+4860,36);
     Libderiv->ABCD[129] = int_stack + 65260;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3960,int_stack+55000,int_stack+5850,36);
     Libderiv->ABCD[119] = int_stack + 3960;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5040,int_stack+7650,int_stack+6840,36);
     Libderiv->ABCD[118] = int_stack + 5040;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+6120,int_stack+8640,int_stack+7200,36);
     Libderiv->ABCD[117] = int_stack + 6120;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+54370,int_stack+55540,int_stack+8190,36);
     Libderiv->ABCD[116] = int_stack + 54370;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+7200,int_stack+9990,int_stack+9180,36);
     Libderiv->ABCD[107] = int_stack + 7200;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8280,int_stack+10980,int_stack+9540,36);
     Libderiv->ABCD[106] = int_stack + 8280;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+9360,int_stack+56080,int_stack+10530,36);
     Libderiv->ABCD[105] = int_stack + 9360;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+10440,int_stack+12330,int_stack+11520,36);
     Libderiv->ABCD[104] = int_stack + 10440;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+55450,int_stack+13320,int_stack+11880,36);
     Libderiv->ABCD[103] = int_stack + 55450;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11520,int_stack+14310,int_stack+12870,36);
     Libderiv->ABCD[95] = int_stack + 11520;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12600,int_stack+14850,int_stack+52930,36);
     Libderiv->ABCD[94] = int_stack + 12600;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+14220,int_stack+15390,int_stack+43420,36);
     Libderiv->ABCD[93] = int_stack + 14220;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+66340,int_stack+15930,int_stack+43780,36);
     Libderiv->ABCD[92] = int_stack + 66340;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+43420,int_stack+16470,int_stack+13860,36);
     Libderiv->ABCD[91] = int_stack + 43420;
 /*--- compute (fp|dd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15300,int_stack+17010,int_stack+46330,36);
     Libderiv->ABCD[90] = int_stack + 15300;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+45700,int_stack+17550,int_stack+45340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[47] = int_stack + 45700;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+16380,int_stack+18540,int_stack+48220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[46] = int_stack + 16380;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+48220,int_stack+19530,int_stack+18090, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[45] = int_stack + 48220;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+17460,int_stack+20520,int_stack+19080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[44] = int_stack + 17460;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+18540,int_stack+21510,int_stack+20070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[43] = int_stack + 18540;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+19620,int_stack+22500,int_stack+21060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[42] = int_stack + 19620;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+20700,int_stack+23040,int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+50230, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[38] = int_stack + 20700;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+22410,int_stack+24030,int_stack+22050, 0.0, zero_stack, 1.0, int_stack+43060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[35] = int_stack + 22410;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+41670,int_stack+25020,int_stack+23580, 0.0, zero_stack, 1.0, int_stack+44980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[34] = int_stack + 41670;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+23490,int_stack+26010,int_stack+24570, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[33] = int_stack + 23490;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+67420,int_stack+27000,int_stack+25560, 0.0, zero_stack, 1.0, int_stack+47320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[32] = int_stack + 67420;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+24570,int_stack+27990,int_stack+26550, 0.0, zero_stack, 1.0, int_stack+47860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[31] = int_stack + 24570;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+25650,int_stack+28980,int_stack+27540, 0.0, zero_stack, 1.0, int_stack+49510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[30] = int_stack + 25650;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+26730,int_stack+29520,int_stack+28530, 0.0, zero_stack, 1.0, int_stack+50230, 1.0, int_stack+51580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[26] = int_stack + 26730;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+27810,int_stack+30510,int_stack+50590, 0.0, zero_stack, 2.0, int_stack+51580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[25] = int_stack + 27810;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+28890,int_stack+31500,int_stack+30060, 1.0, int_stack+43060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[23] = int_stack + 28890;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+29970,int_stack+32490,int_stack+31050, 1.0, int_stack+44980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[22] = int_stack + 29970;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+44500,int_stack+33480,int_stack+32040, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[21] = int_stack + 44500;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+31050,int_stack+34470,int_stack+33030, 1.0, int_stack+47320, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[20] = int_stack + 31050;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+46780,int_stack+35460,int_stack+34020, 1.0, int_stack+47860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[19] = int_stack + 46780;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+32130,int_stack+36450,int_stack+35010, 1.0, int_stack+49510, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[18] = int_stack + 32130;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+33210,int_stack+36990,int_stack+36000, 1.0, int_stack+50230, 0.0, zero_stack, 1.0, int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[14] = int_stack + 33210;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+34290,int_stack+37980,int_stack+51940, 1.0, int_stack+51580, 1.0, int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[13] = int_stack + 34290;
 /*--- compute (fp|dd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+37890,int_stack+38970,int_stack+37530, 2.0, int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,36);
     Libderiv->ABCD[12] = int_stack + 37890;

}
