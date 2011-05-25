#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_fpfd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|fd) integrals */

void d12hrr_order_fpfd(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv2_classes[3][3][143] = int_stack + 3270;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 3370;
 Libderiv->deriv2_classes[3][5][143] = int_stack + 3520;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 3730;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 3880;
 Libderiv->deriv2_classes[4][5][143] = int_stack + 4105;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 4420;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 4520;
 Libderiv->deriv2_classes[3][5][131] = int_stack + 4670;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 4880;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 5030;
 Libderiv->deriv2_classes[4][5][131] = int_stack + 5255;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 5570;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 5670;
 Libderiv->deriv2_classes[3][5][130] = int_stack + 5820;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 6030;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 6180;
 Libderiv->deriv2_classes[4][5][130] = int_stack + 6405;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 6720;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 6820;
 Libderiv->deriv2_classes[3][5][119] = int_stack + 6970;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 7180;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 7330;
 Libderiv->deriv2_classes[4][5][119] = int_stack + 7555;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 7870;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 7970;
 Libderiv->deriv2_classes[3][5][118] = int_stack + 8120;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 8330;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 8480;
 Libderiv->deriv2_classes[4][5][118] = int_stack + 8705;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 9020;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 9120;
 Libderiv->deriv2_classes[3][5][117] = int_stack + 9270;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 9480;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 9630;
 Libderiv->deriv2_classes[4][5][117] = int_stack + 9855;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 10170;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 10270;
 Libderiv->deriv2_classes[3][5][107] = int_stack + 10420;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 10630;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 10780;
 Libderiv->deriv2_classes[4][5][107] = int_stack + 11005;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 11320;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 11420;
 Libderiv->deriv2_classes[3][5][106] = int_stack + 11570;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 11780;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 11930;
 Libderiv->deriv2_classes[4][5][106] = int_stack + 12155;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 12470;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 12570;
 Libderiv->deriv2_classes[3][5][105] = int_stack + 12720;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 12930;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 13080;
 Libderiv->deriv2_classes[4][5][105] = int_stack + 13305;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 13620;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 13720;
 Libderiv->deriv2_classes[3][5][104] = int_stack + 13870;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 14080;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 14230;
 Libderiv->deriv2_classes[4][5][104] = int_stack + 14455;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 14770;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 14870;
 Libderiv->deriv2_classes[3][5][95] = int_stack + 15020;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 15230;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 15380;
 Libderiv->deriv2_classes[4][5][95] = int_stack + 15605;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 15920;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 16020;
 Libderiv->deriv2_classes[3][5][94] = int_stack + 16170;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 16380;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 16530;
 Libderiv->deriv2_classes[4][5][94] = int_stack + 16755;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 17070;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 17170;
 Libderiv->deriv2_classes[3][5][93] = int_stack + 17320;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 17530;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 17680;
 Libderiv->deriv2_classes[4][5][93] = int_stack + 17905;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 18220;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 18320;
 Libderiv->deriv2_classes[3][5][92] = int_stack + 18470;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 18680;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 18830;
 Libderiv->deriv2_classes[4][5][92] = int_stack + 19055;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 19370;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 19470;
 Libderiv->deriv2_classes[3][5][91] = int_stack + 19620;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 19830;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 19980;
 Libderiv->deriv2_classes[4][5][91] = int_stack + 20205;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 20520;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 20620;
 Libderiv->deriv2_classes[3][5][83] = int_stack + 20770;
 Libderiv->deriv_classes[4][3][11] = int_stack + 20980;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 21130;
 Libderiv->deriv_classes[4][4][11] = int_stack + 21280;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 21505;
 Libderiv->deriv2_classes[4][5][83] = int_stack + 21730;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 22045;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 22145;
 Libderiv->deriv2_classes[3][5][82] = int_stack + 22295;
 Libderiv->deriv_classes[4][3][10] = int_stack + 22505;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 22655;
 Libderiv->deriv_classes[4][4][10] = int_stack + 22805;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 23030;
 Libderiv->deriv2_classes[4][5][82] = int_stack + 23255;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 23570;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 23670;
 Libderiv->deriv2_classes[3][5][81] = int_stack + 23820;
 Libderiv->deriv_classes[4][3][9] = int_stack + 24030;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 24180;
 Libderiv->deriv_classes[4][4][9] = int_stack + 24330;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 24555;
 Libderiv->deriv2_classes[4][5][81] = int_stack + 24780;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 25095;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 25195;
 Libderiv->deriv2_classes[3][5][80] = int_stack + 25345;
 Libderiv->deriv_classes[4][3][8] = int_stack + 25555;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 25705;
 Libderiv->deriv_classes[4][4][8] = int_stack + 25855;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 26080;
 Libderiv->deriv2_classes[4][5][80] = int_stack + 26305;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 26620;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 26720;
 Libderiv->deriv2_classes[3][5][79] = int_stack + 26870;
 Libderiv->deriv_classes[4][3][7] = int_stack + 27080;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 27230;
 Libderiv->deriv_classes[4][4][7] = int_stack + 27380;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 27605;
 Libderiv->deriv2_classes[4][5][79] = int_stack + 27830;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 28145;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 28245;
 Libderiv->deriv2_classes[3][5][78] = int_stack + 28395;
 Libderiv->dvrr_classes[4][3] = int_stack + 28605;
 Libderiv->deriv_classes[4][3][6] = int_stack + 28755;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 28905;
 Libderiv->deriv_classes[4][4][6] = int_stack + 29055;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 29280;
 Libderiv->deriv2_classes[4][5][78] = int_stack + 29505;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 29820;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 29920;
 Libderiv->deriv2_classes[3][5][35] = int_stack + 30070;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 30280;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 30430;
 Libderiv->deriv2_classes[4][5][35] = int_stack + 30655;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 30970;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 31070;
 Libderiv->deriv2_classes[3][5][34] = int_stack + 31220;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 31430;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 31580;
 Libderiv->deriv2_classes[4][5][34] = int_stack + 31805;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 32120;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 32220;
 Libderiv->deriv2_classes[3][5][33] = int_stack + 32370;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 32580;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 32730;
 Libderiv->deriv2_classes[4][5][33] = int_stack + 32955;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 33270;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 33370;
 Libderiv->deriv2_classes[3][5][32] = int_stack + 33520;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 33730;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 33880;
 Libderiv->deriv2_classes[4][5][32] = int_stack + 34105;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 34420;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 34520;
 Libderiv->deriv2_classes[3][5][31] = int_stack + 34670;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 34880;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 35030;
 Libderiv->deriv2_classes[4][5][31] = int_stack + 35255;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 35570;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 35670;
 Libderiv->deriv2_classes[3][5][30] = int_stack + 35820;
 Libderiv->deriv_classes[4][3][2] = int_stack + 36030;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 36180;
 Libderiv->deriv_classes[4][4][2] = int_stack + 36330;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 36555;
 Libderiv->deriv2_classes[4][5][30] = int_stack + 36780;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 37095;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 37195;
 Libderiv->deriv2_classes[3][5][26] = int_stack + 37345;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 37555;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 37705;
 Libderiv->deriv2_classes[4][5][26] = int_stack + 37930;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 38245;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 38345;
 Libderiv->deriv2_classes[3][5][23] = int_stack + 38495;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 38705;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 38855;
 Libderiv->deriv2_classes[4][5][23] = int_stack + 39080;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 39395;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 39495;
 Libderiv->deriv2_classes[3][5][22] = int_stack + 39645;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 39855;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 40005;
 Libderiv->deriv2_classes[4][5][22] = int_stack + 40230;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 40545;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 40645;
 Libderiv->deriv2_classes[3][5][21] = int_stack + 40795;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 41005;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 41155;
 Libderiv->deriv2_classes[4][5][21] = int_stack + 41380;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 41695;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 41795;
 Libderiv->deriv2_classes[3][5][20] = int_stack + 41945;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 42155;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 42305;
 Libderiv->deriv2_classes[4][5][20] = int_stack + 42530;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 42845;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 42945;
 Libderiv->deriv2_classes[3][5][19] = int_stack + 43095;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 43305;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 43455;
 Libderiv->deriv2_classes[4][5][19] = int_stack + 43680;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 43995;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 44095;
 Libderiv->deriv2_classes[3][5][18] = int_stack + 44245;
 Libderiv->deriv_classes[4][3][1] = int_stack + 44455;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 44605;
 Libderiv->deriv_classes[4][4][1] = int_stack + 44755;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 44980;
 Libderiv->deriv2_classes[4][5][18] = int_stack + 45205;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 45520;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 45620;
 Libderiv->deriv2_classes[3][5][14] = int_stack + 45770;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 45980;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 46130;
 Libderiv->deriv2_classes[4][5][14] = int_stack + 46355;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 46670;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 46770;
 Libderiv->deriv2_classes[3][5][13] = int_stack + 46920;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 47130;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 47280;
 Libderiv->deriv2_classes[4][5][13] = int_stack + 47505;
 Libderiv->deriv_classes[3][3][11] = int_stack + 47820;
 Libderiv->deriv_classes[3][4][11] = int_stack + 47920;
 Libderiv->deriv_classes[3][5][11] = int_stack + 48070;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 48280;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 48380;
 Libderiv->deriv2_classes[3][5][11] = int_stack + 48530;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 48740;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 48890;
 Libderiv->deriv2_classes[4][5][11] = int_stack + 49115;
 Libderiv->deriv_classes[3][3][10] = int_stack + 49430;
 Libderiv->deriv_classes[3][4][10] = int_stack + 49530;
 Libderiv->deriv_classes[3][5][10] = int_stack + 49680;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 49890;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 49990;
 Libderiv->deriv2_classes[3][5][10] = int_stack + 50140;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 50350;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 50500;
 Libderiv->deriv2_classes[4][5][10] = int_stack + 50725;
 Libderiv->deriv_classes[3][3][9] = int_stack + 51040;
 Libderiv->deriv_classes[3][4][9] = int_stack + 51140;
 Libderiv->deriv_classes[3][5][9] = int_stack + 51290;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 51500;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 51600;
 Libderiv->deriv2_classes[3][5][9] = int_stack + 51750;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 51960;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 52110;
 Libderiv->deriv2_classes[4][5][9] = int_stack + 52335;
 Libderiv->deriv_classes[3][3][8] = int_stack + 52650;
 Libderiv->deriv_classes[3][4][8] = int_stack + 52750;
 Libderiv->deriv_classes[3][5][8] = int_stack + 52900;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 53110;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 53210;
 Libderiv->deriv2_classes[3][5][8] = int_stack + 53360;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 53570;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 53720;
 Libderiv->deriv2_classes[4][5][8] = int_stack + 53945;
 Libderiv->deriv_classes[3][3][7] = int_stack + 54260;
 Libderiv->deriv_classes[3][4][7] = int_stack + 54360;
 Libderiv->deriv_classes[3][5][7] = int_stack + 54510;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 54720;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 54820;
 Libderiv->deriv2_classes[3][5][7] = int_stack + 54970;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 55180;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 55330;
 Libderiv->deriv2_classes[4][5][7] = int_stack + 55555;
 Libderiv->dvrr_classes[3][3] = int_stack + 55870;
 Libderiv->deriv_classes[3][3][6] = int_stack + 55970;
 Libderiv->dvrr_classes[3][4] = int_stack + 56070;
 Libderiv->deriv_classes[3][4][6] = int_stack + 56220;
 Libderiv->deriv_classes[3][5][6] = int_stack + 56370;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 56580;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 56680;
 Libderiv->deriv2_classes[3][5][6] = int_stack + 56830;
 Libderiv->deriv_classes[4][3][0] = int_stack + 57040;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 57190;
 Libderiv->deriv_classes[4][4][0] = int_stack + 57340;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 57565;
 Libderiv->deriv2_classes[4][5][6] = int_stack + 57790;
 Libderiv->deriv_classes[3][3][2] = int_stack + 58105;
 Libderiv->deriv_classes[3][4][2] = int_stack + 58205;
 Libderiv->deriv_classes[3][5][2] = int_stack + 58355;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 58565;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 58665;
 Libderiv->deriv2_classes[3][5][2] = int_stack + 58815;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 59025;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 59175;
 Libderiv->deriv2_classes[4][5][2] = int_stack + 59400;
 Libderiv->deriv_classes[3][3][1] = int_stack + 59715;
 Libderiv->deriv_classes[3][4][1] = int_stack + 59815;
 Libderiv->deriv_classes[3][5][1] = int_stack + 59965;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 60175;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 60275;
 Libderiv->deriv2_classes[3][5][1] = int_stack + 60425;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 60635;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 60785;
 Libderiv->deriv2_classes[4][5][1] = int_stack + 61010;
 Libderiv->deriv_classes[3][3][0] = int_stack + 61325;
 Libderiv->deriv_classes[3][4][0] = int_stack + 61425;
 Libderiv->deriv_classes[3][5][0] = int_stack + 61575;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 61785;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 61885;
 Libderiv->deriv2_classes[3][5][0] = int_stack + 62035;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 62245;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 62395;
 Libderiv->deriv2_classes[4][5][0] = int_stack + 62620;
 memset(int_stack,0,503480);

 Libderiv->dvrr_stack = int_stack + 117310;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_fpfd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+62935,int_stack+56070,int_stack+55870,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+47920,int_stack+47820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55870,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+63535,int_stack+48070,int_stack+47920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56070,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+63985,int_stack+63535,int_stack+63235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62935,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+63535,int_stack+1575,int_stack+28605,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+64585,int_stack+21280,int_stack+20980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28605,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+65035,int_stack+0,int_stack+21280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1575,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+65710,int_stack+65035,int_stack+64585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63535,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+65035,int_stack+49530,int_stack+49430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55870, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+49680,int_stack+49530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56070, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+67060,int_stack+66610,int_stack+65035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62935, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+66610,int_stack+22805,int_stack+22505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28605, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+67660,int_stack+315,int_stack+22805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1575, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+68335,int_stack+67660,int_stack+66610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63535, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+67660,int_stack+51140,int_stack+51040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55870, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+51290,int_stack+51140, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56070, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+69235,int_stack+0,int_stack+67660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62935, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+24330,int_stack+24030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28605, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+69835,int_stack+630,int_stack+24330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+70510,int_stack+69835,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63535, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+69835,int_stack+52750,int_stack+52650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+52900,int_stack+52750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+71410,int_stack+450,int_stack+69835, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+25855,int_stack+25555, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+72010,int_stack+945,int_stack+25855, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+72685,int_stack+72010,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+72010,int_stack+54360,int_stack+54260, 0.0, zero_stack, 1.0, int_stack+55870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+54510,int_stack+54360, 0.0, zero_stack, 1.0, int_stack+56070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+74035,int_stack+73585,int_stack+72010, 0.0, zero_stack, 1.0, int_stack+62935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+27380,int_stack+27080, 0.0, zero_stack, 1.0, int_stack+28605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+74635,int_stack+1260,int_stack+27380, 0.0, zero_stack, 1.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+75310,int_stack+74635,int_stack+73585, 0.0, zero_stack, 1.0, int_stack+63535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+74635,int_stack+56220,int_stack+55970, 1.0, int_stack+55870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+56370,int_stack+56220, 1.0, int_stack+56070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+76210,int_stack+900,int_stack+74635, 1.0, int_stack+62935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+29055,int_stack+28755, 1.0, int_stack+28605, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+76810,int_stack+1800,int_stack+29055, 1.0, int_stack+1575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+77485,int_stack+76810,int_stack+900, 1.0, int_stack+63535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+63535,int_stack+2745,int_stack+56070,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+76810,int_stack+63535,int_stack+62935,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+62935,int_stack+58205,int_stack+58105,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+63535,int_stack+58355,int_stack+58205,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1350,int_stack+63535,int_stack+62935,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+63535,int_stack+36330,int_stack+36030,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+78385,int_stack+2115,int_stack+36330,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+79060,int_stack+78385,int_stack+63535,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+78385,int_stack+59815,int_stack+59715,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1950,int_stack+59965,int_stack+59815,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+79960,int_stack+1950,int_stack+78385,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+44755,int_stack+44455,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+80560,int_stack+2430,int_stack+44755,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+81235,int_stack+80560,int_stack+1950,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+80560,int_stack+61425,int_stack+61325,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2400,int_stack+61575,int_stack+61425,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+82135,int_stack+2400,int_stack+80560,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2400,int_stack+57340,int_stack+57040,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+2955,int_stack+57340,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+83410,int_stack+82735,int_stack+2400,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82735,int_stack+3370,int_stack+3270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+47820,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2850,int_stack+3520,int_stack+3370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+47920,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+84310,int_stack+2850,int_stack+82735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+63235,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82735,int_stack+3880,int_stack+3730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+20980,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2850,int_stack+4105,int_stack+3880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+21280,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+84910,int_stack+2850,int_stack+82735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+64585,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+82735,int_stack+4520,int_stack+4420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47820, 1.0, int_stack+49430,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2850,int_stack+4670,int_stack+4520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47920, 1.0, int_stack+49530,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3300,int_stack+2850,int_stack+82735, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63235, 1.0, int_stack+65035,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+5030,int_stack+4880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20980, 1.0, int_stack+22505,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+5255,int_stack+5030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21280, 1.0, int_stack+22805,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3900,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64585, 1.0, int_stack+66610,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+5670,int_stack+5570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49430, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+5820,int_stack+5670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+49530, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4800,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+65035, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+6180,int_stack+6030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22505, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+6405,int_stack+6180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+22805, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5400,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+66610, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+6820,int_stack+6720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47820, 0.0, zero_stack, 1.0, int_stack+51040,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+6970,int_stack+6820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47920, 0.0, zero_stack, 1.0, int_stack+51140,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6300,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63235, 0.0, zero_stack, 1.0, int_stack+67660,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+7330,int_stack+7180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20980, 0.0, zero_stack, 1.0, int_stack+24030,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+7555,int_stack+7330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21280, 0.0, zero_stack, 1.0, int_stack+24330,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6900,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64585, 0.0, zero_stack, 1.0, int_stack+0,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+7970,int_stack+7870, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49430, 1.0, int_stack+51040, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+8120,int_stack+7970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49530, 1.0, int_stack+51140, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+85810,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65035, 1.0, int_stack+67660, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+8480,int_stack+8330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22505, 1.0, int_stack+24030, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+8705,int_stack+8480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22805, 1.0, int_stack+24330, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7800,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66610, 1.0, int_stack+0, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+9120,int_stack+9020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+51040, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+9270,int_stack+9120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+51140, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8700,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+67660, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+9630,int_stack+9480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24030, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+9855,int_stack+9630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+24330, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+86410,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+10270,int_stack+10170, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47820, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52650,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+10420,int_stack+10270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47920, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52750,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9300,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69835,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+10780,int_stack+10630, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+20980, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25555,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+11005,int_stack+10780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+21280, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+25855,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9900,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+64585, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+450,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+11420,int_stack+11320, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49430, 0.0, zero_stack, 1.0, int_stack+52650, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+11570,int_stack+11420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49530, 0.0, zero_stack, 1.0, int_stack+52750, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10800,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+65035, 0.0, zero_stack, 1.0, int_stack+69835, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+11930,int_stack+11780, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22505, 0.0, zero_stack, 1.0, int_stack+25555, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+12155,int_stack+11930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+22805, 0.0, zero_stack, 1.0, int_stack+25855, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+11400,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+66610, 0.0, zero_stack, 1.0, int_stack+450, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+12570,int_stack+12470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51040, 1.0, int_stack+52650, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+12720,int_stack+12570, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51140, 1.0, int_stack+52750, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12300,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67660, 1.0, int_stack+69835, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+13080,int_stack+12930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24030, 1.0, int_stack+25555, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+13305,int_stack+13080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24330, 1.0, int_stack+25855, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+87310,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+13720,int_stack+13620, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+52650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+13870,int_stack+13720, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+52750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12900,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+69835, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+14230,int_stack+14080, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25555, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+14455,int_stack+14230, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+25855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13500,int_stack+82735,int_stack+2850, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+14870,int_stack+14770, 0.0, zero_stack, 1.0, int_stack+47820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54260,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+15020,int_stack+14870, 0.0, zero_stack, 1.0, int_stack+47920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54360,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14400,int_stack+82735,int_stack+2850, 0.0, zero_stack, 1.0, int_stack+63235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72010,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+15380,int_stack+15230, 0.0, zero_stack, 1.0, int_stack+20980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27080,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+15605,int_stack+15380, 0.0, zero_stack, 1.0, int_stack+21280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27380,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15000,int_stack+82735,int_stack+2850, 0.0, zero_stack, 1.0, int_stack+64585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73585,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+16020,int_stack+15920, 0.0, zero_stack, 1.0, int_stack+49430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54260, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+16170,int_stack+16020, 0.0, zero_stack, 1.0, int_stack+49530, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54360, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+88210,int_stack+82735,int_stack+2850, 0.0, zero_stack, 1.0, int_stack+65035, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+72010, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+16530,int_stack+16380, 0.0, zero_stack, 1.0, int_stack+22505, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27080, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+16755,int_stack+16530, 0.0, zero_stack, 1.0, int_stack+22805, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27380, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15900,int_stack+82735,int_stack+2850, 0.0, zero_stack, 1.0, int_stack+66610, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+73585, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+17170,int_stack+17070, 0.0, zero_stack, 1.0, int_stack+51040, 0.0, zero_stack, 1.0, int_stack+54260, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+17320,int_stack+17170, 0.0, zero_stack, 1.0, int_stack+51140, 0.0, zero_stack, 1.0, int_stack+54360, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16800,int_stack+82735,int_stack+2850, 0.0, zero_stack, 1.0, int_stack+67660, 0.0, zero_stack, 1.0, int_stack+72010, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+17680,int_stack+17530, 0.0, zero_stack, 1.0, int_stack+24030, 0.0, zero_stack, 1.0, int_stack+27080, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+17905,int_stack+17680, 0.0, zero_stack, 1.0, int_stack+24330, 0.0, zero_stack, 1.0, int_stack+27380, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+88810,int_stack+82735,int_stack+2850, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+73585, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+18320,int_stack+18220, 0.0, zero_stack, 1.0, int_stack+52650, 1.0, int_stack+54260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+18470,int_stack+18320, 0.0, zero_stack, 1.0, int_stack+52750, 1.0, int_stack+54360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+17400,int_stack+82735,int_stack+2850, 0.0, zero_stack, 1.0, int_stack+69835, 1.0, int_stack+72010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+18830,int_stack+18680, 0.0, zero_stack, 1.0, int_stack+25555, 1.0, int_stack+27080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+19055,int_stack+18830, 0.0, zero_stack, 1.0, int_stack+25855, 1.0, int_stack+27380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18000,int_stack+82735,int_stack+2850, 0.0, zero_stack, 1.0, int_stack+450, 1.0, int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+19470,int_stack+19370, 0.0, zero_stack, 2.0, int_stack+54260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+19620,int_stack+19470, 0.0, zero_stack, 2.0, int_stack+54360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18900,int_stack+82735,int_stack+2850, 0.0, zero_stack, 2.0, int_stack+72010, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+19980,int_stack+19830, 0.0, zero_stack, 2.0, int_stack+27080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+20205,int_stack+19980, 0.0, zero_stack, 2.0, int_stack+27380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19500,int_stack+82735,int_stack+2850, 0.0, zero_stack, 2.0, int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+20620,int_stack+20520, 1.0, int_stack+47820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55970,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+20770,int_stack+20620, 1.0, int_stack+47920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56220,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+89710,int_stack+82735,int_stack+2850, 1.0, int_stack+63235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74635,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+21505,int_stack+21130, 1.0, int_stack+20980, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28755,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+82735,int_stack+21730,int_stack+21505, 1.0, int_stack+21280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29055,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20400,int_stack+82735,int_stack+2850, 1.0, int_stack+64585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+22145,int_stack+22045, 1.0, int_stack+49430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55970, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+64585,int_stack+22295,int_stack+22145, 1.0, int_stack+49530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56220, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+82735,int_stack+64585,int_stack+63235, 1.0, int_stack+65035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74635, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+23030,int_stack+22655, 1.0, int_stack+22505, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28755, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+64585,int_stack+23255,int_stack+23030, 1.0, int_stack+22805, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29055, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21300,int_stack+64585,int_stack+2850, 1.0, int_stack+66610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+23670,int_stack+23570, 1.0, int_stack+51040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+55970, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+23820,int_stack+23670, 1.0, int_stack+51140, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+56220, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+64585,int_stack+66610,int_stack+63235, 1.0, int_stack+67660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74635, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+66610,int_stack+24555,int_stack+24180, 1.0, int_stack+24030, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28755, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+67660,int_stack+24780,int_stack+24555, 1.0, int_stack+24330, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29055, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22200,int_stack+67660,int_stack+66610, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+25195,int_stack+25095, 1.0, int_stack+52650, 0.0, zero_stack, 1.0, int_stack+55970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+25345,int_stack+25195, 1.0, int_stack+52750, 0.0, zero_stack, 1.0, int_stack+56220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+67660,int_stack+0,int_stack+63235, 1.0, int_stack+69835, 0.0, zero_stack, 1.0, int_stack+74635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+26080,int_stack+25705, 1.0, int_stack+25555, 0.0, zero_stack, 1.0, int_stack+28755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+69835,int_stack+26305,int_stack+26080, 1.0, int_stack+25855, 0.0, zero_stack, 1.0, int_stack+29055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23100,int_stack+69835,int_stack+0, 1.0, int_stack+450, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+26720,int_stack+26620, 1.0, int_stack+54260, 1.0, int_stack+55970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+26870,int_stack+26720, 1.0, int_stack+54360, 1.0, int_stack+56220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+66610,int_stack+63235, 1.0, int_stack+72010, 1.0, int_stack+74635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+66610,int_stack+27605,int_stack+27230, 1.0, int_stack+27080, 1.0, int_stack+28755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+72010,int_stack+27830,int_stack+27605, 1.0, int_stack+27380, 1.0, int_stack+29055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24000,int_stack+72010,int_stack+66610, 1.0, int_stack+73585, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+28245,int_stack+28145, 2.0, int_stack+55970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+28395,int_stack+28245, 2.0, int_stack+56220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+72010,int_stack+73585,int_stack+63235, 2.0, int_stack+74635, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+29280,int_stack+28905, 2.0, int_stack+28755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+74635,int_stack+29505,int_stack+29280, 2.0, int_stack+29055, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+24900,int_stack+74635,int_stack+73585, 2.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+29920,int_stack+29820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58105,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+30070,int_stack+29920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58205,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+74635,int_stack+73585,int_stack+63235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62935,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+30430,int_stack+30280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36030,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+69835,int_stack+30655,int_stack+30430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36330,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+25800,int_stack+69835,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63535,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+31070,int_stack+30970, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58105, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+31220,int_stack+31070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58205, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+69835,int_stack+73585,int_stack+63235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62935, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+31580,int_stack+31430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36030, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+600,int_stack+31805,int_stack+31580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36330, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+26700,int_stack+600,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63535, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+32220,int_stack+32120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58105, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+32370,int_stack+32220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58205, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+600,int_stack+73585,int_stack+63235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62935, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+32730,int_stack+32580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36030, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+27600,int_stack+32955,int_stack+32730, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36330, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+28275,int_stack+27600,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63535, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+33370,int_stack+33270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+33520,int_stack+33370, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+58205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+27600,int_stack+73585,int_stack+63235, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+62935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+33880,int_stack+33730, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+29175,int_stack+34105,int_stack+33880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29850,int_stack+29175,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+34520,int_stack+34420, 0.0, zero_stack, 1.0, int_stack+58105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+34670,int_stack+34520, 0.0, zero_stack, 1.0, int_stack+58205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+29175,int_stack+73585,int_stack+63235, 0.0, zero_stack, 1.0, int_stack+62935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+35030,int_stack+34880, 0.0, zero_stack, 1.0, int_stack+36030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+30750,int_stack+35255,int_stack+35030, 0.0, zero_stack, 1.0, int_stack+36330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+31425,int_stack+30750,int_stack+73585, 0.0, zero_stack, 1.0, int_stack+63535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+63235,int_stack+35670,int_stack+35570, 1.0, int_stack+58105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+35820,int_stack+35670, 1.0, int_stack+58205, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+30750,int_stack+73585,int_stack+63235, 1.0, int_stack+62935, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+36555,int_stack+36180, 1.0, int_stack+36030, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+32325,int_stack+36780,int_stack+36555, 1.0, int_stack+36330, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+33000,int_stack+32325,int_stack+73585, 1.0, int_stack+63535, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+37195,int_stack+37095,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+37345,int_stack+37195,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+32325,int_stack+66610,int_stack+73585,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+37705,int_stack+37555,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+62935,int_stack+37930,int_stack+37705,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+33900,int_stack+62935,int_stack+73585,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+38345,int_stack+38245, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59715,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+38495,int_stack+38345, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59815,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+62935,int_stack+66610,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78385,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+38855,int_stack+38705, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44455,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+34800,int_stack+39080,int_stack+38855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44755,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+35475,int_stack+34800,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1950,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+39495,int_stack+39395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59715, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+39645,int_stack+39495, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59815, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+34800,int_stack+66610,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78385, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+40005,int_stack+39855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44455, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+36375,int_stack+40230,int_stack+40005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44755, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+37050,int_stack+36375,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1950, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+40645,int_stack+40545, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59715, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+40795,int_stack+40645, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59815, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+36375,int_stack+66610,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78385, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+41155,int_stack+41005, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44455, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+37950,int_stack+41380,int_stack+41155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44755, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+38625,int_stack+37950,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1950, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+41795,int_stack+41695, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+41945,int_stack+41795, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+59815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+37950,int_stack+66610,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+78385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+42305,int_stack+42155, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+39525,int_stack+42530,int_stack+42305, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+40200,int_stack+39525,int_stack+73585, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+42945,int_stack+42845, 0.0, zero_stack, 1.0, int_stack+59715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+43095,int_stack+42945, 0.0, zero_stack, 1.0, int_stack+59815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+39525,int_stack+66610,int_stack+73585, 0.0, zero_stack, 1.0, int_stack+78385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+43455,int_stack+43305, 0.0, zero_stack, 1.0, int_stack+44455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+41100,int_stack+43680,int_stack+43455, 0.0, zero_stack, 1.0, int_stack+44755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41775,int_stack+41100,int_stack+73585, 0.0, zero_stack, 1.0, int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+44095,int_stack+43995, 1.0, int_stack+59715, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+66610,int_stack+44245,int_stack+44095, 1.0, int_stack+59815, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+41100,int_stack+66610,int_stack+73585, 1.0, int_stack+78385, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+44980,int_stack+44605, 1.0, int_stack+44455, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+78385,int_stack+45205,int_stack+44980, 1.0, int_stack+44755, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+42675,int_stack+78385,int_stack+73585, 1.0, int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+45620,int_stack+45520,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+45770,int_stack+45620,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+78385,int_stack+73585,int_stack+1950,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+46130,int_stack+45980,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+43575,int_stack+46355,int_stack+46130,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+44250,int_stack+43575,int_stack+1950,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+46770,int_stack+46670,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+46920,int_stack+46770,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+43575,int_stack+73585,int_stack+1950,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+47280,int_stack+47130,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+45150,int_stack+47505,int_stack+47280,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+45825,int_stack+45150,int_stack+1950,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+48380,int_stack+48280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61325,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+48530,int_stack+48380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61425,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+45150,int_stack+73585,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80560,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+48890,int_stack+48740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57040,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+46725,int_stack+49115,int_stack+48890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57340,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+47400,int_stack+46725,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2400,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+49990,int_stack+49890, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61325, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+50140,int_stack+49990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61425, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+46725,int_stack+73585,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80560, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+50500,int_stack+50350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57040, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+48300,int_stack+50725,int_stack+50500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57340, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+48975,int_stack+48300,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2400, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+51600,int_stack+51500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61325, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+51750,int_stack+51600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61425, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+48300,int_stack+73585,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80560, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+52110,int_stack+51960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57040, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+49875,int_stack+52335,int_stack+52110, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57340, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+50550,int_stack+49875,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2400, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+53210,int_stack+53110, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+53360,int_stack+53210, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+61425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+49875,int_stack+73585,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+80560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+53720,int_stack+53570, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+51450,int_stack+53945,int_stack+53720, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+57340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+52125,int_stack+51450,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+54820,int_stack+54720, 0.0, zero_stack, 1.0, int_stack+61325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+54970,int_stack+54820, 0.0, zero_stack, 1.0, int_stack+61425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+51450,int_stack+73585,int_stack+1950, 0.0, zero_stack, 1.0, int_stack+80560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+55330,int_stack+55180, 0.0, zero_stack, 1.0, int_stack+57040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+53025,int_stack+55555,int_stack+55330, 0.0, zero_stack, 1.0, int_stack+57340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+53700,int_stack+53025,int_stack+1950, 0.0, zero_stack, 1.0, int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+56680,int_stack+56580, 1.0, int_stack+61325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+56830,int_stack+56680, 1.0, int_stack+61425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+53025,int_stack+73585,int_stack+1950, 1.0, int_stack+80560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+57565,int_stack+57190, 1.0, int_stack+57040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+80560,int_stack+57790,int_stack+57565, 1.0, int_stack+57340, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+54600,int_stack+80560,int_stack+1950, 1.0, int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+58665,int_stack+58565,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+58815,int_stack+58665,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2250,int_stack+73585,int_stack+1950,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+59175,int_stack+59025,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+80560,int_stack+59400,int_stack+59175,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+55500,int_stack+80560,int_stack+73585,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+60275,int_stack+60175,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+60425,int_stack+60275,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+80560,int_stack+73585,int_stack+1950,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+60785,int_stack+60635,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+56400,int_stack+61010,int_stack+60785,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+57075,int_stack+56400,int_stack+73585,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1950,int_stack+61885,int_stack+61785,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+73585,int_stack+62035,int_stack+61885,10);
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+56400,int_stack+73585,int_stack+1950,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+73585,int_stack+62395,int_stack+62245,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+57975,int_stack+62620,int_stack+62395,15);
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+58650,int_stack+57975,int_stack+73585,15);
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+59550,int_stack+65710,int_stack+63985,60);
     Libderiv->ABCD[11] = int_stack + 59550;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+65185,int_stack+68335,int_stack+67060,60);
     Libderiv->ABCD[10] = int_stack + 65185;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+90310,int_stack+70510,int_stack+69235,60);
     Libderiv->ABCD[9] = int_stack + 90310;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+92110,int_stack+72685,int_stack+71410,60);
     Libderiv->ABCD[8] = int_stack + 92110;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+93910,int_stack+75310,int_stack+74035,60);
     Libderiv->ABCD[7] = int_stack + 93910;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+95710,int_stack+77485,int_stack+76210,60);
     Libderiv->ABCD[6] = int_stack + 95710;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+97510,int_stack+79060,int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[2] = int_stack + 97510;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+99310,int_stack+81235,int_stack+79960, 0.0, zero_stack, 1.0, int_stack+76810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[1] = int_stack + 99310;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+101110,int_stack+83410,int_stack+82135, 1.0, int_stack+76810, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[0] = int_stack + 101110;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+102910,int_stack+84910,int_stack+84310,60);
     Libderiv->ABCD[155] = int_stack + 102910;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+83335,int_stack+3900,int_stack+3300,60);
     Libderiv->ABCD[143] = int_stack + 83335;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+2850,int_stack+5400,int_stack+4800,60);
     Libderiv->ABCD[142] = int_stack + 2850;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+104710,int_stack+6900,int_stack+6300,60);
     Libderiv->ABCD[131] = int_stack + 104710;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+4650,int_stack+7800,int_stack+85810,60);
     Libderiv->ABCD[130] = int_stack + 4650;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+6450,int_stack+86410,int_stack+8700,60);
     Libderiv->ABCD[129] = int_stack + 6450;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+85135,int_stack+9900,int_stack+9300,60);
     Libderiv->ABCD[119] = int_stack + 85135;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8250,int_stack+11400,int_stack+10800,60);
     Libderiv->ABCD[118] = int_stack + 8250;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+10050,int_stack+87310,int_stack+12300,60);
     Libderiv->ABCD[117] = int_stack + 10050;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+106510,int_stack+13500,int_stack+12900,60);
     Libderiv->ABCD[116] = int_stack + 106510;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11850,int_stack+15000,int_stack+14400,60);
     Libderiv->ABCD[107] = int_stack + 11850;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+13650,int_stack+15900,int_stack+88210,60);
     Libderiv->ABCD[106] = int_stack + 13650;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+86935,int_stack+88810,int_stack+16800,60);
     Libderiv->ABCD[105] = int_stack + 86935;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15450,int_stack+18000,int_stack+17400,60);
     Libderiv->ABCD[104] = int_stack + 15450;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+108310,int_stack+19500,int_stack+18900,60);
     Libderiv->ABCD[103] = int_stack + 108310;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+17250,int_stack+20400,int_stack+89710,60);
     Libderiv->ABCD[95] = int_stack + 17250;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+19050,int_stack+21300,int_stack+82735,60);
     Libderiv->ABCD[94] = int_stack + 19050;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+110110,int_stack+22200,int_stack+64585,60);
     Libderiv->ABCD[93] = int_stack + 110110;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+20850,int_stack+23100,int_stack+67660,60);
     Libderiv->ABCD[92] = int_stack + 20850;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+111910,int_stack+24000,int_stack+0,60);
     Libderiv->ABCD[91] = int_stack + 111910;
 /*--- compute (fp|fd) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+22650,int_stack+24900,int_stack+72010,60);
     Libderiv->ABCD[90] = int_stack + 22650;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+72010,int_stack+25800,int_stack+74635, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+63985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[47] = int_stack + 72010;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+24450,int_stack+26700,int_stack+69835, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+67060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[46] = int_stack + 24450;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+113710,int_stack+28275,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+69235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[45] = int_stack + 113710;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+115510,int_stack+29850,int_stack+27600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+71410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[44] = int_stack + 115510;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+26250,int_stack+31425,int_stack+29175, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+74035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[43] = int_stack + 26250;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+28050,int_stack+33000,int_stack+30750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+76210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[42] = int_stack + 28050;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+29850,int_stack+33900,int_stack+32325, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[38] = int_stack + 29850;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+31650,int_stack+35475,int_stack+62935, 0.0, zero_stack, 1.0, int_stack+63985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[35] = int_stack + 31650;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+61350,int_stack+37050,int_stack+34800, 0.0, zero_stack, 1.0, int_stack+67060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[34] = int_stack + 61350;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+33450,int_stack+38625,int_stack+36375, 0.0, zero_stack, 1.0, int_stack+69235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[33] = int_stack + 33450;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+35250,int_stack+40200,int_stack+37950, 0.0, zero_stack, 1.0, int_stack+71410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[32] = int_stack + 35250;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+37050,int_stack+41775,int_stack+39525, 0.0, zero_stack, 1.0, int_stack+74035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[31] = int_stack + 37050;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+38850,int_stack+42675,int_stack+41100, 0.0, zero_stack, 1.0, int_stack+76210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[30] = int_stack + 38850;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+40650,int_stack+44250,int_stack+78385, 0.0, zero_stack, 1.0, int_stack+1350, 1.0, int_stack+79960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[26] = int_stack + 40650;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+76810,int_stack+45825,int_stack+43575, 0.0, zero_stack, 2.0, int_stack+79960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[25] = int_stack + 76810;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+42450,int_stack+47400,int_stack+45150, 1.0, int_stack+63985, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[23] = int_stack + 42450;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+44250,int_stack+48975,int_stack+46725, 1.0, int_stack+67060, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[22] = int_stack + 44250;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+46050,int_stack+50550,int_stack+48300, 1.0, int_stack+69235, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[21] = int_stack + 46050;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+47850,int_stack+52125,int_stack+49875, 1.0, int_stack+71410, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[20] = int_stack + 47850;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+49650,int_stack+53700,int_stack+51450, 1.0, int_stack+74035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[19] = int_stack + 49650;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+73810,int_stack+54600,int_stack+53025, 1.0, int_stack+76210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[18] = int_stack + 73810;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+51450,int_stack+55500,int_stack+2250, 1.0, int_stack+1350, 0.0, zero_stack, 1.0, int_stack+82135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[14] = int_stack + 51450;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+53250,int_stack+57075,int_stack+80560, 1.0, int_stack+79960, 1.0, int_stack+82135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[13] = int_stack + 53250;
 /*--- compute (fp|fd) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+58650,int_stack+56400, 2.0, int_stack+82135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,60);
     Libderiv->ABCD[12] = int_stack + 0;

}
