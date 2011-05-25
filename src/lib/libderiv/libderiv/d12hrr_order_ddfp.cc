#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_ddfp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dd|fp) integrals */

void d12hrr_order_ddfp(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv2_classes[2][3][143] = int_stack + 2325;
 Libderiv->deriv2_classes[2][4][143] = int_stack + 2385;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 2475;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 2575;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 2725;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 2875;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 3100;
 Libderiv->deriv2_classes[2][4][131] = int_stack + 3160;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 3250;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 3350;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 3500;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 3650;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 3875;
 Libderiv->deriv2_classes[2][4][130] = int_stack + 3935;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 4025;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 4125;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 4275;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 4425;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 4650;
 Libderiv->deriv2_classes[2][4][119] = int_stack + 4710;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 4800;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 4900;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 5050;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 5200;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 5425;
 Libderiv->deriv2_classes[2][4][118] = int_stack + 5485;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 5575;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 5675;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 5825;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 5975;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 6200;
 Libderiv->deriv2_classes[2][4][117] = int_stack + 6260;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 6350;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 6450;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 6600;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 6750;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 6975;
 Libderiv->deriv2_classes[2][4][107] = int_stack + 7035;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 7125;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 7225;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 7375;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 7525;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 7750;
 Libderiv->deriv2_classes[2][4][106] = int_stack + 7810;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 7900;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 8000;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 8150;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 8300;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 8525;
 Libderiv->deriv2_classes[2][4][105] = int_stack + 8585;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 8675;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 8775;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 8925;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 9075;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 9300;
 Libderiv->deriv2_classes[2][4][104] = int_stack + 9360;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 9450;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 9550;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 9700;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 9850;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 10075;
 Libderiv->deriv2_classes[2][4][95] = int_stack + 10135;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 10225;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 10325;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 10475;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 10625;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 10850;
 Libderiv->deriv2_classes[2][4][94] = int_stack + 10910;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 11000;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 11100;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 11250;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 11400;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 11625;
 Libderiv->deriv2_classes[2][4][93] = int_stack + 11685;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 11775;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 11875;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 12025;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 12175;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 12400;
 Libderiv->deriv2_classes[2][4][92] = int_stack + 12460;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 12550;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 12650;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 12800;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 12950;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 13175;
 Libderiv->deriv2_classes[2][4][91] = int_stack + 13235;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 13325;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 13425;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 13575;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 13725;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 13950;
 Libderiv->deriv2_classes[2][4][83] = int_stack + 14010;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 14100;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 14200;
 Libderiv->deriv_classes[4][3][11] = int_stack + 14350;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 14500;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 14650;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 14875;
 Libderiv->deriv2_classes[2][4][82] = int_stack + 14935;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 15025;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 15125;
 Libderiv->deriv_classes[4][3][10] = int_stack + 15275;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 15425;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 15575;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 15800;
 Libderiv->deriv2_classes[2][4][81] = int_stack + 15860;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 15950;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 16050;
 Libderiv->deriv_classes[4][3][9] = int_stack + 16200;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 16350;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 16500;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 16725;
 Libderiv->deriv2_classes[2][4][80] = int_stack + 16785;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 16875;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 16975;
 Libderiv->deriv_classes[4][3][8] = int_stack + 17125;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 17275;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 17425;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 17650;
 Libderiv->deriv2_classes[2][4][79] = int_stack + 17710;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 17800;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 17900;
 Libderiv->deriv_classes[4][3][7] = int_stack + 18050;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 18200;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 18350;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 18575;
 Libderiv->deriv2_classes[2][4][78] = int_stack + 18635;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 18725;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 18825;
 Libderiv->deriv_classes[4][3][6] = int_stack + 18975;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 19125;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 19275;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 19500;
 Libderiv->deriv2_classes[2][4][35] = int_stack + 19560;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 19650;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 19750;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 19900;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 20050;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 20275;
 Libderiv->deriv2_classes[2][4][34] = int_stack + 20335;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 20425;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 20525;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 20675;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 20825;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 21050;
 Libderiv->deriv2_classes[2][4][33] = int_stack + 21110;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 21200;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 21300;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 21450;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 21600;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 21825;
 Libderiv->deriv2_classes[2][4][32] = int_stack + 21885;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 21975;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 22075;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 22225;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 22375;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 22600;
 Libderiv->deriv2_classes[2][4][31] = int_stack + 22660;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 22750;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 22850;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 23000;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 23150;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 23375;
 Libderiv->deriv2_classes[2][4][30] = int_stack + 23435;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 23525;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 23625;
 Libderiv->deriv_classes[4][3][2] = int_stack + 23775;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 23925;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 24075;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 24300;
 Libderiv->deriv2_classes[2][4][26] = int_stack + 24360;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 24450;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 24550;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 24700;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 24850;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 25075;
 Libderiv->deriv2_classes[2][4][23] = int_stack + 25135;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 25225;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 25325;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 25475;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 25625;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 25850;
 Libderiv->deriv2_classes[2][4][22] = int_stack + 25910;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 26000;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 26100;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 26250;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 26400;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 26625;
 Libderiv->deriv2_classes[2][4][21] = int_stack + 26685;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 26775;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 26875;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 27025;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 27175;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 27400;
 Libderiv->deriv2_classes[2][4][20] = int_stack + 27460;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 27550;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 27650;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 27800;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 27950;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 28175;
 Libderiv->deriv2_classes[2][4][19] = int_stack + 28235;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 28325;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 28425;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 28575;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 28725;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 28950;
 Libderiv->deriv2_classes[2][4][18] = int_stack + 29010;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 29100;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 29200;
 Libderiv->deriv_classes[4][3][1] = int_stack + 29350;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 29500;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 29650;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 29875;
 Libderiv->deriv2_classes[2][4][14] = int_stack + 29935;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 30025;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 30125;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 30275;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 30425;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 30650;
 Libderiv->deriv2_classes[2][4][13] = int_stack + 30710;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 30800;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 30900;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 31050;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 31200;
 Libderiv->deriv_classes[2][3][11] = int_stack + 31425;
 Libderiv->deriv_classes[2][4][11] = int_stack + 31485;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 31575;
 Libderiv->deriv2_classes[2][4][11] = int_stack + 31635;
 Libderiv->deriv_classes[3][3][11] = int_stack + 31725;
 Libderiv->deriv_classes[3][4][11] = int_stack + 31825;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 31975;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 32075;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 32225;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 32375;
 Libderiv->deriv_classes[2][3][10] = int_stack + 32600;
 Libderiv->deriv_classes[2][4][10] = int_stack + 32660;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 32750;
 Libderiv->deriv2_classes[2][4][10] = int_stack + 32810;
 Libderiv->deriv_classes[3][3][10] = int_stack + 32900;
 Libderiv->deriv_classes[3][4][10] = int_stack + 33000;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 33150;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 33250;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 33400;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 33550;
 Libderiv->deriv_classes[2][3][9] = int_stack + 33775;
 Libderiv->deriv_classes[2][4][9] = int_stack + 33835;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 33925;
 Libderiv->deriv2_classes[2][4][9] = int_stack + 33985;
 Libderiv->deriv_classes[3][3][9] = int_stack + 34075;
 Libderiv->deriv_classes[3][4][9] = int_stack + 34175;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 34325;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 34425;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 34575;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 34725;
 Libderiv->deriv_classes[2][3][8] = int_stack + 34950;
 Libderiv->deriv_classes[2][4][8] = int_stack + 35010;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 35100;
 Libderiv->deriv2_classes[2][4][8] = int_stack + 35160;
 Libderiv->deriv_classes[3][3][8] = int_stack + 35250;
 Libderiv->deriv_classes[3][4][8] = int_stack + 35350;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 35500;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 35600;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 35750;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 35900;
 Libderiv->deriv_classes[2][3][7] = int_stack + 36125;
 Libderiv->deriv_classes[2][4][7] = int_stack + 36185;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 36275;
 Libderiv->deriv2_classes[2][4][7] = int_stack + 36335;
 Libderiv->deriv_classes[3][3][7] = int_stack + 36425;
 Libderiv->deriv_classes[3][4][7] = int_stack + 36525;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 36675;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 36775;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 36925;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 37075;
 Libderiv->deriv_classes[2][3][6] = int_stack + 37300;
 Libderiv->deriv_classes[2][4][6] = int_stack + 37360;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 37450;
 Libderiv->deriv2_classes[2][4][6] = int_stack + 37510;
 Libderiv->dvrr_classes[3][3] = int_stack + 37600;
 Libderiv->deriv_classes[3][3][6] = int_stack + 37700;
 Libderiv->deriv_classes[3][4][6] = int_stack + 37800;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 37950;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 38050;
 Libderiv->deriv_classes[4][3][0] = int_stack + 38200;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 38350;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 38500;
 Libderiv->deriv_classes[2][3][2] = int_stack + 38725;
 Libderiv->deriv_classes[2][4][2] = int_stack + 38785;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 38875;
 Libderiv->deriv2_classes[2][4][2] = int_stack + 38935;
 Libderiv->deriv_classes[3][3][2] = int_stack + 39025;
 Libderiv->deriv_classes[3][4][2] = int_stack + 39125;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 39275;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 39375;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 39525;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 39675;
 Libderiv->deriv_classes[2][3][1] = int_stack + 39900;
 Libderiv->deriv_classes[2][4][1] = int_stack + 39960;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 40050;
 Libderiv->deriv2_classes[2][4][1] = int_stack + 40110;
 Libderiv->deriv_classes[3][3][1] = int_stack + 40200;
 Libderiv->deriv_classes[3][4][1] = int_stack + 40300;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 40450;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 40550;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 40700;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 40850;
 Libderiv->dvrr_classes[2][3] = int_stack + 41075;
 Libderiv->dvrr_classes[2][4] = int_stack + 41135;
 Libderiv->deriv_classes[2][3][0] = int_stack + 41225;
 Libderiv->deriv_classes[2][4][0] = int_stack + 41285;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 41375;
 Libderiv->deriv2_classes[2][4][0] = int_stack + 41435;
 Libderiv->deriv_classes[3][3][0] = int_stack + 41525;
 Libderiv->deriv_classes[3][4][0] = int_stack + 41625;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 41775;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 41875;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 42025;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 42175;
 memset(int_stack,0,339200);

 Libderiv->dvrr_stack = int_stack + 85270;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_ddfp(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42400,int_stack+31485,int_stack+31425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41075,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42580,int_stack+31825,int_stack+31725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37600,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+42880,int_stack+42580,int_stack+42400,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43420,int_stack+0,int_stack+14350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+43870,int_stack+43420,int_stack+42580,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43420,int_stack+32660,int_stack+32600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41075, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44770,int_stack+33000,int_stack+32900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37600, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+45070,int_stack+44770,int_stack+43420,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45610,int_stack+225,int_stack+15275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+46060,int_stack+45610,int_stack+44770,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45610,int_stack+33835,int_stack+33775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41075, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+34175,int_stack+34075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37600, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+46960,int_stack+0,int_stack+45610,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47500,int_stack+450,int_stack+16200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+47950,int_stack+47500,int_stack+0,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47500,int_stack+35010,int_stack+34950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+35350,int_stack+35250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+48850,int_stack+300,int_stack+47500,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49390,int_stack+675,int_stack+17125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+49840,int_stack+49390,int_stack+300,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49390,int_stack+36185,int_stack+36125, 0.0, zero_stack, 1.0, int_stack+41075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+36525,int_stack+36425, 0.0, zero_stack, 1.0, int_stack+37600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+50740,int_stack+600,int_stack+49390,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51280,int_stack+900,int_stack+18050, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+51730,int_stack+51280,int_stack+600,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51280,int_stack+37360,int_stack+37300, 1.0, int_stack+41075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+52630,int_stack+37800,int_stack+37700, 1.0, int_stack+37600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+52930,int_stack+52630,int_stack+51280,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+1275,int_stack+18975, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+53920,int_stack+53470,int_stack+52630,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+41135,int_stack+41075,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1950,int_stack+37600,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+54820,int_stack+900,int_stack+53470,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53650,int_stack+38785,int_stack+38725,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+39125,int_stack+39025,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+55360,int_stack+1200,int_stack+53650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+55900,int_stack+1500,int_stack+23775,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+56350,int_stack+55900,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+55900,int_stack+39960,int_stack+39900,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+57250,int_stack+40300,int_stack+40200,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+57550,int_stack+57250,int_stack+55900, 0.0, zero_stack, 1.0, int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+58090,int_stack+1725,int_stack+29350,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+58540,int_stack+58090,int_stack+57250, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+58090,int_stack+41285,int_stack+41225,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1500,int_stack+41625,int_stack+41525,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+59440,int_stack+1500,int_stack+58090, 1.0, int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+2100,int_stack+38200,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+60430,int_stack+59980,int_stack+1500, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+2385,int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+31425,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+2575,int_stack+2475, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+31725,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1800,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+2875,int_stack+2725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14350,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+61330,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+3160,int_stack+3100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31425, 1.0, int_stack+32600,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+3350,int_stack+3250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31725, 1.0, int_stack+32900,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2340,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+3650,int_stack+3500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14350, 1.0, int_stack+15275,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+2880,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+3935,int_stack+3875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+32600, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+4125,int_stack+4025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+32900, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+62230,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+4425,int_stack+4275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+15275, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+62770,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+4710,int_stack+4650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31425, 0.0, zero_stack, 1.0, int_stack+33775,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+4900,int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31725, 0.0, zero_stack, 1.0, int_stack+34075,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3780,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+5200,int_stack+5050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14350, 0.0, zero_stack, 1.0, int_stack+16200,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+4320,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+5485,int_stack+5425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32600, 1.0, int_stack+33775, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+5675,int_stack+5575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32900, 1.0, int_stack+34075, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5220,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+5975,int_stack+5825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15275, 1.0, int_stack+16200, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+63670,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+6260,int_stack+6200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+33775, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+6450,int_stack+6350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+34075, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+5760,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+6750,int_stack+6600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+16200, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+64570,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+7035,int_stack+6975, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31425, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34950,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+7225,int_stack+7125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+35250,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+6300,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+7525,int_stack+7375, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17125,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+6840,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+7810,int_stack+7750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32600, 0.0, zero_stack, 1.0, int_stack+34950, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+8000,int_stack+7900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32900, 0.0, zero_stack, 1.0, int_stack+35250, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+65470,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+8300,int_stack+8150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15275, 0.0, zero_stack, 1.0, int_stack+17125, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+66010,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+8585,int_stack+8525, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33775, 1.0, int_stack+34950, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+8775,int_stack+8675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34075, 1.0, int_stack+35250, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+7740,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+9075,int_stack+8925, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16200, 1.0, int_stack+17125, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8280,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+9360,int_stack+9300, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+34950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+9550,int_stack+9450, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+35250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+66910,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+9850,int_stack+9700, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+17125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+67450,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+10135,int_stack+10075, 0.0, zero_stack, 1.0, int_stack+31425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36125,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+10325,int_stack+10225, 0.0, zero_stack, 1.0, int_stack+31725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36425,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9180,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+10625,int_stack+10475, 0.0, zero_stack, 1.0, int_stack+14350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18050,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+9720,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+10910,int_stack+10850, 0.0, zero_stack, 1.0, int_stack+32600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36125, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+11100,int_stack+11000, 0.0, zero_stack, 1.0, int_stack+32900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+36425, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+10620,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+11400,int_stack+11250, 0.0, zero_stack, 1.0, int_stack+15275, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18050, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+68350,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+11685,int_stack+11625, 0.0, zero_stack, 1.0, int_stack+33775, 0.0, zero_stack, 1.0, int_stack+36125, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+11875,int_stack+11775, 0.0, zero_stack, 1.0, int_stack+34075, 0.0, zero_stack, 1.0, int_stack+36425, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+11160,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+12175,int_stack+12025, 0.0, zero_stack, 1.0, int_stack+16200, 0.0, zero_stack, 1.0, int_stack+18050, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+69250,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+12460,int_stack+12400, 0.0, zero_stack, 1.0, int_stack+34950, 1.0, int_stack+36125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+12650,int_stack+12550, 0.0, zero_stack, 1.0, int_stack+35250, 1.0, int_stack+36425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+11700,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+12950,int_stack+12800, 0.0, zero_stack, 1.0, int_stack+17125, 1.0, int_stack+18050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12240,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+13235,int_stack+13175, 0.0, zero_stack, 2.0, int_stack+36125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+13425,int_stack+13325, 0.0, zero_stack, 2.0, int_stack+36425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+70150,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+13725,int_stack+13575, 0.0, zero_stack, 2.0, int_stack+18050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+70690,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+14010,int_stack+13950, 1.0, int_stack+31425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37300,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+14200,int_stack+14100, 1.0, int_stack+31725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37700,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+13140,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+14650,int_stack+14500, 1.0, int_stack+14350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18975,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+13680,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+14935,int_stack+14875, 1.0, int_stack+32600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37300, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+15125,int_stack+15025, 1.0, int_stack+32900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37700, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+14580,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+15575,int_stack+15425, 1.0, int_stack+15275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18975, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+71590,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+15860,int_stack+15800, 1.0, int_stack+33775, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37300, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+16050,int_stack+15950, 1.0, int_stack+34075, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37700, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+15120,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+16500,int_stack+16350, 1.0, int_stack+16200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+18975, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+15660,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+16785,int_stack+16725, 1.0, int_stack+34950, 0.0, zero_stack, 1.0, int_stack+37300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+16975,int_stack+16875, 1.0, int_stack+35250, 0.0, zero_stack, 1.0, int_stack+37700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+16560,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+17425,int_stack+17275, 1.0, int_stack+17125, 0.0, zero_stack, 1.0, int_stack+18975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+72490,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+17710,int_stack+17650, 1.0, int_stack+36125, 1.0, int_stack+37300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+17900,int_stack+17800, 1.0, int_stack+36425, 1.0, int_stack+37700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+17100,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+18350,int_stack+18200, 1.0, int_stack+18050, 1.0, int_stack+18975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+17640,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+18635,int_stack+18575, 2.0, int_stack+37300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+18825,int_stack+18725, 2.0, int_stack+37700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+73390,int_stack+900,int_stack+53470,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+19275,int_stack+19125, 2.0, int_stack+18975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+18540,int_stack+59980,int_stack+900,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+19560,int_stack+19500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38725,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+19750,int_stack+19650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39025,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+73930,int_stack+900,int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+20050,int_stack+19900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23775,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+74470,int_stack+59980,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+20335,int_stack+20275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38725, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+20525,int_stack+20425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39025, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+19440,int_stack+900,int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+43420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+20825,int_stack+20675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23775, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+19980,int_stack+59980,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+44770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+21110,int_stack+21050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38725, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+21300,int_stack+21200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39025, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+20880,int_stack+900,int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+21600,int_stack+21450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23775, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+75370,int_stack+59980,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+21885,int_stack+21825, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+22075,int_stack+21975, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+21420,int_stack+900,int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+47500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+22375,int_stack+22225, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+23775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+76270,int_stack+59980,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+22660,int_stack+22600, 0.0, zero_stack, 1.0, int_stack+38725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+22850,int_stack+22750, 0.0, zero_stack, 1.0, int_stack+39025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+21960,int_stack+900,int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+49390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+23150,int_stack+23000, 0.0, zero_stack, 1.0, int_stack+23775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+77170,int_stack+59980,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+23435,int_stack+23375, 1.0, int_stack+38725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+23625,int_stack+23525, 1.0, int_stack+39025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+22500,int_stack+900,int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+51280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+24075,int_stack+23925, 1.0, int_stack+23775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+23040,int_stack+59980,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+24360,int_stack+24300,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+24550,int_stack+24450,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+23940,int_stack+900,int_stack+53470, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+53650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+24850,int_stack+24700,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+78070,int_stack+59980,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+25135,int_stack+25075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39900,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+25325,int_stack+25225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40200,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+24480,int_stack+900,int_stack+53470, 0.0, zero_stack, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+25625,int_stack+25475, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29350,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+78970,int_stack+59980,int_stack+900, 0.0, zero_stack, 1.0, int_stack+42580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+25910,int_stack+25850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39900, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+26100,int_stack+26000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40200, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+25020,int_stack+900,int_stack+53470, 0.0, zero_stack, 1.0, int_stack+43420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+26400,int_stack+26250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29350, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+25560,int_stack+59980,int_stack+900, 0.0, zero_stack, 1.0, int_stack+44770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+26685,int_stack+26625, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39900, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+26875,int_stack+26775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40200, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+26460,int_stack+900,int_stack+53470, 0.0, zero_stack, 1.0, int_stack+45610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+27175,int_stack+27025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29350, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+79870,int_stack+59980,int_stack+900, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+27460,int_stack+27400, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+39900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+27650,int_stack+27550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+40200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+27000,int_stack+900,int_stack+53470, 0.0, zero_stack, 1.0, int_stack+47500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+27950,int_stack+27800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+80770,int_stack+59980,int_stack+900, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+28235,int_stack+28175, 0.0, zero_stack, 1.0, int_stack+39900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+28425,int_stack+28325, 0.0, zero_stack, 1.0, int_stack+40200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+27540,int_stack+900,int_stack+53470, 0.0, zero_stack, 1.0, int_stack+49390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+28725,int_stack+28575, 0.0, zero_stack, 1.0, int_stack+29350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+81670,int_stack+59980,int_stack+900, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+29010,int_stack+28950, 1.0, int_stack+39900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+29200,int_stack+29100, 1.0, int_stack+40200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+28080,int_stack+900,int_stack+53470, 0.0, zero_stack, 1.0, int_stack+51280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+29650,int_stack+29500, 1.0, int_stack+29350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+28620,int_stack+59980,int_stack+900, 0.0, zero_stack, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+29935,int_stack+29875,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+30125,int_stack+30025,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+29520,int_stack+900,int_stack+53470, 0.0, zero_stack, 1.0, int_stack+53650, 1.0, int_stack+55900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+30425,int_stack+30275,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+82570,int_stack+59980,int_stack+900, 0.0, zero_stack, 1.0, int_stack+1200, 1.0, int_stack+57250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+30710,int_stack+30650,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+30900,int_stack+30800,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+30060,int_stack+900,int_stack+53470, 0.0, zero_stack, 2.0, int_stack+55900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+31200,int_stack+31050,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+30600,int_stack+59980,int_stack+900, 0.0, zero_stack, 2.0, int_stack+57250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+31635,int_stack+31575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41225,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+32075,int_stack+31975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41525,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+31500,int_stack+900,int_stack+53470, 1.0, int_stack+42400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+59980,int_stack+32375,int_stack+32225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38200,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+83470,int_stack+59980,int_stack+900, 1.0, int_stack+42580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+32810,int_stack+32750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41225, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+33250,int_stack+33150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41525, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+32040,int_stack+900,int_stack+53470, 1.0, int_stack+43420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+43420,int_stack+33550,int_stack+33400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38200, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+32580,int_stack+43420,int_stack+900, 1.0, int_stack+44770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+33985,int_stack+33925, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41225, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44770,int_stack+34425,int_stack+34325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41525, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+33480,int_stack+44770,int_stack+53470, 1.0, int_stack+45610, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+45610,int_stack+34725,int_stack+34575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38200, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+34020,int_stack+45610,int_stack+44770, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+35160,int_stack+35100, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+35600,int_stack+35500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+41525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+34920,int_stack+0,int_stack+53470, 1.0, int_stack+47500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+47500,int_stack+35900,int_stack+35750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+38200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+84370,int_stack+47500,int_stack+0, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+36335,int_stack+36275, 0.0, zero_stack, 1.0, int_stack+41225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44770,int_stack+36775,int_stack+36675, 0.0, zero_stack, 1.0, int_stack+41525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+44770,int_stack+53470, 1.0, int_stack+49390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+49390,int_stack+37075,int_stack+36925, 0.0, zero_stack, 1.0, int_stack+38200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+35460,int_stack+49390,int_stack+44770, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+37510,int_stack+37450, 1.0, int_stack+41225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+44770,int_stack+38050,int_stack+37950, 1.0, int_stack+41525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+540,int_stack+44770,int_stack+53470, 1.0, int_stack+51280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+51280,int_stack+38500,int_stack+38350, 1.0, int_stack+38200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+36360,int_stack+51280,int_stack+44770, 1.0, int_stack+52630, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+38935,int_stack+38875,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52630,int_stack+39375,int_stack+39275,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+37260,int_stack+52630,int_stack+53470, 1.0, int_stack+53650, 0.0, zero_stack, 1.0, int_stack+58090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+53470,int_stack+39675,int_stack+39525,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+37800,int_stack+53470,int_stack+52630, 1.0, int_stack+1200, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+52630,int_stack+40110,int_stack+40050,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+44770,int_stack+40550,int_stack+40450,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+38700,int_stack+44770,int_stack+52630, 1.0, int_stack+55900, 1.0, int_stack+58090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+55900,int_stack+40850,int_stack+40700,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+39240,int_stack+55900,int_stack+44770, 1.0, int_stack+57250, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+57250,int_stack+41435,int_stack+41375,6);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+44770,int_stack+41875,int_stack+41775,10);
 /*--- compute (dp|fp) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+40140,int_stack+44770,int_stack+57250, 2.0, int_stack+58090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+58090,int_stack+42175,int_stack+42025,15);
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+40680,int_stack+58090,int_stack+44770, 2.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+41580,int_stack+43870,int_stack+42880,30);
     Libderiv->ABCD[11] = int_stack + 41580;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+43420,int_stack+46060,int_stack+45070,30);
     Libderiv->ABCD[10] = int_stack + 43420;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+45610,int_stack+47950,int_stack+46960,30);
     Libderiv->ABCD[9] = int_stack + 45610;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+47500,int_stack+49840,int_stack+48850,30);
     Libderiv->ABCD[8] = int_stack + 47500;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+49390,int_stack+51730,int_stack+50740,30);
     Libderiv->ABCD[7] = int_stack + 49390;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+51280,int_stack+53920,int_stack+52930,30);
     Libderiv->ABCD[6] = int_stack + 51280;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+53470,int_stack+56350,int_stack+55360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+54820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[2] = int_stack + 53470;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+55900,int_stack+58540,int_stack+57550, 0.0, zero_stack, 1.0, int_stack+54820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[1] = int_stack + 55900;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+58090,int_stack+60430,int_stack+59440, 1.0, int_stack+54820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[0] = int_stack + 58090;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+59980,int_stack+61330,int_stack+1800,30);
     Libderiv->ABCD[155] = int_stack + 59980;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+1080,int_stack+2880,int_stack+2340,30);
     Libderiv->ABCD[143] = int_stack + 1080;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+2160,int_stack+62770,int_stack+62230,30);
     Libderiv->ABCD[142] = int_stack + 2160;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+61060,int_stack+4320,int_stack+3780,30);
     Libderiv->ABCD[131] = int_stack + 61060;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+3240,int_stack+63670,int_stack+5220,30);
     Libderiv->ABCD[130] = int_stack + 3240;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+4320,int_stack+64570,int_stack+5760,30);
     Libderiv->ABCD[129] = int_stack + 4320;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+62140,int_stack+6840,int_stack+6300,30);
     Libderiv->ABCD[119] = int_stack + 62140;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+5400,int_stack+66010,int_stack+65470,30);
     Libderiv->ABCD[118] = int_stack + 5400;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+6480,int_stack+8280,int_stack+7740,30);
     Libderiv->ABCD[117] = int_stack + 6480;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+7560,int_stack+67450,int_stack+66910,30);
     Libderiv->ABCD[116] = int_stack + 7560;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+63220,int_stack+9720,int_stack+9180,30);
     Libderiv->ABCD[107] = int_stack + 63220;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+8640,int_stack+68350,int_stack+10620,30);
     Libderiv->ABCD[106] = int_stack + 8640;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+9720,int_stack+69250,int_stack+11160,30);
     Libderiv->ABCD[105] = int_stack + 9720;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+64300,int_stack+12240,int_stack+11700,30);
     Libderiv->ABCD[104] = int_stack + 64300;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+10800,int_stack+70690,int_stack+70150,30);
     Libderiv->ABCD[103] = int_stack + 10800;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+11880,int_stack+13680,int_stack+13140,30);
     Libderiv->ABCD[95] = int_stack + 11880;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+12960,int_stack+71590,int_stack+14580,30);
     Libderiv->ABCD[94] = int_stack + 12960;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+14040,int_stack+15660,int_stack+15120,30);
     Libderiv->ABCD[93] = int_stack + 14040;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+15120,int_stack+72490,int_stack+16560,30);
     Libderiv->ABCD[92] = int_stack + 15120;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+65380,int_stack+17640,int_stack+17100,30);
     Libderiv->ABCD[91] = int_stack + 65380;
 /*--- compute (dd|fp) ---*/
   hrr1_build_dd(Libderiv->AB,int_stack+16200,int_stack+18540,int_stack+73390,30);
     Libderiv->ABCD[90] = int_stack + 16200;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+17280,int_stack+74470,int_stack+73930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+42880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[47] = int_stack + 17280;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+18360,int_stack+19980,int_stack+19440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+45070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[46] = int_stack + 18360;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+19440,int_stack+75370,int_stack+20880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+46960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[45] = int_stack + 19440;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+66460,int_stack+76270,int_stack+21420, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+48850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[44] = int_stack + 66460;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+20520,int_stack+77170,int_stack+21960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+50740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[43] = int_stack + 20520;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+67540,int_stack+23040,int_stack+22500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[42] = int_stack + 67540;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+21600,int_stack+78070,int_stack+23940, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+55360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[38] = int_stack + 21600;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+22680,int_stack+78970,int_stack+24480, 0.0, zero_stack, 1.0, int_stack+42880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[35] = int_stack + 22680;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+23760,int_stack+25560,int_stack+25020, 0.0, zero_stack, 1.0, int_stack+45070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[34] = int_stack + 23760;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+24840,int_stack+79870,int_stack+26460, 0.0, zero_stack, 1.0, int_stack+46960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[33] = int_stack + 24840;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+25920,int_stack+80770,int_stack+27000, 0.0, zero_stack, 1.0, int_stack+48850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[32] = int_stack + 25920;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+68620,int_stack+81670,int_stack+27540, 0.0, zero_stack, 1.0, int_stack+50740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[31] = int_stack + 68620;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+27000,int_stack+28620,int_stack+28080, 0.0, zero_stack, 1.0, int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[30] = int_stack + 27000;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+28080,int_stack+82570,int_stack+29520, 0.0, zero_stack, 1.0, int_stack+55360, 1.0, int_stack+57550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[26] = int_stack + 28080;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+69700,int_stack+30600,int_stack+30060, 0.0, zero_stack, 2.0, int_stack+57550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[25] = int_stack + 69700;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+29160,int_stack+83470,int_stack+31500, 1.0, int_stack+42880, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[23] = int_stack + 29160;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+30240,int_stack+32580,int_stack+32040, 1.0, int_stack+45070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[22] = int_stack + 30240;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+31320,int_stack+34020,int_stack+33480, 1.0, int_stack+46960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[21] = int_stack + 31320;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+32400,int_stack+84370,int_stack+34920, 1.0, int_stack+48850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[20] = int_stack + 32400;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+33480,int_stack+35460,int_stack+0, 1.0, int_stack+50740, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[19] = int_stack + 33480;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+34560,int_stack+36360,int_stack+540, 1.0, int_stack+52930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[18] = int_stack + 34560;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+0,int_stack+37800,int_stack+37260, 1.0, int_stack+55360, 0.0, zero_stack, 1.0, int_stack+59440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[14] = int_stack + 0;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+35640,int_stack+39240,int_stack+38700, 1.0, int_stack+57550, 1.0, int_stack+59440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[13] = int_stack + 35640;
 /*--- compute (dd|fp) ---*/
   d1hrr1_build_dd(Libderiv->AB,int_stack+36720,int_stack+40680,int_stack+40140, 2.0, int_stack+59440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[12] = int_stack + 36720;

}
