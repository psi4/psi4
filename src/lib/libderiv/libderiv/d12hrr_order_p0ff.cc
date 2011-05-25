#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0ff(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|ff) integrals */

void d12hrr_order_p0ff(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][6][11] = int_stack + 0;
 Libderiv->deriv_classes[1][6][10] = int_stack + 84;
 Libderiv->deriv_classes[1][6][9] = int_stack + 168;
 Libderiv->deriv_classes[1][6][8] = int_stack + 252;
 Libderiv->deriv_classes[1][6][7] = int_stack + 336;
 Libderiv->dvrr_classes[1][5] = int_stack + 420;
 Libderiv->deriv_classes[1][6][6] = int_stack + 483;
 Libderiv->deriv_classes[1][6][2] = int_stack + 567;
 Libderiv->deriv_classes[1][6][1] = int_stack + 651;
 Libderiv->deriv_classes[1][6][0] = int_stack + 735;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 819;
 Libderiv->deriv2_classes[1][4][143] = int_stack + 849;
 Libderiv->deriv2_classes[1][5][143] = int_stack + 894;
 Libderiv->deriv2_classes[1][6][143] = int_stack + 957;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 1041;
 Libderiv->deriv2_classes[1][4][131] = int_stack + 1071;
 Libderiv->deriv2_classes[1][5][131] = int_stack + 1116;
 Libderiv->deriv2_classes[1][6][131] = int_stack + 1179;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 1263;
 Libderiv->deriv2_classes[1][4][130] = int_stack + 1293;
 Libderiv->deriv2_classes[1][5][130] = int_stack + 1338;
 Libderiv->deriv2_classes[1][6][130] = int_stack + 1401;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 1485;
 Libderiv->deriv2_classes[1][4][119] = int_stack + 1515;
 Libderiv->deriv2_classes[1][5][119] = int_stack + 1560;
 Libderiv->deriv2_classes[1][6][119] = int_stack + 1623;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 1707;
 Libderiv->deriv2_classes[1][4][118] = int_stack + 1737;
 Libderiv->deriv2_classes[1][5][118] = int_stack + 1782;
 Libderiv->deriv2_classes[1][6][118] = int_stack + 1845;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 1929;
 Libderiv->deriv2_classes[1][4][117] = int_stack + 1959;
 Libderiv->deriv2_classes[1][5][117] = int_stack + 2004;
 Libderiv->deriv2_classes[1][6][117] = int_stack + 2067;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 2151;
 Libderiv->deriv2_classes[1][4][107] = int_stack + 2181;
 Libderiv->deriv2_classes[1][5][107] = int_stack + 2226;
 Libderiv->deriv2_classes[1][6][107] = int_stack + 2289;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 2373;
 Libderiv->deriv2_classes[1][4][106] = int_stack + 2403;
 Libderiv->deriv2_classes[1][5][106] = int_stack + 2448;
 Libderiv->deriv2_classes[1][6][106] = int_stack + 2511;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 2595;
 Libderiv->deriv2_classes[1][4][105] = int_stack + 2625;
 Libderiv->deriv2_classes[1][5][105] = int_stack + 2670;
 Libderiv->deriv2_classes[1][6][105] = int_stack + 2733;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 2817;
 Libderiv->deriv2_classes[1][4][104] = int_stack + 2847;
 Libderiv->deriv2_classes[1][5][104] = int_stack + 2892;
 Libderiv->deriv2_classes[1][6][104] = int_stack + 2955;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 3039;
 Libderiv->deriv2_classes[1][4][95] = int_stack + 3069;
 Libderiv->deriv2_classes[1][5][95] = int_stack + 3114;
 Libderiv->deriv2_classes[1][6][95] = int_stack + 3177;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 3261;
 Libderiv->deriv2_classes[1][4][94] = int_stack + 3291;
 Libderiv->deriv2_classes[1][5][94] = int_stack + 3336;
 Libderiv->deriv2_classes[1][6][94] = int_stack + 3399;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 3483;
 Libderiv->deriv2_classes[1][4][93] = int_stack + 3513;
 Libderiv->deriv2_classes[1][5][93] = int_stack + 3558;
 Libderiv->deriv2_classes[1][6][93] = int_stack + 3621;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 3705;
 Libderiv->deriv2_classes[1][4][92] = int_stack + 3735;
 Libderiv->deriv2_classes[1][5][92] = int_stack + 3780;
 Libderiv->deriv2_classes[1][6][92] = int_stack + 3843;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 3927;
 Libderiv->deriv2_classes[1][4][91] = int_stack + 3957;
 Libderiv->deriv2_classes[1][5][91] = int_stack + 4002;
 Libderiv->deriv2_classes[1][6][91] = int_stack + 4065;
 Libderiv->deriv_classes[1][3][11] = int_stack + 4149;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 4179;
 Libderiv->deriv_classes[1][4][11] = int_stack + 4209;
 Libderiv->deriv2_classes[1][4][83] = int_stack + 4254;
 Libderiv->deriv_classes[1][5][11] = int_stack + 4299;
 Libderiv->deriv2_classes[1][5][83] = int_stack + 4362;
 Libderiv->deriv2_classes[1][6][83] = int_stack + 4425;
 Libderiv->deriv_classes[1][3][10] = int_stack + 4509;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 4539;
 Libderiv->deriv_classes[1][4][10] = int_stack + 4569;
 Libderiv->deriv2_classes[1][4][82] = int_stack + 4614;
 Libderiv->deriv_classes[1][5][10] = int_stack + 4659;
 Libderiv->deriv2_classes[1][5][82] = int_stack + 4722;
 Libderiv->deriv2_classes[1][6][82] = int_stack + 4785;
 Libderiv->deriv_classes[1][3][9] = int_stack + 4869;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 4899;
 Libderiv->deriv_classes[1][4][9] = int_stack + 4929;
 Libderiv->deriv2_classes[1][4][81] = int_stack + 4974;
 Libderiv->deriv_classes[1][5][9] = int_stack + 5019;
 Libderiv->deriv2_classes[1][5][81] = int_stack + 5082;
 Libderiv->deriv2_classes[1][6][81] = int_stack + 5145;
 Libderiv->deriv_classes[1][3][8] = int_stack + 5229;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 5259;
 Libderiv->deriv_classes[1][4][8] = int_stack + 5289;
 Libderiv->deriv2_classes[1][4][80] = int_stack + 5334;
 Libderiv->deriv_classes[1][5][8] = int_stack + 5379;
 Libderiv->deriv2_classes[1][5][80] = int_stack + 5442;
 Libderiv->deriv2_classes[1][6][80] = int_stack + 5505;
 Libderiv->deriv_classes[1][3][7] = int_stack + 5589;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 5619;
 Libderiv->deriv_classes[1][4][7] = int_stack + 5649;
 Libderiv->deriv2_classes[1][4][79] = int_stack + 5694;
 Libderiv->deriv_classes[1][5][7] = int_stack + 5739;
 Libderiv->deriv2_classes[1][5][79] = int_stack + 5802;
 Libderiv->deriv2_classes[1][6][79] = int_stack + 5865;
 Libderiv->dvrr_classes[1][3] = int_stack + 5949;
 Libderiv->deriv_classes[1][3][6] = int_stack + 5979;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 6009;
 Libderiv->dvrr_classes[1][4] = int_stack + 6039;
 Libderiv->deriv_classes[1][4][6] = int_stack + 6084;
 Libderiv->deriv2_classes[1][4][78] = int_stack + 6129;
 Libderiv->deriv_classes[1][5][6] = int_stack + 6174;
 Libderiv->deriv2_classes[1][5][78] = int_stack + 6237;
 Libderiv->deriv2_classes[1][6][78] = int_stack + 6300;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 6384;
 Libderiv->deriv2_classes[1][4][35] = int_stack + 6414;
 Libderiv->deriv2_classes[1][5][35] = int_stack + 6459;
 Libderiv->deriv2_classes[1][6][35] = int_stack + 6522;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 6606;
 Libderiv->deriv2_classes[1][4][34] = int_stack + 6636;
 Libderiv->deriv2_classes[1][5][34] = int_stack + 6681;
 Libderiv->deriv2_classes[1][6][34] = int_stack + 6744;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 6828;
 Libderiv->deriv2_classes[1][4][33] = int_stack + 6858;
 Libderiv->deriv2_classes[1][5][33] = int_stack + 6903;
 Libderiv->deriv2_classes[1][6][33] = int_stack + 6966;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 7050;
 Libderiv->deriv2_classes[1][4][32] = int_stack + 7080;
 Libderiv->deriv2_classes[1][5][32] = int_stack + 7125;
 Libderiv->deriv2_classes[1][6][32] = int_stack + 7188;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 7272;
 Libderiv->deriv2_classes[1][4][31] = int_stack + 7302;
 Libderiv->deriv2_classes[1][5][31] = int_stack + 7347;
 Libderiv->deriv2_classes[1][6][31] = int_stack + 7410;
 Libderiv->deriv_classes[1][3][2] = int_stack + 7494;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 7524;
 Libderiv->deriv_classes[1][4][2] = int_stack + 7554;
 Libderiv->deriv2_classes[1][4][30] = int_stack + 7599;
 Libderiv->deriv_classes[1][5][2] = int_stack + 7644;
 Libderiv->deriv2_classes[1][5][30] = int_stack + 7707;
 Libderiv->deriv2_classes[1][6][30] = int_stack + 7770;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 7854;
 Libderiv->deriv2_classes[1][4][26] = int_stack + 7884;
 Libderiv->deriv2_classes[1][5][26] = int_stack + 7929;
 Libderiv->deriv2_classes[1][6][26] = int_stack + 7992;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 8076;
 Libderiv->deriv2_classes[1][4][23] = int_stack + 8106;
 Libderiv->deriv2_classes[1][5][23] = int_stack + 8151;
 Libderiv->deriv2_classes[1][6][23] = int_stack + 8214;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 8298;
 Libderiv->deriv2_classes[1][4][22] = int_stack + 8328;
 Libderiv->deriv2_classes[1][5][22] = int_stack + 8373;
 Libderiv->deriv2_classes[1][6][22] = int_stack + 8436;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 8520;
 Libderiv->deriv2_classes[1][4][21] = int_stack + 8550;
 Libderiv->deriv2_classes[1][5][21] = int_stack + 8595;
 Libderiv->deriv2_classes[1][6][21] = int_stack + 8658;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 8742;
 Libderiv->deriv2_classes[1][4][20] = int_stack + 8772;
 Libderiv->deriv2_classes[1][5][20] = int_stack + 8817;
 Libderiv->deriv2_classes[1][6][20] = int_stack + 8880;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 8964;
 Libderiv->deriv2_classes[1][4][19] = int_stack + 8994;
 Libderiv->deriv2_classes[1][5][19] = int_stack + 9039;
 Libderiv->deriv2_classes[1][6][19] = int_stack + 9102;
 Libderiv->deriv_classes[1][3][1] = int_stack + 9186;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 9216;
 Libderiv->deriv_classes[1][4][1] = int_stack + 9246;
 Libderiv->deriv2_classes[1][4][18] = int_stack + 9291;
 Libderiv->deriv_classes[1][5][1] = int_stack + 9336;
 Libderiv->deriv2_classes[1][5][18] = int_stack + 9399;
 Libderiv->deriv2_classes[1][6][18] = int_stack + 9462;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 9546;
 Libderiv->deriv2_classes[1][4][14] = int_stack + 9576;
 Libderiv->deriv2_classes[1][5][14] = int_stack + 9621;
 Libderiv->deriv2_classes[1][6][14] = int_stack + 9684;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 9768;
 Libderiv->deriv2_classes[1][4][13] = int_stack + 9798;
 Libderiv->deriv2_classes[1][5][13] = int_stack + 9843;
 Libderiv->deriv2_classes[1][6][13] = int_stack + 9906;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 9990;
 Libderiv->deriv2_classes[1][4][11] = int_stack + 10020;
 Libderiv->deriv2_classes[1][5][11] = int_stack + 10065;
 Libderiv->deriv2_classes[1][6][11] = int_stack + 10128;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 10212;
 Libderiv->deriv2_classes[1][4][10] = int_stack + 10242;
 Libderiv->deriv2_classes[1][5][10] = int_stack + 10287;
 Libderiv->deriv2_classes[1][6][10] = int_stack + 10350;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 10434;
 Libderiv->deriv2_classes[1][4][9] = int_stack + 10464;
 Libderiv->deriv2_classes[1][5][9] = int_stack + 10509;
 Libderiv->deriv2_classes[1][6][9] = int_stack + 10572;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 10656;
 Libderiv->deriv2_classes[1][4][8] = int_stack + 10686;
 Libderiv->deriv2_classes[1][5][8] = int_stack + 10731;
 Libderiv->deriv2_classes[1][6][8] = int_stack + 10794;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 10878;
 Libderiv->deriv2_classes[1][4][7] = int_stack + 10908;
 Libderiv->deriv2_classes[1][5][7] = int_stack + 10953;
 Libderiv->deriv2_classes[1][6][7] = int_stack + 11016;
 Libderiv->deriv_classes[1][3][0] = int_stack + 11100;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 11130;
 Libderiv->deriv_classes[1][4][0] = int_stack + 11160;
 Libderiv->deriv2_classes[1][4][6] = int_stack + 11205;
 Libderiv->deriv_classes[1][5][0] = int_stack + 11250;
 Libderiv->deriv2_classes[1][5][6] = int_stack + 11313;
 Libderiv->deriv2_classes[1][6][6] = int_stack + 11376;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 11460;
 Libderiv->deriv2_classes[1][4][2] = int_stack + 11490;
 Libderiv->deriv2_classes[1][5][2] = int_stack + 11535;
 Libderiv->deriv2_classes[1][6][2] = int_stack + 11598;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 11682;
 Libderiv->deriv2_classes[1][4][1] = int_stack + 11712;
 Libderiv->deriv2_classes[1][5][1] = int_stack + 11757;
 Libderiv->deriv2_classes[1][6][1] = int_stack + 11820;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 11904;
 Libderiv->deriv2_classes[1][4][0] = int_stack + 11934;
 Libderiv->deriv2_classes[1][5][0] = int_stack + 11979;
 Libderiv->deriv2_classes[1][6][0] = int_stack + 12042;
 memset(int_stack,0,97008);

 Libderiv->dvrr_stack = int_stack + 28995;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0ff(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12126,int_stack+6039,int_stack+5949,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12216,int_stack+420,int_stack+6039,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+12351,int_stack+12216,int_stack+12126,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12531,int_stack+4209,int_stack+4149, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5949,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12621,int_stack+4299,int_stack+4209, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6039,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12756,int_stack+12621,int_stack+12531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12126,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12936,int_stack+0,int_stack+4299, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13125,int_stack+12936,int_stack+12621, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12216,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12936,int_stack+4569,int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5949, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13395,int_stack+4659,int_stack+4569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6039, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13530,int_stack+13395,int_stack+12936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12126, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13710,int_stack+84,int_stack+4659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+13899,int_stack+13710,int_stack+13395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12216, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13710,int_stack+4929,int_stack+4869, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5949, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+5019,int_stack+4929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6039, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14169,int_stack+0,int_stack+13710, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12126, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14349,int_stack+168,int_stack+5019, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+14538,int_stack+14349,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12216, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14349,int_stack+5289,int_stack+5229, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14808,int_stack+5379,int_stack+5289, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14943,int_stack+14808,int_stack+14349, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15123,int_stack+252,int_stack+5379, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15312,int_stack+15123,int_stack+14808, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15123,int_stack+5649,int_stack+5589, 0.0, zero_stack, 1.0, int_stack+5949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+135,int_stack+5739,int_stack+5649, 0.0, zero_stack, 1.0, int_stack+6039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15582,int_stack+135,int_stack+15123, 0.0, zero_stack, 1.0, int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15762,int_stack+336,int_stack+5739, 0.0, zero_stack, 1.0, int_stack+420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+15951,int_stack+15762,int_stack+135, 0.0, zero_stack, 1.0, int_stack+12216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15762,int_stack+6084,int_stack+5979, 1.0, int_stack+5949, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+6174,int_stack+6084, 1.0, int_stack+6039, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+16221,int_stack+270,int_stack+15762, 1.0, int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16401,int_stack+483,int_stack+6174, 1.0, int_stack+420, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16590,int_stack+16401,int_stack+270, 1.0, int_stack+12216, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+16401,int_stack+7554,int_stack+7494,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12126,int_stack+7644,int_stack+7554,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+16860,int_stack+12126,int_stack+16401,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+17040,int_stack+567,int_stack+7644,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17229,int_stack+17040,int_stack+12126,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12261,int_stack+9246,int_stack+9186,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17040,int_stack+9336,int_stack+9246,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+405,int_stack+17040,int_stack+12261,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+17499,int_stack+651,int_stack+9336,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+17688,int_stack+17499,int_stack+17040,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+17499,int_stack+11160,int_stack+11100,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+585,int_stack+11250,int_stack+11160,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17958,int_stack+585,int_stack+17499,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+18138,int_stack+735,int_stack+11250,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+18327,int_stack+18138,int_stack+585,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18138,int_stack+849,int_stack+819, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4149,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+894,int_stack+849, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4209,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18732,int_stack+18597,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12531,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18138,int_stack+957,int_stack+894, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4299,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+720,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12621,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18597,int_stack+1071,int_stack+1041, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4149, 1.0, int_stack+4509,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18138,int_stack+1116,int_stack+1071, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4209, 1.0, int_stack+4569,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18912,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12531, 1.0, int_stack+12936,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+19092,int_stack+1179,int_stack+1116, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4299, 1.0, int_stack+4659,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+990,int_stack+19092,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12621, 1.0, int_stack+13395,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18138,int_stack+1293,int_stack+1263, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4509, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+1338,int_stack+1293, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4569, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19092,int_stack+18597,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12936, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18138,int_stack+1401,int_stack+1338, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4659, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19272,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13395, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18597,int_stack+1515,int_stack+1485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4149, 0.0, zero_stack, 1.0, int_stack+4869,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18138,int_stack+1560,int_stack+1515, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4209, 0.0, zero_stack, 1.0, int_stack+4929,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19542,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12531, 0.0, zero_stack, 1.0, int_stack+13710,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+19722,int_stack+1623,int_stack+1560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4299, 0.0, zero_stack, 1.0, int_stack+5019,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1260,int_stack+19722,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12621, 0.0, zero_stack, 1.0, int_stack+0,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18138,int_stack+1737,int_stack+1707, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4509, 1.0, int_stack+4869, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+1782,int_stack+1737, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4569, 1.0, int_stack+4929, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19722,int_stack+18597,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12936, 1.0, int_stack+13710, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18138,int_stack+1845,int_stack+1782, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4659, 1.0, int_stack+5019, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1530,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13395, 1.0, int_stack+0, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18597,int_stack+1959,int_stack+1929, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4869, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18138,int_stack+2004,int_stack+1959, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4929, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+1800,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13710, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+19902,int_stack+2067,int_stack+2004, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5019, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20091,int_stack+19902,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18138,int_stack+2181,int_stack+2151, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4149, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5229,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+2226,int_stack+2181, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4209, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5289,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+19902,int_stack+18597,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12531, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14349,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18138,int_stack+2289,int_stack+2226, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4299, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5379,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20361,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12621, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14808,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18597,int_stack+2403,int_stack+2373, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4509, 0.0, zero_stack, 1.0, int_stack+5229, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18138,int_stack+2448,int_stack+2403, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4569, 0.0, zero_stack, 1.0, int_stack+5289, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20631,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12936, 0.0, zero_stack, 1.0, int_stack+14349, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+20811,int_stack+2511,int_stack+2448, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4659, 0.0, zero_stack, 1.0, int_stack+5379, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1980,int_stack+20811,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13395, 0.0, zero_stack, 1.0, int_stack+14808, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18138,int_stack+2625,int_stack+2595, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4869, 1.0, int_stack+5229, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+2670,int_stack+2625, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4929, 1.0, int_stack+5289, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+20811,int_stack+18597,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13710, 1.0, int_stack+14349, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18138,int_stack+2733,int_stack+2670, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5019, 1.0, int_stack+5379, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2250,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 1.0, int_stack+14808, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18597,int_stack+2847,int_stack+2817, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5229, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18138,int_stack+2892,int_stack+2847, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5289, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2520,int_stack+18138,int_stack+18597, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14349, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2700,int_stack+2955,int_stack+2892, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5379, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+20991,int_stack+2700,int_stack+18138, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14808, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18138,int_stack+3069,int_stack+3039, 0.0, zero_stack, 1.0, int_stack+4149, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5589,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+3114,int_stack+3069, 0.0, zero_stack, 1.0, int_stack+4209, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5649,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2700,int_stack+18597,int_stack+18138, 0.0, zero_stack, 1.0, int_stack+12531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15123,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18138,int_stack+3177,int_stack+3114, 0.0, zero_stack, 1.0, int_stack+4299, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5739,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2880,int_stack+18138,int_stack+18597, 0.0, zero_stack, 1.0, int_stack+12621, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18597,int_stack+3291,int_stack+3261, 0.0, zero_stack, 1.0, int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5589, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18138,int_stack+3336,int_stack+3291, 0.0, zero_stack, 1.0, int_stack+4569, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5649, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3150,int_stack+18138,int_stack+18597, 0.0, zero_stack, 1.0, int_stack+12936, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15123, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21261,int_stack+3399,int_stack+3336, 0.0, zero_stack, 1.0, int_stack+4659, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5739, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21450,int_stack+21261,int_stack+18138, 0.0, zero_stack, 1.0, int_stack+13395, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18138,int_stack+3513,int_stack+3483, 0.0, zero_stack, 1.0, int_stack+4869, 0.0, zero_stack, 1.0, int_stack+5589, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+3558,int_stack+3513, 0.0, zero_stack, 1.0, int_stack+4929, 0.0, zero_stack, 1.0, int_stack+5649, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21261,int_stack+18597,int_stack+18138, 0.0, zero_stack, 1.0, int_stack+13710, 0.0, zero_stack, 1.0, int_stack+15123, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18138,int_stack+3621,int_stack+3558, 0.0, zero_stack, 1.0, int_stack+5019, 0.0, zero_stack, 1.0, int_stack+5739, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3330,int_stack+18138,int_stack+18597, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18597,int_stack+3735,int_stack+3705, 0.0, zero_stack, 1.0, int_stack+5229, 1.0, int_stack+5589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3600,int_stack+3780,int_stack+3735, 0.0, zero_stack, 1.0, int_stack+5289, 1.0, int_stack+5649, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+18138,int_stack+3600,int_stack+18597, 0.0, zero_stack, 1.0, int_stack+14349, 1.0, int_stack+15123, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21720,int_stack+3843,int_stack+3780, 0.0, zero_stack, 1.0, int_stack+5379, 1.0, int_stack+5739, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+21909,int_stack+21720,int_stack+3600, 0.0, zero_stack, 1.0, int_stack+14808, 1.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+3957,int_stack+3927, 0.0, zero_stack, 2.0, int_stack+5589, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+4002,int_stack+3957, 0.0, zero_stack, 2.0, int_stack+5649, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3690,int_stack+18597,int_stack+3600, 0.0, zero_stack, 2.0, int_stack+15123, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+21720,int_stack+4065,int_stack+4002, 0.0, zero_stack, 2.0, int_stack+5739, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3870,int_stack+21720,int_stack+18597, 0.0, zero_stack, 2.0, int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4254,int_stack+4179, 1.0, int_stack+4149, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5979,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+4362,int_stack+4254, 1.0, int_stack+4209, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6084,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+21720,int_stack+18597,int_stack+3600, 1.0, int_stack+12531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15762,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22179,int_stack+4425,int_stack+4362, 1.0, int_stack+4299, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6174,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4140,int_stack+22179,int_stack+18597, 1.0, int_stack+12621, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4614,int_stack+4539, 1.0, int_stack+4509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5979, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+18597,int_stack+4722,int_stack+4614, 1.0, int_stack+4569, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6084, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22179,int_stack+18597,int_stack+3600, 1.0, int_stack+12936, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15762, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12936,int_stack+4785,int_stack+4722, 1.0, int_stack+4659, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6174, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22359,int_stack+12936,int_stack+18597, 1.0, int_stack+13395, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4974,int_stack+4899, 1.0, int_stack+4869, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5979, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13395,int_stack+5082,int_stack+4974, 1.0, int_stack+4929, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6084, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12936,int_stack+13395,int_stack+3600, 1.0, int_stack+13710, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15762, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13710,int_stack+5145,int_stack+5082, 1.0, int_stack+5019, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6174, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+22629,int_stack+13710,int_stack+13395, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+5334,int_stack+5259, 1.0, int_stack+5229, 0.0, zero_stack, 1.0, int_stack+5979, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+5442,int_stack+5334, 1.0, int_stack+5289, 0.0, zero_stack, 1.0, int_stack+6084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+13710,int_stack+0,int_stack+3600, 1.0, int_stack+14349, 0.0, zero_stack, 1.0, int_stack+15762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14349,int_stack+5505,int_stack+5442, 1.0, int_stack+5379, 0.0, zero_stack, 1.0, int_stack+6174, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4410,int_stack+14349,int_stack+0, 1.0, int_stack+14808, 0.0, zero_stack, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+5694,int_stack+5619, 1.0, int_stack+5589, 1.0, int_stack+5979, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14808,int_stack+5802,int_stack+5694, 1.0, int_stack+5649, 1.0, int_stack+6084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+14349,int_stack+14808,int_stack+3600, 1.0, int_stack+15123, 1.0, int_stack+15762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15123,int_stack+5865,int_stack+5802, 1.0, int_stack+5739, 1.0, int_stack+6174, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4680,int_stack+15123,int_stack+14808, 1.0, int_stack+135, 1.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+6129,int_stack+6009, 2.0, int_stack+5979, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14808,int_stack+6237,int_stack+6129, 2.0, int_stack+6084, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15123,int_stack+14808,int_stack+3600, 2.0, int_stack+15762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+15762,int_stack+6300,int_stack+6237, 2.0, int_stack+6174, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+15762,int_stack+14808, 2.0, int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+6414,int_stack+6384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7494,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+6459,int_stack+6414, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7554,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+15762,int_stack+270,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16401,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12531,int_stack+6522,int_stack+6459, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7644,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4950,int_stack+12531,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12126,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+6636,int_stack+6606, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7494, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+6681,int_stack+6636, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7554, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12531,int_stack+270,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16401, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5220,int_stack+6744,int_stack+6681, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7644, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5409,int_stack+5220,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12126, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+6858,int_stack+6828, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7494, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+6903,int_stack+6858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7554, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5220,int_stack+270,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16401, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5679,int_stack+6966,int_stack+6903, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7644, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5868,int_stack+5679,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12126, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+7080,int_stack+7050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+7125,int_stack+7080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5679,int_stack+270,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6138,int_stack+7188,int_stack+7125, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6327,int_stack+6138,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+7302,int_stack+7272, 0.0, zero_stack, 1.0, int_stack+7494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+7347,int_stack+7302, 0.0, zero_stack, 1.0, int_stack+7554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6138,int_stack+270,int_stack+3600, 0.0, zero_stack, 1.0, int_stack+16401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+6597,int_stack+7410,int_stack+7347, 0.0, zero_stack, 1.0, int_stack+7644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6786,int_stack+6597,int_stack+270, 0.0, zero_stack, 1.0, int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+7599,int_stack+7524, 1.0, int_stack+7494, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+7707,int_stack+7599, 1.0, int_stack+7554, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+6597,int_stack+270,int_stack+3600, 1.0, int_stack+16401, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+16401,int_stack+7770,int_stack+7707, 1.0, int_stack+7644, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7056,int_stack+16401,int_stack+270, 1.0, int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+7884,int_stack+7854,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+12126,int_stack+7929,int_stack+7884,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+16401,int_stack+12126,int_stack+3600,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7326,int_stack+7992,int_stack+7929,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+7515,int_stack+7326,int_stack+12126,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+8106,int_stack+8076, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9186,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12126,int_stack+8151,int_stack+8106, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9246,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7326,int_stack+12126,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12261,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7785,int_stack+8214,int_stack+8151, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9336,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7974,int_stack+7785,int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17040,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+8328,int_stack+8298, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9186, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12126,int_stack+8373,int_stack+8328, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9246, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7785,int_stack+12126,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12261, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+22899,int_stack+8436,int_stack+8373, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9336, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23088,int_stack+22899,int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17040, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+8550,int_stack+8520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9186, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12126,int_stack+8595,int_stack+8550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9246, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+22899,int_stack+12126,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12261, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23358,int_stack+8658,int_stack+8595, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9336, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+23547,int_stack+23358,int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17040, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+8772,int_stack+8742, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12126,int_stack+8817,int_stack+8772, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+23358,int_stack+12126,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8244,int_stack+8880,int_stack+8817, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8433,int_stack+8244,int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+8994,int_stack+8964, 0.0, zero_stack, 1.0, int_stack+9186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12126,int_stack+9039,int_stack+8994, 0.0, zero_stack, 1.0, int_stack+9246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8244,int_stack+12126,int_stack+3600, 0.0, zero_stack, 1.0, int_stack+12261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8703,int_stack+9102,int_stack+9039, 0.0, zero_stack, 1.0, int_stack+9336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+8892,int_stack+8703,int_stack+12126, 0.0, zero_stack, 1.0, int_stack+17040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+9291,int_stack+9216, 1.0, int_stack+9186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12126,int_stack+9399,int_stack+9291, 1.0, int_stack+9246, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8703,int_stack+12126,int_stack+3600, 1.0, int_stack+12261, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+23817,int_stack+9462,int_stack+9399, 1.0, int_stack+9336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24006,int_stack+23817,int_stack+12126, 1.0, int_stack+17040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+9576,int_stack+9546,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+9621,int_stack+9576,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17040,int_stack+270,int_stack+3600,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+23817,int_stack+9684,int_stack+9621,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+24276,int_stack+23817,int_stack+270,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+9798,int_stack+9768,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+9843,int_stack+9798,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+23817,int_stack+270,int_stack+3600,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12126,int_stack+9906,int_stack+9843,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+24546,int_stack+12126,int_stack+270,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+10020,int_stack+9990, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11100,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+10065,int_stack+10020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11160,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+12126,int_stack+270,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17499,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9162,int_stack+10128,int_stack+10065, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11250,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9351,int_stack+9162,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+10242,int_stack+10212, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11100, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+10287,int_stack+10242, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11160, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9162,int_stack+270,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17499, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9621,int_stack+10350,int_stack+10287, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11250, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9810,int_stack+9621,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+10464,int_stack+10434, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11100, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+10509,int_stack+10464, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11160, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+9621,int_stack+270,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17499, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10080,int_stack+10572,int_stack+10509, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11250, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10269,int_stack+10080,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+10686,int_stack+10656, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+10731,int_stack+10686, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10080,int_stack+270,int_stack+3600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17499, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10539,int_stack+10794,int_stack+10731, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+24816,int_stack+10539,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+10908,int_stack+10878, 0.0, zero_stack, 1.0, int_stack+11100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+10953,int_stack+10908, 0.0, zero_stack, 1.0, int_stack+11160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10539,int_stack+270,int_stack+3600, 0.0, zero_stack, 1.0, int_stack+17499, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10719,int_stack+11016,int_stack+10953, 0.0, zero_stack, 1.0, int_stack+11250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+25086,int_stack+10719,int_stack+270, 0.0, zero_stack, 1.0, int_stack+585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+11205,int_stack+11130, 1.0, int_stack+11100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+11313,int_stack+11205, 1.0, int_stack+11160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+10719,int_stack+270,int_stack+3600, 1.0, int_stack+17499, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17499,int_stack+11376,int_stack+11313, 1.0, int_stack+11250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10899,int_stack+17499,int_stack+270, 1.0, int_stack+585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+11490,int_stack+11460,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+585,int_stack+11535,int_stack+11490,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+17499,int_stack+585,int_stack+3600,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11169,int_stack+11598,int_stack+11535,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11358,int_stack+11169,int_stack+585,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+11712,int_stack+11682,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+585,int_stack+11757,int_stack+11712,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+11169,int_stack+585,int_stack+3600,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+25356,int_stack+11820,int_stack+11757,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+11628,int_stack+25356,int_stack+585,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+11934,int_stack+11904,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+585,int_stack+11979,int_stack+11934,3);
 /*--- compute (p0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+25356,int_stack+585,int_stack+3600,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+25536,int_stack+12042,int_stack+11979,3);
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+25725,int_stack+25536,int_stack+585,3);
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+25995,int_stack+13125,int_stack+12756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12351,3);
     Libderiv->ABCD[11] = int_stack + 25995;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+26295,int_stack+13899,int_stack+13530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12351, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 26295;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+26595,int_stack+14538,int_stack+14169, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12351, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 26595;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+14529,int_stack+15312,int_stack+14943, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 14529;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+13116,int_stack+15951,int_stack+15582, 0.0, zero_stack, 1.0, int_stack+12351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 13116;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+26895,int_stack+16590,int_stack+16221, 1.0, int_stack+12351, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 26895;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+27195,int_stack+17229,int_stack+16860,3);
     Libderiv->ABCD[2] = int_stack + 27195;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+27495,int_stack+17688,int_stack+405,3);
     Libderiv->ABCD[1] = int_stack + 27495;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+27795,int_stack+18327,int_stack+17958,3);
     Libderiv->ABCD[0] = int_stack + 27795;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+28095,int_stack+720,int_stack+18732, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12756,3);
     Libderiv->ABCD[155] = int_stack + 28095;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+585,int_stack+990,int_stack+18912, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12756, 1.0, int_stack+13530,3);
     Libderiv->ABCD[143] = int_stack + 585;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+885,int_stack+19272,int_stack+19092, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13530, 0.0, zero_stack,3);
     Libderiv->ABCD[142] = int_stack + 885;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+28395,int_stack+1260,int_stack+19542, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12756, 0.0, zero_stack, 1.0, int_stack+14169,3);
     Libderiv->ABCD[131] = int_stack + 28395;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1185,int_stack+1530,int_stack+19722, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13530, 1.0, int_stack+14169, 0.0, zero_stack,3);
     Libderiv->ABCD[130] = int_stack + 1185;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1485,int_stack+20091,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14169, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[129] = int_stack + 1485;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+28695,int_stack+20361,int_stack+19902, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12756, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14943,3);
     Libderiv->ABCD[119] = int_stack + 28695;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18318,int_stack+1980,int_stack+20631, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13530, 0.0, zero_stack, 1.0, int_stack+14943, 0.0, zero_stack,3);
     Libderiv->ABCD[118] = int_stack + 18318;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+1785,int_stack+2250,int_stack+20811, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14169, 1.0, int_stack+14943, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[117] = int_stack + 1785;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2085,int_stack+20991,int_stack+2520, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14943, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[116] = int_stack + 2085;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2385,int_stack+2880,int_stack+2700, 0.0, zero_stack, 1.0, int_stack+12756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15582,3);
     Libderiv->ABCD[107] = int_stack + 2385;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2685,int_stack+21450,int_stack+3150, 0.0, zero_stack, 1.0, int_stack+13530, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15582, 0.0, zero_stack,3);
     Libderiv->ABCD[106] = int_stack + 2685;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+2985,int_stack+3330,int_stack+21261, 0.0, zero_stack, 1.0, int_stack+14169, 0.0, zero_stack, 1.0, int_stack+15582, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[105] = int_stack + 2985;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3285,int_stack+21909,int_stack+18138, 0.0, zero_stack, 1.0, int_stack+14943, 1.0, int_stack+15582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[104] = int_stack + 3285;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18618,int_stack+3870,int_stack+3690, 0.0, zero_stack, 2.0, int_stack+15582, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[103] = int_stack + 18618;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3585,int_stack+4140,int_stack+21720, 1.0, int_stack+12756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16221,3);
     Libderiv->ABCD[95] = int_stack + 3585;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+3885,int_stack+22359,int_stack+22179, 1.0, int_stack+13530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16221, 0.0, zero_stack,3);
     Libderiv->ABCD[94] = int_stack + 3885;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+18918,int_stack+22629,int_stack+12936, 1.0, int_stack+14169, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16221, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[93] = int_stack + 18918;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12711,int_stack+4410,int_stack+13710, 1.0, int_stack+14943, 0.0, zero_stack, 1.0, int_stack+16221, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[92] = int_stack + 12711;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4185,int_stack+4680,int_stack+14349, 1.0, int_stack+15582, 1.0, int_stack+16221, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[91] = int_stack + 4185;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4485,int_stack+0,int_stack+15123, 2.0, int_stack+16221, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[90] = int_stack + 4485;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+0,int_stack+4950,int_stack+15762, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16860,3);
     Libderiv->ABCD[47] = int_stack + 0;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+4785,int_stack+5409,int_stack+12531, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16860, 0.0, zero_stack,3);
     Libderiv->ABCD[46] = int_stack + 4785;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+12306,int_stack+5868,int_stack+5220, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16860, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[45] = int_stack + 12306;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5085,int_stack+6327,int_stack+5679, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[44] = int_stack + 5085;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5385,int_stack+6786,int_stack+6138, 0.0, zero_stack, 1.0, int_stack+16860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[43] = int_stack + 5385;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+5685,int_stack+7056,int_stack+6597, 1.0, int_stack+16860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[42] = int_stack + 5685;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+5985,int_stack+7515,int_stack+16401,3);
     Libderiv->ABCD[38] = int_stack + 5985;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6285,int_stack+7974,int_stack+7326, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+405,3);
     Libderiv->ABCD[35] = int_stack + 6285;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6585,int_stack+23088,int_stack+7785, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+405, 0.0, zero_stack,3);
     Libderiv->ABCD[34] = int_stack + 6585;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+6885,int_stack+23547,int_stack+22899, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+405, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[33] = int_stack + 6885;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7185,int_stack+8433,int_stack+23358, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[32] = int_stack + 7185;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7485,int_stack+8892,int_stack+8244, 0.0, zero_stack, 1.0, int_stack+405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[31] = int_stack + 7485;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+7785,int_stack+24006,int_stack+8703, 1.0, int_stack+405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[30] = int_stack + 7785;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+8085,int_stack+24276,int_stack+17040,3);
     Libderiv->ABCD[26] = int_stack + 8085;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+8385,int_stack+24546,int_stack+23817,3);
     Libderiv->ABCD[25] = int_stack + 8385;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8685,int_stack+9351,int_stack+12126, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17958,3);
     Libderiv->ABCD[23] = int_stack + 8685;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+11898,int_stack+9810,int_stack+9162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17958, 0.0, zero_stack,3);
     Libderiv->ABCD[22] = int_stack + 11898;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+8985,int_stack+10269,int_stack+9621, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17958, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[21] = int_stack + 8985;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9285,int_stack+24816,int_stack+10080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+17958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[20] = int_stack + 9285;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9585,int_stack+25086,int_stack+10539, 0.0, zero_stack, 1.0, int_stack+17958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[19] = int_stack + 9585;
 /*--- compute (p0|ff) ---*/
   d1hrr3_build_ff(Libderiv->CD,int_stack+9885,int_stack+10899,int_stack+10719, 1.0, int_stack+17958, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[18] = int_stack + 9885;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+10185,int_stack+11358,int_stack+17499,3);
     Libderiv->ABCD[14] = int_stack + 10185;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+10485,int_stack+11628,int_stack+11169,3);
     Libderiv->ABCD[13] = int_stack + 10485;
 /*--- compute (p0|ff) ---*/
   hrr3_build_ff(Libderiv->CD,int_stack+10785,int_stack+25725,int_stack+25356,3);
     Libderiv->ABCD[12] = int_stack + 10785;

}
