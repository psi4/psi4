#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_fpfp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|fp) integrals */

void d12hrr_order_fpfp(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv2_classes[3][3][143] = int_stack + 2325;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 2425;
 Libderiv->deriv2_classes[4][3][143] = int_stack + 2575;
 Libderiv->deriv2_classes[4][4][143] = int_stack + 2725;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 2950;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 3050;
 Libderiv->deriv2_classes[4][3][131] = int_stack + 3200;
 Libderiv->deriv2_classes[4][4][131] = int_stack + 3350;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 3575;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 3675;
 Libderiv->deriv2_classes[4][3][130] = int_stack + 3825;
 Libderiv->deriv2_classes[4][4][130] = int_stack + 3975;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 4200;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 4300;
 Libderiv->deriv2_classes[4][3][119] = int_stack + 4450;
 Libderiv->deriv2_classes[4][4][119] = int_stack + 4600;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 4825;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 4925;
 Libderiv->deriv2_classes[4][3][118] = int_stack + 5075;
 Libderiv->deriv2_classes[4][4][118] = int_stack + 5225;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 5450;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 5550;
 Libderiv->deriv2_classes[4][3][117] = int_stack + 5700;
 Libderiv->deriv2_classes[4][4][117] = int_stack + 5850;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 6075;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 6175;
 Libderiv->deriv2_classes[4][3][107] = int_stack + 6325;
 Libderiv->deriv2_classes[4][4][107] = int_stack + 6475;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 6700;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 6800;
 Libderiv->deriv2_classes[4][3][106] = int_stack + 6950;
 Libderiv->deriv2_classes[4][4][106] = int_stack + 7100;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 7325;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 7425;
 Libderiv->deriv2_classes[4][3][105] = int_stack + 7575;
 Libderiv->deriv2_classes[4][4][105] = int_stack + 7725;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 7950;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 8050;
 Libderiv->deriv2_classes[4][3][104] = int_stack + 8200;
 Libderiv->deriv2_classes[4][4][104] = int_stack + 8350;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 8575;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 8675;
 Libderiv->deriv2_classes[4][3][95] = int_stack + 8825;
 Libderiv->deriv2_classes[4][4][95] = int_stack + 8975;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 9200;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 9300;
 Libderiv->deriv2_classes[4][3][94] = int_stack + 9450;
 Libderiv->deriv2_classes[4][4][94] = int_stack + 9600;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 9825;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 9925;
 Libderiv->deriv2_classes[4][3][93] = int_stack + 10075;
 Libderiv->deriv2_classes[4][4][93] = int_stack + 10225;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 10450;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 10550;
 Libderiv->deriv2_classes[4][3][92] = int_stack + 10700;
 Libderiv->deriv2_classes[4][4][92] = int_stack + 10850;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 11075;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 11175;
 Libderiv->deriv2_classes[4][3][91] = int_stack + 11325;
 Libderiv->deriv2_classes[4][4][91] = int_stack + 11475;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 11700;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 11800;
 Libderiv->deriv_classes[4][3][11] = int_stack + 11950;
 Libderiv->deriv2_classes[4][3][83] = int_stack + 12100;
 Libderiv->deriv2_classes[4][4][83] = int_stack + 12250;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 12475;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 12575;
 Libderiv->deriv_classes[4][3][10] = int_stack + 12725;
 Libderiv->deriv2_classes[4][3][82] = int_stack + 12875;
 Libderiv->deriv2_classes[4][4][82] = int_stack + 13025;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 13250;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 13350;
 Libderiv->deriv_classes[4][3][9] = int_stack + 13500;
 Libderiv->deriv2_classes[4][3][81] = int_stack + 13650;
 Libderiv->deriv2_classes[4][4][81] = int_stack + 13800;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 14025;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 14125;
 Libderiv->deriv_classes[4][3][8] = int_stack + 14275;
 Libderiv->deriv2_classes[4][3][80] = int_stack + 14425;
 Libderiv->deriv2_classes[4][4][80] = int_stack + 14575;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 14800;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 14900;
 Libderiv->deriv_classes[4][3][7] = int_stack + 15050;
 Libderiv->deriv2_classes[4][3][79] = int_stack + 15200;
 Libderiv->deriv2_classes[4][4][79] = int_stack + 15350;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 15575;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 15675;
 Libderiv->deriv_classes[4][3][6] = int_stack + 15825;
 Libderiv->deriv2_classes[4][3][78] = int_stack + 15975;
 Libderiv->deriv2_classes[4][4][78] = int_stack + 16125;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 16350;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 16450;
 Libderiv->deriv2_classes[4][3][35] = int_stack + 16600;
 Libderiv->deriv2_classes[4][4][35] = int_stack + 16750;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 16975;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 17075;
 Libderiv->deriv2_classes[4][3][34] = int_stack + 17225;
 Libderiv->deriv2_classes[4][4][34] = int_stack + 17375;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 17600;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 17700;
 Libderiv->deriv2_classes[4][3][33] = int_stack + 17850;
 Libderiv->deriv2_classes[4][4][33] = int_stack + 18000;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 18225;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 18325;
 Libderiv->deriv2_classes[4][3][32] = int_stack + 18475;
 Libderiv->deriv2_classes[4][4][32] = int_stack + 18625;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 18850;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 18950;
 Libderiv->deriv2_classes[4][3][31] = int_stack + 19100;
 Libderiv->deriv2_classes[4][4][31] = int_stack + 19250;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 19475;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 19575;
 Libderiv->deriv_classes[4][3][2] = int_stack + 19725;
 Libderiv->deriv2_classes[4][3][30] = int_stack + 19875;
 Libderiv->deriv2_classes[4][4][30] = int_stack + 20025;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 20250;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 20350;
 Libderiv->deriv2_classes[4][3][26] = int_stack + 20500;
 Libderiv->deriv2_classes[4][4][26] = int_stack + 20650;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 20875;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 20975;
 Libderiv->deriv2_classes[4][3][23] = int_stack + 21125;
 Libderiv->deriv2_classes[4][4][23] = int_stack + 21275;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 21500;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 21600;
 Libderiv->deriv2_classes[4][3][22] = int_stack + 21750;
 Libderiv->deriv2_classes[4][4][22] = int_stack + 21900;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 22125;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 22225;
 Libderiv->deriv2_classes[4][3][21] = int_stack + 22375;
 Libderiv->deriv2_classes[4][4][21] = int_stack + 22525;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 22750;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 22850;
 Libderiv->deriv2_classes[4][3][20] = int_stack + 23000;
 Libderiv->deriv2_classes[4][4][20] = int_stack + 23150;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 23375;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 23475;
 Libderiv->deriv2_classes[4][3][19] = int_stack + 23625;
 Libderiv->deriv2_classes[4][4][19] = int_stack + 23775;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 24000;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 24100;
 Libderiv->deriv_classes[4][3][1] = int_stack + 24250;
 Libderiv->deriv2_classes[4][3][18] = int_stack + 24400;
 Libderiv->deriv2_classes[4][4][18] = int_stack + 24550;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 24775;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 24875;
 Libderiv->deriv2_classes[4][3][14] = int_stack + 25025;
 Libderiv->deriv2_classes[4][4][14] = int_stack + 25175;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 25400;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 25500;
 Libderiv->deriv2_classes[4][3][13] = int_stack + 25650;
 Libderiv->deriv2_classes[4][4][13] = int_stack + 25800;
 Libderiv->deriv_classes[3][3][11] = int_stack + 26025;
 Libderiv->deriv_classes[3][4][11] = int_stack + 26125;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 26275;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 26375;
 Libderiv->deriv2_classes[4][3][11] = int_stack + 26525;
 Libderiv->deriv2_classes[4][4][11] = int_stack + 26675;
 Libderiv->deriv_classes[3][3][10] = int_stack + 26900;
 Libderiv->deriv_classes[3][4][10] = int_stack + 27000;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 27150;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 27250;
 Libderiv->deriv2_classes[4][3][10] = int_stack + 27400;
 Libderiv->deriv2_classes[4][4][10] = int_stack + 27550;
 Libderiv->deriv_classes[3][3][9] = int_stack + 27775;
 Libderiv->deriv_classes[3][4][9] = int_stack + 27875;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 28025;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 28125;
 Libderiv->deriv2_classes[4][3][9] = int_stack + 28275;
 Libderiv->deriv2_classes[4][4][9] = int_stack + 28425;
 Libderiv->deriv_classes[3][3][8] = int_stack + 28650;
 Libderiv->deriv_classes[3][4][8] = int_stack + 28750;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 28900;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 29000;
 Libderiv->deriv2_classes[4][3][8] = int_stack + 29150;
 Libderiv->deriv2_classes[4][4][8] = int_stack + 29300;
 Libderiv->deriv_classes[3][3][7] = int_stack + 29525;
 Libderiv->deriv_classes[3][4][7] = int_stack + 29625;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 29775;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 29875;
 Libderiv->deriv2_classes[4][3][7] = int_stack + 30025;
 Libderiv->deriv2_classes[4][4][7] = int_stack + 30175;
 Libderiv->dvrr_classes[3][3] = int_stack + 30400;
 Libderiv->deriv_classes[3][3][6] = int_stack + 30500;
 Libderiv->deriv_classes[3][4][6] = int_stack + 30600;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 30750;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 30850;
 Libderiv->deriv_classes[4][3][0] = int_stack + 31000;
 Libderiv->deriv2_classes[4][3][6] = int_stack + 31150;
 Libderiv->deriv2_classes[4][4][6] = int_stack + 31300;
 Libderiv->deriv_classes[3][3][2] = int_stack + 31525;
 Libderiv->deriv_classes[3][4][2] = int_stack + 31625;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 31775;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 31875;
 Libderiv->deriv2_classes[4][3][2] = int_stack + 32025;
 Libderiv->deriv2_classes[4][4][2] = int_stack + 32175;
 Libderiv->deriv_classes[3][3][1] = int_stack + 32400;
 Libderiv->deriv_classes[3][4][1] = int_stack + 32500;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 32650;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 32750;
 Libderiv->deriv2_classes[4][3][1] = int_stack + 32900;
 Libderiv->deriv2_classes[4][4][1] = int_stack + 33050;
 Libderiv->deriv_classes[3][3][0] = int_stack + 33275;
 Libderiv->deriv_classes[3][4][0] = int_stack + 33375;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 33525;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 33625;
 Libderiv->deriv2_classes[4][3][0] = int_stack + 33775;
 Libderiv->deriv2_classes[4][4][0] = int_stack + 33925;
 memset(int_stack,0,273200);

 Libderiv->dvrr_stack = int_stack + 53350;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_fpfp(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+34150,int_stack+26125,int_stack+26025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30400,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+34450,int_stack+0,int_stack+11950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+34900,int_stack+27000,int_stack+26900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30400, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+35200,int_stack+225,int_stack+12725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+27875,int_stack+27775, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30400, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+35650,int_stack+450,int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+28750,int_stack+28650, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36100,int_stack+675,int_stack+14275, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+29625,int_stack+29525, 0.0, zero_stack, 1.0, int_stack+30400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36550,int_stack+900,int_stack+15050, 0.0, zero_stack, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37000,int_stack+30600,int_stack+30500, 1.0, int_stack+30400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+37300,int_stack+1275,int_stack+15825, 1.0, int_stack+1125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1950,int_stack+30400,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+31625,int_stack+31525,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+37750,int_stack+1500,int_stack+19725,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+38200,int_stack+32500,int_stack+32400,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+38500,int_stack+1725,int_stack+24250,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1500,int_stack+33375,int_stack+33275,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+38950,int_stack+2100,int_stack+31000,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2425,int_stack+2325, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+26025,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+2725,int_stack+2575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+11950,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2550,int_stack+3050,int_stack+2950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26025, 1.0, int_stack+26900,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+39400,int_stack+3350,int_stack+3200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11950, 1.0, int_stack+12725,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2850,int_stack+3675,int_stack+3575, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+26900, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3150,int_stack+3975,int_stack+3825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+12725, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4300,int_stack+4200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26025, 0.0, zero_stack, 1.0, int_stack+27775,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3900,int_stack+4600,int_stack+4450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11950, 0.0, zero_stack, 1.0, int_stack+13500,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4350,int_stack+4925,int_stack+4825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26900, 1.0, int_stack+27775, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+39850,int_stack+5225,int_stack+5075, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12725, 1.0, int_stack+13500, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4650,int_stack+5550,int_stack+5450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+27775, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4950,int_stack+5850,int_stack+5700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5400,int_stack+6175,int_stack+6075, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26025, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+28650,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5700,int_stack+6475,int_stack+6325, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+11950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14275,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6150,int_stack+6800,int_stack+6700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+26900, 0.0, zero_stack, 1.0, int_stack+28650, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6450,int_stack+7100,int_stack+6950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12725, 0.0, zero_stack, 1.0, int_stack+14275, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6900,int_stack+7425,int_stack+7325, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+27775, 1.0, int_stack+28650, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40300,int_stack+7725,int_stack+7575, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+13500, 1.0, int_stack+14275, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+8050,int_stack+7950, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+28650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7500,int_stack+8350,int_stack+8200, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+14275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7950,int_stack+8675,int_stack+8575, 0.0, zero_stack, 1.0, int_stack+26025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29525,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8250,int_stack+8975,int_stack+8825, 0.0, zero_stack, 1.0, int_stack+11950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15050,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8700,int_stack+9300,int_stack+9200, 0.0, zero_stack, 1.0, int_stack+26900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+29525, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9000,int_stack+9600,int_stack+9450, 0.0, zero_stack, 1.0, int_stack+12725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15050, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9450,int_stack+9925,int_stack+9825, 0.0, zero_stack, 1.0, int_stack+27775, 0.0, zero_stack, 1.0, int_stack+29525, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+40750,int_stack+10225,int_stack+10075, 0.0, zero_stack, 1.0, int_stack+13500, 0.0, zero_stack, 1.0, int_stack+15050, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9750,int_stack+10550,int_stack+10450, 0.0, zero_stack, 1.0, int_stack+28650, 1.0, int_stack+29525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10050,int_stack+10850,int_stack+10700, 0.0, zero_stack, 1.0, int_stack+14275, 1.0, int_stack+15050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10500,int_stack+11175,int_stack+11075, 0.0, zero_stack, 2.0, int_stack+29525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10800,int_stack+11475,int_stack+11325, 0.0, zero_stack, 2.0, int_stack+15050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11250,int_stack+11800,int_stack+11700, 1.0, int_stack+26025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30500,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41200,int_stack+12250,int_stack+12100, 1.0, int_stack+11950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15825,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11550,int_stack+12575,int_stack+12475, 1.0, int_stack+26900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30500, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11850,int_stack+13025,int_stack+12875, 1.0, int_stack+12725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15825, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12300,int_stack+13350,int_stack+13250, 1.0, int_stack+27775, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+30500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12600,int_stack+13800,int_stack+13650, 1.0, int_stack+13500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15825, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13050,int_stack+14125,int_stack+14025, 1.0, int_stack+28650, 0.0, zero_stack, 1.0, int_stack+30500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13350,int_stack+14575,int_stack+14425, 1.0, int_stack+14275, 0.0, zero_stack, 1.0, int_stack+15825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13800,int_stack+14900,int_stack+14800, 1.0, int_stack+29525, 1.0, int_stack+30500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14100,int_stack+15350,int_stack+15200, 1.0, int_stack+15050, 1.0, int_stack+15825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14550,int_stack+15675,int_stack+15575, 2.0, int_stack+30500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14850,int_stack+16125,int_stack+15975, 2.0, int_stack+15825, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15300,int_stack+16450,int_stack+16350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31525,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15600,int_stack+16750,int_stack+16600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19725,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16050,int_stack+17075,int_stack+16975, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31525, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16350,int_stack+17375,int_stack+17225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19725, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16800,int_stack+17700,int_stack+17600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31525, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17100,int_stack+18000,int_stack+17850, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19725, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17550,int_stack+18325,int_stack+18225, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+17850,int_stack+18625,int_stack+18475, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+19725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18300,int_stack+18950,int_stack+18850, 0.0, zero_stack, 1.0, int_stack+31525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+18600,int_stack+19250,int_stack+19100, 0.0, zero_stack, 1.0, int_stack+19725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+19050,int_stack+19575,int_stack+19475, 1.0, int_stack+31525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+41650,int_stack+20025,int_stack+19875, 1.0, int_stack+19725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19350,int_stack+20350,int_stack+20250,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+19650,int_stack+20650,int_stack+20500,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20100,int_stack+20975,int_stack+20875, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32400,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20400,int_stack+21275,int_stack+21125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24250,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+20850,int_stack+21600,int_stack+21500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32400, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21150,int_stack+21900,int_stack+21750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24250, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21600,int_stack+22225,int_stack+22125, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32400, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+21900,int_stack+22525,int_stack+22375, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24250, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22350,int_stack+22850,int_stack+22750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+32400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+42100,int_stack+23150,int_stack+23000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+24250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22650,int_stack+23475,int_stack+23375, 0.0, zero_stack, 1.0, int_stack+32400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+22950,int_stack+23775,int_stack+23625, 0.0, zero_stack, 1.0, int_stack+24250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23400,int_stack+24100,int_stack+24000, 1.0, int_stack+32400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+23700,int_stack+24550,int_stack+24400, 1.0, int_stack+24250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24150,int_stack+24875,int_stack+24775,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24450,int_stack+25175,int_stack+25025,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+24900,int_stack+25500,int_stack+25400,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+25200,int_stack+25800,int_stack+25650,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25650,int_stack+26375,int_stack+26275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33275,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+25950,int_stack+26675,int_stack+26525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31000,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26400,int_stack+27250,int_stack+27150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33275, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+26700,int_stack+27550,int_stack+27400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31000, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27150,int_stack+28125,int_stack+28025, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33275, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27450,int_stack+28425,int_stack+28275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31000, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+27900,int_stack+29000,int_stack+28900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+33275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28200,int_stack+29300,int_stack+29150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+31000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28650,int_stack+29875,int_stack+29775, 0.0, zero_stack, 1.0, int_stack+33275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+28950,int_stack+30175,int_stack+30025, 0.0, zero_stack, 1.0, int_stack+31000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29400,int_stack+30850,int_stack+30750, 1.0, int_stack+33275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+29700,int_stack+31300,int_stack+31150, 1.0, int_stack+31000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+30150,int_stack+31875,int_stack+31775,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+30450,int_stack+32175,int_stack+32025,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+30900,int_stack+32750,int_stack+32650,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+31200,int_stack+33050,int_stack+32900,15);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+31650,int_stack+33625,int_stack+33525,10);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+31950,int_stack+33925,int_stack+33775,15);
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+32400,int_stack+34450,int_stack+34150,30);
     Libderiv->ABCD[11] = int_stack + 32400;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+42550,int_stack+35200,int_stack+34900,30);
     Libderiv->ABCD[10] = int_stack + 42550;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+43450,int_stack+35650,int_stack+0,30);
     Libderiv->ABCD[9] = int_stack + 43450;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+35200,int_stack+36100,int_stack+300,30);
     Libderiv->ABCD[8] = int_stack + 35200;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+44350,int_stack+36550,int_stack+600,30);
     Libderiv->ABCD[7] = int_stack + 44350;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+36100,int_stack+37300,int_stack+37000,30);
     Libderiv->ABCD[6] = int_stack + 36100;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+45250,int_stack+37750,int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[2] = int_stack + 45250;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+37300,int_stack+38500,int_stack+38200, 0.0, zero_stack, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[1] = int_stack + 37300;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+46150,int_stack+38950,int_stack+1500, 1.0, int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[0] = int_stack + 46150;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+38500,int_stack+2100,int_stack+1800,30);
     Libderiv->ABCD[155] = int_stack + 38500;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+47050,int_stack+39400,int_stack+2550,30);
     Libderiv->ABCD[143] = int_stack + 47050;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1800,int_stack+3150,int_stack+2850,30);
     Libderiv->ABCD[142] = int_stack + 1800;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+2700,int_stack+3900,int_stack+3600,30);
     Libderiv->ABCD[131] = int_stack + 2700;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+47950,int_stack+39850,int_stack+4350,30);
     Libderiv->ABCD[130] = int_stack + 47950;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+39400,int_stack+4950,int_stack+4650,30);
     Libderiv->ABCD[129] = int_stack + 39400;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+3600,int_stack+5700,int_stack+5400,30);
     Libderiv->ABCD[119] = int_stack + 3600;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+4500,int_stack+6450,int_stack+6150,30);
     Libderiv->ABCD[118] = int_stack + 4500;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+5400,int_stack+40300,int_stack+6900,30);
     Libderiv->ABCD[117] = int_stack + 5400;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+6300,int_stack+7500,int_stack+7200,30);
     Libderiv->ABCD[116] = int_stack + 6300;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+48850,int_stack+8250,int_stack+7950,30);
     Libderiv->ABCD[107] = int_stack + 48850;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+7200,int_stack+9000,int_stack+8700,30);
     Libderiv->ABCD[106] = int_stack + 7200;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8100,int_stack+40750,int_stack+9450,30);
     Libderiv->ABCD[105] = int_stack + 8100;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+40300,int_stack+10050,int_stack+9750,30);
     Libderiv->ABCD[104] = int_stack + 40300;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+9000,int_stack+10800,int_stack+10500,30);
     Libderiv->ABCD[103] = int_stack + 9000;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+9900,int_stack+41200,int_stack+11250,30);
     Libderiv->ABCD[95] = int_stack + 9900;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+49750,int_stack+11850,int_stack+11550,30);
     Libderiv->ABCD[94] = int_stack + 49750;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+10800,int_stack+12600,int_stack+12300,30);
     Libderiv->ABCD[93] = int_stack + 10800;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11700,int_stack+13350,int_stack+13050,30);
     Libderiv->ABCD[92] = int_stack + 11700;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+12600,int_stack+14100,int_stack+13800,30);
     Libderiv->ABCD[91] = int_stack + 12600;
 /*--- compute (fp|fp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+13500,int_stack+14850,int_stack+14550,30);
     Libderiv->ABCD[90] = int_stack + 13500;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+14400,int_stack+15600,int_stack+15300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[47] = int_stack + 14400;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+50650,int_stack+16350,int_stack+16050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+34900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[46] = int_stack + 50650;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+15300,int_stack+17100,int_stack+16800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[45] = int_stack + 15300;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+16200,int_stack+17850,int_stack+17550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[44] = int_stack + 16200;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+17100,int_stack+18600,int_stack+18300, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[43] = int_stack + 17100;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+18000,int_stack+41650,int_stack+19050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+37000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[42] = int_stack + 18000;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+41200,int_stack+19650,int_stack+19350, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[38] = int_stack + 41200;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+18900,int_stack+20400,int_stack+20100, 0.0, zero_stack, 1.0, int_stack+34150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[35] = int_stack + 18900;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+19800,int_stack+21150,int_stack+20850, 0.0, zero_stack, 1.0, int_stack+34900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[34] = int_stack + 19800;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+20700,int_stack+21900,int_stack+21600, 0.0, zero_stack, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[33] = int_stack + 20700;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+51550,int_stack+42100,int_stack+22350, 0.0, zero_stack, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[32] = int_stack + 51550;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+21600,int_stack+22950,int_stack+22650, 0.0, zero_stack, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[31] = int_stack + 21600;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+22500,int_stack+23700,int_stack+23400, 0.0, zero_stack, 1.0, int_stack+37000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[30] = int_stack + 22500;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+52450,int_stack+24450,int_stack+24150, 0.0, zero_stack, 1.0, int_stack+1200, 1.0, int_stack+38200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[26] = int_stack + 52450;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+23400,int_stack+25200,int_stack+24900, 0.0, zero_stack, 2.0, int_stack+38200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[25] = int_stack + 23400;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+24300,int_stack+25950,int_stack+25650, 1.0, int_stack+34150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[23] = int_stack + 24300;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+25200,int_stack+26700,int_stack+26400, 1.0, int_stack+34900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[22] = int_stack + 25200;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+26100,int_stack+27450,int_stack+27150, 1.0, int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[21] = int_stack + 26100;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+27000,int_stack+28200,int_stack+27900, 1.0, int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[20] = int_stack + 27000;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+33300,int_stack+28950,int_stack+28650, 1.0, int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[19] = int_stack + 33300;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+27900,int_stack+29700,int_stack+29400, 1.0, int_stack+37000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[18] = int_stack + 27900;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+28800,int_stack+30450,int_stack+30150, 1.0, int_stack+1200, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[14] = int_stack + 28800;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+29700,int_stack+31200,int_stack+30900, 1.0, int_stack+38200, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[13] = int_stack + 29700;
 /*--- compute (fp|fp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+30600,int_stack+31950,int_stack+31650, 2.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[12] = int_stack + 30600;

}
