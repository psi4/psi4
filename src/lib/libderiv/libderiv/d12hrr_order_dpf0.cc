#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_dpf0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|f0) integrals */

void d12hrr_order_dpf0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][3][10] = int_stack + 100;
 Libderiv->deriv_classes[3][3][9] = int_stack + 200;
 Libderiv->deriv_classes[3][3][8] = int_stack + 300;
 Libderiv->deriv_classes[3][3][7] = int_stack + 400;
 Libderiv->deriv_classes[3][3][6] = int_stack + 500;
 Libderiv->deriv_classes[3][3][2] = int_stack + 600;
 Libderiv->deriv_classes[3][3][1] = int_stack + 700;
 Libderiv->dvrr_classes[2][3] = int_stack + 800;
 Libderiv->deriv_classes[3][3][0] = int_stack + 860;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 960;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 1020;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 1120;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 1180;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 1280;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 1340;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 1440;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 1500;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 1600;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 1660;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 1760;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 1820;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 1920;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 1980;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 2080;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 2140;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 2240;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 2300;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 2400;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 2460;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 2560;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 2620;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 2720;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 2780;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 2880;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 2940;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 3040;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 3100;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 3200;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 3260;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 3360;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 3420;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 3520;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 3580;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 3680;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 3740;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 3840;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 3900;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 4000;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 4060;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 4160;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 4220;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 4320;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 4380;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 4480;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 4540;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 4640;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 4700;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 4800;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 4860;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 4960;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 5020;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 5120;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 5180;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 5280;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 5340;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 5440;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 5500;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 5600;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 5660;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 5760;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 5820;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 5920;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 5980;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 6080;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 6140;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 6240;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 6300;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 6400;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 6460;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 6560;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 6620;
 Libderiv->deriv_classes[2][3][11] = int_stack + 6720;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 6780;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 6840;
 Libderiv->deriv_classes[2][3][10] = int_stack + 6940;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 7000;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 7060;
 Libderiv->deriv_classes[2][3][9] = int_stack + 7160;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 7220;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 7280;
 Libderiv->deriv_classes[2][3][8] = int_stack + 7380;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 7440;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 7500;
 Libderiv->deriv_classes[2][3][7] = int_stack + 7600;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 7660;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 7720;
 Libderiv->deriv_classes[2][3][6] = int_stack + 7820;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 7880;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 7940;
 Libderiv->deriv_classes[2][3][2] = int_stack + 8040;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 8100;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 8160;
 Libderiv->deriv_classes[2][3][1] = int_stack + 8260;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 8320;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 8380;
 Libderiv->deriv_classes[2][3][0] = int_stack + 8480;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 8540;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 8600;
 memset(int_stack,0,69600);

 Libderiv->dvrr_stack = int_stack + 10320;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_dpf0(Libderiv, Data);
   Data++;
 }

 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8700,int_stack+0,int_stack+6720,10);
     Libderiv->ABCD[11] = int_stack + 8700;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+8880,int_stack+100,int_stack+6940,10);
     Libderiv->ABCD[10] = int_stack + 8880;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+200,int_stack+7160,10);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9060,int_stack+300,int_stack+7380,10);
     Libderiv->ABCD[8] = int_stack + 9060;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+180,int_stack+400,int_stack+7600,10);
     Libderiv->ABCD[7] = int_stack + 180;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9240,int_stack+500,int_stack+7820,10);
     Libderiv->ABCD[6] = int_stack + 9240;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+360,int_stack+600,int_stack+8040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[2] = int_stack + 360;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9420,int_stack+700,int_stack+8260, 0.0, zero_stack, 1.0, int_stack+800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[1] = int_stack + 9420;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+540,int_stack+860,int_stack+8480, 1.0, int_stack+800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[0] = int_stack + 540;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+720,int_stack+1020,int_stack+960,10);
     Libderiv->ABCD[155] = int_stack + 720;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+900,int_stack+1180,int_stack+1120,10);
     Libderiv->ABCD[143] = int_stack + 900;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1080,int_stack+1340,int_stack+1280,10);
     Libderiv->ABCD[142] = int_stack + 1080;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1260,int_stack+1500,int_stack+1440,10);
     Libderiv->ABCD[131] = int_stack + 1260;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9600,int_stack+1660,int_stack+1600,10);
     Libderiv->ABCD[130] = int_stack + 9600;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1440,int_stack+1820,int_stack+1760,10);
     Libderiv->ABCD[129] = int_stack + 1440;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1620,int_stack+1980,int_stack+1920,10);
     Libderiv->ABCD[119] = int_stack + 1620;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1800,int_stack+2140,int_stack+2080,10);
     Libderiv->ABCD[118] = int_stack + 1800;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+1980,int_stack+2300,int_stack+2240,10);
     Libderiv->ABCD[117] = int_stack + 1980;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2160,int_stack+2460,int_stack+2400,10);
     Libderiv->ABCD[116] = int_stack + 2160;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2340,int_stack+2620,int_stack+2560,10);
     Libderiv->ABCD[107] = int_stack + 2340;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2520,int_stack+2780,int_stack+2720,10);
     Libderiv->ABCD[106] = int_stack + 2520;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2700,int_stack+2940,int_stack+2880,10);
     Libderiv->ABCD[105] = int_stack + 2700;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+9780,int_stack+3100,int_stack+3040,10);
     Libderiv->ABCD[104] = int_stack + 9780;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2880,int_stack+3260,int_stack+3200,10);
     Libderiv->ABCD[103] = int_stack + 2880;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3060,int_stack+3420,int_stack+3360,10);
     Libderiv->ABCD[95] = int_stack + 3060;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3240,int_stack+3580,int_stack+3520,10);
     Libderiv->ABCD[94] = int_stack + 3240;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3420,int_stack+3740,int_stack+3680,10);
     Libderiv->ABCD[93] = int_stack + 3420;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3600,int_stack+3900,int_stack+3840,10);
     Libderiv->ABCD[92] = int_stack + 3600;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3780,int_stack+4060,int_stack+4000,10);
     Libderiv->ABCD[91] = int_stack + 3780;
 /*--- compute (dp|f0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+3960,int_stack+4220,int_stack+4160,10);
     Libderiv->ABCD[90] = int_stack + 3960;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+4140,int_stack+4380,int_stack+4320, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[47] = int_stack + 4140;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+9960,int_stack+4540,int_stack+4480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[46] = int_stack + 9960;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+4320,int_stack+4700,int_stack+4640, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[45] = int_stack + 4320;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+4500,int_stack+4860,int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[44] = int_stack + 4500;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+4680,int_stack+5020,int_stack+4960, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[43] = int_stack + 4680;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+4860,int_stack+5180,int_stack+5120, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[42] = int_stack + 4860;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5040,int_stack+5340,int_stack+5280, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+8040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[38] = int_stack + 5040;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5220,int_stack+5500,int_stack+5440, 0.0, zero_stack, 1.0, int_stack+6720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[35] = int_stack + 5220;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5400,int_stack+5660,int_stack+5600, 0.0, zero_stack, 1.0, int_stack+6940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[34] = int_stack + 5400;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5580,int_stack+5820,int_stack+5760, 0.0, zero_stack, 1.0, int_stack+7160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[33] = int_stack + 5580;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+10140,int_stack+5980,int_stack+5920, 0.0, zero_stack, 1.0, int_stack+7380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[32] = int_stack + 10140;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5760,int_stack+6140,int_stack+6080, 0.0, zero_stack, 1.0, int_stack+7600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[31] = int_stack + 5760;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+5940,int_stack+6300,int_stack+6240, 0.0, zero_stack, 1.0, int_stack+7820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[30] = int_stack + 5940;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6120,int_stack+6460,int_stack+6400, 0.0, zero_stack, 1.0, int_stack+8040, 1.0, int_stack+8260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[26] = int_stack + 6120;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6300,int_stack+6620,int_stack+6560, 0.0, zero_stack, 2.0, int_stack+8260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[25] = int_stack + 6300;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6480,int_stack+6840,int_stack+6780, 1.0, int_stack+6720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[23] = int_stack + 6480;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6660,int_stack+7060,int_stack+7000, 1.0, int_stack+6940, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[22] = int_stack + 6660;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+6840,int_stack+7280,int_stack+7220, 1.0, int_stack+7160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[21] = int_stack + 6840;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7020,int_stack+7500,int_stack+7440, 1.0, int_stack+7380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[20] = int_stack + 7020;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7200,int_stack+7720,int_stack+7660, 1.0, int_stack+7600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[19] = int_stack + 7200;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7380,int_stack+7940,int_stack+7880, 1.0, int_stack+7820, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[18] = int_stack + 7380;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7560,int_stack+8160,int_stack+8100, 1.0, int_stack+8040, 0.0, zero_stack, 1.0, int_stack+8480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[14] = int_stack + 7560;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7740,int_stack+8380,int_stack+8320, 1.0, int_stack+8260, 1.0, int_stack+8480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[13] = int_stack + 7740;
 /*--- compute (dp|f0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+7920,int_stack+8600,int_stack+8540, 2.0, int_stack+8480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[12] = int_stack + 7920;

}
