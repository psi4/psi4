#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_f0dp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|dp) integrals */

void d12hrr_order_f0dp(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->dvrr_classes[3][2] = int_stack + 500;
 Libderiv->deriv_classes[3][3][6] = int_stack + 560;
 Libderiv->deriv_classes[3][3][2] = int_stack + 660;
 Libderiv->deriv_classes[3][3][1] = int_stack + 760;
 Libderiv->deriv_classes[3][3][0] = int_stack + 860;
 Libderiv->deriv2_classes[3][2][143] = int_stack + 960;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 1020;
 Libderiv->deriv2_classes[3][2][131] = int_stack + 1120;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 1180;
 Libderiv->deriv2_classes[3][2][130] = int_stack + 1280;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 1340;
 Libderiv->deriv2_classes[3][2][119] = int_stack + 1440;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 1500;
 Libderiv->deriv2_classes[3][2][118] = int_stack + 1600;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 1660;
 Libderiv->deriv2_classes[3][2][117] = int_stack + 1760;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 1820;
 Libderiv->deriv2_classes[3][2][107] = int_stack + 1920;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 1980;
 Libderiv->deriv2_classes[3][2][106] = int_stack + 2080;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 2140;
 Libderiv->deriv2_classes[3][2][105] = int_stack + 2240;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 2300;
 Libderiv->deriv2_classes[3][2][104] = int_stack + 2400;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 2460;
 Libderiv->deriv2_classes[3][2][95] = int_stack + 2560;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 2620;
 Libderiv->deriv2_classes[3][2][94] = int_stack + 2720;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 2780;
 Libderiv->deriv2_classes[3][2][93] = int_stack + 2880;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 2940;
 Libderiv->deriv2_classes[3][2][92] = int_stack + 3040;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 3100;
 Libderiv->deriv2_classes[3][2][91] = int_stack + 3200;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 3260;
 Libderiv->deriv_classes[3][2][11] = int_stack + 3360;
 Libderiv->deriv2_classes[3][2][83] = int_stack + 3420;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 3480;
 Libderiv->deriv_classes[3][2][10] = int_stack + 3580;
 Libderiv->deriv2_classes[3][2][82] = int_stack + 3640;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 3700;
 Libderiv->deriv_classes[3][2][9] = int_stack + 3800;
 Libderiv->deriv2_classes[3][2][81] = int_stack + 3860;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 3920;
 Libderiv->deriv_classes[3][2][8] = int_stack + 4020;
 Libderiv->deriv2_classes[3][2][80] = int_stack + 4080;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 4140;
 Libderiv->deriv_classes[3][2][7] = int_stack + 4240;
 Libderiv->deriv2_classes[3][2][79] = int_stack + 4300;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 4360;
 Libderiv->deriv_classes[3][2][6] = int_stack + 4460;
 Libderiv->deriv2_classes[3][2][78] = int_stack + 4520;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 4580;
 Libderiv->deriv2_classes[3][2][35] = int_stack + 4680;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 4740;
 Libderiv->deriv2_classes[3][2][34] = int_stack + 4840;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 4900;
 Libderiv->deriv2_classes[3][2][33] = int_stack + 5000;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 5060;
 Libderiv->deriv2_classes[3][2][32] = int_stack + 5160;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 5220;
 Libderiv->deriv2_classes[3][2][31] = int_stack + 5320;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 5380;
 Libderiv->deriv_classes[3][2][2] = int_stack + 5480;
 Libderiv->deriv2_classes[3][2][30] = int_stack + 5540;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 5600;
 Libderiv->deriv2_classes[3][2][26] = int_stack + 5700;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 5760;
 Libderiv->deriv2_classes[3][2][23] = int_stack + 5860;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 5920;
 Libderiv->deriv2_classes[3][2][22] = int_stack + 6020;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 6080;
 Libderiv->deriv2_classes[3][2][21] = int_stack + 6180;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 6240;
 Libderiv->deriv2_classes[3][2][20] = int_stack + 6340;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 6400;
 Libderiv->deriv2_classes[3][2][19] = int_stack + 6500;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 6560;
 Libderiv->deriv_classes[3][2][1] = int_stack + 6660;
 Libderiv->deriv2_classes[3][2][18] = int_stack + 6720;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 6780;
 Libderiv->deriv2_classes[3][2][14] = int_stack + 6880;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 6940;
 Libderiv->deriv2_classes[3][2][13] = int_stack + 7040;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 7100;
 Libderiv->deriv2_classes[3][2][11] = int_stack + 7200;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 7260;
 Libderiv->deriv2_classes[3][2][10] = int_stack + 7360;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 7420;
 Libderiv->deriv2_classes[3][2][9] = int_stack + 7520;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 7580;
 Libderiv->deriv2_classes[3][2][8] = int_stack + 7680;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 7740;
 Libderiv->deriv2_classes[3][2][7] = int_stack + 7840;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 7900;
 Libderiv->deriv_classes[3][2][0] = int_stack + 8000;
 Libderiv->deriv2_classes[3][2][6] = int_stack + 8060;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 8120;
 Libderiv->deriv2_classes[3][2][2] = int_stack + 8220;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 8280;
 Libderiv->deriv2_classes[3][2][1] = int_stack + 8380;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 8440;
 Libderiv->deriv2_classes[3][2][0] = int_stack + 8540;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 8600;
 memset(int_stack,0,69600);

 Libderiv->dvrr_stack = int_stack + 9960;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_f0dp(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8700,int_stack+0,int_stack+3360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500,10);
     Libderiv->ABCD[11] = int_stack + 8700;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+8880,int_stack+100,int_stack+3580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 8880;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+200,int_stack+3800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9060,int_stack+300,int_stack+4020, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 9060;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+180,int_stack+400,int_stack+4240, 0.0, zero_stack, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 180;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9240,int_stack+560,int_stack+4460, 1.0, int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 9240;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+360,int_stack+660,int_stack+5480,10);
     Libderiv->ABCD[2] = int_stack + 360;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+760,int_stack+6660,10);
     Libderiv->ABCD[1] = int_stack + 540;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+9420,int_stack+860,int_stack+8000,10);
     Libderiv->ABCD[0] = int_stack + 9420;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+720,int_stack+1020,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3360,10);
     Libderiv->ABCD[155] = int_stack + 720;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+900,int_stack+1180,int_stack+1120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3360, 1.0, int_stack+3580,10);
     Libderiv->ABCD[143] = int_stack + 900;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1340,int_stack+1280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3580, 0.0, zero_stack,10);
     Libderiv->ABCD[142] = int_stack + 1080;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1260,int_stack+1500,int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3360, 0.0, zero_stack, 1.0, int_stack+3800,10);
     Libderiv->ABCD[131] = int_stack + 1260;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9600,int_stack+1660,int_stack+1600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3580, 1.0, int_stack+3800, 0.0, zero_stack,10);
     Libderiv->ABCD[130] = int_stack + 9600;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1440,int_stack+1820,int_stack+1760, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+3800, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[129] = int_stack + 1440;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1620,int_stack+1980,int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4020,10);
     Libderiv->ABCD[119] = int_stack + 1620;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1800,int_stack+2140,int_stack+2080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3580, 0.0, zero_stack, 1.0, int_stack+4020, 0.0, zero_stack,10);
     Libderiv->ABCD[118] = int_stack + 1800;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1980,int_stack+2300,int_stack+2240, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3800, 1.0, int_stack+4020, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[117] = int_stack + 1980;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2160,int_stack+2460,int_stack+2400, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+4020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[116] = int_stack + 2160;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2340,int_stack+2620,int_stack+2560, 0.0, zero_stack, 1.0, int_stack+3360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4240,10);
     Libderiv->ABCD[107] = int_stack + 2340;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2520,int_stack+2780,int_stack+2720, 0.0, zero_stack, 1.0, int_stack+3580, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4240, 0.0, zero_stack,10);
     Libderiv->ABCD[106] = int_stack + 2520;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2700,int_stack+2940,int_stack+2880, 0.0, zero_stack, 1.0, int_stack+3800, 0.0, zero_stack, 1.0, int_stack+4240, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[105] = int_stack + 2700;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+9780,int_stack+3100,int_stack+3040, 0.0, zero_stack, 1.0, int_stack+4020, 1.0, int_stack+4240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[104] = int_stack + 9780;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+2880,int_stack+3260,int_stack+3200, 0.0, zero_stack, 2.0, int_stack+4240, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[103] = int_stack + 2880;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3060,int_stack+3480,int_stack+3420, 1.0, int_stack+3360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4460,10);
     Libderiv->ABCD[95] = int_stack + 3060;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3240,int_stack+3700,int_stack+3640, 1.0, int_stack+3580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4460, 0.0, zero_stack,10);
     Libderiv->ABCD[94] = int_stack + 3240;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3420,int_stack+3920,int_stack+3860, 1.0, int_stack+3800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4460, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[93] = int_stack + 3420;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3600,int_stack+4140,int_stack+4080, 1.0, int_stack+4020, 0.0, zero_stack, 1.0, int_stack+4460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[92] = int_stack + 3600;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3780,int_stack+4360,int_stack+4300, 1.0, int_stack+4240, 1.0, int_stack+4460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[91] = int_stack + 3780;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3960,int_stack+4580,int_stack+4520, 2.0, int_stack+4460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[90] = int_stack + 3960;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4140,int_stack+4740,int_stack+4680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5480,10);
     Libderiv->ABCD[47] = int_stack + 4140;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4320,int_stack+4900,int_stack+4840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5480, 0.0, zero_stack,10);
     Libderiv->ABCD[46] = int_stack + 4320;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4500,int_stack+5060,int_stack+5000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5480, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[45] = int_stack + 4500;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4680,int_stack+5220,int_stack+5160, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[44] = int_stack + 4680;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+4860,int_stack+5380,int_stack+5320, 0.0, zero_stack, 1.0, int_stack+5480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[43] = int_stack + 4860;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5040,int_stack+5600,int_stack+5540, 1.0, int_stack+5480, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[42] = int_stack + 5040;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+5220,int_stack+5760,int_stack+5700,10);
     Libderiv->ABCD[38] = int_stack + 5220;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5400,int_stack+5920,int_stack+5860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6660,10);
     Libderiv->ABCD[35] = int_stack + 5400;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5580,int_stack+6080,int_stack+6020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6660, 0.0, zero_stack,10);
     Libderiv->ABCD[34] = int_stack + 5580;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5760,int_stack+6240,int_stack+6180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6660, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[33] = int_stack + 5760;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+5940,int_stack+6400,int_stack+6340, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[32] = int_stack + 5940;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6120,int_stack+6560,int_stack+6500, 0.0, zero_stack, 1.0, int_stack+6660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[31] = int_stack + 6120;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6300,int_stack+6780,int_stack+6720, 1.0, int_stack+6660, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[30] = int_stack + 6300;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+6480,int_stack+6940,int_stack+6880,10);
     Libderiv->ABCD[26] = int_stack + 6480;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+6660,int_stack+7100,int_stack+7040,10);
     Libderiv->ABCD[25] = int_stack + 6660;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+6840,int_stack+7260,int_stack+7200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8000,10);
     Libderiv->ABCD[23] = int_stack + 6840;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7020,int_stack+7420,int_stack+7360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8000, 0.0, zero_stack,10);
     Libderiv->ABCD[22] = int_stack + 7020;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7200,int_stack+7580,int_stack+7520, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8000, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[21] = int_stack + 7200;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7380,int_stack+7740,int_stack+7680, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[20] = int_stack + 7380;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7560,int_stack+7900,int_stack+7840, 0.0, zero_stack, 1.0, int_stack+8000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[19] = int_stack + 7560;
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+7740,int_stack+8120,int_stack+8060, 1.0, int_stack+8000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[18] = int_stack + 7740;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+7920,int_stack+8280,int_stack+8220,10);
     Libderiv->ABCD[14] = int_stack + 7920;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+8100,int_stack+8440,int_stack+8380,10);
     Libderiv->ABCD[13] = int_stack + 8100;
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+8280,int_stack+8600,int_stack+8540,10);
     Libderiv->ABCD[12] = int_stack + 8280;

}
