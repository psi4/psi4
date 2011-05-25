#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_d0f0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|f0) integrals */

void d12hrr_order_d0f0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][3][10] = int_stack + 60;
 Libderiv->deriv_classes[2][3][9] = int_stack + 120;
 Libderiv->deriv_classes[2][3][8] = int_stack + 180;
 Libderiv->deriv_classes[2][3][7] = int_stack + 240;
 Libderiv->deriv_classes[2][3][6] = int_stack + 300;
 Libderiv->deriv_classes[2][3][2] = int_stack + 360;
 Libderiv->deriv_classes[2][3][1] = int_stack + 420;
 Libderiv->deriv_classes[2][3][0] = int_stack + 480;
 Libderiv->deriv2_classes[2][3][143] = int_stack + 540;
 Libderiv->deriv2_classes[2][3][131] = int_stack + 600;
 Libderiv->deriv2_classes[2][3][130] = int_stack + 660;
 Libderiv->deriv2_classes[2][3][119] = int_stack + 720;
 Libderiv->deriv2_classes[2][3][118] = int_stack + 780;
 Libderiv->deriv2_classes[2][3][117] = int_stack + 840;
 Libderiv->deriv2_classes[2][3][107] = int_stack + 900;
 Libderiv->deriv2_classes[2][3][106] = int_stack + 960;
 Libderiv->deriv2_classes[2][3][105] = int_stack + 1020;
 Libderiv->deriv2_classes[2][3][104] = int_stack + 1080;
 Libderiv->deriv2_classes[2][3][95] = int_stack + 1140;
 Libderiv->deriv2_classes[2][3][94] = int_stack + 1200;
 Libderiv->deriv2_classes[2][3][93] = int_stack + 1260;
 Libderiv->deriv2_classes[2][3][92] = int_stack + 1320;
 Libderiv->deriv2_classes[2][3][91] = int_stack + 1380;
 Libderiv->deriv2_classes[2][3][83] = int_stack + 1440;
 Libderiv->deriv2_classes[2][3][82] = int_stack + 1500;
 Libderiv->deriv2_classes[2][3][81] = int_stack + 1560;
 Libderiv->deriv2_classes[2][3][80] = int_stack + 1620;
 Libderiv->deriv2_classes[2][3][79] = int_stack + 1680;
 Libderiv->deriv2_classes[2][3][78] = int_stack + 1740;
 Libderiv->deriv2_classes[2][3][35] = int_stack + 1800;
 Libderiv->deriv2_classes[2][3][34] = int_stack + 1860;
 Libderiv->deriv2_classes[2][3][33] = int_stack + 1920;
 Libderiv->deriv2_classes[2][3][32] = int_stack + 1980;
 Libderiv->deriv2_classes[2][3][31] = int_stack + 2040;
 Libderiv->deriv2_classes[2][3][30] = int_stack + 2100;
 Libderiv->deriv2_classes[2][3][26] = int_stack + 2160;
 Libderiv->deriv2_classes[2][3][23] = int_stack + 2220;
 Libderiv->deriv2_classes[2][3][22] = int_stack + 2280;
 Libderiv->deriv2_classes[2][3][21] = int_stack + 2340;
 Libderiv->deriv2_classes[2][3][20] = int_stack + 2400;
 Libderiv->deriv2_classes[2][3][19] = int_stack + 2460;
 Libderiv->deriv2_classes[2][3][18] = int_stack + 2520;
 Libderiv->deriv2_classes[2][3][14] = int_stack + 2580;
 Libderiv->deriv2_classes[2][3][13] = int_stack + 2640;
 Libderiv->deriv2_classes[2][3][11] = int_stack + 2700;
 Libderiv->deriv2_classes[2][3][10] = int_stack + 2760;
 Libderiv->deriv2_classes[2][3][9] = int_stack + 2820;
 Libderiv->deriv2_classes[2][3][8] = int_stack + 2880;
 Libderiv->deriv2_classes[2][3][7] = int_stack + 2940;
 Libderiv->deriv2_classes[2][3][6] = int_stack + 3000;
 Libderiv->deriv2_classes[2][3][2] = int_stack + 3060;
 Libderiv->deriv2_classes[2][3][1] = int_stack + 3120;
 Libderiv->deriv2_classes[2][3][0] = int_stack + 3180;
 memset(int_stack,0,25920);

 Libderiv->dvrr_stack = int_stack + 3240;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_d0f0(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[10] = int_stack + 60;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[9] = int_stack + 120;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[8] = int_stack + 180;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[7] = int_stack + 240;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[6] = int_stack + 300;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[2] = int_stack + 360;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[1] = int_stack + 420;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[0] = int_stack + 480;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[155] = int_stack + 540;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[143] = int_stack + 600;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[142] = int_stack + 660;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[131] = int_stack + 720;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[130] = int_stack + 780;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[129] = int_stack + 840;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[119] = int_stack + 900;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[118] = int_stack + 960;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[117] = int_stack + 1020;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[116] = int_stack + 1080;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[107] = int_stack + 1140;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[106] = int_stack + 1200;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[105] = int_stack + 1260;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[104] = int_stack + 1320;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[103] = int_stack + 1380;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[95] = int_stack + 1440;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[94] = int_stack + 1500;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[93] = int_stack + 1560;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[92] = int_stack + 1620;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[91] = int_stack + 1680;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[90] = int_stack + 1740;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[47] = int_stack + 1800;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[46] = int_stack + 1860;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[45] = int_stack + 1920;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[44] = int_stack + 1980;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[43] = int_stack + 2040;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[42] = int_stack + 2100;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[38] = int_stack + 2160;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[35] = int_stack + 2220;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[34] = int_stack + 2280;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[33] = int_stack + 2340;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[32] = int_stack + 2400;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[31] = int_stack + 2460;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[30] = int_stack + 2520;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[26] = int_stack + 2580;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[25] = int_stack + 2640;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[23] = int_stack + 2700;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[22] = int_stack + 2760;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[21] = int_stack + 2820;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[20] = int_stack + 2880;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[19] = int_stack + 2940;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[18] = int_stack + 3000;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[14] = int_stack + 3060;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[13] = int_stack + 3120;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[12] = int_stack + 3180;

}
