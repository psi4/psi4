#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0f0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|f0) integrals */

void d12hrr_order_p0f0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[1][3][10] = int_stack + 30;
 Libderiv->deriv_classes[1][3][9] = int_stack + 60;
 Libderiv->deriv_classes[1][3][8] = int_stack + 90;
 Libderiv->deriv_classes[1][3][7] = int_stack + 120;
 Libderiv->deriv_classes[1][3][6] = int_stack + 150;
 Libderiv->deriv_classes[1][3][2] = int_stack + 180;
 Libderiv->deriv_classes[1][3][1] = int_stack + 210;
 Libderiv->deriv_classes[1][3][0] = int_stack + 240;
 Libderiv->deriv2_classes[1][3][143] = int_stack + 270;
 Libderiv->deriv2_classes[1][3][131] = int_stack + 300;
 Libderiv->deriv2_classes[1][3][130] = int_stack + 330;
 Libderiv->deriv2_classes[1][3][119] = int_stack + 360;
 Libderiv->deriv2_classes[1][3][118] = int_stack + 390;
 Libderiv->deriv2_classes[1][3][117] = int_stack + 420;
 Libderiv->deriv2_classes[1][3][107] = int_stack + 450;
 Libderiv->deriv2_classes[1][3][106] = int_stack + 480;
 Libderiv->deriv2_classes[1][3][105] = int_stack + 510;
 Libderiv->deriv2_classes[1][3][104] = int_stack + 540;
 Libderiv->deriv2_classes[1][3][95] = int_stack + 570;
 Libderiv->deriv2_classes[1][3][94] = int_stack + 600;
 Libderiv->deriv2_classes[1][3][93] = int_stack + 630;
 Libderiv->deriv2_classes[1][3][92] = int_stack + 660;
 Libderiv->deriv2_classes[1][3][91] = int_stack + 690;
 Libderiv->deriv2_classes[1][3][83] = int_stack + 720;
 Libderiv->deriv2_classes[1][3][82] = int_stack + 750;
 Libderiv->deriv2_classes[1][3][81] = int_stack + 780;
 Libderiv->deriv2_classes[1][3][80] = int_stack + 810;
 Libderiv->deriv2_classes[1][3][79] = int_stack + 840;
 Libderiv->deriv2_classes[1][3][78] = int_stack + 870;
 Libderiv->deriv2_classes[1][3][35] = int_stack + 900;
 Libderiv->deriv2_classes[1][3][34] = int_stack + 930;
 Libderiv->deriv2_classes[1][3][33] = int_stack + 960;
 Libderiv->deriv2_classes[1][3][32] = int_stack + 990;
 Libderiv->deriv2_classes[1][3][31] = int_stack + 1020;
 Libderiv->deriv2_classes[1][3][30] = int_stack + 1050;
 Libderiv->deriv2_classes[1][3][26] = int_stack + 1080;
 Libderiv->deriv2_classes[1][3][23] = int_stack + 1110;
 Libderiv->deriv2_classes[1][3][22] = int_stack + 1140;
 Libderiv->deriv2_classes[1][3][21] = int_stack + 1170;
 Libderiv->deriv2_classes[1][3][20] = int_stack + 1200;
 Libderiv->deriv2_classes[1][3][19] = int_stack + 1230;
 Libderiv->deriv2_classes[1][3][18] = int_stack + 1260;
 Libderiv->deriv2_classes[1][3][14] = int_stack + 1290;
 Libderiv->deriv2_classes[1][3][13] = int_stack + 1320;
 Libderiv->deriv2_classes[1][3][11] = int_stack + 1350;
 Libderiv->deriv2_classes[1][3][10] = int_stack + 1380;
 Libderiv->deriv2_classes[1][3][9] = int_stack + 1410;
 Libderiv->deriv2_classes[1][3][8] = int_stack + 1440;
 Libderiv->deriv2_classes[1][3][7] = int_stack + 1470;
 Libderiv->deriv2_classes[1][3][6] = int_stack + 1500;
 Libderiv->deriv2_classes[1][3][2] = int_stack + 1530;
 Libderiv->deriv2_classes[1][3][1] = int_stack + 1560;
 Libderiv->deriv2_classes[1][3][0] = int_stack + 1590;
 memset(int_stack,0,12960);

 Libderiv->dvrr_stack = int_stack + 1620;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0f0(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[10] = int_stack + 30;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[9] = int_stack + 60;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[8] = int_stack + 90;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[7] = int_stack + 120;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[6] = int_stack + 150;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[2] = int_stack + 180;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[1] = int_stack + 210;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[0] = int_stack + 240;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[155] = int_stack + 270;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[143] = int_stack + 300;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[142] = int_stack + 330;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[131] = int_stack + 360;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[130] = int_stack + 390;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[129] = int_stack + 420;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[119] = int_stack + 450;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[118] = int_stack + 480;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[117] = int_stack + 510;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[116] = int_stack + 540;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[107] = int_stack + 570;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[106] = int_stack + 600;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[105] = int_stack + 630;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[104] = int_stack + 660;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[103] = int_stack + 690;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[95] = int_stack + 720;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[94] = int_stack + 750;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[93] = int_stack + 780;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[92] = int_stack + 810;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[91] = int_stack + 840;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[90] = int_stack + 870;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[47] = int_stack + 900;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[46] = int_stack + 930;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[45] = int_stack + 960;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[44] = int_stack + 990;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[43] = int_stack + 1020;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[42] = int_stack + 1050;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[38] = int_stack + 1080;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[35] = int_stack + 1110;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[34] = int_stack + 1140;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[33] = int_stack + 1170;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[32] = int_stack + 1200;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[31] = int_stack + 1230;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[30] = int_stack + 1260;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[26] = int_stack + 1290;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[25] = int_stack + 1320;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[23] = int_stack + 1350;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[22] = int_stack + 1380;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[21] = int_stack + 1410;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[20] = int_stack + 1440;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[19] = int_stack + 1470;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[18] = int_stack + 1500;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[14] = int_stack + 1530;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[13] = int_stack + 1560;
 /*--- compute (p0|f0) ---*/
     Libderiv->ABCD[12] = int_stack + 1590;

}
