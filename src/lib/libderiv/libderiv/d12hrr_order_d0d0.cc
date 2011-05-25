#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_d0d0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|d0) integrals */

void d12hrr_order_d0d0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][2][11] = int_stack + 0;
 Libderiv->deriv_classes[2][2][10] = int_stack + 36;
 Libderiv->deriv_classes[2][2][9] = int_stack + 72;
 Libderiv->deriv_classes[2][2][8] = int_stack + 108;
 Libderiv->deriv_classes[2][2][7] = int_stack + 144;
 Libderiv->deriv_classes[2][2][6] = int_stack + 180;
 Libderiv->deriv_classes[2][2][2] = int_stack + 216;
 Libderiv->deriv_classes[2][2][1] = int_stack + 252;
 Libderiv->deriv_classes[2][2][0] = int_stack + 288;
 Libderiv->deriv2_classes[2][2][143] = int_stack + 324;
 Libderiv->deriv2_classes[2][2][131] = int_stack + 360;
 Libderiv->deriv2_classes[2][2][130] = int_stack + 396;
 Libderiv->deriv2_classes[2][2][119] = int_stack + 432;
 Libderiv->deriv2_classes[2][2][118] = int_stack + 468;
 Libderiv->deriv2_classes[2][2][117] = int_stack + 504;
 Libderiv->deriv2_classes[2][2][107] = int_stack + 540;
 Libderiv->deriv2_classes[2][2][106] = int_stack + 576;
 Libderiv->deriv2_classes[2][2][105] = int_stack + 612;
 Libderiv->deriv2_classes[2][2][104] = int_stack + 648;
 Libderiv->deriv2_classes[2][2][95] = int_stack + 684;
 Libderiv->deriv2_classes[2][2][94] = int_stack + 720;
 Libderiv->deriv2_classes[2][2][93] = int_stack + 756;
 Libderiv->deriv2_classes[2][2][92] = int_stack + 792;
 Libderiv->deriv2_classes[2][2][91] = int_stack + 828;
 Libderiv->deriv2_classes[2][2][83] = int_stack + 864;
 Libderiv->deriv2_classes[2][2][82] = int_stack + 900;
 Libderiv->deriv2_classes[2][2][81] = int_stack + 936;
 Libderiv->deriv2_classes[2][2][80] = int_stack + 972;
 Libderiv->deriv2_classes[2][2][79] = int_stack + 1008;
 Libderiv->deriv2_classes[2][2][78] = int_stack + 1044;
 Libderiv->deriv2_classes[2][2][35] = int_stack + 1080;
 Libderiv->deriv2_classes[2][2][34] = int_stack + 1116;
 Libderiv->deriv2_classes[2][2][33] = int_stack + 1152;
 Libderiv->deriv2_classes[2][2][32] = int_stack + 1188;
 Libderiv->deriv2_classes[2][2][31] = int_stack + 1224;
 Libderiv->deriv2_classes[2][2][30] = int_stack + 1260;
 Libderiv->deriv2_classes[2][2][26] = int_stack + 1296;
 Libderiv->deriv2_classes[2][2][23] = int_stack + 1332;
 Libderiv->deriv2_classes[2][2][22] = int_stack + 1368;
 Libderiv->deriv2_classes[2][2][21] = int_stack + 1404;
 Libderiv->deriv2_classes[2][2][20] = int_stack + 1440;
 Libderiv->deriv2_classes[2][2][19] = int_stack + 1476;
 Libderiv->deriv2_classes[2][2][18] = int_stack + 1512;
 Libderiv->deriv2_classes[2][2][14] = int_stack + 1548;
 Libderiv->deriv2_classes[2][2][13] = int_stack + 1584;
 Libderiv->deriv2_classes[2][2][11] = int_stack + 1620;
 Libderiv->deriv2_classes[2][2][10] = int_stack + 1656;
 Libderiv->deriv2_classes[2][2][9] = int_stack + 1692;
 Libderiv->deriv2_classes[2][2][8] = int_stack + 1728;
 Libderiv->deriv2_classes[2][2][7] = int_stack + 1764;
 Libderiv->deriv2_classes[2][2][6] = int_stack + 1800;
 Libderiv->deriv2_classes[2][2][2] = int_stack + 1836;
 Libderiv->deriv2_classes[2][2][1] = int_stack + 1872;
 Libderiv->deriv2_classes[2][2][0] = int_stack + 1908;
 memset(int_stack,0,15552);

 Libderiv->dvrr_stack = int_stack + 1944;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_d0d0(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[10] = int_stack + 36;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[9] = int_stack + 72;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[8] = int_stack + 108;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[7] = int_stack + 144;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[6] = int_stack + 180;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[2] = int_stack + 216;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[1] = int_stack + 252;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[0] = int_stack + 288;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[155] = int_stack + 324;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[143] = int_stack + 360;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[142] = int_stack + 396;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[131] = int_stack + 432;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[130] = int_stack + 468;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[129] = int_stack + 504;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[119] = int_stack + 540;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[118] = int_stack + 576;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[117] = int_stack + 612;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[116] = int_stack + 648;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[107] = int_stack + 684;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[106] = int_stack + 720;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[105] = int_stack + 756;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[104] = int_stack + 792;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[103] = int_stack + 828;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[95] = int_stack + 864;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[94] = int_stack + 900;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[93] = int_stack + 936;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[92] = int_stack + 972;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[91] = int_stack + 1008;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[90] = int_stack + 1044;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[47] = int_stack + 1080;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[46] = int_stack + 1116;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[45] = int_stack + 1152;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[44] = int_stack + 1188;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[43] = int_stack + 1224;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[42] = int_stack + 1260;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[38] = int_stack + 1296;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[35] = int_stack + 1332;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[34] = int_stack + 1368;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[33] = int_stack + 1404;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[32] = int_stack + 1440;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[31] = int_stack + 1476;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[30] = int_stack + 1512;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[26] = int_stack + 1548;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[25] = int_stack + 1584;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[23] = int_stack + 1620;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[22] = int_stack + 1656;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[21] = int_stack + 1692;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[20] = int_stack + 1728;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[19] = int_stack + 1764;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[18] = int_stack + 1800;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[14] = int_stack + 1836;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[13] = int_stack + 1872;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[12] = int_stack + 1908;

}
