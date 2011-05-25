#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0d0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|d0) integrals */

void d12hrr_order_p0d0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][2][11] = int_stack + 0;
 Libderiv->deriv_classes[1][2][10] = int_stack + 18;
 Libderiv->deriv_classes[1][2][9] = int_stack + 36;
 Libderiv->deriv_classes[1][2][8] = int_stack + 54;
 Libderiv->deriv_classes[1][2][7] = int_stack + 72;
 Libderiv->deriv_classes[1][2][6] = int_stack + 90;
 Libderiv->deriv_classes[1][2][2] = int_stack + 108;
 Libderiv->deriv_classes[1][2][1] = int_stack + 126;
 Libderiv->deriv_classes[1][2][0] = int_stack + 144;
 Libderiv->deriv2_classes[1][2][143] = int_stack + 162;
 Libderiv->deriv2_classes[1][2][131] = int_stack + 180;
 Libderiv->deriv2_classes[1][2][130] = int_stack + 198;
 Libderiv->deriv2_classes[1][2][119] = int_stack + 216;
 Libderiv->deriv2_classes[1][2][118] = int_stack + 234;
 Libderiv->deriv2_classes[1][2][117] = int_stack + 252;
 Libderiv->deriv2_classes[1][2][107] = int_stack + 270;
 Libderiv->deriv2_classes[1][2][106] = int_stack + 288;
 Libderiv->deriv2_classes[1][2][105] = int_stack + 306;
 Libderiv->deriv2_classes[1][2][104] = int_stack + 324;
 Libderiv->deriv2_classes[1][2][95] = int_stack + 342;
 Libderiv->deriv2_classes[1][2][94] = int_stack + 360;
 Libderiv->deriv2_classes[1][2][93] = int_stack + 378;
 Libderiv->deriv2_classes[1][2][92] = int_stack + 396;
 Libderiv->deriv2_classes[1][2][91] = int_stack + 414;
 Libderiv->deriv2_classes[1][2][83] = int_stack + 432;
 Libderiv->deriv2_classes[1][2][82] = int_stack + 450;
 Libderiv->deriv2_classes[1][2][81] = int_stack + 468;
 Libderiv->deriv2_classes[1][2][80] = int_stack + 486;
 Libderiv->deriv2_classes[1][2][79] = int_stack + 504;
 Libderiv->deriv2_classes[1][2][78] = int_stack + 522;
 Libderiv->deriv2_classes[1][2][35] = int_stack + 540;
 Libderiv->deriv2_classes[1][2][34] = int_stack + 558;
 Libderiv->deriv2_classes[1][2][33] = int_stack + 576;
 Libderiv->deriv2_classes[1][2][32] = int_stack + 594;
 Libderiv->deriv2_classes[1][2][31] = int_stack + 612;
 Libderiv->deriv2_classes[1][2][30] = int_stack + 630;
 Libderiv->deriv2_classes[1][2][26] = int_stack + 648;
 Libderiv->deriv2_classes[1][2][23] = int_stack + 666;
 Libderiv->deriv2_classes[1][2][22] = int_stack + 684;
 Libderiv->deriv2_classes[1][2][21] = int_stack + 702;
 Libderiv->deriv2_classes[1][2][20] = int_stack + 720;
 Libderiv->deriv2_classes[1][2][19] = int_stack + 738;
 Libderiv->deriv2_classes[1][2][18] = int_stack + 756;
 Libderiv->deriv2_classes[1][2][14] = int_stack + 774;
 Libderiv->deriv2_classes[1][2][13] = int_stack + 792;
 Libderiv->deriv2_classes[1][2][11] = int_stack + 810;
 Libderiv->deriv2_classes[1][2][10] = int_stack + 828;
 Libderiv->deriv2_classes[1][2][9] = int_stack + 846;
 Libderiv->deriv2_classes[1][2][8] = int_stack + 864;
 Libderiv->deriv2_classes[1][2][7] = int_stack + 882;
 Libderiv->deriv2_classes[1][2][6] = int_stack + 900;
 Libderiv->deriv2_classes[1][2][2] = int_stack + 918;
 Libderiv->deriv2_classes[1][2][1] = int_stack + 936;
 Libderiv->deriv2_classes[1][2][0] = int_stack + 954;
 memset(int_stack,0,7776);

 Libderiv->dvrr_stack = int_stack + 972;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0d0(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[10] = int_stack + 18;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[9] = int_stack + 36;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[8] = int_stack + 54;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[7] = int_stack + 72;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[6] = int_stack + 90;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[2] = int_stack + 108;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[1] = int_stack + 126;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[0] = int_stack + 144;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[155] = int_stack + 162;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[143] = int_stack + 180;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[142] = int_stack + 198;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[131] = int_stack + 216;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[130] = int_stack + 234;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[129] = int_stack + 252;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[119] = int_stack + 270;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[118] = int_stack + 288;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[117] = int_stack + 306;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[116] = int_stack + 324;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[107] = int_stack + 342;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[106] = int_stack + 360;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[105] = int_stack + 378;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[104] = int_stack + 396;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[103] = int_stack + 414;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[95] = int_stack + 432;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[94] = int_stack + 450;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[93] = int_stack + 468;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[92] = int_stack + 486;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[91] = int_stack + 504;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[90] = int_stack + 522;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[47] = int_stack + 540;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[46] = int_stack + 558;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[45] = int_stack + 576;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[44] = int_stack + 594;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[43] = int_stack + 612;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[42] = int_stack + 630;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[38] = int_stack + 648;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[35] = int_stack + 666;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[34] = int_stack + 684;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[33] = int_stack + 702;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[32] = int_stack + 720;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[31] = int_stack + 738;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[30] = int_stack + 756;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[26] = int_stack + 774;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[25] = int_stack + 792;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[23] = int_stack + 810;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[22] = int_stack + 828;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[21] = int_stack + 846;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[20] = int_stack + 864;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[19] = int_stack + 882;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[18] = int_stack + 900;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[14] = int_stack + 918;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[13] = int_stack + 936;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[12] = int_stack + 954;

}
