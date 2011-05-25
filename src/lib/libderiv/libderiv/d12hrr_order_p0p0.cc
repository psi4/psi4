#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_p0p0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|p0) integrals */

void d12hrr_order_p0p0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][1][11] = int_stack + 0;
 Libderiv->deriv_classes[1][1][10] = int_stack + 9;
 Libderiv->deriv_classes[1][1][9] = int_stack + 18;
 Libderiv->deriv_classes[1][1][8] = int_stack + 27;
 Libderiv->deriv_classes[1][1][7] = int_stack + 36;
 Libderiv->deriv_classes[1][1][6] = int_stack + 45;
 Libderiv->deriv_classes[1][1][2] = int_stack + 54;
 Libderiv->deriv_classes[1][1][1] = int_stack + 63;
 Libderiv->deriv_classes[1][1][0] = int_stack + 72;
 Libderiv->deriv2_classes[1][1][143] = int_stack + 81;
 Libderiv->deriv2_classes[1][1][131] = int_stack + 90;
 Libderiv->deriv2_classes[1][1][130] = int_stack + 99;
 Libderiv->deriv2_classes[1][1][119] = int_stack + 108;
 Libderiv->deriv2_classes[1][1][118] = int_stack + 117;
 Libderiv->deriv2_classes[1][1][117] = int_stack + 126;
 Libderiv->deriv2_classes[1][1][107] = int_stack + 135;
 Libderiv->deriv2_classes[1][1][106] = int_stack + 144;
 Libderiv->deriv2_classes[1][1][105] = int_stack + 153;
 Libderiv->deriv2_classes[1][1][104] = int_stack + 162;
 Libderiv->deriv2_classes[1][1][95] = int_stack + 171;
 Libderiv->deriv2_classes[1][1][94] = int_stack + 180;
 Libderiv->deriv2_classes[1][1][93] = int_stack + 189;
 Libderiv->deriv2_classes[1][1][92] = int_stack + 198;
 Libderiv->deriv2_classes[1][1][91] = int_stack + 207;
 Libderiv->deriv2_classes[1][1][83] = int_stack + 216;
 Libderiv->deriv2_classes[1][1][82] = int_stack + 225;
 Libderiv->deriv2_classes[1][1][81] = int_stack + 234;
 Libderiv->deriv2_classes[1][1][80] = int_stack + 243;
 Libderiv->deriv2_classes[1][1][79] = int_stack + 252;
 Libderiv->deriv2_classes[1][1][78] = int_stack + 261;
 Libderiv->deriv2_classes[1][1][35] = int_stack + 270;
 Libderiv->deriv2_classes[1][1][34] = int_stack + 279;
 Libderiv->deriv2_classes[1][1][33] = int_stack + 288;
 Libderiv->deriv2_classes[1][1][32] = int_stack + 297;
 Libderiv->deriv2_classes[1][1][31] = int_stack + 306;
 Libderiv->deriv2_classes[1][1][30] = int_stack + 315;
 Libderiv->deriv2_classes[1][1][26] = int_stack + 324;
 Libderiv->deriv2_classes[1][1][23] = int_stack + 333;
 Libderiv->deriv2_classes[1][1][22] = int_stack + 342;
 Libderiv->deriv2_classes[1][1][21] = int_stack + 351;
 Libderiv->deriv2_classes[1][1][20] = int_stack + 360;
 Libderiv->deriv2_classes[1][1][19] = int_stack + 369;
 Libderiv->deriv2_classes[1][1][18] = int_stack + 378;
 Libderiv->deriv2_classes[1][1][14] = int_stack + 387;
 Libderiv->deriv2_classes[1][1][13] = int_stack + 396;
 Libderiv->deriv2_classes[1][1][11] = int_stack + 405;
 Libderiv->deriv2_classes[1][1][10] = int_stack + 414;
 Libderiv->deriv2_classes[1][1][9] = int_stack + 423;
 Libderiv->deriv2_classes[1][1][8] = int_stack + 432;
 Libderiv->deriv2_classes[1][1][7] = int_stack + 441;
 Libderiv->deriv2_classes[1][1][6] = int_stack + 450;
 Libderiv->deriv2_classes[1][1][2] = int_stack + 459;
 Libderiv->deriv2_classes[1][1][1] = int_stack + 468;
 Libderiv->deriv2_classes[1][1][0] = int_stack + 477;
 memset(int_stack,0,3888);

 Libderiv->dvrr_stack = int_stack + 486;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_p0p0(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[10] = int_stack + 9;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[9] = int_stack + 18;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[8] = int_stack + 27;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[7] = int_stack + 36;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[6] = int_stack + 45;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[2] = int_stack + 54;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[1] = int_stack + 63;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[0] = int_stack + 72;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[155] = int_stack + 81;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[143] = int_stack + 90;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[142] = int_stack + 99;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[131] = int_stack + 108;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[130] = int_stack + 117;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[129] = int_stack + 126;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[119] = int_stack + 135;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[118] = int_stack + 144;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[117] = int_stack + 153;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[116] = int_stack + 162;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[107] = int_stack + 171;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[106] = int_stack + 180;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[105] = int_stack + 189;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[104] = int_stack + 198;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[103] = int_stack + 207;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[95] = int_stack + 216;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[94] = int_stack + 225;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[93] = int_stack + 234;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[92] = int_stack + 243;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[91] = int_stack + 252;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[90] = int_stack + 261;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[47] = int_stack + 270;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[46] = int_stack + 279;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[45] = int_stack + 288;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[44] = int_stack + 297;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[43] = int_stack + 306;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[42] = int_stack + 315;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[38] = int_stack + 324;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[35] = int_stack + 333;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[34] = int_stack + 342;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[33] = int_stack + 351;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[32] = int_stack + 360;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[31] = int_stack + 369;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[30] = int_stack + 378;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[26] = int_stack + 387;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[25] = int_stack + 396;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[23] = int_stack + 405;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[22] = int_stack + 414;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[21] = int_stack + 423;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[20] = int_stack + 432;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[19] = int_stack + 441;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[18] = int_stack + 450;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[14] = int_stack + 459;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[13] = int_stack + 468;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[12] = int_stack + 477;

}
