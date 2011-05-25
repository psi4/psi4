#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00p0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|p0) integrals */

void d12hrr_order_00p0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][1][11] = int_stack + 0;
 Libderiv->deriv_classes[0][1][10] = int_stack + 3;
 Libderiv->deriv_classes[0][1][9] = int_stack + 6;
 Libderiv->deriv_classes[0][1][8] = int_stack + 9;
 Libderiv->deriv_classes[0][1][7] = int_stack + 12;
 Libderiv->deriv_classes[0][1][6] = int_stack + 15;
 Libderiv->deriv_classes[0][1][2] = int_stack + 18;
 Libderiv->deriv_classes[0][1][1] = int_stack + 21;
 Libderiv->deriv_classes[0][1][0] = int_stack + 24;
 Libderiv->deriv2_classes[0][1][143] = int_stack + 27;
 Libderiv->deriv2_classes[0][1][131] = int_stack + 30;
 Libderiv->deriv2_classes[0][1][130] = int_stack + 33;
 Libderiv->deriv2_classes[0][1][119] = int_stack + 36;
 Libderiv->deriv2_classes[0][1][118] = int_stack + 39;
 Libderiv->deriv2_classes[0][1][117] = int_stack + 42;
 Libderiv->deriv2_classes[0][1][107] = int_stack + 45;
 Libderiv->deriv2_classes[0][1][106] = int_stack + 48;
 Libderiv->deriv2_classes[0][1][105] = int_stack + 51;
 Libderiv->deriv2_classes[0][1][104] = int_stack + 54;
 Libderiv->deriv2_classes[0][1][95] = int_stack + 57;
 Libderiv->deriv2_classes[0][1][94] = int_stack + 60;
 Libderiv->deriv2_classes[0][1][93] = int_stack + 63;
 Libderiv->deriv2_classes[0][1][92] = int_stack + 66;
 Libderiv->deriv2_classes[0][1][91] = int_stack + 69;
 Libderiv->deriv2_classes[0][1][83] = int_stack + 72;
 Libderiv->deriv2_classes[0][1][82] = int_stack + 75;
 Libderiv->deriv2_classes[0][1][81] = int_stack + 78;
 Libderiv->deriv2_classes[0][1][80] = int_stack + 81;
 Libderiv->deriv2_classes[0][1][79] = int_stack + 84;
 Libderiv->deriv2_classes[0][1][78] = int_stack + 87;
 Libderiv->deriv2_classes[0][1][35] = int_stack + 90;
 Libderiv->deriv2_classes[0][1][34] = int_stack + 93;
 Libderiv->deriv2_classes[0][1][33] = int_stack + 96;
 Libderiv->deriv2_classes[0][1][32] = int_stack + 99;
 Libderiv->deriv2_classes[0][1][31] = int_stack + 102;
 Libderiv->deriv2_classes[0][1][30] = int_stack + 105;
 Libderiv->deriv2_classes[0][1][26] = int_stack + 108;
 Libderiv->deriv2_classes[0][1][23] = int_stack + 111;
 Libderiv->deriv2_classes[0][1][22] = int_stack + 114;
 Libderiv->deriv2_classes[0][1][21] = int_stack + 117;
 Libderiv->deriv2_classes[0][1][20] = int_stack + 120;
 Libderiv->deriv2_classes[0][1][19] = int_stack + 123;
 Libderiv->deriv2_classes[0][1][18] = int_stack + 126;
 Libderiv->deriv2_classes[0][1][14] = int_stack + 129;
 Libderiv->deriv2_classes[0][1][13] = int_stack + 132;
 Libderiv->deriv2_classes[0][1][11] = int_stack + 135;
 Libderiv->deriv2_classes[0][1][10] = int_stack + 138;
 Libderiv->deriv2_classes[0][1][9] = int_stack + 141;
 Libderiv->deriv2_classes[0][1][8] = int_stack + 144;
 Libderiv->deriv2_classes[0][1][7] = int_stack + 147;
 Libderiv->deriv2_classes[0][1][6] = int_stack + 150;
 Libderiv->deriv2_classes[0][1][2] = int_stack + 153;
 Libderiv->deriv2_classes[0][1][1] = int_stack + 156;
 Libderiv->deriv2_classes[0][1][0] = int_stack + 159;
 memset(int_stack,0,1296);

 Libderiv->dvrr_stack = int_stack + 162;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00p0(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[10] = int_stack + 3;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[9] = int_stack + 6;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[8] = int_stack + 9;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[7] = int_stack + 12;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[6] = int_stack + 15;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[2] = int_stack + 18;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[1] = int_stack + 21;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[0] = int_stack + 24;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[155] = int_stack + 27;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[143] = int_stack + 30;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[142] = int_stack + 33;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[131] = int_stack + 36;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[130] = int_stack + 39;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[129] = int_stack + 42;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[119] = int_stack + 45;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[118] = int_stack + 48;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[117] = int_stack + 51;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[116] = int_stack + 54;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[107] = int_stack + 57;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[106] = int_stack + 60;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[105] = int_stack + 63;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[104] = int_stack + 66;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[103] = int_stack + 69;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[95] = int_stack + 72;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[94] = int_stack + 75;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[93] = int_stack + 78;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[92] = int_stack + 81;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[91] = int_stack + 84;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[90] = int_stack + 87;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[47] = int_stack + 90;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[46] = int_stack + 93;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[45] = int_stack + 96;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[44] = int_stack + 99;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[43] = int_stack + 102;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[42] = int_stack + 105;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[38] = int_stack + 108;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[35] = int_stack + 111;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[34] = int_stack + 114;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[33] = int_stack + 117;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[32] = int_stack + 120;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[31] = int_stack + 123;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[30] = int_stack + 126;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[26] = int_stack + 129;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[25] = int_stack + 132;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[23] = int_stack + 135;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[22] = int_stack + 138;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[21] = int_stack + 141;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[20] = int_stack + 144;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[19] = int_stack + 147;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[18] = int_stack + 150;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[14] = int_stack + 153;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[13] = int_stack + 156;
 /*--- compute (00|p0) ---*/
     Libderiv->ABCD[12] = int_stack + 159;

}
