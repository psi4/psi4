#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_0000(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|00) integrals */

void d12hrr_order_0000(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][0][11] = int_stack + 0;
 Libderiv->deriv_classes[0][0][10] = int_stack + 1;
 Libderiv->deriv_classes[0][0][9] = int_stack + 2;
 Libderiv->deriv_classes[0][0][8] = int_stack + 3;
 Libderiv->deriv_classes[0][0][7] = int_stack + 4;
 Libderiv->deriv_classes[0][0][6] = int_stack + 5;
 Libderiv->deriv_classes[0][0][2] = int_stack + 6;
 Libderiv->deriv_classes[0][0][1] = int_stack + 7;
 Libderiv->deriv_classes[0][0][0] = int_stack + 8;
 Libderiv->deriv2_classes[0][0][143] = int_stack + 9;
 Libderiv->deriv2_classes[0][0][131] = int_stack + 10;
 Libderiv->deriv2_classes[0][0][130] = int_stack + 11;
 Libderiv->deriv2_classes[0][0][119] = int_stack + 12;
 Libderiv->deriv2_classes[0][0][118] = int_stack + 13;
 Libderiv->deriv2_classes[0][0][117] = int_stack + 14;
 Libderiv->deriv2_classes[0][0][107] = int_stack + 15;
 Libderiv->deriv2_classes[0][0][106] = int_stack + 16;
 Libderiv->deriv2_classes[0][0][105] = int_stack + 17;
 Libderiv->deriv2_classes[0][0][104] = int_stack + 18;
 Libderiv->deriv2_classes[0][0][95] = int_stack + 19;
 Libderiv->deriv2_classes[0][0][94] = int_stack + 20;
 Libderiv->deriv2_classes[0][0][93] = int_stack + 21;
 Libderiv->deriv2_classes[0][0][92] = int_stack + 22;
 Libderiv->deriv2_classes[0][0][91] = int_stack + 23;
 Libderiv->deriv2_classes[0][0][83] = int_stack + 24;
 Libderiv->deriv2_classes[0][0][82] = int_stack + 25;
 Libderiv->deriv2_classes[0][0][81] = int_stack + 26;
 Libderiv->deriv2_classes[0][0][80] = int_stack + 27;
 Libderiv->deriv2_classes[0][0][79] = int_stack + 28;
 Libderiv->deriv2_classes[0][0][78] = int_stack + 29;
 Libderiv->deriv2_classes[0][0][35] = int_stack + 30;
 Libderiv->deriv2_classes[0][0][34] = int_stack + 31;
 Libderiv->deriv2_classes[0][0][33] = int_stack + 32;
 Libderiv->deriv2_classes[0][0][32] = int_stack + 33;
 Libderiv->deriv2_classes[0][0][31] = int_stack + 34;
 Libderiv->deriv2_classes[0][0][30] = int_stack + 35;
 Libderiv->deriv2_classes[0][0][26] = int_stack + 36;
 Libderiv->deriv2_classes[0][0][23] = int_stack + 37;
 Libderiv->deriv2_classes[0][0][22] = int_stack + 38;
 Libderiv->deriv2_classes[0][0][21] = int_stack + 39;
 Libderiv->deriv2_classes[0][0][20] = int_stack + 40;
 Libderiv->deriv2_classes[0][0][19] = int_stack + 41;
 Libderiv->deriv2_classes[0][0][18] = int_stack + 42;
 Libderiv->deriv2_classes[0][0][14] = int_stack + 43;
 Libderiv->deriv2_classes[0][0][13] = int_stack + 44;
 Libderiv->deriv2_classes[0][0][11] = int_stack + 45;
 Libderiv->deriv2_classes[0][0][10] = int_stack + 46;
 Libderiv->deriv2_classes[0][0][9] = int_stack + 47;
 Libderiv->deriv2_classes[0][0][8] = int_stack + 48;
 Libderiv->deriv2_classes[0][0][7] = int_stack + 49;
 Libderiv->deriv2_classes[0][0][6] = int_stack + 50;
 Libderiv->deriv2_classes[0][0][2] = int_stack + 51;
 Libderiv->deriv2_classes[0][0][1] = int_stack + 52;
 Libderiv->deriv2_classes[0][0][0] = int_stack + 53;
 memset(int_stack,0,432);

 Libderiv->dvrr_stack = int_stack + 54;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_0000(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|00) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[10] = int_stack + 1;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[9] = int_stack + 2;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[8] = int_stack + 3;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[7] = int_stack + 4;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[6] = int_stack + 5;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[2] = int_stack + 6;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[1] = int_stack + 7;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[0] = int_stack + 8;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[155] = int_stack + 9;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[143] = int_stack + 10;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[142] = int_stack + 11;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[131] = int_stack + 12;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[130] = int_stack + 13;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[129] = int_stack + 14;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[119] = int_stack + 15;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[118] = int_stack + 16;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[117] = int_stack + 17;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[116] = int_stack + 18;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[107] = int_stack + 19;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[106] = int_stack + 20;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[105] = int_stack + 21;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[104] = int_stack + 22;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[103] = int_stack + 23;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[95] = int_stack + 24;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[94] = int_stack + 25;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[93] = int_stack + 26;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[92] = int_stack + 27;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[91] = int_stack + 28;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[90] = int_stack + 29;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[47] = int_stack + 30;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[46] = int_stack + 31;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[45] = int_stack + 32;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[44] = int_stack + 33;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[43] = int_stack + 34;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[42] = int_stack + 35;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[38] = int_stack + 36;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[35] = int_stack + 37;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[34] = int_stack + 38;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[33] = int_stack + 39;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[32] = int_stack + 40;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[31] = int_stack + 41;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[30] = int_stack + 42;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[26] = int_stack + 43;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[25] = int_stack + 44;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[23] = int_stack + 45;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[22] = int_stack + 46;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[21] = int_stack + 47;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[20] = int_stack + 48;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[19] = int_stack + 49;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[18] = int_stack + 50;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[14] = int_stack + 51;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[13] = int_stack + 52;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[12] = int_stack + 53;

}
