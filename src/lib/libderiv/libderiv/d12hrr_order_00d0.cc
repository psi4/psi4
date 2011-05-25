#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00d0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|d0) integrals */

void d12hrr_order_00d0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][2][11] = int_stack + 0;
 Libderiv->deriv_classes[0][2][10] = int_stack + 6;
 Libderiv->deriv_classes[0][2][9] = int_stack + 12;
 Libderiv->deriv_classes[0][2][8] = int_stack + 18;
 Libderiv->deriv_classes[0][2][7] = int_stack + 24;
 Libderiv->deriv_classes[0][2][6] = int_stack + 30;
 Libderiv->deriv_classes[0][2][2] = int_stack + 36;
 Libderiv->deriv_classes[0][2][1] = int_stack + 42;
 Libderiv->deriv_classes[0][2][0] = int_stack + 48;
 Libderiv->deriv2_classes[0][2][143] = int_stack + 54;
 Libderiv->deriv2_classes[0][2][131] = int_stack + 60;
 Libderiv->deriv2_classes[0][2][130] = int_stack + 66;
 Libderiv->deriv2_classes[0][2][119] = int_stack + 72;
 Libderiv->deriv2_classes[0][2][118] = int_stack + 78;
 Libderiv->deriv2_classes[0][2][117] = int_stack + 84;
 Libderiv->deriv2_classes[0][2][107] = int_stack + 90;
 Libderiv->deriv2_classes[0][2][106] = int_stack + 96;
 Libderiv->deriv2_classes[0][2][105] = int_stack + 102;
 Libderiv->deriv2_classes[0][2][104] = int_stack + 108;
 Libderiv->deriv2_classes[0][2][95] = int_stack + 114;
 Libderiv->deriv2_classes[0][2][94] = int_stack + 120;
 Libderiv->deriv2_classes[0][2][93] = int_stack + 126;
 Libderiv->deriv2_classes[0][2][92] = int_stack + 132;
 Libderiv->deriv2_classes[0][2][91] = int_stack + 138;
 Libderiv->deriv2_classes[0][2][83] = int_stack + 144;
 Libderiv->deriv2_classes[0][2][82] = int_stack + 150;
 Libderiv->deriv2_classes[0][2][81] = int_stack + 156;
 Libderiv->deriv2_classes[0][2][80] = int_stack + 162;
 Libderiv->deriv2_classes[0][2][79] = int_stack + 168;
 Libderiv->deriv2_classes[0][2][78] = int_stack + 174;
 Libderiv->deriv2_classes[0][2][35] = int_stack + 180;
 Libderiv->deriv2_classes[0][2][34] = int_stack + 186;
 Libderiv->deriv2_classes[0][2][33] = int_stack + 192;
 Libderiv->deriv2_classes[0][2][32] = int_stack + 198;
 Libderiv->deriv2_classes[0][2][31] = int_stack + 204;
 Libderiv->deriv2_classes[0][2][30] = int_stack + 210;
 Libderiv->deriv2_classes[0][2][26] = int_stack + 216;
 Libderiv->deriv2_classes[0][2][23] = int_stack + 222;
 Libderiv->deriv2_classes[0][2][22] = int_stack + 228;
 Libderiv->deriv2_classes[0][2][21] = int_stack + 234;
 Libderiv->deriv2_classes[0][2][20] = int_stack + 240;
 Libderiv->deriv2_classes[0][2][19] = int_stack + 246;
 Libderiv->deriv2_classes[0][2][18] = int_stack + 252;
 Libderiv->deriv2_classes[0][2][14] = int_stack + 258;
 Libderiv->deriv2_classes[0][2][13] = int_stack + 264;
 Libderiv->deriv2_classes[0][2][11] = int_stack + 270;
 Libderiv->deriv2_classes[0][2][10] = int_stack + 276;
 Libderiv->deriv2_classes[0][2][9] = int_stack + 282;
 Libderiv->deriv2_classes[0][2][8] = int_stack + 288;
 Libderiv->deriv2_classes[0][2][7] = int_stack + 294;
 Libderiv->deriv2_classes[0][2][6] = int_stack + 300;
 Libderiv->deriv2_classes[0][2][2] = int_stack + 306;
 Libderiv->deriv2_classes[0][2][1] = int_stack + 312;
 Libderiv->deriv2_classes[0][2][0] = int_stack + 318;
 memset(int_stack,0,2592);

 Libderiv->dvrr_stack = int_stack + 324;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00d0(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[10] = int_stack + 6;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[9] = int_stack + 12;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[8] = int_stack + 18;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[7] = int_stack + 24;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[6] = int_stack + 30;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[2] = int_stack + 36;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[1] = int_stack + 42;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[0] = int_stack + 48;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[155] = int_stack + 54;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[143] = int_stack + 60;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[142] = int_stack + 66;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[131] = int_stack + 72;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[130] = int_stack + 78;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[129] = int_stack + 84;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[119] = int_stack + 90;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[118] = int_stack + 96;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[117] = int_stack + 102;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[116] = int_stack + 108;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[107] = int_stack + 114;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[106] = int_stack + 120;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[105] = int_stack + 126;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[104] = int_stack + 132;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[103] = int_stack + 138;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[95] = int_stack + 144;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[94] = int_stack + 150;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[93] = int_stack + 156;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[92] = int_stack + 162;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[91] = int_stack + 168;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[90] = int_stack + 174;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[47] = int_stack + 180;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[46] = int_stack + 186;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[45] = int_stack + 192;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[44] = int_stack + 198;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[43] = int_stack + 204;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[42] = int_stack + 210;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[38] = int_stack + 216;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[35] = int_stack + 222;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[34] = int_stack + 228;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[33] = int_stack + 234;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[32] = int_stack + 240;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[31] = int_stack + 246;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[30] = int_stack + 252;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[26] = int_stack + 258;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[25] = int_stack + 264;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[23] = int_stack + 270;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[22] = int_stack + 276;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[21] = int_stack + 282;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[20] = int_stack + 288;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[19] = int_stack + 294;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[18] = int_stack + 300;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[14] = int_stack + 306;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[13] = int_stack + 312;
 /*--- compute (00|d0) ---*/
     Libderiv->ABCD[12] = int_stack + 318;

}
