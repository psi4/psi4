#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_f0f0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|f0) integrals */

void d12hrr_order_f0f0(Libderiv_t *Libderiv, int num_prim_comb)
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
 Libderiv->deriv_classes[3][3][0] = int_stack + 800;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 900;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 1000;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 1100;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 1200;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 1300;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 1400;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 1500;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 1600;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 1700;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 1800;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 1900;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 2000;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 2100;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 2200;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 2300;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 2400;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 2500;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 2600;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 2700;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 2800;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 2900;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 3000;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 3100;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 3200;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 3300;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 3400;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 3500;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 3600;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 3700;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 3800;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 3900;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 4000;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 4100;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 4200;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 4300;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 4400;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 4500;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 4600;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 4700;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 4800;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 4900;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 5000;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 5100;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 5200;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 5300;
 memset(int_stack,0,43200);

 Libderiv->dvrr_stack = int_stack + 5400;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_f0f0(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[10] = int_stack + 100;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[9] = int_stack + 200;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[8] = int_stack + 300;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[7] = int_stack + 400;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[6] = int_stack + 500;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[2] = int_stack + 600;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[1] = int_stack + 700;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[0] = int_stack + 800;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[155] = int_stack + 900;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[143] = int_stack + 1000;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[142] = int_stack + 1100;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[131] = int_stack + 1200;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[130] = int_stack + 1300;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[129] = int_stack + 1400;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[119] = int_stack + 1500;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[118] = int_stack + 1600;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[117] = int_stack + 1700;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[116] = int_stack + 1800;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[107] = int_stack + 1900;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[106] = int_stack + 2000;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[105] = int_stack + 2100;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[104] = int_stack + 2200;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[103] = int_stack + 2300;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[95] = int_stack + 2400;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[94] = int_stack + 2500;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[93] = int_stack + 2600;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[92] = int_stack + 2700;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[91] = int_stack + 2800;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[90] = int_stack + 2900;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[47] = int_stack + 3000;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[46] = int_stack + 3100;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[45] = int_stack + 3200;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[44] = int_stack + 3300;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[43] = int_stack + 3400;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[42] = int_stack + 3500;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[38] = int_stack + 3600;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[35] = int_stack + 3700;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[34] = int_stack + 3800;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[33] = int_stack + 3900;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[32] = int_stack + 4000;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[31] = int_stack + 4100;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[30] = int_stack + 4200;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[26] = int_stack + 4300;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[25] = int_stack + 4400;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[23] = int_stack + 4500;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[22] = int_stack + 4600;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[21] = int_stack + 4700;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[20] = int_stack + 4800;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[19] = int_stack + 4900;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[18] = int_stack + 5000;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[14] = int_stack + 5100;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[13] = int_stack + 5200;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[12] = int_stack + 5300;

}
