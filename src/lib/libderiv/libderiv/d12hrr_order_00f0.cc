#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_00f0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|f0) integrals */

void d12hrr_order_00f0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][3][11] = int_stack + 0;
 Libderiv->deriv_classes[0][3][10] = int_stack + 10;
 Libderiv->deriv_classes[0][3][9] = int_stack + 20;
 Libderiv->deriv_classes[0][3][8] = int_stack + 30;
 Libderiv->deriv_classes[0][3][7] = int_stack + 40;
 Libderiv->deriv_classes[0][3][6] = int_stack + 50;
 Libderiv->deriv_classes[0][3][2] = int_stack + 60;
 Libderiv->deriv_classes[0][3][1] = int_stack + 70;
 Libderiv->deriv_classes[0][3][0] = int_stack + 80;
 Libderiv->deriv2_classes[0][3][143] = int_stack + 90;
 Libderiv->deriv2_classes[0][3][131] = int_stack + 100;
 Libderiv->deriv2_classes[0][3][130] = int_stack + 110;
 Libderiv->deriv2_classes[0][3][119] = int_stack + 120;
 Libderiv->deriv2_classes[0][3][118] = int_stack + 130;
 Libderiv->deriv2_classes[0][3][117] = int_stack + 140;
 Libderiv->deriv2_classes[0][3][107] = int_stack + 150;
 Libderiv->deriv2_classes[0][3][106] = int_stack + 160;
 Libderiv->deriv2_classes[0][3][105] = int_stack + 170;
 Libderiv->deriv2_classes[0][3][104] = int_stack + 180;
 Libderiv->deriv2_classes[0][3][95] = int_stack + 190;
 Libderiv->deriv2_classes[0][3][94] = int_stack + 200;
 Libderiv->deriv2_classes[0][3][93] = int_stack + 210;
 Libderiv->deriv2_classes[0][3][92] = int_stack + 220;
 Libderiv->deriv2_classes[0][3][91] = int_stack + 230;
 Libderiv->deriv2_classes[0][3][83] = int_stack + 240;
 Libderiv->deriv2_classes[0][3][82] = int_stack + 250;
 Libderiv->deriv2_classes[0][3][81] = int_stack + 260;
 Libderiv->deriv2_classes[0][3][80] = int_stack + 270;
 Libderiv->deriv2_classes[0][3][79] = int_stack + 280;
 Libderiv->deriv2_classes[0][3][78] = int_stack + 290;
 Libderiv->deriv2_classes[0][3][35] = int_stack + 300;
 Libderiv->deriv2_classes[0][3][34] = int_stack + 310;
 Libderiv->deriv2_classes[0][3][33] = int_stack + 320;
 Libderiv->deriv2_classes[0][3][32] = int_stack + 330;
 Libderiv->deriv2_classes[0][3][31] = int_stack + 340;
 Libderiv->deriv2_classes[0][3][30] = int_stack + 350;
 Libderiv->deriv2_classes[0][3][26] = int_stack + 360;
 Libderiv->deriv2_classes[0][3][23] = int_stack + 370;
 Libderiv->deriv2_classes[0][3][22] = int_stack + 380;
 Libderiv->deriv2_classes[0][3][21] = int_stack + 390;
 Libderiv->deriv2_classes[0][3][20] = int_stack + 400;
 Libderiv->deriv2_classes[0][3][19] = int_stack + 410;
 Libderiv->deriv2_classes[0][3][18] = int_stack + 420;
 Libderiv->deriv2_classes[0][3][14] = int_stack + 430;
 Libderiv->deriv2_classes[0][3][13] = int_stack + 440;
 Libderiv->deriv2_classes[0][3][11] = int_stack + 450;
 Libderiv->deriv2_classes[0][3][10] = int_stack + 460;
 Libderiv->deriv2_classes[0][3][9] = int_stack + 470;
 Libderiv->deriv2_classes[0][3][8] = int_stack + 480;
 Libderiv->deriv2_classes[0][3][7] = int_stack + 490;
 Libderiv->deriv2_classes[0][3][6] = int_stack + 500;
 Libderiv->deriv2_classes[0][3][2] = int_stack + 510;
 Libderiv->deriv2_classes[0][3][1] = int_stack + 520;
 Libderiv->deriv2_classes[0][3][0] = int_stack + 530;
 memset(int_stack,0,4320);

 Libderiv->dvrr_stack = int_stack + 540;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_00f0(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[10] = int_stack + 10;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[9] = int_stack + 20;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[8] = int_stack + 30;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[7] = int_stack + 40;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[6] = int_stack + 50;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[2] = int_stack + 60;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[1] = int_stack + 70;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[0] = int_stack + 80;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[155] = int_stack + 90;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[143] = int_stack + 100;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[142] = int_stack + 110;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[131] = int_stack + 120;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[130] = int_stack + 130;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[129] = int_stack + 140;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[119] = int_stack + 150;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[118] = int_stack + 160;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[117] = int_stack + 170;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[116] = int_stack + 180;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[107] = int_stack + 190;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[106] = int_stack + 200;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[105] = int_stack + 210;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[104] = int_stack + 220;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[103] = int_stack + 230;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[95] = int_stack + 240;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[94] = int_stack + 250;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[93] = int_stack + 260;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[92] = int_stack + 270;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[91] = int_stack + 280;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[90] = int_stack + 290;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[47] = int_stack + 300;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[46] = int_stack + 310;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[45] = int_stack + 320;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[44] = int_stack + 330;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[43] = int_stack + 340;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[42] = int_stack + 350;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[38] = int_stack + 360;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[35] = int_stack + 370;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[34] = int_stack + 380;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[33] = int_stack + 390;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[32] = int_stack + 400;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[31] = int_stack + 410;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[30] = int_stack + 420;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[26] = int_stack + 430;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[25] = int_stack + 440;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[23] = int_stack + 450;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[22] = int_stack + 460;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[21] = int_stack + 470;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[20] = int_stack + 480;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[19] = int_stack + 490;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[18] = int_stack + 500;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[14] = int_stack + 510;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[13] = int_stack + 520;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[12] = int_stack + 530;

}
