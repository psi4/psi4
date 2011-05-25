#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d12vrr_order_f0fp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|fp) integrals */

void d12hrr_order_f0fp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][10] = int_stack + 150;
 Libderiv->deriv_classes[3][4][9] = int_stack + 300;
 Libderiv->deriv_classes[3][4][8] = int_stack + 450;
 Libderiv->deriv_classes[3][4][7] = int_stack + 600;
 Libderiv->dvrr_classes[3][3] = int_stack + 750;
 Libderiv->deriv_classes[3][4][6] = int_stack + 850;
 Libderiv->deriv_classes[3][4][2] = int_stack + 1000;
 Libderiv->deriv_classes[3][4][1] = int_stack + 1150;
 Libderiv->deriv_classes[3][4][0] = int_stack + 1300;
 Libderiv->deriv2_classes[3][3][143] = int_stack + 1450;
 Libderiv->deriv2_classes[3][4][143] = int_stack + 1550;
 Libderiv->deriv2_classes[3][3][131] = int_stack + 1700;
 Libderiv->deriv2_classes[3][4][131] = int_stack + 1800;
 Libderiv->deriv2_classes[3][3][130] = int_stack + 1950;
 Libderiv->deriv2_classes[3][4][130] = int_stack + 2050;
 Libderiv->deriv2_classes[3][3][119] = int_stack + 2200;
 Libderiv->deriv2_classes[3][4][119] = int_stack + 2300;
 Libderiv->deriv2_classes[3][3][118] = int_stack + 2450;
 Libderiv->deriv2_classes[3][4][118] = int_stack + 2550;
 Libderiv->deriv2_classes[3][3][117] = int_stack + 2700;
 Libderiv->deriv2_classes[3][4][117] = int_stack + 2800;
 Libderiv->deriv2_classes[3][3][107] = int_stack + 2950;
 Libderiv->deriv2_classes[3][4][107] = int_stack + 3050;
 Libderiv->deriv2_classes[3][3][106] = int_stack + 3200;
 Libderiv->deriv2_classes[3][4][106] = int_stack + 3300;
 Libderiv->deriv2_classes[3][3][105] = int_stack + 3450;
 Libderiv->deriv2_classes[3][4][105] = int_stack + 3550;
 Libderiv->deriv2_classes[3][3][104] = int_stack + 3700;
 Libderiv->deriv2_classes[3][4][104] = int_stack + 3800;
 Libderiv->deriv2_classes[3][3][95] = int_stack + 3950;
 Libderiv->deriv2_classes[3][4][95] = int_stack + 4050;
 Libderiv->deriv2_classes[3][3][94] = int_stack + 4200;
 Libderiv->deriv2_classes[3][4][94] = int_stack + 4300;
 Libderiv->deriv2_classes[3][3][93] = int_stack + 4450;
 Libderiv->deriv2_classes[3][4][93] = int_stack + 4550;
 Libderiv->deriv2_classes[3][3][92] = int_stack + 4700;
 Libderiv->deriv2_classes[3][4][92] = int_stack + 4800;
 Libderiv->deriv2_classes[3][3][91] = int_stack + 4950;
 Libderiv->deriv2_classes[3][4][91] = int_stack + 5050;
 Libderiv->deriv_classes[3][3][11] = int_stack + 5200;
 Libderiv->deriv2_classes[3][3][83] = int_stack + 5300;
 Libderiv->deriv2_classes[3][4][83] = int_stack + 5400;
 Libderiv->deriv_classes[3][3][10] = int_stack + 5550;
 Libderiv->deriv2_classes[3][3][82] = int_stack + 5650;
 Libderiv->deriv2_classes[3][4][82] = int_stack + 5750;
 Libderiv->deriv_classes[3][3][9] = int_stack + 5900;
 Libderiv->deriv2_classes[3][3][81] = int_stack + 6000;
 Libderiv->deriv2_classes[3][4][81] = int_stack + 6100;
 Libderiv->deriv_classes[3][3][8] = int_stack + 6250;
 Libderiv->deriv2_classes[3][3][80] = int_stack + 6350;
 Libderiv->deriv2_classes[3][4][80] = int_stack + 6450;
 Libderiv->deriv_classes[3][3][7] = int_stack + 6600;
 Libderiv->deriv2_classes[3][3][79] = int_stack + 6700;
 Libderiv->deriv2_classes[3][4][79] = int_stack + 6800;
 Libderiv->deriv_classes[3][3][6] = int_stack + 6950;
 Libderiv->deriv2_classes[3][3][78] = int_stack + 7050;
 Libderiv->deriv2_classes[3][4][78] = int_stack + 7150;
 Libderiv->deriv2_classes[3][3][35] = int_stack + 7300;
 Libderiv->deriv2_classes[3][4][35] = int_stack + 7400;
 Libderiv->deriv2_classes[3][3][34] = int_stack + 7550;
 Libderiv->deriv2_classes[3][4][34] = int_stack + 7650;
 Libderiv->deriv2_classes[3][3][33] = int_stack + 7800;
 Libderiv->deriv2_classes[3][4][33] = int_stack + 7900;
 Libderiv->deriv2_classes[3][3][32] = int_stack + 8050;
 Libderiv->deriv2_classes[3][4][32] = int_stack + 8150;
 Libderiv->deriv2_classes[3][3][31] = int_stack + 8300;
 Libderiv->deriv2_classes[3][4][31] = int_stack + 8400;
 Libderiv->deriv_classes[3][3][2] = int_stack + 8550;
 Libderiv->deriv2_classes[3][3][30] = int_stack + 8650;
 Libderiv->deriv2_classes[3][4][30] = int_stack + 8750;
 Libderiv->deriv2_classes[3][3][26] = int_stack + 8900;
 Libderiv->deriv2_classes[3][4][26] = int_stack + 9000;
 Libderiv->deriv2_classes[3][3][23] = int_stack + 9150;
 Libderiv->deriv2_classes[3][4][23] = int_stack + 9250;
 Libderiv->deriv2_classes[3][3][22] = int_stack + 9400;
 Libderiv->deriv2_classes[3][4][22] = int_stack + 9500;
 Libderiv->deriv2_classes[3][3][21] = int_stack + 9650;
 Libderiv->deriv2_classes[3][4][21] = int_stack + 9750;
 Libderiv->deriv2_classes[3][3][20] = int_stack + 9900;
 Libderiv->deriv2_classes[3][4][20] = int_stack + 10000;
 Libderiv->deriv2_classes[3][3][19] = int_stack + 10150;
 Libderiv->deriv2_classes[3][4][19] = int_stack + 10250;
 Libderiv->deriv_classes[3][3][1] = int_stack + 10400;
 Libderiv->deriv2_classes[3][3][18] = int_stack + 10500;
 Libderiv->deriv2_classes[3][4][18] = int_stack + 10600;
 Libderiv->deriv2_classes[3][3][14] = int_stack + 10750;
 Libderiv->deriv2_classes[3][4][14] = int_stack + 10850;
 Libderiv->deriv2_classes[3][3][13] = int_stack + 11000;
 Libderiv->deriv2_classes[3][4][13] = int_stack + 11100;
 Libderiv->deriv2_classes[3][3][11] = int_stack + 11250;
 Libderiv->deriv2_classes[3][4][11] = int_stack + 11350;
 Libderiv->deriv2_classes[3][3][10] = int_stack + 11500;
 Libderiv->deriv2_classes[3][4][10] = int_stack + 11600;
 Libderiv->deriv2_classes[3][3][9] = int_stack + 11750;
 Libderiv->deriv2_classes[3][4][9] = int_stack + 11850;
 Libderiv->deriv2_classes[3][3][8] = int_stack + 12000;
 Libderiv->deriv2_classes[3][4][8] = int_stack + 12100;
 Libderiv->deriv2_classes[3][3][7] = int_stack + 12250;
 Libderiv->deriv2_classes[3][4][7] = int_stack + 12350;
 Libderiv->deriv_classes[3][3][0] = int_stack + 12500;
 Libderiv->deriv2_classes[3][3][6] = int_stack + 12600;
 Libderiv->deriv2_classes[3][4][6] = int_stack + 12700;
 Libderiv->deriv2_classes[3][3][2] = int_stack + 12850;
 Libderiv->deriv2_classes[3][4][2] = int_stack + 12950;
 Libderiv->deriv2_classes[3][3][1] = int_stack + 13100;
 Libderiv->deriv2_classes[3][4][1] = int_stack + 13200;
 Libderiv->deriv2_classes[3][3][0] = int_stack + 13350;
 Libderiv->deriv2_classes[3][4][0] = int_stack + 13450;
 memset(int_stack,0,108800);

 Libderiv->dvrr_stack = int_stack + 16600;
 for(i=0;i<num_prim_comb;i++) {
   d12vrr_order_f0fp(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13600,int_stack+0,int_stack+5200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750,10);
     Libderiv->ABCD[11] = int_stack + 13600;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+13900,int_stack+150,int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 13900;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+300,int_stack+5900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14200,int_stack+450,int_stack+6250, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 14200;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+600,int_stack+6600, 0.0, zero_stack, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 300;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+14500,int_stack+850,int_stack+6950, 1.0, int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 14500;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+1000,int_stack+8550,10);
     Libderiv->ABCD[2] = int_stack + 600;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+14800,int_stack+1150,int_stack+10400,10);
     Libderiv->ABCD[1] = int_stack + 14800;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1300,int_stack+12500,10);
     Libderiv->ABCD[0] = int_stack + 900;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15100,int_stack+1550,int_stack+1450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5200,10);
     Libderiv->ABCD[155] = int_stack + 15100;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+1800,int_stack+1700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5200, 1.0, int_stack+5550,10);
     Libderiv->ABCD[143] = int_stack + 1200;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1500,int_stack+2050,int_stack+1950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5550, 0.0, zero_stack,10);
     Libderiv->ABCD[142] = int_stack + 1500;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2300,int_stack+2200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5200, 0.0, zero_stack, 1.0, int_stack+5900,10);
     Libderiv->ABCD[131] = int_stack + 1800;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2100,int_stack+2550,int_stack+2450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5550, 1.0, int_stack+5900, 0.0, zero_stack,10);
     Libderiv->ABCD[130] = int_stack + 2100;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2400,int_stack+2800,int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+5900, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[129] = int_stack + 2400;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15400,int_stack+3050,int_stack+2950, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6250,10);
     Libderiv->ABCD[119] = int_stack + 15400;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+3300,int_stack+3200, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5550, 0.0, zero_stack, 1.0, int_stack+6250, 0.0, zero_stack,10);
     Libderiv->ABCD[118] = int_stack + 2700;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3000,int_stack+3550,int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5900, 1.0, int_stack+6250, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[117] = int_stack + 3000;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3300,int_stack+3800,int_stack+3700, 0.0, zero_stack, 0.0, zero_stack, 2.0, int_stack+6250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[116] = int_stack + 3300;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4050,int_stack+3950, 0.0, zero_stack, 1.0, int_stack+5200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6600,10);
     Libderiv->ABCD[107] = int_stack + 3600;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3900,int_stack+4300,int_stack+4200, 0.0, zero_stack, 1.0, int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6600, 0.0, zero_stack,10);
     Libderiv->ABCD[106] = int_stack + 3900;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+15700,int_stack+4550,int_stack+4450, 0.0, zero_stack, 1.0, int_stack+5900, 0.0, zero_stack, 1.0, int_stack+6600, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[105] = int_stack + 15700;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4200,int_stack+4800,int_stack+4700, 0.0, zero_stack, 1.0, int_stack+6250, 1.0, int_stack+6600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[104] = int_stack + 4200;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4500,int_stack+5050,int_stack+4950, 0.0, zero_stack, 2.0, int_stack+6600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[103] = int_stack + 4500;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4800,int_stack+5400,int_stack+5300, 1.0, int_stack+5200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6950,10);
     Libderiv->ABCD[95] = int_stack + 4800;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5100,int_stack+5750,int_stack+5650, 1.0, int_stack+5550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6950, 0.0, zero_stack,10);
     Libderiv->ABCD[94] = int_stack + 5100;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5400,int_stack+6100,int_stack+6000, 1.0, int_stack+5900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6950, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[93] = int_stack + 5400;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+5700,int_stack+6450,int_stack+6350, 1.0, int_stack+6250, 0.0, zero_stack, 1.0, int_stack+6950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[92] = int_stack + 5700;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6000,int_stack+6800,int_stack+6700, 1.0, int_stack+6600, 1.0, int_stack+6950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[91] = int_stack + 6000;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6300,int_stack+7150,int_stack+7050, 2.0, int_stack+6950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[90] = int_stack + 6300;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6600,int_stack+7400,int_stack+7300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8550,10);
     Libderiv->ABCD[47] = int_stack + 6600;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+6900,int_stack+7650,int_stack+7550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8550, 0.0, zero_stack,10);
     Libderiv->ABCD[46] = int_stack + 6900;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7200,int_stack+7900,int_stack+7800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8550, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[45] = int_stack + 7200;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7500,int_stack+8150,int_stack+8050, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[44] = int_stack + 7500;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7800,int_stack+8400,int_stack+8300, 0.0, zero_stack, 1.0, int_stack+8550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[43] = int_stack + 7800;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8100,int_stack+8750,int_stack+8650, 1.0, int_stack+8550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[42] = int_stack + 8100;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+8400,int_stack+9000,int_stack+8900,10);
     Libderiv->ABCD[38] = int_stack + 8400;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+8700,int_stack+9250,int_stack+9150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10400,10);
     Libderiv->ABCD[35] = int_stack + 8700;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9000,int_stack+9500,int_stack+9400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10400, 0.0, zero_stack,10);
     Libderiv->ABCD[34] = int_stack + 9000;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9300,int_stack+9750,int_stack+9650, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10400, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[33] = int_stack + 9300;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9600,int_stack+10000,int_stack+9900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[32] = int_stack + 9600;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16000,int_stack+10250,int_stack+10150, 0.0, zero_stack, 1.0, int_stack+10400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[31] = int_stack + 16000;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+9900,int_stack+10600,int_stack+10500, 1.0, int_stack+10400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[30] = int_stack + 9900;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10200,int_stack+10850,int_stack+10750,10);
     Libderiv->ABCD[26] = int_stack + 10200;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+10500,int_stack+11100,int_stack+11000,10);
     Libderiv->ABCD[25] = int_stack + 10500;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+10800,int_stack+11350,int_stack+11250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12500,10);
     Libderiv->ABCD[23] = int_stack + 10800;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11100,int_stack+11600,int_stack+11500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12500, 0.0, zero_stack,10);
     Libderiv->ABCD[22] = int_stack + 11100;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11400,int_stack+11850,int_stack+11750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12500, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[21] = int_stack + 11400;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+11700,int_stack+12100,int_stack+12000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+12500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[20] = int_stack + 11700;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+16300,int_stack+12350,int_stack+12250, 0.0, zero_stack, 1.0, int_stack+12500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[19] = int_stack + 16300;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+12000,int_stack+12700,int_stack+12600, 1.0, int_stack+12500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[18] = int_stack + 12000;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12300,int_stack+12950,int_stack+12850,10);
     Libderiv->ABCD[14] = int_stack + 12300;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12600,int_stack+13200,int_stack+13100,10);
     Libderiv->ABCD[13] = int_stack + 12600;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+12900,int_stack+13450,int_stack+13350,10);
     Libderiv->ABCD[12] = int_stack + 12900;

}
