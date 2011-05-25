#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppdp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|dp) integrals */

void d1hrr_order_ppdp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][2][11] = int_stack + 0;
 Libderiv->deriv_classes[1][3][11] = int_stack + 18;
 Libderiv->deriv_classes[2][2][11] = int_stack + 48;
 Libderiv->deriv_classes[2][3][11] = int_stack + 84;
 Libderiv->deriv_classes[1][2][10] = int_stack + 144;
 Libderiv->deriv_classes[1][3][10] = int_stack + 162;
 Libderiv->deriv_classes[2][2][10] = int_stack + 192;
 Libderiv->deriv_classes[2][3][10] = int_stack + 228;
 Libderiv->deriv_classes[1][2][9] = int_stack + 288;
 Libderiv->deriv_classes[1][3][9] = int_stack + 306;
 Libderiv->deriv_classes[2][2][9] = int_stack + 336;
 Libderiv->deriv_classes[2][3][9] = int_stack + 372;
 Libderiv->deriv_classes[1][2][8] = int_stack + 432;
 Libderiv->deriv_classes[1][3][8] = int_stack + 450;
 Libderiv->deriv_classes[2][2][8] = int_stack + 480;
 Libderiv->deriv_classes[2][3][8] = int_stack + 516;
 Libderiv->deriv_classes[1][2][7] = int_stack + 576;
 Libderiv->deriv_classes[1][3][7] = int_stack + 594;
 Libderiv->deriv_classes[2][2][7] = int_stack + 624;
 Libderiv->deriv_classes[2][3][7] = int_stack + 660;
 Libderiv->deriv_classes[1][2][6] = int_stack + 720;
 Libderiv->deriv_classes[1][3][6] = int_stack + 738;
 Libderiv->dvrr_classes[2][2] = int_stack + 768;
 Libderiv->deriv_classes[2][2][6] = int_stack + 804;
 Libderiv->deriv_classes[2][3][6] = int_stack + 840;
 Libderiv->deriv_classes[1][2][2] = int_stack + 900;
 Libderiv->deriv_classes[1][3][2] = int_stack + 918;
 Libderiv->deriv_classes[2][2][2] = int_stack + 948;
 Libderiv->deriv_classes[2][3][2] = int_stack + 984;
 Libderiv->deriv_classes[1][2][1] = int_stack + 1044;
 Libderiv->deriv_classes[1][3][1] = int_stack + 1062;
 Libderiv->deriv_classes[2][2][1] = int_stack + 1092;
 Libderiv->deriv_classes[2][3][1] = int_stack + 1128;
 Libderiv->dvrr_classes[1][2] = int_stack + 1188;
 Libderiv->dvrr_classes[1][3] = int_stack + 1206;
 Libderiv->deriv_classes[1][2][0] = int_stack + 1236;
 Libderiv->deriv_classes[1][3][0] = int_stack + 1254;
 Libderiv->deriv_classes[2][2][0] = int_stack + 1284;
 Libderiv->deriv_classes[2][3][0] = int_stack + 1320;
 memset(int_stack,0,11040);

 Libderiv->dvrr_stack = int_stack + 1812;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppdp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1380,int_stack+18,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1188,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1434,int_stack+84,int_stack+48, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+768,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+162,int_stack+144, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1188, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+228,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+768, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+162,int_stack+306,int_stack+288, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1188, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+372,int_stack+336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+768, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+324,int_stack+450,int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1542,int_stack+516,int_stack+480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+378,int_stack+594,int_stack+576, 0.0, zero_stack, 1.0, int_stack+1188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+432,int_stack+660,int_stack+624, 0.0, zero_stack, 1.0, int_stack+768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+540,int_stack+738,int_stack+720, 1.0, int_stack+1188, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+594,int_stack+840,int_stack+804, 1.0, int_stack+768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+702,int_stack+1206,int_stack+1188,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+756,int_stack+918,int_stack+900,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+810,int_stack+984,int_stack+948,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+918,int_stack+1062,int_stack+1044,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+972,int_stack+1128,int_stack+1092,6);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1080,int_stack+1254,int_stack+1236,3);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1134,int_stack+1320,int_stack+1284,6);
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1650,int_stack+1434,int_stack+1380,18);
     Libderiv->ABCD[11] = int_stack + 1650;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1242,int_stack+54,int_stack+0,18);
     Libderiv->ABCD[10] = int_stack + 1242;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+216,int_stack+162,18);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+162,int_stack+1542,int_stack+324,18);
     Libderiv->ABCD[8] = int_stack + 162;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1404,int_stack+432,int_stack+378,18);
     Libderiv->ABCD[7] = int_stack + 1404;
 /*--- compute (pp|dp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+324,int_stack+594,int_stack+540,18);
     Libderiv->ABCD[6] = int_stack + 324;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+486,int_stack+810,int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[2] = int_stack + 486;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+756,int_stack+972,int_stack+918, 0.0, zero_stack, 1.0, int_stack+702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[1] = int_stack + 756;
 /*--- compute (pp|dp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+918,int_stack+1134,int_stack+1080, 1.0, int_stack+702, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,18);
     Libderiv->ABCD[0] = int_stack + 918;

}
