#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|dd) integrals */

void d1hrr_order_p0dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][2][11] = int_stack + 0;
 Libderiv->deriv_classes[1][3][11] = int_stack + 18;
 Libderiv->deriv_classes[1][4][11] = int_stack + 48;
 Libderiv->deriv_classes[1][2][10] = int_stack + 93;
 Libderiv->deriv_classes[1][3][10] = int_stack + 111;
 Libderiv->deriv_classes[1][4][10] = int_stack + 141;
 Libderiv->deriv_classes[1][2][9] = int_stack + 186;
 Libderiv->deriv_classes[1][3][9] = int_stack + 204;
 Libderiv->deriv_classes[1][4][9] = int_stack + 234;
 Libderiv->deriv_classes[1][2][8] = int_stack + 279;
 Libderiv->deriv_classes[1][3][8] = int_stack + 297;
 Libderiv->deriv_classes[1][4][8] = int_stack + 327;
 Libderiv->deriv_classes[1][2][7] = int_stack + 372;
 Libderiv->deriv_classes[1][3][7] = int_stack + 390;
 Libderiv->deriv_classes[1][4][7] = int_stack + 420;
 Libderiv->dvrr_classes[1][2] = int_stack + 465;
 Libderiv->deriv_classes[1][2][6] = int_stack + 483;
 Libderiv->dvrr_classes[1][3] = int_stack + 501;
 Libderiv->deriv_classes[1][3][6] = int_stack + 531;
 Libderiv->deriv_classes[1][4][6] = int_stack + 561;
 Libderiv->deriv_classes[1][2][2] = int_stack + 606;
 Libderiv->deriv_classes[1][3][2] = int_stack + 624;
 Libderiv->deriv_classes[1][4][2] = int_stack + 654;
 Libderiv->deriv_classes[1][2][1] = int_stack + 699;
 Libderiv->deriv_classes[1][3][1] = int_stack + 717;
 Libderiv->deriv_classes[1][4][1] = int_stack + 747;
 Libderiv->deriv_classes[1][2][0] = int_stack + 792;
 Libderiv->deriv_classes[1][3][0] = int_stack + 810;
 Libderiv->deriv_classes[1][4][0] = int_stack + 840;
 memset(int_stack,0,7080);

 Libderiv->dvrr_stack = int_stack + 1551;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+885,int_stack+501,int_stack+465,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+939,int_stack+18,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+465,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+993,int_stack+48,int_stack+18, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+501,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+111,int_stack+93, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+465, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1083,int_stack+141,int_stack+111, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+501, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+54,int_stack+204,int_stack+186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+465, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+108,int_stack+234,int_stack+204, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+501, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+198,int_stack+297,int_stack+279, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1173,int_stack+327,int_stack+297, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+501, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+252,int_stack+390,int_stack+372, 0.0, zero_stack, 1.0, int_stack+465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1263,int_stack+420,int_stack+390, 0.0, zero_stack, 1.0, int_stack+501, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+306,int_stack+531,int_stack+483, 1.0, int_stack+465, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+561,int_stack+531, 1.0, int_stack+501, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+450,int_stack+624,int_stack+606,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+504,int_stack+654,int_stack+624,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+594,int_stack+717,int_stack+699,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1353,int_stack+747,int_stack+717,3);
 /*--- compute (p0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+648,int_stack+810,int_stack+792,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+702,int_stack+840,int_stack+810,3);
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1443,int_stack+993,int_stack+939, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+885,3);
     Libderiv->ABCD[11] = int_stack + 1443;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+939,int_stack+1083,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+885, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 939;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1047,int_stack+108,int_stack+54, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+885, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 1047;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+1173,int_stack+198, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1155,int_stack+1263,int_stack+252, 0.0, zero_stack, 1.0, int_stack+885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 1155;
 /*--- compute (p0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+108,int_stack+360,int_stack+306, 1.0, int_stack+885, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 108;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+216,int_stack+504,int_stack+450,3);
     Libderiv->ABCD[2] = int_stack + 216;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+324,int_stack+1353,int_stack+594,3);
     Libderiv->ABCD[1] = int_stack + 324;
 /*--- compute (p0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1263,int_stack+702,int_stack+648,3);
     Libderiv->ABCD[0] = int_stack + 1263;

}
