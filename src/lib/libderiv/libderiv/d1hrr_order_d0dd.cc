#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|dd) integrals */

void d1hrr_order_d0dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][2][11] = int_stack + 0;
 Libderiv->deriv_classes[2][3][11] = int_stack + 36;
 Libderiv->deriv_classes[2][4][11] = int_stack + 96;
 Libderiv->deriv_classes[2][2][10] = int_stack + 186;
 Libderiv->deriv_classes[2][3][10] = int_stack + 222;
 Libderiv->deriv_classes[2][4][10] = int_stack + 282;
 Libderiv->deriv_classes[2][2][9] = int_stack + 372;
 Libderiv->deriv_classes[2][3][9] = int_stack + 408;
 Libderiv->deriv_classes[2][4][9] = int_stack + 468;
 Libderiv->deriv_classes[2][2][8] = int_stack + 558;
 Libderiv->deriv_classes[2][3][8] = int_stack + 594;
 Libderiv->deriv_classes[2][4][8] = int_stack + 654;
 Libderiv->deriv_classes[2][2][7] = int_stack + 744;
 Libderiv->deriv_classes[2][3][7] = int_stack + 780;
 Libderiv->deriv_classes[2][4][7] = int_stack + 840;
 Libderiv->dvrr_classes[2][2] = int_stack + 930;
 Libderiv->deriv_classes[2][2][6] = int_stack + 966;
 Libderiv->dvrr_classes[2][3] = int_stack + 1002;
 Libderiv->deriv_classes[2][3][6] = int_stack + 1062;
 Libderiv->deriv_classes[2][4][6] = int_stack + 1122;
 Libderiv->deriv_classes[2][2][2] = int_stack + 1212;
 Libderiv->deriv_classes[2][3][2] = int_stack + 1248;
 Libderiv->deriv_classes[2][4][2] = int_stack + 1308;
 Libderiv->deriv_classes[2][2][1] = int_stack + 1398;
 Libderiv->deriv_classes[2][3][1] = int_stack + 1434;
 Libderiv->deriv_classes[2][4][1] = int_stack + 1494;
 Libderiv->deriv_classes[2][2][0] = int_stack + 1584;
 Libderiv->deriv_classes[2][3][0] = int_stack + 1620;
 Libderiv->deriv_classes[2][4][0] = int_stack + 1680;
 memset(int_stack,0,14160);

 Libderiv->dvrr_stack = int_stack + 3102;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1770,int_stack+1002,int_stack+930,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1878,int_stack+36,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+930,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1986,int_stack+96,int_stack+36, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1002,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+222,int_stack+186, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+930, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2166,int_stack+282,int_stack+222, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1002, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+108,int_stack+408,int_stack+372, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+930, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+216,int_stack+468,int_stack+408, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1002, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+396,int_stack+594,int_stack+558, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2346,int_stack+654,int_stack+594, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+504,int_stack+780,int_stack+744, 0.0, zero_stack, 1.0, int_stack+930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2526,int_stack+840,int_stack+780, 0.0, zero_stack, 1.0, int_stack+1002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+612,int_stack+1062,int_stack+966, 1.0, int_stack+930, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+1122,int_stack+1062, 1.0, int_stack+1002, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+900,int_stack+1248,int_stack+1212,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1008,int_stack+1308,int_stack+1248,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1188,int_stack+1434,int_stack+1398,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2706,int_stack+1494,int_stack+1434,6);
 /*--- compute (d0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1296,int_stack+1620,int_stack+1584,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1404,int_stack+1680,int_stack+1620,6);
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2886,int_stack+1986,int_stack+1878, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1770,6);
     Libderiv->ABCD[11] = int_stack + 2886;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+1878,int_stack+2166,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1770, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 1878;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2094,int_stack+216,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1770, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 2094;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+2346,int_stack+396, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+2310,int_stack+2526,int_stack+504, 0.0, zero_stack, 1.0, int_stack+1770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 2310;
 /*--- compute (d0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+216,int_stack+720,int_stack+612, 1.0, int_stack+1770, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 216;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+432,int_stack+1008,int_stack+900,6);
     Libderiv->ABCD[2] = int_stack + 432;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+648,int_stack+2706,int_stack+1188,6);
     Libderiv->ABCD[1] = int_stack + 648;
 /*--- compute (d0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+2526,int_stack+1404,int_stack+1296,6);
     Libderiv->ABCD[0] = int_stack + 2526;

}
