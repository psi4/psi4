#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|dd) integrals */

void d1hrr_order_00dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][2][11] = int_stack + 0;
 Libderiv->deriv_classes[0][3][11] = int_stack + 6;
 Libderiv->deriv_classes[0][4][11] = int_stack + 16;
 Libderiv->deriv_classes[0][2][10] = int_stack + 31;
 Libderiv->deriv_classes[0][3][10] = int_stack + 37;
 Libderiv->deriv_classes[0][4][10] = int_stack + 47;
 Libderiv->deriv_classes[0][2][9] = int_stack + 62;
 Libderiv->deriv_classes[0][3][9] = int_stack + 68;
 Libderiv->deriv_classes[0][4][9] = int_stack + 78;
 Libderiv->deriv_classes[0][2][8] = int_stack + 93;
 Libderiv->deriv_classes[0][3][8] = int_stack + 99;
 Libderiv->deriv_classes[0][4][8] = int_stack + 109;
 Libderiv->deriv_classes[0][2][7] = int_stack + 124;
 Libderiv->deriv_classes[0][3][7] = int_stack + 130;
 Libderiv->deriv_classes[0][4][7] = int_stack + 140;
 Libderiv->dvrr_classes[0][2] = int_stack + 155;
 Libderiv->deriv_classes[0][2][6] = int_stack + 161;
 Libderiv->dvrr_classes[0][3] = int_stack + 167;
 Libderiv->deriv_classes[0][3][6] = int_stack + 177;
 Libderiv->deriv_classes[0][4][6] = int_stack + 187;
 Libderiv->deriv_classes[0][2][2] = int_stack + 202;
 Libderiv->deriv_classes[0][3][2] = int_stack + 208;
 Libderiv->deriv_classes[0][4][2] = int_stack + 218;
 Libderiv->deriv_classes[0][2][1] = int_stack + 233;
 Libderiv->deriv_classes[0][3][1] = int_stack + 239;
 Libderiv->deriv_classes[0][4][1] = int_stack + 249;
 Libderiv->deriv_classes[0][2][0] = int_stack + 264;
 Libderiv->deriv_classes[0][3][0] = int_stack + 270;
 Libderiv->deriv_classes[0][4][0] = int_stack + 280;
 memset(int_stack,0,2360);

 Libderiv->dvrr_stack = int_stack + 517;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+295,int_stack+167,int_stack+155,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+313,int_stack+6,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+331,int_stack+16,int_stack+6, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+37,int_stack+31, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+361,int_stack+47,int_stack+37, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+18,int_stack+68,int_stack+62, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+36,int_stack+78,int_stack+68, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+66,int_stack+99,int_stack+93, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+391,int_stack+109,int_stack+99, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+167, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+84,int_stack+130,int_stack+124, 0.0, zero_stack, 1.0, int_stack+155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+421,int_stack+140,int_stack+130, 0.0, zero_stack, 1.0, int_stack+167, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+102,int_stack+177,int_stack+161, 1.0, int_stack+155, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+120,int_stack+187,int_stack+177, 1.0, int_stack+167, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+150,int_stack+208,int_stack+202,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+168,int_stack+218,int_stack+208,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+198,int_stack+239,int_stack+233,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+451,int_stack+249,int_stack+239,1);
 /*--- compute (00|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+216,int_stack+270,int_stack+264,1);
 /*--- compute (00|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+234,int_stack+280,int_stack+270,1);
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+481,int_stack+331,int_stack+313, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295,1);
     Libderiv->ABCD[11] = int_stack + 481;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+313,int_stack+361,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295, 0.0, zero_stack,1);
     Libderiv->ABCD[10] = int_stack + 313;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+349,int_stack+36,int_stack+18, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[9] = int_stack + 349;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+391,int_stack+66, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+385,int_stack+421,int_stack+84, 0.0, zero_stack, 1.0, int_stack+295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[7] = int_stack + 385;
 /*--- compute (00|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+36,int_stack+120,int_stack+102, 1.0, int_stack+295, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,1);
     Libderiv->ABCD[6] = int_stack + 36;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+72,int_stack+168,int_stack+150,1);
     Libderiv->ABCD[2] = int_stack + 72;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+108,int_stack+451,int_stack+198,1);
     Libderiv->ABCD[1] = int_stack + 108;
 /*--- compute (00|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+421,int_stack+234,int_stack+216,1);
     Libderiv->ABCD[0] = int_stack + 421;

}
