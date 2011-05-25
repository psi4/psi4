#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0dd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|dd) integrals */

void d1hrr_order_f0dd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][2][11] = int_stack + 0;
 Libderiv->deriv_classes[3][3][11] = int_stack + 60;
 Libderiv->deriv_classes[3][4][11] = int_stack + 160;
 Libderiv->deriv_classes[3][2][10] = int_stack + 310;
 Libderiv->deriv_classes[3][3][10] = int_stack + 370;
 Libderiv->deriv_classes[3][4][10] = int_stack + 470;
 Libderiv->deriv_classes[3][2][9] = int_stack + 620;
 Libderiv->deriv_classes[3][3][9] = int_stack + 680;
 Libderiv->deriv_classes[3][4][9] = int_stack + 780;
 Libderiv->deriv_classes[3][2][8] = int_stack + 930;
 Libderiv->deriv_classes[3][3][8] = int_stack + 990;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1090;
 Libderiv->deriv_classes[3][2][7] = int_stack + 1240;
 Libderiv->deriv_classes[3][3][7] = int_stack + 1300;
 Libderiv->deriv_classes[3][4][7] = int_stack + 1400;
 Libderiv->dvrr_classes[3][2] = int_stack + 1550;
 Libderiv->deriv_classes[3][2][6] = int_stack + 1610;
 Libderiv->dvrr_classes[3][3] = int_stack + 1670;
 Libderiv->deriv_classes[3][3][6] = int_stack + 1770;
 Libderiv->deriv_classes[3][4][6] = int_stack + 1870;
 Libderiv->deriv_classes[3][2][2] = int_stack + 2020;
 Libderiv->deriv_classes[3][3][2] = int_stack + 2080;
 Libderiv->deriv_classes[3][4][2] = int_stack + 2180;
 Libderiv->deriv_classes[3][2][1] = int_stack + 2330;
 Libderiv->deriv_classes[3][3][1] = int_stack + 2390;
 Libderiv->deriv_classes[3][4][1] = int_stack + 2490;
 Libderiv->deriv_classes[3][2][0] = int_stack + 2640;
 Libderiv->deriv_classes[3][3][0] = int_stack + 2700;
 Libderiv->deriv_classes[3][4][0] = int_stack + 2800;
 memset(int_stack,0,23600);

 Libderiv->dvrr_stack = int_stack + 5170;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0dd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2950,int_stack+1670,int_stack+1550,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+3130,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1550,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3310,int_stack+160,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1670,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+0,int_stack+370,int_stack+310, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1550, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3610,int_stack+470,int_stack+370, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1670, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+180,int_stack+680,int_stack+620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1550, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+780,int_stack+680, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1670, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+660,int_stack+990,int_stack+930, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+3910,int_stack+1090,int_stack+990, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+840,int_stack+1300,int_stack+1240, 0.0, zero_stack, 1.0, int_stack+1550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4210,int_stack+1400,int_stack+1300, 0.0, zero_stack, 1.0, int_stack+1670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   d1hrr3_build_dp(Libderiv->CD,int_stack+1020,int_stack+1770,int_stack+1610, 1.0, int_stack+1550, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+1870,int_stack+1770, 1.0, int_stack+1670, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1500,int_stack+2080,int_stack+2020,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1680,int_stack+2180,int_stack+2080,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+1980,int_stack+2390,int_stack+2330,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4510,int_stack+2490,int_stack+2390,10);
 /*--- compute (f0|dp) ---*/
   hrr3_build_dp(Libderiv->CD,int_stack+2160,int_stack+2700,int_stack+2640,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2340,int_stack+2800,int_stack+2700,10);
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+4810,int_stack+3310,int_stack+3130, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2950,10);
     Libderiv->ABCD[11] = int_stack + 4810;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3130,int_stack+3610,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2950, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 3130;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3490,int_stack+360,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2950, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 3490;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+0,int_stack+3910,int_stack+660, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+3850,int_stack+4210,int_stack+840, 0.0, zero_stack, 1.0, int_stack+2950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 3850;
 /*--- compute (f0|dd) ---*/
   d1hrr3_build_dd(Libderiv->CD,int_stack+360,int_stack+1200,int_stack+1020, 1.0, int_stack+2950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 360;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+720,int_stack+1680,int_stack+1500,10);
     Libderiv->ABCD[2] = int_stack + 720;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+1080,int_stack+4510,int_stack+1980,10);
     Libderiv->ABCD[1] = int_stack + 1080;
 /*--- compute (f0|dd) ---*/
   hrr3_build_dd(Libderiv->CD,int_stack+4210,int_stack+2340,int_stack+2160,10);
     Libderiv->ABCD[0] = int_stack + 4210;

}
