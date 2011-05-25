#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_dpg0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (dp|g0) integrals */

void d1hrr_order_dpg0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 90;
 Libderiv->deriv_classes[2][4][10] = int_stack + 240;
 Libderiv->deriv_classes[3][4][10] = int_stack + 330;
 Libderiv->deriv_classes[2][4][9] = int_stack + 480;
 Libderiv->deriv_classes[3][4][9] = int_stack + 570;
 Libderiv->deriv_classes[2][4][8] = int_stack + 720;
 Libderiv->deriv_classes[3][4][8] = int_stack + 810;
 Libderiv->deriv_classes[2][4][7] = int_stack + 960;
 Libderiv->deriv_classes[3][4][7] = int_stack + 1050;
 Libderiv->deriv_classes[2][4][6] = int_stack + 1200;
 Libderiv->deriv_classes[3][4][6] = int_stack + 1290;
 Libderiv->deriv_classes[2][4][2] = int_stack + 1440;
 Libderiv->deriv_classes[3][4][2] = int_stack + 1530;
 Libderiv->deriv_classes[2][4][1] = int_stack + 1680;
 Libderiv->deriv_classes[3][4][1] = int_stack + 1770;
 Libderiv->dvrr_classes[2][4] = int_stack + 1920;
 Libderiv->deriv_classes[2][4][0] = int_stack + 2010;
 Libderiv->deriv_classes[3][4][0] = int_stack + 2100;
 memset(int_stack,0,18000);

 Libderiv->dvrr_stack = int_stack + 2790;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_dpg0(Libderiv, Data);
   Data++;
 }

 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2250,int_stack+90,int_stack+0,15);
     Libderiv->ABCD[11] = int_stack + 2250;
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+2520,int_stack+330,int_stack+240,15);
     Libderiv->ABCD[10] = int_stack + 2520;
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+0,int_stack+570,int_stack+480,15);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+270,int_stack+810,int_stack+720,15);
     Libderiv->ABCD[8] = int_stack + 270;
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+540,int_stack+1050,int_stack+960,15);
     Libderiv->ABCD[7] = int_stack + 540;
 /*--- compute (dp|g0) ---*/
   hrr1_build_dp(Libderiv->AB,int_stack+810,int_stack+1290,int_stack+1200,15);
     Libderiv->ABCD[6] = int_stack + 810;
 /*--- compute (dp|g0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1080,int_stack+1530,int_stack+1440, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[2] = int_stack + 1080;
 /*--- compute (dp|g0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1350,int_stack+1770,int_stack+1680, 0.0, zero_stack, 1.0, int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[1] = int_stack + 1350;
 /*--- compute (dp|g0) ---*/
   d1hrr1_build_dp(Libderiv->AB,int_stack+1620,int_stack+2100,int_stack+2010, 1.0, int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[0] = int_stack + 1620;

}
