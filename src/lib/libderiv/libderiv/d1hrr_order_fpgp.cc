#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_fpgp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (fp|gp) integrals */

void d1hrr_order_fpgp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[4][4][11] = int_stack + 360;
 Libderiv->deriv_classes[4][5][11] = int_stack + 585;
 Libderiv->deriv_classes[3][4][10] = int_stack + 900;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1050;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1260;
 Libderiv->deriv_classes[4][5][10] = int_stack + 1485;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1800;
 Libderiv->deriv_classes[3][5][9] = int_stack + 1950;
 Libderiv->deriv_classes[4][4][9] = int_stack + 2160;
 Libderiv->deriv_classes[4][5][9] = int_stack + 2385;
 Libderiv->deriv_classes[3][4][8] = int_stack + 2700;
 Libderiv->deriv_classes[3][5][8] = int_stack + 2850;
 Libderiv->deriv_classes[4][4][8] = int_stack + 3060;
 Libderiv->deriv_classes[4][5][8] = int_stack + 3285;
 Libderiv->deriv_classes[3][4][7] = int_stack + 3600;
 Libderiv->deriv_classes[3][5][7] = int_stack + 3750;
 Libderiv->deriv_classes[4][4][7] = int_stack + 3960;
 Libderiv->deriv_classes[4][5][7] = int_stack + 4185;
 Libderiv->deriv_classes[3][4][6] = int_stack + 4500;
 Libderiv->deriv_classes[3][5][6] = int_stack + 4650;
 Libderiv->dvrr_classes[4][4] = int_stack + 4860;
 Libderiv->deriv_classes[4][4][6] = int_stack + 5085;
 Libderiv->deriv_classes[4][5][6] = int_stack + 5310;
 Libderiv->deriv_classes[3][4][2] = int_stack + 5625;
 Libderiv->deriv_classes[3][5][2] = int_stack + 5775;
 Libderiv->deriv_classes[4][4][2] = int_stack + 5985;
 Libderiv->deriv_classes[4][5][2] = int_stack + 6210;
 Libderiv->deriv_classes[3][4][1] = int_stack + 6525;
 Libderiv->deriv_classes[3][5][1] = int_stack + 6675;
 Libderiv->deriv_classes[4][4][1] = int_stack + 6885;
 Libderiv->deriv_classes[4][5][1] = int_stack + 7110;
 Libderiv->dvrr_classes[3][4] = int_stack + 7425;
 Libderiv->dvrr_classes[3][5] = int_stack + 7575;
 Libderiv->deriv_classes[3][4][0] = int_stack + 7785;
 Libderiv->deriv_classes[3][5][0] = int_stack + 7935;
 Libderiv->deriv_classes[4][4][0] = int_stack + 8145;
 Libderiv->deriv_classes[4][5][0] = int_stack + 8370;
 memset(int_stack,0,69480);

 Libderiv->dvrr_stack = int_stack + 15885;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_fpgp(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8685,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7425,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9135,int_stack+585,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4860,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1050,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7425, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+1485,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4860, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1125,int_stack+1950,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7425, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9810,int_stack+2385,int_stack+2160, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4860, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1575,int_stack+2850,int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2025,int_stack+3285,int_stack+3060, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2700,int_stack+3750,int_stack+3600, 0.0, zero_stack, 1.0, int_stack+7425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3150,int_stack+4185,int_stack+3960, 0.0, zero_stack, 1.0, int_stack+4860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3825,int_stack+4650,int_stack+4500, 1.0, int_stack+7425, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+10485,int_stack+5310,int_stack+5085, 1.0, int_stack+4860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4275,int_stack+7575,int_stack+7425,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4725,int_stack+5775,int_stack+5625,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5175,int_stack+6210,int_stack+5985,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5850,int_stack+6675,int_stack+6525,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11160,int_stack+7110,int_stack+6885,15);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6300,int_stack+7935,int_stack+7785,10);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6750,int_stack+8370,int_stack+8145,15);
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+11835,int_stack+9135,int_stack+8685,45);
     Libderiv->ABCD[11] = int_stack + 11835;
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+7425,int_stack+450,int_stack+0,45);
     Libderiv->ABCD[10] = int_stack + 7425;
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+13185,int_stack+9810,int_stack+1125,45);
     Libderiv->ABCD[9] = int_stack + 13185;
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+0,int_stack+2025,int_stack+1575,45);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+1350,int_stack+3150,int_stack+2700,45);
     Libderiv->ABCD[7] = int_stack + 1350;
 /*--- compute (fp|gp) ---*/
   hrr1_build_fp(Libderiv->AB,int_stack+8775,int_stack+10485,int_stack+3825,45);
     Libderiv->ABCD[6] = int_stack + 8775;
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+2700,int_stack+5175,int_stack+4725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[2] = int_stack + 2700;
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+14535,int_stack+11160,int_stack+5850, 0.0, zero_stack, 1.0, int_stack+4275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[1] = int_stack + 14535;
 /*--- compute (fp|gp) ---*/
   d1hrr1_build_fp(Libderiv->AB,int_stack+4725,int_stack+6750,int_stack+6300, 1.0, int_stack+4275, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[0] = int_stack + 4725;

}
