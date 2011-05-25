#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_g0gf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (g0|gf) integrals */

void d1hrr_order_g0gf(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 540;
 Libderiv->deriv_classes[4][7][11] = int_stack + 960;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1500;
 Libderiv->deriv_classes[4][5][10] = int_stack + 1725;
 Libderiv->deriv_classes[4][6][10] = int_stack + 2040;
 Libderiv->deriv_classes[4][7][10] = int_stack + 2460;
 Libderiv->deriv_classes[4][4][9] = int_stack + 3000;
 Libderiv->deriv_classes[4][5][9] = int_stack + 3225;
 Libderiv->deriv_classes[4][6][9] = int_stack + 3540;
 Libderiv->deriv_classes[4][7][9] = int_stack + 3960;
 Libderiv->deriv_classes[4][4][8] = int_stack + 4500;
 Libderiv->deriv_classes[4][5][8] = int_stack + 4725;
 Libderiv->deriv_classes[4][6][8] = int_stack + 5040;
 Libderiv->deriv_classes[4][7][8] = int_stack + 5460;
 Libderiv->deriv_classes[4][4][7] = int_stack + 6000;
 Libderiv->deriv_classes[4][5][7] = int_stack + 6225;
 Libderiv->deriv_classes[4][6][7] = int_stack + 6540;
 Libderiv->deriv_classes[4][7][7] = int_stack + 6960;
 Libderiv->dvrr_classes[4][4] = int_stack + 7500;
 Libderiv->deriv_classes[4][4][6] = int_stack + 7725;
 Libderiv->dvrr_classes[4][5] = int_stack + 7950;
 Libderiv->deriv_classes[4][5][6] = int_stack + 8265;
 Libderiv->dvrr_classes[4][6] = int_stack + 8580;
 Libderiv->deriv_classes[4][6][6] = int_stack + 9000;
 Libderiv->deriv_classes[4][7][6] = int_stack + 9420;
 Libderiv->deriv_classes[4][4][2] = int_stack + 9960;
 Libderiv->deriv_classes[4][5][2] = int_stack + 10185;
 Libderiv->deriv_classes[4][6][2] = int_stack + 10500;
 Libderiv->deriv_classes[4][7][2] = int_stack + 10920;
 Libderiv->deriv_classes[4][4][1] = int_stack + 11460;
 Libderiv->deriv_classes[4][5][1] = int_stack + 11685;
 Libderiv->deriv_classes[4][6][1] = int_stack + 12000;
 Libderiv->deriv_classes[4][7][1] = int_stack + 12420;
 Libderiv->deriv_classes[4][4][0] = int_stack + 12960;
 Libderiv->deriv_classes[4][5][0] = int_stack + 13185;
 Libderiv->deriv_classes[4][6][0] = int_stack + 13500;
 Libderiv->deriv_classes[4][7][0] = int_stack + 13920;
 memset(int_stack,0,115680);

 Libderiv->dvrr_stack = int_stack + 37410;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_g0gf(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+14460,int_stack+7950,int_stack+7500,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15135,int_stack+8580,int_stack+7950,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+16080,int_stack+15135,int_stack+14460,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17430,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7500,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18105,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7950,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+19050,int_stack+18105,int_stack+17430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14460,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+960,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8580,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+21660,int_stack+20400,int_stack+18105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15135,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20400,int_stack+1725,int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7500, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17430,int_stack+2040,int_stack+1725, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7950, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+17430,int_stack+20400, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14460, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+2460,int_stack+2040, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8580, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+23550,int_stack+20400,int_stack+17430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15135, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17430,int_stack+3225,int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7500, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18105,int_stack+3540,int_stack+3225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7950, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1350,int_stack+18105,int_stack+17430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14460, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+3960,int_stack+3540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8580, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+25440,int_stack+20400,int_stack+18105, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15135, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20400,int_stack+4725,int_stack+4500, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17430,int_stack+5040,int_stack+4725, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2700,int_stack+17430,int_stack+20400, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+14460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+5460,int_stack+5040, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+8580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+4050,int_stack+20400,int_stack+17430, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+15135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+17430,int_stack+6225,int_stack+6000, 0.0, zero_stack, 1.0, int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+18105,int_stack+6540,int_stack+6225, 0.0, zero_stack, 1.0, int_stack+7950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+27330,int_stack+18105,int_stack+17430, 0.0, zero_stack, 1.0, int_stack+14460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+6960,int_stack+6540, 0.0, zero_stack, 1.0, int_stack+8580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+28680,int_stack+20400,int_stack+18105, 0.0, zero_stack, 1.0, int_stack+15135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+20400,int_stack+8265,int_stack+7725, 1.0, int_stack+7500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+17430,int_stack+9000,int_stack+8265, 1.0, int_stack+7950, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+5940,int_stack+17430,int_stack+20400, 1.0, int_stack+14460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+9420,int_stack+9000, 1.0, int_stack+8580, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+7290,int_stack+20400,int_stack+17430, 1.0, int_stack+15135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17430,int_stack+10185,int_stack+9960,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+18105,int_stack+10500,int_stack+10185,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+14460,int_stack+18105,int_stack+17430,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+10920,int_stack+10500,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+9180,int_stack+20400,int_stack+18105,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+20400,int_stack+11685,int_stack+11460,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+17430,int_stack+12000,int_stack+11685,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+30570,int_stack+17430,int_stack+20400,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+12420,int_stack+12000,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+11070,int_stack+20400,int_stack+17430,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+17430,int_stack+13185,int_stack+12960,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+18105,int_stack+13500,int_stack+13185,15);
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+31920,int_stack+18105,int_stack+17430,15);
 /*--- compute (g0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+20400,int_stack+13920,int_stack+13500,15);
 /*--- compute (g0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+33270,int_stack+20400,int_stack+18105,15);
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+35160,int_stack+21660,int_stack+19050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16080,15);
     Libderiv->ABCD[11] = int_stack + 35160;
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+17430,int_stack+23550,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16080, 0.0, zero_stack,15);
     Libderiv->ABCD[10] = int_stack + 17430;
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+19680,int_stack+25440,int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16080, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[9] = int_stack + 19680;
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+4050,int_stack+2700, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+16080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+2250,int_stack+28680,int_stack+27330, 0.0, zero_stack, 1.0, int_stack+16080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[7] = int_stack + 2250;
 /*--- compute (g0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+21930,int_stack+7290,int_stack+5940, 1.0, int_stack+16080, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[6] = int_stack + 21930;
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+4500,int_stack+9180,int_stack+14460,15);
     Libderiv->ABCD[2] = int_stack + 4500;
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+6750,int_stack+11070,int_stack+30570,15);
     Libderiv->ABCD[1] = int_stack + 6750;
 /*--- compute (g0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+9000,int_stack+33270,int_stack+31920,15);
     Libderiv->ABCD[0] = int_stack + 9000;

}
