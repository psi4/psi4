#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0gd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|gd) integrals */

void d1hrr_order_f0gd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][6][11] = int_stack + 360;
 Libderiv->deriv_classes[3][4][10] = int_stack + 640;
 Libderiv->deriv_classes[3][5][10] = int_stack + 790;
 Libderiv->deriv_classes[3][6][10] = int_stack + 1000;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1280;
 Libderiv->deriv_classes[3][5][9] = int_stack + 1430;
 Libderiv->deriv_classes[3][6][9] = int_stack + 1640;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1920;
 Libderiv->deriv_classes[3][5][8] = int_stack + 2070;
 Libderiv->deriv_classes[3][6][8] = int_stack + 2280;
 Libderiv->deriv_classes[3][4][7] = int_stack + 2560;
 Libderiv->deriv_classes[3][5][7] = int_stack + 2710;
 Libderiv->deriv_classes[3][6][7] = int_stack + 2920;
 Libderiv->dvrr_classes[3][4] = int_stack + 3200;
 Libderiv->deriv_classes[3][4][6] = int_stack + 3350;
 Libderiv->dvrr_classes[3][5] = int_stack + 3500;
 Libderiv->deriv_classes[3][5][6] = int_stack + 3710;
 Libderiv->deriv_classes[3][6][6] = int_stack + 3920;
 Libderiv->deriv_classes[3][4][2] = int_stack + 4200;
 Libderiv->deriv_classes[3][5][2] = int_stack + 4350;
 Libderiv->deriv_classes[3][6][2] = int_stack + 4560;
 Libderiv->deriv_classes[3][4][1] = int_stack + 4840;
 Libderiv->deriv_classes[3][5][1] = int_stack + 4990;
 Libderiv->deriv_classes[3][6][1] = int_stack + 5200;
 Libderiv->deriv_classes[3][4][0] = int_stack + 5480;
 Libderiv->deriv_classes[3][5][0] = int_stack + 5630;
 Libderiv->deriv_classes[3][6][0] = int_stack + 5840;
 memset(int_stack,0,48960);

 Libderiv->dvrr_stack = int_stack + 12600;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0gd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6120,int_stack+3500,int_stack+3200,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6570,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3200,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7020,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3500,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+790,int_stack+640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3200, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+7650,int_stack+1000,int_stack+790, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3500, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+450,int_stack+1430,int_stack+1280, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3200, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8280,int_stack+1640,int_stack+1430, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+2070,int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+1350,int_stack+2280,int_stack+2070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1980,int_stack+2710,int_stack+2560, 0.0, zero_stack, 1.0, int_stack+3200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+8910,int_stack+2920,int_stack+2710, 0.0, zero_stack, 1.0, int_stack+3500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2430,int_stack+3710,int_stack+3350, 1.0, int_stack+3200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+9540,int_stack+3920,int_stack+3710, 1.0, int_stack+3500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2880,int_stack+4350,int_stack+4200,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3330,int_stack+4560,int_stack+4350,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3960,int_stack+4990,int_stack+4840,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10170,int_stack+5200,int_stack+4990,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4410,int_stack+5630,int_stack+5480,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4860,int_stack+5840,int_stack+5630,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+10800,int_stack+7020,int_stack+6570, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6120,10);
     Libderiv->ABCD[11] = int_stack + 10800;
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6570,int_stack+7650,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6120, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 6570;
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11700,int_stack+8280,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6120, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 11700;
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+1350,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+900,int_stack+8910,int_stack+1980, 0.0, zero_stack, 1.0, int_stack+6120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 900;
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7470,int_stack+9540,int_stack+2430, 1.0, int_stack+6120, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 7470;
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1800,int_stack+3330,int_stack+2880,10);
     Libderiv->ABCD[2] = int_stack + 1800;
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2700,int_stack+10170,int_stack+3960,10);
     Libderiv->ABCD[1] = int_stack + 2700;
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+8370,int_stack+4860,int_stack+4410,10);
     Libderiv->ABCD[0] = int_stack + 8370;

}
