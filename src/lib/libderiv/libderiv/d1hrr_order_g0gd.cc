#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_g0gd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (g0|gd) integrals */

void d1hrr_order_g0gd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[4][6][11] = int_stack + 540;
 Libderiv->deriv_classes[4][4][10] = int_stack + 960;
 Libderiv->deriv_classes[4][5][10] = int_stack + 1185;
 Libderiv->deriv_classes[4][6][10] = int_stack + 1500;
 Libderiv->deriv_classes[4][4][9] = int_stack + 1920;
 Libderiv->deriv_classes[4][5][9] = int_stack + 2145;
 Libderiv->deriv_classes[4][6][9] = int_stack + 2460;
 Libderiv->deriv_classes[4][4][8] = int_stack + 2880;
 Libderiv->deriv_classes[4][5][8] = int_stack + 3105;
 Libderiv->deriv_classes[4][6][8] = int_stack + 3420;
 Libderiv->deriv_classes[4][4][7] = int_stack + 3840;
 Libderiv->deriv_classes[4][5][7] = int_stack + 4065;
 Libderiv->deriv_classes[4][6][7] = int_stack + 4380;
 Libderiv->dvrr_classes[4][4] = int_stack + 4800;
 Libderiv->deriv_classes[4][4][6] = int_stack + 5025;
 Libderiv->dvrr_classes[4][5] = int_stack + 5250;
 Libderiv->deriv_classes[4][5][6] = int_stack + 5565;
 Libderiv->deriv_classes[4][6][6] = int_stack + 5880;
 Libderiv->deriv_classes[4][4][2] = int_stack + 6300;
 Libderiv->deriv_classes[4][5][2] = int_stack + 6525;
 Libderiv->deriv_classes[4][6][2] = int_stack + 6840;
 Libderiv->deriv_classes[4][4][1] = int_stack + 7260;
 Libderiv->deriv_classes[4][5][1] = int_stack + 7485;
 Libderiv->deriv_classes[4][6][1] = int_stack + 7800;
 Libderiv->deriv_classes[4][4][0] = int_stack + 8220;
 Libderiv->deriv_classes[4][5][0] = int_stack + 8445;
 Libderiv->deriv_classes[4][6][0] = int_stack + 8760;
 memset(int_stack,0,73440);

 Libderiv->dvrr_stack = int_stack + 18900;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_g0gd(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9180,int_stack+5250,int_stack+4800,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9855,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+10530,int_stack+540,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5250,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1185,int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11475,int_stack+1500,int_stack+1185, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5250, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+675,int_stack+2145,int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12420,int_stack+2460,int_stack+2145, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5250, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1350,int_stack+3105,int_stack+2880, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2025,int_stack+3420,int_stack+3105, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2970,int_stack+4065,int_stack+3840, 0.0, zero_stack, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+13365,int_stack+4380,int_stack+4065, 0.0, zero_stack, 1.0, int_stack+5250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3645,int_stack+5565,int_stack+5025, 1.0, int_stack+4800, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+14310,int_stack+5880,int_stack+5565, 1.0, int_stack+5250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4320,int_stack+6525,int_stack+6300,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+4995,int_stack+6840,int_stack+6525,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+5940,int_stack+7485,int_stack+7260,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+15255,int_stack+7800,int_stack+7485,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6615,int_stack+8445,int_stack+8220,15);
 /*--- compute (g0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+7290,int_stack+8760,int_stack+8445,15);
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+16200,int_stack+10530,int_stack+9855, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9180,15);
     Libderiv->ABCD[11] = int_stack + 16200;
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+9855,int_stack+11475,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9180, 0.0, zero_stack,15);
     Libderiv->ABCD[10] = int_stack + 9855;
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+17550,int_stack+12420,int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9180, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[9] = int_stack + 17550;
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+2025,int_stack+1350, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1350,int_stack+13365,int_stack+2970, 0.0, zero_stack, 1.0, int_stack+9180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[7] = int_stack + 1350;
 /*--- compute (g0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+11205,int_stack+14310,int_stack+3645, 1.0, int_stack+9180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[6] = int_stack + 11205;
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2700,int_stack+4995,int_stack+4320,15);
     Libderiv->ABCD[2] = int_stack + 2700;
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+4050,int_stack+15255,int_stack+5940,15);
     Libderiv->ABCD[1] = int_stack + 4050;
 /*--- compute (g0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+12555,int_stack+7290,int_stack+6615,15);
     Libderiv->ABCD[0] = int_stack + 12555;

}
