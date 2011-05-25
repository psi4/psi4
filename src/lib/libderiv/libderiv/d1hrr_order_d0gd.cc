#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0gd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|gd) integrals */

void d1hrr_order_d0gd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][5][11] = int_stack + 90;
 Libderiv->deriv_classes[2][6][11] = int_stack + 216;
 Libderiv->deriv_classes[2][4][10] = int_stack + 384;
 Libderiv->deriv_classes[2][5][10] = int_stack + 474;
 Libderiv->deriv_classes[2][6][10] = int_stack + 600;
 Libderiv->deriv_classes[2][4][9] = int_stack + 768;
 Libderiv->deriv_classes[2][5][9] = int_stack + 858;
 Libderiv->deriv_classes[2][6][9] = int_stack + 984;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1152;
 Libderiv->deriv_classes[2][5][8] = int_stack + 1242;
 Libderiv->deriv_classes[2][6][8] = int_stack + 1368;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1536;
 Libderiv->deriv_classes[2][5][7] = int_stack + 1626;
 Libderiv->deriv_classes[2][6][7] = int_stack + 1752;
 Libderiv->dvrr_classes[2][4] = int_stack + 1920;
 Libderiv->deriv_classes[2][4][6] = int_stack + 2010;
 Libderiv->dvrr_classes[2][5] = int_stack + 2100;
 Libderiv->deriv_classes[2][5][6] = int_stack + 2226;
 Libderiv->deriv_classes[2][6][6] = int_stack + 2352;
 Libderiv->deriv_classes[2][4][2] = int_stack + 2520;
 Libderiv->deriv_classes[2][5][2] = int_stack + 2610;
 Libderiv->deriv_classes[2][6][2] = int_stack + 2736;
 Libderiv->deriv_classes[2][4][1] = int_stack + 2904;
 Libderiv->deriv_classes[2][5][1] = int_stack + 2994;
 Libderiv->deriv_classes[2][6][1] = int_stack + 3120;
 Libderiv->deriv_classes[2][4][0] = int_stack + 3288;
 Libderiv->deriv_classes[2][5][0] = int_stack + 3378;
 Libderiv->deriv_classes[2][6][0] = int_stack + 3504;
 memset(int_stack,0,29376);

 Libderiv->dvrr_stack = int_stack + 7560;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0gd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3672,int_stack+2100,int_stack+1920,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3942,int_stack+90,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1920,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4212,int_stack+216,int_stack+90, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+474,int_stack+384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1920, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4590,int_stack+600,int_stack+474, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+858,int_stack+768, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1920, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+4968,int_stack+984,int_stack+858, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+1242,int_stack+1152, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+810,int_stack+1368,int_stack+1242, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1188,int_stack+1626,int_stack+1536, 0.0, zero_stack, 1.0, int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5346,int_stack+1752,int_stack+1626, 0.0, zero_stack, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1458,int_stack+2226,int_stack+2010, 1.0, int_stack+1920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+5724,int_stack+2352,int_stack+2226, 1.0, int_stack+2100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1728,int_stack+2610,int_stack+2520,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1998,int_stack+2736,int_stack+2610,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2376,int_stack+2994,int_stack+2904,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+6102,int_stack+3120,int_stack+2994,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2646,int_stack+3378,int_stack+3288,6);
 /*--- compute (d0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+2916,int_stack+3504,int_stack+3378,6);
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+6480,int_stack+4212,int_stack+3942, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3672,6);
     Libderiv->ABCD[11] = int_stack + 6480;
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3942,int_stack+4590,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3672, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 3942;
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+7020,int_stack+4968,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3672, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 7020;
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+810,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+540,int_stack+5346,int_stack+1188, 0.0, zero_stack, 1.0, int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 540;
 /*--- compute (d0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+4482,int_stack+5724,int_stack+1458, 1.0, int_stack+3672, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 4482;
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1080,int_stack+1998,int_stack+1728,6);
     Libderiv->ABCD[2] = int_stack + 1080;
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+1620,int_stack+6102,int_stack+2376,6);
     Libderiv->ABCD[1] = int_stack + 1620;
 /*--- compute (d0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+5022,int_stack+2916,int_stack+2646,6);
     Libderiv->ABCD[0] = int_stack + 5022;

}
