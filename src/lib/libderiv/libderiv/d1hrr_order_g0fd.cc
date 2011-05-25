#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_g0fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (g0|fd) integrals */

void d1hrr_order_g0fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][3][11] = int_stack + 0;
 Libderiv->deriv_classes[4][4][11] = int_stack + 150;
 Libderiv->deriv_classes[4][5][11] = int_stack + 375;
 Libderiv->deriv_classes[4][3][10] = int_stack + 690;
 Libderiv->deriv_classes[4][4][10] = int_stack + 840;
 Libderiv->deriv_classes[4][5][10] = int_stack + 1065;
 Libderiv->deriv_classes[4][3][9] = int_stack + 1380;
 Libderiv->deriv_classes[4][4][9] = int_stack + 1530;
 Libderiv->deriv_classes[4][5][9] = int_stack + 1755;
 Libderiv->deriv_classes[4][3][8] = int_stack + 2070;
 Libderiv->deriv_classes[4][4][8] = int_stack + 2220;
 Libderiv->deriv_classes[4][5][8] = int_stack + 2445;
 Libderiv->deriv_classes[4][3][7] = int_stack + 2760;
 Libderiv->deriv_classes[4][4][7] = int_stack + 2910;
 Libderiv->deriv_classes[4][5][7] = int_stack + 3135;
 Libderiv->dvrr_classes[4][3] = int_stack + 3450;
 Libderiv->deriv_classes[4][3][6] = int_stack + 3600;
 Libderiv->dvrr_classes[4][4] = int_stack + 3750;
 Libderiv->deriv_classes[4][4][6] = int_stack + 3975;
 Libderiv->deriv_classes[4][5][6] = int_stack + 4200;
 Libderiv->deriv_classes[4][3][2] = int_stack + 4515;
 Libderiv->deriv_classes[4][4][2] = int_stack + 4665;
 Libderiv->deriv_classes[4][5][2] = int_stack + 4890;
 Libderiv->deriv_classes[4][3][1] = int_stack + 5205;
 Libderiv->deriv_classes[4][4][1] = int_stack + 5355;
 Libderiv->deriv_classes[4][5][1] = int_stack + 5580;
 Libderiv->deriv_classes[4][3][0] = int_stack + 5895;
 Libderiv->deriv_classes[4][4][0] = int_stack + 6045;
 Libderiv->deriv_classes[4][5][0] = int_stack + 6270;
 memset(int_stack,0,52680);

 Libderiv->dvrr_stack = int_stack + 11535;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_g0fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+6585,int_stack+3750,int_stack+3450,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+7035,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3450,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+7485,int_stack+375,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+840,int_stack+690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3450, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8160,int_stack+1065,int_stack+840, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+1530,int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3450, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+8835,int_stack+1755,int_stack+1530, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+2220,int_stack+2070, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1350,int_stack+2445,int_stack+2220, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2025,int_stack+2910,int_stack+2760, 0.0, zero_stack, 1.0, int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+9510,int_stack+3135,int_stack+2910, 0.0, zero_stack, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2475,int_stack+3975,int_stack+3600, 1.0, int_stack+3450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2925,int_stack+4200,int_stack+3975, 1.0, int_stack+3750, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3600,int_stack+4665,int_stack+4515,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10185,int_stack+4890,int_stack+4665,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4050,int_stack+5355,int_stack+5205,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4500,int_stack+5580,int_stack+5355,15);
 /*--- compute (g0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+5175,int_stack+6045,int_stack+5895,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+10860,int_stack+6270,int_stack+6045,15);
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5625,int_stack+7485,int_stack+7035, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6585,15);
     Libderiv->ABCD[11] = int_stack + 5625;
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7035,int_stack+8160,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6585, 0.0, zero_stack,15);
     Libderiv->ABCD[10] = int_stack + 7035;
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+7935,int_stack+8835,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6585, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[9] = int_stack + 7935;
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+1350,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+900,int_stack+9510,int_stack+2025, 0.0, zero_stack, 1.0, int_stack+6585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[7] = int_stack + 900;
 /*--- compute (g0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+8835,int_stack+2925,int_stack+2475, 1.0, int_stack+6585, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
     Libderiv->ABCD[6] = int_stack + 8835;
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1800,int_stack+10185,int_stack+3600,15);
     Libderiv->ABCD[2] = int_stack + 1800;
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2700,int_stack+4500,int_stack+4050,15);
     Libderiv->ABCD[1] = int_stack + 2700;
 /*--- compute (g0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+3600,int_stack+10860,int_stack+5175,15);
     Libderiv->ABCD[0] = int_stack + 3600;

}
