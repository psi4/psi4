#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|fd) integrals */

void d1hrr_order_f0fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 100;
 Libderiv->deriv_classes[3][5][11] = int_stack + 250;
 Libderiv->deriv_classes[3][3][10] = int_stack + 460;
 Libderiv->deriv_classes[3][4][10] = int_stack + 560;
 Libderiv->deriv_classes[3][5][10] = int_stack + 710;
 Libderiv->deriv_classes[3][3][9] = int_stack + 920;
 Libderiv->deriv_classes[3][4][9] = int_stack + 1020;
 Libderiv->deriv_classes[3][5][9] = int_stack + 1170;
 Libderiv->deriv_classes[3][3][8] = int_stack + 1380;
 Libderiv->deriv_classes[3][4][8] = int_stack + 1480;
 Libderiv->deriv_classes[3][5][8] = int_stack + 1630;
 Libderiv->deriv_classes[3][3][7] = int_stack + 1840;
 Libderiv->deriv_classes[3][4][7] = int_stack + 1940;
 Libderiv->deriv_classes[3][5][7] = int_stack + 2090;
 Libderiv->dvrr_classes[3][3] = int_stack + 2300;
 Libderiv->deriv_classes[3][3][6] = int_stack + 2400;
 Libderiv->dvrr_classes[3][4] = int_stack + 2500;
 Libderiv->deriv_classes[3][4][6] = int_stack + 2650;
 Libderiv->deriv_classes[3][5][6] = int_stack + 2800;
 Libderiv->deriv_classes[3][3][2] = int_stack + 3010;
 Libderiv->deriv_classes[3][4][2] = int_stack + 3110;
 Libderiv->deriv_classes[3][5][2] = int_stack + 3260;
 Libderiv->deriv_classes[3][3][1] = int_stack + 3470;
 Libderiv->deriv_classes[3][4][1] = int_stack + 3570;
 Libderiv->deriv_classes[3][5][1] = int_stack + 3720;
 Libderiv->deriv_classes[3][3][0] = int_stack + 3930;
 Libderiv->deriv_classes[3][4][0] = int_stack + 4030;
 Libderiv->deriv_classes[3][5][0] = int_stack + 4180;
 memset(int_stack,0,35120);

 Libderiv->dvrr_stack = int_stack + 7690;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+4390,int_stack+2500,int_stack+2300,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+4690,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2300,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4990,int_stack+250,int_stack+100, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2500,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+560,int_stack+460, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2300, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5440,int_stack+710,int_stack+560, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2500, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+1020,int_stack+920, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2300, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5890,int_stack+1170,int_stack+1020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2500, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+1480,int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+900,int_stack+1630,int_stack+1480, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1350,int_stack+1940,int_stack+1840, 0.0, zero_stack, 1.0, int_stack+2300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+6340,int_stack+2090,int_stack+1940, 0.0, zero_stack, 1.0, int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+1650,int_stack+2650,int_stack+2400, 1.0, int_stack+2300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1950,int_stack+2800,int_stack+2650, 1.0, int_stack+2500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2400,int_stack+3110,int_stack+3010,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6790,int_stack+3260,int_stack+3110,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2700,int_stack+3570,int_stack+3470,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+3000,int_stack+3720,int_stack+3570,10);
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+3450,int_stack+4030,int_stack+3930,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7240,int_stack+4180,int_stack+4030,10);
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3750,int_stack+4990,int_stack+4690, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4390,10);
     Libderiv->ABCD[11] = int_stack + 3750;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+4690,int_stack+5440,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4390, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 4690;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5290,int_stack+5890,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4390, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 5290;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+900,int_stack+600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+4390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+600,int_stack+6340,int_stack+1350, 0.0, zero_stack, 1.0, int_stack+4390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 600;
 /*--- compute (f0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+5890,int_stack+1950,int_stack+1650, 1.0, int_stack+4390, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 5890;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1200,int_stack+6790,int_stack+2400,10);
     Libderiv->ABCD[2] = int_stack + 1200;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1800,int_stack+3000,int_stack+2700,10);
     Libderiv->ABCD[1] = int_stack + 1800;
 /*--- compute (f0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+2400,int_stack+7240,int_stack+3450,10);
     Libderiv->ABCD[0] = int_stack + 2400;

}
