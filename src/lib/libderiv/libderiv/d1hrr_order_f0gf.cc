#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0gf(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|gf) integrals */

void d1hrr_order_f0gf(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][4][11] = int_stack + 0;
 Libderiv->deriv_classes[3][5][11] = int_stack + 150;
 Libderiv->deriv_classes[3][6][11] = int_stack + 360;
 Libderiv->deriv_classes[3][7][11] = int_stack + 640;
 Libderiv->deriv_classes[3][4][10] = int_stack + 1000;
 Libderiv->deriv_classes[3][5][10] = int_stack + 1150;
 Libderiv->deriv_classes[3][6][10] = int_stack + 1360;
 Libderiv->deriv_classes[3][7][10] = int_stack + 1640;
 Libderiv->deriv_classes[3][4][9] = int_stack + 2000;
 Libderiv->deriv_classes[3][5][9] = int_stack + 2150;
 Libderiv->deriv_classes[3][6][9] = int_stack + 2360;
 Libderiv->deriv_classes[3][7][9] = int_stack + 2640;
 Libderiv->deriv_classes[3][4][8] = int_stack + 3000;
 Libderiv->deriv_classes[3][5][8] = int_stack + 3150;
 Libderiv->deriv_classes[3][6][8] = int_stack + 3360;
 Libderiv->deriv_classes[3][7][8] = int_stack + 3640;
 Libderiv->deriv_classes[3][4][7] = int_stack + 4000;
 Libderiv->deriv_classes[3][5][7] = int_stack + 4150;
 Libderiv->deriv_classes[3][6][7] = int_stack + 4360;
 Libderiv->deriv_classes[3][7][7] = int_stack + 4640;
 Libderiv->dvrr_classes[3][4] = int_stack + 5000;
 Libderiv->deriv_classes[3][4][6] = int_stack + 5150;
 Libderiv->dvrr_classes[3][5] = int_stack + 5300;
 Libderiv->deriv_classes[3][5][6] = int_stack + 5510;
 Libderiv->dvrr_classes[3][6] = int_stack + 5720;
 Libderiv->deriv_classes[3][6][6] = int_stack + 6000;
 Libderiv->deriv_classes[3][7][6] = int_stack + 6280;
 Libderiv->deriv_classes[3][4][2] = int_stack + 6640;
 Libderiv->deriv_classes[3][5][2] = int_stack + 6790;
 Libderiv->deriv_classes[3][6][2] = int_stack + 7000;
 Libderiv->deriv_classes[3][7][2] = int_stack + 7280;
 Libderiv->deriv_classes[3][4][1] = int_stack + 7640;
 Libderiv->deriv_classes[3][5][1] = int_stack + 7790;
 Libderiv->deriv_classes[3][6][1] = int_stack + 8000;
 Libderiv->deriv_classes[3][7][1] = int_stack + 8280;
 Libderiv->deriv_classes[3][4][0] = int_stack + 8640;
 Libderiv->deriv_classes[3][5][0] = int_stack + 8790;
 Libderiv->deriv_classes[3][6][0] = int_stack + 9000;
 Libderiv->deriv_classes[3][7][0] = int_stack + 9280;
 memset(int_stack,0,77120);

 Libderiv->dvrr_stack = int_stack + 24940;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0gf(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9640,int_stack+5300,int_stack+5000,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+10090,int_stack+5720,int_stack+5300,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+10720,int_stack+10090,int_stack+9640,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11620,int_stack+150,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5000,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12070,int_stack+360,int_stack+150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5300,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+12700,int_stack+12070,int_stack+11620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9640,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+640,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5720,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+14440,int_stack+13600,int_stack+12070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10090,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13600,int_stack+1150,int_stack+1000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5000, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11620,int_stack+1360,int_stack+1150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5300, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+11620,int_stack+13600, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9640, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+1640,int_stack+1360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5720, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+15700,int_stack+13600,int_stack+11620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10090, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11620,int_stack+2150,int_stack+2000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5000, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12070,int_stack+2360,int_stack+2150, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5300, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+900,int_stack+12070,int_stack+11620, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9640, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+2640,int_stack+2360, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5720, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+16960,int_stack+13600,int_stack+12070, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10090, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13600,int_stack+3150,int_stack+3000, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11620,int_stack+3360,int_stack+3150, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1800,int_stack+11620,int_stack+13600, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+9640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+3640,int_stack+3360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+5720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+2700,int_stack+13600,int_stack+11620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+11620,int_stack+4150,int_stack+4000, 0.0, zero_stack, 1.0, int_stack+5000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+12070,int_stack+4360,int_stack+4150, 0.0, zero_stack, 1.0, int_stack+5300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+18220,int_stack+12070,int_stack+11620, 0.0, zero_stack, 1.0, int_stack+9640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+4640,int_stack+4360, 0.0, zero_stack, 1.0, int_stack+5720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+19120,int_stack+13600,int_stack+12070, 0.0, zero_stack, 1.0, int_stack+10090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13600,int_stack+5510,int_stack+5150, 1.0, int_stack+5000, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+11620,int_stack+6000,int_stack+5510, 1.0, int_stack+5300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3960,int_stack+11620,int_stack+13600, 1.0, int_stack+9640, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|ip) ---*/
   d1hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+6280,int_stack+6000, 1.0, int_stack+5720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|hd) ---*/
   d1hrr3_build_hd(Libderiv->CD,int_stack+4860,int_stack+13600,int_stack+11620, 1.0, int_stack+10090, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11620,int_stack+6790,int_stack+6640,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12070,int_stack+7000,int_stack+6790,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+9640,int_stack+12070,int_stack+11620,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+7280,int_stack+7000,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+6120,int_stack+13600,int_stack+12070,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+13600,int_stack+7790,int_stack+7640,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+11620,int_stack+8000,int_stack+7790,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+20380,int_stack+11620,int_stack+13600,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+8280,int_stack+8000,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+7380,int_stack+13600,int_stack+11620,10);
 /*--- compute (f0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+11620,int_stack+8790,int_stack+8640,10);
 /*--- compute (f0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+12070,int_stack+9000,int_stack+8790,10);
 /*--- compute (f0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+21280,int_stack+12070,int_stack+11620,10);
 /*--- compute (f0|ip) ---*/
   hrr3_build_ip(Libderiv->CD,int_stack+13600,int_stack+9280,int_stack+9000,10);
 /*--- compute (f0|hd) ---*/
   hrr3_build_hd(Libderiv->CD,int_stack+22180,int_stack+13600,int_stack+12070,10);
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+23440,int_stack+14440,int_stack+12700, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10720,10);
     Libderiv->ABCD[11] = int_stack + 23440;
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+11620,int_stack+15700,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10720, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 11620;
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+13120,int_stack+16960,int_stack+900, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10720, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 13120;
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+0,int_stack+2700,int_stack+1800, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+1500,int_stack+19120,int_stack+18220, 0.0, zero_stack, 1.0, int_stack+10720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 1500;
 /*--- compute (f0|gf) ---*/
   d1hrr3_build_gf(Libderiv->CD,int_stack+14620,int_stack+4860,int_stack+3960, 1.0, int_stack+10720, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 14620;
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+3000,int_stack+6120,int_stack+9640,10);
     Libderiv->ABCD[2] = int_stack + 3000;
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+4500,int_stack+7380,int_stack+20380,10);
     Libderiv->ABCD[1] = int_stack + 4500;
 /*--- compute (f0|gf) ---*/
   hrr3_build_gf(Libderiv->CD,int_stack+6000,int_stack+22180,int_stack+21280,10);
     Libderiv->ABCD[0] = int_stack + 6000;

}
