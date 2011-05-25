#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0fd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|fd) integrals */

void d1hrr_order_d0fd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][11] = int_stack + 60;
 Libderiv->deriv_classes[2][5][11] = int_stack + 150;
 Libderiv->deriv_classes[2][3][10] = int_stack + 276;
 Libderiv->deriv_classes[2][4][10] = int_stack + 336;
 Libderiv->deriv_classes[2][5][10] = int_stack + 426;
 Libderiv->deriv_classes[2][3][9] = int_stack + 552;
 Libderiv->deriv_classes[2][4][9] = int_stack + 612;
 Libderiv->deriv_classes[2][5][9] = int_stack + 702;
 Libderiv->deriv_classes[2][3][8] = int_stack + 828;
 Libderiv->deriv_classes[2][4][8] = int_stack + 888;
 Libderiv->deriv_classes[2][5][8] = int_stack + 978;
 Libderiv->deriv_classes[2][3][7] = int_stack + 1104;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1164;
 Libderiv->deriv_classes[2][5][7] = int_stack + 1254;
 Libderiv->dvrr_classes[2][3] = int_stack + 1380;
 Libderiv->deriv_classes[2][3][6] = int_stack + 1440;
 Libderiv->dvrr_classes[2][4] = int_stack + 1500;
 Libderiv->deriv_classes[2][4][6] = int_stack + 1590;
 Libderiv->deriv_classes[2][5][6] = int_stack + 1680;
 Libderiv->deriv_classes[2][3][2] = int_stack + 1806;
 Libderiv->deriv_classes[2][4][2] = int_stack + 1866;
 Libderiv->deriv_classes[2][5][2] = int_stack + 1956;
 Libderiv->deriv_classes[2][3][1] = int_stack + 2082;
 Libderiv->deriv_classes[2][4][1] = int_stack + 2142;
 Libderiv->deriv_classes[2][5][1] = int_stack + 2232;
 Libderiv->deriv_classes[2][3][0] = int_stack + 2358;
 Libderiv->deriv_classes[2][4][0] = int_stack + 2418;
 Libderiv->deriv_classes[2][5][0] = int_stack + 2508;
 memset(int_stack,0,21072);

 Libderiv->dvrr_stack = int_stack + 4614;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0fd(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2634,int_stack+1500,int_stack+1380,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2814,int_stack+60,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1380,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2994,int_stack+150,int_stack+60, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+336,int_stack+276, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1380, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3264,int_stack+426,int_stack+336, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+180,int_stack+612,int_stack+552, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1380, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3534,int_stack+702,int_stack+612, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+888,int_stack+828, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+978,int_stack+888, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+810,int_stack+1164,int_stack+1104, 0.0, zero_stack, 1.0, int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3804,int_stack+1254,int_stack+1164, 0.0, zero_stack, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+990,int_stack+1590,int_stack+1440, 1.0, int_stack+1380, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1170,int_stack+1680,int_stack+1590, 1.0, int_stack+1500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1440,int_stack+1866,int_stack+1806,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4074,int_stack+1956,int_stack+1866,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1620,int_stack+2142,int_stack+2082,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1800,int_stack+2232,int_stack+2142,6);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2070,int_stack+2418,int_stack+2358,6);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+4344,int_stack+2508,int_stack+2418,6);
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2250,int_stack+2994,int_stack+2814, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2634,6);
     Libderiv->ABCD[11] = int_stack + 2250;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+2814,int_stack+3264,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2634, 0.0, zero_stack,6);
     Libderiv->ABCD[10] = int_stack + 2814;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3174,int_stack+3534,int_stack+180, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2634, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[9] = int_stack + 3174;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+0,int_stack+540,int_stack+360, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+360,int_stack+3804,int_stack+810, 0.0, zero_stack, 1.0, int_stack+2634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[7] = int_stack + 360;
 /*--- compute (d0|fd) ---*/
   d1hrr3_build_fd(Libderiv->CD,int_stack+3534,int_stack+1170,int_stack+990, 1.0, int_stack+2634, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
     Libderiv->ABCD[6] = int_stack + 3534;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+720,int_stack+4074,int_stack+1440,6);
     Libderiv->ABCD[2] = int_stack + 720;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1080,int_stack+1800,int_stack+1620,6);
     Libderiv->ABCD[1] = int_stack + 1080;
 /*--- compute (d0|fd) ---*/
   hrr3_build_fd(Libderiv->CD,int_stack+1440,int_stack+4344,int_stack+2070,6);
     Libderiv->ABCD[0] = int_stack + 1440;

}
