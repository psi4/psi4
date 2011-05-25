#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppfp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|fp) integrals */

void d1hrr_order_ppfp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][3][11] = int_stack + 0;
 Libderiv->deriv_classes[1][4][11] = int_stack + 30;
 Libderiv->deriv_classes[2][3][11] = int_stack + 75;
 Libderiv->deriv_classes[2][4][11] = int_stack + 135;
 Libderiv->deriv_classes[1][3][10] = int_stack + 225;
 Libderiv->deriv_classes[1][4][10] = int_stack + 255;
 Libderiv->deriv_classes[2][3][10] = int_stack + 300;
 Libderiv->deriv_classes[2][4][10] = int_stack + 360;
 Libderiv->deriv_classes[1][3][9] = int_stack + 450;
 Libderiv->deriv_classes[1][4][9] = int_stack + 480;
 Libderiv->deriv_classes[2][3][9] = int_stack + 525;
 Libderiv->deriv_classes[2][4][9] = int_stack + 585;
 Libderiv->deriv_classes[1][3][8] = int_stack + 675;
 Libderiv->deriv_classes[1][4][8] = int_stack + 705;
 Libderiv->deriv_classes[2][3][8] = int_stack + 750;
 Libderiv->deriv_classes[2][4][8] = int_stack + 810;
 Libderiv->deriv_classes[1][3][7] = int_stack + 900;
 Libderiv->deriv_classes[1][4][7] = int_stack + 930;
 Libderiv->deriv_classes[2][3][7] = int_stack + 975;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1035;
 Libderiv->deriv_classes[1][3][6] = int_stack + 1125;
 Libderiv->deriv_classes[1][4][6] = int_stack + 1155;
 Libderiv->dvrr_classes[2][3] = int_stack + 1200;
 Libderiv->deriv_classes[2][3][6] = int_stack + 1260;
 Libderiv->deriv_classes[2][4][6] = int_stack + 1320;
 Libderiv->deriv_classes[1][3][2] = int_stack + 1410;
 Libderiv->deriv_classes[1][4][2] = int_stack + 1440;
 Libderiv->deriv_classes[2][3][2] = int_stack + 1485;
 Libderiv->deriv_classes[2][4][2] = int_stack + 1545;
 Libderiv->deriv_classes[1][3][1] = int_stack + 1635;
 Libderiv->deriv_classes[1][4][1] = int_stack + 1665;
 Libderiv->deriv_classes[2][3][1] = int_stack + 1710;
 Libderiv->deriv_classes[2][4][1] = int_stack + 1770;
 Libderiv->dvrr_classes[1][3] = int_stack + 1860;
 Libderiv->dvrr_classes[1][4] = int_stack + 1890;
 Libderiv->deriv_classes[1][3][0] = int_stack + 1935;
 Libderiv->deriv_classes[1][4][0] = int_stack + 1965;
 Libderiv->deriv_classes[2][3][0] = int_stack + 2010;
 Libderiv->deriv_classes[2][4][0] = int_stack + 2070;
 memset(int_stack,0,17280);

 Libderiv->dvrr_stack = int_stack + 3060;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppfp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2160,int_stack+30,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1860,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2250,int_stack+135,int_stack+75, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+255,int_stack+225, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1860, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+90,int_stack+360,int_stack+300, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+270,int_stack+480,int_stack+450, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1860, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2430,int_stack+585,int_stack+525, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+360,int_stack+705,int_stack+675, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+450,int_stack+810,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+630,int_stack+930,int_stack+900, 0.0, zero_stack, 1.0, int_stack+1860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+720,int_stack+1035,int_stack+975, 0.0, zero_stack, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1155,int_stack+1125, 1.0, int_stack+1860, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+990,int_stack+1320,int_stack+1260, 1.0, int_stack+1200, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1170,int_stack+1890,int_stack+1860,3);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1260,int_stack+1440,int_stack+1410,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+2610,int_stack+1545,int_stack+1485,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1350,int_stack+1665,int_stack+1635,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1440,int_stack+1770,int_stack+1710,6);
 /*--- compute (p0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1620,int_stack+1965,int_stack+1935,3);
 /*--- compute (d0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1710,int_stack+2070,int_stack+2010,6);
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+1890,int_stack+2250,int_stack+2160,30);
     Libderiv->ABCD[11] = int_stack + 1890;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2160,int_stack+90,int_stack+0,30);
     Libderiv->ABCD[10] = int_stack + 2160;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+2430,int_stack+270,30);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2790,int_stack+450,int_stack+360,30);
     Libderiv->ABCD[8] = int_stack + 2790;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+270,int_stack+720,int_stack+630,30);
     Libderiv->ABCD[7] = int_stack + 270;
 /*--- compute (pp|fp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+540,int_stack+990,int_stack+900,30);
     Libderiv->ABCD[6] = int_stack + 540;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+810,int_stack+2610,int_stack+1260, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[2] = int_stack + 810;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2430,int_stack+1440,int_stack+1350, 0.0, zero_stack, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[1] = int_stack + 2430;
 /*--- compute (pp|fp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1260,int_stack+1710,int_stack+1620, 1.0, int_stack+1170, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,30);
     Libderiv->ABCD[0] = int_stack + 1260;

}
