#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0gd(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|gd) integrals */

void d1hrr_order_p0gd(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[1][5][11] = int_stack + 45;
 Libderiv->deriv_classes[1][6][11] = int_stack + 108;
 Libderiv->deriv_classes[1][4][10] = int_stack + 192;
 Libderiv->deriv_classes[1][5][10] = int_stack + 237;
 Libderiv->deriv_classes[1][6][10] = int_stack + 300;
 Libderiv->deriv_classes[1][4][9] = int_stack + 384;
 Libderiv->deriv_classes[1][5][9] = int_stack + 429;
 Libderiv->deriv_classes[1][6][9] = int_stack + 492;
 Libderiv->deriv_classes[1][4][8] = int_stack + 576;
 Libderiv->deriv_classes[1][5][8] = int_stack + 621;
 Libderiv->deriv_classes[1][6][8] = int_stack + 684;
 Libderiv->deriv_classes[1][4][7] = int_stack + 768;
 Libderiv->deriv_classes[1][5][7] = int_stack + 813;
 Libderiv->deriv_classes[1][6][7] = int_stack + 876;
 Libderiv->dvrr_classes[1][4] = int_stack + 960;
 Libderiv->deriv_classes[1][4][6] = int_stack + 1005;
 Libderiv->dvrr_classes[1][5] = int_stack + 1050;
 Libderiv->deriv_classes[1][5][6] = int_stack + 1113;
 Libderiv->deriv_classes[1][6][6] = int_stack + 1176;
 Libderiv->deriv_classes[1][4][2] = int_stack + 1260;
 Libderiv->deriv_classes[1][5][2] = int_stack + 1305;
 Libderiv->deriv_classes[1][6][2] = int_stack + 1368;
 Libderiv->deriv_classes[1][4][1] = int_stack + 1452;
 Libderiv->deriv_classes[1][5][1] = int_stack + 1497;
 Libderiv->deriv_classes[1][6][1] = int_stack + 1560;
 Libderiv->deriv_classes[1][4][0] = int_stack + 1644;
 Libderiv->deriv_classes[1][5][0] = int_stack + 1689;
 Libderiv->deriv_classes[1][6][0] = int_stack + 1752;
 memset(int_stack,0,14688);

 Libderiv->dvrr_stack = int_stack + 3780;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0gd(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1836,int_stack+1050,int_stack+960,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1971,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+960,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2106,int_stack+108,int_stack+45, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+237,int_stack+192, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+960, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2295,int_stack+300,int_stack+237, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+135,int_stack+429,int_stack+384, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+960, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2484,int_stack+492,int_stack+429, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+270,int_stack+621,int_stack+576, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+405,int_stack+684,int_stack+621, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+594,int_stack+813,int_stack+768, 0.0, zero_stack, 1.0, int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2673,int_stack+876,int_stack+813, 0.0, zero_stack, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+729,int_stack+1113,int_stack+1005, 1.0, int_stack+960, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|hp) ---*/
   d1hrr3_build_hp(Libderiv->CD,int_stack+2862,int_stack+1176,int_stack+1113, 1.0, int_stack+1050, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+864,int_stack+1305,int_stack+1260,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+999,int_stack+1368,int_stack+1305,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1188,int_stack+1497,int_stack+1452,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+3051,int_stack+1560,int_stack+1497,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1323,int_stack+1689,int_stack+1644,3);
 /*--- compute (p0|hp) ---*/
   hrr3_build_hp(Libderiv->CD,int_stack+1458,int_stack+1752,int_stack+1689,3);
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3240,int_stack+2106,int_stack+1971, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1836,3);
     Libderiv->ABCD[11] = int_stack + 3240;
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+1971,int_stack+2295,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1836, 0.0, zero_stack,3);
     Libderiv->ABCD[10] = int_stack + 1971;
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+3510,int_stack+2484,int_stack+135, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1836, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[9] = int_stack + 3510;
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+0,int_stack+405,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+270,int_stack+2673,int_stack+594, 0.0, zero_stack, 1.0, int_stack+1836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[7] = int_stack + 270;
 /*--- compute (p0|gd) ---*/
   d1hrr3_build_gd(Libderiv->CD,int_stack+2241,int_stack+2862,int_stack+729, 1.0, int_stack+1836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
     Libderiv->ABCD[6] = int_stack + 2241;
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+540,int_stack+999,int_stack+864,3);
     Libderiv->ABCD[2] = int_stack + 540;
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+810,int_stack+3051,int_stack+1188,3);
     Libderiv->ABCD[1] = int_stack + 810;
 /*--- compute (p0|gd) ---*/
   hrr3_build_gd(Libderiv->CD,int_stack+2511,int_stack+1458,int_stack+1323,3);
     Libderiv->ABCD[0] = int_stack + 2511;

}
