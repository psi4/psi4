#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_ppgp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|gp) integrals */

void d1hrr_order_ppgp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][4][11] = int_stack + 0;
 Libderiv->deriv_classes[1][5][11] = int_stack + 45;
 Libderiv->deriv_classes[2][4][11] = int_stack + 108;
 Libderiv->deriv_classes[2][5][11] = int_stack + 198;
 Libderiv->deriv_classes[1][4][10] = int_stack + 324;
 Libderiv->deriv_classes[1][5][10] = int_stack + 369;
 Libderiv->deriv_classes[2][4][10] = int_stack + 432;
 Libderiv->deriv_classes[2][5][10] = int_stack + 522;
 Libderiv->deriv_classes[1][4][9] = int_stack + 648;
 Libderiv->deriv_classes[1][5][9] = int_stack + 693;
 Libderiv->deriv_classes[2][4][9] = int_stack + 756;
 Libderiv->deriv_classes[2][5][9] = int_stack + 846;
 Libderiv->deriv_classes[1][4][8] = int_stack + 972;
 Libderiv->deriv_classes[1][5][8] = int_stack + 1017;
 Libderiv->deriv_classes[2][4][8] = int_stack + 1080;
 Libderiv->deriv_classes[2][5][8] = int_stack + 1170;
 Libderiv->deriv_classes[1][4][7] = int_stack + 1296;
 Libderiv->deriv_classes[1][5][7] = int_stack + 1341;
 Libderiv->deriv_classes[2][4][7] = int_stack + 1404;
 Libderiv->deriv_classes[2][5][7] = int_stack + 1494;
 Libderiv->deriv_classes[1][4][6] = int_stack + 1620;
 Libderiv->deriv_classes[1][5][6] = int_stack + 1665;
 Libderiv->dvrr_classes[2][4] = int_stack + 1728;
 Libderiv->deriv_classes[2][4][6] = int_stack + 1818;
 Libderiv->deriv_classes[2][5][6] = int_stack + 1908;
 Libderiv->deriv_classes[1][4][2] = int_stack + 2034;
 Libderiv->deriv_classes[1][5][2] = int_stack + 2079;
 Libderiv->deriv_classes[2][4][2] = int_stack + 2142;
 Libderiv->deriv_classes[2][5][2] = int_stack + 2232;
 Libderiv->deriv_classes[1][4][1] = int_stack + 2358;
 Libderiv->deriv_classes[1][5][1] = int_stack + 2403;
 Libderiv->deriv_classes[2][4][1] = int_stack + 2466;
 Libderiv->deriv_classes[2][5][1] = int_stack + 2556;
 Libderiv->dvrr_classes[1][4] = int_stack + 2682;
 Libderiv->dvrr_classes[1][5] = int_stack + 2727;
 Libderiv->deriv_classes[1][4][0] = int_stack + 2790;
 Libderiv->deriv_classes[1][5][0] = int_stack + 2835;
 Libderiv->deriv_classes[2][4][0] = int_stack + 2898;
 Libderiv->deriv_classes[2][5][0] = int_stack + 2988;
 memset(int_stack,0,24912);

 Libderiv->dvrr_stack = int_stack + 4464;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_ppgp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3114,int_stack+45,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2682,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3249,int_stack+198,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1728,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+369,int_stack+324, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2682, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+135,int_stack+522,int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1728, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+405,int_stack+693,int_stack+648, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2682, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3519,int_stack+846,int_stack+756, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1728, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+540,int_stack+1017,int_stack+972, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+2682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+675,int_stack+1170,int_stack+1080, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+945,int_stack+1341,int_stack+1296, 0.0, zero_stack, 1.0, int_stack+2682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1080,int_stack+1494,int_stack+1404, 0.0, zero_stack, 1.0, int_stack+1728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1350,int_stack+1665,int_stack+1620, 1.0, int_stack+2682, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3789,int_stack+1908,int_stack+1818, 1.0, int_stack+1728, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1485,int_stack+2727,int_stack+2682,3);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1620,int_stack+2079,int_stack+2034,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+1755,int_stack+2232,int_stack+2142,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2025,int_stack+2403,int_stack+2358,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2160,int_stack+2556,int_stack+2466,6);
 /*--- compute (p0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2430,int_stack+2835,int_stack+2790,3);
 /*--- compute (d0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+2565,int_stack+2988,int_stack+2898,6);
 /*--- compute (pp|gp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+4059,int_stack+3249,int_stack+3114,45);
     Libderiv->ABCD[11] = int_stack + 4059;
 /*--- compute (pp|gp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+2835,int_stack+135,int_stack+0,45);
     Libderiv->ABCD[10] = int_stack + 2835;
 /*--- compute (pp|gp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+3519,int_stack+405,45);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (pp|gp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+3240,int_stack+675,int_stack+540,45);
     Libderiv->ABCD[8] = int_stack + 3240;
 /*--- compute (pp|gp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+405,int_stack+1080,int_stack+945,45);
     Libderiv->ABCD[7] = int_stack + 405;
 /*--- compute (pp|gp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+810,int_stack+3789,int_stack+1350,45);
     Libderiv->ABCD[6] = int_stack + 810;
 /*--- compute (pp|gp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+3645,int_stack+1755,int_stack+1620, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[2] = int_stack + 3645;
 /*--- compute (pp|gp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+1620,int_stack+2160,int_stack+2025, 0.0, zero_stack, 1.0, int_stack+1485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[1] = int_stack + 1620;
 /*--- compute (pp|gp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+2025,int_stack+2565,int_stack+2430, 1.0, int_stack+1485, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[0] = int_stack + 2025;

}
