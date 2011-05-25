#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_pppp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (pp|pp) integrals */

void d1hrr_order_pppp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][1][11] = int_stack + 0;
 Libderiv->deriv_classes[1][2][11] = int_stack + 9;
 Libderiv->deriv_classes[2][1][11] = int_stack + 27;
 Libderiv->deriv_classes[2][2][11] = int_stack + 45;
 Libderiv->deriv_classes[1][1][10] = int_stack + 81;
 Libderiv->deriv_classes[1][2][10] = int_stack + 90;
 Libderiv->deriv_classes[2][1][10] = int_stack + 108;
 Libderiv->deriv_classes[2][2][10] = int_stack + 126;
 Libderiv->deriv_classes[1][1][9] = int_stack + 162;
 Libderiv->deriv_classes[1][2][9] = int_stack + 171;
 Libderiv->deriv_classes[2][1][9] = int_stack + 189;
 Libderiv->deriv_classes[2][2][9] = int_stack + 207;
 Libderiv->deriv_classes[1][1][8] = int_stack + 243;
 Libderiv->deriv_classes[1][2][8] = int_stack + 252;
 Libderiv->deriv_classes[2][1][8] = int_stack + 270;
 Libderiv->deriv_classes[2][2][8] = int_stack + 288;
 Libderiv->deriv_classes[1][1][7] = int_stack + 324;
 Libderiv->deriv_classes[1][2][7] = int_stack + 333;
 Libderiv->deriv_classes[2][1][7] = int_stack + 351;
 Libderiv->deriv_classes[2][2][7] = int_stack + 369;
 Libderiv->deriv_classes[1][1][6] = int_stack + 405;
 Libderiv->deriv_classes[1][2][6] = int_stack + 414;
 Libderiv->dvrr_classes[2][1] = int_stack + 432;
 Libderiv->deriv_classes[2][1][6] = int_stack + 450;
 Libderiv->deriv_classes[2][2][6] = int_stack + 468;
 Libderiv->deriv_classes[1][1][2] = int_stack + 504;
 Libderiv->deriv_classes[1][2][2] = int_stack + 513;
 Libderiv->deriv_classes[2][1][2] = int_stack + 531;
 Libderiv->deriv_classes[2][2][2] = int_stack + 549;
 Libderiv->deriv_classes[1][1][1] = int_stack + 585;
 Libderiv->deriv_classes[1][2][1] = int_stack + 594;
 Libderiv->deriv_classes[2][1][1] = int_stack + 612;
 Libderiv->deriv_classes[2][2][1] = int_stack + 630;
 Libderiv->dvrr_classes[1][1] = int_stack + 666;
 Libderiv->dvrr_classes[1][2] = int_stack + 675;
 Libderiv->deriv_classes[1][1][0] = int_stack + 693;
 Libderiv->deriv_classes[1][2][0] = int_stack + 702;
 Libderiv->deriv_classes[2][1][0] = int_stack + 720;
 Libderiv->deriv_classes[2][2][0] = int_stack + 738;
 memset(int_stack,0,6192);

 Libderiv->dvrr_stack = int_stack + 936;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_pppp(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+774,int_stack+9,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+666,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+801,int_stack+45,int_stack+27, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+0,int_stack+90,int_stack+81, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+666, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+27,int_stack+126,int_stack+108, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+81,int_stack+171,int_stack+162, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+666, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+108,int_stack+207,int_stack+189, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+162,int_stack+252,int_stack+243, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+189,int_stack+288,int_stack+270, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+243,int_stack+333,int_stack+324, 0.0, zero_stack, 1.0, int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+270,int_stack+369,int_stack+351, 0.0, zero_stack, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+324,int_stack+414,int_stack+405, 1.0, int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,3);
 /*--- compute (d0|pp) ---*/
   d1hrr3_build_pp(Libderiv->CD,int_stack+351,int_stack+468,int_stack+450, 1.0, int_stack+432, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+405,int_stack+675,int_stack+666,3);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+666,int_stack+513,int_stack+504,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+432,int_stack+549,int_stack+531,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+486,int_stack+594,int_stack+585,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+513,int_stack+630,int_stack+612,6);
 /*--- compute (p0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+567,int_stack+702,int_stack+693,3);
 /*--- compute (d0|pp) ---*/
   hrr3_build_pp(Libderiv->CD,int_stack+594,int_stack+738,int_stack+720,6);
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+693,int_stack+801,int_stack+774,9);
     Libderiv->ABCD[11] = int_stack + 693;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+774,int_stack+27,int_stack+0,9);
     Libderiv->ABCD[10] = int_stack + 774;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+0,int_stack+108,int_stack+81,9);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+81,int_stack+189,int_stack+162,9);
     Libderiv->ABCD[8] = int_stack + 81;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+162,int_stack+270,int_stack+243,9);
     Libderiv->ABCD[7] = int_stack + 162;
 /*--- compute (pp|pp) ---*/
   hrr1_build_pp(Libderiv->AB,int_stack+243,int_stack+351,int_stack+324,9);
     Libderiv->ABCD[6] = int_stack + 243;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+324,int_stack+432,int_stack+666, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[2] = int_stack + 324;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+855,int_stack+513,int_stack+486, 0.0, zero_stack, 1.0, int_stack+405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[1] = int_stack + 855;
 /*--- compute (pp|pp) ---*/
   d1hrr1_build_pp(Libderiv->AB,int_stack+432,int_stack+594,int_stack+567, 1.0, int_stack+405, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,9);
     Libderiv->ABCD[0] = int_stack + 432;

}
