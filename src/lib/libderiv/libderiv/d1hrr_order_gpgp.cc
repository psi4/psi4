#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_gpgp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (gp|gp) integrals */

void d1hrr_order_gpgp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][5][11] = int_stack + 225;
 Libderiv->deriv_classes[5][4][11] = int_stack + 540;
 Libderiv->deriv_classes[5][5][11] = int_stack + 855;
 Libderiv->deriv_classes[4][4][10] = int_stack + 1296;
 Libderiv->deriv_classes[4][5][10] = int_stack + 1521;
 Libderiv->deriv_classes[5][4][10] = int_stack + 1836;
 Libderiv->deriv_classes[5][5][10] = int_stack + 2151;
 Libderiv->deriv_classes[4][4][9] = int_stack + 2592;
 Libderiv->deriv_classes[4][5][9] = int_stack + 2817;
 Libderiv->deriv_classes[5][4][9] = int_stack + 3132;
 Libderiv->deriv_classes[5][5][9] = int_stack + 3447;
 Libderiv->deriv_classes[4][4][8] = int_stack + 3888;
 Libderiv->deriv_classes[4][5][8] = int_stack + 4113;
 Libderiv->deriv_classes[5][4][8] = int_stack + 4428;
 Libderiv->deriv_classes[5][5][8] = int_stack + 4743;
 Libderiv->deriv_classes[4][4][7] = int_stack + 5184;
 Libderiv->deriv_classes[4][5][7] = int_stack + 5409;
 Libderiv->deriv_classes[5][4][7] = int_stack + 5724;
 Libderiv->deriv_classes[5][5][7] = int_stack + 6039;
 Libderiv->deriv_classes[4][4][6] = int_stack + 6480;
 Libderiv->deriv_classes[4][5][6] = int_stack + 6705;
 Libderiv->dvrr_classes[5][4] = int_stack + 7020;
 Libderiv->deriv_classes[5][4][6] = int_stack + 7335;
 Libderiv->deriv_classes[5][5][6] = int_stack + 7650;
 Libderiv->deriv_classes[4][4][2] = int_stack + 8091;
 Libderiv->deriv_classes[4][5][2] = int_stack + 8316;
 Libderiv->deriv_classes[5][4][2] = int_stack + 8631;
 Libderiv->deriv_classes[5][5][2] = int_stack + 8946;
 Libderiv->deriv_classes[4][4][1] = int_stack + 9387;
 Libderiv->deriv_classes[4][5][1] = int_stack + 9612;
 Libderiv->deriv_classes[5][4][1] = int_stack + 9927;
 Libderiv->deriv_classes[5][5][1] = int_stack + 10242;
 Libderiv->dvrr_classes[4][4] = int_stack + 10683;
 Libderiv->dvrr_classes[4][5] = int_stack + 10908;
 Libderiv->deriv_classes[4][4][0] = int_stack + 11223;
 Libderiv->deriv_classes[4][5][0] = int_stack + 11448;
 Libderiv->deriv_classes[5][4][0] = int_stack + 11763;
 Libderiv->deriv_classes[5][5][0] = int_stack + 12078;
 memset(int_stack,0,100152);

 Libderiv->dvrr_stack = int_stack + 23049;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_gpgp(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+12519,int_stack+225,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10683,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+13194,int_stack+855,int_stack+540, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7020,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+0,int_stack+1521,int_stack+1296, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10683, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+675,int_stack+2151,int_stack+1836, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7020, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+1620,int_stack+2817,int_stack+2592, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10683, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+14139,int_stack+3447,int_stack+3132, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7020, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2295,int_stack+4113,int_stack+3888, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+10683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+2970,int_stack+4743,int_stack+4428, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+7020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+3915,int_stack+5409,int_stack+5184, 0.0, zero_stack, 1.0, int_stack+10683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+4590,int_stack+6039,int_stack+5724, 0.0, zero_stack, 1.0, int_stack+7020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+5535,int_stack+6705,int_stack+6480, 1.0, int_stack+10683, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,15);
 /*--- compute (h0|gp) ---*/
   d1hrr3_build_gp(Libderiv->CD,int_stack+15084,int_stack+7650,int_stack+7335, 1.0, int_stack+7020, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6210,int_stack+10908,int_stack+10683,15);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+6885,int_stack+8316,int_stack+8091,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+7560,int_stack+8946,int_stack+8631,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+8505,int_stack+9612,int_stack+9387,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+16029,int_stack+10242,int_stack+9927,21);
 /*--- compute (g0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9180,int_stack+11448,int_stack+11223,15);
 /*--- compute (h0|gp) ---*/
   hrr3_build_gp(Libderiv->CD,int_stack+9855,int_stack+12078,int_stack+11763,21);
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+16974,int_stack+13194,int_stack+12519,45);
     Libderiv->ABCD[11] = int_stack + 16974;
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+10800,int_stack+675,int_stack+0,45);
     Libderiv->ABCD[10] = int_stack + 10800;
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+18999,int_stack+14139,int_stack+1620,45);
     Libderiv->ABCD[9] = int_stack + 18999;
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+0,int_stack+2970,int_stack+2295,45);
     Libderiv->ABCD[8] = int_stack + 0;
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+12825,int_stack+4590,int_stack+3915,45);
     Libderiv->ABCD[7] = int_stack + 12825;
 /*--- compute (gp|gp) ---*/
   hrr1_build_gp(Libderiv->AB,int_stack+2025,int_stack+15084,int_stack+5535,45);
     Libderiv->ABCD[6] = int_stack + 2025;
 /*--- compute (gp|gp) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+4050,int_stack+7560,int_stack+6885, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[2] = int_stack + 4050;
 /*--- compute (gp|gp) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+21024,int_stack+16029,int_stack+8505, 0.0, zero_stack, 1.0, int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[1] = int_stack + 21024;
 /*--- compute (gp|gp) ---*/
   d1hrr1_build_gp(Libderiv->AB,int_stack+6885,int_stack+9855,int_stack+9180, 1.0, int_stack+6210, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,45);
     Libderiv->ABCD[0] = int_stack + 6885;

}
