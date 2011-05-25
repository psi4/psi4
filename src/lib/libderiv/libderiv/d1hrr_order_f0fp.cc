#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0fp(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|fp) integrals */

void d1hrr_order_f0fp(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][4][11] = int_stack + 100;
 Libderiv->deriv_classes[3][3][10] = int_stack + 250;
 Libderiv->deriv_classes[3][4][10] = int_stack + 350;
 Libderiv->deriv_classes[3][3][9] = int_stack + 500;
 Libderiv->deriv_classes[3][4][9] = int_stack + 600;
 Libderiv->deriv_classes[3][3][8] = int_stack + 750;
 Libderiv->deriv_classes[3][4][8] = int_stack + 850;
 Libderiv->deriv_classes[3][3][7] = int_stack + 1000;
 Libderiv->deriv_classes[3][4][7] = int_stack + 1100;
 Libderiv->dvrr_classes[3][3] = int_stack + 1250;
 Libderiv->deriv_classes[3][3][6] = int_stack + 1350;
 Libderiv->deriv_classes[3][4][6] = int_stack + 1450;
 Libderiv->deriv_classes[3][3][2] = int_stack + 1600;
 Libderiv->deriv_classes[3][4][2] = int_stack + 1700;
 Libderiv->deriv_classes[3][3][1] = int_stack + 1850;
 Libderiv->deriv_classes[3][4][1] = int_stack + 1950;
 Libderiv->deriv_classes[3][3][0] = int_stack + 2100;
 Libderiv->deriv_classes[3][4][0] = int_stack + 2200;
 memset(int_stack,0,18800);

 Libderiv->dvrr_stack = int_stack + 2950;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0fp(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2350,int_stack+100,int_stack+0, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1250,10);
     Libderiv->ABCD[11] = int_stack + 2350;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+2650,int_stack+350,int_stack+250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1250, 0.0, zero_stack,10);
     Libderiv->ABCD[10] = int_stack + 2650;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+0,int_stack+600,int_stack+500, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1250, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[9] = int_stack + 0;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+300,int_stack+850,int_stack+750, 0.0, zero_stack, 0.0, zero_stack, 1.0, int_stack+1250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[8] = int_stack + 300;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+600,int_stack+1100,int_stack+1000, 0.0, zero_stack, 1.0, int_stack+1250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[7] = int_stack + 600;
 /*--- compute (f0|fp) ---*/
   d1hrr3_build_fp(Libderiv->CD,int_stack+900,int_stack+1450,int_stack+1350, 1.0, int_stack+1250, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack, 0.0, zero_stack,10);
     Libderiv->ABCD[6] = int_stack + 900;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1200,int_stack+1700,int_stack+1600,10);
     Libderiv->ABCD[2] = int_stack + 1200;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1500,int_stack+1950,int_stack+1850,10);
     Libderiv->ABCD[1] = int_stack + 1500;
 /*--- compute (f0|fp) ---*/
   hrr3_build_fp(Libderiv->CD,int_stack+1800,int_stack+2200,int_stack+2100,10);
     Libderiv->ABCD[0] = int_stack + 1800;

}
