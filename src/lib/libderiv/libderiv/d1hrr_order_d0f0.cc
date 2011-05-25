#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0f0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|f0) integrals */

void d1hrr_order_d0f0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][3][11] = int_stack + 0;
 Libderiv->deriv_classes[2][3][10] = int_stack + 60;
 Libderiv->deriv_classes[2][3][9] = int_stack + 120;
 Libderiv->deriv_classes[2][3][8] = int_stack + 180;
 Libderiv->deriv_classes[2][3][7] = int_stack + 240;
 Libderiv->deriv_classes[2][3][6] = int_stack + 300;
 Libderiv->deriv_classes[2][3][2] = int_stack + 360;
 Libderiv->deriv_classes[2][3][1] = int_stack + 420;
 Libderiv->deriv_classes[2][3][0] = int_stack + 480;
 memset(int_stack,0,4320);

 Libderiv->dvrr_stack = int_stack + 540;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0f0(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[10] = int_stack + 60;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[9] = int_stack + 120;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[8] = int_stack + 180;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[7] = int_stack + 240;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[6] = int_stack + 300;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[2] = int_stack + 360;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[1] = int_stack + 420;
 /*--- compute (d0|f0) ---*/
     Libderiv->ABCD[0] = int_stack + 480;

}
