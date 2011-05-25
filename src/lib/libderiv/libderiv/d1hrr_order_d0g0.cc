#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0g0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|g0) integrals */

void d1hrr_order_d0g0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][4][11] = int_stack + 0;
 Libderiv->deriv_classes[2][4][10] = int_stack + 90;
 Libderiv->deriv_classes[2][4][9] = int_stack + 180;
 Libderiv->deriv_classes[2][4][8] = int_stack + 270;
 Libderiv->deriv_classes[2][4][7] = int_stack + 360;
 Libderiv->deriv_classes[2][4][6] = int_stack + 450;
 Libderiv->deriv_classes[2][4][2] = int_stack + 540;
 Libderiv->deriv_classes[2][4][1] = int_stack + 630;
 Libderiv->deriv_classes[2][4][0] = int_stack + 720;
 memset(int_stack,0,6480);

 Libderiv->dvrr_stack = int_stack + 810;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0g0(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[10] = int_stack + 90;
 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[9] = int_stack + 180;
 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[8] = int_stack + 270;
 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[7] = int_stack + 360;
 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[6] = int_stack + 450;
 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[2] = int_stack + 540;
 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[1] = int_stack + 630;
 /*--- compute (d0|g0) ---*/
     Libderiv->ABCD[0] = int_stack + 720;

}
