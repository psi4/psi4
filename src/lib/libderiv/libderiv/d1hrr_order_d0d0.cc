#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_d0d0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (d0|d0) integrals */

void d1hrr_order_d0d0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[2][2][11] = int_stack + 0;
 Libderiv->deriv_classes[2][2][10] = int_stack + 36;
 Libderiv->deriv_classes[2][2][9] = int_stack + 72;
 Libderiv->deriv_classes[2][2][8] = int_stack + 108;
 Libderiv->deriv_classes[2][2][7] = int_stack + 144;
 Libderiv->deriv_classes[2][2][6] = int_stack + 180;
 Libderiv->deriv_classes[2][2][2] = int_stack + 216;
 Libderiv->deriv_classes[2][2][1] = int_stack + 252;
 Libderiv->deriv_classes[2][2][0] = int_stack + 288;
 memset(int_stack,0,2592);

 Libderiv->dvrr_stack = int_stack + 324;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_d0d0(Libderiv, Data);
   Data++;
 }

 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[10] = int_stack + 36;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[9] = int_stack + 72;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[8] = int_stack + 108;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[7] = int_stack + 144;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[6] = int_stack + 180;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[2] = int_stack + 216;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[1] = int_stack + 252;
 /*--- compute (d0|d0) ---*/
     Libderiv->ABCD[0] = int_stack + 288;

}
