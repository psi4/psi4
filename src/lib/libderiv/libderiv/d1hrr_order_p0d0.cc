#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0d0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|d0) integrals */

void d1hrr_order_p0d0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][2][11] = int_stack + 0;
 Libderiv->deriv_classes[1][2][10] = int_stack + 18;
 Libderiv->deriv_classes[1][2][9] = int_stack + 36;
 Libderiv->deriv_classes[1][2][8] = int_stack + 54;
 Libderiv->deriv_classes[1][2][7] = int_stack + 72;
 Libderiv->deriv_classes[1][2][6] = int_stack + 90;
 Libderiv->deriv_classes[1][2][2] = int_stack + 108;
 Libderiv->deriv_classes[1][2][1] = int_stack + 126;
 Libderiv->deriv_classes[1][2][0] = int_stack + 144;
 memset(int_stack,0,1296);

 Libderiv->dvrr_stack = int_stack + 162;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0d0(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[10] = int_stack + 18;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[9] = int_stack + 36;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[8] = int_stack + 54;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[7] = int_stack + 72;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[6] = int_stack + 90;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[2] = int_stack + 108;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[1] = int_stack + 126;
 /*--- compute (p0|d0) ---*/
     Libderiv->ABCD[0] = int_stack + 144;

}
