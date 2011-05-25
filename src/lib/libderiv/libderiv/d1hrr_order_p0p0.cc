#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_p0p0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (p0|p0) integrals */

void d1hrr_order_p0p0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[1][1][11] = int_stack + 0;
 Libderiv->deriv_classes[1][1][10] = int_stack + 9;
 Libderiv->deriv_classes[1][1][9] = int_stack + 18;
 Libderiv->deriv_classes[1][1][8] = int_stack + 27;
 Libderiv->deriv_classes[1][1][7] = int_stack + 36;
 Libderiv->deriv_classes[1][1][6] = int_stack + 45;
 Libderiv->deriv_classes[1][1][2] = int_stack + 54;
 Libderiv->deriv_classes[1][1][1] = int_stack + 63;
 Libderiv->deriv_classes[1][1][0] = int_stack + 72;
 memset(int_stack,0,648);

 Libderiv->dvrr_stack = int_stack + 81;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_p0p0(Libderiv, Data);
   Data++;
 }

 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[10] = int_stack + 9;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[9] = int_stack + 18;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[8] = int_stack + 27;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[7] = int_stack + 36;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[6] = int_stack + 45;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[2] = int_stack + 54;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[1] = int_stack + 63;
 /*--- compute (p0|p0) ---*/
     Libderiv->ABCD[0] = int_stack + 72;

}
