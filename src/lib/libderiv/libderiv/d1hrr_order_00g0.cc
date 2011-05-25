#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00g0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|g0) integrals */

void d1hrr_order_00g0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][4][11] = int_stack + 0;
 Libderiv->deriv_classes[0][4][10] = int_stack + 15;
 Libderiv->deriv_classes[0][4][9] = int_stack + 30;
 Libderiv->deriv_classes[0][4][8] = int_stack + 45;
 Libderiv->deriv_classes[0][4][7] = int_stack + 60;
 Libderiv->deriv_classes[0][4][6] = int_stack + 75;
 Libderiv->deriv_classes[0][4][2] = int_stack + 90;
 Libderiv->deriv_classes[0][4][1] = int_stack + 105;
 Libderiv->deriv_classes[0][4][0] = int_stack + 120;
 memset(int_stack,0,1080);

 Libderiv->dvrr_stack = int_stack + 135;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00g0(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[10] = int_stack + 15;
 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[9] = int_stack + 30;
 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[8] = int_stack + 45;
 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[7] = int_stack + 60;
 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[6] = int_stack + 75;
 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[2] = int_stack + 90;
 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[1] = int_stack + 105;
 /*--- compute (00|g0) ---*/
     Libderiv->ABCD[0] = int_stack + 120;

}
