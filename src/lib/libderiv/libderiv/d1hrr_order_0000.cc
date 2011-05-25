#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_0000(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|00) integrals */

void d1hrr_order_0000(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][0][11] = int_stack + 0;
 Libderiv->deriv_classes[0][0][10] = int_stack + 1;
 Libderiv->deriv_classes[0][0][9] = int_stack + 2;
 Libderiv->deriv_classes[0][0][8] = int_stack + 3;
 Libderiv->deriv_classes[0][0][7] = int_stack + 4;
 Libderiv->deriv_classes[0][0][6] = int_stack + 5;
 Libderiv->deriv_classes[0][0][2] = int_stack + 6;
 Libderiv->deriv_classes[0][0][1] = int_stack + 7;
 Libderiv->deriv_classes[0][0][0] = int_stack + 8;
 memset(int_stack,0,72);

 Libderiv->dvrr_stack = int_stack + 9;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_0000(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|00) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[10] = int_stack + 1;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[9] = int_stack + 2;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[8] = int_stack + 3;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[7] = int_stack + 4;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[6] = int_stack + 5;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[2] = int_stack + 6;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[1] = int_stack + 7;
 /*--- compute (00|00) ---*/
     Libderiv->ABCD[0] = int_stack + 8;

}
