#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_00f0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (00|f0) integrals */

void d1hrr_order_00f0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[0][3][11] = int_stack + 0;
 Libderiv->deriv_classes[0][3][10] = int_stack + 10;
 Libderiv->deriv_classes[0][3][9] = int_stack + 20;
 Libderiv->deriv_classes[0][3][8] = int_stack + 30;
 Libderiv->deriv_classes[0][3][7] = int_stack + 40;
 Libderiv->deriv_classes[0][3][6] = int_stack + 50;
 Libderiv->deriv_classes[0][3][2] = int_stack + 60;
 Libderiv->deriv_classes[0][3][1] = int_stack + 70;
 Libderiv->deriv_classes[0][3][0] = int_stack + 80;
 memset(int_stack,0,720);

 Libderiv->dvrr_stack = int_stack + 90;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_00f0(Libderiv, Data);
   Data++;
 }

 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[10] = int_stack + 10;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[9] = int_stack + 20;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[8] = int_stack + 30;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[7] = int_stack + 40;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[6] = int_stack + 50;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[2] = int_stack + 60;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[1] = int_stack + 70;
 /*--- compute (00|f0) ---*/
     Libderiv->ABCD[0] = int_stack + 80;

}
