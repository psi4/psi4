#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_f0f0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (f0|f0) integrals */

void d1hrr_order_f0f0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[3][3][11] = int_stack + 0;
 Libderiv->deriv_classes[3][3][10] = int_stack + 100;
 Libderiv->deriv_classes[3][3][9] = int_stack + 200;
 Libderiv->deriv_classes[3][3][8] = int_stack + 300;
 Libderiv->deriv_classes[3][3][7] = int_stack + 400;
 Libderiv->deriv_classes[3][3][6] = int_stack + 500;
 Libderiv->deriv_classes[3][3][2] = int_stack + 600;
 Libderiv->deriv_classes[3][3][1] = int_stack + 700;
 Libderiv->deriv_classes[3][3][0] = int_stack + 800;
 memset(int_stack,0,7200);

 Libderiv->dvrr_stack = int_stack + 900;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_f0f0(Libderiv, Data);
   Data++;
 }

 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[10] = int_stack + 100;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[9] = int_stack + 200;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[8] = int_stack + 300;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[7] = int_stack + 400;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[6] = int_stack + 500;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[2] = int_stack + 600;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[1] = int_stack + 700;
 /*--- compute (f0|f0) ---*/
     Libderiv->ABCD[0] = int_stack + 800;

}
