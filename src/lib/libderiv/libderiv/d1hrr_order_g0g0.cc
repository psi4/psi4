#include <stdio.h>
#include <string.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/hrr_header.h>

#include "d1hrr_header.h"

extern void d1vrr_order_g0g0(Libderiv_t *, prim_data *);

  /* Computes derivatives of (g0|g0) integrals */

void d1hrr_order_g0g0(Libderiv_t *Libderiv, int num_prim_comb)
{
 prim_data *Data = Libderiv->PrimQuartet;
 double *int_stack = Libderiv->int_stack;
 double *zero_stack = Libderiv->zero_stack;
 int i,j;
 double tmp, *target;

 Libderiv->deriv_classes[4][4][11] = int_stack + 0;
 Libderiv->deriv_classes[4][4][10] = int_stack + 225;
 Libderiv->deriv_classes[4][4][9] = int_stack + 450;
 Libderiv->deriv_classes[4][4][8] = int_stack + 675;
 Libderiv->deriv_classes[4][4][7] = int_stack + 900;
 Libderiv->deriv_classes[4][4][6] = int_stack + 1125;
 Libderiv->deriv_classes[4][4][2] = int_stack + 1350;
 Libderiv->deriv_classes[4][4][1] = int_stack + 1575;
 Libderiv->deriv_classes[4][4][0] = int_stack + 1800;
 memset(int_stack,0,16200);

 Libderiv->dvrr_stack = int_stack + 2025;
 for(i=0;i<num_prim_comb;i++) {
   d1vrr_order_g0g0(Libderiv, Data);
   Data++;
 }

 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[11] = int_stack + 0;
 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[10] = int_stack + 225;
 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[9] = int_stack + 450;
 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[8] = int_stack + 675;
 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[7] = int_stack + 900;
 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[6] = int_stack + 1125;
 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[2] = int_stack + 1350;
 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[1] = int_stack + 1575;
 /*--- compute (g0|g0) ---*/
     Libderiv->ABCD[0] = int_stack + 1800;

}
