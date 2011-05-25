#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|00) integrals */

void d1vrr_order_0000(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 0 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0p(Libderiv->CD,dvrr_stack+3,dvrr_stack+0,Data->F,1);


 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+6, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1,dvrr_stack+9, dvrr_stack+3, NULL);
 tmp = dvrr_stack + 9;
 target_ptr = Libderiv->deriv_classes[0][0][11];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1,dvrr_stack+10, dvrr_stack+3, NULL);
 tmp = dvrr_stack + 10;
 target_ptr = Libderiv->deriv_classes[0][0][10];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1,dvrr_stack+11, dvrr_stack+3, NULL);
 tmp = dvrr_stack + 11;
 target_ptr = Libderiv->deriv_classes[0][0][9];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_0(Data,1,1,dvrr_stack+3, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 3;
 target_ptr = Libderiv->deriv_classes[0][0][8];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+4, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 4;
 target_ptr = Libderiv->deriv_classes[0][0][7];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+5, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 5;
 target_ptr = Libderiv->deriv_classes[0][0][6];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,1,dvrr_stack+0, dvrr_stack+6, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[0][0][2];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+1, dvrr_stack+6, NULL);
 tmp = dvrr_stack + 1;
 target_ptr = Libderiv->deriv_classes[0][0][1];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+2, dvrr_stack+6, NULL);
 tmp = dvrr_stack + 2;
 target_ptr = Libderiv->deriv_classes[0][0][0];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];


}

