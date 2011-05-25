#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (p0|p0) integrals */

void d1vrr_order_p0p0(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+1);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+15, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+18, dvrr_stack+0, dvrr_stack+15, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+24, dvrr_stack+3, dvrr_stack+0, Data->F+0, Data->F+1, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+30, dvrr_stack+24, dvrr_stack+18, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+48,dvrr_stack+30,dvrr_stack+6,3);


 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+18, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+21, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+75, dvrr_stack+0, dvrr_stack+15, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+84, dvrr_stack+6, dvrr_stack+75, dvrr_stack+3, dvrr_stack+0, dvrr_stack+21);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+21, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv_classes[1][1][11];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+75, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 75;
 target_ptr = Libderiv->deriv_classes[1][1][10];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+6, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 6;
 target_ptr = Libderiv->deriv_classes[1][1][9];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+48, dvrr_stack+30, dvrr_stack+18);
 tmp = dvrr_stack + 48;
 target_ptr = Libderiv->deriv_classes[1][1][8];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+57, dvrr_stack+30, dvrr_stack+18);
 tmp = dvrr_stack + 57;
 target_ptr = Libderiv->deriv_classes[1][1][7];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+66, dvrr_stack+30, dvrr_stack+18);
 tmp = dvrr_stack + 66;
 target_ptr = Libderiv->deriv_classes[1][1][6];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+30, dvrr_stack+84, dvrr_stack+3);
 tmp = dvrr_stack + 30;
 target_ptr = Libderiv->deriv_classes[1][1][2];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+39, dvrr_stack+84, dvrr_stack+3);
 tmp = dvrr_stack + 39;
 target_ptr = Libderiv->deriv_classes[1][1][1];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+102, dvrr_stack+84, dvrr_stack+3);
 tmp = dvrr_stack + 102;
 target_ptr = Libderiv->deriv_classes[1][1][0];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];


}

