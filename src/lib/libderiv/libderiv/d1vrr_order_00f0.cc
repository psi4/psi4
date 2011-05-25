#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|f0) integrals */

void d1vrr_order_00f0(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+9, dvrr_stack+0, dvrr_stack+6, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+9, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+31, dvrr_stack+6, dvrr_stack+3, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+37, dvrr_stack+9, dvrr_stack+31, dvrr_stack+0, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+47, dvrr_stack+21, dvrr_stack+37, dvrr_stack+15, dvrr_stack+9, NULL);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+62,dvrr_stack+47,dvrr_stack+21,1);


 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+92, dvrr_stack+21, dvrr_stack+37, NULL, NULL, dvrr_stack+9);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+21, dvrr_stack+62, NULL);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+31, dvrr_stack+62, NULL);
 tmp = dvrr_stack + 31;
 target_ptr = Libderiv->deriv_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+0, dvrr_stack+62, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+62, dvrr_stack+47, dvrr_stack+15);
 tmp = dvrr_stack + 62;
 target_ptr = Libderiv->deriv_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+72, dvrr_stack+47, dvrr_stack+15);
 tmp = dvrr_stack + 72;
 target_ptr = Libderiv->deriv_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+82, dvrr_stack+47, dvrr_stack+15);
 tmp = dvrr_stack + 82;
 target_ptr = Libderiv->deriv_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+41, dvrr_stack+92, NULL);
 tmp = dvrr_stack + 41;
 target_ptr = Libderiv->deriv_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+51, dvrr_stack+92, NULL);
 tmp = dvrr_stack + 51;
 target_ptr = Libderiv->deriv_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+10, dvrr_stack+92, NULL);
 tmp = dvrr_stack + 10;
 target_ptr = Libderiv->deriv_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];


}

