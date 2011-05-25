#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (p0|pp) integrals */

void d1vrr_order_p0pp(Libderiv_t *Libderiv, prim_data *Data)
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
 tmp = dvrr_stack + 6;
 target_ptr = Libderiv->dvrr_classes[1][1];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

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


 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+75, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+78, dvrr_stack+15, dvrr_stack+75, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+84, dvrr_stack+18, dvrr_stack+78, dvrr_stack+0, dvrr_stack+15, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+94, dvrr_stack+24, dvrr_stack+18, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+104, dvrr_stack+94, dvrr_stack+84, NULL, NULL, dvrr_stack+18);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+134,dvrr_stack+104,dvrr_stack+30,3);


 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+75, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+84, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+87, dvrr_stack+0, dvrr_stack+15, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+188, dvrr_stack+6, dvrr_stack+87, dvrr_stack+3, dvrr_stack+0, dvrr_stack+84);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+206, dvrr_stack+18, dvrr_stack+78, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+224, dvrr_stack+30, dvrr_stack+206, dvrr_stack+24, dvrr_stack+18, dvrr_stack+87);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+15, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 15;
 target_ptr = Libderiv->deriv_classes[1][1][11];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+206, dvrr_stack+134, NULL);
 tmp = dvrr_stack + 206;
 target_ptr = Libderiv->deriv_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+78, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 78;
 target_ptr = Libderiv->deriv_classes[1][1][10];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+260, dvrr_stack+134, NULL);
 tmp = dvrr_stack + 260;
 target_ptr = Libderiv->deriv_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+87, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 87;
 target_ptr = Libderiv->deriv_classes[1][1][9];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+48, dvrr_stack+134, NULL);
 tmp = dvrr_stack + 48;
 target_ptr = Libderiv->deriv_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+66, dvrr_stack+30, dvrr_stack+75);
 tmp = dvrr_stack + 66;
 target_ptr = Libderiv->deriv_classes[1][1][8];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+134, dvrr_stack+104, dvrr_stack+6);
 tmp = dvrr_stack + 134;
 target_ptr = Libderiv->deriv_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+152, dvrr_stack+30, dvrr_stack+75);
 tmp = dvrr_stack + 152;
 target_ptr = Libderiv->deriv_classes[1][1][7];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+161, dvrr_stack+104, dvrr_stack+6);
 tmp = dvrr_stack + 161;
 target_ptr = Libderiv->deriv_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+179, dvrr_stack+30, dvrr_stack+75);
 tmp = dvrr_stack + 179;
 target_ptr = Libderiv->deriv_classes[1][1][6];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+30, dvrr_stack+104, dvrr_stack+6);
 tmp = dvrr_stack + 30;
 target_ptr = Libderiv->deriv_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+6, dvrr_stack+188, dvrr_stack+3);
 tmp = dvrr_stack + 6;
 target_ptr = Libderiv->deriv_classes[1][1][2];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+96, dvrr_stack+224, dvrr_stack+24);
 tmp = dvrr_stack + 96;
 target_ptr = Libderiv->deriv_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+114, dvrr_stack+188, dvrr_stack+3);
 tmp = dvrr_stack + 114;
 target_ptr = Libderiv->deriv_classes[1][1][1];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+278, dvrr_stack+224, dvrr_stack+24);
 tmp = dvrr_stack + 278;
 target_ptr = Libderiv->deriv_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+123, dvrr_stack+188, dvrr_stack+3);
 tmp = dvrr_stack + 123;
 target_ptr = Libderiv->deriv_classes[1][1][0];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+188, dvrr_stack+224, dvrr_stack+24);
 tmp = dvrr_stack + 188;
 target_ptr = Libderiv->deriv_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];


}

