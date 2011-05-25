#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|pp) integrals */

void d12vrr_order_00pp(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+0, dvrr_stack+3, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+0, Data->F+1, NULL, NULL, NULL);
 tmp = dvrr_stack + 12;
 target_ptr = Libderiv->dvrr_classes[0][1];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+12, dvrr_stack+0, Data->F+0, Data->F+1, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+6, NULL, NULL, dvrr_stack+0);

 /* compute (0 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+39,dvrr_stack+15,dvrr_stack+12,1);


 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+48, dvrr_stack+15, dvrr_stack+6, dvrr_stack+12, dvrr_stack+0, NULL);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+58,dvrr_stack+48,dvrr_stack+15,1);


 /* compute (0 0 | 1 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pd(Libderiv->CD,dvrr_stack+76,dvrr_stack+58,dvrr_stack+39,1);


 /* compute (0 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,3,dvrr_stack+94, dvrr_stack+76, dvrr_stack+12);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+103, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+106, dvrr_stack+3, dvrr_stack+103, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+112, dvrr_stack+6, dvrr_stack+106, dvrr_stack+0, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+122, dvrr_stack+48, dvrr_stack+112, dvrr_stack+15, dvrr_stack+6, NULL);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+137,dvrr_stack+122,dvrr_stack+48,1);


 /* compute (0 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+167,dvrr_stack+137,dvrr_stack+58,1);


 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,6,dvrr_stack+203, dvrr_stack+167, dvrr_stack+15);

 /* compute (0 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,3,dvrr_stack+221, dvrr_stack+76, dvrr_stack+12);

 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,6,dvrr_stack+230, dvrr_stack+167, dvrr_stack+15);

 /* compute (0 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,3,dvrr_stack+248, dvrr_stack+76, dvrr_stack+12);

 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,6,dvrr_stack+76, dvrr_stack+167, dvrr_stack+15);

 /* compute (0 0 | 0 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0p(Libderiv->CD,dvrr_stack+103,dvrr_stack+12,Data->F,1);


 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1,dvrr_stack+167, dvrr_stack+103, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+168, dvrr_stack+58, NULL);
 tmp = dvrr_stack + 168;
 target_ptr = Libderiv->deriv_classes[0][2][11];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+174, dvrr_stack+39, NULL);
 tmp = dvrr_stack + 174;
 target_ptr = Libderiv->deriv_classes[0][1][11];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+177, dvrr_stack+137, NULL);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1,dvrr_stack+187, dvrr_stack+103, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+188, dvrr_stack+58, NULL);
 tmp = dvrr_stack + 188;
 target_ptr = Libderiv->deriv_classes[0][2][10];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+194, dvrr_stack+39, NULL);
 tmp = dvrr_stack + 194;
 target_ptr = Libderiv->deriv_classes[0][1][10];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+257, dvrr_stack+137, NULL);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1,dvrr_stack+197, dvrr_stack+103, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+267, dvrr_stack+58, NULL);
 tmp = dvrr_stack + 267;
 target_ptr = Libderiv->deriv_classes[0][2][9];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+103, dvrr_stack+39, NULL);
 tmp = dvrr_stack + 103;
 target_ptr = Libderiv->deriv_classes[0][1][9];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+58, dvrr_stack+137, NULL);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_0(Data,1,1,dvrr_stack+137, dvrr_stack+12, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+138, dvrr_stack+48, dvrr_stack+12);
 tmp = dvrr_stack + 138;
 target_ptr = Libderiv->deriv_classes[0][2][8];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+144, dvrr_stack+15, Data->F+0);
 tmp = dvrr_stack + 144;
 target_ptr = Libderiv->deriv_classes[0][1][8];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+147, dvrr_stack+122, dvrr_stack+15);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+157, dvrr_stack+12, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+158, dvrr_stack+48, dvrr_stack+12);
 tmp = dvrr_stack + 158;
 target_ptr = Libderiv->deriv_classes[0][2][7];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+164, dvrr_stack+15, Data->F+0);
 tmp = dvrr_stack + 164;
 target_ptr = Libderiv->deriv_classes[0][1][7];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+273, dvrr_stack+122, dvrr_stack+15);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+39, dvrr_stack+12, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+40, dvrr_stack+48, dvrr_stack+12);
 tmp = dvrr_stack + 40;
 target_ptr = Libderiv->deriv_classes[0][2][6];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+68, dvrr_stack+15, Data->F+0);
 tmp = dvrr_stack + 68;
 target_ptr = Libderiv->deriv_classes[0][1][6];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+283, dvrr_stack+122, dvrr_stack+15);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+122, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+1);

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+293,dvrr_stack+21,dvrr_stack+122,3);


 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+320, dvrr_stack+293, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+329, dvrr_stack+48, dvrr_stack+112, NULL, NULL, dvrr_stack+6);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+359,dvrr_stack+329,dvrr_stack+21,3);


 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+413, dvrr_stack+359, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+112, dvrr_stack+293, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+431, dvrr_stack+359, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+46, dvrr_stack+293, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+293, dvrr_stack+359, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+55, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+311, dvrr_stack+21, dvrr_stack+55);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+359, dvrr_stack+329, dvrr_stack+122);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+377, dvrr_stack+21, dvrr_stack+55);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+386, dvrr_stack+329, dvrr_stack+122);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+404, dvrr_stack+21, dvrr_stack+55);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+449, dvrr_stack+329, dvrr_stack+122);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+55, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+329, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+338, dvrr_stack+122, dvrr_stack+329, dvrr_stack+12, dvrr_stack+0, dvrr_stack+55);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+467, dvrr_stack+338, dvrr_stack+12);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+476, dvrr_stack+6, dvrr_stack+106, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+494, dvrr_stack+21, dvrr_stack+476, dvrr_stack+15, dvrr_stack+6, dvrr_stack+329);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+476, dvrr_stack+494, dvrr_stack+15);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+329, dvrr_stack+338, dvrr_stack+12);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+530, dvrr_stack+494, dvrr_stack+15);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+0, dvrr_stack+338, dvrr_stack+12);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+338, dvrr_stack+494, dvrr_stack+15);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+106, dvrr_stack+21, NULL);
 tmp = dvrr_stack + 106;
 target_ptr = Libderiv->deriv_classes[0][2][2];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+131, dvrr_stack+21, NULL);
 tmp = dvrr_stack + 131;
 target_ptr = Libderiv->deriv_classes[0][2][1];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+494, dvrr_stack+21, NULL);
 tmp = dvrr_stack + 494;
 target_ptr = Libderiv->deriv_classes[0][2][0];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,3,dvrr_stack+356, dvrr_stack+94, NULL);
 tmp = dvrr_stack + 356;
 target_ptr = Libderiv->deriv2_classes[0][1][143];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,6,dvrr_stack+500, dvrr_stack+203, NULL);
 tmp = dvrr_stack + 500;
 target_ptr = Libderiv->deriv2_classes[0][2][143];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,3,dvrr_stack+55, dvrr_stack+94, NULL);
 tmp = dvrr_stack + 55;
 target_ptr = Libderiv->deriv2_classes[0][1][131];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,6,dvrr_stack+506, dvrr_stack+203, NULL);
 tmp = dvrr_stack + 506;
 target_ptr = Libderiv->deriv2_classes[0][2][131];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+512, dvrr_stack+221, NULL);
 tmp = dvrr_stack + 512;
 target_ptr = Libderiv->deriv2_classes[0][1][130];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+515, dvrr_stack+230, NULL);
 tmp = dvrr_stack + 515;
 target_ptr = Libderiv->deriv2_classes[0][2][130];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,3,dvrr_stack+521, dvrr_stack+94, NULL);
 tmp = dvrr_stack + 521;
 target_ptr = Libderiv->deriv2_classes[0][1][119];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,6,dvrr_stack+524, dvrr_stack+203, NULL);
 tmp = dvrr_stack + 524;
 target_ptr = Libderiv->deriv2_classes[0][2][119];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+94, dvrr_stack+221, NULL);
 tmp = dvrr_stack + 94;
 target_ptr = Libderiv->deriv2_classes[0][1][118];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+97, dvrr_stack+230, NULL);
 tmp = dvrr_stack + 97;
 target_ptr = Libderiv->deriv2_classes[0][2][118];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+9, dvrr_stack+248, NULL);
 tmp = dvrr_stack + 9;
 target_ptr = Libderiv->deriv2_classes[0][1][117];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+12, dvrr_stack+76, NULL);
 tmp = dvrr_stack + 12;
 target_ptr = Libderiv->deriv2_classes[0][2][117];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+18, dvrr_stack+168, dvrr_stack+167);
 tmp = dvrr_stack + 18;
 target_ptr = Libderiv->deriv2_classes[0][1][107];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+21, dvrr_stack+177, dvrr_stack+174);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv2_classes[0][2][107];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+27, dvrr_stack+188, dvrr_stack+187);
 tmp = dvrr_stack + 27;
 target_ptr = Libderiv->deriv2_classes[0][1][106];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+30, dvrr_stack+257, dvrr_stack+194);
 tmp = dvrr_stack + 30;
 target_ptr = Libderiv->deriv2_classes[0][2][106];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+36, dvrr_stack+267, dvrr_stack+197);
 tmp = dvrr_stack + 36;
 target_ptr = Libderiv->deriv2_classes[0][1][105];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+71, dvrr_stack+58, dvrr_stack+103);
 tmp = dvrr_stack + 71;
 target_ptr = Libderiv->deriv2_classes[0][2][105];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+77, dvrr_stack+138, dvrr_stack+137);
 tmp = dvrr_stack + 77;
 target_ptr = Libderiv->deriv2_classes[0][1][104];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+80, dvrr_stack+147, dvrr_stack+144);
 tmp = dvrr_stack + 80;
 target_ptr = Libderiv->deriv2_classes[0][2][104];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+86, dvrr_stack+168, dvrr_stack+167);
 tmp = dvrr_stack + 86;
 target_ptr = Libderiv->deriv2_classes[0][1][95];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+198, dvrr_stack+177, dvrr_stack+174);
 tmp = dvrr_stack + 198;
 target_ptr = Libderiv->deriv2_classes[0][2][95];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+89, dvrr_stack+188, dvrr_stack+187);
 tmp = dvrr_stack + 89;
 target_ptr = Libderiv->deriv2_classes[0][1][94];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+204, dvrr_stack+257, dvrr_stack+194);
 tmp = dvrr_stack + 204;
 target_ptr = Libderiv->deriv2_classes[0][2][94];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+210, dvrr_stack+267, dvrr_stack+197);
 tmp = dvrr_stack + 210;
 target_ptr = Libderiv->deriv2_classes[0][1][93];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+213, dvrr_stack+58, dvrr_stack+103);
 tmp = dvrr_stack + 213;
 target_ptr = Libderiv->deriv2_classes[0][2][93];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+219, dvrr_stack+138, dvrr_stack+137);
 tmp = dvrr_stack + 219;
 target_ptr = Libderiv->deriv2_classes[0][1][92];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+222, dvrr_stack+147, dvrr_stack+144);
 tmp = dvrr_stack + 222;
 target_ptr = Libderiv->deriv2_classes[0][2][92];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+228, dvrr_stack+158, dvrr_stack+157);
 tmp = dvrr_stack + 228;
 target_ptr = Libderiv->deriv2_classes[0][1][91];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+231, dvrr_stack+273, dvrr_stack+164);
 tmp = dvrr_stack + 231;
 target_ptr = Libderiv->deriv2_classes[0][2][91];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+237, dvrr_stack+168, dvrr_stack+167);
 tmp = dvrr_stack + 237;
 target_ptr = Libderiv->deriv2_classes[0][1][83];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+167, dvrr_stack+177, dvrr_stack+174);
 tmp = dvrr_stack + 167;
 target_ptr = Libderiv->deriv2_classes[0][2][83];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+173, dvrr_stack+188, dvrr_stack+187);
 tmp = dvrr_stack + 173;
 target_ptr = Libderiv->deriv2_classes[0][1][82];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+176, dvrr_stack+257, dvrr_stack+194);
 tmp = dvrr_stack + 176;
 target_ptr = Libderiv->deriv2_classes[0][2][82];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+182, dvrr_stack+267, dvrr_stack+197);
 tmp = dvrr_stack + 182;
 target_ptr = Libderiv->deriv2_classes[0][1][81];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+185, dvrr_stack+58, dvrr_stack+103);
 tmp = dvrr_stack + 185;
 target_ptr = Libderiv->deriv2_classes[0][2][81];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+103, dvrr_stack+138, dvrr_stack+137);
 tmp = dvrr_stack + 103;
 target_ptr = Libderiv->deriv2_classes[0][1][80];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+137, dvrr_stack+147, dvrr_stack+144);
 tmp = dvrr_stack + 137;
 target_ptr = Libderiv->deriv2_classes[0][2][80];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+143, dvrr_stack+158, dvrr_stack+157);
 tmp = dvrr_stack + 143;
 target_ptr = Libderiv->deriv2_classes[0][1][79];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+146, dvrr_stack+273, dvrr_stack+164);
 tmp = dvrr_stack + 146;
 target_ptr = Libderiv->deriv2_classes[0][2][79];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+152, dvrr_stack+40, dvrr_stack+39);
 tmp = dvrr_stack + 152;
 target_ptr = Libderiv->deriv2_classes[0][1][78];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+39, dvrr_stack+283, dvrr_stack+68);
 tmp = dvrr_stack + 39;
 target_ptr = Libderiv->deriv2_classes[0][2][78];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,3,dvrr_stack+155, dvrr_stack+320, NULL);
 tmp = dvrr_stack + 155;
 target_ptr = Libderiv->deriv2_classes[0][1][35];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,6,dvrr_stack+158, dvrr_stack+413, NULL);
 tmp = dvrr_stack + 158;
 target_ptr = Libderiv->deriv2_classes[0][2][35];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+164, dvrr_stack+112, NULL);
 tmp = dvrr_stack + 164;
 target_ptr = Libderiv->deriv2_classes[0][1][34];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+58, dvrr_stack+431, NULL);
 tmp = dvrr_stack + 58;
 target_ptr = Libderiv->deriv2_classes[0][2][34];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+64, dvrr_stack+46, NULL);
 tmp = dvrr_stack + 64;
 target_ptr = Libderiv->deriv2_classes[0][1][33];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+191, dvrr_stack+293, NULL);
 tmp = dvrr_stack + 191;
 target_ptr = Libderiv->deriv2_classes[0][2][33];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+67, dvrr_stack+311, NULL);
 tmp = dvrr_stack + 67;
 target_ptr = Libderiv->deriv2_classes[0][1][32];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+240, dvrr_stack+359, NULL);
 tmp = dvrr_stack + 240;
 target_ptr = Libderiv->deriv2_classes[0][2][32];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+246, dvrr_stack+377, NULL);
 tmp = dvrr_stack + 246;
 target_ptr = Libderiv->deriv2_classes[0][1][31];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+249, dvrr_stack+386, NULL);
 tmp = dvrr_stack + 249;
 target_ptr = Libderiv->deriv2_classes[0][2][31];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+255, dvrr_stack+122, NULL);
 tmp = dvrr_stack + 255;
 target_ptr = Libderiv->deriv_classes[0][1][2];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+258, dvrr_stack+404, NULL);
 tmp = dvrr_stack + 258;
 target_ptr = Libderiv->deriv2_classes[0][1][30];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+261, dvrr_stack+449, NULL);
 tmp = dvrr_stack + 261;
 target_ptr = Libderiv->deriv2_classes[0][2][30];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+267, dvrr_stack+467, NULL);
 tmp = dvrr_stack + 267;
 target_ptr = Libderiv->deriv2_classes[0][1][26];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+270, dvrr_stack+476, NULL);
 tmp = dvrr_stack + 270;
 target_ptr = Libderiv->deriv2_classes[0][2][26];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,3,dvrr_stack+276, dvrr_stack+320, NULL);
 tmp = dvrr_stack + 276;
 target_ptr = Libderiv->deriv2_classes[0][1][23];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,6,dvrr_stack+279, dvrr_stack+413, NULL);
 tmp = dvrr_stack + 279;
 target_ptr = Libderiv->deriv2_classes[0][2][23];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+285, dvrr_stack+112, NULL);
 tmp = dvrr_stack + 285;
 target_ptr = Libderiv->deriv2_classes[0][1][22];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+548, dvrr_stack+431, NULL);
 tmp = dvrr_stack + 548;
 target_ptr = Libderiv->deriv2_classes[0][2][22];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+288, dvrr_stack+46, NULL);
 tmp = dvrr_stack + 288;
 target_ptr = Libderiv->deriv2_classes[0][1][21];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+554, dvrr_stack+293, NULL);
 tmp = dvrr_stack + 554;
 target_ptr = Libderiv->deriv2_classes[0][2][21];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+560, dvrr_stack+311, NULL);
 tmp = dvrr_stack + 560;
 target_ptr = Libderiv->deriv2_classes[0][1][20];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+563, dvrr_stack+359, NULL);
 tmp = dvrr_stack + 563;
 target_ptr = Libderiv->deriv2_classes[0][2][20];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+569, dvrr_stack+377, NULL);
 tmp = dvrr_stack + 569;
 target_ptr = Libderiv->deriv2_classes[0][1][19];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+572, dvrr_stack+386, NULL);
 tmp = dvrr_stack + 572;
 target_ptr = Libderiv->deriv2_classes[0][2][19];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+578, dvrr_stack+122, NULL);
 tmp = dvrr_stack + 578;
 target_ptr = Libderiv->deriv_classes[0][1][1];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+581, dvrr_stack+404, NULL);
 tmp = dvrr_stack + 581;
 target_ptr = Libderiv->deriv2_classes[0][1][18];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+584, dvrr_stack+449, NULL);
 tmp = dvrr_stack + 584;
 target_ptr = Libderiv->deriv2_classes[0][2][18];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+590, dvrr_stack+467, NULL);
 tmp = dvrr_stack + 590;
 target_ptr = Libderiv->deriv2_classes[0][1][14];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+593, dvrr_stack+476, NULL);
 tmp = dvrr_stack + 593;
 target_ptr = Libderiv->deriv2_classes[0][2][14];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+599, dvrr_stack+329, NULL);
 tmp = dvrr_stack + 599;
 target_ptr = Libderiv->deriv2_classes[0][1][13];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+602, dvrr_stack+530, NULL);
 tmp = dvrr_stack + 602;
 target_ptr = Libderiv->deriv2_classes[0][2][13];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,3,dvrr_stack+608, dvrr_stack+320, NULL);
 tmp = dvrr_stack + 608;
 target_ptr = Libderiv->deriv2_classes[0][1][11];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,6,dvrr_stack+320, dvrr_stack+413, NULL);
 tmp = dvrr_stack + 320;
 target_ptr = Libderiv->deriv2_classes[0][2][11];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+326, dvrr_stack+112, NULL);
 tmp = dvrr_stack + 326;
 target_ptr = Libderiv->deriv2_classes[0][1][10];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+112, dvrr_stack+431, NULL);
 tmp = dvrr_stack + 112;
 target_ptr = Libderiv->deriv2_classes[0][2][10];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+118, dvrr_stack+46, NULL);
 tmp = dvrr_stack + 118;
 target_ptr = Libderiv->deriv2_classes[0][1][9];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+413, dvrr_stack+293, NULL);
 tmp = dvrr_stack + 413;
 target_ptr = Libderiv->deriv2_classes[0][2][9];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+419, dvrr_stack+311, NULL);
 tmp = dvrr_stack + 419;
 target_ptr = Libderiv->deriv2_classes[0][1][8];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+422, dvrr_stack+359, NULL);
 tmp = dvrr_stack + 422;
 target_ptr = Libderiv->deriv2_classes[0][2][8];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+359, dvrr_stack+377, NULL);
 tmp = dvrr_stack + 359;
 target_ptr = Libderiv->deriv2_classes[0][1][7];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+362, dvrr_stack+386, NULL);
 tmp = dvrr_stack + 362;
 target_ptr = Libderiv->deriv2_classes[0][2][7];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+368, dvrr_stack+122, NULL);
 tmp = dvrr_stack + 368;
 target_ptr = Libderiv->deriv_classes[0][1][0];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+371, dvrr_stack+404, NULL);
 tmp = dvrr_stack + 371;
 target_ptr = Libderiv->deriv2_classes[0][1][6];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+374, dvrr_stack+449, NULL);
 tmp = dvrr_stack + 374;
 target_ptr = Libderiv->deriv2_classes[0][2][6];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+380, dvrr_stack+467, NULL);
 tmp = dvrr_stack + 380;
 target_ptr = Libderiv->deriv2_classes[0][1][2];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+383, dvrr_stack+476, NULL);
 tmp = dvrr_stack + 383;
 target_ptr = Libderiv->deriv2_classes[0][2][2];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+389, dvrr_stack+329, NULL);
 tmp = dvrr_stack + 389;
 target_ptr = Libderiv->deriv2_classes[0][1][1];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+329, dvrr_stack+530, NULL);
 tmp = dvrr_stack + 329;
 target_ptr = Libderiv->deriv2_classes[0][2][1];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+335, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 335;
 target_ptr = Libderiv->deriv2_classes[0][1][0];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+0, dvrr_stack+338, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[0][2][0];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];


}

