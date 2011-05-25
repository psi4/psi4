#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|dp) integrals */

void d12vrr_order_00dp(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+0, dvrr_stack+12, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+6, dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+34, dvrr_stack+31, dvrr_stack+3, Data->F+0, Data->F+1, NULL);
 tmp = dvrr_stack + 34;
 target_ptr = Libderiv->dvrr_classes[0][2];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+34, dvrr_stack+6, dvrr_stack+31, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+50, dvrr_stack+40, dvrr_stack+21, NULL, NULL, dvrr_stack+6);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+80,dvrr_stack+40,dvrr_stack+34,1);


 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+98, dvrr_stack+40, dvrr_stack+21, dvrr_stack+34, dvrr_stack+6, NULL);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+113,dvrr_stack+98,dvrr_stack+40,1);


 /* compute (0 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+143,dvrr_stack+113,dvrr_stack+80,1);


 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,6,dvrr_stack+179, dvrr_stack+143, dvrr_stack+34);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+197, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+200, dvrr_stack+12, dvrr_stack+197, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+206, dvrr_stack+15, dvrr_stack+200, dvrr_stack+0, dvrr_stack+12, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+216, dvrr_stack+21, dvrr_stack+206, dvrr_stack+6, dvrr_stack+15, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+231, dvrr_stack+98, dvrr_stack+216, dvrr_stack+40, dvrr_stack+21, NULL);

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+252,dvrr_stack+231,dvrr_stack+98,1);


 /* compute (0 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+297,dvrr_stack+252,dvrr_stack+113,1);


 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,10,dvrr_stack+357, dvrr_stack+297, dvrr_stack+40);

 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,6,dvrr_stack+387, dvrr_stack+143, dvrr_stack+34);

 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,10,dvrr_stack+405, dvrr_stack+297, dvrr_stack+40);

 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,6,dvrr_stack+435, dvrr_stack+143, dvrr_stack+34);

 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,10,dvrr_stack+143, dvrr_stack+297, dvrr_stack+40);

 /* compute (0 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+197,dvrr_stack+34,dvrr_stack+31,1);


 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+12, dvrr_stack+197, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+297, dvrr_stack+113, NULL);
 tmp = dvrr_stack + 297;
 target_ptr = Libderiv->deriv_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+173, dvrr_stack+80, NULL);
 tmp = dvrr_stack + 173;
 target_ptr = Libderiv->deriv_classes[0][2][11];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+307, dvrr_stack+252, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+322, dvrr_stack+197, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+325, dvrr_stack+113, NULL);
 tmp = dvrr_stack + 325;
 target_ptr = Libderiv->deriv_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+335, dvrr_stack+80, NULL);
 tmp = dvrr_stack + 335;
 target_ptr = Libderiv->deriv_classes[0][2][10];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+341, dvrr_stack+252, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+453, dvrr_stack+197, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+456, dvrr_stack+113, NULL);
 tmp = dvrr_stack + 456;
 target_ptr = Libderiv->deriv_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+113, dvrr_stack+80, NULL);
 tmp = dvrr_stack + 113;
 target_ptr = Libderiv->deriv_classes[0][2][9];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+80, dvrr_stack+252, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+95, dvrr_stack+34, Data->F+0);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+252, dvrr_stack+98, dvrr_stack+34);
 tmp = dvrr_stack + 252;
 target_ptr = Libderiv->deriv_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+262, dvrr_stack+40, dvrr_stack+31);
 tmp = dvrr_stack + 262;
 target_ptr = Libderiv->deriv_classes[0][2][8];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+268, dvrr_stack+231, dvrr_stack+40);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+283, dvrr_stack+34, Data->F+0);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+286, dvrr_stack+98, dvrr_stack+34);
 tmp = dvrr_stack + 286;
 target_ptr = Libderiv->deriv_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+119, dvrr_stack+40, dvrr_stack+31);
 tmp = dvrr_stack + 119;
 target_ptr = Libderiv->deriv_classes[0][2][7];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+125, dvrr_stack+231, dvrr_stack+40);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+140, dvrr_stack+34, Data->F+0);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+466, dvrr_stack+98, dvrr_stack+34);
 tmp = dvrr_stack + 466;
 target_ptr = Libderiv->deriv_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+197, dvrr_stack+40, dvrr_stack+31);
 tmp = dvrr_stack + 197;
 target_ptr = Libderiv->deriv_classes[0][2][6];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+476, dvrr_stack+231, dvrr_stack+40);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+231, dvrr_stack+34, dvrr_stack+6, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+491,dvrr_stack+50,dvrr_stack+231,3);


 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+545, dvrr_stack+491, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+563, dvrr_stack+98, dvrr_stack+216, NULL, NULL, dvrr_stack+21);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+608,dvrr_stack+563,dvrr_stack+50,3);


 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+698, dvrr_stack+608, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+728, dvrr_stack+491, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+746, dvrr_stack+608, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+776, dvrr_stack+491, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+491, dvrr_stack+608, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+608, dvrr_stack+31, dvrr_stack+3, NULL, NULL, Data->F+1);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+617, dvrr_stack+50, dvrr_stack+608);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+635, dvrr_stack+563, dvrr_stack+231);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+665, dvrr_stack+50, dvrr_stack+608);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+794, dvrr_stack+563, dvrr_stack+231);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+521, dvrr_stack+50, dvrr_stack+608);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+824, dvrr_stack+563, dvrr_stack+231);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+563, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+572, dvrr_stack+6, dvrr_stack+15, NULL, NULL, dvrr_stack+0);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+854, dvrr_stack+231, dvrr_stack+572, dvrr_stack+34, dvrr_stack+6, dvrr_stack+563);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+590, dvrr_stack+854, dvrr_stack+34);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+890, dvrr_stack+21, dvrr_stack+206, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+920, dvrr_stack+50, dvrr_stack+890, dvrr_stack+40, dvrr_stack+21, dvrr_stack+572);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+890, dvrr_stack+920, dvrr_stack+40);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+15, dvrr_stack+854, dvrr_stack+34);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+980, dvrr_stack+920, dvrr_stack+40);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+563, dvrr_stack+854, dvrr_stack+34);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+854, dvrr_stack+920, dvrr_stack+40);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+920, dvrr_stack+50, NULL);
 tmp = dvrr_stack + 920;
 target_ptr = Libderiv->deriv_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+930, dvrr_stack+50, NULL);
 tmp = dvrr_stack + 930;
 target_ptr = Libderiv->deriv_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+940, dvrr_stack+50, NULL);
 tmp = dvrr_stack + 940;
 target_ptr = Libderiv->deriv_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,6,dvrr_stack+884, dvrr_stack+179, NULL);
 tmp = dvrr_stack + 884;
 target_ptr = Libderiv->deriv2_classes[0][2][143];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,10,dvrr_stack+950, dvrr_stack+357, NULL);
 tmp = dvrr_stack + 950;
 target_ptr = Libderiv->deriv2_classes[0][3][143];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,6,dvrr_stack+539, dvrr_stack+179, NULL);
 tmp = dvrr_stack + 539;
 target_ptr = Libderiv->deriv2_classes[0][2][131];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,10,dvrr_stack+960, dvrr_stack+357, NULL);
 tmp = dvrr_stack + 960;
 target_ptr = Libderiv->deriv2_classes[0][3][131];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+970, dvrr_stack+387, NULL);
 tmp = dvrr_stack + 970;
 target_ptr = Libderiv->deriv2_classes[0][2][130];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+33, dvrr_stack+405, NULL);
 tmp = dvrr_stack + 33;
 target_ptr = Libderiv->deriv2_classes[0][3][130];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,6,dvrr_stack+43, dvrr_stack+179, NULL);
 tmp = dvrr_stack + 43;
 target_ptr = Libderiv->deriv2_classes[0][2][119];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,10,dvrr_stack+179, dvrr_stack+357, NULL);
 tmp = dvrr_stack + 179;
 target_ptr = Libderiv->deriv2_classes[0][3][119];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+189, dvrr_stack+387, NULL);
 tmp = dvrr_stack + 189;
 target_ptr = Libderiv->deriv2_classes[0][2][118];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+49, dvrr_stack+405, NULL);
 tmp = dvrr_stack + 49;
 target_ptr = Libderiv->deriv2_classes[0][3][118];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+59, dvrr_stack+435, NULL);
 tmp = dvrr_stack + 59;
 target_ptr = Libderiv->deriv2_classes[0][2][117];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+65, dvrr_stack+143, NULL);
 tmp = dvrr_stack + 65;
 target_ptr = Libderiv->deriv2_classes[0][3][117];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+143, dvrr_stack+297, dvrr_stack+12);
 tmp = dvrr_stack + 143;
 target_ptr = Libderiv->deriv2_classes[0][2][107];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+149, dvrr_stack+307, dvrr_stack+173);
 tmp = dvrr_stack + 149;
 target_ptr = Libderiv->deriv2_classes[0][3][107];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+159, dvrr_stack+325, dvrr_stack+322);
 tmp = dvrr_stack + 159;
 target_ptr = Libderiv->deriv2_classes[0][2][106];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+0, dvrr_stack+341, dvrr_stack+335);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[0][3][106];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+165, dvrr_stack+456, dvrr_stack+453);
 tmp = dvrr_stack + 165;
 target_ptr = Libderiv->deriv2_classes[0][2][105];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+683, dvrr_stack+80, dvrr_stack+113);
 tmp = dvrr_stack + 683;
 target_ptr = Libderiv->deriv2_classes[0][3][105];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+581, dvrr_stack+252, dvrr_stack+95);
 tmp = dvrr_stack + 581;
 target_ptr = Libderiv->deriv2_classes[0][2][104];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+98, dvrr_stack+268, dvrr_stack+262);
 tmp = dvrr_stack + 98;
 target_ptr = Libderiv->deriv2_classes[0][3][104];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+608, dvrr_stack+297, dvrr_stack+12);
 tmp = dvrr_stack + 608;
 target_ptr = Libderiv->deriv2_classes[0][2][95];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+203, dvrr_stack+307, dvrr_stack+173);
 tmp = dvrr_stack + 203;
 target_ptr = Libderiv->deriv2_classes[0][3][95];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+213, dvrr_stack+325, dvrr_stack+322);
 tmp = dvrr_stack + 213;
 target_ptr = Libderiv->deriv2_classes[0][2][94];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+219, dvrr_stack+341, dvrr_stack+335);
 tmp = dvrr_stack + 219;
 target_ptr = Libderiv->deriv2_classes[0][3][94];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+356, dvrr_stack+456, dvrr_stack+453);
 tmp = dvrr_stack + 356;
 target_ptr = Libderiv->deriv2_classes[0][2][93];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+362, dvrr_stack+80, dvrr_stack+113);
 tmp = dvrr_stack + 362;
 target_ptr = Libderiv->deriv2_classes[0][3][93];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+372, dvrr_stack+252, dvrr_stack+95);
 tmp = dvrr_stack + 372;
 target_ptr = Libderiv->deriv2_classes[0][2][92];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+378, dvrr_stack+268, dvrr_stack+262);
 tmp = dvrr_stack + 378;
 target_ptr = Libderiv->deriv2_classes[0][3][92];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+388, dvrr_stack+286, dvrr_stack+283);
 tmp = dvrr_stack + 388;
 target_ptr = Libderiv->deriv2_classes[0][2][91];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+394, dvrr_stack+125, dvrr_stack+119);
 tmp = dvrr_stack + 394;
 target_ptr = Libderiv->deriv2_classes[0][3][91];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+404, dvrr_stack+297, dvrr_stack+12);
 tmp = dvrr_stack + 404;
 target_ptr = Libderiv->deriv2_classes[0][2][83];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+296, dvrr_stack+307, dvrr_stack+173);
 tmp = dvrr_stack + 296;
 target_ptr = Libderiv->deriv2_classes[0][3][83];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+171, dvrr_stack+325, dvrr_stack+322);
 tmp = dvrr_stack + 171;
 target_ptr = Libderiv->deriv2_classes[0][2][82];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+306, dvrr_stack+341, dvrr_stack+335);
 tmp = dvrr_stack + 306;
 target_ptr = Libderiv->deriv2_classes[0][3][82];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+316, dvrr_stack+456, dvrr_stack+453);
 tmp = dvrr_stack + 316;
 target_ptr = Libderiv->deriv2_classes[0][2][81];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+322, dvrr_stack+80, dvrr_stack+113);
 tmp = dvrr_stack + 322;
 target_ptr = Libderiv->deriv2_classes[0][3][81];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+75, dvrr_stack+252, dvrr_stack+95);
 tmp = dvrr_stack + 75;
 target_ptr = Libderiv->deriv2_classes[0][2][80];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+81, dvrr_stack+268, dvrr_stack+262);
 tmp = dvrr_stack + 81;
 target_ptr = Libderiv->deriv2_classes[0][3][80];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+91, dvrr_stack+286, dvrr_stack+283);
 tmp = dvrr_stack + 91;
 target_ptr = Libderiv->deriv2_classes[0][2][79];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+108, dvrr_stack+125, dvrr_stack+119);
 tmp = dvrr_stack + 108;
 target_ptr = Libderiv->deriv2_classes[0][3][79];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+118, dvrr_stack+466, dvrr_stack+140);
 tmp = dvrr_stack + 118;
 target_ptr = Libderiv->deriv2_classes[0][2][78];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+124, dvrr_stack+476, dvrr_stack+197);
 tmp = dvrr_stack + 124;
 target_ptr = Libderiv->deriv2_classes[0][3][78];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,6,dvrr_stack+195, dvrr_stack+545, NULL);
 tmp = dvrr_stack + 195;
 target_ptr = Libderiv->deriv2_classes[0][2][35];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,10,dvrr_stack+249, dvrr_stack+698, NULL);
 tmp = dvrr_stack + 249;
 target_ptr = Libderiv->deriv2_classes[0][3][35];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+134, dvrr_stack+728, NULL);
 tmp = dvrr_stack + 134;
 target_ptr = Libderiv->deriv2_classes[0][2][34];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+259, dvrr_stack+746, NULL);
 tmp = dvrr_stack + 259;
 target_ptr = Libderiv->deriv2_classes[0][3][34];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+269, dvrr_stack+776, NULL);
 tmp = dvrr_stack + 269;
 target_ptr = Libderiv->deriv2_classes[0][2][33];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+275, dvrr_stack+491, NULL);
 tmp = dvrr_stack + 275;
 target_ptr = Libderiv->deriv2_classes[0][3][33];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+285, dvrr_stack+617, NULL);
 tmp = dvrr_stack + 285;
 target_ptr = Libderiv->deriv2_classes[0][2][32];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+332, dvrr_stack+635, NULL);
 tmp = dvrr_stack + 332;
 target_ptr = Libderiv->deriv2_classes[0][3][32];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+342, dvrr_stack+665, NULL);
 tmp = dvrr_stack + 342;
 target_ptr = Libderiv->deriv2_classes[0][2][31];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+410, dvrr_stack+794, NULL);
 tmp = dvrr_stack + 410;
 target_ptr = Libderiv->deriv2_classes[0][3][31];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+348, dvrr_stack+231, NULL);
 tmp = dvrr_stack + 348;
 target_ptr = Libderiv->deriv_classes[0][2][2];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+420, dvrr_stack+521, NULL);
 tmp = dvrr_stack + 420;
 target_ptr = Libderiv->deriv2_classes[0][2][30];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+426, dvrr_stack+824, NULL);
 tmp = dvrr_stack + 426;
 target_ptr = Libderiv->deriv2_classes[0][3][30];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+436, dvrr_stack+590, NULL);
 tmp = dvrr_stack + 436;
 target_ptr = Libderiv->deriv2_classes[0][2][26];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+442, dvrr_stack+890, NULL);
 tmp = dvrr_stack + 442;
 target_ptr = Libderiv->deriv2_classes[0][3][26];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,6,dvrr_stack+452, dvrr_stack+545, NULL);
 tmp = dvrr_stack + 452;
 target_ptr = Libderiv->deriv2_classes[0][2][23];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,10,dvrr_stack+458, dvrr_stack+698, NULL);
 tmp = dvrr_stack + 458;
 target_ptr = Libderiv->deriv2_classes[0][3][23];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+468, dvrr_stack+728, NULL);
 tmp = dvrr_stack + 468;
 target_ptr = Libderiv->deriv2_classes[0][2][22];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+474, dvrr_stack+746, NULL);
 tmp = dvrr_stack + 474;
 target_ptr = Libderiv->deriv2_classes[0][3][22];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+484, dvrr_stack+776, NULL);
 tmp = dvrr_stack + 484;
 target_ptr = Libderiv->deriv2_classes[0][2][21];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1010, dvrr_stack+491, NULL);
 tmp = dvrr_stack + 1010;
 target_ptr = Libderiv->deriv2_classes[0][3][21];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1020, dvrr_stack+617, NULL);
 tmp = dvrr_stack + 1020;
 target_ptr = Libderiv->deriv2_classes[0][2][20];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1026, dvrr_stack+635, NULL);
 tmp = dvrr_stack + 1026;
 target_ptr = Libderiv->deriv2_classes[0][3][20];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1036, dvrr_stack+665, NULL);
 tmp = dvrr_stack + 1036;
 target_ptr = Libderiv->deriv2_classes[0][2][19];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1042, dvrr_stack+794, NULL);
 tmp = dvrr_stack + 1042;
 target_ptr = Libderiv->deriv2_classes[0][3][19];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1052, dvrr_stack+231, NULL);
 tmp = dvrr_stack + 1052;
 target_ptr = Libderiv->deriv_classes[0][2][1];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1058, dvrr_stack+521, NULL);
 tmp = dvrr_stack + 1058;
 target_ptr = Libderiv->deriv2_classes[0][2][18];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1064, dvrr_stack+824, NULL);
 tmp = dvrr_stack + 1064;
 target_ptr = Libderiv->deriv2_classes[0][3][18];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1074, dvrr_stack+590, NULL);
 tmp = dvrr_stack + 1074;
 target_ptr = Libderiv->deriv2_classes[0][2][14];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1080, dvrr_stack+890, NULL);
 tmp = dvrr_stack + 1080;
 target_ptr = Libderiv->deriv2_classes[0][3][14];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1090, dvrr_stack+15, NULL);
 tmp = dvrr_stack + 1090;
 target_ptr = Libderiv->deriv2_classes[0][2][13];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1096, dvrr_stack+980, NULL);
 tmp = dvrr_stack + 1096;
 target_ptr = Libderiv->deriv2_classes[0][3][13];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,6,dvrr_stack+1106, dvrr_stack+545, NULL);
 tmp = dvrr_stack + 1106;
 target_ptr = Libderiv->deriv2_classes[0][2][11];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,10,dvrr_stack+545, dvrr_stack+698, NULL);
 tmp = dvrr_stack + 545;
 target_ptr = Libderiv->deriv2_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+555, dvrr_stack+728, NULL);
 tmp = dvrr_stack + 555;
 target_ptr = Libderiv->deriv2_classes[0][2][10];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+693, dvrr_stack+746, NULL);
 tmp = dvrr_stack + 693;
 target_ptr = Libderiv->deriv2_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+703, dvrr_stack+776, NULL);
 tmp = dvrr_stack + 703;
 target_ptr = Libderiv->deriv2_classes[0][2][9];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+709, dvrr_stack+491, NULL);
 tmp = dvrr_stack + 709;
 target_ptr = Libderiv->deriv2_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+719, dvrr_stack+617, NULL);
 tmp = dvrr_stack + 719;
 target_ptr = Libderiv->deriv2_classes[0][2][8];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+614, dvrr_stack+635, NULL);
 tmp = dvrr_stack + 614;
 target_ptr = Libderiv->deriv2_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+624, dvrr_stack+665, NULL);
 tmp = dvrr_stack + 624;
 target_ptr = Libderiv->deriv2_classes[0][2][7];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+630, dvrr_stack+794, NULL);
 tmp = dvrr_stack + 630;
 target_ptr = Libderiv->deriv2_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+640, dvrr_stack+231, NULL);
 tmp = dvrr_stack + 640;
 target_ptr = Libderiv->deriv_classes[0][2][0];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+646, dvrr_stack+521, NULL);
 tmp = dvrr_stack + 646;
 target_ptr = Libderiv->deriv2_classes[0][2][6];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+652, dvrr_stack+824, NULL);
 tmp = dvrr_stack + 652;
 target_ptr = Libderiv->deriv2_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+662, dvrr_stack+590, NULL);
 tmp = dvrr_stack + 662;
 target_ptr = Libderiv->deriv2_classes[0][2][2];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+587, dvrr_stack+890, NULL);
 tmp = dvrr_stack + 587;
 target_ptr = Libderiv->deriv2_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+890, dvrr_stack+15, NULL);
 tmp = dvrr_stack + 890;
 target_ptr = Libderiv->deriv2_classes[0][2][1];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+896, dvrr_stack+980, NULL);
 tmp = dvrr_stack + 896;
 target_ptr = Libderiv->deriv2_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+906, dvrr_stack+563, NULL);
 tmp = dvrr_stack + 906;
 target_ptr = Libderiv->deriv2_classes[0][2][0];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+561, dvrr_stack+854, NULL);
 tmp = dvrr_stack + 561;
 target_ptr = Libderiv->deriv2_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];


}

