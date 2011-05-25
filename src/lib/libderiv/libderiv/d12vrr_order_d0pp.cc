#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (d0|pp) integrals */

void d12vrr_order_d0pp(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+12, dvrr_stack+0, Data->F+0, Data->F+1, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+6, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+39, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+42, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+45, dvrr_stack+3, dvrr_stack+42, NULL, NULL, Data->F+3);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+54, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+63, dvrr_stack+54, dvrr_stack+45, dvrr_stack+0, dvrr_stack+3, dvrr_stack+39);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+81, dvrr_stack+3, dvrr_stack+42, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+87, dvrr_stack+6, dvrr_stack+81, NULL, NULL, dvrr_stack+3);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+105, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+108, dvrr_stack+42, dvrr_stack+105, Data->F+3, Data->F+4, NULL);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+114, dvrr_stack+81, dvrr_stack+108, NULL, NULL, dvrr_stack+42);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+132, dvrr_stack+87, dvrr_stack+114, dvrr_stack+6, dvrr_stack+81, dvrr_stack+45);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+168, dvrr_stack+21, dvrr_stack+87, dvrr_stack+15, dvrr_stack+6, dvrr_stack+54);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+204, dvrr_stack+168, dvrr_stack+132, dvrr_stack+21, dvrr_stack+87, dvrr_stack+63);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+264, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+267, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+1);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+276, dvrr_stack+267, dvrr_stack+54, dvrr_stack+12, dvrr_stack+0, dvrr_stack+264);
 tmp = dvrr_stack + 276;
 target_ptr = Libderiv->dvrr_classes[2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+294,dvrr_stack+168,dvrr_stack+276,6);


 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+348, dvrr_stack+6, dvrr_stack+81, dvrr_stack+0, dvrr_stack+3, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+358, dvrr_stack+15, dvrr_stack+6, dvrr_stack+12, dvrr_stack+0, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+368, dvrr_stack+81, dvrr_stack+108, dvrr_stack+3, dvrr_stack+42, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+378, dvrr_stack+348, dvrr_stack+368, NULL, NULL, dvrr_stack+81);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+408, dvrr_stack+358, dvrr_stack+348, NULL, NULL, dvrr_stack+6);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+438, dvrr_stack+408, dvrr_stack+378, dvrr_stack+358, dvrr_stack+348, dvrr_stack+87);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+498,dvrr_stack+438,dvrr_stack+168,6);


 /* compute (2 0 | 1 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pd(Libderiv->CD,dvrr_stack+606,dvrr_stack+498,dvrr_stack+294,6);


 /* compute (2 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,18,dvrr_stack+714, dvrr_stack+606, dvrr_stack+276);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+768, dvrr_stack+348, dvrr_stack+368, dvrr_stack+6, dvrr_stack+81, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+783, dvrr_stack+358, dvrr_stack+348, dvrr_stack+15, dvrr_stack+6, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+105, dvrr_stack+0, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+358, dvrr_stack+108, dvrr_stack+6, dvrr_stack+42, dvrr_stack+105, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+798, dvrr_stack+368, dvrr_stack+358, dvrr_stack+81, dvrr_stack+108, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+813, dvrr_stack+768, dvrr_stack+798, NULL, NULL, dvrr_stack+368);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+858, dvrr_stack+783, dvrr_stack+768, NULL, NULL, dvrr_stack+348);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+903, dvrr_stack+858, dvrr_stack+813, dvrr_stack+783, dvrr_stack+768, dvrr_stack+378);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+993,dvrr_stack+903,dvrr_stack+438,6);


 /* compute (2 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+1173,dvrr_stack+993,dvrr_stack+498,6);


 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,36,dvrr_stack+768, dvrr_stack+1173, dvrr_stack+168);

 /* compute (2 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,18,dvrr_stack+1389, dvrr_stack+606, dvrr_stack+276);

 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,36,dvrr_stack+1443, dvrr_stack+1173, dvrr_stack+168);

 /* compute (2 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,18,dvrr_stack+1551, dvrr_stack+606, dvrr_stack+276);

 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,36,dvrr_stack+606, dvrr_stack+1173, dvrr_stack+168);

 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+1173, dvrr_stack+0, dvrr_stack+264, Data->F+0, Data->F+1, NULL);

 /* compute (2 0 | 0 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0p(Libderiv->CD,dvrr_stack+1179,dvrr_stack+276,dvrr_stack+1173,6);


 /* compute (2 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+1197, dvrr_stack+1179, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+1203, dvrr_stack+498, NULL);
 tmp = dvrr_stack + 1203;
 target_ptr = Libderiv->deriv_classes[2][2][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+1239, dvrr_stack+294, NULL);
 tmp = dvrr_stack + 1239;
 target_ptr = Libderiv->deriv_classes[2][1][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+1257, dvrr_stack+993, NULL);

 /* compute (2 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+1317, dvrr_stack+1179, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+1323, dvrr_stack+498, NULL);
 tmp = dvrr_stack + 1323;
 target_ptr = Libderiv->deriv_classes[2][2][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+1359, dvrr_stack+294, NULL);
 tmp = dvrr_stack + 1359;
 target_ptr = Libderiv->deriv_classes[2][1][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+1605, dvrr_stack+993, NULL);

 /* compute (2 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+1377, dvrr_stack+1179, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+1665, dvrr_stack+498, NULL);
 tmp = dvrr_stack + 1665;
 target_ptr = Libderiv->deriv_classes[2][2][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+1179, dvrr_stack+294, NULL);
 tmp = dvrr_stack + 1179;
 target_ptr = Libderiv->deriv_classes[2][1][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+498, dvrr_stack+993, NULL);

 /* compute (2 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_0(Data,6,1,dvrr_stack+1383, dvrr_stack+276, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+993, dvrr_stack+438, dvrr_stack+276);
 tmp = dvrr_stack + 993;
 target_ptr = Libderiv->deriv_classes[2][2][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+1029, dvrr_stack+168, dvrr_stack+1173);
 tmp = dvrr_stack + 1029;
 target_ptr = Libderiv->deriv_classes[2][1][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+1047, dvrr_stack+903, dvrr_stack+168);

 /* compute (2 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_0(Data,6,1,dvrr_stack+1107, dvrr_stack+276, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+1113, dvrr_stack+438, dvrr_stack+276);
 tmp = dvrr_stack + 1113;
 target_ptr = Libderiv->deriv_classes[2][2][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+1149, dvrr_stack+168, dvrr_stack+1173);
 tmp = dvrr_stack + 1149;
 target_ptr = Libderiv->deriv_classes[2][1][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+1701, dvrr_stack+903, dvrr_stack+168);

 /* compute (2 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_0(Data,6,1,dvrr_stack+1167, dvrr_stack+276, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+294, dvrr_stack+438, dvrr_stack+276);
 tmp = dvrr_stack + 294;
 target_ptr = Libderiv->deriv_classes[2][2][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+330, dvrr_stack+168, dvrr_stack+1173);
 tmp = dvrr_stack + 330;
 target_ptr = Libderiv->deriv_classes[2][1][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+1761, dvrr_stack+903, dvrr_stack+168);

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+558,dvrr_stack+21,dvrr_stack+267,3);


 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+585, dvrr_stack+558, NULL);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+594, dvrr_stack+264, dvrr_stack+39, Data->F+1, Data->F+2, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+876, dvrr_stack+276, dvrr_stack+63, dvrr_stack+267, dvrr_stack+54, dvrr_stack+594);

 /* compute (3 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+1821,dvrr_stack+204,dvrr_stack+876,10);


 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+906, dvrr_stack+1821, NULL);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+936,dvrr_stack+408,dvrr_stack+21,3);


 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+1911, dvrr_stack+936, NULL);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1929, dvrr_stack+368, dvrr_stack+358, NULL, NULL, dvrr_stack+108);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1959, dvrr_stack+378, dvrr_stack+1929, dvrr_stack+348, dvrr_stack+368, dvrr_stack+114);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+2019, dvrr_stack+438, dvrr_stack+1959, dvrr_stack+408, dvrr_stack+378, dvrr_stack+132);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+2119,dvrr_stack+2019,dvrr_stack+204,10);


 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+438, dvrr_stack+2119, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+348, dvrr_stack+558, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+357, dvrr_stack+1821, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+387, dvrr_stack+936, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+1929, dvrr_stack+2119, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+1989, dvrr_stack+558, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+2299, dvrr_stack+1821, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+1821, dvrr_stack+936, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+1839, dvrr_stack+2119, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+2119, dvrr_stack+21, dvrr_stack+0);

 /* compute (3 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+2128, dvrr_stack+1173, dvrr_stack+594, dvrr_stack+0, dvrr_stack+264, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,10,1,dvrr_stack+2138, dvrr_stack+204, dvrr_stack+2128);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+2168, dvrr_stack+408, dvrr_stack+267);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+2186, dvrr_stack+2019, dvrr_stack+876);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+2246, dvrr_stack+21, dvrr_stack+0);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,10,1,dvrr_stack+2255, dvrr_stack+204, dvrr_stack+2128);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+936, dvrr_stack+408, dvrr_stack+267);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+2329, dvrr_stack+2019, dvrr_stack+876);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+2285, dvrr_stack+21, dvrr_stack+0);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,10,1,dvrr_stack+954, dvrr_stack+204, dvrr_stack+2128);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+558, dvrr_stack+408, dvrr_stack+267);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+2389, dvrr_stack+2019, dvrr_stack+876);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+984, dvrr_stack+276, dvrr_stack+12);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+1173, dvrr_stack+39, dvrr_stack+0, Data->F+2, Data->F+3, NULL);

 /* compute (3 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+2128, dvrr_stack+594, dvrr_stack+1173, dvrr_stack+264, dvrr_stack+39, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+576, dvrr_stack+42, dvrr_stack+105, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+405, dvrr_stack+45, dvrr_stack+576, dvrr_stack+3, dvrr_stack+42, dvrr_stack+0);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+1998, dvrr_stack+63, dvrr_stack+405, dvrr_stack+54, dvrr_stack+45, dvrr_stack+1173);

 /* compute (4 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+2028, dvrr_stack+876, dvrr_stack+1998, dvrr_stack+276, dvrr_stack+63, dvrr_stack+2128);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,3,dvrr_stack+39, dvrr_stack+2028, dvrr_stack+276);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+2073, dvrr_stack+168, dvrr_stack+15);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2091, dvrr_stack+108, dvrr_stack+6, NULL, NULL, dvrr_stack+105);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+2449, dvrr_stack+114, dvrr_stack+2091, dvrr_stack+81, dvrr_stack+108, dvrr_stack+576);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+2485, dvrr_stack+132, dvrr_stack+2449, dvrr_stack+87, dvrr_stack+114, dvrr_stack+405);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+2545, dvrr_stack+204, dvrr_stack+2485, dvrr_stack+168, dvrr_stack+132, dvrr_stack+1998);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+2449, dvrr_stack+2545, dvrr_stack+168);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+576, dvrr_stack+276, dvrr_stack+12);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,3,dvrr_stack+1998, dvrr_stack+2028, dvrr_stack+276);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+405, dvrr_stack+168, dvrr_stack+15);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+69, dvrr_stack+2545, dvrr_stack+168);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+423, dvrr_stack+276, dvrr_stack+12);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,3,dvrr_stack+2509, dvrr_stack+2028, dvrr_stack+276);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+276, dvrr_stack+168, dvrr_stack+15);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+2635, dvrr_stack+2545, dvrr_stack+168);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+2028, dvrr_stack+204, dvrr_stack+21);
 tmp = dvrr_stack + 2028;
 target_ptr = Libderiv->deriv_classes[2][2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+2539, dvrr_stack+204, dvrr_stack+21);
 tmp = dvrr_stack + 2539;
 target_ptr = Libderiv->deriv_classes[2][2][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2575, dvrr_stack+204, dvrr_stack+21);
 tmp = dvrr_stack + 2575;
 target_ptr = Libderiv->deriv_classes[2][2][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,18,dvrr_stack+2611, dvrr_stack+714, NULL);
 tmp = dvrr_stack + 2611;
 target_ptr = Libderiv->deriv2_classes[2][1][143];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,36,dvrr_stack+0, dvrr_stack+768, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[2][2][143];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,18,dvrr_stack+2091, dvrr_stack+714, NULL);
 tmp = dvrr_stack + 2091;
 target_ptr = Libderiv->deriv2_classes[2][1][131];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,36,dvrr_stack+129, dvrr_stack+768, NULL);
 tmp = dvrr_stack + 129;
 target_ptr = Libderiv->deriv2_classes[2][2][131];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+165, dvrr_stack+1389, NULL);
 tmp = dvrr_stack + 165;
 target_ptr = Libderiv->deriv2_classes[2][1][130];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+183, dvrr_stack+1443, NULL);
 tmp = dvrr_stack + 183;
 target_ptr = Libderiv->deriv2_classes[2][2][130];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,18,dvrr_stack+219, dvrr_stack+714, NULL);
 tmp = dvrr_stack + 219;
 target_ptr = Libderiv->deriv2_classes[2][1][119];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,36,dvrr_stack+714, dvrr_stack+768, NULL);
 tmp = dvrr_stack + 714;
 target_ptr = Libderiv->deriv2_classes[2][2][119];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+750, dvrr_stack+1389, NULL);
 tmp = dvrr_stack + 750;
 target_ptr = Libderiv->deriv2_classes[2][1][118];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+1389, dvrr_stack+1443, NULL);
 tmp = dvrr_stack + 1389;
 target_ptr = Libderiv->deriv2_classes[2][2][118];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+1425, dvrr_stack+1551, NULL);
 tmp = dvrr_stack + 1425;
 target_ptr = Libderiv->deriv2_classes[2][1][117];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+1443, dvrr_stack+606, NULL);
 tmp = dvrr_stack + 1443;
 target_ptr = Libderiv->deriv2_classes[2][2][117];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+1479, dvrr_stack+1203, dvrr_stack+1197);
 tmp = dvrr_stack + 1479;
 target_ptr = Libderiv->deriv2_classes[2][1][107];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+1497, dvrr_stack+1257, dvrr_stack+1239);
 tmp = dvrr_stack + 1497;
 target_ptr = Libderiv->deriv2_classes[2][2][107];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+1533, dvrr_stack+1323, dvrr_stack+1317);
 tmp = dvrr_stack + 1533;
 target_ptr = Libderiv->deriv2_classes[2][1][106];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+1551, dvrr_stack+1605, dvrr_stack+1359);
 tmp = dvrr_stack + 1551;
 target_ptr = Libderiv->deriv2_classes[2][2][106];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+1587, dvrr_stack+1665, dvrr_stack+1377);
 tmp = dvrr_stack + 1587;
 target_ptr = Libderiv->deriv2_classes[2][1][105];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+768, dvrr_stack+498, dvrr_stack+1179);
 tmp = dvrr_stack + 768;
 target_ptr = Libderiv->deriv2_classes[2][2][105];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+804, dvrr_stack+993, dvrr_stack+1383);
 tmp = dvrr_stack + 804;
 target_ptr = Libderiv->deriv2_classes[2][1][104];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+822, dvrr_stack+1047, dvrr_stack+1029);
 tmp = dvrr_stack + 822;
 target_ptr = Libderiv->deriv2_classes[2][2][104];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+858, dvrr_stack+1203, dvrr_stack+1197);
 tmp = dvrr_stack + 858;
 target_ptr = Libderiv->deriv2_classes[2][1][95];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+594, dvrr_stack+1257, dvrr_stack+1239);
 tmp = dvrr_stack + 594;
 target_ptr = Libderiv->deriv2_classes[2][2][95];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+237, dvrr_stack+1323, dvrr_stack+1317);
 tmp = dvrr_stack + 237;
 target_ptr = Libderiv->deriv2_classes[2][1][94];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+630, dvrr_stack+1605, dvrr_stack+1359);
 tmp = dvrr_stack + 630;
 target_ptr = Libderiv->deriv2_classes[2][2][94];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+666, dvrr_stack+1665, dvrr_stack+1377);
 tmp = dvrr_stack + 666;
 target_ptr = Libderiv->deriv2_classes[2][1][93];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+2695, dvrr_stack+498, dvrr_stack+1179);
 tmp = dvrr_stack + 2695;
 target_ptr = Libderiv->deriv2_classes[2][2][93];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+684, dvrr_stack+993, dvrr_stack+1383);
 tmp = dvrr_stack + 684;
 target_ptr = Libderiv->deriv2_classes[2][1][92];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+2731, dvrr_stack+1047, dvrr_stack+1029);
 tmp = dvrr_stack + 2731;
 target_ptr = Libderiv->deriv2_classes[2][2][92];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+2767, dvrr_stack+1113, dvrr_stack+1107);
 tmp = dvrr_stack + 2767;
 target_ptr = Libderiv->deriv2_classes[2][1][91];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+2785, dvrr_stack+1701, dvrr_stack+1149);
 tmp = dvrr_stack + 2785;
 target_ptr = Libderiv->deriv2_classes[2][2][91];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+2821, dvrr_stack+1203, dvrr_stack+1197);
 tmp = dvrr_stack + 2821;
 target_ptr = Libderiv->deriv2_classes[2][1][83];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1197, dvrr_stack+1257, dvrr_stack+1239);
 tmp = dvrr_stack + 1197;
 target_ptr = Libderiv->deriv2_classes[2][2][83];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+1233, dvrr_stack+1323, dvrr_stack+1317);
 tmp = dvrr_stack + 1233;
 target_ptr = Libderiv->deriv2_classes[2][1][82];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1251, dvrr_stack+1605, dvrr_stack+1359);
 tmp = dvrr_stack + 1251;
 target_ptr = Libderiv->deriv2_classes[2][2][82];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+1605, dvrr_stack+1665, dvrr_stack+1377);
 tmp = dvrr_stack + 1605;
 target_ptr = Libderiv->deriv2_classes[2][1][81];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1623, dvrr_stack+498, dvrr_stack+1179);
 tmp = dvrr_stack + 1623;
 target_ptr = Libderiv->deriv2_classes[2][2][81];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+498, dvrr_stack+993, dvrr_stack+1383);
 tmp = dvrr_stack + 498;
 target_ptr = Libderiv->deriv2_classes[2][1][80];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+993, dvrr_stack+1047, dvrr_stack+1029);
 tmp = dvrr_stack + 993;
 target_ptr = Libderiv->deriv2_classes[2][2][80];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+1029, dvrr_stack+1113, dvrr_stack+1107);
 tmp = dvrr_stack + 1029;
 target_ptr = Libderiv->deriv2_classes[2][1][79];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1047, dvrr_stack+1701, dvrr_stack+1149);
 tmp = dvrr_stack + 1047;
 target_ptr = Libderiv->deriv2_classes[2][2][79];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+1083, dvrr_stack+294, dvrr_stack+1167);
 tmp = dvrr_stack + 1083;
 target_ptr = Libderiv->deriv2_classes[2][1][78];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+294, dvrr_stack+1761, dvrr_stack+330);
 tmp = dvrr_stack + 294;
 target_ptr = Libderiv->deriv2_classes[2][2][78];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,3,dvrr_stack+330, dvrr_stack+906, dvrr_stack+585);
 tmp = dvrr_stack + 330;
 target_ptr = Libderiv->deriv2_classes[2][1][35];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1101, dvrr_stack+438, dvrr_stack+1911);
 tmp = dvrr_stack + 1101;
 target_ptr = Libderiv->deriv2_classes[2][2][35];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+1137, dvrr_stack+357, dvrr_stack+348);
 tmp = dvrr_stack + 1137;
 target_ptr = Libderiv->deriv2_classes[2][1][34];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1155, dvrr_stack+1929, dvrr_stack+387);
 tmp = dvrr_stack + 1155;
 target_ptr = Libderiv->deriv2_classes[2][2][34];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+516, dvrr_stack+2299, dvrr_stack+1989);
 tmp = dvrr_stack + 516;
 target_ptr = Libderiv->deriv2_classes[2][1][33];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1659, dvrr_stack+1839, dvrr_stack+1821);
 tmp = dvrr_stack + 1659;
 target_ptr = Libderiv->deriv2_classes[2][2][33];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+534, dvrr_stack+2138, dvrr_stack+2119);
 tmp = dvrr_stack + 534;
 target_ptr = Libderiv->deriv2_classes[2][1][32];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1695, dvrr_stack+2186, dvrr_stack+2168);
 tmp = dvrr_stack + 1695;
 target_ptr = Libderiv->deriv2_classes[2][2][32];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+1731, dvrr_stack+2255, dvrr_stack+2246);
 tmp = dvrr_stack + 1731;
 target_ptr = Libderiv->deriv2_classes[2][1][31];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1749, dvrr_stack+2329, dvrr_stack+936);
 tmp = dvrr_stack + 1749;
 target_ptr = Libderiv->deriv2_classes[2][2][31];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+1785, dvrr_stack+876, dvrr_stack+267);
 tmp = dvrr_stack + 1785;
 target_ptr = Libderiv->deriv_classes[2][1][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+1803, dvrr_stack+954, dvrr_stack+2285);
 tmp = dvrr_stack + 1803;
 target_ptr = Libderiv->deriv2_classes[2][1][30];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1287, dvrr_stack+2389, dvrr_stack+558);
 tmp = dvrr_stack + 1287;
 target_ptr = Libderiv->deriv2_classes[2][2][30];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+1323, dvrr_stack+39, dvrr_stack+984);
 tmp = dvrr_stack + 1323;
 target_ptr = Libderiv->deriv2_classes[2][1][26];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1341, dvrr_stack+2449, dvrr_stack+2073);
 tmp = dvrr_stack + 1341;
 target_ptr = Libderiv->deriv2_classes[2][2][26];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,3,dvrr_stack+2839, dvrr_stack+906, dvrr_stack+585);
 tmp = dvrr_stack + 2839;
 target_ptr = Libderiv->deriv2_classes[2][1][23];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,6,dvrr_stack+2857, dvrr_stack+438, dvrr_stack+1911);
 tmp = dvrr_stack + 2857;
 target_ptr = Libderiv->deriv2_classes[2][2][23];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+2893, dvrr_stack+357, dvrr_stack+348);
 tmp = dvrr_stack + 2893;
 target_ptr = Libderiv->deriv2_classes[2][1][22];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+2911, dvrr_stack+1929, dvrr_stack+387);
 tmp = dvrr_stack + 2911;
 target_ptr = Libderiv->deriv2_classes[2][2][22];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+2947, dvrr_stack+2299, dvrr_stack+1989);
 tmp = dvrr_stack + 2947;
 target_ptr = Libderiv->deriv2_classes[2][1][21];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+2965, dvrr_stack+1839, dvrr_stack+1821);
 tmp = dvrr_stack + 2965;
 target_ptr = Libderiv->deriv2_classes[2][2][21];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+3001, dvrr_stack+2138, dvrr_stack+2119);
 tmp = dvrr_stack + 3001;
 target_ptr = Libderiv->deriv2_classes[2][1][20];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3019, dvrr_stack+2186, dvrr_stack+2168);
 tmp = dvrr_stack + 3019;
 target_ptr = Libderiv->deriv2_classes[2][2][20];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+3055, dvrr_stack+2255, dvrr_stack+2246);
 tmp = dvrr_stack + 3055;
 target_ptr = Libderiv->deriv2_classes[2][1][19];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3073, dvrr_stack+2329, dvrr_stack+936);
 tmp = dvrr_stack + 3073;
 target_ptr = Libderiv->deriv2_classes[2][2][19];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+3109, dvrr_stack+876, dvrr_stack+267);
 tmp = dvrr_stack + 3109;
 target_ptr = Libderiv->deriv_classes[2][1][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+3127, dvrr_stack+954, dvrr_stack+2285);
 tmp = dvrr_stack + 3127;
 target_ptr = Libderiv->deriv2_classes[2][1][18];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3145, dvrr_stack+2389, dvrr_stack+558);
 tmp = dvrr_stack + 3145;
 target_ptr = Libderiv->deriv2_classes[2][2][18];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+3181, dvrr_stack+39, dvrr_stack+984);
 tmp = dvrr_stack + 3181;
 target_ptr = Libderiv->deriv2_classes[2][1][14];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3199, dvrr_stack+2449, dvrr_stack+2073);
 tmp = dvrr_stack + 3199;
 target_ptr = Libderiv->deriv2_classes[2][2][14];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+3235, dvrr_stack+1998, dvrr_stack+576);
 tmp = dvrr_stack + 3235;
 target_ptr = Libderiv->deriv2_classes[2][1][13];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3253, dvrr_stack+69, dvrr_stack+405);
 tmp = dvrr_stack + 3253;
 target_ptr = Libderiv->deriv2_classes[2][2][13];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,3,dvrr_stack+3289, dvrr_stack+906, dvrr_stack+585);
 tmp = dvrr_stack + 3289;
 target_ptr = Libderiv->deriv2_classes[2][1][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,6,dvrr_stack+3307, dvrr_stack+438, dvrr_stack+1911);
 tmp = dvrr_stack + 3307;
 target_ptr = Libderiv->deriv2_classes[2][2][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+906, dvrr_stack+357, dvrr_stack+348);
 tmp = dvrr_stack + 906;
 target_ptr = Libderiv->deriv2_classes[2][1][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+348, dvrr_stack+1929, dvrr_stack+387);
 tmp = dvrr_stack + 348;
 target_ptr = Libderiv->deriv2_classes[2][2][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+384, dvrr_stack+2299, dvrr_stack+1989);
 tmp = dvrr_stack + 384;
 target_ptr = Libderiv->deriv2_classes[2][1][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+432, dvrr_stack+1839, dvrr_stack+1821);
 tmp = dvrr_stack + 432;
 target_ptr = Libderiv->deriv2_classes[2][2][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+1821, dvrr_stack+2138, dvrr_stack+2119);
 tmp = dvrr_stack + 1821;
 target_ptr = Libderiv->deriv2_classes[2][1][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1839, dvrr_stack+2186, dvrr_stack+2168);
 tmp = dvrr_stack + 1839;
 target_ptr = Libderiv->deriv2_classes[2][2][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+1875, dvrr_stack+2255, dvrr_stack+2246);
 tmp = dvrr_stack + 1875;
 target_ptr = Libderiv->deriv2_classes[2][1][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1893, dvrr_stack+2329, dvrr_stack+936);
 tmp = dvrr_stack + 1893;
 target_ptr = Libderiv->deriv2_classes[2][2][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+1929, dvrr_stack+876, dvrr_stack+267);
 tmp = dvrr_stack + 1929;
 target_ptr = Libderiv->deriv_classes[2][1][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+876, dvrr_stack+954, dvrr_stack+2285);
 tmp = dvrr_stack + 876;
 target_ptr = Libderiv->deriv2_classes[2][1][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1947, dvrr_stack+2389, dvrr_stack+558);
 tmp = dvrr_stack + 1947;
 target_ptr = Libderiv->deriv2_classes[2][2][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+924, dvrr_stack+39, dvrr_stack+984);
 tmp = dvrr_stack + 924;
 target_ptr = Libderiv->deriv2_classes[2][1][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+942, dvrr_stack+2449, dvrr_stack+2073);
 tmp = dvrr_stack + 942;
 target_ptr = Libderiv->deriv2_classes[2][2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+552, dvrr_stack+1998, dvrr_stack+576);
 tmp = dvrr_stack + 552;
 target_ptr = Libderiv->deriv2_classes[2][1][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1983, dvrr_stack+69, dvrr_stack+405);
 tmp = dvrr_stack + 1983;
 target_ptr = Libderiv->deriv2_classes[2][2][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+402, dvrr_stack+2509, dvrr_stack+423);
 tmp = dvrr_stack + 402;
 target_ptr = Libderiv->deriv2_classes[2][1][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2109, dvrr_stack+2635, dvrr_stack+276);
 tmp = dvrr_stack + 2109;
 target_ptr = Libderiv->deriv2_classes[2][2][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];


}

