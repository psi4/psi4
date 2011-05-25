#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|ff) integrals */

void d12vrr_order_00ff(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+9, dvrr_stack+0, dvrr_stack+6, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+9, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+34, dvrr_stack+31, dvrr_stack+3, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+34, dvrr_stack+15, dvrr_stack+31, dvrr_stack+3, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+50, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+53, dvrr_stack+6, dvrr_stack+50, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+59, dvrr_stack+9, dvrr_stack+53, dvrr_stack+0, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+69, dvrr_stack+21, dvrr_stack+59, dvrr_stack+15, dvrr_stack+9, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+84, dvrr_stack+40, dvrr_stack+21, dvrr_stack+34, dvrr_stack+15, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+99, dvrr_stack+84, dvrr_stack+69, dvrr_stack+40, dvrr_stack+21, NULL);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+120, dvrr_stack+50, dvrr_stack+0, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+126, dvrr_stack+53, dvrr_stack+120, dvrr_stack+6, dvrr_stack+50, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+136, dvrr_stack+59, dvrr_stack+126, dvrr_stack+9, dvrr_stack+53, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+151, dvrr_stack+69, dvrr_stack+136, dvrr_stack+21, dvrr_stack+59, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+172, dvrr_stack+99, dvrr_stack+151, dvrr_stack+84, dvrr_stack+69, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+9, dvrr_stack+6, dvrr_stack+31, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+200, dvrr_stack+9, dvrr_stack+34, dvrr_stack+6, dvrr_stack+31, NULL);
 tmp = dvrr_stack + 200;
 target_ptr = Libderiv->dvrr_classes[0][3];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+210, dvrr_stack+200, dvrr_stack+40, dvrr_stack+9, dvrr_stack+34, NULL);
 tmp = dvrr_stack + 210;
 target_ptr = Libderiv->dvrr_classes[0][4];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+225, dvrr_stack+210, dvrr_stack+84, dvrr_stack+200, dvrr_stack+40, NULL);
 tmp = dvrr_stack + 225;
 target_ptr = Libderiv->dvrr_classes[0][5];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+246, dvrr_stack+225, dvrr_stack+99, dvrr_stack+210, dvrr_stack+84, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+274, dvrr_stack+246, dvrr_stack+172, NULL, NULL, dvrr_stack+99);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+358,dvrr_stack+210,dvrr_stack+200,1);


 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+388,dvrr_stack+225,dvrr_stack+210,1);


 /* compute (0 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+433,dvrr_stack+388,dvrr_stack+358,1);


 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,10,dvrr_stack+493, dvrr_stack+433, dvrr_stack+200);

 /* compute (0 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+523,dvrr_stack+246,dvrr_stack+225,1);


 /* compute (0 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+586,dvrr_stack+523,dvrr_stack+388,1);


 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,15,dvrr_stack+676, dvrr_stack+586, dvrr_stack+210);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+721, dvrr_stack+246, dvrr_stack+172, dvrr_stack+225, dvrr_stack+99, NULL);

 /* compute (0 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+757,dvrr_stack+721,dvrr_stack+246,1);


 /* compute (0 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+841,dvrr_stack+757,dvrr_stack+523,1);


 /* compute (0 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,21,dvrr_stack+967, dvrr_stack+841, dvrr_stack+225);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1030, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1033, dvrr_stack+0, dvrr_stack+1030, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1039, dvrr_stack+120, dvrr_stack+1033, dvrr_stack+50, dvrr_stack+0, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1049, dvrr_stack+126, dvrr_stack+1039, dvrr_stack+53, dvrr_stack+120, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1064, dvrr_stack+136, dvrr_stack+1049, dvrr_stack+59, dvrr_stack+126, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1030, dvrr_stack+151, dvrr_stack+1064, dvrr_stack+69, dvrr_stack+136, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1058, dvrr_stack+172, dvrr_stack+1030, dvrr_stack+99, dvrr_stack+151, NULL);

 /* compute (0 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1094, dvrr_stack+721, dvrr_stack+1058, dvrr_stack+246, dvrr_stack+172, NULL);

 /* compute (0 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+1139,dvrr_stack+1094,dvrr_stack+721,1);


 /* compute (0 0 | 6 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_id(Libderiv->CD,dvrr_stack+1247,dvrr_stack+1139,dvrr_stack+757,1);


 /* compute (0 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,28,dvrr_stack+1415, dvrr_stack+1247, dvrr_stack+246);

 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,10,dvrr_stack+120, dvrr_stack+433, dvrr_stack+200);

 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,15,dvrr_stack+1499, dvrr_stack+586, dvrr_stack+210);

 /* compute (0 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,21,dvrr_stack+1544, dvrr_stack+841, dvrr_stack+225);

 /* compute (0 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,28,dvrr_stack+1607, dvrr_stack+1247, dvrr_stack+246);

 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,10,dvrr_stack+1691, dvrr_stack+433, dvrr_stack+200);

 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,15,dvrr_stack+433, dvrr_stack+586, dvrr_stack+210);

 /* compute (0 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,21,dvrr_stack+586, dvrr_stack+841, dvrr_stack+225);

 /* compute (0 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,28,dvrr_stack+841, dvrr_stack+1247, dvrr_stack+246);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+1247,dvrr_stack+200,dvrr_stack+9,1);


 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+1265, dvrr_stack+1247, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+478, dvrr_stack+388, NULL);
 tmp = dvrr_stack + 478;
 target_ptr = Libderiv->deriv_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+1271, dvrr_stack+358, NULL);
 tmp = dvrr_stack + 1271;
 target_ptr = Libderiv->deriv_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,21,dvrr_stack+1281, dvrr_stack+523, NULL);
 tmp = dvrr_stack + 1281;
 target_ptr = Libderiv->deriv_classes[0][5][11];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,28,dvrr_stack+1302, dvrr_stack+757, NULL);
 tmp = dvrr_stack + 1302;
 target_ptr = Libderiv->deriv_classes[0][6][11];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+1330, dvrr_stack+1139, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+1366, dvrr_stack+1247, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+1372, dvrr_stack+388, NULL);
 tmp = dvrr_stack + 1372;
 target_ptr = Libderiv->deriv_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+1387, dvrr_stack+358, NULL);
 tmp = dvrr_stack + 1387;
 target_ptr = Libderiv->deriv_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+925, dvrr_stack+523, NULL);
 tmp = dvrr_stack + 925;
 target_ptr = Libderiv->deriv_classes[0][5][10];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,28,dvrr_stack+1721, dvrr_stack+757, NULL);
 tmp = dvrr_stack + 1721;
 target_ptr = Libderiv->deriv_classes[0][6][10];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+1749, dvrr_stack+1139, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+1397, dvrr_stack+1247, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+1247, dvrr_stack+388, NULL);
 tmp = dvrr_stack + 1247;
 target_ptr = Libderiv->deriv_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+388, dvrr_stack+358, NULL);
 tmp = dvrr_stack + 388;
 target_ptr = Libderiv->deriv_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+946, dvrr_stack+523, NULL);
 tmp = dvrr_stack + 946;
 target_ptr = Libderiv->deriv_classes[0][5][9];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,28,dvrr_stack+523, dvrr_stack+757, NULL);
 tmp = dvrr_stack + 523;
 target_ptr = Libderiv->deriv_classes[0][6][9];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+757, dvrr_stack+1139, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+1139, dvrr_stack+200, dvrr_stack+6);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+1145, dvrr_stack+225, dvrr_stack+200);
 tmp = dvrr_stack + 1145;
 target_ptr = Libderiv->deriv_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+1160, dvrr_stack+210, dvrr_stack+9);
 tmp = dvrr_stack + 1160;
 target_ptr = Libderiv->deriv_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+1170, dvrr_stack+246, dvrr_stack+210);
 tmp = dvrr_stack + 1170;
 target_ptr = Libderiv->deriv_classes[0][5][8];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,1,1,dvrr_stack+1191, dvrr_stack+721, dvrr_stack+225);
 tmp = dvrr_stack + 1191;
 target_ptr = Libderiv->deriv_classes[0][6][8];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,1,1,dvrr_stack+793, dvrr_stack+1094, dvrr_stack+246);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+1219, dvrr_stack+200, dvrr_stack+6);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+1225, dvrr_stack+225, dvrr_stack+200);
 tmp = dvrr_stack + 1225;
 target_ptr = Libderiv->deriv_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+829, dvrr_stack+210, dvrr_stack+9);
 tmp = dvrr_stack + 829;
 target_ptr = Libderiv->deriv_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+551, dvrr_stack+246, dvrr_stack+210);
 tmp = dvrr_stack + 551;
 target_ptr = Libderiv->deriv_classes[0][5][7];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+358, dvrr_stack+721, dvrr_stack+225);
 tmp = dvrr_stack + 358;
 target_ptr = Libderiv->deriv_classes[0][6][7];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,1,1,dvrr_stack+1785, dvrr_stack+1094, dvrr_stack+246);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+1240, dvrr_stack+200, dvrr_stack+6);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+398, dvrr_stack+225, dvrr_stack+200);
 tmp = dvrr_stack + 398;
 target_ptr = Libderiv->deriv_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+572, dvrr_stack+210, dvrr_stack+9);
 tmp = dvrr_stack + 572;
 target_ptr = Libderiv->deriv_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+649, dvrr_stack+246, dvrr_stack+210);
 tmp = dvrr_stack + 649;
 target_ptr = Libderiv->deriv_classes[0][5][6];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+1821, dvrr_stack+721, dvrr_stack+225);
 tmp = dvrr_stack + 1821;
 target_ptr = Libderiv->deriv_classes[0][6][6];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,1,1,dvrr_stack+1849, dvrr_stack+1094, dvrr_stack+246);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1094, dvrr_stack+200, dvrr_stack+40, NULL, NULL, dvrr_stack+34);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1885, dvrr_stack+210, dvrr_stack+84, NULL, NULL, dvrr_stack+40);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+1930,dvrr_stack+1885,dvrr_stack+1094,3);


 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+2020, dvrr_stack+1930, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2050, dvrr_stack+225, dvrr_stack+99, NULL, NULL, dvrr_stack+84);

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+2113,dvrr_stack+2050,dvrr_stack+1885,3);


 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+2248, dvrr_stack+2113, NULL);

 /* compute (1 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+2293,dvrr_stack+274,dvrr_stack+2050,3);


 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,63,dvrr_stack+2482, dvrr_stack+2293, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2545, dvrr_stack+721, dvrr_stack+1058, NULL, NULL, dvrr_stack+172);

 /* compute (1 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+2653,dvrr_stack+2545,dvrr_stack+274,3);


 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,84,dvrr_stack+2905, dvrr_stack+2653, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+1058, dvrr_stack+1930, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+2989, dvrr_stack+2113, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+3034, dvrr_stack+2293, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,84,dvrr_stack+3097, dvrr_stack+2653, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+721, dvrr_stack+1930, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+1930, dvrr_stack+2113, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+2113, dvrr_stack+2293, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,84,dvrr_stack+2293, dvrr_stack+2653, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2653, dvrr_stack+9, dvrr_stack+34, NULL, NULL, dvrr_stack+31);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2671, dvrr_stack+1885, dvrr_stack+2653);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+1975, dvrr_stack+2050, dvrr_stack+1094);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+2701, dvrr_stack+274, dvrr_stack+1885);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+2764, dvrr_stack+2545, dvrr_stack+2050);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+2848, dvrr_stack+1885, dvrr_stack+2653);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+2377, dvrr_stack+2050, dvrr_stack+1094);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+2176, dvrr_stack+274, dvrr_stack+1885);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+3181, dvrr_stack+2545, dvrr_stack+2050);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2422, dvrr_stack+1885, dvrr_stack+2653);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+3265, dvrr_stack+2050, dvrr_stack+1094);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+3310, dvrr_stack+274, dvrr_stack+1885);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+3373, dvrr_stack+2545, dvrr_stack+2050);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2545, dvrr_stack+34, dvrr_stack+15, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+2452, dvrr_stack+40, dvrr_stack+21, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+2563, dvrr_stack+1094, dvrr_stack+2452, dvrr_stack+200, dvrr_stack+40, dvrr_stack+2545);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+2623, dvrr_stack+2563, dvrr_stack+200);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3457, dvrr_stack+84, dvrr_stack+69, NULL, NULL, dvrr_stack+21);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3502, dvrr_stack+1885, dvrr_stack+3457, dvrr_stack+210, dvrr_stack+84, dvrr_stack+2452);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+0, dvrr_stack+3502, dvrr_stack+210);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3592, dvrr_stack+99, dvrr_stack+151, NULL, NULL, dvrr_stack+69);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3655, dvrr_stack+2050, dvrr_stack+3592, dvrr_stack+225, dvrr_stack+99, dvrr_stack+3457);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+45, dvrr_stack+3655, dvrr_stack+225);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3781, dvrr_stack+172, dvrr_stack+1030, NULL, NULL, dvrr_stack+151);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3865, dvrr_stack+274, dvrr_stack+3781, dvrr_stack+246, dvrr_stack+172, dvrr_stack+3592);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+3781, dvrr_stack+3865, dvrr_stack+246);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+2452, dvrr_stack+2563, dvrr_stack+200);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+3457, dvrr_stack+3502, dvrr_stack+210);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+3592, dvrr_stack+3655, dvrr_stack+225);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+4033, dvrr_stack+3865, dvrr_stack+246);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+150, dvrr_stack+2563, dvrr_stack+200);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+2545, dvrr_stack+3502, dvrr_stack+210);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+3502, dvrr_stack+3655, dvrr_stack+225);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+3655, dvrr_stack+3865, dvrr_stack+246);

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+1030, dvrr_stack+274, NULL);
 tmp = dvrr_stack + 1030;
 target_ptr = Libderiv->deriv_classes[0][6][2];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+3865, dvrr_stack+274, NULL);
 tmp = dvrr_stack + 3865;
 target_ptr = Libderiv->deriv_classes[0][6][1];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+3893, dvrr_stack+274, NULL);
 tmp = dvrr_stack + 3893;
 target_ptr = Libderiv->deriv_classes[0][6][0];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,10,dvrr_stack+3921, dvrr_stack+493, NULL);
 tmp = dvrr_stack + 3921;
 target_ptr = Libderiv->deriv2_classes[0][3][143];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,15,dvrr_stack+1124, dvrr_stack+676, NULL);
 tmp = dvrr_stack + 1124;
 target_ptr = Libderiv->deriv2_classes[0][4][143];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,21,dvrr_stack+3931, dvrr_stack+967, NULL);
 tmp = dvrr_stack + 3931;
 target_ptr = Libderiv->deriv2_classes[0][5][143];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,28,dvrr_stack+3952, dvrr_stack+1415, NULL);
 tmp = dvrr_stack + 3952;
 target_ptr = Libderiv->deriv2_classes[0][6][143];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,10,dvrr_stack+3980, dvrr_stack+493, NULL);
 tmp = dvrr_stack + 3980;
 target_ptr = Libderiv->deriv2_classes[0][3][131];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,15,dvrr_stack+3990, dvrr_stack+676, NULL);
 tmp = dvrr_stack + 3990;
 target_ptr = Libderiv->deriv2_classes[0][4][131];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,21,dvrr_stack+4005, dvrr_stack+967, NULL);
 tmp = dvrr_stack + 4005;
 target_ptr = Libderiv->deriv2_classes[0][5][131];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,28,dvrr_stack+3739, dvrr_stack+1415, NULL);
 tmp = dvrr_stack + 3739;
 target_ptr = Libderiv->deriv2_classes[0][6][131];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+3767, dvrr_stack+120, NULL);
 tmp = dvrr_stack + 3767;
 target_ptr = Libderiv->deriv2_classes[0][3][130];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+3565, dvrr_stack+1499, NULL);
 tmp = dvrr_stack + 3565;
 target_ptr = Libderiv->deriv2_classes[0][4][130];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+2590, dvrr_stack+1544, NULL);
 tmp = dvrr_stack + 2590;
 target_ptr = Libderiv->deriv2_classes[0][5][130];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,28,dvrr_stack+180, dvrr_stack+1607, NULL);
 tmp = dvrr_stack + 180;
 target_ptr = Libderiv->deriv2_classes[0][6][130];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,10,dvrr_stack+3580, dvrr_stack+493, NULL);
 tmp = dvrr_stack + 3580;
 target_ptr = Libderiv->deriv2_classes[0][3][119];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,15,dvrr_stack+493, dvrr_stack+676, NULL);
 tmp = dvrr_stack + 493;
 target_ptr = Libderiv->deriv2_classes[0][4][119];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,21,dvrr_stack+2878, dvrr_stack+967, NULL);
 tmp = dvrr_stack + 2878;
 target_ptr = Libderiv->deriv2_classes[0][5][119];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,28,dvrr_stack+967, dvrr_stack+1415, NULL);
 tmp = dvrr_stack + 967;
 target_ptr = Libderiv->deriv2_classes[0][6][119];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+995, dvrr_stack+120, NULL);
 tmp = dvrr_stack + 995;
 target_ptr = Libderiv->deriv2_classes[0][3][118];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+508, dvrr_stack+1499, NULL);
 tmp = dvrr_stack + 508;
 target_ptr = Libderiv->deriv2_classes[0][4][118];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+1005, dvrr_stack+1544, NULL);
 tmp = dvrr_stack + 1005;
 target_ptr = Libderiv->deriv2_classes[0][5][118];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,28,dvrr_stack+1403, dvrr_stack+1607, NULL);
 tmp = dvrr_stack + 1403;
 target_ptr = Libderiv->deriv2_classes[0][6][118];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+2611, dvrr_stack+1691, NULL);
 tmp = dvrr_stack + 2611;
 target_ptr = Libderiv->deriv2_classes[0][3][117];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+2653, dvrr_stack+433, NULL);
 tmp = dvrr_stack + 2653;
 target_ptr = Libderiv->deriv2_classes[0][4][117];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+413, dvrr_stack+586, NULL);
 tmp = dvrr_stack + 413;
 target_ptr = Libderiv->deriv2_classes[0][5][117];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,28,dvrr_stack+582, dvrr_stack+841, NULL);
 tmp = dvrr_stack + 582;
 target_ptr = Libderiv->deriv2_classes[0][6][117];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+839, dvrr_stack+478, dvrr_stack+1265);
 tmp = dvrr_stack + 839;
 target_ptr = Libderiv->deriv2_classes[0][3][107];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+849, dvrr_stack+1281, dvrr_stack+1271);
 tmp = dvrr_stack + 849;
 target_ptr = Libderiv->deriv2_classes[0][4][107];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+864, dvrr_stack+1302, dvrr_stack+478);
 tmp = dvrr_stack + 864;
 target_ptr = Libderiv->deriv2_classes[0][5][107];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_i(Data,1,1,dvrr_stack+885, dvrr_stack+1330, dvrr_stack+1281);
 tmp = dvrr_stack + 885;
 target_ptr = Libderiv->deriv2_classes[0][6][107];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+913, dvrr_stack+1372, dvrr_stack+1366);
 tmp = dvrr_stack + 913;
 target_ptr = Libderiv->deriv2_classes[0][3][106];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+610, dvrr_stack+925, dvrr_stack+1387);
 tmp = dvrr_stack + 610;
 target_ptr = Libderiv->deriv2_classes[0][4][106];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+625, dvrr_stack+1721, dvrr_stack+1372);
 tmp = dvrr_stack + 625;
 target_ptr = Libderiv->deriv2_classes[0][5][106];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_i(Data,1,1,dvrr_stack+434, dvrr_stack+1749, dvrr_stack+925);
 tmp = dvrr_stack + 434;
 target_ptr = Libderiv->deriv2_classes[0][6][106];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+462, dvrr_stack+1247, dvrr_stack+1397);
 tmp = dvrr_stack + 462;
 target_ptr = Libderiv->deriv2_classes[0][3][105];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+1431, dvrr_stack+946, dvrr_stack+388);
 tmp = dvrr_stack + 1431;
 target_ptr = Libderiv->deriv2_classes[0][4][105];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+1446, dvrr_stack+523, dvrr_stack+1247);
 tmp = dvrr_stack + 1446;
 target_ptr = Libderiv->deriv2_classes[0][5][105];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_i(Data,1,1,dvrr_stack+1467, dvrr_stack+757, dvrr_stack+946);
 tmp = dvrr_stack + 1467;
 target_ptr = Libderiv->deriv2_classes[0][6][105];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+1495, dvrr_stack+1145, dvrr_stack+1139);
 tmp = dvrr_stack + 1495;
 target_ptr = Libderiv->deriv2_classes[0][3][104];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+1505, dvrr_stack+1170, dvrr_stack+1160);
 tmp = dvrr_stack + 1505;
 target_ptr = Libderiv->deriv2_classes[0][4][104];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+1520, dvrr_stack+1191, dvrr_stack+1145);
 tmp = dvrr_stack + 1520;
 target_ptr = Libderiv->deriv2_classes[0][5][104];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_i(Data,1,1,dvrr_stack+1541, dvrr_stack+793, dvrr_stack+1170);
 tmp = dvrr_stack + 1541;
 target_ptr = Libderiv->deriv2_classes[0][6][104];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+1569, dvrr_stack+478, dvrr_stack+1265);
 tmp = dvrr_stack + 1569;
 target_ptr = Libderiv->deriv2_classes[0][3][95];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+1579, dvrr_stack+1281, dvrr_stack+1271);
 tmp = dvrr_stack + 1579;
 target_ptr = Libderiv->deriv2_classes[0][4][95];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+1594, dvrr_stack+1302, dvrr_stack+478);
 tmp = dvrr_stack + 1594;
 target_ptr = Libderiv->deriv2_classes[0][5][95];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+1615, dvrr_stack+1330, dvrr_stack+1281);
 tmp = dvrr_stack + 1615;
 target_ptr = Libderiv->deriv2_classes[0][6][95];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+1643, dvrr_stack+1372, dvrr_stack+1366);
 tmp = dvrr_stack + 1643;
 target_ptr = Libderiv->deriv2_classes[0][3][94];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+1653, dvrr_stack+925, dvrr_stack+1387);
 tmp = dvrr_stack + 1653;
 target_ptr = Libderiv->deriv2_classes[0][4][94];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+1668, dvrr_stack+1721, dvrr_stack+1372);
 tmp = dvrr_stack + 1668;
 target_ptr = Libderiv->deriv2_classes[0][5][94];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+1689, dvrr_stack+1749, dvrr_stack+925);
 tmp = dvrr_stack + 1689;
 target_ptr = Libderiv->deriv2_classes[0][6][94];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+670, dvrr_stack+1247, dvrr_stack+1397);
 tmp = dvrr_stack + 670;
 target_ptr = Libderiv->deriv2_classes[0][3][93];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+680, dvrr_stack+946, dvrr_stack+388);
 tmp = dvrr_stack + 680;
 target_ptr = Libderiv->deriv2_classes[0][4][93];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+695, dvrr_stack+523, dvrr_stack+1247);
 tmp = dvrr_stack + 695;
 target_ptr = Libderiv->deriv2_classes[0][5][93];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+208, dvrr_stack+757, dvrr_stack+946);
 tmp = dvrr_stack + 208;
 target_ptr = Libderiv->deriv2_classes[0][6][93];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+236, dvrr_stack+1145, dvrr_stack+1139);
 tmp = dvrr_stack + 236;
 target_ptr = Libderiv->deriv2_classes[0][3][92];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+246, dvrr_stack+1170, dvrr_stack+1160);
 tmp = dvrr_stack + 246;
 target_ptr = Libderiv->deriv2_classes[0][4][92];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+261, dvrr_stack+1191, dvrr_stack+1145);
 tmp = dvrr_stack + 261;
 target_ptr = Libderiv->deriv2_classes[0][5][92];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+282, dvrr_stack+793, dvrr_stack+1170);
 tmp = dvrr_stack + 282;
 target_ptr = Libderiv->deriv2_classes[0][6][92];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+310, dvrr_stack+1225, dvrr_stack+1219);
 tmp = dvrr_stack + 310;
 target_ptr = Libderiv->deriv2_classes[0][3][91];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+320, dvrr_stack+551, dvrr_stack+829);
 tmp = dvrr_stack + 320;
 target_ptr = Libderiv->deriv2_classes[0][4][91];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+335, dvrr_stack+358, dvrr_stack+1225);
 tmp = dvrr_stack + 335;
 target_ptr = Libderiv->deriv2_classes[0][5][91];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+108, dvrr_stack+1785, dvrr_stack+551);
 tmp = dvrr_stack + 108;
 target_ptr = Libderiv->deriv2_classes[0][6][91];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+136, dvrr_stack+478, dvrr_stack+1265);
 tmp = dvrr_stack + 136;
 target_ptr = Libderiv->deriv2_classes[0][3][83];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+4117, dvrr_stack+1281, dvrr_stack+1271);
 tmp = dvrr_stack + 4117;
 target_ptr = Libderiv->deriv2_classes[0][4][83];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+4132, dvrr_stack+1302, dvrr_stack+478);
 tmp = dvrr_stack + 4132;
 target_ptr = Libderiv->deriv2_classes[0][5][83];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+1302, dvrr_stack+1330, dvrr_stack+1281);
 tmp = dvrr_stack + 1302;
 target_ptr = Libderiv->deriv2_classes[0][6][83];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1330, dvrr_stack+1372, dvrr_stack+1366);
 tmp = dvrr_stack + 1330;
 target_ptr = Libderiv->deriv2_classes[0][3][82];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+1340, dvrr_stack+925, dvrr_stack+1387);
 tmp = dvrr_stack + 1340;
 target_ptr = Libderiv->deriv2_classes[0][4][82];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+472, dvrr_stack+1721, dvrr_stack+1372);
 tmp = dvrr_stack + 472;
 target_ptr = Libderiv->deriv2_classes[0][5][82];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+1355, dvrr_stack+1749, dvrr_stack+925);
 tmp = dvrr_stack + 1355;
 target_ptr = Libderiv->deriv2_classes[0][6][82];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1383, dvrr_stack+1247, dvrr_stack+1397);
 tmp = dvrr_stack + 1383;
 target_ptr = Libderiv->deriv2_classes[0][3][81];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+923, dvrr_stack+946, dvrr_stack+388);
 tmp = dvrr_stack + 923;
 target_ptr = Libderiv->deriv2_classes[0][4][81];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+1262, dvrr_stack+523, dvrr_stack+1247);
 tmp = dvrr_stack + 1262;
 target_ptr = Libderiv->deriv2_classes[0][5][81];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+523, dvrr_stack+757, dvrr_stack+946);
 tmp = dvrr_stack + 523;
 target_ptr = Libderiv->deriv2_classes[0][6][81];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1393, dvrr_stack+1145, dvrr_stack+1139);
 tmp = dvrr_stack + 1393;
 target_ptr = Libderiv->deriv2_classes[0][3][80];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+751, dvrr_stack+1170, dvrr_stack+1160);
 tmp = dvrr_stack + 751;
 target_ptr = Libderiv->deriv2_classes[0][4][80];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+766, dvrr_stack+1191, dvrr_stack+1145);
 tmp = dvrr_stack + 766;
 target_ptr = Libderiv->deriv2_classes[0][5][80];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+1191, dvrr_stack+793, dvrr_stack+1170);
 tmp = dvrr_stack + 1191;
 target_ptr = Libderiv->deriv2_classes[0][6][80];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1139, dvrr_stack+1225, dvrr_stack+1219);
 tmp = dvrr_stack + 1139;
 target_ptr = Libderiv->deriv2_classes[0][3][79];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+1149, dvrr_stack+551, dvrr_stack+829);
 tmp = dvrr_stack + 1149;
 target_ptr = Libderiv->deriv2_classes[0][4][79];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+1164, dvrr_stack+358, dvrr_stack+1225);
 tmp = dvrr_stack + 1164;
 target_ptr = Libderiv->deriv2_classes[0][5][79];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+787, dvrr_stack+1785, dvrr_stack+551);
 tmp = dvrr_stack + 787;
 target_ptr = Libderiv->deriv2_classes[0][6][79];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+551, dvrr_stack+398, dvrr_stack+1240);
 tmp = dvrr_stack + 551;
 target_ptr = Libderiv->deriv2_classes[0][3][78];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+1219, dvrr_stack+649, dvrr_stack+572);
 tmp = dvrr_stack + 1219;
 target_ptr = Libderiv->deriv2_classes[0][4][78];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+561, dvrr_stack+1821, dvrr_stack+398);
 tmp = dvrr_stack + 561;
 target_ptr = Libderiv->deriv2_classes[0][5][78];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+1234, dvrr_stack+1849, dvrr_stack+649);
 tmp = dvrr_stack + 1234;
 target_ptr = Libderiv->deriv2_classes[0][6][78];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,10,dvrr_stack+815, dvrr_stack+2020, NULL);
 tmp = dvrr_stack + 815;
 target_ptr = Libderiv->deriv2_classes[0][3][35];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,15,dvrr_stack+938, dvrr_stack+2248, NULL);
 tmp = dvrr_stack + 938;
 target_ptr = Libderiv->deriv2_classes[0][4][35];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,21,dvrr_stack+646, dvrr_stack+2482, NULL);
 tmp = dvrr_stack + 646;
 target_ptr = Libderiv->deriv2_classes[0][5][35];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,28,dvrr_stack+1717, dvrr_stack+2905, NULL);
 tmp = dvrr_stack + 1717;
 target_ptr = Libderiv->deriv2_classes[0][6][35];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+825, dvrr_stack+1058, NULL);
 tmp = dvrr_stack + 825;
 target_ptr = Libderiv->deriv2_classes[0][3][34];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1283, dvrr_stack+2989, NULL);
 tmp = dvrr_stack + 1283;
 target_ptr = Libderiv->deriv2_classes[0][4][34];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+1745, dvrr_stack+3034, NULL);
 tmp = dvrr_stack + 1745;
 target_ptr = Libderiv->deriv2_classes[0][5][34];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+1766, dvrr_stack+3097, NULL);
 tmp = dvrr_stack + 1766;
 target_ptr = Libderiv->deriv2_classes[0][6][34];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+953, dvrr_stack+721, NULL);
 tmp = dvrr_stack + 953;
 target_ptr = Libderiv->deriv2_classes[0][3][33];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1794, dvrr_stack+1930, NULL);
 tmp = dvrr_stack + 1794;
 target_ptr = Libderiv->deriv2_classes[0][4][33];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+1809, dvrr_stack+2113, NULL);
 tmp = dvrr_stack + 1809;
 target_ptr = Libderiv->deriv2_classes[0][5][33];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+1830, dvrr_stack+2293, NULL);
 tmp = dvrr_stack + 1830;
 target_ptr = Libderiv->deriv2_classes[0][6][33];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1858, dvrr_stack+2671, NULL);
 tmp = dvrr_stack + 1858;
 target_ptr = Libderiv->deriv2_classes[0][3][32];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1868, dvrr_stack+1975, NULL);
 tmp = dvrr_stack + 1868;
 target_ptr = Libderiv->deriv2_classes[0][4][32];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+356, dvrr_stack+2701, NULL);
 tmp = dvrr_stack + 356;
 target_ptr = Libderiv->deriv2_classes[0][5][32];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+377, dvrr_stack+2764, NULL);
 tmp = dvrr_stack + 377;
 target_ptr = Libderiv->deriv2_classes[0][6][32];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+4153, dvrr_stack+2848, NULL);
 tmp = dvrr_stack + 4153;
 target_ptr = Libderiv->deriv2_classes[0][3][31];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+4163, dvrr_stack+2377, NULL);
 tmp = dvrr_stack + 4163;
 target_ptr = Libderiv->deriv2_classes[0][4][31];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+4178, dvrr_stack+2176, NULL);
 tmp = dvrr_stack + 4178;
 target_ptr = Libderiv->deriv2_classes[0][5][31];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+4199, dvrr_stack+3181, NULL);
 tmp = dvrr_stack + 4199;
 target_ptr = Libderiv->deriv2_classes[0][6][31];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+4227, dvrr_stack+1094, NULL);
 tmp = dvrr_stack + 4227;
 target_ptr = Libderiv->deriv_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+4237, dvrr_stack+2422, NULL);
 tmp = dvrr_stack + 4237;
 target_ptr = Libderiv->deriv2_classes[0][3][30];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+4247, dvrr_stack+1885, NULL);
 tmp = dvrr_stack + 4247;
 target_ptr = Libderiv->deriv_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+4262, dvrr_stack+3265, NULL);
 tmp = dvrr_stack + 4262;
 target_ptr = Libderiv->deriv2_classes[0][4][30];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+4277, dvrr_stack+2050, NULL);
 tmp = dvrr_stack + 4277;
 target_ptr = Libderiv->deriv_classes[0][5][2];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+4298, dvrr_stack+3310, NULL);
 tmp = dvrr_stack + 4298;
 target_ptr = Libderiv->deriv2_classes[0][5][30];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+4319, dvrr_stack+3373, NULL);
 tmp = dvrr_stack + 4319;
 target_ptr = Libderiv->deriv2_classes[0][6][30];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+4347, dvrr_stack+2623, NULL);
 tmp = dvrr_stack + 4347;
 target_ptr = Libderiv->deriv2_classes[0][3][26];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+4357, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 4357;
 target_ptr = Libderiv->deriv2_classes[0][4][26];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+4372, dvrr_stack+45, NULL);
 tmp = dvrr_stack + 4372;
 target_ptr = Libderiv->deriv2_classes[0][5][26];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+4393, dvrr_stack+3781, NULL);
 tmp = dvrr_stack + 4393;
 target_ptr = Libderiv->deriv2_classes[0][6][26];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,10,dvrr_stack+4421, dvrr_stack+2020, NULL);
 tmp = dvrr_stack + 4421;
 target_ptr = Libderiv->deriv2_classes[0][3][23];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,15,dvrr_stack+4431, dvrr_stack+2248, NULL);
 tmp = dvrr_stack + 4431;
 target_ptr = Libderiv->deriv2_classes[0][4][23];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,21,dvrr_stack+4446, dvrr_stack+2482, NULL);
 tmp = dvrr_stack + 4446;
 target_ptr = Libderiv->deriv2_classes[0][5][23];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,28,dvrr_stack+4467, dvrr_stack+2905, NULL);
 tmp = dvrr_stack + 4467;
 target_ptr = Libderiv->deriv2_classes[0][6][23];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4495, dvrr_stack+1058, NULL);
 tmp = dvrr_stack + 4495;
 target_ptr = Libderiv->deriv2_classes[0][3][22];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+4505, dvrr_stack+2989, NULL);
 tmp = dvrr_stack + 4505;
 target_ptr = Libderiv->deriv2_classes[0][4][22];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+4520, dvrr_stack+3034, NULL);
 tmp = dvrr_stack + 4520;
 target_ptr = Libderiv->deriv2_classes[0][5][22];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+4541, dvrr_stack+3097, NULL);
 tmp = dvrr_stack + 4541;
 target_ptr = Libderiv->deriv2_classes[0][6][22];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4569, dvrr_stack+721, NULL);
 tmp = dvrr_stack + 4569;
 target_ptr = Libderiv->deriv2_classes[0][3][21];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+4579, dvrr_stack+1930, NULL);
 tmp = dvrr_stack + 4579;
 target_ptr = Libderiv->deriv2_classes[0][4][21];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+4594, dvrr_stack+2113, NULL);
 tmp = dvrr_stack + 4594;
 target_ptr = Libderiv->deriv2_classes[0][5][21];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+4615, dvrr_stack+2293, NULL);
 tmp = dvrr_stack + 4615;
 target_ptr = Libderiv->deriv2_classes[0][6][21];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4643, dvrr_stack+2671, NULL);
 tmp = dvrr_stack + 4643;
 target_ptr = Libderiv->deriv2_classes[0][3][20];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+4653, dvrr_stack+1975, NULL);
 tmp = dvrr_stack + 4653;
 target_ptr = Libderiv->deriv2_classes[0][4][20];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+4668, dvrr_stack+2701, NULL);
 tmp = dvrr_stack + 4668;
 target_ptr = Libderiv->deriv2_classes[0][5][20];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+4689, dvrr_stack+2764, NULL);
 tmp = dvrr_stack + 4689;
 target_ptr = Libderiv->deriv2_classes[0][6][20];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4717, dvrr_stack+2848, NULL);
 tmp = dvrr_stack + 4717;
 target_ptr = Libderiv->deriv2_classes[0][3][19];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+4727, dvrr_stack+2377, NULL);
 tmp = dvrr_stack + 4727;
 target_ptr = Libderiv->deriv2_classes[0][4][19];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+4742, dvrr_stack+2176, NULL);
 tmp = dvrr_stack + 4742;
 target_ptr = Libderiv->deriv2_classes[0][5][19];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+4763, dvrr_stack+3181, NULL);
 tmp = dvrr_stack + 4763;
 target_ptr = Libderiv->deriv2_classes[0][6][19];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4791, dvrr_stack+1094, NULL);
 tmp = dvrr_stack + 4791;
 target_ptr = Libderiv->deriv_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4801, dvrr_stack+2422, NULL);
 tmp = dvrr_stack + 4801;
 target_ptr = Libderiv->deriv2_classes[0][3][18];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+4811, dvrr_stack+1885, NULL);
 tmp = dvrr_stack + 4811;
 target_ptr = Libderiv->deriv_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+4826, dvrr_stack+3265, NULL);
 tmp = dvrr_stack + 4826;
 target_ptr = Libderiv->deriv2_classes[0][4][18];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+4841, dvrr_stack+2050, NULL);
 tmp = dvrr_stack + 4841;
 target_ptr = Libderiv->deriv_classes[0][5][1];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+4862, dvrr_stack+3310, NULL);
 tmp = dvrr_stack + 4862;
 target_ptr = Libderiv->deriv2_classes[0][5][18];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+4883, dvrr_stack+3373, NULL);
 tmp = dvrr_stack + 4883;
 target_ptr = Libderiv->deriv2_classes[0][6][18];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4911, dvrr_stack+2623, NULL);
 tmp = dvrr_stack + 4911;
 target_ptr = Libderiv->deriv2_classes[0][3][14];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+4921, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 4921;
 target_ptr = Libderiv->deriv2_classes[0][4][14];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+4936, dvrr_stack+45, NULL);
 tmp = dvrr_stack + 4936;
 target_ptr = Libderiv->deriv2_classes[0][5][14];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+4957, dvrr_stack+3781, NULL);
 tmp = dvrr_stack + 4957;
 target_ptr = Libderiv->deriv2_classes[0][6][14];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4985, dvrr_stack+2452, NULL);
 tmp = dvrr_stack + 4985;
 target_ptr = Libderiv->deriv2_classes[0][3][13];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+4995, dvrr_stack+3457, NULL);
 tmp = dvrr_stack + 4995;
 target_ptr = Libderiv->deriv2_classes[0][4][13];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+5010, dvrr_stack+3592, NULL);
 tmp = dvrr_stack + 5010;
 target_ptr = Libderiv->deriv2_classes[0][5][13];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+5031, dvrr_stack+4033, NULL);
 tmp = dvrr_stack + 5031;
 target_ptr = Libderiv->deriv2_classes[0][6][13];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,10,dvrr_stack+5059, dvrr_stack+2020, NULL);
 tmp = dvrr_stack + 5059;
 target_ptr = Libderiv->deriv2_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,15,dvrr_stack+2020, dvrr_stack+2248, NULL);
 tmp = dvrr_stack + 2020;
 target_ptr = Libderiv->deriv2_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,21,dvrr_stack+2239, dvrr_stack+2482, NULL);
 tmp = dvrr_stack + 2239;
 target_ptr = Libderiv->deriv2_classes[0][5][11];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,28,dvrr_stack+2482, dvrr_stack+2905, NULL);
 tmp = dvrr_stack + 2482;
 target_ptr = Libderiv->deriv2_classes[0][6][11];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+2510, dvrr_stack+1058, NULL);
 tmp = dvrr_stack + 2510;
 target_ptr = Libderiv->deriv2_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+2035, dvrr_stack+2989, NULL);
 tmp = dvrr_stack + 2035;
 target_ptr = Libderiv->deriv2_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+1058, dvrr_stack+3034, NULL);
 tmp = dvrr_stack + 1058;
 target_ptr = Libderiv->deriv2_classes[0][5][10];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+2899, dvrr_stack+3097, NULL);
 tmp = dvrr_stack + 2899;
 target_ptr = Libderiv->deriv2_classes[0][6][10];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1079, dvrr_stack+721, NULL);
 tmp = dvrr_stack + 1079;
 target_ptr = Libderiv->deriv2_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+2520, dvrr_stack+1930, NULL);
 tmp = dvrr_stack + 2520;
 target_ptr = Libderiv->deriv2_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+1930, dvrr_stack+2113, NULL);
 tmp = dvrr_stack + 1930;
 target_ptr = Libderiv->deriv2_classes[0][5][9];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+2113, dvrr_stack+2293, NULL);
 tmp = dvrr_stack + 2113;
 target_ptr = Libderiv->deriv2_classes[0][6][9];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+2535, dvrr_stack+2671, NULL);
 tmp = dvrr_stack + 2535;
 target_ptr = Libderiv->deriv2_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+2141, dvrr_stack+1975, NULL);
 tmp = dvrr_stack + 2141;
 target_ptr = Libderiv->deriv2_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+1951, dvrr_stack+2701, NULL);
 tmp = dvrr_stack + 1951;
 target_ptr = Libderiv->deriv2_classes[0][5][8];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+1972, dvrr_stack+2764, NULL);
 tmp = dvrr_stack + 1972;
 target_ptr = Libderiv->deriv2_classes[0][6][8];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+2156, dvrr_stack+2848, NULL);
 tmp = dvrr_stack + 2156;
 target_ptr = Libderiv->deriv2_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+2000, dvrr_stack+2377, NULL);
 tmp = dvrr_stack + 2000;
 target_ptr = Libderiv->deriv2_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+2668, dvrr_stack+2176, NULL);
 tmp = dvrr_stack + 2668;
 target_ptr = Libderiv->deriv2_classes[0][5][7];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+2166, dvrr_stack+3181, NULL);
 tmp = dvrr_stack + 2166;
 target_ptr = Libderiv->deriv2_classes[0][6][7];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+2194, dvrr_stack+1094, NULL);
 tmp = dvrr_stack + 2194;
 target_ptr = Libderiv->deriv_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+2204, dvrr_stack+2422, NULL);
 tmp = dvrr_stack + 2204;
 target_ptr = Libderiv->deriv2_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+2214, dvrr_stack+1885, NULL);
 tmp = dvrr_stack + 2214;
 target_ptr = Libderiv->deriv_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1089, dvrr_stack+3265, NULL);
 tmp = dvrr_stack + 1089;
 target_ptr = Libderiv->deriv2_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+2689, dvrr_stack+2050, NULL);
 tmp = dvrr_stack + 2689;
 target_ptr = Libderiv->deriv_classes[0][5][0];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+2050, dvrr_stack+3310, NULL);
 tmp = dvrr_stack + 2050;
 target_ptr = Libderiv->deriv2_classes[0][5][6];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+2071, dvrr_stack+3373, NULL);
 tmp = dvrr_stack + 2071;
 target_ptr = Libderiv->deriv2_classes[0][6][6];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+2229, dvrr_stack+2623, NULL);
 tmp = dvrr_stack + 2229;
 target_ptr = Libderiv->deriv2_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1104, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 1104;
 target_ptr = Libderiv->deriv2_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+0, dvrr_stack+45, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[0][5][2];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+21, dvrr_stack+3781, NULL);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv2_classes[0][6][2];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+49, dvrr_stack+2452, NULL);
 tmp = dvrr_stack + 49;
 target_ptr = Libderiv->deriv2_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+59, dvrr_stack+3457, NULL);
 tmp = dvrr_stack + 59;
 target_ptr = Libderiv->deriv2_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+74, dvrr_stack+3592, NULL);
 tmp = dvrr_stack + 74;
 target_ptr = Libderiv->deriv2_classes[0][5][1];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+3777, dvrr_stack+4033, NULL);
 tmp = dvrr_stack + 3777;
 target_ptr = Libderiv->deriv2_classes[0][6][1];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+95, dvrr_stack+150, NULL);
 tmp = dvrr_stack + 95;
 target_ptr = Libderiv->deriv2_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+4026, dvrr_stack+2545, NULL);
 tmp = dvrr_stack + 4026;
 target_ptr = Libderiv->deriv2_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+2545, dvrr_stack+3502, NULL);
 tmp = dvrr_stack + 2545;
 target_ptr = Libderiv->deriv2_classes[0][5][0];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+4041, dvrr_stack+3655, NULL);
 tmp = dvrr_stack + 4041;
 target_ptr = Libderiv->deriv2_classes[0][6][0];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];


}

