#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (dp|fd) integrals */

void d1vrr_order_dpfd(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+0, dvrr_stack+3, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+12, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+6, NULL, NULL, dvrr_stack+0);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+39, dvrr_stack+15, dvrr_stack+6, dvrr_stack+12, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+49, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+52, dvrr_stack+49, dvrr_stack+12, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+58, dvrr_stack+52, dvrr_stack+15, dvrr_stack+49, dvrr_stack+12, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+49, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+68, dvrr_stack+3, dvrr_stack+49, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+74, dvrr_stack+6, dvrr_stack+68, dvrr_stack+0, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+84, dvrr_stack+39, dvrr_stack+74, NULL, NULL, dvrr_stack+6);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+114, dvrr_stack+58, dvrr_stack+39, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+144, dvrr_stack+114, dvrr_stack+84, dvrr_stack+58, dvrr_stack+39, dvrr_stack+21);
 tmp = dvrr_stack + 144;
 target_ptr = Libderiv->dvrr_classes[2][3];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+204, dvrr_stack+39, dvrr_stack+74, dvrr_stack+15, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+219, dvrr_stack+58, dvrr_stack+39, dvrr_stack+52, dvrr_stack+15, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+234, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+237, dvrr_stack+49, dvrr_stack+234, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+243, dvrr_stack+68, dvrr_stack+237, dvrr_stack+3, dvrr_stack+49, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+253, dvrr_stack+74, dvrr_stack+243, dvrr_stack+6, dvrr_stack+68, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+268, dvrr_stack+204, dvrr_stack+253, NULL, NULL, dvrr_stack+74);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+313, dvrr_stack+219, dvrr_stack+204, NULL, NULL, dvrr_stack+39);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+358, dvrr_stack+313, dvrr_stack+268, dvrr_stack+219, dvrr_stack+204, dvrr_stack+84);
 tmp = dvrr_stack + 358;
 target_ptr = Libderiv->dvrr_classes[2][4];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+448,dvrr_stack+358,dvrr_stack+144,6);


 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+628, dvrr_stack+204, dvrr_stack+253, dvrr_stack+39, dvrr_stack+74, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+649, dvrr_stack+219, dvrr_stack+204, dvrr_stack+58, dvrr_stack+39, NULL);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+58, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+61, dvrr_stack+234, dvrr_stack+58, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+670, dvrr_stack+237, dvrr_stack+61, dvrr_stack+49, dvrr_stack+234, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+680, dvrr_stack+243, dvrr_stack+670, dvrr_stack+68, dvrr_stack+237, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+695, dvrr_stack+253, dvrr_stack+680, dvrr_stack+74, dvrr_stack+243, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+716, dvrr_stack+628, dvrr_stack+695, NULL, NULL, dvrr_stack+253);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+779, dvrr_stack+649, dvrr_stack+628, NULL, NULL, dvrr_stack+204);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+842, dvrr_stack+779, dvrr_stack+716, dvrr_stack+649, dvrr_stack+628, dvrr_stack+268);
 tmp = dvrr_stack + 842;
 target_ptr = Libderiv->dvrr_classes[2][5];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+968,dvrr_stack+842,dvrr_stack+358,6);


 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1238, dvrr_stack+628, dvrr_stack+695, dvrr_stack+204, dvrr_stack+253, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1266, dvrr_stack+649, dvrr_stack+628, dvrr_stack+219, dvrr_stack+204, NULL);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+219, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+222, dvrr_stack+58, dvrr_stack+219, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+649, dvrr_stack+61, dvrr_stack+222, dvrr_stack+234, dvrr_stack+58, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1294, dvrr_stack+670, dvrr_stack+649, dvrr_stack+237, dvrr_stack+61, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1309, dvrr_stack+680, dvrr_stack+1294, dvrr_stack+243, dvrr_stack+670, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1330, dvrr_stack+695, dvrr_stack+1309, dvrr_stack+253, dvrr_stack+680, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1358, dvrr_stack+1238, dvrr_stack+1330, NULL, NULL, dvrr_stack+695);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1442, dvrr_stack+1266, dvrr_stack+1238, NULL, NULL, dvrr_stack+628);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1526, dvrr_stack+1442, dvrr_stack+1358, dvrr_stack+1266, dvrr_stack+1238, dvrr_stack+716);

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+1694,dvrr_stack+1526,dvrr_stack+842,6);


 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+228, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+3);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1266, dvrr_stack+6, dvrr_stack+68, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+2072, dvrr_stack+21, dvrr_stack+1266, dvrr_stack+15, dvrr_stack+6, dvrr_stack+228);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+2108, dvrr_stack+74, dvrr_stack+243, NULL, NULL, dvrr_stack+68);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+2138, dvrr_stack+84, dvrr_stack+2108, dvrr_stack+39, dvrr_stack+74, dvrr_stack+1266);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+2198, dvrr_stack+144, dvrr_stack+2138, dvrr_stack+114, dvrr_stack+84, dvrr_stack+2072);
 tmp = dvrr_stack + 2198;
 target_ptr = Libderiv->dvrr_classes[3][3];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2298, dvrr_stack+253, dvrr_stack+680, NULL, NULL, dvrr_stack+243);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2343, dvrr_stack+268, dvrr_stack+2298, dvrr_stack+204, dvrr_stack+253, dvrr_stack+2108);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2433, dvrr_stack+358, dvrr_stack+2343, dvrr_stack+313, dvrr_stack+268, dvrr_stack+2138);
 tmp = dvrr_stack + 2433;
 target_ptr = Libderiv->dvrr_classes[3][4];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+2583,dvrr_stack+2433,dvrr_stack+2198,10);


 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2883, dvrr_stack+695, dvrr_stack+1309, NULL, NULL, dvrr_stack+680);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2946, dvrr_stack+716, dvrr_stack+2883, dvrr_stack+628, dvrr_stack+695, dvrr_stack+2298);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3072, dvrr_stack+842, dvrr_stack+2946, dvrr_stack+779, dvrr_stack+716, dvrr_stack+2343);

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+3282,dvrr_stack+3072,dvrr_stack+2433,10);


 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+628, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+631, dvrr_stack+219, dvrr_stack+628, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+39, dvrr_stack+222, dvrr_stack+631, dvrr_stack+58, dvrr_stack+219, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+628, dvrr_stack+649, dvrr_stack+39, dvrr_stack+61, dvrr_stack+222, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+204, dvrr_stack+1294, dvrr_stack+628, dvrr_stack+670, dvrr_stack+649, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+628, dvrr_stack+1309, dvrr_stack+204, dvrr_stack+680, dvrr_stack+1294, NULL);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3732, dvrr_stack+1330, dvrr_stack+628, NULL, NULL, dvrr_stack+1309);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3816, dvrr_stack+1358, dvrr_stack+3732, dvrr_stack+1238, dvrr_stack+1330, dvrr_stack+2883);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3984, dvrr_stack+1526, dvrr_stack+3816, dvrr_stack+1442, dvrr_stack+1358, dvrr_stack+2946);

 /* compute (3 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+4264,dvrr_stack+3984,dvrr_stack+3072,10);


 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1330, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1339, dvrr_stack+52, dvrr_stack+15, NULL, NULL, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1357, dvrr_stack+1339, dvrr_stack+21, dvrr_stack+52, dvrr_stack+15, dvrr_stack+1330);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+225, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+1393, dvrr_stack+1330, dvrr_stack+228, dvrr_stack+12, dvrr_stack+0, dvrr_stack+225);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+1411, dvrr_stack+1357, dvrr_stack+2072, dvrr_stack+1339, dvrr_stack+21, dvrr_stack+1393);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+225, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+12, dvrr_stack+3, dvrr_stack+49, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+1393, dvrr_stack+228, dvrr_stack+12, dvrr_stack+0, dvrr_stack+3, dvrr_stack+225);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1330, dvrr_stack+68, dvrr_stack+237, NULL, NULL, dvrr_stack+49);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1471, dvrr_stack+1266, dvrr_stack+1330, dvrr_stack+6, dvrr_stack+68, dvrr_stack+12);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+3732, dvrr_stack+2072, dvrr_stack+1471, dvrr_stack+21, dvrr_stack+1266, dvrr_stack+1393);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+2072, dvrr_stack+243, dvrr_stack+670, NULL, NULL, dvrr_stack+237);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+0, dvrr_stack+2108, dvrr_stack+2072, dvrr_stack+74, dvrr_stack+243, dvrr_stack+1330);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+3792, dvrr_stack+2138, dvrr_stack+0, dvrr_stack+84, dvrr_stack+2108, dvrr_stack+1471);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+4894, dvrr_stack+2198, dvrr_stack+3792, dvrr_stack+144, dvrr_stack+2138, dvrr_stack+3732);

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3732, dvrr_stack+680, dvrr_stack+1294, NULL, NULL, dvrr_stack+670);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2102, dvrr_stack+2298, dvrr_stack+3732, dvrr_stack+253, dvrr_stack+680, dvrr_stack+2072);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5044, dvrr_stack+2343, dvrr_stack+2102, dvrr_stack+268, dvrr_stack+2298, dvrr_stack+0);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5194, dvrr_stack+2433, dvrr_stack+5044, dvrr_stack+358, dvrr_stack+2343, dvrr_stack+3792);

 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+1309, dvrr_stack+204, NULL, NULL, dvrr_stack+1294);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2298, dvrr_stack+2883, dvrr_stack+0, dvrr_stack+695, dvrr_stack+1309, dvrr_stack+3732);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3732, dvrr_stack+2946, dvrr_stack+2298, dvrr_stack+716, dvrr_stack+2883, dvrr_stack+2102);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+5419, dvrr_stack+3072, dvrr_stack+3732, dvrr_stack+842, dvrr_stack+2946, dvrr_stack+5044);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+5044, dvrr_stack+448, NULL);
 tmp = dvrr_stack + 5044;
 target_ptr = Libderiv->deriv_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+5104, dvrr_stack+968, NULL);
 tmp = dvrr_stack + 5104;
 target_ptr = Libderiv->deriv_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+2072, dvrr_stack+1694, NULL);
 tmp = dvrr_stack + 2072;
 target_ptr = Libderiv->deriv_classes[2][5][11];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+3732, dvrr_stack+2583, NULL);
 tmp = dvrr_stack + 3732;
 target_ptr = Libderiv->deriv_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+3832, dvrr_stack+3282, NULL);
 tmp = dvrr_stack + 3832;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+5734, dvrr_stack+4264, NULL);
 tmp = dvrr_stack + 5734;
 target_ptr = Libderiv->deriv_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+2883, dvrr_stack+448, NULL);
 tmp = dvrr_stack + 2883;
 target_ptr = Libderiv->deriv_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+2943, dvrr_stack+968, NULL);
 tmp = dvrr_stack + 2943;
 target_ptr = Libderiv->deriv_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+2298, dvrr_stack+1694, NULL);
 tmp = dvrr_stack + 2298;
 target_ptr = Libderiv->deriv_classes[2][5][10];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+0, dvrr_stack+2583, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+628, dvrr_stack+3282, NULL);
 tmp = dvrr_stack + 628;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+5944, dvrr_stack+4264, NULL);
 tmp = dvrr_stack + 5944;
 target_ptr = Libderiv->deriv_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+204, dvrr_stack+448, NULL);
 tmp = dvrr_stack + 204;
 target_ptr = Libderiv->deriv_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+448, dvrr_stack+968, NULL);
 tmp = dvrr_stack + 448;
 target_ptr = Libderiv->deriv_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+968, dvrr_stack+1694, NULL);
 tmp = dvrr_stack + 968;
 target_ptr = Libderiv->deriv_classes[2][5][9];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+1694, dvrr_stack+2583, NULL);
 tmp = dvrr_stack + 1694;
 target_ptr = Libderiv->deriv_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+2583, dvrr_stack+3282, NULL);
 tmp = dvrr_stack + 2583;
 target_ptr = Libderiv->deriv_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+3282, dvrr_stack+4264, NULL);
 tmp = dvrr_stack + 3282;
 target_ptr = Libderiv->deriv_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+4264, dvrr_stack+358, dvrr_stack+1357);
 tmp = dvrr_stack + 4264;
 target_ptr = Libderiv->deriv_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+538, dvrr_stack+842, dvrr_stack+144);
 tmp = dvrr_stack + 538;
 target_ptr = Libderiv->deriv_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+4324, dvrr_stack+1526, dvrr_stack+358);
 tmp = dvrr_stack + 4324;
 target_ptr = Libderiv->deriv_classes[2][5][8];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+4450, dvrr_stack+2433, dvrr_stack+1411);
 tmp = dvrr_stack + 4450;
 target_ptr = Libderiv->deriv_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+2733, dvrr_stack+3072, dvrr_stack+2198);
 tmp = dvrr_stack + 2733;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+4550, dvrr_stack+3984, dvrr_stack+2433);
 tmp = dvrr_stack + 4550;
 target_ptr = Libderiv->deriv_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+4760, dvrr_stack+358, dvrr_stack+1357);
 tmp = dvrr_stack + 4760;
 target_ptr = Libderiv->deriv_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+3492, dvrr_stack+842, dvrr_stack+144);
 tmp = dvrr_stack + 3492;
 target_ptr = Libderiv->deriv_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+3582, dvrr_stack+1526, dvrr_stack+358);
 tmp = dvrr_stack + 3582;
 target_ptr = Libderiv->deriv_classes[2][5][7];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+1794, dvrr_stack+2433, dvrr_stack+1411);
 tmp = dvrr_stack + 1794;
 target_ptr = Libderiv->deriv_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+1894, dvrr_stack+3072, dvrr_stack+2198);
 tmp = dvrr_stack + 1894;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+1094, dvrr_stack+3984, dvrr_stack+2433);
 tmp = dvrr_stack + 1094;
 target_ptr = Libderiv->deriv_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+4820, dvrr_stack+358, dvrr_stack+1357);
 tmp = dvrr_stack + 4820;
 target_ptr = Libderiv->deriv_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+1304, dvrr_stack+842, dvrr_stack+144);
 tmp = dvrr_stack + 1304;
 target_ptr = Libderiv->deriv_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+6154, dvrr_stack+1526, dvrr_stack+358);
 tmp = dvrr_stack + 6154;
 target_ptr = Libderiv->deriv_classes[2][5][6];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+1471, dvrr_stack+2433, dvrr_stack+1411);
 tmp = dvrr_stack + 1471;
 target_ptr = Libderiv->deriv_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+6280, dvrr_stack+3072, dvrr_stack+2198);
 tmp = dvrr_stack + 6280;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+6430, dvrr_stack+3984, dvrr_stack+2433);
 tmp = dvrr_stack + 6430;
 target_ptr = Libderiv->deriv_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+1394, dvrr_stack+2198, dvrr_stack+114);
 tmp = dvrr_stack + 1394;
 target_ptr = Libderiv->deriv_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+3982, dvrr_stack+2433, dvrr_stack+313);
 tmp = dvrr_stack + 3982;
 target_ptr = Libderiv->deriv_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+4072, dvrr_stack+3072, dvrr_stack+779);
 tmp = dvrr_stack + 4072;
 target_ptr = Libderiv->deriv_classes[2][5][2];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+1571, dvrr_stack+4894, dvrr_stack+144);
 tmp = dvrr_stack + 1571;
 target_ptr = Libderiv->deriv_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+6640, dvrr_stack+5194, dvrr_stack+358);
 tmp = dvrr_stack + 6640;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+6790, dvrr_stack+5419, dvrr_stack+842);
 tmp = dvrr_stack + 6790;
 target_ptr = Libderiv->deriv_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+4198, dvrr_stack+2198, dvrr_stack+114);
 tmp = dvrr_stack + 4198;
 target_ptr = Libderiv->deriv_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+7000, dvrr_stack+2433, dvrr_stack+313);
 tmp = dvrr_stack + 7000;
 target_ptr = Libderiv->deriv_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+7090, dvrr_stack+3072, dvrr_stack+779);
 tmp = dvrr_stack + 7090;
 target_ptr = Libderiv->deriv_classes[2][5][1];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+7216, dvrr_stack+4894, dvrr_stack+144);
 tmp = dvrr_stack + 7216;
 target_ptr = Libderiv->deriv_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+7316, dvrr_stack+5194, dvrr_stack+358);
 tmp = dvrr_stack + 7316;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+7466, dvrr_stack+5419, dvrr_stack+842);
 tmp = dvrr_stack + 7466;
 target_ptr = Libderiv->deriv_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+7676, dvrr_stack+2198, dvrr_stack+114);
 tmp = dvrr_stack + 7676;
 target_ptr = Libderiv->deriv_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+2198, dvrr_stack+2433, dvrr_stack+313);
 tmp = dvrr_stack + 2198;
 target_ptr = Libderiv->deriv_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+2424, dvrr_stack+3072, dvrr_stack+779);
 tmp = dvrr_stack + 2424;
 target_ptr = Libderiv->deriv_classes[2][5][0];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+3033, dvrr_stack+4894, dvrr_stack+144);
 tmp = dvrr_stack + 3033;
 target_ptr = Libderiv->deriv_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+4880, dvrr_stack+5194, dvrr_stack+358);
 tmp = dvrr_stack + 4880;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+5194, dvrr_stack+5419, dvrr_stack+842);
 tmp = dvrr_stack + 5194;
 target_ptr = Libderiv->deriv_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];


}

