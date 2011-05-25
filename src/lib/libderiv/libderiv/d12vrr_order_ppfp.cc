#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (pp|fp) integrals */

void d12vrr_order_ppfp(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+34, dvrr_stack+6, dvrr_stack+31, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+50, dvrr_stack+40, dvrr_stack+21, NULL, NULL, dvrr_stack+6);
 tmp = dvrr_stack + 50;
 target_ptr = Libderiv->dvrr_classes[1][3];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+80, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+83, dvrr_stack+12, dvrr_stack+80, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+89, dvrr_stack+15, dvrr_stack+83, dvrr_stack+0, dvrr_stack+12, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+99, dvrr_stack+21, dvrr_stack+89, dvrr_stack+6, dvrr_stack+15, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+114, dvrr_stack+40, dvrr_stack+21, dvrr_stack+34, dvrr_stack+6, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+129, dvrr_stack+114, dvrr_stack+99, NULL, NULL, dvrr_stack+21);
 tmp = dvrr_stack + 129;
 target_ptr = Libderiv->dvrr_classes[1][4];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+174,dvrr_stack+129,dvrr_stack+50,3);


 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+264, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+267, dvrr_stack+80, dvrr_stack+264, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+273, dvrr_stack+83, dvrr_stack+267, dvrr_stack+12, dvrr_stack+80, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+283, dvrr_stack+89, dvrr_stack+273, dvrr_stack+15, dvrr_stack+83, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+298, dvrr_stack+99, dvrr_stack+283, dvrr_stack+21, dvrr_stack+89, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+319, dvrr_stack+114, dvrr_stack+99, dvrr_stack+40, dvrr_stack+21, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+340, dvrr_stack+319, dvrr_stack+298, NULL, NULL, dvrr_stack+99);

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+403,dvrr_stack+340,dvrr_stack+129,3);


 /* compute (1 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+538,dvrr_stack+403,dvrr_stack+174,3);


 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,30,dvrr_stack+718, dvrr_stack+538, dvrr_stack+50);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+808, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+811, dvrr_stack+264, dvrr_stack+808, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+817, dvrr_stack+267, dvrr_stack+811, dvrr_stack+80, dvrr_stack+264, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+827, dvrr_stack+273, dvrr_stack+817, dvrr_stack+83, dvrr_stack+267, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+842, dvrr_stack+283, dvrr_stack+827, dvrr_stack+89, dvrr_stack+273, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+863, dvrr_stack+298, dvrr_stack+842, dvrr_stack+99, dvrr_stack+283, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+891, dvrr_stack+319, dvrr_stack+298, dvrr_stack+114, dvrr_stack+99, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+919, dvrr_stack+891, dvrr_stack+863, NULL, NULL, dvrr_stack+298);

 /* compute (1 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+1003,dvrr_stack+919,dvrr_stack+340,3);


 /* compute (1 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+1192,dvrr_stack+1003,dvrr_stack+403,3);


 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,45,dvrr_stack+1462, dvrr_stack+1192, dvrr_stack+129);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1597, dvrr_stack+6, dvrr_stack+15, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1615, dvrr_stack+21, dvrr_stack+89, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1645, dvrr_stack+50, dvrr_stack+1615, dvrr_stack+40, dvrr_stack+21, dvrr_stack+1597);
 tmp = dvrr_stack + 1645;
 target_ptr = Libderiv->dvrr_classes[2][3];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1705, dvrr_stack+99, dvrr_stack+283, NULL, NULL, dvrr_stack+89);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1750, dvrr_stack+129, dvrr_stack+1705, dvrr_stack+114, dvrr_stack+99, dvrr_stack+1615);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+1840,dvrr_stack+1750,dvrr_stack+1645,6);


 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2020, dvrr_stack+298, dvrr_stack+842, NULL, NULL, dvrr_stack+283);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2083, dvrr_stack+340, dvrr_stack+2020, dvrr_stack+319, dvrr_stack+298, dvrr_stack+1705);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+2209,dvrr_stack+2083,dvrr_stack+1750,6);


 /* compute (2 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+2479,dvrr_stack+2209,dvrr_stack+1840,6);


 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,60,dvrr_stack+2839, dvrr_stack+2479, dvrr_stack+1645);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3019, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+3022, dvrr_stack+808, dvrr_stack+3019, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+3028, dvrr_stack+811, dvrr_stack+3022, dvrr_stack+264, dvrr_stack+808, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3038, dvrr_stack+817, dvrr_stack+3028, dvrr_stack+267, dvrr_stack+811, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3053, dvrr_stack+827, dvrr_stack+3038, dvrr_stack+273, dvrr_stack+817, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3074, dvrr_stack+842, dvrr_stack+3053, dvrr_stack+283, dvrr_stack+827, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3102, dvrr_stack+863, dvrr_stack+3074, NULL, NULL, dvrr_stack+842);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3186, dvrr_stack+919, dvrr_stack+3102, dvrr_stack+891, dvrr_stack+863, dvrr_stack+2020);

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+3354,dvrr_stack+3186,dvrr_stack+2083,6);


 /* compute (2 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+3732,dvrr_stack+3354,dvrr_stack+2209,6);


 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,90,dvrr_stack+4272, dvrr_stack+3732, dvrr_stack+1750);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,30,dvrr_stack+3074, dvrr_stack+538, dvrr_stack+50);

 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,45,dvrr_stack+4542, dvrr_stack+1192, dvrr_stack+129);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,60,dvrr_stack+4677, dvrr_stack+2479, dvrr_stack+1645);

 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,90,dvrr_stack+4857, dvrr_stack+3732, dvrr_stack+1750);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,30,dvrr_stack+5127, dvrr_stack+538, dvrr_stack+50);

 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,45,dvrr_stack+538, dvrr_stack+1192, dvrr_stack+129);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,60,dvrr_stack+1192, dvrr_stack+2479, dvrr_stack+1645);

 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,90,dvrr_stack+2479, dvrr_stack+3732, dvrr_stack+1750);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+3732, dvrr_stack+34, dvrr_stack+6, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+3750,dvrr_stack+50,dvrr_stack+3732,3);


 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+3804, dvrr_stack+3750, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+673, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 673;
 target_ptr = Libderiv->deriv_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+3822, dvrr_stack+174, NULL);
 tmp = dvrr_stack + 3822;
 target_ptr = Libderiv->deriv_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,63,dvrr_stack+3852, dvrr_stack+1003, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+808, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3915, dvrr_stack+3732, dvrr_stack+1597, dvrr_stack+34, dvrr_stack+6, dvrr_stack+808);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+3951,dvrr_stack+1645,dvrr_stack+3915,6);


 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+4059, dvrr_stack+3951, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+2749, dvrr_stack+2209, NULL);
 tmp = dvrr_stack + 2749;
 target_ptr = Libderiv->deriv_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+4095, dvrr_stack+1840, NULL);
 tmp = dvrr_stack + 4095;
 target_ptr = Libderiv->deriv_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+5217, dvrr_stack+3354, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+4155, dvrr_stack+3750, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+4173, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 4173;
 target_ptr = Libderiv->deriv_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+4218, dvrr_stack+174, NULL);
 tmp = dvrr_stack + 4218;
 target_ptr = Libderiv->deriv_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+1372, dvrr_stack+1003, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+863, dvrr_stack+3951, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+5343, dvrr_stack+2209, NULL);
 tmp = dvrr_stack + 5343;
 target_ptr = Libderiv->deriv_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+5433, dvrr_stack+1840, NULL);
 tmp = dvrr_stack + 5433;
 target_ptr = Libderiv->deriv_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+5493, dvrr_stack+3354, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+4248, dvrr_stack+3750, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+3750, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 3750;
 target_ptr = Libderiv->deriv_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+403, dvrr_stack+174, NULL);
 tmp = dvrr_stack + 403;
 target_ptr = Libderiv->deriv_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+174, dvrr_stack+1003, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+1003, dvrr_stack+3951, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+3951, dvrr_stack+2209, NULL);
 tmp = dvrr_stack + 3951;
 target_ptr = Libderiv->deriv_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+2209, dvrr_stack+1840, NULL);
 tmp = dvrr_stack + 2209;
 target_ptr = Libderiv->deriv_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+1840, dvrr_stack+3354, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3795, dvrr_stack+31, dvrr_stack+3, NULL, NULL, Data->F+1);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+4041, dvrr_stack+50, dvrr_stack+3795);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+3354, dvrr_stack+340, dvrr_stack+50);
 tmp = dvrr_stack + 3354;
 target_ptr = Libderiv->deriv_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+237, dvrr_stack+129, dvrr_stack+3732);
 tmp = dvrr_stack + 237;
 target_ptr = Libderiv->deriv_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+3399, dvrr_stack+919, dvrr_stack+129);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3462, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+3465, dvrr_stack+3795, dvrr_stack+808, dvrr_stack+31, dvrr_stack+3, dvrr_stack+3462);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+3483, dvrr_stack+1645, dvrr_stack+3465);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+3519, dvrr_stack+2083, dvrr_stack+1645);
 tmp = dvrr_stack + 3519;
 target_ptr = Libderiv->deriv_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+3609, dvrr_stack+1750, dvrr_stack+3915);
 tmp = dvrr_stack + 3609;
 target_ptr = Libderiv->deriv_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+2269, dvrr_stack+3186, dvrr_stack+1750);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+3669, dvrr_stack+50, dvrr_stack+3795);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+3687, dvrr_stack+340, dvrr_stack+50);
 tmp = dvrr_stack + 3687;
 target_ptr = Libderiv->deriv_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1966, dvrr_stack+129, dvrr_stack+3732);
 tmp = dvrr_stack + 1966;
 target_ptr = Libderiv->deriv_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+2395, dvrr_stack+919, dvrr_stack+129);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+1039, dvrr_stack+1645, dvrr_stack+3465);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+1075, dvrr_stack+2083, dvrr_stack+1645);
 tmp = dvrr_stack + 1075;
 target_ptr = Libderiv->deriv_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+433, dvrr_stack+1750, dvrr_stack+3915);
 tmp = dvrr_stack + 433;
 target_ptr = Libderiv->deriv_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+5619, dvrr_stack+3186, dvrr_stack+1750);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+1996, dvrr_stack+50, dvrr_stack+3795);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+493, dvrr_stack+340, dvrr_stack+50);
 tmp = dvrr_stack + 493;
 target_ptr = Libderiv->deriv_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+5745, dvrr_stack+129, dvrr_stack+3732);
 tmp = dvrr_stack + 5745;
 target_ptr = Libderiv->deriv_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+5775, dvrr_stack+919, dvrr_stack+129);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+899, dvrr_stack+1645, dvrr_stack+3465);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+5838, dvrr_stack+2083, dvrr_stack+1645);
 tmp = dvrr_stack + 5838;
 target_ptr = Libderiv->deriv_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+935, dvrr_stack+1750, dvrr_stack+3915);
 tmp = dvrr_stack + 935;
 target_ptr = Libderiv->deriv_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+5928, dvrr_stack+3186, dvrr_stack+1750);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+3164,dvrr_stack+114,dvrr_stack+40,1);


 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+3462, dvrr_stack+3164, NULL);

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+3194,dvrr_stack+319,dvrr_stack+114,1);


 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+2458, dvrr_stack+3194, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3795, dvrr_stack+0, dvrr_stack+12, NULL, NULL, Data->F+3);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1165, dvrr_stack+15, dvrr_stack+83, NULL, NULL, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3239, dvrr_stack+1597, dvrr_stack+1165, dvrr_stack+6, dvrr_stack+15, dvrr_stack+3795);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+3275, dvrr_stack+89, dvrr_stack+273, NULL, NULL, dvrr_stack+83);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+6054, dvrr_stack+1615, dvrr_stack+3275, dvrr_stack+21, dvrr_stack+89, dvrr_stack+1165);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+6114, dvrr_stack+1645, dvrr_stack+6054, dvrr_stack+50, dvrr_stack+1615, dvrr_stack+3239);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3305, dvrr_stack+283, dvrr_stack+827, NULL, NULL, dvrr_stack+273);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6214, dvrr_stack+1705, dvrr_stack+3305, dvrr_stack+99, dvrr_stack+283, dvrr_stack+3275);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6304, dvrr_stack+1750, dvrr_stack+6214, dvrr_stack+129, dvrr_stack+1705, dvrr_stack+6054);

 /* compute (3 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+6454,dvrr_stack+6304,dvrr_stack+6114,10);


 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+6754, dvrr_stack+6454, NULL);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6854, dvrr_stack+842, dvrr_stack+3053, NULL, NULL, dvrr_stack+827);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6917, dvrr_stack+2020, dvrr_stack+6854, dvrr_stack+298, dvrr_stack+842, dvrr_stack+3305);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+7043, dvrr_stack+2083, dvrr_stack+6917, dvrr_stack+340, dvrr_stack+2020, dvrr_stack+6214);

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+7253,dvrr_stack+7043,dvrr_stack+6304,10);


 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+6854, dvrr_stack+7253, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+340, dvrr_stack+3164, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+99, dvrr_stack+3194, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+2014, dvrr_stack+6454, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+7703, dvrr_stack+7253, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+350, dvrr_stack+3164, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+3164, dvrr_stack+3194, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+7853, dvrr_stack+6454, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+6454, dvrr_stack+7253, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+7253, dvrr_stack+114, dvrr_stack+34);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+7263, dvrr_stack+319, dvrr_stack+40);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+7278, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+7281, dvrr_stack+808, dvrr_stack+3795, dvrr_stack+3, dvrr_stack+0, dvrr_stack+7278);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+3179, dvrr_stack+3915, dvrr_stack+3239, dvrr_stack+3732, dvrr_stack+1597, dvrr_stack+7281);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+7278, dvrr_stack+6304, dvrr_stack+3179);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+6604, dvrr_stack+7043, dvrr_stack+6114);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+3732, dvrr_stack+114, dvrr_stack+34);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+3915, dvrr_stack+319, dvrr_stack+40);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+7378, dvrr_stack+6304, dvrr_stack+3179);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+7478, dvrr_stack+7043, dvrr_stack+6114);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+3930, dvrr_stack+114, dvrr_stack+34);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+7628, dvrr_stack+319, dvrr_stack+40);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+7953, dvrr_stack+6304, dvrr_stack+3179);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+8053, dvrr_stack+7043, dvrr_stack+6114);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+3179, dvrr_stack+50, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+7643, dvrr_stack+6114, dvrr_stack+50);
 tmp = dvrr_stack + 7643;
 target_ptr = Libderiv->deriv_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+3189, dvrr_stack+129, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+7004, dvrr_stack+6304, dvrr_stack+129);
 tmp = dvrr_stack + 7004;
 target_ptr = Libderiv->deriv_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+3204, dvrr_stack+1645, dvrr_stack+40);
 tmp = dvrr_stack + 3204;
 target_ptr = Libderiv->deriv_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3234, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3, dvrr_stack+12, dvrr_stack+80, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+360, dvrr_stack+3795, dvrr_stack+3, dvrr_stack+0, dvrr_stack+12, dvrr_stack+3234);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+378, dvrr_stack+83, dvrr_stack+267, NULL, NULL, dvrr_stack+80);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+298, dvrr_stack+1165, dvrr_stack+378, dvrr_stack+15, dvrr_stack+83, dvrr_stack+3);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+7094, dvrr_stack+3239, dvrr_stack+298, dvrr_stack+1597, dvrr_stack+1165, dvrr_stack+360);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+3234, dvrr_stack+273, dvrr_stack+817, NULL, NULL, dvrr_stack+267);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+7154, dvrr_stack+3275, dvrr_stack+3234, dvrr_stack+89, dvrr_stack+273, dvrr_stack+378);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+8203, dvrr_stack+6054, dvrr_stack+7154, dvrr_stack+1615, dvrr_stack+3275, dvrr_stack+298);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+8303, dvrr_stack+6114, dvrr_stack+8203, dvrr_stack+1645, dvrr_stack+6054, dvrr_stack+7094);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+8453, dvrr_stack+8303, dvrr_stack+1645);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+7094, dvrr_stack+1750, dvrr_stack+114);
 tmp = dvrr_stack + 7094;
 target_ptr = Libderiv->deriv_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6054, dvrr_stack+827, dvrr_stack+3038, NULL, NULL, dvrr_stack+817);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2114, dvrr_stack+3305, dvrr_stack+6054, dvrr_stack+283, dvrr_stack+827, dvrr_stack+3234);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8553, dvrr_stack+6214, dvrr_stack+2114, dvrr_stack+1705, dvrr_stack+3305, dvrr_stack+7154);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8703, dvrr_stack+6304, dvrr_stack+8553, dvrr_stack+1750, dvrr_stack+6214, dvrr_stack+8203);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+8553, dvrr_stack+8703, dvrr_stack+1750);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+8203, dvrr_stack+50, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+6054, dvrr_stack+6114, dvrr_stack+50);
 tmp = dvrr_stack + 6054;
 target_ptr = Libderiv->deriv_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+8213, dvrr_stack+129, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+6214, dvrr_stack+6304, dvrr_stack+129);
 tmp = dvrr_stack + 6214;
 target_ptr = Libderiv->deriv_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+8228, dvrr_stack+1645, dvrr_stack+40);
 tmp = dvrr_stack + 8228;
 target_ptr = Libderiv->deriv_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+3234, dvrr_stack+8303, dvrr_stack+1645);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+8258, dvrr_stack+1750, dvrr_stack+114);
 tmp = dvrr_stack + 8258;
 target_ptr = Libderiv->deriv_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+8928, dvrr_stack+8703, dvrr_stack+1750);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1705, dvrr_stack+50, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+2114, dvrr_stack+6114, dvrr_stack+50);
 tmp = dvrr_stack + 2114;
 target_ptr = Libderiv->deriv_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+50, dvrr_stack+129, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+6114, dvrr_stack+6304, dvrr_stack+129);
 tmp = dvrr_stack + 6114;
 target_ptr = Libderiv->deriv_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+129, dvrr_stack+1645, dvrr_stack+40);
 tmp = dvrr_stack + 129;
 target_ptr = Libderiv->deriv_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+6304, dvrr_stack+8303, dvrr_stack+1645);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+8303, dvrr_stack+1750, dvrr_stack+114);
 tmp = dvrr_stack + 8303;
 target_ptr = Libderiv->deriv_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+9078, dvrr_stack+8703, dvrr_stack+1750);

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,30,dvrr_stack+8703, dvrr_stack+718, NULL);
 tmp = dvrr_stack + 8703;
 target_ptr = Libderiv->deriv2_classes[1][3][143];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,45,dvrr_stack+8733, dvrr_stack+1462, NULL);
 tmp = dvrr_stack + 8733;
 target_ptr = Libderiv->deriv2_classes[1][4][143];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,60,dvrr_stack+8778, dvrr_stack+2839, NULL);
 tmp = dvrr_stack + 8778;
 target_ptr = Libderiv->deriv2_classes[2][3][143];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,90,dvrr_stack+8838, dvrr_stack+4272, NULL);
 tmp = dvrr_stack + 8838;
 target_ptr = Libderiv->deriv2_classes[2][4][143];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,30,dvrr_stack+8348, dvrr_stack+718, NULL);
 tmp = dvrr_stack + 8348;
 target_ptr = Libderiv->deriv2_classes[1][3][131];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,45,dvrr_stack+8378, dvrr_stack+1462, NULL);
 tmp = dvrr_stack + 8378;
 target_ptr = Libderiv->deriv2_classes[1][4][131];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,60,dvrr_stack+1715, dvrr_stack+2839, NULL);
 tmp = dvrr_stack + 1715;
 target_ptr = Libderiv->deriv2_classes[2][3][131];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,90,dvrr_stack+7139, dvrr_stack+4272, NULL);
 tmp = dvrr_stack + 7139;
 target_ptr = Libderiv->deriv2_classes[2][4][131];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+8423, dvrr_stack+3074, NULL);
 tmp = dvrr_stack + 8423;
 target_ptr = Libderiv->deriv2_classes[1][3][130];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+6404, dvrr_stack+4542, NULL);
 tmp = dvrr_stack + 6404;
 target_ptr = Libderiv->deriv2_classes[1][4][130];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+1775, dvrr_stack+4677, NULL);
 tmp = dvrr_stack + 1775;
 target_ptr = Libderiv->deriv2_classes[2][3][130];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+1597, dvrr_stack+4857, NULL);
 tmp = dvrr_stack + 1597;
 target_ptr = Libderiv->deriv2_classes[2][4][130];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,30,dvrr_stack+65, dvrr_stack+718, NULL);
 tmp = dvrr_stack + 65;
 target_ptr = Libderiv->deriv2_classes[1][3][119];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,45,dvrr_stack+718, dvrr_stack+1462, NULL);
 tmp = dvrr_stack + 718;
 target_ptr = Libderiv->deriv2_classes[1][4][119];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,60,dvrr_stack+763, dvrr_stack+2839, NULL);
 tmp = dvrr_stack + 763;
 target_ptr = Libderiv->deriv2_classes[2][3][119];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,90,dvrr_stack+2839, dvrr_stack+4272, NULL);
 tmp = dvrr_stack + 2839;
 target_ptr = Libderiv->deriv2_classes[2][4][119];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+2929, dvrr_stack+3074, NULL);
 tmp = dvrr_stack + 2929;
 target_ptr = Libderiv->deriv2_classes[1][3][118];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+2959, dvrr_stack+4542, NULL);
 tmp = dvrr_stack + 2959;
 target_ptr = Libderiv->deriv2_classes[1][4][118];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+3004, dvrr_stack+4677, NULL);
 tmp = dvrr_stack + 3004;
 target_ptr = Libderiv->deriv2_classes[2][3][118];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+3064, dvrr_stack+4857, NULL);
 tmp = dvrr_stack + 3064;
 target_ptr = Libderiv->deriv2_classes[2][4][118];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+823, dvrr_stack+5127, NULL);
 tmp = dvrr_stack + 823;
 target_ptr = Libderiv->deriv2_classes[1][3][117];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+267, dvrr_stack+538, NULL);
 tmp = dvrr_stack + 267;
 target_ptr = Libderiv->deriv2_classes[1][4][117];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+538, dvrr_stack+1192, NULL);
 tmp = dvrr_stack + 538;
 target_ptr = Libderiv->deriv2_classes[2][3][117];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+1165, dvrr_stack+2479, NULL);
 tmp = dvrr_stack + 1165;
 target_ptr = Libderiv->deriv2_classes[2][4][117];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+598, dvrr_stack+673, dvrr_stack+3804);
 tmp = dvrr_stack + 598;
 target_ptr = Libderiv->deriv2_classes[1][3][107];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+628, dvrr_stack+3852, dvrr_stack+3822);
 tmp = dvrr_stack + 628;
 target_ptr = Libderiv->deriv2_classes[1][4][107];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+1255, dvrr_stack+2749, dvrr_stack+4059);
 tmp = dvrr_stack + 1255;
 target_ptr = Libderiv->deriv2_classes[2][3][107];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+2473, dvrr_stack+5217, dvrr_stack+4095);
 tmp = dvrr_stack + 2473;
 target_ptr = Libderiv->deriv2_classes[2][4][107];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2174, dvrr_stack+4173, dvrr_stack+4155);
 tmp = dvrr_stack + 2174;
 target_ptr = Libderiv->deriv2_classes[1][3][106];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+1315, dvrr_stack+1372, dvrr_stack+4218);
 tmp = dvrr_stack + 1315;
 target_ptr = Libderiv->deriv2_classes[1][4][106];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+2563, dvrr_stack+5343, dvrr_stack+863);
 tmp = dvrr_stack + 2563;
 target_ptr = Libderiv->deriv2_classes[2][3][106];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+2623, dvrr_stack+5493, dvrr_stack+5433);
 tmp = dvrr_stack + 2623;
 target_ptr = Libderiv->deriv2_classes[2][4][106];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+360, dvrr_stack+3750, dvrr_stack+4248);
 tmp = dvrr_stack + 360;
 target_ptr = Libderiv->deriv2_classes[1][3][105];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+0, dvrr_stack+174, dvrr_stack+403);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[1][4][105];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+4266, dvrr_stack+3951, dvrr_stack+1003);
 tmp = dvrr_stack + 4266;
 target_ptr = Libderiv->deriv2_classes[2][3][105];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+4326, dvrr_stack+1840, dvrr_stack+2209);
 tmp = dvrr_stack + 4326;
 target_ptr = Libderiv->deriv2_classes[2][4][105];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2713, dvrr_stack+3354, dvrr_stack+4041);
 tmp = dvrr_stack + 2713;
 target_ptr = Libderiv->deriv2_classes[1][3][104];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+4416, dvrr_stack+3399, dvrr_stack+237);
 tmp = dvrr_stack + 4416;
 target_ptr = Libderiv->deriv2_classes[1][4][104];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+4461, dvrr_stack+3519, dvrr_stack+3483);
 tmp = dvrr_stack + 4461;
 target_ptr = Libderiv->deriv2_classes[2][3][104];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+4521, dvrr_stack+2269, dvrr_stack+3609);
 tmp = dvrr_stack + 4521;
 target_ptr = Libderiv->deriv2_classes[2][4][104];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+4611, dvrr_stack+673, dvrr_stack+3804);
 tmp = dvrr_stack + 4611;
 target_ptr = Libderiv->deriv2_classes[1][3][95];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+4641, dvrr_stack+3852, dvrr_stack+3822);
 tmp = dvrr_stack + 4641;
 target_ptr = Libderiv->deriv2_classes[1][4][95];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+4686, dvrr_stack+2749, dvrr_stack+4059);
 tmp = dvrr_stack + 4686;
 target_ptr = Libderiv->deriv2_classes[2][3][95];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+4746, dvrr_stack+5217, dvrr_stack+4095);
 tmp = dvrr_stack + 4746;
 target_ptr = Libderiv->deriv2_classes[2][4][95];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+4836, dvrr_stack+4173, dvrr_stack+4155);
 tmp = dvrr_stack + 4836;
 target_ptr = Libderiv->deriv2_classes[1][3][94];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+4866, dvrr_stack+1372, dvrr_stack+4218);
 tmp = dvrr_stack + 4866;
 target_ptr = Libderiv->deriv2_classes[1][4][94];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+4911, dvrr_stack+5343, dvrr_stack+863);
 tmp = dvrr_stack + 4911;
 target_ptr = Libderiv->deriv2_classes[2][3][94];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+4971, dvrr_stack+5493, dvrr_stack+5433);
 tmp = dvrr_stack + 4971;
 target_ptr = Libderiv->deriv2_classes[2][4][94];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+5061, dvrr_stack+3750, dvrr_stack+4248);
 tmp = dvrr_stack + 5061;
 target_ptr = Libderiv->deriv2_classes[1][3][93];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+5091, dvrr_stack+174, dvrr_stack+403);
 tmp = dvrr_stack + 5091;
 target_ptr = Libderiv->deriv2_classes[1][4][93];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+5136, dvrr_stack+3951, dvrr_stack+1003);
 tmp = dvrr_stack + 5136;
 target_ptr = Libderiv->deriv2_classes[2][3][93];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+1435, dvrr_stack+1840, dvrr_stack+2209);
 tmp = dvrr_stack + 1435;
 target_ptr = Libderiv->deriv2_classes[2][4][93];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1525, dvrr_stack+3354, dvrr_stack+4041);
 tmp = dvrr_stack + 1525;
 target_ptr = Libderiv->deriv2_classes[1][3][92];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+9228, dvrr_stack+3399, dvrr_stack+237);
 tmp = dvrr_stack + 9228;
 target_ptr = Libderiv->deriv2_classes[1][4][92];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+9273, dvrr_stack+3519, dvrr_stack+3483);
 tmp = dvrr_stack + 9273;
 target_ptr = Libderiv->deriv2_classes[2][3][92];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+9333, dvrr_stack+2269, dvrr_stack+3609);
 tmp = dvrr_stack + 9333;
 target_ptr = Libderiv->deriv2_classes[2][4][92];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1555, dvrr_stack+3687, dvrr_stack+3669);
 tmp = dvrr_stack + 1555;
 target_ptr = Libderiv->deriv2_classes[1][3][91];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+9423, dvrr_stack+2395, dvrr_stack+1966);
 tmp = dvrr_stack + 9423;
 target_ptr = Libderiv->deriv2_classes[1][4][91];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+9468, dvrr_stack+1075, dvrr_stack+1039);
 tmp = dvrr_stack + 9468;
 target_ptr = Libderiv->deriv2_classes[2][3][91];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+9528, dvrr_stack+5619, dvrr_stack+433);
 tmp = dvrr_stack + 9528;
 target_ptr = Libderiv->deriv2_classes[2][4][91];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+9618, dvrr_stack+673, dvrr_stack+3804);
 tmp = dvrr_stack + 9618;
 target_ptr = Libderiv->deriv2_classes[1][3][83];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+9648, dvrr_stack+3852, dvrr_stack+3822);
 tmp = dvrr_stack + 9648;
 target_ptr = Libderiv->deriv2_classes[1][4][83];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+3852, dvrr_stack+2749, dvrr_stack+4059);
 tmp = dvrr_stack + 3852;
 target_ptr = Libderiv->deriv2_classes[2][3][83];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+9693, dvrr_stack+5217, dvrr_stack+4095);
 tmp = dvrr_stack + 9693;
 target_ptr = Libderiv->deriv2_classes[2][4][83];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+4059, dvrr_stack+4173, dvrr_stack+4155);
 tmp = dvrr_stack + 4059;
 target_ptr = Libderiv->deriv2_classes[1][3][82];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+5196, dvrr_stack+1372, dvrr_stack+4218);
 tmp = dvrr_stack + 5196;
 target_ptr = Libderiv->deriv2_classes[1][4][82];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+1360, dvrr_stack+5343, dvrr_stack+863);
 tmp = dvrr_stack + 1360;
 target_ptr = Libderiv->deriv2_classes[2][3][82];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+5241, dvrr_stack+5493, dvrr_stack+5433);
 tmp = dvrr_stack + 5241;
 target_ptr = Libderiv->deriv2_classes[2][4][82];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+5493, dvrr_stack+3750, dvrr_stack+4248);
 tmp = dvrr_stack + 5493;
 target_ptr = Libderiv->deriv2_classes[1][3][81];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+5523, dvrr_stack+174, dvrr_stack+403);
 tmp = dvrr_stack + 5523;
 target_ptr = Libderiv->deriv2_classes[1][4][81];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+159, dvrr_stack+3951, dvrr_stack+1003);
 tmp = dvrr_stack + 159;
 target_ptr = Libderiv->deriv2_classes[2][3][81];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+9783, dvrr_stack+1840, dvrr_stack+2209);
 tmp = dvrr_stack + 9783;
 target_ptr = Libderiv->deriv2_classes[2][4][81];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+5568, dvrr_stack+3354, dvrr_stack+4041);
 tmp = dvrr_stack + 5568;
 target_ptr = Libderiv->deriv2_classes[1][3][80];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+853, dvrr_stack+3399, dvrr_stack+237);
 tmp = dvrr_stack + 853;
 target_ptr = Libderiv->deriv2_classes[1][4][80];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+3399, dvrr_stack+3519, dvrr_stack+3483);
 tmp = dvrr_stack + 3399;
 target_ptr = Libderiv->deriv2_classes[2][3][80];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+1835, dvrr_stack+2269, dvrr_stack+3609);
 tmp = dvrr_stack + 1835;
 target_ptr = Libderiv->deriv2_classes[2][4][80];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2269, dvrr_stack+3687, dvrr_stack+3669);
 tmp = dvrr_stack + 2269;
 target_ptr = Libderiv->deriv2_classes[1][3][79];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+2299, dvrr_stack+2395, dvrr_stack+1966);
 tmp = dvrr_stack + 2299;
 target_ptr = Libderiv->deriv2_classes[1][4][79];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+2344, dvrr_stack+1075, dvrr_stack+1039);
 tmp = dvrr_stack + 2344;
 target_ptr = Libderiv->deriv2_classes[2][3][79];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+9873, dvrr_stack+5619, dvrr_stack+433);
 tmp = dvrr_stack + 9873;
 target_ptr = Libderiv->deriv2_classes[2][4][79];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2404, dvrr_stack+493, dvrr_stack+1996);
 tmp = dvrr_stack + 2404;
 target_ptr = Libderiv->deriv2_classes[1][3][78];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+5598, dvrr_stack+5775, dvrr_stack+5745);
 tmp = dvrr_stack + 5598;
 target_ptr = Libderiv->deriv2_classes[1][4][78];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+5775, dvrr_stack+5838, dvrr_stack+899);
 tmp = dvrr_stack + 5775;
 target_ptr = Libderiv->deriv2_classes[2][3][78];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+5643, dvrr_stack+5928, dvrr_stack+935);
 tmp = dvrr_stack + 5643;
 target_ptr = Libderiv->deriv2_classes[2][4][78];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,10,dvrr_stack+5928, dvrr_stack+4095, dvrr_stack+3462);
 tmp = dvrr_stack + 5928;
 target_ptr = Libderiv->deriv2_classes[1][3][35];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,15,dvrr_stack+5958, dvrr_stack+2749, dvrr_stack+2458);
 tmp = dvrr_stack + 5958;
 target_ptr = Libderiv->deriv2_classes[1][4][35];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,10,dvrr_stack+995, dvrr_stack+6754, dvrr_stack+3822);
 tmp = dvrr_stack + 995;
 target_ptr = Libderiv->deriv2_classes[2][3][35];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,15,dvrr_stack+9963, dvrr_stack+6854, dvrr_stack+673);
 tmp = dvrr_stack + 9963;
 target_ptr = Libderiv->deriv2_classes[2][4][35];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+6003, dvrr_stack+5433, dvrr_stack+340);
 tmp = dvrr_stack + 6003;
 target_ptr = Libderiv->deriv2_classes[1][3][34];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+3472, dvrr_stack+5343, dvrr_stack+99);
 tmp = dvrr_stack + 3472;
 target_ptr = Libderiv->deriv2_classes[1][4][34];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+10053, dvrr_stack+2014, dvrr_stack+4218);
 tmp = dvrr_stack + 10053;
 target_ptr = Libderiv->deriv2_classes[2][3][34];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+10113, dvrr_stack+7703, dvrr_stack+4173);
 tmp = dvrr_stack + 10113;
 target_ptr = Libderiv->deriv2_classes[2][4][34];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+898, dvrr_stack+2209, dvrr_stack+350);
 tmp = dvrr_stack + 898;
 target_ptr = Libderiv->deriv2_classes[1][3][33];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+10203, dvrr_stack+3951, dvrr_stack+3164);
 tmp = dvrr_stack + 10203;
 target_ptr = Libderiv->deriv2_classes[1][4][33];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+10248, dvrr_stack+7853, dvrr_stack+403);
 tmp = dvrr_stack + 10248;
 target_ptr = Libderiv->deriv2_classes[2][3][33];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+10308, dvrr_stack+6454, dvrr_stack+3750);
 tmp = dvrr_stack + 10308;
 target_ptr = Libderiv->deriv2_classes[2][4][33];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+1925, dvrr_stack+3609, dvrr_stack+7253);
 tmp = dvrr_stack + 1925;
 target_ptr = Libderiv->deriv2_classes[1][3][32];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+10398, dvrr_stack+3519, dvrr_stack+7263);
 tmp = dvrr_stack + 10398;
 target_ptr = Libderiv->deriv2_classes[1][4][32];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+10443, dvrr_stack+7278, dvrr_stack+237);
 tmp = dvrr_stack + 10443;
 target_ptr = Libderiv->deriv2_classes[2][3][32];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+10503, dvrr_stack+6604, dvrr_stack+3354);
 tmp = dvrr_stack + 10503;
 target_ptr = Libderiv->deriv2_classes[2][4][32];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+10593, dvrr_stack+433, dvrr_stack+3732);
 tmp = dvrr_stack + 10593;
 target_ptr = Libderiv->deriv2_classes[1][3][31];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+10623, dvrr_stack+1075, dvrr_stack+3915);
 tmp = dvrr_stack + 10623;
 target_ptr = Libderiv->deriv2_classes[1][4][31];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+10668, dvrr_stack+7378, dvrr_stack+1966);
 tmp = dvrr_stack + 10668;
 target_ptr = Libderiv->deriv2_classes[2][3][31];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+10728, dvrr_stack+7478, dvrr_stack+3687);
 tmp = dvrr_stack + 10728;
 target_ptr = Libderiv->deriv2_classes[2][4][31];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+10818, dvrr_stack+935, dvrr_stack+3930);
 tmp = dvrr_stack + 10818;
 target_ptr = Libderiv->deriv2_classes[1][3][30];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+10848, dvrr_stack+5838, dvrr_stack+7628);
 tmp = dvrr_stack + 10848;
 target_ptr = Libderiv->deriv2_classes[1][4][30];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+10893, dvrr_stack+7953, dvrr_stack+5745);
 tmp = dvrr_stack + 10893;
 target_ptr = Libderiv->deriv2_classes[2][3][30];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+10953, dvrr_stack+8053, dvrr_stack+493);
 tmp = dvrr_stack + 10953;
 target_ptr = Libderiv->deriv2_classes[2][4][30];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+11043, dvrr_stack+7643, dvrr_stack+3179);
 tmp = dvrr_stack + 11043;
 target_ptr = Libderiv->deriv2_classes[1][3][26];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+11073, dvrr_stack+7004, dvrr_stack+3189);
 tmp = dvrr_stack + 11073;
 target_ptr = Libderiv->deriv2_classes[1][4][26];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+11118, dvrr_stack+8453, dvrr_stack+3204);
 tmp = dvrr_stack + 11118;
 target_ptr = Libderiv->deriv2_classes[2][3][26];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+11178, dvrr_stack+8553, dvrr_stack+7094);
 tmp = dvrr_stack + 11178;
 target_ptr = Libderiv->deriv2_classes[2][4][26];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,10,dvrr_stack+11268, dvrr_stack+4095, dvrr_stack+3462);
 tmp = dvrr_stack + 11268;
 target_ptr = Libderiv->deriv2_classes[1][3][23];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,15,dvrr_stack+11298, dvrr_stack+2749, dvrr_stack+2458);
 tmp = dvrr_stack + 11298;
 target_ptr = Libderiv->deriv2_classes[1][4][23];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,10,dvrr_stack+11343, dvrr_stack+6754, dvrr_stack+3822);
 tmp = dvrr_stack + 11343;
 target_ptr = Libderiv->deriv2_classes[2][3][23];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,15,dvrr_stack+11403, dvrr_stack+6854, dvrr_stack+673);
 tmp = dvrr_stack + 11403;
 target_ptr = Libderiv->deriv2_classes[2][4][23];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+11493, dvrr_stack+5433, dvrr_stack+340);
 tmp = dvrr_stack + 11493;
 target_ptr = Libderiv->deriv2_classes[1][3][22];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+11523, dvrr_stack+5343, dvrr_stack+99);
 tmp = dvrr_stack + 11523;
 target_ptr = Libderiv->deriv2_classes[1][4][22];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+11568, dvrr_stack+2014, dvrr_stack+4218);
 tmp = dvrr_stack + 11568;
 target_ptr = Libderiv->deriv2_classes[2][3][22];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+11628, dvrr_stack+7703, dvrr_stack+4173);
 tmp = dvrr_stack + 11628;
 target_ptr = Libderiv->deriv2_classes[2][4][22];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+11718, dvrr_stack+2209, dvrr_stack+350);
 tmp = dvrr_stack + 11718;
 target_ptr = Libderiv->deriv2_classes[1][3][21];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+11748, dvrr_stack+3951, dvrr_stack+3164);
 tmp = dvrr_stack + 11748;
 target_ptr = Libderiv->deriv2_classes[1][4][21];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+11793, dvrr_stack+7853, dvrr_stack+403);
 tmp = dvrr_stack + 11793;
 target_ptr = Libderiv->deriv2_classes[2][3][21];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+11853, dvrr_stack+6454, dvrr_stack+3750);
 tmp = dvrr_stack + 11853;
 target_ptr = Libderiv->deriv2_classes[2][4][21];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+11943, dvrr_stack+3609, dvrr_stack+7253);
 tmp = dvrr_stack + 11943;
 target_ptr = Libderiv->deriv2_classes[1][3][20];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+11973, dvrr_stack+3519, dvrr_stack+7263);
 tmp = dvrr_stack + 11973;
 target_ptr = Libderiv->deriv2_classes[1][4][20];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+12018, dvrr_stack+7278, dvrr_stack+237);
 tmp = dvrr_stack + 12018;
 target_ptr = Libderiv->deriv2_classes[2][3][20];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+12078, dvrr_stack+6604, dvrr_stack+3354);
 tmp = dvrr_stack + 12078;
 target_ptr = Libderiv->deriv2_classes[2][4][20];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+12168, dvrr_stack+433, dvrr_stack+3732);
 tmp = dvrr_stack + 12168;
 target_ptr = Libderiv->deriv2_classes[1][3][19];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+12198, dvrr_stack+1075, dvrr_stack+3915);
 tmp = dvrr_stack + 12198;
 target_ptr = Libderiv->deriv2_classes[1][4][19];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+12243, dvrr_stack+7378, dvrr_stack+1966);
 tmp = dvrr_stack + 12243;
 target_ptr = Libderiv->deriv2_classes[2][3][19];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+12303, dvrr_stack+7478, dvrr_stack+3687);
 tmp = dvrr_stack + 12303;
 target_ptr = Libderiv->deriv2_classes[2][4][19];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+12393, dvrr_stack+935, dvrr_stack+3930);
 tmp = dvrr_stack + 12393;
 target_ptr = Libderiv->deriv2_classes[1][3][18];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+12423, dvrr_stack+5838, dvrr_stack+7628);
 tmp = dvrr_stack + 12423;
 target_ptr = Libderiv->deriv2_classes[1][4][18];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+12468, dvrr_stack+7953, dvrr_stack+5745);
 tmp = dvrr_stack + 12468;
 target_ptr = Libderiv->deriv2_classes[2][3][18];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+12528, dvrr_stack+8053, dvrr_stack+493);
 tmp = dvrr_stack + 12528;
 target_ptr = Libderiv->deriv2_classes[2][4][18];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+12618, dvrr_stack+7643, dvrr_stack+3179);
 tmp = dvrr_stack + 12618;
 target_ptr = Libderiv->deriv2_classes[1][3][14];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+12648, dvrr_stack+7004, dvrr_stack+3189);
 tmp = dvrr_stack + 12648;
 target_ptr = Libderiv->deriv2_classes[1][4][14];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+12693, dvrr_stack+8453, dvrr_stack+3204);
 tmp = dvrr_stack + 12693;
 target_ptr = Libderiv->deriv2_classes[2][3][14];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+12753, dvrr_stack+8553, dvrr_stack+7094);
 tmp = dvrr_stack + 12753;
 target_ptr = Libderiv->deriv2_classes[2][4][14];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+12843, dvrr_stack+6054, dvrr_stack+8203);
 tmp = dvrr_stack + 12843;
 target_ptr = Libderiv->deriv2_classes[1][3][13];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+12873, dvrr_stack+6214, dvrr_stack+8213);
 tmp = dvrr_stack + 12873;
 target_ptr = Libderiv->deriv2_classes[1][4][13];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+12918, dvrr_stack+3234, dvrr_stack+8228);
 tmp = dvrr_stack + 12918;
 target_ptr = Libderiv->deriv2_classes[2][3][13];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+12978, dvrr_stack+8928, dvrr_stack+8258);
 tmp = dvrr_stack + 12978;
 target_ptr = Libderiv->deriv2_classes[2][4][13];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,10,dvrr_stack+13068, dvrr_stack+4095, dvrr_stack+3462);
 tmp = dvrr_stack + 13068;
 target_ptr = Libderiv->deriv2_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,15,dvrr_stack+4089, dvrr_stack+2749, dvrr_stack+2458);
 tmp = dvrr_stack + 4089;
 target_ptr = Libderiv->deriv2_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,10,dvrr_stack+2743, dvrr_stack+6754, dvrr_stack+3822);
 tmp = dvrr_stack + 2743;
 target_ptr = Libderiv->deriv2_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,15,dvrr_stack+6754, dvrr_stack+6854, dvrr_stack+673);
 tmp = dvrr_stack + 6754;
 target_ptr = Libderiv->deriv2_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+673, dvrr_stack+5433, dvrr_stack+340);
 tmp = dvrr_stack + 673;
 target_ptr = Libderiv->deriv2_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+5433, dvrr_stack+5343, dvrr_stack+99);
 tmp = dvrr_stack + 5433;
 target_ptr = Libderiv->deriv2_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+6844, dvrr_stack+2014, dvrr_stack+4218);
 tmp = dvrr_stack + 6844;
 target_ptr = Libderiv->deriv2_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+6904, dvrr_stack+7703, dvrr_stack+4173);
 tmp = dvrr_stack + 6904;
 target_ptr = Libderiv->deriv2_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+7703, dvrr_stack+2209, dvrr_stack+350);
 tmp = dvrr_stack + 7703;
 target_ptr = Libderiv->deriv2_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+7733, dvrr_stack+3951, dvrr_stack+3164);
 tmp = dvrr_stack + 7733;
 target_ptr = Libderiv->deriv2_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+7778, dvrr_stack+7853, dvrr_stack+403);
 tmp = dvrr_stack + 7778;
 target_ptr = Libderiv->deriv2_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+7838, dvrr_stack+6454, dvrr_stack+3750);
 tmp = dvrr_stack + 7838;
 target_ptr = Libderiv->deriv2_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+1996, dvrr_stack+3609, dvrr_stack+7253);
 tmp = dvrr_stack + 1996;
 target_ptr = Libderiv->deriv2_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3609, dvrr_stack+3519, dvrr_stack+7263);
 tmp = dvrr_stack + 3609;
 target_ptr = Libderiv->deriv2_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+2026, dvrr_stack+7278, dvrr_stack+237);
 tmp = dvrr_stack + 2026;
 target_ptr = Libderiv->deriv2_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+4134, dvrr_stack+6604, dvrr_stack+3354);
 tmp = dvrr_stack + 4134;
 target_ptr = Libderiv->deriv2_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3654, dvrr_stack+433, dvrr_stack+3732);
 tmp = dvrr_stack + 3654;
 target_ptr = Libderiv->deriv2_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3732, dvrr_stack+1075, dvrr_stack+3915);
 tmp = dvrr_stack + 3732;
 target_ptr = Libderiv->deriv2_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+3777, dvrr_stack+7378, dvrr_stack+1966);
 tmp = dvrr_stack + 3777;
 target_ptr = Libderiv->deriv2_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+6449, dvrr_stack+7478, dvrr_stack+3687);
 tmp = dvrr_stack + 6449;
 target_ptr = Libderiv->deriv2_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3684, dvrr_stack+935, dvrr_stack+3930);
 tmp = dvrr_stack + 3684;
 target_ptr = Libderiv->deriv2_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3912, dvrr_stack+5838, dvrr_stack+7628);
 tmp = dvrr_stack + 3912;
 target_ptr = Libderiv->deriv2_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+5835, dvrr_stack+7953, dvrr_stack+5745);
 tmp = dvrr_stack + 5835;
 target_ptr = Libderiv->deriv2_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+7928, dvrr_stack+8053, dvrr_stack+493);
 tmp = dvrr_stack + 7928;
 target_ptr = Libderiv->deriv2_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+8018, dvrr_stack+7643, dvrr_stack+3179);
 tmp = dvrr_stack + 8018;
 target_ptr = Libderiv->deriv2_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+8048, dvrr_stack+7004, dvrr_stack+3189);
 tmp = dvrr_stack + 8048;
 target_ptr = Libderiv->deriv2_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+8093, dvrr_stack+8453, dvrr_stack+3204);
 tmp = dvrr_stack + 8093;
 target_ptr = Libderiv->deriv2_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+8453, dvrr_stack+8553, dvrr_stack+7094);
 tmp = dvrr_stack + 8453;
 target_ptr = Libderiv->deriv2_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+8543, dvrr_stack+6054, dvrr_stack+8203);
 tmp = dvrr_stack + 8543;
 target_ptr = Libderiv->deriv2_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+8573, dvrr_stack+6214, dvrr_stack+8213);
 tmp = dvrr_stack + 8573;
 target_ptr = Libderiv->deriv2_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+8618, dvrr_stack+3234, dvrr_stack+8228);
 tmp = dvrr_stack + 8618;
 target_ptr = Libderiv->deriv2_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+8153, dvrr_stack+8928, dvrr_stack+8258);
 tmp = dvrr_stack + 8153;
 target_ptr = Libderiv->deriv2_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+8928, dvrr_stack+2114, dvrr_stack+1705);
 tmp = dvrr_stack + 8928;
 target_ptr = Libderiv->deriv2_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+8958, dvrr_stack+6114, dvrr_stack+50);
 tmp = dvrr_stack + 8958;
 target_ptr = Libderiv->deriv2_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+8243, dvrr_stack+6304, dvrr_stack+129);
 tmp = dvrr_stack + 8243;
 target_ptr = Libderiv->deriv2_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+6994, dvrr_stack+9078, dvrr_stack+8303);
 tmp = dvrr_stack + 6994;
 target_ptr = Libderiv->deriv2_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];


}

