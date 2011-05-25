#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (pp|gf) integrals */

void d1vrr_order_ppgf(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+9, dvrr_stack+0, dvrr_stack+6, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+9, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+34, dvrr_stack+6, dvrr_stack+31, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+9, dvrr_stack+34, dvrr_stack+0, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+50, dvrr_stack+21, dvrr_stack+40, dvrr_stack+15, dvrr_stack+9, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+65, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+68, dvrr_stack+65, dvrr_stack+3, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+74, dvrr_stack+68, dvrr_stack+15, dvrr_stack+65, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+84, dvrr_stack+74, dvrr_stack+21, dvrr_stack+68, dvrr_stack+15, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+99, dvrr_stack+84, dvrr_stack+50, NULL, NULL, dvrr_stack+21);
 tmp = dvrr_stack + 99;
 target_ptr = Libderiv->dvrr_classes[1][4];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+65, dvrr_stack+31, dvrr_stack+3, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+144, dvrr_stack+34, dvrr_stack+65, dvrr_stack+6, dvrr_stack+31, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+154, dvrr_stack+40, dvrr_stack+144, dvrr_stack+9, dvrr_stack+34, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+169, dvrr_stack+50, dvrr_stack+154, dvrr_stack+21, dvrr_stack+40, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+190, dvrr_stack+84, dvrr_stack+50, dvrr_stack+74, dvrr_stack+21, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+211, dvrr_stack+190, dvrr_stack+169, NULL, NULL, dvrr_stack+50);
 tmp = dvrr_stack + 211;
 target_ptr = Libderiv->dvrr_classes[1][5];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+274,dvrr_stack+211,dvrr_stack+99,3);


 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+71, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+409, dvrr_stack+3, dvrr_stack+71, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+415, dvrr_stack+65, dvrr_stack+409, dvrr_stack+31, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+425, dvrr_stack+144, dvrr_stack+415, dvrr_stack+34, dvrr_stack+65, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+440, dvrr_stack+154, dvrr_stack+425, dvrr_stack+40, dvrr_stack+144, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+461, dvrr_stack+169, dvrr_stack+440, dvrr_stack+50, dvrr_stack+154, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+489, dvrr_stack+190, dvrr_stack+169, dvrr_stack+84, dvrr_stack+50, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+517, dvrr_stack+489, dvrr_stack+461, NULL, NULL, dvrr_stack+169);
 tmp = dvrr_stack + 517;
 target_ptr = Libderiv->dvrr_classes[1][6];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+601,dvrr_stack+517,dvrr_stack+211,3);


 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+790, dvrr_stack+71, dvrr_stack+31, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+796, dvrr_stack+409, dvrr_stack+790, dvrr_stack+3, dvrr_stack+71, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+806, dvrr_stack+415, dvrr_stack+796, dvrr_stack+65, dvrr_stack+409, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+821, dvrr_stack+425, dvrr_stack+806, dvrr_stack+144, dvrr_stack+415, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+842, dvrr_stack+440, dvrr_stack+821, dvrr_stack+154, dvrr_stack+425, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+870, dvrr_stack+461, dvrr_stack+842, dvrr_stack+169, dvrr_stack+440, NULL);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+906, dvrr_stack+489, dvrr_stack+461, dvrr_stack+190, dvrr_stack+169, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+942, dvrr_stack+906, dvrr_stack+870, NULL, NULL, dvrr_stack+461);
 tmp = dvrr_stack + 942;
 target_ptr = Libderiv->dvrr_classes[1][7];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+1050,dvrr_stack+942,dvrr_stack+517,3);


 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+65, dvrr_stack+31, dvrr_stack+3, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1302, dvrr_stack+790, dvrr_stack+65, dvrr_stack+71, dvrr_stack+31, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1312, dvrr_stack+796, dvrr_stack+1302, dvrr_stack+409, dvrr_stack+790, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1327, dvrr_stack+806, dvrr_stack+1312, dvrr_stack+415, dvrr_stack+796, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1348, dvrr_stack+821, dvrr_stack+1327, dvrr_stack+425, dvrr_stack+806, NULL);

 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1376, dvrr_stack+842, dvrr_stack+1348, dvrr_stack+440, dvrr_stack+821, NULL);

 /* compute (0 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1412, dvrr_stack+870, dvrr_stack+1376, dvrr_stack+461, dvrr_stack+842, NULL);

 /* compute (0 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1457, dvrr_stack+906, dvrr_stack+870, dvrr_stack+489, dvrr_stack+461, NULL);

 /* compute (1 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1502, dvrr_stack+1457, dvrr_stack+1412, NULL, NULL, dvrr_stack+870);

 /* compute (1 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+1637,dvrr_stack+1502,dvrr_stack+942,3);


 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1961, dvrr_stack+21, dvrr_stack+40, NULL, NULL, dvrr_stack+9);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1991, dvrr_stack+50, dvrr_stack+154, NULL, NULL, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2036, dvrr_stack+99, dvrr_stack+1991, dvrr_stack+84, dvrr_stack+50, dvrr_stack+1961);
 tmp = dvrr_stack + 2036;
 target_ptr = Libderiv->dvrr_classes[2][4];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2126, dvrr_stack+169, dvrr_stack+440, NULL, NULL, dvrr_stack+154);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2189, dvrr_stack+211, dvrr_stack+2126, dvrr_stack+190, dvrr_stack+169, dvrr_stack+1991);
 tmp = dvrr_stack + 2189;
 target_ptr = Libderiv->dvrr_classes[2][5];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+2315,dvrr_stack+2189,dvrr_stack+2036,6);


 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2585, dvrr_stack+461, dvrr_stack+842, NULL, NULL, dvrr_stack+440);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2669, dvrr_stack+517, dvrr_stack+2585, dvrr_stack+489, dvrr_stack+461, dvrr_stack+2126);
 tmp = dvrr_stack + 2669;
 target_ptr = Libderiv->dvrr_classes[2][6];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+2837,dvrr_stack+2669,dvrr_stack+2189,6);


 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3215, dvrr_stack+870, dvrr_stack+1376, NULL, NULL, dvrr_stack+842);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3323, dvrr_stack+942, dvrr_stack+3215, dvrr_stack+906, dvrr_stack+870, dvrr_stack+2585);

 /* compute (2 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+3539,dvrr_stack+3323,dvrr_stack+2669,6);


 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+71, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+409, dvrr_stack+3, dvrr_stack+71, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+415, dvrr_stack+65, dvrr_stack+409, dvrr_stack+31, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4043, dvrr_stack+1302, dvrr_stack+415, dvrr_stack+790, dvrr_stack+65, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4058, dvrr_stack+1312, dvrr_stack+4043, dvrr_stack+796, dvrr_stack+1302, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4079, dvrr_stack+1327, dvrr_stack+4058, dvrr_stack+806, dvrr_stack+1312, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4043, dvrr_stack+1348, dvrr_stack+4079, dvrr_stack+821, dvrr_stack+1327, NULL);

 /* compute (0 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+4079, dvrr_stack+1376, dvrr_stack+4043, dvrr_stack+842, dvrr_stack+1348, NULL);

 /* compute (1 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+4124, dvrr_stack+1412, dvrr_stack+4079, NULL, NULL, dvrr_stack+1376);

 /* compute (2 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+4259, dvrr_stack+1502, dvrr_stack+4124, dvrr_stack+1457, dvrr_stack+1412, dvrr_stack+3215);

 /* compute (2 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+4529,dvrr_stack+4259,dvrr_stack+3323,6);


 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1412, dvrr_stack+74, dvrr_stack+21, NULL, NULL, dvrr_stack+15);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1442, dvrr_stack+15, dvrr_stack+9, NULL, NULL, dvrr_stack+0);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4079, dvrr_stack+1412, dvrr_stack+1961, dvrr_stack+74, dvrr_stack+21, dvrr_stack+1442);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1442, dvrr_stack+9, dvrr_stack+34, NULL, NULL, dvrr_stack+6);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1460, dvrr_stack+40, dvrr_stack+144, NULL, NULL, dvrr_stack+34);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4139, dvrr_stack+1961, dvrr_stack+1460, dvrr_stack+21, dvrr_stack+40, dvrr_stack+1442);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+154, dvrr_stack+425, NULL, NULL, dvrr_stack+144);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5177, dvrr_stack+1991, dvrr_stack+0, dvrr_stack+50, dvrr_stack+154, dvrr_stack+1460);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5267, dvrr_stack+2036, dvrr_stack+5177, dvrr_stack+99, dvrr_stack+1991, dvrr_stack+4139);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4139, dvrr_stack+440, dvrr_stack+821, NULL, NULL, dvrr_stack+425);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+5417, dvrr_stack+2126, dvrr_stack+4139, dvrr_stack+169, dvrr_stack+440, dvrr_stack+0);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+5543, dvrr_stack+2189, dvrr_stack+5417, dvrr_stack+211, dvrr_stack+2126, dvrr_stack+5177);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+842, dvrr_stack+1348, NULL, NULL, dvrr_stack+821);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+5753, dvrr_stack+2585, dvrr_stack+0, dvrr_stack+461, dvrr_stack+842, dvrr_stack+4139);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+5921, dvrr_stack+2669, dvrr_stack+5753, dvrr_stack+517, dvrr_stack+2585, dvrr_stack+5417);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+5417, dvrr_stack+1376, dvrr_stack+4043, NULL, NULL, dvrr_stack+1348);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6201, dvrr_stack+3215, dvrr_stack+5417, dvrr_stack+870, dvrr_stack+1376, dvrr_stack+0);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6417, dvrr_stack+3323, dvrr_stack+6201, dvrr_stack+942, dvrr_stack+3215, dvrr_stack+5753);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+5753, dvrr_stack+274, NULL);
 tmp = dvrr_stack + 5753;
 target_ptr = Libderiv->deriv_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,63,dvrr_stack+2126, dvrr_stack+601, NULL);
 tmp = dvrr_stack + 2126;
 target_ptr = Libderiv->deriv_classes[1][5][11];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,84,dvrr_stack+0, dvrr_stack+1050, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[1][6][11];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,108,dvrr_stack+3215, dvrr_stack+1637, NULL);
 tmp = dvrr_stack + 3215;
 target_ptr = Libderiv->deriv_classes[1][7][11];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+5177, dvrr_stack+2315, NULL);
 tmp = dvrr_stack + 5177;
 target_ptr = Libderiv->deriv_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+5417, dvrr_stack+2837, NULL);
 tmp = dvrr_stack + 5417;
 target_ptr = Libderiv->deriv_classes[2][5][11];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,168,dvrr_stack+6201, dvrr_stack+3539, NULL);
 tmp = dvrr_stack + 6201;
 target_ptr = Libderiv->deriv_classes[2][6][11];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,216,dvrr_stack+6777, dvrr_stack+4529, NULL);
 tmp = dvrr_stack + 6777;
 target_ptr = Libderiv->deriv_classes[2][7][11];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+5798, dvrr_stack+274, NULL);
 tmp = dvrr_stack + 5798;
 target_ptr = Libderiv->deriv_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+5843, dvrr_stack+601, NULL);
 tmp = dvrr_stack + 5843;
 target_ptr = Libderiv->deriv_classes[1][5][10];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,84,dvrr_stack+2585, dvrr_stack+1050, NULL);
 tmp = dvrr_stack + 2585;
 target_ptr = Libderiv->deriv_classes[1][6][10];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,108,dvrr_stack+4139, dvrr_stack+1637, NULL);
 tmp = dvrr_stack + 4139;
 target_ptr = Libderiv->deriv_classes[1][7][10];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+1302, dvrr_stack+2315, NULL);
 tmp = dvrr_stack + 1302;
 target_ptr = Libderiv->deriv_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+6993, dvrr_stack+2837, NULL);
 tmp = dvrr_stack + 6993;
 target_ptr = Libderiv->deriv_classes[2][5][10];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,168,dvrr_stack+7119, dvrr_stack+3539, NULL);
 tmp = dvrr_stack + 7119;
 target_ptr = Libderiv->deriv_classes[2][6][10];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,216,dvrr_stack+7287, dvrr_stack+4529, NULL);
 tmp = dvrr_stack + 7287;
 target_ptr = Libderiv->deriv_classes[2][7][10];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+6369, dvrr_stack+274, NULL);
 tmp = dvrr_stack + 6369;
 target_ptr = Libderiv->deriv_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+274, dvrr_stack+601, NULL);
 tmp = dvrr_stack + 274;
 target_ptr = Libderiv->deriv_classes[1][5][9];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,84,dvrr_stack+601, dvrr_stack+1050, NULL);
 tmp = dvrr_stack + 601;
 target_ptr = Libderiv->deriv_classes[1][6][9];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,108,dvrr_stack+1050, dvrr_stack+1637, NULL);
 tmp = dvrr_stack + 1050;
 target_ptr = Libderiv->deriv_classes[1][7][9];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+1637, dvrr_stack+2315, NULL);
 tmp = dvrr_stack + 1637;
 target_ptr = Libderiv->deriv_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+2315, dvrr_stack+2837, NULL);
 tmp = dvrr_stack + 2315;
 target_ptr = Libderiv->deriv_classes[2][5][9];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,168,dvrr_stack+2837, dvrr_stack+3539, NULL);
 tmp = dvrr_stack + 2837;
 target_ptr = Libderiv->deriv_classes[2][6][9];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,216,dvrr_stack+3539, dvrr_stack+4529, NULL);
 tmp = dvrr_stack + 3539;
 target_ptr = Libderiv->deriv_classes[2][7][9];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+4529, dvrr_stack+211, dvrr_stack+1412);
 tmp = dvrr_stack + 4529;
 target_ptr = Libderiv->deriv_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+4574, dvrr_stack+517, dvrr_stack+99);
 tmp = dvrr_stack + 4574;
 target_ptr = Libderiv->deriv_classes[1][5][8];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+4637, dvrr_stack+942, dvrr_stack+211);
 tmp = dvrr_stack + 4637;
 target_ptr = Libderiv->deriv_classes[1][6][8];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,3,1,dvrr_stack+4721, dvrr_stack+1502, dvrr_stack+517);
 tmp = dvrr_stack + 4721;
 target_ptr = Libderiv->deriv_classes[1][7][8];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+4829, dvrr_stack+2189, dvrr_stack+4079);
 tmp = dvrr_stack + 4829;
 target_ptr = Libderiv->deriv_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+4919, dvrr_stack+2669, dvrr_stack+2036);
 tmp = dvrr_stack + 4919;
 target_ptr = Libderiv->deriv_classes[2][5][8];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+3755, dvrr_stack+3323, dvrr_stack+2189);
 tmp = dvrr_stack + 3755;
 target_ptr = Libderiv->deriv_classes[2][6][8];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,6,1,dvrr_stack+1727, dvrr_stack+4259, dvrr_stack+2669);
 tmp = dvrr_stack + 1727;
 target_ptr = Libderiv->deriv_classes[2][7][8];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+5045, dvrr_stack+211, dvrr_stack+1412);
 tmp = dvrr_stack + 5045;
 target_ptr = Libderiv->deriv_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+5090, dvrr_stack+517, dvrr_stack+99);
 tmp = dvrr_stack + 5090;
 target_ptr = Libderiv->deriv_classes[1][5][7];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+3923, dvrr_stack+942, dvrr_stack+211);
 tmp = dvrr_stack + 3923;
 target_ptr = Libderiv->deriv_classes[1][6][7];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,3,1,dvrr_stack+3005, dvrr_stack+1502, dvrr_stack+517);
 tmp = dvrr_stack + 3005;
 target_ptr = Libderiv->deriv_classes[1][7][7];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+3113, dvrr_stack+2189, dvrr_stack+4079);
 tmp = dvrr_stack + 3113;
 target_ptr = Libderiv->deriv_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+2441, dvrr_stack+2669, dvrr_stack+2036);
 tmp = dvrr_stack + 2441;
 target_ptr = Libderiv->deriv_classes[2][5][7];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+685, dvrr_stack+3323, dvrr_stack+2189);
 tmp = dvrr_stack + 685;
 target_ptr = Libderiv->deriv_classes[2][6][7];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,6,1,dvrr_stack+7503, dvrr_stack+4259, dvrr_stack+2669);
 tmp = dvrr_stack + 7503;
 target_ptr = Libderiv->deriv_classes[2][7][7];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+4007, dvrr_stack+211, dvrr_stack+1412);
 tmp = dvrr_stack + 4007;
 target_ptr = Libderiv->deriv_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+1943, dvrr_stack+517, dvrr_stack+99);
 tmp = dvrr_stack + 1943;
 target_ptr = Libderiv->deriv_classes[1][5][6];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+1158, dvrr_stack+942, dvrr_stack+211);
 tmp = dvrr_stack + 1158;
 target_ptr = Libderiv->deriv_classes[1][6][6];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,3,1,dvrr_stack+337, dvrr_stack+1502, dvrr_stack+517);
 tmp = dvrr_stack + 337;
 target_ptr = Libderiv->deriv_classes[1][7][6];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+7719, dvrr_stack+2189, dvrr_stack+4079);
 tmp = dvrr_stack + 7719;
 target_ptr = Libderiv->deriv_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+7809, dvrr_stack+2669, dvrr_stack+2036);
 tmp = dvrr_stack + 7809;
 target_ptr = Libderiv->deriv_classes[2][5][6];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+1392, dvrr_stack+3323, dvrr_stack+2189);
 tmp = dvrr_stack + 1392;
 target_ptr = Libderiv->deriv_classes[2][6][6];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,6,1,dvrr_stack+7935, dvrr_stack+4259, dvrr_stack+2669);
 tmp = dvrr_stack + 7935;
 target_ptr = Libderiv->deriv_classes[2][7][6];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+4052, dvrr_stack+2036, dvrr_stack+84);
 tmp = dvrr_stack + 4052;
 target_ptr = Libderiv->deriv_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+4247, dvrr_stack+2189, dvrr_stack+190);
 tmp = dvrr_stack + 4247;
 target_ptr = Libderiv->deriv_classes[1][5][2];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+4310, dvrr_stack+2669, dvrr_stack+489);
 tmp = dvrr_stack + 4310;
 target_ptr = Libderiv->deriv_classes[1][6][2];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,36,dvrr_stack+4394, dvrr_stack+3323, dvrr_stack+906);
 tmp = dvrr_stack + 4394;
 target_ptr = Libderiv->deriv_classes[1][7][2];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+8151, dvrr_stack+5267, dvrr_stack+99);
 tmp = dvrr_stack + 8151;
 target_ptr = Libderiv->deriv_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+8241, dvrr_stack+5543, dvrr_stack+211);
 tmp = dvrr_stack + 8241;
 target_ptr = Libderiv->deriv_classes[2][5][2];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+8367, dvrr_stack+5921, dvrr_stack+517);
 tmp = dvrr_stack + 8367;
 target_ptr = Libderiv->deriv_classes[2][6][2];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,36,dvrr_stack+8535, dvrr_stack+6417, dvrr_stack+942);
 tmp = dvrr_stack + 8535;
 target_ptr = Libderiv->deriv_classes[2][7][2];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+1242, dvrr_stack+2036, dvrr_stack+84);
 tmp = dvrr_stack + 1242;
 target_ptr = Libderiv->deriv_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+8751, dvrr_stack+2189, dvrr_stack+190);
 tmp = dvrr_stack + 8751;
 target_ptr = Libderiv->deriv_classes[1][5][1];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+8814, dvrr_stack+2669, dvrr_stack+489);
 tmp = dvrr_stack + 8814;
 target_ptr = Libderiv->deriv_classes[1][6][1];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,36,dvrr_stack+8898, dvrr_stack+3323, dvrr_stack+906);
 tmp = dvrr_stack + 8898;
 target_ptr = Libderiv->deriv_classes[1][7][1];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+9006, dvrr_stack+5267, dvrr_stack+99);
 tmp = dvrr_stack + 9006;
 target_ptr = Libderiv->deriv_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+9096, dvrr_stack+5543, dvrr_stack+211);
 tmp = dvrr_stack + 9096;
 target_ptr = Libderiv->deriv_classes[2][5][1];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+9222, dvrr_stack+5921, dvrr_stack+517);
 tmp = dvrr_stack + 9222;
 target_ptr = Libderiv->deriv_classes[2][6][1];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,36,dvrr_stack+9390, dvrr_stack+6417, dvrr_stack+942);
 tmp = dvrr_stack + 9390;
 target_ptr = Libderiv->deriv_classes[2][7][1];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+853, dvrr_stack+2036, dvrr_stack+84);
 tmp = dvrr_stack + 853;
 target_ptr = Libderiv->deriv_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+2006, dvrr_stack+2189, dvrr_stack+190);
 tmp = dvrr_stack + 2006;
 target_ptr = Libderiv->deriv_classes[1][5][0];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+2189, dvrr_stack+2669, dvrr_stack+489);
 tmp = dvrr_stack + 2189;
 target_ptr = Libderiv->deriv_classes[1][6][0];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,36,dvrr_stack+2669, dvrr_stack+3323, dvrr_stack+906);
 tmp = dvrr_stack + 2669;
 target_ptr = Libderiv->deriv_classes[1][7][0];
 for(i=0;i<108;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+3323, dvrr_stack+5267, dvrr_stack+99);
 tmp = dvrr_stack + 3323;
 target_ptr = Libderiv->deriv_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+3413, dvrr_stack+5543, dvrr_stack+211);
 tmp = dvrr_stack + 3413;
 target_ptr = Libderiv->deriv_classes[2][5][0];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+5543, dvrr_stack+5921, dvrr_stack+517);
 tmp = dvrr_stack + 5543;
 target_ptr = Libderiv->deriv_classes[2][6][0];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,36,dvrr_stack+5906, dvrr_stack+6417, dvrr_stack+942);
 tmp = dvrr_stack + 5906;
 target_ptr = Libderiv->deriv_classes[2][7][0];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];


}

