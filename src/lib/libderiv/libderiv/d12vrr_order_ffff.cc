#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (ff|ff) integrals */

void d12vrr_order_ffff(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+3);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+21, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+24, dvrr_stack+21, dvrr_stack+3, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+30, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+33, dvrr_stack+0, dvrr_stack+30, Data->F+3, Data->F+4, NULL);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+39, dvrr_stack+15, dvrr_stack+33, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+57, dvrr_stack+24, dvrr_stack+15, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+75, dvrr_stack+57, dvrr_stack+39, dvrr_stack+24, dvrr_stack+15, dvrr_stack+6);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+15, dvrr_stack+33, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+121, dvrr_stack+24, dvrr_stack+15, dvrr_stack+21, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+131, dvrr_stack+121, dvrr_stack+111, NULL, NULL, dvrr_stack+15);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+161, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+164, dvrr_stack+161, dvrr_stack+21, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+170, dvrr_stack+164, dvrr_stack+24, dvrr_stack+161, dvrr_stack+21, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+180, dvrr_stack+170, dvrr_stack+121, NULL, NULL, dvrr_stack+24);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+210, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+213, dvrr_stack+30, dvrr_stack+210, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+219, dvrr_stack+33, dvrr_stack+213, dvrr_stack+0, dvrr_stack+30, NULL);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+229, dvrr_stack+111, dvrr_stack+219, NULL, NULL, dvrr_stack+33);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+259, dvrr_stack+131, dvrr_stack+229, dvrr_stack+121, dvrr_stack+111, dvrr_stack+39);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+319, dvrr_stack+180, dvrr_stack+131, dvrr_stack+170, dvrr_stack+121, dvrr_stack+57);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+379, dvrr_stack+319, dvrr_stack+259, dvrr_stack+180, dvrr_stack+131, dvrr_stack+75);
 tmp = dvrr_stack + 379;
 target_ptr = Libderiv->dvrr_classes[3][3];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+479, dvrr_stack+111, dvrr_stack+219, dvrr_stack+15, dvrr_stack+33, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+494, dvrr_stack+121, dvrr_stack+111, dvrr_stack+24, dvrr_stack+15, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+509, dvrr_stack+494, dvrr_stack+479, NULL, NULL, dvrr_stack+111);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+554, dvrr_stack+170, dvrr_stack+121, dvrr_stack+164, dvrr_stack+24, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+569, dvrr_stack+554, dvrr_stack+494, NULL, NULL, dvrr_stack+121);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+614, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+617, dvrr_stack+210, dvrr_stack+614, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+623, dvrr_stack+213, dvrr_stack+617, dvrr_stack+30, dvrr_stack+210, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+633, dvrr_stack+219, dvrr_stack+623, dvrr_stack+33, dvrr_stack+213, NULL);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+648, dvrr_stack+479, dvrr_stack+633, NULL, NULL, dvrr_stack+219);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+693, dvrr_stack+509, dvrr_stack+648, dvrr_stack+494, dvrr_stack+479, dvrr_stack+229);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+783, dvrr_stack+569, dvrr_stack+509, dvrr_stack+554, dvrr_stack+494, dvrr_stack+131);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+873, dvrr_stack+783, dvrr_stack+693, dvrr_stack+569, dvrr_stack+509, dvrr_stack+259);
 tmp = dvrr_stack + 873;
 target_ptr = Libderiv->dvrr_classes[3][4];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+1023,dvrr_stack+873,dvrr_stack+379,10);


 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1323, dvrr_stack+479, dvrr_stack+633, dvrr_stack+111, dvrr_stack+219, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1344, dvrr_stack+494, dvrr_stack+479, dvrr_stack+121, dvrr_stack+111, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1365, dvrr_stack+1344, dvrr_stack+1323, NULL, NULL, dvrr_stack+479);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1428, dvrr_stack+554, dvrr_stack+494, dvrr_stack+170, dvrr_stack+121, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1449, dvrr_stack+1428, dvrr_stack+1344, NULL, NULL, dvrr_stack+494);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+121, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+124, dvrr_stack+614, dvrr_stack+121, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+170, dvrr_stack+617, dvrr_stack+124, dvrr_stack+210, dvrr_stack+614, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1512, dvrr_stack+623, dvrr_stack+170, dvrr_stack+213, dvrr_stack+617, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1527, dvrr_stack+633, dvrr_stack+1512, dvrr_stack+219, dvrr_stack+623, NULL);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1548, dvrr_stack+1323, dvrr_stack+1527, NULL, NULL, dvrr_stack+633);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1611, dvrr_stack+1365, dvrr_stack+1548, dvrr_stack+1344, dvrr_stack+1323, dvrr_stack+648);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1737, dvrr_stack+1449, dvrr_stack+1365, dvrr_stack+1428, dvrr_stack+1344, dvrr_stack+509);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1863, dvrr_stack+1737, dvrr_stack+1611, dvrr_stack+1449, dvrr_stack+1365, dvrr_stack+693);
 tmp = dvrr_stack + 1863;
 target_ptr = Libderiv->dvrr_classes[3][5];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+2073,dvrr_stack+1863,dvrr_stack+873,10);


 /* compute (3 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+2523,dvrr_stack+2073,dvrr_stack+1023,10);


 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,100,dvrr_stack+3123, dvrr_stack+2523, dvrr_stack+379);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3423, dvrr_stack+1323, dvrr_stack+1527, dvrr_stack+479, dvrr_stack+633, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3451, dvrr_stack+1344, dvrr_stack+1323, dvrr_stack+494, dvrr_stack+479, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3479, dvrr_stack+3451, dvrr_stack+3423, NULL, NULL, dvrr_stack+1323);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3563, dvrr_stack+1428, dvrr_stack+1344, dvrr_stack+554, dvrr_stack+494, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3591, dvrr_stack+3563, dvrr_stack+3451, NULL, NULL, dvrr_stack+1344);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+494, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+497, dvrr_stack+121, dvrr_stack+494, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+554, dvrr_stack+124, dvrr_stack+497, dvrr_stack+614, dvrr_stack+121, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3675, dvrr_stack+170, dvrr_stack+554, dvrr_stack+617, dvrr_stack+124, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3690, dvrr_stack+1512, dvrr_stack+3675, dvrr_stack+623, dvrr_stack+170, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3711, dvrr_stack+1527, dvrr_stack+3690, dvrr_stack+633, dvrr_stack+1512, NULL);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3739, dvrr_stack+3423, dvrr_stack+3711, NULL, NULL, dvrr_stack+1527);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3823, dvrr_stack+3479, dvrr_stack+3739, dvrr_stack+3451, dvrr_stack+3423, dvrr_stack+1548);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3991, dvrr_stack+3591, dvrr_stack+3479, dvrr_stack+3563, dvrr_stack+3451, dvrr_stack+1365);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4159, dvrr_stack+3991, dvrr_stack+3823, dvrr_stack+3591, dvrr_stack+3479, dvrr_stack+1611);
 tmp = dvrr_stack + 4159;
 target_ptr = Libderiv->dvrr_classes[3][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+4439,dvrr_stack+4159,dvrr_stack+1863,10);


 /* compute (3 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+5069,dvrr_stack+4439,dvrr_stack+2073,10);


 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+5969, dvrr_stack+5069, dvrr_stack+873);

 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6419, dvrr_stack+3423, dvrr_stack+3711, dvrr_stack+1323, dvrr_stack+1527, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6455, dvrr_stack+3451, dvrr_stack+3423, dvrr_stack+1344, dvrr_stack+1323, NULL);

 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6491, dvrr_stack+6455, dvrr_stack+6419, NULL, NULL, dvrr_stack+3423);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6599, dvrr_stack+3563, dvrr_stack+3451, dvrr_stack+1428, dvrr_stack+1344, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6635, dvrr_stack+6599, dvrr_stack+6455, NULL, NULL, dvrr_stack+3451);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1344, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+503, dvrr_stack+494, dvrr_stack+1344, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1347, dvrr_stack+497, dvrr_stack+503, dvrr_stack+121, dvrr_stack+494, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1428, dvrr_stack+554, dvrr_stack+1347, dvrr_stack+124, dvrr_stack+497, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6743, dvrr_stack+3675, dvrr_stack+1428, dvrr_stack+170, dvrr_stack+554, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6764, dvrr_stack+3690, dvrr_stack+6743, dvrr_stack+1512, dvrr_stack+3675, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6792, dvrr_stack+3711, dvrr_stack+6764, dvrr_stack+1527, dvrr_stack+3690, NULL);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6828, dvrr_stack+6419, dvrr_stack+6792, NULL, NULL, dvrr_stack+3711);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6936, dvrr_stack+6491, dvrr_stack+6828, dvrr_stack+6455, dvrr_stack+6419, dvrr_stack+3739);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7152, dvrr_stack+6635, dvrr_stack+6491, dvrr_stack+6599, dvrr_stack+6455, dvrr_stack+3479);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7368, dvrr_stack+7152, dvrr_stack+6936, dvrr_stack+6635, dvrr_stack+6491, dvrr_stack+3823);

 /* compute (3 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+7728,dvrr_stack+7368,dvrr_stack+4159,10);


 /* compute (3 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+8568,dvrr_stack+7728,dvrr_stack+4439,10);


 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,210,dvrr_stack+9828, dvrr_stack+8568, dvrr_stack+1863);

 /* compute (0 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6635, dvrr_stack+6419, dvrr_stack+6792, dvrr_stack+3423, dvrr_stack+3711, NULL);

 /* compute (0 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6680, dvrr_stack+6455, dvrr_stack+6419, dvrr_stack+3451, dvrr_stack+3423, NULL);

 /* compute (1 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10458, dvrr_stack+6680, dvrr_stack+6635, NULL, NULL, dvrr_stack+6419);

 /* compute (0 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10593, dvrr_stack+6599, dvrr_stack+6455, dvrr_stack+3563, dvrr_stack+3451, NULL);

 /* compute (1 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10638, dvrr_stack+10593, dvrr_stack+6680, NULL, NULL, dvrr_stack+6455);

 /* compute (0 0 | 1 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6455, Data->F+10, Data->F+11, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1443, dvrr_stack+1344, dvrr_stack+6455, Data->F+9, Data->F+10, NULL);

 /* compute (0 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+6458, dvrr_stack+503, dvrr_stack+1443, dvrr_stack+494, dvrr_stack+1344, NULL);

 /* compute (0 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6468, dvrr_stack+1347, dvrr_stack+6458, dvrr_stack+497, dvrr_stack+503, NULL);

 /* compute (0 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3451, dvrr_stack+1428, dvrr_stack+6468, dvrr_stack+554, dvrr_stack+1347, NULL);

 /* compute (0 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3563, dvrr_stack+6743, dvrr_stack+3451, dvrr_stack+3675, dvrr_stack+1428, NULL);

 /* compute (0 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6599, dvrr_stack+6764, dvrr_stack+3563, dvrr_stack+3690, dvrr_stack+6743, NULL);

 /* compute (0 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10773, dvrr_stack+6792, dvrr_stack+6599, dvrr_stack+3711, dvrr_stack+6764, NULL);

 /* compute (1 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10818, dvrr_stack+6635, dvrr_stack+10773, NULL, NULL, dvrr_stack+6792);

 /* compute (2 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10953, dvrr_stack+10458, dvrr_stack+10818, dvrr_stack+6680, dvrr_stack+6635, dvrr_stack+6828);

 /* compute (2 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+11223, dvrr_stack+10638, dvrr_stack+10458, dvrr_stack+10593, dvrr_stack+6680, dvrr_stack+6491);

 /* compute (3 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+11493, dvrr_stack+11223, dvrr_stack+10953, dvrr_stack+10638, dvrr_stack+10458, dvrr_stack+6936);

 /* compute (3 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+11943,dvrr_stack+11493,dvrr_stack+7368,10);


 /* compute (3 0 | 6 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_id(Libderiv->CD,dvrr_stack+13023,dvrr_stack+11943,dvrr_stack+7728,10);


 /* compute (3 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,280,dvrr_stack+14703, dvrr_stack+13023, dvrr_stack+4159);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+6680, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+6683, dvrr_stack+0, dvrr_stack+30, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+6692, dvrr_stack+6, dvrr_stack+6683, dvrr_stack+3, dvrr_stack+0, dvrr_stack+6680);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+6710, dvrr_stack+33, dvrr_stack+213, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+10593, dvrr_stack+39, dvrr_stack+6710, dvrr_stack+15, dvrr_stack+33, dvrr_stack+6683);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+10629, dvrr_stack+75, dvrr_stack+10593, dvrr_stack+57, dvrr_stack+39, dvrr_stack+6692);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+10689, dvrr_stack+219, dvrr_stack+623, NULL, NULL, dvrr_stack+213);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+15543, dvrr_stack+229, dvrr_stack+10689, dvrr_stack+111, dvrr_stack+219, dvrr_stack+6710);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+15603, dvrr_stack+259, dvrr_stack+15543, dvrr_stack+131, dvrr_stack+229, dvrr_stack+10593);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+15703, dvrr_stack+379, dvrr_stack+15603, dvrr_stack+319, dvrr_stack+259, dvrr_stack+10629);
 tmp = dvrr_stack + 15703;
 target_ptr = Libderiv->dvrr_classes[4][3];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10719, dvrr_stack+633, dvrr_stack+1512, NULL, NULL, dvrr_stack+623);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15853, dvrr_stack+648, dvrr_stack+10719, dvrr_stack+479, dvrr_stack+633, dvrr_stack+10689);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15943, dvrr_stack+693, dvrr_stack+15853, dvrr_stack+509, dvrr_stack+648, dvrr_stack+15543);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+16093, dvrr_stack+873, dvrr_stack+15943, dvrr_stack+783, dvrr_stack+693, dvrr_stack+15603);
 tmp = dvrr_stack + 16093;
 target_ptr = Libderiv->dvrr_classes[4][4];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+16318,dvrr_stack+16093,dvrr_stack+15703,15);


 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+16768, dvrr_stack+1527, dvrr_stack+3690, NULL, NULL, dvrr_stack+1512);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+16831, dvrr_stack+1548, dvrr_stack+16768, dvrr_stack+1323, dvrr_stack+1527, dvrr_stack+10719);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+16957, dvrr_stack+1611, dvrr_stack+16831, dvrr_stack+1365, dvrr_stack+1548, dvrr_stack+15853);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+17167, dvrr_stack+1863, dvrr_stack+16957, dvrr_stack+1737, dvrr_stack+1611, dvrr_stack+15943);
 tmp = dvrr_stack + 17167;
 target_ptr = Libderiv->dvrr_classes[4][5];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+17482,dvrr_stack+17167,dvrr_stack+16093,15);


 /* compute (4 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+18157,dvrr_stack+17482,dvrr_stack+16318,15);


 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+19057, dvrr_stack+18157, dvrr_stack+15703);

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+19507, dvrr_stack+3711, dvrr_stack+6764, NULL, NULL, dvrr_stack+3690);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+19591, dvrr_stack+3739, dvrr_stack+19507, dvrr_stack+3423, dvrr_stack+3711, dvrr_stack+16768);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+19759, dvrr_stack+3823, dvrr_stack+19591, dvrr_stack+3479, dvrr_stack+3739, dvrr_stack+16831);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+20039, dvrr_stack+4159, dvrr_stack+19759, dvrr_stack+3991, dvrr_stack+3823, dvrr_stack+16957);
 tmp = dvrr_stack + 20039;
 target_ptr = Libderiv->dvrr_classes[4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+20459,dvrr_stack+20039,dvrr_stack+17167,15);


 /* compute (4 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+21404,dvrr_stack+20459,dvrr_stack+17482,15);


 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,225,dvrr_stack+22754, dvrr_stack+21404, dvrr_stack+16093);

 /* compute (1 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+23429, dvrr_stack+6792, dvrr_stack+6599, NULL, NULL, dvrr_stack+6764);

 /* compute (2 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+23537, dvrr_stack+6828, dvrr_stack+23429, dvrr_stack+6419, dvrr_stack+6792, dvrr_stack+19507);

 /* compute (3 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+23753, dvrr_stack+6936, dvrr_stack+23537, dvrr_stack+6491, dvrr_stack+6828, dvrr_stack+19591);

 /* compute (4 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+24113, dvrr_stack+7368, dvrr_stack+23753, dvrr_stack+7152, dvrr_stack+6936, dvrr_stack+19759);

 /* compute (4 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+24653,dvrr_stack+24113,dvrr_stack+20039,15);


 /* compute (4 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+25913,dvrr_stack+24653,dvrr_stack+20459,15);


 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,315,dvrr_stack+27803, dvrr_stack+25913, dvrr_stack+17167);

 /* compute (0 0 | 1 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6419, Data->F+11, Data->F+12, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+6455, dvrr_stack+6419, Data->F+10, Data->F+11, NULL);

 /* compute (0 0 | 3 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+1443, dvrr_stack+15, dvrr_stack+1344, dvrr_stack+6455, NULL);

 /* compute (0 0 | 4 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+479, dvrr_stack+6458, dvrr_stack+111, dvrr_stack+503, dvrr_stack+1443, NULL);

 /* compute (0 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6422, dvrr_stack+6468, dvrr_stack+479, dvrr_stack+1347, dvrr_stack+6458, NULL);

 /* compute (0 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3423, dvrr_stack+3451, dvrr_stack+6422, dvrr_stack+1428, dvrr_stack+6468, NULL);

 /* compute (0 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+509, dvrr_stack+3563, dvrr_stack+3423, dvrr_stack+6743, dvrr_stack+3451, NULL);

 /* compute (0 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6483, dvrr_stack+6599, dvrr_stack+509, dvrr_stack+6764, dvrr_stack+3563, NULL);

 /* compute (1 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+28748, dvrr_stack+10773, dvrr_stack+6483, NULL, NULL, dvrr_stack+6599);

 /* compute (2 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+28883, dvrr_stack+10818, dvrr_stack+28748, dvrr_stack+6635, dvrr_stack+10773, dvrr_stack+23429);

 /* compute (3 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+29153, dvrr_stack+10953, dvrr_stack+28883, dvrr_stack+10458, dvrr_stack+10818, dvrr_stack+23537);

 /* compute (4 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+29603, dvrr_stack+11493, dvrr_stack+29153, dvrr_stack+11223, dvrr_stack+10953, dvrr_stack+23753);

 /* compute (4 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+30278,dvrr_stack+29603,dvrr_stack+24113,15);


 /* compute (4 0 | 6 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_id(Libderiv->CD,dvrr_stack+31898,dvrr_stack+30278,dvrr_stack+24653,15);


 /* compute (4 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,420,dvrr_stack+34418, dvrr_stack+31898, dvrr_stack+20039);

 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+11223, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+11226, dvrr_stack+6680, dvrr_stack+11223, Data->F+3, Data->F+4, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+545, dvrr_stack+30, dvrr_stack+210, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+11232, dvrr_stack+6683, dvrr_stack+545, dvrr_stack+0, dvrr_stack+30, dvrr_stack+11223);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+11250, dvrr_stack+6692, dvrr_stack+11232, dvrr_stack+6, dvrr_stack+6683, dvrr_stack+11226);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+11280, dvrr_stack+213, dvrr_stack+617, NULL, NULL, dvrr_stack+210);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+11298, dvrr_stack+6710, dvrr_stack+11280, dvrr_stack+33, dvrr_stack+213, dvrr_stack+545);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+11334, dvrr_stack+10593, dvrr_stack+11298, dvrr_stack+39, dvrr_stack+6710, dvrr_stack+11232);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+11394, dvrr_stack+10629, dvrr_stack+11334, dvrr_stack+75, dvrr_stack+10593, dvrr_stack+11250);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+10458, dvrr_stack+623, dvrr_stack+170, NULL, NULL, dvrr_stack+617);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+10488, dvrr_stack+10689, dvrr_stack+10458, dvrr_stack+219, dvrr_stack+623, dvrr_stack+11280);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+35678, dvrr_stack+15543, dvrr_stack+10488, dvrr_stack+229, dvrr_stack+10689, dvrr_stack+11298);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+35778, dvrr_stack+15603, dvrr_stack+35678, dvrr_stack+259, dvrr_stack+15543, dvrr_stack+11334);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+35928, dvrr_stack+15703, dvrr_stack+35778, dvrr_stack+379, dvrr_stack+15603, dvrr_stack+11394);
 tmp = dvrr_stack + 35928;
 target_ptr = Libderiv->dvrr_classes[5][3];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10548, dvrr_stack+1512, dvrr_stack+3675, NULL, NULL, dvrr_stack+170);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+219, dvrr_stack+10719, dvrr_stack+10548, dvrr_stack+633, dvrr_stack+1512, dvrr_stack+10458);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+36138, dvrr_stack+15853, dvrr_stack+219, dvrr_stack+648, dvrr_stack+10719, dvrr_stack+10488);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+36288, dvrr_stack+15943, dvrr_stack+36138, dvrr_stack+693, dvrr_stack+15853, dvrr_stack+35678);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+36513, dvrr_stack+16093, dvrr_stack+36288, dvrr_stack+873, dvrr_stack+15943, dvrr_stack+35778);
 tmp = dvrr_stack + 36513;
 target_ptr = Libderiv->dvrr_classes[5][4];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+36828,dvrr_stack+36513,dvrr_stack+35928,21);


 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+633, dvrr_stack+3690, dvrr_stack+6743, NULL, NULL, dvrr_stack+3675);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+37458, dvrr_stack+16768, dvrr_stack+633, dvrr_stack+1527, dvrr_stack+3690, dvrr_stack+10548);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+37584, dvrr_stack+16831, dvrr_stack+37458, dvrr_stack+1548, dvrr_stack+16768, dvrr_stack+219);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+37794, dvrr_stack+16957, dvrr_stack+37584, dvrr_stack+1611, dvrr_stack+16831, dvrr_stack+36138);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+38109, dvrr_stack+17167, dvrr_stack+37794, dvrr_stack+1863, dvrr_stack+16957, dvrr_stack+36288);
 tmp = dvrr_stack + 38109;
 target_ptr = Libderiv->dvrr_classes[5][5];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+38550,dvrr_stack+38109,dvrr_stack+36513,21);


 /* compute (5 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+39495,dvrr_stack+38550,dvrr_stack+36828,21);


 /* compute (5 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,210,dvrr_stack+40755, dvrr_stack+39495, dvrr_stack+35928);

 /* compute (1 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1527, dvrr_stack+6764, dvrr_stack+3563, NULL, NULL, dvrr_stack+6743);

 /* compute (2 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+41385, dvrr_stack+19507, dvrr_stack+1527, dvrr_stack+3711, dvrr_stack+6764, dvrr_stack+633);

 /* compute (3 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+41553, dvrr_stack+19591, dvrr_stack+41385, dvrr_stack+3739, dvrr_stack+19507, dvrr_stack+37458);

 /* compute (4 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+41833, dvrr_stack+19759, dvrr_stack+41553, dvrr_stack+3823, dvrr_stack+19591, dvrr_stack+37584);

 /* compute (5 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+42253, dvrr_stack+20039, dvrr_stack+41833, dvrr_stack+4159, dvrr_stack+19759, dvrr_stack+37794);
 tmp = dvrr_stack + 42253;
 target_ptr = Libderiv->dvrr_classes[5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+42841,dvrr_stack+42253,dvrr_stack+38109,21);


 /* compute (5 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+44164,dvrr_stack+42841,dvrr_stack+38550,21);


 /* compute (5 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,315,dvrr_stack+46054, dvrr_stack+44164, dvrr_stack+36513);

 /* compute (1 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3711, dvrr_stack+6599, dvrr_stack+509, NULL, NULL, dvrr_stack+3563);

 /* compute (2 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+46999, dvrr_stack+23429, dvrr_stack+3711, dvrr_stack+6792, dvrr_stack+6599, dvrr_stack+1527);

 /* compute (3 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+47215, dvrr_stack+23537, dvrr_stack+46999, dvrr_stack+6828, dvrr_stack+23429, dvrr_stack+41385);

 /* compute (4 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+47575, dvrr_stack+23753, dvrr_stack+47215, dvrr_stack+6936, dvrr_stack+23537, dvrr_stack+41553);

 /* compute (5 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+48115, dvrr_stack+24113, dvrr_stack+47575, dvrr_stack+7368, dvrr_stack+23753, dvrr_stack+41833);

 /* compute (5 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+48871,dvrr_stack+48115,dvrr_stack+42253,21);


 /* compute (5 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+50635,dvrr_stack+48871,dvrr_stack+42841,21);


 /* compute (5 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,441,dvrr_stack+53281, dvrr_stack+50635, dvrr_stack+38109);

 /* compute (0 0 | 1 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+12, Data->F+13, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6792, dvrr_stack+6419, dvrr_stack+0, Data->F+11, Data->F+12, NULL);

 /* compute (0 0 | 3 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+309, dvrr_stack+15, dvrr_stack+6792, dvrr_stack+6455, dvrr_stack+6419, NULL);

 /* compute (0 0 | 4 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6443, dvrr_stack+111, dvrr_stack+309, dvrr_stack+1443, dvrr_stack+15, NULL);

 /* compute (0 0 | 5 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6798, dvrr_stack+479, dvrr_stack+6443, dvrr_stack+6458, dvrr_stack+111, NULL);

 /* compute (0 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6819, dvrr_stack+6422, dvrr_stack+6798, dvrr_stack+6468, dvrr_stack+479, NULL);

 /* compute (0 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6847, dvrr_stack+3423, dvrr_stack+6819, dvrr_stack+3451, dvrr_stack+6422, NULL);

 /* compute (0 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6635, dvrr_stack+509, dvrr_stack+6847, dvrr_stack+3563, dvrr_stack+3423, NULL);

 /* compute (1 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6883, dvrr_stack+6483, dvrr_stack+6635, NULL, NULL, dvrr_stack+509);

 /* compute (2 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+54604, dvrr_stack+28748, dvrr_stack+6883, dvrr_stack+10773, dvrr_stack+6483, dvrr_stack+3711);

 /* compute (3 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+54874, dvrr_stack+28883, dvrr_stack+54604, dvrr_stack+10818, dvrr_stack+28748, dvrr_stack+46999);

 /* compute (4 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+55324, dvrr_stack+29153, dvrr_stack+54874, dvrr_stack+10953, dvrr_stack+28883, dvrr_stack+47215);

 /* compute (5 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+55999, dvrr_stack+29603, dvrr_stack+55324, dvrr_stack+11493, dvrr_stack+29153, dvrr_stack+47575);

 /* compute (5 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+56944,dvrr_stack+55999,dvrr_stack+48115,21);


 /* compute (5 0 | 6 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_id(Libderiv->CD,dvrr_stack+59212,dvrr_stack+56944,dvrr_stack+48871,21);


 /* compute (5 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,588,dvrr_stack+62740, dvrr_stack+59212, dvrr_stack+42253);

 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+1443, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+7018, dvrr_stack+11223, dvrr_stack+1443, Data->F+4, Data->F+5, NULL);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+7024, dvrr_stack+11226, dvrr_stack+7018, dvrr_stack+6680, dvrr_stack+11223, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+11484, dvrr_stack+210, dvrr_stack+614, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+7034, dvrr_stack+545, dvrr_stack+11484, dvrr_stack+30, dvrr_stack+210, dvrr_stack+1443);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+7052, dvrr_stack+11232, dvrr_stack+7034, dvrr_stack+6683, dvrr_stack+545, dvrr_stack+7018);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+7082, dvrr_stack+11250, dvrr_stack+7052, dvrr_stack+6692, dvrr_stack+11232, dvrr_stack+7024);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+30, dvrr_stack+617, dvrr_stack+124, NULL, NULL, dvrr_stack+614);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3819, dvrr_stack+11280, dvrr_stack+30, dvrr_stack+213, dvrr_stack+617, dvrr_stack+11484);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+3855, dvrr_stack+11298, dvrr_stack+3819, dvrr_stack+6710, dvrr_stack+11280, dvrr_stack+7034);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+1611, dvrr_stack+11334, dvrr_stack+3855, dvrr_stack+10593, dvrr_stack+11298, dvrr_stack+7052);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+10764, dvrr_stack+11394, dvrr_stack+1611, dvrr_stack+10629, dvrr_stack+11334, dvrr_stack+7082);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+10593, dvrr_stack+170, dvrr_stack+554, NULL, NULL, dvrr_stack+124);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+3915, dvrr_stack+10458, dvrr_stack+10593, dvrr_stack+623, dvrr_stack+170, dvrr_stack+30);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+10890, dvrr_stack+10488, dvrr_stack+3915, dvrr_stack+10689, dvrr_stack+10458, dvrr_stack+3819);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10990, dvrr_stack+35678, dvrr_stack+10890, dvrr_stack+15543, dvrr_stack+10488, dvrr_stack+3855);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+64504, dvrr_stack+35778, dvrr_stack+10990, dvrr_stack+15603, dvrr_stack+35678, dvrr_stack+1611);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+64714, dvrr_stack+35928, dvrr_stack+64504, dvrr_stack+15703, dvrr_stack+35778, dvrr_stack+10764);
 tmp = dvrr_stack + 64714;
 target_ptr = Libderiv->dvrr_classes[6][3];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15543, dvrr_stack+3675, dvrr_stack+1428, NULL, NULL, dvrr_stack+554);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15588, dvrr_stack+10548, dvrr_stack+15543, dvrr_stack+1512, dvrr_stack+3675, dvrr_stack+10593);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+64994, dvrr_stack+219, dvrr_stack+15588, dvrr_stack+10719, dvrr_stack+10548, dvrr_stack+3915);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+65144, dvrr_stack+36138, dvrr_stack+64994, dvrr_stack+15853, dvrr_stack+219, dvrr_stack+10890);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+65369, dvrr_stack+36288, dvrr_stack+65144, dvrr_stack+15943, dvrr_stack+36138, dvrr_stack+10990);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+65684, dvrr_stack+36513, dvrr_stack+65369, dvrr_stack+16093, dvrr_stack+36288, dvrr_stack+64504);
 tmp = dvrr_stack + 65684;
 target_ptr = Libderiv->dvrr_classes[6][4];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+66104,dvrr_stack+65684,dvrr_stack+64714,28);


 /* compute (1 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15853, dvrr_stack+6743, dvrr_stack+3451, NULL, NULL, dvrr_stack+1428);

 /* compute (2 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15916, dvrr_stack+633, dvrr_stack+15853, dvrr_stack+3690, dvrr_stack+6743, dvrr_stack+15543);

 /* compute (3 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+66944, dvrr_stack+37458, dvrr_stack+15916, dvrr_stack+16768, dvrr_stack+633, dvrr_stack+15588);

 /* compute (4 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+67154, dvrr_stack+37584, dvrr_stack+66944, dvrr_stack+16831, dvrr_stack+37458, dvrr_stack+64994);

 /* compute (5 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+67469, dvrr_stack+37794, dvrr_stack+67154, dvrr_stack+16957, dvrr_stack+37584, dvrr_stack+65144);

 /* compute (6 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+67910, dvrr_stack+38109, dvrr_stack+67469, dvrr_stack+17167, dvrr_stack+37794, dvrr_stack+65369);
 tmp = dvrr_stack + 67910;
 target_ptr = Libderiv->dvrr_classes[6][5];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+68498,dvrr_stack+67910,dvrr_stack+65684,28);


 /* compute (6 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+69758,dvrr_stack+68498,dvrr_stack+66104,28);


 /* compute (6 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,280,dvrr_stack+71438, dvrr_stack+69758, dvrr_stack+64714);

 /* compute (1 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16768, dvrr_stack+3563, dvrr_stack+3423, NULL, NULL, dvrr_stack+3451);

 /* compute (2 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16852, dvrr_stack+1527, dvrr_stack+16768, dvrr_stack+6764, dvrr_stack+3563, dvrr_stack+15853);

 /* compute (3 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+72278, dvrr_stack+41385, dvrr_stack+16852, dvrr_stack+19507, dvrr_stack+1527, dvrr_stack+15916);

 /* compute (4 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+72558, dvrr_stack+41553, dvrr_stack+72278, dvrr_stack+19591, dvrr_stack+41385, dvrr_stack+66944);

 /* compute (5 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+72978, dvrr_stack+41833, dvrr_stack+72558, dvrr_stack+19759, dvrr_stack+41553, dvrr_stack+67154);

 /* compute (6 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+73566, dvrr_stack+42253, dvrr_stack+72978, dvrr_stack+20039, dvrr_stack+41833, dvrr_stack+67469);

 /* compute (6 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+74350,dvrr_stack+73566,dvrr_stack+67910,28);


 /* compute (6 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+76114,dvrr_stack+74350,dvrr_stack+68498,28);


 /* compute (6 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,420,dvrr_stack+78634, dvrr_stack+76114, dvrr_stack+65684);

 /* compute (1 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+19507, dvrr_stack+509, dvrr_stack+6847, NULL, NULL, dvrr_stack+3423);

 /* compute (2 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+19615, dvrr_stack+3711, dvrr_stack+19507, dvrr_stack+6599, dvrr_stack+509, dvrr_stack+16768);

 /* compute (3 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+79894, dvrr_stack+46999, dvrr_stack+19615, dvrr_stack+23429, dvrr_stack+3711, dvrr_stack+16852);

 /* compute (4 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+80254, dvrr_stack+47215, dvrr_stack+79894, dvrr_stack+23537, dvrr_stack+46999, dvrr_stack+72278);

 /* compute (5 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+80794, dvrr_stack+47575, dvrr_stack+80254, dvrr_stack+23753, dvrr_stack+47215, dvrr_stack+72558);

 /* compute (6 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+81550, dvrr_stack+48115, dvrr_stack+80794, dvrr_stack+24113, dvrr_stack+47575, dvrr_stack+72978);

 /* compute (6 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+82558,dvrr_stack+81550,dvrr_stack+73566,28);


 /* compute (6 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+84910,dvrr_stack+82558,dvrr_stack+74350,28);


 /* compute (6 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,588,dvrr_stack+88438, dvrr_stack+84910, dvrr_stack+67910);

 /* compute (0 0 | 1 0) m=13 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1446, Data->F+13, Data->F+14, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+10623, dvrr_stack+0, dvrr_stack+1446, Data->F+12, Data->F+13, NULL);

 /* compute (0 0 | 3 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+623, dvrr_stack+6792, dvrr_stack+10623, dvrr_stack+6419, dvrr_stack+0, NULL);

 /* compute (0 0 | 4 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1512, dvrr_stack+309, dvrr_stack+623, dvrr_stack+15, dvrr_stack+6792, NULL);

 /* compute (0 0 | 5 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3690, dvrr_stack+6443, dvrr_stack+1512, dvrr_stack+111, dvrr_stack+309, NULL);

 /* compute (0 0 | 6 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+23429, dvrr_stack+6798, dvrr_stack+3690, dvrr_stack+479, dvrr_stack+6443, NULL);

 /* compute (0 0 | 7 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1701, dvrr_stack+6819, dvrr_stack+23429, dvrr_stack+6422, dvrr_stack+6798, NULL);

 /* compute (0 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+23457, dvrr_stack+6847, dvrr_stack+1701, dvrr_stack+3423, dvrr_stack+6819, NULL);

 /* compute (1 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+23502, dvrr_stack+6635, dvrr_stack+23457, NULL, NULL, dvrr_stack+6847);

 /* compute (2 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+23637, dvrr_stack+6883, dvrr_stack+23502, dvrr_stack+6483, dvrr_stack+6635, dvrr_stack+19507);

 /* compute (3 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+90202, dvrr_stack+54604, dvrr_stack+23637, dvrr_stack+28748, dvrr_stack+6883, dvrr_stack+19615);

 /* compute (4 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+90652, dvrr_stack+54874, dvrr_stack+90202, dvrr_stack+28883, dvrr_stack+54604, dvrr_stack+79894);

 /* compute (5 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+91327, dvrr_stack+55324, dvrr_stack+90652, dvrr_stack+29153, dvrr_stack+54874, dvrr_stack+80254);

 /* compute (6 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+92272, dvrr_stack+55999, dvrr_stack+91327, dvrr_stack+29603, dvrr_stack+55324, dvrr_stack+80794);

 /* compute (6 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+93532,dvrr_stack+92272,dvrr_stack+81550,28);


 /* compute (6 0 | 6 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_id(Libderiv->CD,dvrr_stack+96556,dvrr_stack+93532,dvrr_stack+82558,28);


 /* compute (6 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,784,dvrr_stack+101260, dvrr_stack+96556, dvrr_stack+73566);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,100,dvrr_stack+54604, dvrr_stack+2523, dvrr_stack+379);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+54904, dvrr_stack+5069, dvrr_stack+873);

 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,210,dvrr_stack+55354, dvrr_stack+8568, dvrr_stack+1863);

 /* compute (3 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,280,dvrr_stack+90202, dvrr_stack+13023, dvrr_stack+4159);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+91042, dvrr_stack+18157, dvrr_stack+15703);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,225,dvrr_stack+91492, dvrr_stack+21404, dvrr_stack+16093);

 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,315,dvrr_stack+103612, dvrr_stack+25913, dvrr_stack+17167);

 /* compute (4 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,420,dvrr_stack+104557, dvrr_stack+31898, dvrr_stack+20039);

 /* compute (5 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,210,dvrr_stack+28748, dvrr_stack+39495, dvrr_stack+35928);

 /* compute (5 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,315,dvrr_stack+105817, dvrr_stack+44164, dvrr_stack+36513);

 /* compute (5 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,441,dvrr_stack+106762, dvrr_stack+50635, dvrr_stack+38109);

 /* compute (5 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,588,dvrr_stack+108085, dvrr_stack+59212, dvrr_stack+42253);

 /* compute (6 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,280,dvrr_stack+109849, dvrr_stack+69758, dvrr_stack+64714);

 /* compute (6 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,420,dvrr_stack+110689, dvrr_stack+76114, dvrr_stack+65684);

 /* compute (6 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,588,dvrr_stack+111949, dvrr_stack+84910, dvrr_stack+67910);

 /* compute (6 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,784,dvrr_stack+113713, dvrr_stack+96556, dvrr_stack+73566);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,100,dvrr_stack+23457, dvrr_stack+2523, dvrr_stack+379);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+2523, dvrr_stack+5069, dvrr_stack+873);

 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,210,dvrr_stack+5069, dvrr_stack+8568, dvrr_stack+1863);

 /* compute (3 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,280,dvrr_stack+8568, dvrr_stack+13023, dvrr_stack+4159);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+13023, dvrr_stack+18157, dvrr_stack+15703);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,225,dvrr_stack+18157, dvrr_stack+21404, dvrr_stack+16093);

 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,315,dvrr_stack+21404, dvrr_stack+25913, dvrr_stack+17167);

 /* compute (4 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,420,dvrr_stack+25913, dvrr_stack+31898, dvrr_stack+20039);

 /* compute (5 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,210,dvrr_stack+27173, dvrr_stack+39495, dvrr_stack+35928);

 /* compute (5 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,315,dvrr_stack+39495, dvrr_stack+44164, dvrr_stack+36513);

 /* compute (5 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,441,dvrr_stack+44164, dvrr_stack+50635, dvrr_stack+38109);

 /* compute (5 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,588,dvrr_stack+50635, dvrr_stack+59212, dvrr_stack+42253);

 /* compute (6 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,280,dvrr_stack+59212, dvrr_stack+69758, dvrr_stack+64714);

 /* compute (6 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,420,dvrr_stack+69758, dvrr_stack+76114, dvrr_stack+65684);

 /* compute (6 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,588,dvrr_stack+76114, dvrr_stack+84910, dvrr_stack+67910);

 /* compute (6 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,784,dvrr_stack+84910, dvrr_stack+96556, dvrr_stack+73566);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+6683, dvrr_stack+21, dvrr_stack+3, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+96556, dvrr_stack+6683, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+0);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+96574, dvrr_stack+164, dvrr_stack+24, NULL, NULL, dvrr_stack+21);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+96592, dvrr_stack+96574, dvrr_stack+57, dvrr_stack+164, dvrr_stack+24, dvrr_stack+6683);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+96628, dvrr_stack+96592, dvrr_stack+75, dvrr_stack+96574, dvrr_stack+57, dvrr_stack+96556);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+96688,dvrr_stack+379,dvrr_stack+96628,10);


 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+96868, dvrr_stack+96688, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+2973, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 2973;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+96928, dvrr_stack+1023, NULL);
 tmp = dvrr_stack + 96928;
 target_ptr = Libderiv->deriv_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+97028, dvrr_stack+4439, NULL);
 tmp = dvrr_stack + 97028;
 target_ptr = Libderiv->deriv_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,280,dvrr_stack+97238, dvrr_stack+7728, NULL);
 tmp = dvrr_stack + 97238;
 target_ptr = Libderiv->deriv_classes[3][6][11];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,360,dvrr_stack+97518, dvrr_stack+11943, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+24, dvrr_stack+0, dvrr_stack+6680, Data->F+2, Data->F+3, NULL);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+97878, dvrr_stack+96556, dvrr_stack+6692, dvrr_stack+6683, dvrr_stack+6, dvrr_stack+24);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+97908, dvrr_stack+96628, dvrr_stack+10629, dvrr_stack+96592, dvrr_stack+75, dvrr_stack+97878);

 /* compute (4 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+5699,dvrr_stack+15703,dvrr_stack+97908,15);


 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+97998, dvrr_stack+5699, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+18832, dvrr_stack+17482, NULL);
 tmp = dvrr_stack + 18832;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+98088, dvrr_stack+16318, NULL);
 tmp = dvrr_stack + 98088;
 target_ptr = Libderiv->deriv_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+40440, dvrr_stack+20459, NULL);
 tmp = dvrr_stack + 40440;
 target_ptr = Libderiv->deriv_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+71018, dvrr_stack+24653, NULL);
 tmp = dvrr_stack + 71018;
 target_ptr = Libderiv->deriv_classes[4][6][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,540,dvrr_stack+98238, dvrr_stack+30278, NULL);

 /* compute (3 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+309, dvrr_stack+24, dvrr_stack+11226, dvrr_stack+0, dvrr_stack+6680, NULL);

 /* compute (4 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+98778, dvrr_stack+97878, dvrr_stack+11250, dvrr_stack+96556, dvrr_stack+6692, dvrr_stack+309);

 /* compute (5 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+98823, dvrr_stack+97908, dvrr_stack+11394, dvrr_stack+96628, dvrr_stack+10629, dvrr_stack+98778);

 /* compute (5 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+98949,dvrr_stack+35928,dvrr_stack+98823,21);


 /* compute (5 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+99327, dvrr_stack+98949, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+99453, dvrr_stack+38550, NULL);
 tmp = dvrr_stack + 99453;
 target_ptr = Libderiv->deriv_classes[5][4][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+99768, dvrr_stack+36828, NULL);
 tmp = dvrr_stack + 99768;
 target_ptr = Libderiv->deriv_classes[5][3][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,441,dvrr_stack+99978, dvrr_stack+42841, NULL);
 tmp = dvrr_stack + 99978;
 target_ptr = Libderiv->deriv_classes[5][5][11];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+100419, dvrr_stack+48871, NULL);
 tmp = dvrr_stack + 100419;
 target_ptr = Libderiv->deriv_classes[5][6][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,756,dvrr_stack+77878, dvrr_stack+56944, NULL);

 /* compute (4 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+55984, dvrr_stack+309, dvrr_stack+7024, dvrr_stack+24, dvrr_stack+11226, NULL);

 /* compute (5 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+101007, dvrr_stack+98778, dvrr_stack+7082, dvrr_stack+97878, dvrr_stack+11250, dvrr_stack+55984);

 /* compute (6 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+101070, dvrr_stack+98823, dvrr_stack+10764, dvrr_stack+97908, dvrr_stack+11394, dvrr_stack+101007);

 /* compute (6 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+87262,dvrr_stack+64714,dvrr_stack+101070,28);


 /* compute (6 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,168,dvrr_stack+87766, dvrr_stack+87262, NULL);

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+9408, dvrr_stack+68498, NULL);
 tmp = dvrr_stack + 9408;
 target_ptr = Libderiv->deriv_classes[6][4][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,280,dvrr_stack+87934, dvrr_stack+66104, NULL);
 tmp = dvrr_stack + 87934;
 target_ptr = Libderiv->deriv_classes[6][3][11];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+60052, dvrr_stack+74350, NULL);
 tmp = dvrr_stack + 60052;
 target_ptr = Libderiv->deriv_classes[6][5][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,784,dvrr_stack+60640, dvrr_stack+82558, NULL);
 tmp = dvrr_stack + 60640;
 target_ptr = Libderiv->deriv_classes[6][6][11];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1008,dvrr_stack+61424, dvrr_stack+93532, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+88214, dvrr_stack+96688, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+88274, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 88274;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+62432, dvrr_stack+1023, NULL);
 tmp = dvrr_stack + 62432;
 target_ptr = Libderiv->deriv_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+52399, dvrr_stack+4439, NULL);
 tmp = dvrr_stack + 52399;
 target_ptr = Libderiv->deriv_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+52609, dvrr_stack+7728, NULL);
 tmp = dvrr_stack + 52609;
 target_ptr = Libderiv->deriv_classes[3][6][10];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,360,dvrr_stack+52889, dvrr_stack+11943, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+62532, dvrr_stack+5699, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+29378, dvrr_stack+17482, NULL);
 tmp = dvrr_stack + 29378;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+45487, dvrr_stack+16318, NULL);
 tmp = dvrr_stack + 45487;
 target_ptr = Libderiv->deriv_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+45637, dvrr_stack+20459, NULL);
 tmp = dvrr_stack + 45637;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+31898, dvrr_stack+24653, NULL);
 tmp = dvrr_stack + 31898;
 target_ptr = Libderiv->deriv_classes[4][6][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,540,dvrr_stack+32318, dvrr_stack+30278, NULL);

 /* compute (5 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+32858, dvrr_stack+98949, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+32984, dvrr_stack+38550, NULL);
 tmp = dvrr_stack + 32984;
 target_ptr = Libderiv->deriv_classes[5][4][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+33299, dvrr_stack+36828, NULL);
 tmp = dvrr_stack + 33299;
 target_ptr = Libderiv->deriv_classes[5][3][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+33509, dvrr_stack+42841, NULL);
 tmp = dvrr_stack + 33509;
 target_ptr = Libderiv->deriv_classes[5][5][10];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+13473, dvrr_stack+48871, NULL);
 tmp = dvrr_stack + 13473;
 target_ptr = Libderiv->deriv_classes[5][6][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,756,dvrr_stack+116065, dvrr_stack+56944, NULL);

 /* compute (6 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,168,dvrr_stack+33950, dvrr_stack+87262, NULL);

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+14061, dvrr_stack+68498, NULL);
 tmp = dvrr_stack + 14061;
 target_ptr = Libderiv->deriv_classes[6][4][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+34118, dvrr_stack+66104, NULL);
 tmp = dvrr_stack + 34118;
 target_ptr = Libderiv->deriv_classes[6][3][10];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+116821, dvrr_stack+74350, NULL);
 tmp = dvrr_stack + 116821;
 target_ptr = Libderiv->deriv_classes[6][5][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,784,dvrr_stack+117409, dvrr_stack+82558, NULL);
 tmp = dvrr_stack + 117409;
 target_ptr = Libderiv->deriv_classes[6][6][10];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1008,dvrr_stack+118193, dvrr_stack+93532, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+62622, dvrr_stack+96688, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+96688, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 96688;
 target_ptr = Libderiv->deriv_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+2073, dvrr_stack+1023, NULL);
 tmp = dvrr_stack + 2073;
 target_ptr = Libderiv->deriv_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+1023, dvrr_stack+4439, NULL);
 tmp = dvrr_stack + 1023;
 target_ptr = Libderiv->deriv_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+4439, dvrr_stack+7728, NULL);
 tmp = dvrr_stack + 4439;
 target_ptr = Libderiv->deriv_classes[3][6][9];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,360,dvrr_stack+7728, dvrr_stack+11943, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+11943, dvrr_stack+5699, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+5699, dvrr_stack+17482, NULL);
 tmp = dvrr_stack + 5699;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+17482, dvrr_stack+16318, NULL);
 tmp = dvrr_stack + 17482;
 target_ptr = Libderiv->deriv_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+16318, dvrr_stack+20459, NULL);
 tmp = dvrr_stack + 16318;
 target_ptr = Libderiv->deriv_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+20459, dvrr_stack+24653, NULL);
 tmp = dvrr_stack + 20459;
 target_ptr = Libderiv->deriv_classes[4][6][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,540,dvrr_stack+24653, dvrr_stack+30278, NULL);

 /* compute (5 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+30278, dvrr_stack+98949, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+98949, dvrr_stack+38550, NULL);
 tmp = dvrr_stack + 98949;
 target_ptr = Libderiv->deriv_classes[5][4][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+38550, dvrr_stack+36828, NULL);
 tmp = dvrr_stack + 38550;
 target_ptr = Libderiv->deriv_classes[5][3][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+36828, dvrr_stack+42841, NULL);
 tmp = dvrr_stack + 36828;
 target_ptr = Libderiv->deriv_classes[5][5][9];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+42841, dvrr_stack+48871, NULL);
 tmp = dvrr_stack + 42841;
 target_ptr = Libderiv->deriv_classes[5][6][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,756,dvrr_stack+48871, dvrr_stack+56944, NULL);

 /* compute (6 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,168,dvrr_stack+56944, dvrr_stack+87262, NULL);

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+87262, dvrr_stack+68498, NULL);
 tmp = dvrr_stack + 87262;
 target_ptr = Libderiv->deriv_classes[6][4][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+68498, dvrr_stack+66104, NULL);
 tmp = dvrr_stack + 68498;
 target_ptr = Libderiv->deriv_classes[6][3][9];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+66104, dvrr_stack+74350, NULL);
 tmp = dvrr_stack + 66104;
 target_ptr = Libderiv->deriv_classes[6][5][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,784,dvrr_stack+74350, dvrr_stack+82558, NULL);
 tmp = dvrr_stack + 74350;
 target_ptr = Libderiv->deriv_classes[6][6][9];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1008,dvrr_stack+49627, dvrr_stack+93532, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+6419, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+164, dvrr_stack+6419, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+93532, dvrr_stack+161, dvrr_stack+21, NULL, NULL, Data->F+1);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+96574, dvrr_stack+93532, dvrr_stack+6683, dvrr_stack+161, dvrr_stack+21, dvrr_stack+6419);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+96838, dvrr_stack+96574, dvrr_stack+96556, dvrr_stack+93532, dvrr_stack+6683, dvrr_stack+164);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+93532, dvrr_stack+379, dvrr_stack+96838);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+93592, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 93592;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+93742, dvrr_stack+873, dvrr_stack+96628);
 tmp = dvrr_stack + 93742;
 target_ptr = Libderiv->deriv_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+93842, dvrr_stack+4159, dvrr_stack+873);
 tmp = dvrr_stack + 93842;
 target_ptr = Libderiv->deriv_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+94052, dvrr_stack+7368, dvrr_stack+1863);
 tmp = dvrr_stack + 94052;
 target_ptr = Libderiv->deriv_classes[3][6][8];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,10,1,dvrr_stack+94332, dvrr_stack+11493, dvrr_stack+4159);

 /* compute (3 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+623, dvrr_stack+164, dvrr_stack+24, dvrr_stack+6419, dvrr_stack+0, NULL);

 /* compute (4 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+5924, dvrr_stack+96838, dvrr_stack+97878, dvrr_stack+96574, dvrr_stack+96556, dvrr_stack+623);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,15,1,dvrr_stack+94692, dvrr_stack+15703, dvrr_stack+5924);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+94782, dvrr_stack+17167, dvrr_stack+15703);
 tmp = dvrr_stack + 94782;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+95007, dvrr_stack+16093, dvrr_stack+97908);
 tmp = dvrr_stack + 95007;
 target_ptr = Libderiv->deriv_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+95157, dvrr_stack+20039, dvrr_stack+16093);
 tmp = dvrr_stack + 95157;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+95472, dvrr_stack+24113, dvrr_stack+17167);
 tmp = dvrr_stack + 95472;
 target_ptr = Libderiv->deriv_classes[4][6][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,15,1,dvrr_stack+95892, dvrr_stack+29603, dvrr_stack+20039);

 /* compute (4 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+6443, dvrr_stack+623, dvrr_stack+309, dvrr_stack+164, dvrr_stack+24, NULL);

 /* compute (5 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+99264, dvrr_stack+5924, dvrr_stack+98778, dvrr_stack+96838, dvrr_stack+97878, dvrr_stack+6443);

 /* compute (5 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,21,1,dvrr_stack+96432, dvrr_stack+35928, dvrr_stack+99264);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+82558, dvrr_stack+38109, dvrr_stack+35928);
 tmp = dvrr_stack + 82558;
 target_ptr = Libderiv->deriv_classes[5][4][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+82873, dvrr_stack+36513, dvrr_stack+98823);
 tmp = dvrr_stack + 82873;
 target_ptr = Libderiv->deriv_classes[5][3][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+83083, dvrr_stack+42253, dvrr_stack+36513);
 tmp = dvrr_stack + 83083;
 target_ptr = Libderiv->deriv_classes[5][5][8];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+83524, dvrr_stack+48115, dvrr_stack+38109);
 tmp = dvrr_stack + 83524;
 target_ptr = Libderiv->deriv_classes[5][6][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,21,1,dvrr_stack+84112, dvrr_stack+55999, dvrr_stack+42253);

 /* compute (5 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+3690, dvrr_stack+6443, dvrr_stack+55984, dvrr_stack+623, dvrr_stack+309, NULL);

 /* compute (6 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+87682, dvrr_stack+99264, dvrr_stack+101007, dvrr_stack+5924, dvrr_stack+98778, dvrr_stack+3690);

 /* compute (6 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,28,1,dvrr_stack+75134, dvrr_stack+64714, dvrr_stack+87682);

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+75302, dvrr_stack+67910, dvrr_stack+64714);
 tmp = dvrr_stack + 75302;
 target_ptr = Libderiv->deriv_classes[6][4][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,28,1,dvrr_stack+75722, dvrr_stack+65684, dvrr_stack+101070);
 tmp = dvrr_stack + 75722;
 target_ptr = Libderiv->deriv_classes[6][3][8];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+68778, dvrr_stack+73566, dvrr_stack+65684);
 tmp = dvrr_stack + 68778;
 target_ptr = Libderiv->deriv_classes[6][5][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+57112, dvrr_stack+81550, dvrr_stack+67910);
 tmp = dvrr_stack + 57112;
 target_ptr = Libderiv->deriv_classes[6][6][8];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,28,1,dvrr_stack+57896, dvrr_stack+92272, dvrr_stack+73566);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+76002, dvrr_stack+379, dvrr_stack+96838);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+66692, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 66692;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+66842, dvrr_stack+873, dvrr_stack+96628);
 tmp = dvrr_stack + 66842;
 target_ptr = Libderiv->deriv_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+69366, dvrr_stack+4159, dvrr_stack+873);
 tmp = dvrr_stack + 69366;
 target_ptr = Libderiv->deriv_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+58904, dvrr_stack+7368, dvrr_stack+1863);
 tmp = dvrr_stack + 58904;
 target_ptr = Libderiv->deriv_classes[3][6][7];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,10,1,dvrr_stack+43429, dvrr_stack+11493, dvrr_stack+4159);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,15,1,dvrr_stack+69576, dvrr_stack+15703, dvrr_stack+5924);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+43789, dvrr_stack+17167, dvrr_stack+15703);
 tmp = dvrr_stack + 43789;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+44014, dvrr_stack+16093, dvrr_stack+97908);
 tmp = dvrr_stack + 44014;
 target_ptr = Libderiv->deriv_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+38760, dvrr_stack+20039, dvrr_stack+16093);
 tmp = dvrr_stack + 38760;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+39075, dvrr_stack+24113, dvrr_stack+17167);
 tmp = dvrr_stack + 39075;
 target_ptr = Libderiv->deriv_classes[4][6][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,15,1,dvrr_stack+30404, dvrr_stack+29603, dvrr_stack+20039);

 /* compute (5 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,21,1,dvrr_stack+37269, dvrr_stack+35928, dvrr_stack+99264);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+30944, dvrr_stack+38109, dvrr_stack+35928);
 tmp = dvrr_stack + 30944;
 target_ptr = Libderiv->deriv_classes[5][4][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+31259, dvrr_stack+36513, dvrr_stack+98823);
 tmp = dvrr_stack + 31259;
 target_ptr = Libderiv->deriv_classes[5][3][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+25193, dvrr_stack+42253, dvrr_stack+36513);
 tmp = dvrr_stack + 25193;
 target_ptr = Libderiv->deriv_classes[5][5][7];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+12033, dvrr_stack+48115, dvrr_stack+38109);
 tmp = dvrr_stack + 12033;
 target_ptr = Libderiv->deriv_classes[5][6][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,21,1,dvrr_stack+119201, dvrr_stack+55999, dvrr_stack+42253);

 /* compute (6 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,28,1,dvrr_stack+31469, dvrr_stack+64714, dvrr_stack+87682);

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+20879, dvrr_stack+67910, dvrr_stack+64714);
 tmp = dvrr_stack + 20879;
 target_ptr = Libderiv->deriv_classes[6][4][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,28,1,dvrr_stack+17632, dvrr_stack+65684, dvrr_stack+101070);
 tmp = dvrr_stack + 17632;
 target_ptr = Libderiv->deriv_classes[6][3][7];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+119957, dvrr_stack+73566, dvrr_stack+65684);
 tmp = dvrr_stack + 119957;
 target_ptr = Libderiv->deriv_classes[6][5][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+120545, dvrr_stack+81550, dvrr_stack+67910);
 tmp = dvrr_stack + 120545;
 target_ptr = Libderiv->deriv_classes[6][6][7];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,28,1,dvrr_stack+121329, dvrr_stack+92272, dvrr_stack+73566);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+69666, dvrr_stack+379, dvrr_stack+96838);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+31637, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 31637;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+31787, dvrr_stack+873, dvrr_stack+96628);
 tmp = dvrr_stack + 31787;
 target_ptr = Libderiv->deriv_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+25634, dvrr_stack+4159, dvrr_stack+873);
 tmp = dvrr_stack + 25634;
 target_ptr = Libderiv->deriv_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+12621, dvrr_stack+7368, dvrr_stack+1863);
 tmp = dvrr_stack + 12621;
 target_ptr = Libderiv->deriv_classes[3][6][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,10,1,dvrr_stack+7368, dvrr_stack+11493, dvrr_stack+4159);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,15,1,dvrr_stack+11493, dvrr_stack+15703, dvrr_stack+5924);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+11583, dvrr_stack+17167, dvrr_stack+15703);
 tmp = dvrr_stack + 11583;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+17912, dvrr_stack+16093, dvrr_stack+97908);
 tmp = dvrr_stack + 17912;
 target_ptr = Libderiv->deriv_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+8088, dvrr_stack+20039, dvrr_stack+16093);
 tmp = dvrr_stack + 8088;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+122337, dvrr_stack+24113, dvrr_stack+17167);
 tmp = dvrr_stack + 122337;
 target_ptr = Libderiv->deriv_classes[4][6][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,15,1,dvrr_stack+23757, dvrr_stack+29603, dvrr_stack+20039);

 /* compute (5 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,21,1,dvrr_stack+29603, dvrr_stack+35928, dvrr_stack+99264);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+29729, dvrr_stack+38109, dvrr_stack+35928);
 tmp = dvrr_stack + 29729;
 target_ptr = Libderiv->deriv_classes[5][4][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+30044, dvrr_stack+36513, dvrr_stack+98823);
 tmp = dvrr_stack + 30044;
 target_ptr = Libderiv->deriv_classes[5][3][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+122757, dvrr_stack+42253, dvrr_stack+36513);
 tmp = dvrr_stack + 122757;
 target_ptr = Libderiv->deriv_classes[5][5][6];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+123198, dvrr_stack+48115, dvrr_stack+38109);
 tmp = dvrr_stack + 123198;
 target_ptr = Libderiv->deriv_classes[5][6][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,21,1,dvrr_stack+123786, dvrr_stack+55999, dvrr_stack+42253);

 /* compute (6 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,28,1,dvrr_stack+55999, dvrr_stack+64714, dvrr_stack+87682);

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+56167, dvrr_stack+67910, dvrr_stack+64714);
 tmp = dvrr_stack + 56167;
 target_ptr = Libderiv->deriv_classes[6][4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,28,1,dvrr_stack+56587, dvrr_stack+65684, dvrr_stack+101070);
 tmp = dvrr_stack + 56587;
 target_ptr = Libderiv->deriv_classes[6][3][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+124542, dvrr_stack+73566, dvrr_stack+65684);
 tmp = dvrr_stack + 124542;
 target_ptr = Libderiv->deriv_classes[6][5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+125130, dvrr_stack+81550, dvrr_stack+67910);
 tmp = dvrr_stack + 125130;
 target_ptr = Libderiv->deriv_classes[6][6][6];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,28,1,dvrr_stack+125914, dvrr_stack+92272, dvrr_stack+73566);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+4719,dvrr_stack+783,dvrr_stack+319,6);


 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+96628, dvrr_stack+4719, NULL);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+2173,dvrr_stack+1737,dvrr_stack+783,6);


 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+11808, dvrr_stack+2173, NULL);

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+22349,dvrr_stack+3991,dvrr_stack+1737,6);


 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+16633, dvrr_stack+22349, NULL);

 /* compute (2 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+92167,dvrr_stack+7152,dvrr_stack+3991,6);


 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,168,dvrr_stack+4899, dvrr_stack+92167, NULL);

 /* compute (1 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+6419, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+213, dvrr_stack+1443, dvrr_stack+6419, Data->F+5, Data->F+6, NULL);

 /* compute (3 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+623, dvrr_stack+7018, dvrr_stack+213, dvrr_stack+11223, dvrr_stack+1443, NULL);

 /* compute (4 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+6443, dvrr_stack+7024, dvrr_stack+623, dvrr_stack+11226, dvrr_stack+7018, NULL);

 /* compute (1 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+16759, dvrr_stack+614, dvrr_stack+121, NULL, NULL, Data->F+7);

 /* compute (2 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+87682, dvrr_stack+11484, dvrr_stack+16759, dvrr_stack+210, dvrr_stack+614, dvrr_stack+6419);

 /* compute (3 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+96838, dvrr_stack+7034, dvrr_stack+87682, dvrr_stack+545, dvrr_stack+11484, dvrr_stack+213);

 /* compute (4 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+5924, dvrr_stack+7052, dvrr_stack+96838, dvrr_stack+11232, dvrr_stack+7034, dvrr_stack+623);

 /* compute (5 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+99264, dvrr_stack+7082, dvrr_stack+5924, dvrr_stack+11250, dvrr_stack+7052, dvrr_stack+6443);

 /* compute (1 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+87700, dvrr_stack+124, dvrr_stack+497, NULL, NULL, dvrr_stack+121);

 /* compute (2 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+87718, dvrr_stack+30, dvrr_stack+87700, dvrr_stack+617, dvrr_stack+124, dvrr_stack+16759);

 /* compute (3 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+56867, dvrr_stack+3819, dvrr_stack+87718, dvrr_stack+11280, dvrr_stack+30, dvrr_stack+87682);

 /* compute (4 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+97878, dvrr_stack+3855, dvrr_stack+56867, dvrr_stack+11298, dvrr_stack+3819, dvrr_stack+96838);

 /* compute (5 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+8403, dvrr_stack+1611, dvrr_stack+97878, dvrr_stack+11334, dvrr_stack+3855, dvrr_stack+5924);

 /* compute (6 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+14481, dvrr_stack+10764, dvrr_stack+8403, dvrr_stack+11394, dvrr_stack+1611, dvrr_stack+99264);

 /* compute (1 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+97968, dvrr_stack+554, dvrr_stack+1347, NULL, NULL, dvrr_stack+497);

 /* compute (2 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+37395, dvrr_stack+10593, dvrr_stack+97968, dvrr_stack+170, dvrr_stack+554, dvrr_stack+87700);

 /* compute (3 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+21299, dvrr_stack+3915, dvrr_stack+37395, dvrr_stack+10458, dvrr_stack+10593, dvrr_stack+87718);

 /* compute (4 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+92671, dvrr_stack+10890, dvrr_stack+21299, dvrr_stack+10488, dvrr_stack+3915, dvrr_stack+56867);

 /* compute (5 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+92821, dvrr_stack+10990, dvrr_stack+92671, dvrr_stack+35678, dvrr_stack+10890, dvrr_stack+97878);

 /* compute (6 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+93031, dvrr_stack+64504, dvrr_stack+92821, dvrr_stack+35778, dvrr_stack+10990, dvrr_stack+8403);

 /* compute (7 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+126922, dvrr_stack+64714, dvrr_stack+93031, dvrr_stack+35928, dvrr_stack+64504, dvrr_stack+14481);

 /* compute (1 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11898, dvrr_stack+1428, dvrr_stack+6468, NULL, NULL, dvrr_stack+1347);

 /* compute (2 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10458, dvrr_stack+15543, dvrr_stack+11898, dvrr_stack+3675, dvrr_stack+1428, dvrr_stack+97968);

 /* compute (3 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+35678, dvrr_stack+15588, dvrr_stack+10458, dvrr_stack+10548, dvrr_stack+15543, dvrr_stack+37395);

 /* compute (4 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6483, dvrr_stack+64994, dvrr_stack+35678, dvrr_stack+219, dvrr_stack+15588, dvrr_stack+21299);

 /* compute (5 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+24297, dvrr_stack+65144, dvrr_stack+6483, dvrr_stack+36138, dvrr_stack+64994, dvrr_stack+92671);

 /* compute (6 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+127282, dvrr_stack+65369, dvrr_stack+24297, dvrr_stack+36288, dvrr_stack+65144, dvrr_stack+92821);

 /* compute (7 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+127702, dvrr_stack+65684, dvrr_stack+127282, dvrr_stack+36513, dvrr_stack+65369, dvrr_stack+93031);

 /* compute (7 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+128242,dvrr_stack+127702,dvrr_stack+126922,36);


 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,360,dvrr_stack+36138, dvrr_stack+128242, NULL);

 /* compute (1 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+219, dvrr_stack+3451, dvrr_stack+6422, NULL, NULL, dvrr_stack+6468);

 /* compute (2 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+93311, dvrr_stack+15853, dvrr_stack+219, dvrr_stack+6743, dvrr_stack+3451, dvrr_stack+11898);

 /* compute (3 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+11140, dvrr_stack+15916, dvrr_stack+93311, dvrr_stack+633, dvrr_stack+15853, dvrr_stack+10458);

 /* compute (4 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+129322, dvrr_stack+66944, dvrr_stack+11140, dvrr_stack+37458, dvrr_stack+15916, dvrr_stack+35678);

 /* compute (5 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+129637, dvrr_stack+67154, dvrr_stack+129322, dvrr_stack+37584, dvrr_stack+66944, dvrr_stack+6483);

 /* compute (6 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+130078, dvrr_stack+67469, dvrr_stack+129637, dvrr_stack+37794, dvrr_stack+67154, dvrr_stack+24297);

 /* compute (7 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+130666, dvrr_stack+67910, dvrr_stack+130078, dvrr_stack+38109, dvrr_stack+67469, dvrr_stack+127282);

 /* compute (7 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+131422,dvrr_stack+130666,dvrr_stack+127702,36);


 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,540,dvrr_stack+37455, dvrr_stack+131422, NULL);

 /* compute (1 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+633, dvrr_stack+3423, dvrr_stack+6819, NULL, NULL, dvrr_stack+6422);

 /* compute (2 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+19831, dvrr_stack+16768, dvrr_stack+633, dvrr_stack+3563, dvrr_stack+3423, dvrr_stack+219);

 /* compute (3 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+133042, dvrr_stack+16852, dvrr_stack+19831, dvrr_stack+1527, dvrr_stack+16768, dvrr_stack+93311);

 /* compute (4 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+133322, dvrr_stack+72278, dvrr_stack+133042, dvrr_stack+41385, dvrr_stack+16852, dvrr_stack+11140);

 /* compute (5 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+133742, dvrr_stack+72558, dvrr_stack+133322, dvrr_stack+41553, dvrr_stack+72278, dvrr_stack+129322);

 /* compute (6 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+134330, dvrr_stack+72978, dvrr_stack+133742, dvrr_stack+41833, dvrr_stack+72558, dvrr_stack+129637);

 /* compute (7 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+135114, dvrr_stack+73566, dvrr_stack+134330, dvrr_stack+42253, dvrr_stack+72978, dvrr_stack+130078);

 /* compute (7 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+136122,dvrr_stack+135114,dvrr_stack+130666,36);


 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,756,dvrr_stack+41385, dvrr_stack+136122, NULL);

 /* compute (1 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+42141, dvrr_stack+6847, dvrr_stack+1701, NULL, NULL, dvrr_stack+6819);

 /* compute (2 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+138390, dvrr_stack+19507, dvrr_stack+42141, dvrr_stack+509, dvrr_stack+6847, dvrr_stack+633);

 /* compute (3 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+138606, dvrr_stack+19615, dvrr_stack+138390, dvrr_stack+3711, dvrr_stack+19507, dvrr_stack+19831);

 /* compute (4 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+138966, dvrr_stack+79894, dvrr_stack+138606, dvrr_stack+46999, dvrr_stack+19615, dvrr_stack+133042);

 /* compute (5 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+139506, dvrr_stack+80254, dvrr_stack+138966, dvrr_stack+47215, dvrr_stack+79894, dvrr_stack+133322);

 /* compute (6 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+138390, dvrr_stack+80794, dvrr_stack+139506, dvrr_stack+47575, dvrr_stack+80254, dvrr_stack+133742);

 /* compute (7 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+139398, dvrr_stack+81550, dvrr_stack+138390, dvrr_stack+48115, dvrr_stack+80794, dvrr_stack+134330);

 /* compute (7 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+140694,dvrr_stack+139398,dvrr_stack+135114,36);


 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1008,dvrr_stack+138390, dvrr_stack+140694, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+79894, dvrr_stack+4719, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+6708, dvrr_stack+2173, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+79954, dvrr_stack+22349, NULL);

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,168,dvrr_stack+80080, dvrr_stack+92167, NULL);

 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,360,dvrr_stack+80248, dvrr_stack+128242, NULL);

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,540,dvrr_stack+80608, dvrr_stack+131422, NULL);

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,756,dvrr_stack+81148, dvrr_stack+136122, NULL);

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1008,dvrr_stack+46999, dvrr_stack+140694, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+81904, dvrr_stack+4719, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+4719, dvrr_stack+2173, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+2173, dvrr_stack+22349, NULL);

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,168,dvrr_stack+22349, dvrr_stack+92167, NULL);

 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,360,dvrr_stack+92167, dvrr_stack+128242, NULL);

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,540,dvrr_stack+128242, dvrr_stack+131422, NULL);

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,756,dvrr_stack+131422, dvrr_stack+136122, NULL);

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1008,dvrr_stack+136122, dvrr_stack+140694, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+140694, dvrr_stack+783, dvrr_stack+96592);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+4809, dvrr_stack+1737, dvrr_stack+319);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+140754, dvrr_stack+3991, dvrr_stack+783);

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+140880, dvrr_stack+7152, dvrr_stack+1737);

 /* compute (5 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+141048, dvrr_stack+55984, dvrr_stack+6443, dvrr_stack+309, dvrr_stack+7024, NULL);

 /* compute (6 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+141069, dvrr_stack+101007, dvrr_stack+99264, dvrr_stack+98778, dvrr_stack+7082, dvrr_stack+141048);

 /* compute (7 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+141153, dvrr_stack+101070, dvrr_stack+14481, dvrr_stack+98823, dvrr_stack+10764, dvrr_stack+141069);

 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,36,1,dvrr_stack+141369, dvrr_stack+127702, dvrr_stack+141153);

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,36,1,dvrr_stack+128782, dvrr_stack+130666, dvrr_stack+126922);

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,36,1,dvrr_stack+141729, dvrr_stack+135114, dvrr_stack+127702);

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,36,1,dvrr_stack+142485, dvrr_stack+139398, dvrr_stack+130666);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+141048, dvrr_stack+783, dvrr_stack+96592);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+98778, dvrr_stack+1737, dvrr_stack+319);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+101007, dvrr_stack+3991, dvrr_stack+783);

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+143493, dvrr_stack+7152, dvrr_stack+1737);

 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,36,1,dvrr_stack+137130, dvrr_stack+127702, dvrr_stack+141153);

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,36,1,dvrr_stack+137490, dvrr_stack+130666, dvrr_stack+126922);

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,36,1,dvrr_stack+132178, dvrr_stack+135114, dvrr_stack+127702);

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,36,1,dvrr_stack+143661, dvrr_stack+139398, dvrr_stack+130666);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+98868, dvrr_stack+783, dvrr_stack+96592);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+101133, dvrr_stack+1737, dvrr_stack+319);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+144669, dvrr_stack+3991, dvrr_stack+783);

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+144795, dvrr_stack+7152, dvrr_stack+1737);

 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,36,1,dvrr_stack+138030, dvrr_stack+127702, dvrr_stack+141153);

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,36,1,dvrr_stack+81964, dvrr_stack+130666, dvrr_stack+126922);

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,36,1,dvrr_stack+48007, dvrr_stack+135114, dvrr_stack+127702);

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,36,1,dvrr_stack+144963, dvrr_stack+139398, dvrr_stack+130666);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+139398, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+139458, dvrr_stack+35928, dvrr_stack+379);
 tmp = dvrr_stack + 139458;
 target_ptr = Libderiv->deriv_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+139608, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+139698, dvrr_stack+36513, dvrr_stack+873);
 tmp = dvrr_stack + 139698;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+139923, dvrr_stack+1863, dvrr_stack+1449);

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+140049, dvrr_stack+38109, dvrr_stack+1863);
 tmp = dvrr_stack + 140049;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+140364, dvrr_stack+4159, dvrr_stack+3591);

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+145971, dvrr_stack+42253, dvrr_stack+4159);
 tmp = dvrr_stack + 145971;
 target_ptr = Libderiv->deriv_classes[4][6][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+35828, dvrr_stack+15703, dvrr_stack+319);
 tmp = dvrr_stack + 35828;
 target_ptr = Libderiv->deriv_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+141108, dvrr_stack+64714, dvrr_stack+15703);
 tmp = dvrr_stack + 141108;
 target_ptr = Libderiv->deriv_classes[5][3][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+140532, dvrr_stack+16093, dvrr_stack+783);
 tmp = dvrr_stack + 140532;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+146391, dvrr_stack+65684, dvrr_stack+16093);
 tmp = dvrr_stack + 146391;
 target_ptr = Libderiv->deriv_classes[5][4][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+146706, dvrr_stack+17167, dvrr_stack+1737);
 tmp = dvrr_stack + 146706;
 target_ptr = Libderiv->deriv_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+146916, dvrr_stack+67910, dvrr_stack+17167);
 tmp = dvrr_stack + 146916;
 target_ptr = Libderiv->deriv_classes[5][5][2];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+147357, dvrr_stack+20039, dvrr_stack+3991);
 tmp = dvrr_stack + 147357;
 target_ptr = Libderiv->deriv_classes[3][6][2];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+147637, dvrr_stack+73566, dvrr_stack+20039);
 tmp = dvrr_stack + 147637;
 target_ptr = Libderiv->deriv_classes[5][6][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,10,dvrr_stack+148225, dvrr_stack+126922, dvrr_stack+35928);
 tmp = dvrr_stack + 148225;
 target_ptr = Libderiv->deriv_classes[6][3][2];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+148505, dvrr_stack+127702, dvrr_stack+36513);
 tmp = dvrr_stack + 148505;
 target_ptr = Libderiv->deriv_classes[6][4][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+148925, dvrr_stack+130666, dvrr_stack+38109);
 tmp = dvrr_stack + 148925;
 target_ptr = Libderiv->deriv_classes[6][5][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+149513, dvrr_stack+135114, dvrr_stack+42253);
 tmp = dvrr_stack + 149513;
 target_ptr = Libderiv->deriv_classes[6][6][2];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+210, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+617, dvrr_stack+6419, dvrr_stack+210, Data->F+6, Data->F+7, NULL);

 /* compute (3 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+140682, dvrr_stack+213, dvrr_stack+617, dvrr_stack+1443, dvrr_stack+6419, NULL);

 /* compute (4 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+55984, dvrr_stack+623, dvrr_stack+140682, dvrr_stack+7018, dvrr_stack+213, NULL);

 /* compute (5 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+98928, dvrr_stack+6443, dvrr_stack+55984, dvrr_stack+7024, dvrr_stack+623, NULL);

 /* compute (1 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+623, dvrr_stack+121, dvrr_stack+494, NULL, NULL, Data->F+8);

 /* compute (2 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+141318, dvrr_stack+16759, dvrr_stack+623, dvrr_stack+614, dvrr_stack+121, dvrr_stack+210);

 /* compute (3 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+0, dvrr_stack+87682, dvrr_stack+141318, dvrr_stack+11484, dvrr_stack+16759, dvrr_stack+617);

 /* compute (4 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+509, dvrr_stack+96838, dvrr_stack+0, dvrr_stack+7034, dvrr_stack+87682, dvrr_stack+140682);

 /* compute (5 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+150297, dvrr_stack+5924, dvrr_stack+509, dvrr_stack+7052, dvrr_stack+96838, dvrr_stack+55984);

 /* compute (6 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+150360, dvrr_stack+99264, dvrr_stack+150297, dvrr_stack+7082, dvrr_stack+5924, dvrr_stack+98928);

 /* compute (1 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+87682, dvrr_stack+497, dvrr_stack+503, NULL, NULL, dvrr_stack+494);

 /* compute (2 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1701, dvrr_stack+87700, dvrr_stack+87682, dvrr_stack+124, dvrr_stack+497, dvrr_stack+623);

 /* compute (3 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+99264, dvrr_stack+87718, dvrr_stack+1701, dvrr_stack+30, dvrr_stack+87700, dvrr_stack+141318);

 /* compute (4 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+30, dvrr_stack+56867, dvrr_stack+99264, dvrr_stack+3819, dvrr_stack+87718, dvrr_stack+0);

 /* compute (5 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+150444, dvrr_stack+97878, dvrr_stack+30, dvrr_stack+3855, dvrr_stack+56867, dvrr_stack+509);

 /* compute (6 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+150570, dvrr_stack+8403, dvrr_stack+150444, dvrr_stack+1611, dvrr_stack+97878, dvrr_stack+150297);

 /* compute (7 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+150738, dvrr_stack+14481, dvrr_stack+150570, dvrr_stack+10764, dvrr_stack+8403, dvrr_stack+150360);

 /* compute (1 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+0, dvrr_stack+1347, dvrr_stack+6458, NULL, NULL, dvrr_stack+503);

 /* compute (2 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+120, dvrr_stack+97968, dvrr_stack+0, dvrr_stack+554, dvrr_stack+1347, dvrr_stack+87682);

 /* compute (3 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+8403, dvrr_stack+37395, dvrr_stack+120, dvrr_stack+10593, dvrr_stack+97968, dvrr_stack+1701);

 /* compute (4 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+14481, dvrr_stack+21299, dvrr_stack+8403, dvrr_stack+3915, dvrr_stack+37395, dvrr_stack+99264);

 /* compute (5 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+22517, dvrr_stack+92671, dvrr_stack+14481, dvrr_stack+10890, dvrr_stack+21299, dvrr_stack+30);

 /* compute (6 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+19507, dvrr_stack+92821, dvrr_stack+22517, dvrr_stack+10990, dvrr_stack+92671, dvrr_stack+150444);

 /* compute (7 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+6847, dvrr_stack+93031, dvrr_stack+19507, dvrr_stack+64504, dvrr_stack+92821, dvrr_stack+150570);

 /* compute (8 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+92527, dvrr_stack+126922, dvrr_stack+6847, dvrr_stack+64714, dvrr_stack+93031, dvrr_stack+150738);

 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,10,dvrr_stack+150297, dvrr_stack+92527, dvrr_stack+64714);

 /* compute (1 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5924, dvrr_stack+6468, dvrr_stack+479, NULL, NULL, dvrr_stack+6458);

 /* compute (2 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+30, dvrr_stack+11898, dvrr_stack+5924, dvrr_stack+1428, dvrr_stack+6468, dvrr_stack+0);

 /* compute (3 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+64504, dvrr_stack+10458, dvrr_stack+30, dvrr_stack+15543, dvrr_stack+11898, dvrr_stack+120);

 /* compute (4 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1512, dvrr_stack+35678, dvrr_stack+64504, dvrr_stack+15588, dvrr_stack+10458, dvrr_stack+8403);

 /* compute (5 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10458, dvrr_stack+6483, dvrr_stack+1512, dvrr_stack+64994, dvrr_stack+35678, dvrr_stack+14481);

 /* compute (6 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+150657, dvrr_stack+24297, dvrr_stack+10458, dvrr_stack+65144, dvrr_stack+6483, dvrr_stack+22517);

 /* compute (7 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+151077, dvrr_stack+127282, dvrr_stack+150657, dvrr_stack+65369, dvrr_stack+24297, dvrr_stack+19507);

 /* compute (8 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+64994, dvrr_stack+127702, dvrr_stack+151077, dvrr_stack+65684, dvrr_stack+127282, dvrr_stack+6847);

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,15,dvrr_stack+151617, dvrr_stack+64994, dvrr_stack+65684);

 /* compute (1 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+99264, dvrr_stack+6422, dvrr_stack+6798, NULL, NULL, dvrr_stack+479);

 /* compute (2 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6847, dvrr_stack+219, dvrr_stack+99264, dvrr_stack+3451, dvrr_stack+6422, dvrr_stack+5924);

 /* compute (3 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6973, dvrr_stack+93311, dvrr_stack+6847, dvrr_stack+15853, dvrr_stack+219, dvrr_stack+30);

 /* compute (4 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+127282, dvrr_stack+11140, dvrr_stack+6973, dvrr_stack+15916, dvrr_stack+93311, dvrr_stack+64504);

 /* compute (5 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+152157, dvrr_stack+129322, dvrr_stack+127282, dvrr_stack+66944, dvrr_stack+11140, dvrr_stack+1512);

 /* compute (6 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+10773, dvrr_stack+129637, dvrr_stack+152157, dvrr_stack+67154, dvrr_stack+129322, dvrr_stack+10458);

 /* compute (7 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+152598, dvrr_stack+130078, dvrr_stack+10773, dvrr_stack+67469, dvrr_stack+129637, dvrr_stack+150657);

 /* compute (8 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+66942, dvrr_stack+130666, dvrr_stack+152598, dvrr_stack+67910, dvrr_stack+130078, dvrr_stack+151077);

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,21,dvrr_stack+150657, dvrr_stack+66942, dvrr_stack+67910);

 /* compute (1 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+87682, dvrr_stack+6819, dvrr_stack+23429, NULL, NULL, dvrr_stack+6798);

 /* compute (2 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+151413, dvrr_stack+633, dvrr_stack+87682, dvrr_stack+3423, dvrr_stack+6819, dvrr_stack+99264);

 /* compute (3 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10458, dvrr_stack+19831, dvrr_stack+151413, dvrr_stack+16768, dvrr_stack+633, dvrr_stack+6847);

 /* compute (4 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+129322, dvrr_stack+133042, dvrr_stack+10458, dvrr_stack+16852, dvrr_stack+19831, dvrr_stack+6973);

 /* compute (5 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+129742, dvrr_stack+133322, dvrr_stack+129322, dvrr_stack+72278, dvrr_stack+133042, dvrr_stack+127282);

 /* compute (6 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+153354, dvrr_stack+133742, dvrr_stack+129742, dvrr_stack+72558, dvrr_stack+133322, dvrr_stack+152157);

 /* compute (7 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+129322, dvrr_stack+134330, dvrr_stack+153354, dvrr_stack+72978, dvrr_stack+133742, dvrr_stack+10773);

 /* compute (8 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+153354, dvrr_stack+135114, dvrr_stack+129322, dvrr_stack+73566, dvrr_stack+134330, dvrr_stack+152598);

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,28,dvrr_stack+129322, dvrr_stack+153354, dvrr_stack+73566);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+37395, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+35678, dvrr_stack+35928, dvrr_stack+379);
 tmp = dvrr_stack + 35678;
 target_ptr = Libderiv->deriv_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+479, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+1512, dvrr_stack+36513, dvrr_stack+873);
 tmp = dvrr_stack + 1512;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+130330, dvrr_stack+1863, dvrr_stack+1449);

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+154614, dvrr_stack+38109, dvrr_stack+1863);
 tmp = dvrr_stack + 154614;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+3423, dvrr_stack+4159, dvrr_stack+3591);

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+127282, dvrr_stack+42253, dvrr_stack+4159);
 tmp = dvrr_stack + 127282;
 target_ptr = Libderiv->deriv_classes[4][6][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+130456, dvrr_stack+15703, dvrr_stack+319);
 tmp = dvrr_stack + 130456;
 target_ptr = Libderiv->deriv_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+64504, dvrr_stack+64714, dvrr_stack+15703);
 tmp = dvrr_stack + 64504;
 target_ptr = Libderiv->deriv_classes[5][3][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+152157, dvrr_stack+16093, dvrr_stack+783);
 tmp = dvrr_stack + 152157;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+152307, dvrr_stack+65684, dvrr_stack+16093);
 tmp = dvrr_stack + 152307;
 target_ptr = Libderiv->deriv_classes[5][4][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+152622, dvrr_stack+17167, dvrr_stack+1737);
 tmp = dvrr_stack + 152622;
 target_ptr = Libderiv->deriv_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+152832, dvrr_stack+67910, dvrr_stack+17167);
 tmp = dvrr_stack + 152832;
 target_ptr = Libderiv->deriv_classes[5][5][1];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+72278, dvrr_stack+20039, dvrr_stack+3991);
 tmp = dvrr_stack + 72278;
 target_ptr = Libderiv->deriv_classes[3][6][1];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+72558, dvrr_stack+73566, dvrr_stack+20039);
 tmp = dvrr_stack + 72558;
 target_ptr = Libderiv->deriv_classes[5][6][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,10,dvrr_stack+73146, dvrr_stack+126922, dvrr_stack+35928);
 tmp = dvrr_stack + 73146;
 target_ptr = Libderiv->deriv_classes[6][3][1];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+10458, dvrr_stack+127702, dvrr_stack+36513);
 tmp = dvrr_stack + 10458;
 target_ptr = Libderiv->deriv_classes[6][4][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+10878, dvrr_stack+130666, dvrr_stack+38109);
 tmp = dvrr_stack + 10878;
 target_ptr = Libderiv->deriv_classes[6][5][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+132934, dvrr_stack+135114, dvrr_stack+42253);
 tmp = dvrr_stack + 132934;
 target_ptr = Libderiv->deriv_classes[6][6][1];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,10,dvrr_stack+6798, dvrr_stack+92527, dvrr_stack+64714);

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,15,dvrr_stack+133718, dvrr_stack+64994, dvrr_stack+65684);

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,21,dvrr_stack+134258, dvrr_stack+66942, dvrr_stack+67910);

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,28,dvrr_stack+154929, dvrr_stack+153354, dvrr_stack+73566);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+130556, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+151413, dvrr_stack+35928, dvrr_stack+379);
 tmp = dvrr_stack + 151413;
 target_ptr = Libderiv->deriv_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+379, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+15853, dvrr_stack+36513, dvrr_stack+873);
 tmp = dvrr_stack + 15853;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+873, dvrr_stack+1863, dvrr_stack+1449);

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+19507, dvrr_stack+38109, dvrr_stack+1863);
 tmp = dvrr_stack + 19507;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+1863, dvrr_stack+4159, dvrr_stack+3591);

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+92977, dvrr_stack+42253, dvrr_stack+4159);
 tmp = dvrr_stack + 92977;
 target_ptr = Libderiv->deriv_classes[4][6][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+135014, dvrr_stack+15703, dvrr_stack+319);
 tmp = dvrr_stack + 135014;
 target_ptr = Libderiv->deriv_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+7158, dvrr_stack+64714, dvrr_stack+15703);
 tmp = dvrr_stack + 7158;
 target_ptr = Libderiv->deriv_classes[5][3][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+4159, dvrr_stack+16093, dvrr_stack+783);
 tmp = dvrr_stack + 4159;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+3591, dvrr_stack+65684, dvrr_stack+16093);
 tmp = dvrr_stack + 3591;
 target_ptr = Libderiv->deriv_classes[5][4][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+569, dvrr_stack+17167, dvrr_stack+1737);
 tmp = dvrr_stack + 569;
 target_ptr = Libderiv->deriv_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+155937, dvrr_stack+67910, dvrr_stack+17167);
 tmp = dvrr_stack + 155937;
 target_ptr = Libderiv->deriv_classes[5][5][0];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+156378, dvrr_stack+20039, dvrr_stack+3991);
 tmp = dvrr_stack + 156378;
 target_ptr = Libderiv->deriv_classes[3][6][0];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+16759, dvrr_stack+73566, dvrr_stack+20039);
 tmp = dvrr_stack + 16759;
 target_ptr = Libderiv->deriv_classes[5][6][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+156658, dvrr_stack+126922, dvrr_stack+35928);
 tmp = dvrr_stack + 156658;
 target_ptr = Libderiv->deriv_classes[6][3][0];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+19822, dvrr_stack+127702, dvrr_stack+36513);
 tmp = dvrr_stack + 19822;
 target_ptr = Libderiv->deriv_classes[6][4][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+156938, dvrr_stack+130666, dvrr_stack+38109);
 tmp = dvrr_stack + 156938;
 target_ptr = Libderiv->deriv_classes[6][5][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+130616, dvrr_stack+135114, dvrr_stack+42253);
 tmp = dvrr_stack + 130616;
 target_ptr = Libderiv->deriv_classes[6][6][0];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,10,dvrr_stack+126922, dvrr_stack+92527, dvrr_stack+64714);

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,15,dvrr_stack+127702, dvrr_stack+64994, dvrr_stack+65684);

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,21,dvrr_stack+64714, dvrr_stack+66942, dvrr_stack+67910);

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,28,dvrr_stack+135114, dvrr_stack+153354, dvrr_stack+73566);

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,100,dvrr_stack+66942, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 66942;
 target_ptr = Libderiv->deriv2_classes[3][3][143];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+67042, dvrr_stack+5969, NULL);
 tmp = dvrr_stack + 67042;
 target_ptr = Libderiv->deriv2_classes[3][4][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,210,dvrr_stack+35928, dvrr_stack+9828, NULL);
 tmp = dvrr_stack + 35928;
 target_ptr = Libderiv->deriv2_classes[3][5][143];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,280,dvrr_stack+67192, dvrr_stack+14703, NULL);
 tmp = dvrr_stack + 67192;
 target_ptr = Libderiv->deriv2_classes[3][6][143];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+67472, dvrr_stack+19057, NULL);
 tmp = dvrr_stack + 67472;
 target_ptr = Libderiv->deriv2_classes[4][3][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,225,dvrr_stack+67622, dvrr_stack+22754, NULL);
 tmp = dvrr_stack + 67622;
 target_ptr = Libderiv->deriv2_classes[4][4][143];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,315,dvrr_stack+67847, dvrr_stack+27803, NULL);
 tmp = dvrr_stack + 67847;
 target_ptr = Libderiv->deriv2_classes[4][5][143];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,420,dvrr_stack+65470, dvrr_stack+34418, NULL);
 tmp = dvrr_stack + 65470;
 target_ptr = Libderiv->deriv2_classes[4][6][143];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,210,dvrr_stack+68162, dvrr_stack+40755, NULL);
 tmp = dvrr_stack + 68162;
 target_ptr = Libderiv->deriv2_classes[5][3][143];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,315,dvrr_stack+92527, dvrr_stack+46054, NULL);
 tmp = dvrr_stack + 92527;
 target_ptr = Libderiv->deriv2_classes[5][4][143];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,441,dvrr_stack+157526, dvrr_stack+53281, NULL);
 tmp = dvrr_stack + 157526;
 target_ptr = Libderiv->deriv2_classes[5][5][143];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,588,dvrr_stack+153273, dvrr_stack+62740, NULL);
 tmp = dvrr_stack + 153273;
 target_ptr = Libderiv->deriv2_classes[5][6][143];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,280,dvrr_stack+153861, dvrr_stack+71438, NULL);
 tmp = dvrr_stack + 153861;
 target_ptr = Libderiv->deriv2_classes[6][3][143];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,420,dvrr_stack+154141, dvrr_stack+78634, NULL);
 tmp = dvrr_stack + 154141;
 target_ptr = Libderiv->deriv2_classes[6][4][143];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,588,dvrr_stack+73426, dvrr_stack+88438, NULL);
 tmp = dvrr_stack + 73426;
 target_ptr = Libderiv->deriv2_classes[6][5][143];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,784,dvrr_stack+157967, dvrr_stack+101260, NULL);
 tmp = dvrr_stack + 157967;
 target_ptr = Libderiv->deriv2_classes[6][6][143];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,100,dvrr_stack+68372, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 68372;
 target_ptr = Libderiv->deriv2_classes[3][3][131];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+65890, dvrr_stack+5969, NULL);
 tmp = dvrr_stack + 65890;
 target_ptr = Libderiv->deriv2_classes[3][4][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,210,dvrr_stack+3906, dvrr_stack+9828, NULL);
 tmp = dvrr_stack + 3906;
 target_ptr = Libderiv->deriv2_classes[3][5][131];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,280,dvrr_stack+74014, dvrr_stack+14703, NULL);
 tmp = dvrr_stack + 74014;
 target_ptr = Libderiv->deriv2_classes[3][6][131];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+158751, dvrr_stack+19057, NULL);
 tmp = dvrr_stack + 158751;
 target_ptr = Libderiv->deriv2_classes[4][3][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,225,dvrr_stack+16078, dvrr_stack+22754, NULL);
 tmp = dvrr_stack + 16078;
 target_ptr = Libderiv->deriv2_classes[4][4][131];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,315,dvrr_stack+24297, dvrr_stack+27803, NULL);
 tmp = dvrr_stack + 24297;
 target_ptr = Libderiv->deriv2_classes[4][5][131];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,420,dvrr_stack+42141, dvrr_stack+34418, NULL);
 tmp = dvrr_stack + 42141;
 target_ptr = Libderiv->deriv2_classes[4][6][131];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,210,dvrr_stack+20242, dvrr_stack+40755, NULL);
 tmp = dvrr_stack + 20242;
 target_ptr = Libderiv->deriv2_classes[5][3][131];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,315,dvrr_stack+0, dvrr_stack+46054, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[5][4][131];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,441,dvrr_stack+37995, dvrr_stack+53281, NULL);
 tmp = dvrr_stack + 37995;
 target_ptr = Libderiv->deriv2_classes[5][5][131];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,588,dvrr_stack+158901, dvrr_stack+62740, NULL);
 tmp = dvrr_stack + 158901;
 target_ptr = Libderiv->deriv2_classes[5][6][131];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,280,dvrr_stack+42561, dvrr_stack+71438, NULL);
 tmp = dvrr_stack + 42561;
 target_ptr = Libderiv->deriv2_classes[6][3][131];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,420,dvrr_stack+159489, dvrr_stack+78634, NULL);
 tmp = dvrr_stack + 159489;
 target_ptr = Libderiv->deriv2_classes[6][4][131];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,588,dvrr_stack+159909, dvrr_stack+88438, NULL);
 tmp = dvrr_stack + 159909;
 target_ptr = Libderiv->deriv2_classes[6][5][131];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,784,dvrr_stack+160497, dvrr_stack+101260, NULL);
 tmp = dvrr_stack + 160497;
 target_ptr = Libderiv->deriv2_classes[6][6][131];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+92842, dvrr_stack+54604, NULL);
 tmp = dvrr_stack + 92842;
 target_ptr = Libderiv->deriv2_classes[3][3][130];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+161281, dvrr_stack+54904, NULL);
 tmp = dvrr_stack + 161281;
 target_ptr = Libderiv->deriv2_classes[3][4][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+161431, dvrr_stack+55354, NULL);
 tmp = dvrr_stack + 161431;
 target_ptr = Libderiv->deriv2_classes[3][5][130];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+161641, dvrr_stack+90202, NULL);
 tmp = dvrr_stack + 161641;
 target_ptr = Libderiv->deriv2_classes[3][6][130];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+22517, dvrr_stack+91042, NULL);
 tmp = dvrr_stack + 22517;
 target_ptr = Libderiv->deriv2_classes[4][3][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+15543, dvrr_stack+91492, NULL);
 tmp = dvrr_stack + 15543;
 target_ptr = Libderiv->deriv2_classes[4][4][130];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+36498, dvrr_stack+103612, NULL);
 tmp = dvrr_stack + 36498;
 target_ptr = Libderiv->deriv2_classes[4][5][130];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+161921, dvrr_stack+104557, NULL);
 tmp = dvrr_stack + 161921;
 target_ptr = Libderiv->deriv2_classes[4][6][130];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+162341, dvrr_stack+28748, NULL);
 tmp = dvrr_stack + 162341;
 target_ptr = Libderiv->deriv2_classes[5][3][130];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+162551, dvrr_stack+105817, NULL);
 tmp = dvrr_stack + 162551;
 target_ptr = Libderiv->deriv2_classes[5][4][130];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+162866, dvrr_stack+106762, NULL);
 tmp = dvrr_stack + 162866;
 target_ptr = Libderiv->deriv2_classes[5][5][130];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+163307, dvrr_stack+108085, NULL);
 tmp = dvrr_stack + 163307;
 target_ptr = Libderiv->deriv2_classes[5][6][130];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+6419, dvrr_stack+109849, NULL);
 tmp = dvrr_stack + 6419;
 target_ptr = Libderiv->deriv2_classes[6][3][130];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+163895, dvrr_stack+110689, NULL);
 tmp = dvrr_stack + 163895;
 target_ptr = Libderiv->deriv2_classes[6][4][130];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+164315, dvrr_stack+111949, NULL);
 tmp = dvrr_stack + 164315;
 target_ptr = Libderiv->deriv2_classes[6][5][130];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,784,dvrr_stack+164903, dvrr_stack+113713, NULL);
 tmp = dvrr_stack + 164903;
 target_ptr = Libderiv->deriv2_classes[6][6][130];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,100,dvrr_stack+1737, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 1737;
 target_ptr = Libderiv->deriv2_classes[3][3][119];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+3123, dvrr_stack+5969, NULL);
 tmp = dvrr_stack + 3123;
 target_ptr = Libderiv->deriv2_classes[3][4][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,210,dvrr_stack+165687, dvrr_stack+9828, NULL);
 tmp = dvrr_stack + 165687;
 target_ptr = Libderiv->deriv2_classes[3][5][119];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,280,dvrr_stack+9828, dvrr_stack+14703, NULL);
 tmp = dvrr_stack + 9828;
 target_ptr = Libderiv->deriv2_classes[3][6][119];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+3273, dvrr_stack+19057, NULL);
 tmp = dvrr_stack + 3273;
 target_ptr = Libderiv->deriv2_classes[4][3][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,225,dvrr_stack+19057, dvrr_stack+22754, NULL);
 tmp = dvrr_stack + 19057;
 target_ptr = Libderiv->deriv2_classes[4][4][119];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,315,dvrr_stack+10108, dvrr_stack+27803, NULL);
 tmp = dvrr_stack + 10108;
 target_ptr = Libderiv->deriv2_classes[4][5][119];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,420,dvrr_stack+27803, dvrr_stack+34418, NULL);
 tmp = dvrr_stack + 27803;
 target_ptr = Libderiv->deriv2_classes[4][6][119];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,210,dvrr_stack+28223, dvrr_stack+40755, NULL);
 tmp = dvrr_stack + 28223;
 target_ptr = Libderiv->deriv2_classes[5][3][119];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,315,dvrr_stack+28433, dvrr_stack+46054, NULL);
 tmp = dvrr_stack + 28433;
 target_ptr = Libderiv->deriv2_classes[5][4][119];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,441,dvrr_stack+40755, dvrr_stack+53281, NULL);
 tmp = dvrr_stack + 40755;
 target_ptr = Libderiv->deriv2_classes[5][5][119];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,588,dvrr_stack+22667, dvrr_stack+62740, NULL);
 tmp = dvrr_stack + 22667;
 target_ptr = Libderiv->deriv2_classes[5][6][119];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,280,dvrr_stack+5924, dvrr_stack+71438, NULL);
 tmp = dvrr_stack + 5924;
 target_ptr = Libderiv->deriv2_classes[6][3][119];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,420,dvrr_stack+71438, dvrr_stack+78634, NULL);
 tmp = dvrr_stack + 71438;
 target_ptr = Libderiv->deriv2_classes[6][4][119];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,588,dvrr_stack+78634, dvrr_stack+88438, NULL);
 tmp = dvrr_stack + 78634;
 target_ptr = Libderiv->deriv2_classes[6][5][119];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,784,dvrr_stack+14481, dvrr_stack+101260, NULL);
 tmp = dvrr_stack + 14481;
 target_ptr = Libderiv->deriv2_classes[6][6][119];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+79222, dvrr_stack+54604, NULL);
 tmp = dvrr_stack + 79222;
 target_ptr = Libderiv->deriv2_classes[3][3][118];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+79322, dvrr_stack+54904, NULL);
 tmp = dvrr_stack + 79322;
 target_ptr = Libderiv->deriv2_classes[3][4][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+79472, dvrr_stack+55354, NULL);
 tmp = dvrr_stack + 79472;
 target_ptr = Libderiv->deriv2_classes[3][5][118];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+71858, dvrr_stack+90202, NULL);
 tmp = dvrr_stack + 71858;
 target_ptr = Libderiv->deriv2_classes[3][6][118];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+79682, dvrr_stack+91042, NULL);
 tmp = dvrr_stack + 79682;
 target_ptr = Libderiv->deriv2_classes[4][3][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+19282, dvrr_stack+91492, NULL);
 tmp = dvrr_stack + 19282;
 target_ptr = Libderiv->deriv2_classes[4][4][118];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+101223, dvrr_stack+103612, NULL);
 tmp = dvrr_stack + 101223;
 target_ptr = Libderiv->deriv2_classes[4][5][118];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+101538, dvrr_stack+104557, NULL);
 tmp = dvrr_stack + 101538;
 target_ptr = Libderiv->deriv2_classes[4][6][118];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+6204, dvrr_stack+28748, NULL);
 tmp = dvrr_stack + 6204;
 target_ptr = Libderiv->deriv2_classes[5][3][118];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+28748, dvrr_stack+105817, NULL);
 tmp = dvrr_stack + 28748;
 target_ptr = Libderiv->deriv2_classes[5][4][118];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+101958, dvrr_stack+106762, NULL);
 tmp = dvrr_stack + 101958;
 target_ptr = Libderiv->deriv2_classes[5][5][118];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+102399, dvrr_stack+108085, NULL);
 tmp = dvrr_stack + 102399;
 target_ptr = Libderiv->deriv2_classes[5][6][118];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+29063, dvrr_stack+109849, NULL);
 tmp = dvrr_stack + 29063;
 target_ptr = Libderiv->deriv2_classes[6][3][118];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+102987, dvrr_stack+110689, NULL);
 tmp = dvrr_stack + 102987;
 target_ptr = Libderiv->deriv2_classes[6][4][118];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+103407, dvrr_stack+111949, NULL);
 tmp = dvrr_stack + 103407;
 target_ptr = Libderiv->deriv2_classes[6][5][118];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,784,dvrr_stack+103995, dvrr_stack+113713, NULL);
 tmp = dvrr_stack + 103995;
 target_ptr = Libderiv->deriv2_classes[6][6][118];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+72138, dvrr_stack+23457, NULL);
 tmp = dvrr_stack + 72138;
 target_ptr = Libderiv->deriv2_classes[3][3][117];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+41196, dvrr_stack+2523, NULL);
 tmp = dvrr_stack + 41196;
 target_ptr = Libderiv->deriv2_classes[3][4][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+23255, dvrr_stack+5069, NULL);
 tmp = dvrr_stack + 23255;
 target_ptr = Libderiv->deriv2_classes[3][5][117];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+23465, dvrr_stack+8568, NULL);
 tmp = dvrr_stack + 23465;
 target_ptr = Libderiv->deriv2_classes[3][6][117];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+15265, dvrr_stack+13023, NULL);
 tmp = dvrr_stack + 15265;
 target_ptr = Libderiv->deriv2_classes[4][3][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+8403, dvrr_stack+18157, NULL);
 tmp = dvrr_stack + 8403;
 target_ptr = Libderiv->deriv2_classes[4][4][117];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+8628, dvrr_stack+21404, NULL);
 tmp = dvrr_stack + 8628;
 target_ptr = Libderiv->deriv2_classes[4][5][117];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+8943, dvrr_stack+25913, NULL);
 tmp = dvrr_stack + 8943;
 target_ptr = Libderiv->deriv2_classes[4][6][117];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+21299, dvrr_stack+27173, NULL);
 tmp = dvrr_stack + 21299;
 target_ptr = Libderiv->deriv2_classes[5][3][117];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+21509, dvrr_stack+39495, NULL);
 tmp = dvrr_stack + 21509;
 target_ptr = Libderiv->deriv2_classes[5][4][117];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+39495, dvrr_stack+44164, NULL);
 tmp = dvrr_stack + 39495;
 target_ptr = Libderiv->deriv2_classes[5][5][117];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+44164, dvrr_stack+50635, NULL);
 tmp = dvrr_stack + 44164;
 target_ptr = Libderiv->deriv2_classes[5][6][117];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+50635, dvrr_stack+59212, NULL);
 tmp = dvrr_stack + 50635;
 target_ptr = Libderiv->deriv2_classes[6][3][117];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+50915, dvrr_stack+69758, NULL);
 tmp = dvrr_stack + 50915;
 target_ptr = Libderiv->deriv2_classes[6][4][117];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+51335, dvrr_stack+76114, NULL);
 tmp = dvrr_stack + 51335;
 target_ptr = Libderiv->deriv2_classes[6][5][117];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,784,dvrr_stack+104779, dvrr_stack+84910, NULL);
 tmp = dvrr_stack + 104779;
 target_ptr = Libderiv->deriv2_classes[6][6][117];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+51923, dvrr_stack+2973, dvrr_stack+96868);
 tmp = dvrr_stack + 51923;
 target_ptr = Libderiv->deriv2_classes[3][3][107];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+52023, dvrr_stack+97028, dvrr_stack+96928);
 tmp = dvrr_stack + 52023;
 target_ptr = Libderiv->deriv2_classes[3][4][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+52173, dvrr_stack+97238, dvrr_stack+2973);
 tmp = dvrr_stack + 52173;
 target_ptr = Libderiv->deriv2_classes[3][5][107];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+44752, dvrr_stack+97518, dvrr_stack+97028);
 tmp = dvrr_stack + 44752;
 target_ptr = Libderiv->deriv2_classes[3][6][107];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+45032, dvrr_stack+18832, dvrr_stack+97998);
 tmp = dvrr_stack + 45032;
 target_ptr = Libderiv->deriv2_classes[4][3][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+45182, dvrr_stack+40440, dvrr_stack+98088);
 tmp = dvrr_stack + 45182;
 target_ptr = Libderiv->deriv2_classes[4][4][107];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+39936, dvrr_stack+71018, dvrr_stack+18832);
 tmp = dvrr_stack + 39936;
 target_ptr = Libderiv->deriv2_classes[4][5][107];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+21824, dvrr_stack+98238, dvrr_stack+40440);
 tmp = dvrr_stack + 21824;
 target_ptr = Libderiv->deriv2_classes[4][6][107];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+105563, dvrr_stack+99453, dvrr_stack+99327);
 tmp = dvrr_stack + 105563;
 target_ptr = Libderiv->deriv2_classes[5][3][107];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+105773, dvrr_stack+99978, dvrr_stack+99768);
 tmp = dvrr_stack + 105773;
 target_ptr = Libderiv->deriv2_classes[5][4][107];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+106088, dvrr_stack+100419, dvrr_stack+99453);
 tmp = dvrr_stack + 106088;
 target_ptr = Libderiv->deriv2_classes[5][5][107];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+106529, dvrr_stack+77878, dvrr_stack+99978);
 tmp = dvrr_stack + 106529;
 target_ptr = Libderiv->deriv2_classes[5][6][107];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,28,1,dvrr_stack+107117, dvrr_stack+9408, dvrr_stack+87766);
 tmp = dvrr_stack + 107117;
 target_ptr = Libderiv->deriv2_classes[6][3][107];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+107397, dvrr_stack+60052, dvrr_stack+87934);
 tmp = dvrr_stack + 107397;
 target_ptr = Libderiv->deriv2_classes[6][4][107];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+107817, dvrr_stack+60640, dvrr_stack+9408);
 tmp = dvrr_stack + 107817;
 target_ptr = Libderiv->deriv2_classes[6][5][107];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+108405, dvrr_stack+61424, dvrr_stack+60052);
 tmp = dvrr_stack + 108405;
 target_ptr = Libderiv->deriv2_classes[6][6][107];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+40251, dvrr_stack+88274, dvrr_stack+88214);
 tmp = dvrr_stack + 40251;
 target_ptr = Libderiv->deriv2_classes[3][3][106];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+109189, dvrr_stack+52399, dvrr_stack+62432);
 tmp = dvrr_stack + 109189;
 target_ptr = Libderiv->deriv2_classes[3][4][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+109339, dvrr_stack+52609, dvrr_stack+88274);
 tmp = dvrr_stack + 109339;
 target_ptr = Libderiv->deriv2_classes[3][5][106];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+109549, dvrr_stack+52889, dvrr_stack+52399);
 tmp = dvrr_stack + 109549;
 target_ptr = Libderiv->deriv2_classes[3][6][106];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+109829, dvrr_stack+29378, dvrr_stack+62532);
 tmp = dvrr_stack + 109829;
 target_ptr = Libderiv->deriv2_classes[4][3][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+109979, dvrr_stack+45637, dvrr_stack+45487);
 tmp = dvrr_stack + 109979;
 target_ptr = Libderiv->deriv2_classes[4][4][106];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+110204, dvrr_stack+31898, dvrr_stack+29378);
 tmp = dvrr_stack + 110204;
 target_ptr = Libderiv->deriv2_classes[4][5][106];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+110519, dvrr_stack+32318, dvrr_stack+45637);
 tmp = dvrr_stack + 110519;
 target_ptr = Libderiv->deriv2_classes[4][6][106];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+110939, dvrr_stack+32984, dvrr_stack+32858);
 tmp = dvrr_stack + 110939;
 target_ptr = Libderiv->deriv2_classes[5][3][106];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+111149, dvrr_stack+33509, dvrr_stack+33299);
 tmp = dvrr_stack + 111149;
 target_ptr = Libderiv->deriv2_classes[5][4][106];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+111464, dvrr_stack+13473, dvrr_stack+32984);
 tmp = dvrr_stack + 111464;
 target_ptr = Libderiv->deriv2_classes[5][5][106];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+111905, dvrr_stack+116065, dvrr_stack+33509);
 tmp = dvrr_stack + 111905;
 target_ptr = Libderiv->deriv2_classes[5][6][106];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,28,1,dvrr_stack+112493, dvrr_stack+14061, dvrr_stack+33950);
 tmp = dvrr_stack + 112493;
 target_ptr = Libderiv->deriv2_classes[6][3][106];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+112773, dvrr_stack+116821, dvrr_stack+34118);
 tmp = dvrr_stack + 112773;
 target_ptr = Libderiv->deriv2_classes[6][4][106];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+113193, dvrr_stack+117409, dvrr_stack+14061);
 tmp = dvrr_stack + 113193;
 target_ptr = Libderiv->deriv2_classes[6][5][106];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+113781, dvrr_stack+118193, dvrr_stack+116821);
 tmp = dvrr_stack + 113781;
 target_ptr = Libderiv->deriv2_classes[6][6][106];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+4309, dvrr_stack+96688, dvrr_stack+62622);
 tmp = dvrr_stack + 4309;
 target_ptr = Libderiv->deriv2_classes[3][3][105];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+114565, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 114565;
 target_ptr = Libderiv->deriv2_classes[3][4][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+114715, dvrr_stack+4439, dvrr_stack+96688);
 tmp = dvrr_stack + 114715;
 target_ptr = Libderiv->deriv2_classes[3][5][105];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+114925, dvrr_stack+7728, dvrr_stack+1023);
 tmp = dvrr_stack + 114925;
 target_ptr = Libderiv->deriv2_classes[3][6][105];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+115205, dvrr_stack+5699, dvrr_stack+11943);
 tmp = dvrr_stack + 115205;
 target_ptr = Libderiv->deriv2_classes[4][3][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+115355, dvrr_stack+16318, dvrr_stack+17482);
 tmp = dvrr_stack + 115355;
 target_ptr = Libderiv->deriv2_classes[4][4][105];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+115580, dvrr_stack+20459, dvrr_stack+5699);
 tmp = dvrr_stack + 115580;
 target_ptr = Libderiv->deriv2_classes[4][5][105];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+2299, dvrr_stack+24653, dvrr_stack+16318);
 tmp = dvrr_stack + 2299;
 target_ptr = Libderiv->deriv2_classes[4][6][105];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+2719, dvrr_stack+98949, dvrr_stack+30278);
 tmp = dvrr_stack + 2719;
 target_ptr = Libderiv->deriv2_classes[5][3][105];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+84868, dvrr_stack+36828, dvrr_stack+38550);
 tmp = dvrr_stack + 84868;
 target_ptr = Libderiv->deriv2_classes[5][4][105];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+85183, dvrr_stack+42841, dvrr_stack+98949);
 tmp = dvrr_stack + 85183;
 target_ptr = Libderiv->deriv2_classes[5][5][105];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+85624, dvrr_stack+48871, dvrr_stack+36828);
 tmp = dvrr_stack + 85624;
 target_ptr = Libderiv->deriv2_classes[5][6][105];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,28,1,dvrr_stack+86212, dvrr_stack+87262, dvrr_stack+56944);
 tmp = dvrr_stack + 86212;
 target_ptr = Libderiv->deriv2_classes[6][3][105];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+86492, dvrr_stack+66104, dvrr_stack+68498);
 tmp = dvrr_stack + 86492;
 target_ptr = Libderiv->deriv2_classes[6][4][105];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+76062, dvrr_stack+74350, dvrr_stack+87262);
 tmp = dvrr_stack + 76062;
 target_ptr = Libderiv->deriv2_classes[6][5][105];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+76650, dvrr_stack+49627, dvrr_stack+66104);
 tmp = dvrr_stack + 76650;
 target_ptr = Libderiv->deriv2_classes[6][6][105];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+165897, dvrr_stack+93592, dvrr_stack+93532);
 tmp = dvrr_stack + 165897;
 target_ptr = Libderiv->deriv2_classes[3][3][104];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+115895, dvrr_stack+93842, dvrr_stack+93742);
 tmp = dvrr_stack + 115895;
 target_ptr = Libderiv->deriv2_classes[3][4][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+86912, dvrr_stack+94052, dvrr_stack+93592);
 tmp = dvrr_stack + 86912;
 target_ptr = Libderiv->deriv2_classes[3][5][104];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+77434, dvrr_stack+94332, dvrr_stack+93842);
 tmp = dvrr_stack + 77434;
 target_ptr = Libderiv->deriv2_classes[3][6][104];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+77714, dvrr_stack+94782, dvrr_stack+94692);
 tmp = dvrr_stack + 77714;
 target_ptr = Libderiv->deriv2_classes[4][3][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+69726, dvrr_stack+95157, dvrr_stack+95007);
 tmp = dvrr_stack + 69726;
 target_ptr = Libderiv->deriv2_classes[4][4][104];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+69951, dvrr_stack+95472, dvrr_stack+94782);
 tmp = dvrr_stack + 69951;
 target_ptr = Libderiv->deriv2_classes[4][5][104];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+70266, dvrr_stack+95892, dvrr_stack+95157);
 tmp = dvrr_stack + 70266;
 target_ptr = Libderiv->deriv2_classes[4][6][104];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+70686, dvrr_stack+82558, dvrr_stack+96432);
 tmp = dvrr_stack + 70686;
 target_ptr = Libderiv->deriv2_classes[5][3][104];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+59184, dvrr_stack+83083, dvrr_stack+82873);
 tmp = dvrr_stack + 59184;
 target_ptr = Libderiv->deriv2_classes[5][4][104];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+59499, dvrr_stack+83524, dvrr_stack+82558);
 tmp = dvrr_stack + 59499;
 target_ptr = Libderiv->deriv2_classes[5][5][104];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+25844, dvrr_stack+84112, dvrr_stack+83083);
 tmp = dvrr_stack + 25844;
 target_ptr = Libderiv->deriv2_classes[5][6][104];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,28,1,dvrr_stack+26432, dvrr_stack+75302, dvrr_stack+75134);
 tmp = dvrr_stack + 26432;
 target_ptr = Libderiv->deriv2_classes[6][3][104];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+26712, dvrr_stack+68778, dvrr_stack+75722);
 tmp = dvrr_stack + 26712;
 target_ptr = Libderiv->deriv2_classes[6][4][104];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+27132, dvrr_stack+57112, dvrr_stack+75302);
 tmp = dvrr_stack + 27132;
 target_ptr = Libderiv->deriv2_classes[6][5][104];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+88424, dvrr_stack+57896, dvrr_stack+68778);
 tmp = dvrr_stack + 88424;
 target_ptr = Libderiv->deriv2_classes[6][6][104];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+15415, dvrr_stack+2973, dvrr_stack+96868);
 tmp = dvrr_stack + 15415;
 target_ptr = Libderiv->deriv2_classes[3][3][95];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+18062, dvrr_stack+97028, dvrr_stack+96928);
 tmp = dvrr_stack + 18062;
 target_ptr = Libderiv->deriv2_classes[3][4][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+18212, dvrr_stack+97238, dvrr_stack+2973);
 tmp = dvrr_stack + 18212;
 target_ptr = Libderiv->deriv2_classes[3][5][95];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+18422, dvrr_stack+97518, dvrr_stack+97028);
 tmp = dvrr_stack + 18422;
 target_ptr = Libderiv->deriv2_classes[3][6][95];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+12901, dvrr_stack+18832, dvrr_stack+97998);
 tmp = dvrr_stack + 12901;
 target_ptr = Libderiv->deriv2_classes[4][3][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+13051, dvrr_stack+40440, dvrr_stack+98088);
 tmp = dvrr_stack + 13051;
 target_ptr = Libderiv->deriv2_classes[4][4][95];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+5067, dvrr_stack+71018, dvrr_stack+18832);
 tmp = dvrr_stack + 5067;
 target_ptr = Libderiv->deriv2_classes[4][5][95];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+89208, dvrr_stack+98238, dvrr_stack+40440);
 tmp = dvrr_stack + 89208;
 target_ptr = Libderiv->deriv2_classes[4][6][95];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+5382, dvrr_stack+99453, dvrr_stack+99327);
 tmp = dvrr_stack + 5382;
 target_ptr = Libderiv->deriv2_classes[5][3][95];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+89628, dvrr_stack+99978, dvrr_stack+99768);
 tmp = dvrr_stack + 89628;
 target_ptr = Libderiv->deriv2_classes[5][4][95];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+89943, dvrr_stack+100419, dvrr_stack+99453);
 tmp = dvrr_stack + 89943;
 target_ptr = Libderiv->deriv2_classes[5][5][95];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+90384, dvrr_stack+77878, dvrr_stack+99978);
 tmp = dvrr_stack + 90384;
 target_ptr = Libderiv->deriv2_classes[5][6][95];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,28,1,dvrr_stack+90972, dvrr_stack+9408, dvrr_stack+87766);
 tmp = dvrr_stack + 90972;
 target_ptr = Libderiv->deriv2_classes[6][3][95];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+91252, dvrr_stack+60052, dvrr_stack+87934);
 tmp = dvrr_stack + 91252;
 target_ptr = Libderiv->deriv2_classes[6][4][95];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+62682, dvrr_stack+60640, dvrr_stack+9408);
 tmp = dvrr_stack + 62682;
 target_ptr = Libderiv->deriv2_classes[6][5][95];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+63270, dvrr_stack+61424, dvrr_stack+60052);
 tmp = dvrr_stack + 63270;
 target_ptr = Libderiv->deriv2_classes[6][6][95];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+22244, dvrr_stack+88274, dvrr_stack+88214);
 tmp = dvrr_stack + 22244;
 target_ptr = Libderiv->deriv2_classes[3][3][94];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+13276, dvrr_stack+52399, dvrr_stack+62432);
 tmp = dvrr_stack + 13276;
 target_ptr = Libderiv->deriv2_classes[3][4][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+1233, dvrr_stack+52609, dvrr_stack+88274);
 tmp = dvrr_stack + 1233;
 target_ptr = Libderiv->deriv2_classes[3][5][94];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+91672, dvrr_stack+52889, dvrr_stack+52399);
 tmp = dvrr_stack + 91672;
 target_ptr = Libderiv->deriv2_classes[3][6][94];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+91952, dvrr_stack+29378, dvrr_stack+62532);
 tmp = dvrr_stack + 91952;
 target_ptr = Libderiv->deriv2_classes[4][3][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+64054, dvrr_stack+45637, dvrr_stack+45487);
 tmp = dvrr_stack + 64054;
 target_ptr = Libderiv->deriv2_classes[4][4][94];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+53249, dvrr_stack+31898, dvrr_stack+29378);
 tmp = dvrr_stack + 53249;
 target_ptr = Libderiv->deriv2_classes[4][5][94];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+53564, dvrr_stack+32318, dvrr_stack+45637);
 tmp = dvrr_stack + 53564;
 target_ptr = Libderiv->deriv2_classes[4][6][94];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+64279, dvrr_stack+32984, dvrr_stack+32858);
 tmp = dvrr_stack + 64279;
 target_ptr = Libderiv->deriv2_classes[5][3][94];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+53984, dvrr_stack+33509, dvrr_stack+33299);
 tmp = dvrr_stack + 53984;
 target_ptr = Libderiv->deriv2_classes[5][4][94];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+54299, dvrr_stack+13473, dvrr_stack+32984);
 tmp = dvrr_stack + 54299;
 target_ptr = Libderiv->deriv2_classes[5][5][94];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+54740, dvrr_stack+116065, dvrr_stack+33509);
 tmp = dvrr_stack + 54740;
 target_ptr = Libderiv->deriv2_classes[5][6][94];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,28,1,dvrr_stack+55328, dvrr_stack+14061, dvrr_stack+33950);
 tmp = dvrr_stack + 55328;
 target_ptr = Libderiv->deriv2_classes[6][3][94];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+45952, dvrr_stack+116821, dvrr_stack+34118);
 tmp = dvrr_stack + 45952;
 target_ptr = Libderiv->deriv2_classes[6][4][94];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+46372, dvrr_stack+117409, dvrr_stack+14061);
 tmp = dvrr_stack + 46372;
 target_ptr = Libderiv->deriv2_classes[6][5][94];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+34398, dvrr_stack+118193, dvrr_stack+116821);
 tmp = dvrr_stack + 34398;
 target_ptr = Libderiv->deriv2_classes[6][6][94];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+97878, dvrr_stack+96688, dvrr_stack+62622);
 tmp = dvrr_stack + 97878;
 target_ptr = Libderiv->deriv2_classes[3][3][93];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+55608, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 55608;
 target_ptr = Libderiv->deriv2_classes[3][4][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+55758, dvrr_stack+4439, dvrr_stack+96688);
 tmp = dvrr_stack + 55758;
 target_ptr = Libderiv->deriv2_classes[3][5][93];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+35182, dvrr_stack+7728, dvrr_stack+1023);
 tmp = dvrr_stack + 35182;
 target_ptr = Libderiv->deriv2_classes[3][6][93];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+35462, dvrr_stack+5699, dvrr_stack+11943);
 tmp = dvrr_stack + 35462;
 target_ptr = Libderiv->deriv2_classes[4][3][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+165997, dvrr_stack+16318, dvrr_stack+17482);
 tmp = dvrr_stack + 165997;
 target_ptr = Libderiv->deriv2_classes[4][4][93];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+166222, dvrr_stack+20459, dvrr_stack+5699);
 tmp = dvrr_stack + 166222;
 target_ptr = Libderiv->deriv2_classes[4][5][93];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+166537, dvrr_stack+24653, dvrr_stack+16318);
 tmp = dvrr_stack + 166537;
 target_ptr = Libderiv->deriv2_classes[4][6][93];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+166957, dvrr_stack+98949, dvrr_stack+30278);
 tmp = dvrr_stack + 166957;
 target_ptr = Libderiv->deriv2_classes[5][3][93];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+167167, dvrr_stack+36828, dvrr_stack+38550);
 tmp = dvrr_stack + 167167;
 target_ptr = Libderiv->deriv2_classes[5][4][93];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+167482, dvrr_stack+42841, dvrr_stack+98949);
 tmp = dvrr_stack + 167482;
 target_ptr = Libderiv->deriv2_classes[5][5][93];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+167923, dvrr_stack+48871, dvrr_stack+36828);
 tmp = dvrr_stack + 167923;
 target_ptr = Libderiv->deriv2_classes[5][6][93];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,28,1,dvrr_stack+168511, dvrr_stack+87262, dvrr_stack+56944);
 tmp = dvrr_stack + 168511;
 target_ptr = Libderiv->deriv2_classes[6][3][93];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+168791, dvrr_stack+66104, dvrr_stack+68498);
 tmp = dvrr_stack + 168791;
 target_ptr = Libderiv->deriv2_classes[6][4][93];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+169211, dvrr_stack+74350, dvrr_stack+87262);
 tmp = dvrr_stack + 169211;
 target_ptr = Libderiv->deriv2_classes[6][5][93];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+169799, dvrr_stack+49627, dvrr_stack+66104);
 tmp = dvrr_stack + 169799;
 target_ptr = Libderiv->deriv2_classes[6][6][93];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+170583, dvrr_stack+93592, dvrr_stack+93532);
 tmp = dvrr_stack + 170583;
 target_ptr = Libderiv->deriv2_classes[3][3][92];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+170683, dvrr_stack+93842, dvrr_stack+93742);
 tmp = dvrr_stack + 170683;
 target_ptr = Libderiv->deriv2_classes[3][4][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+170833, dvrr_stack+94052, dvrr_stack+93592);
 tmp = dvrr_stack + 170833;
 target_ptr = Libderiv->deriv2_classes[3][5][92];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+171043, dvrr_stack+94332, dvrr_stack+93842);
 tmp = dvrr_stack + 171043;
 target_ptr = Libderiv->deriv2_classes[3][6][92];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+171323, dvrr_stack+94782, dvrr_stack+94692);
 tmp = dvrr_stack + 171323;
 target_ptr = Libderiv->deriv2_classes[4][3][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+171473, dvrr_stack+95157, dvrr_stack+95007);
 tmp = dvrr_stack + 171473;
 target_ptr = Libderiv->deriv2_classes[4][4][92];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+171698, dvrr_stack+95472, dvrr_stack+94782);
 tmp = dvrr_stack + 171698;
 target_ptr = Libderiv->deriv2_classes[4][5][92];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+172013, dvrr_stack+95892, dvrr_stack+95157);
 tmp = dvrr_stack + 172013;
 target_ptr = Libderiv->deriv2_classes[4][6][92];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+172433, dvrr_stack+82558, dvrr_stack+96432);
 tmp = dvrr_stack + 172433;
 target_ptr = Libderiv->deriv2_classes[5][3][92];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+172643, dvrr_stack+83083, dvrr_stack+82873);
 tmp = dvrr_stack + 172643;
 target_ptr = Libderiv->deriv2_classes[5][4][92];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+172958, dvrr_stack+83524, dvrr_stack+82558);
 tmp = dvrr_stack + 172958;
 target_ptr = Libderiv->deriv2_classes[5][5][92];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+173399, dvrr_stack+84112, dvrr_stack+83083);
 tmp = dvrr_stack + 173399;
 target_ptr = Libderiv->deriv2_classes[5][6][92];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,28,1,dvrr_stack+173987, dvrr_stack+75302, dvrr_stack+75134);
 tmp = dvrr_stack + 173987;
 target_ptr = Libderiv->deriv2_classes[6][3][92];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+174267, dvrr_stack+68778, dvrr_stack+75722);
 tmp = dvrr_stack + 174267;
 target_ptr = Libderiv->deriv2_classes[6][4][92];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+174687, dvrr_stack+57112, dvrr_stack+75302);
 tmp = dvrr_stack + 174687;
 target_ptr = Libderiv->deriv2_classes[6][5][92];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+175275, dvrr_stack+57896, dvrr_stack+68778);
 tmp = dvrr_stack + 175275;
 target_ptr = Libderiv->deriv2_classes[6][6][92];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+176059, dvrr_stack+66692, dvrr_stack+76002);
 tmp = dvrr_stack + 176059;
 target_ptr = Libderiv->deriv2_classes[3][3][91];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+176159, dvrr_stack+69366, dvrr_stack+66842);
 tmp = dvrr_stack + 176159;
 target_ptr = Libderiv->deriv2_classes[3][4][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+176309, dvrr_stack+58904, dvrr_stack+66692);
 tmp = dvrr_stack + 176309;
 target_ptr = Libderiv->deriv2_classes[3][5][91];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+176519, dvrr_stack+43429, dvrr_stack+69366);
 tmp = dvrr_stack + 176519;
 target_ptr = Libderiv->deriv2_classes[3][6][91];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+176799, dvrr_stack+43789, dvrr_stack+69576);
 tmp = dvrr_stack + 176799;
 target_ptr = Libderiv->deriv2_classes[4][3][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+176949, dvrr_stack+38760, dvrr_stack+44014);
 tmp = dvrr_stack + 176949;
 target_ptr = Libderiv->deriv2_classes[4][4][91];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+177174, dvrr_stack+39075, dvrr_stack+43789);
 tmp = dvrr_stack + 177174;
 target_ptr = Libderiv->deriv2_classes[4][5][91];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+177489, dvrr_stack+30404, dvrr_stack+38760);
 tmp = dvrr_stack + 177489;
 target_ptr = Libderiv->deriv2_classes[4][6][91];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+177909, dvrr_stack+30944, dvrr_stack+37269);
 tmp = dvrr_stack + 177909;
 target_ptr = Libderiv->deriv2_classes[5][3][91];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+178119, dvrr_stack+25193, dvrr_stack+31259);
 tmp = dvrr_stack + 178119;
 target_ptr = Libderiv->deriv2_classes[5][4][91];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+178434, dvrr_stack+12033, dvrr_stack+30944);
 tmp = dvrr_stack + 178434;
 target_ptr = Libderiv->deriv2_classes[5][5][91];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+178875, dvrr_stack+119201, dvrr_stack+25193);
 tmp = dvrr_stack + 178875;
 target_ptr = Libderiv->deriv2_classes[5][6][91];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,28,1,dvrr_stack+179463, dvrr_stack+20879, dvrr_stack+31469);
 tmp = dvrr_stack + 179463;
 target_ptr = Libderiv->deriv2_classes[6][3][91];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+179743, dvrr_stack+119957, dvrr_stack+17632);
 tmp = dvrr_stack + 179743;
 target_ptr = Libderiv->deriv2_classes[6][4][91];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+180163, dvrr_stack+120545, dvrr_stack+20879);
 tmp = dvrr_stack + 180163;
 target_ptr = Libderiv->deriv2_classes[6][5][91];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+180751, dvrr_stack+121329, dvrr_stack+119957);
 tmp = dvrr_stack + 180751;
 target_ptr = Libderiv->deriv2_classes[6][6][91];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+181535, dvrr_stack+2973, dvrr_stack+96868);
 tmp = dvrr_stack + 181535;
 target_ptr = Libderiv->deriv2_classes[3][3][83];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+181635, dvrr_stack+97028, dvrr_stack+96928);
 tmp = dvrr_stack + 181635;
 target_ptr = Libderiv->deriv2_classes[3][4][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+181785, dvrr_stack+97238, dvrr_stack+2973);
 tmp = dvrr_stack + 181785;
 target_ptr = Libderiv->deriv2_classes[3][5][83];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+181995, dvrr_stack+97518, dvrr_stack+97028);
 tmp = dvrr_stack + 181995;
 target_ptr = Libderiv->deriv2_classes[3][6][83];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+97518, dvrr_stack+18832, dvrr_stack+97998);
 tmp = dvrr_stack + 97518;
 target_ptr = Libderiv->deriv2_classes[4][3][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+182275, dvrr_stack+40440, dvrr_stack+98088);
 tmp = dvrr_stack + 182275;
 target_ptr = Libderiv->deriv2_classes[4][4][83];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+182500, dvrr_stack+71018, dvrr_stack+18832);
 tmp = dvrr_stack + 182500;
 target_ptr = Libderiv->deriv2_classes[4][5][83];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+182815, dvrr_stack+98238, dvrr_stack+40440);
 tmp = dvrr_stack + 182815;
 target_ptr = Libderiv->deriv2_classes[4][6][83];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+97668, dvrr_stack+99453, dvrr_stack+99327);
 tmp = dvrr_stack + 97668;
 target_ptr = Libderiv->deriv2_classes[5][3][83];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+98238, dvrr_stack+99978, dvrr_stack+99768);
 tmp = dvrr_stack + 98238;
 target_ptr = Libderiv->deriv2_classes[5][4][83];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+183235, dvrr_stack+100419, dvrr_stack+99453);
 tmp = dvrr_stack + 183235;
 target_ptr = Libderiv->deriv2_classes[5][5][83];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+183676, dvrr_stack+77878, dvrr_stack+99978);
 tmp = dvrr_stack + 183676;
 target_ptr = Libderiv->deriv2_classes[5][6][83];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,28,1,dvrr_stack+184264, dvrr_stack+9408, dvrr_stack+87766);
 tmp = dvrr_stack + 184264;
 target_ptr = Libderiv->deriv2_classes[6][3][83];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+184544, dvrr_stack+60052, dvrr_stack+87934);
 tmp = dvrr_stack + 184544;
 target_ptr = Libderiv->deriv2_classes[6][4][83];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+77864, dvrr_stack+60640, dvrr_stack+9408);
 tmp = dvrr_stack + 77864;
 target_ptr = Libderiv->deriv2_classes[6][5][83];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+184964, dvrr_stack+61424, dvrr_stack+60052);
 tmp = dvrr_stack + 184964;
 target_ptr = Libderiv->deriv2_classes[6][6][83];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+61424, dvrr_stack+88274, dvrr_stack+88214);
 tmp = dvrr_stack + 61424;
 target_ptr = Libderiv->deriv2_classes[3][3][82];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+61524, dvrr_stack+52399, dvrr_stack+62432);
 tmp = dvrr_stack + 61524;
 target_ptr = Libderiv->deriv2_classes[3][4][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+61674, dvrr_stack+52609, dvrr_stack+88274);
 tmp = dvrr_stack + 61674;
 target_ptr = Libderiv->deriv2_classes[3][5][82];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+61884, dvrr_stack+52889, dvrr_stack+52399);
 tmp = dvrr_stack + 61884;
 target_ptr = Libderiv->deriv2_classes[3][6][82];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+52889, dvrr_stack+29378, dvrr_stack+62532);
 tmp = dvrr_stack + 52889;
 target_ptr = Libderiv->deriv2_classes[4][3][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+98553, dvrr_stack+45637, dvrr_stack+45487);
 tmp = dvrr_stack + 98553;
 target_ptr = Libderiv->deriv2_classes[4][4][82];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+185748, dvrr_stack+31898, dvrr_stack+29378);
 tmp = dvrr_stack + 185748;
 target_ptr = Libderiv->deriv2_classes[4][5][82];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+186063, dvrr_stack+32318, dvrr_stack+45637);
 tmp = dvrr_stack + 186063;
 target_ptr = Libderiv->deriv2_classes[4][6][82];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+53039, dvrr_stack+32984, dvrr_stack+32858);
 tmp = dvrr_stack + 53039;
 target_ptr = Libderiv->deriv2_classes[5][3][82];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+32318, dvrr_stack+33509, dvrr_stack+33299);
 tmp = dvrr_stack + 32318;
 target_ptr = Libderiv->deriv2_classes[5][4][82];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+186483, dvrr_stack+13473, dvrr_stack+32984);
 tmp = dvrr_stack + 186483;
 target_ptr = Libderiv->deriv2_classes[5][5][82];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+186924, dvrr_stack+116065, dvrr_stack+33509);
 tmp = dvrr_stack + 186924;
 target_ptr = Libderiv->deriv2_classes[5][6][82];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,28,1,dvrr_stack+32633, dvrr_stack+14061, dvrr_stack+33950);
 tmp = dvrr_stack + 32633;
 target_ptr = Libderiv->deriv2_classes[6][3][82];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+187512, dvrr_stack+116821, dvrr_stack+34118);
 tmp = dvrr_stack + 187512;
 target_ptr = Libderiv->deriv2_classes[6][4][82];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+116045, dvrr_stack+117409, dvrr_stack+14061);
 tmp = dvrr_stack + 116045;
 target_ptr = Libderiv->deriv2_classes[6][5][82];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+187932, dvrr_stack+118193, dvrr_stack+116821);
 tmp = dvrr_stack + 187932;
 target_ptr = Libderiv->deriv2_classes[6][6][82];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+118193, dvrr_stack+96688, dvrr_stack+62622);
 tmp = dvrr_stack + 118193;
 target_ptr = Libderiv->deriv2_classes[3][3][81];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+62532, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 62532;
 target_ptr = Libderiv->deriv2_classes[3][4][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+118293, dvrr_stack+4439, dvrr_stack+96688);
 tmp = dvrr_stack + 118293;
 target_ptr = Libderiv->deriv2_classes[3][5][81];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+118503, dvrr_stack+7728, dvrr_stack+1023);
 tmp = dvrr_stack + 118503;
 target_ptr = Libderiv->deriv2_classes[3][6][81];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+7728, dvrr_stack+5699, dvrr_stack+11943);
 tmp = dvrr_stack + 7728;
 target_ptr = Libderiv->deriv2_classes[4][3][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+118783, dvrr_stack+16318, dvrr_stack+17482);
 tmp = dvrr_stack + 118783;
 target_ptr = Libderiv->deriv2_classes[4][4][81];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+188716, dvrr_stack+20459, dvrr_stack+5699);
 tmp = dvrr_stack + 188716;
 target_ptr = Libderiv->deriv2_classes[4][5][81];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+189031, dvrr_stack+24653, dvrr_stack+16318);
 tmp = dvrr_stack + 189031;
 target_ptr = Libderiv->deriv2_classes[4][6][81];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+7878, dvrr_stack+98949, dvrr_stack+30278);
 tmp = dvrr_stack + 7878;
 target_ptr = Libderiv->deriv2_classes[5][3][81];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+189451, dvrr_stack+36828, dvrr_stack+38550);
 tmp = dvrr_stack + 189451;
 target_ptr = Libderiv->deriv2_classes[5][4][81];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+24612, dvrr_stack+42841, dvrr_stack+98949);
 tmp = dvrr_stack + 24612;
 target_ptr = Libderiv->deriv2_classes[5][5][81];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+189766, dvrr_stack+48871, dvrr_stack+36828);
 tmp = dvrr_stack + 189766;
 target_ptr = Libderiv->deriv2_classes[5][6][81];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,28,1,dvrr_stack+190354, dvrr_stack+87262, dvrr_stack+56944);
 tmp = dvrr_stack + 190354;
 target_ptr = Libderiv->deriv2_classes[6][3][81];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+48763, dvrr_stack+66104, dvrr_stack+68498);
 tmp = dvrr_stack + 48763;
 target_ptr = Libderiv->deriv2_classes[6][4][81];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+190634, dvrr_stack+74350, dvrr_stack+87262);
 tmp = dvrr_stack + 190634;
 target_ptr = Libderiv->deriv2_classes[6][5][81];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+191222, dvrr_stack+49627, dvrr_stack+66104);
 tmp = dvrr_stack + 191222;
 target_ptr = Libderiv->deriv2_classes[6][6][81];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+119008, dvrr_stack+93592, dvrr_stack+93532);
 tmp = dvrr_stack + 119008;
 target_ptr = Libderiv->deriv2_classes[3][3][80];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+30254, dvrr_stack+93842, dvrr_stack+93742);
 tmp = dvrr_stack + 30254;
 target_ptr = Libderiv->deriv2_classes[3][4][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+62164, dvrr_stack+94052, dvrr_stack+93592);
 tmp = dvrr_stack + 62164;
 target_ptr = Libderiv->deriv2_classes[3][5][80];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+192006, dvrr_stack+94332, dvrr_stack+93842);
 tmp = dvrr_stack + 192006;
 target_ptr = Libderiv->deriv2_classes[3][6][80];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+94332, dvrr_stack+94782, dvrr_stack+94692);
 tmp = dvrr_stack + 94332;
 target_ptr = Libderiv->deriv2_classes[4][3][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+94482, dvrr_stack+95157, dvrr_stack+95007);
 tmp = dvrr_stack + 94482;
 target_ptr = Libderiv->deriv2_classes[4][4][80];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+192286, dvrr_stack+95472, dvrr_stack+94782);
 tmp = dvrr_stack + 192286;
 target_ptr = Libderiv->deriv2_classes[4][5][80];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+49183, dvrr_stack+95892, dvrr_stack+95157);
 tmp = dvrr_stack + 49183;
 target_ptr = Libderiv->deriv2_classes[4][6][80];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+95892, dvrr_stack+82558, dvrr_stack+96432);
 tmp = dvrr_stack + 95892;
 target_ptr = Libderiv->deriv2_classes[5][3][80];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+96102, dvrr_stack+83083, dvrr_stack+82873);
 tmp = dvrr_stack + 96102;
 target_ptr = Libderiv->deriv2_classes[5][4][80];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+49603, dvrr_stack+83524, dvrr_stack+82558);
 tmp = dvrr_stack + 49603;
 target_ptr = Libderiv->deriv2_classes[5][5][80];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+50044, dvrr_stack+84112, dvrr_stack+83083);
 tmp = dvrr_stack + 50044;
 target_ptr = Libderiv->deriv2_classes[5][6][80];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,28,1,dvrr_stack+84112, dvrr_stack+75302, dvrr_stack+75134);
 tmp = dvrr_stack + 84112;
 target_ptr = Libderiv->deriv2_classes[6][3][80];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+84392, dvrr_stack+68778, dvrr_stack+75722);
 tmp = dvrr_stack + 84392;
 target_ptr = Libderiv->deriv2_classes[6][4][80];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+192601, dvrr_stack+57112, dvrr_stack+75302);
 tmp = dvrr_stack + 192601;
 target_ptr = Libderiv->deriv2_classes[6][5][80];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+193189, dvrr_stack+57896, dvrr_stack+68778);
 tmp = dvrr_stack + 193189;
 target_ptr = Libderiv->deriv2_classes[6][6][80];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+57896, dvrr_stack+66692, dvrr_stack+76002);
 tmp = dvrr_stack + 57896;
 target_ptr = Libderiv->deriv2_classes[3][3][79];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+57996, dvrr_stack+69366, dvrr_stack+66842);
 tmp = dvrr_stack + 57996;
 target_ptr = Libderiv->deriv2_classes[3][4][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+58146, dvrr_stack+58904, dvrr_stack+66692);
 tmp = dvrr_stack + 58146;
 target_ptr = Libderiv->deriv2_classes[3][5][79];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+58356, dvrr_stack+43429, dvrr_stack+69366);
 tmp = dvrr_stack + 58356;
 target_ptr = Libderiv->deriv2_classes[3][6][79];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+43429, dvrr_stack+43789, dvrr_stack+69576);
 tmp = dvrr_stack + 43429;
 target_ptr = Libderiv->deriv2_classes[4][3][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+58636, dvrr_stack+38760, dvrr_stack+44014);
 tmp = dvrr_stack + 58636;
 target_ptr = Libderiv->deriv2_classes[4][4][79];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+193973, dvrr_stack+39075, dvrr_stack+43789);
 tmp = dvrr_stack + 193973;
 target_ptr = Libderiv->deriv2_classes[4][5][79];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+194288, dvrr_stack+30404, dvrr_stack+38760);
 tmp = dvrr_stack + 194288;
 target_ptr = Libderiv->deriv2_classes[4][6][79];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+43579, dvrr_stack+30944, dvrr_stack+37269);
 tmp = dvrr_stack + 43579;
 target_ptr = Libderiv->deriv2_classes[5][3][79];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+30404, dvrr_stack+25193, dvrr_stack+31259);
 tmp = dvrr_stack + 30404;
 target_ptr = Libderiv->deriv2_classes[5][4][79];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+194708, dvrr_stack+12033, dvrr_stack+30944);
 tmp = dvrr_stack + 194708;
 target_ptr = Libderiv->deriv2_classes[5][5][79];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+195149, dvrr_stack+119201, dvrr_stack+25193);
 tmp = dvrr_stack + 195149;
 target_ptr = Libderiv->deriv2_classes[5][6][79];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,28,1,dvrr_stack+119108, dvrr_stack+20879, dvrr_stack+31469);
 tmp = dvrr_stack + 119108;
 target_ptr = Libderiv->deriv2_classes[6][3][79];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+119388, dvrr_stack+119957, dvrr_stack+17632);
 tmp = dvrr_stack + 119388;
 target_ptr = Libderiv->deriv2_classes[6][4][79];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+195737, dvrr_stack+120545, dvrr_stack+20879);
 tmp = dvrr_stack + 195737;
 target_ptr = Libderiv->deriv2_classes[6][5][79];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+196325, dvrr_stack+121329, dvrr_stack+119957);
 tmp = dvrr_stack + 196325;
 target_ptr = Libderiv->deriv2_classes[6][6][79];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+121329, dvrr_stack+31637, dvrr_stack+69666);
 tmp = dvrr_stack + 121329;
 target_ptr = Libderiv->deriv2_classes[3][3][78];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+69576, dvrr_stack+25634, dvrr_stack+31787);
 tmp = dvrr_stack + 69576;
 target_ptr = Libderiv->deriv2_classes[3][4][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+121429, dvrr_stack+12621, dvrr_stack+31637);
 tmp = dvrr_stack + 121429;
 target_ptr = Libderiv->deriv2_classes[3][5][78];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+121639, dvrr_stack+7368, dvrr_stack+25634);
 tmp = dvrr_stack + 121639;
 target_ptr = Libderiv->deriv2_classes[3][6][78];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+7368, dvrr_stack+11583, dvrr_stack+11493);
 tmp = dvrr_stack + 7368;
 target_ptr = Libderiv->deriv2_classes[4][3][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+30719, dvrr_stack+8088, dvrr_stack+17912);
 tmp = dvrr_stack + 30719;
 target_ptr = Libderiv->deriv2_classes[4][4][78];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+121919, dvrr_stack+122337, dvrr_stack+11583);
 tmp = dvrr_stack + 121919;
 target_ptr = Libderiv->deriv2_classes[4][5][78];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+197109, dvrr_stack+23757, dvrr_stack+8088);
 tmp = dvrr_stack + 197109;
 target_ptr = Libderiv->deriv2_classes[4][6][78];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+7518, dvrr_stack+29729, dvrr_stack+29603);
 tmp = dvrr_stack + 7518;
 target_ptr = Libderiv->deriv2_classes[5][3][78];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+197529, dvrr_stack+122757, dvrr_stack+30044);
 tmp = dvrr_stack + 197529;
 target_ptr = Libderiv->deriv2_classes[5][4][78];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+23745, dvrr_stack+123198, dvrr_stack+29729);
 tmp = dvrr_stack + 23745;
 target_ptr = Libderiv->deriv2_classes[5][5][78];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+197844, dvrr_stack+123786, dvrr_stack+122757);
 tmp = dvrr_stack + 197844;
 target_ptr = Libderiv->deriv2_classes[5][6][78];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,28,1,dvrr_stack+123786, dvrr_stack+56167, dvrr_stack+55999);
 tmp = dvrr_stack + 123786;
 target_ptr = Libderiv->deriv2_classes[6][3][78];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+124066, dvrr_stack+124542, dvrr_stack+56587);
 tmp = dvrr_stack + 124066;
 target_ptr = Libderiv->deriv2_classes[6][4][78];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+198432, dvrr_stack+125130, dvrr_stack+56167);
 tmp = dvrr_stack + 198432;
 target_ptr = Libderiv->deriv2_classes[6][5][78];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+199020, dvrr_stack+125914, dvrr_stack+124542);
 tmp = dvrr_stack + 199020;
 target_ptr = Libderiv->deriv2_classes[6][6][78];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,10,dvrr_stack+125914, dvrr_stack+98088, dvrr_stack+96628);
 tmp = dvrr_stack + 125914;
 target_ptr = Libderiv->deriv2_classes[3][3][35];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,15,dvrr_stack+126014, dvrr_stack+18832, dvrr_stack+11808);
 tmp = dvrr_stack + 126014;
 target_ptr = Libderiv->deriv2_classes[3][4][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,21,dvrr_stack+126164, dvrr_stack+40440, dvrr_stack+16633);
 tmp = dvrr_stack + 126164;
 target_ptr = Libderiv->deriv2_classes[3][5][35];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,28,dvrr_stack+126374, dvrr_stack+71018, dvrr_stack+4899);
 tmp = dvrr_stack + 126374;
 target_ptr = Libderiv->deriv2_classes[3][6][35];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,10,dvrr_stack+126654, dvrr_stack+99768, dvrr_stack+96928);
 tmp = dvrr_stack + 126654;
 target_ptr = Libderiv->deriv2_classes[4][3][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,15,dvrr_stack+87682, dvrr_stack+99453, dvrr_stack+2973);
 tmp = dvrr_stack + 87682;
 target_ptr = Libderiv->deriv2_classes[4][4][35];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,21,dvrr_stack+199804, dvrr_stack+99978, dvrr_stack+97028);
 tmp = dvrr_stack + 199804;
 target_ptr = Libderiv->deriv2_classes[4][5][35];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,28,dvrr_stack+200119, dvrr_stack+100419, dvrr_stack+97238);
 tmp = dvrr_stack + 200119;
 target_ptr = Libderiv->deriv2_classes[4][6][35];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_h(Data,10,dvrr_stack+96417, dvrr_stack+87934, dvrr_stack+98088);
 tmp = dvrr_stack + 96417;
 target_ptr = Libderiv->deriv2_classes[5][3][35];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_h(Data,15,dvrr_stack+200539, dvrr_stack+9408, dvrr_stack+18832);
 tmp = dvrr_stack + 200539;
 target_ptr = Libderiv->deriv2_classes[5][4][35];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_h(Data,21,dvrr_stack+200854, dvrr_stack+60052, dvrr_stack+40440);
 tmp = dvrr_stack + 200854;
 target_ptr = Libderiv->deriv2_classes[5][5][35];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_h(Data,28,dvrr_stack+201295, dvrr_stack+60640, dvrr_stack+71018);
 tmp = dvrr_stack + 201295;
 target_ptr = Libderiv->deriv2_classes[5][6][35];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_i(Data,10,dvrr_stack+201883, dvrr_stack+36138, dvrr_stack+99768);
 tmp = dvrr_stack + 201883;
 target_ptr = Libderiv->deriv2_classes[6][3][35];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_i(Data,15,dvrr_stack+202163, dvrr_stack+37455, dvrr_stack+99453);
 tmp = dvrr_stack + 202163;
 target_ptr = Libderiv->deriv2_classes[6][4][35];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_i(Data,21,dvrr_stack+202583, dvrr_stack+41385, dvrr_stack+99978);
 tmp = dvrr_stack + 202583;
 target_ptr = Libderiv->deriv2_classes[6][5][35];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_i(Data,28,dvrr_stack+203171, dvrr_stack+138390, dvrr_stack+100419);
 tmp = dvrr_stack + 203171;
 target_ptr = Libderiv->deriv2_classes[6][6][35];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+126804, dvrr_stack+45487, dvrr_stack+79894);
 tmp = dvrr_stack + 126804;
 target_ptr = Libderiv->deriv2_classes[3][3][34];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+31469, dvrr_stack+29378, dvrr_stack+6708);
 tmp = dvrr_stack + 31469;
 target_ptr = Libderiv->deriv2_classes[3][4][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+56867, dvrr_stack+45637, dvrr_stack+79954);
 tmp = dvrr_stack + 56867;
 target_ptr = Libderiv->deriv2_classes[3][5][34];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+203955, dvrr_stack+31898, dvrr_stack+80080);
 tmp = dvrr_stack + 203955;
 target_ptr = Libderiv->deriv2_classes[3][6][34];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+75134, dvrr_stack+33299, dvrr_stack+62432);
 tmp = dvrr_stack + 75134;
 target_ptr = Libderiv->deriv2_classes[4][3][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+204235, dvrr_stack+32984, dvrr_stack+88274);
 tmp = dvrr_stack + 204235;
 target_ptr = Libderiv->deriv2_classes[4][4][34];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+204460, dvrr_stack+33509, dvrr_stack+52399);
 tmp = dvrr_stack + 204460;
 target_ptr = Libderiv->deriv2_classes[4][5][34];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+204775, dvrr_stack+13473, dvrr_stack+52609);
 tmp = dvrr_stack + 204775;
 target_ptr = Libderiv->deriv2_classes[4][6][34];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+205195, dvrr_stack+34118, dvrr_stack+45487);
 tmp = dvrr_stack + 205195;
 target_ptr = Libderiv->deriv2_classes[5][3][34];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+205405, dvrr_stack+14061, dvrr_stack+29378);
 tmp = dvrr_stack + 205405;
 target_ptr = Libderiv->deriv2_classes[5][4][34];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+205720, dvrr_stack+116821, dvrr_stack+45637);
 tmp = dvrr_stack + 205720;
 target_ptr = Libderiv->deriv2_classes[5][5][34];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+206161, dvrr_stack+117409, dvrr_stack+31898);
 tmp = dvrr_stack + 206161;
 target_ptr = Libderiv->deriv2_classes[5][6][34];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_i(Data,10,dvrr_stack+206749, dvrr_stack+80248, dvrr_stack+33299);
 tmp = dvrr_stack + 206749;
 target_ptr = Libderiv->deriv2_classes[6][3][34];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+207029, dvrr_stack+80608, dvrr_stack+32984);
 tmp = dvrr_stack + 207029;
 target_ptr = Libderiv->deriv2_classes[6][4][34];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+207449, dvrr_stack+81148, dvrr_stack+33509);
 tmp = dvrr_stack + 207449;
 target_ptr = Libderiv->deriv2_classes[6][5][34];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+208037, dvrr_stack+46999, dvrr_stack+13473);
 tmp = dvrr_stack + 208037;
 target_ptr = Libderiv->deriv2_classes[6][6][34];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+29603, dvrr_stack+17482, dvrr_stack+81904);
 tmp = dvrr_stack + 29603;
 target_ptr = Libderiv->deriv2_classes[3][3][33];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+33950, dvrr_stack+5699, dvrr_stack+4719);
 tmp = dvrr_stack + 33950;
 target_ptr = Libderiv->deriv2_classes[3][4][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+208821, dvrr_stack+16318, dvrr_stack+2173);
 tmp = dvrr_stack + 208821;
 target_ptr = Libderiv->deriv2_classes[3][5][33];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+209031, dvrr_stack+20459, dvrr_stack+22349);
 tmp = dvrr_stack + 209031;
 target_ptr = Libderiv->deriv2_classes[3][6][33];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+209311, dvrr_stack+38550, dvrr_stack+2073);
 tmp = dvrr_stack + 209311;
 target_ptr = Libderiv->deriv2_classes[4][3][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+209461, dvrr_stack+98949, dvrr_stack+96688);
 tmp = dvrr_stack + 209461;
 target_ptr = Libderiv->deriv2_classes[4][4][33];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+209686, dvrr_stack+36828, dvrr_stack+1023);
 tmp = dvrr_stack + 209686;
 target_ptr = Libderiv->deriv2_classes[4][5][33];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+210001, dvrr_stack+42841, dvrr_stack+4439);
 tmp = dvrr_stack + 210001;
 target_ptr = Libderiv->deriv2_classes[4][6][33];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+210421, dvrr_stack+68498, dvrr_stack+17482);
 tmp = dvrr_stack + 210421;
 target_ptr = Libderiv->deriv2_classes[5][3][33];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+210631, dvrr_stack+87262, dvrr_stack+5699);
 tmp = dvrr_stack + 210631;
 target_ptr = Libderiv->deriv2_classes[5][4][33];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+210946, dvrr_stack+66104, dvrr_stack+16318);
 tmp = dvrr_stack + 210946;
 target_ptr = Libderiv->deriv2_classes[5][5][33];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+211387, dvrr_stack+74350, dvrr_stack+20459);
 tmp = dvrr_stack + 211387;
 target_ptr = Libderiv->deriv2_classes[5][6][33];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_i(Data,10,dvrr_stack+211975, dvrr_stack+92167, dvrr_stack+38550);
 tmp = dvrr_stack + 211975;
 target_ptr = Libderiv->deriv2_classes[6][3][33];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+212255, dvrr_stack+128242, dvrr_stack+98949);
 tmp = dvrr_stack + 212255;
 target_ptr = Libderiv->deriv2_classes[6][4][33];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+212675, dvrr_stack+131422, dvrr_stack+36828);
 tmp = dvrr_stack + 212675;
 target_ptr = Libderiv->deriv2_classes[6][5][33];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+213263, dvrr_stack+136122, dvrr_stack+42841);
 tmp = dvrr_stack + 213263;
 target_ptr = Libderiv->deriv2_classes[6][6][33];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+122234, dvrr_stack+95007, dvrr_stack+140694);
 tmp = dvrr_stack + 122234;
 target_ptr = Libderiv->deriv2_classes[3][3][32];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+214047, dvrr_stack+94782, dvrr_stack+4809);
 tmp = dvrr_stack + 214047;
 target_ptr = Libderiv->deriv2_classes[3][4][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+214197, dvrr_stack+95157, dvrr_stack+140754);
 tmp = dvrr_stack + 214197;
 target_ptr = Libderiv->deriv2_classes[3][5][32];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+214407, dvrr_stack+95472, dvrr_stack+140880);
 tmp = dvrr_stack + 214407;
 target_ptr = Libderiv->deriv2_classes[3][6][32];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+214687, dvrr_stack+82873, dvrr_stack+93742);
 tmp = dvrr_stack + 214687;
 target_ptr = Libderiv->deriv2_classes[4][3][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+214837, dvrr_stack+82558, dvrr_stack+93592);
 tmp = dvrr_stack + 214837;
 target_ptr = Libderiv->deriv2_classes[4][4][32];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+215062, dvrr_stack+83083, dvrr_stack+93842);
 tmp = dvrr_stack + 215062;
 target_ptr = Libderiv->deriv2_classes[4][5][32];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+215377, dvrr_stack+83524, dvrr_stack+94052);
 tmp = dvrr_stack + 215377;
 target_ptr = Libderiv->deriv2_classes[4][6][32];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+215797, dvrr_stack+75722, dvrr_stack+95007);
 tmp = dvrr_stack + 215797;
 target_ptr = Libderiv->deriv2_classes[5][3][32];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+216007, dvrr_stack+75302, dvrr_stack+94782);
 tmp = dvrr_stack + 216007;
 target_ptr = Libderiv->deriv2_classes[5][4][32];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+216322, dvrr_stack+68778, dvrr_stack+95157);
 tmp = dvrr_stack + 216322;
 target_ptr = Libderiv->deriv2_classes[5][5][32];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+216763, dvrr_stack+57112, dvrr_stack+95472);
 tmp = dvrr_stack + 216763;
 target_ptr = Libderiv->deriv2_classes[5][6][32];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_i(Data,10,dvrr_stack+217351, dvrr_stack+141369, dvrr_stack+82873);
 tmp = dvrr_stack + 217351;
 target_ptr = Libderiv->deriv2_classes[6][3][32];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+217631, dvrr_stack+128782, dvrr_stack+82558);
 tmp = dvrr_stack + 217631;
 target_ptr = Libderiv->deriv2_classes[6][4][32];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+218051, dvrr_stack+141729, dvrr_stack+83083);
 tmp = dvrr_stack + 218051;
 target_ptr = Libderiv->deriv2_classes[6][5][32];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+218639, dvrr_stack+142485, dvrr_stack+83524);
 tmp = dvrr_stack + 218639;
 target_ptr = Libderiv->deriv2_classes[6][6][32];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+37269, dvrr_stack+44014, dvrr_stack+141048);
 tmp = dvrr_stack + 37269;
 target_ptr = Libderiv->deriv2_classes[3][3][31];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+99264, dvrr_stack+43789, dvrr_stack+98778);
 tmp = dvrr_stack + 99264;
 target_ptr = Libderiv->deriv2_classes[3][4][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+219423, dvrr_stack+38760, dvrr_stack+101007);
 tmp = dvrr_stack + 219423;
 target_ptr = Libderiv->deriv2_classes[3][5][31];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+219633, dvrr_stack+39075, dvrr_stack+143493);
 tmp = dvrr_stack + 219633;
 target_ptr = Libderiv->deriv2_classes[3][6][31];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+219913, dvrr_stack+31259, dvrr_stack+66842);
 tmp = dvrr_stack + 219913;
 target_ptr = Libderiv->deriv2_classes[4][3][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+220063, dvrr_stack+30944, dvrr_stack+66692);
 tmp = dvrr_stack + 220063;
 target_ptr = Libderiv->deriv2_classes[4][4][31];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+220288, dvrr_stack+25193, dvrr_stack+69366);
 tmp = dvrr_stack + 220288;
 target_ptr = Libderiv->deriv2_classes[4][5][31];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+220603, dvrr_stack+12033, dvrr_stack+58904);
 tmp = dvrr_stack + 220603;
 target_ptr = Libderiv->deriv2_classes[4][6][31];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+221023, dvrr_stack+17632, dvrr_stack+44014);
 tmp = dvrr_stack + 221023;
 target_ptr = Libderiv->deriv2_classes[5][3][31];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+221233, dvrr_stack+20879, dvrr_stack+43789);
 tmp = dvrr_stack + 221233;
 target_ptr = Libderiv->deriv2_classes[5][4][31];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+221548, dvrr_stack+119957, dvrr_stack+38760);
 tmp = dvrr_stack + 221548;
 target_ptr = Libderiv->deriv2_classes[5][5][31];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+221989, dvrr_stack+120545, dvrr_stack+39075);
 tmp = dvrr_stack + 221989;
 target_ptr = Libderiv->deriv2_classes[5][6][31];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_i(Data,10,dvrr_stack+222577, dvrr_stack+137130, dvrr_stack+31259);
 tmp = dvrr_stack + 222577;
 target_ptr = Libderiv->deriv2_classes[6][3][31];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+222857, dvrr_stack+137490, dvrr_stack+30944);
 tmp = dvrr_stack + 222857;
 target_ptr = Libderiv->deriv2_classes[6][4][31];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+223277, dvrr_stack+132178, dvrr_stack+25193);
 tmp = dvrr_stack + 223277;
 target_ptr = Libderiv->deriv2_classes[6][5][31];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+223865, dvrr_stack+143661, dvrr_stack+12033);
 tmp = dvrr_stack + 223865;
 target_ptr = Libderiv->deriv2_classes[6][6][31];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+119808, dvrr_stack+17912, dvrr_stack+98868);
 tmp = dvrr_stack + 119808;
 target_ptr = Libderiv->deriv2_classes[3][3][30];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+224649, dvrr_stack+11583, dvrr_stack+101133);
 tmp = dvrr_stack + 224649;
 target_ptr = Libderiv->deriv2_classes[3][4][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+224799, dvrr_stack+8088, dvrr_stack+144669);
 tmp = dvrr_stack + 224799;
 target_ptr = Libderiv->deriv2_classes[3][5][30];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+225009, dvrr_stack+122337, dvrr_stack+144795);
 tmp = dvrr_stack + 225009;
 target_ptr = Libderiv->deriv2_classes[3][6][30];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+225289, dvrr_stack+30044, dvrr_stack+31787);
 tmp = dvrr_stack + 225289;
 target_ptr = Libderiv->deriv2_classes[4][3][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+225439, dvrr_stack+29729, dvrr_stack+31637);
 tmp = dvrr_stack + 225439;
 target_ptr = Libderiv->deriv2_classes[4][4][30];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+225664, dvrr_stack+122757, dvrr_stack+25634);
 tmp = dvrr_stack + 225664;
 target_ptr = Libderiv->deriv2_classes[4][5][30];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+225979, dvrr_stack+123198, dvrr_stack+12621);
 tmp = dvrr_stack + 225979;
 target_ptr = Libderiv->deriv2_classes[4][6][30];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+226399, dvrr_stack+56587, dvrr_stack+17912);
 tmp = dvrr_stack + 226399;
 target_ptr = Libderiv->deriv2_classes[5][3][30];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+226609, dvrr_stack+56167, dvrr_stack+11583);
 tmp = dvrr_stack + 226609;
 target_ptr = Libderiv->deriv2_classes[5][4][30];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+226924, dvrr_stack+124542, dvrr_stack+8088);
 tmp = dvrr_stack + 226924;
 target_ptr = Libderiv->deriv2_classes[5][5][30];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+227365, dvrr_stack+125130, dvrr_stack+122337);
 tmp = dvrr_stack + 227365;
 target_ptr = Libderiv->deriv2_classes[5][6][30];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_i(Data,10,dvrr_stack+227953, dvrr_stack+138030, dvrr_stack+30044);
 tmp = dvrr_stack + 227953;
 target_ptr = Libderiv->deriv2_classes[6][3][30];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+228233, dvrr_stack+81964, dvrr_stack+29729);
 tmp = dvrr_stack + 228233;
 target_ptr = Libderiv->deriv2_classes[6][4][30];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+228653, dvrr_stack+48007, dvrr_stack+122757);
 tmp = dvrr_stack + 228653;
 target_ptr = Libderiv->deriv2_classes[6][5][30];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+229241, dvrr_stack+144963, dvrr_stack+123198);
 tmp = dvrr_stack + 229241;
 target_ptr = Libderiv->deriv2_classes[6][6][30];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+230025, dvrr_stack+139458, dvrr_stack+139398);
 tmp = dvrr_stack + 230025;
 target_ptr = Libderiv->deriv2_classes[3][3][26];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+230125, dvrr_stack+139698, dvrr_stack+139608);
 tmp = dvrr_stack + 230125;
 target_ptr = Libderiv->deriv2_classes[3][4][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+230275, dvrr_stack+140049, dvrr_stack+139923);
 tmp = dvrr_stack + 230275;
 target_ptr = Libderiv->deriv2_classes[3][5][26];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+230485, dvrr_stack+145971, dvrr_stack+140364);
 tmp = dvrr_stack + 230485;
 target_ptr = Libderiv->deriv2_classes[3][6][26];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+230765, dvrr_stack+141108, dvrr_stack+35828);
 tmp = dvrr_stack + 230765;
 target_ptr = Libderiv->deriv2_classes[4][3][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+230915, dvrr_stack+146391, dvrr_stack+140532);
 tmp = dvrr_stack + 230915;
 target_ptr = Libderiv->deriv2_classes[4][4][26];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+231140, dvrr_stack+146916, dvrr_stack+146706);
 tmp = dvrr_stack + 231140;
 target_ptr = Libderiv->deriv2_classes[4][5][26];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+231455, dvrr_stack+147637, dvrr_stack+147357);
 tmp = dvrr_stack + 231455;
 target_ptr = Libderiv->deriv2_classes[4][6][26];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+231875, dvrr_stack+148225, dvrr_stack+139458);
 tmp = dvrr_stack + 231875;
 target_ptr = Libderiv->deriv2_classes[5][3][26];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+232085, dvrr_stack+148505, dvrr_stack+139698);
 tmp = dvrr_stack + 232085;
 target_ptr = Libderiv->deriv2_classes[5][4][26];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+232400, dvrr_stack+148925, dvrr_stack+140049);
 tmp = dvrr_stack + 232400;
 target_ptr = Libderiv->deriv2_classes[5][5][26];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+232841, dvrr_stack+149513, dvrr_stack+145971);
 tmp = dvrr_stack + 232841;
 target_ptr = Libderiv->deriv2_classes[5][6][26];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,10,dvrr_stack+233429, dvrr_stack+150297, dvrr_stack+141108);
 tmp = dvrr_stack + 233429;
 target_ptr = Libderiv->deriv2_classes[6][3][26];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+233709, dvrr_stack+151617, dvrr_stack+146391);
 tmp = dvrr_stack + 233709;
 target_ptr = Libderiv->deriv2_classes[6][4][26];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+234129, dvrr_stack+150657, dvrr_stack+146916);
 tmp = dvrr_stack + 234129;
 target_ptr = Libderiv->deriv2_classes[6][5][26];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+234717, dvrr_stack+129322, dvrr_stack+147637);
 tmp = dvrr_stack + 234717;
 target_ptr = Libderiv->deriv2_classes[6][6][26];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,10,dvrr_stack+235501, dvrr_stack+98088, dvrr_stack+96628);
 tmp = dvrr_stack + 235501;
 target_ptr = Libderiv->deriv2_classes[3][3][23];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,15,dvrr_stack+235601, dvrr_stack+18832, dvrr_stack+11808);
 tmp = dvrr_stack + 235601;
 target_ptr = Libderiv->deriv2_classes[3][4][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,21,dvrr_stack+235751, dvrr_stack+40440, dvrr_stack+16633);
 tmp = dvrr_stack + 235751;
 target_ptr = Libderiv->deriv2_classes[3][5][23];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,28,dvrr_stack+235961, dvrr_stack+71018, dvrr_stack+4899);
 tmp = dvrr_stack + 235961;
 target_ptr = Libderiv->deriv2_classes[3][6][23];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,10,dvrr_stack+236241, dvrr_stack+99768, dvrr_stack+96928);
 tmp = dvrr_stack + 236241;
 target_ptr = Libderiv->deriv2_classes[4][3][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,15,dvrr_stack+236391, dvrr_stack+99453, dvrr_stack+2973);
 tmp = dvrr_stack + 236391;
 target_ptr = Libderiv->deriv2_classes[4][4][23];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,21,dvrr_stack+236616, dvrr_stack+99978, dvrr_stack+97028);
 tmp = dvrr_stack + 236616;
 target_ptr = Libderiv->deriv2_classes[4][5][23];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,28,dvrr_stack+236931, dvrr_stack+100419, dvrr_stack+97238);
 tmp = dvrr_stack + 236931;
 target_ptr = Libderiv->deriv2_classes[4][6][23];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_h(Data,10,dvrr_stack+237351, dvrr_stack+87934, dvrr_stack+98088);
 tmp = dvrr_stack + 237351;
 target_ptr = Libderiv->deriv2_classes[5][3][23];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_h(Data,15,dvrr_stack+237561, dvrr_stack+9408, dvrr_stack+18832);
 tmp = dvrr_stack + 237561;
 target_ptr = Libderiv->deriv2_classes[5][4][23];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_h(Data,21,dvrr_stack+237876, dvrr_stack+60052, dvrr_stack+40440);
 tmp = dvrr_stack + 237876;
 target_ptr = Libderiv->deriv2_classes[5][5][23];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_h(Data,28,dvrr_stack+238317, dvrr_stack+60640, dvrr_stack+71018);
 tmp = dvrr_stack + 238317;
 target_ptr = Libderiv->deriv2_classes[5][6][23];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_i(Data,10,dvrr_stack+238905, dvrr_stack+36138, dvrr_stack+99768);
 tmp = dvrr_stack + 238905;
 target_ptr = Libderiv->deriv2_classes[6][3][23];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_i(Data,15,dvrr_stack+239185, dvrr_stack+37455, dvrr_stack+99453);
 tmp = dvrr_stack + 239185;
 target_ptr = Libderiv->deriv2_classes[6][4][23];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_i(Data,21,dvrr_stack+239605, dvrr_stack+41385, dvrr_stack+99978);
 tmp = dvrr_stack + 239605;
 target_ptr = Libderiv->deriv2_classes[6][5][23];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_i(Data,28,dvrr_stack+240193, dvrr_stack+138390, dvrr_stack+100419);
 tmp = dvrr_stack + 240193;
 target_ptr = Libderiv->deriv2_classes[6][6][23];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+11466, dvrr_stack+45487, dvrr_stack+79894);
 tmp = dvrr_stack + 11466;
 target_ptr = Libderiv->deriv2_classes[3][3][22];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+116633, dvrr_stack+29378, dvrr_stack+6708);
 tmp = dvrr_stack + 116633;
 target_ptr = Libderiv->deriv2_classes[3][4][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+240977, dvrr_stack+45637, dvrr_stack+79954);
 tmp = dvrr_stack + 240977;
 target_ptr = Libderiv->deriv2_classes[3][5][22];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+241187, dvrr_stack+31898, dvrr_stack+80080);
 tmp = dvrr_stack + 241187;
 target_ptr = Libderiv->deriv2_classes[3][6][22];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+241467, dvrr_stack+33299, dvrr_stack+62432);
 tmp = dvrr_stack + 241467;
 target_ptr = Libderiv->deriv2_classes[4][3][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+241617, dvrr_stack+32984, dvrr_stack+88274);
 tmp = dvrr_stack + 241617;
 target_ptr = Libderiv->deriv2_classes[4][4][22];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+241842, dvrr_stack+33509, dvrr_stack+52399);
 tmp = dvrr_stack + 241842;
 target_ptr = Libderiv->deriv2_classes[4][5][22];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+242157, dvrr_stack+13473, dvrr_stack+52609);
 tmp = dvrr_stack + 242157;
 target_ptr = Libderiv->deriv2_classes[4][6][22];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+242577, dvrr_stack+34118, dvrr_stack+45487);
 tmp = dvrr_stack + 242577;
 target_ptr = Libderiv->deriv2_classes[5][3][22];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+242787, dvrr_stack+14061, dvrr_stack+29378);
 tmp = dvrr_stack + 242787;
 target_ptr = Libderiv->deriv2_classes[5][4][22];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+243102, dvrr_stack+116821, dvrr_stack+45637);
 tmp = dvrr_stack + 243102;
 target_ptr = Libderiv->deriv2_classes[5][5][22];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+243543, dvrr_stack+117409, dvrr_stack+31898);
 tmp = dvrr_stack + 243543;
 target_ptr = Libderiv->deriv2_classes[5][6][22];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_i(Data,10,dvrr_stack+244131, dvrr_stack+80248, dvrr_stack+33299);
 tmp = dvrr_stack + 244131;
 target_ptr = Libderiv->deriv2_classes[6][3][22];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+244411, dvrr_stack+80608, dvrr_stack+32984);
 tmp = dvrr_stack + 244411;
 target_ptr = Libderiv->deriv2_classes[6][4][22];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+244831, dvrr_stack+81148, dvrr_stack+33509);
 tmp = dvrr_stack + 244831;
 target_ptr = Libderiv->deriv2_classes[6][5][22];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+245419, dvrr_stack+46999, dvrr_stack+13473);
 tmp = dvrr_stack + 245419;
 target_ptr = Libderiv->deriv2_classes[6][6][22];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+246203, dvrr_stack+17482, dvrr_stack+81904);
 tmp = dvrr_stack + 246203;
 target_ptr = Libderiv->deriv2_classes[3][3][21];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+246303, dvrr_stack+5699, dvrr_stack+4719);
 tmp = dvrr_stack + 246303;
 target_ptr = Libderiv->deriv2_classes[3][4][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+246453, dvrr_stack+16318, dvrr_stack+2173);
 tmp = dvrr_stack + 246453;
 target_ptr = Libderiv->deriv2_classes[3][5][21];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+246663, dvrr_stack+20459, dvrr_stack+22349);
 tmp = dvrr_stack + 246663;
 target_ptr = Libderiv->deriv2_classes[3][6][21];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+93397, dvrr_stack+38550, dvrr_stack+2073);
 tmp = dvrr_stack + 93397;
 target_ptr = Libderiv->deriv2_classes[4][3][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+246943, dvrr_stack+98949, dvrr_stack+96688);
 tmp = dvrr_stack + 246943;
 target_ptr = Libderiv->deriv2_classes[4][4][21];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+247168, dvrr_stack+36828, dvrr_stack+1023);
 tmp = dvrr_stack + 247168;
 target_ptr = Libderiv->deriv2_classes[4][5][21];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+247483, dvrr_stack+42841, dvrr_stack+4439);
 tmp = dvrr_stack + 247483;
 target_ptr = Libderiv->deriv2_classes[4][6][21];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+247903, dvrr_stack+68498, dvrr_stack+17482);
 tmp = dvrr_stack + 247903;
 target_ptr = Libderiv->deriv2_classes[5][3][21];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+248113, dvrr_stack+87262, dvrr_stack+5699);
 tmp = dvrr_stack + 248113;
 target_ptr = Libderiv->deriv2_classes[5][4][21];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+248428, dvrr_stack+66104, dvrr_stack+16318);
 tmp = dvrr_stack + 248428;
 target_ptr = Libderiv->deriv2_classes[5][5][21];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+248869, dvrr_stack+74350, dvrr_stack+20459);
 tmp = dvrr_stack + 248869;
 target_ptr = Libderiv->deriv2_classes[5][6][21];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_i(Data,10,dvrr_stack+249457, dvrr_stack+92167, dvrr_stack+38550);
 tmp = dvrr_stack + 249457;
 target_ptr = Libderiv->deriv2_classes[6][3][21];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+249737, dvrr_stack+128242, dvrr_stack+98949);
 tmp = dvrr_stack + 249737;
 target_ptr = Libderiv->deriv2_classes[6][4][21];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+250157, dvrr_stack+131422, dvrr_stack+36828);
 tmp = dvrr_stack + 250157;
 target_ptr = Libderiv->deriv2_classes[6][5][21];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+250745, dvrr_stack+136122, dvrr_stack+42841);
 tmp = dvrr_stack + 250745;
 target_ptr = Libderiv->deriv2_classes[6][6][21];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+251529, dvrr_stack+95007, dvrr_stack+140694);
 tmp = dvrr_stack + 251529;
 target_ptr = Libderiv->deriv2_classes[3][3][20];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+251629, dvrr_stack+94782, dvrr_stack+4809);
 tmp = dvrr_stack + 251629;
 target_ptr = Libderiv->deriv2_classes[3][4][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+251779, dvrr_stack+95157, dvrr_stack+140754);
 tmp = dvrr_stack + 251779;
 target_ptr = Libderiv->deriv2_classes[3][5][20];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+251989, dvrr_stack+95472, dvrr_stack+140880);
 tmp = dvrr_stack + 251989;
 target_ptr = Libderiv->deriv2_classes[3][6][20];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+252269, dvrr_stack+82873, dvrr_stack+93742);
 tmp = dvrr_stack + 252269;
 target_ptr = Libderiv->deriv2_classes[4][3][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+252419, dvrr_stack+82558, dvrr_stack+93592);
 tmp = dvrr_stack + 252419;
 target_ptr = Libderiv->deriv2_classes[4][4][20];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+252644, dvrr_stack+83083, dvrr_stack+93842);
 tmp = dvrr_stack + 252644;
 target_ptr = Libderiv->deriv2_classes[4][5][20];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+252959, dvrr_stack+83524, dvrr_stack+94052);
 tmp = dvrr_stack + 252959;
 target_ptr = Libderiv->deriv2_classes[4][6][20];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+253379, dvrr_stack+75722, dvrr_stack+95007);
 tmp = dvrr_stack + 253379;
 target_ptr = Libderiv->deriv2_classes[5][3][20];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+253589, dvrr_stack+75302, dvrr_stack+94782);
 tmp = dvrr_stack + 253589;
 target_ptr = Libderiv->deriv2_classes[5][4][20];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+253904, dvrr_stack+68778, dvrr_stack+95157);
 tmp = dvrr_stack + 253904;
 target_ptr = Libderiv->deriv2_classes[5][5][20];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+254345, dvrr_stack+57112, dvrr_stack+95472);
 tmp = dvrr_stack + 254345;
 target_ptr = Libderiv->deriv2_classes[5][6][20];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_i(Data,10,dvrr_stack+254933, dvrr_stack+141369, dvrr_stack+82873);
 tmp = dvrr_stack + 254933;
 target_ptr = Libderiv->deriv2_classes[6][3][20];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+255213, dvrr_stack+128782, dvrr_stack+82558);
 tmp = dvrr_stack + 255213;
 target_ptr = Libderiv->deriv2_classes[6][4][20];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+255633, dvrr_stack+141729, dvrr_stack+83083);
 tmp = dvrr_stack + 255633;
 target_ptr = Libderiv->deriv2_classes[6][5][20];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+256221, dvrr_stack+142485, dvrr_stack+83524);
 tmp = dvrr_stack + 256221;
 target_ptr = Libderiv->deriv2_classes[6][6][20];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+257005, dvrr_stack+44014, dvrr_stack+141048);
 tmp = dvrr_stack + 257005;
 target_ptr = Libderiv->deriv2_classes[3][3][19];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+257105, dvrr_stack+43789, dvrr_stack+98778);
 tmp = dvrr_stack + 257105;
 target_ptr = Libderiv->deriv2_classes[3][4][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+257255, dvrr_stack+38760, dvrr_stack+101007);
 tmp = dvrr_stack + 257255;
 target_ptr = Libderiv->deriv2_classes[3][5][19];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+257465, dvrr_stack+39075, dvrr_stack+143493);
 tmp = dvrr_stack + 257465;
 target_ptr = Libderiv->deriv2_classes[3][6][19];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+257745, dvrr_stack+31259, dvrr_stack+66842);
 tmp = dvrr_stack + 257745;
 target_ptr = Libderiv->deriv2_classes[4][3][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+257895, dvrr_stack+30944, dvrr_stack+66692);
 tmp = dvrr_stack + 257895;
 target_ptr = Libderiv->deriv2_classes[4][4][19];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+258120, dvrr_stack+25193, dvrr_stack+69366);
 tmp = dvrr_stack + 258120;
 target_ptr = Libderiv->deriv2_classes[4][5][19];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+258435, dvrr_stack+12033, dvrr_stack+58904);
 tmp = dvrr_stack + 258435;
 target_ptr = Libderiv->deriv2_classes[4][6][19];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+258855, dvrr_stack+17632, dvrr_stack+44014);
 tmp = dvrr_stack + 258855;
 target_ptr = Libderiv->deriv2_classes[5][3][19];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+259065, dvrr_stack+20879, dvrr_stack+43789);
 tmp = dvrr_stack + 259065;
 target_ptr = Libderiv->deriv2_classes[5][4][19];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+259380, dvrr_stack+119957, dvrr_stack+38760);
 tmp = dvrr_stack + 259380;
 target_ptr = Libderiv->deriv2_classes[5][5][19];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+259821, dvrr_stack+120545, dvrr_stack+39075);
 tmp = dvrr_stack + 259821;
 target_ptr = Libderiv->deriv2_classes[5][6][19];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_i(Data,10,dvrr_stack+260409, dvrr_stack+137130, dvrr_stack+31259);
 tmp = dvrr_stack + 260409;
 target_ptr = Libderiv->deriv2_classes[6][3][19];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+260689, dvrr_stack+137490, dvrr_stack+30944);
 tmp = dvrr_stack + 260689;
 target_ptr = Libderiv->deriv2_classes[6][4][19];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+261109, dvrr_stack+132178, dvrr_stack+25193);
 tmp = dvrr_stack + 261109;
 target_ptr = Libderiv->deriv2_classes[6][5][19];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+261697, dvrr_stack+143661, dvrr_stack+12033);
 tmp = dvrr_stack + 261697;
 target_ptr = Libderiv->deriv2_classes[6][6][19];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+262481, dvrr_stack+17912, dvrr_stack+98868);
 tmp = dvrr_stack + 262481;
 target_ptr = Libderiv->deriv2_classes[3][3][18];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+262581, dvrr_stack+11583, dvrr_stack+101133);
 tmp = dvrr_stack + 262581;
 target_ptr = Libderiv->deriv2_classes[3][4][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+262731, dvrr_stack+8088, dvrr_stack+144669);
 tmp = dvrr_stack + 262731;
 target_ptr = Libderiv->deriv2_classes[3][5][18];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+262941, dvrr_stack+122337, dvrr_stack+144795);
 tmp = dvrr_stack + 262941;
 target_ptr = Libderiv->deriv2_classes[3][6][18];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+263221, dvrr_stack+30044, dvrr_stack+31787);
 tmp = dvrr_stack + 263221;
 target_ptr = Libderiv->deriv2_classes[4][3][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+263371, dvrr_stack+29729, dvrr_stack+31637);
 tmp = dvrr_stack + 263371;
 target_ptr = Libderiv->deriv2_classes[4][4][18];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+263596, dvrr_stack+122757, dvrr_stack+25634);
 tmp = dvrr_stack + 263596;
 target_ptr = Libderiv->deriv2_classes[4][5][18];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+263911, dvrr_stack+123198, dvrr_stack+12621);
 tmp = dvrr_stack + 263911;
 target_ptr = Libderiv->deriv2_classes[4][6][18];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+264331, dvrr_stack+56587, dvrr_stack+17912);
 tmp = dvrr_stack + 264331;
 target_ptr = Libderiv->deriv2_classes[5][3][18];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+264541, dvrr_stack+56167, dvrr_stack+11583);
 tmp = dvrr_stack + 264541;
 target_ptr = Libderiv->deriv2_classes[5][4][18];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+264856, dvrr_stack+124542, dvrr_stack+8088);
 tmp = dvrr_stack + 264856;
 target_ptr = Libderiv->deriv2_classes[5][5][18];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+265297, dvrr_stack+125130, dvrr_stack+122337);
 tmp = dvrr_stack + 265297;
 target_ptr = Libderiv->deriv2_classes[5][6][18];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_i(Data,10,dvrr_stack+265885, dvrr_stack+138030, dvrr_stack+30044);
 tmp = dvrr_stack + 265885;
 target_ptr = Libderiv->deriv2_classes[6][3][18];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+266165, dvrr_stack+81964, dvrr_stack+29729);
 tmp = dvrr_stack + 266165;
 target_ptr = Libderiv->deriv2_classes[6][4][18];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+266585, dvrr_stack+48007, dvrr_stack+122757);
 tmp = dvrr_stack + 266585;
 target_ptr = Libderiv->deriv2_classes[6][5][18];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+267173, dvrr_stack+144963, dvrr_stack+123198);
 tmp = dvrr_stack + 267173;
 target_ptr = Libderiv->deriv2_classes[6][6][18];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+25053, dvrr_stack+139458, dvrr_stack+139398);
 tmp = dvrr_stack + 25053;
 target_ptr = Libderiv->deriv2_classes[3][3][14];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+78452, dvrr_stack+139698, dvrr_stack+139608);
 tmp = dvrr_stack + 78452;
 target_ptr = Libderiv->deriv2_classes[3][4][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+267957, dvrr_stack+140049, dvrr_stack+139923);
 tmp = dvrr_stack + 267957;
 target_ptr = Libderiv->deriv2_classes[3][5][14];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+268167, dvrr_stack+145971, dvrr_stack+140364);
 tmp = dvrr_stack + 268167;
 target_ptr = Libderiv->deriv2_classes[3][6][14];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+268447, dvrr_stack+141108, dvrr_stack+35828);
 tmp = dvrr_stack + 268447;
 target_ptr = Libderiv->deriv2_classes[4][3][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+268597, dvrr_stack+146391, dvrr_stack+140532);
 tmp = dvrr_stack + 268597;
 target_ptr = Libderiv->deriv2_classes[4][4][14];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+268822, dvrr_stack+146916, dvrr_stack+146706);
 tmp = dvrr_stack + 268822;
 target_ptr = Libderiv->deriv2_classes[4][5][14];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+269137, dvrr_stack+147637, dvrr_stack+147357);
 tmp = dvrr_stack + 269137;
 target_ptr = Libderiv->deriv2_classes[4][6][14];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+269557, dvrr_stack+148225, dvrr_stack+139458);
 tmp = dvrr_stack + 269557;
 target_ptr = Libderiv->deriv2_classes[5][3][14];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+269767, dvrr_stack+148505, dvrr_stack+139698);
 tmp = dvrr_stack + 269767;
 target_ptr = Libderiv->deriv2_classes[5][4][14];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+270082, dvrr_stack+148925, dvrr_stack+140049);
 tmp = dvrr_stack + 270082;
 target_ptr = Libderiv->deriv2_classes[5][5][14];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+270523, dvrr_stack+149513, dvrr_stack+145971);
 tmp = dvrr_stack + 270523;
 target_ptr = Libderiv->deriv2_classes[5][6][14];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,10,dvrr_stack+271111, dvrr_stack+150297, dvrr_stack+141108);
 tmp = dvrr_stack + 271111;
 target_ptr = Libderiv->deriv2_classes[6][3][14];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+271391, dvrr_stack+151617, dvrr_stack+146391);
 tmp = dvrr_stack + 271391;
 target_ptr = Libderiv->deriv2_classes[6][4][14];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+271811, dvrr_stack+150657, dvrr_stack+146916);
 tmp = dvrr_stack + 271811;
 target_ptr = Libderiv->deriv2_classes[6][5][14];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+272399, dvrr_stack+129322, dvrr_stack+147637);
 tmp = dvrr_stack + 272399;
 target_ptr = Libderiv->deriv2_classes[6][6][14];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+273183, dvrr_stack+35678, dvrr_stack+37395);
 tmp = dvrr_stack + 273183;
 target_ptr = Libderiv->deriv2_classes[3][3][13];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+273283, dvrr_stack+1512, dvrr_stack+479);
 tmp = dvrr_stack + 273283;
 target_ptr = Libderiv->deriv2_classes[3][4][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+273433, dvrr_stack+154614, dvrr_stack+130330);
 tmp = dvrr_stack + 273433;
 target_ptr = Libderiv->deriv2_classes[3][5][13];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+273643, dvrr_stack+127282, dvrr_stack+3423);
 tmp = dvrr_stack + 273643;
 target_ptr = Libderiv->deriv2_classes[3][6][13];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+55968, dvrr_stack+64504, dvrr_stack+130456);
 tmp = dvrr_stack + 55968;
 target_ptr = Libderiv->deriv2_classes[4][3][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+273923, dvrr_stack+152307, dvrr_stack+152157);
 tmp = dvrr_stack + 273923;
 target_ptr = Libderiv->deriv2_classes[4][4][13];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+274148, dvrr_stack+152832, dvrr_stack+152622);
 tmp = dvrr_stack + 274148;
 target_ptr = Libderiv->deriv2_classes[4][5][13];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+274463, dvrr_stack+72558, dvrr_stack+72278);
 tmp = dvrr_stack + 274463;
 target_ptr = Libderiv->deriv2_classes[4][6][13];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+274883, dvrr_stack+73146, dvrr_stack+35678);
 tmp = dvrr_stack + 274883;
 target_ptr = Libderiv->deriv2_classes[5][3][13];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+275093, dvrr_stack+10458, dvrr_stack+1512);
 tmp = dvrr_stack + 275093;
 target_ptr = Libderiv->deriv2_classes[5][4][13];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+275408, dvrr_stack+10878, dvrr_stack+154614);
 tmp = dvrr_stack + 275408;
 target_ptr = Libderiv->deriv2_classes[5][5][13];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+275849, dvrr_stack+132934, dvrr_stack+127282);
 tmp = dvrr_stack + 275849;
 target_ptr = Libderiv->deriv2_classes[5][6][13];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,10,dvrr_stack+276437, dvrr_stack+6798, dvrr_stack+64504);
 tmp = dvrr_stack + 276437;
 target_ptr = Libderiv->deriv2_classes[6][3][13];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+276717, dvrr_stack+133718, dvrr_stack+152307);
 tmp = dvrr_stack + 276717;
 target_ptr = Libderiv->deriv2_classes[6][4][13];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+277137, dvrr_stack+134258, dvrr_stack+152832);
 tmp = dvrr_stack + 277137;
 target_ptr = Libderiv->deriv2_classes[6][5][13];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+277725, dvrr_stack+154929, dvrr_stack+72558);
 tmp = dvrr_stack + 277725;
 target_ptr = Libderiv->deriv2_classes[6][6][13];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,10,dvrr_stack+278509, dvrr_stack+98088, dvrr_stack+96628);
 tmp = dvrr_stack + 278509;
 target_ptr = Libderiv->deriv2_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,15,dvrr_stack+278609, dvrr_stack+18832, dvrr_stack+11808);
 tmp = dvrr_stack + 278609;
 target_ptr = Libderiv->deriv2_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,21,dvrr_stack+11808, dvrr_stack+40440, dvrr_stack+16633);
 tmp = dvrr_stack + 11808;
 target_ptr = Libderiv->deriv2_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,28,dvrr_stack+278759, dvrr_stack+71018, dvrr_stack+4899);
 tmp = dvrr_stack + 278759;
 target_ptr = Libderiv->deriv2_classes[3][6][11];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,10,dvrr_stack+4899, dvrr_stack+99768, dvrr_stack+96928);
 tmp = dvrr_stack + 4899;
 target_ptr = Libderiv->deriv2_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,15,dvrr_stack+279039, dvrr_stack+99453, dvrr_stack+2973);
 tmp = dvrr_stack + 279039;
 target_ptr = Libderiv->deriv2_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,21,dvrr_stack+279264, dvrr_stack+99978, dvrr_stack+97028);
 tmp = dvrr_stack + 279264;
 target_ptr = Libderiv->deriv2_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,28,dvrr_stack+279579, dvrr_stack+100419, dvrr_stack+97238);
 tmp = dvrr_stack + 279579;
 target_ptr = Libderiv->deriv2_classes[4][6][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_h(Data,10,dvrr_stack+96838, dvrr_stack+87934, dvrr_stack+98088);
 tmp = dvrr_stack + 96838;
 target_ptr = Libderiv->deriv2_classes[5][3][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_h(Data,15,dvrr_stack+87907, dvrr_stack+9408, dvrr_stack+18832);
 tmp = dvrr_stack + 87907;
 target_ptr = Libderiv->deriv2_classes[5][4][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_h(Data,21,dvrr_stack+9363, dvrr_stack+60052, dvrr_stack+40440);
 tmp = dvrr_stack + 9363;
 target_ptr = Libderiv->deriv2_classes[5][5][11];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_h(Data,28,dvrr_stack+59940, dvrr_stack+60640, dvrr_stack+71018);
 tmp = dvrr_stack + 59940;
 target_ptr = Libderiv->deriv2_classes[5][6][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_i(Data,10,dvrr_stack+40351, dvrr_stack+36138, dvrr_stack+99768);
 tmp = dvrr_stack + 40351;
 target_ptr = Libderiv->deriv2_classes[6][3][11];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_i(Data,15,dvrr_stack+97048, dvrr_stack+37455, dvrr_stack+99453);
 tmp = dvrr_stack + 97048;
 target_ptr = Libderiv->deriv2_classes[6][4][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_i(Data,21,dvrr_stack+60528, dvrr_stack+41385, dvrr_stack+99978);
 tmp = dvrr_stack + 60528;
 target_ptr = Libderiv->deriv2_classes[6][5][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_i(Data,28,dvrr_stack+41346, dvrr_stack+138390, dvrr_stack+100419);
 tmp = dvrr_stack + 41346;
 target_ptr = Libderiv->deriv2_classes[6][6][11];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+138390, dvrr_stack+45487, dvrr_stack+79894);
 tmp = dvrr_stack + 138390;
 target_ptr = Libderiv->deriv2_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+138490, dvrr_stack+29378, dvrr_stack+6708);
 tmp = dvrr_stack + 138490;
 target_ptr = Libderiv->deriv2_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+138640, dvrr_stack+45637, dvrr_stack+79954);
 tmp = dvrr_stack + 138640;
 target_ptr = Libderiv->deriv2_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+138850, dvrr_stack+31898, dvrr_stack+80080);
 tmp = dvrr_stack + 138850;
 target_ptr = Libderiv->deriv2_classes[3][6][10];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+139130, dvrr_stack+33299, dvrr_stack+62432);
 tmp = dvrr_stack + 139130;
 target_ptr = Libderiv->deriv2_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+37455, dvrr_stack+32984, dvrr_stack+88274);
 tmp = dvrr_stack + 37455;
 target_ptr = Libderiv->deriv2_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+37680, dvrr_stack+33509, dvrr_stack+52399);
 tmp = dvrr_stack + 37680;
 target_ptr = Libderiv->deriv2_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+99414, dvrr_stack+13473, dvrr_stack+52609);
 tmp = dvrr_stack + 99414;
 target_ptr = Libderiv->deriv2_classes[4][6][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+36138, dvrr_stack+34118, dvrr_stack+45487);
 tmp = dvrr_stack + 36138;
 target_ptr = Libderiv->deriv2_classes[5][3][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+52383, dvrr_stack+14061, dvrr_stack+29378);
 tmp = dvrr_stack + 52383;
 target_ptr = Libderiv->deriv2_classes[5][4][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+99834, dvrr_stack+116821, dvrr_stack+45637);
 tmp = dvrr_stack + 99834;
 target_ptr = Libderiv->deriv2_classes[5][5][10];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+100275, dvrr_stack+117409, dvrr_stack+31898);
 tmp = dvrr_stack + 100275;
 target_ptr = Libderiv->deriv2_classes[5][6][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+14061, dvrr_stack+80248, dvrr_stack+33299);
 tmp = dvrr_stack + 14061;
 target_ptr = Libderiv->deriv2_classes[6][3][10];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+45407, dvrr_stack+80608, dvrr_stack+32984);
 tmp = dvrr_stack + 45407;
 target_ptr = Libderiv->deriv2_classes[6][4][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+32913, dvrr_stack+81148, dvrr_stack+33509);
 tmp = dvrr_stack + 32913;
 target_ptr = Libderiv->deriv2_classes[6][5][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+79832, dvrr_stack+46999, dvrr_stack+13473);
 tmp = dvrr_stack + 79832;
 target_ptr = Libderiv->deriv2_classes[6][6][10];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+14341, dvrr_stack+17482, dvrr_stack+81904);
 tmp = dvrr_stack + 14341;
 target_ptr = Libderiv->deriv2_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+36348, dvrr_stack+5699, dvrr_stack+4719);
 tmp = dvrr_stack + 36348;
 target_ptr = Libderiv->deriv2_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+34100, dvrr_stack+16318, dvrr_stack+2173);
 tmp = dvrr_stack + 34100;
 target_ptr = Libderiv->deriv2_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+33501, dvrr_stack+20459, dvrr_stack+22349);
 tmp = dvrr_stack + 33501;
 target_ptr = Libderiv->deriv2_classes[3][6][9];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+33781, dvrr_stack+38550, dvrr_stack+2073);
 tmp = dvrr_stack + 33781;
 target_ptr = Libderiv->deriv2_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+29343, dvrr_stack+98949, dvrr_stack+96688);
 tmp = dvrr_stack + 29343;
 target_ptr = Libderiv->deriv2_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+80616, dvrr_stack+36828, dvrr_stack+1023);
 tmp = dvrr_stack + 80616;
 target_ptr = Libderiv->deriv2_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+80931, dvrr_stack+42841, dvrr_stack+4439);
 tmp = dvrr_stack + 80931;
 target_ptr = Libderiv->deriv2_classes[4][6][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+96627, dvrr_stack+68498, dvrr_stack+17482);
 tmp = dvrr_stack + 96627;
 target_ptr = Libderiv->deriv2_classes[5][3][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+81351, dvrr_stack+87262, dvrr_stack+5699);
 tmp = dvrr_stack + 81351;
 target_ptr = Libderiv->deriv2_classes[5][4][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+116783, dvrr_stack+66104, dvrr_stack+16318);
 tmp = dvrr_stack + 116783;
 target_ptr = Libderiv->deriv2_classes[5][5][9];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+66040, dvrr_stack+74350, dvrr_stack+20459);
 tmp = dvrr_stack + 66040;
 target_ptr = Libderiv->deriv2_classes[5][6][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+81666, dvrr_stack+92167, dvrr_stack+38550);
 tmp = dvrr_stack + 81666;
 target_ptr = Libderiv->deriv2_classes[6][3][9];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+74294, dvrr_stack+128242, dvrr_stack+98949);
 tmp = dvrr_stack + 74294;
 target_ptr = Libderiv->deriv2_classes[6][4][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+117224, dvrr_stack+131422, dvrr_stack+36828);
 tmp = dvrr_stack + 117224;
 target_ptr = Libderiv->deriv2_classes[6][5][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+46960, dvrr_stack+136122, dvrr_stack+42841);
 tmp = dvrr_stack + 46960;
 target_ptr = Libderiv->deriv2_classes[6][6][9];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+42841, dvrr_stack+95007, dvrr_stack+140694);
 tmp = dvrr_stack + 42841;
 target_ptr = Libderiv->deriv2_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+42941, dvrr_stack+94782, dvrr_stack+4809);
 tmp = dvrr_stack + 42941;
 target_ptr = Libderiv->deriv2_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+43091, dvrr_stack+95157, dvrr_stack+140754);
 tmp = dvrr_stack + 43091;
 target_ptr = Libderiv->deriv2_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+136122, dvrr_stack+95472, dvrr_stack+140880);
 tmp = dvrr_stack + 136122;
 target_ptr = Libderiv->deriv2_classes[3][6][8];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+136402, dvrr_stack+82873, dvrr_stack+93742);
 tmp = dvrr_stack + 136402;
 target_ptr = Libderiv->deriv2_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+136552, dvrr_stack+82558, dvrr_stack+93592);
 tmp = dvrr_stack + 136552;
 target_ptr = Libderiv->deriv2_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+136777, dvrr_stack+83083, dvrr_stack+93842);
 tmp = dvrr_stack + 136777;
 target_ptr = Libderiv->deriv2_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+74714, dvrr_stack+83524, dvrr_stack+94052);
 tmp = dvrr_stack + 74714;
 target_ptr = Libderiv->deriv2_classes[4][6][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+128242, dvrr_stack+75722, dvrr_stack+95007);
 tmp = dvrr_stack + 128242;
 target_ptr = Libderiv->deriv2_classes[5][3][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+75722, dvrr_stack+75302, dvrr_stack+94782);
 tmp = dvrr_stack + 75722;
 target_ptr = Libderiv->deriv2_classes[5][4][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+94707, dvrr_stack+68778, dvrr_stack+95157);
 tmp = dvrr_stack + 94707;
 target_ptr = Libderiv->deriv2_classes[5][5][8];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+68472, dvrr_stack+57112, dvrr_stack+95472);
 tmp = dvrr_stack + 68472;
 target_ptr = Libderiv->deriv2_classes[5][6][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+128452, dvrr_stack+141369, dvrr_stack+82873);
 tmp = dvrr_stack + 128452;
 target_ptr = Libderiv->deriv2_classes[6][3][8];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+75284, dvrr_stack+128782, dvrr_stack+82558);
 tmp = dvrr_stack + 75284;
 target_ptr = Libderiv->deriv2_classes[6][4][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+128732, dvrr_stack+141729, dvrr_stack+83083);
 tmp = dvrr_stack + 128732;
 target_ptr = Libderiv->deriv2_classes[6][5][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+57077, dvrr_stack+142485, dvrr_stack+83524);
 tmp = dvrr_stack + 57077;
 target_ptr = Libderiv->deriv2_classes[6][6][8];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+43301, dvrr_stack+44014, dvrr_stack+141048);
 tmp = dvrr_stack + 43301;
 target_ptr = Libderiv->deriv2_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+95148, dvrr_stack+43789, dvrr_stack+98778);
 tmp = dvrr_stack + 95148;
 target_ptr = Libderiv->deriv2_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+95298, dvrr_stack+38760, dvrr_stack+101007);
 tmp = dvrr_stack + 95298;
 target_ptr = Libderiv->deriv2_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+95508, dvrr_stack+39075, dvrr_stack+143493);
 tmp = dvrr_stack + 95508;
 target_ptr = Libderiv->deriv2_classes[3][6][7];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+62374, dvrr_stack+31259, dvrr_stack+66842);
 tmp = dvrr_stack + 62374;
 target_ptr = Libderiv->deriv2_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+69060, dvrr_stack+30944, dvrr_stack+66692);
 tmp = dvrr_stack + 69060;
 target_ptr = Libderiv->deriv2_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+4409, dvrr_stack+25193, dvrr_stack+69366);
 tmp = dvrr_stack + 4409;
 target_ptr = Libderiv->deriv2_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+131400, dvrr_stack+12033, dvrr_stack+58904);
 tmp = dvrr_stack + 131400;
 target_ptr = Libderiv->deriv2_classes[4][6][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+58861, dvrr_stack+17632, dvrr_stack+44014);
 tmp = dvrr_stack + 58861;
 target_ptr = Libderiv->deriv2_classes[5][3][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+131820, dvrr_stack+20879, dvrr_stack+43789);
 tmp = dvrr_stack + 131820;
 target_ptr = Libderiv->deriv2_classes[5][4][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+16303, dvrr_stack+119957, dvrr_stack+38760);
 tmp = dvrr_stack + 16303;
 target_ptr = Libderiv->deriv2_classes[5][5][7];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+119908, dvrr_stack+120545, dvrr_stack+39075);
 tmp = dvrr_stack + 119908;
 target_ptr = Libderiv->deriv2_classes[5][6][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+43789, dvrr_stack+137130, dvrr_stack+31259);
 tmp = dvrr_stack + 43789;
 target_ptr = Libderiv->deriv2_classes[6][3][7];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+120496, dvrr_stack+137490, dvrr_stack+30944);
 tmp = dvrr_stack + 120496;
 target_ptr = Libderiv->deriv2_classes[6][4][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+137092, dvrr_stack+132178, dvrr_stack+25193);
 tmp = dvrr_stack + 137092;
 target_ptr = Libderiv->deriv2_classes[6][5][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+132135, dvrr_stack+143661, dvrr_stack+12033);
 tmp = dvrr_stack + 132135;
 target_ptr = Libderiv->deriv2_classes[6][6][7];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+30944, dvrr_stack+17912, dvrr_stack+98868);
 tmp = dvrr_stack + 30944;
 target_ptr = Libderiv->deriv2_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+31044, dvrr_stack+11583, dvrr_stack+101133);
 tmp = dvrr_stack + 31044;
 target_ptr = Libderiv->deriv2_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+31194, dvrr_stack+8088, dvrr_stack+144669);
 tmp = dvrr_stack + 31194;
 target_ptr = Libderiv->deriv2_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+98778, dvrr_stack+122337, dvrr_stack+144795);
 tmp = dvrr_stack + 98778;
 target_ptr = Libderiv->deriv2_classes[3][6][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+99058, dvrr_stack+30044, dvrr_stack+31787);
 tmp = dvrr_stack + 99058;
 target_ptr = Libderiv->deriv2_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+31787, dvrr_stack+29729, dvrr_stack+31637);
 tmp = dvrr_stack + 31787;
 target_ptr = Libderiv->deriv2_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+137680, dvrr_stack+122757, dvrr_stack+25634);
 tmp = dvrr_stack + 137680;
 target_ptr = Libderiv->deriv2_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+12018, dvrr_stack+123198, dvrr_stack+12621);
 tmp = dvrr_stack + 12018;
 target_ptr = Libderiv->deriv2_classes[4][6][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+32012, dvrr_stack+56587, dvrr_stack+17912);
 tmp = dvrr_stack + 32012;
 target_ptr = Libderiv->deriv2_classes[5][3][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+12438, dvrr_stack+56167, dvrr_stack+11583);
 tmp = dvrr_stack + 12438;
 target_ptr = Libderiv->deriv2_classes[5][4][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+20452, dvrr_stack+124542, dvrr_stack+8088);
 tmp = dvrr_stack + 20452;
 target_ptr = Libderiv->deriv2_classes[5][5][6];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+124486, dvrr_stack+125130, dvrr_stack+122337);
 tmp = dvrr_stack + 124486;
 target_ptr = Libderiv->deriv2_classes[5][6][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+8088, dvrr_stack+138030, dvrr_stack+30044);
 tmp = dvrr_stack + 8088;
 target_ptr = Libderiv->deriv2_classes[6][3][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+125074, dvrr_stack+81964, dvrr_stack+29729);
 tmp = dvrr_stack + 125074;
 target_ptr = Libderiv->deriv2_classes[6][4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+81946, dvrr_stack+48007, dvrr_stack+122757);
 tmp = dvrr_stack + 81946;
 target_ptr = Libderiv->deriv2_classes[6][5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+122334, dvrr_stack+144963, dvrr_stack+123198);
 tmp = dvrr_stack + 122334;
 target_ptr = Libderiv->deriv2_classes[6][6][6];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+137995, dvrr_stack+139458, dvrr_stack+139398);
 tmp = dvrr_stack + 137995;
 target_ptr = Libderiv->deriv2_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+138095, dvrr_stack+139698, dvrr_stack+139608);
 tmp = dvrr_stack + 138095;
 target_ptr = Libderiv->deriv2_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+125494, dvrr_stack+140049, dvrr_stack+139923);
 tmp = dvrr_stack + 125494;
 target_ptr = Libderiv->deriv2_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+29703, dvrr_stack+145971, dvrr_stack+140364);
 tmp = dvrr_stack + 29703;
 target_ptr = Libderiv->deriv2_classes[3][6][2];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+140364, dvrr_stack+141108, dvrr_stack+35828);
 tmp = dvrr_stack + 140364;
 target_ptr = Libderiv->deriv2_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+29983, dvrr_stack+146391, dvrr_stack+140532);
 tmp = dvrr_stack + 29983;
 target_ptr = Libderiv->deriv2_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+140514, dvrr_stack+146916, dvrr_stack+146706);
 tmp = dvrr_stack + 140514;
 target_ptr = Libderiv->deriv2_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+123118, dvrr_stack+147637, dvrr_stack+147357);
 tmp = dvrr_stack + 123118;
 target_ptr = Libderiv->deriv2_classes[4][6][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+146706, dvrr_stack+148225, dvrr_stack+139458);
 tmp = dvrr_stack + 146706;
 target_ptr = Libderiv->deriv2_classes[5][3][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+139280, dvrr_stack+148505, dvrr_stack+139698);
 tmp = dvrr_stack + 139280;
 target_ptr = Libderiv->deriv2_classes[5][4][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+148225, dvrr_stack+148925, dvrr_stack+140049);
 tmp = dvrr_stack + 148225;
 target_ptr = Libderiv->deriv2_classes[5][5][2];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+148666, dvrr_stack+149513, dvrr_stack+145971);
 tmp = dvrr_stack + 148666;
 target_ptr = Libderiv->deriv2_classes[5][6][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+147357, dvrr_stack+150297, dvrr_stack+141108);
 tmp = dvrr_stack + 147357;
 target_ptr = Libderiv->deriv2_classes[6][3][2];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+149254, dvrr_stack+151617, dvrr_stack+146391);
 tmp = dvrr_stack + 149254;
 target_ptr = Libderiv->deriv2_classes[6][4][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+149674, dvrr_stack+150657, dvrr_stack+146916);
 tmp = dvrr_stack + 149674;
 target_ptr = Libderiv->deriv2_classes[6][5][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+150262, dvrr_stack+129322, dvrr_stack+147637);
 tmp = dvrr_stack + 150262;
 target_ptr = Libderiv->deriv2_classes[6][6][2];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+35828, dvrr_stack+35678, dvrr_stack+37395);
 tmp = dvrr_stack + 35828;
 target_ptr = Libderiv->deriv2_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+147637, dvrr_stack+1512, dvrr_stack+479);
 tmp = dvrr_stack + 147637;
 target_ptr = Libderiv->deriv2_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+125704, dvrr_stack+154614, dvrr_stack+130330);
 tmp = dvrr_stack + 125704;
 target_ptr = Libderiv->deriv2_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+147787, dvrr_stack+127282, dvrr_stack+3423);
 tmp = dvrr_stack + 147787;
 target_ptr = Libderiv->deriv2_classes[3][6][1];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+3423, dvrr_stack+64504, dvrr_stack+130456);
 tmp = dvrr_stack + 3423;
 target_ptr = Libderiv->deriv2_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+146916, dvrr_stack+152307, dvrr_stack+152157);
 tmp = dvrr_stack + 146916;
 target_ptr = Libderiv->deriv2_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+151046, dvrr_stack+152832, dvrr_stack+152622);
 tmp = dvrr_stack + 151046;
 target_ptr = Libderiv->deriv2_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+140829, dvrr_stack+72558, dvrr_stack+72278);
 tmp = dvrr_stack + 140829;
 target_ptr = Libderiv->deriv2_classes[4][6][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+152622, dvrr_stack+73146, dvrr_stack+35678);
 tmp = dvrr_stack + 152622;
 target_ptr = Libderiv->deriv2_classes[5][3][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+141249, dvrr_stack+10458, dvrr_stack+1512);
 tmp = dvrr_stack + 141249;
 target_ptr = Libderiv->deriv2_classes[5][4][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+141564, dvrr_stack+10878, dvrr_stack+154614);
 tmp = dvrr_stack + 141564;
 target_ptr = Libderiv->deriv2_classes[5][5][1];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+142005, dvrr_stack+132934, dvrr_stack+127282);
 tmp = dvrr_stack + 142005;
 target_ptr = Libderiv->deriv2_classes[5][6][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+73146, dvrr_stack+6798, dvrr_stack+64504);
 tmp = dvrr_stack + 73146;
 target_ptr = Libderiv->deriv2_classes[6][3][1];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+127282, dvrr_stack+133718, dvrr_stack+152307);
 tmp = dvrr_stack + 127282;
 target_ptr = Libderiv->deriv2_classes[6][4][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+142593, dvrr_stack+134258, dvrr_stack+152832);
 tmp = dvrr_stack + 142593;
 target_ptr = Libderiv->deriv2_classes[6][5][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+143181, dvrr_stack+154929, dvrr_stack+72558);
 tmp = dvrr_stack + 143181;
 target_ptr = Libderiv->deriv2_classes[6][6][1];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+469, dvrr_stack+151413, dvrr_stack+130556);
 tmp = dvrr_stack + 469;
 target_ptr = Libderiv->deriv2_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+152832, dvrr_stack+15853, dvrr_stack+379);
 tmp = dvrr_stack + 152832;
 target_ptr = Libderiv->deriv2_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+152982, dvrr_stack+19507, dvrr_stack+873);
 tmp = dvrr_stack + 152982;
 target_ptr = Libderiv->deriv2_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+143965, dvrr_stack+92977, dvrr_stack+1863);
 tmp = dvrr_stack + 143965;
 target_ptr = Libderiv->deriv2_classes[3][6][0];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+148067, dvrr_stack+7158, dvrr_stack+135014);
 tmp = dvrr_stack + 148067;
 target_ptr = Libderiv->deriv2_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+64489, dvrr_stack+3591, dvrr_stack+4159);
 tmp = dvrr_stack + 64489;
 target_ptr = Libderiv->deriv2_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+144245, dvrr_stack+155937, dvrr_stack+569);
 tmp = dvrr_stack + 144245;
 target_ptr = Libderiv->deriv2_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+569, dvrr_stack+16759, dvrr_stack+156378);
 tmp = dvrr_stack + 569;
 target_ptr = Libderiv->deriv2_classes[4][6][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+156378, dvrr_stack+156658, dvrr_stack+151413);
 tmp = dvrr_stack + 156378;
 target_ptr = Libderiv->deriv2_classes[5][3][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+156588, dvrr_stack+19822, dvrr_stack+15853);
 tmp = dvrr_stack + 156588;
 target_ptr = Libderiv->deriv2_classes[5][4][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+151361, dvrr_stack+156938, dvrr_stack+19507);
 tmp = dvrr_stack + 151361;
 target_ptr = Libderiv->deriv2_classes[5][5][0];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+19507, dvrr_stack+130616, dvrr_stack+92977);
 tmp = dvrr_stack + 19507;
 target_ptr = Libderiv->deriv2_classes[5][6][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,10,dvrr_stack+156903, dvrr_stack+126922, dvrr_stack+7158);
 tmp = dvrr_stack + 156903;
 target_ptr = Libderiv->deriv2_classes[6][3][0];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+151802, dvrr_stack+127702, dvrr_stack+3591);
 tmp = dvrr_stack + 151802;
 target_ptr = Libderiv->deriv2_classes[6][4][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+144560, dvrr_stack+64714, dvrr_stack+155937);
 tmp = dvrr_stack + 144560;
 target_ptr = Libderiv->deriv2_classes[6][5][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+145148, dvrr_stack+135114, dvrr_stack+16759);
 tmp = dvrr_stack + 145148;
 target_ptr = Libderiv->deriv2_classes[6][6][0];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];


}

