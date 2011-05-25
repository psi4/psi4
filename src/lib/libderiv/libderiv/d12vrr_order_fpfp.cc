#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (fp|fp) integrals */

void d12vrr_order_fpfp(Libderiv_t *Libderiv, prim_data *Data)
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
 vrr_build_xxxx(am,Data,dvrr_stack+1449, dvrr_stack+1323, dvrr_stack+1527, dvrr_stack+479, dvrr_stack+633, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1477, dvrr_stack+1344, dvrr_stack+1323, dvrr_stack+494, dvrr_stack+479, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3423, dvrr_stack+1477, dvrr_stack+1449, NULL, NULL, dvrr_stack+1323);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3507, dvrr_stack+1428, dvrr_stack+1344, dvrr_stack+554, dvrr_stack+494, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3535, dvrr_stack+3507, dvrr_stack+1477, NULL, NULL, dvrr_stack+1344);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1344, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1347, dvrr_stack+121, dvrr_stack+1344, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1353, dvrr_stack+124, dvrr_stack+1347, dvrr_stack+614, dvrr_stack+121, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+494, dvrr_stack+170, dvrr_stack+1353, dvrr_stack+617, dvrr_stack+124, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1428, dvrr_stack+1512, dvrr_stack+494, dvrr_stack+623, dvrr_stack+170, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3619, dvrr_stack+1527, dvrr_stack+1428, dvrr_stack+633, dvrr_stack+1512, NULL);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3647, dvrr_stack+1449, dvrr_stack+3619, NULL, NULL, dvrr_stack+1527);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3731, dvrr_stack+3423, dvrr_stack+3647, dvrr_stack+1477, dvrr_stack+1449, dvrr_stack+1548);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3899, dvrr_stack+3535, dvrr_stack+3423, dvrr_stack+3507, dvrr_stack+1477, dvrr_stack+1365);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4067, dvrr_stack+3899, dvrr_stack+3731, dvrr_stack+3535, dvrr_stack+3423, dvrr_stack+1611);

 /* compute (3 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+4347,dvrr_stack+4067,dvrr_stack+1863,10);


 /* compute (3 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+4977,dvrr_stack+4347,dvrr_stack+2073,10);


 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+5877, dvrr_stack+4977, dvrr_stack+873);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+1477, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1480, dvrr_stack+0, dvrr_stack+30, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+1489, dvrr_stack+6, dvrr_stack+1480, dvrr_stack+3, dvrr_stack+0, dvrr_stack+1477);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+3507, dvrr_stack+33, dvrr_stack+213, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3525, dvrr_stack+39, dvrr_stack+3507, dvrr_stack+15, dvrr_stack+33, dvrr_stack+1480);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+6327, dvrr_stack+75, dvrr_stack+3525, dvrr_stack+57, dvrr_stack+39, dvrr_stack+1489);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+3561, dvrr_stack+219, dvrr_stack+623, NULL, NULL, dvrr_stack+213);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+6387, dvrr_stack+229, dvrr_stack+3561, dvrr_stack+111, dvrr_stack+219, dvrr_stack+3507);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+6447, dvrr_stack+259, dvrr_stack+6387, dvrr_stack+131, dvrr_stack+229, dvrr_stack+3525);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+6547, dvrr_stack+379, dvrr_stack+6447, dvrr_stack+319, dvrr_stack+259, dvrr_stack+6327);
 tmp = dvrr_stack + 6547;
 target_ptr = Libderiv->dvrr_classes[4][3];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6697, dvrr_stack+633, dvrr_stack+1512, NULL, NULL, dvrr_stack+623);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6742, dvrr_stack+648, dvrr_stack+6697, dvrr_stack+479, dvrr_stack+633, dvrr_stack+3561);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6832, dvrr_stack+693, dvrr_stack+6742, dvrr_stack+509, dvrr_stack+648, dvrr_stack+6387);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6982, dvrr_stack+873, dvrr_stack+6832, dvrr_stack+783, dvrr_stack+693, dvrr_stack+6447);

 /* compute (4 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+7207,dvrr_stack+6982,dvrr_stack+6547,15);


 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+7657, dvrr_stack+1527, dvrr_stack+1428, NULL, NULL, dvrr_stack+1512);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+7720, dvrr_stack+1548, dvrr_stack+7657, dvrr_stack+1323, dvrr_stack+1527, dvrr_stack+6697);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+7846, dvrr_stack+1611, dvrr_stack+7720, dvrr_stack+1365, dvrr_stack+1548, dvrr_stack+6742);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+8056, dvrr_stack+1863, dvrr_stack+7846, dvrr_stack+1737, dvrr_stack+1611, dvrr_stack+6832);

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+8371,dvrr_stack+8056,dvrr_stack+6982,15);


 /* compute (4 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+9046,dvrr_stack+8371,dvrr_stack+7207,15);


 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+9946, dvrr_stack+9046, dvrr_stack+6547);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1323, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+1344, dvrr_stack+1323, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+1347, dvrr_stack+15, dvrr_stack+121, dvrr_stack+1344, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+479, dvrr_stack+1353, dvrr_stack+111, dvrr_stack+124, dvrr_stack+1347, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1323, dvrr_stack+494, dvrr_stack+479, dvrr_stack+170, dvrr_stack+1353, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3591, dvrr_stack+1428, dvrr_stack+1323, dvrr_stack+1512, dvrr_stack+494, NULL);

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10396, dvrr_stack+3619, dvrr_stack+3591, NULL, NULL, dvrr_stack+1428);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10480, dvrr_stack+3647, dvrr_stack+10396, dvrr_stack+1449, dvrr_stack+3619, dvrr_stack+7657);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10648, dvrr_stack+3731, dvrr_stack+10480, dvrr_stack+3423, dvrr_stack+3647, dvrr_stack+7720);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10928, dvrr_stack+4067, dvrr_stack+10648, dvrr_stack+3899, dvrr_stack+3731, dvrr_stack+7846);

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+11348,dvrr_stack+10928,dvrr_stack+8056,15);


 /* compute (4 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+12293,dvrr_stack+11348,dvrr_stack+8371,15);


 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,225,dvrr_stack+13643, dvrr_stack+12293, dvrr_stack+6982);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,100,dvrr_stack+10396, dvrr_stack+2523, dvrr_stack+379);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+3591, dvrr_stack+4977, dvrr_stack+873);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+14318, dvrr_stack+9046, dvrr_stack+6547);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,225,dvrr_stack+14768, dvrr_stack+12293, dvrr_stack+6982);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,100,dvrr_stack+15443, dvrr_stack+2523, dvrr_stack+379);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+2523, dvrr_stack+4977, dvrr_stack+873);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+4977, dvrr_stack+9046, dvrr_stack+6547);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,225,dvrr_stack+9046, dvrr_stack+12293, dvrr_stack+6982);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+12293, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1344, dvrr_stack+21, dvrr_stack+3, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+12296, dvrr_stack+1344, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+12293);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+12314, dvrr_stack+164, dvrr_stack+24, NULL, NULL, dvrr_stack+21);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+12332, dvrr_stack+12314, dvrr_stack+57, dvrr_stack+164, dvrr_stack+24, dvrr_stack+1344);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+509, dvrr_stack+12332, dvrr_stack+75, dvrr_stack+12314, dvrr_stack+57, dvrr_stack+12296);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+12368,dvrr_stack+379,dvrr_stack+509,10);


 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+12548, dvrr_stack+12368, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+2973, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 2973;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+12608, dvrr_stack+1023, NULL);
 tmp = dvrr_stack + 12608;
 target_ptr = Libderiv->deriv_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+12708, dvrr_stack+4347, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+24, dvrr_stack+12293, dvrr_stack+1477, Data->F+2, Data->F+3, NULL);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+12918, dvrr_stack+12296, dvrr_stack+1489, dvrr_stack+1344, dvrr_stack+6, dvrr_stack+24);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+12948, dvrr_stack+509, dvrr_stack+6327, dvrr_stack+12332, dvrr_stack+75, dvrr_stack+12918);

 /* compute (4 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+13038,dvrr_stack+6547,dvrr_stack+12948,15);


 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+13308, dvrr_stack+13038, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+9721, dvrr_stack+8371, NULL);
 tmp = dvrr_stack + 9721;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+13398, dvrr_stack+7207, NULL);
 tmp = dvrr_stack + 13398;
 target_ptr = Libderiv->deriv_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+5427, dvrr_stack+11348, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+13548, dvrr_stack+12368, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+10696, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 10696;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+5742, dvrr_stack+1023, NULL);
 tmp = dvrr_stack + 5742;
 target_ptr = Libderiv->deriv_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+15743, dvrr_stack+4347, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+15953, dvrr_stack+13038, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+16043, dvrr_stack+8371, NULL);
 tmp = dvrr_stack + 16043;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+16268, dvrr_stack+7207, NULL);
 tmp = dvrr_stack + 16268;
 target_ptr = Libderiv->deriv_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+16418, dvrr_stack+11348, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+3423, dvrr_stack+12368, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+12368, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 12368;
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
 deriv_build_DX_0(Data,210,dvrr_stack+1023, dvrr_stack+4347, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+1233, dvrr_stack+13038, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+13038, dvrr_stack+8371, NULL);
 tmp = dvrr_stack + 13038;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+8371, dvrr_stack+7207, NULL);
 tmp = dvrr_stack + 8371;
 target_ptr = Libderiv->deriv_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+7207, dvrr_stack+11348, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+164, dvrr_stack+3, dvrr_stack+12293, Data->F+1, Data->F+2, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+11348, dvrr_stack+161, dvrr_stack+21, NULL, NULL, Data->F+1);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+57, dvrr_stack+11348, dvrr_stack+1344, dvrr_stack+161, dvrr_stack+21, dvrr_stack+3);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+12518, dvrr_stack+57, dvrr_stack+12296, dvrr_stack+11348, dvrr_stack+1344, dvrr_stack+164);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+11348, dvrr_stack+379, dvrr_stack+12518);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+11408, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 11408;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+11558, dvrr_stack+873, dvrr_stack+509);
 tmp = dvrr_stack + 11558;
 target_ptr = Libderiv->deriv_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+11658, dvrr_stack+4067, dvrr_stack+873);

 /* compute (3 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+11868, dvrr_stack+164, dvrr_stack+24, dvrr_stack+3, dvrr_stack+12293, NULL);

 /* compute (4 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+13263, dvrr_stack+12518, dvrr_stack+12918, dvrr_stack+57, dvrr_stack+12296, dvrr_stack+11868);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,15,1,dvrr_stack+11868, dvrr_stack+6547, dvrr_stack+13263);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+11958, dvrr_stack+8056, dvrr_stack+6547);
 tmp = dvrr_stack + 11958;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+8521, dvrr_stack+6982, dvrr_stack+12948);
 tmp = dvrr_stack + 8521;
 target_ptr = Libderiv->deriv_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+8671, dvrr_stack+10928, dvrr_stack+6982);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+8986, dvrr_stack+379, dvrr_stack+12518);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+4347, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 4347;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+12183, dvrr_stack+873, dvrr_stack+509);
 tmp = dvrr_stack + 12183;
 target_ptr = Libderiv->deriv_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+4497, dvrr_stack+4067, dvrr_stack+873);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,15,1,dvrr_stack+7522, dvrr_stack+6547, dvrr_stack+13263);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+4707, dvrr_stack+8056, dvrr_stack+6547);
 tmp = dvrr_stack + 4707;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+2173, dvrr_stack+6982, dvrr_stack+12948);
 tmp = dvrr_stack + 2173;
 target_ptr = Libderiv->deriv_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+16733, dvrr_stack+10928, dvrr_stack+6982);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+2323, dvrr_stack+379, dvrr_stack+12518);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+17048, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 17048;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+2383, dvrr_stack+873, dvrr_stack+509);
 tmp = dvrr_stack + 2383;
 target_ptr = Libderiv->deriv_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+17198, dvrr_stack+4067, dvrr_stack+873);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,15,1,dvrr_stack+4041, dvrr_stack+6547, dvrr_stack+13263);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+17408, dvrr_stack+8056, dvrr_stack+6547);
 tmp = dvrr_stack + 17408;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+4131, dvrr_stack+6982, dvrr_stack+12948);
 tmp = dvrr_stack + 4131;
 target_ptr = Libderiv->deriv_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+17633, dvrr_stack+10928, dvrr_stack+6982);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+10846,dvrr_stack+783,dvrr_stack+319,6);


 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+11026, dvrr_stack+10846, NULL);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+17948,dvrr_stack+1737,dvrr_stack+783,6);


 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+11086, dvrr_stack+17948, NULL);

 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+13263, dvrr_stack+1477, dvrr_stack+3, Data->F+3, Data->F+4, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1344, dvrr_stack+30, dvrr_stack+210, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+57, dvrr_stack+1480, dvrr_stack+1344, dvrr_stack+0, dvrr_stack+30, dvrr_stack+3);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+12518, dvrr_stack+1489, dvrr_stack+57, dvrr_stack+6, dvrr_stack+1480, dvrr_stack+13263);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+6, dvrr_stack+213, dvrr_stack+617, NULL, NULL, dvrr_stack+210);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+13269, dvrr_stack+3507, dvrr_stack+6, dvrr_stack+33, dvrr_stack+213, dvrr_stack+1344);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+11176, dvrr_stack+3525, dvrr_stack+13269, dvrr_stack+39, dvrr_stack+3507, dvrr_stack+57);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+11236, dvrr_stack+6327, dvrr_stack+11176, dvrr_stack+75, dvrr_stack+3525, dvrr_stack+12518);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+75, dvrr_stack+623, dvrr_stack+170, NULL, NULL, dvrr_stack+617);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4281, dvrr_stack+3561, dvrr_stack+75, dvrr_stack+219, dvrr_stack+623, dvrr_stack+6);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+18218, dvrr_stack+6387, dvrr_stack+4281, dvrr_stack+229, dvrr_stack+3561, dvrr_stack+13269);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+18318, dvrr_stack+6447, dvrr_stack+18218, dvrr_stack+259, dvrr_stack+6387, dvrr_stack+11176);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+18468, dvrr_stack+6547, dvrr_stack+18318, dvrr_stack+379, dvrr_stack+6447, dvrr_stack+11236);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+7612, dvrr_stack+1512, dvrr_stack+494, NULL, NULL, dvrr_stack+170);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+219, dvrr_stack+6697, dvrr_stack+7612, dvrr_stack+633, dvrr_stack+1512, dvrr_stack+75);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+18678, dvrr_stack+6742, dvrr_stack+219, dvrr_stack+648, dvrr_stack+6697, dvrr_stack+4281);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+18828, dvrr_stack+6832, dvrr_stack+18678, dvrr_stack+693, dvrr_stack+6742, dvrr_stack+18218);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+19053, dvrr_stack+6982, dvrr_stack+18828, dvrr_stack+873, dvrr_stack+6832, dvrr_stack+18318);

 /* compute (5 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+19368,dvrr_stack+19053,dvrr_stack+18468,21);


 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+19998, dvrr_stack+19368, NULL);

 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+633, dvrr_stack+1428, dvrr_stack+1323, NULL, NULL, dvrr_stack+494);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+20208, dvrr_stack+7657, dvrr_stack+633, dvrr_stack+1527, dvrr_stack+1428, dvrr_stack+7612);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+20334, dvrr_stack+7720, dvrr_stack+20208, dvrr_stack+1548, dvrr_stack+7657, dvrr_stack+219);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+20544, dvrr_stack+7846, dvrr_stack+20334, dvrr_stack+1611, dvrr_stack+7720, dvrr_stack+18678);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+20859, dvrr_stack+8056, dvrr_stack+20544, dvrr_stack+1863, dvrr_stack+7846, dvrr_stack+18828);

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+21300,dvrr_stack+20859,dvrr_stack+19053,21);


 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+7657, dvrr_stack+21300, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+1863, dvrr_stack+10846, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+1923, dvrr_stack+17948, NULL);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+1527, dvrr_stack+19368, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+7972, dvrr_stack+21300, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+2013, dvrr_stack+10846, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+10846, dvrr_stack+17948, NULL);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+17948, dvrr_stack+19368, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+19368, dvrr_stack+21300, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+18158, dvrr_stack+783, dvrr_stack+12332);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+10936, dvrr_stack+1737, dvrr_stack+319);

 /* compute (3 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+309, dvrr_stack+24, dvrr_stack+13263, dvrr_stack+12293, dvrr_stack+1477, NULL);

 /* compute (4 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+4932, dvrr_stack+12918, dvrr_stack+12518, dvrr_stack+12296, dvrr_stack+1489, dvrr_stack+309);

 /* compute (5 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+21300, dvrr_stack+12948, dvrr_stack+11236, dvrr_stack+509, dvrr_stack+6327, dvrr_stack+4932);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+21426, dvrr_stack+19053, dvrr_stack+21300);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+19683, dvrr_stack+20859, dvrr_stack+18468);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+509, dvrr_stack+783, dvrr_stack+12332);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+12918, dvrr_stack+1737, dvrr_stack+319);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+21636, dvrr_stack+19053, dvrr_stack+21300);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+21846, dvrr_stack+20859, dvrr_stack+18468);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+22161, dvrr_stack+783, dvrr_stack+12332);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+22221, dvrr_stack+1737, dvrr_stack+319);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+22311, dvrr_stack+19053, dvrr_stack+21300);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+22521, dvrr_stack+20859, dvrr_stack+18468);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+1737, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+633, dvrr_stack+18468, dvrr_stack+379);
 tmp = dvrr_stack + 633;
 target_ptr = Libderiv->deriv_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+22836, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+20208, dvrr_stack+19053, dvrr_stack+873);
 tmp = dvrr_stack + 20208;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+20433, dvrr_stack+6547, dvrr_stack+319);
 tmp = dvrr_stack + 20433;
 target_ptr = Libderiv->deriv_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+24, dvrr_stack+3, dvrr_stack+0, Data->F+4, Data->F+5, NULL);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+309, dvrr_stack+13263, dvrr_stack+24, dvrr_stack+1477, dvrr_stack+3, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1797, dvrr_stack+210, dvrr_stack+614, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+1806, dvrr_stack+1344, dvrr_stack+1797, dvrr_stack+30, dvrr_stack+210, dvrr_stack+0);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+13008, dvrr_stack+57, dvrr_stack+1806, dvrr_stack+1480, dvrr_stack+1344, dvrr_stack+24);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+4932, dvrr_stack+12518, dvrr_stack+13008, dvrr_stack+1489, dvrr_stack+57, dvrr_stack+309);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+12518, dvrr_stack+617, dvrr_stack+124, NULL, NULL, dvrr_stack+614);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+24, dvrr_stack+6, dvrr_stack+12518, dvrr_stack+213, dvrr_stack+617, dvrr_stack+1797);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+22926, dvrr_stack+13269, dvrr_stack+24, dvrr_stack+3507, dvrr_stack+6, dvrr_stack+1806);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+20533, dvrr_stack+11176, dvrr_stack+22926, dvrr_stack+3525, dvrr_stack+13269, dvrr_stack+13008);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+20623, dvrr_stack+11236, dvrr_stack+20533, dvrr_stack+6327, dvrr_stack+11176, dvrr_stack+4932);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+13008, dvrr_stack+170, dvrr_stack+1353, NULL, NULL, dvrr_stack+124);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+6327, dvrr_stack+75, dvrr_stack+13008, dvrr_stack+623, dvrr_stack+170, dvrr_stack+12518);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+11176, dvrr_stack+4281, dvrr_stack+6327, dvrr_stack+3561, dvrr_stack+75, dvrr_stack+24);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+18218, dvrr_stack+11176, dvrr_stack+6387, dvrr_stack+4281, dvrr_stack+22926);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+20749, dvrr_stack+18318, dvrr_stack+0, dvrr_stack+6447, dvrr_stack+18218, dvrr_stack+20533);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+20959, dvrr_stack+18468, dvrr_stack+20749, dvrr_stack+6547, dvrr_stack+18318, dvrr_stack+20623);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+20533, dvrr_stack+20959, dvrr_stack+6547);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+18218, dvrr_stack+6982, dvrr_stack+783);
 tmp = dvrr_stack + 18218;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4932, dvrr_stack+494, dvrr_stack+479, NULL, NULL, dvrr_stack+1353);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+18368, dvrr_stack+7612, dvrr_stack+4932, dvrr_stack+1512, dvrr_stack+494, dvrr_stack+13008);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6387, dvrr_stack+219, dvrr_stack+18368, dvrr_stack+6697, dvrr_stack+7612, dvrr_stack+6327);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+22926, dvrr_stack+18678, dvrr_stack+6387, dvrr_stack+6742, dvrr_stack+219, dvrr_stack+11176);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+23151, dvrr_stack+18828, dvrr_stack+22926, dvrr_stack+6832, dvrr_stack+18678, dvrr_stack+0);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+23466, dvrr_stack+19053, dvrr_stack+23151, dvrr_stack+6982, dvrr_stack+18828, dvrr_stack+20749);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+18678, dvrr_stack+23466, dvrr_stack+6982);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+18993, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+0, dvrr_stack+18468, dvrr_stack+379);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+22926, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+23016, dvrr_stack+19053, dvrr_stack+873);
 tmp = dvrr_stack + 23016;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+18368, dvrr_stack+6547, dvrr_stack+319);
 tmp = dvrr_stack + 18368;
 target_ptr = Libderiv->deriv_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+23241, dvrr_stack+20959, dvrr_stack+6547);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+11176, dvrr_stack+6982, dvrr_stack+783);
 tmp = dvrr_stack + 11176;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+23886, dvrr_stack+23466, dvrr_stack+6982);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+6327, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+150, dvrr_stack+18468, dvrr_stack+379);
 tmp = dvrr_stack + 150;
 target_ptr = Libderiv->deriv_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+379, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+6697, dvrr_stack+19053, dvrr_stack+873);
 tmp = dvrr_stack + 6697;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+873, dvrr_stack+6547, dvrr_stack+319);
 tmp = dvrr_stack + 873;
 target_ptr = Libderiv->deriv_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+18468, dvrr_stack+20959, dvrr_stack+6547);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+19053, dvrr_stack+6982, dvrr_stack+783);
 tmp = dvrr_stack + 19053;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+20743, dvrr_stack+23466, dvrr_stack+6982);

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,100,dvrr_stack+19203, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 19203;
 target_ptr = Libderiv->deriv2_classes[3][3][143];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+23451, dvrr_stack+5877, NULL);
 tmp = dvrr_stack + 23451;
 target_ptr = Libderiv->deriv2_classes[3][4][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+23601, dvrr_stack+9946, NULL);
 tmp = dvrr_stack + 23601;
 target_ptr = Libderiv->deriv2_classes[4][3][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,225,dvrr_stack+6387, dvrr_stack+13643, NULL);
 tmp = dvrr_stack + 6387;
 target_ptr = Libderiv->deriv2_classes[4][4][143];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,100,dvrr_stack+23751, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 23751;
 target_ptr = Libderiv->deriv2_classes[3][3][131];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+6922, dvrr_stack+5877, NULL);
 tmp = dvrr_stack + 6922;
 target_ptr = Libderiv->deriv2_classes[3][4][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+21058, dvrr_stack+9946, NULL);
 tmp = dvrr_stack + 21058;
 target_ptr = Libderiv->deriv2_classes[4][3][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,225,dvrr_stack+24201, dvrr_stack+13643, NULL);
 tmp = dvrr_stack + 24201;
 target_ptr = Libderiv->deriv2_classes[4][4][131];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+7072, dvrr_stack+10396, NULL);
 tmp = dvrr_stack + 7072;
 target_ptr = Libderiv->deriv2_classes[3][3][130];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+21208, dvrr_stack+3591, NULL);
 tmp = dvrr_stack + 21208;
 target_ptr = Libderiv->deriv2_classes[3][4][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+24426, dvrr_stack+14318, NULL);
 tmp = dvrr_stack + 24426;
 target_ptr = Libderiv->deriv2_classes[4][3][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+24576, dvrr_stack+14768, NULL);
 tmp = dvrr_stack + 24576;
 target_ptr = Libderiv->deriv2_classes[4][4][130];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,100,dvrr_stack+24801, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 24801;
 target_ptr = Libderiv->deriv2_classes[3][3][119];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+3123, dvrr_stack+5877, NULL);
 tmp = dvrr_stack + 3123;
 target_ptr = Libderiv->deriv2_classes[3][4][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+3273, dvrr_stack+9946, NULL);
 tmp = dvrr_stack + 3273;
 target_ptr = Libderiv->deriv2_classes[4][3][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,225,dvrr_stack+9946, dvrr_stack+13643, NULL);
 tmp = dvrr_stack + 9946;
 target_ptr = Libderiv->deriv2_classes[4][4][119];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+24901, dvrr_stack+10396, NULL);
 tmp = dvrr_stack + 24901;
 target_ptr = Libderiv->deriv2_classes[3][3][118];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+10171, dvrr_stack+3591, NULL);
 tmp = dvrr_stack + 10171;
 target_ptr = Libderiv->deriv2_classes[3][4][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+10321, dvrr_stack+14318, NULL);
 tmp = dvrr_stack + 10321;
 target_ptr = Libderiv->deriv2_classes[4][3][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+10471, dvrr_stack+14768, NULL);
 tmp = dvrr_stack + 10471;
 target_ptr = Libderiv->deriv2_classes[4][4][118];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+25001, dvrr_stack+15443, NULL);
 tmp = dvrr_stack + 25001;
 target_ptr = Libderiv->deriv2_classes[3][3][117];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+25101, dvrr_stack+2523, NULL);
 tmp = dvrr_stack + 25101;
 target_ptr = Libderiv->deriv2_classes[3][4][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+25251, dvrr_stack+4977, NULL);
 tmp = dvrr_stack + 25251;
 target_ptr = Libderiv->deriv2_classes[4][3][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+4932, dvrr_stack+9046, NULL);
 tmp = dvrr_stack + 4932;
 target_ptr = Libderiv->deriv2_classes[4][4][117];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+25401, dvrr_stack+2973, dvrr_stack+12548);
 tmp = dvrr_stack + 25401;
 target_ptr = Libderiv->deriv2_classes[3][3][107];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+9046, dvrr_stack+12708, dvrr_stack+12608);
 tmp = dvrr_stack + 9046;
 target_ptr = Libderiv->deriv2_classes[3][4][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+9196, dvrr_stack+9721, dvrr_stack+13308);
 tmp = dvrr_stack + 9196;
 target_ptr = Libderiv->deriv2_classes[4][3][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+9346, dvrr_stack+5427, dvrr_stack+13398);
 tmp = dvrr_stack + 9346;
 target_ptr = Libderiv->deriv2_classes[4][4][107];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+9571, dvrr_stack+10696, dvrr_stack+13548);
 tmp = dvrr_stack + 9571;
 target_ptr = Libderiv->deriv2_classes[3][3][106];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+5157, dvrr_stack+15743, dvrr_stack+5742);
 tmp = dvrr_stack + 5157;
 target_ptr = Libderiv->deriv2_classes[3][4][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+25501, dvrr_stack+16043, dvrr_stack+15953);
 tmp = dvrr_stack + 25501;
 target_ptr = Libderiv->deriv2_classes[4][3][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+25651, dvrr_stack+16418, dvrr_stack+16268);
 tmp = dvrr_stack + 25651;
 target_ptr = Libderiv->deriv2_classes[4][4][106];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+5307, dvrr_stack+12368, dvrr_stack+3423);
 tmp = dvrr_stack + 5307;
 target_ptr = Libderiv->deriv2_classes[3][3][105];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+1323, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 1323;
 target_ptr = Libderiv->deriv2_classes[3][4][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+2483, dvrr_stack+13038, dvrr_stack+1233);
 tmp = dvrr_stack + 2483;
 target_ptr = Libderiv->deriv2_classes[4][3][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+2633, dvrr_stack+7207, dvrr_stack+8371);
 tmp = dvrr_stack + 2633;
 target_ptr = Libderiv->deriv2_classes[4][4][105];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+25876, dvrr_stack+11408, dvrr_stack+11348);
 tmp = dvrr_stack + 25876;
 target_ptr = Libderiv->deriv2_classes[3][3][104];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+13608, dvrr_stack+11658, dvrr_stack+11558);
 tmp = dvrr_stack + 13608;
 target_ptr = Libderiv->deriv2_classes[3][4][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+13758, dvrr_stack+11958, dvrr_stack+11868);
 tmp = dvrr_stack + 13758;
 target_ptr = Libderiv->deriv2_classes[4][3][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+13908, dvrr_stack+8671, dvrr_stack+8521);
 tmp = dvrr_stack + 13908;
 target_ptr = Libderiv->deriv2_classes[4][4][104];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+2858, dvrr_stack+2973, dvrr_stack+12548);
 tmp = dvrr_stack + 2858;
 target_ptr = Libderiv->deriv2_classes[3][3][95];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+14133, dvrr_stack+12708, dvrr_stack+12608);
 tmp = dvrr_stack + 14133;
 target_ptr = Libderiv->deriv2_classes[3][4][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+14283, dvrr_stack+9721, dvrr_stack+13308);
 tmp = dvrr_stack + 14283;
 target_ptr = Libderiv->deriv2_classes[4][3][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+14433, dvrr_stack+5427, dvrr_stack+13398);
 tmp = dvrr_stack + 14433;
 target_ptr = Libderiv->deriv2_classes[4][4][95];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+14658, dvrr_stack+10696, dvrr_stack+13548);
 tmp = dvrr_stack + 14658;
 target_ptr = Libderiv->deriv2_classes[3][3][94];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+14758, dvrr_stack+15743, dvrr_stack+5742);
 tmp = dvrr_stack + 14758;
 target_ptr = Libderiv->deriv2_classes[3][4][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+14908, dvrr_stack+16043, dvrr_stack+15953);
 tmp = dvrr_stack + 14908;
 target_ptr = Libderiv->deriv2_classes[4][3][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+15058, dvrr_stack+16418, dvrr_stack+16268);
 tmp = dvrr_stack + 15058;
 target_ptr = Libderiv->deriv2_classes[4][4][94];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+15283, dvrr_stack+12368, dvrr_stack+3423);
 tmp = dvrr_stack + 15283;
 target_ptr = Libderiv->deriv2_classes[3][3][93];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+15383, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 15383;
 target_ptr = Libderiv->deriv2_classes[3][4][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+15533, dvrr_stack+13038, dvrr_stack+1233);
 tmp = dvrr_stack + 15533;
 target_ptr = Libderiv->deriv2_classes[4][3][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+5842, dvrr_stack+7207, dvrr_stack+8371);
 tmp = dvrr_stack + 5842;
 target_ptr = Libderiv->deriv2_classes[4][4][93];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+6067, dvrr_stack+11408, dvrr_stack+11348);
 tmp = dvrr_stack + 6067;
 target_ptr = Libderiv->deriv2_classes[3][3][92];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+6167, dvrr_stack+11658, dvrr_stack+11558);
 tmp = dvrr_stack + 6167;
 target_ptr = Libderiv->deriv2_classes[3][4][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+3483, dvrr_stack+11958, dvrr_stack+11868);
 tmp = dvrr_stack + 3483;
 target_ptr = Libderiv->deriv2_classes[4][3][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+3633, dvrr_stack+8671, dvrr_stack+8521);
 tmp = dvrr_stack + 3633;
 target_ptr = Libderiv->deriv2_classes[4][4][92];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+3858, dvrr_stack+4347, dvrr_stack+8986);
 tmp = dvrr_stack + 3858;
 target_ptr = Libderiv->deriv2_classes[3][3][91];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+25976, dvrr_stack+4497, dvrr_stack+12183);
 tmp = dvrr_stack + 25976;
 target_ptr = Libderiv->deriv2_classes[3][4][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+26126, dvrr_stack+4707, dvrr_stack+7522);
 tmp = dvrr_stack + 26126;
 target_ptr = Libderiv->deriv2_classes[4][3][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+26276, dvrr_stack+16733, dvrr_stack+2173);
 tmp = dvrr_stack + 26276;
 target_ptr = Libderiv->deriv2_classes[4][4][91];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+26501, dvrr_stack+2973, dvrr_stack+12548);
 tmp = dvrr_stack + 26501;
 target_ptr = Libderiv->deriv2_classes[3][3][83];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+26601, dvrr_stack+12708, dvrr_stack+12608);
 tmp = dvrr_stack + 26601;
 target_ptr = Libderiv->deriv2_classes[3][4][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+12708, dvrr_stack+9721, dvrr_stack+13308);
 tmp = dvrr_stack + 12708;
 target_ptr = Libderiv->deriv2_classes[4][3][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+26751, dvrr_stack+5427, dvrr_stack+13398);
 tmp = dvrr_stack + 26751;
 target_ptr = Libderiv->deriv2_classes[4][4][83];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+5407, dvrr_stack+10696, dvrr_stack+13548);
 tmp = dvrr_stack + 5407;
 target_ptr = Libderiv->deriv2_classes[3][3][82];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+5507, dvrr_stack+15743, dvrr_stack+5742);
 tmp = dvrr_stack + 5507;
 target_ptr = Libderiv->deriv2_classes[3][4][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+15683, dvrr_stack+16043, dvrr_stack+15953);
 tmp = dvrr_stack + 15683;
 target_ptr = Libderiv->deriv2_classes[4][3][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+26976, dvrr_stack+16418, dvrr_stack+16268);
 tmp = dvrr_stack + 26976;
 target_ptr = Libderiv->deriv2_classes[4][4][82];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+16418, dvrr_stack+12368, dvrr_stack+3423);
 tmp = dvrr_stack + 16418;
 target_ptr = Libderiv->deriv2_classes[3][3][81];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+16518, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 16518;
 target_ptr = Libderiv->deriv2_classes[3][4][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+973, dvrr_stack+13038, dvrr_stack+1233);
 tmp = dvrr_stack + 973;
 target_ptr = Libderiv->deriv2_classes[4][3][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+27201, dvrr_stack+7207, dvrr_stack+8371);
 tmp = dvrr_stack + 27201;
 target_ptr = Libderiv->deriv2_classes[4][4][81];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+1123, dvrr_stack+11408, dvrr_stack+11348);
 tmp = dvrr_stack + 1123;
 target_ptr = Libderiv->deriv2_classes[3][3][80];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+7172, dvrr_stack+11658, dvrr_stack+11558);
 tmp = dvrr_stack + 7172;
 target_ptr = Libderiv->deriv2_classes[3][4][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+11658, dvrr_stack+11958, dvrr_stack+11868);
 tmp = dvrr_stack + 11658;
 target_ptr = Libderiv->deriv2_classes[4][3][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+27426, dvrr_stack+8671, dvrr_stack+8521);
 tmp = dvrr_stack + 27426;
 target_ptr = Libderiv->deriv2_classes[4][4][80];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+1223, dvrr_stack+4347, dvrr_stack+8986);
 tmp = dvrr_stack + 1223;
 target_ptr = Libderiv->deriv2_classes[3][3][79];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+11808, dvrr_stack+4497, dvrr_stack+12183);
 tmp = dvrr_stack + 11808;
 target_ptr = Libderiv->deriv2_classes[3][4][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+4497, dvrr_stack+4707, dvrr_stack+7522);
 tmp = dvrr_stack + 4497;
 target_ptr = Libderiv->deriv2_classes[4][3][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+8671, dvrr_stack+16733, dvrr_stack+2173);
 tmp = dvrr_stack + 8671;
 target_ptr = Libderiv->deriv2_classes[4][4][79];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+8896, dvrr_stack+17048, dvrr_stack+2323);
 tmp = dvrr_stack + 8896;
 target_ptr = Libderiv->deriv2_classes[3][3][78];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+16668, dvrr_stack+17198, dvrr_stack+2383);
 tmp = dvrr_stack + 16668;
 target_ptr = Libderiv->deriv2_classes[3][4][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+17198, dvrr_stack+17408, dvrr_stack+4041);
 tmp = dvrr_stack + 17198;
 target_ptr = Libderiv->deriv2_classes[4][3][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+16818, dvrr_stack+17633, dvrr_stack+4131);
 tmp = dvrr_stack + 16818;
 target_ptr = Libderiv->deriv2_classes[4][4][78];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,10,dvrr_stack+17633, dvrr_stack+13398, dvrr_stack+11026);
 tmp = dvrr_stack + 17633;
 target_ptr = Libderiv->deriv2_classes[3][3][35];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,15,dvrr_stack+17733, dvrr_stack+9721, dvrr_stack+11086);
 tmp = dvrr_stack + 17733;
 target_ptr = Libderiv->deriv2_classes[3][4][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,10,dvrr_stack+7322, dvrr_stack+19998, dvrr_stack+12608);
 tmp = dvrr_stack + 7322;
 target_ptr = Libderiv->deriv2_classes[4][3][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,15,dvrr_stack+27651, dvrr_stack+7657, dvrr_stack+2973);
 tmp = dvrr_stack + 27651;
 target_ptr = Libderiv->deriv2_classes[4][4][35];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+7472, dvrr_stack+16268, dvrr_stack+1863);
 tmp = dvrr_stack + 7472;
 target_ptr = Libderiv->deriv2_classes[3][3][34];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+15833, dvrr_stack+16043, dvrr_stack+1923);
 tmp = dvrr_stack + 15833;
 target_ptr = Libderiv->deriv2_classes[3][4][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+3958, dvrr_stack+1527, dvrr_stack+5742);
 tmp = dvrr_stack + 3958;
 target_ptr = Libderiv->deriv2_classes[4][3][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+27876, dvrr_stack+7972, dvrr_stack+10696);
 tmp = dvrr_stack + 27876;
 target_ptr = Libderiv->deriv2_classes[4][4][34];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+28101, dvrr_stack+8371, dvrr_stack+2013);
 tmp = dvrr_stack + 28101;
 target_ptr = Libderiv->deriv2_classes[3][3][33];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+28201, dvrr_stack+13038, dvrr_stack+10846);
 tmp = dvrr_stack + 28201;
 target_ptr = Libderiv->deriv2_classes[3][4][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+28351, dvrr_stack+17948, dvrr_stack+2073);
 tmp = dvrr_stack + 28351;
 target_ptr = Libderiv->deriv2_classes[4][3][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+28501, dvrr_stack+19368, dvrr_stack+12368);
 tmp = dvrr_stack + 28501;
 target_ptr = Libderiv->deriv2_classes[4][4][33];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+28726, dvrr_stack+8521, dvrr_stack+18158);
 tmp = dvrr_stack + 28726;
 target_ptr = Libderiv->deriv2_classes[3][3][32];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+28826, dvrr_stack+11958, dvrr_stack+10936);
 tmp = dvrr_stack + 28826;
 target_ptr = Libderiv->deriv2_classes[3][4][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+28976, dvrr_stack+21426, dvrr_stack+11558);
 tmp = dvrr_stack + 28976;
 target_ptr = Libderiv->deriv2_classes[4][3][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+29126, dvrr_stack+19683, dvrr_stack+11408);
 tmp = dvrr_stack + 29126;
 target_ptr = Libderiv->deriv2_classes[4][4][32];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+29351, dvrr_stack+2173, dvrr_stack+509);
 tmp = dvrr_stack + 29351;
 target_ptr = Libderiv->deriv2_classes[3][3][31];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+29451, dvrr_stack+4707, dvrr_stack+12918);
 tmp = dvrr_stack + 29451;
 target_ptr = Libderiv->deriv2_classes[3][4][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+29601, dvrr_stack+21636, dvrr_stack+12183);
 tmp = dvrr_stack + 29601;
 target_ptr = Libderiv->deriv2_classes[4][3][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+29751, dvrr_stack+21846, dvrr_stack+4347);
 tmp = dvrr_stack + 29751;
 target_ptr = Libderiv->deriv2_classes[4][4][31];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+13263, dvrr_stack+4131, dvrr_stack+22161);
 tmp = dvrr_stack + 13263;
 target_ptr = Libderiv->deriv2_classes[3][3][30];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+29976, dvrr_stack+17408, dvrr_stack+22221);
 tmp = dvrr_stack + 29976;
 target_ptr = Libderiv->deriv2_classes[3][4][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+30126, dvrr_stack+22311, dvrr_stack+2383);
 tmp = dvrr_stack + 30126;
 target_ptr = Libderiv->deriv2_classes[4][3][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+30276, dvrr_stack+22521, dvrr_stack+17048);
 tmp = dvrr_stack + 30276;
 target_ptr = Libderiv->deriv2_classes[4][4][30];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+30501, dvrr_stack+633, dvrr_stack+1737);
 tmp = dvrr_stack + 30501;
 target_ptr = Libderiv->deriv2_classes[3][3][26];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+30601, dvrr_stack+20208, dvrr_stack+22836);
 tmp = dvrr_stack + 30601;
 target_ptr = Libderiv->deriv2_classes[3][4][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+30751, dvrr_stack+20533, dvrr_stack+20433);
 tmp = dvrr_stack + 30751;
 target_ptr = Libderiv->deriv2_classes[4][3][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+30901, dvrr_stack+18678, dvrr_stack+18218);
 tmp = dvrr_stack + 30901;
 target_ptr = Libderiv->deriv2_classes[4][4][26];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,10,dvrr_stack+31126, dvrr_stack+13398, dvrr_stack+11026);
 tmp = dvrr_stack + 31126;
 target_ptr = Libderiv->deriv2_classes[3][3][23];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,15,dvrr_stack+31226, dvrr_stack+9721, dvrr_stack+11086);
 tmp = dvrr_stack + 31226;
 target_ptr = Libderiv->deriv2_classes[3][4][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,10,dvrr_stack+31376, dvrr_stack+19998, dvrr_stack+12608);
 tmp = dvrr_stack + 31376;
 target_ptr = Libderiv->deriv2_classes[4][3][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,15,dvrr_stack+31526, dvrr_stack+7657, dvrr_stack+2973);
 tmp = dvrr_stack + 31526;
 target_ptr = Libderiv->deriv2_classes[4][4][23];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+31751, dvrr_stack+16268, dvrr_stack+1863);
 tmp = dvrr_stack + 31751;
 target_ptr = Libderiv->deriv2_classes[3][3][22];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+31851, dvrr_stack+16043, dvrr_stack+1923);
 tmp = dvrr_stack + 31851;
 target_ptr = Libderiv->deriv2_classes[3][4][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+32001, dvrr_stack+1527, dvrr_stack+5742);
 tmp = dvrr_stack + 32001;
 target_ptr = Libderiv->deriv2_classes[4][3][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+32151, dvrr_stack+7972, dvrr_stack+10696);
 tmp = dvrr_stack + 32151;
 target_ptr = Libderiv->deriv2_classes[4][4][22];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+32376, dvrr_stack+8371, dvrr_stack+2013);
 tmp = dvrr_stack + 32376;
 target_ptr = Libderiv->deriv2_classes[3][3][21];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+32476, dvrr_stack+13038, dvrr_stack+10846);
 tmp = dvrr_stack + 32476;
 target_ptr = Libderiv->deriv2_classes[3][4][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+32626, dvrr_stack+17948, dvrr_stack+2073);
 tmp = dvrr_stack + 32626;
 target_ptr = Libderiv->deriv2_classes[4][3][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+32776, dvrr_stack+19368, dvrr_stack+12368);
 tmp = dvrr_stack + 32776;
 target_ptr = Libderiv->deriv2_classes[4][4][21];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+33001, dvrr_stack+8521, dvrr_stack+18158);
 tmp = dvrr_stack + 33001;
 target_ptr = Libderiv->deriv2_classes[3][3][20];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+33101, dvrr_stack+11958, dvrr_stack+10936);
 tmp = dvrr_stack + 33101;
 target_ptr = Libderiv->deriv2_classes[3][4][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+33251, dvrr_stack+21426, dvrr_stack+11558);
 tmp = dvrr_stack + 33251;
 target_ptr = Libderiv->deriv2_classes[4][3][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+33401, dvrr_stack+19683, dvrr_stack+11408);
 tmp = dvrr_stack + 33401;
 target_ptr = Libderiv->deriv2_classes[4][4][20];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+33626, dvrr_stack+2173, dvrr_stack+509);
 tmp = dvrr_stack + 33626;
 target_ptr = Libderiv->deriv2_classes[3][3][19];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+33726, dvrr_stack+4707, dvrr_stack+12918);
 tmp = dvrr_stack + 33726;
 target_ptr = Libderiv->deriv2_classes[3][4][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+33876, dvrr_stack+21636, dvrr_stack+12183);
 tmp = dvrr_stack + 33876;
 target_ptr = Libderiv->deriv2_classes[4][3][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+34026, dvrr_stack+21846, dvrr_stack+4347);
 tmp = dvrr_stack + 34026;
 target_ptr = Libderiv->deriv2_classes[4][4][19];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+34251, dvrr_stack+4131, dvrr_stack+22161);
 tmp = dvrr_stack + 34251;
 target_ptr = Libderiv->deriv2_classes[3][3][18];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+34351, dvrr_stack+17408, dvrr_stack+22221);
 tmp = dvrr_stack + 34351;
 target_ptr = Libderiv->deriv2_classes[3][4][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+34501, dvrr_stack+22311, dvrr_stack+2383);
 tmp = dvrr_stack + 34501;
 target_ptr = Libderiv->deriv2_classes[4][3][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+34651, dvrr_stack+22521, dvrr_stack+17048);
 tmp = dvrr_stack + 34651;
 target_ptr = Libderiv->deriv2_classes[4][4][18];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+34876, dvrr_stack+633, dvrr_stack+1737);
 tmp = dvrr_stack + 34876;
 target_ptr = Libderiv->deriv2_classes[3][3][14];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+34976, dvrr_stack+20208, dvrr_stack+22836);
 tmp = dvrr_stack + 34976;
 target_ptr = Libderiv->deriv2_classes[3][4][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+35126, dvrr_stack+20533, dvrr_stack+20433);
 tmp = dvrr_stack + 35126;
 target_ptr = Libderiv->deriv2_classes[4][3][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+35276, dvrr_stack+18678, dvrr_stack+18218);
 tmp = dvrr_stack + 35276;
 target_ptr = Libderiv->deriv2_classes[4][4][14];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+35501, dvrr_stack+0, dvrr_stack+18993);
 tmp = dvrr_stack + 35501;
 target_ptr = Libderiv->deriv2_classes[3][3][13];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+35601, dvrr_stack+23016, dvrr_stack+22926);
 tmp = dvrr_stack + 35601;
 target_ptr = Libderiv->deriv2_classes[3][4][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+35751, dvrr_stack+23241, dvrr_stack+18368);
 tmp = dvrr_stack + 35751;
 target_ptr = Libderiv->deriv2_classes[4][3][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+35901, dvrr_stack+23886, dvrr_stack+11176);
 tmp = dvrr_stack + 35901;
 target_ptr = Libderiv->deriv2_classes[4][4][13];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,10,dvrr_stack+36126, dvrr_stack+13398, dvrr_stack+11026);
 tmp = dvrr_stack + 36126;
 target_ptr = Libderiv->deriv2_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,15,dvrr_stack+36226, dvrr_stack+9721, dvrr_stack+11086);
 tmp = dvrr_stack + 36226;
 target_ptr = Libderiv->deriv2_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,10,dvrr_stack+11026, dvrr_stack+19998, dvrr_stack+12608);
 tmp = dvrr_stack + 11026;
 target_ptr = Libderiv->deriv2_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,15,dvrr_stack+9671, dvrr_stack+7657, dvrr_stack+2973);
 tmp = dvrr_stack + 9671;
 target_ptr = Libderiv->deriv2_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+19998, dvrr_stack+16268, dvrr_stack+1863);
 tmp = dvrr_stack + 19998;
 target_ptr = Libderiv->deriv2_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+16268, dvrr_stack+16043, dvrr_stack+1923);
 tmp = dvrr_stack + 16268;
 target_ptr = Libderiv->deriv2_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+7572, dvrr_stack+1527, dvrr_stack+5742);
 tmp = dvrr_stack + 7572;
 target_ptr = Libderiv->deriv2_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+7722, dvrr_stack+7972, dvrr_stack+10696);
 tmp = dvrr_stack + 7722;
 target_ptr = Libderiv->deriv2_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+10696, dvrr_stack+8371, dvrr_stack+2013);
 tmp = dvrr_stack + 10696;
 target_ptr = Libderiv->deriv2_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+7947, dvrr_stack+13038, dvrr_stack+10846);
 tmp = dvrr_stack + 7947;
 target_ptr = Libderiv->deriv2_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+8097, dvrr_stack+17948, dvrr_stack+2073);
 tmp = dvrr_stack + 8097;
 target_ptr = Libderiv->deriv2_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+17883, dvrr_stack+19368, dvrr_stack+12368);
 tmp = dvrr_stack + 17883;
 target_ptr = Libderiv->deriv2_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+10796, dvrr_stack+8521, dvrr_stack+18158);
 tmp = dvrr_stack + 10796;
 target_ptr = Libderiv->deriv2_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+19303, dvrr_stack+11958, dvrr_stack+10936);
 tmp = dvrr_stack + 19303;
 target_ptr = Libderiv->deriv2_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+11958, dvrr_stack+21426, dvrr_stack+11558);
 tmp = dvrr_stack + 11958;
 target_ptr = Libderiv->deriv2_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+19453, dvrr_stack+19683, dvrr_stack+11408);
 tmp = dvrr_stack + 19453;
 target_ptr = Libderiv->deriv2_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+10896, dvrr_stack+2173, dvrr_stack+509);
 tmp = dvrr_stack + 10896;
 target_ptr = Libderiv->deriv2_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+19678, dvrr_stack+4707, dvrr_stack+12918);
 tmp = dvrr_stack + 19678;
 target_ptr = Libderiv->deriv2_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+4647, dvrr_stack+21636, dvrr_stack+12183);
 tmp = dvrr_stack + 4647;
 target_ptr = Libderiv->deriv2_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+12108, dvrr_stack+21846, dvrr_stack+4347);
 tmp = dvrr_stack + 12108;
 target_ptr = Libderiv->deriv2_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+12333, dvrr_stack+4131, dvrr_stack+22161);
 tmp = dvrr_stack + 12333;
 target_ptr = Libderiv->deriv2_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+12433, dvrr_stack+17408, dvrr_stack+22221);
 tmp = dvrr_stack + 12433;
 target_ptr = Libderiv->deriv2_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+17348, dvrr_stack+22311, dvrr_stack+2383);
 tmp = dvrr_stack + 17348;
 target_ptr = Libderiv->deriv2_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+12858, dvrr_stack+22521, dvrr_stack+17048);
 tmp = dvrr_stack + 12858;
 target_ptr = Libderiv->deriv2_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+12583, dvrr_stack+633, dvrr_stack+1737);
 tmp = dvrr_stack + 12583;
 target_ptr = Libderiv->deriv2_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+17043, dvrr_stack+20208, dvrr_stack+22836);
 tmp = dvrr_stack + 17043;
 target_ptr = Libderiv->deriv2_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+20098, dvrr_stack+20533, dvrr_stack+20433);
 tmp = dvrr_stack + 20098;
 target_ptr = Libderiv->deriv2_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+20248, dvrr_stack+18678, dvrr_stack+18218);
 tmp = dvrr_stack + 20248;
 target_ptr = Libderiv->deriv2_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+18678, dvrr_stack+0, dvrr_stack+18993);
 tmp = dvrr_stack + 18678;
 target_ptr = Libderiv->deriv2_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+0, dvrr_stack+23016, dvrr_stack+22926);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+18778, dvrr_stack+23241, dvrr_stack+18368);
 tmp = dvrr_stack + 18778;
 target_ptr = Libderiv->deriv2_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+20473, dvrr_stack+23886, dvrr_stack+11176);
 tmp = dvrr_stack + 20473;
 target_ptr = Libderiv->deriv2_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+11176, dvrr_stack+150, dvrr_stack+6327);
 tmp = dvrr_stack + 11176;
 target_ptr = Libderiv->deriv2_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+150, dvrr_stack+6697, dvrr_stack+379);
 tmp = dvrr_stack + 150;
 target_ptr = Libderiv->deriv2_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+300, dvrr_stack+18468, dvrr_stack+873);
 tmp = dvrr_stack + 300;
 target_ptr = Libderiv->deriv2_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+450, dvrr_stack+20743, dvrr_stack+19053);
 tmp = dvrr_stack + 450;
 target_ptr = Libderiv->deriv2_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];


}

