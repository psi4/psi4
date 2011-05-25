#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (fp|fd) integrals */

void d12vrr_order_fpfd(Libderiv_t *Libderiv, prim_data *Data)
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
 vrr_build_xxxx(am,Data,dvrr_stack+3591, dvrr_stack+3423, dvrr_stack+3711, dvrr_stack+1323, dvrr_stack+1527, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3627, dvrr_stack+3451, dvrr_stack+3423, dvrr_stack+1344, dvrr_stack+1323, NULL);

 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6419, dvrr_stack+3627, dvrr_stack+3591, NULL, NULL, dvrr_stack+3423);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6527, dvrr_stack+3563, dvrr_stack+3451, dvrr_stack+1428, dvrr_stack+1344, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6563, dvrr_stack+6527, dvrr_stack+3627, NULL, NULL, dvrr_stack+3451);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3451, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+503, dvrr_stack+494, dvrr_stack+3451, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+3454, dvrr_stack+497, dvrr_stack+503, dvrr_stack+121, dvrr_stack+494, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3464, dvrr_stack+554, dvrr_stack+3454, dvrr_stack+124, dvrr_stack+497, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1344, dvrr_stack+3675, dvrr_stack+3464, dvrr_stack+170, dvrr_stack+554, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3563, dvrr_stack+3690, dvrr_stack+1344, dvrr_stack+1512, dvrr_stack+3675, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6671, dvrr_stack+3711, dvrr_stack+3563, dvrr_stack+1527, dvrr_stack+3690, NULL);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6707, dvrr_stack+3591, dvrr_stack+6671, NULL, NULL, dvrr_stack+3711);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6815, dvrr_stack+6419, dvrr_stack+6707, dvrr_stack+3627, dvrr_stack+3591, dvrr_stack+3739);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7031, dvrr_stack+6563, dvrr_stack+6419, dvrr_stack+6527, dvrr_stack+3627, dvrr_stack+3479);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7247, dvrr_stack+7031, dvrr_stack+6815, dvrr_stack+6563, dvrr_stack+6419, dvrr_stack+3823);

 /* compute (3 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+7607,dvrr_stack+7247,dvrr_stack+4159,10);


 /* compute (3 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+8447,dvrr_stack+7607,dvrr_stack+4439,10);


 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,210,dvrr_stack+9707, dvrr_stack+8447, dvrr_stack+1863);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+121, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3627, dvrr_stack+0, dvrr_stack+30, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+3636, dvrr_stack+6, dvrr_stack+3627, dvrr_stack+3, dvrr_stack+0, dvrr_stack+121);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+3654, dvrr_stack+33, dvrr_stack+213, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+6527, dvrr_stack+39, dvrr_stack+3654, dvrr_stack+15, dvrr_stack+33, dvrr_stack+3627);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+6563, dvrr_stack+75, dvrr_stack+6527, dvrr_stack+57, dvrr_stack+39, dvrr_stack+3636);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+6623, dvrr_stack+219, dvrr_stack+623, NULL, NULL, dvrr_stack+213);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+10337, dvrr_stack+229, dvrr_stack+6623, dvrr_stack+111, dvrr_stack+219, dvrr_stack+3654);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+10397, dvrr_stack+259, dvrr_stack+10337, dvrr_stack+131, dvrr_stack+229, dvrr_stack+6527);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10497, dvrr_stack+379, dvrr_stack+10397, dvrr_stack+319, dvrr_stack+259, dvrr_stack+6563);
 tmp = dvrr_stack + 10497;
 target_ptr = Libderiv->dvrr_classes[4][3];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10647, dvrr_stack+633, dvrr_stack+1512, NULL, NULL, dvrr_stack+623);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10692, dvrr_stack+648, dvrr_stack+10647, dvrr_stack+479, dvrr_stack+633, dvrr_stack+6623);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10782, dvrr_stack+693, dvrr_stack+10692, dvrr_stack+509, dvrr_stack+648, dvrr_stack+10337);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10932, dvrr_stack+873, dvrr_stack+10782, dvrr_stack+783, dvrr_stack+693, dvrr_stack+10397);
 tmp = dvrr_stack + 10932;
 target_ptr = Libderiv->dvrr_classes[4][4];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+11157,dvrr_stack+10932,dvrr_stack+10497,15);


 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+11607, dvrr_stack+1527, dvrr_stack+3690, NULL, NULL, dvrr_stack+1512);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+11670, dvrr_stack+1548, dvrr_stack+11607, dvrr_stack+1323, dvrr_stack+1527, dvrr_stack+10647);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+11796, dvrr_stack+1611, dvrr_stack+11670, dvrr_stack+1365, dvrr_stack+1548, dvrr_stack+10692);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+12006, dvrr_stack+1863, dvrr_stack+11796, dvrr_stack+1737, dvrr_stack+1611, dvrr_stack+10782);

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+12321,dvrr_stack+12006,dvrr_stack+10932,15);


 /* compute (4 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+12996,dvrr_stack+12321,dvrr_stack+11157,15);


 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+13896, dvrr_stack+12996, dvrr_stack+10497);

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1365, dvrr_stack+3711, dvrr_stack+3563, NULL, NULL, dvrr_stack+3690);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14346, dvrr_stack+3739, dvrr_stack+1365, dvrr_stack+3423, dvrr_stack+3711, dvrr_stack+11607);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14514, dvrr_stack+3823, dvrr_stack+14346, dvrr_stack+3479, dvrr_stack+3739, dvrr_stack+11670);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14794, dvrr_stack+4159, dvrr_stack+14514, dvrr_stack+3991, dvrr_stack+3823, dvrr_stack+11796);

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+15214,dvrr_stack+14794,dvrr_stack+12006,15);


 /* compute (4 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+16159,dvrr_stack+15214,dvrr_stack+12321,15);


 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,225,dvrr_stack+17509, dvrr_stack+16159, dvrr_stack+10932);

 /* compute (0 0 | 1 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3672, Data->F+10, Data->F+11, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3451, dvrr_stack+3672, Data->F+9, Data->F+10, NULL);

 /* compute (0 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+503, dvrr_stack+15, dvrr_stack+494, dvrr_stack+3451, NULL);

 /* compute (0 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3479, dvrr_stack+3454, dvrr_stack+111, dvrr_stack+497, dvrr_stack+503, NULL);

 /* compute (0 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1323, dvrr_stack+3464, dvrr_stack+3479, dvrr_stack+554, dvrr_stack+3454, NULL);

 /* compute (0 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3479, dvrr_stack+1344, dvrr_stack+1323, dvrr_stack+3675, dvrr_stack+3464, NULL);

 /* compute (0 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3507, dvrr_stack+3563, dvrr_stack+3479, dvrr_stack+3690, dvrr_stack+1344, NULL);

 /* compute (1 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18184, dvrr_stack+6671, dvrr_stack+3507, NULL, NULL, dvrr_stack+3563);

 /* compute (2 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18292, dvrr_stack+6707, dvrr_stack+18184, dvrr_stack+3591, dvrr_stack+6671, dvrr_stack+1365);

 /* compute (3 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18508, dvrr_stack+6815, dvrr_stack+18292, dvrr_stack+6419, dvrr_stack+6707, dvrr_stack+14346);

 /* compute (4 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18868, dvrr_stack+7247, dvrr_stack+18508, dvrr_stack+7031, dvrr_stack+6815, dvrr_stack+14514);

 /* compute (4 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+19408,dvrr_stack+18868,dvrr_stack+14794,15);


 /* compute (4 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+20668,dvrr_stack+19408,dvrr_stack+15214,15);


 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,315,dvrr_stack+22558, dvrr_stack+20668, dvrr_stack+12006);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,100,dvrr_stack+18184, dvrr_stack+2523, dvrr_stack+379);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+6653, dvrr_stack+5069, dvrr_stack+873);

 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,210,dvrr_stack+23503, dvrr_stack+8447, dvrr_stack+1863);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+24133, dvrr_stack+12996, dvrr_stack+10497);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,225,dvrr_stack+24583, dvrr_stack+16159, dvrr_stack+10932);

 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,315,dvrr_stack+25258, dvrr_stack+20668, dvrr_stack+12006);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,100,dvrr_stack+18484, dvrr_stack+2523, dvrr_stack+379);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+2523, dvrr_stack+5069, dvrr_stack+873);

 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,210,dvrr_stack+5069, dvrr_stack+8447, dvrr_stack+1863);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+8447, dvrr_stack+12996, dvrr_stack+10497);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,225,dvrr_stack+12996, dvrr_stack+16159, dvrr_stack+10932);

 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,315,dvrr_stack+16159, dvrr_stack+20668, dvrr_stack+12006);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3672, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+20668, dvrr_stack+21, dvrr_stack+3, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+20677, dvrr_stack+20668, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+3672);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+20695, dvrr_stack+164, dvrr_stack+24, NULL, NULL, dvrr_stack+21);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3591, dvrr_stack+20695, dvrr_stack+57, dvrr_stack+164, dvrr_stack+24, dvrr_stack+20668);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+20713, dvrr_stack+3591, dvrr_stack+75, dvrr_stack+20695, dvrr_stack+57, dvrr_stack+20677);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+20773,dvrr_stack+379,dvrr_stack+20713,10);


 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+20953, dvrr_stack+20773, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+2973, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 2973;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+21013, dvrr_stack+1023, NULL);
 tmp = dvrr_stack + 21013;
 target_ptr = Libderiv->deriv_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+21113, dvrr_stack+4439, NULL);
 tmp = dvrr_stack + 21113;
 target_ptr = Libderiv->deriv_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,280,dvrr_stack+21323, dvrr_stack+7607, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+24, dvrr_stack+3672, dvrr_stack+121, Data->F+2, Data->F+3, NULL);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+21603, dvrr_stack+20677, dvrr_stack+3636, dvrr_stack+20668, dvrr_stack+6, dvrr_stack+24);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+21633, dvrr_stack+20713, dvrr_stack+6563, dvrr_stack+3591, dvrr_stack+75, dvrr_stack+21603);

 /* compute (4 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+5699,dvrr_stack+10497,dvrr_stack+21633,15);


 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+21723, dvrr_stack+5699, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+13671, dvrr_stack+12321, NULL);
 tmp = dvrr_stack + 13671;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+21813, dvrr_stack+11157, NULL);
 tmp = dvrr_stack + 21813;
 target_ptr = Libderiv->deriv_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+21963, dvrr_stack+15214, NULL);
 tmp = dvrr_stack + 21963;
 target_ptr = Libderiv->deriv_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+8897, dvrr_stack+19408, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+22278, dvrr_stack+20773, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+22338, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 22338;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+17104, dvrr_stack+1023, NULL);
 tmp = dvrr_stack + 17104;
 target_ptr = Libderiv->deriv_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+17204, dvrr_stack+4439, NULL);
 tmp = dvrr_stack + 17204;
 target_ptr = Libderiv->deriv_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+9317, dvrr_stack+7607, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+17414, dvrr_stack+5699, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+26203, dvrr_stack+12321, NULL);
 tmp = dvrr_stack + 26203;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+26428, dvrr_stack+11157, NULL);
 tmp = dvrr_stack + 26428;
 target_ptr = Libderiv->deriv_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+26578, dvrr_stack+15214, NULL);
 tmp = dvrr_stack + 26578;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+26893, dvrr_stack+19408, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+22488, dvrr_stack+20773, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+20773, dvrr_stack+2073, NULL);
 tmp = dvrr_stack + 20773;
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
 deriv_build_DX_0(Data,280,dvrr_stack+4439, dvrr_stack+7607, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+1233, dvrr_stack+5699, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+5699, dvrr_stack+12321, NULL);
 tmp = dvrr_stack + 5699;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+12321, dvrr_stack+11157, NULL);
 tmp = dvrr_stack + 12321;
 target_ptr = Libderiv->deriv_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+11157, dvrr_stack+15214, NULL);
 tmp = dvrr_stack + 11157;
 target_ptr = Libderiv->deriv_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+15214, dvrr_stack+19408, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+164, dvrr_stack+3, dvrr_stack+3672, Data->F+1, Data->F+2, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+19408, dvrr_stack+161, dvrr_stack+21, NULL, NULL, Data->F+1);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+57, dvrr_stack+19408, dvrr_stack+20668, dvrr_stack+161, dvrr_stack+21, dvrr_stack+3);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+20923, dvrr_stack+57, dvrr_stack+20677, dvrr_stack+19408, dvrr_stack+20668, dvrr_stack+164);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+19408, dvrr_stack+379, dvrr_stack+20923);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+19468, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 19468;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+19618, dvrr_stack+873, dvrr_stack+20713);
 tmp = dvrr_stack + 19618;
 target_ptr = Libderiv->deriv_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+19718, dvrr_stack+4159, dvrr_stack+873);
 tmp = dvrr_stack + 19718;
 target_ptr = Libderiv->deriv_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+19928, dvrr_stack+7247, dvrr_stack+1863);

 /* compute (3 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+22548, dvrr_stack+164, dvrr_stack+24, dvrr_stack+3, dvrr_stack+3672, NULL);

 /* compute (4 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+5924, dvrr_stack+20923, dvrr_stack+21603, dvrr_stack+57, dvrr_stack+20677, dvrr_stack+22548);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,15,1,dvrr_stack+20208, dvrr_stack+10497, dvrr_stack+5924);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+20298, dvrr_stack+12006, dvrr_stack+10497);
 tmp = dvrr_stack + 20298;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+20523, dvrr_stack+10932, dvrr_stack+21633);
 tmp = dvrr_stack + 20523;
 target_ptr = Libderiv->deriv_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+15634, dvrr_stack+14794, dvrr_stack+10932);
 tmp = dvrr_stack + 15634;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+12471, dvrr_stack+18868, dvrr_stack+12006);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+15949, dvrr_stack+379, dvrr_stack+20923);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+16009, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 16009;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+11472, dvrr_stack+873, dvrr_stack+20713);
 tmp = dvrr_stack + 11472;
 target_ptr = Libderiv->deriv_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+7607, dvrr_stack+4159, dvrr_stack+873);
 tmp = dvrr_stack + 7607;
 target_ptr = Libderiv->deriv_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+7817, dvrr_stack+7247, dvrr_stack+1863);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,15,1,dvrr_stack+12891, dvrr_stack+10497, dvrr_stack+5924);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+8097, dvrr_stack+12006, dvrr_stack+10497);
 tmp = dvrr_stack + 8097;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+4719, dvrr_stack+10932, dvrr_stack+21633);
 tmp = dvrr_stack + 4719;
 target_ptr = Libderiv->deriv_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+2173, dvrr_stack+14794, dvrr_stack+10932);
 tmp = dvrr_stack + 2173;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+27313, dvrr_stack+18868, dvrr_stack+12006);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+8322, dvrr_stack+379, dvrr_stack+20923);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+4869, dvrr_stack+1863, dvrr_stack+379);
 tmp = dvrr_stack + 4869;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+9597, dvrr_stack+873, dvrr_stack+20713);
 tmp = dvrr_stack + 9597;
 target_ptr = Libderiv->deriv_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+27733, dvrr_stack+4159, dvrr_stack+873);
 tmp = dvrr_stack + 27733;
 target_ptr = Libderiv->deriv_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+27943, dvrr_stack+7247, dvrr_stack+1863);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,15,1,dvrr_stack+6419, dvrr_stack+10497, dvrr_stack+5924);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+7103, dvrr_stack+12006, dvrr_stack+10497);
 tmp = dvrr_stack + 7103;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+7328, dvrr_stack+10932, dvrr_stack+21633);
 tmp = dvrr_stack + 7328;
 target_ptr = Libderiv->deriv_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+28223, dvrr_stack+14794, dvrr_stack+10932);
 tmp = dvrr_stack + 28223;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+28538, dvrr_stack+18868, dvrr_stack+12006);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+18784,dvrr_stack+783,dvrr_stack+319,6);


 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+8382, dvrr_stack+18784, NULL);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+18964,dvrr_stack+1737,dvrr_stack+783,6);


 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+19234, dvrr_stack+18964, NULL);

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+28958,dvrr_stack+3991,dvrr_stack+1737,6);


 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+7478, dvrr_stack+28958, NULL);

 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+5924, dvrr_stack+121, dvrr_stack+3, Data->F+3, Data->F+4, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+15, dvrr_stack+30, dvrr_stack+210, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+57, dvrr_stack+3627, dvrr_stack+15, dvrr_stack+0, dvrr_stack+30, dvrr_stack+3);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+20923, dvrr_stack+3636, dvrr_stack+57, dvrr_stack+6, dvrr_stack+3627, dvrr_stack+5924);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+20695, dvrr_stack+213, dvrr_stack+617, NULL, NULL, dvrr_stack+210);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+5930, dvrr_stack+3654, dvrr_stack+20695, dvrr_stack+33, dvrr_stack+213, dvrr_stack+15);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+19324, dvrr_stack+6527, dvrr_stack+5930, dvrr_stack+39, dvrr_stack+3654, dvrr_stack+57);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+29336, dvrr_stack+6563, dvrr_stack+19324, dvrr_stack+75, dvrr_stack+6527, dvrr_stack+20923);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+75, dvrr_stack+623, dvrr_stack+170, NULL, NULL, dvrr_stack+617);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+479, dvrr_stack+6623, dvrr_stack+75, dvrr_stack+219, dvrr_stack+623, dvrr_stack+20695);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+29426, dvrr_stack+10337, dvrr_stack+479, dvrr_stack+229, dvrr_stack+6623, dvrr_stack+5930);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+29526, dvrr_stack+10397, dvrr_stack+29426, dvrr_stack+259, dvrr_stack+10337, dvrr_stack+19324);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+29676, dvrr_stack+10497, dvrr_stack+29526, dvrr_stack+379, dvrr_stack+10397, dvrr_stack+29336);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+219, dvrr_stack+1512, dvrr_stack+3675, NULL, NULL, dvrr_stack+170);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+29886, dvrr_stack+10647, dvrr_stack+219, dvrr_stack+633, dvrr_stack+1512, dvrr_stack+75);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+29976, dvrr_stack+10692, dvrr_stack+29886, dvrr_stack+648, dvrr_stack+10647, dvrr_stack+479);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+30126, dvrr_stack+10782, dvrr_stack+29976, dvrr_stack+693, dvrr_stack+10692, dvrr_stack+29426);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+30351, dvrr_stack+10932, dvrr_stack+30126, dvrr_stack+873, dvrr_stack+10782, dvrr_stack+29526);

 /* compute (5 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+30666,dvrr_stack+30351,dvrr_stack+29676,21);


 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+31296, dvrr_stack+30666, NULL);

 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+633, dvrr_stack+3690, dvrr_stack+1344, NULL, NULL, dvrr_stack+3675);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+31506, dvrr_stack+11607, dvrr_stack+633, dvrr_stack+1527, dvrr_stack+3690, dvrr_stack+219);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+31632, dvrr_stack+11670, dvrr_stack+31506, dvrr_stack+1548, dvrr_stack+11607, dvrr_stack+29886);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+31842, dvrr_stack+11796, dvrr_stack+31632, dvrr_stack+1611, dvrr_stack+11670, dvrr_stack+29976);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+32157, dvrr_stack+12006, dvrr_stack+31842, dvrr_stack+1863, dvrr_stack+11796, dvrr_stack+30126);

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+32598,dvrr_stack+32157,dvrr_stack+30351,21);


 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+33543, dvrr_stack+32598, NULL);

 /* compute (1 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1527, dvrr_stack+3563, dvrr_stack+3479, NULL, NULL, dvrr_stack+1344);

 /* compute (2 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+33858, dvrr_stack+1365, dvrr_stack+1527, dvrr_stack+3711, dvrr_stack+3563, dvrr_stack+633);

 /* compute (3 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+34026, dvrr_stack+14346, dvrr_stack+33858, dvrr_stack+3739, dvrr_stack+1365, dvrr_stack+31506);

 /* compute (4 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+34306, dvrr_stack+14514, dvrr_stack+34026, dvrr_stack+3823, dvrr_stack+14346, dvrr_stack+31632);

 /* compute (5 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+34726, dvrr_stack+14794, dvrr_stack+34306, dvrr_stack+4159, dvrr_stack+14514, dvrr_stack+31842);

 /* compute (5 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+35314,dvrr_stack+34726,dvrr_stack+32157,21);


 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,441,dvrr_stack+14346, dvrr_stack+35314, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+4159, dvrr_stack+18784, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+4219, dvrr_stack+18964, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+4309, dvrr_stack+28958, NULL);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+1527, dvrr_stack+30666, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+14787, dvrr_stack+32598, NULL);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+33858, dvrr_stack+35314, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+15102, dvrr_stack+18784, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+18784, dvrr_stack+18964, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+18874, dvrr_stack+28958, NULL);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+28958, dvrr_stack+30666, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+30666, dvrr_stack+32598, NULL);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+32598, dvrr_stack+35314, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+35314, dvrr_stack+783, dvrr_stack+3591);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+35374, dvrr_stack+1737, dvrr_stack+319);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+35464, dvrr_stack+3991, dvrr_stack+783);

 /* compute (3 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+22548, dvrr_stack+24, dvrr_stack+5924, dvrr_stack+3672, dvrr_stack+121, NULL);

 /* compute (4 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+35590, dvrr_stack+21603, dvrr_stack+20923, dvrr_stack+20677, dvrr_stack+3636, dvrr_stack+22548);

 /* compute (5 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+35635, dvrr_stack+21633, dvrr_stack+29336, dvrr_stack+20713, dvrr_stack+6563, dvrr_stack+35590);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+35761, dvrr_stack+30351, dvrr_stack+35635);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+30981, dvrr_stack+32157, dvrr_stack+29676);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+35971, dvrr_stack+34726, dvrr_stack+30351);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+20713, dvrr_stack+783, dvrr_stack+3591);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+21603, dvrr_stack+1737, dvrr_stack+319);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+36412, dvrr_stack+3991, dvrr_stack+783);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+36538, dvrr_stack+30351, dvrr_stack+35635);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+33039, dvrr_stack+32157, dvrr_stack+29676);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+36748, dvrr_stack+34726, dvrr_stack+30351);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+37189, dvrr_stack+783, dvrr_stack+3591);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+37249, dvrr_stack+1737, dvrr_stack+319);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+37339, dvrr_stack+3991, dvrr_stack+783);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+37465, dvrr_stack+30351, dvrr_stack+35635);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+37675, dvrr_stack+32157, dvrr_stack+29676);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+3711, dvrr_stack+34726, dvrr_stack+30351);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+35590, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+33354, dvrr_stack+29676, dvrr_stack+379);
 tmp = dvrr_stack + 33354;
 target_ptr = Libderiv->deriv_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+35650, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+19000, dvrr_stack+30351, dvrr_stack+873);
 tmp = dvrr_stack + 19000;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+29168, dvrr_stack+1863, dvrr_stack+1449);

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+34299, dvrr_stack+32157, dvrr_stack+1863);
 tmp = dvrr_stack + 34299;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+34614, dvrr_stack+10497, dvrr_stack+319);
 tmp = dvrr_stack + 34614;
 target_ptr = Libderiv->deriv_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3672, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+24, dvrr_stack+3, dvrr_stack+3672, Data->F+4, Data->F+5, NULL);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+22548, dvrr_stack+5924, dvrr_stack+24, dvrr_stack+121, dvrr_stack+3, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+19225, dvrr_stack+210, dvrr_stack+614, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+6509, dvrr_stack+15, dvrr_stack+19225, dvrr_stack+30, dvrr_stack+210, dvrr_stack+3672);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+21693, dvrr_stack+57, dvrr_stack+6509, dvrr_stack+3627, dvrr_stack+15, dvrr_stack+24);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+15162, dvrr_stack+20923, dvrr_stack+21693, dvrr_stack+3636, dvrr_stack+57, dvrr_stack+22548);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+20923, dvrr_stack+617, dvrr_stack+124, NULL, NULL, dvrr_stack+614);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+33504, dvrr_stack+20695, dvrr_stack+20923, dvrr_stack+213, dvrr_stack+617, dvrr_stack+19225);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+1365, dvrr_stack+5930, dvrr_stack+33504, dvrr_stack+3654, dvrr_stack+20695, dvrr_stack+6509);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+34714, dvrr_stack+19324, dvrr_stack+1365, dvrr_stack+6527, dvrr_stack+5930, dvrr_stack+21693);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+34804, dvrr_stack+29336, dvrr_stack+34714, dvrr_stack+6563, dvrr_stack+19324, dvrr_stack+15162);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+21693, dvrr_stack+170, dvrr_stack+554, NULL, NULL, dvrr_stack+124);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+19324, dvrr_stack+75, dvrr_stack+21693, dvrr_stack+623, dvrr_stack+170, dvrr_stack+20923);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+6509, dvrr_stack+479, dvrr_stack+19324, dvrr_stack+6623, dvrr_stack+75, dvrr_stack+33504);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+34930, dvrr_stack+29426, dvrr_stack+6509, dvrr_stack+10337, dvrr_stack+479, dvrr_stack+1365);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+35080, dvrr_stack+29526, dvrr_stack+34930, dvrr_stack+10397, dvrr_stack+29426, dvrr_stack+34714);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+37990, dvrr_stack+29676, dvrr_stack+35080, dvrr_stack+10497, dvrr_stack+29526, dvrr_stack+34804);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+34714, dvrr_stack+37990, dvrr_stack+10497);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+10337, dvrr_stack+10932, dvrr_stack+783);
 tmp = dvrr_stack + 10337;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5924, dvrr_stack+3675, dvrr_stack+3464, NULL, NULL, dvrr_stack+554);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+479, dvrr_stack+219, dvrr_stack+5924, dvrr_stack+1512, dvrr_stack+3675, dvrr_stack+21693);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+38270, dvrr_stack+29886, dvrr_stack+479, dvrr_stack+10647, dvrr_stack+219, dvrr_stack+19324);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+38420, dvrr_stack+29976, dvrr_stack+38270, dvrr_stack+10692, dvrr_stack+29886, dvrr_stack+6509);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+38645, dvrr_stack+30126, dvrr_stack+38420, dvrr_stack+10782, dvrr_stack+29976, dvrr_stack+34930);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+38960, dvrr_stack+30351, dvrr_stack+38645, dvrr_stack+10932, dvrr_stack+30126, dvrr_stack+35080);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+29886, dvrr_stack+38960, dvrr_stack+10932);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+10647, dvrr_stack+12006, dvrr_stack+1737);
 tmp = dvrr_stack + 10647;
 target_ptr = Libderiv->deriv_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6509, dvrr_stack+1344, dvrr_stack+1323, NULL, NULL, dvrr_stack+3464);

 /* compute (2 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+30201, dvrr_stack+633, dvrr_stack+6509, dvrr_stack+3690, dvrr_stack+1344, dvrr_stack+5924);

 /* compute (3 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+34924, dvrr_stack+31506, dvrr_stack+30201, dvrr_stack+11607, dvrr_stack+633, dvrr_stack+479);

 /* compute (4 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+39380, dvrr_stack+31632, dvrr_stack+34924, dvrr_stack+11670, dvrr_stack+31506, dvrr_stack+38270);

 /* compute (5 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+39695, dvrr_stack+31842, dvrr_stack+39380, dvrr_stack+11796, dvrr_stack+31632, dvrr_stack+38420);

 /* compute (6 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40136, dvrr_stack+32157, dvrr_stack+39695, dvrr_stack+12006, dvrr_stack+31842, dvrr_stack+38645);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+39380, dvrr_stack+40136, dvrr_stack+12006);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+39821, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+30201, dvrr_stack+29676, dvrr_stack+379);
 tmp = dvrr_stack + 30201;
 target_ptr = Libderiv->deriv_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+479, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+39881, dvrr_stack+30351, dvrr_stack+873);
 tmp = dvrr_stack + 39881;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+1323, dvrr_stack+1863, dvrr_stack+1449);

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+38270, dvrr_stack+32157, dvrr_stack+1863);
 tmp = dvrr_stack + 38270;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+38585, dvrr_stack+10497, dvrr_stack+319);
 tmp = dvrr_stack + 38585;
 target_ptr = Libderiv->deriv_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+38685, dvrr_stack+37990, dvrr_stack+10497);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+31506, dvrr_stack+10932, dvrr_stack+783);
 tmp = dvrr_stack + 31506;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+31656, dvrr_stack+38960, dvrr_stack+10932);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+34924, dvrr_stack+12006, dvrr_stack+1737);
 tmp = dvrr_stack + 34924;
 target_ptr = Libderiv->deriv_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+40724, dvrr_stack+40136, dvrr_stack+12006);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+38895, dvrr_stack+379, dvrr_stack+180);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+31971, dvrr_stack+29676, dvrr_stack+379);
 tmp = dvrr_stack + 31971;
 target_ptr = Libderiv->deriv_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+379, dvrr_stack+873, dvrr_stack+569);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+41165, dvrr_stack+30351, dvrr_stack+873);
 tmp = dvrr_stack + 41165;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+873, dvrr_stack+1863, dvrr_stack+1449);

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+30351, dvrr_stack+32157, dvrr_stack+1863);
 tmp = dvrr_stack + 30351;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+1863, dvrr_stack+10497, dvrr_stack+319);
 tmp = dvrr_stack + 1863;
 target_ptr = Libderiv->deriv_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+569, dvrr_stack+37990, dvrr_stack+10497);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+37990, dvrr_stack+10932, dvrr_stack+783);
 tmp = dvrr_stack + 37990;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+32121, dvrr_stack+38960, dvrr_stack+10932);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+38955, dvrr_stack+12006, dvrr_stack+1737);
 tmp = dvrr_stack + 38955;
 target_ptr = Libderiv->deriv_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+41390, dvrr_stack+40136, dvrr_stack+12006);

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,100,dvrr_stack+1737, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 1737;
 target_ptr = Libderiv->deriv2_classes[3][3][143];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+40106, dvrr_stack+5969, NULL);
 tmp = dvrr_stack + 40106;
 target_ptr = Libderiv->deriv2_classes[3][4][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,210,dvrr_stack+40256, dvrr_stack+9707, NULL);
 tmp = dvrr_stack + 40256;
 target_ptr = Libderiv->deriv2_classes[3][5][143];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+40466, dvrr_stack+13896, NULL);
 tmp = dvrr_stack + 40466;
 target_ptr = Libderiv->deriv2_classes[4][3][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,225,dvrr_stack+10857, dvrr_stack+17509, NULL);
 tmp = dvrr_stack + 10857;
 target_ptr = Libderiv->deriv2_classes[4][4][143];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,315,dvrr_stack+29294, dvrr_stack+22558, NULL);
 tmp = dvrr_stack + 29294;
 target_ptr = Libderiv->deriv2_classes[4][5][143];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,100,dvrr_stack+38140, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 38140;
 target_ptr = Libderiv->deriv2_classes[3][3][131];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+39165, dvrr_stack+5969, NULL);
 tmp = dvrr_stack + 39165;
 target_ptr = Libderiv->deriv2_classes[3][4][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,210,dvrr_stack+29609, dvrr_stack+9707, NULL);
 tmp = dvrr_stack + 29609;
 target_ptr = Libderiv->deriv2_classes[3][5][131];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+32436, dvrr_stack+13896, NULL);
 tmp = dvrr_stack + 32436;
 target_ptr = Libderiv->deriv2_classes[4][3][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,225,dvrr_stack+0, dvrr_stack+17509, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[4][4][131];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,315,dvrr_stack+11572, dvrr_stack+22558, NULL);
 tmp = dvrr_stack + 11572;
 target_ptr = Libderiv->deriv2_classes[4][5][131];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+1963, dvrr_stack+18184, NULL);
 tmp = dvrr_stack + 1963;
 target_ptr = Libderiv->deriv2_classes[3][3][130];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+35134, dvrr_stack+6653, NULL);
 tmp = dvrr_stack + 35134;
 target_ptr = Libderiv->deriv2_classes[3][4][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+11887, dvrr_stack+23503, NULL);
 tmp = dvrr_stack + 11887;
 target_ptr = Libderiv->deriv2_classes[3][5][130];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+10487, dvrr_stack+24133, NULL);
 tmp = dvrr_stack + 10487;
 target_ptr = Libderiv->deriv2_classes[4][3][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+3423, dvrr_stack+24583, NULL);
 tmp = dvrr_stack + 3423;
 target_ptr = Libderiv->deriv2_classes[4][4][130];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+41831, dvrr_stack+25258, NULL);
 tmp = dvrr_stack + 41831;
 target_ptr = Libderiv->deriv2_classes[4][5][130];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,100,dvrr_stack+40616, dvrr_stack+3123, NULL);
 tmp = dvrr_stack + 40616;
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
 deriv_build_DX_0(Data,210,dvrr_stack+5924, dvrr_stack+9707, NULL);
 tmp = dvrr_stack + 5924;
 target_ptr = Libderiv->deriv2_classes[3][5][119];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+3273, dvrr_stack+13896, NULL);
 tmp = dvrr_stack + 3273;
 target_ptr = Libderiv->deriv2_classes[4][3][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,225,dvrr_stack+13896, dvrr_stack+17509, NULL);
 tmp = dvrr_stack + 13896;
 target_ptr = Libderiv->deriv2_classes[4][4][119];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,315,dvrr_stack+42146, dvrr_stack+22558, NULL);
 tmp = dvrr_stack + 42146;
 target_ptr = Libderiv->deriv2_classes[4][5][119];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+14121, dvrr_stack+18184, NULL);
 tmp = dvrr_stack + 14121;
 target_ptr = Libderiv->deriv2_classes[3][3][118];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+6134, dvrr_stack+6653, NULL);
 tmp = dvrr_stack + 6134;
 target_ptr = Libderiv->deriv2_classes[3][4][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+6509, dvrr_stack+23503, NULL);
 tmp = dvrr_stack + 6509;
 target_ptr = Libderiv->deriv2_classes[3][5][118];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+6719, dvrr_stack+24133, NULL);
 tmp = dvrr_stack + 6719;
 target_ptr = Libderiv->deriv2_classes[4][3][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+6869, dvrr_stack+24583, NULL);
 tmp = dvrr_stack + 6869;
 target_ptr = Libderiv->deriv2_classes[4][4][118];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+22548, dvrr_stack+25258, NULL);
 tmp = dvrr_stack + 22548;
 target_ptr = Libderiv->deriv2_classes[4][5][118];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+14221, dvrr_stack+18484, NULL);
 tmp = dvrr_stack + 14221;
 target_ptr = Libderiv->deriv2_classes[3][3][117];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+22863, dvrr_stack+2523, NULL);
 tmp = dvrr_stack + 22863;
 target_ptr = Libderiv->deriv2_classes[3][4][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+23013, dvrr_stack+5069, NULL);
 tmp = dvrr_stack + 23013;
 target_ptr = Libderiv->deriv2_classes[3][5][117];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+23223, dvrr_stack+8447, NULL);
 tmp = dvrr_stack + 23223;
 target_ptr = Libderiv->deriv2_classes[4][3][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+23373, dvrr_stack+12996, NULL);
 tmp = dvrr_stack + 23373;
 target_ptr = Libderiv->deriv2_classes[4][4][117];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+23598, dvrr_stack+16159, NULL);
 tmp = dvrr_stack + 23598;
 target_ptr = Libderiv->deriv2_classes[4][5][117];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+16159, dvrr_stack+2973, dvrr_stack+20953);
 tmp = dvrr_stack + 16159;
 target_ptr = Libderiv->deriv2_classes[3][3][107];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+16259, dvrr_stack+21113, dvrr_stack+21013);
 tmp = dvrr_stack + 16259;
 target_ptr = Libderiv->deriv2_classes[3][4][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+16409, dvrr_stack+21323, dvrr_stack+2973);
 tmp = dvrr_stack + 16409;
 target_ptr = Libderiv->deriv2_classes[3][5][107];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+16619, dvrr_stack+13671, dvrr_stack+21723);
 tmp = dvrr_stack + 16619;
 target_ptr = Libderiv->deriv2_classes[4][3][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+16769, dvrr_stack+21963, dvrr_stack+21813);
 tmp = dvrr_stack + 16769;
 target_ptr = Libderiv->deriv2_classes[4][4][107];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+23913, dvrr_stack+8897, dvrr_stack+13671);
 tmp = dvrr_stack + 23913;
 target_ptr = Libderiv->deriv2_classes[4][5][107];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+16994, dvrr_stack+22338, dvrr_stack+22278);
 tmp = dvrr_stack + 16994;
 target_ptr = Libderiv->deriv2_classes[3][3][106];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+24228, dvrr_stack+17204, dvrr_stack+17104);
 tmp = dvrr_stack + 24228;
 target_ptr = Libderiv->deriv2_classes[3][4][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+24378, dvrr_stack+9317, dvrr_stack+22338);
 tmp = dvrr_stack + 24378;
 target_ptr = Libderiv->deriv2_classes[3][5][106];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+24588, dvrr_stack+26203, dvrr_stack+17414);
 tmp = dvrr_stack + 24588;
 target_ptr = Libderiv->deriv2_classes[4][3][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+24738, dvrr_stack+26578, dvrr_stack+26428);
 tmp = dvrr_stack + 24738;
 target_ptr = Libderiv->deriv2_classes[4][4][106];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+24963, dvrr_stack+26893, dvrr_stack+26203);
 tmp = dvrr_stack + 24963;
 target_ptr = Libderiv->deriv2_classes[4][5][106];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+6284, dvrr_stack+20773, dvrr_stack+22488);
 tmp = dvrr_stack + 6284;
 target_ptr = Libderiv->deriv2_classes[3][3][105];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+25278, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 25278;
 target_ptr = Libderiv->deriv2_classes[3][4][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+25428, dvrr_stack+4439, dvrr_stack+20773);
 tmp = dvrr_stack + 25428;
 target_ptr = Libderiv->deriv2_classes[3][5][105];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+25638, dvrr_stack+5699, dvrr_stack+1233);
 tmp = dvrr_stack + 25638;
 target_ptr = Libderiv->deriv2_classes[4][3][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+25788, dvrr_stack+11157, dvrr_stack+12321);
 tmp = dvrr_stack + 25788;
 target_ptr = Libderiv->deriv2_classes[4][4][105];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+42461, dvrr_stack+15214, dvrr_stack+5699);
 tmp = dvrr_stack + 42461;
 target_ptr = Libderiv->deriv2_classes[4][5][105];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+26013, dvrr_stack+19468, dvrr_stack+19408);
 tmp = dvrr_stack + 26013;
 target_ptr = Libderiv->deriv2_classes[3][3][104];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+42776, dvrr_stack+19718, dvrr_stack+19618);
 tmp = dvrr_stack + 42776;
 target_ptr = Libderiv->deriv2_classes[3][4][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+12097, dvrr_stack+19928, dvrr_stack+19468);
 tmp = dvrr_stack + 12097;
 target_ptr = Libderiv->deriv2_classes[3][5][104];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+225, dvrr_stack+20298, dvrr_stack+20208);
 tmp = dvrr_stack + 225;
 target_ptr = Libderiv->deriv2_classes[4][3][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+12981, dvrr_stack+15634, dvrr_stack+20523);
 tmp = dvrr_stack + 12981;
 target_ptr = Libderiv->deriv2_classes[4][4][104];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+13206, dvrr_stack+12471, dvrr_stack+20298);
 tmp = dvrr_stack + 13206;
 target_ptr = Libderiv->deriv2_classes[4][5][104];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+13521, dvrr_stack+2973, dvrr_stack+20953);
 tmp = dvrr_stack + 13521;
 target_ptr = Libderiv->deriv2_classes[3][3][95];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+8442, dvrr_stack+21113, dvrr_stack+21013);
 tmp = dvrr_stack + 8442;
 target_ptr = Libderiv->deriv2_classes[3][4][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+8592, dvrr_stack+21323, dvrr_stack+2973);
 tmp = dvrr_stack + 8592;
 target_ptr = Libderiv->deriv2_classes[3][5][95];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+5019, dvrr_stack+13671, dvrr_stack+21723);
 tmp = dvrr_stack + 5019;
 target_ptr = Libderiv->deriv2_classes[4][3][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+5169, dvrr_stack+21963, dvrr_stack+21813);
 tmp = dvrr_stack + 5169;
 target_ptr = Libderiv->deriv2_classes[4][4][95];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+2488, dvrr_stack+8897, dvrr_stack+13671);
 tmp = dvrr_stack + 2488;
 target_ptr = Libderiv->deriv2_classes[4][5][95];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+5394, dvrr_stack+22338, dvrr_stack+22278);
 tmp = dvrr_stack + 5394;
 target_ptr = Libderiv->deriv2_classes[3][3][94];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+5494, dvrr_stack+17204, dvrr_stack+17104);
 tmp = dvrr_stack + 5494;
 target_ptr = Libderiv->deriv2_classes[3][4][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+17504, dvrr_stack+9317, dvrr_stack+22338);
 tmp = dvrr_stack + 17504;
 target_ptr = Libderiv->deriv2_classes[3][5][94];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+2803, dvrr_stack+26203, dvrr_stack+17414);
 tmp = dvrr_stack + 2803;
 target_ptr = Libderiv->deriv2_classes[4][3][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+17714, dvrr_stack+26578, dvrr_stack+26428);
 tmp = dvrr_stack + 17714;
 target_ptr = Libderiv->deriv2_classes[4][4][94];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+17939, dvrr_stack+26893, dvrr_stack+26203);
 tmp = dvrr_stack + 17939;
 target_ptr = Libderiv->deriv2_classes[4][5][94];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+18254, dvrr_stack+20773, dvrr_stack+22488);
 tmp = dvrr_stack + 18254;
 target_ptr = Libderiv->deriv2_classes[3][3][93];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+18354, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 18354;
 target_ptr = Libderiv->deriv2_classes[3][4][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+18504, dvrr_stack+4439, dvrr_stack+20773);
 tmp = dvrr_stack + 18504;
 target_ptr = Libderiv->deriv2_classes[3][5][93];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+9697, dvrr_stack+5699, dvrr_stack+1233);
 tmp = dvrr_stack + 9697;
 target_ptr = Libderiv->deriv2_classes[4][3][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+9847, dvrr_stack+11157, dvrr_stack+12321);
 tmp = dvrr_stack + 9847;
 target_ptr = Libderiv->deriv2_classes[4][4][93];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+42926, dvrr_stack+15214, dvrr_stack+5699);
 tmp = dvrr_stack + 42926;
 target_ptr = Libderiv->deriv2_classes[4][5][93];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+43241, dvrr_stack+19468, dvrr_stack+19408);
 tmp = dvrr_stack + 43241;
 target_ptr = Libderiv->deriv2_classes[3][3][92];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+43341, dvrr_stack+19718, dvrr_stack+19618);
 tmp = dvrr_stack + 43341;
 target_ptr = Libderiv->deriv2_classes[3][4][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+43491, dvrr_stack+19928, dvrr_stack+19468);
 tmp = dvrr_stack + 43491;
 target_ptr = Libderiv->deriv2_classes[3][5][92];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+43701, dvrr_stack+20298, dvrr_stack+20208);
 tmp = dvrr_stack + 43701;
 target_ptr = Libderiv->deriv2_classes[4][3][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+10072, dvrr_stack+15634, dvrr_stack+20523);
 tmp = dvrr_stack + 10072;
 target_ptr = Libderiv->deriv2_classes[4][4][92];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+43851, dvrr_stack+12471, dvrr_stack+20298);
 tmp = dvrr_stack + 43851;
 target_ptr = Libderiv->deriv2_classes[4][5][92];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+44166, dvrr_stack+16009, dvrr_stack+15949);
 tmp = dvrr_stack + 44166;
 target_ptr = Libderiv->deriv2_classes[3][3][91];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+44266, dvrr_stack+7607, dvrr_stack+11472);
 tmp = dvrr_stack + 44266;
 target_ptr = Libderiv->deriv2_classes[3][4][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+44416, dvrr_stack+7817, dvrr_stack+16009);
 tmp = dvrr_stack + 44416;
 target_ptr = Libderiv->deriv2_classes[3][5][91];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+44626, dvrr_stack+8097, dvrr_stack+12891);
 tmp = dvrr_stack + 44626;
 target_ptr = Libderiv->deriv2_classes[4][3][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+44776, dvrr_stack+2173, dvrr_stack+4719);
 tmp = dvrr_stack + 44776;
 target_ptr = Libderiv->deriv2_classes[4][4][91];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+45001, dvrr_stack+27313, dvrr_stack+8097);
 tmp = dvrr_stack + 45001;
 target_ptr = Libderiv->deriv2_classes[4][5][91];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+45316, dvrr_stack+2973, dvrr_stack+20953);
 tmp = dvrr_stack + 45316;
 target_ptr = Libderiv->deriv2_classes[3][3][83];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+45416, dvrr_stack+21113, dvrr_stack+21013);
 tmp = dvrr_stack + 45416;
 target_ptr = Libderiv->deriv2_classes[3][4][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+45566, dvrr_stack+21323, dvrr_stack+2973);
 tmp = dvrr_stack + 45566;
 target_ptr = Libderiv->deriv2_classes[3][5][83];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+21323, dvrr_stack+13671, dvrr_stack+21723);
 tmp = dvrr_stack + 21323;
 target_ptr = Libderiv->deriv2_classes[4][3][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+45776, dvrr_stack+21963, dvrr_stack+21813);
 tmp = dvrr_stack + 45776;
 target_ptr = Libderiv->deriv2_classes[4][4][83];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+46001, dvrr_stack+8897, dvrr_stack+13671);
 tmp = dvrr_stack + 46001;
 target_ptr = Libderiv->deriv2_classes[4][5][83];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+21473, dvrr_stack+22338, dvrr_stack+22278);
 tmp = dvrr_stack + 21473;
 target_ptr = Libderiv->deriv2_classes[3][3][82];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+46316, dvrr_stack+17204, dvrr_stack+17104);
 tmp = dvrr_stack + 46316;
 target_ptr = Libderiv->deriv2_classes[3][4][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+8802, dvrr_stack+9317, dvrr_stack+22338);
 tmp = dvrr_stack + 8802;
 target_ptr = Libderiv->deriv2_classes[3][5][82];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+9012, dvrr_stack+26203, dvrr_stack+17414);
 tmp = dvrr_stack + 9012;
 target_ptr = Libderiv->deriv2_classes[4][3][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+9162, dvrr_stack+26578, dvrr_stack+26428);
 tmp = dvrr_stack + 9162;
 target_ptr = Libderiv->deriv2_classes[4][4][82];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+46466, dvrr_stack+26893, dvrr_stack+26203);
 tmp = dvrr_stack + 46466;
 target_ptr = Libderiv->deriv2_classes[4][5][82];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+26893, dvrr_stack+20773, dvrr_stack+22488);
 tmp = dvrr_stack + 26893;
 target_ptr = Libderiv->deriv2_classes[3][3][81];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+26993, dvrr_stack+1023, dvrr_stack+2073);
 tmp = dvrr_stack + 26993;
 target_ptr = Libderiv->deriv2_classes[3][4][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+9387, dvrr_stack+4439, dvrr_stack+20773);
 tmp = dvrr_stack + 9387;
 target_ptr = Libderiv->deriv2_classes[3][5][81];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+27143, dvrr_stack+5699, dvrr_stack+1233);
 tmp = dvrr_stack + 27143;
 target_ptr = Libderiv->deriv2_classes[4][3][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+46781, dvrr_stack+11157, dvrr_stack+12321);
 tmp = dvrr_stack + 46781;
 target_ptr = Libderiv->deriv2_classes[4][4][81];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+47006, dvrr_stack+15214, dvrr_stack+5699);
 tmp = dvrr_stack + 47006;
 target_ptr = Libderiv->deriv2_classes[4][5][81];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+21693, dvrr_stack+19468, dvrr_stack+19408);
 tmp = dvrr_stack + 21693;
 target_ptr = Libderiv->deriv2_classes[3][3][80];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+15162, dvrr_stack+19718, dvrr_stack+19618);
 tmp = dvrr_stack + 15162;
 target_ptr = Libderiv->deriv2_classes[3][4][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+15312, dvrr_stack+19928, dvrr_stack+19468);
 tmp = dvrr_stack + 15312;
 target_ptr = Libderiv->deriv2_classes[3][5][80];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+19928, dvrr_stack+20298, dvrr_stack+20208);
 tmp = dvrr_stack + 19928;
 target_ptr = Libderiv->deriv2_classes[4][3][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+4435, dvrr_stack+15634, dvrr_stack+20523);
 tmp = dvrr_stack + 4435;
 target_ptr = Libderiv->deriv2_classes[4][4][80];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+47321, dvrr_stack+12471, dvrr_stack+20298);
 tmp = dvrr_stack + 47321;
 target_ptr = Libderiv->deriv2_classes[4][5][80];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+12471, dvrr_stack+16009, dvrr_stack+15949);
 tmp = dvrr_stack + 12471;
 target_ptr = Libderiv->deriv2_classes[3][3][79];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+12571, dvrr_stack+7607, dvrr_stack+11472);
 tmp = dvrr_stack + 12571;
 target_ptr = Libderiv->deriv2_classes[3][4][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+20078, dvrr_stack+7817, dvrr_stack+16009);
 tmp = dvrr_stack + 20078;
 target_ptr = Libderiv->deriv2_classes[3][5][79];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+7817, dvrr_stack+8097, dvrr_stack+12891);
 tmp = dvrr_stack + 7817;
 target_ptr = Libderiv->deriv2_classes[4][3][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+12721, dvrr_stack+2173, dvrr_stack+4719);
 tmp = dvrr_stack + 12721;
 target_ptr = Libderiv->deriv2_classes[4][4][79];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+47636, dvrr_stack+27313, dvrr_stack+8097);
 tmp = dvrr_stack + 47636;
 target_ptr = Libderiv->deriv2_classes[4][5][79];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+7967, dvrr_stack+4869, dvrr_stack+8322);
 tmp = dvrr_stack + 7967;
 target_ptr = Libderiv->deriv2_classes[3][3][78];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+27293, dvrr_stack+27733, dvrr_stack+9597);
 tmp = dvrr_stack + 27293;
 target_ptr = Libderiv->deriv2_classes[3][4][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+27443, dvrr_stack+27943, dvrr_stack+4869);
 tmp = dvrr_stack + 27443;
 target_ptr = Libderiv->deriv2_classes[3][5][78];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+27943, dvrr_stack+7103, dvrr_stack+6419);
 tmp = dvrr_stack + 27943;
 target_ptr = Libderiv->deriv2_classes[4][3][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+47951, dvrr_stack+28223, dvrr_stack+7328);
 tmp = dvrr_stack + 47951;
 target_ptr = Libderiv->deriv2_classes[4][4][78];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+48176, dvrr_stack+28538, dvrr_stack+7103);
 tmp = dvrr_stack + 48176;
 target_ptr = Libderiv->deriv2_classes[4][5][78];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,10,dvrr_stack+28538, dvrr_stack+21813, dvrr_stack+8382);
 tmp = dvrr_stack + 28538;
 target_ptr = Libderiv->deriv2_classes[3][3][35];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,15,dvrr_stack+28638, dvrr_stack+13671, dvrr_stack+19234);
 tmp = dvrr_stack + 28638;
 target_ptr = Libderiv->deriv2_classes[3][4][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,21,dvrr_stack+48491, dvrr_stack+21963, dvrr_stack+7478);
 tmp = dvrr_stack + 48491;
 target_ptr = Libderiv->deriv2_classes[3][5][35];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,10,dvrr_stack+28788, dvrr_stack+31296, dvrr_stack+21013);
 tmp = dvrr_stack + 28788;
 target_ptr = Libderiv->deriv2_classes[4][3][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,15,dvrr_stack+48701, dvrr_stack+33543, dvrr_stack+2973);
 tmp = dvrr_stack + 48701;
 target_ptr = Libderiv->deriv2_classes[4][4][35];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,21,dvrr_stack+48926, dvrr_stack+14346, dvrr_stack+21113);
 tmp = dvrr_stack + 48926;
 target_ptr = Libderiv->deriv2_classes[4][5][35];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+28093, dvrr_stack+26428, dvrr_stack+4159);
 tmp = dvrr_stack + 28093;
 target_ptr = Libderiv->deriv2_classes[3][3][34];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+49241, dvrr_stack+26203, dvrr_stack+4219);
 tmp = dvrr_stack + 49241;
 target_ptr = Libderiv->deriv2_classes[3][4][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+49391, dvrr_stack+26578, dvrr_stack+4309);
 tmp = dvrr_stack + 49391;
 target_ptr = Libderiv->deriv2_classes[3][5][34];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+49601, dvrr_stack+1527, dvrr_stack+17104);
 tmp = dvrr_stack + 49601;
 target_ptr = Libderiv->deriv2_classes[4][3][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+49751, dvrr_stack+14787, dvrr_stack+22338);
 tmp = dvrr_stack + 49751;
 target_ptr = Libderiv->deriv2_classes[4][4][34];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+49976, dvrr_stack+33858, dvrr_stack+17204);
 tmp = dvrr_stack + 49976;
 target_ptr = Libderiv->deriv2_classes[4][5][34];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+6384, dvrr_stack+12321, dvrr_stack+15102);
 tmp = dvrr_stack + 6384;
 target_ptr = Libderiv->deriv2_classes[3][3][33];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+50291, dvrr_stack+5699, dvrr_stack+18784);
 tmp = dvrr_stack + 50291;
 target_ptr = Libderiv->deriv2_classes[3][4][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+50441, dvrr_stack+11157, dvrr_stack+18874);
 tmp = dvrr_stack + 50441;
 target_ptr = Libderiv->deriv2_classes[3][5][33];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+50651, dvrr_stack+28958, dvrr_stack+2073);
 tmp = dvrr_stack + 50651;
 target_ptr = Libderiv->deriv2_classes[4][3][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+50801, dvrr_stack+30666, dvrr_stack+20773);
 tmp = dvrr_stack + 50801;
 target_ptr = Libderiv->deriv2_classes[4][4][33];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+51026, dvrr_stack+32598, dvrr_stack+1023);
 tmp = dvrr_stack + 51026;
 target_ptr = Libderiv->deriv2_classes[4][5][33];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+19324, dvrr_stack+20523, dvrr_stack+35314);
 tmp = dvrr_stack + 19324;
 target_ptr = Libderiv->deriv2_classes[3][3][32];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+51341, dvrr_stack+20298, dvrr_stack+35374);
 tmp = dvrr_stack + 51341;
 target_ptr = Libderiv->deriv2_classes[3][4][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+51491, dvrr_stack+15634, dvrr_stack+35464);
 tmp = dvrr_stack + 51491;
 target_ptr = Libderiv->deriv2_classes[3][5][32];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+51701, dvrr_stack+35761, dvrr_stack+19618);
 tmp = dvrr_stack + 51701;
 target_ptr = Libderiv->deriv2_classes[4][3][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+51851, dvrr_stack+30981, dvrr_stack+19468);
 tmp = dvrr_stack + 51851;
 target_ptr = Libderiv->deriv2_classes[4][4][32];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+52076, dvrr_stack+35971, dvrr_stack+19718);
 tmp = dvrr_stack + 52076;
 target_ptr = Libderiv->deriv2_classes[4][5][32];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+15522, dvrr_stack+4719, dvrr_stack+20713);
 tmp = dvrr_stack + 15522;
 target_ptr = Libderiv->deriv2_classes[3][3][31];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+52391, dvrr_stack+8097, dvrr_stack+21603);
 tmp = dvrr_stack + 52391;
 target_ptr = Libderiv->deriv2_classes[3][4][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+52541, dvrr_stack+2173, dvrr_stack+36412);
 tmp = dvrr_stack + 52541;
 target_ptr = Libderiv->deriv2_classes[3][5][31];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+52751, dvrr_stack+36538, dvrr_stack+11472);
 tmp = dvrr_stack + 52751;
 target_ptr = Libderiv->deriv2_classes[4][3][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+52901, dvrr_stack+33039, dvrr_stack+16009);
 tmp = dvrr_stack + 52901;
 target_ptr = Libderiv->deriv2_classes[4][4][31];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+53126, dvrr_stack+36748, dvrr_stack+7607);
 tmp = dvrr_stack + 53126;
 target_ptr = Libderiv->deriv2_classes[4][5][31];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+53441, dvrr_stack+7328, dvrr_stack+37189);
 tmp = dvrr_stack + 53441;
 target_ptr = Libderiv->deriv2_classes[3][3][30];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+53541, dvrr_stack+7103, dvrr_stack+37249);
 tmp = dvrr_stack + 53541;
 target_ptr = Libderiv->deriv2_classes[3][4][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+53691, dvrr_stack+28223, dvrr_stack+37339);
 tmp = dvrr_stack + 53691;
 target_ptr = Libderiv->deriv2_classes[3][5][30];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+53901, dvrr_stack+37465, dvrr_stack+9597);
 tmp = dvrr_stack + 53901;
 target_ptr = Libderiv->deriv2_classes[4][3][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+54051, dvrr_stack+37675, dvrr_stack+4869);
 tmp = dvrr_stack + 54051;
 target_ptr = Libderiv->deriv2_classes[4][4][30];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+54276, dvrr_stack+3711, dvrr_stack+27733);
 tmp = dvrr_stack + 54276;
 target_ptr = Libderiv->deriv2_classes[4][5][30];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+54591, dvrr_stack+33354, dvrr_stack+35590);
 tmp = dvrr_stack + 54591;
 target_ptr = Libderiv->deriv2_classes[3][3][26];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+54691, dvrr_stack+19000, dvrr_stack+35650);
 tmp = dvrr_stack + 54691;
 target_ptr = Libderiv->deriv2_classes[3][4][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+54841, dvrr_stack+34299, dvrr_stack+29168);
 tmp = dvrr_stack + 54841;
 target_ptr = Libderiv->deriv2_classes[3][5][26];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+55051, dvrr_stack+34714, dvrr_stack+34614);
 tmp = dvrr_stack + 55051;
 target_ptr = Libderiv->deriv2_classes[4][3][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+55201, dvrr_stack+29886, dvrr_stack+10337);
 tmp = dvrr_stack + 55201;
 target_ptr = Libderiv->deriv2_classes[4][4][26];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+55426, dvrr_stack+39380, dvrr_stack+10647);
 tmp = dvrr_stack + 55426;
 target_ptr = Libderiv->deriv2_classes[4][5][26];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,10,dvrr_stack+55741, dvrr_stack+21813, dvrr_stack+8382);
 tmp = dvrr_stack + 55741;
 target_ptr = Libderiv->deriv2_classes[3][3][23];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,15,dvrr_stack+55841, dvrr_stack+13671, dvrr_stack+19234);
 tmp = dvrr_stack + 55841;
 target_ptr = Libderiv->deriv2_classes[3][4][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,21,dvrr_stack+55991, dvrr_stack+21963, dvrr_stack+7478);
 tmp = dvrr_stack + 55991;
 target_ptr = Libderiv->deriv2_classes[3][5][23];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,10,dvrr_stack+56201, dvrr_stack+31296, dvrr_stack+21013);
 tmp = dvrr_stack + 56201;
 target_ptr = Libderiv->deriv2_classes[4][3][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,15,dvrr_stack+56351, dvrr_stack+33543, dvrr_stack+2973);
 tmp = dvrr_stack + 56351;
 target_ptr = Libderiv->deriv2_classes[4][4][23];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,21,dvrr_stack+56576, dvrr_stack+14346, dvrr_stack+21113);
 tmp = dvrr_stack + 56576;
 target_ptr = Libderiv->deriv2_classes[4][5][23];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+56891, dvrr_stack+26428, dvrr_stack+4159);
 tmp = dvrr_stack + 56891;
 target_ptr = Libderiv->deriv2_classes[3][3][22];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+56991, dvrr_stack+26203, dvrr_stack+4219);
 tmp = dvrr_stack + 56991;
 target_ptr = Libderiv->deriv2_classes[3][4][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+57141, dvrr_stack+26578, dvrr_stack+4309);
 tmp = dvrr_stack + 57141;
 target_ptr = Libderiv->deriv2_classes[3][5][22];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+57351, dvrr_stack+1527, dvrr_stack+17104);
 tmp = dvrr_stack + 57351;
 target_ptr = Libderiv->deriv2_classes[4][3][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+57501, dvrr_stack+14787, dvrr_stack+22338);
 tmp = dvrr_stack + 57501;
 target_ptr = Libderiv->deriv2_classes[4][4][22];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+57726, dvrr_stack+33858, dvrr_stack+17204);
 tmp = dvrr_stack + 57726;
 target_ptr = Libderiv->deriv2_classes[4][5][22];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+58041, dvrr_stack+12321, dvrr_stack+15102);
 tmp = dvrr_stack + 58041;
 target_ptr = Libderiv->deriv2_classes[3][3][21];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+58141, dvrr_stack+5699, dvrr_stack+18784);
 tmp = dvrr_stack + 58141;
 target_ptr = Libderiv->deriv2_classes[3][4][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+58291, dvrr_stack+11157, dvrr_stack+18874);
 tmp = dvrr_stack + 58291;
 target_ptr = Libderiv->deriv2_classes[3][5][21];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+58501, dvrr_stack+28958, dvrr_stack+2073);
 tmp = dvrr_stack + 58501;
 target_ptr = Libderiv->deriv2_classes[4][3][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+58651, dvrr_stack+30666, dvrr_stack+20773);
 tmp = dvrr_stack + 58651;
 target_ptr = Libderiv->deriv2_classes[4][4][21];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+58876, dvrr_stack+32598, dvrr_stack+1023);
 tmp = dvrr_stack + 58876;
 target_ptr = Libderiv->deriv2_classes[4][5][21];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+59191, dvrr_stack+20523, dvrr_stack+35314);
 tmp = dvrr_stack + 59191;
 target_ptr = Libderiv->deriv2_classes[3][3][20];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+59291, dvrr_stack+20298, dvrr_stack+35374);
 tmp = dvrr_stack + 59291;
 target_ptr = Libderiv->deriv2_classes[3][4][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+59441, dvrr_stack+15634, dvrr_stack+35464);
 tmp = dvrr_stack + 59441;
 target_ptr = Libderiv->deriv2_classes[3][5][20];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+59651, dvrr_stack+35761, dvrr_stack+19618);
 tmp = dvrr_stack + 59651;
 target_ptr = Libderiv->deriv2_classes[4][3][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+59801, dvrr_stack+30981, dvrr_stack+19468);
 tmp = dvrr_stack + 59801;
 target_ptr = Libderiv->deriv2_classes[4][4][20];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+60026, dvrr_stack+35971, dvrr_stack+19718);
 tmp = dvrr_stack + 60026;
 target_ptr = Libderiv->deriv2_classes[4][5][20];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+60341, dvrr_stack+4719, dvrr_stack+20713);
 tmp = dvrr_stack + 60341;
 target_ptr = Libderiv->deriv2_classes[3][3][19];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+60441, dvrr_stack+8097, dvrr_stack+21603);
 tmp = dvrr_stack + 60441;
 target_ptr = Libderiv->deriv2_classes[3][4][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+60591, dvrr_stack+2173, dvrr_stack+36412);
 tmp = dvrr_stack + 60591;
 target_ptr = Libderiv->deriv2_classes[3][5][19];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+60801, dvrr_stack+36538, dvrr_stack+11472);
 tmp = dvrr_stack + 60801;
 target_ptr = Libderiv->deriv2_classes[4][3][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+60951, dvrr_stack+33039, dvrr_stack+16009);
 tmp = dvrr_stack + 60951;
 target_ptr = Libderiv->deriv2_classes[4][4][19];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+61176, dvrr_stack+36748, dvrr_stack+7607);
 tmp = dvrr_stack + 61176;
 target_ptr = Libderiv->deriv2_classes[4][5][19];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+61491, dvrr_stack+7328, dvrr_stack+37189);
 tmp = dvrr_stack + 61491;
 target_ptr = Libderiv->deriv2_classes[3][3][18];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+61591, dvrr_stack+7103, dvrr_stack+37249);
 tmp = dvrr_stack + 61591;
 target_ptr = Libderiv->deriv2_classes[3][4][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+61741, dvrr_stack+28223, dvrr_stack+37339);
 tmp = dvrr_stack + 61741;
 target_ptr = Libderiv->deriv2_classes[3][5][18];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+61951, dvrr_stack+37465, dvrr_stack+9597);
 tmp = dvrr_stack + 61951;
 target_ptr = Libderiv->deriv2_classes[4][3][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+62101, dvrr_stack+37675, dvrr_stack+4869);
 tmp = dvrr_stack + 62101;
 target_ptr = Libderiv->deriv2_classes[4][4][18];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+62326, dvrr_stack+3711, dvrr_stack+27733);
 tmp = dvrr_stack + 62326;
 target_ptr = Libderiv->deriv2_classes[4][5][18];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+62641, dvrr_stack+33354, dvrr_stack+35590);
 tmp = dvrr_stack + 62641;
 target_ptr = Libderiv->deriv2_classes[3][3][14];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+62741, dvrr_stack+19000, dvrr_stack+35650);
 tmp = dvrr_stack + 62741;
 target_ptr = Libderiv->deriv2_classes[3][4][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+62891, dvrr_stack+34299, dvrr_stack+29168);
 tmp = dvrr_stack + 62891;
 target_ptr = Libderiv->deriv2_classes[3][5][14];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+63101, dvrr_stack+34714, dvrr_stack+34614);
 tmp = dvrr_stack + 63101;
 target_ptr = Libderiv->deriv2_classes[4][3][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+63251, dvrr_stack+29886, dvrr_stack+10337);
 tmp = dvrr_stack + 63251;
 target_ptr = Libderiv->deriv2_classes[4][4][14];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+63476, dvrr_stack+39380, dvrr_stack+10647);
 tmp = dvrr_stack + 63476;
 target_ptr = Libderiv->deriv2_classes[4][5][14];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+63791, dvrr_stack+30201, dvrr_stack+39821);
 tmp = dvrr_stack + 63791;
 target_ptr = Libderiv->deriv2_classes[3][3][13];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+63891, dvrr_stack+39881, dvrr_stack+479);
 tmp = dvrr_stack + 63891;
 target_ptr = Libderiv->deriv2_classes[3][4][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+64041, dvrr_stack+38270, dvrr_stack+1323);
 tmp = dvrr_stack + 64041;
 target_ptr = Libderiv->deriv2_classes[3][5][13];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+64251, dvrr_stack+38685, dvrr_stack+38585);
 tmp = dvrr_stack + 64251;
 target_ptr = Libderiv->deriv2_classes[4][3][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+64401, dvrr_stack+31656, dvrr_stack+31506);
 tmp = dvrr_stack + 64401;
 target_ptr = Libderiv->deriv2_classes[4][4][13];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+64626, dvrr_stack+40724, dvrr_stack+34924);
 tmp = dvrr_stack + 64626;
 target_ptr = Libderiv->deriv2_classes[4][5][13];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,10,dvrr_stack+64941, dvrr_stack+21813, dvrr_stack+8382);
 tmp = dvrr_stack + 64941;
 target_ptr = Libderiv->deriv2_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,15,dvrr_stack+21793, dvrr_stack+13671, dvrr_stack+19234);
 tmp = dvrr_stack + 21793;
 target_ptr = Libderiv->deriv2_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,21,dvrr_stack+65041, dvrr_stack+21963, dvrr_stack+7478);
 tmp = dvrr_stack + 65041;
 target_ptr = Libderiv->deriv2_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,10,dvrr_stack+21943, dvrr_stack+31296, dvrr_stack+21013);
 tmp = dvrr_stack + 21943;
 target_ptr = Libderiv->deriv2_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,15,dvrr_stack+22093, dvrr_stack+33543, dvrr_stack+2973);
 tmp = dvrr_stack + 22093;
 target_ptr = Libderiv->deriv2_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,21,dvrr_stack+33504, dvrr_stack+14346, dvrr_stack+21113);
 tmp = dvrr_stack + 33504;
 target_ptr = Libderiv->deriv2_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+31296, dvrr_stack+26428, dvrr_stack+4159);
 tmp = dvrr_stack + 31296;
 target_ptr = Libderiv->deriv2_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+26428, dvrr_stack+26203, dvrr_stack+4219);
 tmp = dvrr_stack + 26428;
 target_ptr = Libderiv->deriv2_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+14321, dvrr_stack+26578, dvrr_stack+4309);
 tmp = dvrr_stack + 14321;
 target_ptr = Libderiv->deriv2_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+26578, dvrr_stack+1527, dvrr_stack+17104);
 tmp = dvrr_stack + 26578;
 target_ptr = Libderiv->deriv2_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+14531, dvrr_stack+14787, dvrr_stack+22338);
 tmp = dvrr_stack + 14531;
 target_ptr = Libderiv->deriv2_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+26113, dvrr_stack+33858, dvrr_stack+17204);
 tmp = dvrr_stack + 26113;
 target_ptr = Libderiv->deriv2_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+26728, dvrr_stack+12321, dvrr_stack+15102);
 tmp = dvrr_stack + 26728;
 target_ptr = Libderiv->deriv2_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+17094, dvrr_stack+5699, dvrr_stack+18784);
 tmp = dvrr_stack + 17094;
 target_ptr = Libderiv->deriv2_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+17244, dvrr_stack+11157, dvrr_stack+18874);
 tmp = dvrr_stack + 17244;
 target_ptr = Libderiv->deriv2_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+14756, dvrr_stack+28958, dvrr_stack+2073);
 tmp = dvrr_stack + 14756;
 target_ptr = Libderiv->deriv2_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+28938, dvrr_stack+30666, dvrr_stack+20773);
 tmp = dvrr_stack + 28938;
 target_ptr = Libderiv->deriv2_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+30666, dvrr_stack+32598, dvrr_stack+1023);
 tmp = dvrr_stack + 30666;
 target_ptr = Libderiv->deriv2_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+20773, dvrr_stack+20523, dvrr_stack+35314);
 tmp = dvrr_stack + 20773;
 target_ptr = Libderiv->deriv2_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+20523, dvrr_stack+20298, dvrr_stack+35374);
 tmp = dvrr_stack + 20523;
 target_ptr = Libderiv->deriv2_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+20873, dvrr_stack+15634, dvrr_stack+35464);
 tmp = dvrr_stack + 20873;
 target_ptr = Libderiv->deriv2_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+21083, dvrr_stack+35761, dvrr_stack+19618);
 tmp = dvrr_stack + 21083;
 target_ptr = Libderiv->deriv2_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+20288, dvrr_stack+30981, dvrr_stack+19468);
 tmp = dvrr_stack + 20288;
 target_ptr = Libderiv->deriv2_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+30981, dvrr_stack+35971, dvrr_stack+19718);
 tmp = dvrr_stack + 30981;
 target_ptr = Libderiv->deriv2_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+31396, dvrr_stack+4719, dvrr_stack+20713);
 tmp = dvrr_stack + 31396;
 target_ptr = Libderiv->deriv2_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+14906, dvrr_stack+8097, dvrr_stack+21603);
 tmp = dvrr_stack + 14906;
 target_ptr = Libderiv->deriv2_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+8067, dvrr_stack+2173, dvrr_stack+36412);
 tmp = dvrr_stack + 8067;
 target_ptr = Libderiv->deriv2_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+8277, dvrr_stack+36538, dvrr_stack+11472);
 tmp = dvrr_stack + 8277;
 target_ptr = Libderiv->deriv2_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+2063, dvrr_stack+33039, dvrr_stack+16009);
 tmp = dvrr_stack + 2063;
 target_ptr = Libderiv->deriv2_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+999, dvrr_stack+36748, dvrr_stack+7607);
 tmp = dvrr_stack + 999;
 target_ptr = Libderiv->deriv2_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+20673, dvrr_stack+7328, dvrr_stack+37189);
 tmp = dvrr_stack + 20673;
 target_ptr = Libderiv->deriv2_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+7328, dvrr_stack+7103, dvrr_stack+37249);
 tmp = dvrr_stack + 7328;
 target_ptr = Libderiv->deriv2_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+7478, dvrr_stack+28223, dvrr_stack+37339);
 tmp = dvrr_stack + 7478;
 target_ptr = Libderiv->deriv2_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+28193, dvrr_stack+37465, dvrr_stack+9597);
 tmp = dvrr_stack + 28193;
 target_ptr = Libderiv->deriv2_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+1449, dvrr_stack+37675, dvrr_stack+4869);
 tmp = dvrr_stack + 1449;
 target_ptr = Libderiv->deriv2_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+32586, dvrr_stack+3711, dvrr_stack+27733);
 tmp = dvrr_stack + 32586;
 target_ptr = Libderiv->deriv2_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+9597, dvrr_stack+33354, dvrr_stack+35590);
 tmp = dvrr_stack + 9597;
 target_ptr = Libderiv->deriv2_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+28343, dvrr_stack+19000, dvrr_stack+35650);
 tmp = dvrr_stack + 28343;
 target_ptr = Libderiv->deriv2_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+27653, dvrr_stack+34299, dvrr_stack+29168);
 tmp = dvrr_stack + 27653;
 target_ptr = Libderiv->deriv2_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+2288, dvrr_stack+34714, dvrr_stack+34614);
 tmp = dvrr_stack + 2288;
 target_ptr = Libderiv->deriv2_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+32901, dvrr_stack+29886, dvrr_stack+10337);
 tmp = dvrr_stack + 32901;
 target_ptr = Libderiv->deriv2_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+33126, dvrr_stack+39380, dvrr_stack+10647);
 tmp = dvrr_stack + 33126;
 target_ptr = Libderiv->deriv2_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+7688, dvrr_stack+30201, dvrr_stack+39821);
 tmp = dvrr_stack + 7688;
 target_ptr = Libderiv->deriv2_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+39315, dvrr_stack+39881, dvrr_stack+479);
 tmp = dvrr_stack + 39315;
 target_ptr = Libderiv->deriv2_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+39465, dvrr_stack+38270, dvrr_stack+1323);
 tmp = dvrr_stack + 39465;
 target_ptr = Libderiv->deriv2_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+38240, dvrr_stack+38685, dvrr_stack+38585);
 tmp = dvrr_stack + 38240;
 target_ptr = Libderiv->deriv2_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+38390, dvrr_stack+31656, dvrr_stack+31506);
 tmp = dvrr_stack + 38390;
 target_ptr = Libderiv->deriv2_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+31496, dvrr_stack+40724, dvrr_stack+34924);
 tmp = dvrr_stack + 31496;
 target_ptr = Libderiv->deriv2_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+469, dvrr_stack+31971, dvrr_stack+38895);
 tmp = dvrr_stack + 469;
 target_ptr = Libderiv->deriv2_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+31811, dvrr_stack+41165, dvrr_stack+379);
 tmp = dvrr_stack + 31811;
 target_ptr = Libderiv->deriv2_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+38615, dvrr_stack+30351, dvrr_stack+873);
 tmp = dvrr_stack + 38615;
 target_ptr = Libderiv->deriv2_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+31961, dvrr_stack+569, dvrr_stack+1863);
 tmp = dvrr_stack + 31961;
 target_ptr = Libderiv->deriv2_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+569, dvrr_stack+32121, dvrr_stack+37990);
 tmp = dvrr_stack + 569;
 target_ptr = Libderiv->deriv2_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+32111, dvrr_stack+41390, dvrr_stack+38955);
 tmp = dvrr_stack + 32111;
 target_ptr = Libderiv->deriv2_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];


}

