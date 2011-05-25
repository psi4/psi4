#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (gd|gg) integrals */

void d1vrr_order_gdgg(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+4);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+21, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+24, dvrr_stack+21, dvrr_stack+3, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+30, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+33, dvrr_stack+0, dvrr_stack+30, Data->F+4, Data->F+5, NULL);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+39, dvrr_stack+15, dvrr_stack+33, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+57, dvrr_stack+24, dvrr_stack+15, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+75, dvrr_stack+57, dvrr_stack+39, dvrr_stack+24, dvrr_stack+15, dvrr_stack+6);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+15, dvrr_stack+33, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+121, dvrr_stack+24, dvrr_stack+15, dvrr_stack+21, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+131, dvrr_stack+121, dvrr_stack+111, NULL, NULL, dvrr_stack+15);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+161, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+164, dvrr_stack+161, dvrr_stack+21, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+170, dvrr_stack+164, dvrr_stack+24, dvrr_stack+161, dvrr_stack+21, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+180, dvrr_stack+170, dvrr_stack+121, NULL, NULL, dvrr_stack+24);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+210, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+213, dvrr_stack+30, dvrr_stack+210, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+219, dvrr_stack+33, dvrr_stack+213, dvrr_stack+0, dvrr_stack+30, NULL);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+229, dvrr_stack+111, dvrr_stack+219, NULL, NULL, dvrr_stack+33);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+259, dvrr_stack+131, dvrr_stack+229, dvrr_stack+121, dvrr_stack+111, dvrr_stack+39);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+319, dvrr_stack+180, dvrr_stack+131, dvrr_stack+170, dvrr_stack+121, dvrr_stack+57);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+379, dvrr_stack+319, dvrr_stack+259, dvrr_stack+180, dvrr_stack+131, dvrr_stack+75);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+479, dvrr_stack+121, dvrr_stack+111, dvrr_stack+24, dvrr_stack+15, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+494, dvrr_stack+170, dvrr_stack+121, dvrr_stack+164, dvrr_stack+24, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+509, dvrr_stack+111, dvrr_stack+219, dvrr_stack+15, dvrr_stack+33, NULL);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+524, dvrr_stack+479, dvrr_stack+509, NULL, NULL, dvrr_stack+111);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+569, dvrr_stack+494, dvrr_stack+479, NULL, NULL, dvrr_stack+121);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+614, dvrr_stack+569, dvrr_stack+524, dvrr_stack+494, dvrr_stack+479, dvrr_stack+131);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+704, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+707, dvrr_stack+704, dvrr_stack+161, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+713, dvrr_stack+707, dvrr_stack+164, dvrr_stack+704, dvrr_stack+161, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+723, dvrr_stack+713, dvrr_stack+170, dvrr_stack+707, dvrr_stack+164, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+738, dvrr_stack+723, dvrr_stack+494, NULL, NULL, dvrr_stack+170);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+783, dvrr_stack+738, dvrr_stack+569, dvrr_stack+723, dvrr_stack+494, dvrr_stack+180);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+161, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+704, dvrr_stack+210, dvrr_stack+161, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+873, dvrr_stack+213, dvrr_stack+704, dvrr_stack+30, dvrr_stack+210, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+883, dvrr_stack+219, dvrr_stack+873, dvrr_stack+33, dvrr_stack+213, NULL);

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+898, dvrr_stack+509, dvrr_stack+883, NULL, NULL, dvrr_stack+219);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+943, dvrr_stack+524, dvrr_stack+898, dvrr_stack+479, dvrr_stack+509, dvrr_stack+229);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1033, dvrr_stack+614, dvrr_stack+943, dvrr_stack+569, dvrr_stack+524, dvrr_stack+259);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1183, dvrr_stack+783, dvrr_stack+614, dvrr_stack+738, dvrr_stack+569, dvrr_stack+319);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1333, dvrr_stack+1183, dvrr_stack+1033, dvrr_stack+783, dvrr_stack+614, dvrr_stack+379);
 tmp = dvrr_stack + 1333;
 target_ptr = Libderiv->dvrr_classes[4][4];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+738, dvrr_stack+479, dvrr_stack+509, dvrr_stack+121, dvrr_stack+111, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+759, dvrr_stack+494, dvrr_stack+479, dvrr_stack+170, dvrr_stack+121, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+780, dvrr_stack+509, dvrr_stack+883, dvrr_stack+111, dvrr_stack+219, NULL);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+801, dvrr_stack+738, dvrr_stack+780, NULL, NULL, dvrr_stack+509);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1558, dvrr_stack+759, dvrr_stack+738, NULL, NULL, dvrr_stack+479);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1621, dvrr_stack+1558, dvrr_stack+801, dvrr_stack+759, dvrr_stack+738, dvrr_stack+524);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1747, dvrr_stack+723, dvrr_stack+494, dvrr_stack+713, dvrr_stack+170, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+1747, dvrr_stack+759, NULL, NULL, dvrr_stack+494);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1831, dvrr_stack+1768, dvrr_stack+1558, dvrr_stack+1747, dvrr_stack+759, dvrr_stack+569);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+710, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+569, dvrr_stack+161, dvrr_stack+710, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+121, dvrr_stack+704, dvrr_stack+569, dvrr_stack+210, dvrr_stack+161, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+575, dvrr_stack+873, dvrr_stack+121, dvrr_stack+213, dvrr_stack+704, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+590, dvrr_stack+883, dvrr_stack+575, dvrr_stack+219, dvrr_stack+873, NULL);

 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1957, dvrr_stack+780, dvrr_stack+590, NULL, NULL, dvrr_stack+883);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2020, dvrr_stack+801, dvrr_stack+1957, dvrr_stack+738, dvrr_stack+780, dvrr_stack+898);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2146, dvrr_stack+1621, dvrr_stack+2020, dvrr_stack+1558, dvrr_stack+801, dvrr_stack+943);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2356, dvrr_stack+1831, dvrr_stack+1621, dvrr_stack+1768, dvrr_stack+1558, dvrr_stack+614);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2566, dvrr_stack+2356, dvrr_stack+2146, dvrr_stack+1831, dvrr_stack+1621, dvrr_stack+1033);
 tmp = dvrr_stack + 2566;
 target_ptr = Libderiv->dvrr_classes[4][5];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+2881,dvrr_stack+2566,dvrr_stack+1333,15);


 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+738, dvrr_stack+780, dvrr_stack+479, dvrr_stack+509, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1796, dvrr_stack+759, dvrr_stack+738, dvrr_stack+494, dvrr_stack+479, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1824, dvrr_stack+780, dvrr_stack+590, dvrr_stack+509, dvrr_stack+883, NULL);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1852, dvrr_stack+1768, dvrr_stack+1824, NULL, NULL, dvrr_stack+780);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3556, dvrr_stack+1796, dvrr_stack+1768, NULL, NULL, dvrr_stack+738);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3640, dvrr_stack+3556, dvrr_stack+1852, dvrr_stack+1796, dvrr_stack+1768, dvrr_stack+801);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3808, dvrr_stack+1747, dvrr_stack+759, dvrr_stack+723, dvrr_stack+494, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3836, dvrr_stack+3808, dvrr_stack+1796, NULL, NULL, dvrr_stack+759);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3920, dvrr_stack+3836, dvrr_stack+3556, dvrr_stack+3808, dvrr_stack+1796, dvrr_stack+1558);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+611, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1558, dvrr_stack+710, dvrr_stack+611, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1564, dvrr_stack+569, dvrr_stack+1558, dvrr_stack+161, dvrr_stack+710, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+723, dvrr_stack+121, dvrr_stack+1564, dvrr_stack+704, dvrr_stack+569, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1936, dvrr_stack+575, dvrr_stack+723, dvrr_stack+873, dvrr_stack+121, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1574, dvrr_stack+590, dvrr_stack+1936, dvrr_stack+883, dvrr_stack+575, NULL);

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4088, dvrr_stack+1824, dvrr_stack+1574, NULL, NULL, dvrr_stack+590);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4172, dvrr_stack+1852, dvrr_stack+4088, dvrr_stack+1768, dvrr_stack+1824, dvrr_stack+1957);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4340, dvrr_stack+3640, dvrr_stack+4172, dvrr_stack+3556, dvrr_stack+1852, dvrr_stack+2020);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4620, dvrr_stack+3920, dvrr_stack+3640, dvrr_stack+3836, dvrr_stack+3556, dvrr_stack+1621);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4900, dvrr_stack+4620, dvrr_stack+4340, dvrr_stack+3920, dvrr_stack+3640, dvrr_stack+2146);
 tmp = dvrr_stack + 4900;
 target_ptr = Libderiv->dvrr_classes[4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+5320,dvrr_stack+4900,dvrr_stack+2566,15);


 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3836, dvrr_stack+1768, dvrr_stack+1824, dvrr_stack+738, dvrr_stack+780, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3872, dvrr_stack+1796, dvrr_stack+1768, dvrr_stack+759, dvrr_stack+738, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3908, dvrr_stack+1824, dvrr_stack+1574, dvrr_stack+780, dvrr_stack+590, NULL);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3944, dvrr_stack+3836, dvrr_stack+3908, NULL, NULL, dvrr_stack+1824);

 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6265, dvrr_stack+3872, dvrr_stack+3836, NULL, NULL, dvrr_stack+1768);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6373, dvrr_stack+6265, dvrr_stack+3944, dvrr_stack+3872, dvrr_stack+3836, dvrr_stack+1852);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4052, dvrr_stack+3808, dvrr_stack+1796, dvrr_stack+1747, dvrr_stack+759, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6589, dvrr_stack+4052, dvrr_stack+3872, NULL, NULL, dvrr_stack+1796);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6697, dvrr_stack+6589, dvrr_stack+6265, dvrr_stack+4052, dvrr_stack+3872, dvrr_stack+3556);

 /* compute (0 0 | 1 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3556, Data->F+10, Data->F+11, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+3559, dvrr_stack+611, dvrr_stack+3556, Data->F+9, Data->F+10, NULL);

 /* compute (0 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+3565, dvrr_stack+1558, dvrr_stack+3559, dvrr_stack+710, dvrr_stack+611, NULL);

 /* compute (0 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3575, dvrr_stack+1564, dvrr_stack+3565, dvrr_stack+569, dvrr_stack+1558, NULL);

 /* compute (0 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1747, dvrr_stack+723, dvrr_stack+3575, dvrr_stack+121, dvrr_stack+1564, NULL);

 /* compute (0 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3590, dvrr_stack+1936, dvrr_stack+1747, dvrr_stack+575, dvrr_stack+723, NULL);

 /* compute (0 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+738, dvrr_stack+1574, dvrr_stack+3590, dvrr_stack+590, dvrr_stack+1936, NULL);

 /* compute (1 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6913, dvrr_stack+3908, dvrr_stack+738, NULL, NULL, dvrr_stack+1574);

 /* compute (2 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7021, dvrr_stack+3944, dvrr_stack+6913, dvrr_stack+3836, dvrr_stack+3908, dvrr_stack+4088);

 /* compute (3 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7237, dvrr_stack+6373, dvrr_stack+7021, dvrr_stack+6265, dvrr_stack+3944, dvrr_stack+4172);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7597, dvrr_stack+6697, dvrr_stack+6373, dvrr_stack+6589, dvrr_stack+6265, dvrr_stack+3640);

 /* compute (4 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7957, dvrr_stack+7597, dvrr_stack+7237, dvrr_stack+6697, dvrr_stack+6373, dvrr_stack+4340);
 tmp = dvrr_stack + 7957;
 target_ptr = Libderiv->dvrr_classes[4][7];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+8497,dvrr_stack+7957,dvrr_stack+4900,15);


 /* compute (0 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6589, dvrr_stack+3836, dvrr_stack+3908, dvrr_stack+1768, dvrr_stack+1824, NULL);

 /* compute (0 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6634, dvrr_stack+3872, dvrr_stack+3836, dvrr_stack+1796, dvrr_stack+1768, NULL);

 /* compute (0 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6679, dvrr_stack+3908, dvrr_stack+738, dvrr_stack+1824, dvrr_stack+1574, NULL);

 /* compute (1 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6724, dvrr_stack+6589, dvrr_stack+6679, NULL, NULL, dvrr_stack+3908);

 /* compute (1 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9757, dvrr_stack+6634, dvrr_stack+6589, NULL, NULL, dvrr_stack+3836);

 /* compute (2 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9892, dvrr_stack+9757, dvrr_stack+6724, dvrr_stack+6634, dvrr_stack+6589, dvrr_stack+3944);

 /* compute (0 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6859, dvrr_stack+4052, dvrr_stack+3872, dvrr_stack+3808, dvrr_stack+1796, NULL);

 /* compute (1 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10162, dvrr_stack+6859, dvrr_stack+6634, NULL, NULL, dvrr_stack+3872);

 /* compute (2 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10297, dvrr_stack+10162, dvrr_stack+9757, dvrr_stack+6859, dvrr_stack+6634, dvrr_stack+6265);

 /* compute (0 0 | 1 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+710, Data->F+11, Data->F+12, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+774, dvrr_stack+3556, dvrr_stack+710, Data->F+10, Data->F+11, NULL);

 /* compute (0 0 | 3 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+6265, dvrr_stack+3559, dvrr_stack+774, dvrr_stack+611, dvrr_stack+3556, NULL);

 /* compute (0 0 | 4 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6275, dvrr_stack+3565, dvrr_stack+6265, dvrr_stack+1558, dvrr_stack+3559, NULL);

 /* compute (0 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6290, dvrr_stack+3575, dvrr_stack+6275, dvrr_stack+1564, dvrr_stack+3565, NULL);

 /* compute (0 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3808, dvrr_stack+1747, dvrr_stack+6290, dvrr_stack+723, dvrr_stack+3575, NULL);

 /* compute (0 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6311, dvrr_stack+3590, dvrr_stack+3808, dvrr_stack+1936, dvrr_stack+1747, NULL);

 /* compute (0 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+738, dvrr_stack+6311, dvrr_stack+1574, dvrr_stack+3590, NULL);

 /* compute (1 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10567, dvrr_stack+6679, dvrr_stack+1768, NULL, NULL, dvrr_stack+738);

 /* compute (2 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10702, dvrr_stack+6724, dvrr_stack+10567, dvrr_stack+6589, dvrr_stack+6679, dvrr_stack+6913);

 /* compute (3 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10972, dvrr_stack+9892, dvrr_stack+10702, dvrr_stack+9757, dvrr_stack+6724, dvrr_stack+7021);

 /* compute (3 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+11422, dvrr_stack+10297, dvrr_stack+9892, dvrr_stack+10162, dvrr_stack+9757, dvrr_stack+6373);

 /* compute (4 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+11872, dvrr_stack+11422, dvrr_stack+10972, dvrr_stack+10297, dvrr_stack+9892, dvrr_stack+7237);
 tmp = dvrr_stack + 11872;
 target_ptr = Libderiv->dvrr_classes[4][8];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+12547,dvrr_stack+11872,dvrr_stack+7957,15);


 /* compute (0 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10162, dvrr_stack+6589, dvrr_stack+6679, dvrr_stack+3836, dvrr_stack+3908, NULL);

 /* compute (0 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10217, dvrr_stack+6634, dvrr_stack+6589, dvrr_stack+3872, dvrr_stack+3836, NULL);

 /* compute (0 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10272, dvrr_stack+6679, dvrr_stack+1768, dvrr_stack+3908, dvrr_stack+738, NULL);

 /* compute (1 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10327, dvrr_stack+10162, dvrr_stack+10272, NULL, NULL, dvrr_stack+6679);

 /* compute (1 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+14167, dvrr_stack+10217, dvrr_stack+10162, NULL, NULL, dvrr_stack+6589);

 /* compute (2 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+14332, dvrr_stack+14167, dvrr_stack+10327, dvrr_stack+10217, dvrr_stack+10162, dvrr_stack+6724);

 /* compute (0 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10492, dvrr_stack+6859, dvrr_stack+6634, dvrr_stack+4052, dvrr_stack+3872, NULL);

 /* compute (1 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+14662, dvrr_stack+10492, dvrr_stack+10217, NULL, NULL, dvrr_stack+6634);

 /* compute (2 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+14827, dvrr_stack+14662, dvrr_stack+14167, dvrr_stack+10492, dvrr_stack+10217, dvrr_stack+9757);

 /* compute (0 0 | 1 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+611, Data->F+12, Data->F+13, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1558, dvrr_stack+710, dvrr_stack+611, Data->F+11, Data->F+12, NULL);

 /* compute (0 0 | 3 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+9757, dvrr_stack+774, dvrr_stack+1558, dvrr_stack+3556, dvrr_stack+710, NULL);

 /* compute (0 0 | 4 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+9767, dvrr_stack+6265, dvrr_stack+9757, dvrr_stack+3559, dvrr_stack+774, NULL);

 /* compute (0 0 | 5 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+9782, dvrr_stack+6275, dvrr_stack+9767, dvrr_stack+3565, dvrr_stack+6265, NULL);

 /* compute (0 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+9803, dvrr_stack+6290, dvrr_stack+9782, dvrr_stack+3575, dvrr_stack+6275, NULL);

 /* compute (0 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4052, dvrr_stack+3808, dvrr_stack+9803, dvrr_stack+1747, dvrr_stack+6290, NULL);

 /* compute (0 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9831, dvrr_stack+6311, dvrr_stack+4052, dvrr_stack+3590, dvrr_stack+3808, NULL);

 /* compute (0 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10217, dvrr_stack+1768, dvrr_stack+9831, dvrr_stack+738, dvrr_stack+6311, NULL);

 /* compute (1 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+15157, dvrr_stack+10272, dvrr_stack+10217, NULL, NULL, dvrr_stack+1768);

 /* compute (2 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+15322, dvrr_stack+10327, dvrr_stack+15157, dvrr_stack+10162, dvrr_stack+10272, dvrr_stack+10567);

 /* compute (3 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+15652, dvrr_stack+14332, dvrr_stack+15322, dvrr_stack+14167, dvrr_stack+10327, dvrr_stack+10702);

 /* compute (3 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+16202, dvrr_stack+14827, dvrr_stack+14332, dvrr_stack+14662, dvrr_stack+14167, dvrr_stack+9892);

 /* compute (4 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+16752, dvrr_stack+16202, dvrr_stack+15652, dvrr_stack+14827, dvrr_stack+14332, dvrr_stack+10972);

 /* compute (4 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+17577,dvrr_stack+16752,dvrr_stack+11872,15);


 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+14167, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+864, dvrr_stack+0, dvrr_stack+30, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+14170, dvrr_stack+6, dvrr_stack+864, dvrr_stack+3, dvrr_stack+0, dvrr_stack+14167);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+14188, dvrr_stack+33, dvrr_stack+213, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+14206, dvrr_stack+39, dvrr_stack+14188, dvrr_stack+15, dvrr_stack+33, dvrr_stack+864);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+14242, dvrr_stack+75, dvrr_stack+14206, dvrr_stack+57, dvrr_stack+39, dvrr_stack+14170);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+14302, dvrr_stack+219, dvrr_stack+873, NULL, NULL, dvrr_stack+213);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+14662, dvrr_stack+229, dvrr_stack+14302, dvrr_stack+111, dvrr_stack+219, dvrr_stack+14188);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+14722, dvrr_stack+259, dvrr_stack+14662, dvrr_stack+131, dvrr_stack+229, dvrr_stack+14206);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+14822, dvrr_stack+379, dvrr_stack+14722, dvrr_stack+319, dvrr_stack+259, dvrr_stack+14242);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14972, dvrr_stack+883, dvrr_stack+575, NULL, NULL, dvrr_stack+873);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6589, dvrr_stack+898, dvrr_stack+14972, dvrr_stack+509, dvrr_stack+883, dvrr_stack+14302);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+19602, dvrr_stack+943, dvrr_stack+6589, dvrr_stack+524, dvrr_stack+898, dvrr_stack+14662);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+19752, dvrr_stack+1033, dvrr_stack+19602, dvrr_stack+614, dvrr_stack+943, dvrr_stack+14722);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+19977, dvrr_stack+1333, dvrr_stack+19752, dvrr_stack+1183, dvrr_stack+1033, dvrr_stack+14822);
 tmp = dvrr_stack + 19977;
 target_ptr = Libderiv->dvrr_classes[5][4];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+614, dvrr_stack+590, dvrr_stack+1936, NULL, NULL, dvrr_stack+575);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15017, dvrr_stack+1957, dvrr_stack+614, dvrr_stack+780, dvrr_stack+590, dvrr_stack+14972);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+20292, dvrr_stack+2020, dvrr_stack+15017, dvrr_stack+801, dvrr_stack+1957, dvrr_stack+6589);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+20502, dvrr_stack+2146, dvrr_stack+20292, dvrr_stack+1621, dvrr_stack+2020, dvrr_stack+19602);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+20817, dvrr_stack+2566, dvrr_stack+20502, dvrr_stack+2356, dvrr_stack+2146, dvrr_stack+19752);
 tmp = dvrr_stack + 20817;
 target_ptr = Libderiv->dvrr_classes[5][5];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+21258,dvrr_stack+20817,dvrr_stack+19977,21);


 /* compute (1 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+780, dvrr_stack+1574, dvrr_stack+3590, NULL, NULL, dvrr_stack+1936);

 /* compute (2 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+22203, dvrr_stack+4088, dvrr_stack+780, dvrr_stack+1824, dvrr_stack+1574, dvrr_stack+614);

 /* compute (3 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+22371, dvrr_stack+4172, dvrr_stack+22203, dvrr_stack+1852, dvrr_stack+4088, dvrr_stack+15017);

 /* compute (4 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+22651, dvrr_stack+4340, dvrr_stack+22371, dvrr_stack+3640, dvrr_stack+4172, dvrr_stack+20292);

 /* compute (5 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+23071, dvrr_stack+4900, dvrr_stack+22651, dvrr_stack+4620, dvrr_stack+4340, dvrr_stack+20502);
 tmp = dvrr_stack + 23071;
 target_ptr = Libderiv->dvrr_classes[5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+23659,dvrr_stack+23071,dvrr_stack+20817,21);


 /* compute (1 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1813, dvrr_stack+738, dvrr_stack+6311, NULL, NULL, dvrr_stack+3590);

 /* compute (2 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+24982, dvrr_stack+6913, dvrr_stack+1813, dvrr_stack+3908, dvrr_stack+738, dvrr_stack+780);

 /* compute (3 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+25198, dvrr_stack+7021, dvrr_stack+24982, dvrr_stack+3944, dvrr_stack+6913, dvrr_stack+22203);

 /* compute (4 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+25558, dvrr_stack+7237, dvrr_stack+25198, dvrr_stack+6373, dvrr_stack+7021, dvrr_stack+22371);

 /* compute (5 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+26098, dvrr_stack+7957, dvrr_stack+25558, dvrr_stack+7597, dvrr_stack+7237, dvrr_stack+22651);
 tmp = dvrr_stack + 26098;
 target_ptr = Libderiv->dvrr_classes[5][7];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+26854,dvrr_stack+26098,dvrr_stack+23071,21);


 /* compute (1 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3836, dvrr_stack+1768, dvrr_stack+9831, NULL, NULL, dvrr_stack+6311);

 /* compute (2 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+28618, dvrr_stack+10567, dvrr_stack+3836, dvrr_stack+6679, dvrr_stack+1768, dvrr_stack+1813);

 /* compute (3 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+28888, dvrr_stack+10702, dvrr_stack+28618, dvrr_stack+6724, dvrr_stack+10567, dvrr_stack+24982);

 /* compute (4 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+29338, dvrr_stack+10972, dvrr_stack+28888, dvrr_stack+9892, dvrr_stack+10702, dvrr_stack+25198);

 /* compute (5 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+30013, dvrr_stack+11872, dvrr_stack+29338, dvrr_stack+11422, dvrr_stack+10972, dvrr_stack+25558);
 tmp = dvrr_stack + 30013;
 target_ptr = Libderiv->dvrr_classes[5][8];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+30958,dvrr_stack+30013,dvrr_stack+26098,21);


 /* compute (0 0 | 1 0) m=13 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6679, Data->F+13, Data->F+14, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+611, dvrr_stack+6679, Data->F+12, Data->F+13, NULL);

 /* compute (0 0 | 3 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+1558, dvrr_stack+15, dvrr_stack+710, dvrr_stack+611, NULL);

 /* compute (0 0 | 4 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1921, dvrr_stack+9757, dvrr_stack+111, dvrr_stack+774, dvrr_stack+1558, NULL);

 /* compute (0 0 | 5 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6682, dvrr_stack+9767, dvrr_stack+1921, dvrr_stack+6265, dvrr_stack+9757, NULL);

 /* compute (0 0 | 6 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6703, dvrr_stack+9782, dvrr_stack+6682, dvrr_stack+6275, dvrr_stack+9767, NULL);

 /* compute (0 0 | 7 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6731, dvrr_stack+9803, dvrr_stack+6703, dvrr_stack+6290, dvrr_stack+9782, NULL);

 /* compute (0 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6767, dvrr_stack+4052, dvrr_stack+6731, dvrr_stack+3808, dvrr_stack+9803, NULL);

 /* compute (0 0 | 9 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+6812, dvrr_stack+9831, dvrr_stack+6767, dvrr_stack+6311, dvrr_stack+4052, NULL);

 /* compute (1 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+9876, dvrr_stack+10217, dvrr_stack+6812, NULL, NULL, dvrr_stack+9831);

 /* compute (2 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+33226, dvrr_stack+15157, dvrr_stack+9876, dvrr_stack+10272, dvrr_stack+10217, dvrr_stack+3836);

 /* compute (3 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+33556, dvrr_stack+15322, dvrr_stack+33226, dvrr_stack+10327, dvrr_stack+15157, dvrr_stack+28618);

 /* compute (4 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+34106, dvrr_stack+15652, dvrr_stack+33556, dvrr_stack+14332, dvrr_stack+15322, dvrr_stack+28888);

 /* compute (5 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+34931, dvrr_stack+16752, dvrr_stack+34106, dvrr_stack+16202, dvrr_stack+15652, dvrr_stack+29338);

 /* compute (5 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+36086,dvrr_stack+34931,dvrr_stack+30013,21);


 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+710, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+774, dvrr_stack+14167, dvrr_stack+710, Data->F+4, Data->F+5, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+16202, dvrr_stack+30, dvrr_stack+210, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+16211, dvrr_stack+864, dvrr_stack+16202, dvrr_stack+0, dvrr_stack+30, dvrr_stack+710);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+131, dvrr_stack+14170, dvrr_stack+16211, dvrr_stack+6, dvrr_stack+864, dvrr_stack+774);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+16229, dvrr_stack+213, dvrr_stack+704, NULL, NULL, dvrr_stack+210);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+16247, dvrr_stack+14188, dvrr_stack+16229, dvrr_stack+33, dvrr_stack+213, dvrr_stack+16202);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+16283, dvrr_stack+14206, dvrr_stack+16247, dvrr_stack+39, dvrr_stack+14188, dvrr_stack+16211);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+479, dvrr_stack+14242, dvrr_stack+16283, dvrr_stack+75, dvrr_stack+14206, dvrr_stack+131);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+16343, dvrr_stack+873, dvrr_stack+121, NULL, NULL, dvrr_stack+704);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+16373, dvrr_stack+14302, dvrr_stack+16343, dvrr_stack+219, dvrr_stack+873, dvrr_stack+16229);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+16433, dvrr_stack+14662, dvrr_stack+16373, dvrr_stack+229, dvrr_stack+14302, dvrr_stack+16247);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+16533, dvrr_stack+14722, dvrr_stack+16433, dvrr_stack+259, dvrr_stack+14662, dvrr_stack+16283);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+14332, dvrr_stack+14822, dvrr_stack+16533, dvrr_stack+379, dvrr_stack+14722, dvrr_stack+479);

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+219, dvrr_stack+575, dvrr_stack+723, NULL, NULL, dvrr_stack+121);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14542, dvrr_stack+14972, dvrr_stack+219, dvrr_stack+883, dvrr_stack+575, dvrr_stack+16343);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10272, dvrr_stack+6589, dvrr_stack+14542, dvrr_stack+898, dvrr_stack+14972, dvrr_stack+16373);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6347, dvrr_stack+19602, dvrr_stack+10272, dvrr_stack+943, dvrr_stack+6589, dvrr_stack+16433);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+38921, dvrr_stack+19752, dvrr_stack+6347, dvrr_stack+1033, dvrr_stack+19602, dvrr_stack+16533);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+39236, dvrr_stack+19977, dvrr_stack+38921, dvrr_stack+1333, dvrr_stack+19752, dvrr_stack+14332);
 tmp = dvrr_stack + 39236;
 target_ptr = Libderiv->dvrr_classes[6][4];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+883, dvrr_stack+1936, dvrr_stack+1747, NULL, NULL, dvrr_stack+723);

 /* compute (2 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+946, dvrr_stack+614, dvrr_stack+883, dvrr_stack+590, dvrr_stack+1936, dvrr_stack+219);

 /* compute (3 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+39656, dvrr_stack+15017, dvrr_stack+946, dvrr_stack+1957, dvrr_stack+614, dvrr_stack+14542);

 /* compute (4 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+39866, dvrr_stack+20292, dvrr_stack+39656, dvrr_stack+2020, dvrr_stack+15017, dvrr_stack+10272);

 /* compute (5 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40181, dvrr_stack+20502, dvrr_stack+39866, dvrr_stack+2146, dvrr_stack+20292, dvrr_stack+6347);

 /* compute (6 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40622, dvrr_stack+20817, dvrr_stack+40181, dvrr_stack+2566, dvrr_stack+20502, dvrr_stack+38921);
 tmp = dvrr_stack + 40622;
 target_ptr = Libderiv->dvrr_classes[6][5];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+41210,dvrr_stack+40622,dvrr_stack+39236,28);


 /* compute (1 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1957, dvrr_stack+3590, dvrr_stack+3808, NULL, NULL, dvrr_stack+1747);

 /* compute (2 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2041, dvrr_stack+780, dvrr_stack+1957, dvrr_stack+1574, dvrr_stack+3590, dvrr_stack+883);

 /* compute (3 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+42470, dvrr_stack+22203, dvrr_stack+2041, dvrr_stack+4088, dvrr_stack+780, dvrr_stack+946);

 /* compute (4 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+42750, dvrr_stack+22371, dvrr_stack+42470, dvrr_stack+4172, dvrr_stack+22203, dvrr_stack+39656);

 /* compute (5 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+43170, dvrr_stack+22651, dvrr_stack+42750, dvrr_stack+4340, dvrr_stack+22371, dvrr_stack+39866);

 /* compute (6 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+43758, dvrr_stack+23071, dvrr_stack+43170, dvrr_stack+4900, dvrr_stack+22651, dvrr_stack+40181);
 tmp = dvrr_stack + 43758;
 target_ptr = Libderiv->dvrr_classes[6][6];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+44542,dvrr_stack+43758,dvrr_stack+40622,28);


 /* compute (1 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4088, dvrr_stack+6311, dvrr_stack+4052, NULL, NULL, dvrr_stack+3808);

 /* compute (2 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4196, dvrr_stack+1813, dvrr_stack+4088, dvrr_stack+738, dvrr_stack+6311, dvrr_stack+1957);

 /* compute (3 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+46306, dvrr_stack+24982, dvrr_stack+4196, dvrr_stack+6913, dvrr_stack+1813, dvrr_stack+2041);

 /* compute (4 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+46666, dvrr_stack+25198, dvrr_stack+46306, dvrr_stack+7021, dvrr_stack+24982, dvrr_stack+42470);

 /* compute (5 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+47206, dvrr_stack+25558, dvrr_stack+46666, dvrr_stack+7237, dvrr_stack+25198, dvrr_stack+42750);

 /* compute (6 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+47962, dvrr_stack+26098, dvrr_stack+47206, dvrr_stack+7957, dvrr_stack+25558, dvrr_stack+43170);
 tmp = dvrr_stack + 47962;
 target_ptr = Libderiv->dvrr_classes[6][7];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+48970,dvrr_stack+47962,dvrr_stack+43758,28);


 /* compute (1 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+4412, dvrr_stack+9831, dvrr_stack+6767, NULL, NULL, dvrr_stack+4052);

 /* compute (2 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6867, dvrr_stack+3836, dvrr_stack+4412, dvrr_stack+1768, dvrr_stack+9831, dvrr_stack+4088);

 /* compute (3 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+7137, dvrr_stack+28618, dvrr_stack+6867, dvrr_stack+10567, dvrr_stack+3836, dvrr_stack+4196);

 /* compute (4 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+51322, dvrr_stack+28888, dvrr_stack+7137, dvrr_stack+10702, dvrr_stack+28618, dvrr_stack+46306);

 /* compute (5 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+51997, dvrr_stack+29338, dvrr_stack+51322, dvrr_stack+10972, dvrr_stack+28888, dvrr_stack+46666);

 /* compute (6 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+52942, dvrr_stack+30013, dvrr_stack+51997, dvrr_stack+11872, dvrr_stack+29338, dvrr_stack+47206);

 /* compute (6 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+54202,dvrr_stack+52942,dvrr_stack+47962,28);


 /* compute (0 0 | 1 0) m=14 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+14, Data->F+15, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=13 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1768, dvrr_stack+6679, dvrr_stack+0, Data->F+13, Data->F+14, NULL);

 /* compute (0 0 | 3 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+7587, dvrr_stack+15, dvrr_stack+1768, dvrr_stack+611, dvrr_stack+6679, NULL);

 /* compute (0 0 | 4 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+111, dvrr_stack+7587, dvrr_stack+1558, dvrr_stack+15, NULL);

 /* compute (0 0 | 5 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1783, dvrr_stack+1921, dvrr_stack+1768, dvrr_stack+9757, dvrr_stack+111, NULL);

 /* compute (0 0 | 6 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+738, dvrr_stack+6682, dvrr_stack+1783, dvrr_stack+9767, dvrr_stack+1921, NULL);

 /* compute (0 0 | 7 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+6703, dvrr_stack+738, dvrr_stack+9782, dvrr_stack+6682, NULL);

 /* compute (0 0 | 8 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9757, dvrr_stack+6731, dvrr_stack+1768, dvrr_stack+9803, dvrr_stack+6703, NULL);

 /* compute (0 0 | 9 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+264, dvrr_stack+6767, dvrr_stack+9757, dvrr_stack+4052, dvrr_stack+6731, NULL);

 /* compute (1 0 | 9 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+1574, dvrr_stack+6812, dvrr_stack+264, NULL, NULL, dvrr_stack+6767);

 /* compute (2 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10422, dvrr_stack+9876, dvrr_stack+1574, dvrr_stack+10217, dvrr_stack+6812, dvrr_stack+4412);

 /* compute (3 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10752, dvrr_stack+33226, dvrr_stack+10422, dvrr_stack+15157, dvrr_stack+9876, dvrr_stack+6867);

 /* compute (4 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+57226, dvrr_stack+33556, dvrr_stack+10752, dvrr_stack+15322, dvrr_stack+33226, dvrr_stack+7137);

 /* compute (5 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+58051, dvrr_stack+34106, dvrr_stack+57226, dvrr_stack+15652, dvrr_stack+33556, dvrr_stack+51322);

 /* compute (6 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+59206, dvrr_stack+34931, dvrr_stack+58051, dvrr_stack+16752, dvrr_stack+34106, dvrr_stack+51997);

 /* compute (6 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+60746,dvrr_stack+59206,dvrr_stack+52942,28);


 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+57226, dvrr_stack+21, dvrr_stack+3, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+57235, dvrr_stack+57226, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+0);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+57253, dvrr_stack+164, dvrr_stack+24, NULL, NULL, dvrr_stack+21);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+738, dvrr_stack+57253, dvrr_stack+57, dvrr_stack+164, dvrr_stack+24, dvrr_stack+57226);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+57271, dvrr_stack+738, dvrr_stack+75, dvrr_stack+57253, dvrr_stack+57, dvrr_stack+57235);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+14632, dvrr_stack+713, dvrr_stack+170, NULL, NULL, dvrr_stack+164);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+57331, dvrr_stack+14632, dvrr_stack+180, dvrr_stack+713, dvrr_stack+170, dvrr_stack+57253);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+57391, dvrr_stack+57331, dvrr_stack+319, dvrr_stack+14632, dvrr_stack+180, dvrr_stack+738);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+57491, dvrr_stack+57391, dvrr_stack+379, dvrr_stack+57331, dvrr_stack+319, dvrr_stack+57271);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+1558, dvrr_stack+0, dvrr_stack+14167, Data->F+3, Data->F+4, NULL);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+14632, dvrr_stack+57235, dvrr_stack+14170, dvrr_stack+57226, dvrr_stack+6, dvrr_stack+1558);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+57641, dvrr_stack+57271, dvrr_stack+14242, dvrr_stack+738, dvrr_stack+75, dvrr_stack+14632);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+57731, dvrr_stack+57491, dvrr_stack+14822, dvrr_stack+57391, dvrr_stack+379, dvrr_stack+57641);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+713, dvrr_stack+1558, dvrr_stack+774, dvrr_stack+0, dvrr_stack+14167, NULL);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+14632, dvrr_stack+131, dvrr_stack+57235, dvrr_stack+14170, dvrr_stack+713);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+57331, dvrr_stack+57641, dvrr_stack+479, dvrr_stack+57271, dvrr_stack+14242, dvrr_stack+1768);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+57941, dvrr_stack+57731, dvrr_stack+14332, dvrr_stack+57491, dvrr_stack+14822, dvrr_stack+57331);

 /* compute (1 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+1768, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+1558, dvrr_stack+710, dvrr_stack+1768, Data->F+5, Data->F+6, NULL);

 /* compute (3 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+713, dvrr_stack+774, dvrr_stack+1558, dvrr_stack+14167, dvrr_stack+710, NULL);

 /* compute (1 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1771, dvrr_stack+210, dvrr_stack+161, NULL, NULL, Data->F+7);

 /* compute (2 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+1780, dvrr_stack+16202, dvrr_stack+1771, dvrr_stack+30, dvrr_stack+210, dvrr_stack+1768);

 /* compute (3 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+14632, dvrr_stack+16211, dvrr_stack+1780, dvrr_stack+864, dvrr_stack+16202, dvrr_stack+1558);

 /* compute (4 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+57641, dvrr_stack+131, dvrr_stack+14632, dvrr_stack+14170, dvrr_stack+16211, dvrr_stack+713);

 /* compute (1 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+131, dvrr_stack+704, dvrr_stack+569, NULL, NULL, dvrr_stack+161);

 /* compute (2 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+149, dvrr_stack+16229, dvrr_stack+131, dvrr_stack+213, dvrr_stack+704, dvrr_stack+1771);

 /* compute (3 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+0, dvrr_stack+16247, dvrr_stack+149, dvrr_stack+14188, dvrr_stack+16229, dvrr_stack+1780);

 /* compute (4 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+57226, dvrr_stack+16283, dvrr_stack+0, dvrr_stack+14206, dvrr_stack+16247, dvrr_stack+14632);

 /* compute (5 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+57316, dvrr_stack+479, dvrr_stack+57226, dvrr_stack+14242, dvrr_stack+16283, dvrr_stack+57641);

 /* compute (1 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+14632, dvrr_stack+121, dvrr_stack+1564, NULL, NULL, dvrr_stack+569);

 /* compute (2 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+57641, dvrr_stack+16343, dvrr_stack+14632, dvrr_stack+873, dvrr_stack+121, dvrr_stack+131);

 /* compute (3 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+14167, dvrr_stack+16373, dvrr_stack+57641, dvrr_stack+14302, dvrr_stack+16343, dvrr_stack+149);

 /* compute (4 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+60, dvrr_stack+16433, dvrr_stack+14167, dvrr_stack+14662, dvrr_stack+16373, dvrr_stack+0);

 /* compute (5 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+58221, dvrr_stack+16533, dvrr_stack+60, dvrr_stack+14722, dvrr_stack+16433, dvrr_stack+57226);

 /* compute (6 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+58431, dvrr_stack+14332, dvrr_stack+58221, dvrr_stack+14822, dvrr_stack+16533, dvrr_stack+57316);

 /* compute (1 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+723, dvrr_stack+3575, NULL, NULL, dvrr_stack+1564);

 /* compute (2 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+57226, dvrr_stack+219, dvrr_stack+1768, dvrr_stack+575, dvrr_stack+723, dvrr_stack+14632);

 /* compute (3 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14632, dvrr_stack+14542, dvrr_stack+57226, dvrr_stack+14972, dvrr_stack+219, dvrr_stack+57641);

 /* compute (4 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14782, dvrr_stack+10272, dvrr_stack+14632, dvrr_stack+6589, dvrr_stack+14542, dvrr_stack+14167);

 /* compute (5 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14167, dvrr_stack+6347, dvrr_stack+14782, dvrr_stack+19602, dvrr_stack+10272, dvrr_stack+60);

 /* compute (6 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+38921, dvrr_stack+14167, dvrr_stack+19752, dvrr_stack+6347, dvrr_stack+58221);

 /* compute (7 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+33226, dvrr_stack+39236, dvrr_stack+0, dvrr_stack+19977, dvrr_stack+38921, dvrr_stack+58431);

 /* compute (1 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+38921, dvrr_stack+1747, dvrr_stack+6290, NULL, NULL, dvrr_stack+3575);

 /* compute (2 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+38984, dvrr_stack+883, dvrr_stack+38921, dvrr_stack+1936, dvrr_stack+1747, dvrr_stack+1768);

 /* compute (3 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+58221, dvrr_stack+946, dvrr_stack+38984, dvrr_stack+614, dvrr_stack+883, dvrr_stack+57226);

 /* compute (4 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+58431, dvrr_stack+39656, dvrr_stack+58221, dvrr_stack+15017, dvrr_stack+946, dvrr_stack+14632);

 /* compute (5 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+58746, dvrr_stack+39866, dvrr_stack+58431, dvrr_stack+20292, dvrr_stack+39656, dvrr_stack+14782);

 /* compute (6 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+14482, dvrr_stack+40181, dvrr_stack+58746, dvrr_stack+20502, dvrr_stack+39866, dvrr_stack+14167);

 /* compute (7 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15070, dvrr_stack+40622, dvrr_stack+14482, dvrr_stack+20817, dvrr_stack+40181, dvrr_stack+0);

 /* compute (1 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+3808, dvrr_stack+9803, NULL, NULL, dvrr_stack+6290);

 /* compute (2 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+84, dvrr_stack+1957, dvrr_stack+0, dvrr_stack+3590, dvrr_stack+3808, dvrr_stack+38921);

 /* compute (3 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3556, dvrr_stack+2041, dvrr_stack+84, dvrr_stack+780, dvrr_stack+1957, dvrr_stack+38984);

 /* compute (4 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+252, dvrr_stack+42470, dvrr_stack+3556, dvrr_stack+22203, dvrr_stack+2041, dvrr_stack+58221);

 /* compute (5 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+39656, dvrr_stack+42750, dvrr_stack+252, dvrr_stack+22371, dvrr_stack+42470, dvrr_stack+58431);

 /* compute (6 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+15826, dvrr_stack+43170, dvrr_stack+39656, dvrr_stack+22651, dvrr_stack+42750, dvrr_stack+58746);

 /* compute (7 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+33766, dvrr_stack+43758, dvrr_stack+15826, dvrr_stack+23071, dvrr_stack+43170, dvrr_stack+14482);

 /* compute (1 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+42470, dvrr_stack+4052, dvrr_stack+6731, NULL, NULL, dvrr_stack+9803);

 /* compute (2 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+42578, dvrr_stack+4088, dvrr_stack+42470, dvrr_stack+6311, dvrr_stack+4052, dvrr_stack+0);

 /* compute (3 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+42794, dvrr_stack+4196, dvrr_stack+42578, dvrr_stack+1813, dvrr_stack+4088, dvrr_stack+84);

 /* compute (4 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+43154, dvrr_stack+46306, dvrr_stack+42794, dvrr_stack+24982, dvrr_stack+4196, dvrr_stack+3556);

 /* compute (5 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+58221, dvrr_stack+46666, dvrr_stack+43154, dvrr_stack+25198, dvrr_stack+46306, dvrr_stack+252);

 /* compute (6 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+47206, dvrr_stack+58221, dvrr_stack+25558, dvrr_stack+46666, dvrr_stack+39656);

 /* compute (7 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+9876, dvrr_stack+47962, dvrr_stack+0, dvrr_stack+26098, dvrr_stack+47206, dvrr_stack+15826);

 /* compute (1 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+15826, dvrr_stack+6767, dvrr_stack+9757, NULL, NULL, dvrr_stack+6731);

 /* compute (2 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+15961, dvrr_stack+4412, dvrr_stack+15826, dvrr_stack+9831, dvrr_stack+6767, dvrr_stack+42470);

 /* compute (3 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+16231, dvrr_stack+6867, dvrr_stack+15961, dvrr_stack+3836, dvrr_stack+4412, dvrr_stack+42578);

 /* compute (4 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+39656, dvrr_stack+7137, dvrr_stack+16231, dvrr_stack+28618, dvrr_stack+6867, dvrr_stack+42794);

 /* compute (5 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+46306, dvrr_stack+51322, dvrr_stack+39656, dvrr_stack+28888, dvrr_stack+7137, dvrr_stack+43154);

 /* compute (6 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+42470, dvrr_stack+51997, dvrr_stack+46306, dvrr_stack+29338, dvrr_stack+51322, dvrr_stack+58221);

 /* compute (7 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+46306, dvrr_stack+52942, dvrr_stack+42470, dvrr_stack+30013, dvrr_stack+51997, dvrr_stack+0);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+0, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+38921, dvrr_stack+5320, NULL);
 tmp = dvrr_stack + 38921;
 target_ptr = Libderiv->deriv_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+225, dvrr_stack+8497, NULL);
 tmp = dvrr_stack + 225;
 target_ptr = Libderiv->deriv_classes[4][6][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,540,dvrr_stack+42470, dvrr_stack+12547, NULL);
 tmp = dvrr_stack + 42470;
 target_ptr = Libderiv->deriv_classes[4][7][11];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,675,dvrr_stack+43010, dvrr_stack+17577, NULL);
 tmp = dvrr_stack + 43010;
 target_ptr = Libderiv->deriv_classes[4][8][11];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+645, dvrr_stack+21258, NULL);
 tmp = dvrr_stack + 645;
 target_ptr = Libderiv->deriv_classes[5][4][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,441,dvrr_stack+58221, dvrr_stack+23659, NULL);
 tmp = dvrr_stack + 58221;
 target_ptr = Libderiv->deriv_classes[5][5][11];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+51322, dvrr_stack+26854, NULL);
 tmp = dvrr_stack + 51322;
 target_ptr = Libderiv->deriv_classes[5][6][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,756,dvrr_stack+51910, dvrr_stack+30958, NULL);
 tmp = dvrr_stack + 51910;
 target_ptr = Libderiv->deriv_classes[5][7][11];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,945,dvrr_stack+39656, dvrr_stack+36086, NULL);
 tmp = dvrr_stack + 39656;
 target_ptr = Libderiv->deriv_classes[5][8][11];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+58662, dvrr_stack+41210, NULL);
 tmp = dvrr_stack + 58662;
 target_ptr = Libderiv->deriv_classes[6][4][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+28618, dvrr_stack+44542, NULL);
 tmp = dvrr_stack + 28618;
 target_ptr = Libderiv->deriv_classes[6][5][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,784,dvrr_stack+29206, dvrr_stack+48970, NULL);
 tmp = dvrr_stack + 29206;
 target_ptr = Libderiv->deriv_classes[6][6][11];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1008,dvrr_stack+3556, dvrr_stack+54202, NULL);
 tmp = dvrr_stack + 3556;
 target_ptr = Libderiv->deriv_classes[6][7][11];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1260,dvrr_stack+6265, dvrr_stack+60746, NULL);
 tmp = dvrr_stack + 6265;
 target_ptr = Libderiv->deriv_classes[6][8][11];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+52666, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 52666;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+15826, dvrr_stack+5320, NULL);
 tmp = dvrr_stack + 15826;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+16141, dvrr_stack+8497, NULL);
 tmp = dvrr_stack + 16141;
 target_ptr = Libderiv->deriv_classes[4][6][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,540,dvrr_stack+24982, dvrr_stack+12547, NULL);
 tmp = dvrr_stack + 24982;
 target_ptr = Libderiv->deriv_classes[4][7][10];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,675,dvrr_stack+22203, dvrr_stack+17577, NULL);
 tmp = dvrr_stack + 22203;
 target_ptr = Libderiv->deriv_classes[4][8][10];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+25522, dvrr_stack+21258, NULL);
 tmp = dvrr_stack + 25522;
 target_ptr = Libderiv->deriv_classes[5][4][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+14167, dvrr_stack+23659, NULL);
 tmp = dvrr_stack + 14167;
 target_ptr = Libderiv->deriv_classes[5][5][10];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+1558, dvrr_stack+26854, NULL);
 tmp = dvrr_stack + 1558;
 target_ptr = Libderiv->deriv_classes[5][6][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,756,dvrr_stack+64526, dvrr_stack+30958, NULL);
 tmp = dvrr_stack + 64526;
 target_ptr = Libderiv->deriv_classes[5][7][10];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,945,dvrr_stack+65282, dvrr_stack+36086, NULL);
 tmp = dvrr_stack + 65282;
 target_ptr = Libderiv->deriv_classes[5][8][10];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+14608, dvrr_stack+41210, NULL);
 tmp = dvrr_stack + 14608;
 target_ptr = Libderiv->deriv_classes[6][4][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+66227, dvrr_stack+44542, NULL);
 tmp = dvrr_stack + 66227;
 target_ptr = Libderiv->deriv_classes[6][5][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,784,dvrr_stack+66815, dvrr_stack+48970, NULL);
 tmp = dvrr_stack + 66815;
 target_ptr = Libderiv->deriv_classes[6][6][10];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1008,dvrr_stack+67599, dvrr_stack+54202, NULL);
 tmp = dvrr_stack + 67599;
 target_ptr = Libderiv->deriv_classes[6][7][10];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1260,dvrr_stack+68607, dvrr_stack+60746, NULL);
 tmp = dvrr_stack + 68607;
 target_ptr = Libderiv->deriv_classes[6][8][10];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+25837, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 25837;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+2881, dvrr_stack+5320, NULL);
 tmp = dvrr_stack + 2881;
 target_ptr = Libderiv->deriv_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+5320, dvrr_stack+8497, NULL);
 tmp = dvrr_stack + 5320;
 target_ptr = Libderiv->deriv_classes[4][6][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,540,dvrr_stack+8497, dvrr_stack+12547, NULL);
 tmp = dvrr_stack + 8497;
 target_ptr = Libderiv->deriv_classes[4][7][9];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,675,dvrr_stack+12547, dvrr_stack+17577, NULL);
 tmp = dvrr_stack + 12547;
 target_ptr = Libderiv->deriv_classes[4][8][9];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+17577, dvrr_stack+21258, NULL);
 tmp = dvrr_stack + 17577;
 target_ptr = Libderiv->deriv_classes[5][4][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+21258, dvrr_stack+23659, NULL);
 tmp = dvrr_stack + 21258;
 target_ptr = Libderiv->deriv_classes[5][5][9];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+23659, dvrr_stack+26854, NULL);
 tmp = dvrr_stack + 23659;
 target_ptr = Libderiv->deriv_classes[5][6][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,756,dvrr_stack+26854, dvrr_stack+30958, NULL);
 tmp = dvrr_stack + 26854;
 target_ptr = Libderiv->deriv_classes[5][7][9];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,945,dvrr_stack+13222, dvrr_stack+36086, NULL);
 tmp = dvrr_stack + 13222;
 target_ptr = Libderiv->deriv_classes[5][8][9];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+36086, dvrr_stack+41210, NULL);
 tmp = dvrr_stack + 36086;
 target_ptr = Libderiv->deriv_classes[6][4][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+41210, dvrr_stack+44542, NULL);
 tmp = dvrr_stack + 41210;
 target_ptr = Libderiv->deriv_classes[6][5][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,784,dvrr_stack+44542, dvrr_stack+48970, NULL);
 tmp = dvrr_stack + 44542;
 target_ptr = Libderiv->deriv_classes[6][6][9];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1008,dvrr_stack+27610, dvrr_stack+54202, NULL);
 tmp = dvrr_stack + 27610;
 target_ptr = Libderiv->deriv_classes[6][7][9];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1260,dvrr_stack+54202, dvrr_stack+60746, NULL);
 tmp = dvrr_stack + 54202;
 target_ptr = Libderiv->deriv_classes[6][8][9];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+60746, dvrr_stack+2566, dvrr_stack+57491);
 tmp = dvrr_stack + 60746;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+60971, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 60971;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+61286, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 61286;
 target_ptr = Libderiv->deriv_classes[4][6][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,15,1,dvrr_stack+61706, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 61706;
 target_ptr = Libderiv->deriv_classes[4][7][8];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,15,1,dvrr_stack+62246, dvrr_stack+16752, dvrr_stack+7957);
 tmp = dvrr_stack + 62246;
 target_ptr = Libderiv->deriv_classes[4][8][8];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+62921, dvrr_stack+20817, dvrr_stack+57731);
 tmp = dvrr_stack + 62921;
 target_ptr = Libderiv->deriv_classes[5][4][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+63236, dvrr_stack+23071, dvrr_stack+19977);
 tmp = dvrr_stack + 63236;
 target_ptr = Libderiv->deriv_classes[5][5][8];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+63677, dvrr_stack+26098, dvrr_stack+20817);
 tmp = dvrr_stack + 63677;
 target_ptr = Libderiv->deriv_classes[5][6][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,21,1,dvrr_stack+55462, dvrr_stack+30013, dvrr_stack+23071);
 tmp = dvrr_stack + 55462;
 target_ptr = Libderiv->deriv_classes[5][7][8];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,21,1,dvrr_stack+56218, dvrr_stack+34931, dvrr_stack+26098);
 tmp = dvrr_stack + 56218;
 target_ptr = Libderiv->deriv_classes[5][8][8];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+48970, dvrr_stack+40622, dvrr_stack+57941);
 tmp = dvrr_stack + 48970;
 target_ptr = Libderiv->deriv_classes[6][4][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+49390, dvrr_stack+43758, dvrr_stack+39236);
 tmp = dvrr_stack + 49390;
 target_ptr = Libderiv->deriv_classes[6][5][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+49978, dvrr_stack+47962, dvrr_stack+40622);
 tmp = dvrr_stack + 49978;
 target_ptr = Libderiv->deriv_classes[6][6][8];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,28,1,dvrr_stack+36506, dvrr_stack+52942, dvrr_stack+43758);
 tmp = dvrr_stack + 36506;
 target_ptr = Libderiv->deriv_classes[6][7][8];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,28,1,dvrr_stack+37514, dvrr_stack+59206, dvrr_stack+47962);
 tmp = dvrr_stack + 37514;
 target_ptr = Libderiv->deriv_classes[6][8][8];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+64265, dvrr_stack+2566, dvrr_stack+57491);
 tmp = dvrr_stack + 64265;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+57163, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 57163;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+50762, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 50762;
 target_ptr = Libderiv->deriv_classes[4][6][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,15,1,dvrr_stack+45326, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 45326;
 target_ptr = Libderiv->deriv_classes[4][7][7];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,15,1,dvrr_stack+30958, dvrr_stack+16752, dvrr_stack+7957);
 tmp = dvrr_stack + 30958;
 target_ptr = Libderiv->deriv_classes[4][8][7];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+45866, dvrr_stack+20817, dvrr_stack+57731);
 tmp = dvrr_stack + 45866;
 target_ptr = Libderiv->deriv_classes[5][4][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+41798, dvrr_stack+23071, dvrr_stack+19977);
 tmp = dvrr_stack + 41798;
 target_ptr = Libderiv->deriv_classes[5][5][7];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+31633, dvrr_stack+26098, dvrr_stack+20817);
 tmp = dvrr_stack + 31633;
 target_ptr = Libderiv->deriv_classes[5][6][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,21,1,dvrr_stack+32221, dvrr_stack+30013, dvrr_stack+23071);
 tmp = dvrr_stack + 32221;
 target_ptr = Libderiv->deriv_classes[5][7][7];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,21,1,dvrr_stack+17892, dvrr_stack+34931, dvrr_stack+26098);
 tmp = dvrr_stack + 17892;
 target_ptr = Libderiv->deriv_classes[5][8][7];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+24247, dvrr_stack+40622, dvrr_stack+57941);
 tmp = dvrr_stack + 24247;
 target_ptr = Libderiv->deriv_classes[6][4][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+18837, dvrr_stack+43758, dvrr_stack+39236);
 tmp = dvrr_stack + 18837;
 target_ptr = Libderiv->deriv_classes[6][5][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+9037, dvrr_stack+47962, dvrr_stack+40622);
 tmp = dvrr_stack + 9037;
 target_ptr = Libderiv->deriv_classes[6][6][7];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,28,1,dvrr_stack+69867, dvrr_stack+52942, dvrr_stack+43758);
 tmp = dvrr_stack + 69867;
 target_ptr = Libderiv->deriv_classes[6][7][7];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,28,1,dvrr_stack+70875, dvrr_stack+59206, dvrr_stack+47962);
 tmp = dvrr_stack + 70875;
 target_ptr = Libderiv->deriv_classes[6][8][7];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+42239, dvrr_stack+2566, dvrr_stack+57491);
 tmp = dvrr_stack + 42239;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+24667, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 24667;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+21699, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 21699;
 target_ptr = Libderiv->deriv_classes[4][6][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,15,1,dvrr_stack+19425, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 19425;
 target_ptr = Libderiv->deriv_classes[4][7][6];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,15,1,dvrr_stack+72135, dvrr_stack+16752, dvrr_stack+7957);
 tmp = dvrr_stack + 72135;
 target_ptr = Libderiv->deriv_classes[4][8][6];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+5740, dvrr_stack+20817, dvrr_stack+57731);
 tmp = dvrr_stack + 5740;
 target_ptr = Libderiv->deriv_classes[5][4][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+57478, dvrr_stack+23071, dvrr_stack+19977);
 tmp = dvrr_stack + 57478;
 target_ptr = Libderiv->deriv_classes[5][5][6];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+16561, dvrr_stack+26098, dvrr_stack+20817);
 tmp = dvrr_stack + 16561;
 target_ptr = Libderiv->deriv_classes[5][6][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,21,1,dvrr_stack+72810, dvrr_stack+30013, dvrr_stack+23071);
 tmp = dvrr_stack + 72810;
 target_ptr = Libderiv->deriv_classes[5][7][6];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,21,1,dvrr_stack+73566, dvrr_stack+34931, dvrr_stack+26098);
 tmp = dvrr_stack + 73566;
 target_ptr = Libderiv->deriv_classes[5][8][6];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+17149, dvrr_stack+40622, dvrr_stack+57941);
 tmp = dvrr_stack + 17149;
 target_ptr = Libderiv->deriv_classes[6][4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+34774, dvrr_stack+43758, dvrr_stack+39236);
 tmp = dvrr_stack + 34774;
 target_ptr = Libderiv->deriv_classes[6][5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+74511, dvrr_stack+47962, dvrr_stack+40622);
 tmp = dvrr_stack + 74511;
 target_ptr = Libderiv->deriv_classes[6][6][6];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,28,1,dvrr_stack+75295, dvrr_stack+52942, dvrr_stack+43758);
 tmp = dvrr_stack + 75295;
 target_ptr = Libderiv->deriv_classes[6][7][6];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,28,1,dvrr_stack+76303, dvrr_stack+59206, dvrr_stack+47962);
 tmp = dvrr_stack + 76303;
 target_ptr = Libderiv->deriv_classes[6][8][6];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+57919, dvrr_stack+19977, dvrr_stack+1183);
 tmp = dvrr_stack + 57919;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+3196, dvrr_stack+20817, dvrr_stack+2356);
 tmp = dvrr_stack + 3196;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+59082, dvrr_stack+23071, dvrr_stack+4620);
 tmp = dvrr_stack + 59082;
 target_ptr = Libderiv->deriv_classes[4][6][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,36,dvrr_stack+59502, dvrr_stack+26098, dvrr_stack+7597);
 tmp = dvrr_stack + 59502;
 target_ptr = Libderiv->deriv_classes[4][7][2];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,45,dvrr_stack+60042, dvrr_stack+30013, dvrr_stack+11422);
 tmp = dvrr_stack + 60042;
 target_ptr = Libderiv->deriv_classes[4][8][2];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+20292, dvrr_stack+39236, dvrr_stack+1333);
 tmp = dvrr_stack + 20292;
 target_ptr = Libderiv->deriv_classes[5][4][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+35362, dvrr_stack+40622, dvrr_stack+2566);
 tmp = dvrr_stack + 35362;
 target_ptr = Libderiv->deriv_classes[5][5][2];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+77563, dvrr_stack+43758, dvrr_stack+4900);
 tmp = dvrr_stack + 77563;
 target_ptr = Libderiv->deriv_classes[5][6][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,36,dvrr_stack+78151, dvrr_stack+47962, dvrr_stack+7957);
 tmp = dvrr_stack + 78151;
 target_ptr = Libderiv->deriv_classes[5][7][2];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,45,dvrr_stack+78907, dvrr_stack+52942, dvrr_stack+11872);
 tmp = dvrr_stack + 78907;
 target_ptr = Libderiv->deriv_classes[5][8][2];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+79852, dvrr_stack+33226, dvrr_stack+19977);
 tmp = dvrr_stack + 79852;
 target_ptr = Libderiv->deriv_classes[6][4][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+80272, dvrr_stack+15070, dvrr_stack+20817);
 tmp = dvrr_stack + 80272;
 target_ptr = Libderiv->deriv_classes[6][5][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+80860, dvrr_stack+33766, dvrr_stack+23071);
 tmp = dvrr_stack + 80860;
 target_ptr = Libderiv->deriv_classes[6][6][2];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,36,dvrr_stack+81644, dvrr_stack+9876, dvrr_stack+26098);
 tmp = dvrr_stack + 81644;
 target_ptr = Libderiv->deriv_classes[6][7][2];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,45,dvrr_stack+82652, dvrr_stack+46306, dvrr_stack+30013);
 tmp = dvrr_stack + 82652;
 target_ptr = Libderiv->deriv_classes[6][8][2];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+32977, dvrr_stack+19977, dvrr_stack+1183);
 tmp = dvrr_stack + 32977;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+83912, dvrr_stack+20817, dvrr_stack+2356);
 tmp = dvrr_stack + 83912;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+84227, dvrr_stack+23071, dvrr_stack+4620);
 tmp = dvrr_stack + 84227;
 target_ptr = Libderiv->deriv_classes[4][6][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,36,dvrr_stack+84647, dvrr_stack+26098, dvrr_stack+7597);
 tmp = dvrr_stack + 84647;
 target_ptr = Libderiv->deriv_classes[4][7][1];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,45,dvrr_stack+85187, dvrr_stack+30013, dvrr_stack+11422);
 tmp = dvrr_stack + 85187;
 target_ptr = Libderiv->deriv_classes[4][8][1];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+85862, dvrr_stack+39236, dvrr_stack+1333);
 tmp = dvrr_stack + 85862;
 target_ptr = Libderiv->deriv_classes[5][4][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+86177, dvrr_stack+40622, dvrr_stack+2566);
 tmp = dvrr_stack + 86177;
 target_ptr = Libderiv->deriv_classes[5][5][1];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+86618, dvrr_stack+43758, dvrr_stack+4900);
 tmp = dvrr_stack + 86618;
 target_ptr = Libderiv->deriv_classes[5][6][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,36,dvrr_stack+87206, dvrr_stack+47962, dvrr_stack+7957);
 tmp = dvrr_stack + 87206;
 target_ptr = Libderiv->deriv_classes[5][7][1];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,45,dvrr_stack+87962, dvrr_stack+52942, dvrr_stack+11872);
 tmp = dvrr_stack + 87962;
 target_ptr = Libderiv->deriv_classes[5][8][1];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+88907, dvrr_stack+33226, dvrr_stack+19977);
 tmp = dvrr_stack + 88907;
 target_ptr = Libderiv->deriv_classes[6][4][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+89327, dvrr_stack+15070, dvrr_stack+20817);
 tmp = dvrr_stack + 89327;
 target_ptr = Libderiv->deriv_classes[6][5][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+89915, dvrr_stack+33766, dvrr_stack+23071);
 tmp = dvrr_stack + 89915;
 target_ptr = Libderiv->deriv_classes[6][6][1];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,36,dvrr_stack+90699, dvrr_stack+9876, dvrr_stack+26098);
 tmp = dvrr_stack + 90699;
 target_ptr = Libderiv->deriv_classes[6][7][1];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,45,dvrr_stack+91707, dvrr_stack+46306, dvrr_stack+30013);
 tmp = dvrr_stack + 91707;
 target_ptr = Libderiv->deriv_classes[6][8][1];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+35803, dvrr_stack+19977, dvrr_stack+1183);
 tmp = dvrr_stack + 35803;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+960, dvrr_stack+20817, dvrr_stack+2356);
 tmp = dvrr_stack + 960;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+2146, dvrr_stack+23071, dvrr_stack+4620);
 tmp = dvrr_stack + 2146;
 target_ptr = Libderiv->deriv_classes[4][6][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,36,dvrr_stack+92967, dvrr_stack+26098, dvrr_stack+7597);
 tmp = dvrr_stack + 92967;
 target_ptr = Libderiv->deriv_classes[4][7][0];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,45,dvrr_stack+93507, dvrr_stack+30013, dvrr_stack+11422);
 tmp = dvrr_stack + 93507;
 target_ptr = Libderiv->deriv_classes[4][8][0];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+4564, dvrr_stack+39236, dvrr_stack+1333);
 tmp = dvrr_stack + 4564;
 target_ptr = Libderiv->deriv_classes[5][4][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+11172, dvrr_stack+40622, dvrr_stack+2566);
 tmp = dvrr_stack + 11172;
 target_ptr = Libderiv->deriv_classes[5][5][0];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+40601, dvrr_stack+43758, dvrr_stack+4900);
 tmp = dvrr_stack + 40601;
 target_ptr = Libderiv->deriv_classes[5][6][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,36,dvrr_stack+43685, dvrr_stack+47962, dvrr_stack+7957);
 tmp = dvrr_stack + 43685;
 target_ptr = Libderiv->deriv_classes[5][7][0];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,45,dvrr_stack+47926, dvrr_stack+52942, dvrr_stack+11872);
 tmp = dvrr_stack + 47926;
 target_ptr = Libderiv->deriv_classes[5][8][0];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+39236, dvrr_stack+33226, dvrr_stack+19977);
 tmp = dvrr_stack + 39236;
 target_ptr = Libderiv->deriv_classes[6][4][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+52891, dvrr_stack+15070, dvrr_stack+20817);
 tmp = dvrr_stack + 52891;
 target_ptr = Libderiv->deriv_classes[6][5][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+15028, dvrr_stack+33766, dvrr_stack+23071);
 tmp = dvrr_stack + 15028;
 target_ptr = Libderiv->deriv_classes[6][6][0];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,36,dvrr_stack+33202, dvrr_stack+9876, dvrr_stack+26098);
 tmp = dvrr_stack + 33202;
 target_ptr = Libderiv->deriv_classes[6][7][0];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,45,dvrr_stack+9821, dvrr_stack+46306, dvrr_stack+30013);
 tmp = dvrr_stack + 9821;
 target_ptr = Libderiv->deriv_classes[6][8][0];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];


}

