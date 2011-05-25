#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (gg|gg) integrals */

void d1vrr_order_gggg(Libderiv_t *Libderiv, prim_data *Data)
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
 _BUILD_00p0(Data,dvrr_stack+6265, Data->F+11, Data->F+12, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+774, dvrr_stack+3556, dvrr_stack+6265, Data->F+10, Data->F+11, NULL);

 /* compute (0 0 | 3 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+6268, dvrr_stack+3559, dvrr_stack+774, dvrr_stack+611, dvrr_stack+3556, NULL);

 /* compute (0 0 | 4 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6278, dvrr_stack+3565, dvrr_stack+6268, dvrr_stack+1558, dvrr_stack+3559, NULL);

 /* compute (0 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6293, dvrr_stack+3575, dvrr_stack+6278, dvrr_stack+1564, dvrr_stack+3565, NULL);

 /* compute (0 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3808, dvrr_stack+1747, dvrr_stack+6293, dvrr_stack+723, dvrr_stack+3575, NULL);

 /* compute (0 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6314, dvrr_stack+3590, dvrr_stack+3808, dvrr_stack+1936, dvrr_stack+1747, NULL);

 /* compute (0 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+738, dvrr_stack+6314, dvrr_stack+1574, dvrr_stack+3590, NULL);

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
 _BUILD_00p0(Data,dvrr_stack+9757, Data->F+12, Data->F+13, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+9760, dvrr_stack+6265, dvrr_stack+9757, Data->F+11, Data->F+12, NULL);

 /* compute (0 0 | 3 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+9766, dvrr_stack+774, dvrr_stack+9760, dvrr_stack+3556, dvrr_stack+6265, NULL);

 /* compute (0 0 | 4 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+9776, dvrr_stack+6268, dvrr_stack+9766, dvrr_stack+3559, dvrr_stack+774, NULL);

 /* compute (0 0 | 5 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+9791, dvrr_stack+6278, dvrr_stack+9776, dvrr_stack+3565, dvrr_stack+6268, NULL);

 /* compute (0 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+9812, dvrr_stack+6293, dvrr_stack+9791, dvrr_stack+3575, dvrr_stack+6278, NULL);

 /* compute (0 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4052, dvrr_stack+3808, dvrr_stack+9812, dvrr_stack+1747, dvrr_stack+6293, NULL);

 /* compute (0 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9840, dvrr_stack+6314, dvrr_stack+4052, dvrr_stack+3590, dvrr_stack+3808, NULL);

 /* compute (0 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10217, dvrr_stack+1768, dvrr_stack+9840, dvrr_stack+738, dvrr_stack+6314, NULL);

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
 _BUILD_p000(Data,dvrr_stack+3556, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+864, dvrr_stack+0, dvrr_stack+30, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+14167, dvrr_stack+6, dvrr_stack+864, dvrr_stack+3, dvrr_stack+0, dvrr_stack+3556);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+14185, dvrr_stack+33, dvrr_stack+213, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+14203, dvrr_stack+39, dvrr_stack+14185, dvrr_stack+15, dvrr_stack+33, dvrr_stack+864);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+14239, dvrr_stack+75, dvrr_stack+14203, dvrr_stack+57, dvrr_stack+39, dvrr_stack+14167);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+479, dvrr_stack+219, dvrr_stack+873, NULL, NULL, dvrr_stack+213);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+14662, dvrr_stack+229, dvrr_stack+479, dvrr_stack+111, dvrr_stack+219, dvrr_stack+14185);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+14722, dvrr_stack+259, dvrr_stack+14662, dvrr_stack+131, dvrr_stack+229, dvrr_stack+14203);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+14822, dvrr_stack+379, dvrr_stack+14722, dvrr_stack+319, dvrr_stack+259, dvrr_stack+14239);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14972, dvrr_stack+883, dvrr_stack+575, NULL, NULL, dvrr_stack+873);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6589, dvrr_stack+898, dvrr_stack+14972, dvrr_stack+509, dvrr_stack+883, dvrr_stack+479);

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
 vrr_build_xxxx(am,Data,dvrr_stack+1813, dvrr_stack+738, dvrr_stack+6314, NULL, NULL, dvrr_stack+3590);

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
 vrr_build_xxxx(am,Data,dvrr_stack+3836, dvrr_stack+1768, dvrr_stack+9840, NULL, NULL, dvrr_stack+6314);

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
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+9757, dvrr_stack+6679, Data->F+12, Data->F+13, NULL);

 /* compute (0 0 | 3 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+9760, dvrr_stack+15, dvrr_stack+6265, dvrr_stack+9757, NULL);

 /* compute (0 0 | 4 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1921, dvrr_stack+9766, dvrr_stack+111, dvrr_stack+774, dvrr_stack+9760, NULL);

 /* compute (0 0 | 5 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6682, dvrr_stack+9776, dvrr_stack+1921, dvrr_stack+6268, dvrr_stack+9766, NULL);

 /* compute (0 0 | 6 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6703, dvrr_stack+9791, dvrr_stack+6682, dvrr_stack+6278, dvrr_stack+9776, NULL);

 /* compute (0 0 | 7 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6731, dvrr_stack+9812, dvrr_stack+6703, dvrr_stack+6293, dvrr_stack+9791, NULL);

 /* compute (0 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6767, dvrr_stack+4052, dvrr_stack+6731, dvrr_stack+3808, dvrr_stack+9812, NULL);

 /* compute (0 0 | 9 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+6812, dvrr_stack+9840, dvrr_stack+6767, dvrr_stack+6314, dvrr_stack+4052, NULL);

 /* compute (1 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+9885, dvrr_stack+10217, dvrr_stack+6812, NULL, NULL, dvrr_stack+9840);

 /* compute (2 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+33226, dvrr_stack+15157, dvrr_stack+9885, dvrr_stack+10272, dvrr_stack+10217, dvrr_stack+3836);

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
 _BUILD_p000(Data,dvrr_stack+6265, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+774, dvrr_stack+3556, dvrr_stack+6265, Data->F+4, Data->F+5, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+16202, dvrr_stack+30, dvrr_stack+210, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+16211, dvrr_stack+864, dvrr_stack+16202, dvrr_stack+0, dvrr_stack+30, dvrr_stack+6265);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+131, dvrr_stack+14167, dvrr_stack+16211, dvrr_stack+6, dvrr_stack+864, dvrr_stack+774);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+16229, dvrr_stack+213, dvrr_stack+704, NULL, NULL, dvrr_stack+210);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+16247, dvrr_stack+14185, dvrr_stack+16229, dvrr_stack+33, dvrr_stack+213, dvrr_stack+16202);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+509, dvrr_stack+14203, dvrr_stack+16247, dvrr_stack+39, dvrr_stack+14185, dvrr_stack+16211);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+16283, dvrr_stack+14239, dvrr_stack+509, dvrr_stack+75, dvrr_stack+14203, dvrr_stack+131);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+16373, dvrr_stack+873, dvrr_stack+121, NULL, NULL, dvrr_stack+704);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+16403, dvrr_stack+479, dvrr_stack+16373, dvrr_stack+219, dvrr_stack+873, dvrr_stack+16229);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+16463, dvrr_stack+14662, dvrr_stack+16403, dvrr_stack+229, dvrr_stack+479, dvrr_stack+16247);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+16563, dvrr_stack+14722, dvrr_stack+16463, dvrr_stack+259, dvrr_stack+14662, dvrr_stack+509);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10272, dvrr_stack+14822, dvrr_stack+16563, dvrr_stack+379, dvrr_stack+14722, dvrr_stack+16283);

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+219, dvrr_stack+575, dvrr_stack+723, NULL, NULL, dvrr_stack+121);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14299, dvrr_stack+14972, dvrr_stack+219, dvrr_stack+883, dvrr_stack+575, dvrr_stack+16373);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14389, dvrr_stack+6589, dvrr_stack+14299, dvrr_stack+898, dvrr_stack+14972, dvrr_stack+16403);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6350, dvrr_stack+19602, dvrr_stack+14389, dvrr_stack+943, dvrr_stack+6589, dvrr_stack+16463);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+38921, dvrr_stack+19752, dvrr_stack+6350, dvrr_stack+1033, dvrr_stack+19602, dvrr_stack+16563);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+39236, dvrr_stack+19977, dvrr_stack+38921, dvrr_stack+1333, dvrr_stack+19752, dvrr_stack+10272);
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
 vrr_build_xxxx(am,Data,dvrr_stack+39656, dvrr_stack+15017, dvrr_stack+946, dvrr_stack+1957, dvrr_stack+614, dvrr_stack+14299);

 /* compute (4 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+39866, dvrr_stack+20292, dvrr_stack+39656, dvrr_stack+2020, dvrr_stack+15017, dvrr_stack+14389);

 /* compute (5 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40181, dvrr_stack+20502, dvrr_stack+39866, dvrr_stack+2146, dvrr_stack+20292, dvrr_stack+6350);

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
 vrr_build_xxxx(am,Data,dvrr_stack+4088, dvrr_stack+6314, dvrr_stack+4052, NULL, NULL, dvrr_stack+3808);

 /* compute (2 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4196, dvrr_stack+1813, dvrr_stack+4088, dvrr_stack+738, dvrr_stack+6314, dvrr_stack+1957);

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
 vrr_build_xxxx(am,Data,dvrr_stack+4412, dvrr_stack+9840, dvrr_stack+6767, NULL, NULL, dvrr_stack+4052);

 /* compute (2 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6867, dvrr_stack+3836, dvrr_stack+4412, dvrr_stack+1768, dvrr_stack+9840, dvrr_stack+4088);

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
 tmp = dvrr_stack + 52942;
 target_ptr = Libderiv->dvrr_classes[6][8];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

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
 _BUILD_00f0(Data,dvrr_stack+7587, dvrr_stack+15, dvrr_stack+1768, dvrr_stack+9757, dvrr_stack+6679, NULL);

 /* compute (0 0 | 4 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1774, dvrr_stack+111, dvrr_stack+7587, dvrr_stack+9760, dvrr_stack+15, NULL);

 /* compute (0 0 | 5 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+590, dvrr_stack+1921, dvrr_stack+1774, dvrr_stack+9766, dvrr_stack+111, NULL);

 /* compute (0 0 | 6 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+738, dvrr_stack+6682, dvrr_stack+590, dvrr_stack+9776, dvrr_stack+1921, NULL);

 /* compute (0 0 | 7 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4547, dvrr_stack+6703, dvrr_stack+738, dvrr_stack+9791, dvrr_stack+6682, NULL);

 /* compute (0 0 | 8 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1574, dvrr_stack+6731, dvrr_stack+4547, dvrr_stack+9812, dvrr_stack+6703, NULL);

 /* compute (0 0 | 9 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+264, dvrr_stack+6767, dvrr_stack+1574, dvrr_stack+4052, dvrr_stack+6731, NULL);

 /* compute (1 0 | 9 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10482, dvrr_stack+6812, dvrr_stack+264, NULL, NULL, dvrr_stack+6767);

 /* compute (2 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+10647, dvrr_stack+9885, dvrr_stack+10482, dvrr_stack+10217, dvrr_stack+6812, dvrr_stack+4412);

 /* compute (3 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+57226, dvrr_stack+33226, dvrr_stack+10647, dvrr_stack+15157, dvrr_stack+9885, dvrr_stack+6867);

 /* compute (4 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+57776, dvrr_stack+33556, dvrr_stack+57226, dvrr_stack+15322, dvrr_stack+33226, dvrr_stack+7137);

 /* compute (5 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+58601, dvrr_stack+34106, dvrr_stack+57776, dvrr_stack+15652, dvrr_stack+33556, dvrr_stack+51322);

 /* compute (6 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+59756, dvrr_stack+34931, dvrr_stack+58601, dvrr_stack+16752, dvrr_stack+34106, dvrr_stack+51997);

 /* compute (6 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+61296,dvrr_stack+59756,dvrr_stack+52942,28);


 /* compute (1 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+9757, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+9760, dvrr_stack+6265, dvrr_stack+9757, Data->F+5, Data->F+6, NULL);

 /* compute (3 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+9766, dvrr_stack+774, dvrr_stack+9760, dvrr_stack+3556, dvrr_stack+6265, NULL);

 /* compute (1 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1789, dvrr_stack+210, dvrr_stack+161, NULL, NULL, Data->F+7);

 /* compute (2 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+4583, dvrr_stack+16202, dvrr_stack+1789, dvrr_stack+30, dvrr_stack+210, dvrr_stack+9757);

 /* compute (3 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+1619, dvrr_stack+16211, dvrr_stack+4583, dvrr_stack+864, dvrr_stack+16202, dvrr_stack+9760);

 /* compute (4 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+1649, dvrr_stack+131, dvrr_stack+1619, dvrr_stack+14167, dvrr_stack+16211, dvrr_stack+9766);

 /* compute (1 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+30, dvrr_stack+704, dvrr_stack+569, NULL, NULL, dvrr_stack+161);

 /* compute (2 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1694, dvrr_stack+16229, dvrr_stack+30, dvrr_stack+213, dvrr_stack+704, dvrr_stack+1789);

 /* compute (3 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+2209, dvrr_stack+16247, dvrr_stack+1694, dvrr_stack+14185, dvrr_stack+16229, dvrr_stack+4583);

 /* compute (4 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+1072, dvrr_stack+509, dvrr_stack+2209, dvrr_stack+14203, dvrr_stack+16247, dvrr_stack+1619);

 /* compute (5 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+10977, dvrr_stack+16283, dvrr_stack+1072, dvrr_stack+14239, dvrr_stack+509, dvrr_stack+1649);

 /* compute (1 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+14185, dvrr_stack+121, dvrr_stack+1564, NULL, NULL, dvrr_stack+569);

 /* compute (2 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+2269, dvrr_stack+16373, dvrr_stack+14185, dvrr_stack+873, dvrr_stack+121, dvrr_stack+30);

 /* compute (3 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+11103, dvrr_stack+16403, dvrr_stack+2269, dvrr_stack+479, dvrr_stack+16373, dvrr_stack+1694);

 /* compute (4 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+11203, dvrr_stack+16463, dvrr_stack+11103, dvrr_stack+14662, dvrr_stack+16403, dvrr_stack+2209);

 /* compute (5 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+15143, dvrr_stack+16563, dvrr_stack+11203, dvrr_stack+14722, dvrr_stack+16463, dvrr_stack+1072);

 /* compute (6 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+14539, dvrr_stack+10272, dvrr_stack+15143, dvrr_stack+14822, dvrr_stack+16563, dvrr_stack+10977);

 /* compute (1 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11353, dvrr_stack+723, dvrr_stack+3575, NULL, NULL, dvrr_stack+1564);

 /* compute (2 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15353, dvrr_stack+219, dvrr_stack+11353, dvrr_stack+575, dvrr_stack+723, dvrr_stack+14185);

 /* compute (3 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15443, dvrr_stack+14299, dvrr_stack+15353, dvrr_stack+14972, dvrr_stack+219, dvrr_stack+2269);

 /* compute (4 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15593, dvrr_stack+14389, dvrr_stack+15443, dvrr_stack+6589, dvrr_stack+14299, dvrr_stack+11103);

 /* compute (5 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15818, dvrr_stack+6350, dvrr_stack+15593, dvrr_stack+19602, dvrr_stack+14389, dvrr_stack+11203);

 /* compute (6 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+65076, dvrr_stack+38921, dvrr_stack+15818, dvrr_stack+19752, dvrr_stack+6350, dvrr_stack+15143);

 /* compute (7 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+65496, dvrr_stack+39236, dvrr_stack+65076, dvrr_stack+19977, dvrr_stack+38921, dvrr_stack+14539);
 tmp = dvrr_stack + 65496;
 target_ptr = Libderiv->dvrr_classes[7][4];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+19602, dvrr_stack+1747, dvrr_stack+6293, NULL, NULL, dvrr_stack+3575);

 /* compute (2 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+19665, dvrr_stack+883, dvrr_stack+19602, dvrr_stack+1936, dvrr_stack+1747, dvrr_stack+11353);

 /* compute (3 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+10050, dvrr_stack+946, dvrr_stack+19665, dvrr_stack+614, dvrr_stack+883, dvrr_stack+15353);

 /* compute (4 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+66036, dvrr_stack+39656, dvrr_stack+10050, dvrr_stack+15017, dvrr_stack+946, dvrr_stack+15443);

 /* compute (5 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+66351, dvrr_stack+39866, dvrr_stack+66036, dvrr_stack+20292, dvrr_stack+39656, dvrr_stack+15593);

 /* compute (6 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+66792, dvrr_stack+40181, dvrr_stack+66351, dvrr_stack+20502, dvrr_stack+39866, dvrr_stack+15818);

 /* compute (7 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+67380, dvrr_stack+40622, dvrr_stack+66792, dvrr_stack+20817, dvrr_stack+40181, dvrr_stack+65076);
 tmp = dvrr_stack + 67380;
 target_ptr = Libderiv->dvrr_classes[7][5];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+68136,dvrr_stack+67380,dvrr_stack+65496,36);


 /* compute (1 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+20292, dvrr_stack+3808, dvrr_stack+9812, NULL, NULL, dvrr_stack+6293);

 /* compute (2 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+20376, dvrr_stack+1957, dvrr_stack+20292, dvrr_stack+3590, dvrr_stack+3808, dvrr_stack+19602);

 /* compute (3 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+69756, dvrr_stack+2041, dvrr_stack+20376, dvrr_stack+780, dvrr_stack+1957, dvrr_stack+19665);

 /* compute (4 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+70036, dvrr_stack+42470, dvrr_stack+69756, dvrr_stack+22203, dvrr_stack+2041, dvrr_stack+10050);

 /* compute (5 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+70456, dvrr_stack+42750, dvrr_stack+70036, dvrr_stack+22371, dvrr_stack+42470, dvrr_stack+66036);

 /* compute (6 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+71044, dvrr_stack+43170, dvrr_stack+70456, dvrr_stack+22651, dvrr_stack+42750, dvrr_stack+66351);

 /* compute (7 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+71828, dvrr_stack+43758, dvrr_stack+71044, dvrr_stack+23071, dvrr_stack+43170, dvrr_stack+66792);
 tmp = dvrr_stack + 71828;
 target_ptr = Libderiv->dvrr_classes[7][6];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+72836,dvrr_stack+71828,dvrr_stack+67380,36);


 /* compute (1 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+22203, dvrr_stack+4052, dvrr_stack+6731, NULL, NULL, dvrr_stack+9812);

 /* compute (2 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+22311, dvrr_stack+4088, dvrr_stack+22203, dvrr_stack+6314, dvrr_stack+4052, dvrr_stack+20292);

 /* compute (3 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+22527, dvrr_stack+4196, dvrr_stack+22311, dvrr_stack+1813, dvrr_stack+4088, dvrr_stack+20376);

 /* compute (4 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+75104, dvrr_stack+46306, dvrr_stack+22527, dvrr_stack+24982, dvrr_stack+4196, dvrr_stack+69756);

 /* compute (5 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+75644, dvrr_stack+46666, dvrr_stack+75104, dvrr_stack+25198, dvrr_stack+46306, dvrr_stack+70036);

 /* compute (6 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+76400, dvrr_stack+47206, dvrr_stack+75644, dvrr_stack+25558, dvrr_stack+46666, dvrr_stack+70456);

 /* compute (7 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+77408, dvrr_stack+47962, dvrr_stack+76400, dvrr_stack+26098, dvrr_stack+47206, dvrr_stack+71044);
 tmp = dvrr_stack + 77408;
 target_ptr = Libderiv->dvrr_classes[7][7];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+78704,dvrr_stack+77408,dvrr_stack+71828,36);


 /* compute (1 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+24982, dvrr_stack+6767, dvrr_stack+1574, NULL, NULL, dvrr_stack+6731);

 /* compute (2 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+25117, dvrr_stack+4412, dvrr_stack+24982, dvrr_stack+9840, dvrr_stack+6767, dvrr_stack+22203);

 /* compute (3 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+25387, dvrr_stack+6867, dvrr_stack+25117, dvrr_stack+3836, dvrr_stack+4412, dvrr_stack+22311);

 /* compute (4 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+81728, dvrr_stack+7137, dvrr_stack+25387, dvrr_stack+28618, dvrr_stack+6867, dvrr_stack+22527);

 /* compute (5 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+82403, dvrr_stack+51322, dvrr_stack+81728, dvrr_stack+28888, dvrr_stack+7137, dvrr_stack+75104);

 /* compute (6 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+83348, dvrr_stack+51997, dvrr_stack+82403, dvrr_stack+29338, dvrr_stack+51322, dvrr_stack+75644);

 /* compute (7 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+84608, dvrr_stack+52942, dvrr_stack+83348, dvrr_stack+30013, dvrr_stack+51997, dvrr_stack+76400);
 tmp = dvrr_stack + 84608;
 target_ptr = Libderiv->dvrr_classes[7][8];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+86228,dvrr_stack+84608,dvrr_stack+77408,36);


 /* compute (0 0 | 1 0) m=15 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+14819, Data->F+15, Data->F+16, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=14 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+213, dvrr_stack+0, dvrr_stack+14819, Data->F+14, Data->F+15, NULL);

 /* compute (0 0 | 3 0) m=13 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+28618, dvrr_stack+1768, dvrr_stack+213, dvrr_stack+6679, dvrr_stack+0, NULL);

 /* compute (0 0 | 4 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+575, dvrr_stack+7587, dvrr_stack+28618, dvrr_stack+15, dvrr_stack+1768, NULL);

 /* compute (0 0 | 5 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1936, dvrr_stack+1774, dvrr_stack+575, dvrr_stack+111, dvrr_stack+7587, NULL);

 /* compute (0 0 | 6 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+28628, dvrr_stack+590, dvrr_stack+1936, dvrr_stack+1921, dvrr_stack+1774, NULL);

 /* compute (0 0 | 7 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6314, dvrr_stack+738, dvrr_stack+28628, dvrr_stack+6682, dvrr_stack+590, NULL);

 /* compute (0 0 | 8 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9840, dvrr_stack+4547, dvrr_stack+6314, dvrr_stack+6703, dvrr_stack+738, NULL);

 /* compute (0 0 | 9 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+28656, dvrr_stack+1574, dvrr_stack+9840, dvrr_stack+6731, dvrr_stack+4547, NULL);

 /* compute (1 0 | 9 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+28711, dvrr_stack+264, dvrr_stack+28656, NULL, NULL, dvrr_stack+1574);

 /* compute (2 0 | 9 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+28876, dvrr_stack+10482, dvrr_stack+28711, dvrr_stack+6812, dvrr_stack+264, dvrr_stack+24982);

 /* compute (3 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+29206, dvrr_stack+10647, dvrr_stack+28876, dvrr_stack+9885, dvrr_stack+10482, dvrr_stack+25117);

 /* compute (4 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+90116, dvrr_stack+57226, dvrr_stack+29206, dvrr_stack+33226, dvrr_stack+10647, dvrr_stack+25387);

 /* compute (5 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+90941, dvrr_stack+57776, dvrr_stack+90116, dvrr_stack+33556, dvrr_stack+57226, dvrr_stack+81728);

 /* compute (6 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+92096, dvrr_stack+58601, dvrr_stack+90941, dvrr_stack+34106, dvrr_stack+57776, dvrr_stack+82403);

 /* compute (7 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+93636, dvrr_stack+59756, dvrr_stack+92096, dvrr_stack+34931, dvrr_stack+58601, dvrr_stack+83348);

 /* compute (7 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+95616,dvrr_stack+93636,dvrr_stack+84608,36);


 /* compute (1 0 | 0 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+33226, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+15, dvrr_stack+9757, dvrr_stack+33226, Data->F+6, Data->F+7, NULL);

 /* compute (3 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+111, dvrr_stack+9760, dvrr_stack+15, dvrr_stack+6265, dvrr_stack+9757, NULL);

 /* compute (4 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+33229, dvrr_stack+9766, dvrr_stack+111, dvrr_stack+774, dvrr_stack+9760, NULL);

 /* compute (1 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+48, dvrr_stack+161, dvrr_stack+710, NULL, NULL, Data->F+8);

 /* compute (2 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+33244, dvrr_stack+1789, dvrr_stack+48, dvrr_stack+210, dvrr_stack+161, dvrr_stack+33226);

 /* compute (3 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+479, dvrr_stack+4583, dvrr_stack+33244, dvrr_stack+16202, dvrr_stack+1789, dvrr_stack+15);

 /* compute (4 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+33262, dvrr_stack+1619, dvrr_stack+479, dvrr_stack+16211, dvrr_stack+4583, dvrr_stack+111);

 /* compute (5 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+33307, dvrr_stack+1649, dvrr_stack+33262, dvrr_stack+131, dvrr_stack+1619, dvrr_stack+33229);

 /* compute (1 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+33370, dvrr_stack+569, dvrr_stack+1558, NULL, NULL, dvrr_stack+710);

 /* compute (2 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+33388, dvrr_stack+30, dvrr_stack+33370, dvrr_stack+704, dvrr_stack+569, dvrr_stack+48);

 /* compute (3 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+33424, dvrr_stack+1694, dvrr_stack+33388, dvrr_stack+16229, dvrr_stack+30, dvrr_stack+33244);

 /* compute (4 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+33484, dvrr_stack+2209, dvrr_stack+33424, dvrr_stack+16247, dvrr_stack+1694, dvrr_stack+479);

 /* compute (5 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+33574, dvrr_stack+1072, dvrr_stack+33484, dvrr_stack+509, dvrr_stack+2209, dvrr_stack+33262);

 /* compute (6 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+33700, dvrr_stack+10977, dvrr_stack+33574, dvrr_stack+16283, dvrr_stack+1072, dvrr_stack+33307);

 /* compute (1 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+509, dvrr_stack+1564, dvrr_stack+3565, NULL, NULL, dvrr_stack+1558);

 /* compute (2 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+33868, dvrr_stack+14185, dvrr_stack+509, dvrr_stack+121, dvrr_stack+1564, dvrr_stack+33370);

 /* compute (3 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+33928, dvrr_stack+2269, dvrr_stack+33868, dvrr_stack+16373, dvrr_stack+14185, dvrr_stack+33388);

 /* compute (4 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+16133, dvrr_stack+11103, dvrr_stack+33928, dvrr_stack+16403, dvrr_stack+2269, dvrr_stack+33424);

 /* compute (5 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+34028, dvrr_stack+11203, dvrr_stack+16133, dvrr_stack+16463, dvrr_stack+11103, dvrr_stack+33484);

 /* compute (6 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+34238, dvrr_stack+15143, dvrr_stack+34028, dvrr_stack+16563, dvrr_stack+11203, dvrr_stack+33574);

 /* compute (7 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+16373, dvrr_stack+14539, dvrr_stack+34238, dvrr_stack+10272, dvrr_stack+15143, dvrr_stack+33700);

 /* compute (1 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+34518, dvrr_stack+3575, dvrr_stack+6278, NULL, NULL, dvrr_stack+3565);

 /* compute (2 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+34563, dvrr_stack+11353, dvrr_stack+34518, dvrr_stack+723, dvrr_stack+3575, dvrr_stack+509);

 /* compute (3 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+34653, dvrr_stack+15353, dvrr_stack+34563, dvrr_stack+219, dvrr_stack+11353, dvrr_stack+33868);

 /* compute (4 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+29756, dvrr_stack+15443, dvrr_stack+34653, dvrr_stack+14299, dvrr_stack+15353, dvrr_stack+33928);

 /* compute (5 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+100476, dvrr_stack+15593, dvrr_stack+29756, dvrr_stack+14389, dvrr_stack+15443, dvrr_stack+16133);

 /* compute (6 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+100791, dvrr_stack+15818, dvrr_stack+100476, dvrr_stack+6350, dvrr_stack+15593, dvrr_stack+34028);

 /* compute (7 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+101211, dvrr_stack+65076, dvrr_stack+100791, dvrr_stack+38921, dvrr_stack+15818, dvrr_stack+34238);

 /* compute (8 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+101751, dvrr_stack+65496, dvrr_stack+101211, dvrr_stack+39236, dvrr_stack+65076, dvrr_stack+16373);
 tmp = dvrr_stack + 101751;
 target_ptr = Libderiv->dvrr_classes[8][4];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+38921, dvrr_stack+6293, dvrr_stack+9791, NULL, NULL, dvrr_stack+6278);

 /* compute (2 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+38984, dvrr_stack+19602, dvrr_stack+38921, dvrr_stack+1747, dvrr_stack+6293, dvrr_stack+34518);

 /* compute (3 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6350, dvrr_stack+19665, dvrr_stack+38984, dvrr_stack+883, dvrr_stack+19602, dvrr_stack+34563);

 /* compute (4 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+102426, dvrr_stack+10050, dvrr_stack+6350, dvrr_stack+946, dvrr_stack+19665, dvrr_stack+34653);

 /* compute (5 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+102741, dvrr_stack+66036, dvrr_stack+102426, dvrr_stack+39656, dvrr_stack+10050, dvrr_stack+29756);

 /* compute (6 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+103182, dvrr_stack+66351, dvrr_stack+102741, dvrr_stack+39866, dvrr_stack+66036, dvrr_stack+100476);

 /* compute (7 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+103770, dvrr_stack+66792, dvrr_stack+103182, dvrr_stack+40181, dvrr_stack+66351, dvrr_stack+100791);

 /* compute (8 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+39656, dvrr_stack+67380, dvrr_stack+103770, dvrr_stack+40622, dvrr_stack+66792, dvrr_stack+101211);
 tmp = dvrr_stack + 39656;
 target_ptr = Libderiv->dvrr_classes[8][5];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+104526,dvrr_stack+39656,dvrr_stack+101751,45);


 /* compute (1 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+39110, dvrr_stack+9812, dvrr_stack+6703, NULL, NULL, dvrr_stack+9791);

 /* compute (2 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14299, dvrr_stack+20292, dvrr_stack+39110, dvrr_stack+3808, dvrr_stack+9812, dvrr_stack+38921);

 /* compute (3 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+780, dvrr_stack+20376, dvrr_stack+14299, dvrr_stack+1957, dvrr_stack+20292, dvrr_stack+38984);

 /* compute (4 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3590, dvrr_stack+69756, dvrr_stack+780, dvrr_stack+2041, dvrr_stack+20376, dvrr_stack+6350);

 /* compute (5 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+106551, dvrr_stack+70036, dvrr_stack+3590, dvrr_stack+42470, dvrr_stack+69756, dvrr_stack+102426);

 /* compute (6 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+107139, dvrr_stack+70456, dvrr_stack+106551, dvrr_stack+42750, dvrr_stack+70036, dvrr_stack+102741);

 /* compute (7 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+107923, dvrr_stack+71044, dvrr_stack+107139, dvrr_stack+43170, dvrr_stack+70456, dvrr_stack+103182);

 /* compute (8 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+42470, dvrr_stack+71828, dvrr_stack+107923, dvrr_stack+43758, dvrr_stack+71044, dvrr_stack+103770);
 tmp = dvrr_stack + 42470;
 target_ptr = Libderiv->dvrr_classes[8][6];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+108931,dvrr_stack+42470,dvrr_stack+39656,45);


 /* compute (1 0 | 7 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1957, dvrr_stack+6731, dvrr_stack+4547, NULL, NULL, dvrr_stack+6703);

 /* compute (2 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+25837, dvrr_stack+22203, dvrr_stack+1957, dvrr_stack+4052, dvrr_stack+6731, dvrr_stack+39110);

 /* compute (3 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+111766, dvrr_stack+22311, dvrr_stack+25837, dvrr_stack+4088, dvrr_stack+22203, dvrr_stack+14299);

 /* compute (4 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+112126, dvrr_stack+22527, dvrr_stack+111766, dvrr_stack+4196, dvrr_stack+22311, dvrr_stack+780);

 /* compute (5 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+112666, dvrr_stack+75104, dvrr_stack+112126, dvrr_stack+46306, dvrr_stack+22527, dvrr_stack+3590);

 /* compute (6 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+113422, dvrr_stack+75644, dvrr_stack+112666, dvrr_stack+46666, dvrr_stack+75104, dvrr_stack+106551);

 /* compute (7 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+114430, dvrr_stack+76400, dvrr_stack+113422, dvrr_stack+47206, dvrr_stack+75644, dvrr_stack+107139);

 /* compute (8 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+46306, dvrr_stack+77408, dvrr_stack+114430, dvrr_stack+47962, dvrr_stack+76400, dvrr_stack+107923);
 tmp = dvrr_stack + 46306;
 target_ptr = Libderiv->dvrr_classes[8][7];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+115726,dvrr_stack+46306,dvrr_stack+42470,45);


 /* compute (1 0 | 8 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+2065, dvrr_stack+1574, dvrr_stack+9840, NULL, NULL, dvrr_stack+4547);

 /* compute (2 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+4010, dvrr_stack+24982, dvrr_stack+2065, dvrr_stack+6767, dvrr_stack+1574, dvrr_stack+1957);

 /* compute (3 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+119506, dvrr_stack+25117, dvrr_stack+4010, dvrr_stack+4412, dvrr_stack+24982, dvrr_stack+25837);

 /* compute (4 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+119956, dvrr_stack+25387, dvrr_stack+119506, dvrr_stack+6867, dvrr_stack+25117, dvrr_stack+111766);

 /* compute (5 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+120631, dvrr_stack+81728, dvrr_stack+119956, dvrr_stack+7137, dvrr_stack+25387, dvrr_stack+112126);

 /* compute (6 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+121576, dvrr_stack+82403, dvrr_stack+120631, dvrr_stack+51322, dvrr_stack+81728, dvrr_stack+112666);

 /* compute (7 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+122836, dvrr_stack+83348, dvrr_stack+121576, dvrr_stack+51997, dvrr_stack+82403, dvrr_stack+113422);

 /* compute (8 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+124456, dvrr_stack+84608, dvrr_stack+122836, dvrr_stack+52942, dvrr_stack+83348, dvrr_stack+114430);

 /* compute (8 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+126481,dvrr_stack+124456,dvrr_stack+46306,45);


 /* compute (0 0 | 1 0) m=16 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+210, Data->F+16, Data->F+17, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=15 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+51322, dvrr_stack+14819, dvrr_stack+210, Data->F+15, Data->F+16, NULL);

 /* compute (0 0 | 3 0) m=14 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+121, dvrr_stack+213, dvrr_stack+51322, dvrr_stack+0, dvrr_stack+14819, NULL);

 /* compute (0 0 | 4 0) m=13 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+723, dvrr_stack+28618, dvrr_stack+121, dvrr_stack+1768, dvrr_stack+213, NULL);

 /* compute (0 0 | 5 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40601, dvrr_stack+575, dvrr_stack+723, dvrr_stack+7587, dvrr_stack+28618, NULL);

 /* compute (0 0 | 6 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+43730, dvrr_stack+1936, dvrr_stack+40601, dvrr_stack+1774, dvrr_stack+575, NULL);

 /* compute (0 0 | 7 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+47926, dvrr_stack+28628, dvrr_stack+43730, dvrr_stack+590, dvrr_stack+1936, NULL);

 /* compute (0 0 | 8 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+26053, dvrr_stack+6314, dvrr_stack+47926, dvrr_stack+738, dvrr_stack+28628, NULL);

 /* compute (0 0 | 9 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+51322, dvrr_stack+9840, dvrr_stack+26053, dvrr_stack+4547, dvrr_stack+6314, NULL);

 /* compute (1 0 | 9 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+9885, dvrr_stack+28656, dvrr_stack+51322, NULL, NULL, dvrr_stack+9840);

 /* compute (2 0 | 9 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+51322, dvrr_stack+28711, dvrr_stack+9885, dvrr_stack+264, dvrr_stack+28656, dvrr_stack+2065);

 /* compute (3 0 | 9 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+51652, dvrr_stack+28876, dvrr_stack+51322, dvrr_stack+10482, dvrr_stack+28711, dvrr_stack+4010);

 /* compute (4 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+6767, dvrr_stack+29206, dvrr_stack+51652, dvrr_stack+10647, dvrr_stack+28876, dvrr_stack+119506);

 /* compute (5 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+51322, dvrr_stack+90116, dvrr_stack+6767, dvrr_stack+57226, dvrr_stack+29206, dvrr_stack+119956);

 /* compute (6 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+131341, dvrr_stack+90941, dvrr_stack+51322, dvrr_stack+57776, dvrr_stack+90116, dvrr_stack+120631);

 /* compute (7 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+132881, dvrr_stack+92096, dvrr_stack+131341, dvrr_stack+58601, dvrr_stack+90941, dvrr_stack+121576);

 /* compute (8 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+57226, dvrr_stack+93636, dvrr_stack+132881, dvrr_stack+59756, dvrr_stack+92096, dvrr_stack+122836);

 /* compute (8 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+131341,dvrr_stack+57226,dvrr_stack+124456,45);


 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+14819, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+2200, dvrr_stack+21, dvrr_stack+3, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+137416, dvrr_stack+2200, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+14819);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+137434, dvrr_stack+164, dvrr_stack+24, NULL, NULL, dvrr_stack+21);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+47926, dvrr_stack+137434, dvrr_stack+57, dvrr_stack+164, dvrr_stack+24, dvrr_stack+2200);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+137452, dvrr_stack+47926, dvrr_stack+75, dvrr_stack+137434, dvrr_stack+57, dvrr_stack+137416);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+539, dvrr_stack+713, dvrr_stack+170, NULL, NULL, dvrr_stack+164);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+137512, dvrr_stack+539, dvrr_stack+180, dvrr_stack+713, dvrr_stack+170, dvrr_stack+137434);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+137572, dvrr_stack+137512, dvrr_stack+319, dvrr_stack+539, dvrr_stack+180, dvrr_stack+47926);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+164, dvrr_stack+137572, dvrr_stack+379, dvrr_stack+137512, dvrr_stack+319, dvrr_stack+137452);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+0, dvrr_stack+14819, dvrr_stack+3556, Data->F+3, Data->F+4, NULL);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+539, dvrr_stack+137416, dvrr_stack+14167, dvrr_stack+2200, dvrr_stack+6, dvrr_stack+0);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+137672, dvrr_stack+137452, dvrr_stack+14239, dvrr_stack+47926, dvrr_stack+75, dvrr_stack+539);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+137762, dvrr_stack+164, dvrr_stack+14822, dvrr_stack+137572, dvrr_stack+379, dvrr_stack+137672);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+121, dvrr_stack+0, dvrr_stack+774, dvrr_stack+14819, dvrr_stack+3556, NULL);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+137512, dvrr_stack+539, dvrr_stack+131, dvrr_stack+137416, dvrr_stack+14167, dvrr_stack+121);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+314, dvrr_stack+137672, dvrr_stack+16283, dvrr_stack+137452, dvrr_stack+14239, dvrr_stack+137512);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+90116, dvrr_stack+137762, dvrr_stack+10272, dvrr_stack+164, dvrr_stack+14822, dvrr_stack+314);

 /* compute (4 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+14167, dvrr_stack+121, dvrr_stack+9766, dvrr_stack+0, dvrr_stack+774, NULL);

 /* compute (5 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+137416, dvrr_stack+137512, dvrr_stack+1649, dvrr_stack+539, dvrr_stack+131, dvrr_stack+14167);

 /* compute (6 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+14819, dvrr_stack+314, dvrr_stack+10977, dvrr_stack+137672, dvrr_stack+16283, dvrr_stack+137416);

 /* compute (7 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+90396, dvrr_stack+90116, dvrr_stack+14539, dvrr_stack+137762, dvrr_stack+10272, dvrr_stack+14819);

 /* compute (5 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+40601, dvrr_stack+14167, dvrr_stack+33229, dvrr_stack+121, dvrr_stack+9766, NULL);

 /* compute (6 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+14215, dvrr_stack+137416, dvrr_stack+33307, dvrr_stack+137512, dvrr_stack+1649, dvrr_stack+40601);

 /* compute (7 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+137416, dvrr_stack+14819, dvrr_stack+33700, dvrr_stack+314, dvrr_stack+10977, dvrr_stack+14215);

 /* compute (8 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+90756, dvrr_stack+90396, dvrr_stack+16373, dvrr_stack+90116, dvrr_stack+14539, dvrr_stack+137416);

 /* compute (1 0 | 0 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3556, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+137416, dvrr_stack+33226, dvrr_stack+3556, Data->F+7, Data->F+8, NULL);

 /* compute (3 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+137422, dvrr_stack+15, dvrr_stack+137416, dvrr_stack+9757, dvrr_stack+33226, NULL);

 /* compute (4 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+111, dvrr_stack+137422, dvrr_stack+9760, dvrr_stack+15, NULL);

 /* compute (5 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+40601, dvrr_stack+33229, dvrr_stack+0, dvrr_stack+9766, dvrr_stack+111, NULL);

 /* compute (1 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+2200, dvrr_stack+710, dvrr_stack+611, NULL, NULL, Data->F+9);

 /* compute (2 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+33226, dvrr_stack+48, dvrr_stack+2200, dvrr_stack+161, dvrr_stack+710, dvrr_stack+3556);

 /* compute (3 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+539, dvrr_stack+33244, dvrr_stack+33226, dvrr_stack+1789, dvrr_stack+48, dvrr_stack+137416);

 /* compute (4 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+48, dvrr_stack+479, dvrr_stack+539, dvrr_stack+4583, dvrr_stack+33244, dvrr_stack+137422);

 /* compute (5 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+137416, dvrr_stack+33262, dvrr_stack+48, dvrr_stack+1619, dvrr_stack+479, dvrr_stack+0);

 /* compute (6 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+14215, dvrr_stack+33307, dvrr_stack+137416, dvrr_stack+1649, dvrr_stack+33262, dvrr_stack+40601);

 /* compute (1 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+14167, dvrr_stack+1558, dvrr_stack+3559, NULL, NULL, dvrr_stack+611);

 /* compute (2 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+47926, dvrr_stack+33370, dvrr_stack+14167, dvrr_stack+569, dvrr_stack+1558, dvrr_stack+2200);

 /* compute (3 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+569, dvrr_stack+33388, dvrr_stack+47926, dvrr_stack+30, dvrr_stack+33370, dvrr_stack+33226);

 /* compute (4 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+16283, dvrr_stack+33424, dvrr_stack+569, dvrr_stack+1694, dvrr_stack+33388, dvrr_stack+539);

 /* compute (5 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+33226, dvrr_stack+33484, dvrr_stack+16283, dvrr_stack+2209, dvrr_stack+33424, dvrr_stack+48);

 /* compute (6 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+1619, dvrr_stack+33574, dvrr_stack+33226, dvrr_stack+1072, dvrr_stack+33484, dvrr_stack+137416);

 /* compute (7 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+137416, dvrr_stack+33700, dvrr_stack+1619, dvrr_stack+10977, dvrr_stack+33574, dvrr_stack+14215);

 /* compute (1 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+539, dvrr_stack+3565, dvrr_stack+6268, NULL, NULL, dvrr_stack+3559);

 /* compute (2 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+14215, dvrr_stack+509, dvrr_stack+539, dvrr_stack+1564, dvrr_stack+3565, dvrr_stack+14167);

 /* compute (3 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+137632, dvrr_stack+33868, dvrr_stack+14215, dvrr_stack+14185, dvrr_stack+509, dvrr_stack+47926);

 /* compute (4 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+33352, dvrr_stack+33928, dvrr_stack+137632, dvrr_stack+2269, dvrr_stack+33868, dvrr_stack+569);

 /* compute (5 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+33502, dvrr_stack+16133, dvrr_stack+33352, dvrr_stack+11103, dvrr_stack+33928, dvrr_stack+16283);

 /* compute (6 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+33712, dvrr_stack+34028, dvrr_stack+33502, dvrr_stack+11203, dvrr_stack+16133, dvrr_stack+33226);

 /* compute (7 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+91206, dvrr_stack+34238, dvrr_stack+33712, dvrr_stack+15143, dvrr_stack+34028, dvrr_stack+1619);

 /* compute (8 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+14819, dvrr_stack+16373, dvrr_stack+91206, dvrr_stack+14539, dvrr_stack+34238, dvrr_stack+137416);

 /* compute (1 0 | 4 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+137416, dvrr_stack+6278, dvrr_stack+9776, NULL, NULL, dvrr_stack+6268);

 /* compute (2 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+137461, dvrr_stack+34518, dvrr_stack+137416, dvrr_stack+3575, dvrr_stack+6278, dvrr_stack+539);

 /* compute (3 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1619, dvrr_stack+34563, dvrr_stack+137461, dvrr_stack+11353, dvrr_stack+34518, dvrr_stack+14215);

 /* compute (4 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+16133, dvrr_stack+34653, dvrr_stack+1619, dvrr_stack+15353, dvrr_stack+34563, dvrr_stack+137632);

 /* compute (5 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+16358, dvrr_stack+29756, dvrr_stack+16133, dvrr_stack+15443, dvrr_stack+34653, dvrr_stack+33352);

 /* compute (6 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+33992, dvrr_stack+100476, dvrr_stack+16358, dvrr_stack+15593, dvrr_stack+29756, dvrr_stack+33502);

 /* compute (7 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15269, dvrr_stack+100791, dvrr_stack+33992, dvrr_stack+15818, dvrr_stack+100476, dvrr_stack+33712);

 /* compute (8 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+33226, dvrr_stack+101211, dvrr_stack+15269, dvrr_stack+65076, dvrr_stack+100791, dvrr_stack+91206);

 /* compute (9 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 9;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+91206, dvrr_stack+101751, dvrr_stack+33226, dvrr_stack+65496, dvrr_stack+101211, dvrr_stack+14819);

 /* compute (1 0 | 5 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+92031, dvrr_stack+9791, dvrr_stack+6682, NULL, NULL, dvrr_stack+9776);

 /* compute (2 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+92094, dvrr_stack+38921, dvrr_stack+92031, dvrr_stack+6293, dvrr_stack+9791, dvrr_stack+137416);

 /* compute (3 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+92220, dvrr_stack+38984, dvrr_stack+92094, dvrr_stack+19602, dvrr_stack+38921, dvrr_stack+137461);

 /* compute (4 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+137416, dvrr_stack+6350, dvrr_stack+92220, dvrr_stack+19665, dvrr_stack+38984, dvrr_stack+1619);

 /* compute (5 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+92430, dvrr_stack+102426, dvrr_stack+137416, dvrr_stack+10050, dvrr_stack+6350, dvrr_stack+16133);

 /* compute (6 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+92871, dvrr_stack+102741, dvrr_stack+92430, dvrr_stack+66036, dvrr_stack+102426, dvrr_stack+16358);

 /* compute (7 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+100476, dvrr_stack+103182, dvrr_stack+92871, dvrr_stack+66351, dvrr_stack+102741, dvrr_stack+33992);

 /* compute (8 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+33901, dvrr_stack+103770, dvrr_stack+100476, dvrr_stack+66792, dvrr_stack+103182, dvrr_stack+15269);

 /* compute (9 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 9;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+102426, dvrr_stack+39656, dvrr_stack+33901, dvrr_stack+67380, dvrr_stack+103770, dvrr_stack+33226);

 /* compute (1 0 | 6 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+33226, dvrr_stack+6703, dvrr_stack+738, NULL, NULL, dvrr_stack+6682);

 /* compute (2 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+33310, dvrr_stack+39110, dvrr_stack+33226, dvrr_stack+9812, dvrr_stack+6703, dvrr_stack+92031);

 /* compute (3 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+33478, dvrr_stack+14299, dvrr_stack+33310, dvrr_stack+20292, dvrr_stack+39110, dvrr_stack+92094);

 /* compute (4 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+65076, dvrr_stack+780, dvrr_stack+33478, dvrr_stack+20376, dvrr_stack+14299, dvrr_stack+92220);

 /* compute (5 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+103581, dvrr_stack+3590, dvrr_stack+65076, dvrr_stack+69756, dvrr_stack+780, dvrr_stack+137416);

 /* compute (6 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+66036, dvrr_stack+106551, dvrr_stack+103581, dvrr_stack+70036, dvrr_stack+3590, dvrr_stack+92430);

 /* compute (7 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14167, dvrr_stack+107139, dvrr_stack+66036, dvrr_stack+70456, dvrr_stack+106551, dvrr_stack+92871);

 /* compute (8 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+69756, dvrr_stack+107923, dvrr_stack+14167, dvrr_stack+71044, dvrr_stack+107139, dvrr_stack+100476);

 /* compute (9 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 9;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+92031, dvrr_stack+42470, dvrr_stack+69756, dvrr_stack+71828, dvrr_stack+107923, dvrr_stack+33901);

 /* compute (1 0 | 7 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+100476, dvrr_stack+4547, dvrr_stack+6314, NULL, NULL, dvrr_stack+738);

 /* compute (2 0 | 7 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+100584, dvrr_stack+1957, dvrr_stack+100476, dvrr_stack+6731, dvrr_stack+4547, dvrr_stack+33226);

 /* compute (3 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+100800, dvrr_stack+25837, dvrr_stack+100584, dvrr_stack+22203, dvrr_stack+1957, dvrr_stack+33310);

 /* compute (4 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+101160, dvrr_stack+111766, dvrr_stack+100800, dvrr_stack+22311, dvrr_stack+25837, dvrr_stack+33478);

 /* compute (5 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+33226, dvrr_stack+112126, dvrr_stack+101160, dvrr_stack+22527, dvrr_stack+111766, dvrr_stack+65076);

 /* compute (6 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+106551, dvrr_stack+112666, dvrr_stack+33226, dvrr_stack+75104, dvrr_stack+112126, dvrr_stack+103581);

 /* compute (7 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+107559, dvrr_stack+113422, dvrr_stack+106551, dvrr_stack+75644, dvrr_stack+112666, dvrr_stack+66036);

 /* compute (8 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+51322, dvrr_stack+114430, dvrr_stack+107559, dvrr_stack+76400, dvrr_stack+113422, dvrr_stack+14167);

 /* compute (9 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 9;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+14167, dvrr_stack+46306, dvrr_stack+51322, dvrr_stack+77408, dvrr_stack+114430, dvrr_stack+69756);

 /* compute (1 0 | 8 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+69756, dvrr_stack+9840, dvrr_stack+26053, NULL, NULL, dvrr_stack+6314);

 /* compute (2 0 | 8 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+69891, dvrr_stack+2065, dvrr_stack+69756, dvrr_stack+1574, dvrr_stack+9840, dvrr_stack+100476);

 /* compute (3 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+70161, dvrr_stack+4010, dvrr_stack+69891, dvrr_stack+24982, dvrr_stack+2065, dvrr_stack+100584);

 /* compute (4 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+70611, dvrr_stack+119506, dvrr_stack+70161, dvrr_stack+25117, dvrr_stack+4010, dvrr_stack+100800);

 /* compute (5 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+103581, dvrr_stack+119956, dvrr_stack+70611, dvrr_stack+25387, dvrr_stack+119506, dvrr_stack+101160);

 /* compute (6 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+100476, dvrr_stack+120631, dvrr_stack+103581, dvrr_stack+81728, dvrr_stack+119956, dvrr_stack+33226);

 /* compute (7 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+33226, dvrr_stack+121576, dvrr_stack+100476, dvrr_stack+82403, dvrr_stack+120631, dvrr_stack+106551);

 /* compute (8 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+119506, dvrr_stack+122836, dvrr_stack+33226, dvrr_stack+83348, dvrr_stack+121576, dvrr_stack+107559);

 /* compute (9 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 9;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+81728, dvrr_stack+124456, dvrr_stack+119506, dvrr_stack+84608, dvrr_stack+122836, dvrr_stack+51322);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+51322, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 51322;
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
 deriv_build_DZ_0(Data,420,dvrr_stack+65076, dvrr_stack+8497, NULL);
 tmp = dvrr_stack + 65076;
 target_ptr = Libderiv->deriv_classes[4][6][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,540,dvrr_stack+51547, dvrr_stack+12547, NULL);
 tmp = dvrr_stack + 51547;
 target_ptr = Libderiv->deriv_classes[4][7][11];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,675,dvrr_stack+52087, dvrr_stack+17577, NULL);
 tmp = dvrr_stack + 52087;
 target_ptr = Libderiv->deriv_classes[4][8][11];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+119506, dvrr_stack+21258, NULL);
 tmp = dvrr_stack + 119506;
 target_ptr = Libderiv->deriv_classes[5][4][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,441,dvrr_stack+119821, dvrr_stack+23659, NULL);
 tmp = dvrr_stack + 119821;
 target_ptr = Libderiv->deriv_classes[5][5][11];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+120262, dvrr_stack+26854, NULL);
 tmp = dvrr_stack + 120262;
 target_ptr = Libderiv->deriv_classes[5][6][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,756,dvrr_stack+120850, dvrr_stack+30958, NULL);
 tmp = dvrr_stack + 120850;
 target_ptr = Libderiv->deriv_classes[5][7][11];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,945,dvrr_stack+103581, dvrr_stack+36086, NULL);
 tmp = dvrr_stack + 103581;
 target_ptr = Libderiv->deriv_classes[5][8][11];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+121606, dvrr_stack+41210, NULL);
 tmp = dvrr_stack + 121606;
 target_ptr = Libderiv->deriv_classes[6][4][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+122026, dvrr_stack+44542, NULL);
 tmp = dvrr_stack + 122026;
 target_ptr = Libderiv->deriv_classes[6][5][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,784,dvrr_stack+122614, dvrr_stack+48970, NULL);
 tmp = dvrr_stack + 122614;
 target_ptr = Libderiv->deriv_classes[6][6][11];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1008,dvrr_stack+123398, dvrr_stack+54202, NULL);
 tmp = dvrr_stack + 123398;
 target_ptr = Libderiv->deriv_classes[6][7][11];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1260,dvrr_stack+33226, dvrr_stack+61296, NULL);
 tmp = dvrr_stack + 33226;
 target_ptr = Libderiv->deriv_classes[6][8][11];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,540,dvrr_stack+106551, dvrr_stack+68136, NULL);
 tmp = dvrr_stack + 106551;
 target_ptr = Libderiv->deriv_classes[7][4][11];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,756,dvrr_stack+107091, dvrr_stack+72836, NULL);
 tmp = dvrr_stack + 107091;
 target_ptr = Libderiv->deriv_classes[7][5][11];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1008,dvrr_stack+107847, dvrr_stack+78704, NULL);
 tmp = dvrr_stack + 107847;
 target_ptr = Libderiv->deriv_classes[7][6][11];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1296,dvrr_stack+69756, dvrr_stack+86228, NULL);
 tmp = dvrr_stack + 69756;
 target_ptr = Libderiv->deriv_classes[7][7][11];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1620,dvrr_stack+75104, dvrr_stack+95616, NULL);
 tmp = dvrr_stack + 75104;
 target_ptr = Libderiv->deriv_classes[7][8][11];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,675,dvrr_stack+100476, dvrr_stack+104526, NULL);
 tmp = dvrr_stack + 100476;
 target_ptr = Libderiv->deriv_classes[8][4][11];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,945,dvrr_stack+24982, dvrr_stack+108931, NULL);
 tmp = dvrr_stack + 24982;
 target_ptr = Libderiv->deriv_classes[8][5][11];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1260,dvrr_stack+66036, dvrr_stack+115726, NULL);
 tmp = dvrr_stack + 66036;
 target_ptr = Libderiv->deriv_classes[8][6][11];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1620,dvrr_stack+111766, dvrr_stack+126481, NULL);
 tmp = dvrr_stack + 111766;
 target_ptr = Libderiv->deriv_classes[8][7][11];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,2025,dvrr_stack+113386, dvrr_stack+131341, NULL);
 tmp = dvrr_stack + 113386;
 target_ptr = Libderiv->deriv_classes[8][8][11];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+34486, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 34486;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+115411, dvrr_stack+5320, NULL);
 tmp = dvrr_stack + 115411;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+101151, dvrr_stack+8497, NULL);
 tmp = dvrr_stack + 101151;
 target_ptr = Libderiv->deriv_classes[4][6][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,540,dvrr_stack+71052, dvrr_stack+12547, NULL);
 tmp = dvrr_stack + 71052;
 target_ptr = Libderiv->deriv_classes[4][7][10];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,675,dvrr_stack+76724, dvrr_stack+17577, NULL);
 tmp = dvrr_stack + 76724;
 target_ptr = Libderiv->deriv_classes[4][8][10];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+84203, dvrr_stack+21258, NULL);
 tmp = dvrr_stack + 84203;
 target_ptr = Libderiv->deriv_classes[5][4][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+16147, dvrr_stack+23659, NULL);
 tmp = dvrr_stack + 16147;
 target_ptr = Libderiv->deriv_classes[5][5][10];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+22203, dvrr_stack+26854, NULL);
 tmp = dvrr_stack + 22203;
 target_ptr = Libderiv->deriv_classes[5][6][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,756,dvrr_stack+1558, dvrr_stack+30958, NULL);
 tmp = dvrr_stack + 1558;
 target_ptr = Libderiv->deriv_classes[5][7][10];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,945,dvrr_stack+3556, dvrr_stack+36086, NULL);
 tmp = dvrr_stack + 3556;
 target_ptr = Libderiv->deriv_classes[5][8][10];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+20292, dvrr_stack+41210, NULL);
 tmp = dvrr_stack + 20292;
 target_ptr = Libderiv->deriv_classes[6][4][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+9757, dvrr_stack+44542, NULL);
 tmp = dvrr_stack + 9757;
 target_ptr = Libderiv->deriv_classes[6][5][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,784,dvrr_stack+10345, dvrr_stack+48970, NULL);
 tmp = dvrr_stack + 10345;
 target_ptr = Libderiv->deriv_classes[6][6][10];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1008,dvrr_stack+28618, dvrr_stack+54202, NULL);
 tmp = dvrr_stack + 28618;
 target_ptr = Libderiv->deriv_classes[6][7][10];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1260,dvrr_stack+6265, dvrr_stack+61296, NULL);
 tmp = dvrr_stack + 6265;
 target_ptr = Libderiv->deriv_classes[6][8][10];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,540,dvrr_stack+314, dvrr_stack+68136, NULL);
 tmp = dvrr_stack + 314;
 target_ptr = Libderiv->deriv_classes[7][4][10];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,756,dvrr_stack+137972, dvrr_stack+72836, NULL);
 tmp = dvrr_stack + 137972;
 target_ptr = Libderiv->deriv_classes[7][5][10];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1008,dvrr_stack+138728, dvrr_stack+78704, NULL);
 tmp = dvrr_stack + 138728;
 target_ptr = Libderiv->deriv_classes[7][6][10];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1296,dvrr_stack+139736, dvrr_stack+86228, NULL);
 tmp = dvrr_stack + 139736;
 target_ptr = Libderiv->deriv_classes[7][7][10];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1620,dvrr_stack+141032, dvrr_stack+95616, NULL);
 tmp = dvrr_stack + 141032;
 target_ptr = Libderiv->deriv_classes[7][8][10];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,675,dvrr_stack+142652, dvrr_stack+104526, NULL);
 tmp = dvrr_stack + 142652;
 target_ptr = Libderiv->deriv_classes[8][4][10];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,945,dvrr_stack+143327, dvrr_stack+108931, NULL);
 tmp = dvrr_stack + 143327;
 target_ptr = Libderiv->deriv_classes[8][5][10];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1260,dvrr_stack+144272, dvrr_stack+115726, NULL);
 tmp = dvrr_stack + 144272;
 target_ptr = Libderiv->deriv_classes[8][6][10];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1620,dvrr_stack+145532, dvrr_stack+126481, NULL);
 tmp = dvrr_stack + 145532;
 target_ptr = Libderiv->deriv_classes[8][7][10];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,2025,dvrr_stack+147152, dvrr_stack+131341, NULL);
 tmp = dvrr_stack + 147152;
 target_ptr = Libderiv->deriv_classes[8][8][10];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+71592, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 71592;
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
 deriv_build_DX_0(Data,1260,dvrr_stack+54202, dvrr_stack+61296, NULL);
 tmp = dvrr_stack + 54202;
 target_ptr = Libderiv->deriv_classes[6][8][9];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,540,dvrr_stack+61296, dvrr_stack+68136, NULL);
 tmp = dvrr_stack + 61296;
 target_ptr = Libderiv->deriv_classes[7][4][9];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,756,dvrr_stack+68136, dvrr_stack+72836, NULL);
 tmp = dvrr_stack + 68136;
 target_ptr = Libderiv->deriv_classes[7][5][9];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1008,dvrr_stack+72836, dvrr_stack+78704, NULL);
 tmp = dvrr_stack + 72836;
 target_ptr = Libderiv->deriv_classes[7][6][9];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1296,dvrr_stack+78704, dvrr_stack+86228, NULL);
 tmp = dvrr_stack + 78704;
 target_ptr = Libderiv->deriv_classes[7][7][9];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1620,dvrr_stack+86228, dvrr_stack+95616, NULL);
 tmp = dvrr_stack + 86228;
 target_ptr = Libderiv->deriv_classes[7][8][9];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,675,dvrr_stack+95616, dvrr_stack+104526, NULL);
 tmp = dvrr_stack + 95616;
 target_ptr = Libderiv->deriv_classes[8][4][9];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,945,dvrr_stack+104526, dvrr_stack+108931, NULL);
 tmp = dvrr_stack + 104526;
 target_ptr = Libderiv->deriv_classes[8][5][9];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1260,dvrr_stack+73844, dvrr_stack+115726, NULL);
 tmp = dvrr_stack + 73844;
 target_ptr = Libderiv->deriv_classes[8][6][9];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1620,dvrr_stack+115726, dvrr_stack+126481, NULL);
 tmp = dvrr_stack + 115726;
 target_ptr = Libderiv->deriv_classes[8][7][9];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,2025,dvrr_stack+126481, dvrr_stack+131341, NULL);
 tmp = dvrr_stack + 126481;
 target_ptr = Libderiv->deriv_classes[8][8][9];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+128506, dvrr_stack+2566, dvrr_stack+164);
 tmp = dvrr_stack + 128506;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+128731, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 128731;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+129046, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 129046;
 target_ptr = Libderiv->deriv_classes[4][6][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,15,1,dvrr_stack+129466, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 129466;
 target_ptr = Libderiv->deriv_classes[4][7][8];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,15,1,dvrr_stack+130006, dvrr_stack+16752, dvrr_stack+7957);
 tmp = dvrr_stack + 130006;
 target_ptr = Libderiv->deriv_classes[4][8][8];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+130681, dvrr_stack+20817, dvrr_stack+137762);
 tmp = dvrr_stack + 130681;
 target_ptr = Libderiv->deriv_classes[5][4][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+130996, dvrr_stack+23071, dvrr_stack+19977);
 tmp = dvrr_stack + 130996;
 target_ptr = Libderiv->deriv_classes[5][5][8];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+131437, dvrr_stack+26098, dvrr_stack+20817);
 tmp = dvrr_stack + 131437;
 target_ptr = Libderiv->deriv_classes[5][6][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,21,1,dvrr_stack+132025, dvrr_stack+30013, dvrr_stack+23071);
 tmp = dvrr_stack + 132025;
 target_ptr = Libderiv->deriv_classes[5][7][8];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,21,1,dvrr_stack+132781, dvrr_stack+34931, dvrr_stack+26098);
 tmp = dvrr_stack + 132781;
 target_ptr = Libderiv->deriv_classes[5][8][8];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+133726, dvrr_stack+40622, dvrr_stack+90116);
 tmp = dvrr_stack + 133726;
 target_ptr = Libderiv->deriv_classes[6][4][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+134146, dvrr_stack+43758, dvrr_stack+39236);
 tmp = dvrr_stack + 134146;
 target_ptr = Libderiv->deriv_classes[6][5][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+134734, dvrr_stack+47962, dvrr_stack+40622);
 tmp = dvrr_stack + 134734;
 target_ptr = Libderiv->deriv_classes[6][6][8];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,28,1,dvrr_stack+135518, dvrr_stack+52942, dvrr_stack+43758);
 tmp = dvrr_stack + 135518;
 target_ptr = Libderiv->deriv_classes[6][7][8];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,28,1,dvrr_stack+117346, dvrr_stack+59756, dvrr_stack+47962);
 tmp = dvrr_stack + 117346;
 target_ptr = Libderiv->deriv_classes[6][8][8];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,36,1,dvrr_stack+136526, dvrr_stack+67380, dvrr_stack+90396);
 tmp = dvrr_stack + 136526;
 target_ptr = Libderiv->deriv_classes[7][4][8];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,36,1,dvrr_stack+118606, dvrr_stack+71828, dvrr_stack+65496);
 tmp = dvrr_stack + 118606;
 target_ptr = Libderiv->deriv_classes[7][5][8];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,36,1,dvrr_stack+105471, dvrr_stack+77408, dvrr_stack+67380);
 tmp = dvrr_stack + 105471;
 target_ptr = Libderiv->deriv_classes[7][6][8];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,36,1,dvrr_stack+96291, dvrr_stack+84608, dvrr_stack+71828);
 tmp = dvrr_stack + 96291;
 target_ptr = Libderiv->deriv_classes[7][7][8];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,36,1,dvrr_stack+97587, dvrr_stack+93636, dvrr_stack+77408);
 tmp = dvrr_stack + 97587;
 target_ptr = Libderiv->deriv_classes[7][8][8];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,45,1,dvrr_stack+137066, dvrr_stack+39656, dvrr_stack+90756);
 tmp = dvrr_stack + 137066;
 target_ptr = Libderiv->deriv_classes[8][4][8];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,45,1,dvrr_stack+99207, dvrr_stack+42470, dvrr_stack+101751);
 tmp = dvrr_stack + 99207;
 target_ptr = Libderiv->deriv_classes[8][5][8];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,45,1,dvrr_stack+87848, dvrr_stack+46306, dvrr_stack+39656);
 tmp = dvrr_stack + 87848;
 target_ptr = Libderiv->deriv_classes[8][6][8];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,45,1,dvrr_stack+80000, dvrr_stack+124456, dvrr_stack+42470);
 tmp = dvrr_stack + 80000;
 target_ptr = Libderiv->deriv_classes[8][7][8];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,45,1,dvrr_stack+61836, dvrr_stack+57226, dvrr_stack+46306);
 tmp = dvrr_stack + 61836;
 target_ptr = Libderiv->deriv_classes[8][8][8];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+100152, dvrr_stack+2566, dvrr_stack+164);
 tmp = dvrr_stack + 100152;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+89108, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 89108;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+89423, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 89423;
 target_ptr = Libderiv->deriv_classes[4][6][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,15,1,dvrr_stack+68892, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 68892;
 target_ptr = Libderiv->deriv_classes[4][7][7];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,15,1,dvrr_stack+63861, dvrr_stack+16752, dvrr_stack+7957);
 tmp = dvrr_stack + 63861;
 target_ptr = Libderiv->deriv_classes[4][8][7];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+69432, dvrr_stack+20817, dvrr_stack+137762);
 tmp = dvrr_stack + 69432;
 target_ptr = Libderiv->deriv_classes[5][4][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+64536, dvrr_stack+23071, dvrr_stack+19977);
 tmp = dvrr_stack + 64536;
 target_ptr = Libderiv->deriv_classes[5][5][7];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+55462, dvrr_stack+26098, dvrr_stack+20817);
 tmp = dvrr_stack + 55462;
 target_ptr = Libderiv->deriv_classes[5][6][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,21,1,dvrr_stack+56050, dvrr_stack+30013, dvrr_stack+23071);
 tmp = dvrr_stack + 56050;
 target_ptr = Libderiv->deriv_classes[5][7][7];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,21,1,dvrr_stack+48970, dvrr_stack+34931, dvrr_stack+26098);
 tmp = dvrr_stack + 48970;
 target_ptr = Libderiv->deriv_classes[5][8][7];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+56806, dvrr_stack+40622, dvrr_stack+90116);
 tmp = dvrr_stack + 56806;
 target_ptr = Libderiv->deriv_classes[6][4][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+49915, dvrr_stack+43758, dvrr_stack+39236);
 tmp = dvrr_stack + 49915;
 target_ptr = Libderiv->deriv_classes[6][5][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+50503, dvrr_stack+47962, dvrr_stack+40622);
 tmp = dvrr_stack + 50503;
 target_ptr = Libderiv->deriv_classes[6][6][7];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,28,1,dvrr_stack+36506, dvrr_stack+52942, dvrr_stack+43758);
 tmp = dvrr_stack + 36506;
 target_ptr = Libderiv->deriv_classes[6][7][7];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,28,1,dvrr_stack+37514, dvrr_stack+59756, dvrr_stack+47962);
 tmp = dvrr_stack + 37514;
 target_ptr = Libderiv->deriv_classes[6][8][7];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,36,1,dvrr_stack+45326, dvrr_stack+67380, dvrr_stack+90396);
 tmp = dvrr_stack + 45326;
 target_ptr = Libderiv->deriv_classes[7][4][7];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,36,1,dvrr_stack+30958, dvrr_stack+71828, dvrr_stack+65496);
 tmp = dvrr_stack + 30958;
 target_ptr = Libderiv->deriv_classes[7][5][7];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,36,1,dvrr_stack+31714, dvrr_stack+77408, dvrr_stack+67380);
 tmp = dvrr_stack + 31714;
 target_ptr = Libderiv->deriv_classes[7][6][7];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,36,1,dvrr_stack+17892, dvrr_stack+84608, dvrr_stack+71828);
 tmp = dvrr_stack + 17892;
 target_ptr = Libderiv->deriv_classes[7][7][7];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,36,1,dvrr_stack+108855, dvrr_stack+93636, dvrr_stack+77408);
 tmp = dvrr_stack + 108855;
 target_ptr = Libderiv->deriv_classes[7][8][7];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,45,1,dvrr_stack+24247, dvrr_stack+39656, dvrr_stack+90756);
 tmp = dvrr_stack + 24247;
 target_ptr = Libderiv->deriv_classes[8][4][7];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,45,1,dvrr_stack+110475, dvrr_stack+42470, dvrr_stack+101751);
 tmp = dvrr_stack + 110475;
 target_ptr = Libderiv->deriv_classes[8][5][7];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,45,1,dvrr_stack+149177, dvrr_stack+46306, dvrr_stack+39656);
 tmp = dvrr_stack + 149177;
 target_ptr = Libderiv->deriv_classes[8][6][7];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,45,1,dvrr_stack+150437, dvrr_stack+124456, dvrr_stack+42470);
 tmp = dvrr_stack + 150437;
 target_ptr = Libderiv->deriv_classes[8][7][7];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,45,1,dvrr_stack+152057, dvrr_stack+57226, dvrr_stack+46306);
 tmp = dvrr_stack + 152057;
 target_ptr = Libderiv->deriv_classes[8][8][7];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+89843, dvrr_stack+2566, dvrr_stack+164);
 tmp = dvrr_stack + 89843;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+45866, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 45866;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+41798, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 41798;
 target_ptr = Libderiv->deriv_classes[4][6][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,15,1,dvrr_stack+19188, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 19188;
 target_ptr = Libderiv->deriv_classes[4][7][6];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,15,1,dvrr_stack+9037, dvrr_stack+16752, dvrr_stack+7957);
 tmp = dvrr_stack + 9037;
 target_ptr = Libderiv->deriv_classes[4][8][6];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+32722, dvrr_stack+20817, dvrr_stack+137762);
 tmp = dvrr_stack + 32722;
 target_ptr = Libderiv->deriv_classes[5][4][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+21699, dvrr_stack+23071, dvrr_stack+19977);
 tmp = dvrr_stack + 21699;
 target_ptr = Libderiv->deriv_classes[5][5][6];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+16588, dvrr_stack+26098, dvrr_stack+20817);
 tmp = dvrr_stack + 16588;
 target_ptr = Libderiv->deriv_classes[5][6][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,21,1,dvrr_stack+154082, dvrr_stack+30013, dvrr_stack+23071);
 tmp = dvrr_stack + 154082;
 target_ptr = Libderiv->deriv_classes[5][7][6];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,21,1,dvrr_stack+154838, dvrr_stack+34931, dvrr_stack+26098);
 tmp = dvrr_stack + 154838;
 target_ptr = Libderiv->deriv_classes[5][8][6];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+5740, dvrr_stack+40622, dvrr_stack+90116);
 tmp = dvrr_stack + 5740;
 target_ptr = Libderiv->deriv_classes[6][4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+34711, dvrr_stack+43758, dvrr_stack+39236);
 tmp = dvrr_stack + 34711;
 target_ptr = Libderiv->deriv_classes[6][5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+35299, dvrr_stack+47962, dvrr_stack+40622);
 tmp = dvrr_stack + 35299;
 target_ptr = Libderiv->deriv_classes[6][6][6];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,28,1,dvrr_stack+155783, dvrr_stack+52942, dvrr_stack+43758);
 tmp = dvrr_stack + 155783;
 target_ptr = Libderiv->deriv_classes[6][7][6];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,28,1,dvrr_stack+156791, dvrr_stack+59756, dvrr_stack+47962);
 tmp = dvrr_stack + 156791;
 target_ptr = Libderiv->deriv_classes[6][8][6];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,36,1,dvrr_stack+158051, dvrr_stack+67380, dvrr_stack+90396);
 tmp = dvrr_stack + 158051;
 target_ptr = Libderiv->deriv_classes[7][4][6];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,36,1,dvrr_stack+59701, dvrr_stack+71828, dvrr_stack+65496);
 tmp = dvrr_stack + 59701;
 target_ptr = Libderiv->deriv_classes[7][5][6];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,36,1,dvrr_stack+158591, dvrr_stack+77408, dvrr_stack+67380);
 tmp = dvrr_stack + 158591;
 target_ptr = Libderiv->deriv_classes[7][6][6];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,36,1,dvrr_stack+159599, dvrr_stack+84608, dvrr_stack+71828);
 tmp = dvrr_stack + 159599;
 target_ptr = Libderiv->deriv_classes[7][7][6];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,36,1,dvrr_stack+160895, dvrr_stack+93636, dvrr_stack+77408);
 tmp = dvrr_stack + 160895;
 target_ptr = Libderiv->deriv_classes[7][8][6];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,45,1,dvrr_stack+90068, dvrr_stack+39656, dvrr_stack+90756);
 tmp = dvrr_stack + 90068;
 target_ptr = Libderiv->deriv_classes[8][4][6];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,45,1,dvrr_stack+93571, dvrr_stack+42470, dvrr_stack+101751);
 tmp = dvrr_stack + 93571;
 target_ptr = Libderiv->deriv_classes[8][5][6];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,45,1,dvrr_stack+162515, dvrr_stack+46306, dvrr_stack+39656);
 tmp = dvrr_stack + 162515;
 target_ptr = Libderiv->deriv_classes[8][6][6];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,45,1,dvrr_stack+163775, dvrr_stack+124456, dvrr_stack+42470);
 tmp = dvrr_stack + 163775;
 target_ptr = Libderiv->deriv_classes[8][7][6];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,45,1,dvrr_stack+165395, dvrr_stack+57226, dvrr_stack+46306);
 tmp = dvrr_stack + 165395;
 target_ptr = Libderiv->deriv_classes[8][8][6];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+57226, dvrr_stack+19977, dvrr_stack+1183);
 tmp = dvrr_stack + 57226;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+57451, dvrr_stack+20817, dvrr_stack+2356);
 tmp = dvrr_stack + 57451;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+57766, dvrr_stack+23071, dvrr_stack+4620);
 tmp = dvrr_stack + 57766;
 target_ptr = Libderiv->deriv_classes[4][6][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,36,dvrr_stack+58186, dvrr_stack+26098, dvrr_stack+7597);
 tmp = dvrr_stack + 58186;
 target_ptr = Libderiv->deriv_classes[4][7][2];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,45,dvrr_stack+58726, dvrr_stack+30013, dvrr_stack+11422);
 tmp = dvrr_stack + 58726;
 target_ptr = Libderiv->deriv_classes[4][8][2];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+90743, dvrr_stack+39236, dvrr_stack+1333);
 tmp = dvrr_stack + 90743;
 target_ptr = Libderiv->deriv_classes[5][4][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+94516, dvrr_stack+40622, dvrr_stack+2566);
 tmp = dvrr_stack + 94516;
 target_ptr = Libderiv->deriv_classes[5][5][2];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+94957, dvrr_stack+43758, dvrr_stack+4900);
 tmp = dvrr_stack + 94957;
 target_ptr = Libderiv->deriv_classes[5][6][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,36,dvrr_stack+60457, dvrr_stack+47962, dvrr_stack+7957);
 tmp = dvrr_stack + 60457;
 target_ptr = Libderiv->deriv_classes[5][7][2];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,45,dvrr_stack+167420, dvrr_stack+52942, dvrr_stack+11872);
 tmp = dvrr_stack + 167420;
 target_ptr = Libderiv->deriv_classes[5][8][2];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+168365, dvrr_stack+65496, dvrr_stack+19977);
 tmp = dvrr_stack + 168365;
 target_ptr = Libderiv->deriv_classes[6][4][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+168785, dvrr_stack+67380, dvrr_stack+20817);
 tmp = dvrr_stack + 168785;
 target_ptr = Libderiv->deriv_classes[6][5][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+169373, dvrr_stack+71828, dvrr_stack+23071);
 tmp = dvrr_stack + 169373;
 target_ptr = Libderiv->deriv_classes[6][6][2];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,36,dvrr_stack+170157, dvrr_stack+77408, dvrr_stack+26098);
 tmp = dvrr_stack + 170157;
 target_ptr = Libderiv->deriv_classes[6][7][2];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,45,dvrr_stack+171165, dvrr_stack+84608, dvrr_stack+30013);
 tmp = dvrr_stack + 171165;
 target_ptr = Libderiv->deriv_classes[6][8][2];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,15,dvrr_stack+172425, dvrr_stack+101751, dvrr_stack+39236);
 tmp = dvrr_stack + 172425;
 target_ptr = Libderiv->deriv_classes[7][4][2];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,21,dvrr_stack+172965, dvrr_stack+39656, dvrr_stack+40622);
 tmp = dvrr_stack + 172965;
 target_ptr = Libderiv->deriv_classes[7][5][2];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,28,dvrr_stack+173721, dvrr_stack+42470, dvrr_stack+43758);
 tmp = dvrr_stack + 173721;
 target_ptr = Libderiv->deriv_classes[7][6][2];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,36,dvrr_stack+174729, dvrr_stack+46306, dvrr_stack+47962);
 tmp = dvrr_stack + 174729;
 target_ptr = Libderiv->deriv_classes[7][7][2];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,45,dvrr_stack+176025, dvrr_stack+124456, dvrr_stack+52942);
 tmp = dvrr_stack + 176025;
 target_ptr = Libderiv->deriv_classes[7][8][2];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_l(Data,15,dvrr_stack+177645, dvrr_stack+91206, dvrr_stack+65496);
 tmp = dvrr_stack + 177645;
 target_ptr = Libderiv->deriv_classes[8][4][2];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_l(Data,21,dvrr_stack+178320, dvrr_stack+102426, dvrr_stack+67380);
 tmp = dvrr_stack + 178320;
 target_ptr = Libderiv->deriv_classes[8][5][2];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_l(Data,28,dvrr_stack+179265, dvrr_stack+92031, dvrr_stack+71828);
 tmp = dvrr_stack + 179265;
 target_ptr = Libderiv->deriv_classes[8][6][2];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_l(Data,36,dvrr_stack+180525, dvrr_stack+14167, dvrr_stack+77408);
 tmp = dvrr_stack + 180525;
 target_ptr = Libderiv->deriv_classes[8][7][2];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_l(Data,45,dvrr_stack+182145, dvrr_stack+81728, dvrr_stack+84608);
 tmp = dvrr_stack + 182145;
 target_ptr = Libderiv->deriv_classes[8][8][2];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+59401, dvrr_stack+19977, dvrr_stack+1183);
 tmp = dvrr_stack + 59401;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+3196, dvrr_stack+20817, dvrr_stack+2356);
 tmp = dvrr_stack + 3196;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+184170, dvrr_stack+23071, dvrr_stack+4620);
 tmp = dvrr_stack + 184170;
 target_ptr = Libderiv->deriv_classes[4][6][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,36,dvrr_stack+184590, dvrr_stack+26098, dvrr_stack+7597);
 tmp = dvrr_stack + 184590;
 target_ptr = Libderiv->deriv_classes[4][7][1];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,45,dvrr_stack+185130, dvrr_stack+30013, dvrr_stack+11422);
 tmp = dvrr_stack + 185130;
 target_ptr = Libderiv->deriv_classes[4][8][1];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+111420, dvrr_stack+39236, dvrr_stack+1333);
 tmp = dvrr_stack + 111420;
 target_ptr = Libderiv->deriv_classes[5][4][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+185805, dvrr_stack+40622, dvrr_stack+2566);
 tmp = dvrr_stack + 185805;
 target_ptr = Libderiv->deriv_classes[5][5][1];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+186246, dvrr_stack+43758, dvrr_stack+4900);
 tmp = dvrr_stack + 186246;
 target_ptr = Libderiv->deriv_classes[5][6][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,36,dvrr_stack+186834, dvrr_stack+47962, dvrr_stack+7957);
 tmp = dvrr_stack + 186834;
 target_ptr = Libderiv->deriv_classes[5][7][1];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,45,dvrr_stack+187590, dvrr_stack+52942, dvrr_stack+11872);
 tmp = dvrr_stack + 187590;
 target_ptr = Libderiv->deriv_classes[5][8][1];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+188535, dvrr_stack+65496, dvrr_stack+19977);
 tmp = dvrr_stack + 188535;
 target_ptr = Libderiv->deriv_classes[6][4][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+188955, dvrr_stack+67380, dvrr_stack+20817);
 tmp = dvrr_stack + 188955;
 target_ptr = Libderiv->deriv_classes[6][5][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+189543, dvrr_stack+71828, dvrr_stack+23071);
 tmp = dvrr_stack + 189543;
 target_ptr = Libderiv->deriv_classes[6][6][1];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,36,dvrr_stack+190327, dvrr_stack+77408, dvrr_stack+26098);
 tmp = dvrr_stack + 190327;
 target_ptr = Libderiv->deriv_classes[6][7][1];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,45,dvrr_stack+191335, dvrr_stack+84608, dvrr_stack+30013);
 tmp = dvrr_stack + 191335;
 target_ptr = Libderiv->deriv_classes[6][8][1];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,15,dvrr_stack+192595, dvrr_stack+101751, dvrr_stack+39236);
 tmp = dvrr_stack + 192595;
 target_ptr = Libderiv->deriv_classes[7][4][1];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,21,dvrr_stack+193135, dvrr_stack+39656, dvrr_stack+40622);
 tmp = dvrr_stack + 193135;
 target_ptr = Libderiv->deriv_classes[7][5][1];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,28,dvrr_stack+193891, dvrr_stack+42470, dvrr_stack+43758);
 tmp = dvrr_stack + 193891;
 target_ptr = Libderiv->deriv_classes[7][6][1];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,36,dvrr_stack+194899, dvrr_stack+46306, dvrr_stack+47962);
 tmp = dvrr_stack + 194899;
 target_ptr = Libderiv->deriv_classes[7][7][1];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,45,dvrr_stack+196195, dvrr_stack+124456, dvrr_stack+52942);
 tmp = dvrr_stack + 196195;
 target_ptr = Libderiv->deriv_classes[7][8][1];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_l(Data,15,dvrr_stack+197815, dvrr_stack+91206, dvrr_stack+65496);
 tmp = dvrr_stack + 197815;
 target_ptr = Libderiv->deriv_classes[8][4][1];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_l(Data,21,dvrr_stack+198490, dvrr_stack+102426, dvrr_stack+67380);
 tmp = dvrr_stack + 198490;
 target_ptr = Libderiv->deriv_classes[8][5][1];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_l(Data,28,dvrr_stack+199435, dvrr_stack+92031, dvrr_stack+71828);
 tmp = dvrr_stack + 199435;
 target_ptr = Libderiv->deriv_classes[8][6][1];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_l(Data,36,dvrr_stack+200695, dvrr_stack+14167, dvrr_stack+77408);
 tmp = dvrr_stack + 200695;
 target_ptr = Libderiv->deriv_classes[8][7][1];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_l(Data,45,dvrr_stack+202315, dvrr_stack+81728, dvrr_stack+84608);
 tmp = dvrr_stack + 202315;
 target_ptr = Libderiv->deriv_classes[8][8][1];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+137741, dvrr_stack+19977, dvrr_stack+1183);
 tmp = dvrr_stack + 137741;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+17176, dvrr_stack+20817, dvrr_stack+2356);
 tmp = dvrr_stack + 17176;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+854, dvrr_stack+23071, dvrr_stack+4620);
 tmp = dvrr_stack + 854;
 target_ptr = Libderiv->deriv_classes[4][6][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,36,dvrr_stack+204340, dvrr_stack+26098, dvrr_stack+7597);
 tmp = dvrr_stack + 204340;
 target_ptr = Libderiv->deriv_classes[4][7][0];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,45,dvrr_stack+204880, dvrr_stack+30013, dvrr_stack+11422);
 tmp = dvrr_stack + 204880;
 target_ptr = Libderiv->deriv_classes[4][8][0];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+4501, dvrr_stack+39236, dvrr_stack+1333);
 tmp = dvrr_stack + 4501;
 target_ptr = Libderiv->deriv_classes[5][4][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+11129, dvrr_stack+40622, dvrr_stack+2566);
 tmp = dvrr_stack + 11129;
 target_ptr = Libderiv->deriv_classes[5][5][0];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+205555, dvrr_stack+43758, dvrr_stack+4900);
 tmp = dvrr_stack + 205555;
 target_ptr = Libderiv->deriv_classes[5][6][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,36,dvrr_stack+206143, dvrr_stack+47962, dvrr_stack+7957);
 tmp = dvrr_stack + 206143;
 target_ptr = Libderiv->deriv_classes[5][7][0];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,45,dvrr_stack+7525, dvrr_stack+52942, dvrr_stack+11872);
 tmp = dvrr_stack + 7525;
 target_ptr = Libderiv->deriv_classes[5][8][0];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+2314, dvrr_stack+65496, dvrr_stack+19977);
 tmp = dvrr_stack + 2314;
 target_ptr = Libderiv->deriv_classes[6][4][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+11570, dvrr_stack+67380, dvrr_stack+20817);
 tmp = dvrr_stack + 11570;
 target_ptr = Libderiv->deriv_classes[6][5][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+206899, dvrr_stack+71828, dvrr_stack+23071);
 tmp = dvrr_stack + 206899;
 target_ptr = Libderiv->deriv_classes[6][6][0];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,36,dvrr_stack+207683, dvrr_stack+77408, dvrr_stack+26098);
 tmp = dvrr_stack + 207683;
 target_ptr = Libderiv->deriv_classes[6][7][0];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,45,dvrr_stack+208691, dvrr_stack+84608, dvrr_stack+30013);
 tmp = dvrr_stack + 208691;
 target_ptr = Libderiv->deriv_classes[6][8][0];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,15,dvrr_stack+19728, dvrr_stack+101751, dvrr_stack+39236);
 tmp = dvrr_stack + 19728;
 target_ptr = Libderiv->deriv_classes[7][4][0];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,21,dvrr_stack+101571, dvrr_stack+39656, dvrr_stack+40622);
 tmp = dvrr_stack + 101571;
 target_ptr = Libderiv->deriv_classes[7][5][0];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,28,dvrr_stack+39236, dvrr_stack+42470, dvrr_stack+43758);
 tmp = dvrr_stack + 39236;
 target_ptr = Libderiv->deriv_classes[7][6][0];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,36,dvrr_stack+42218, dvrr_stack+46306, dvrr_stack+47962);
 tmp = dvrr_stack + 42218;
 target_ptr = Libderiv->deriv_classes[7][7][0];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,45,dvrr_stack+46181, dvrr_stack+124456, dvrr_stack+52942);
 tmp = dvrr_stack + 46181;
 target_ptr = Libderiv->deriv_classes[7][8][0];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_l(Data,15,dvrr_stack+40244, dvrr_stack+91206, dvrr_stack+65496);
 tmp = dvrr_stack + 40244;
 target_ptr = Libderiv->deriv_classes[8][4][0];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_l(Data,21,dvrr_stack+91058, dvrr_stack+102426, dvrr_stack+67380);
 tmp = dvrr_stack + 91058;
 target_ptr = Libderiv->deriv_classes[8][5][0];
 for(i=0;i<945;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_l(Data,28,dvrr_stack+52762, dvrr_stack+92031, dvrr_stack+71828);
 tmp = dvrr_stack + 52762;
 target_ptr = Libderiv->deriv_classes[8][6][0];
 for(i=0;i<1260;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_l(Data,36,dvrr_stack+124406, dvrr_stack+14167, dvrr_stack+77408);
 tmp = dvrr_stack + 124406;
 target_ptr = Libderiv->deriv_classes[8][7][0];
 for(i=0;i<1620;i++)
   target_ptr[i] += tmp[i];

 /* compute (8 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_l(Data,45,dvrr_stack+209951, dvrr_stack+81728, dvrr_stack+84608);
 tmp = dvrr_stack + 209951;
 target_ptr = Libderiv->deriv_classes[8][8][0];
 for(i=0;i<2025;i++)
   target_ptr[i] += tmp[i];


}

