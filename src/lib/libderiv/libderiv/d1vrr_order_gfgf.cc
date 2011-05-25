#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (gf|gf) integrals */

void d1vrr_order_gfgf(Libderiv_t *Libderiv, prim_data *Data)
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
 vrr_build_xxxx(am,Data,dvrr_stack+6314, dvrr_stack+1747, dvrr_stack+6293, dvrr_stack+723, dvrr_stack+3575, NULL);

 /* compute (0 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4052, dvrr_stack+3590, dvrr_stack+6314, dvrr_stack+1936, dvrr_stack+1747, NULL);

 /* compute (0 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6634, dvrr_stack+738, dvrr_stack+4052, dvrr_stack+1574, dvrr_stack+3590, NULL);

 /* compute (1 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10567, dvrr_stack+6679, dvrr_stack+6634, NULL, NULL, dvrr_stack+738);

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

 /* compute (4 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+12547,dvrr_stack+11872,dvrr_stack+7957,15);


 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+611, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+864, dvrr_stack+0, dvrr_stack+30, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+9757, dvrr_stack+6, dvrr_stack+864, dvrr_stack+3, dvrr_stack+0, dvrr_stack+611);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+9775, dvrr_stack+33, dvrr_stack+213, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+9793, dvrr_stack+39, dvrr_stack+9775, dvrr_stack+15, dvrr_stack+33, dvrr_stack+864);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+9829, dvrr_stack+75, dvrr_stack+9793, dvrr_stack+57, dvrr_stack+39, dvrr_stack+9757);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+479, dvrr_stack+219, dvrr_stack+873, NULL, NULL, dvrr_stack+213);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+10162, dvrr_stack+229, dvrr_stack+479, dvrr_stack+111, dvrr_stack+219, dvrr_stack+9775);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+3808, dvrr_stack+259, dvrr_stack+10162, dvrr_stack+131, dvrr_stack+229, dvrr_stack+9793);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10222, dvrr_stack+379, dvrr_stack+3808, dvrr_stack+319, dvrr_stack+259, dvrr_stack+9829);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6589, dvrr_stack+883, dvrr_stack+575, NULL, NULL, dvrr_stack+873);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10372, dvrr_stack+898, dvrr_stack+6589, dvrr_stack+509, dvrr_stack+883, dvrr_stack+479);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14167, dvrr_stack+943, dvrr_stack+10372, dvrr_stack+524, dvrr_stack+898, dvrr_stack+10162);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14317, dvrr_stack+1033, dvrr_stack+14167, dvrr_stack+614, dvrr_stack+943, dvrr_stack+3808);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14542, dvrr_stack+1333, dvrr_stack+14317, dvrr_stack+1183, dvrr_stack+1033, dvrr_stack+10222);
 tmp = dvrr_stack + 14542;
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
 vrr_build_xxxx(am,Data,dvrr_stack+14857, dvrr_stack+1957, dvrr_stack+614, dvrr_stack+780, dvrr_stack+590, dvrr_stack+6589);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+14983, dvrr_stack+2020, dvrr_stack+14857, dvrr_stack+801, dvrr_stack+1957, dvrr_stack+10372);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15193, dvrr_stack+2146, dvrr_stack+14983, dvrr_stack+1621, dvrr_stack+2020, dvrr_stack+14167);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15508, dvrr_stack+2566, dvrr_stack+15193, dvrr_stack+2356, dvrr_stack+2146, dvrr_stack+14317);
 tmp = dvrr_stack + 15508;
 target_ptr = Libderiv->dvrr_classes[5][5];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+15949,dvrr_stack+15508,dvrr_stack+14542,21);


 /* compute (1 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+780, dvrr_stack+1574, dvrr_stack+3590, NULL, NULL, dvrr_stack+1936);

 /* compute (2 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16894, dvrr_stack+4088, dvrr_stack+780, dvrr_stack+1824, dvrr_stack+1574, dvrr_stack+614);

 /* compute (3 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+17062, dvrr_stack+4172, dvrr_stack+16894, dvrr_stack+1852, dvrr_stack+4088, dvrr_stack+14857);

 /* compute (4 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+17342, dvrr_stack+4340, dvrr_stack+17062, dvrr_stack+3640, dvrr_stack+4172, dvrr_stack+14983);

 /* compute (5 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+17762, dvrr_stack+4900, dvrr_stack+17342, dvrr_stack+4620, dvrr_stack+4340, dvrr_stack+15193);
 tmp = dvrr_stack + 17762;
 target_ptr = Libderiv->dvrr_classes[5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+18350,dvrr_stack+17762,dvrr_stack+15508,21);


 /* compute (1 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1768, dvrr_stack+738, dvrr_stack+4052, NULL, NULL, dvrr_stack+3590);

 /* compute (2 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+19673, dvrr_stack+6913, dvrr_stack+1768, dvrr_stack+3908, dvrr_stack+738, dvrr_stack+780);

 /* compute (3 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+19889, dvrr_stack+7021, dvrr_stack+19673, dvrr_stack+3944, dvrr_stack+6913, dvrr_stack+16894);

 /* compute (4 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+20249, dvrr_stack+7237, dvrr_stack+19889, dvrr_stack+6373, dvrr_stack+7021, dvrr_stack+17062);

 /* compute (5 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+20789, dvrr_stack+7957, dvrr_stack+20249, dvrr_stack+7597, dvrr_stack+7237, dvrr_stack+17342);
 tmp = dvrr_stack + 20789;
 target_ptr = Libderiv->dvrr_classes[5][7];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+21545,dvrr_stack+20789,dvrr_stack+17762,21);


 /* compute (0 0 | 1 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+9889, Data->F+12, Data->F+13, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+6265, dvrr_stack+9889, Data->F+11, Data->F+12, NULL);

 /* compute (0 0 | 3 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+774, dvrr_stack+15, dvrr_stack+3556, dvrr_stack+6265, NULL);

 /* compute (0 0 | 4 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3908, dvrr_stack+6268, dvrr_stack+111, dvrr_stack+3559, dvrr_stack+774, NULL);

 /* compute (0 0 | 5 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3923, dvrr_stack+6278, dvrr_stack+3908, dvrr_stack+3565, dvrr_stack+6268, NULL);

 /* compute (0 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3944, dvrr_stack+6293, dvrr_stack+3923, dvrr_stack+3575, dvrr_stack+6278, NULL);

 /* compute (0 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3972, dvrr_stack+6314, dvrr_stack+3944, dvrr_stack+1747, dvrr_stack+6293, NULL);

 /* compute (0 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+509, dvrr_stack+4052, dvrr_stack+3972, dvrr_stack+3590, dvrr_stack+6314, NULL);

 /* compute (1 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6342, dvrr_stack+6634, dvrr_stack+509, NULL, NULL, dvrr_stack+4052);

 /* compute (2 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+23309, dvrr_stack+10567, dvrr_stack+6342, dvrr_stack+6679, dvrr_stack+6634, dvrr_stack+1768);

 /* compute (3 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+23579, dvrr_stack+10702, dvrr_stack+23309, dvrr_stack+6724, dvrr_stack+10567, dvrr_stack+19673);

 /* compute (4 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+24029, dvrr_stack+10972, dvrr_stack+23579, dvrr_stack+9892, dvrr_stack+10702, dvrr_stack+19889);

 /* compute (5 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+24704, dvrr_stack+11872, dvrr_stack+24029, dvrr_stack+11422, dvrr_stack+10972, dvrr_stack+20249);

 /* compute (5 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+25649,dvrr_stack+24704,dvrr_stack+20789,21);


 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+11422, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+11425, dvrr_stack+611, dvrr_stack+11422, Data->F+4, Data->F+5, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3556, dvrr_stack+30, dvrr_stack+210, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+11431, dvrr_stack+864, dvrr_stack+3556, dvrr_stack+0, dvrr_stack+30, dvrr_stack+11422);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+131, dvrr_stack+9757, dvrr_stack+11431, dvrr_stack+6, dvrr_stack+864, dvrr_stack+11425);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+11449, dvrr_stack+213, dvrr_stack+704, NULL, NULL, dvrr_stack+210);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+11467, dvrr_stack+9775, dvrr_stack+11449, dvrr_stack+33, dvrr_stack+213, dvrr_stack+3556);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+1876, dvrr_stack+9793, dvrr_stack+11467, dvrr_stack+39, dvrr_stack+9775, dvrr_stack+11431);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+11503, dvrr_stack+9829, dvrr_stack+1876, dvrr_stack+75, dvrr_stack+9793, dvrr_stack+131);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+11593, dvrr_stack+873, dvrr_stack+121, NULL, NULL, dvrr_stack+704);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+11623, dvrr_stack+479, dvrr_stack+11593, dvrr_stack+219, dvrr_stack+873, dvrr_stack+11449);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+11683, dvrr_stack+10162, dvrr_stack+11623, dvrr_stack+229, dvrr_stack+479, dvrr_stack+11467);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+9892, dvrr_stack+3808, dvrr_stack+11683, dvrr_stack+259, dvrr_stack+10162, dvrr_stack+1876);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+6679, dvrr_stack+10222, dvrr_stack+9892, dvrr_stack+379, dvrr_stack+3808, dvrr_stack+11503);

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+219, dvrr_stack+575, dvrr_stack+723, NULL, NULL, dvrr_stack+121);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10042, dvrr_stack+6589, dvrr_stack+219, dvrr_stack+883, dvrr_stack+575, dvrr_stack+11593);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3618, dvrr_stack+10372, dvrr_stack+10042, dvrr_stack+898, dvrr_stack+6589, dvrr_stack+11623);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+27917, dvrr_stack+14167, dvrr_stack+3618, dvrr_stack+943, dvrr_stack+10372, dvrr_stack+11683);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+28142, dvrr_stack+14317, dvrr_stack+27917, dvrr_stack+1033, dvrr_stack+14167, dvrr_stack+9892);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+28457, dvrr_stack+14542, dvrr_stack+28142, dvrr_stack+1333, dvrr_stack+14317, dvrr_stack+6679);
 tmp = dvrr_stack + 28457;
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
 vrr_build_xxxx(am,Data,dvrr_stack+28877, dvrr_stack+14857, dvrr_stack+946, dvrr_stack+1957, dvrr_stack+614, dvrr_stack+10042);

 /* compute (4 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+29087, dvrr_stack+14983, dvrr_stack+28877, dvrr_stack+2020, dvrr_stack+14857, dvrr_stack+3618);

 /* compute (5 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+29402, dvrr_stack+15193, dvrr_stack+29087, dvrr_stack+2146, dvrr_stack+14983, dvrr_stack+27917);

 /* compute (6 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+29843, dvrr_stack+15508, dvrr_stack+29402, dvrr_stack+2566, dvrr_stack+15193, dvrr_stack+28142);
 tmp = dvrr_stack + 29843;
 target_ptr = Libderiv->dvrr_classes[6][5];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+30431,dvrr_stack+29843,dvrr_stack+28457,28);


 /* compute (1 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1957, dvrr_stack+3590, dvrr_stack+6314, NULL, NULL, dvrr_stack+1747);

 /* compute (2 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2041, dvrr_stack+780, dvrr_stack+1957, dvrr_stack+1574, dvrr_stack+3590, dvrr_stack+883);

 /* compute (3 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+31691, dvrr_stack+16894, dvrr_stack+2041, dvrr_stack+4088, dvrr_stack+780, dvrr_stack+946);

 /* compute (4 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+31971, dvrr_stack+17062, dvrr_stack+31691, dvrr_stack+4172, dvrr_stack+16894, dvrr_stack+28877);

 /* compute (5 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+32391, dvrr_stack+17342, dvrr_stack+31971, dvrr_stack+4340, dvrr_stack+17062, dvrr_stack+29087);

 /* compute (6 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+32979, dvrr_stack+17762, dvrr_stack+32391, dvrr_stack+4900, dvrr_stack+17342, dvrr_stack+29402);
 tmp = dvrr_stack + 32979;
 target_ptr = Libderiv->dvrr_classes[6][6];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+33763,dvrr_stack+32979,dvrr_stack+29843,28);


 /* compute (1 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4088, dvrr_stack+4052, dvrr_stack+3972, NULL, NULL, dvrr_stack+6314);

 /* compute (2 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4196, dvrr_stack+1768, dvrr_stack+4088, dvrr_stack+738, dvrr_stack+4052, dvrr_stack+1957);

 /* compute (3 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+35527, dvrr_stack+19673, dvrr_stack+4196, dvrr_stack+6913, dvrr_stack+1768, dvrr_stack+2041);

 /* compute (4 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+35887, dvrr_stack+19889, dvrr_stack+35527, dvrr_stack+7021, dvrr_stack+19673, dvrr_stack+31691);

 /* compute (5 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+36427, dvrr_stack+20249, dvrr_stack+35887, dvrr_stack+7237, dvrr_stack+19889, dvrr_stack+31971);

 /* compute (6 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+37183, dvrr_stack+20789, dvrr_stack+36427, dvrr_stack+7957, dvrr_stack+20249, dvrr_stack+32391);
 tmp = dvrr_stack + 37183;
 target_ptr = Libderiv->dvrr_classes[6][7];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+38191,dvrr_stack+37183,dvrr_stack+32979,28);


 /* compute (0 0 | 1 0) m=13 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+13, Data->F+14, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+738, dvrr_stack+9889, dvrr_stack+0, Data->F+12, Data->F+13, NULL);

 /* compute (0 0 | 3 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+744, dvrr_stack+15, dvrr_stack+738, dvrr_stack+6265, dvrr_stack+9889, NULL);

 /* compute (0 0 | 4 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+554, dvrr_stack+111, dvrr_stack+744, dvrr_stack+774, dvrr_stack+15, NULL);

 /* compute (0 0 | 5 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+590, dvrr_stack+3908, dvrr_stack+554, dvrr_stack+6268, dvrr_stack+111, NULL);

 /* compute (0 0 | 6 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4412, dvrr_stack+3923, dvrr_stack+590, dvrr_stack+6278, dvrr_stack+3908, NULL);

 /* compute (0 0 | 7 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4440, dvrr_stack+3944, dvrr_stack+4412, dvrr_stack+6293, dvrr_stack+3923, NULL);

 /* compute (0 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+4476, dvrr_stack+3972, dvrr_stack+4440, dvrr_stack+6314, dvrr_stack+3944, NULL);

 /* compute (1 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1574, dvrr_stack+509, dvrr_stack+4476, NULL, NULL, dvrr_stack+3972);

 /* compute (2 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6889, dvrr_stack+6342, dvrr_stack+1574, dvrr_stack+6634, dvrr_stack+509, dvrr_stack+4088);

 /* compute (3 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+40543, dvrr_stack+23309, dvrr_stack+6889, dvrr_stack+10567, dvrr_stack+6342, dvrr_stack+4196);

 /* compute (4 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+40993, dvrr_stack+23579, dvrr_stack+40543, dvrr_stack+10702, dvrr_stack+23309, dvrr_stack+35527);

 /* compute (5 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+41668, dvrr_stack+24029, dvrr_stack+40993, dvrr_stack+10972, dvrr_stack+23579, dvrr_stack+35887);

 /* compute (6 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+42613, dvrr_stack+24704, dvrr_stack+41668, dvrr_stack+11872, dvrr_stack+24029, dvrr_stack+36427);

 /* compute (6 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+43873,dvrr_stack+42613,dvrr_stack+37183,28);


 /* compute (1 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+6634, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+6637, dvrr_stack+11422, dvrr_stack+6634, Data->F+5, Data->F+6, NULL);

 /* compute (3 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+6643, dvrr_stack+11425, dvrr_stack+6637, dvrr_stack+611, dvrr_stack+11422, NULL);

 /* compute (1 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+6653, dvrr_stack+210, dvrr_stack+161, NULL, NULL, Data->F+7);

 /* compute (2 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+754, dvrr_stack+3556, dvrr_stack+6653, dvrr_stack+30, dvrr_stack+210, dvrr_stack+6634);

 /* compute (3 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+10132, dvrr_stack+11431, dvrr_stack+754, dvrr_stack+864, dvrr_stack+3556, dvrr_stack+6637);

 /* compute (4 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+4521, dvrr_stack+131, dvrr_stack+10132, dvrr_stack+9757, dvrr_stack+11431, dvrr_stack+6643);

 /* compute (1 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+30, dvrr_stack+704, dvrr_stack+569, NULL, NULL, dvrr_stack+161);

 /* compute (2 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+4566, dvrr_stack+11449, dvrr_stack+30, dvrr_stack+213, dvrr_stack+704, dvrr_stack+6653);

 /* compute (3 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+2209, dvrr_stack+11467, dvrr_stack+4566, dvrr_stack+9775, dvrr_stack+11449, dvrr_stack+754);

 /* compute (4 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+1072, dvrr_stack+1876, dvrr_stack+2209, dvrr_stack+9793, dvrr_stack+11467, dvrr_stack+10132);

 /* compute (5 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+7159, dvrr_stack+11503, dvrr_stack+1072, dvrr_stack+9829, dvrr_stack+1876, dvrr_stack+4521);

 /* compute (1 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+9775, dvrr_stack+121, dvrr_stack+1564, NULL, NULL, dvrr_stack+569);

 /* compute (2 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+2269, dvrr_stack+11593, dvrr_stack+9775, dvrr_stack+873, dvrr_stack+121, dvrr_stack+30);

 /* compute (3 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+7285, dvrr_stack+11623, dvrr_stack+2269, dvrr_stack+479, dvrr_stack+11593, dvrr_stack+4566);

 /* compute (4 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+7385, dvrr_stack+11683, dvrr_stack+7285, dvrr_stack+10162, dvrr_stack+11623, dvrr_stack+2209);

 /* compute (5 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10462, dvrr_stack+9892, dvrr_stack+7385, dvrr_stack+3808, dvrr_stack+11683, dvrr_stack+1072);

 /* compute (6 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10672, dvrr_stack+6679, dvrr_stack+10462, dvrr_stack+10222, dvrr_stack+9892, dvrr_stack+7159);

 /* compute (1 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10162, dvrr_stack+723, dvrr_stack+3575, NULL, NULL, dvrr_stack+1564);

 /* compute (2 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10952, dvrr_stack+219, dvrr_stack+10162, dvrr_stack+575, dvrr_stack+723, dvrr_stack+9775);

 /* compute (3 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11042, dvrr_stack+10042, dvrr_stack+10952, dvrr_stack+6589, dvrr_stack+219, dvrr_stack+2269);

 /* compute (4 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11192, dvrr_stack+3618, dvrr_stack+11042, dvrr_stack+10372, dvrr_stack+10042, dvrr_stack+7285);

 /* compute (5 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+46897, dvrr_stack+27917, dvrr_stack+11192, dvrr_stack+14167, dvrr_stack+3618, dvrr_stack+7385);

 /* compute (6 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+47212, dvrr_stack+28142, dvrr_stack+46897, dvrr_stack+14317, dvrr_stack+27917, dvrr_stack+10462);

 /* compute (7 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+47632, dvrr_stack+28457, dvrr_stack+47212, dvrr_stack+14542, dvrr_stack+28142, dvrr_stack+10672);
 tmp = dvrr_stack + 47632;
 target_ptr = Libderiv->dvrr_classes[7][4];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+14167, dvrr_stack+1747, dvrr_stack+6293, NULL, NULL, dvrr_stack+3575);

 /* compute (2 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+14230, dvrr_stack+883, dvrr_stack+14167, dvrr_stack+1936, dvrr_stack+1747, dvrr_stack+10162);

 /* compute (3 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+48172, dvrr_stack+946, dvrr_stack+14230, dvrr_stack+614, dvrr_stack+883, dvrr_stack+10952);

 /* compute (4 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+48382, dvrr_stack+28877, dvrr_stack+48172, dvrr_stack+14857, dvrr_stack+946, dvrr_stack+11042);

 /* compute (5 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+48697, dvrr_stack+29087, dvrr_stack+48382, dvrr_stack+14983, dvrr_stack+28877, dvrr_stack+11192);

 /* compute (6 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+49138, dvrr_stack+29402, dvrr_stack+48697, dvrr_stack+15193, dvrr_stack+29087, dvrr_stack+46897);

 /* compute (7 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+49726, dvrr_stack+29843, dvrr_stack+49138, dvrr_stack+15508, dvrr_stack+29402, dvrr_stack+47212);
 tmp = dvrr_stack + 49726;
 target_ptr = Libderiv->dvrr_classes[7][5];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+50482,dvrr_stack+49726,dvrr_stack+47632,36);


 /* compute (1 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14857, dvrr_stack+6314, dvrr_stack+3944, NULL, NULL, dvrr_stack+6293);

 /* compute (2 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14941, dvrr_stack+1957, dvrr_stack+14857, dvrr_stack+3590, dvrr_stack+6314, dvrr_stack+14167);

 /* compute (3 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+15109, dvrr_stack+2041, dvrr_stack+14941, dvrr_stack+780, dvrr_stack+1957, dvrr_stack+14230);

 /* compute (4 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+52102, dvrr_stack+31691, dvrr_stack+15109, dvrr_stack+16894, dvrr_stack+2041, dvrr_stack+48172);

 /* compute (5 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+52522, dvrr_stack+31971, dvrr_stack+52102, dvrr_stack+17062, dvrr_stack+31691, dvrr_stack+48382);

 /* compute (6 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+53110, dvrr_stack+32391, dvrr_stack+52522, dvrr_stack+17342, dvrr_stack+31971, dvrr_stack+48697);

 /* compute (7 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+53894, dvrr_stack+32979, dvrr_stack+53110, dvrr_stack+17762, dvrr_stack+32391, dvrr_stack+49138);
 tmp = dvrr_stack + 53894;
 target_ptr = Libderiv->dvrr_classes[7][6];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+54902,dvrr_stack+53894,dvrr_stack+49726,36);


 /* compute (1 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+16894, dvrr_stack+3972, dvrr_stack+4440, NULL, NULL, dvrr_stack+3944);

 /* compute (2 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+17002, dvrr_stack+4088, dvrr_stack+16894, dvrr_stack+4052, dvrr_stack+3972, dvrr_stack+14857);

 /* compute (3 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+17218, dvrr_stack+4196, dvrr_stack+17002, dvrr_stack+1768, dvrr_stack+4088, dvrr_stack+14941);

 /* compute (4 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+57170, dvrr_stack+35527, dvrr_stack+17218, dvrr_stack+19673, dvrr_stack+4196, dvrr_stack+15109);

 /* compute (5 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+57710, dvrr_stack+35887, dvrr_stack+57170, dvrr_stack+19889, dvrr_stack+35527, dvrr_stack+52102);

 /* compute (6 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+58466, dvrr_stack+36427, dvrr_stack+57710, dvrr_stack+20249, dvrr_stack+35887, dvrr_stack+52522);

 /* compute (7 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+59474, dvrr_stack+37183, dvrr_stack+58466, dvrr_stack+20789, dvrr_stack+36427, dvrr_stack+53110);

 /* compute (7 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+60770,dvrr_stack+59474,dvrr_stack+53894,36);


 /* compute (0 0 | 1 0) m=14 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+19673, Data->F+14, Data->F+15, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=13 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+213, dvrr_stack+0, dvrr_stack+19673, Data->F+13, Data->F+14, NULL);

 /* compute (0 0 | 3 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+19673, dvrr_stack+738, dvrr_stack+213, dvrr_stack+9889, dvrr_stack+0, NULL);

 /* compute (0 0 | 4 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+575, dvrr_stack+744, dvrr_stack+19673, dvrr_stack+15, dvrr_stack+738, NULL);

 /* compute (0 0 | 5 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1936, dvrr_stack+554, dvrr_stack+575, dvrr_stack+111, dvrr_stack+744, NULL);

 /* compute (0 0 | 6 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3590, dvrr_stack+590, dvrr_stack+1936, dvrr_stack+3908, dvrr_stack+554, NULL);

 /* compute (0 0 | 7 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+19673, dvrr_stack+4412, dvrr_stack+3590, dvrr_stack+3923, dvrr_stack+590, NULL);

 /* compute (0 0 | 8 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+19709, dvrr_stack+4440, dvrr_stack+19673, dvrr_stack+3944, dvrr_stack+4412, NULL);

 /* compute (1 0 | 8 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+19754, dvrr_stack+4476, dvrr_stack+19709, NULL, NULL, dvrr_stack+4440);

 /* compute (2 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+19889, dvrr_stack+1574, dvrr_stack+19754, dvrr_stack+509, dvrr_stack+4476, dvrr_stack+16894);

 /* compute (3 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+20159, dvrr_stack+6889, dvrr_stack+19889, dvrr_stack+6342, dvrr_stack+1574, dvrr_stack+17002);

 /* compute (4 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+63794, dvrr_stack+40543, dvrr_stack+20159, dvrr_stack+23309, dvrr_stack+6889, dvrr_stack+17218);

 /* compute (5 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+19709, dvrr_stack+40993, dvrr_stack+63794, dvrr_stack+23579, dvrr_stack+40543, dvrr_stack+57170);

 /* compute (6 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+63794, dvrr_stack+41668, dvrr_stack+19709, dvrr_stack+24029, dvrr_stack+40993, dvrr_stack+57710);

 /* compute (7 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+65054, dvrr_stack+42613, dvrr_stack+63794, dvrr_stack+24704, dvrr_stack+41668, dvrr_stack+58466);

 /* compute (7 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+66674,dvrr_stack+65054,dvrr_stack+59474,36);


 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+48, dvrr_stack+21, dvrr_stack+3, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+4602, dvrr_stack+48, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+0);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+63794, dvrr_stack+164, dvrr_stack+24, NULL, NULL, dvrr_stack+21);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+575, dvrr_stack+63794, dvrr_stack+57, dvrr_stack+164, dvrr_stack+24, dvrr_stack+48);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+63812, dvrr_stack+575, dvrr_stack+75, dvrr_stack+63794, dvrr_stack+57, dvrr_stack+4602);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+63872, dvrr_stack+713, dvrr_stack+170, NULL, NULL, dvrr_stack+164);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+63902, dvrr_stack+63872, dvrr_stack+180, dvrr_stack+713, dvrr_stack+170, dvrr_stack+63794);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+63962, dvrr_stack+63902, dvrr_stack+319, dvrr_stack+63872, dvrr_stack+180, dvrr_stack+575);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+64062, dvrr_stack+63962, dvrr_stack+379, dvrr_stack+63902, dvrr_stack+319, dvrr_stack+63812);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+213, dvrr_stack+0, dvrr_stack+611, Data->F+3, Data->F+4, NULL);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+63872, dvrr_stack+4602, dvrr_stack+9757, dvrr_stack+48, dvrr_stack+6, dvrr_stack+213);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+614, dvrr_stack+63812, dvrr_stack+9829, dvrr_stack+575, dvrr_stack+75, dvrr_stack+63872);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+64212, dvrr_stack+64062, dvrr_stack+10222, dvrr_stack+63962, dvrr_stack+379, dvrr_stack+614);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+713, dvrr_stack+213, dvrr_stack+11425, dvrr_stack+0, dvrr_stack+611, NULL);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+4476, dvrr_stack+63872, dvrr_stack+131, dvrr_stack+4602, dvrr_stack+9757, dvrr_stack+713);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+63902, dvrr_stack+614, dvrr_stack+11503, dvrr_stack+63812, dvrr_stack+9829, dvrr_stack+4476);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+64422, dvrr_stack+64212, dvrr_stack+6679, dvrr_stack+64062, dvrr_stack+10222, dvrr_stack+63902);

 /* compute (4 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+9757, dvrr_stack+713, dvrr_stack+6643, dvrr_stack+213, dvrr_stack+11425, NULL);

 /* compute (5 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+48, dvrr_stack+4476, dvrr_stack+4521, dvrr_stack+63872, dvrr_stack+131, dvrr_stack+9757);

 /* compute (6 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+64702, dvrr_stack+63902, dvrr_stack+7159, dvrr_stack+614, dvrr_stack+11503, dvrr_stack+48);

 /* compute (7 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+19709, dvrr_stack+64422, dvrr_stack+10672, dvrr_stack+64212, dvrr_stack+6679, dvrr_stack+64702);

 /* compute (1 0 | 0 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+64702, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+213, dvrr_stack+6634, dvrr_stack+64702, Data->F+6, Data->F+7, NULL);

 /* compute (3 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+713, dvrr_stack+6637, dvrr_stack+213, dvrr_stack+11422, dvrr_stack+6634, NULL);

 /* compute (4 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 0;
 vrr_build_xxxx(am,Data,dvrr_stack+64705, dvrr_stack+6643, dvrr_stack+713, dvrr_stack+11425, dvrr_stack+6637, NULL);

 /* compute (1 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+64720, dvrr_stack+161, dvrr_stack+710, NULL, NULL, Data->F+8);

 /* compute (2 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+9757, dvrr_stack+6653, dvrr_stack+64720, dvrr_stack+210, dvrr_stack+161, dvrr_stack+64702);

 /* compute (3 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+0, dvrr_stack+754, dvrr_stack+9757, dvrr_stack+3556, dvrr_stack+6653, dvrr_stack+213);

 /* compute (4 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+4476, dvrr_stack+10132, dvrr_stack+0, dvrr_stack+11431, dvrr_stack+754, dvrr_stack+713);

 /* compute (5 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+64729, dvrr_stack+4521, dvrr_stack+4476, dvrr_stack+131, dvrr_stack+10132, dvrr_stack+64705);

 /* compute (1 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+64702, dvrr_stack+569, dvrr_stack+1558, NULL, NULL, dvrr_stack+710);

 /* compute (2 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+131, dvrr_stack+30, dvrr_stack+64702, dvrr_stack+704, dvrr_stack+569, dvrr_stack+64720);

 /* compute (3 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+64792, dvrr_stack+4566, dvrr_stack+131, dvrr_stack+11449, dvrr_stack+30, dvrr_stack+9757);

 /* compute (4 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+30, dvrr_stack+2209, dvrr_stack+64792, dvrr_stack+11467, dvrr_stack+4566, dvrr_stack+0);

 /* compute (5 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+64852, dvrr_stack+1072, dvrr_stack+30, dvrr_stack+1876, dvrr_stack+2209, dvrr_stack+4476);

 /* compute (6 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+63794, dvrr_stack+7159, dvrr_stack+64852, dvrr_stack+11503, dvrr_stack+1072, dvrr_stack+64729);

 /* compute (1 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+0, dvrr_stack+1564, dvrr_stack+3565, NULL, NULL, dvrr_stack+1558);

 /* compute (2 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+2209, dvrr_stack+9775, dvrr_stack+0, dvrr_stack+121, dvrr_stack+1564, dvrr_stack+64702);

 /* compute (3 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+63962, dvrr_stack+2269, dvrr_stack+2209, dvrr_stack+11593, dvrr_stack+9775, dvrr_stack+131);

 /* compute (4 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+1558, dvrr_stack+7285, dvrr_stack+63962, dvrr_stack+11623, dvrr_stack+2269, dvrr_stack+64792);

 /* compute (5 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+20069, dvrr_stack+7385, dvrr_stack+1558, dvrr_stack+11683, dvrr_stack+7285, dvrr_stack+30);

 /* compute (6 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+20279, dvrr_stack+10462, dvrr_stack+20069, dvrr_stack+9892, dvrr_stack+7385, dvrr_stack+64852);

 /* compute (7 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+40543, dvrr_stack+10672, dvrr_stack+20279, dvrr_stack+6679, dvrr_stack+10462, dvrr_stack+63794);

 /* compute (1 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+63794, dvrr_stack+3575, dvrr_stack+6278, NULL, NULL, dvrr_stack+3565);

 /* compute (2 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+63839, dvrr_stack+10162, dvrr_stack+63794, dvrr_stack+723, dvrr_stack+3575, dvrr_stack+0);

 /* compute (3 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+10952, dvrr_stack+63839, dvrr_stack+219, dvrr_stack+10162, dvrr_stack+2209);

 /* compute (4 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+150, dvrr_stack+11042, dvrr_stack+0, dvrr_stack+10042, dvrr_stack+10952, dvrr_stack+63962);

 /* compute (5 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+375, dvrr_stack+11192, dvrr_stack+150, dvrr_stack+3618, dvrr_stack+11042, dvrr_stack+1558);

 /* compute (6 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+9757, dvrr_stack+46897, dvrr_stack+375, dvrr_stack+27917, dvrr_stack+11192, dvrr_stack+20069);

 /* compute (7 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10177, dvrr_stack+47212, dvrr_stack+9757, dvrr_stack+28142, dvrr_stack+46897, dvrr_stack+20279);

 /* compute (8 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+20069, dvrr_stack+47632, dvrr_stack+10177, dvrr_stack+28457, dvrr_stack+47212, dvrr_stack+40543);

 /* compute (1 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40543, dvrr_stack+6293, dvrr_stack+3923, NULL, NULL, dvrr_stack+6278);

 /* compute (2 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40606, dvrr_stack+14167, dvrr_stack+40543, dvrr_stack+1747, dvrr_stack+6293, dvrr_stack+63794);

 /* compute (3 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40732, dvrr_stack+14230, dvrr_stack+40606, dvrr_stack+883, dvrr_stack+14167, dvrr_stack+63839);

 /* compute (4 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40942, dvrr_stack+48172, dvrr_stack+40732, dvrr_stack+946, dvrr_stack+14230, dvrr_stack+0);

 /* compute (5 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+41257, dvrr_stack+48382, dvrr_stack+40942, dvrr_stack+28877, dvrr_stack+48172, dvrr_stack+150);

 /* compute (6 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+41698, dvrr_stack+48697, dvrr_stack+41257, dvrr_stack+29087, dvrr_stack+48382, dvrr_stack+375);

 /* compute (7 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+49138, dvrr_stack+41698, dvrr_stack+29402, dvrr_stack+48697, dvrr_stack+9757);

 /* compute (8 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+48172, dvrr_stack+49726, dvrr_stack+0, dvrr_stack+29843, dvrr_stack+49138, dvrr_stack+10177);

 /* compute (1 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+9757, dvrr_stack+3944, dvrr_stack+4412, NULL, NULL, dvrr_stack+3923);

 /* compute (2 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+9841, dvrr_stack+14857, dvrr_stack+9757, dvrr_stack+6314, dvrr_stack+3944, dvrr_stack+40543);

 /* compute (3 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10009, dvrr_stack+14941, dvrr_stack+9841, dvrr_stack+1957, dvrr_stack+14857, dvrr_stack+40606);

 /* compute (4 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10289, dvrr_stack+15109, dvrr_stack+10009, dvrr_stack+2041, dvrr_stack+14941, dvrr_stack+40732);

 /* compute (5 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10709, dvrr_stack+52102, dvrr_stack+10289, dvrr_stack+31691, dvrr_stack+15109, dvrr_stack+40942);

 /* compute (6 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+28877, dvrr_stack+52522, dvrr_stack+10709, dvrr_stack+31971, dvrr_stack+52102, dvrr_stack+41257);

 /* compute (7 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+40543, dvrr_stack+53110, dvrr_stack+28877, dvrr_stack+32391, dvrr_stack+52522, dvrr_stack+41698);

 /* compute (8 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+31691, dvrr_stack+53894, dvrr_stack+40543, dvrr_stack+32979, dvrr_stack+53110, dvrr_stack+0);

 /* compute (1 0 | 7 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+4440, dvrr_stack+19673, NULL, NULL, dvrr_stack+4412);

 /* compute (2 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+108, dvrr_stack+16894, dvrr_stack+0, dvrr_stack+3972, dvrr_stack+4440, dvrr_stack+9757);

 /* compute (3 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+324, dvrr_stack+17002, dvrr_stack+108, dvrr_stack+4088, dvrr_stack+16894, dvrr_stack+9841);

 /* compute (4 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+27917, dvrr_stack+17218, dvrr_stack+324, dvrr_stack+4196, dvrr_stack+17002, dvrr_stack+10009);

 /* compute (5 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+57170, dvrr_stack+27917, dvrr_stack+35527, dvrr_stack+17218, dvrr_stack+10289);

 /* compute (6 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+52102, dvrr_stack+57710, dvrr_stack+0, dvrr_stack+35887, dvrr_stack+57170, dvrr_stack+10709);

 /* compute (7 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+9757, dvrr_stack+58466, dvrr_stack+52102, dvrr_stack+36427, dvrr_stack+57710, dvrr_stack+28877);

 /* compute (8 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 8;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+52102, dvrr_stack+59474, dvrr_stack+9757, dvrr_stack+37183, dvrr_stack+58466, dvrr_stack+40543);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+40543, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 40543;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+40768, dvrr_stack+5320, NULL);
 tmp = dvrr_stack + 40768;
 target_ptr = Libderiv->deriv_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+41083, dvrr_stack+8497, NULL);
 tmp = dvrr_stack + 41083;
 target_ptr = Libderiv->deriv_classes[4][6][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,540,dvrr_stack+27917, dvrr_stack+12547, NULL);
 tmp = dvrr_stack + 27917;
 target_ptr = Libderiv->deriv_classes[4][7][11];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+41503, dvrr_stack+15949, NULL);
 tmp = dvrr_stack + 41503;
 target_ptr = Libderiv->deriv_classes[5][4][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,441,dvrr_stack+41818, dvrr_stack+18350, NULL);
 tmp = dvrr_stack + 41818;
 target_ptr = Libderiv->deriv_classes[5][5][11];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+9757, dvrr_stack+21545, NULL);
 tmp = dvrr_stack + 9757;
 target_ptr = Libderiv->deriv_classes[5][6][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,756,dvrr_stack+10345, dvrr_stack+25649, NULL);
 tmp = dvrr_stack + 10345;
 target_ptr = Libderiv->deriv_classes[5][7][11];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+11101, dvrr_stack+30431, NULL);
 tmp = dvrr_stack + 11101;
 target_ptr = Libderiv->deriv_classes[6][4][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+28877, dvrr_stack+33763, NULL);
 tmp = dvrr_stack + 28877;
 target_ptr = Libderiv->deriv_classes[6][5][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,784,dvrr_stack+57170, dvrr_stack+38191, NULL);
 tmp = dvrr_stack + 57170;
 target_ptr = Libderiv->deriv_classes[6][6][11];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1008,dvrr_stack+57954, dvrr_stack+43873, NULL);
 tmp = dvrr_stack + 57954;
 target_ptr = Libderiv->deriv_classes[6][7][11];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,540,dvrr_stack+0, dvrr_stack+50482, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[7][4][11];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,756,dvrr_stack+35527, dvrr_stack+54902, NULL);
 tmp = dvrr_stack + 35527;
 target_ptr = Libderiv->deriv_classes[7][5][11];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1008,dvrr_stack+3556, dvrr_stack+60770, NULL);
 tmp = dvrr_stack + 3556;
 target_ptr = Libderiv->deriv_classes[7][6][11];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1296,dvrr_stack+23309, dvrr_stack+66674, NULL);
 tmp = dvrr_stack + 23309;
 target_ptr = Libderiv->deriv_classes[7][7][11];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+42259, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 42259;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+11521, dvrr_stack+5320, NULL);
 tmp = dvrr_stack + 11521;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+58962, dvrr_stack+8497, NULL);
 tmp = dvrr_stack + 58962;
 target_ptr = Libderiv->deriv_classes[4][6][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,540,dvrr_stack+540, dvrr_stack+12547, NULL);
 tmp = dvrr_stack + 540;
 target_ptr = Libderiv->deriv_classes[4][7][10];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+29465, dvrr_stack+15949, NULL);
 tmp = dvrr_stack + 29465;
 target_ptr = Libderiv->deriv_classes[5][4][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+36283, dvrr_stack+18350, NULL);
 tmp = dvrr_stack + 36283;
 target_ptr = Libderiv->deriv_classes[5][5][10];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+16894, dvrr_stack+21545, NULL);
 tmp = dvrr_stack + 16894;
 target_ptr = Libderiv->deriv_classes[5][6][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,756,dvrr_stack+1558, dvrr_stack+25649, NULL);
 tmp = dvrr_stack + 1558;
 target_ptr = Libderiv->deriv_classes[5][7][10];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+36724, dvrr_stack+30431, NULL);
 tmp = dvrr_stack + 36724;
 target_ptr = Libderiv->deriv_classes[6][4][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+14857, dvrr_stack+33763, NULL);
 tmp = dvrr_stack + 14857;
 target_ptr = Libderiv->deriv_classes[6][5][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,784,dvrr_stack+6265, dvrr_stack+38191, NULL);
 tmp = dvrr_stack + 6265;
 target_ptr = Libderiv->deriv_classes[6][6][10];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1008,dvrr_stack+70562, dvrr_stack+43873, NULL);
 tmp = dvrr_stack + 70562;
 target_ptr = Libderiv->deriv_classes[6][7][10];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,540,dvrr_stack+49117, dvrr_stack+50482, NULL);
 tmp = dvrr_stack + 49117;
 target_ptr = Libderiv->deriv_classes[7][4][10];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,756,dvrr_stack+71570, dvrr_stack+54902, NULL);
 tmp = dvrr_stack + 71570;
 target_ptr = Libderiv->deriv_classes[7][5][10];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1008,dvrr_stack+72326, dvrr_stack+60770, NULL);
 tmp = dvrr_stack + 72326;
 target_ptr = Libderiv->deriv_classes[7][6][10];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1296,dvrr_stack+73334, dvrr_stack+66674, NULL);
 tmp = dvrr_stack + 73334;
 target_ptr = Libderiv->deriv_classes[7][7][10];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+17482, dvrr_stack+2881, NULL);
 tmp = dvrr_stack + 17482;
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

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+12547, dvrr_stack+15949, NULL);
 tmp = dvrr_stack + 12547;
 target_ptr = Libderiv->deriv_classes[5][4][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+15949, dvrr_stack+18350, NULL);
 tmp = dvrr_stack + 15949;
 target_ptr = Libderiv->deriv_classes[5][5][9];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+18350, dvrr_stack+21545, NULL);
 tmp = dvrr_stack + 18350;
 target_ptr = Libderiv->deriv_classes[5][6][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,756,dvrr_stack+21545, dvrr_stack+25649, NULL);
 tmp = dvrr_stack + 21545;
 target_ptr = Libderiv->deriv_classes[5][7][9];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+25649, dvrr_stack+30431, NULL);
 tmp = dvrr_stack + 25649;
 target_ptr = Libderiv->deriv_classes[6][4][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+30431, dvrr_stack+33763, NULL);
 tmp = dvrr_stack + 30431;
 target_ptr = Libderiv->deriv_classes[6][5][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,784,dvrr_stack+33763, dvrr_stack+38191, NULL);
 tmp = dvrr_stack + 33763;
 target_ptr = Libderiv->deriv_classes[6][6][9];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1008,dvrr_stack+22301, dvrr_stack+43873, NULL);
 tmp = dvrr_stack + 22301;
 target_ptr = Libderiv->deriv_classes[6][7][9];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,540,dvrr_stack+43873, dvrr_stack+50482, NULL);
 tmp = dvrr_stack + 43873;
 target_ptr = Libderiv->deriv_classes[7][4][9];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,756,dvrr_stack+50482, dvrr_stack+54902, NULL);
 tmp = dvrr_stack + 50482;
 target_ptr = Libderiv->deriv_classes[7][5][9];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1008,dvrr_stack+54902, dvrr_stack+60770, NULL);
 tmp = dvrr_stack + 54902;
 target_ptr = Libderiv->deriv_classes[7][6][9];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1296,dvrr_stack+60770, dvrr_stack+66674, NULL);
 tmp = dvrr_stack + 60770;
 target_ptr = Libderiv->deriv_classes[7][7][9];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+66674, dvrr_stack+2566, dvrr_stack+64062);
 tmp = dvrr_stack + 66674;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+66899, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 66899;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+67214, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 67214;
 target_ptr = Libderiv->deriv_classes[4][6][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,15,1,dvrr_stack+67634, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 67634;
 target_ptr = Libderiv->deriv_classes[4][7][8];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+68174, dvrr_stack+15508, dvrr_stack+64212);
 tmp = dvrr_stack + 68174;
 target_ptr = Libderiv->deriv_classes[5][4][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+68489, dvrr_stack+17762, dvrr_stack+14542);
 tmp = dvrr_stack + 68489;
 target_ptr = Libderiv->deriv_classes[5][5][8];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+68930, dvrr_stack+20789, dvrr_stack+15508);
 tmp = dvrr_stack + 68930;
 target_ptr = Libderiv->deriv_classes[5][6][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,21,1,dvrr_stack+69518, dvrr_stack+24704, dvrr_stack+17762);
 tmp = dvrr_stack + 69518;
 target_ptr = Libderiv->deriv_classes[5][7][8];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+62066, dvrr_stack+29843, dvrr_stack+64422);
 tmp = dvrr_stack + 62066;
 target_ptr = Libderiv->deriv_classes[6][4][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+62486, dvrr_stack+32979, dvrr_stack+28457);
 tmp = dvrr_stack + 62486;
 target_ptr = Libderiv->deriv_classes[6][5][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+63074, dvrr_stack+37183, dvrr_stack+29843);
 tmp = dvrr_stack + 63074;
 target_ptr = Libderiv->deriv_classes[6][6][8];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,28,1,dvrr_stack+55910, dvrr_stack+42613, dvrr_stack+32979);
 tmp = dvrr_stack + 55910;
 target_ptr = Libderiv->deriv_classes[6][7][8];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,36,1,dvrr_stack+51238, dvrr_stack+49726, dvrr_stack+19709);
 tmp = dvrr_stack + 51238;
 target_ptr = Libderiv->deriv_classes[7][4][8];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,36,1,dvrr_stack+44413, dvrr_stack+53894, dvrr_stack+47632);
 tmp = dvrr_stack + 44413;
 target_ptr = Libderiv->deriv_classes[7][5][8];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,36,1,dvrr_stack+45169, dvrr_stack+59474, dvrr_stack+49726);
 tmp = dvrr_stack + 45169;
 target_ptr = Libderiv->deriv_classes[7][6][8];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,36,1,dvrr_stack+46177, dvrr_stack+65054, dvrr_stack+53894);
 tmp = dvrr_stack + 46177;
 target_ptr = Libderiv->deriv_classes[7][7][8];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+70274, dvrr_stack+2566, dvrr_stack+64062);
 tmp = dvrr_stack + 70274;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+51778, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 51778;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+38191, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 38191;
 target_ptr = Libderiv->deriv_classes[4][6][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,15,1,dvrr_stack+38611, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 38611;
 target_ptr = Libderiv->deriv_classes[4][7][7];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+39151, dvrr_stack+15508, dvrr_stack+64212);
 tmp = dvrr_stack + 39151;
 target_ptr = Libderiv->deriv_classes[5][4][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+39466, dvrr_stack+17762, dvrr_stack+14542);
 tmp = dvrr_stack + 39466;
 target_ptr = Libderiv->deriv_classes[5][5][7];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+39907, dvrr_stack+20789, dvrr_stack+15508);
 tmp = dvrr_stack + 39907;
 target_ptr = Libderiv->deriv_classes[5][6][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,21,1,dvrr_stack+34547, dvrr_stack+24704, dvrr_stack+17762);
 tmp = dvrr_stack + 34547;
 target_ptr = Libderiv->deriv_classes[5][7][7];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+31019, dvrr_stack+29843, dvrr_stack+64422);
 tmp = dvrr_stack + 31019;
 target_ptr = Libderiv->deriv_classes[6][4][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+26069, dvrr_stack+32979, dvrr_stack+28457);
 tmp = dvrr_stack + 26069;
 target_ptr = Libderiv->deriv_classes[6][5][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+26657, dvrr_stack+37183, dvrr_stack+29843);
 tmp = dvrr_stack + 26657;
 target_ptr = Libderiv->deriv_classes[6][6][7];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,28,1,dvrr_stack+12862, dvrr_stack+42613, dvrr_stack+32979);
 tmp = dvrr_stack + 12862;
 target_ptr = Libderiv->deriv_classes[6][7][7];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,36,1,dvrr_stack+18938, dvrr_stack+49726, dvrr_stack+19709);
 tmp = dvrr_stack + 18938;
 target_ptr = Libderiv->deriv_classes[7][4][7];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,36,1,dvrr_stack+74630, dvrr_stack+53894, dvrr_stack+47632);
 tmp = dvrr_stack + 74630;
 target_ptr = Libderiv->deriv_classes[7][5][7];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,36,1,dvrr_stack+75386, dvrr_stack+59474, dvrr_stack+49726);
 tmp = dvrr_stack + 75386;
 target_ptr = Libderiv->deriv_classes[7][6][7];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,36,1,dvrr_stack+76394, dvrr_stack+65054, dvrr_stack+53894);
 tmp = dvrr_stack + 76394;
 target_ptr = Libderiv->deriv_classes[7][7][7];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+56918, dvrr_stack+2566, dvrr_stack+64062);
 tmp = dvrr_stack + 56918;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+63858, dvrr_stack+4900, dvrr_stack+1333);
 tmp = dvrr_stack + 63858;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+27441, dvrr_stack+7957, dvrr_stack+2566);
 tmp = dvrr_stack + 27441;
 target_ptr = Libderiv->deriv_classes[4][6][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,15,1,dvrr_stack+13870, dvrr_stack+11872, dvrr_stack+4900);
 tmp = dvrr_stack + 13870;
 target_ptr = Libderiv->deriv_classes[4][7][6];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+16390, dvrr_stack+15508, dvrr_stack+64212);
 tmp = dvrr_stack + 16390;
 target_ptr = Libderiv->deriv_classes[5][4][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+9037, dvrr_stack+17762, dvrr_stack+14542);
 tmp = dvrr_stack + 9037;
 target_ptr = Libderiv->deriv_classes[5][5][6];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+11836, dvrr_stack+20789, dvrr_stack+15508);
 tmp = dvrr_stack + 11836;
 target_ptr = Libderiv->deriv_classes[5][6][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,21,1,dvrr_stack+77690, dvrr_stack+24704, dvrr_stack+17762);
 tmp = dvrr_stack + 77690;
 target_ptr = Libderiv->deriv_classes[5][7][6];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+5740, dvrr_stack+29843, dvrr_stack+64422);
 tmp = dvrr_stack + 5740;
 target_ptr = Libderiv->deriv_classes[6][4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+64173, dvrr_stack+32979, dvrr_stack+28457);
 tmp = dvrr_stack + 64173;
 target_ptr = Libderiv->deriv_classes[6][5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+24605, dvrr_stack+37183, dvrr_stack+29843);
 tmp = dvrr_stack + 24605;
 target_ptr = Libderiv->deriv_classes[6][6][6];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,28,1,dvrr_stack+78446, dvrr_stack+42613, dvrr_stack+32979);
 tmp = dvrr_stack + 78446;
 target_ptr = Libderiv->deriv_classes[6][7][6];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,36,1,dvrr_stack+42484, dvrr_stack+49726, dvrr_stack+19709);
 tmp = dvrr_stack + 42484;
 target_ptr = Libderiv->deriv_classes[7][4][6];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,36,1,dvrr_stack+43024, dvrr_stack+53894, dvrr_stack+47632);
 tmp = dvrr_stack + 43024;
 target_ptr = Libderiv->deriv_classes[7][5][6];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,36,1,dvrr_stack+79454, dvrr_stack+59474, dvrr_stack+49726);
 tmp = dvrr_stack + 79454;
 target_ptr = Libderiv->deriv_classes[7][6][6];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,36,1,dvrr_stack+80462, dvrr_stack+65054, dvrr_stack+53894);
 tmp = dvrr_stack + 80462;
 target_ptr = Libderiv->deriv_classes[7][7][6];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+64761, dvrr_stack+14542, dvrr_stack+1183);
 tmp = dvrr_stack + 64761;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+64986, dvrr_stack+15508, dvrr_stack+2356);
 tmp = dvrr_stack + 64986;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+65301, dvrr_stack+17762, dvrr_stack+4620);
 tmp = dvrr_stack + 65301;
 target_ptr = Libderiv->deriv_classes[4][6][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,36,dvrr_stack+65721, dvrr_stack+20789, dvrr_stack+7597);
 tmp = dvrr_stack + 65721;
 target_ptr = Libderiv->deriv_classes[4][7][2];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+66261, dvrr_stack+28457, dvrr_stack+1333);
 tmp = dvrr_stack + 66261;
 target_ptr = Libderiv->deriv_classes[5][4][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+19478, dvrr_stack+29843, dvrr_stack+2566);
 tmp = dvrr_stack + 19478;
 target_ptr = Libderiv->deriv_classes[5][5][2];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+81758, dvrr_stack+32979, dvrr_stack+4900);
 tmp = dvrr_stack + 81758;
 target_ptr = Libderiv->deriv_classes[5][6][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,36,dvrr_stack+82346, dvrr_stack+37183, dvrr_stack+7957);
 tmp = dvrr_stack + 82346;
 target_ptr = Libderiv->deriv_classes[5][7][2];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+83102, dvrr_stack+47632, dvrr_stack+14542);
 tmp = dvrr_stack + 83102;
 target_ptr = Libderiv->deriv_classes[6][4][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+83522, dvrr_stack+49726, dvrr_stack+15508);
 tmp = dvrr_stack + 83522;
 target_ptr = Libderiv->deriv_classes[6][5][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+84110, dvrr_stack+53894, dvrr_stack+17762);
 tmp = dvrr_stack + 84110;
 target_ptr = Libderiv->deriv_classes[6][6][2];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,36,dvrr_stack+84894, dvrr_stack+59474, dvrr_stack+20789);
 tmp = dvrr_stack + 84894;
 target_ptr = Libderiv->deriv_classes[6][7][2];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,15,dvrr_stack+7049, dvrr_stack+20069, dvrr_stack+28457);
 tmp = dvrr_stack + 7049;
 target_ptr = Libderiv->deriv_classes[7][4][2];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,21,dvrr_stack+85902, dvrr_stack+48172, dvrr_stack+29843);
 tmp = dvrr_stack + 85902;
 target_ptr = Libderiv->deriv_classes[7][5][2];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,28,dvrr_stack+86658, dvrr_stack+31691, dvrr_stack+32979);
 tmp = dvrr_stack + 86658;
 target_ptr = Libderiv->deriv_classes[7][6][2];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_k(Data,36,dvrr_stack+87666, dvrr_stack+52102, dvrr_stack+37183);
 tmp = dvrr_stack + 87666;
 target_ptr = Libderiv->deriv_classes[7][7][2];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+31439, dvrr_stack+14542, dvrr_stack+1183);
 tmp = dvrr_stack + 31439;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+3196, dvrr_stack+15508, dvrr_stack+2356);
 tmp = dvrr_stack + 3196;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+88962, dvrr_stack+17762, dvrr_stack+4620);
 tmp = dvrr_stack + 88962;
 target_ptr = Libderiv->deriv_classes[4][6][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,36,dvrr_stack+89382, dvrr_stack+20789, dvrr_stack+7597);
 tmp = dvrr_stack + 89382;
 target_ptr = Libderiv->deriv_classes[4][7][1];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+89922, dvrr_stack+28457, dvrr_stack+1333);
 tmp = dvrr_stack + 89922;
 target_ptr = Libderiv->deriv_classes[5][4][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+90237, dvrr_stack+29843, dvrr_stack+2566);
 tmp = dvrr_stack + 90237;
 target_ptr = Libderiv->deriv_classes[5][5][1];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+90678, dvrr_stack+32979, dvrr_stack+4900);
 tmp = dvrr_stack + 90678;
 target_ptr = Libderiv->deriv_classes[5][6][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,36,dvrr_stack+91266, dvrr_stack+37183, dvrr_stack+7957);
 tmp = dvrr_stack + 91266;
 target_ptr = Libderiv->deriv_classes[5][7][1];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+92022, dvrr_stack+47632, dvrr_stack+14542);
 tmp = dvrr_stack + 92022;
 target_ptr = Libderiv->deriv_classes[6][4][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+92442, dvrr_stack+49726, dvrr_stack+15508);
 tmp = dvrr_stack + 92442;
 target_ptr = Libderiv->deriv_classes[6][5][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+93030, dvrr_stack+53894, dvrr_stack+17762);
 tmp = dvrr_stack + 93030;
 target_ptr = Libderiv->deriv_classes[6][6][1];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,36,dvrr_stack+93814, dvrr_stack+59474, dvrr_stack+20789);
 tmp = dvrr_stack + 93814;
 target_ptr = Libderiv->deriv_classes[6][7][1];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,15,dvrr_stack+94822, dvrr_stack+20069, dvrr_stack+28457);
 tmp = dvrr_stack + 94822;
 target_ptr = Libderiv->deriv_classes[7][4][1];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,21,dvrr_stack+95362, dvrr_stack+48172, dvrr_stack+29843);
 tmp = dvrr_stack + 95362;
 target_ptr = Libderiv->deriv_classes[7][5][1];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,28,dvrr_stack+96118, dvrr_stack+31691, dvrr_stack+32979);
 tmp = dvrr_stack + 96118;
 target_ptr = Libderiv->deriv_classes[7][6][1];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_k(Data,36,dvrr_stack+97126, dvrr_stack+52102, dvrr_stack+37183);
 tmp = dvrr_stack + 97126;
 target_ptr = Libderiv->deriv_classes[7][7][1];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+9478, dvrr_stack+14542, dvrr_stack+1183);
 tmp = dvrr_stack + 9478;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+98422, dvrr_stack+15508, dvrr_stack+2356);
 tmp = dvrr_stack + 98422;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+98737, dvrr_stack+17762, dvrr_stack+4620);
 tmp = dvrr_stack + 98737;
 target_ptr = Libderiv->deriv_classes[4][6][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,36,dvrr_stack+99157, dvrr_stack+20789, dvrr_stack+7597);
 tmp = dvrr_stack + 99157;
 target_ptr = Libderiv->deriv_classes[4][7][0];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+4564, dvrr_stack+28457, dvrr_stack+1333);
 tmp = dvrr_stack + 4564;
 target_ptr = Libderiv->deriv_classes[5][4][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+1080, dvrr_stack+29843, dvrr_stack+2566);
 tmp = dvrr_stack + 1080;
 target_ptr = Libderiv->deriv_classes[5][5][0];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+99697, dvrr_stack+32979, dvrr_stack+4900);
 tmp = dvrr_stack + 99697;
 target_ptr = Libderiv->deriv_classes[5][6][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,36,dvrr_stack+100285, dvrr_stack+37183, dvrr_stack+7957);
 tmp = dvrr_stack + 100285;
 target_ptr = Libderiv->deriv_classes[5][7][0];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+2314, dvrr_stack+47632, dvrr_stack+14542);
 tmp = dvrr_stack + 2314;
 target_ptr = Libderiv->deriv_classes[6][4][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+47473, dvrr_stack+49726, dvrr_stack+15508);
 tmp = dvrr_stack + 47473;
 target_ptr = Libderiv->deriv_classes[6][5][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+49657, dvrr_stack+53894, dvrr_stack+17762);
 tmp = dvrr_stack + 49657;
 target_ptr = Libderiv->deriv_classes[6][6][0];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,36,dvrr_stack+53722, dvrr_stack+59474, dvrr_stack+20789);
 tmp = dvrr_stack + 53722;
 target_ptr = Libderiv->deriv_classes[6][7][0];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,15,dvrr_stack+59382, dvrr_stack+20069, dvrr_stack+28457);
 tmp = dvrr_stack + 59382;
 target_ptr = Libderiv->deriv_classes[7][4][0];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,21,dvrr_stack+19919, dvrr_stack+48172, dvrr_stack+29843);
 tmp = dvrr_stack + 19919;
 target_ptr = Libderiv->deriv_classes[7][5][0];
 for(i=0;i<756;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,28,dvrr_stack+48061, dvrr_stack+31691, dvrr_stack+32979);
 tmp = dvrr_stack + 48061;
 target_ptr = Libderiv->deriv_classes[7][6][0];
 for(i=0;i<1008;i++)
   target_ptr[i] += tmp[i];

 /* compute (7 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_k(Data,36,dvrr_stack+31664, dvrr_stack+52102, dvrr_stack+37183);
 tmp = dvrr_stack + 31664;
 target_ptr = Libderiv->deriv_classes[7][7][0];
 for(i=0;i<1296;i++)
   target_ptr[i] += tmp[i];


}

