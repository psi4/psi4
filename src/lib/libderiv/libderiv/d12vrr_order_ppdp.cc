#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (pp|dp) integrals */

void d12vrr_order_ppdp(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+0, dvrr_stack+3, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+12, dvrr_stack+0, Data->F+0, Data->F+1, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+6, NULL, NULL, dvrr_stack+0);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->dvrr_classes[1][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+39, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+42, dvrr_stack+3, dvrr_stack+39, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+48, dvrr_stack+6, dvrr_stack+42, dvrr_stack+0, dvrr_stack+3, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+58, dvrr_stack+15, dvrr_stack+6, dvrr_stack+12, dvrr_stack+0, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+68, dvrr_stack+58, dvrr_stack+48, NULL, NULL, dvrr_stack+6);
 tmp = dvrr_stack + 68;
 target_ptr = Libderiv->dvrr_classes[1][3];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+98,dvrr_stack+68,dvrr_stack+21,3);


 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+152, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+155, dvrr_stack+39, dvrr_stack+152, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+161, dvrr_stack+42, dvrr_stack+155, dvrr_stack+3, dvrr_stack+39, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+171, dvrr_stack+48, dvrr_stack+161, dvrr_stack+6, dvrr_stack+42, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+186, dvrr_stack+58, dvrr_stack+48, dvrr_stack+15, dvrr_stack+6, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+201, dvrr_stack+186, dvrr_stack+171, NULL, NULL, dvrr_stack+48);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+246,dvrr_stack+201,dvrr_stack+68,3);


 /* compute (1 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+336,dvrr_stack+246,dvrr_stack+98,3);


 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,18,dvrr_stack+444, dvrr_stack+336, dvrr_stack+21);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+498, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+501, dvrr_stack+152, dvrr_stack+498, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+507, dvrr_stack+155, dvrr_stack+501, dvrr_stack+39, dvrr_stack+152, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+517, dvrr_stack+161, dvrr_stack+507, dvrr_stack+42, dvrr_stack+155, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+532, dvrr_stack+171, dvrr_stack+517, dvrr_stack+48, dvrr_stack+161, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+553, dvrr_stack+186, dvrr_stack+171, dvrr_stack+58, dvrr_stack+48, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+574, dvrr_stack+553, dvrr_stack+532, NULL, NULL, dvrr_stack+171);

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+637,dvrr_stack+574,dvrr_stack+201,3);


 /* compute (1 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+772,dvrr_stack+637,dvrr_stack+246,3);


 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,30,dvrr_stack+952, dvrr_stack+772, dvrr_stack+68);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1042, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+2);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1051, dvrr_stack+6, dvrr_stack+42, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1069, dvrr_stack+21, dvrr_stack+1051, dvrr_stack+15, dvrr_stack+6, dvrr_stack+1042);
 tmp = dvrr_stack + 1069;
 target_ptr = Libderiv->dvrr_classes[2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1105, dvrr_stack+48, dvrr_stack+161, NULL, NULL, dvrr_stack+42);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1135, dvrr_stack+68, dvrr_stack+1105, dvrr_stack+58, dvrr_stack+48, dvrr_stack+1051);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+1195,dvrr_stack+1135,dvrr_stack+1069,6);


 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1303, dvrr_stack+171, dvrr_stack+517, NULL, NULL, dvrr_stack+161);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1348, dvrr_stack+201, dvrr_stack+1303, dvrr_stack+186, dvrr_stack+171, dvrr_stack+1105);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+1438,dvrr_stack+1348,dvrr_stack+1135,6);


 /* compute (2 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+1618,dvrr_stack+1438,dvrr_stack+1195,6);


 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,36,dvrr_stack+1834, dvrr_stack+1618, dvrr_stack+1069);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1942, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1945, dvrr_stack+498, dvrr_stack+1942, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1951, dvrr_stack+501, dvrr_stack+1945, dvrr_stack+152, dvrr_stack+498, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1961, dvrr_stack+507, dvrr_stack+1951, dvrr_stack+155, dvrr_stack+501, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1976, dvrr_stack+517, dvrr_stack+1961, dvrr_stack+161, dvrr_stack+507, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1997, dvrr_stack+532, dvrr_stack+1976, NULL, NULL, dvrr_stack+517);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2060, dvrr_stack+574, dvrr_stack+1997, dvrr_stack+553, dvrr_stack+532, dvrr_stack+1303);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+2186,dvrr_stack+2060,dvrr_stack+1348,6);


 /* compute (2 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+2456,dvrr_stack+2186,dvrr_stack+1438,6);


 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,60,dvrr_stack+2816, dvrr_stack+2456, dvrr_stack+1135);

 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,18,dvrr_stack+1976, dvrr_stack+336, dvrr_stack+21);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,30,dvrr_stack+2996, dvrr_stack+772, dvrr_stack+68);

 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,36,dvrr_stack+3086, dvrr_stack+1618, dvrr_stack+1069);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,60,dvrr_stack+3194, dvrr_stack+2456, dvrr_stack+1135);

 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,18,dvrr_stack+3374, dvrr_stack+336, dvrr_stack+21);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,30,dvrr_stack+336, dvrr_stack+772, dvrr_stack+68);

 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,36,dvrr_stack+772, dvrr_stack+1618, dvrr_stack+1069);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,60,dvrr_stack+1618, dvrr_stack+2456, dvrr_stack+1135);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1942, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+1);

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+2456,dvrr_stack+21,dvrr_stack+1942,3);


 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+2483, dvrr_stack+2456, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+2030, dvrr_stack+246, NULL);
 tmp = dvrr_stack + 2030;
 target_ptr = Libderiv->deriv_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+426, dvrr_stack+98, NULL);
 tmp = dvrr_stack + 426;
 target_ptr = Libderiv->deriv_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+2492, dvrr_stack+637, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+498, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+2537, dvrr_stack+1942, dvrr_stack+1042, dvrr_stack+12, dvrr_stack+0, dvrr_stack+498);

 /* compute (2 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+2555,dvrr_stack+1069,dvrr_stack+2537,6);


 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+2609, dvrr_stack+2555, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+2627, dvrr_stack+1438, NULL);
 tmp = dvrr_stack + 2627;
 target_ptr = Libderiv->deriv_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+1798, dvrr_stack+1195, NULL);
 tmp = dvrr_stack + 1798;
 target_ptr = Libderiv->deriv_classes[2][2][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+2687, dvrr_stack+2186, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+2777, dvrr_stack+2456, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+2786, dvrr_stack+246, NULL);
 tmp = dvrr_stack + 2786;
 target_ptr = Libderiv->deriv_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+880, dvrr_stack+98, NULL);
 tmp = dvrr_stack + 880;
 target_ptr = Libderiv->deriv_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+898, dvrr_stack+637, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+532, dvrr_stack+2555, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+3428, dvrr_stack+1438, NULL);
 tmp = dvrr_stack + 3428;
 target_ptr = Libderiv->deriv_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+3488, dvrr_stack+1195, NULL);
 tmp = dvrr_stack + 3488;
 target_ptr = Libderiv->deriv_classes[2][2][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+3524, dvrr_stack+2186, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+943, dvrr_stack+2456, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+3614, dvrr_stack+246, NULL);
 tmp = dvrr_stack + 3614;
 target_ptr = Libderiv->deriv_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+246, dvrr_stack+98, NULL);
 tmp = dvrr_stack + 246;
 target_ptr = Libderiv->deriv_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+98, dvrr_stack+637, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+637, dvrr_stack+2555, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+655, dvrr_stack+1438, NULL);
 tmp = dvrr_stack + 655;
 target_ptr = Libderiv->deriv_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+1438, dvrr_stack+1195, NULL);
 tmp = dvrr_stack + 1438;
 target_ptr = Libderiv->deriv_classes[2][2][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+1195, dvrr_stack+2186, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+2186, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+143, dvrr_stack+21, dvrr_stack+2186);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2189, dvrr_stack+201, dvrr_stack+21);
 tmp = dvrr_stack + 2189;
 target_ptr = Libderiv->deriv_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+1285, dvrr_stack+68, dvrr_stack+1942);
 tmp = dvrr_stack + 1285;
 target_ptr = Libderiv->deriv_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+2219, dvrr_stack+574, dvrr_stack+68);

 /* compute (2 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+2264, dvrr_stack+2186, dvrr_stack+498, Data->F+0, Data->F+1, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+2270, dvrr_stack+1069, dvrr_stack+2264);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+2288, dvrr_stack+1348, dvrr_stack+1069);
 tmp = dvrr_stack + 2288;
 target_ptr = Libderiv->deriv_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+2348, dvrr_stack+1135, dvrr_stack+2537);
 tmp = dvrr_stack + 2348;
 target_ptr = Libderiv->deriv_classes[2][2][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+2384, dvrr_stack+2060, dvrr_stack+1135);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+2474, dvrr_stack+21, dvrr_stack+2186);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1474, dvrr_stack+201, dvrr_stack+21);
 tmp = dvrr_stack + 1474;
 target_ptr = Libderiv->deriv_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+1504, dvrr_stack+68, dvrr_stack+1942);
 tmp = dvrr_stack + 1504;
 target_ptr = Libderiv->deriv_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+1522, dvrr_stack+574, dvrr_stack+68);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+1567, dvrr_stack+1069, dvrr_stack+2264);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+264, dvrr_stack+1348, dvrr_stack+1069);
 tmp = dvrr_stack + 264;
 target_ptr = Libderiv->deriv_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+2555, dvrr_stack+1135, dvrr_stack+2537);
 tmp = dvrr_stack + 2555;
 target_ptr = Libderiv->deriv_classes[2][2][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+3644, dvrr_stack+2060, dvrr_stack+1135);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+1585, dvrr_stack+21, dvrr_stack+2186);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+715, dvrr_stack+201, dvrr_stack+21);
 tmp = dvrr_stack + 715;
 target_ptr = Libderiv->deriv_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+2591, dvrr_stack+68, dvrr_stack+1942);
 tmp = dvrr_stack + 2591;
 target_ptr = Libderiv->deriv_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+3734, dvrr_stack+574, dvrr_stack+68);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+1594, dvrr_stack+1069, dvrr_stack+2264);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+550, dvrr_stack+1348, dvrr_stack+1069);
 tmp = dvrr_stack + 550;
 target_ptr = Libderiv->deriv_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+3779, dvrr_stack+1135, dvrr_stack+2537);
 tmp = dvrr_stack + 3779;
 target_ptr = Libderiv->deriv_classes[2][2][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+3815, dvrr_stack+2060, dvrr_stack+1135);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+2060,dvrr_stack+58,dvrr_stack+15,1);


 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+2264, dvrr_stack+2060, NULL);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+2078,dvrr_stack+186,dvrr_stack+58,1);


 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+2108, dvrr_stack+2078, NULL);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+2118, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+2121, dvrr_stack+3, dvrr_stack+39, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+2130, dvrr_stack+1042, dvrr_stack+2121, dvrr_stack+0, dvrr_stack+3, dvrr_stack+2118);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2148, dvrr_stack+42, dvrr_stack+155, NULL, NULL, dvrr_stack+39);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3905, dvrr_stack+1051, dvrr_stack+2148, dvrr_stack+6, dvrr_stack+42, dvrr_stack+2121);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+3941, dvrr_stack+1069, dvrr_stack+3905, dvrr_stack+21, dvrr_stack+1051, dvrr_stack+2130);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+4001, dvrr_stack+161, dvrr_stack+507, NULL, NULL, dvrr_stack+155);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4031, dvrr_stack+1105, dvrr_stack+4001, dvrr_stack+48, dvrr_stack+161, dvrr_stack+2148);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+4091, dvrr_stack+1135, dvrr_stack+4031, dvrr_stack+68, dvrr_stack+1105, dvrr_stack+3905);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+4191,dvrr_stack+4091,dvrr_stack+3941,10);


 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+4371, dvrr_stack+4191, NULL);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4431, dvrr_stack+517, dvrr_stack+1961, NULL, NULL, dvrr_stack+507);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4476, dvrr_stack+1303, dvrr_stack+4431, dvrr_stack+171, dvrr_stack+517, dvrr_stack+4001);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4566, dvrr_stack+1348, dvrr_stack+4476, dvrr_stack+201, dvrr_stack+1303, dvrr_stack+4031);

 /* compute (3 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+4716,dvrr_stack+4566,dvrr_stack+4091,10);


 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+1303, dvrr_stack+4716, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+6, dvrr_stack+2060, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+48, dvrr_stack+2078, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+4431, dvrr_stack+4191, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+5016, dvrr_stack+4716, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+1612, dvrr_stack+2060, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+2060, dvrr_stack+2078, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+4491, dvrr_stack+4191, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+4191, dvrr_stack+4716, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+4716, dvrr_stack+58, dvrr_stack+12);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+4722, dvrr_stack+186, dvrr_stack+15);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+4732, dvrr_stack+498, dvrr_stack+2118, Data->F+1, Data->F+2, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+4738, dvrr_stack+2537, dvrr_stack+2130, dvrr_stack+1942, dvrr_stack+1042, dvrr_stack+4732);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+4768, dvrr_stack+4091, dvrr_stack+4738);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+4828, dvrr_stack+4566, dvrr_stack+3941);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+4732, dvrr_stack+58, dvrr_stack+12);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+2537, dvrr_stack+186, dvrr_stack+15);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+4928, dvrr_stack+4091, dvrr_stack+4738);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+5116, dvrr_stack+4566, dvrr_stack+3941);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+1942, dvrr_stack+58, dvrr_stack+12);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+4988, dvrr_stack+186, dvrr_stack+15);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+4291, dvrr_stack+4091, dvrr_stack+4738);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+5216, dvrr_stack+4566, dvrr_stack+3941);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+4738, dvrr_stack+21, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+2070, dvrr_stack+3941, dvrr_stack+21);
 tmp = dvrr_stack + 2070;
 target_ptr = Libderiv->deriv_classes[2][2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+4744, dvrr_stack+68, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+171, dvrr_stack+4091, dvrr_stack+68);
 tmp = dvrr_stack + 171;
 target_ptr = Libderiv->deriv_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+4998, dvrr_stack+1069, dvrr_stack+15);
 tmp = dvrr_stack + 4998;
 target_ptr = Libderiv->deriv_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+12, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+4754, dvrr_stack+2118, dvrr_stack+12, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+4351, dvrr_stack+39, dvrr_stack+152, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+1403, dvrr_stack+2121, dvrr_stack+4351, dvrr_stack+3, dvrr_stack+39, dvrr_stack+12);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+4551, dvrr_stack+2130, dvrr_stack+1403, dvrr_stack+1042, dvrr_stack+2121, dvrr_stack+4754);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2118, dvrr_stack+155, dvrr_stack+501, NULL, NULL, dvrr_stack+152);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+4581, dvrr_stack+2148, dvrr_stack+2118, dvrr_stack+42, dvrr_stack+155, dvrr_stack+4351);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+4617, dvrr_stack+3905, dvrr_stack+4581, dvrr_stack+1051, dvrr_stack+2148, dvrr_stack+1403);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+5316, dvrr_stack+3941, dvrr_stack+4617, dvrr_stack+1069, dvrr_stack+3905, dvrr_stack+4551);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+5406, dvrr_stack+5316, dvrr_stack+1069);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+4551, dvrr_stack+1135, dvrr_stack+58);
 tmp = dvrr_stack + 4551;
 target_ptr = Libderiv->deriv_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+3905, dvrr_stack+507, dvrr_stack+1951, NULL, NULL, dvrr_stack+501);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+5466, dvrr_stack+4001, dvrr_stack+3905, dvrr_stack+161, dvrr_stack+507, dvrr_stack+2118);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+5526, dvrr_stack+4031, dvrr_stack+5466, dvrr_stack+1105, dvrr_stack+4001, dvrr_stack+4581);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+5626, dvrr_stack+4091, dvrr_stack+5526, dvrr_stack+1135, dvrr_stack+4031, dvrr_stack+4617);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+4581, dvrr_stack+5626, dvrr_stack+1135);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+0, dvrr_stack+21, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3905, dvrr_stack+3941, dvrr_stack+21);
 tmp = dvrr_stack + 3905;
 target_ptr = Libderiv->deriv_classes[2][2][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4681, dvrr_stack+68, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+4001, dvrr_stack+4091, dvrr_stack+68);
 tmp = dvrr_stack + 4001;
 target_ptr = Libderiv->deriv_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+4691, dvrr_stack+1069, dvrr_stack+15);
 tmp = dvrr_stack + 4691;
 target_ptr = Libderiv->deriv_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+5466, dvrr_stack+5316, dvrr_stack+1069);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+4061, dvrr_stack+1135, dvrr_stack+58);
 tmp = dvrr_stack + 4061;
 target_ptr = Libderiv->deriv_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+5526, dvrr_stack+5626, dvrr_stack+1135);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+4709, dvrr_stack+21, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2118, dvrr_stack+3941, dvrr_stack+21);
 tmp = dvrr_stack + 2118;
 target_ptr = Libderiv->deriv_classes[2][2][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+21, dvrr_stack+68, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+3941, dvrr_stack+4091, dvrr_stack+68);
 tmp = dvrr_stack + 3941;
 target_ptr = Libderiv->deriv_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+68, dvrr_stack+1069, dvrr_stack+15);
 tmp = dvrr_stack + 68;
 target_ptr = Libderiv->deriv_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+4091, dvrr_stack+5316, dvrr_stack+1069);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5316, dvrr_stack+1135, dvrr_stack+58);
 tmp = dvrr_stack + 5316;
 target_ptr = Libderiv->deriv_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+5776, dvrr_stack+5626, dvrr_stack+1135);

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,18,dvrr_stack+5626, dvrr_stack+444, NULL);
 tmp = dvrr_stack + 5626;
 target_ptr = Libderiv->deriv2_classes[1][2][143];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,30,dvrr_stack+5644, dvrr_stack+952, NULL);
 tmp = dvrr_stack + 5644;
 target_ptr = Libderiv->deriv2_classes[1][3][143];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,36,dvrr_stack+5674, dvrr_stack+1834, NULL);
 tmp = dvrr_stack + 5674;
 target_ptr = Libderiv->deriv2_classes[2][2][143];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,60,dvrr_stack+5346, dvrr_stack+2816, NULL);
 tmp = dvrr_stack + 5346;
 target_ptr = Libderiv->deriv2_classes[2][3][143];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,18,dvrr_stack+5710, dvrr_stack+444, NULL);
 tmp = dvrr_stack + 5710;
 target_ptr = Libderiv->deriv2_classes[1][2][131];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,30,dvrr_stack+5728, dvrr_stack+952, NULL);
 tmp = dvrr_stack + 5728;
 target_ptr = Libderiv->deriv2_classes[1][3][131];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,36,dvrr_stack+4151, dvrr_stack+1834, NULL);
 tmp = dvrr_stack + 4151;
 target_ptr = Libderiv->deriv2_classes[2][2][131];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,60,dvrr_stack+1042, dvrr_stack+2816, NULL);
 tmp = dvrr_stack + 1042;
 target_ptr = Libderiv->deriv2_classes[2][3][131];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+5758, dvrr_stack+1976, NULL);
 tmp = dvrr_stack + 5758;
 target_ptr = Libderiv->deriv2_classes[1][2][130];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+2154, dvrr_stack+2996, NULL);
 tmp = dvrr_stack + 2154;
 target_ptr = Libderiv->deriv2_classes[1][3][130];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+1102, dvrr_stack+3086, NULL);
 tmp = dvrr_stack + 1102;
 target_ptr = Libderiv->deriv2_classes[2][2][130];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+5876, dvrr_stack+3194, NULL);
 tmp = dvrr_stack + 5876;
 target_ptr = Libderiv->deriv2_classes[2][3][130];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,18,dvrr_stack+1403, dvrr_stack+444, NULL);
 tmp = dvrr_stack + 1403;
 target_ptr = Libderiv->deriv2_classes[1][2][119];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,30,dvrr_stack+444, dvrr_stack+952, NULL);
 tmp = dvrr_stack + 444;
 target_ptr = Libderiv->deriv2_classes[1][3][119];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,36,dvrr_stack+952, dvrr_stack+1834, NULL);
 tmp = dvrr_stack + 952;
 target_ptr = Libderiv->deriv2_classes[2][2][119];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,60,dvrr_stack+1834, dvrr_stack+2816, NULL);
 tmp = dvrr_stack + 1834;
 target_ptr = Libderiv->deriv2_classes[2][3][119];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+2816, dvrr_stack+1976, NULL);
 tmp = dvrr_stack + 2816;
 target_ptr = Libderiv->deriv2_classes[1][2][118];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+2834, dvrr_stack+2996, NULL);
 tmp = dvrr_stack + 2834;
 target_ptr = Libderiv->deriv2_classes[1][3][118];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+2864, dvrr_stack+3086, NULL);
 tmp = dvrr_stack + 2864;
 target_ptr = Libderiv->deriv2_classes[2][2][118];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+2900, dvrr_stack+3194, NULL);
 tmp = dvrr_stack + 2900;
 target_ptr = Libderiv->deriv2_classes[2][3][118];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+2960, dvrr_stack+3374, NULL);
 tmp = dvrr_stack + 2960;
 target_ptr = Libderiv->deriv2_classes[1][2][117];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+2978, dvrr_stack+336, NULL);
 tmp = dvrr_stack + 2978;
 target_ptr = Libderiv->deriv2_classes[1][3][117];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+3008, dvrr_stack+772, NULL);
 tmp = dvrr_stack + 3008;
 target_ptr = Libderiv->deriv2_classes[2][2][117];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+3044, dvrr_stack+1618, NULL);
 tmp = dvrr_stack + 3044;
 target_ptr = Libderiv->deriv2_classes[2][3][117];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+1618, dvrr_stack+2030, dvrr_stack+2483);
 tmp = dvrr_stack + 1618;
 target_ptr = Libderiv->deriv2_classes[1][2][107];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+1636, dvrr_stack+2492, dvrr_stack+426);
 tmp = dvrr_stack + 1636;
 target_ptr = Libderiv->deriv2_classes[1][3][107];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+1666, dvrr_stack+2627, dvrr_stack+2609);
 tmp = dvrr_stack + 1666;
 target_ptr = Libderiv->deriv2_classes[2][2][107];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+1702, dvrr_stack+2687, dvrr_stack+1798);
 tmp = dvrr_stack + 1702;
 target_ptr = Libderiv->deriv2_classes[2][3][107];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+1762, dvrr_stack+2786, dvrr_stack+2777);
 tmp = dvrr_stack + 1762;
 target_ptr = Libderiv->deriv2_classes[1][2][106];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+3104, dvrr_stack+898, dvrr_stack+880);
 tmp = dvrr_stack + 3104;
 target_ptr = Libderiv->deriv2_classes[1][3][106];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+3134, dvrr_stack+3428, dvrr_stack+532);
 tmp = dvrr_stack + 3134;
 target_ptr = Libderiv->deriv2_classes[2][2][106];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+3170, dvrr_stack+3524, dvrr_stack+3488);
 tmp = dvrr_stack + 3170;
 target_ptr = Libderiv->deriv2_classes[2][3][106];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+1780, dvrr_stack+3614, dvrr_stack+943);
 tmp = dvrr_stack + 1780;
 target_ptr = Libderiv->deriv2_classes[1][2][105];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+3230, dvrr_stack+98, dvrr_stack+246);
 tmp = dvrr_stack + 3230;
 target_ptr = Libderiv->deriv2_classes[1][3][105];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+3260, dvrr_stack+655, dvrr_stack+637);
 tmp = dvrr_stack + 3260;
 target_ptr = Libderiv->deriv2_classes[2][2][105];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+3296, dvrr_stack+1195, dvrr_stack+1438);
 tmp = dvrr_stack + 3296;
 target_ptr = Libderiv->deriv2_classes[2][3][105];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+3356, dvrr_stack+2189, dvrr_stack+143);
 tmp = dvrr_stack + 3356;
 target_ptr = Libderiv->deriv2_classes[1][2][104];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+3374, dvrr_stack+2219, dvrr_stack+1285);
 tmp = dvrr_stack + 3374;
 target_ptr = Libderiv->deriv2_classes[1][3][104];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+1894, dvrr_stack+2288, dvrr_stack+2270);
 tmp = dvrr_stack + 1894;
 target_ptr = Libderiv->deriv2_classes[2][2][104];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+1948, dvrr_stack+2384, dvrr_stack+2348);
 tmp = dvrr_stack + 1948;
 target_ptr = Libderiv->deriv2_classes[2][3][104];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+3404, dvrr_stack+2030, dvrr_stack+2483);
 tmp = dvrr_stack + 3404;
 target_ptr = Libderiv->deriv2_classes[1][2][95];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+988, dvrr_stack+2492, dvrr_stack+426);
 tmp = dvrr_stack + 988;
 target_ptr = Libderiv->deriv2_classes[1][3][95];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+474, dvrr_stack+2627, dvrr_stack+2609);
 tmp = dvrr_stack + 474;
 target_ptr = Libderiv->deriv2_classes[2][2][95];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+745, dvrr_stack+2687, dvrr_stack+1798);
 tmp = dvrr_stack + 745;
 target_ptr = Libderiv->deriv2_classes[2][3][95];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+1018, dvrr_stack+2786, dvrr_stack+2777);
 tmp = dvrr_stack + 1018;
 target_ptr = Libderiv->deriv2_classes[1][2][94];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1138, dvrr_stack+898, dvrr_stack+880);
 tmp = dvrr_stack + 1138;
 target_ptr = Libderiv->deriv2_classes[1][3][94];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+805, dvrr_stack+3428, dvrr_stack+532);
 tmp = dvrr_stack + 805;
 target_ptr = Libderiv->deriv2_classes[2][2][94];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+324, dvrr_stack+3524, dvrr_stack+3488);
 tmp = dvrr_stack + 324;
 target_ptr = Libderiv->deriv2_classes[2][3][94];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+510, dvrr_stack+3614, dvrr_stack+943);
 tmp = dvrr_stack + 510;
 target_ptr = Libderiv->deriv2_classes[1][2][93];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+841, dvrr_stack+98, dvrr_stack+246);
 tmp = dvrr_stack + 841;
 target_ptr = Libderiv->deriv2_classes[1][3][93];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+384, dvrr_stack+655, dvrr_stack+637);
 tmp = dvrr_stack + 384;
 target_ptr = Libderiv->deriv2_classes[2][2][93];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+5936, dvrr_stack+1195, dvrr_stack+1438);
 tmp = dvrr_stack + 5936;
 target_ptr = Libderiv->deriv2_classes[2][3][93];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+4351, dvrr_stack+2189, dvrr_stack+143);
 tmp = dvrr_stack + 4351;
 target_ptr = Libderiv->deriv2_classes[1][2][92];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+5996, dvrr_stack+2219, dvrr_stack+1285);
 tmp = dvrr_stack + 5996;
 target_ptr = Libderiv->deriv2_classes[1][3][92];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+6026, dvrr_stack+2288, dvrr_stack+2270);
 tmp = dvrr_stack + 6026;
 target_ptr = Libderiv->deriv2_classes[2][2][92];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+6062, dvrr_stack+2384, dvrr_stack+2348);
 tmp = dvrr_stack + 6062;
 target_ptr = Libderiv->deriv2_classes[2][3][92];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+152, dvrr_stack+1474, dvrr_stack+2474);
 tmp = dvrr_stack + 152;
 target_ptr = Libderiv->deriv2_classes[1][2][91];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+6122, dvrr_stack+1522, dvrr_stack+1504);
 tmp = dvrr_stack + 6122;
 target_ptr = Libderiv->deriv2_classes[1][3][91];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+6152, dvrr_stack+264, dvrr_stack+1567);
 tmp = dvrr_stack + 6152;
 target_ptr = Libderiv->deriv2_classes[2][2][91];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+6188, dvrr_stack+3644, dvrr_stack+2555);
 tmp = dvrr_stack + 6188;
 target_ptr = Libderiv->deriv2_classes[2][3][91];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+1168, dvrr_stack+2030, dvrr_stack+2483);
 tmp = dvrr_stack + 1168;
 target_ptr = Libderiv->deriv2_classes[1][2][83];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+6248, dvrr_stack+2492, dvrr_stack+426);
 tmp = dvrr_stack + 6248;
 target_ptr = Libderiv->deriv2_classes[1][3][83];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+2483, dvrr_stack+2627, dvrr_stack+2609);
 tmp = dvrr_stack + 2483;
 target_ptr = Libderiv->deriv2_classes[2][2][83];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+6278, dvrr_stack+2687, dvrr_stack+1798);
 tmp = dvrr_stack + 6278;
 target_ptr = Libderiv->deriv2_classes[2][3][83];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+2609, dvrr_stack+2786, dvrr_stack+2777);
 tmp = dvrr_stack + 2609;
 target_ptr = Libderiv->deriv2_classes[1][2][82];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2687, dvrr_stack+898, dvrr_stack+880);
 tmp = dvrr_stack + 2687;
 target_ptr = Libderiv->deriv2_classes[1][3][82];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+898, dvrr_stack+3428, dvrr_stack+532);
 tmp = dvrr_stack + 898;
 target_ptr = Libderiv->deriv2_classes[2][2][82];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+2717, dvrr_stack+3524, dvrr_stack+3488);
 tmp = dvrr_stack + 2717;
 target_ptr = Libderiv->deriv2_classes[2][3][82];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+2519, dvrr_stack+3614, dvrr_stack+943);
 tmp = dvrr_stack + 2519;
 target_ptr = Libderiv->deriv2_classes[1][2][81];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+3524, dvrr_stack+98, dvrr_stack+246);
 tmp = dvrr_stack + 3524;
 target_ptr = Libderiv->deriv2_classes[1][3][81];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+3554, dvrr_stack+655, dvrr_stack+637);
 tmp = dvrr_stack + 3554;
 target_ptr = Libderiv->deriv2_classes[2][2][81];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+6338, dvrr_stack+1195, dvrr_stack+1438);
 tmp = dvrr_stack + 6338;
 target_ptr = Libderiv->deriv2_classes[2][3][81];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+934, dvrr_stack+2189, dvrr_stack+143);
 tmp = dvrr_stack + 934;
 target_ptr = Libderiv->deriv2_classes[1][2][80];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+86, dvrr_stack+2219, dvrr_stack+1285);
 tmp = dvrr_stack + 86;
 target_ptr = Libderiv->deriv2_classes[1][3][80];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+116, dvrr_stack+2288, dvrr_stack+2270);
 tmp = dvrr_stack + 116;
 target_ptr = Libderiv->deriv2_classes[2][2][80];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+1186, dvrr_stack+2384, dvrr_stack+2348);
 tmp = dvrr_stack + 1186;
 target_ptr = Libderiv->deriv2_classes[2][3][80];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+2270, dvrr_stack+1474, dvrr_stack+2474);
 tmp = dvrr_stack + 2270;
 target_ptr = Libderiv->deriv2_classes[1][2][79];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2384, dvrr_stack+1522, dvrr_stack+1504);
 tmp = dvrr_stack + 2384;
 target_ptr = Libderiv->deriv2_classes[1][3][79];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1522, dvrr_stack+264, dvrr_stack+1567);
 tmp = dvrr_stack + 1522;
 target_ptr = Libderiv->deriv2_classes[2][2][79];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+2414, dvrr_stack+3644, dvrr_stack+2555);
 tmp = dvrr_stack + 2414;
 target_ptr = Libderiv->deriv2_classes[2][3][79];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+3644, dvrr_stack+715, dvrr_stack+1585);
 tmp = dvrr_stack + 3644;
 target_ptr = Libderiv->deriv2_classes[1][2][78];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+3662, dvrr_stack+3734, dvrr_stack+2591);
 tmp = dvrr_stack + 3662;
 target_ptr = Libderiv->deriv2_classes[1][3][78];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1558, dvrr_stack+550, dvrr_stack+1594);
 tmp = dvrr_stack + 1558;
 target_ptr = Libderiv->deriv2_classes[2][2][78];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+3692, dvrr_stack+3815, dvrr_stack+3779);
 tmp = dvrr_stack + 3692;
 target_ptr = Libderiv->deriv2_classes[2][3][78];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,6,dvrr_stack+1594, dvrr_stack+1798, dvrr_stack+2264);
 tmp = dvrr_stack + 1594;
 target_ptr = Libderiv->deriv2_classes[1][2][35];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,10,dvrr_stack+3815, dvrr_stack+2627, dvrr_stack+2108);
 tmp = dvrr_stack + 3815;
 target_ptr = Libderiv->deriv2_classes[1][3][35];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,6,dvrr_stack+3845, dvrr_stack+4371, dvrr_stack+426);
 tmp = dvrr_stack + 3845;
 target_ptr = Libderiv->deriv2_classes[2][2][35];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6398, dvrr_stack+1303, dvrr_stack+2030);
 tmp = dvrr_stack + 6398;
 target_ptr = Libderiv->deriv2_classes[2][3][35];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+3881, dvrr_stack+3488, dvrr_stack+6);
 tmp = dvrr_stack + 3881;
 target_ptr = Libderiv->deriv2_classes[1][2][34];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+2219, dvrr_stack+3428, dvrr_stack+48);
 tmp = dvrr_stack + 2219;
 target_ptr = Libderiv->deriv2_classes[1][3][34];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1246, dvrr_stack+4431, dvrr_stack+880);
 tmp = dvrr_stack + 1246;
 target_ptr = Libderiv->deriv2_classes[2][2][34];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6458, dvrr_stack+5016, dvrr_stack+2786);
 tmp = dvrr_stack + 6458;
 target_ptr = Libderiv->deriv2_classes[2][3][34];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+3752, dvrr_stack+1438, dvrr_stack+1612);
 tmp = dvrr_stack + 3752;
 target_ptr = Libderiv->deriv2_classes[1][2][33];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+610, dvrr_stack+655, dvrr_stack+2060);
 tmp = dvrr_stack + 610;
 target_ptr = Libderiv->deriv2_classes[1][3][33];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+6518, dvrr_stack+4491, dvrr_stack+246);
 tmp = dvrr_stack + 6518;
 target_ptr = Libderiv->deriv2_classes[2][2][33];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6554, dvrr_stack+4191, dvrr_stack+3614);
 tmp = dvrr_stack + 6554;
 target_ptr = Libderiv->deriv2_classes[2][3][33];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+3590, dvrr_stack+2348, dvrr_stack+4716);
 tmp = dvrr_stack + 3590;
 target_ptr = Libderiv->deriv2_classes[1][2][32];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+6614, dvrr_stack+2288, dvrr_stack+4722);
 tmp = dvrr_stack + 6614;
 target_ptr = Libderiv->deriv2_classes[1][3][32];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+6644, dvrr_stack+4768, dvrr_stack+1285);
 tmp = dvrr_stack + 6644;
 target_ptr = Libderiv->deriv2_classes[2][2][32];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6680, dvrr_stack+4828, dvrr_stack+2189);
 tmp = dvrr_stack + 6680;
 target_ptr = Libderiv->deriv2_classes[2][3][32];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+528, dvrr_stack+2555, dvrr_stack+4732);
 tmp = dvrr_stack + 528;
 target_ptr = Libderiv->deriv2_classes[1][2][31];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+6740, dvrr_stack+264, dvrr_stack+2537);
 tmp = dvrr_stack + 6740;
 target_ptr = Libderiv->deriv2_classes[1][3][31];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+6770, dvrr_stack+4928, dvrr_stack+1504);
 tmp = dvrr_stack + 6770;
 target_ptr = Libderiv->deriv2_classes[2][2][31];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6806, dvrr_stack+5116, dvrr_stack+1474);
 tmp = dvrr_stack + 6806;
 target_ptr = Libderiv->deriv2_classes[2][3][31];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+2008, dvrr_stack+3779, dvrr_stack+1942);
 tmp = dvrr_stack + 2008;
 target_ptr = Libderiv->deriv2_classes[1][2][30];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+6866, dvrr_stack+550, dvrr_stack+4988);
 tmp = dvrr_stack + 6866;
 target_ptr = Libderiv->deriv2_classes[1][3][30];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+6896, dvrr_stack+4291, dvrr_stack+2591);
 tmp = dvrr_stack + 6896;
 target_ptr = Libderiv->deriv2_classes[2][2][30];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6932, dvrr_stack+5216, dvrr_stack+715);
 tmp = dvrr_stack + 6932;
 target_ptr = Libderiv->deriv2_classes[2][3][30];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+6992, dvrr_stack+2070, dvrr_stack+4738);
 tmp = dvrr_stack + 6992;
 target_ptr = Libderiv->deriv2_classes[1][2][26];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+7010, dvrr_stack+171, dvrr_stack+4744);
 tmp = dvrr_stack + 7010;
 target_ptr = Libderiv->deriv2_classes[1][3][26];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+7040, dvrr_stack+5406, dvrr_stack+4998);
 tmp = dvrr_stack + 7040;
 target_ptr = Libderiv->deriv2_classes[2][2][26];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+7076, dvrr_stack+4581, dvrr_stack+4551);
 tmp = dvrr_stack + 7076;
 target_ptr = Libderiv->deriv2_classes[2][3][26];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,6,dvrr_stack+7136, dvrr_stack+1798, dvrr_stack+2264);
 tmp = dvrr_stack + 7136;
 target_ptr = Libderiv->deriv2_classes[1][2][23];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,10,dvrr_stack+7154, dvrr_stack+2627, dvrr_stack+2108);
 tmp = dvrr_stack + 7154;
 target_ptr = Libderiv->deriv2_classes[1][3][23];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,6,dvrr_stack+7184, dvrr_stack+4371, dvrr_stack+426);
 tmp = dvrr_stack + 7184;
 target_ptr = Libderiv->deriv2_classes[2][2][23];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,10,dvrr_stack+7220, dvrr_stack+1303, dvrr_stack+2030);
 tmp = dvrr_stack + 7220;
 target_ptr = Libderiv->deriv2_classes[2][3][23];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+7280, dvrr_stack+3488, dvrr_stack+6);
 tmp = dvrr_stack + 7280;
 target_ptr = Libderiv->deriv2_classes[1][2][22];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+7298, dvrr_stack+3428, dvrr_stack+48);
 tmp = dvrr_stack + 7298;
 target_ptr = Libderiv->deriv2_classes[1][3][22];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+7328, dvrr_stack+4431, dvrr_stack+880);
 tmp = dvrr_stack + 7328;
 target_ptr = Libderiv->deriv2_classes[2][2][22];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+7364, dvrr_stack+5016, dvrr_stack+2786);
 tmp = dvrr_stack + 7364;
 target_ptr = Libderiv->deriv2_classes[2][3][22];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+7424, dvrr_stack+1438, dvrr_stack+1612);
 tmp = dvrr_stack + 7424;
 target_ptr = Libderiv->deriv2_classes[1][2][21];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+7442, dvrr_stack+655, dvrr_stack+2060);
 tmp = dvrr_stack + 7442;
 target_ptr = Libderiv->deriv2_classes[1][3][21];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+7472, dvrr_stack+4491, dvrr_stack+246);
 tmp = dvrr_stack + 7472;
 target_ptr = Libderiv->deriv2_classes[2][2][21];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+7508, dvrr_stack+4191, dvrr_stack+3614);
 tmp = dvrr_stack + 7508;
 target_ptr = Libderiv->deriv2_classes[2][3][21];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+7568, dvrr_stack+2348, dvrr_stack+4716);
 tmp = dvrr_stack + 7568;
 target_ptr = Libderiv->deriv2_classes[1][2][20];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+7586, dvrr_stack+2288, dvrr_stack+4722);
 tmp = dvrr_stack + 7586;
 target_ptr = Libderiv->deriv2_classes[1][3][20];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+7616, dvrr_stack+4768, dvrr_stack+1285);
 tmp = dvrr_stack + 7616;
 target_ptr = Libderiv->deriv2_classes[2][2][20];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+7652, dvrr_stack+4828, dvrr_stack+2189);
 tmp = dvrr_stack + 7652;
 target_ptr = Libderiv->deriv2_classes[2][3][20];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+7712, dvrr_stack+2555, dvrr_stack+4732);
 tmp = dvrr_stack + 7712;
 target_ptr = Libderiv->deriv2_classes[1][2][19];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+7730, dvrr_stack+264, dvrr_stack+2537);
 tmp = dvrr_stack + 7730;
 target_ptr = Libderiv->deriv2_classes[1][3][19];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+7760, dvrr_stack+4928, dvrr_stack+1504);
 tmp = dvrr_stack + 7760;
 target_ptr = Libderiv->deriv2_classes[2][2][19];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+7796, dvrr_stack+5116, dvrr_stack+1474);
 tmp = dvrr_stack + 7796;
 target_ptr = Libderiv->deriv2_classes[2][3][19];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+7856, dvrr_stack+3779, dvrr_stack+1942);
 tmp = dvrr_stack + 7856;
 target_ptr = Libderiv->deriv2_classes[1][2][18];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+7874, dvrr_stack+550, dvrr_stack+4988);
 tmp = dvrr_stack + 7874;
 target_ptr = Libderiv->deriv2_classes[1][3][18];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+7904, dvrr_stack+4291, dvrr_stack+2591);
 tmp = dvrr_stack + 7904;
 target_ptr = Libderiv->deriv2_classes[2][2][18];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+7940, dvrr_stack+5216, dvrr_stack+715);
 tmp = dvrr_stack + 7940;
 target_ptr = Libderiv->deriv2_classes[2][3][18];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+8000, dvrr_stack+2070, dvrr_stack+4738);
 tmp = dvrr_stack + 8000;
 target_ptr = Libderiv->deriv2_classes[1][2][14];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+8018, dvrr_stack+171, dvrr_stack+4744);
 tmp = dvrr_stack + 8018;
 target_ptr = Libderiv->deriv2_classes[1][3][14];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+8048, dvrr_stack+5406, dvrr_stack+4998);
 tmp = dvrr_stack + 8048;
 target_ptr = Libderiv->deriv2_classes[2][2][14];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+8084, dvrr_stack+4581, dvrr_stack+4551);
 tmp = dvrr_stack + 8084;
 target_ptr = Libderiv->deriv2_classes[2][3][14];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+8144, dvrr_stack+3905, dvrr_stack+0);
 tmp = dvrr_stack + 8144;
 target_ptr = Libderiv->deriv2_classes[1][2][13];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+8162, dvrr_stack+4001, dvrr_stack+4681);
 tmp = dvrr_stack + 8162;
 target_ptr = Libderiv->deriv2_classes[1][3][13];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+8192, dvrr_stack+5466, dvrr_stack+4691);
 tmp = dvrr_stack + 8192;
 target_ptr = Libderiv->deriv2_classes[2][2][13];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+8228, dvrr_stack+5526, dvrr_stack+4061);
 tmp = dvrr_stack + 8228;
 target_ptr = Libderiv->deriv2_classes[2][3][13];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,6,dvrr_stack+8288, dvrr_stack+1798, dvrr_stack+2264);
 tmp = dvrr_stack + 8288;
 target_ptr = Libderiv->deriv2_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,10,dvrr_stack+1798, dvrr_stack+2627, dvrr_stack+2108);
 tmp = dvrr_stack + 1798;
 target_ptr = Libderiv->deriv2_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,6,dvrr_stack+2627, dvrr_stack+4371, dvrr_stack+426);
 tmp = dvrr_stack + 2627;
 target_ptr = Libderiv->deriv2_classes[2][2][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,10,dvrr_stack+4369, dvrr_stack+1303, dvrr_stack+2030);
 tmp = dvrr_stack + 4369;
 target_ptr = Libderiv->deriv2_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1303, dvrr_stack+3488, dvrr_stack+6);
 tmp = dvrr_stack + 1303;
 target_ptr = Libderiv->deriv2_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3488, dvrr_stack+3428, dvrr_stack+48);
 tmp = dvrr_stack + 3488;
 target_ptr = Libderiv->deriv2_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1321, dvrr_stack+4431, dvrr_stack+880);
 tmp = dvrr_stack + 1321;
 target_ptr = Libderiv->deriv2_classes[2][2][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+3422, dvrr_stack+5016, dvrr_stack+2786);
 tmp = dvrr_stack + 3422;
 target_ptr = Libderiv->deriv2_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+5016, dvrr_stack+1438, dvrr_stack+1612);
 tmp = dvrr_stack + 5016;
 target_ptr = Libderiv->deriv2_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5034, dvrr_stack+655, dvrr_stack+2060);
 tmp = dvrr_stack + 5034;
 target_ptr = Libderiv->deriv2_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+5064, dvrr_stack+4491, dvrr_stack+246);
 tmp = dvrr_stack + 5064;
 target_ptr = Libderiv->deriv2_classes[2][2][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+4429, dvrr_stack+4191, dvrr_stack+3614);
 tmp = dvrr_stack + 4429;
 target_ptr = Libderiv->deriv2_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1357, dvrr_stack+2348, dvrr_stack+4716);
 tmp = dvrr_stack + 1357;
 target_ptr = Libderiv->deriv2_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+2348, dvrr_stack+2288, dvrr_stack+4722);
 tmp = dvrr_stack + 2348;
 target_ptr = Libderiv->deriv2_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+3608, dvrr_stack+4768, dvrr_stack+1285);
 tmp = dvrr_stack + 3608;
 target_ptr = Libderiv->deriv2_classes[2][2][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+2288, dvrr_stack+4828, dvrr_stack+2189);
 tmp = dvrr_stack + 2288;
 target_ptr = Libderiv->deriv2_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1375, dvrr_stack+2555, dvrr_stack+4732);
 tmp = dvrr_stack + 1375;
 target_ptr = Libderiv->deriv2_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+2777, dvrr_stack+264, dvrr_stack+2537);
 tmp = dvrr_stack + 2777;
 target_ptr = Libderiv->deriv2_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2537, dvrr_stack+4928, dvrr_stack+1504);
 tmp = dvrr_stack + 2537;
 target_ptr = Libderiv->deriv2_classes[2][2][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+4187, dvrr_stack+5116, dvrr_stack+1474);
 tmp = dvrr_stack + 4187;
 target_ptr = Libderiv->deriv2_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2573, dvrr_stack+3779, dvrr_stack+1942);
 tmp = dvrr_stack + 2573;
 target_ptr = Libderiv->deriv2_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5100, dvrr_stack+550, dvrr_stack+4988);
 tmp = dvrr_stack + 5100;
 target_ptr = Libderiv->deriv2_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+5130, dvrr_stack+4291, dvrr_stack+2591);
 tmp = dvrr_stack + 5130;
 target_ptr = Libderiv->deriv2_classes[2][2][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+546, dvrr_stack+5216, dvrr_stack+715);
 tmp = dvrr_stack + 546;
 target_ptr = Libderiv->deriv2_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2591, dvrr_stack+2070, dvrr_stack+4738);
 tmp = dvrr_stack + 2591;
 target_ptr = Libderiv->deriv2_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5166, dvrr_stack+171, dvrr_stack+4744);
 tmp = dvrr_stack + 5166;
 target_ptr = Libderiv->deriv2_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+5196, dvrr_stack+5406, dvrr_stack+4998);
 tmp = dvrr_stack + 5196;
 target_ptr = Libderiv->deriv2_classes[2][2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+5406, dvrr_stack+4581, dvrr_stack+4551);
 tmp = dvrr_stack + 5406;
 target_ptr = Libderiv->deriv2_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1930, dvrr_stack+3905, dvrr_stack+0);
 tmp = dvrr_stack + 1930;
 target_ptr = Libderiv->deriv2_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5232, dvrr_stack+4001, dvrr_stack+4681);
 tmp = dvrr_stack + 5232;
 target_ptr = Libderiv->deriv2_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+4001, dvrr_stack+5466, dvrr_stack+4691);
 tmp = dvrr_stack + 4001;
 target_ptr = Libderiv->deriv2_classes[2][2][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+5466, dvrr_stack+5526, dvrr_stack+4061);
 tmp = dvrr_stack + 5466;
 target_ptr = Libderiv->deriv2_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+5526, dvrr_stack+2118, dvrr_stack+4709);
 tmp = dvrr_stack + 5526;
 target_ptr = Libderiv->deriv2_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5544, dvrr_stack+3941, dvrr_stack+21);
 tmp = dvrr_stack + 5544;
 target_ptr = Libderiv->deriv2_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+5574, dvrr_stack+4091, dvrr_stack+68);
 tmp = dvrr_stack + 5574;
 target_ptr = Libderiv->deriv2_classes[2][2][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+4037, dvrr_stack+5776, dvrr_stack+5316);
 tmp = dvrr_stack + 4037;
 target_ptr = Libderiv->deriv2_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];


}

