#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (pp|d0) integrals */

void d12vrr_order_ppd0(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+498, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+2);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+507, dvrr_stack+6, dvrr_stack+42, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+525, dvrr_stack+21, dvrr_stack+507, dvrr_stack+15, dvrr_stack+6, dvrr_stack+498);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+561, dvrr_stack+48, dvrr_stack+161, NULL, NULL, dvrr_stack+42);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+591, dvrr_stack+68, dvrr_stack+561, dvrr_stack+58, dvrr_stack+48, dvrr_stack+507);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+651,dvrr_stack+591,dvrr_stack+525,6);


 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+759, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+762, dvrr_stack+152, dvrr_stack+759, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+768, dvrr_stack+155, dvrr_stack+762, dvrr_stack+39, dvrr_stack+152, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+778, dvrr_stack+161, dvrr_stack+768, dvrr_stack+42, dvrr_stack+155, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+793, dvrr_stack+171, dvrr_stack+778, NULL, NULL, dvrr_stack+161);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+838, dvrr_stack+201, dvrr_stack+793, dvrr_stack+186, dvrr_stack+171, dvrr_stack+561);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+928,dvrr_stack+838,dvrr_stack+591,6);


 /* compute (2 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+1108,dvrr_stack+928,dvrr_stack+651,6);


 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,36,dvrr_stack+1324, dvrr_stack+1108, dvrr_stack+525);

 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,18,dvrr_stack+778, dvrr_stack+336, dvrr_stack+21);

 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,36,dvrr_stack+1432, dvrr_stack+1108, dvrr_stack+525);

 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,18,dvrr_stack+1540, dvrr_stack+336, dvrr_stack+21);

 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,36,dvrr_stack+336, dvrr_stack+1108, dvrr_stack+525);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1108, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+1);

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+1117,dvrr_stack+21,dvrr_stack+1108,3);


 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+1144, dvrr_stack+1117, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+171, dvrr_stack+246, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+759, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+1153, dvrr_stack+1108, dvrr_stack+498, dvrr_stack+12, dvrr_stack+0, dvrr_stack+759);

 /* compute (2 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+1171,dvrr_stack+525,dvrr_stack+1153,6);


 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+1225, dvrr_stack+1171, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+1243, dvrr_stack+928, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+1303, dvrr_stack+1117, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+1594, dvrr_stack+246, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+1624, dvrr_stack+1171, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+1642, dvrr_stack+928, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+1312, dvrr_stack+1117, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+1702, dvrr_stack+246, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+246, dvrr_stack+1171, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+264, dvrr_stack+928, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+1321, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+928, dvrr_stack+21, dvrr_stack+1321);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+937, dvrr_stack+201, dvrr_stack+21);

 /* compute (2 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+832, dvrr_stack+1321, dvrr_stack+759, Data->F+0, Data->F+1, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+967, dvrr_stack+525, dvrr_stack+832);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+985, dvrr_stack+838, dvrr_stack+525);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+1045, dvrr_stack+21, dvrr_stack+1321);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1054, dvrr_stack+201, dvrr_stack+21);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+1084, dvrr_stack+525, dvrr_stack+832);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+1732, dvrr_stack+838, dvrr_stack+525);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+1171, dvrr_stack+21, dvrr_stack+1321);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+1180, dvrr_stack+201, dvrr_stack+21);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+201, dvrr_stack+525, dvrr_stack+832);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+1792, dvrr_stack+838, dvrr_stack+525);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+832,dvrr_stack+58,dvrr_stack+15,1);


 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+1102, dvrr_stack+832, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+850, dvrr_stack+651, NULL);
 tmp = dvrr_stack + 850;
 target_ptr = Libderiv->deriv_classes[2][2][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+886, dvrr_stack+98, NULL);
 tmp = dvrr_stack + 886;
 target_ptr = Libderiv->deriv_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+1321, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+904, dvrr_stack+3, dvrr_stack+39, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+219, dvrr_stack+498, dvrr_stack+904, dvrr_stack+0, dvrr_stack+3, dvrr_stack+1321);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1117, dvrr_stack+42, dvrr_stack+155, NULL, NULL, dvrr_stack+39);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1852, dvrr_stack+507, dvrr_stack+1117, dvrr_stack+6, dvrr_stack+42, dvrr_stack+904);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+1888, dvrr_stack+525, dvrr_stack+1852, dvrr_stack+21, dvrr_stack+507, dvrr_stack+219);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1948, dvrr_stack+161, dvrr_stack+768, NULL, NULL, dvrr_stack+155);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1978, dvrr_stack+561, dvrr_stack+1948, dvrr_stack+48, dvrr_stack+161, dvrr_stack+1117);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+2038, dvrr_stack+591, dvrr_stack+1978, dvrr_stack+68, dvrr_stack+561, dvrr_stack+1852);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+2138,dvrr_stack+2038,dvrr_stack+1888,10);


 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+1948, dvrr_stack+2138, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+6, dvrr_stack+832, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+2318, dvrr_stack+651, NULL);
 tmp = dvrr_stack + 2318;
 target_ptr = Libderiv->deriv_classes[2][2][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+561, dvrr_stack+98, NULL);
 tmp = dvrr_stack + 561;
 target_ptr = Libderiv->deriv_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+2354, dvrr_stack+2138, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+579, dvrr_stack+832, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+2414, dvrr_stack+651, NULL);
 tmp = dvrr_stack + 2414;
 target_ptr = Libderiv->deriv_classes[2][2][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+832, dvrr_stack+98, NULL);
 tmp = dvrr_stack + 832;
 target_ptr = Libderiv->deriv_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+651, dvrr_stack+2138, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+585, dvrr_stack+58, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+2138, dvrr_stack+591, dvrr_stack+1153);
 tmp = dvrr_stack + 2138;
 target_ptr = Libderiv->deriv_classes[2][2][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+2174, dvrr_stack+68, dvrr_stack+1108);
 tmp = dvrr_stack + 2174;
 target_ptr = Libderiv->deriv_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+2192, dvrr_stack+759, dvrr_stack+1321, Data->F+1, Data->F+2, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+2008, dvrr_stack+1153, dvrr_stack+219, dvrr_stack+1108, dvrr_stack+498, dvrr_stack+2192);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+2192, dvrr_stack+2038, dvrr_stack+2008);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+2252, dvrr_stack+58, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+2258, dvrr_stack+591, dvrr_stack+1153);
 tmp = dvrr_stack + 2258;
 target_ptr = Libderiv->deriv_classes[2][2][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+2294, dvrr_stack+68, dvrr_stack+1108);
 tmp = dvrr_stack + 2294;
 target_ptr = Libderiv->deriv_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+2450, dvrr_stack+2038, dvrr_stack+2008);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+2312, dvrr_stack+58, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+98, dvrr_stack+591, dvrr_stack+1153);
 tmp = dvrr_stack + 98;
 target_ptr = Libderiv->deriv_classes[2][2][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+1153, dvrr_stack+68, dvrr_stack+1108);
 tmp = dvrr_stack + 1153;
 target_ptr = Libderiv->deriv_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+591, dvrr_stack+2038, dvrr_stack+2008);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+2008, dvrr_stack+21, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+2014, dvrr_stack+1888, dvrr_stack+21);
 tmp = dvrr_stack + 2014;
 target_ptr = Libderiv->deriv_classes[2][2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+134, dvrr_stack+525, dvrr_stack+15);
 tmp = dvrr_stack + 134;
 target_ptr = Libderiv->deriv_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+12, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+2050, dvrr_stack+1321, dvrr_stack+12, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1108, dvrr_stack+39, dvrr_stack+152, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+2056, dvrr_stack+904, dvrr_stack+1108, dvrr_stack+3, dvrr_stack+39, dvrr_stack+12);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+2074, dvrr_stack+219, dvrr_stack+2056, dvrr_stack+498, dvrr_stack+904, dvrr_stack+2050);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+904, dvrr_stack+155, dvrr_stack+762, NULL, NULL, dvrr_stack+152);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+711, dvrr_stack+1117, dvrr_stack+904, dvrr_stack+42, dvrr_stack+155, dvrr_stack+1108);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+2510, dvrr_stack+1852, dvrr_stack+711, dvrr_stack+507, dvrr_stack+1117, dvrr_stack+2056);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+2570, dvrr_stack+1888, dvrr_stack+2510, dvrr_stack+525, dvrr_stack+1852, dvrr_stack+2074);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+2510, dvrr_stack+2570, dvrr_stack+525);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+0, dvrr_stack+21, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+1852, dvrr_stack+1888, dvrr_stack+21);
 tmp = dvrr_stack + 1852;
 target_ptr = Libderiv->deriv_classes[2][2][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+711, dvrr_stack+525, dvrr_stack+15);
 tmp = dvrr_stack + 711;
 target_ptr = Libderiv->deriv_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+2050, dvrr_stack+2570, dvrr_stack+525);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+729, dvrr_stack+21, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1108, dvrr_stack+1888, dvrr_stack+21);
 tmp = dvrr_stack + 1108;
 target_ptr = Libderiv->deriv_classes[2][2][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+21, dvrr_stack+525, dvrr_stack+15);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+1888, dvrr_stack+2570, dvrr_stack+525);

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,18,dvrr_stack+2570, dvrr_stack+444, NULL);
 tmp = dvrr_stack + 2570;
 target_ptr = Libderiv->deriv2_classes[1][2][143];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,36,dvrr_stack+2588, dvrr_stack+1324, NULL);
 tmp = dvrr_stack + 2588;
 target_ptr = Libderiv->deriv2_classes[2][2][143];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,18,dvrr_stack+2624, dvrr_stack+444, NULL);
 tmp = dvrr_stack + 2624;
 target_ptr = Libderiv->deriv2_classes[1][2][131];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,36,dvrr_stack+2642, dvrr_stack+1324, NULL);
 tmp = dvrr_stack + 2642;
 target_ptr = Libderiv->deriv2_classes[2][2][131];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+2678, dvrr_stack+778, NULL);
 tmp = dvrr_stack + 2678;
 target_ptr = Libderiv->deriv2_classes[1][2][130];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+2696, dvrr_stack+1432, NULL);
 tmp = dvrr_stack + 2696;
 target_ptr = Libderiv->deriv2_classes[2][2][130];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,18,dvrr_stack+2732, dvrr_stack+444, NULL);
 tmp = dvrr_stack + 2732;
 target_ptr = Libderiv->deriv2_classes[1][2][119];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,36,dvrr_stack+444, dvrr_stack+1324, NULL);
 tmp = dvrr_stack + 444;
 target_ptr = Libderiv->deriv2_classes[2][2][119];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+480, dvrr_stack+778, NULL);
 tmp = dvrr_stack + 480;
 target_ptr = Libderiv->deriv2_classes[1][2][118];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+498, dvrr_stack+1432, NULL);
 tmp = dvrr_stack + 498;
 target_ptr = Libderiv->deriv2_classes[2][2][118];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+534, dvrr_stack+1540, NULL);
 tmp = dvrr_stack + 534;
 target_ptr = Libderiv->deriv2_classes[1][2][117];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+2750, dvrr_stack+336, NULL);
 tmp = dvrr_stack + 2750;
 target_ptr = Libderiv->deriv2_classes[2][2][117];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+2786, dvrr_stack+171, dvrr_stack+1144);
 tmp = dvrr_stack + 2786;
 target_ptr = Libderiv->deriv2_classes[1][2][107];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+2804, dvrr_stack+1243, dvrr_stack+1225);
 tmp = dvrr_stack + 2804;
 target_ptr = Libderiv->deriv2_classes[2][2][107];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+2840, dvrr_stack+1594, dvrr_stack+1303);
 tmp = dvrr_stack + 2840;
 target_ptr = Libderiv->deriv2_classes[1][2][106];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+2858, dvrr_stack+1642, dvrr_stack+1624);
 tmp = dvrr_stack + 2858;
 target_ptr = Libderiv->deriv2_classes[2][2][106];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+2894, dvrr_stack+1702, dvrr_stack+1312);
 tmp = dvrr_stack + 2894;
 target_ptr = Libderiv->deriv2_classes[1][2][105];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+2912, dvrr_stack+264, dvrr_stack+246);
 tmp = dvrr_stack + 2912;
 target_ptr = Libderiv->deriv2_classes[2][2][105];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+2948, dvrr_stack+937, dvrr_stack+928);
 tmp = dvrr_stack + 2948;
 target_ptr = Libderiv->deriv2_classes[1][2][104];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+39, dvrr_stack+985, dvrr_stack+967);
 tmp = dvrr_stack + 39;
 target_ptr = Libderiv->deriv2_classes[2][2][104];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+2966, dvrr_stack+171, dvrr_stack+1144);
 tmp = dvrr_stack + 2966;
 target_ptr = Libderiv->deriv2_classes[1][2][95];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+735, dvrr_stack+1243, dvrr_stack+1225);
 tmp = dvrr_stack + 735;
 target_ptr = Libderiv->deriv2_classes[2][2][95];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+75, dvrr_stack+1594, dvrr_stack+1303);
 tmp = dvrr_stack + 75;
 target_ptr = Libderiv->deriv2_classes[1][2][94];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+771, dvrr_stack+1642, dvrr_stack+1624);
 tmp = dvrr_stack + 771;
 target_ptr = Libderiv->deriv2_classes[2][2][94];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+807, dvrr_stack+1702, dvrr_stack+1312);
 tmp = dvrr_stack + 807;
 target_ptr = Libderiv->deriv2_classes[1][2][93];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+1321, dvrr_stack+264, dvrr_stack+246);
 tmp = dvrr_stack + 1321;
 target_ptr = Libderiv->deriv2_classes[2][2][93];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+904, dvrr_stack+937, dvrr_stack+928);
 tmp = dvrr_stack + 904;
 target_ptr = Libderiv->deriv2_classes[1][2][92];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+1357, dvrr_stack+985, dvrr_stack+967);
 tmp = dvrr_stack + 1357;
 target_ptr = Libderiv->deriv2_classes[2][2][92];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+152, dvrr_stack+1054, dvrr_stack+1045);
 tmp = dvrr_stack + 152;
 target_ptr = Libderiv->deriv2_classes[1][2][91];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+1393, dvrr_stack+1732, dvrr_stack+1084);
 tmp = dvrr_stack + 1393;
 target_ptr = Libderiv->deriv2_classes[2][2][91];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+2110, dvrr_stack+171, dvrr_stack+1144);
 tmp = dvrr_stack + 2110;
 target_ptr = Libderiv->deriv2_classes[1][2][83];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1429, dvrr_stack+1243, dvrr_stack+1225);
 tmp = dvrr_stack + 1429;
 target_ptr = Libderiv->deriv2_classes[2][2][83];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+170, dvrr_stack+1594, dvrr_stack+1303);
 tmp = dvrr_stack + 170;
 target_ptr = Libderiv->deriv2_classes[1][2][82];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1465, dvrr_stack+1642, dvrr_stack+1624);
 tmp = dvrr_stack + 1465;
 target_ptr = Libderiv->deriv2_classes[2][2][82];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+219, dvrr_stack+1702, dvrr_stack+1312);
 tmp = dvrr_stack + 219;
 target_ptr = Libderiv->deriv2_classes[1][2][81];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1501, dvrr_stack+264, dvrr_stack+246);
 tmp = dvrr_stack + 1501;
 target_ptr = Libderiv->deriv2_classes[2][2][81];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+237, dvrr_stack+937, dvrr_stack+928);
 tmp = dvrr_stack + 237;
 target_ptr = Libderiv->deriv2_classes[1][2][80];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+922, dvrr_stack+985, dvrr_stack+967);
 tmp = dvrr_stack + 922;
 target_ptr = Libderiv->deriv2_classes[2][2][80];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+958, dvrr_stack+1054, dvrr_stack+1045);
 tmp = dvrr_stack + 958;
 target_ptr = Libderiv->deriv2_classes[1][2][79];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+976, dvrr_stack+1732, dvrr_stack+1084);
 tmp = dvrr_stack + 976;
 target_ptr = Libderiv->deriv2_classes[2][2][79];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+1012, dvrr_stack+1180, dvrr_stack+1171);
 tmp = dvrr_stack + 1012;
 target_ptr = Libderiv->deriv2_classes[1][2][78];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1171, dvrr_stack+1792, dvrr_stack+201);
 tmp = dvrr_stack + 1171;
 target_ptr = Libderiv->deriv2_classes[2][2][78];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,6,dvrr_stack+1207, dvrr_stack+850, dvrr_stack+1102);
 tmp = dvrr_stack + 1207;
 target_ptr = Libderiv->deriv2_classes[1][2][35];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1225, dvrr_stack+1948, dvrr_stack+886);
 tmp = dvrr_stack + 1225;
 target_ptr = Libderiv->deriv2_classes[2][2][35];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+1261, dvrr_stack+2318, dvrr_stack+6);
 tmp = dvrr_stack + 1261;
 target_ptr = Libderiv->deriv2_classes[1][2][34];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1279, dvrr_stack+2354, dvrr_stack+561);
 tmp = dvrr_stack + 1279;
 target_ptr = Libderiv->deriv2_classes[2][2][34];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+1030, dvrr_stack+2414, dvrr_stack+579);
 tmp = dvrr_stack + 1030;
 target_ptr = Libderiv->deriv2_classes[1][2][33];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1048, dvrr_stack+651, dvrr_stack+832);
 tmp = dvrr_stack + 1048;
 target_ptr = Libderiv->deriv2_classes[2][2][33];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+1084, dvrr_stack+2138, dvrr_stack+585);
 tmp = dvrr_stack + 1084;
 target_ptr = Libderiv->deriv2_classes[1][2][32];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+255, dvrr_stack+2192, dvrr_stack+2174);
 tmp = dvrr_stack + 255;
 target_ptr = Libderiv->deriv2_classes[2][2][32];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+188, dvrr_stack+2258, dvrr_stack+2252);
 tmp = dvrr_stack + 188;
 target_ptr = Libderiv->deriv2_classes[1][2][31];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+291, dvrr_stack+2450, dvrr_stack+2294);
 tmp = dvrr_stack + 291;
 target_ptr = Libderiv->deriv2_classes[2][2][31];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+327, dvrr_stack+98, dvrr_stack+2312);
 tmp = dvrr_stack + 327;
 target_ptr = Libderiv->deriv2_classes[1][2][30];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+345, dvrr_stack+591, dvrr_stack+1153);
 tmp = dvrr_stack + 345;
 target_ptr = Libderiv->deriv2_classes[2][2][30];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+381, dvrr_stack+2014, dvrr_stack+2008);
 tmp = dvrr_stack + 381;
 target_ptr = Libderiv->deriv2_classes[1][2][26];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+399, dvrr_stack+2510, dvrr_stack+134);
 tmp = dvrr_stack + 399;
 target_ptr = Libderiv->deriv2_classes[2][2][26];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,6,dvrr_stack+1537, dvrr_stack+850, dvrr_stack+1102);
 tmp = dvrr_stack + 1537;
 target_ptr = Libderiv->deriv2_classes[1][2][23];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,6,dvrr_stack+1555, dvrr_stack+1948, dvrr_stack+886);
 tmp = dvrr_stack + 1555;
 target_ptr = Libderiv->deriv2_classes[2][2][23];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1591, dvrr_stack+2318, dvrr_stack+6);
 tmp = dvrr_stack + 1591;
 target_ptr = Libderiv->deriv2_classes[1][2][22];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+1609, dvrr_stack+2354, dvrr_stack+561);
 tmp = dvrr_stack + 1609;
 target_ptr = Libderiv->deriv2_classes[2][2][22];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1645, dvrr_stack+2414, dvrr_stack+579);
 tmp = dvrr_stack + 1645;
 target_ptr = Libderiv->deriv2_classes[1][2][21];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+1663, dvrr_stack+651, dvrr_stack+832);
 tmp = dvrr_stack + 1663;
 target_ptr = Libderiv->deriv2_classes[2][2][21];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1699, dvrr_stack+2138, dvrr_stack+585);
 tmp = dvrr_stack + 1699;
 target_ptr = Libderiv->deriv2_classes[1][2][20];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+1717, dvrr_stack+2192, dvrr_stack+2174);
 tmp = dvrr_stack + 1717;
 target_ptr = Libderiv->deriv2_classes[2][2][20];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1753, dvrr_stack+2258, dvrr_stack+2252);
 tmp = dvrr_stack + 1753;
 target_ptr = Libderiv->deriv2_classes[1][2][19];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+1771, dvrr_stack+2450, dvrr_stack+2294);
 tmp = dvrr_stack + 1771;
 target_ptr = Libderiv->deriv2_classes[2][2][19];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1807, dvrr_stack+98, dvrr_stack+2312);
 tmp = dvrr_stack + 1807;
 target_ptr = Libderiv->deriv2_classes[1][2][18];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+2984, dvrr_stack+591, dvrr_stack+1153);
 tmp = dvrr_stack + 2984;
 target_ptr = Libderiv->deriv2_classes[2][2][18];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+3020, dvrr_stack+2014, dvrr_stack+2008);
 tmp = dvrr_stack + 3020;
 target_ptr = Libderiv->deriv2_classes[1][2][14];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3038, dvrr_stack+2510, dvrr_stack+134);
 tmp = dvrr_stack + 3038;
 target_ptr = Libderiv->deriv2_classes[2][2][14];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+3074, dvrr_stack+1852, dvrr_stack+0);
 tmp = dvrr_stack + 3074;
 target_ptr = Libderiv->deriv2_classes[1][2][13];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3092, dvrr_stack+2050, dvrr_stack+711);
 tmp = dvrr_stack + 3092;
 target_ptr = Libderiv->deriv2_classes[2][2][13];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,6,dvrr_stack+3128, dvrr_stack+850, dvrr_stack+1102);
 tmp = dvrr_stack + 3128;
 target_ptr = Libderiv->deriv2_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,6,dvrr_stack+850, dvrr_stack+1948, dvrr_stack+886);
 tmp = dvrr_stack + 850;
 target_ptr = Libderiv->deriv2_classes[2][2][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+886, dvrr_stack+2318, dvrr_stack+6);
 tmp = dvrr_stack + 886;
 target_ptr = Libderiv->deriv2_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2318, dvrr_stack+2354, dvrr_stack+561);
 tmp = dvrr_stack + 2318;
 target_ptr = Libderiv->deriv2_classes[2][2][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2354, dvrr_stack+2414, dvrr_stack+579);
 tmp = dvrr_stack + 2354;
 target_ptr = Libderiv->deriv2_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2372, dvrr_stack+651, dvrr_stack+832);
 tmp = dvrr_stack + 2372;
 target_ptr = Libderiv->deriv2_classes[2][2][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+651, dvrr_stack+2138, dvrr_stack+585);
 tmp = dvrr_stack + 651;
 target_ptr = Libderiv->deriv2_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+669, dvrr_stack+2192, dvrr_stack+2174);
 tmp = dvrr_stack + 669;
 target_ptr = Libderiv->deriv2_classes[2][2][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2408, dvrr_stack+2258, dvrr_stack+2252);
 tmp = dvrr_stack + 2408;
 target_ptr = Libderiv->deriv2_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1948, dvrr_stack+2450, dvrr_stack+2294);
 tmp = dvrr_stack + 1948;
 target_ptr = Libderiv->deriv2_classes[2][2][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2426, dvrr_stack+98, dvrr_stack+2312);
 tmp = dvrr_stack + 2426;
 target_ptr = Libderiv->deriv2_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2444, dvrr_stack+591, dvrr_stack+1153);
 tmp = dvrr_stack + 2444;
 target_ptr = Libderiv->deriv2_classes[2][2][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2480, dvrr_stack+2014, dvrr_stack+2008);
 tmp = dvrr_stack + 2480;
 target_ptr = Libderiv->deriv2_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1984, dvrr_stack+2510, dvrr_stack+134);
 tmp = dvrr_stack + 1984;
 target_ptr = Libderiv->deriv2_classes[2][2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2498, dvrr_stack+1852, dvrr_stack+0);
 tmp = dvrr_stack + 2498;
 target_ptr = Libderiv->deriv2_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2516, dvrr_stack+2050, dvrr_stack+711);
 tmp = dvrr_stack + 2516;
 target_ptr = Libderiv->deriv2_classes[2][2][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2552, dvrr_stack+1108, dvrr_stack+729);
 tmp = dvrr_stack + 2552;
 target_ptr = Libderiv->deriv2_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+2020, dvrr_stack+1888, dvrr_stack+21);
 tmp = dvrr_stack + 2020;
 target_ptr = Libderiv->deriv2_classes[2][2][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];


}

