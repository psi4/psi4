#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (p0|pp) integrals */

void d12vrr_order_p0pp(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+15, dvrr_stack+0, dvrr_stack+12, NULL, NULL, Data->F+2);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+24, dvrr_stack+0, dvrr_stack+12, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+30, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+33, dvrr_stack+12, dvrr_stack+30, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+39, dvrr_stack+24, dvrr_stack+33, NULL, NULL, dvrr_stack+12);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+57, dvrr_stack+6, dvrr_stack+24, NULL, NULL, dvrr_stack+0);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+75, dvrr_stack+57, dvrr_stack+39, dvrr_stack+6, dvrr_stack+24, dvrr_stack+15);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+111, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+1);
 tmp = dvrr_stack + 111;
 target_ptr = Libderiv->dvrr_classes[1][1];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+120,dvrr_stack+57,dvrr_stack+111,3);


 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+147, dvrr_stack+24, dvrr_stack+33, dvrr_stack+0, dvrr_stack+12, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+157, dvrr_stack+6, dvrr_stack+24, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+167, dvrr_stack+157, dvrr_stack+147, NULL, NULL, dvrr_stack+24);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+197,dvrr_stack+167,dvrr_stack+57,3);


 /* compute (1 0 | 1 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pd(Libderiv->CD,dvrr_stack+251,dvrr_stack+197,dvrr_stack+120,3);


 /* compute (1 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,9,dvrr_stack+305, dvrr_stack+251, dvrr_stack+111);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+332, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+335, dvrr_stack+30, dvrr_stack+332, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+341, dvrr_stack+33, dvrr_stack+335, dvrr_stack+12, dvrr_stack+30, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+351, dvrr_stack+147, dvrr_stack+341, dvrr_stack+24, dvrr_stack+33, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+366, dvrr_stack+157, dvrr_stack+147, dvrr_stack+6, dvrr_stack+24, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+381, dvrr_stack+366, dvrr_stack+351, NULL, NULL, dvrr_stack+147);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+426,dvrr_stack+381,dvrr_stack+167,3);


 /* compute (1 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+516,dvrr_stack+426,dvrr_stack+197,3);


 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,18,dvrr_stack+624, dvrr_stack+516, dvrr_stack+57);

 /* compute (1 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,9,dvrr_stack+351, dvrr_stack+251, dvrr_stack+111);

 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,18,dvrr_stack+678, dvrr_stack+516, dvrr_stack+57);

 /* compute (1 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,9,dvrr_stack+732, dvrr_stack+251, dvrr_stack+111);

 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,18,dvrr_stack+251, dvrr_stack+516, dvrr_stack+57);

 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+378, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 0 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0p(Libderiv->CD,dvrr_stack+516,dvrr_stack+111,dvrr_stack+378,3);


 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+332, dvrr_stack+516, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+525, dvrr_stack+197, NULL);
 tmp = dvrr_stack + 525;
 target_ptr = Libderiv->deriv_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+543, dvrr_stack+120, NULL);
 tmp = dvrr_stack + 543;
 target_ptr = Libderiv->deriv_classes[1][1][11];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+552, dvrr_stack+426, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+582, dvrr_stack+516, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+585, dvrr_stack+197, NULL);
 tmp = dvrr_stack + 585;
 target_ptr = Libderiv->deriv_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+603, dvrr_stack+120, NULL);
 tmp = dvrr_stack + 603;
 target_ptr = Libderiv->deriv_classes[1][1][10];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+759, dvrr_stack+426, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+612, dvrr_stack+516, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+789, dvrr_stack+197, NULL);
 tmp = dvrr_stack + 789;
 target_ptr = Libderiv->deriv_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+516, dvrr_stack+120, NULL);
 tmp = dvrr_stack + 516;
 target_ptr = Libderiv->deriv_classes[1][1][9];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+197, dvrr_stack+426, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_0(Data,3,1,dvrr_stack+426, dvrr_stack+111, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+429, dvrr_stack+167, dvrr_stack+111);
 tmp = dvrr_stack + 429;
 target_ptr = Libderiv->deriv_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+615, dvrr_stack+57, dvrr_stack+378);
 tmp = dvrr_stack + 615;
 target_ptr = Libderiv->deriv_classes[1][1][8];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+447, dvrr_stack+381, dvrr_stack+57);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_0(Data,3,1,dvrr_stack+477, dvrr_stack+111, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+480, dvrr_stack+167, dvrr_stack+111);
 tmp = dvrr_stack + 480;
 target_ptr = Libderiv->deriv_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+498, dvrr_stack+57, dvrr_stack+378);
 tmp = dvrr_stack + 498;
 target_ptr = Libderiv->deriv_classes[1][1][7];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+807, dvrr_stack+381, dvrr_stack+57);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_0(Data,3,1,dvrr_stack+507, dvrr_stack+111, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+120, dvrr_stack+167, dvrr_stack+111);
 tmp = dvrr_stack + 120;
 target_ptr = Libderiv->deriv_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+138, dvrr_stack+57, dvrr_stack+378);
 tmp = dvrr_stack + 138;
 target_ptr = Libderiv->deriv_classes[1][1][6];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+837, dvrr_stack+381, dvrr_stack+57);

 /* compute (0 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+381,dvrr_stack+6,dvrr_stack+3,1);


 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+390, dvrr_stack+381, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+393, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+396, dvrr_stack+111, dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, dvrr_stack+393);

 /* compute (2 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+867,dvrr_stack+75,dvrr_stack+396,6);


 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+227, dvrr_stack+867, NULL);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+921,dvrr_stack+157,dvrr_stack+6,1);


 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+510, dvrr_stack+921, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+939, dvrr_stack+147, dvrr_stack+341, NULL, NULL, dvrr_stack+33);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+969, dvrr_stack+167, dvrr_stack+939, dvrr_stack+157, dvrr_stack+147, dvrr_stack+39);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+1029,dvrr_stack+969,dvrr_stack+75,6);


 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+1137, dvrr_stack+1029, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+147, dvrr_stack+381, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+939, dvrr_stack+867, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+245, dvrr_stack+921, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+1173, dvrr_stack+1029, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+150, dvrr_stack+381, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+167, dvrr_stack+867, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+867, dvrr_stack+921, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+873, dvrr_stack+1029, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+1029, dvrr_stack+6, Data->F+0);

 /* compute (2 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+1032, dvrr_stack+378, dvrr_stack+393, Data->F+0, Data->F+1, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+1038, dvrr_stack+75, dvrr_stack+1032);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+378, dvrr_stack+157, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+1056, dvrr_stack+969, dvrr_stack+396);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+384, dvrr_stack+6, Data->F+0);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+1092, dvrr_stack+75, dvrr_stack+1032);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+1110, dvrr_stack+157, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+1209, dvrr_stack+969, dvrr_stack+396);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+387, dvrr_stack+6, Data->F+0);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+1116, dvrr_stack+75, dvrr_stack+1032);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+1032, dvrr_stack+157, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1245, dvrr_stack+969, dvrr_stack+396);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+1134, dvrr_stack+111, NULL);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+909, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+912, dvrr_stack+393, dvrr_stack+909, Data->F+1, Data->F+2, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+918, dvrr_stack+12, dvrr_stack+30, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+957, dvrr_stack+15, dvrr_stack+918, dvrr_stack+0, dvrr_stack+12, dvrr_stack+909);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+975, dvrr_stack+396, dvrr_stack+957, dvrr_stack+111, dvrr_stack+15, dvrr_stack+912);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+1005, dvrr_stack+975, dvrr_stack+111);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+1023, dvrr_stack+57, NULL);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1281, dvrr_stack+33, dvrr_stack+335, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1299, dvrr_stack+39, dvrr_stack+1281, dvrr_stack+24, dvrr_stack+33, dvrr_stack+918);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+1335, dvrr_stack+75, dvrr_stack+1299, dvrr_stack+57, dvrr_stack+39, dvrr_stack+957);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+1281, dvrr_stack+1335, dvrr_stack+57);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+0, dvrr_stack+111, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+957, dvrr_stack+975, dvrr_stack+111);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1317, dvrr_stack+57, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+12, dvrr_stack+1335, dvrr_stack+57);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+393, dvrr_stack+111, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+909, dvrr_stack+975, dvrr_stack+111);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+111, dvrr_stack+57, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1395, dvrr_stack+1335, dvrr_stack+57);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+975, dvrr_stack+75, dvrr_stack+6);
 tmp = dvrr_stack + 975;
 target_ptr = Libderiv->deriv_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1323, dvrr_stack+75, dvrr_stack+6);
 tmp = dvrr_stack + 1323;
 target_ptr = Libderiv->deriv_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1341, dvrr_stack+75, dvrr_stack+6);
 tmp = dvrr_stack + 1341;
 target_ptr = Libderiv->deriv_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,9,dvrr_stack+993, dvrr_stack+305, NULL);
 tmp = dvrr_stack + 993;
 target_ptr = Libderiv->deriv2_classes[1][1][143];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,18,dvrr_stack+1359, dvrr_stack+624, NULL);
 tmp = dvrr_stack + 1359;
 target_ptr = Libderiv->deriv2_classes[1][2][143];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,9,dvrr_stack+1377, dvrr_stack+305, NULL);
 tmp = dvrr_stack + 1377;
 target_ptr = Libderiv->deriv2_classes[1][1][131];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,18,dvrr_stack+48, dvrr_stack+624, NULL);
 tmp = dvrr_stack + 48;
 target_ptr = Libderiv->deriv2_classes[1][2][131];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+1386, dvrr_stack+351, NULL);
 tmp = dvrr_stack + 1386;
 target_ptr = Libderiv->deriv2_classes[1][1][130];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+66, dvrr_stack+678, NULL);
 tmp = dvrr_stack + 66;
 target_ptr = Libderiv->deriv2_classes[1][2][130];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,9,dvrr_stack+335, dvrr_stack+305, NULL);
 tmp = dvrr_stack + 335;
 target_ptr = Libderiv->deriv2_classes[1][1][119];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,18,dvrr_stack+305, dvrr_stack+624, NULL);
 tmp = dvrr_stack + 305;
 target_ptr = Libderiv->deriv2_classes[1][2][119];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+323, dvrr_stack+351, NULL);
 tmp = dvrr_stack + 323;
 target_ptr = Libderiv->deriv2_classes[1][1][118];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+624, dvrr_stack+678, NULL);
 tmp = dvrr_stack + 624;
 target_ptr = Libderiv->deriv2_classes[1][2][118];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+642, dvrr_stack+732, NULL);
 tmp = dvrr_stack + 642;
 target_ptr = Libderiv->deriv2_classes[1][1][117];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+651, dvrr_stack+251, NULL);
 tmp = dvrr_stack + 651;
 target_ptr = Libderiv->deriv2_classes[1][2][117];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+251, dvrr_stack+525, dvrr_stack+332);
 tmp = dvrr_stack + 251;
 target_ptr = Libderiv->deriv2_classes[1][1][107];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+260, dvrr_stack+552, dvrr_stack+543);
 tmp = dvrr_stack + 260;
 target_ptr = Libderiv->deriv2_classes[1][2][107];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+278, dvrr_stack+585, dvrr_stack+582);
 tmp = dvrr_stack + 278;
 target_ptr = Libderiv->deriv2_classes[1][1][106];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+287, dvrr_stack+759, dvrr_stack+603);
 tmp = dvrr_stack + 287;
 target_ptr = Libderiv->deriv2_classes[1][2][106];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+669, dvrr_stack+789, dvrr_stack+612);
 tmp = dvrr_stack + 669;
 target_ptr = Libderiv->deriv2_classes[1][1][105];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+678, dvrr_stack+197, dvrr_stack+516);
 tmp = dvrr_stack + 678;
 target_ptr = Libderiv->deriv2_classes[1][2][105];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+696, dvrr_stack+429, dvrr_stack+426);
 tmp = dvrr_stack + 696;
 target_ptr = Libderiv->deriv2_classes[1][1][104];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+705, dvrr_stack+447, dvrr_stack+615);
 tmp = dvrr_stack + 705;
 target_ptr = Libderiv->deriv2_classes[1][2][104];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+723, dvrr_stack+525, dvrr_stack+332);
 tmp = dvrr_stack + 723;
 target_ptr = Libderiv->deriv2_classes[1][1][95];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+732, dvrr_stack+552, dvrr_stack+543);
 tmp = dvrr_stack + 732;
 target_ptr = Libderiv->deriv2_classes[1][2][95];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+750, dvrr_stack+585, dvrr_stack+582);
 tmp = dvrr_stack + 750;
 target_ptr = Libderiv->deriv2_classes[1][1][94];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+344, dvrr_stack+759, dvrr_stack+603);
 tmp = dvrr_stack + 344;
 target_ptr = Libderiv->deriv2_classes[1][2][94];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+362, dvrr_stack+789, dvrr_stack+612);
 tmp = dvrr_stack + 362;
 target_ptr = Libderiv->deriv2_classes[1][1][93];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+84, dvrr_stack+197, dvrr_stack+516);
 tmp = dvrr_stack + 84;
 target_ptr = Libderiv->deriv2_classes[1][2][93];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+102, dvrr_stack+429, dvrr_stack+426);
 tmp = dvrr_stack + 102;
 target_ptr = Libderiv->deriv2_classes[1][1][92];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+1431, dvrr_stack+447, dvrr_stack+615);
 tmp = dvrr_stack + 1431;
 target_ptr = Libderiv->deriv2_classes[1][2][92];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+927, dvrr_stack+480, dvrr_stack+477);
 tmp = dvrr_stack + 927;
 target_ptr = Libderiv->deriv2_classes[1][1][91];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+1449, dvrr_stack+807, dvrr_stack+498);
 tmp = dvrr_stack + 1449;
 target_ptr = Libderiv->deriv2_classes[1][2][91];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+153, dvrr_stack+525, dvrr_stack+332);
 tmp = dvrr_stack + 153;
 target_ptr = Libderiv->deriv2_classes[1][1][83];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+525, dvrr_stack+552, dvrr_stack+543);
 tmp = dvrr_stack + 525;
 target_ptr = Libderiv->deriv2_classes[1][2][83];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+543, dvrr_stack+585, dvrr_stack+582);
 tmp = dvrr_stack + 543;
 target_ptr = Libderiv->deriv2_classes[1][1][82];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+552, dvrr_stack+759, dvrr_stack+603);
 tmp = dvrr_stack + 552;
 target_ptr = Libderiv->deriv2_classes[1][2][82];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+759, dvrr_stack+789, dvrr_stack+612);
 tmp = dvrr_stack + 759;
 target_ptr = Libderiv->deriv2_classes[1][1][81];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+768, dvrr_stack+197, dvrr_stack+516);
 tmp = dvrr_stack + 768;
 target_ptr = Libderiv->deriv2_classes[1][2][81];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+516, dvrr_stack+429, dvrr_stack+426);
 tmp = dvrr_stack + 516;
 target_ptr = Libderiv->deriv2_classes[1][1][80];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+786, dvrr_stack+447, dvrr_stack+615);
 tmp = dvrr_stack + 786;
 target_ptr = Libderiv->deriv2_classes[1][2][80];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+570, dvrr_stack+480, dvrr_stack+477);
 tmp = dvrr_stack + 570;
 target_ptr = Libderiv->deriv2_classes[1][1][79];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+579, dvrr_stack+807, dvrr_stack+498);
 tmp = dvrr_stack + 579;
 target_ptr = Libderiv->deriv2_classes[1][2][79];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+804, dvrr_stack+120, dvrr_stack+507);
 tmp = dvrr_stack + 804;
 target_ptr = Libderiv->deriv2_classes[1][1][78];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+813, dvrr_stack+837, dvrr_stack+138);
 tmp = dvrr_stack + 813;
 target_ptr = Libderiv->deriv2_classes[1][2][78];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,3,dvrr_stack+831, dvrr_stack+227, dvrr_stack+390);
 tmp = dvrr_stack + 831;
 target_ptr = Libderiv->deriv2_classes[1][1][35];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,6,dvrr_stack+840, dvrr_stack+1137, dvrr_stack+510);
 tmp = dvrr_stack + 840;
 target_ptr = Libderiv->deriv2_classes[1][2][35];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+858, dvrr_stack+939, dvrr_stack+147);
 tmp = dvrr_stack + 858;
 target_ptr = Libderiv->deriv2_classes[1][1][34];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+597, dvrr_stack+1173, dvrr_stack+245);
 tmp = dvrr_stack + 597;
 target_ptr = Libderiv->deriv2_classes[1][2][34];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+615, dvrr_stack+167, dvrr_stack+150);
 tmp = dvrr_stack + 615;
 target_ptr = Libderiv->deriv2_classes[1][1][33];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+117, dvrr_stack+873, dvrr_stack+867);
 tmp = dvrr_stack + 117;
 target_ptr = Libderiv->deriv2_classes[1][2][33];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+135, dvrr_stack+1038, dvrr_stack+1029);
 tmp = dvrr_stack + 135;
 target_ptr = Libderiv->deriv2_classes[1][1][32];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+185, dvrr_stack+1056, dvrr_stack+378);
 tmp = dvrr_stack + 185;
 target_ptr = Libderiv->deriv2_classes[1][2][32];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+203, dvrr_stack+1092, dvrr_stack+384);
 tmp = dvrr_stack + 203;
 target_ptr = Libderiv->deriv2_classes[1][1][31];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+414, dvrr_stack+1209, dvrr_stack+1110);
 tmp = dvrr_stack + 414;
 target_ptr = Libderiv->deriv2_classes[1][2][31];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+212, dvrr_stack+396, dvrr_stack+3);
 tmp = dvrr_stack + 212;
 target_ptr = Libderiv->deriv_classes[1][1][2];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+432, dvrr_stack+1116, dvrr_stack+387);
 tmp = dvrr_stack + 432;
 target_ptr = Libderiv->deriv2_classes[1][1][30];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+441, dvrr_stack+1245, dvrr_stack+1032);
 tmp = dvrr_stack + 441;
 target_ptr = Libderiv->deriv2_classes[1][2][30];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+459, dvrr_stack+1005, dvrr_stack+1134);
 tmp = dvrr_stack + 459;
 target_ptr = Libderiv->deriv2_classes[1][1][26];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+468, dvrr_stack+1281, dvrr_stack+1023);
 tmp = dvrr_stack + 468;
 target_ptr = Libderiv->deriv2_classes[1][2][26];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,3,dvrr_stack+486, dvrr_stack+227, dvrr_stack+390);
 tmp = dvrr_stack + 486;
 target_ptr = Libderiv->deriv2_classes[1][1][23];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,6,dvrr_stack+1467, dvrr_stack+1137, dvrr_stack+510);
 tmp = dvrr_stack + 1467;
 target_ptr = Libderiv->deriv2_classes[1][2][23];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+495, dvrr_stack+939, dvrr_stack+147);
 tmp = dvrr_stack + 495;
 target_ptr = Libderiv->deriv2_classes[1][1][22];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1485, dvrr_stack+1173, dvrr_stack+245);
 tmp = dvrr_stack + 1485;
 target_ptr = Libderiv->deriv2_classes[1][2][22];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+1503, dvrr_stack+167, dvrr_stack+150);
 tmp = dvrr_stack + 1503;
 target_ptr = Libderiv->deriv2_classes[1][1][21];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1512, dvrr_stack+873, dvrr_stack+867);
 tmp = dvrr_stack + 1512;
 target_ptr = Libderiv->deriv2_classes[1][2][21];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+1530, dvrr_stack+1038, dvrr_stack+1029);
 tmp = dvrr_stack + 1530;
 target_ptr = Libderiv->deriv2_classes[1][1][20];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1539, dvrr_stack+1056, dvrr_stack+378);
 tmp = dvrr_stack + 1539;
 target_ptr = Libderiv->deriv2_classes[1][2][20];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+1557, dvrr_stack+1092, dvrr_stack+384);
 tmp = dvrr_stack + 1557;
 target_ptr = Libderiv->deriv2_classes[1][1][19];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1566, dvrr_stack+1209, dvrr_stack+1110);
 tmp = dvrr_stack + 1566;
 target_ptr = Libderiv->deriv2_classes[1][2][19];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+1584, dvrr_stack+396, dvrr_stack+3);
 tmp = dvrr_stack + 1584;
 target_ptr = Libderiv->deriv_classes[1][1][1];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+1593, dvrr_stack+1116, dvrr_stack+387);
 tmp = dvrr_stack + 1593;
 target_ptr = Libderiv->deriv2_classes[1][1][18];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1602, dvrr_stack+1245, dvrr_stack+1032);
 tmp = dvrr_stack + 1602;
 target_ptr = Libderiv->deriv2_classes[1][2][18];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+1620, dvrr_stack+1005, dvrr_stack+1134);
 tmp = dvrr_stack + 1620;
 target_ptr = Libderiv->deriv2_classes[1][1][14];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1629, dvrr_stack+1281, dvrr_stack+1023);
 tmp = dvrr_stack + 1629;
 target_ptr = Libderiv->deriv2_classes[1][2][14];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+1647, dvrr_stack+957, dvrr_stack+0);
 tmp = dvrr_stack + 1647;
 target_ptr = Libderiv->deriv2_classes[1][1][13];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1656, dvrr_stack+12, dvrr_stack+1317);
 tmp = dvrr_stack + 1656;
 target_ptr = Libderiv->deriv2_classes[1][2][13];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,3,dvrr_stack+1674, dvrr_stack+227, dvrr_stack+390);
 tmp = dvrr_stack + 1674;
 target_ptr = Libderiv->deriv2_classes[1][1][11];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,6,dvrr_stack+221, dvrr_stack+1137, dvrr_stack+510);
 tmp = dvrr_stack + 221;
 target_ptr = Libderiv->deriv2_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+1137, dvrr_stack+939, dvrr_stack+147);
 tmp = dvrr_stack + 1137;
 target_ptr = Libderiv->deriv2_classes[1][1][10];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1146, dvrr_stack+1173, dvrr_stack+245);
 tmp = dvrr_stack + 1146;
 target_ptr = Libderiv->deriv2_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+1164, dvrr_stack+167, dvrr_stack+150);
 tmp = dvrr_stack + 1164;
 target_ptr = Libderiv->deriv2_classes[1][1][9];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1173, dvrr_stack+873, dvrr_stack+867);
 tmp = dvrr_stack + 1173;
 target_ptr = Libderiv->deriv2_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+144, dvrr_stack+1038, dvrr_stack+1029);
 tmp = dvrr_stack + 144;
 target_ptr = Libderiv->deriv2_classes[1][1][8];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1038, dvrr_stack+1056, dvrr_stack+378);
 tmp = dvrr_stack + 1038;
 target_ptr = Libderiv->deriv2_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+1056, dvrr_stack+1092, dvrr_stack+384);
 tmp = dvrr_stack + 1056;
 target_ptr = Libderiv->deriv2_classes[1][1][7];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1191, dvrr_stack+1209, dvrr_stack+1110);
 tmp = dvrr_stack + 1191;
 target_ptr = Libderiv->deriv2_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+1209, dvrr_stack+396, dvrr_stack+3);
 tmp = dvrr_stack + 1209;
 target_ptr = Libderiv->deriv_classes[1][1][0];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+3, dvrr_stack+1116, dvrr_stack+387);
 tmp = dvrr_stack + 3;
 target_ptr = Libderiv->deriv2_classes[1][1][6];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+396, dvrr_stack+1245, dvrr_stack+1032);
 tmp = dvrr_stack + 396;
 target_ptr = Libderiv->deriv2_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+1029, dvrr_stack+1005, dvrr_stack+1134);
 tmp = dvrr_stack + 1029;
 target_ptr = Libderiv->deriv2_classes[1][1][2];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1218, dvrr_stack+1281, dvrr_stack+1023);
 tmp = dvrr_stack + 1218;
 target_ptr = Libderiv->deriv2_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+1236, dvrr_stack+957, dvrr_stack+0);
 tmp = dvrr_stack + 1236;
 target_ptr = Libderiv->deriv2_classes[1][1][1];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1245, dvrr_stack+12, dvrr_stack+1317);
 tmp = dvrr_stack + 1245;
 target_ptr = Libderiv->deriv2_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+12, dvrr_stack+909, dvrr_stack+393);
 tmp = dvrr_stack + 12;
 target_ptr = Libderiv->deriv2_classes[1][1][0];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+21, dvrr_stack+1395, dvrr_stack+111);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv2_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];


}

