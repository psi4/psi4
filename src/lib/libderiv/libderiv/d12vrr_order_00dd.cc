#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|dd) integrals */

void d12vrr_order_00dd(Libderiv_t *Libderiv, prim_data *Data)
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
 tmp = dvrr_stack + 68;
 target_ptr = Libderiv->dvrr_classes[0][2];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+74, dvrr_stack+68, dvrr_stack+15, dvrr_stack+65, dvrr_stack+3, NULL);
 tmp = dvrr_stack + 74;
 target_ptr = Libderiv->dvrr_classes[0][3];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+84, dvrr_stack+74, dvrr_stack+21, dvrr_stack+68, dvrr_stack+15, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+99, dvrr_stack+84, dvrr_stack+50, NULL, NULL, dvrr_stack+21);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+144,dvrr_stack+74,dvrr_stack+68,1);


 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+162,dvrr_stack+84,dvrr_stack+74,1);


 /* compute (0 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+192,dvrr_stack+162,dvrr_stack+144,1);


 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,6,dvrr_stack+228, dvrr_stack+192, dvrr_stack+68);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+246, dvrr_stack+84, dvrr_stack+50, dvrr_stack+74, dvrr_stack+21, NULL);

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+267,dvrr_stack+246,dvrr_stack+84,1);


 /* compute (0 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+312,dvrr_stack+267,dvrr_stack+162,1);


 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,10,dvrr_stack+372, dvrr_stack+312, dvrr_stack+74);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+402, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+405, dvrr_stack+31, dvrr_stack+402, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+411, dvrr_stack+34, dvrr_stack+405, dvrr_stack+6, dvrr_stack+31, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+421, dvrr_stack+40, dvrr_stack+411, dvrr_stack+9, dvrr_stack+34, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+436, dvrr_stack+50, dvrr_stack+421, dvrr_stack+21, dvrr_stack+40, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+457, dvrr_stack+246, dvrr_stack+436, dvrr_stack+84, dvrr_stack+50, NULL);

 /* compute (0 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+485,dvrr_stack+457,dvrr_stack+246,1);


 /* compute (0 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+548,dvrr_stack+485,dvrr_stack+267,1);


 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,15,dvrr_stack+638, dvrr_stack+548, dvrr_stack+84);

 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,6,dvrr_stack+402, dvrr_stack+192, dvrr_stack+68);

 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,10,dvrr_stack+683, dvrr_stack+312, dvrr_stack+74);

 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,15,dvrr_stack+713, dvrr_stack+548, dvrr_stack+84);

 /* compute (0 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,6,dvrr_stack+758, dvrr_stack+192, dvrr_stack+68);

 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,10,dvrr_stack+192, dvrr_stack+312, dvrr_stack+74);

 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,15,dvrr_stack+312, dvrr_stack+548, dvrr_stack+84);

 /* compute (0 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+31,dvrr_stack+68,dvrr_stack+65,1);


 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+6, dvrr_stack+31, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+548, dvrr_stack+162, NULL);
 tmp = dvrr_stack + 548;
 target_ptr = Libderiv->deriv_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+222, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 222;
 target_ptr = Libderiv->deriv_classes[0][2][11];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+357, dvrr_stack+267, NULL);
 tmp = dvrr_stack + 357;
 target_ptr = Libderiv->deriv_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,21,dvrr_stack+558, dvrr_stack+485, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+579, dvrr_stack+31, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+582, dvrr_stack+162, NULL);
 tmp = dvrr_stack + 582;
 target_ptr = Libderiv->deriv_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+592, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 592;
 target_ptr = Libderiv->deriv_classes[0][2][10];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+598, dvrr_stack+267, NULL);
 tmp = dvrr_stack + 598;
 target_ptr = Libderiv->deriv_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+613, dvrr_stack+485, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+634, dvrr_stack+31, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+776, dvrr_stack+162, NULL);
 tmp = dvrr_stack + 776;
 target_ptr = Libderiv->deriv_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+162, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 162;
 target_ptr = Libderiv->deriv_classes[0][2][9];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+144, dvrr_stack+267, NULL);
 tmp = dvrr_stack + 144;
 target_ptr = Libderiv->deriv_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+267, dvrr_stack+485, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+159, dvrr_stack+68, Data->F+0);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+485, dvrr_stack+84, dvrr_stack+68);
 tmp = dvrr_stack + 485;
 target_ptr = Libderiv->deriv_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+495, dvrr_stack+74, dvrr_stack+65);
 tmp = dvrr_stack + 495;
 target_ptr = Libderiv->deriv_classes[0][2][8];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+501, dvrr_stack+246, dvrr_stack+74);
 tmp = dvrr_stack + 501;
 target_ptr = Libderiv->deriv_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+516, dvrr_stack+457, dvrr_stack+84);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+537, dvrr_stack+68, Data->F+0);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+288, dvrr_stack+84, dvrr_stack+68);
 tmp = dvrr_stack + 288;
 target_ptr = Libderiv->deriv_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+540, dvrr_stack+74, dvrr_stack+65);
 tmp = dvrr_stack + 540;
 target_ptr = Libderiv->deriv_classes[0][2][7];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+168, dvrr_stack+246, dvrr_stack+74);
 tmp = dvrr_stack + 168;
 target_ptr = Libderiv->deriv_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+786, dvrr_stack+457, dvrr_stack+84);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+298, dvrr_stack+68, Data->F+0);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+301, dvrr_stack+84, dvrr_stack+68);
 tmp = dvrr_stack + 301;
 target_ptr = Libderiv->deriv_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+183, dvrr_stack+74, dvrr_stack+65);
 tmp = dvrr_stack + 183;
 target_ptr = Libderiv->deriv_classes[0][2][6];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+807, dvrr_stack+246, dvrr_stack+74);
 tmp = dvrr_stack + 807;
 target_ptr = Libderiv->deriv_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+822, dvrr_stack+457, dvrr_stack+84);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+457, dvrr_stack+68, dvrr_stack+15, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+843, dvrr_stack+74, dvrr_stack+21, NULL, NULL, dvrr_stack+15);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+873,dvrr_stack+843,dvrr_stack+457,3);


 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+927, dvrr_stack+873, NULL);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+945,dvrr_stack+99,dvrr_stack+843,3);


 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+1035, dvrr_stack+945, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1065, dvrr_stack+246, dvrr_stack+436, NULL, NULL, dvrr_stack+50);

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+1128,dvrr_stack+1065,dvrr_stack+99,3);


 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+1263, dvrr_stack+1128, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+436, dvrr_stack+873, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+1308, dvrr_stack+945, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+1338, dvrr_stack+1128, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+246, dvrr_stack+873, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+873, dvrr_stack+945, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+945, dvrr_stack+1128, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+31, dvrr_stack+65, dvrr_stack+3, NULL, NULL, Data->F+1);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+1128, dvrr_stack+843, dvrr_stack+31);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+1146, dvrr_stack+99, dvrr_stack+457);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+990, dvrr_stack+1065, dvrr_stack+843);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+1176, dvrr_stack+843, dvrr_stack+31);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1194, dvrr_stack+99, dvrr_stack+457);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+1383, dvrr_stack+1065, dvrr_stack+843);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+1224, dvrr_stack+843, dvrr_stack+31);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+1428, dvrr_stack+99, dvrr_stack+457);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+1458, dvrr_stack+1065, dvrr_stack+843);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+31, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1065, dvrr_stack+15, dvrr_stack+9, NULL, NULL, dvrr_stack+0);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1083, dvrr_stack+457, dvrr_stack+1065, dvrr_stack+68, dvrr_stack+15, dvrr_stack+31);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+1242, dvrr_stack+1083, dvrr_stack+68);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1503, dvrr_stack+21, dvrr_stack+40, NULL, NULL, dvrr_stack+9);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1533, dvrr_stack+843, dvrr_stack+1503, dvrr_stack+74, dvrr_stack+21, dvrr_stack+1065);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+9, dvrr_stack+1533, dvrr_stack+74);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1593, dvrr_stack+50, dvrr_stack+421, NULL, NULL, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1638, dvrr_stack+99, dvrr_stack+1593, dvrr_stack+84, dvrr_stack+50, dvrr_stack+1503);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1593, dvrr_stack+1638, dvrr_stack+84);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+1065, dvrr_stack+1083, dvrr_stack+68);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+1503, dvrr_stack+1533, dvrr_stack+74);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+1728, dvrr_stack+1638, dvrr_stack+84);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+39, dvrr_stack+1083, dvrr_stack+68);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+1083, dvrr_stack+1533, dvrr_stack+74);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+1533, dvrr_stack+1638, dvrr_stack+84);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1578, dvrr_stack+99, NULL);
 tmp = dvrr_stack + 1578;
 target_ptr = Libderiv->deriv_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+1113, dvrr_stack+99, NULL);
 tmp = dvrr_stack + 1113;
 target_ptr = Libderiv->deriv_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1638, dvrr_stack+99, NULL);
 tmp = dvrr_stack + 1638;
 target_ptr = Libderiv->deriv_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,6,dvrr_stack+0, dvrr_stack+228, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[0][2][143];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,10,dvrr_stack+475, dvrr_stack+372, NULL);
 tmp = dvrr_stack + 475;
 target_ptr = Libderiv->deriv2_classes[0][3][143];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,15,dvrr_stack+1653, dvrr_stack+638, NULL);
 tmp = dvrr_stack + 1653;
 target_ptr = Libderiv->deriv2_classes[0][4][143];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,6,dvrr_stack+1668, dvrr_stack+228, NULL);
 tmp = dvrr_stack + 1668;
 target_ptr = Libderiv->deriv2_classes[0][2][131];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,10,dvrr_stack+1674, dvrr_stack+372, NULL);
 tmp = dvrr_stack + 1674;
 target_ptr = Libderiv->deriv2_classes[0][3][131];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,15,dvrr_stack+1684, dvrr_stack+638, NULL);
 tmp = dvrr_stack + 1684;
 target_ptr = Libderiv->deriv2_classes[0][4][131];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+1699, dvrr_stack+402, NULL);
 tmp = dvrr_stack + 1699;
 target_ptr = Libderiv->deriv2_classes[0][2][130];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+1705, dvrr_stack+683, NULL);
 tmp = dvrr_stack + 1705;
 target_ptr = Libderiv->deriv2_classes[0][3][130];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+57, dvrr_stack+713, NULL);
 tmp = dvrr_stack + 57;
 target_ptr = Libderiv->deriv2_classes[0][4][130];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,6,dvrr_stack+1715, dvrr_stack+228, NULL);
 tmp = dvrr_stack + 1715;
 target_ptr = Libderiv->deriv2_classes[0][2][119];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,10,dvrr_stack+228, dvrr_stack+372, NULL);
 tmp = dvrr_stack + 228;
 target_ptr = Libderiv->deriv2_classes[0][3][119];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,15,dvrr_stack+372, dvrr_stack+638, NULL);
 tmp = dvrr_stack + 372;
 target_ptr = Libderiv->deriv2_classes[0][4][119];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+387, dvrr_stack+402, NULL);
 tmp = dvrr_stack + 387;
 target_ptr = Libderiv->deriv2_classes[0][2][118];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+393, dvrr_stack+683, NULL);
 tmp = dvrr_stack + 393;
 target_ptr = Libderiv->deriv2_classes[0][3][118];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+403, dvrr_stack+713, NULL);
 tmp = dvrr_stack + 403;
 target_ptr = Libderiv->deriv2_classes[0][4][118];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+418, dvrr_stack+758, NULL);
 tmp = dvrr_stack + 418;
 target_ptr = Libderiv->deriv2_classes[0][2][117];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+424, dvrr_stack+192, NULL);
 tmp = dvrr_stack + 424;
 target_ptr = Libderiv->deriv2_classes[0][3][117];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+72, dvrr_stack+312, NULL);
 tmp = dvrr_stack + 72;
 target_ptr = Libderiv->deriv2_classes[0][4][117];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+238, dvrr_stack+548, dvrr_stack+6);
 tmp = dvrr_stack + 238;
 target_ptr = Libderiv->deriv2_classes[0][2][107];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+87, dvrr_stack+357, dvrr_stack+222);
 tmp = dvrr_stack + 87;
 target_ptr = Libderiv->deriv2_classes[0][3][107];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+97, dvrr_stack+558, dvrr_stack+548);
 tmp = dvrr_stack + 97;
 target_ptr = Libderiv->deriv2_classes[0][4][107];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+1721, dvrr_stack+582, dvrr_stack+579);
 tmp = dvrr_stack + 1721;
 target_ptr = Libderiv->deriv2_classes[0][2][106];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+112, dvrr_stack+598, dvrr_stack+592);
 tmp = dvrr_stack + 112;
 target_ptr = Libderiv->deriv2_classes[0][3][106];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+122, dvrr_stack+613, dvrr_stack+582);
 tmp = dvrr_stack + 122;
 target_ptr = Libderiv->deriv2_classes[0][4][106];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+137, dvrr_stack+776, dvrr_stack+634);
 tmp = dvrr_stack + 137;
 target_ptr = Libderiv->deriv2_classes[0][2][105];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+903, dvrr_stack+144, dvrr_stack+162);
 tmp = dvrr_stack + 903;
 target_ptr = Libderiv->deriv2_classes[0][3][105];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+311, dvrr_stack+267, dvrr_stack+776);
 tmp = dvrr_stack + 311;
 target_ptr = Libderiv->deriv2_classes[0][4][105];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+913, dvrr_stack+485, dvrr_stack+159);
 tmp = dvrr_stack + 913;
 target_ptr = Libderiv->deriv2_classes[0][2][104];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+326, dvrr_stack+501, dvrr_stack+495);
 tmp = dvrr_stack + 326;
 target_ptr = Libderiv->deriv2_classes[0][3][104];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+336, dvrr_stack+516, dvrr_stack+485);
 tmp = dvrr_stack + 336;
 target_ptr = Libderiv->deriv2_classes[0][4][104];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+351, dvrr_stack+548, dvrr_stack+6);
 tmp = dvrr_stack + 351;
 target_ptr = Libderiv->deriv2_classes[0][2][95];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+189, dvrr_stack+357, dvrr_stack+222);
 tmp = dvrr_stack + 189;
 target_ptr = Libderiv->deriv2_classes[0][3][95];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+199, dvrr_stack+558, dvrr_stack+548);
 tmp = dvrr_stack + 199;
 target_ptr = Libderiv->deriv2_classes[0][4][95];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+919, dvrr_stack+582, dvrr_stack+579);
 tmp = dvrr_stack + 919;
 target_ptr = Libderiv->deriv2_classes[0][2][94];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+637, dvrr_stack+598, dvrr_stack+592);
 tmp = dvrr_stack + 637;
 target_ptr = Libderiv->deriv2_classes[0][3][94];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+647, dvrr_stack+613, dvrr_stack+582);
 tmp = dvrr_stack + 647;
 target_ptr = Libderiv->deriv2_classes[0][4][94];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+214, dvrr_stack+776, dvrr_stack+634);
 tmp = dvrr_stack + 214;
 target_ptr = Libderiv->deriv2_classes[0][2][93];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+662, dvrr_stack+144, dvrr_stack+162);
 tmp = dvrr_stack + 662;
 target_ptr = Libderiv->deriv2_classes[0][3][93];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+672, dvrr_stack+267, dvrr_stack+776);
 tmp = dvrr_stack + 672;
 target_ptr = Libderiv->deriv2_classes[0][4][93];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+687, dvrr_stack+485, dvrr_stack+159);
 tmp = dvrr_stack + 687;
 target_ptr = Libderiv->deriv2_classes[0][2][92];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+693, dvrr_stack+501, dvrr_stack+495);
 tmp = dvrr_stack + 693;
 target_ptr = Libderiv->deriv2_classes[0][3][92];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+703, dvrr_stack+516, dvrr_stack+485);
 tmp = dvrr_stack + 703;
 target_ptr = Libderiv->deriv2_classes[0][4][92];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+718, dvrr_stack+288, dvrr_stack+537);
 tmp = dvrr_stack + 718;
 target_ptr = Libderiv->deriv2_classes[0][2][91];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+724, dvrr_stack+168, dvrr_stack+540);
 tmp = dvrr_stack + 724;
 target_ptr = Libderiv->deriv2_classes[0][3][91];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+734, dvrr_stack+786, dvrr_stack+288);
 tmp = dvrr_stack + 734;
 target_ptr = Libderiv->deriv2_classes[0][4][91];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+749, dvrr_stack+548, dvrr_stack+6);
 tmp = dvrr_stack + 749;
 target_ptr = Libderiv->deriv2_classes[0][2][83];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+755, dvrr_stack+357, dvrr_stack+222);
 tmp = dvrr_stack + 755;
 target_ptr = Libderiv->deriv2_classes[0][3][83];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+357, dvrr_stack+558, dvrr_stack+548);
 tmp = dvrr_stack + 357;
 target_ptr = Libderiv->deriv2_classes[0][4][83];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+546, dvrr_stack+582, dvrr_stack+579);
 tmp = dvrr_stack + 546;
 target_ptr = Libderiv->deriv2_classes[0][2][82];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+552, dvrr_stack+598, dvrr_stack+592);
 tmp = dvrr_stack + 552;
 target_ptr = Libderiv->deriv2_classes[0][3][82];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+592, dvrr_stack+613, dvrr_stack+582);
 tmp = dvrr_stack + 592;
 target_ptr = Libderiv->deriv2_classes[0][4][82];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+607, dvrr_stack+776, dvrr_stack+634);
 tmp = dvrr_stack + 607;
 target_ptr = Libderiv->deriv2_classes[0][2][81];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+613, dvrr_stack+144, dvrr_stack+162);
 tmp = dvrr_stack + 613;
 target_ptr = Libderiv->deriv2_classes[0][3][81];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+143, dvrr_stack+267, dvrr_stack+776);
 tmp = dvrr_stack + 143;
 target_ptr = Libderiv->deriv2_classes[0][4][81];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+162, dvrr_stack+485, dvrr_stack+159);
 tmp = dvrr_stack + 162;
 target_ptr = Libderiv->deriv2_classes[0][2][80];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+623, dvrr_stack+501, dvrr_stack+495);
 tmp = dvrr_stack + 623;
 target_ptr = Libderiv->deriv2_classes[0][3][80];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+495, dvrr_stack+516, dvrr_stack+485);
 tmp = dvrr_stack + 495;
 target_ptr = Libderiv->deriv2_classes[0][4][80];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+485, dvrr_stack+288, dvrr_stack+537);
 tmp = dvrr_stack + 485;
 target_ptr = Libderiv->deriv2_classes[0][2][79];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+510, dvrr_stack+168, dvrr_stack+540);
 tmp = dvrr_stack + 510;
 target_ptr = Libderiv->deriv2_classes[0][3][79];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+168, dvrr_stack+786, dvrr_stack+288);
 tmp = dvrr_stack + 168;
 target_ptr = Libderiv->deriv2_classes[0][4][79];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+520, dvrr_stack+301, dvrr_stack+298);
 tmp = dvrr_stack + 520;
 target_ptr = Libderiv->deriv2_classes[0][2][78];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+526, dvrr_stack+807, dvrr_stack+183);
 tmp = dvrr_stack + 526;
 target_ptr = Libderiv->deriv2_classes[0][3][78];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+264, dvrr_stack+822, dvrr_stack+301);
 tmp = dvrr_stack + 264;
 target_ptr = Libderiv->deriv2_classes[0][4][78];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,6,dvrr_stack+183, dvrr_stack+927, NULL);
 tmp = dvrr_stack + 183;
 target_ptr = Libderiv->deriv2_classes[0][2][35];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,10,dvrr_stack+536, dvrr_stack+1035, NULL);
 tmp = dvrr_stack + 536;
 target_ptr = Libderiv->deriv2_classes[0][3][35];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,15,dvrr_stack+279, dvrr_stack+1263, NULL);
 tmp = dvrr_stack + 279;
 target_ptr = Libderiv->deriv2_classes[0][4][35];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+294, dvrr_stack+436, NULL);
 tmp = dvrr_stack + 294;
 target_ptr = Libderiv->deriv2_classes[0][2][34];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+300, dvrr_stack+1308, NULL);
 tmp = dvrr_stack + 300;
 target_ptr = Libderiv->deriv2_classes[0][3][34];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+562, dvrr_stack+1338, NULL);
 tmp = dvrr_stack + 562;
 target_ptr = Libderiv->deriv2_classes[0][4][34];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+577, dvrr_stack+246, NULL);
 tmp = dvrr_stack + 577;
 target_ptr = Libderiv->deriv2_classes[0][2][33];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+765, dvrr_stack+873, NULL);
 tmp = dvrr_stack + 765;
 target_ptr = Libderiv->deriv2_classes[0][3][33];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+775, dvrr_stack+945, NULL);
 tmp = dvrr_stack + 775;
 target_ptr = Libderiv->deriv2_classes[0][4][33];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+583, dvrr_stack+1128, NULL);
 tmp = dvrr_stack + 583;
 target_ptr = Libderiv->deriv2_classes[0][2][32];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+790, dvrr_stack+1146, NULL);
 tmp = dvrr_stack + 790;
 target_ptr = Libderiv->deriv2_classes[0][3][32];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+800, dvrr_stack+990, NULL);
 tmp = dvrr_stack + 800;
 target_ptr = Libderiv->deriv2_classes[0][4][32];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+220, dvrr_stack+1176, NULL);
 tmp = dvrr_stack + 220;
 target_ptr = Libderiv->deriv2_classes[0][2][31];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+815, dvrr_stack+1194, NULL);
 tmp = dvrr_stack + 815;
 target_ptr = Libderiv->deriv2_classes[0][3][31];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+825, dvrr_stack+1383, NULL);
 tmp = dvrr_stack + 825;
 target_ptr = Libderiv->deriv2_classes[0][4][31];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+1773, dvrr_stack+457, NULL);
 tmp = dvrr_stack + 1773;
 target_ptr = Libderiv->deriv_classes[0][2][2];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+1779, dvrr_stack+1224, NULL);
 tmp = dvrr_stack + 1779;
 target_ptr = Libderiv->deriv2_classes[0][2][30];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1785, dvrr_stack+843, NULL);
 tmp = dvrr_stack + 1785;
 target_ptr = Libderiv->deriv_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1795, dvrr_stack+1428, NULL);
 tmp = dvrr_stack + 1795;
 target_ptr = Libderiv->deriv2_classes[0][3][30];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1805, dvrr_stack+1458, NULL);
 tmp = dvrr_stack + 1805;
 target_ptr = Libderiv->deriv2_classes[0][4][30];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+1820, dvrr_stack+1242, NULL);
 tmp = dvrr_stack + 1820;
 target_ptr = Libderiv->deriv2_classes[0][2][26];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1826, dvrr_stack+9, NULL);
 tmp = dvrr_stack + 1826;
 target_ptr = Libderiv->deriv2_classes[0][3][26];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1836, dvrr_stack+1593, NULL);
 tmp = dvrr_stack + 1836;
 target_ptr = Libderiv->deriv2_classes[0][4][26];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,6,dvrr_stack+1851, dvrr_stack+927, NULL);
 tmp = dvrr_stack + 1851;
 target_ptr = Libderiv->deriv2_classes[0][2][23];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,10,dvrr_stack+1857, dvrr_stack+1035, NULL);
 tmp = dvrr_stack + 1857;
 target_ptr = Libderiv->deriv2_classes[0][3][23];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,15,dvrr_stack+1867, dvrr_stack+1263, NULL);
 tmp = dvrr_stack + 1867;
 target_ptr = Libderiv->deriv2_classes[0][4][23];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1882, dvrr_stack+436, NULL);
 tmp = dvrr_stack + 1882;
 target_ptr = Libderiv->deriv2_classes[0][2][22];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1888, dvrr_stack+1308, NULL);
 tmp = dvrr_stack + 1888;
 target_ptr = Libderiv->deriv2_classes[0][3][22];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+1898, dvrr_stack+1338, NULL);
 tmp = dvrr_stack + 1898;
 target_ptr = Libderiv->deriv2_classes[0][4][22];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1913, dvrr_stack+246, NULL);
 tmp = dvrr_stack + 1913;
 target_ptr = Libderiv->deriv2_classes[0][2][21];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1919, dvrr_stack+873, NULL);
 tmp = dvrr_stack + 1919;
 target_ptr = Libderiv->deriv2_classes[0][3][21];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+1929, dvrr_stack+945, NULL);
 tmp = dvrr_stack + 1929;
 target_ptr = Libderiv->deriv2_classes[0][4][21];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1944, dvrr_stack+1128, NULL);
 tmp = dvrr_stack + 1944;
 target_ptr = Libderiv->deriv2_classes[0][2][20];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1950, dvrr_stack+1146, NULL);
 tmp = dvrr_stack + 1950;
 target_ptr = Libderiv->deriv2_classes[0][3][20];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+1960, dvrr_stack+990, NULL);
 tmp = dvrr_stack + 1960;
 target_ptr = Libderiv->deriv2_classes[0][4][20];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+1975, dvrr_stack+1176, NULL);
 tmp = dvrr_stack + 1975;
 target_ptr = Libderiv->deriv2_classes[0][2][19];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+1981, dvrr_stack+1194, NULL);
 tmp = dvrr_stack + 1981;
 target_ptr = Libderiv->deriv2_classes[0][3][19];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+1991, dvrr_stack+1383, NULL);
 tmp = dvrr_stack + 1991;
 target_ptr = Libderiv->deriv2_classes[0][4][19];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+2006, dvrr_stack+457, NULL);
 tmp = dvrr_stack + 2006;
 target_ptr = Libderiv->deriv_classes[0][2][1];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+2012, dvrr_stack+1224, NULL);
 tmp = dvrr_stack + 2012;
 target_ptr = Libderiv->deriv2_classes[0][2][18];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+2018, dvrr_stack+843, NULL);
 tmp = dvrr_stack + 2018;
 target_ptr = Libderiv->deriv_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+2028, dvrr_stack+1428, NULL);
 tmp = dvrr_stack + 2028;
 target_ptr = Libderiv->deriv2_classes[0][3][18];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+2038, dvrr_stack+1458, NULL);
 tmp = dvrr_stack + 2038;
 target_ptr = Libderiv->deriv2_classes[0][4][18];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+2053, dvrr_stack+1242, NULL);
 tmp = dvrr_stack + 2053;
 target_ptr = Libderiv->deriv2_classes[0][2][14];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+2059, dvrr_stack+9, NULL);
 tmp = dvrr_stack + 2059;
 target_ptr = Libderiv->deriv2_classes[0][3][14];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+2069, dvrr_stack+1593, NULL);
 tmp = dvrr_stack + 2069;
 target_ptr = Libderiv->deriv2_classes[0][4][14];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+2084, dvrr_stack+1065, NULL);
 tmp = dvrr_stack + 2084;
 target_ptr = Libderiv->deriv2_classes[0][2][13];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+2090, dvrr_stack+1503, NULL);
 tmp = dvrr_stack + 2090;
 target_ptr = Libderiv->deriv2_classes[0][3][13];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+2100, dvrr_stack+1728, NULL);
 tmp = dvrr_stack + 2100;
 target_ptr = Libderiv->deriv2_classes[0][4][13];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,6,dvrr_stack+2115, dvrr_stack+927, NULL);
 tmp = dvrr_stack + 2115;
 target_ptr = Libderiv->deriv2_classes[0][2][11];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,10,dvrr_stack+925, dvrr_stack+1035, NULL);
 tmp = dvrr_stack + 925;
 target_ptr = Libderiv->deriv2_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,15,dvrr_stack+1035, dvrr_stack+1263, NULL);
 tmp = dvrr_stack + 1035;
 target_ptr = Libderiv->deriv2_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+1050, dvrr_stack+436, NULL);
 tmp = dvrr_stack + 1050;
 target_ptr = Libderiv->deriv2_classes[0][2][10];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+935, dvrr_stack+1308, NULL);
 tmp = dvrr_stack + 935;
 target_ptr = Libderiv->deriv2_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+434, dvrr_stack+1338, NULL);
 tmp = dvrr_stack + 434;
 target_ptr = Libderiv->deriv2_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+1056, dvrr_stack+246, NULL);
 tmp = dvrr_stack + 1056;
 target_ptr = Libderiv->deriv2_classes[0][2][9];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+244, dvrr_stack+873, NULL);
 tmp = dvrr_stack + 244;
 target_ptr = Libderiv->deriv2_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+873, dvrr_stack+945, NULL);
 tmp = dvrr_stack + 873;
 target_ptr = Libderiv->deriv2_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+945, dvrr_stack+1128, NULL);
 tmp = dvrr_stack + 945;
 target_ptr = Libderiv->deriv2_classes[0][2][8];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+254, dvrr_stack+1146, NULL);
 tmp = dvrr_stack + 254;
 target_ptr = Libderiv->deriv2_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+888, dvrr_stack+990, NULL);
 tmp = dvrr_stack + 888;
 target_ptr = Libderiv->deriv2_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+1128, dvrr_stack+1176, NULL);
 tmp = dvrr_stack + 1128;
 target_ptr = Libderiv->deriv2_classes[0][2][7];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1134, dvrr_stack+1194, NULL);
 tmp = dvrr_stack + 1134;
 target_ptr = Libderiv->deriv2_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1144, dvrr_stack+1383, NULL);
 tmp = dvrr_stack + 1144;
 target_ptr = Libderiv->deriv2_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+1159, dvrr_stack+457, NULL);
 tmp = dvrr_stack + 1159;
 target_ptr = Libderiv->deriv_classes[0][2][0];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+1165, dvrr_stack+1224, NULL);
 tmp = dvrr_stack + 1165;
 target_ptr = Libderiv->deriv2_classes[0][2][6];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1171, dvrr_stack+843, NULL);
 tmp = dvrr_stack + 1171;
 target_ptr = Libderiv->deriv_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1181, dvrr_stack+1428, NULL);
 tmp = dvrr_stack + 1181;
 target_ptr = Libderiv->deriv2_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1191, dvrr_stack+1458, NULL);
 tmp = dvrr_stack + 1191;
 target_ptr = Libderiv->deriv2_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+1206, dvrr_stack+1242, NULL);
 tmp = dvrr_stack + 1206;
 target_ptr = Libderiv->deriv2_classes[0][2][2];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1212, dvrr_stack+9, NULL);
 tmp = dvrr_stack + 1212;
 target_ptr = Libderiv->deriv2_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1222, dvrr_stack+1593, NULL);
 tmp = dvrr_stack + 1222;
 target_ptr = Libderiv->deriv2_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+1593, dvrr_stack+1065, NULL);
 tmp = dvrr_stack + 1593;
 target_ptr = Libderiv->deriv2_classes[0][2][1];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1599, dvrr_stack+1503, NULL);
 tmp = dvrr_stack + 1599;
 target_ptr = Libderiv->deriv2_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1609, dvrr_stack+1728, NULL);
 tmp = dvrr_stack + 1609;
 target_ptr = Libderiv->deriv2_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+1624, dvrr_stack+39, NULL);
 tmp = dvrr_stack + 1624;
 target_ptr = Libderiv->deriv2_classes[0][2][0];
 for(i=0;i<6;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1237, dvrr_stack+1083, NULL);
 tmp = dvrr_stack + 1237;
 target_ptr = Libderiv->deriv2_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1247, dvrr_stack+1533, NULL);
 tmp = dvrr_stack + 1247;
 target_ptr = Libderiv->deriv2_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];


}

