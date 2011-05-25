#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (p0|fd) integrals */

void d12vrr_order_p0fd(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+34, dvrr_stack+31, dvrr_stack+3, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+34, dvrr_stack+15, dvrr_stack+31, dvrr_stack+3, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+50, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+53, dvrr_stack+6, dvrr_stack+50, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+59, dvrr_stack+9, dvrr_stack+53, dvrr_stack+0, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+69, dvrr_stack+21, dvrr_stack+59, dvrr_stack+15, dvrr_stack+9, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+84, dvrr_stack+40, dvrr_stack+21, dvrr_stack+34, dvrr_stack+15, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+99, dvrr_stack+84, dvrr_stack+69, dvrr_stack+40, dvrr_stack+21, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+120, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+123, dvrr_stack+50, dvrr_stack+120, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+129, dvrr_stack+53, dvrr_stack+123, dvrr_stack+6, dvrr_stack+50, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+139, dvrr_stack+59, dvrr_stack+129, dvrr_stack+9, dvrr_stack+53, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+154, dvrr_stack+69, dvrr_stack+139, NULL, NULL, dvrr_stack+59);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+199, dvrr_stack+69, dvrr_stack+139, dvrr_stack+21, dvrr_stack+59, NULL);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+220, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+223, dvrr_stack+120, dvrr_stack+220, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+229, dvrr_stack+123, dvrr_stack+223, dvrr_stack+50, dvrr_stack+120, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+239, dvrr_stack+129, dvrr_stack+229, dvrr_stack+53, dvrr_stack+123, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+254, dvrr_stack+139, dvrr_stack+239, dvrr_stack+59, dvrr_stack+129, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+275, dvrr_stack+199, dvrr_stack+254, NULL, NULL, dvrr_stack+139);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+338, dvrr_stack+99, dvrr_stack+199, NULL, NULL, dvrr_stack+69);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+401, dvrr_stack+338, dvrr_stack+275, dvrr_stack+99, dvrr_stack+199, dvrr_stack+154);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+527, dvrr_stack+40, dvrr_stack+21, NULL, NULL, dvrr_stack+15);
 tmp = dvrr_stack + 527;
 target_ptr = Libderiv->dvrr_classes[1][3];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+557, dvrr_stack+84, dvrr_stack+69, NULL, NULL, dvrr_stack+21);
 tmp = dvrr_stack + 557;
 target_ptr = Libderiv->dvrr_classes[1][4];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+602,dvrr_stack+557,dvrr_stack+527,3);


 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+692,dvrr_stack+338,dvrr_stack+557,3);


 /* compute (1 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+827,dvrr_stack+692,dvrr_stack+602,3);


 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,30,dvrr_stack+1007, dvrr_stack+827, dvrr_stack+527);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1097, dvrr_stack+199, dvrr_stack+254, dvrr_stack+69, dvrr_stack+139, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1125, dvrr_stack+99, dvrr_stack+199, dvrr_stack+84, dvrr_stack+69, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1153, dvrr_stack+1125, dvrr_stack+1097, NULL, NULL, dvrr_stack+199);

 /* compute (1 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+1237,dvrr_stack+1153,dvrr_stack+338,3);


 /* compute (1 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+1426,dvrr_stack+1237,dvrr_stack+692,3);


 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,45,dvrr_stack+1696, dvrr_stack+1426, dvrr_stack+557);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+50, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1831, dvrr_stack+220, dvrr_stack+50, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1837, dvrr_stack+223, dvrr_stack+1831, dvrr_stack+120, dvrr_stack+220, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1847, dvrr_stack+229, dvrr_stack+1837, dvrr_stack+123, dvrr_stack+223, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1862, dvrr_stack+239, dvrr_stack+1847, dvrr_stack+129, dvrr_stack+229, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1831, dvrr_stack+254, dvrr_stack+1862, dvrr_stack+139, dvrr_stack+239, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1883, dvrr_stack+1097, dvrr_stack+1831, dvrr_stack+199, dvrr_stack+254, NULL);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1919, dvrr_stack+1125, dvrr_stack+1097, dvrr_stack+99, dvrr_stack+199, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1955, dvrr_stack+1919, dvrr_stack+1883, NULL, NULL, dvrr_stack+1097);

 /* compute (1 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+2063,dvrr_stack+1955,dvrr_stack+1153,3);


 /* compute (1 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+2315,dvrr_stack+2063,dvrr_stack+1237,3);


 /* compute (1 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,63,dvrr_stack+2693, dvrr_stack+2315, dvrr_stack+338);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,30,dvrr_stack+2882, dvrr_stack+827, dvrr_stack+527);

 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,45,dvrr_stack+2972, dvrr_stack+1426, dvrr_stack+557);

 /* compute (1 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,63,dvrr_stack+3107, dvrr_stack+2315, dvrr_stack+338);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,30,dvrr_stack+3296, dvrr_stack+827, dvrr_stack+527);

 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,45,dvrr_stack+827, dvrr_stack+1426, dvrr_stack+557);

 /* compute (1 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,63,dvrr_stack+1426, dvrr_stack+2315, dvrr_stack+338);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2315, dvrr_stack+34, dvrr_stack+15, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+2333,dvrr_stack+527,dvrr_stack+2315,3);


 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+2387, dvrr_stack+2333, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+962, dvrr_stack+692, NULL);
 tmp = dvrr_stack + 962;
 target_ptr = Libderiv->deriv_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+2405, dvrr_stack+602, NULL);
 tmp = dvrr_stack + 2405;
 target_ptr = Libderiv->deriv_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,63,dvrr_stack+2435, dvrr_stack+1237, NULL);
 tmp = dvrr_stack + 2435;
 target_ptr = Libderiv->deriv_classes[1][5][11];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,84,dvrr_stack+2498, dvrr_stack+2063, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+2582, dvrr_stack+2333, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+2600, dvrr_stack+692, NULL);
 tmp = dvrr_stack + 2600;
 target_ptr = Libderiv->deriv_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+2645, dvrr_stack+602, NULL);
 tmp = dvrr_stack + 2645;
 target_ptr = Libderiv->deriv_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+1615, dvrr_stack+1237, NULL);
 tmp = dvrr_stack + 1615;
 target_ptr = Libderiv->deriv_classes[1][5][10];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,84,dvrr_stack+3386, dvrr_stack+2063, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+2675, dvrr_stack+2333, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+2333, dvrr_stack+692, NULL);
 tmp = dvrr_stack + 2333;
 target_ptr = Libderiv->deriv_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+692, dvrr_stack+602, NULL);
 tmp = dvrr_stack + 692;
 target_ptr = Libderiv->deriv_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+602, dvrr_stack+1237, NULL);
 tmp = dvrr_stack + 602;
 target_ptr = Libderiv->deriv_classes[1][5][9];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,84,dvrr_stack+1237, dvrr_stack+2063, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+2378, dvrr_stack+31, dvrr_stack+3, NULL, NULL, Data->F+1);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+1678, dvrr_stack+527, dvrr_stack+2378);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+2063, dvrr_stack+338, dvrr_stack+527);
 tmp = dvrr_stack + 2063;
 target_ptr = Libderiv->deriv_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2108, dvrr_stack+557, dvrr_stack+2315);
 tmp = dvrr_stack + 2108;
 target_ptr = Libderiv->deriv_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+2138, dvrr_stack+1153, dvrr_stack+557);
 tmp = dvrr_stack + 2138;
 target_ptr = Libderiv->deriv_classes[1][5][8];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+2201, dvrr_stack+1955, dvrr_stack+338);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+2285, dvrr_stack+527, dvrr_stack+2378);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+1321, dvrr_stack+338, dvrr_stack+527);
 tmp = dvrr_stack + 1321;
 target_ptr = Libderiv->deriv_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1366, dvrr_stack+557, dvrr_stack+2315);
 tmp = dvrr_stack + 1366;
 target_ptr = Libderiv->deriv_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+722, dvrr_stack+1153, dvrr_stack+557);
 tmp = dvrr_stack + 722;
 target_ptr = Libderiv->deriv_classes[1][5][7];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+3470, dvrr_stack+1955, dvrr_stack+338);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+1396, dvrr_stack+527, dvrr_stack+2378);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+1883, dvrr_stack+338, dvrr_stack+527);
 tmp = dvrr_stack + 1883;
 target_ptr = Libderiv->deriv_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+785, dvrr_stack+557, dvrr_stack+2315);
 tmp = dvrr_stack + 785;
 target_ptr = Libderiv->deriv_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+3554, dvrr_stack+1153, dvrr_stack+557);
 tmp = dvrr_stack + 3554;
 target_ptr = Libderiv->deriv_classes[1][5][6];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+3617, dvrr_stack+1955, dvrr_stack+338);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+1928,dvrr_stack+84,dvrr_stack+40,1);


 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+2303, dvrr_stack+1928, NULL);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+665, dvrr_stack+15, dvrr_stack+9, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1958, dvrr_stack+21, dvrr_stack+59, NULL, NULL, dvrr_stack+9);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1988, dvrr_stack+527, dvrr_stack+1958, dvrr_stack+40, dvrr_stack+21, dvrr_stack+665);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3701, dvrr_stack+557, dvrr_stack+154, dvrr_stack+84, dvrr_stack+69, dvrr_stack+1958);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+3791,dvrr_stack+3701,dvrr_stack+1988,6);


 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+3971, dvrr_stack+3791, NULL);

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+4031,dvrr_stack+99,dvrr_stack+84,1);


 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+2048, dvrr_stack+4031, NULL);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+4076,dvrr_stack+401,dvrr_stack+3701,6);


 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+4346, dvrr_stack+4076, NULL);

 /* compute (0 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+4436,dvrr_stack+1125,dvrr_stack+99,1);


 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,21,dvrr_stack+4499, dvrr_stack+4436, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4520, dvrr_stack+1097, dvrr_stack+1831, NULL, NULL, dvrr_stack+254);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4604, dvrr_stack+1153, dvrr_stack+4520, dvrr_stack+1125, dvrr_stack+1097, dvrr_stack+275);

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+4772,dvrr_stack+4604,dvrr_stack+401,6);


 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+5150, dvrr_stack+4772, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+1097, dvrr_stack+1928, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+4520, dvrr_stack+3791, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+1107, dvrr_stack+4031, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+5276, dvrr_stack+4076, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+4580, dvrr_stack+4436, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+5366, dvrr_stack+4772, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+1153, dvrr_stack+1928, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+1163, dvrr_stack+3791, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+3791, dvrr_stack+4031, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+3806, dvrr_stack+4076, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+4031, dvrr_stack+4436, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+4052, dvrr_stack+4772, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+4772, dvrr_stack+84, dvrr_stack+34);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+2378, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+4782, dvrr_stack+2315, dvrr_stack+665, dvrr_stack+34, dvrr_stack+15, dvrr_stack+2378);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+4818, dvrr_stack+3701, dvrr_stack+4782);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+4878, dvrr_stack+99, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+4893, dvrr_stack+401, dvrr_stack+1988);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+4983, dvrr_stack+1125, dvrr_stack+84);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+5004, dvrr_stack+4604, dvrr_stack+3701);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+5130, dvrr_stack+84, dvrr_stack+34);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+4436, dvrr_stack+3701, dvrr_stack+4782);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+4178, dvrr_stack+99, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+4193, dvrr_stack+401, dvrr_stack+1988);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+4283, dvrr_stack+1125, dvrr_stack+84);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+5492, dvrr_stack+4604, dvrr_stack+3701);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+5140, dvrr_stack+84, dvrr_stack+34);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+3896, dvrr_stack+3701, dvrr_stack+4782);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+3956, dvrr_stack+99, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+5618, dvrr_stack+401, dvrr_stack+1988);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+4782, dvrr_stack+1125, dvrr_stack+84);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+5708, dvrr_stack+4604, dvrr_stack+3701);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+4803, dvrr_stack+527, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+2378, dvrr_stack+0, dvrr_stack+6, NULL, NULL, Data->F+3);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+4304, dvrr_stack+9, dvrr_stack+53, NULL, NULL, dvrr_stack+6);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+4601, dvrr_stack+665, dvrr_stack+4304, dvrr_stack+15, dvrr_stack+9, dvrr_stack+2378);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1928, dvrr_stack+59, dvrr_stack+129, NULL, NULL, dvrr_stack+53);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4637, dvrr_stack+1958, dvrr_stack+1928, dvrr_stack+21, dvrr_stack+59, dvrr_stack+4304);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+5834, dvrr_stack+1988, dvrr_stack+4637, dvrr_stack+527, dvrr_stack+1958, dvrr_stack+4601);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+4697, dvrr_stack+5834, dvrr_stack+527);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+4757, dvrr_stack+557, NULL);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5934, dvrr_stack+139, dvrr_stack+239, NULL, NULL, dvrr_stack+129);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5979, dvrr_stack+154, dvrr_stack+5934, dvrr_stack+69, dvrr_stack+139, dvrr_stack+1928);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6069, dvrr_stack+3701, dvrr_stack+5979, dvrr_stack+557, dvrr_stack+154, dvrr_stack+4637);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+4601, dvrr_stack+6069, dvrr_stack+557);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+1928, dvrr_stack+338, NULL);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+120, dvrr_stack+254, dvrr_stack+1862, NULL, NULL, dvrr_stack+239);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6219, dvrr_stack+275, dvrr_stack+120, dvrr_stack+199, dvrr_stack+254, dvrr_stack+5934);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6345, dvrr_stack+401, dvrr_stack+6219, dvrr_stack+338, dvrr_stack+275, dvrr_stack+5979);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+6219, dvrr_stack+6345, dvrr_stack+338);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+5934, dvrr_stack+527, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+5944, dvrr_stack+5834, dvrr_stack+527);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+6004, dvrr_stack+557, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+120, dvrr_stack+6069, dvrr_stack+557);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+6019, dvrr_stack+338, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+210, dvrr_stack+6345, dvrr_stack+338);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+6040, dvrr_stack+527, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+6555, dvrr_stack+5834, dvrr_stack+527);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+527, dvrr_stack+557, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+5834, dvrr_stack+6069, dvrr_stack+557);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+542, dvrr_stack+338, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+6050, dvrr_stack+6345, dvrr_stack+338);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+6345, dvrr_stack+401, dvrr_stack+99);
 tmp = dvrr_stack + 6345;
 target_ptr = Libderiv->deriv_classes[1][5][2];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+6408, dvrr_stack+401, dvrr_stack+99);
 tmp = dvrr_stack + 6408;
 target_ptr = Libderiv->deriv_classes[1][5][1];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+6471, dvrr_stack+401, dvrr_stack+99);
 tmp = dvrr_stack + 6471;
 target_ptr = Libderiv->deriv_classes[1][5][0];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,30,dvrr_stack+563, dvrr_stack+1007, NULL);
 tmp = dvrr_stack + 563;
 target_ptr = Libderiv->deriv2_classes[1][3][143];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,45,dvrr_stack+336, dvrr_stack+1696, NULL);
 tmp = dvrr_stack + 336;
 target_ptr = Libderiv->deriv2_classes[1][4][143];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,63,dvrr_stack+381, dvrr_stack+2693, NULL);
 tmp = dvrr_stack + 381;
 target_ptr = Libderiv->deriv2_classes[1][5][143];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,30,dvrr_stack+6176, dvrr_stack+1007, NULL);
 tmp = dvrr_stack + 6176;
 target_ptr = Libderiv->deriv2_classes[1][3][131];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,45,dvrr_stack+444, dvrr_stack+1696, NULL);
 tmp = dvrr_stack + 444;
 target_ptr = Libderiv->deriv2_classes[1][4][131];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,63,dvrr_stack+6615, dvrr_stack+2693, NULL);
 tmp = dvrr_stack + 6615;
 target_ptr = Libderiv->deriv2_classes[1][5][131];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+489, dvrr_stack+2882, NULL);
 tmp = dvrr_stack + 489;
 target_ptr = Libderiv->deriv2_classes[1][3][130];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+1831, dvrr_stack+2972, NULL);
 tmp = dvrr_stack + 1831;
 target_ptr = Libderiv->deriv2_classes[1][4][130];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+6678, dvrr_stack+3107, NULL);
 tmp = dvrr_stack + 6678;
 target_ptr = Libderiv->deriv2_classes[1][5][130];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,30,dvrr_stack+1949, dvrr_stack+1007, NULL);
 tmp = dvrr_stack + 1949;
 target_ptr = Libderiv->deriv2_classes[1][3][119];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,45,dvrr_stack+1007, dvrr_stack+1696, NULL);
 tmp = dvrr_stack + 1007;
 target_ptr = Libderiv->deriv2_classes[1][4][119];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,63,dvrr_stack+1696, dvrr_stack+2693, NULL);
 tmp = dvrr_stack + 1696;
 target_ptr = Libderiv->deriv2_classes[1][5][119];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+2693, dvrr_stack+2882, NULL);
 tmp = dvrr_stack + 2693;
 target_ptr = Libderiv->deriv2_classes[1][3][118];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+1052, dvrr_stack+2972, NULL);
 tmp = dvrr_stack + 1052;
 target_ptr = Libderiv->deriv2_classes[1][4][118];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+2723, dvrr_stack+3107, NULL);
 tmp = dvrr_stack + 2723;
 target_ptr = Libderiv->deriv2_classes[1][5][118];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+2786, dvrr_stack+3296, NULL);
 tmp = dvrr_stack + 2786;
 target_ptr = Libderiv->deriv2_classes[1][3][117];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+2816, dvrr_stack+827, NULL);
 tmp = dvrr_stack + 2816;
 target_ptr = Libderiv->deriv2_classes[1][4][117];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+2861, dvrr_stack+1426, NULL);
 tmp = dvrr_stack + 2861;
 target_ptr = Libderiv->deriv2_classes[1][5][117];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2924, dvrr_stack+962, dvrr_stack+2387);
 tmp = dvrr_stack + 2924;
 target_ptr = Libderiv->deriv2_classes[1][3][107];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+2954, dvrr_stack+2435, dvrr_stack+2405);
 tmp = dvrr_stack + 2954;
 target_ptr = Libderiv->deriv2_classes[1][4][107];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+2999, dvrr_stack+2498, dvrr_stack+962);
 tmp = dvrr_stack + 2999;
 target_ptr = Libderiv->deriv2_classes[1][5][107];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+3062, dvrr_stack+2600, dvrr_stack+2582);
 tmp = dvrr_stack + 3062;
 target_ptr = Libderiv->deriv2_classes[1][3][106];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+3092, dvrr_stack+1615, dvrr_stack+2645);
 tmp = dvrr_stack + 3092;
 target_ptr = Libderiv->deriv2_classes[1][4][106];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+3137, dvrr_stack+3386, dvrr_stack+2600);
 tmp = dvrr_stack + 3137;
 target_ptr = Libderiv->deriv2_classes[1][5][106];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+3200, dvrr_stack+2333, dvrr_stack+2675);
 tmp = dvrr_stack + 3200;
 target_ptr = Libderiv->deriv2_classes[1][3][105];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+3230, dvrr_stack+602, dvrr_stack+692);
 tmp = dvrr_stack + 3230;
 target_ptr = Libderiv->deriv2_classes[1][4][105];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+3275, dvrr_stack+1237, dvrr_stack+2333);
 tmp = dvrr_stack + 3275;
 target_ptr = Libderiv->deriv2_classes[1][5][105];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+3338, dvrr_stack+2063, dvrr_stack+1678);
 tmp = dvrr_stack + 3338;
 target_ptr = Libderiv->deriv2_classes[1][3][104];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+1759, dvrr_stack+2138, dvrr_stack+2108);
 tmp = dvrr_stack + 1759;
 target_ptr = Libderiv->deriv2_classes[1][4][104];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+1414, dvrr_stack+2201, dvrr_stack+2063);
 tmp = dvrr_stack + 1414;
 target_ptr = Libderiv->deriv2_classes[1][5][104];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+4304, dvrr_stack+962, dvrr_stack+2387);
 tmp = dvrr_stack + 4304;
 target_ptr = Libderiv->deriv2_classes[1][3][95];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+1477, dvrr_stack+2435, dvrr_stack+2405);
 tmp = dvrr_stack + 1477;
 target_ptr = Libderiv->deriv2_classes[1][4][95];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+1522, dvrr_stack+2498, dvrr_stack+962);
 tmp = dvrr_stack + 1522;
 target_ptr = Libderiv->deriv2_classes[1][5][95];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1585, dvrr_stack+2600, dvrr_stack+2582);
 tmp = dvrr_stack + 1585;
 target_ptr = Libderiv->deriv2_classes[1][3][94];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+815, dvrr_stack+1615, dvrr_stack+2645);
 tmp = dvrr_stack + 815;
 target_ptr = Libderiv->deriv2_classes[1][4][94];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+860, dvrr_stack+3386, dvrr_stack+2600);
 tmp = dvrr_stack + 860;
 target_ptr = Libderiv->deriv2_classes[1][5][94];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+0, dvrr_stack+2333, dvrr_stack+2675);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[1][3][93];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+6741, dvrr_stack+602, dvrr_stack+692);
 tmp = dvrr_stack + 6741;
 target_ptr = Libderiv->deriv2_classes[1][4][93];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+6786, dvrr_stack+1237, dvrr_stack+2333);
 tmp = dvrr_stack + 6786;
 target_ptr = Libderiv->deriv2_classes[1][5][93];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1122, dvrr_stack+2063, dvrr_stack+1678);
 tmp = dvrr_stack + 1122;
 target_ptr = Libderiv->deriv2_classes[1][3][92];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+6849, dvrr_stack+2138, dvrr_stack+2108);
 tmp = dvrr_stack + 6849;
 target_ptr = Libderiv->deriv2_classes[1][4][92];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+6894, dvrr_stack+2201, dvrr_stack+2063);
 tmp = dvrr_stack + 6894;
 target_ptr = Libderiv->deriv2_classes[1][5][92];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+923, dvrr_stack+1321, dvrr_stack+2285);
 tmp = dvrr_stack + 923;
 target_ptr = Libderiv->deriv2_classes[1][3][91];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+6957, dvrr_stack+722, dvrr_stack+1366);
 tmp = dvrr_stack + 6957;
 target_ptr = Libderiv->deriv2_classes[1][4][91];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+7002, dvrr_stack+3470, dvrr_stack+1321);
 tmp = dvrr_stack + 7002;
 target_ptr = Libderiv->deriv2_classes[1][5][91];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+50, dvrr_stack+962, dvrr_stack+2387);
 tmp = dvrr_stack + 50;
 target_ptr = Libderiv->deriv2_classes[1][3][83];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+7065, dvrr_stack+2435, dvrr_stack+2405);
 tmp = dvrr_stack + 7065;
 target_ptr = Libderiv->deriv2_classes[1][4][83];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+2378, dvrr_stack+2498, dvrr_stack+962);
 tmp = dvrr_stack + 2378;
 target_ptr = Libderiv->deriv2_classes[1][5][83];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2441, dvrr_stack+2600, dvrr_stack+2582);
 tmp = dvrr_stack + 2441;
 target_ptr = Libderiv->deriv2_classes[1][3][82];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+2471, dvrr_stack+1615, dvrr_stack+2645);
 tmp = dvrr_stack + 2471;
 target_ptr = Libderiv->deriv2_classes[1][4][82];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+1615, dvrr_stack+3386, dvrr_stack+2600);
 tmp = dvrr_stack + 1615;
 target_ptr = Libderiv->deriv2_classes[1][5][82];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+3368, dvrr_stack+2333, dvrr_stack+2675);
 tmp = dvrr_stack + 3368;
 target_ptr = Libderiv->deriv2_classes[1][3][81];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+3398, dvrr_stack+602, dvrr_stack+692);
 tmp = dvrr_stack + 3398;
 target_ptr = Libderiv->deriv2_classes[1][4][81];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+593, dvrr_stack+1237, dvrr_stack+2333);
 tmp = dvrr_stack + 593;
 target_ptr = Libderiv->deriv2_classes[1][5][81];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+656, dvrr_stack+2063, dvrr_stack+1678);
 tmp = dvrr_stack + 656;
 target_ptr = Libderiv->deriv2_classes[1][3][80];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+2516, dvrr_stack+2138, dvrr_stack+2108);
 tmp = dvrr_stack + 2516;
 target_ptr = Libderiv->deriv2_classes[1][4][80];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+2108, dvrr_stack+2201, dvrr_stack+2063);
 tmp = dvrr_stack + 2108;
 target_ptr = Libderiv->deriv2_classes[1][5][80];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2063, dvrr_stack+1321, dvrr_stack+2285);
 tmp = dvrr_stack + 2063;
 target_ptr = Libderiv->deriv2_classes[1][3][79];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+2171, dvrr_stack+722, dvrr_stack+1366);
 tmp = dvrr_stack + 2171;
 target_ptr = Libderiv->deriv2_classes[1][4][79];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+2216, dvrr_stack+3470, dvrr_stack+1321);
 tmp = dvrr_stack + 2216;
 target_ptr = Libderiv->deriv2_classes[1][5][79];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+3443, dvrr_stack+1883, dvrr_stack+1396);
 tmp = dvrr_stack + 3443;
 target_ptr = Libderiv->deriv2_classes[1][3][78];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+3473, dvrr_stack+3554, dvrr_stack+785);
 tmp = dvrr_stack + 3473;
 target_ptr = Libderiv->deriv2_classes[1][4][78];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+3518, dvrr_stack+3617, dvrr_stack+1883);
 tmp = dvrr_stack + 3518;
 target_ptr = Libderiv->deriv2_classes[1][5][78];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,10,dvrr_stack+3581, dvrr_stack+3971, dvrr_stack+2303);
 tmp = dvrr_stack + 3581;
 target_ptr = Libderiv->deriv2_classes[1][3][35];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,15,dvrr_stack+3611, dvrr_stack+4346, dvrr_stack+2048);
 tmp = dvrr_stack + 3611;
 target_ptr = Libderiv->deriv2_classes[1][4][35];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,21,dvrr_stack+686, dvrr_stack+5150, dvrr_stack+4499);
 tmp = dvrr_stack + 686;
 target_ptr = Libderiv->deriv2_classes[1][5][35];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+3656, dvrr_stack+4520, dvrr_stack+1097);
 tmp = dvrr_stack + 3656;
 target_ptr = Libderiv->deriv2_classes[1][3][34];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+749, dvrr_stack+5276, dvrr_stack+1107);
 tmp = dvrr_stack + 749;
 target_ptr = Libderiv->deriv2_classes[1][4][34];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+2561, dvrr_stack+5366, dvrr_stack+4580);
 tmp = dvrr_stack + 2561;
 target_ptr = Libderiv->deriv2_classes[1][5][34];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+2624, dvrr_stack+1163, dvrr_stack+1153);
 tmp = dvrr_stack + 2624;
 target_ptr = Libderiv->deriv2_classes[1][3][33];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1223, dvrr_stack+3806, dvrr_stack+3791);
 tmp = dvrr_stack + 1223;
 target_ptr = Libderiv->deriv2_classes[1][4][33];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+1268, dvrr_stack+4052, dvrr_stack+4031);
 tmp = dvrr_stack + 1268;
 target_ptr = Libderiv->deriv2_classes[1][5][33];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+2654, dvrr_stack+4818, dvrr_stack+4772);
 tmp = dvrr_stack + 2654;
 target_ptr = Libderiv->deriv2_classes[1][3][32];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1331, dvrr_stack+4893, dvrr_stack+4878);
 tmp = dvrr_stack + 1331;
 target_ptr = Libderiv->deriv2_classes[1][4][32];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+2313, dvrr_stack+5004, dvrr_stack+4983);
 tmp = dvrr_stack + 2313;
 target_ptr = Libderiv->deriv2_classes[1][5][32];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+1376, dvrr_stack+4436, dvrr_stack+5130);
 tmp = dvrr_stack + 1376;
 target_ptr = Libderiv->deriv2_classes[1][3][31];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1876, dvrr_stack+4193, dvrr_stack+4178);
 tmp = dvrr_stack + 1876;
 target_ptr = Libderiv->deriv2_classes[1][4][31];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+7110, dvrr_stack+5492, dvrr_stack+4283);
 tmp = dvrr_stack + 7110;
 target_ptr = Libderiv->deriv2_classes[1][5][31];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+953, dvrr_stack+1988, dvrr_stack+40);
 tmp = dvrr_stack + 953;
 target_ptr = Libderiv->deriv_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+7173, dvrr_stack+3896, dvrr_stack+5140);
 tmp = dvrr_stack + 7173;
 target_ptr = Libderiv->deriv2_classes[1][3][30];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+7203, dvrr_stack+3701, dvrr_stack+84);
 tmp = dvrr_stack + 7203;
 target_ptr = Libderiv->deriv_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+7248, dvrr_stack+5618, dvrr_stack+3956);
 tmp = dvrr_stack + 7248;
 target_ptr = Libderiv->deriv2_classes[1][4][30];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+7293, dvrr_stack+5708, dvrr_stack+4782);
 tmp = dvrr_stack + 7293;
 target_ptr = Libderiv->deriv2_classes[1][5][30];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+7356, dvrr_stack+4697, dvrr_stack+4803);
 tmp = dvrr_stack + 7356;
 target_ptr = Libderiv->deriv2_classes[1][3][26];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+7386, dvrr_stack+4601, dvrr_stack+4757);
 tmp = dvrr_stack + 7386;
 target_ptr = Libderiv->deriv2_classes[1][4][26];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+7431, dvrr_stack+6219, dvrr_stack+1928);
 tmp = dvrr_stack + 7431;
 target_ptr = Libderiv->deriv2_classes[1][5][26];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,10,dvrr_stack+7494, dvrr_stack+3971, dvrr_stack+2303);
 tmp = dvrr_stack + 7494;
 target_ptr = Libderiv->deriv2_classes[1][3][23];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,15,dvrr_stack+7524, dvrr_stack+4346, dvrr_stack+2048);
 tmp = dvrr_stack + 7524;
 target_ptr = Libderiv->deriv2_classes[1][4][23];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,21,dvrr_stack+7569, dvrr_stack+5150, dvrr_stack+4499);
 tmp = dvrr_stack + 7569;
 target_ptr = Libderiv->deriv2_classes[1][5][23];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+7632, dvrr_stack+4520, dvrr_stack+1097);
 tmp = dvrr_stack + 7632;
 target_ptr = Libderiv->deriv2_classes[1][3][22];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+7662, dvrr_stack+5276, dvrr_stack+1107);
 tmp = dvrr_stack + 7662;
 target_ptr = Libderiv->deriv2_classes[1][4][22];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+7707, dvrr_stack+5366, dvrr_stack+4580);
 tmp = dvrr_stack + 7707;
 target_ptr = Libderiv->deriv2_classes[1][5][22];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+7770, dvrr_stack+1163, dvrr_stack+1153);
 tmp = dvrr_stack + 7770;
 target_ptr = Libderiv->deriv2_classes[1][3][21];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+7800, dvrr_stack+3806, dvrr_stack+3791);
 tmp = dvrr_stack + 7800;
 target_ptr = Libderiv->deriv2_classes[1][4][21];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+7845, dvrr_stack+4052, dvrr_stack+4031);
 tmp = dvrr_stack + 7845;
 target_ptr = Libderiv->deriv2_classes[1][5][21];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+7908, dvrr_stack+4818, dvrr_stack+4772);
 tmp = dvrr_stack + 7908;
 target_ptr = Libderiv->deriv2_classes[1][3][20];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+7938, dvrr_stack+4893, dvrr_stack+4878);
 tmp = dvrr_stack + 7938;
 target_ptr = Libderiv->deriv2_classes[1][4][20];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+7983, dvrr_stack+5004, dvrr_stack+4983);
 tmp = dvrr_stack + 7983;
 target_ptr = Libderiv->deriv2_classes[1][5][20];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+8046, dvrr_stack+4436, dvrr_stack+5130);
 tmp = dvrr_stack + 8046;
 target_ptr = Libderiv->deriv2_classes[1][3][19];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+8076, dvrr_stack+4193, dvrr_stack+4178);
 tmp = dvrr_stack + 8076;
 target_ptr = Libderiv->deriv2_classes[1][4][19];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+8121, dvrr_stack+5492, dvrr_stack+4283);
 tmp = dvrr_stack + 8121;
 target_ptr = Libderiv->deriv2_classes[1][5][19];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+8184, dvrr_stack+1988, dvrr_stack+40);
 tmp = dvrr_stack + 8184;
 target_ptr = Libderiv->deriv_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+8214, dvrr_stack+3896, dvrr_stack+5140);
 tmp = dvrr_stack + 8214;
 target_ptr = Libderiv->deriv2_classes[1][3][18];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+8244, dvrr_stack+3701, dvrr_stack+84);
 tmp = dvrr_stack + 8244;
 target_ptr = Libderiv->deriv_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+8289, dvrr_stack+5618, dvrr_stack+3956);
 tmp = dvrr_stack + 8289;
 target_ptr = Libderiv->deriv2_classes[1][4][18];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+8334, dvrr_stack+5708, dvrr_stack+4782);
 tmp = dvrr_stack + 8334;
 target_ptr = Libderiv->deriv2_classes[1][5][18];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+8397, dvrr_stack+4697, dvrr_stack+4803);
 tmp = dvrr_stack + 8397;
 target_ptr = Libderiv->deriv2_classes[1][3][14];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+8427, dvrr_stack+4601, dvrr_stack+4757);
 tmp = dvrr_stack + 8427;
 target_ptr = Libderiv->deriv2_classes[1][4][14];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+8472, dvrr_stack+6219, dvrr_stack+1928);
 tmp = dvrr_stack + 8472;
 target_ptr = Libderiv->deriv2_classes[1][5][14];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+8535, dvrr_stack+5944, dvrr_stack+5934);
 tmp = dvrr_stack + 8535;
 target_ptr = Libderiv->deriv2_classes[1][3][13];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+8565, dvrr_stack+120, dvrr_stack+6004);
 tmp = dvrr_stack + 8565;
 target_ptr = Libderiv->deriv2_classes[1][4][13];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+8610, dvrr_stack+210, dvrr_stack+6019);
 tmp = dvrr_stack + 8610;
 target_ptr = Libderiv->deriv2_classes[1][5][13];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,10,dvrr_stack+8673, dvrr_stack+3971, dvrr_stack+2303);
 tmp = dvrr_stack + 8673;
 target_ptr = Libderiv->deriv2_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,15,dvrr_stack+3971, dvrr_stack+4346, dvrr_stack+2048);
 tmp = dvrr_stack + 3971;
 target_ptr = Libderiv->deriv2_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,21,dvrr_stack+4334, dvrr_stack+5150, dvrr_stack+4499);
 tmp = dvrr_stack + 4334;
 target_ptr = Libderiv->deriv2_classes[1][5][11];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5150, dvrr_stack+4520, dvrr_stack+1097);
 tmp = dvrr_stack + 5150;
 target_ptr = Libderiv->deriv2_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+5180, dvrr_stack+5276, dvrr_stack+1107);
 tmp = dvrr_stack + 5180;
 target_ptr = Libderiv->deriv2_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+5225, dvrr_stack+5366, dvrr_stack+4580);
 tmp = dvrr_stack + 5225;
 target_ptr = Libderiv->deriv2_classes[1][5][10];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5288, dvrr_stack+1163, dvrr_stack+1153);
 tmp = dvrr_stack + 5288;
 target_ptr = Libderiv->deriv2_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+5318, dvrr_stack+3806, dvrr_stack+3791);
 tmp = dvrr_stack + 5318;
 target_ptr = Libderiv->deriv2_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+3791, dvrr_stack+4052, dvrr_stack+4031);
 tmp = dvrr_stack + 3791;
 target_ptr = Libderiv->deriv2_classes[1][5][9];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3854, dvrr_stack+4818, dvrr_stack+4772);
 tmp = dvrr_stack + 3854;
 target_ptr = Libderiv->deriv2_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+5363, dvrr_stack+4893, dvrr_stack+4878);
 tmp = dvrr_stack + 5363;
 target_ptr = Libderiv->deriv2_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+5408, dvrr_stack+5004, dvrr_stack+4983);
 tmp = dvrr_stack + 5408;
 target_ptr = Libderiv->deriv2_classes[1][5][8];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+4016, dvrr_stack+4436, dvrr_stack+5130);
 tmp = dvrr_stack + 4016;
 target_ptr = Libderiv->deriv2_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+4046, dvrr_stack+4193, dvrr_stack+4178);
 tmp = dvrr_stack + 4046;
 target_ptr = Libderiv->deriv2_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+4091, dvrr_stack+5492, dvrr_stack+4283);
 tmp = dvrr_stack + 4091;
 target_ptr = Libderiv->deriv2_classes[1][5][7];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5471, dvrr_stack+1988, dvrr_stack+40);
 tmp = dvrr_stack + 5471;
 target_ptr = Libderiv->deriv_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5501, dvrr_stack+3896, dvrr_stack+5140);
 tmp = dvrr_stack + 5501;
 target_ptr = Libderiv->deriv2_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3884, dvrr_stack+3701, dvrr_stack+84);
 tmp = dvrr_stack + 3884;
 target_ptr = Libderiv->deriv_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+5531, dvrr_stack+5618, dvrr_stack+3956);
 tmp = dvrr_stack + 5531;
 target_ptr = Libderiv->deriv2_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+5576, dvrr_stack+5708, dvrr_stack+4782);
 tmp = dvrr_stack + 5576;
 target_ptr = Libderiv->deriv2_classes[1][5][6];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+4772, dvrr_stack+4697, dvrr_stack+4803);
 tmp = dvrr_stack + 4772;
 target_ptr = Libderiv->deriv2_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+4802, dvrr_stack+4601, dvrr_stack+4757);
 tmp = dvrr_stack + 4802;
 target_ptr = Libderiv->deriv2_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+4847, dvrr_stack+6219, dvrr_stack+1928);
 tmp = dvrr_stack + 4847;
 target_ptr = Libderiv->deriv2_classes[1][5][2];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+4910, dvrr_stack+5944, dvrr_stack+5934);
 tmp = dvrr_stack + 4910;
 target_ptr = Libderiv->deriv2_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+4940, dvrr_stack+120, dvrr_stack+6004);
 tmp = dvrr_stack + 4940;
 target_ptr = Libderiv->deriv2_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+4985, dvrr_stack+210, dvrr_stack+6019);
 tmp = dvrr_stack + 4985;
 target_ptr = Libderiv->deriv2_classes[1][5][1];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+5048, dvrr_stack+6555, dvrr_stack+6040);
 tmp = dvrr_stack + 5048;
 target_ptr = Libderiv->deriv2_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+5078, dvrr_stack+5834, dvrr_stack+527);
 tmp = dvrr_stack + 5078;
 target_ptr = Libderiv->deriv2_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+5639, dvrr_stack+6050, dvrr_stack+542);
 tmp = dvrr_stack + 5639;
 target_ptr = Libderiv->deriv2_classes[1][5][0];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];


}

