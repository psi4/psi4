#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (dp|dp) integrals */

void d12vrr_order_dpdp(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+21, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+24, dvrr_stack+21, dvrr_stack+3, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+30, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+33, dvrr_stack+0, dvrr_stack+30, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+39, dvrr_stack+15, dvrr_stack+33, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+57, dvrr_stack+24, dvrr_stack+15, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+75, dvrr_stack+57, dvrr_stack+39, dvrr_stack+24, dvrr_stack+15, dvrr_stack+6);
 tmp = dvrr_stack + 75;
 target_ptr = Libderiv->dvrr_classes[2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+15, dvrr_stack+33, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+121, dvrr_stack+24, dvrr_stack+15, dvrr_stack+21, dvrr_stack+3, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+131, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+134, dvrr_stack+30, dvrr_stack+131, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+140, dvrr_stack+33, dvrr_stack+134, dvrr_stack+0, dvrr_stack+30, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+150, dvrr_stack+111, dvrr_stack+140, NULL, NULL, dvrr_stack+33);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+180, dvrr_stack+121, dvrr_stack+111, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+210, dvrr_stack+180, dvrr_stack+150, dvrr_stack+121, dvrr_stack+111, dvrr_stack+39);
 tmp = dvrr_stack + 210;
 target_ptr = Libderiv->dvrr_classes[2][3];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+270,dvrr_stack+210,dvrr_stack+75,6);


 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+378, dvrr_stack+111, dvrr_stack+140, dvrr_stack+15, dvrr_stack+33, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+393, dvrr_stack+121, dvrr_stack+111, dvrr_stack+24, dvrr_stack+15, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+408, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+411, dvrr_stack+131, dvrr_stack+408, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+417, dvrr_stack+134, dvrr_stack+411, dvrr_stack+30, dvrr_stack+131, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+427, dvrr_stack+140, dvrr_stack+417, dvrr_stack+33, dvrr_stack+134, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+442, dvrr_stack+378, dvrr_stack+427, NULL, NULL, dvrr_stack+140);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+487, dvrr_stack+393, dvrr_stack+378, NULL, NULL, dvrr_stack+111);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+532, dvrr_stack+487, dvrr_stack+442, dvrr_stack+393, dvrr_stack+378, dvrr_stack+150);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+622,dvrr_stack+532,dvrr_stack+210,6);


 /* compute (2 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+802,dvrr_stack+622,dvrr_stack+270,6);


 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,36,dvrr_stack+1018, dvrr_stack+802, dvrr_stack+75);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1126, dvrr_stack+378, dvrr_stack+427, dvrr_stack+111, dvrr_stack+140, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1147, dvrr_stack+393, dvrr_stack+378, dvrr_stack+121, dvrr_stack+111, NULL);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+393, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+396, dvrr_stack+408, dvrr_stack+393, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1168, dvrr_stack+411, dvrr_stack+396, dvrr_stack+131, dvrr_stack+408, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1178, dvrr_stack+417, dvrr_stack+1168, dvrr_stack+134, dvrr_stack+411, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1193, dvrr_stack+427, dvrr_stack+1178, dvrr_stack+140, dvrr_stack+417, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1214, dvrr_stack+1126, dvrr_stack+1193, NULL, NULL, dvrr_stack+427);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1277, dvrr_stack+1147, dvrr_stack+1126, NULL, NULL, dvrr_stack+378);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1340, dvrr_stack+1277, dvrr_stack+1214, dvrr_stack+1147, dvrr_stack+1126, dvrr_stack+442);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+1466,dvrr_stack+1340,dvrr_stack+532,6);


 /* compute (2 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+1736,dvrr_stack+1466,dvrr_stack+622,6);


 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,60,dvrr_stack+2096, dvrr_stack+1736, dvrr_stack+210);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+1147, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1150, dvrr_stack+0, dvrr_stack+30, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+2276, dvrr_stack+6, dvrr_stack+1150, dvrr_stack+3, dvrr_stack+0, dvrr_stack+1147);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2294, dvrr_stack+33, dvrr_stack+134, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+2312, dvrr_stack+39, dvrr_stack+2294, dvrr_stack+15, dvrr_stack+33, dvrr_stack+1150);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+2348, dvrr_stack+75, dvrr_stack+2312, dvrr_stack+57, dvrr_stack+39, dvrr_stack+2276);
 tmp = dvrr_stack + 2348;
 target_ptr = Libderiv->dvrr_classes[3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+2408, dvrr_stack+140, dvrr_stack+417, NULL, NULL, dvrr_stack+134);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+2438, dvrr_stack+150, dvrr_stack+2408, dvrr_stack+111, dvrr_stack+140, dvrr_stack+2294);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+2498, dvrr_stack+210, dvrr_stack+2438, dvrr_stack+180, dvrr_stack+150, dvrr_stack+2312);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+2598,dvrr_stack+2498,dvrr_stack+2348,10);


 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2778, dvrr_stack+427, dvrr_stack+1178, NULL, NULL, dvrr_stack+417);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2823, dvrr_stack+442, dvrr_stack+2778, dvrr_stack+378, dvrr_stack+427, dvrr_stack+2408);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2913, dvrr_stack+532, dvrr_stack+2823, dvrr_stack+487, dvrr_stack+442, dvrr_stack+2438);

 /* compute (3 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+3063,dvrr_stack+2913,dvrr_stack+2498,10);


 /* compute (3 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+3363,dvrr_stack+3063,dvrr_stack+2598,10);


 /* compute (3 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,60,dvrr_stack+3723, dvrr_stack+3363, dvrr_stack+2348);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+378, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+393, dvrr_stack+378, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+111, dvrr_stack+396, dvrr_stack+15, dvrr_stack+408, dvrr_stack+393, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+378, dvrr_stack+1168, dvrr_stack+111, dvrr_stack+411, dvrr_stack+396, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3903, dvrr_stack+1178, dvrr_stack+378, dvrr_stack+417, dvrr_stack+1168, NULL);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3924, dvrr_stack+1193, dvrr_stack+3903, NULL, NULL, dvrr_stack+1178);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3987, dvrr_stack+1214, dvrr_stack+3924, dvrr_stack+1126, dvrr_stack+1193, dvrr_stack+2778);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4113, dvrr_stack+1340, dvrr_stack+3987, dvrr_stack+1277, dvrr_stack+1214, dvrr_stack+2823);

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+4323,dvrr_stack+4113,dvrr_stack+2913,10);


 /* compute (3 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+4773,dvrr_stack+4323,dvrr_stack+3063,10);


 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,100,dvrr_stack+5373, dvrr_stack+4773, dvrr_stack+2498);

 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,36,dvrr_stack+1193, dvrr_stack+802, dvrr_stack+75);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,60,dvrr_stack+3903, dvrr_stack+1736, dvrr_stack+210);

 /* compute (3 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,60,dvrr_stack+5673, dvrr_stack+3363, dvrr_stack+2348);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,100,dvrr_stack+5853, dvrr_stack+4773, dvrr_stack+2498);

 /* compute (2 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,36,dvrr_stack+6153, dvrr_stack+802, dvrr_stack+75);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,60,dvrr_stack+802, dvrr_stack+1736, dvrr_stack+210);

 /* compute (3 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,60,dvrr_stack+1736, dvrr_stack+3363, dvrr_stack+2348);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,100,dvrr_stack+3363, dvrr_stack+4773, dvrr_stack+2498);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+393, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1159, dvrr_stack+21, dvrr_stack+3, NULL, NULL, Data->F+1);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+4773, dvrr_stack+1159, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+393);

 /* compute (2 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+4791,dvrr_stack+75,dvrr_stack+4773,6);


 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+4845, dvrr_stack+4791, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+3663, dvrr_stack+622, NULL);
 tmp = dvrr_stack + 3663;
 target_ptr = Libderiv->deriv_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+982, dvrr_stack+270, NULL);
 tmp = dvrr_stack + 982;
 target_ptr = Libderiv->deriv_classes[2][2][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+4863, dvrr_stack+1466, NULL);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+402, dvrr_stack+393, dvrr_stack+1147, Data->F+1, Data->F+2, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+4083, dvrr_stack+4773, dvrr_stack+2276, dvrr_stack+1159, dvrr_stack+6, dvrr_stack+402);

 /* compute (3 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+4953,dvrr_stack+2348,dvrr_stack+4083,10);


 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+5043, dvrr_stack+4953, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+5073, dvrr_stack+3063, NULL);
 tmp = dvrr_stack + 5073;
 target_ptr = Libderiv->deriv_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+5173, dvrr_stack+2598, NULL);
 tmp = dvrr_stack + 5173;
 target_ptr = Libderiv->deriv_classes[3][2][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+1916, dvrr_stack+4323, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+5233, dvrr_stack+4791, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+5251, dvrr_stack+622, NULL);
 tmp = dvrr_stack + 5251;
 target_ptr = Libderiv->deriv_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+5311, dvrr_stack+270, NULL);
 tmp = dvrr_stack + 5311;
 target_ptr = Libderiv->deriv_classes[2][2][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+6261, dvrr_stack+1466, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+2066, dvrr_stack+4953, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+6351, dvrr_stack+3063, NULL);
 tmp = dvrr_stack + 6351;
 target_ptr = Libderiv->deriv_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+6451, dvrr_stack+2598, NULL);
 tmp = dvrr_stack + 6451;
 target_ptr = Libderiv->deriv_classes[3][2][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+6511, dvrr_stack+4323, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+5347, dvrr_stack+4791, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+6661, dvrr_stack+622, NULL);
 tmp = dvrr_stack + 6661;
 target_ptr = Libderiv->deriv_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+622, dvrr_stack+270, NULL);
 tmp = dvrr_stack + 622;
 target_ptr = Libderiv->deriv_classes[2][2][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+270, dvrr_stack+1466, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+1466, dvrr_stack+4953, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+1496, dvrr_stack+3063, NULL);
 tmp = dvrr_stack + 1496;
 target_ptr = Libderiv->deriv_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+3063, dvrr_stack+2598, NULL);
 tmp = dvrr_stack + 3063;
 target_ptr = Libderiv->deriv_classes[3][2][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+2598, dvrr_stack+4323, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+4323, dvrr_stack+3, dvrr_stack+393, Data->F+0, Data->F+1, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+360, dvrr_stack+75, dvrr_stack+4323);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+4329, dvrr_stack+532, dvrr_stack+75);
 tmp = dvrr_stack + 4329;
 target_ptr = Libderiv->deriv_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+4389, dvrr_stack+210, dvrr_stack+4773);
 tmp = dvrr_stack + 4389;
 target_ptr = Libderiv->deriv_classes[2][2][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+4953, dvrr_stack+1340, dvrr_stack+210);

 /* compute (3 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+4425, dvrr_stack+4323, dvrr_stack+402, dvrr_stack+3, dvrr_stack+393, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,10,1,dvrr_stack+2748, dvrr_stack+2348, dvrr_stack+4425);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+4435, dvrr_stack+2913, dvrr_stack+2348);
 tmp = dvrr_stack + 4435;
 target_ptr = Libderiv->deriv_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+4535, dvrr_stack+2498, dvrr_stack+4083);
 tmp = dvrr_stack + 4535;
 target_ptr = Libderiv->deriv_classes[3][2][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+4595, dvrr_stack+4113, dvrr_stack+2498);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+4745, dvrr_stack+75, dvrr_stack+4323);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+3123, dvrr_stack+532, dvrr_stack+75);
 tmp = dvrr_stack + 3123;
 target_ptr = Libderiv->deriv_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+3183, dvrr_stack+210, dvrr_stack+4773);
 tmp = dvrr_stack + 3183;
 target_ptr = Libderiv->deriv_classes[2][2][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+3219, dvrr_stack+1340, dvrr_stack+210);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,10,1,dvrr_stack+3309, dvrr_stack+2348, dvrr_stack+4425);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+1596, dvrr_stack+2913, dvrr_stack+2348);
 tmp = dvrr_stack + 1596;
 target_ptr = Libderiv->deriv_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+658, dvrr_stack+2498, dvrr_stack+4083);
 tmp = dvrr_stack + 658;
 target_ptr = Libderiv->deriv_classes[3][2][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+6721, dvrr_stack+4113, dvrr_stack+2498);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+3339, dvrr_stack+75, dvrr_stack+4323);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+718, dvrr_stack+532, dvrr_stack+75);
 tmp = dvrr_stack + 718;
 target_ptr = Libderiv->deriv_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+1696, dvrr_stack+210, dvrr_stack+4773);
 tmp = dvrr_stack + 1696;
 target_ptr = Libderiv->deriv_classes[2][2][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+6871, dvrr_stack+1340, dvrr_stack+210);

 /* compute (3 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,10,1,dvrr_stack+4791, dvrr_stack+2348, dvrr_stack+4425);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+1301, dvrr_stack+2913, dvrr_stack+2348);
 tmp = dvrr_stack + 1301;
 target_ptr = Libderiv->deriv_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+1401, dvrr_stack+2498, dvrr_stack+4083);
 tmp = dvrr_stack + 1401;
 target_ptr = Libderiv->deriv_classes[3][2][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+6961, dvrr_stack+4113, dvrr_stack+2498);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+4113,dvrr_stack+180,dvrr_stack+57,3);


 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+4167, dvrr_stack+4113, NULL);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+4185,dvrr_stack+487,dvrr_stack+180,3);


 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+4275, dvrr_stack+4185, NULL);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+3357, dvrr_stack+1147, dvrr_stack+3, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+15, dvrr_stack+30, dvrr_stack+131, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+4305, dvrr_stack+1150, dvrr_stack+15, dvrr_stack+0, dvrr_stack+30, dvrr_stack+3);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+7111, dvrr_stack+2276, dvrr_stack+4305, dvrr_stack+6, dvrr_stack+1150, dvrr_stack+3357);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+778, dvrr_stack+134, dvrr_stack+411, NULL, NULL, dvrr_stack+131);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+7141, dvrr_stack+2294, dvrr_stack+778, dvrr_stack+33, dvrr_stack+134, dvrr_stack+15);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+7177, dvrr_stack+2312, dvrr_stack+7141, dvrr_stack+39, dvrr_stack+2294, dvrr_stack+4305);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+7237, dvrr_stack+2348, dvrr_stack+7177, dvrr_stack+75, dvrr_stack+2312, dvrr_stack+7111);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+7327, dvrr_stack+417, dvrr_stack+1168, NULL, NULL, dvrr_stack+411);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+7357, dvrr_stack+2408, dvrr_stack+7327, dvrr_stack+140, dvrr_stack+417, dvrr_stack+778);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+7417, dvrr_stack+2438, dvrr_stack+7357, dvrr_stack+150, dvrr_stack+2408, dvrr_stack+7141);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+7517, dvrr_stack+2498, dvrr_stack+7417, dvrr_stack+210, dvrr_stack+2438, dvrr_stack+7177);

 /* compute (4 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+7667,dvrr_stack+7517,dvrr_stack+7237,15);


 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+7937, dvrr_stack+7667, NULL);

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8027, dvrr_stack+1178, dvrr_stack+378, NULL, NULL, dvrr_stack+1168);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8072, dvrr_stack+2778, dvrr_stack+8027, dvrr_stack+427, dvrr_stack+1178, dvrr_stack+7327);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8162, dvrr_stack+2823, dvrr_stack+8072, dvrr_stack+442, dvrr_stack+2778, dvrr_stack+7357);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8312, dvrr_stack+2913, dvrr_stack+8162, dvrr_stack+532, dvrr_stack+2823, dvrr_stack+7417);

 /* compute (4 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+8537,dvrr_stack+8312,dvrr_stack+7517,15);


 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+2778, dvrr_stack+8537, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+532, dvrr_stack+4113, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+550, dvrr_stack+4185, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+2928, dvrr_stack+7667, NULL);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+8027, dvrr_stack+8537, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+580, dvrr_stack+4113, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+4113, dvrr_stack+4185, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+4185, dvrr_stack+7667, NULL);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+7667, dvrr_stack+8537, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+8537, dvrr_stack+180, dvrr_stack+1159);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+8555, dvrr_stack+487, dvrr_stack+57);

 /* compute (3 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+4425, dvrr_stack+402, dvrr_stack+3357, dvrr_stack+393, dvrr_stack+1147, NULL);

 /* compute (4 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+3018, dvrr_stack+4083, dvrr_stack+7111, dvrr_stack+4773, dvrr_stack+2276, dvrr_stack+4425);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,15,1,dvrr_stack+8585, dvrr_stack+7517, dvrr_stack+3018);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+8675, dvrr_stack+8312, dvrr_stack+7237);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+378, dvrr_stack+180, dvrr_stack+1159);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+4083, dvrr_stack+487, dvrr_stack+57);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,15,1,dvrr_stack+8825, dvrr_stack+7517, dvrr_stack+3018);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+8915, dvrr_stack+8312, dvrr_stack+7237);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+9065, dvrr_stack+180, dvrr_stack+1159);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+9083, dvrr_stack+487, dvrr_stack+57);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,15,1,dvrr_stack+9113, dvrr_stack+7517, dvrr_stack+3018);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+9203, dvrr_stack+8312, dvrr_stack+7237);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+3018, dvrr_stack+75, dvrr_stack+24);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+9353, dvrr_stack+7237, dvrr_stack+75);
 tmp = dvrr_stack + 9353;
 target_ptr = Libderiv->deriv_classes[3][2][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+9413, dvrr_stack+210, dvrr_stack+121);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+9443, dvrr_stack+7517, dvrr_stack+210);
 tmp = dvrr_stack + 9443;
 target_ptr = Libderiv->deriv_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+9543, dvrr_stack+2348, dvrr_stack+57);
 tmp = dvrr_stack + 9543;
 target_ptr = Libderiv->deriv_classes[2][2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+402, dvrr_stack+3, dvrr_stack+0, Data->F+3, Data->F+4, NULL);

 /* compute (3 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+4425, dvrr_stack+3357, dvrr_stack+402, dvrr_stack+1147, dvrr_stack+3, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1159, dvrr_stack+131, dvrr_stack+408, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+3036, dvrr_stack+15, dvrr_stack+1159, dvrr_stack+30, dvrr_stack+131, dvrr_stack+0);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+9579, dvrr_stack+4305, dvrr_stack+3036, dvrr_stack+1150, dvrr_stack+15, dvrr_stack+402);

 /* compute (4 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+9609, dvrr_stack+7111, dvrr_stack+9579, dvrr_stack+2276, dvrr_stack+4305, dvrr_stack+4425);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2276, dvrr_stack+411, dvrr_stack+396, NULL, NULL, dvrr_stack+408);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+9654, dvrr_stack+778, dvrr_stack+2276, dvrr_stack+134, dvrr_stack+411, dvrr_stack+1159);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+9690, dvrr_stack+7141, dvrr_stack+9654, dvrr_stack+2294, dvrr_stack+778, dvrr_stack+3036);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+9750, dvrr_stack+7177, dvrr_stack+9690, dvrr_stack+2312, dvrr_stack+7141, dvrr_stack+9579);

 /* compute (5 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+9840, dvrr_stack+7237, dvrr_stack+9750, dvrr_stack+2348, dvrr_stack+7177, dvrr_stack+9609);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,6,dvrr_stack+7111, dvrr_stack+9840, dvrr_stack+2348);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+9579, dvrr_stack+2498, dvrr_stack+180);
 tmp = dvrr_stack + 9579;
 target_ptr = Libderiv->deriv_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+2294, dvrr_stack+1168, dvrr_stack+111, NULL, NULL, dvrr_stack+396);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+7817, dvrr_stack+7327, dvrr_stack+2294, dvrr_stack+417, dvrr_stack+1168, dvrr_stack+2276);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+396, dvrr_stack+7357, dvrr_stack+7817, dvrr_stack+2408, dvrr_stack+7327, dvrr_stack+9654);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+8177, dvrr_stack+7417, dvrr_stack+396, dvrr_stack+2438, dvrr_stack+7357, dvrr_stack+9690);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+8327, dvrr_stack+7517, dvrr_stack+8177, dvrr_stack+2498, dvrr_stack+7417, dvrr_stack+9750);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+8177, dvrr_stack+8327, dvrr_stack+2498);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+396, dvrr_stack+75, dvrr_stack+24);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+414, dvrr_stack+7237, dvrr_stack+75);
 tmp = dvrr_stack + 414;
 target_ptr = Libderiv->deriv_classes[3][2][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+474, dvrr_stack+210, dvrr_stack+121);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+7327, dvrr_stack+7517, dvrr_stack+210);
 tmp = dvrr_stack + 7327;
 target_ptr = Libderiv->deriv_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+7201, dvrr_stack+2348, dvrr_stack+57);
 tmp = dvrr_stack + 7201;
 target_ptr = Libderiv->deriv_classes[2][2][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,6,dvrr_stack+7427, dvrr_stack+9840, dvrr_stack+2348);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+2408, dvrr_stack+2498, dvrr_stack+180);
 tmp = dvrr_stack + 2408;
 target_ptr = Libderiv->deriv_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+9639, dvrr_stack+8327, dvrr_stack+2498);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+504, dvrr_stack+75, dvrr_stack+24);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+7817, dvrr_stack+7237, dvrr_stack+75);
 tmp = dvrr_stack + 7817;
 target_ptr = Libderiv->deriv_classes[3][2][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+2468, dvrr_stack+210, dvrr_stack+121);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+75, dvrr_stack+7517, dvrr_stack+210);
 tmp = dvrr_stack + 75;
 target_ptr = Libderiv->deriv_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+210, dvrr_stack+2348, dvrr_stack+57);
 tmp = dvrr_stack + 210;
 target_ptr = Libderiv->deriv_classes[2][2][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,6,dvrr_stack+7237, dvrr_stack+9840, dvrr_stack+2348);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+7877, dvrr_stack+2498, dvrr_stack+180);
 tmp = dvrr_stack + 7877;
 target_ptr = Libderiv->deriv_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+7517, dvrr_stack+8327, dvrr_stack+2498);

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,36,dvrr_stack+2498, dvrr_stack+1018, NULL);
 tmp = dvrr_stack + 2498;
 target_ptr = Libderiv->deriv2_classes[2][2][143];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,60,dvrr_stack+2534, dvrr_stack+2096, NULL);
 tmp = dvrr_stack + 2534;
 target_ptr = Libderiv->deriv2_classes[2][3][143];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,60,dvrr_stack+8327, dvrr_stack+3723, NULL);
 tmp = dvrr_stack + 8327;
 target_ptr = Libderiv->deriv2_classes[3][2][143];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,100,dvrr_stack+8387, dvrr_stack+5373, NULL);
 tmp = dvrr_stack + 8387;
 target_ptr = Libderiv->deriv2_classes[3][3][143];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,36,dvrr_stack+8487, dvrr_stack+1018, NULL);
 tmp = dvrr_stack + 8487;
 target_ptr = Libderiv->deriv2_classes[2][2][131];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,60,dvrr_stack+2276, dvrr_stack+2096, NULL);
 tmp = dvrr_stack + 2276;
 target_ptr = Libderiv->deriv2_classes[2][3][131];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,60,dvrr_stack+2336, dvrr_stack+3723, NULL);
 tmp = dvrr_stack + 2336;
 target_ptr = Libderiv->deriv2_classes[3][2][131];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,100,dvrr_stack+9789, dvrr_stack+5373, NULL);
 tmp = dvrr_stack + 9789;
 target_ptr = Libderiv->deriv2_classes[3][3][131];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+9889, dvrr_stack+1193, NULL);
 tmp = dvrr_stack + 9889;
 target_ptr = Libderiv->deriv2_classes[2][2][130];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+9925, dvrr_stack+3903, NULL);
 tmp = dvrr_stack + 9925;
 target_ptr = Libderiv->deriv2_classes[2][3][130];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+0, dvrr_stack+5673, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[3][2][130];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+9985, dvrr_stack+5853, NULL);
 tmp = dvrr_stack + 9985;
 target_ptr = Libderiv->deriv2_classes[3][3][130];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,36,dvrr_stack+10085, dvrr_stack+1018, NULL);
 tmp = dvrr_stack + 10085;
 target_ptr = Libderiv->deriv2_classes[2][2][119];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,60,dvrr_stack+1018, dvrr_stack+2096, NULL);
 tmp = dvrr_stack + 1018;
 target_ptr = Libderiv->deriv2_classes[2][3][119];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,60,dvrr_stack+2096, dvrr_stack+3723, NULL);
 tmp = dvrr_stack + 2096;
 target_ptr = Libderiv->deriv2_classes[3][2][119];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,100,dvrr_stack+3723, dvrr_stack+5373, NULL);
 tmp = dvrr_stack + 3723;
 target_ptr = Libderiv->deriv2_classes[3][3][119];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+3823, dvrr_stack+1193, NULL);
 tmp = dvrr_stack + 3823;
 target_ptr = Libderiv->deriv2_classes[2][2][118];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+2156, dvrr_stack+3903, NULL);
 tmp = dvrr_stack + 2156;
 target_ptr = Libderiv->deriv2_classes[2][3][118];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+2216, dvrr_stack+5673, NULL);
 tmp = dvrr_stack + 2216;
 target_ptr = Libderiv->deriv2_classes[3][2][118];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+3859, dvrr_stack+5853, NULL);
 tmp = dvrr_stack + 3859;
 target_ptr = Libderiv->deriv2_classes[3][3][118];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+3959, dvrr_stack+6153, NULL);
 tmp = dvrr_stack + 3959;
 target_ptr = Libderiv->deriv2_classes[2][2][117];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+3995, dvrr_stack+802, NULL);
 tmp = dvrr_stack + 3995;
 target_ptr = Libderiv->deriv2_classes[2][3][117];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+1078, dvrr_stack+1736, NULL);
 tmp = dvrr_stack + 1078;
 target_ptr = Libderiv->deriv2_classes[3][2][117];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+1138, dvrr_stack+3363, NULL);
 tmp = dvrr_stack + 1138;
 target_ptr = Libderiv->deriv2_classes[3][3][117];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+1238, dvrr_stack+3663, dvrr_stack+4845);
 tmp = dvrr_stack + 1238;
 target_ptr = Libderiv->deriv2_classes[2][2][107];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+10121, dvrr_stack+4863, dvrr_stack+982);
 tmp = dvrr_stack + 10121;
 target_ptr = Libderiv->deriv2_classes[2][3][107];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+10181, dvrr_stack+5073, dvrr_stack+5043);
 tmp = dvrr_stack + 10181;
 target_ptr = Libderiv->deriv2_classes[3][2][107];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+10241, dvrr_stack+1916, dvrr_stack+5173);
 tmp = dvrr_stack + 10241;
 target_ptr = Libderiv->deriv2_classes[3][3][107];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+10341, dvrr_stack+5251, dvrr_stack+5233);
 tmp = dvrr_stack + 10341;
 target_ptr = Libderiv->deriv2_classes[2][2][106];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+10377, dvrr_stack+6261, dvrr_stack+5311);
 tmp = dvrr_stack + 10377;
 target_ptr = Libderiv->deriv2_classes[2][3][106];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+10437, dvrr_stack+6351, dvrr_stack+2066);
 tmp = dvrr_stack + 10437;
 target_ptr = Libderiv->deriv2_classes[3][2][106];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+10497, dvrr_stack+6511, dvrr_stack+6451);
 tmp = dvrr_stack + 10497;
 target_ptr = Libderiv->deriv2_classes[3][3][106];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+10597, dvrr_stack+6661, dvrr_stack+5347);
 tmp = dvrr_stack + 10597;
 target_ptr = Libderiv->deriv2_classes[2][2][105];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+10633, dvrr_stack+270, dvrr_stack+622);
 tmp = dvrr_stack + 10633;
 target_ptr = Libderiv->deriv2_classes[2][3][105];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+10693, dvrr_stack+1496, dvrr_stack+1466);
 tmp = dvrr_stack + 10693;
 target_ptr = Libderiv->deriv2_classes[3][2][105];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+10753, dvrr_stack+2598, dvrr_stack+3063);
 tmp = dvrr_stack + 10753;
 target_ptr = Libderiv->deriv2_classes[3][3][105];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+10853, dvrr_stack+4329, dvrr_stack+360);
 tmp = dvrr_stack + 10853;
 target_ptr = Libderiv->deriv2_classes[2][2][104];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+10889, dvrr_stack+4953, dvrr_stack+4389);
 tmp = dvrr_stack + 10889;
 target_ptr = Libderiv->deriv2_classes[2][3][104];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+778, dvrr_stack+4435, dvrr_stack+2748);
 tmp = dvrr_stack + 778;
 target_ptr = Libderiv->deriv2_classes[3][2][104];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+838, dvrr_stack+4595, dvrr_stack+4535);
 tmp = dvrr_stack + 838;
 target_ptr = Libderiv->deriv2_classes[3][3][104];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+10949, dvrr_stack+3663, dvrr_stack+4845);
 tmp = dvrr_stack + 10949;
 target_ptr = Libderiv->deriv2_classes[2][2][95];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+3357, dvrr_stack+4863, dvrr_stack+982);
 tmp = dvrr_stack + 3357;
 target_ptr = Libderiv->deriv2_classes[2][3][95];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+3417, dvrr_stack+5073, dvrr_stack+5043);
 tmp = dvrr_stack + 3417;
 target_ptr = Libderiv->deriv2_classes[3][2][95];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+3477, dvrr_stack+1916, dvrr_stack+5173);
 tmp = dvrr_stack + 3477;
 target_ptr = Libderiv->deriv2_classes[3][3][95];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+938, dvrr_stack+5251, dvrr_stack+5233);
 tmp = dvrr_stack + 938;
 target_ptr = Libderiv->deriv2_classes[2][2][94];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+3577, dvrr_stack+6261, dvrr_stack+5311);
 tmp = dvrr_stack + 3577;
 target_ptr = Libderiv->deriv2_classes[2][3][94];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+1732, dvrr_stack+6351, dvrr_stack+2066);
 tmp = dvrr_stack + 1732;
 target_ptr = Libderiv->deriv2_classes[3][2][94];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+1792, dvrr_stack+6511, dvrr_stack+6451);
 tmp = dvrr_stack + 1792;
 target_ptr = Libderiv->deriv2_classes[3][3][94];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+5365, dvrr_stack+6661, dvrr_stack+5347);
 tmp = dvrr_stack + 5365;
 target_ptr = Libderiv->deriv2_classes[2][2][93];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+5401, dvrr_stack+270, dvrr_stack+622);
 tmp = dvrr_stack + 5401;
 target_ptr = Libderiv->deriv2_classes[2][3][93];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+5461, dvrr_stack+1496, dvrr_stack+1466);
 tmp = dvrr_stack + 5461;
 target_ptr = Libderiv->deriv2_classes[3][2][93];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+5521, dvrr_stack+2598, dvrr_stack+3063);
 tmp = dvrr_stack + 5521;
 target_ptr = Libderiv->deriv2_classes[3][3][93];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+5621, dvrr_stack+4329, dvrr_stack+360);
 tmp = dvrr_stack + 5621;
 target_ptr = Libderiv->deriv2_classes[2][2][92];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+5657, dvrr_stack+4953, dvrr_stack+4389);
 tmp = dvrr_stack + 5657;
 target_ptr = Libderiv->deriv2_classes[2][3][92];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+5717, dvrr_stack+4435, dvrr_stack+2748);
 tmp = dvrr_stack + 5717;
 target_ptr = Libderiv->deriv2_classes[3][2][92];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+5777, dvrr_stack+4595, dvrr_stack+4535);
 tmp = dvrr_stack + 5777;
 target_ptr = Libderiv->deriv2_classes[3][3][92];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+5877, dvrr_stack+3123, dvrr_stack+4745);
 tmp = dvrr_stack + 5877;
 target_ptr = Libderiv->deriv2_classes[2][2][91];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+5913, dvrr_stack+3219, dvrr_stack+3183);
 tmp = dvrr_stack + 5913;
 target_ptr = Libderiv->deriv2_classes[2][3][91];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+5973, dvrr_stack+1596, dvrr_stack+3309);
 tmp = dvrr_stack + 5973;
 target_ptr = Libderiv->deriv2_classes[3][2][91];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+6033, dvrr_stack+6721, dvrr_stack+658);
 tmp = dvrr_stack + 6033;
 target_ptr = Libderiv->deriv2_classes[3][3][91];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+6133, dvrr_stack+3663, dvrr_stack+4845);
 tmp = dvrr_stack + 6133;
 target_ptr = Libderiv->deriv2_classes[2][2][83];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+6169, dvrr_stack+4863, dvrr_stack+982);
 tmp = dvrr_stack + 6169;
 target_ptr = Libderiv->deriv2_classes[2][3][83];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+4821, dvrr_stack+5073, dvrr_stack+5043);
 tmp = dvrr_stack + 4821;
 target_ptr = Libderiv->deriv2_classes[3][2][83];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+10985, dvrr_stack+1916, dvrr_stack+5173);
 tmp = dvrr_stack + 10985;
 target_ptr = Libderiv->deriv2_classes[3][3][83];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+11085, dvrr_stack+5251, dvrr_stack+5233);
 tmp = dvrr_stack + 11085;
 target_ptr = Libderiv->deriv2_classes[2][2][82];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+11121, dvrr_stack+6261, dvrr_stack+5311);
 tmp = dvrr_stack + 11121;
 target_ptr = Libderiv->deriv2_classes[2][3][82];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+11181, dvrr_stack+6351, dvrr_stack+2066);
 tmp = dvrr_stack + 11181;
 target_ptr = Libderiv->deriv2_classes[3][2][82];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+11241, dvrr_stack+6511, dvrr_stack+6451);
 tmp = dvrr_stack + 11241;
 target_ptr = Libderiv->deriv2_classes[3][3][82];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+6511, dvrr_stack+6661, dvrr_stack+5347);
 tmp = dvrr_stack + 6511;
 target_ptr = Libderiv->deriv2_classes[2][2][81];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+6547, dvrr_stack+270, dvrr_stack+622);
 tmp = dvrr_stack + 6547;
 target_ptr = Libderiv->deriv2_classes[2][3][81];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+246, dvrr_stack+1496, dvrr_stack+1466);
 tmp = dvrr_stack + 246;
 target_ptr = Libderiv->deriv2_classes[3][2][81];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+11341, dvrr_stack+2598, dvrr_stack+3063);
 tmp = dvrr_stack + 11341;
 target_ptr = Libderiv->deriv2_classes[3][3][81];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+6607, dvrr_stack+4329, dvrr_stack+360);
 tmp = dvrr_stack + 6607;
 target_ptr = Libderiv->deriv2_classes[2][2][80];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+2594, dvrr_stack+4953, dvrr_stack+4389);
 tmp = dvrr_stack + 2594;
 target_ptr = Libderiv->deriv2_classes[2][3][80];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+2654, dvrr_stack+4435, dvrr_stack+2748);
 tmp = dvrr_stack + 2654;
 target_ptr = Libderiv->deriv2_classes[3][2][80];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+11441, dvrr_stack+4595, dvrr_stack+4535);
 tmp = dvrr_stack + 11441;
 target_ptr = Libderiv->deriv2_classes[3][3][80];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+4595, dvrr_stack+3123, dvrr_stack+4745);
 tmp = dvrr_stack + 4595;
 target_ptr = Libderiv->deriv2_classes[2][2][79];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+4631, dvrr_stack+3219, dvrr_stack+3183);
 tmp = dvrr_stack + 4631;
 target_ptr = Libderiv->deriv2_classes[2][3][79];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+3219, dvrr_stack+1596, dvrr_stack+3309);
 tmp = dvrr_stack + 3219;
 target_ptr = Libderiv->deriv2_classes[3][2][79];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+4691, dvrr_stack+6721, dvrr_stack+658);
 tmp = dvrr_stack + 4691;
 target_ptr = Libderiv->deriv2_classes[3][3][79];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+6721, dvrr_stack+718, dvrr_stack+3339);
 tmp = dvrr_stack + 6721;
 target_ptr = Libderiv->deriv2_classes[2][2][78];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+6757, dvrr_stack+6871, dvrr_stack+1696);
 tmp = dvrr_stack + 6757;
 target_ptr = Libderiv->deriv2_classes[2][3][78];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+6817, dvrr_stack+1301, dvrr_stack+4791);
 tmp = dvrr_stack + 6817;
 target_ptr = Libderiv->deriv2_classes[3][2][78];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+11541, dvrr_stack+6961, dvrr_stack+1401);
 tmp = dvrr_stack + 11541;
 target_ptr = Libderiv->deriv2_classes[3][3][78];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,6,dvrr_stack+6877, dvrr_stack+5173, dvrr_stack+4167);
 tmp = dvrr_stack + 6877;
 target_ptr = Libderiv->deriv2_classes[2][2][35];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6913, dvrr_stack+5073, dvrr_stack+4275);
 tmp = dvrr_stack + 6913;
 target_ptr = Libderiv->deriv2_classes[2][3][35];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,6,dvrr_stack+6973, dvrr_stack+7937, dvrr_stack+982);
 tmp = dvrr_stack + 6973;
 target_ptr = Libderiv->deriv2_classes[3][2][35];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,10,dvrr_stack+11641, dvrr_stack+2778, dvrr_stack+3663);
 tmp = dvrr_stack + 11641;
 target_ptr = Libderiv->deriv2_classes[3][3][35];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+7033, dvrr_stack+6451, dvrr_stack+532);
 tmp = dvrr_stack + 7033;
 target_ptr = Libderiv->deriv2_classes[2][2][34];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+3279, dvrr_stack+6351, dvrr_stack+550);
 tmp = dvrr_stack + 3279;
 target_ptr = Libderiv->deriv2_classes[2][3][34];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+2714, dvrr_stack+2928, dvrr_stack+5311);
 tmp = dvrr_stack + 2714;
 target_ptr = Libderiv->deriv2_classes[3][2][34];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+11741, dvrr_stack+8027, dvrr_stack+5251);
 tmp = dvrr_stack + 11741;
 target_ptr = Libderiv->deriv2_classes[3][3][34];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+7069, dvrr_stack+3063, dvrr_stack+580);
 tmp = dvrr_stack + 7069;
 target_ptr = Libderiv->deriv2_classes[2][2][33];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+306, dvrr_stack+1496, dvrr_stack+4113);
 tmp = dvrr_stack + 306;
 target_ptr = Libderiv->deriv2_classes[2][3][33];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+11841, dvrr_stack+4185, dvrr_stack+622);
 tmp = dvrr_stack + 11841;
 target_ptr = Libderiv->deriv2_classes[3][2][33];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+11901, dvrr_stack+7667, dvrr_stack+6661);
 tmp = dvrr_stack + 11901;
 target_ptr = Libderiv->deriv2_classes[3][3][33];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+12001, dvrr_stack+4535, dvrr_stack+8537);
 tmp = dvrr_stack + 12001;
 target_ptr = Libderiv->deriv2_classes[2][2][32];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+12037, dvrr_stack+4435, dvrr_stack+8555);
 tmp = dvrr_stack + 12037;
 target_ptr = Libderiv->deriv2_classes[2][3][32];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+12097, dvrr_stack+8585, dvrr_stack+4389);
 tmp = dvrr_stack + 12097;
 target_ptr = Libderiv->deriv2_classes[3][2][32];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+12157, dvrr_stack+8675, dvrr_stack+4329);
 tmp = dvrr_stack + 12157;
 target_ptr = Libderiv->deriv2_classes[3][3][32];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+12257, dvrr_stack+658, dvrr_stack+378);
 tmp = dvrr_stack + 12257;
 target_ptr = Libderiv->deriv2_classes[2][2][31];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+12293, dvrr_stack+1596, dvrr_stack+4083);
 tmp = dvrr_stack + 12293;
 target_ptr = Libderiv->deriv2_classes[2][3][31];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+12353, dvrr_stack+8825, dvrr_stack+3183);
 tmp = dvrr_stack + 12353;
 target_ptr = Libderiv->deriv2_classes[3][2][31];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+1892, dvrr_stack+8915, dvrr_stack+3123);
 tmp = dvrr_stack + 1892;
 target_ptr = Libderiv->deriv2_classes[3][3][31];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+12413, dvrr_stack+1401, dvrr_stack+9065);
 tmp = dvrr_stack + 12413;
 target_ptr = Libderiv->deriv2_classes[2][2][30];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+1992, dvrr_stack+1301, dvrr_stack+9083);
 tmp = dvrr_stack + 1992;
 target_ptr = Libderiv->deriv2_classes[2][3][30];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+4881, dvrr_stack+9113, dvrr_stack+1696);
 tmp = dvrr_stack + 4881;
 target_ptr = Libderiv->deriv2_classes[3][2][30];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+4941, dvrr_stack+9203, dvrr_stack+718);
 tmp = dvrr_stack + 4941;
 target_ptr = Libderiv->deriv2_classes[3][3][30];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+12449, dvrr_stack+9353, dvrr_stack+3018);
 tmp = dvrr_stack + 12449;
 target_ptr = Libderiv->deriv2_classes[2][2][26];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6229, dvrr_stack+9443, dvrr_stack+9413);
 tmp = dvrr_stack + 6229;
 target_ptr = Libderiv->deriv2_classes[2][3][26];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,6,dvrr_stack+6289, dvrr_stack+7111, dvrr_stack+9543);
 tmp = dvrr_stack + 6289;
 target_ptr = Libderiv->deriv2_classes[3][2][26];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+12485, dvrr_stack+8177, dvrr_stack+9579);
 tmp = dvrr_stack + 12485;
 target_ptr = Libderiv->deriv2_classes[3][3][26];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,6,dvrr_stack+12585, dvrr_stack+5173, dvrr_stack+4167);
 tmp = dvrr_stack + 12585;
 target_ptr = Libderiv->deriv2_classes[2][2][23];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,10,dvrr_stack+12621, dvrr_stack+5073, dvrr_stack+4275);
 tmp = dvrr_stack + 12621;
 target_ptr = Libderiv->deriv2_classes[2][3][23];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,6,dvrr_stack+12681, dvrr_stack+7937, dvrr_stack+982);
 tmp = dvrr_stack + 12681;
 target_ptr = Libderiv->deriv2_classes[3][2][23];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,10,dvrr_stack+12741, dvrr_stack+2778, dvrr_stack+3663);
 tmp = dvrr_stack + 12741;
 target_ptr = Libderiv->deriv2_classes[3][3][23];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+12841, dvrr_stack+6451, dvrr_stack+532);
 tmp = dvrr_stack + 12841;
 target_ptr = Libderiv->deriv2_classes[2][2][22];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+12877, dvrr_stack+6351, dvrr_stack+550);
 tmp = dvrr_stack + 12877;
 target_ptr = Libderiv->deriv2_classes[2][3][22];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+12937, dvrr_stack+2928, dvrr_stack+5311);
 tmp = dvrr_stack + 12937;
 target_ptr = Libderiv->deriv2_classes[3][2][22];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+12997, dvrr_stack+8027, dvrr_stack+5251);
 tmp = dvrr_stack + 12997;
 target_ptr = Libderiv->deriv2_classes[3][3][22];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+13097, dvrr_stack+3063, dvrr_stack+580);
 tmp = dvrr_stack + 13097;
 target_ptr = Libderiv->deriv2_classes[2][2][21];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+13133, dvrr_stack+1496, dvrr_stack+4113);
 tmp = dvrr_stack + 13133;
 target_ptr = Libderiv->deriv2_classes[2][3][21];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+13193, dvrr_stack+4185, dvrr_stack+622);
 tmp = dvrr_stack + 13193;
 target_ptr = Libderiv->deriv2_classes[3][2][21];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+13253, dvrr_stack+7667, dvrr_stack+6661);
 tmp = dvrr_stack + 13253;
 target_ptr = Libderiv->deriv2_classes[3][3][21];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+13353, dvrr_stack+4535, dvrr_stack+8537);
 tmp = dvrr_stack + 13353;
 target_ptr = Libderiv->deriv2_classes[2][2][20];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+13389, dvrr_stack+4435, dvrr_stack+8555);
 tmp = dvrr_stack + 13389;
 target_ptr = Libderiv->deriv2_classes[2][3][20];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+13449, dvrr_stack+8585, dvrr_stack+4389);
 tmp = dvrr_stack + 13449;
 target_ptr = Libderiv->deriv2_classes[3][2][20];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+13509, dvrr_stack+8675, dvrr_stack+4329);
 tmp = dvrr_stack + 13509;
 target_ptr = Libderiv->deriv2_classes[3][3][20];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+13609, dvrr_stack+658, dvrr_stack+378);
 tmp = dvrr_stack + 13609;
 target_ptr = Libderiv->deriv2_classes[2][2][19];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+13645, dvrr_stack+1596, dvrr_stack+4083);
 tmp = dvrr_stack + 13645;
 target_ptr = Libderiv->deriv2_classes[2][3][19];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+13705, dvrr_stack+8825, dvrr_stack+3183);
 tmp = dvrr_stack + 13705;
 target_ptr = Libderiv->deriv2_classes[3][2][19];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+13765, dvrr_stack+8915, dvrr_stack+3123);
 tmp = dvrr_stack + 13765;
 target_ptr = Libderiv->deriv2_classes[3][3][19];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+13865, dvrr_stack+1401, dvrr_stack+9065);
 tmp = dvrr_stack + 13865;
 target_ptr = Libderiv->deriv2_classes[2][2][18];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+13901, dvrr_stack+1301, dvrr_stack+9083);
 tmp = dvrr_stack + 13901;
 target_ptr = Libderiv->deriv2_classes[2][3][18];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+13961, dvrr_stack+9113, dvrr_stack+1696);
 tmp = dvrr_stack + 13961;
 target_ptr = Libderiv->deriv2_classes[3][2][18];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+14021, dvrr_stack+9203, dvrr_stack+718);
 tmp = dvrr_stack + 14021;
 target_ptr = Libderiv->deriv2_classes[3][3][18];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+14121, dvrr_stack+9353, dvrr_stack+3018);
 tmp = dvrr_stack + 14121;
 target_ptr = Libderiv->deriv2_classes[2][2][14];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+14157, dvrr_stack+9443, dvrr_stack+9413);
 tmp = dvrr_stack + 14157;
 target_ptr = Libderiv->deriv2_classes[2][3][14];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+14217, dvrr_stack+7111, dvrr_stack+9543);
 tmp = dvrr_stack + 14217;
 target_ptr = Libderiv->deriv2_classes[3][2][14];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+14277, dvrr_stack+8177, dvrr_stack+9579);
 tmp = dvrr_stack + 14277;
 target_ptr = Libderiv->deriv2_classes[3][3][14];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+14377, dvrr_stack+414, dvrr_stack+396);
 tmp = dvrr_stack + 14377;
 target_ptr = Libderiv->deriv2_classes[2][2][13];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+14413, dvrr_stack+7327, dvrr_stack+474);
 tmp = dvrr_stack + 14413;
 target_ptr = Libderiv->deriv2_classes[2][3][13];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,6,dvrr_stack+14473, dvrr_stack+7427, dvrr_stack+7201);
 tmp = dvrr_stack + 14473;
 target_ptr = Libderiv->deriv2_classes[3][2][13];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+14533, dvrr_stack+9639, dvrr_stack+2408);
 tmp = dvrr_stack + 14533;
 target_ptr = Libderiv->deriv2_classes[3][3][13];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,6,dvrr_stack+14633, dvrr_stack+5173, dvrr_stack+4167);
 tmp = dvrr_stack + 14633;
 target_ptr = Libderiv->deriv2_classes[2][2][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,10,dvrr_stack+5173, dvrr_stack+5073, dvrr_stack+4275);
 tmp = dvrr_stack + 5173;
 target_ptr = Libderiv->deriv2_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,6,dvrr_stack+14669, dvrr_stack+7937, dvrr_stack+982);
 tmp = dvrr_stack + 14669;
 target_ptr = Libderiv->deriv2_classes[3][2][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,10,dvrr_stack+14729, dvrr_stack+2778, dvrr_stack+3663);
 tmp = dvrr_stack + 14729;
 target_ptr = Libderiv->deriv2_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+7937, dvrr_stack+6451, dvrr_stack+532);
 tmp = dvrr_stack + 7937;
 target_ptr = Libderiv->deriv2_classes[2][2][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+6451, dvrr_stack+6351, dvrr_stack+550);
 tmp = dvrr_stack + 6451;
 target_ptr = Libderiv->deriv2_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+2774, dvrr_stack+2928, dvrr_stack+5311);
 tmp = dvrr_stack + 2774;
 target_ptr = Libderiv->deriv2_classes[3][2][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+2834, dvrr_stack+8027, dvrr_stack+5251);
 tmp = dvrr_stack + 2834;
 target_ptr = Libderiv->deriv2_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+7973, dvrr_stack+3063, dvrr_stack+580);
 tmp = dvrr_stack + 7973;
 target_ptr = Libderiv->deriv2_classes[2][2][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+8009, dvrr_stack+1496, dvrr_stack+4113);
 tmp = dvrr_stack + 8009;
 target_ptr = Libderiv->deriv2_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+4113, dvrr_stack+4185, dvrr_stack+622);
 tmp = dvrr_stack + 4113;
 target_ptr = Libderiv->deriv2_classes[3][2][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+4173, dvrr_stack+7667, dvrr_stack+6661);
 tmp = dvrr_stack + 4173;
 target_ptr = Libderiv->deriv2_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+7667, dvrr_stack+4535, dvrr_stack+8537);
 tmp = dvrr_stack + 7667;
 target_ptr = Libderiv->deriv2_classes[2][2][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+4535, dvrr_stack+4435, dvrr_stack+8555);
 tmp = dvrr_stack + 4535;
 target_ptr = Libderiv->deriv2_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+7703, dvrr_stack+8585, dvrr_stack+4389);
 tmp = dvrr_stack + 7703;
 target_ptr = Libderiv->deriv2_classes[3][2][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+4389, dvrr_stack+8675, dvrr_stack+4329);
 tmp = dvrr_stack + 4389;
 target_ptr = Libderiv->deriv2_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+4489, dvrr_stack+658, dvrr_stack+378);
 tmp = dvrr_stack + 4489;
 target_ptr = Libderiv->deriv2_classes[2][2][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+4273, dvrr_stack+1596, dvrr_stack+4083);
 tmp = dvrr_stack + 4273;
 target_ptr = Libderiv->deriv2_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+8069, dvrr_stack+8825, dvrr_stack+3183);
 tmp = dvrr_stack + 8069;
 target_ptr = Libderiv->deriv2_classes[3][2][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+5233, dvrr_stack+8915, dvrr_stack+3123);
 tmp = dvrr_stack + 5233;
 target_ptr = Libderiv->deriv2_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+7763, dvrr_stack+1401, dvrr_stack+9065);
 tmp = dvrr_stack + 7763;
 target_ptr = Libderiv->deriv2_classes[2][2][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+1401, dvrr_stack+1301, dvrr_stack+9083);
 tmp = dvrr_stack + 1401;
 target_ptr = Libderiv->deriv2_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+1461, dvrr_stack+9113, dvrr_stack+1696);
 tmp = dvrr_stack + 1461;
 target_ptr = Libderiv->deriv2_classes[3][2][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+1521, dvrr_stack+9203, dvrr_stack+718);
 tmp = dvrr_stack + 1521;
 target_ptr = Libderiv->deriv2_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+1621, dvrr_stack+9353, dvrr_stack+3018);
 tmp = dvrr_stack + 1621;
 target_ptr = Libderiv->deriv2_classes[2][2][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+1657, dvrr_stack+9443, dvrr_stack+9413);
 tmp = dvrr_stack + 1657;
 target_ptr = Libderiv->deriv2_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+6643, dvrr_stack+7111, dvrr_stack+9543);
 tmp = dvrr_stack + 6643;
 target_ptr = Libderiv->deriv2_classes[3][2][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+1274, dvrr_stack+8177, dvrr_stack+9579);
 tmp = dvrr_stack + 1274;
 target_ptr = Libderiv->deriv2_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+4333, dvrr_stack+414, dvrr_stack+396);
 tmp = dvrr_stack + 4333;
 target_ptr = Libderiv->deriv2_classes[2][2][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+8129, dvrr_stack+7327, dvrr_stack+474);
 tmp = dvrr_stack + 8129;
 target_ptr = Libderiv->deriv2_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+7327, dvrr_stack+7427, dvrr_stack+7201);
 tmp = dvrr_stack + 7327;
 target_ptr = Libderiv->deriv2_classes[3][2][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+7387, dvrr_stack+9639, dvrr_stack+2408);
 tmp = dvrr_stack + 7387;
 target_ptr = Libderiv->deriv2_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+8189, dvrr_stack+7817, dvrr_stack+504);
 tmp = dvrr_stack + 8189;
 target_ptr = Libderiv->deriv2_classes[2][2][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+7799, dvrr_stack+75, dvrr_stack+2468);
 tmp = dvrr_stack + 7799;
 target_ptr = Libderiv->deriv2_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,6,dvrr_stack+8225, dvrr_stack+7237, dvrr_stack+210);
 tmp = dvrr_stack + 8225;
 target_ptr = Libderiv->deriv2_classes[3][2][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+7105, dvrr_stack+7517, dvrr_stack+7877);
 tmp = dvrr_stack + 7105;
 target_ptr = Libderiv->deriv2_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];


}

