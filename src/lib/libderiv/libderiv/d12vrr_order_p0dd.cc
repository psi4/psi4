#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (p0|dd) integrals */

void d12vrr_order_p0dd(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+12, dvrr_stack+3, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+21, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+24, dvrr_stack+0, dvrr_stack+21, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+30, dvrr_stack+6, dvrr_stack+24, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+15, dvrr_stack+6, dvrr_stack+12, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+50, dvrr_stack+40, dvrr_stack+30, dvrr_stack+15, dvrr_stack+6, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+65, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+68, dvrr_stack+21, dvrr_stack+65, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+74, dvrr_stack+24, dvrr_stack+68, dvrr_stack+0, dvrr_stack+21, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+84, dvrr_stack+30, dvrr_stack+74, NULL, NULL, dvrr_stack+24);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+114, dvrr_stack+30, dvrr_stack+74, dvrr_stack+6, dvrr_stack+24, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+129, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+132, dvrr_stack+65, dvrr_stack+129, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+138, dvrr_stack+68, dvrr_stack+132, dvrr_stack+21, dvrr_stack+65, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+148, dvrr_stack+74, dvrr_stack+138, dvrr_stack+24, dvrr_stack+68, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+163, dvrr_stack+114, dvrr_stack+148, NULL, NULL, dvrr_stack+74);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+208, dvrr_stack+50, dvrr_stack+114, NULL, NULL, dvrr_stack+30);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+253, dvrr_stack+208, dvrr_stack+163, dvrr_stack+50, dvrr_stack+114, dvrr_stack+84);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+343, dvrr_stack+15, dvrr_stack+6, NULL, NULL, dvrr_stack+3);
 tmp = dvrr_stack + 343;
 target_ptr = Libderiv->dvrr_classes[1][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+361, dvrr_stack+40, dvrr_stack+30, NULL, NULL, dvrr_stack+6);
 tmp = dvrr_stack + 361;
 target_ptr = Libderiv->dvrr_classes[1][3];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+391,dvrr_stack+361,dvrr_stack+343,3);


 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+445,dvrr_stack+208,dvrr_stack+361,3);


 /* compute (1 0 | 2 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dd(Libderiv->CD,dvrr_stack+535,dvrr_stack+445,dvrr_stack+391,3);


 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,18,dvrr_stack+643, dvrr_stack+535, dvrr_stack+343);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+697, dvrr_stack+114, dvrr_stack+148, dvrr_stack+30, dvrr_stack+74, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+718, dvrr_stack+50, dvrr_stack+114, dvrr_stack+40, dvrr_stack+30, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+739, dvrr_stack+718, dvrr_stack+697, NULL, NULL, dvrr_stack+114);

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+802,dvrr_stack+739,dvrr_stack+208,3);


 /* compute (1 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+937,dvrr_stack+802,dvrr_stack+445,3);


 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,30,dvrr_stack+1117, dvrr_stack+937, dvrr_stack+361);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1207, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1210, dvrr_stack+129, dvrr_stack+1207, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1216, dvrr_stack+132, dvrr_stack+1210, dvrr_stack+65, dvrr_stack+129, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1226, dvrr_stack+138, dvrr_stack+1216, dvrr_stack+68, dvrr_stack+132, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1241, dvrr_stack+148, dvrr_stack+1226, dvrr_stack+74, dvrr_stack+138, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1262, dvrr_stack+697, dvrr_stack+1241, dvrr_stack+114, dvrr_stack+148, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1290, dvrr_stack+718, dvrr_stack+697, dvrr_stack+50, dvrr_stack+114, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1318, dvrr_stack+1290, dvrr_stack+1262, NULL, NULL, dvrr_stack+697);

 /* compute (1 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+1402,dvrr_stack+1318,dvrr_stack+739,3);


 /* compute (1 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+1591,dvrr_stack+1402,dvrr_stack+802,3);


 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,45,dvrr_stack+1861, dvrr_stack+1591, dvrr_stack+208);

 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,18,dvrr_stack+1262, dvrr_stack+535, dvrr_stack+343);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,30,dvrr_stack+1996, dvrr_stack+937, dvrr_stack+361);

 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,45,dvrr_stack+2086, dvrr_stack+1591, dvrr_stack+208);

 /* compute (1 0 | 2 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,18,dvrr_stack+2221, dvrr_stack+535, dvrr_stack+343);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,30,dvrr_stack+535, dvrr_stack+937, dvrr_stack+361);

 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,45,dvrr_stack+937, dvrr_stack+1591, dvrr_stack+208);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+129, dvrr_stack+12, dvrr_stack+3, NULL, NULL, Data->F+1);

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+1591,dvrr_stack+343,dvrr_stack+129,3);


 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+1618, dvrr_stack+1591, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+1627, dvrr_stack+445, NULL);
 tmp = dvrr_stack + 1627;
 target_ptr = Libderiv->deriv_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+625, dvrr_stack+391, NULL);
 tmp = dvrr_stack + 625;
 target_ptr = Libderiv->deriv_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+1072, dvrr_stack+802, NULL);
 tmp = dvrr_stack + 1072;
 target_ptr = Libderiv->deriv_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,63,dvrr_stack+1657, dvrr_stack+1402, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+1720, dvrr_stack+1591, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+1729, dvrr_stack+445, NULL);
 tmp = dvrr_stack + 1729;
 target_ptr = Libderiv->deriv_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+1759, dvrr_stack+391, NULL);
 tmp = dvrr_stack + 1759;
 target_ptr = Libderiv->deriv_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+1777, dvrr_stack+802, NULL);
 tmp = dvrr_stack + 1777;
 target_ptr = Libderiv->deriv_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+2275, dvrr_stack+1402, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+1822, dvrr_stack+1591, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+1831, dvrr_stack+445, NULL);
 tmp = dvrr_stack + 1831;
 target_ptr = Libderiv->deriv_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+445, dvrr_stack+391, NULL);
 tmp = dvrr_stack + 445;
 target_ptr = Libderiv->deriv_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+391, dvrr_stack+802, NULL);
 tmp = dvrr_stack + 391;
 target_ptr = Libderiv->deriv_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+802, dvrr_stack+1402, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+65, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+436, dvrr_stack+343, dvrr_stack+65);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+1402, dvrr_stack+208, dvrr_stack+343);
 tmp = dvrr_stack + 1402;
 target_ptr = Libderiv->deriv_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+1432, dvrr_stack+361, dvrr_stack+129);
 tmp = dvrr_stack + 1432;
 target_ptr = Libderiv->deriv_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+1450, dvrr_stack+739, dvrr_stack+361);
 tmp = dvrr_stack + 1450;
 target_ptr = Libderiv->deriv_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+1495, dvrr_stack+1318, dvrr_stack+208);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+1558, dvrr_stack+343, dvrr_stack+65);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1567, dvrr_stack+208, dvrr_stack+343);
 tmp = dvrr_stack + 1567;
 target_ptr = Libderiv->deriv_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+1597, dvrr_stack+361, dvrr_stack+129);
 tmp = dvrr_stack + 1597;
 target_ptr = Libderiv->deriv_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+865, dvrr_stack+739, dvrr_stack+361);
 tmp = dvrr_stack + 865;
 target_ptr = Libderiv->deriv_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+463, dvrr_stack+1318, dvrr_stack+208);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+526, dvrr_stack+343, dvrr_stack+65);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2338, dvrr_stack+208, dvrr_stack+343);
 tmp = dvrr_stack + 2338;
 target_ptr = Libderiv->deriv_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+910, dvrr_stack+361, dvrr_stack+129);
 tmp = dvrr_stack + 910;
 target_ptr = Libderiv->deriv_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+2368, dvrr_stack+739, dvrr_stack+361);
 tmp = dvrr_stack + 2368;
 target_ptr = Libderiv->deriv_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+2413, dvrr_stack+1318, dvrr_stack+208);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+1316,dvrr_stack+40,dvrr_stack+15,1);


 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+928, dvrr_stack+1316, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1334, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1343, dvrr_stack+6, dvrr_stack+24, NULL, NULL, dvrr_stack+0);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+1361, dvrr_stack+343, dvrr_stack+1343, dvrr_stack+15, dvrr_stack+6, dvrr_stack+1334);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+2476, dvrr_stack+361, dvrr_stack+84, dvrr_stack+40, dvrr_stack+30, dvrr_stack+1343);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+2536,dvrr_stack+2476,dvrr_stack+1361,6);


 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+2644, dvrr_stack+2536, NULL);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+2680,dvrr_stack+50,dvrr_stack+40,1);


 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+1207, dvrr_stack+2680, NULL);

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+2710,dvrr_stack+253,dvrr_stack+2476,6);


 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+2890, dvrr_stack+2710, NULL);

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+2950,dvrr_stack+718,dvrr_stack+50,1);


 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+2995, dvrr_stack+2950, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3010, dvrr_stack+697, dvrr_stack+1241, NULL, NULL, dvrr_stack+148);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3073, dvrr_stack+739, dvrr_stack+3010, dvrr_stack+718, dvrr_stack+697, dvrr_stack+163);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+3199,dvrr_stack+3073,dvrr_stack+253,6);


 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+3469, dvrr_stack+3199, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+697, dvrr_stack+1316, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+3010, dvrr_stack+2536, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+703, dvrr_stack+2680, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+739, dvrr_stack+2710, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+3046, dvrr_stack+2950, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+3559, dvrr_stack+3199, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+3061, dvrr_stack+1316, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+3649, dvrr_stack+2536, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+2536, dvrr_stack+2680, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+2546, dvrr_stack+2710, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+2680, dvrr_stack+2950, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+2695, dvrr_stack+3199, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+3067, dvrr_stack+40, dvrr_stack+12);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+799, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+1316, dvrr_stack+129, dvrr_stack+1334, dvrr_stack+12, dvrr_stack+3, dvrr_stack+799);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+3199, dvrr_stack+2476, dvrr_stack+1316);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+3235, dvrr_stack+50, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+3245, dvrr_stack+253, dvrr_stack+1361);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+3305, dvrr_stack+718, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+3320, dvrr_stack+3073, dvrr_stack+2476);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+129, dvrr_stack+40, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+3410, dvrr_stack+2476, dvrr_stack+1316);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+3446, dvrr_stack+50, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+2785, dvrr_stack+253, dvrr_stack+1361);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+2950, dvrr_stack+718, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+3685, dvrr_stack+3073, dvrr_stack+2476);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+3456, dvrr_stack+40, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+2845, dvrr_stack+2476, dvrr_stack+1316);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1316, dvrr_stack+50, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+3775, dvrr_stack+253, dvrr_stack+1361);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+2965, dvrr_stack+718, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+3835, dvrr_stack+3073, dvrr_stack+2476);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,6,dvrr_stack+3073, dvrr_stack+343, NULL);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+12, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+2881, dvrr_stack+0, dvrr_stack+21, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+3079, dvrr_stack+1334, dvrr_stack+2881, dvrr_stack+3, dvrr_stack+0, dvrr_stack+12);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+3097, dvrr_stack+24, dvrr_stack+68, NULL, NULL, dvrr_stack+21);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3115, dvrr_stack+1343, dvrr_stack+3097, dvrr_stack+6, dvrr_stack+24, dvrr_stack+2881);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+3925, dvrr_stack+1361, dvrr_stack+3115, dvrr_stack+343, dvrr_stack+1343, dvrr_stack+3079);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,6,dvrr_stack+3151, dvrr_stack+3925, dvrr_stack+343);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+3079, dvrr_stack+361, NULL);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1326, dvrr_stack+74, dvrr_stack+138, NULL, NULL, dvrr_stack+68);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+3985, dvrr_stack+84, dvrr_stack+1326, dvrr_stack+30, dvrr_stack+74, dvrr_stack+3097);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+4045, dvrr_stack+2476, dvrr_stack+3985, dvrr_stack+361, dvrr_stack+84, dvrr_stack+3115);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+3089, dvrr_stack+4045, dvrr_stack+361);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+0, dvrr_stack+208, NULL);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+65, dvrr_stack+148, dvrr_stack+1226, NULL, NULL, dvrr_stack+138);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4145, dvrr_stack+163, dvrr_stack+65, dvrr_stack+114, dvrr_stack+148, dvrr_stack+1326);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4235, dvrr_stack+253, dvrr_stack+4145, dvrr_stack+208, dvrr_stack+163, dvrr_stack+3985);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+4145, dvrr_stack+4235, dvrr_stack+208);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,6,dvrr_stack+3985, dvrr_stack+343, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,6,dvrr_stack+3991, dvrr_stack+3925, dvrr_stack+343);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+4027, dvrr_stack+361, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+65, dvrr_stack+4045, dvrr_stack+361);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+2980, dvrr_stack+208, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+4385, dvrr_stack+4235, dvrr_stack+208);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,6,dvrr_stack+4037, dvrr_stack+343, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,6,dvrr_stack+135, dvrr_stack+3925, dvrr_stack+343);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+343, dvrr_stack+361, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+3925, dvrr_stack+4045, dvrr_stack+361);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+353, dvrr_stack+208, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+4043, dvrr_stack+4235, dvrr_stack+208);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1217, dvrr_stack+253, dvrr_stack+50);
 tmp = dvrr_stack + 1217;
 target_ptr = Libderiv->deriv_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+4235, dvrr_stack+253, dvrr_stack+50);
 tmp = dvrr_stack + 4235;
 target_ptr = Libderiv->deriv_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+4280, dvrr_stack+253, dvrr_stack+50);
 tmp = dvrr_stack + 4280;
 target_ptr = Libderiv->deriv_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,18,dvrr_stack+4325, dvrr_stack+643, NULL);
 tmp = dvrr_stack + 4325;
 target_ptr = Libderiv->deriv2_classes[1][2][143];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,30,dvrr_stack+4343, dvrr_stack+1117, NULL);
 tmp = dvrr_stack + 4343;
 target_ptr = Libderiv->deriv2_classes[1][3][143];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,45,dvrr_stack+171, dvrr_stack+1861, NULL);
 tmp = dvrr_stack + 171;
 target_ptr = Libderiv->deriv2_classes[1][4][143];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,18,dvrr_stack+368, dvrr_stack+643, NULL);
 tmp = dvrr_stack + 368;
 target_ptr = Libderiv->deriv2_classes[1][2][131];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,30,dvrr_stack+1326, dvrr_stack+1117, NULL);
 tmp = dvrr_stack + 1326;
 target_ptr = Libderiv->deriv2_classes[1][3][131];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,45,dvrr_stack+216, dvrr_stack+1861, NULL);
 tmp = dvrr_stack + 216;
 target_ptr = Libderiv->deriv2_classes[1][4][131];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+21, dvrr_stack+1262, NULL);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv2_classes[1][2][130];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+261, dvrr_stack+1996, NULL);
 tmp = dvrr_stack + 261;
 target_ptr = Libderiv->deriv2_classes[1][3][130];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+291, dvrr_stack+2086, NULL);
 tmp = dvrr_stack + 291;
 target_ptr = Libderiv->deriv2_classes[1][4][130];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,18,dvrr_stack+2606, dvrr_stack+643, NULL);
 tmp = dvrr_stack + 2606;
 target_ptr = Libderiv->deriv2_classes[1][2][119];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,30,dvrr_stack+643, dvrr_stack+1117, NULL);
 tmp = dvrr_stack + 643;
 target_ptr = Libderiv->deriv2_classes[1][3][119];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,45,dvrr_stack+1117, dvrr_stack+1861, NULL);
 tmp = dvrr_stack + 1117;
 target_ptr = Libderiv->deriv2_classes[1][4][119];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+1861, dvrr_stack+1262, NULL);
 tmp = dvrr_stack + 1861;
 target_ptr = Libderiv->deriv2_classes[1][2][118];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+1262, dvrr_stack+1996, NULL);
 tmp = dvrr_stack + 1262;
 target_ptr = Libderiv->deriv2_classes[1][3][118];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+1162, dvrr_stack+2086, NULL);
 tmp = dvrr_stack + 1162;
 target_ptr = Libderiv->deriv2_classes[1][4][118];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+1292, dvrr_stack+2221, NULL);
 tmp = dvrr_stack + 1292;
 target_ptr = Libderiv->deriv2_classes[1][2][117];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+1879, dvrr_stack+535, NULL);
 tmp = dvrr_stack + 1879;
 target_ptr = Libderiv->deriv2_classes[1][3][117];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+535, dvrr_stack+937, NULL);
 tmp = dvrr_stack + 535;
 target_ptr = Libderiv->deriv2_classes[1][4][117];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+580, dvrr_stack+1627, dvrr_stack+1618);
 tmp = dvrr_stack + 580;
 target_ptr = Libderiv->deriv2_classes[1][2][107];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+1909, dvrr_stack+1072, dvrr_stack+625);
 tmp = dvrr_stack + 1909;
 target_ptr = Libderiv->deriv2_classes[1][3][107];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+1939, dvrr_stack+1657, dvrr_stack+1627);
 tmp = dvrr_stack + 1939;
 target_ptr = Libderiv->deriv2_classes[1][4][107];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+598, dvrr_stack+1729, dvrr_stack+1720);
 tmp = dvrr_stack + 598;
 target_ptr = Libderiv->deriv2_classes[1][2][106];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+1984, dvrr_stack+1777, dvrr_stack+1759);
 tmp = dvrr_stack + 1984;
 target_ptr = Libderiv->deriv2_classes[1][3][106];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+2014, dvrr_stack+2275, dvrr_stack+1729);
 tmp = dvrr_stack + 2014;
 target_ptr = Libderiv->deriv2_classes[1][4][106];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+2059, dvrr_stack+1831, dvrr_stack+1822);
 tmp = dvrr_stack + 2059;
 target_ptr = Libderiv->deriv2_classes[1][2][105];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2077, dvrr_stack+391, dvrr_stack+445);
 tmp = dvrr_stack + 2077;
 target_ptr = Libderiv->deriv2_classes[1][3][105];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+2107, dvrr_stack+802, dvrr_stack+1831);
 tmp = dvrr_stack + 2107;
 target_ptr = Libderiv->deriv2_classes[1][4][105];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+2152, dvrr_stack+1402, dvrr_stack+436);
 tmp = dvrr_stack + 2152;
 target_ptr = Libderiv->deriv2_classes[1][2][104];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2170, dvrr_stack+1450, dvrr_stack+1432);
 tmp = dvrr_stack + 2170;
 target_ptr = Libderiv->deriv2_classes[1][3][104];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+2200, dvrr_stack+1495, dvrr_stack+1402);
 tmp = dvrr_stack + 2200;
 target_ptr = Libderiv->deriv2_classes[1][4][104];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+2245, dvrr_stack+1627, dvrr_stack+1618);
 tmp = dvrr_stack + 2245;
 target_ptr = Libderiv->deriv2_classes[1][2][95];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+934, dvrr_stack+1072, dvrr_stack+625);
 tmp = dvrr_stack + 934;
 target_ptr = Libderiv->deriv2_classes[1][3][95];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+964, dvrr_stack+1657, dvrr_stack+1627);
 tmp = dvrr_stack + 964;
 target_ptr = Libderiv->deriv2_classes[1][4][95];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+673, dvrr_stack+1729, dvrr_stack+1720);
 tmp = dvrr_stack + 673;
 target_ptr = Libderiv->deriv2_classes[1][2][94];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1009, dvrr_stack+1777, dvrr_stack+1759);
 tmp = dvrr_stack + 1009;
 target_ptr = Libderiv->deriv2_classes[1][3][94];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+4475, dvrr_stack+2275, dvrr_stack+1729);
 tmp = dvrr_stack + 4475;
 target_ptr = Libderiv->deriv2_classes[1][4][94];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+2624, dvrr_stack+1831, dvrr_stack+1822);
 tmp = dvrr_stack + 2624;
 target_ptr = Libderiv->deriv2_classes[1][2][93];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1039, dvrr_stack+391, dvrr_stack+445);
 tmp = dvrr_stack + 1039;
 target_ptr = Libderiv->deriv2_classes[1][3][93];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+4520, dvrr_stack+802, dvrr_stack+1831);
 tmp = dvrr_stack + 4520;
 target_ptr = Libderiv->deriv2_classes[1][4][93];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+713, dvrr_stack+1402, dvrr_stack+436);
 tmp = dvrr_stack + 713;
 target_ptr = Libderiv->deriv2_classes[1][2][92];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+4565, dvrr_stack+1450, dvrr_stack+1432);
 tmp = dvrr_stack + 4565;
 target_ptr = Libderiv->deriv2_classes[1][3][92];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+4595, dvrr_stack+1495, dvrr_stack+1402);
 tmp = dvrr_stack + 4595;
 target_ptr = Libderiv->deriv2_classes[1][4][92];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+4640, dvrr_stack+1567, dvrr_stack+1558);
 tmp = dvrr_stack + 4640;
 target_ptr = Libderiv->deriv2_classes[1][2][91];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+4658, dvrr_stack+865, dvrr_stack+1597);
 tmp = dvrr_stack + 4658;
 target_ptr = Libderiv->deriv2_classes[1][3][91];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+4688, dvrr_stack+463, dvrr_stack+1567);
 tmp = dvrr_stack + 4688;
 target_ptr = Libderiv->deriv2_classes[1][4][91];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+4733, dvrr_stack+1627, dvrr_stack+1618);
 tmp = dvrr_stack + 4733;
 target_ptr = Libderiv->deriv2_classes[1][2][83];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+4751, dvrr_stack+1072, dvrr_stack+625);
 tmp = dvrr_stack + 4751;
 target_ptr = Libderiv->deriv2_classes[1][3][83];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+1069, dvrr_stack+1657, dvrr_stack+1627);
 tmp = dvrr_stack + 1069;
 target_ptr = Libderiv->deriv2_classes[1][4][83];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+616, dvrr_stack+1729, dvrr_stack+1720);
 tmp = dvrr_stack + 616;
 target_ptr = Libderiv->deriv2_classes[1][2][82];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+1615, dvrr_stack+1777, dvrr_stack+1759);
 tmp = dvrr_stack + 1615;
 target_ptr = Libderiv->deriv2_classes[1][3][82];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+1759, dvrr_stack+2275, dvrr_stack+1729);
 tmp = dvrr_stack + 1759;
 target_ptr = Libderiv->deriv2_classes[1][4][82];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+1804, dvrr_stack+1831, dvrr_stack+1822);
 tmp = dvrr_stack + 1804;
 target_ptr = Libderiv->deriv2_classes[1][2][81];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+2263, dvrr_stack+391, dvrr_stack+445);
 tmp = dvrr_stack + 2263;
 target_ptr = Libderiv->deriv2_classes[1][3][81];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+2293, dvrr_stack+802, dvrr_stack+1831);
 tmp = dvrr_stack + 2293;
 target_ptr = Libderiv->deriv2_classes[1][4][81];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+445, dvrr_stack+1402, dvrr_stack+436);
 tmp = dvrr_stack + 445;
 target_ptr = Libderiv->deriv2_classes[1][2][80];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+1822, dvrr_stack+1450, dvrr_stack+1432);
 tmp = dvrr_stack + 1822;
 target_ptr = Libderiv->deriv2_classes[1][3][80];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+1432, dvrr_stack+1495, dvrr_stack+1402);
 tmp = dvrr_stack + 1432;
 target_ptr = Libderiv->deriv2_classes[1][4][80];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+1477, dvrr_stack+1567, dvrr_stack+1558);
 tmp = dvrr_stack + 1477;
 target_ptr = Libderiv->deriv2_classes[1][2][79];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+1495, dvrr_stack+865, dvrr_stack+1597);
 tmp = dvrr_stack + 1495;
 target_ptr = Libderiv->deriv2_classes[1][3][79];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+386, dvrr_stack+463, dvrr_stack+1567);
 tmp = dvrr_stack + 386;
 target_ptr = Libderiv->deriv2_classes[1][4][79];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+463, dvrr_stack+2338, dvrr_stack+526);
 tmp = dvrr_stack + 463;
 target_ptr = Libderiv->deriv2_classes[1][2][78];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+481, dvrr_stack+2368, dvrr_stack+910);
 tmp = dvrr_stack + 481;
 target_ptr = Libderiv->deriv2_classes[1][3][78];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+2368, dvrr_stack+2413, dvrr_stack+2338);
 tmp = dvrr_stack + 2368;
 target_ptr = Libderiv->deriv2_classes[1][4][78];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,6,dvrr_stack+2338, dvrr_stack+2644, dvrr_stack+928);
 tmp = dvrr_stack + 2338;
 target_ptr = Libderiv->deriv2_classes[1][2][35];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,10,dvrr_stack+2413, dvrr_stack+2890, dvrr_stack+1207);
 tmp = dvrr_stack + 2413;
 target_ptr = Libderiv->deriv2_classes[1][3][35];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1525, dvrr_stack+3469, dvrr_stack+2995);
 tmp = dvrr_stack + 1525;
 target_ptr = Libderiv->deriv2_classes[1][4][35];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+2443, dvrr_stack+3010, dvrr_stack+697);
 tmp = dvrr_stack + 2443;
 target_ptr = Libderiv->deriv2_classes[1][2][34];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+1570, dvrr_stack+739, dvrr_stack+703);
 tmp = dvrr_stack + 1570;
 target_ptr = Libderiv->deriv2_classes[1][3][34];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+799, dvrr_stack+3559, dvrr_stack+3046);
 tmp = dvrr_stack + 799;
 target_ptr = Libderiv->deriv2_classes[1][4][34];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+511, dvrr_stack+3649, dvrr_stack+3061);
 tmp = dvrr_stack + 511;
 target_ptr = Libderiv->deriv2_classes[1][2][33];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+844, dvrr_stack+2546, dvrr_stack+2536);
 tmp = dvrr_stack + 844;
 target_ptr = Libderiv->deriv2_classes[1][3][33];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+874, dvrr_stack+2695, dvrr_stack+2680);
 tmp = dvrr_stack + 874;
 target_ptr = Libderiv->deriv2_classes[1][4][33];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+1645, dvrr_stack+3199, dvrr_stack+3067);
 tmp = dvrr_stack + 1645;
 target_ptr = Libderiv->deriv2_classes[1][2][32];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+1663, dvrr_stack+3245, dvrr_stack+3235);
 tmp = dvrr_stack + 1663;
 target_ptr = Libderiv->deriv2_classes[1][3][32];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1693, dvrr_stack+3320, dvrr_stack+3305);
 tmp = dvrr_stack + 1693;
 target_ptr = Libderiv->deriv2_classes[1][4][32];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+1738, dvrr_stack+3410, dvrr_stack+129);
 tmp = dvrr_stack + 1738;
 target_ptr = Libderiv->deriv2_classes[1][2][31];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+1397, dvrr_stack+2785, dvrr_stack+3446);
 tmp = dvrr_stack + 1397;
 target_ptr = Libderiv->deriv2_classes[1][3][31];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+4781, dvrr_stack+3685, dvrr_stack+2950);
 tmp = dvrr_stack + 4781;
 target_ptr = Libderiv->deriv2_classes[1][4][31];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+4826, dvrr_stack+1361, dvrr_stack+15);
 tmp = dvrr_stack + 4826;
 target_ptr = Libderiv->deriv_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+4844, dvrr_stack+2845, dvrr_stack+3456);
 tmp = dvrr_stack + 4844;
 target_ptr = Libderiv->deriv2_classes[1][2][30];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+4862, dvrr_stack+2476, dvrr_stack+40);
 tmp = dvrr_stack + 4862;
 target_ptr = Libderiv->deriv_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+4892, dvrr_stack+3775, dvrr_stack+1316);
 tmp = dvrr_stack + 4892;
 target_ptr = Libderiv->deriv2_classes[1][3][30];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+4922, dvrr_stack+3835, dvrr_stack+2965);
 tmp = dvrr_stack + 4922;
 target_ptr = Libderiv->deriv2_classes[1][4][30];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,6,dvrr_stack+4967, dvrr_stack+3151, dvrr_stack+3073);
 tmp = dvrr_stack + 4967;
 target_ptr = Libderiv->deriv2_classes[1][2][26];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+4985, dvrr_stack+3089, dvrr_stack+3079);
 tmp = dvrr_stack + 4985;
 target_ptr = Libderiv->deriv2_classes[1][3][26];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+5015, dvrr_stack+4145, dvrr_stack+0);
 tmp = dvrr_stack + 5015;
 target_ptr = Libderiv->deriv2_classes[1][4][26];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,6,dvrr_stack+5060, dvrr_stack+2644, dvrr_stack+928);
 tmp = dvrr_stack + 5060;
 target_ptr = Libderiv->deriv2_classes[1][2][23];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,10,dvrr_stack+5078, dvrr_stack+2890, dvrr_stack+1207);
 tmp = dvrr_stack + 5078;
 target_ptr = Libderiv->deriv2_classes[1][3][23];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,15,dvrr_stack+5108, dvrr_stack+3469, dvrr_stack+2995);
 tmp = dvrr_stack + 5108;
 target_ptr = Libderiv->deriv2_classes[1][4][23];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+5153, dvrr_stack+3010, dvrr_stack+697);
 tmp = dvrr_stack + 5153;
 target_ptr = Libderiv->deriv2_classes[1][2][22];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+5171, dvrr_stack+739, dvrr_stack+703);
 tmp = dvrr_stack + 5171;
 target_ptr = Libderiv->deriv2_classes[1][3][22];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+5201, dvrr_stack+3559, dvrr_stack+3046);
 tmp = dvrr_stack + 5201;
 target_ptr = Libderiv->deriv2_classes[1][4][22];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+5246, dvrr_stack+3649, dvrr_stack+3061);
 tmp = dvrr_stack + 5246;
 target_ptr = Libderiv->deriv2_classes[1][2][21];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+5264, dvrr_stack+2546, dvrr_stack+2536);
 tmp = dvrr_stack + 5264;
 target_ptr = Libderiv->deriv2_classes[1][3][21];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+5294, dvrr_stack+2695, dvrr_stack+2680);
 tmp = dvrr_stack + 5294;
 target_ptr = Libderiv->deriv2_classes[1][4][21];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+5339, dvrr_stack+3199, dvrr_stack+3067);
 tmp = dvrr_stack + 5339;
 target_ptr = Libderiv->deriv2_classes[1][2][20];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+5357, dvrr_stack+3245, dvrr_stack+3235);
 tmp = dvrr_stack + 5357;
 target_ptr = Libderiv->deriv2_classes[1][3][20];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+5387, dvrr_stack+3320, dvrr_stack+3305);
 tmp = dvrr_stack + 5387;
 target_ptr = Libderiv->deriv2_classes[1][4][20];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+5432, dvrr_stack+3410, dvrr_stack+129);
 tmp = dvrr_stack + 5432;
 target_ptr = Libderiv->deriv2_classes[1][2][19];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+5450, dvrr_stack+2785, dvrr_stack+3446);
 tmp = dvrr_stack + 5450;
 target_ptr = Libderiv->deriv2_classes[1][3][19];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+5480, dvrr_stack+3685, dvrr_stack+2950);
 tmp = dvrr_stack + 5480;
 target_ptr = Libderiv->deriv2_classes[1][4][19];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+5525, dvrr_stack+1361, dvrr_stack+15);
 tmp = dvrr_stack + 5525;
 target_ptr = Libderiv->deriv_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+5543, dvrr_stack+2845, dvrr_stack+3456);
 tmp = dvrr_stack + 5543;
 target_ptr = Libderiv->deriv2_classes[1][2][18];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+5561, dvrr_stack+2476, dvrr_stack+40);
 tmp = dvrr_stack + 5561;
 target_ptr = Libderiv->deriv_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+5591, dvrr_stack+3775, dvrr_stack+1316);
 tmp = dvrr_stack + 5591;
 target_ptr = Libderiv->deriv2_classes[1][3][18];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+5621, dvrr_stack+3835, dvrr_stack+2965);
 tmp = dvrr_stack + 5621;
 target_ptr = Libderiv->deriv2_classes[1][4][18];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+5666, dvrr_stack+3151, dvrr_stack+3073);
 tmp = dvrr_stack + 5666;
 target_ptr = Libderiv->deriv2_classes[1][2][14];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+5684, dvrr_stack+3089, dvrr_stack+3079);
 tmp = dvrr_stack + 5684;
 target_ptr = Libderiv->deriv2_classes[1][3][14];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+5714, dvrr_stack+4145, dvrr_stack+0);
 tmp = dvrr_stack + 5714;
 target_ptr = Libderiv->deriv2_classes[1][4][14];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,6,dvrr_stack+5759, dvrr_stack+3991, dvrr_stack+3985);
 tmp = dvrr_stack + 5759;
 target_ptr = Libderiv->deriv2_classes[1][2][13];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+5777, dvrr_stack+65, dvrr_stack+4027);
 tmp = dvrr_stack + 5777;
 target_ptr = Libderiv->deriv2_classes[1][3][13];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+5807, dvrr_stack+4385, dvrr_stack+2980);
 tmp = dvrr_stack + 5807;
 target_ptr = Libderiv->deriv2_classes[1][4][13];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,6,dvrr_stack+5852, dvrr_stack+2644, dvrr_stack+928);
 tmp = dvrr_stack + 5852;
 target_ptr = Libderiv->deriv2_classes[1][2][11];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,10,dvrr_stack+2642, dvrr_stack+2890, dvrr_stack+1207);
 tmp = dvrr_stack + 2642;
 target_ptr = Libderiv->deriv2_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,15,dvrr_stack+2881, dvrr_stack+3469, dvrr_stack+2995);
 tmp = dvrr_stack + 2881;
 target_ptr = Libderiv->deriv2_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2926, dvrr_stack+3010, dvrr_stack+697);
 tmp = dvrr_stack + 2926;
 target_ptr = Libderiv->deriv2_classes[1][2][10];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+2995, dvrr_stack+739, dvrr_stack+703);
 tmp = dvrr_stack + 2995;
 target_ptr = Libderiv->deriv2_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3462, dvrr_stack+3559, dvrr_stack+3046);
 tmp = dvrr_stack + 3462;
 target_ptr = Libderiv->deriv2_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+3025, dvrr_stack+3649, dvrr_stack+3061);
 tmp = dvrr_stack + 3025;
 target_ptr = Libderiv->deriv2_classes[1][2][9];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3507, dvrr_stack+2546, dvrr_stack+2536);
 tmp = dvrr_stack + 3507;
 target_ptr = Libderiv->deriv2_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+2536, dvrr_stack+2695, dvrr_stack+2680);
 tmp = dvrr_stack + 2536;
 target_ptr = Libderiv->deriv2_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2581, dvrr_stack+3199, dvrr_stack+3067);
 tmp = dvrr_stack + 2581;
 target_ptr = Libderiv->deriv2_classes[1][2][8];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3043, dvrr_stack+3245, dvrr_stack+3235);
 tmp = dvrr_stack + 3043;
 target_ptr = Libderiv->deriv2_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3187, dvrr_stack+3320, dvrr_stack+3305);
 tmp = dvrr_stack + 3187;
 target_ptr = Libderiv->deriv2_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+691, dvrr_stack+3410, dvrr_stack+129);
 tmp = dvrr_stack + 691;
 target_ptr = Libderiv->deriv2_classes[1][2][7];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3232, dvrr_stack+2785, dvrr_stack+3446);
 tmp = dvrr_stack + 3232;
 target_ptr = Libderiv->deriv2_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3262, dvrr_stack+3685, dvrr_stack+2950);
 tmp = dvrr_stack + 3262;
 target_ptr = Libderiv->deriv2_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2944, dvrr_stack+1361, dvrr_stack+15);
 tmp = dvrr_stack + 2944;
 target_ptr = Libderiv->deriv_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+1356, dvrr_stack+2845, dvrr_stack+3456);
 tmp = dvrr_stack + 1356;
 target_ptr = Libderiv->deriv2_classes[1][2][6];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3307, dvrr_stack+2476, dvrr_stack+40);
 tmp = dvrr_stack + 3307;
 target_ptr = Libderiv->deriv_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+2461, dvrr_stack+3775, dvrr_stack+1316);
 tmp = dvrr_stack + 2461;
 target_ptr = Libderiv->deriv2_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+2491, dvrr_stack+3835, dvrr_stack+2965);
 tmp = dvrr_stack + 2491;
 target_ptr = Libderiv->deriv2_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+2962, dvrr_stack+3151, dvrr_stack+3073);
 tmp = dvrr_stack + 2962;
 target_ptr = Libderiv->deriv2_classes[1][2][2];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3149, dvrr_stack+3089, dvrr_stack+3079);
 tmp = dvrr_stack + 3149;
 target_ptr = Libderiv->deriv2_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3073, dvrr_stack+4145, dvrr_stack+0);
 tmp = dvrr_stack + 3073;
 target_ptr = Libderiv->deriv2_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+0, dvrr_stack+3991, dvrr_stack+3985);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[1][2][1];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3985, dvrr_stack+65, dvrr_stack+4027);
 tmp = dvrr_stack + 3985;
 target_ptr = Libderiv->deriv2_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+4133, dvrr_stack+4385, dvrr_stack+2980);
 tmp = dvrr_stack + 4133;
 target_ptr = Libderiv->deriv2_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 2 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,6,dvrr_stack+4015, dvrr_stack+135, dvrr_stack+4037);
 tmp = dvrr_stack + 4015;
 target_ptr = Libderiv->deriv2_classes[1][2][0];
 for(i=0;i<18;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3118, dvrr_stack+3925, dvrr_stack+343);
 tmp = dvrr_stack + 3118;
 target_ptr = Libderiv->deriv2_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+4373, dvrr_stack+4043, dvrr_stack+353);
 tmp = dvrr_stack + 4373;
 target_ptr = Libderiv->deriv2_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];


}

