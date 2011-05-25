#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|gg) integrals */

void d1vrr_order_00gg(Libderiv_t *Libderiv, prim_data *Data)
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
 tmp = dvrr_stack + 50;
 target_ptr = Libderiv->dvrr_classes[0][4];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+12, dvrr_stack+21, dvrr_stack+3, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+65, dvrr_stack+24, dvrr_stack+12, dvrr_stack+0, dvrr_stack+21, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+75, dvrr_stack+30, dvrr_stack+65, dvrr_stack+6, dvrr_stack+24, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+90, dvrr_stack+50, dvrr_stack+75, dvrr_stack+40, dvrr_stack+30, NULL);
 tmp = dvrr_stack + 90;
 target_ptr = Libderiv->dvrr_classes[0][5];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+111,dvrr_stack+90,dvrr_stack+50,1);


 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+156, dvrr_stack+12, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+166, dvrr_stack+65, dvrr_stack+156, dvrr_stack+24, dvrr_stack+12, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+181, dvrr_stack+75, dvrr_stack+166, dvrr_stack+30, dvrr_stack+65, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+202, dvrr_stack+90, dvrr_stack+181, dvrr_stack+50, dvrr_stack+75, NULL);
 tmp = dvrr_stack + 202;
 target_ptr = Libderiv->dvrr_classes[0][6];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+230,dvrr_stack+202,dvrr_stack+90,1);


 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+18, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+21, dvrr_stack+0, dvrr_stack+18, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+293, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+303, dvrr_stack+156, dvrr_stack+293, dvrr_stack+12, dvrr_stack+6, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+318, dvrr_stack+166, dvrr_stack+303, dvrr_stack+65, dvrr_stack+156, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+339, dvrr_stack+181, dvrr_stack+318, dvrr_stack+75, dvrr_stack+166, NULL);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+367, dvrr_stack+202, dvrr_stack+339, dvrr_stack+90, dvrr_stack+181, NULL);
 tmp = dvrr_stack + 367;
 target_ptr = Libderiv->dvrr_classes[0][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+403,dvrr_stack+367,dvrr_stack+202,1);


 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+12, dvrr_stack+18, dvrr_stack+3, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+65, dvrr_stack+21, dvrr_stack+12, dvrr_stack+0, dvrr_stack+18, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+487, dvrr_stack+293, dvrr_stack+65, dvrr_stack+6, dvrr_stack+21, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+502, dvrr_stack+303, dvrr_stack+487, dvrr_stack+156, dvrr_stack+293, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+523, dvrr_stack+318, dvrr_stack+502, dvrr_stack+166, dvrr_stack+303, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+551, dvrr_stack+339, dvrr_stack+523, dvrr_stack+181, dvrr_stack+318, NULL);

 /* compute (0 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+587, dvrr_stack+367, dvrr_stack+551, dvrr_stack+202, dvrr_stack+339, NULL);

 /* compute (0 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+632,dvrr_stack+587,dvrr_stack+367,1);


 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+156, dvrr_stack+12, dvrr_stack+6, dvrr_stack+18, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+166, dvrr_stack+65, dvrr_stack+156, dvrr_stack+21, dvrr_stack+12, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+487, dvrr_stack+166, dvrr_stack+293, dvrr_stack+65, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+740, dvrr_stack+502, dvrr_stack+0, dvrr_stack+303, dvrr_stack+487, NULL);

 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+768, dvrr_stack+523, dvrr_stack+740, dvrr_stack+318, dvrr_stack+502, NULL);

 /* compute (0 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+293, dvrr_stack+551, dvrr_stack+768, dvrr_stack+339, dvrr_stack+523, NULL);

 /* compute (0 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+740, dvrr_stack+587, dvrr_stack+293, dvrr_stack+367, dvrr_stack+551, NULL);

 /* compute (0 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+795,dvrr_stack+740,dvrr_stack+587,1);


 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+930, dvrr_stack+50, dvrr_stack+75, NULL, NULL, dvrr_stack+30);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+487, dvrr_stack+90, dvrr_stack+181, NULL, NULL, dvrr_stack+75);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+975, dvrr_stack+202, dvrr_stack+339, NULL, NULL, dvrr_stack+181);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1059, dvrr_stack+367, dvrr_stack+551, NULL, NULL, dvrr_stack+339);

 /* compute (1 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1167, dvrr_stack+587, dvrr_stack+293, NULL, NULL, dvrr_stack+551);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+293, dvrr_stack+111, NULL);
 tmp = dvrr_stack + 293;
 target_ptr = Libderiv->deriv_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,21,dvrr_stack+308, dvrr_stack+230, NULL);
 tmp = dvrr_stack + 308;
 target_ptr = Libderiv->deriv_classes[0][5][11];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,28,dvrr_stack+329, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 329;
 target_ptr = Libderiv->deriv_classes[0][6][11];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+1302, dvrr_stack+632, NULL);
 tmp = dvrr_stack + 1302;
 target_ptr = Libderiv->deriv_classes[0][7][11];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+1338, dvrr_stack+795, NULL);
 tmp = dvrr_stack + 1338;
 target_ptr = Libderiv->deriv_classes[0][8][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+1383, dvrr_stack+111, NULL);
 tmp = dvrr_stack + 1383;
 target_ptr = Libderiv->deriv_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+1398, dvrr_stack+230, NULL);
 tmp = dvrr_stack + 1398;
 target_ptr = Libderiv->deriv_classes[0][5][10];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,28,dvrr_stack+1419, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 1419;
 target_ptr = Libderiv->deriv_classes[0][6][10];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+1447, dvrr_stack+632, NULL);
 tmp = dvrr_stack + 1447;
 target_ptr = Libderiv->deriv_classes[0][7][10];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+1483, dvrr_stack+795, NULL);
 tmp = dvrr_stack + 1483;
 target_ptr = Libderiv->deriv_classes[0][8][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+1528, dvrr_stack+111, NULL);
 tmp = dvrr_stack + 1528;
 target_ptr = Libderiv->deriv_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+111, dvrr_stack+230, NULL);
 tmp = dvrr_stack + 111;
 target_ptr = Libderiv->deriv_classes[0][5][9];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,28,dvrr_stack+230, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 230;
 target_ptr = Libderiv->deriv_classes[0][6][9];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+403, dvrr_stack+632, NULL);
 tmp = dvrr_stack + 403;
 target_ptr = Libderiv->deriv_classes[0][7][9];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+632, dvrr_stack+795, NULL);
 tmp = dvrr_stack + 632;
 target_ptr = Libderiv->deriv_classes[0][8][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+795, dvrr_stack+90, dvrr_stack+40);
 tmp = dvrr_stack + 795;
 target_ptr = Libderiv->deriv_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+810, dvrr_stack+202, dvrr_stack+50);
 tmp = dvrr_stack + 810;
 target_ptr = Libderiv->deriv_classes[0][5][8];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,1,1,dvrr_stack+831, dvrr_stack+367, dvrr_stack+90);
 tmp = dvrr_stack + 831;
 target_ptr = Libderiv->deriv_classes[0][6][8];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,1,1,dvrr_stack+859, dvrr_stack+587, dvrr_stack+202);
 tmp = dvrr_stack + 859;
 target_ptr = Libderiv->deriv_classes[0][7][8];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,1,1,dvrr_stack+677, dvrr_stack+740, dvrr_stack+367);
 tmp = dvrr_stack + 677;
 target_ptr = Libderiv->deriv_classes[0][8][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+895, dvrr_stack+90, dvrr_stack+40);
 tmp = dvrr_stack + 895;
 target_ptr = Libderiv->deriv_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+439, dvrr_stack+202, dvrr_stack+50);
 tmp = dvrr_stack + 439;
 target_ptr = Libderiv->deriv_classes[0][5][7];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+258, dvrr_stack+367, dvrr_stack+90);
 tmp = dvrr_stack + 258;
 target_ptr = Libderiv->deriv_classes[0][6][7];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,1,1,dvrr_stack+132, dvrr_stack+587, dvrr_stack+202);
 tmp = dvrr_stack + 132;
 target_ptr = Libderiv->deriv_classes[0][7][7];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,1,1,dvrr_stack+1543, dvrr_stack+740, dvrr_stack+367);
 tmp = dvrr_stack + 1543;
 target_ptr = Libderiv->deriv_classes[0][8][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+910, dvrr_stack+90, dvrr_stack+40);
 tmp = dvrr_stack + 910;
 target_ptr = Libderiv->deriv_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+460, dvrr_stack+202, dvrr_stack+50);
 tmp = dvrr_stack + 460;
 target_ptr = Libderiv->deriv_classes[0][5][6];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+168, dvrr_stack+367, dvrr_stack+90);
 tmp = dvrr_stack + 168;
 target_ptr = Libderiv->deriv_classes[0][6][6];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,1,1,dvrr_stack+1588, dvrr_stack+587, dvrr_stack+202);
 tmp = dvrr_stack + 1588;
 target_ptr = Libderiv->deriv_classes[0][7][6];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,1,1,dvrr_stack+1624, dvrr_stack+740, dvrr_stack+367);
 tmp = dvrr_stack + 1624;
 target_ptr = Libderiv->deriv_classes[0][8][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+722, dvrr_stack+930, NULL);
 tmp = dvrr_stack + 722;
 target_ptr = Libderiv->deriv_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+737, dvrr_stack+487, NULL);
 tmp = dvrr_stack + 737;
 target_ptr = Libderiv->deriv_classes[0][5][2];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+758, dvrr_stack+975, NULL);
 tmp = dvrr_stack + 758;
 target_ptr = Libderiv->deriv_classes[0][6][2];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,36,dvrr_stack+357, dvrr_stack+1059, NULL);
 tmp = dvrr_stack + 357;
 target_ptr = Libderiv->deriv_classes[0][7][2];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,45,dvrr_stack+1669, dvrr_stack+1167, NULL);
 tmp = dvrr_stack + 1669;
 target_ptr = Libderiv->deriv_classes[0][8][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+196, dvrr_stack+930, NULL);
 tmp = dvrr_stack + 196;
 target_ptr = Libderiv->deriv_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+1714, dvrr_stack+487, NULL);
 tmp = dvrr_stack + 1714;
 target_ptr = Libderiv->deriv_classes[0][5][1];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+1735, dvrr_stack+975, NULL);
 tmp = dvrr_stack + 1735;
 target_ptr = Libderiv->deriv_classes[0][6][1];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,36,dvrr_stack+1763, dvrr_stack+1059, NULL);
 tmp = dvrr_stack + 1763;
 target_ptr = Libderiv->deriv_classes[0][7][1];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,45,dvrr_stack+1799, dvrr_stack+1167, NULL);
 tmp = dvrr_stack + 1799;
 target_ptr = Libderiv->deriv_classes[0][8][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+211, dvrr_stack+930, NULL);
 tmp = dvrr_stack + 211;
 target_ptr = Libderiv->deriv_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+925, dvrr_stack+487, NULL);
 tmp = dvrr_stack + 925;
 target_ptr = Libderiv->deriv_classes[0][5][0];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+946, dvrr_stack+975, NULL);
 tmp = dvrr_stack + 946;
 target_ptr = Libderiv->deriv_classes[0][6][0];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,36,dvrr_stack+974, dvrr_stack+1059, NULL);
 tmp = dvrr_stack + 974;
 target_ptr = Libderiv->deriv_classes[0][7][0];
 for(i=0;i<36;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,45,dvrr_stack+1010, dvrr_stack+1167, NULL);
 tmp = dvrr_stack + 1010;
 target_ptr = Libderiv->deriv_classes[0][8][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];


}

