#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (p0|ff) integrals */

void d1vrr_order_p0ff(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+0, dvrr_stack+12, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+6, dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+34, dvrr_stack+31, dvrr_stack+3, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+34, dvrr_stack+6, dvrr_stack+31, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+50, dvrr_stack+40, dvrr_stack+21, NULL, NULL, dvrr_stack+6);
 tmp = dvrr_stack + 50;
 target_ptr = Libderiv->dvrr_classes[1][3];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+80, dvrr_stack+12, dvrr_stack+31, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+86, dvrr_stack+15, dvrr_stack+80, dvrr_stack+0, dvrr_stack+12, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+96, dvrr_stack+21, dvrr_stack+86, dvrr_stack+6, dvrr_stack+15, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+111, dvrr_stack+40, dvrr_stack+21, dvrr_stack+34, dvrr_stack+6, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+126, dvrr_stack+111, dvrr_stack+96, NULL, NULL, dvrr_stack+21);
 tmp = dvrr_stack + 126;
 target_ptr = Libderiv->dvrr_classes[1][4];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+171,dvrr_stack+126,dvrr_stack+50,3);


 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+261, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+264, dvrr_stack+31, dvrr_stack+261, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+270, dvrr_stack+80, dvrr_stack+264, dvrr_stack+12, dvrr_stack+31, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+280, dvrr_stack+86, dvrr_stack+270, dvrr_stack+15, dvrr_stack+80, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+295, dvrr_stack+96, dvrr_stack+280, dvrr_stack+21, dvrr_stack+86, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+316, dvrr_stack+111, dvrr_stack+96, dvrr_stack+40, dvrr_stack+21, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+337, dvrr_stack+316, dvrr_stack+295, NULL, NULL, dvrr_stack+96);
 tmp = dvrr_stack + 337;
 target_ptr = Libderiv->dvrr_classes[1][5];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+400,dvrr_stack+337,dvrr_stack+126,3);


 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+535, dvrr_stack+261, dvrr_stack+12, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+541, dvrr_stack+264, dvrr_stack+535, dvrr_stack+31, dvrr_stack+261, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+551, dvrr_stack+270, dvrr_stack+541, dvrr_stack+80, dvrr_stack+264, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+566, dvrr_stack+280, dvrr_stack+551, dvrr_stack+86, dvrr_stack+270, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+587, dvrr_stack+295, dvrr_stack+566, dvrr_stack+96, dvrr_stack+280, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+615, dvrr_stack+316, dvrr_stack+295, dvrr_stack+111, dvrr_stack+96, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+643, dvrr_stack+615, dvrr_stack+587, NULL, NULL, dvrr_stack+295);

 /* compute (1 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+727,dvrr_stack+643,dvrr_stack+337,3);


 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+80, dvrr_stack+12, dvrr_stack+31, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+916, dvrr_stack+535, dvrr_stack+80, dvrr_stack+261, dvrr_stack+12, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+926, dvrr_stack+541, dvrr_stack+916, dvrr_stack+264, dvrr_stack+535, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+941, dvrr_stack+551, dvrr_stack+926, dvrr_stack+270, dvrr_stack+541, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+962, dvrr_stack+566, dvrr_stack+941, dvrr_stack+280, dvrr_stack+551, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+916, dvrr_stack+587, dvrr_stack+962, dvrr_stack+295, dvrr_stack+566, NULL);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+990, dvrr_stack+615, dvrr_stack+587, dvrr_stack+316, dvrr_stack+295, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1026, dvrr_stack+990, dvrr_stack+916, NULL, NULL, dvrr_stack+587);

 /* compute (1 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+1134,dvrr_stack+1026,dvrr_stack+643,3);


 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+916, dvrr_stack+34, dvrr_stack+6, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+934, dvrr_stack+6, dvrr_stack+15, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+990, dvrr_stack+21, dvrr_stack+86, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1386, dvrr_stack+50, dvrr_stack+990, dvrr_stack+40, dvrr_stack+21, dvrr_stack+934);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1446, dvrr_stack+96, dvrr_stack+280, NULL, NULL, dvrr_stack+86);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1491, dvrr_stack+126, dvrr_stack+1446, dvrr_stack+111, dvrr_stack+96, dvrr_stack+990);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1581, dvrr_stack+295, dvrr_stack+566, NULL, NULL, dvrr_stack+280);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1644, dvrr_stack+337, dvrr_stack+1581, dvrr_stack+316, dvrr_stack+295, dvrr_stack+1446);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1770, dvrr_stack+587, dvrr_stack+962, NULL, NULL, dvrr_stack+566);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1854, dvrr_stack+643, dvrr_stack+1770, dvrr_stack+615, dvrr_stack+587, dvrr_stack+1581);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+1581, dvrr_stack+171, NULL);
 tmp = dvrr_stack + 1581;
 target_ptr = Libderiv->deriv_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+1446, dvrr_stack+400, NULL);
 tmp = dvrr_stack + 1446;
 target_ptr = Libderiv->deriv_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,63,dvrr_stack+1770, dvrr_stack+727, NULL);
 tmp = dvrr_stack + 1770;
 target_ptr = Libderiv->deriv_classes[1][5][11];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,84,dvrr_stack+934, dvrr_stack+1134, NULL);
 tmp = dvrr_stack + 934;
 target_ptr = Libderiv->deriv_classes[1][6][11];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+1611, dvrr_stack+171, NULL);
 tmp = dvrr_stack + 1611;
 target_ptr = Libderiv->deriv_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+535, dvrr_stack+400, NULL);
 tmp = dvrr_stack + 535;
 target_ptr = Libderiv->deriv_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+2022, dvrr_stack+727, NULL);
 tmp = dvrr_stack + 2022;
 target_ptr = Libderiv->deriv_classes[1][5][10];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,84,dvrr_stack+2085, dvrr_stack+1134, NULL);
 tmp = dvrr_stack + 2085;
 target_ptr = Libderiv->deriv_classes[1][6][10];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+0, dvrr_stack+171, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+171, dvrr_stack+400, NULL);
 tmp = dvrr_stack + 171;
 target_ptr = Libderiv->deriv_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+400, dvrr_stack+727, NULL);
 tmp = dvrr_stack + 400;
 target_ptr = Libderiv->deriv_classes[1][5][9];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,84,dvrr_stack+727, dvrr_stack+1134, NULL);
 tmp = dvrr_stack + 727;
 target_ptr = Libderiv->deriv_classes[1][6][9];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+1134, dvrr_stack+126, dvrr_stack+916);
 tmp = dvrr_stack + 1134;
 target_ptr = Libderiv->deriv_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+1164, dvrr_stack+337, dvrr_stack+50);
 tmp = dvrr_stack + 1164;
 target_ptr = Libderiv->deriv_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+1209, dvrr_stack+643, dvrr_stack+126);
 tmp = dvrr_stack + 1209;
 target_ptr = Libderiv->deriv_classes[1][5][8];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+1272, dvrr_stack+1026, dvrr_stack+337);
 tmp = dvrr_stack + 1272;
 target_ptr = Libderiv->deriv_classes[1][6][8];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1356, dvrr_stack+126, dvrr_stack+916);
 tmp = dvrr_stack + 1356;
 target_ptr = Libderiv->deriv_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+811, dvrr_stack+337, dvrr_stack+50);
 tmp = dvrr_stack + 811;
 target_ptr = Libderiv->deriv_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+463, dvrr_stack+643, dvrr_stack+126);
 tmp = dvrr_stack + 463;
 target_ptr = Libderiv->deriv_classes[1][5][7];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+216, dvrr_stack+1026, dvrr_stack+337);
 tmp = dvrr_stack + 216;
 target_ptr = Libderiv->deriv_classes[1][6][7];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+856, dvrr_stack+126, dvrr_stack+916);
 tmp = dvrr_stack + 856;
 target_ptr = Libderiv->deriv_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+886, dvrr_stack+337, dvrr_stack+50);
 tmp = dvrr_stack + 886;
 target_ptr = Libderiv->deriv_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+2169, dvrr_stack+643, dvrr_stack+126);
 tmp = dvrr_stack + 2169;
 target_ptr = Libderiv->deriv_classes[1][5][6];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+643, dvrr_stack+1026, dvrr_stack+337);
 tmp = dvrr_stack + 643;
 target_ptr = Libderiv->deriv_classes[1][6][6];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+337, dvrr_stack+1386, dvrr_stack+40);
 tmp = dvrr_stack + 337;
 target_ptr = Libderiv->deriv_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+126, dvrr_stack+1491, dvrr_stack+111);
 tmp = dvrr_stack + 126;
 target_ptr = Libderiv->deriv_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+1018, dvrr_stack+1644, dvrr_stack+316);
 tmp = dvrr_stack + 1018;
 target_ptr = Libderiv->deriv_classes[1][5][2];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+2232, dvrr_stack+1854, dvrr_stack+615);
 tmp = dvrr_stack + 2232;
 target_ptr = Libderiv->deriv_classes[1][6][2];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+367, dvrr_stack+1386, dvrr_stack+40);
 tmp = dvrr_stack + 367;
 target_ptr = Libderiv->deriv_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+50, dvrr_stack+1491, dvrr_stack+111);
 tmp = dvrr_stack + 50;
 target_ptr = Libderiv->deriv_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+2316, dvrr_stack+1644, dvrr_stack+316);
 tmp = dvrr_stack + 2316;
 target_ptr = Libderiv->deriv_classes[1][5][1];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+2379, dvrr_stack+1854, dvrr_stack+615);
 tmp = dvrr_stack + 2379;
 target_ptr = Libderiv->deriv_classes[1][6][1];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+1081, dvrr_stack+1386, dvrr_stack+40);
 tmp = dvrr_stack + 1081;
 target_ptr = Libderiv->deriv_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+1386, dvrr_stack+1491, dvrr_stack+111);
 tmp = dvrr_stack + 1386;
 target_ptr = Libderiv->deriv_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+1491, dvrr_stack+1644, dvrr_stack+316);
 tmp = dvrr_stack + 1491;
 target_ptr = Libderiv->deriv_classes[1][5][0];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+1641, dvrr_stack+1854, dvrr_stack+615);
 tmp = dvrr_stack + 1641;
 target_ptr = Libderiv->deriv_classes[1][6][0];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];


}

