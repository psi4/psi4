#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (pp|ff) integrals */

void d12vrr_order_ppff(Libderiv_t *Libderiv, prim_data *Data)
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
 _BUILD_00p0(Data,dvrr_stack+80, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+83, dvrr_stack+12, dvrr_stack+80, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+89, dvrr_stack+15, dvrr_stack+83, dvrr_stack+0, dvrr_stack+12, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+99, dvrr_stack+21, dvrr_stack+89, dvrr_stack+6, dvrr_stack+15, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+114, dvrr_stack+40, dvrr_stack+21, dvrr_stack+34, dvrr_stack+6, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+129, dvrr_stack+114, dvrr_stack+99, NULL, NULL, dvrr_stack+21);
 tmp = dvrr_stack + 129;
 target_ptr = Libderiv->dvrr_classes[1][4];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+174,dvrr_stack+129,dvrr_stack+50,3);


 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+264, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+267, dvrr_stack+80, dvrr_stack+264, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+273, dvrr_stack+83, dvrr_stack+267, dvrr_stack+12, dvrr_stack+80, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+283, dvrr_stack+89, dvrr_stack+273, dvrr_stack+15, dvrr_stack+83, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+298, dvrr_stack+99, dvrr_stack+283, dvrr_stack+21, dvrr_stack+89, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+319, dvrr_stack+114, dvrr_stack+99, dvrr_stack+40, dvrr_stack+21, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+340, dvrr_stack+319, dvrr_stack+298, NULL, NULL, dvrr_stack+99);
 tmp = dvrr_stack + 340;
 target_ptr = Libderiv->dvrr_classes[1][5];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+403,dvrr_stack+340,dvrr_stack+129,3);


 /* compute (1 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+538,dvrr_stack+403,dvrr_stack+174,3);


 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,30,dvrr_stack+718, dvrr_stack+538, dvrr_stack+50);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+808, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+811, dvrr_stack+264, dvrr_stack+808, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+817, dvrr_stack+267, dvrr_stack+811, dvrr_stack+80, dvrr_stack+264, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+827, dvrr_stack+273, dvrr_stack+817, dvrr_stack+83, dvrr_stack+267, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+842, dvrr_stack+283, dvrr_stack+827, dvrr_stack+89, dvrr_stack+273, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+863, dvrr_stack+298, dvrr_stack+842, dvrr_stack+99, dvrr_stack+283, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+891, dvrr_stack+319, dvrr_stack+298, dvrr_stack+114, dvrr_stack+99, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+919, dvrr_stack+891, dvrr_stack+863, NULL, NULL, dvrr_stack+298);
 tmp = dvrr_stack + 919;
 target_ptr = Libderiv->dvrr_classes[1][6];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+1003,dvrr_stack+919,dvrr_stack+340,3);


 /* compute (1 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+1192,dvrr_stack+1003,dvrr_stack+403,3);


 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,45,dvrr_stack+1462, dvrr_stack+1192, dvrr_stack+129);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1597, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1600, dvrr_stack+808, dvrr_stack+1597, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1606, dvrr_stack+811, dvrr_stack+1600, dvrr_stack+264, dvrr_stack+808, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1616, dvrr_stack+817, dvrr_stack+1606, dvrr_stack+267, dvrr_stack+811, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1631, dvrr_stack+827, dvrr_stack+1616, dvrr_stack+273, dvrr_stack+817, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1652, dvrr_stack+842, dvrr_stack+1631, dvrr_stack+283, dvrr_stack+827, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1680, dvrr_stack+863, dvrr_stack+1652, dvrr_stack+298, dvrr_stack+842, NULL);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1716, dvrr_stack+891, dvrr_stack+863, dvrr_stack+319, dvrr_stack+298, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1752, dvrr_stack+1716, dvrr_stack+1680, NULL, NULL, dvrr_stack+863);

 /* compute (1 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+1860,dvrr_stack+1752,dvrr_stack+919,3);


 /* compute (1 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+2112,dvrr_stack+1860,dvrr_stack+1003,3);


 /* compute (1 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,63,dvrr_stack+2490, dvrr_stack+2112, dvrr_stack+340);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+264, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+2679, dvrr_stack+1597, dvrr_stack+264, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+2685, dvrr_stack+1600, dvrr_stack+2679, dvrr_stack+808, dvrr_stack+1597, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2695, dvrr_stack+1606, dvrr_stack+2685, dvrr_stack+811, dvrr_stack+1600, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2710, dvrr_stack+1616, dvrr_stack+2695, dvrr_stack+817, dvrr_stack+1606, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2731, dvrr_stack+1631, dvrr_stack+2710, dvrr_stack+827, dvrr_stack+1616, NULL);

 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2759, dvrr_stack+1652, dvrr_stack+2731, dvrr_stack+842, dvrr_stack+1631, NULL);

 /* compute (0 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+2795, dvrr_stack+1680, dvrr_stack+2759, dvrr_stack+863, dvrr_stack+1652, NULL);

 /* compute (0 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+2840, dvrr_stack+1716, dvrr_stack+1680, dvrr_stack+891, dvrr_stack+863, NULL);

 /* compute (1 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+2885, dvrr_stack+2840, dvrr_stack+2795, NULL, NULL, dvrr_stack+1680);

 /* compute (1 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+3020,dvrr_stack+2885,dvrr_stack+1752,3);


 /* compute (1 0 | 6 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_id(Libderiv->CD,dvrr_stack+3344,dvrr_stack+3020,dvrr_stack+1860,3);


 /* compute (1 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,84,dvrr_stack+3848, dvrr_stack+3344, dvrr_stack+919);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+4100, dvrr_stack+6, dvrr_stack+15, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+4118, dvrr_stack+21, dvrr_stack+89, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4148, dvrr_stack+50, dvrr_stack+4118, dvrr_stack+40, dvrr_stack+21, dvrr_stack+4100);
 tmp = dvrr_stack + 4148;
 target_ptr = Libderiv->dvrr_classes[2][3];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4208, dvrr_stack+99, dvrr_stack+283, NULL, NULL, dvrr_stack+89);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4253, dvrr_stack+129, dvrr_stack+4208, dvrr_stack+114, dvrr_stack+99, dvrr_stack+4118);
 tmp = dvrr_stack + 4253;
 target_ptr = Libderiv->dvrr_classes[2][4];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+4343,dvrr_stack+4253,dvrr_stack+4148,6);


 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4523, dvrr_stack+298, dvrr_stack+842, NULL, NULL, dvrr_stack+283);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4586, dvrr_stack+340, dvrr_stack+4523, dvrr_stack+319, dvrr_stack+298, dvrr_stack+4208);
 tmp = dvrr_stack + 4586;
 target_ptr = Libderiv->dvrr_classes[2][5];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+4712,dvrr_stack+4586,dvrr_stack+4253,6);


 /* compute (2 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+4982,dvrr_stack+4712,dvrr_stack+4343,6);


 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,60,dvrr_stack+5342, dvrr_stack+4982, dvrr_stack+4148);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+5522, dvrr_stack+863, dvrr_stack+1652, NULL, NULL, dvrr_stack+842);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+5606, dvrr_stack+919, dvrr_stack+5522, dvrr_stack+891, dvrr_stack+863, dvrr_stack+4523);

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+5774,dvrr_stack+5606,dvrr_stack+4586,6);


 /* compute (2 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+6152,dvrr_stack+5774,dvrr_stack+4712,6);


 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,90,dvrr_stack+6692, dvrr_stack+6152, dvrr_stack+4253);

 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+6962, dvrr_stack+1680, dvrr_stack+2759, NULL, NULL, dvrr_stack+1652);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7070, dvrr_stack+1752, dvrr_stack+6962, dvrr_stack+1716, dvrr_stack+1680, dvrr_stack+5522);

 /* compute (2 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+7286,dvrr_stack+7070,dvrr_stack+5606,6);


 /* compute (2 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+7790,dvrr_stack+7286,dvrr_stack+5774,6);


 /* compute (2 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,126,dvrr_stack+8546, dvrr_stack+7790, dvrr_stack+4586);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+808, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+811, dvrr_stack+264, dvrr_stack+808, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+8924, dvrr_stack+2679, dvrr_stack+811, dvrr_stack+1597, dvrr_stack+264, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8934, dvrr_stack+2685, dvrr_stack+8924, dvrr_stack+1600, dvrr_stack+2679, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+8949, dvrr_stack+2695, dvrr_stack+8934, dvrr_stack+1606, dvrr_stack+2685, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+8970, dvrr_stack+2710, dvrr_stack+8949, dvrr_stack+1616, dvrr_stack+2695, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+8924, dvrr_stack+2731, dvrr_stack+8970, dvrr_stack+1631, dvrr_stack+2710, NULL);

 /* compute (0 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+8998, dvrr_stack+2759, dvrr_stack+8924, dvrr_stack+1652, dvrr_stack+2731, NULL);

 /* compute (1 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9043, dvrr_stack+2795, dvrr_stack+8998, NULL, NULL, dvrr_stack+2759);

 /* compute (2 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9178, dvrr_stack+2885, dvrr_stack+9043, dvrr_stack+2840, dvrr_stack+2795, dvrr_stack+6962);

 /* compute (2 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+9448,dvrr_stack+9178,dvrr_stack+7070,6);


 /* compute (2 0 | 6 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_id(Libderiv->CD,dvrr_stack+10096,dvrr_stack+9448,dvrr_stack+7286,6);


 /* compute (2 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,168,dvrr_stack+11104, dvrr_stack+10096, dvrr_stack+5606);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,30,dvrr_stack+2795, dvrr_stack+538, dvrr_stack+50);

 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,45,dvrr_stack+8998, dvrr_stack+1192, dvrr_stack+129);

 /* compute (1 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,63,dvrr_stack+11608, dvrr_stack+2112, dvrr_stack+340);

 /* compute (1 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,84,dvrr_stack+11797, dvrr_stack+3344, dvrr_stack+919);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,60,dvrr_stack+12049, dvrr_stack+4982, dvrr_stack+4148);

 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,90,dvrr_stack+12229, dvrr_stack+6152, dvrr_stack+4253);

 /* compute (2 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,126,dvrr_stack+12499, dvrr_stack+7790, dvrr_stack+4586);

 /* compute (2 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,168,dvrr_stack+12877, dvrr_stack+10096, dvrr_stack+5606);

 /* compute (1 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,30,dvrr_stack+13381, dvrr_stack+538, dvrr_stack+50);

 /* compute (1 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,45,dvrr_stack+538, dvrr_stack+1192, dvrr_stack+129);

 /* compute (1 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,63,dvrr_stack+1192, dvrr_stack+2112, dvrr_stack+340);

 /* compute (1 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,84,dvrr_stack+2112, dvrr_stack+3344, dvrr_stack+919);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,60,dvrr_stack+3344, dvrr_stack+4982, dvrr_stack+4148);

 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,90,dvrr_stack+4982, dvrr_stack+6152, dvrr_stack+4253);

 /* compute (2 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,126,dvrr_stack+6152, dvrr_stack+7790, dvrr_stack+4586);

 /* compute (2 0 | 6 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,168,dvrr_stack+7790, dvrr_stack+10096, dvrr_stack+5606);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+10096, dvrr_stack+34, dvrr_stack+6, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+10114,dvrr_stack+50,dvrr_stack+10096,3);


 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+10168, dvrr_stack+10114, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+673, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 673;
 target_ptr = Libderiv->deriv_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+10186, dvrr_stack+174, NULL);
 tmp = dvrr_stack + 10186;
 target_ptr = Libderiv->deriv_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,63,dvrr_stack+10216, dvrr_stack+1003, NULL);
 tmp = dvrr_stack + 10216;
 target_ptr = Libderiv->deriv_classes[1][5][11];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,84,dvrr_stack+10279, dvrr_stack+1860, NULL);
 tmp = dvrr_stack + 10279;
 target_ptr = Libderiv->deriv_classes[1][6][11];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,108,dvrr_stack+10363, dvrr_stack+3020, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+808, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+10471, dvrr_stack+10096, dvrr_stack+4100, dvrr_stack+34, dvrr_stack+6, dvrr_stack+808);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+10507,dvrr_stack+4148,dvrr_stack+10471,6);


 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+10615, dvrr_stack+10507, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+5252, dvrr_stack+4712, NULL);
 tmp = dvrr_stack + 5252;
 target_ptr = Libderiv->deriv_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+10651, dvrr_stack+4343, NULL);
 tmp = dvrr_stack + 10651;
 target_ptr = Libderiv->deriv_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+2364, dvrr_stack+5774, NULL);
 tmp = dvrr_stack + 2364;
 target_ptr = Libderiv->deriv_classes[2][5][11];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,168,dvrr_stack+10711, dvrr_stack+7286, NULL);
 tmp = dvrr_stack + 10711;
 target_ptr = Libderiv->deriv_classes[2][6][11];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,216,dvrr_stack+10879, dvrr_stack+9448, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+8294, dvrr_stack+10114, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+9133, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 9133;
 target_ptr = Libderiv->deriv_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+8312, dvrr_stack+174, NULL);
 tmp = dvrr_stack + 8312;
 target_ptr = Libderiv->deriv_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+8342, dvrr_stack+1003, NULL);
 tmp = dvrr_stack + 8342;
 target_ptr = Libderiv->deriv_classes[1][5][10];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,84,dvrr_stack+8405, dvrr_stack+1860, NULL);
 tmp = dvrr_stack + 8405;
 target_ptr = Libderiv->deriv_classes[1][6][10];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,108,dvrr_stack+6530, dvrr_stack+3020, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+8489, dvrr_stack+10507, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+3524, dvrr_stack+4712, NULL);
 tmp = dvrr_stack + 3524;
 target_ptr = Libderiv->deriv_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+3614, dvrr_stack+4343, NULL);
 tmp = dvrr_stack + 3614;
 target_ptr = Libderiv->deriv_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+3674, dvrr_stack+5774, NULL);
 tmp = dvrr_stack + 3674;
 target_ptr = Libderiv->deriv_classes[2][5][10];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,168,dvrr_stack+13471, dvrr_stack+7286, NULL);
 tmp = dvrr_stack + 13471;
 target_ptr = Libderiv->deriv_classes[2][6][10];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,216,dvrr_stack+13639, dvrr_stack+9448, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+8525, dvrr_stack+10114, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+10114, dvrr_stack+403, NULL);
 tmp = dvrr_stack + 10114;
 target_ptr = Libderiv->deriv_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+403, dvrr_stack+174, NULL);
 tmp = dvrr_stack + 403;
 target_ptr = Libderiv->deriv_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+174, dvrr_stack+1003, NULL);
 tmp = dvrr_stack + 174;
 target_ptr = Libderiv->deriv_classes[1][5][9];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,84,dvrr_stack+1003, dvrr_stack+1860, NULL);
 tmp = dvrr_stack + 1003;
 target_ptr = Libderiv->deriv_classes[1][6][9];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,108,dvrr_stack+1860, dvrr_stack+3020, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+3020, dvrr_stack+10507, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+10507, dvrr_stack+4712, NULL);
 tmp = dvrr_stack + 10507;
 target_ptr = Libderiv->deriv_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+4712, dvrr_stack+4343, NULL);
 tmp = dvrr_stack + 4712;
 target_ptr = Libderiv->deriv_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+4343, dvrr_stack+5774, NULL);
 tmp = dvrr_stack + 4343;
 target_ptr = Libderiv->deriv_classes[2][5][9];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,168,dvrr_stack+5774, dvrr_stack+7286, NULL);
 tmp = dvrr_stack + 5774;
 target_ptr = Libderiv->deriv_classes[2][6][9];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,216,dvrr_stack+7286, dvrr_stack+9448, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+10159, dvrr_stack+31, dvrr_stack+3, NULL, NULL, Data->F+1);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+10597, dvrr_stack+50, dvrr_stack+10159);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+9448, dvrr_stack+340, dvrr_stack+50);
 tmp = dvrr_stack + 9448;
 target_ptr = Libderiv->deriv_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+237, dvrr_stack+129, dvrr_stack+10096);
 tmp = dvrr_stack + 237;
 target_ptr = Libderiv->deriv_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+9493, dvrr_stack+919, dvrr_stack+129);
 tmp = dvrr_stack + 9493;
 target_ptr = Libderiv->deriv_classes[1][5][8];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+9556, dvrr_stack+1752, dvrr_stack+340);
 tmp = dvrr_stack + 9556;
 target_ptr = Libderiv->deriv_classes[1][6][8];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,3,1,dvrr_stack+9640, dvrr_stack+2885, dvrr_stack+919);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+8543, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+9748, dvrr_stack+10159, dvrr_stack+808, dvrr_stack+31, dvrr_stack+3, dvrr_stack+8543);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+9766, dvrr_stack+4148, dvrr_stack+9748);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+9802, dvrr_stack+4586, dvrr_stack+4148);
 tmp = dvrr_stack + 9802;
 target_ptr = Libderiv->deriv_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+9892, dvrr_stack+4253, dvrr_stack+10471);
 tmp = dvrr_stack + 9892;
 target_ptr = Libderiv->deriv_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+9952, dvrr_stack+5606, dvrr_stack+4253);
 tmp = dvrr_stack + 9952;
 target_ptr = Libderiv->deriv_classes[2][5][8];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+7502, dvrr_stack+7070, dvrr_stack+4586);
 tmp = dvrr_stack + 7502;
 target_ptr = Libderiv->deriv_classes[2][6][8];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,6,1,dvrr_stack+3056, dvrr_stack+9178, dvrr_stack+5606);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+10078, dvrr_stack+50, dvrr_stack+10159);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+7670, dvrr_stack+340, dvrr_stack+50);
 tmp = dvrr_stack + 7670;
 target_ptr = Libderiv->deriv_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+7715, dvrr_stack+129, dvrr_stack+10096);
 tmp = dvrr_stack + 7715;
 target_ptr = Libderiv->deriv_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+5942, dvrr_stack+919, dvrr_stack+129);
 tmp = dvrr_stack + 5942;
 target_ptr = Libderiv->deriv_classes[1][5][7];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+6005, dvrr_stack+1752, dvrr_stack+340);
 tmp = dvrr_stack + 6005;
 target_ptr = Libderiv->deriv_classes[1][6][7];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,3,1,dvrr_stack+4772, dvrr_stack+2885, dvrr_stack+919);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+7745, dvrr_stack+4148, dvrr_stack+9748);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+4880, dvrr_stack+4586, dvrr_stack+4148);
 tmp = dvrr_stack + 4880;
 target_ptr = Libderiv->deriv_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+6089, dvrr_stack+4253, dvrr_stack+10471);
 tmp = dvrr_stack + 6089;
 target_ptr = Libderiv->deriv_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+1968, dvrr_stack+5606, dvrr_stack+4253);
 tmp = dvrr_stack + 1968;
 target_ptr = Libderiv->deriv_classes[2][5][7];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+13855, dvrr_stack+7070, dvrr_stack+4586);
 tmp = dvrr_stack + 13855;
 target_ptr = Libderiv->deriv_classes[2][6][7];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,6,1,dvrr_stack+14023, dvrr_stack+9178, dvrr_stack+5606);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+2094, dvrr_stack+50, dvrr_stack+10159);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+4469, dvrr_stack+340, dvrr_stack+50);
 tmp = dvrr_stack + 4469;
 target_ptr = Libderiv->deriv_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+3272, dvrr_stack+129, dvrr_stack+10096);
 tmp = dvrr_stack + 3272;
 target_ptr = Libderiv->deriv_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+1087, dvrr_stack+919, dvrr_stack+129);
 tmp = dvrr_stack + 1087;
 target_ptr = Libderiv->deriv_classes[1][5][6];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+433, dvrr_stack+1752, dvrr_stack+340);
 tmp = dvrr_stack + 433;
 target_ptr = Libderiv->deriv_classes[1][6][6];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,3,1,dvrr_stack+14239, dvrr_stack+2885, dvrr_stack+919);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+2885, dvrr_stack+4148, dvrr_stack+9748);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+2921, dvrr_stack+4586, dvrr_stack+4148);
 tmp = dvrr_stack + 2921;
 target_ptr = Libderiv->deriv_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+1381, dvrr_stack+4253, dvrr_stack+10471);
 tmp = dvrr_stack + 1381;
 target_ptr = Libderiv->deriv_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+14347, dvrr_stack+5606, dvrr_stack+4253);
 tmp = dvrr_stack + 14347;
 target_ptr = Libderiv->deriv_classes[2][5][6];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+14473, dvrr_stack+7070, dvrr_stack+4586);
 tmp = dvrr_stack + 14473;
 target_ptr = Libderiv->deriv_classes[2][6][6];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,6,1,dvrr_stack+14641, dvrr_stack+9178, dvrr_stack+5606);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+9178,dvrr_stack+114,dvrr_stack+40,1);


 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+8960, dvrr_stack+9178, NULL);

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+9208,dvrr_stack+319,dvrr_stack+114,1);


 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+9253, dvrr_stack+9208, NULL);

 /* compute (0 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+9268,dvrr_stack+891,dvrr_stack+319,1);


 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,21,dvrr_stack+517, dvrr_stack+9268, NULL);

 /* compute (0 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+9331,dvrr_stack+1716,dvrr_stack+891,1);


 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,28,dvrr_stack+9415, dvrr_stack+9331, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3011, dvrr_stack+0, dvrr_stack+12, NULL, NULL, Data->F+3);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+9748, dvrr_stack+15, dvrr_stack+83, NULL, NULL, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3302, dvrr_stack+4100, dvrr_stack+9748, dvrr_stack+6, dvrr_stack+15, dvrr_stack+3011);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+1150, dvrr_stack+89, dvrr_stack+273, NULL, NULL, dvrr_stack+83);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+14857, dvrr_stack+4118, dvrr_stack+1150, dvrr_stack+21, dvrr_stack+89, dvrr_stack+9748);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+14917, dvrr_stack+4148, dvrr_stack+14857, dvrr_stack+50, dvrr_stack+4118, dvrr_stack+3302);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6638, dvrr_stack+283, dvrr_stack+827, NULL, NULL, dvrr_stack+273);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15017, dvrr_stack+4208, dvrr_stack+6638, dvrr_stack+99, dvrr_stack+283, dvrr_stack+1150);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+15107, dvrr_stack+4253, dvrr_stack+15017, dvrr_stack+129, dvrr_stack+4208, dvrr_stack+14857);

 /* compute (3 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+15257,dvrr_stack+15107,dvrr_stack+14917,10);


 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+15557, dvrr_stack+15257, NULL);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15657, dvrr_stack+842, dvrr_stack+1631, NULL, NULL, dvrr_stack+827);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15720, dvrr_stack+4523, dvrr_stack+15657, dvrr_stack+298, dvrr_stack+842, dvrr_stack+6638);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15846, dvrr_stack+4586, dvrr_stack+15720, dvrr_stack+340, dvrr_stack+4523, dvrr_stack+15017);

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+16056,dvrr_stack+15846,dvrr_stack+15107,10);


 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+16506, dvrr_stack+16056, NULL);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16656, dvrr_stack+1652, dvrr_stack+2731, NULL, NULL, dvrr_stack+1631);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16740, dvrr_stack+5522, dvrr_stack+16656, dvrr_stack+863, dvrr_stack+1652, dvrr_stack+15657);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16908, dvrr_stack+5606, dvrr_stack+16740, dvrr_stack+919, dvrr_stack+5522, dvrr_stack+15720);

 /* compute (3 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+17188,dvrr_stack+16908,dvrr_stack+15846,10);


 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+17818, dvrr_stack+17188, NULL);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18028, dvrr_stack+2759, dvrr_stack+8924, NULL, NULL, dvrr_stack+2731);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18136, dvrr_stack+6962, dvrr_stack+18028, dvrr_stack+1680, dvrr_stack+2759, dvrr_stack+16656);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18352, dvrr_stack+7070, dvrr_stack+18136, dvrr_stack+1752, dvrr_stack+6962, dvrr_stack+16740);

 /* compute (3 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+18712,dvrr_stack+18352,dvrr_stack+16908,10);


 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,280,dvrr_stack+6962, dvrr_stack+18712, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+7242, dvrr_stack+9178, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+99, dvrr_stack+9208, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+298, dvrr_stack+9268, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,28,dvrr_stack+863, dvrr_stack+9331, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+1752, dvrr_stack+15257, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+18028, dvrr_stack+16056, NULL);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+19552, dvrr_stack+17188, NULL);

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+19762, dvrr_stack+18712, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+7252, dvrr_stack+9178, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+9178, dvrr_stack+9208, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+1441, dvrr_stack+9268, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,28,dvrr_stack+9268, dvrr_stack+9331, NULL);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+9296, dvrr_stack+15257, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+15257, dvrr_stack+16056, NULL);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+16056, dvrr_stack+17188, NULL);

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+17188, dvrr_stack+18712, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+18712, dvrr_stack+114, dvrr_stack+34);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+18722, dvrr_stack+319, dvrr_stack+40);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+18737, dvrr_stack+891, dvrr_stack+114);

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,1,1,dvrr_stack+18758, dvrr_stack+1716, dvrr_stack+319);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+8543, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+18786, dvrr_stack+808, dvrr_stack+3011, dvrr_stack+3, dvrr_stack+0, dvrr_stack+8543);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+9193, dvrr_stack+10471, dvrr_stack+3302, dvrr_stack+10096, dvrr_stack+4100, dvrr_stack+18786);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+18786, dvrr_stack+15107, dvrr_stack+9193);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+15407, dvrr_stack+15846, dvrr_stack+14917);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+18886, dvrr_stack+16908, dvrr_stack+15107);

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+19096, dvrr_stack+18352, dvrr_stack+15846);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+19376, dvrr_stack+114, dvrr_stack+34);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+19386, dvrr_stack+319, dvrr_stack+40);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+19401, dvrr_stack+891, dvrr_stack+114);

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+19422, dvrr_stack+1716, dvrr_stack+319);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+19450, dvrr_stack+15107, dvrr_stack+9193);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+17468, dvrr_stack+15846, dvrr_stack+14917);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+16266, dvrr_stack+16908, dvrr_stack+15107);

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+20042, dvrr_stack+18352, dvrr_stack+15846);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+10096, dvrr_stack+114, dvrr_stack+34);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+10471, dvrr_stack+319, dvrr_stack+40);

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+10486, dvrr_stack+891, dvrr_stack+114);

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+17618, dvrr_stack+1716, dvrr_stack+319);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+17646, dvrr_stack+15107, dvrr_stack+9193);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+18178, dvrr_stack+15846, dvrr_stack+14917);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+20322, dvrr_stack+16908, dvrr_stack+15107);

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+20532, dvrr_stack+18352, dvrr_stack+15846);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+9193, dvrr_stack+50, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+17746, dvrr_stack+14917, dvrr_stack+50);
 tmp = dvrr_stack + 17746;
 target_ptr = Libderiv->deriv_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+9203, dvrr_stack+129, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+18328, dvrr_stack+15107, dvrr_stack+129);
 tmp = dvrr_stack + 18328;
 target_ptr = Libderiv->deriv_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+9218, dvrr_stack+340, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+18418, dvrr_stack+15846, dvrr_stack+340);
 tmp = dvrr_stack + 18418;
 target_ptr = Libderiv->deriv_classes[2][5][2];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+16476, dvrr_stack+919, NULL);

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+18544, dvrr_stack+16908, dvrr_stack+919);
 tmp = dvrr_stack + 18544;
 target_ptr = Libderiv->deriv_classes[2][6][2];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+2759, dvrr_stack+4148, dvrr_stack+40);
 tmp = dvrr_stack + 2759;
 target_ptr = Libderiv->deriv_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+8543, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3, dvrr_stack+12, dvrr_stack+80, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+9396, dvrr_stack+3011, dvrr_stack+3, dvrr_stack+0, dvrr_stack+12, dvrr_stack+8543);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+7262, dvrr_stack+83, dvrr_stack+267, NULL, NULL, dvrr_stack+80);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+8924, dvrr_stack+9748, dvrr_stack+7262, dvrr_stack+15, dvrr_stack+83, dvrr_stack+3);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+1680, dvrr_stack+3302, dvrr_stack+8924, dvrr_stack+4100, dvrr_stack+9748, dvrr_stack+9396);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+3302, dvrr_stack+273, dvrr_stack+817, NULL, NULL, dvrr_stack+267);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+20812, dvrr_stack+1150, dvrr_stack+3302, dvrr_stack+89, dvrr_stack+273, dvrr_stack+7262);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+20872, dvrr_stack+14857, dvrr_stack+20812, dvrr_stack+4118, dvrr_stack+1150, dvrr_stack+8924);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+20972, dvrr_stack+14917, dvrr_stack+20872, dvrr_stack+4148, dvrr_stack+14857, dvrr_stack+1680);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+21122, dvrr_stack+20972, dvrr_stack+4148);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1680, dvrr_stack+4253, dvrr_stack+114);
 tmp = dvrr_stack + 1680;
 target_ptr = Libderiv->deriv_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14857, dvrr_stack+827, dvrr_stack+1616, NULL, NULL, dvrr_stack+817);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21222, dvrr_stack+6638, dvrr_stack+14857, dvrr_stack+283, dvrr_stack+827, dvrr_stack+3302);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21312, dvrr_stack+15017, dvrr_stack+21222, dvrr_stack+4208, dvrr_stack+6638, dvrr_stack+20812);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21462, dvrr_stack+15107, dvrr_stack+21312, dvrr_stack+4253, dvrr_stack+15017, dvrr_stack+20872);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+20812, dvrr_stack+21462, dvrr_stack+4253);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+15017, dvrr_stack+4586, dvrr_stack+319);
 tmp = dvrr_stack + 15017;
 target_ptr = Libderiv->deriv_classes[1][5][2];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+21687, dvrr_stack+1631, dvrr_stack+2710, NULL, NULL, dvrr_stack+1616);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+21750, dvrr_stack+15657, dvrr_stack+21687, dvrr_stack+842, dvrr_stack+1631, dvrr_stack+14857);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+21876, dvrr_stack+15720, dvrr_stack+21750, dvrr_stack+4523, dvrr_stack+15657, dvrr_stack+21222);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+22086, dvrr_stack+15846, dvrr_stack+21876, dvrr_stack+4586, dvrr_stack+15720, dvrr_stack+21312);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+21222, dvrr_stack+22086, dvrr_stack+4586);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+15657, dvrr_stack+5606, dvrr_stack+891);
 tmp = dvrr_stack + 15657;
 target_ptr = Libderiv->deriv_classes[1][6][2];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+15741, dvrr_stack+2731, dvrr_stack+8970, NULL, NULL, dvrr_stack+2710);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+22401, dvrr_stack+16656, dvrr_stack+15741, dvrr_stack+1652, dvrr_stack+2731, dvrr_stack+21687);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+22569, dvrr_stack+16740, dvrr_stack+22401, dvrr_stack+5522, dvrr_stack+16656, dvrr_stack+21750);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+22849, dvrr_stack+16908, dvrr_stack+22569, dvrr_stack+5606, dvrr_stack+16740, dvrr_stack+21876);

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+22401, dvrr_stack+22849, dvrr_stack+5606);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+20962, dvrr_stack+50, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+14857, dvrr_stack+14917, dvrr_stack+50);
 tmp = dvrr_stack + 14857;
 target_ptr = Libderiv->deriv_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+16656, dvrr_stack+129, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+16671, dvrr_stack+15107, dvrr_stack+129);
 tmp = dvrr_stack + 16671;
 target_ptr = Libderiv->deriv_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+16761, dvrr_stack+340, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+16782, dvrr_stack+15846, dvrr_stack+340);
 tmp = dvrr_stack + 16782;
 target_ptr = Libderiv->deriv_classes[2][5][1];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+8970, dvrr_stack+919, NULL);

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+22681, dvrr_stack+16908, dvrr_stack+919);
 tmp = dvrr_stack + 22681;
 target_ptr = Libderiv->deriv_classes[2][6][1];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+21432, dvrr_stack+4148, dvrr_stack+40);
 tmp = dvrr_stack + 21432;
 target_ptr = Libderiv->deriv_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+21687, dvrr_stack+20972, dvrr_stack+4148);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+4208, dvrr_stack+4253, dvrr_stack+114);
 tmp = dvrr_stack + 4208;
 target_ptr = Libderiv->deriv_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+21787, dvrr_stack+21462, dvrr_stack+4253);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+5522, dvrr_stack+4586, dvrr_stack+319);
 tmp = dvrr_stack + 5522;
 target_ptr = Libderiv->deriv_classes[1][5][1];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+23269, dvrr_stack+22086, dvrr_stack+4586);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+21937, dvrr_stack+5606, dvrr_stack+891);
 tmp = dvrr_stack + 21937;
 target_ptr = Libderiv->deriv_classes[1][6][1];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+23479, dvrr_stack+22849, dvrr_stack+5606);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+5585, dvrr_stack+50, NULL);

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+22021, dvrr_stack+14917, dvrr_stack+50);
 tmp = dvrr_stack + 22021;
 target_ptr = Libderiv->deriv_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+50, dvrr_stack+129, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+14917, dvrr_stack+15107, dvrr_stack+129);
 tmp = dvrr_stack + 14917;
 target_ptr = Libderiv->deriv_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+129, dvrr_stack+340, NULL);

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+15080, dvrr_stack+15846, dvrr_stack+340);
 tmp = dvrr_stack + 15080;
 target_ptr = Libderiv->deriv_classes[2][5][0];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+340, dvrr_stack+919, NULL);

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+15741, dvrr_stack+16908, dvrr_stack+919);
 tmp = dvrr_stack + 15741;
 target_ptr = Libderiv->deriv_classes[2][6][0];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+919, dvrr_stack+4148, dvrr_stack+40);
 tmp = dvrr_stack + 919;
 target_ptr = Libderiv->deriv_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+16908, dvrr_stack+20972, dvrr_stack+4148);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+20972, dvrr_stack+4253, dvrr_stack+114);
 tmp = dvrr_stack + 20972;
 target_ptr = Libderiv->deriv_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+17008, dvrr_stack+21462, dvrr_stack+4253);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+4253, dvrr_stack+4586, dvrr_stack+319);
 tmp = dvrr_stack + 4253;
 target_ptr = Libderiv->deriv_classes[1][5][0];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+21462, dvrr_stack+22086, dvrr_stack+4586);

 /* compute (1 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+21017, dvrr_stack+5606, dvrr_stack+891);
 tmp = dvrr_stack + 21017;
 target_ptr = Libderiv->deriv_classes[1][6][0];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+22081, dvrr_stack+22849, dvrr_stack+5606);

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,30,dvrr_stack+17158, dvrr_stack+718, NULL);
 tmp = dvrr_stack + 17158;
 target_ptr = Libderiv->deriv2_classes[1][3][143];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,45,dvrr_stack+22849, dvrr_stack+1462, NULL);
 tmp = dvrr_stack + 22849;
 target_ptr = Libderiv->deriv2_classes[1][4][143];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,63,dvrr_stack+22894, dvrr_stack+2490, NULL);
 tmp = dvrr_stack + 22894;
 target_ptr = Libderiv->deriv2_classes[1][5][143];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,84,dvrr_stack+22957, dvrr_stack+3848, NULL);
 tmp = dvrr_stack + 22957;
 target_ptr = Libderiv->deriv2_classes[1][6][143];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,60,dvrr_stack+23041, dvrr_stack+5342, NULL);
 tmp = dvrr_stack + 23041;
 target_ptr = Libderiv->deriv2_classes[2][3][143];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,90,dvrr_stack+23101, dvrr_stack+6692, NULL);
 tmp = dvrr_stack + 23101;
 target_ptr = Libderiv->deriv2_classes[2][4][143];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,126,dvrr_stack+5595, dvrr_stack+8546, NULL);
 tmp = dvrr_stack + 5595;
 target_ptr = Libderiv->deriv2_classes[2][5][143];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,168,dvrr_stack+4514, dvrr_stack+11104, NULL);
 tmp = dvrr_stack + 4514;
 target_ptr = Libderiv->deriv2_classes[2][6][143];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,30,dvrr_stack+4682, dvrr_stack+718, NULL);
 tmp = dvrr_stack + 4682;
 target_ptr = Libderiv->deriv2_classes[1][3][131];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,45,dvrr_stack+23191, dvrr_stack+1462, NULL);
 tmp = dvrr_stack + 23191;
 target_ptr = Libderiv->deriv2_classes[1][4][131];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,63,dvrr_stack+15909, dvrr_stack+2490, NULL);
 tmp = dvrr_stack + 15909;
 target_ptr = Libderiv->deriv2_classes[1][5][131];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,84,dvrr_stack+15972, dvrr_stack+3848, NULL);
 tmp = dvrr_stack + 15972;
 target_ptr = Libderiv->deriv2_classes[1][6][131];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,60,dvrr_stack+4100, dvrr_stack+5342, NULL);
 tmp = dvrr_stack + 4100;
 target_ptr = Libderiv->deriv2_classes[2][3][131];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,90,dvrr_stack+23759, dvrr_stack+6692, NULL);
 tmp = dvrr_stack + 23759;
 target_ptr = Libderiv->deriv2_classes[2][4][131];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,126,dvrr_stack+23849, dvrr_stack+8546, NULL);
 tmp = dvrr_stack + 23849;
 target_ptr = Libderiv->deriv2_classes[2][5][131];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,168,dvrr_stack+23975, dvrr_stack+11104, NULL);
 tmp = dvrr_stack + 23975;
 target_ptr = Libderiv->deriv2_classes[2][6][131];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+23236, dvrr_stack+2795, NULL);
 tmp = dvrr_stack + 23236;
 target_ptr = Libderiv->deriv2_classes[1][3][130];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+949, dvrr_stack+8998, NULL);
 tmp = dvrr_stack + 949;
 target_ptr = Libderiv->deriv2_classes[1][4][130];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+2679, dvrr_stack+11608, NULL);
 tmp = dvrr_stack + 2679;
 target_ptr = Libderiv->deriv2_classes[1][5][130];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,84,dvrr_stack+24143, dvrr_stack+11797, NULL);
 tmp = dvrr_stack + 24143;
 target_ptr = Libderiv->deriv2_classes[1][6][130];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+1597, dvrr_stack+12049, NULL);
 tmp = dvrr_stack + 1597;
 target_ptr = Libderiv->deriv2_classes[2][3][130];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+24227, dvrr_stack+12229, NULL);
 tmp = dvrr_stack + 24227;
 target_ptr = Libderiv->deriv2_classes[2][4][130];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+24317, dvrr_stack+12499, NULL);
 tmp = dvrr_stack + 24317;
 target_ptr = Libderiv->deriv2_classes[2][5][130];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,168,dvrr_stack+24443, dvrr_stack+12877, NULL);
 tmp = dvrr_stack + 24443;
 target_ptr = Libderiv->deriv2_classes[2][6][130];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,30,dvrr_stack+368, dvrr_stack+718, NULL);
 tmp = dvrr_stack + 368;
 target_ptr = Libderiv->deriv2_classes[1][3][119];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,45,dvrr_stack+718, dvrr_stack+1462, NULL);
 tmp = dvrr_stack + 718;
 target_ptr = Libderiv->deriv2_classes[1][4][119];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,63,dvrr_stack+1462, dvrr_stack+2490, NULL);
 tmp = dvrr_stack + 1462;
 target_ptr = Libderiv->deriv2_classes[1][5][119];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,84,dvrr_stack+2490, dvrr_stack+3848, NULL);
 tmp = dvrr_stack + 2490;
 target_ptr = Libderiv->deriv2_classes[1][6][119];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,60,dvrr_stack+2574, dvrr_stack+5342, NULL);
 tmp = dvrr_stack + 2574;
 target_ptr = Libderiv->deriv2_classes[2][3][119];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,90,dvrr_stack+5342, dvrr_stack+6692, NULL);
 tmp = dvrr_stack + 5342;
 target_ptr = Libderiv->deriv2_classes[2][4][119];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,126,dvrr_stack+6638, dvrr_stack+8546, NULL);
 tmp = dvrr_stack + 6638;
 target_ptr = Libderiv->deriv2_classes[2][5][119];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,168,dvrr_stack+6764, dvrr_stack+11104, NULL);
 tmp = dvrr_stack + 6764;
 target_ptr = Libderiv->deriv2_classes[2][6][119];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+6932, dvrr_stack+2795, NULL);
 tmp = dvrr_stack + 6932;
 target_ptr = Libderiv->deriv2_classes[1][3][118];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+2634, dvrr_stack+8998, NULL);
 tmp = dvrr_stack + 2634;
 target_ptr = Libderiv->deriv2_classes[1][4][118];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+8998, dvrr_stack+11608, NULL);
 tmp = dvrr_stack + 8998;
 target_ptr = Libderiv->deriv2_classes[1][5][118];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,84,dvrr_stack+5432, dvrr_stack+11797, NULL);
 tmp = dvrr_stack + 5432;
 target_ptr = Libderiv->deriv2_classes[1][6][118];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+9061, dvrr_stack+12049, NULL);
 tmp = dvrr_stack + 9061;
 target_ptr = Libderiv->deriv2_classes[2][3][118];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+763, dvrr_stack+12229, NULL);
 tmp = dvrr_stack + 763;
 target_ptr = Libderiv->deriv2_classes[2][4][118];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+8543, dvrr_stack+12499, NULL);
 tmp = dvrr_stack + 8543;
 target_ptr = Libderiv->deriv2_classes[2][5][118];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,168,dvrr_stack+8669, dvrr_stack+12877, NULL);
 tmp = dvrr_stack + 8669;
 target_ptr = Libderiv->deriv2_classes[2][6][118];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+1525, dvrr_stack+13381, NULL);
 tmp = dvrr_stack + 1525;
 target_ptr = Libderiv->deriv2_classes[1][3][117];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+5721, dvrr_stack+538, NULL);
 tmp = dvrr_stack + 5721;
 target_ptr = Libderiv->deriv2_classes[1][4][117];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+538, dvrr_stack+1192, NULL);
 tmp = dvrr_stack + 538;
 target_ptr = Libderiv->deriv2_classes[1][5][117];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,84,dvrr_stack+1150, dvrr_stack+2112, NULL);
 tmp = dvrr_stack + 1150;
 target_ptr = Libderiv->deriv2_classes[1][6][117];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+2112, dvrr_stack+3344, NULL);
 tmp = dvrr_stack + 2112;
 target_ptr = Libderiv->deriv2_classes[2][3][117];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+2172, dvrr_stack+4982, NULL);
 tmp = dvrr_stack + 2172;
 target_ptr = Libderiv->deriv2_classes[2][4][117];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+3302, dvrr_stack+6152, NULL);
 tmp = dvrr_stack + 3302;
 target_ptr = Libderiv->deriv2_classes[2][5][117];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,168,dvrr_stack+6149, dvrr_stack+7790, NULL);
 tmp = dvrr_stack + 6149;
 target_ptr = Libderiv->deriv2_classes[2][6][117];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+2262, dvrr_stack+673, dvrr_stack+10168);
 tmp = dvrr_stack + 2262;
 target_ptr = Libderiv->deriv2_classes[1][3][107];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+2292, dvrr_stack+10216, dvrr_stack+10186);
 tmp = dvrr_stack + 2292;
 target_ptr = Libderiv->deriv2_classes[1][4][107];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+601, dvrr_stack+10279, dvrr_stack+673);
 tmp = dvrr_stack + 601;
 target_ptr = Libderiv->deriv2_classes[1][5][107];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+3428, dvrr_stack+10363, dvrr_stack+10216);
 tmp = dvrr_stack + 3428;
 target_ptr = Libderiv->deriv2_classes[1][6][107];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+1234, dvrr_stack+5252, dvrr_stack+10615);
 tmp = dvrr_stack + 1234;
 target_ptr = Libderiv->deriv2_classes[2][3][107];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+8837, dvrr_stack+2364, dvrr_stack+10651);
 tmp = dvrr_stack + 8837;
 target_ptr = Libderiv->deriv2_classes[2][4][107];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+7781, dvrr_stack+10711, dvrr_stack+5252);
 tmp = dvrr_stack + 7781;
 target_ptr = Libderiv->deriv2_classes[2][5][107];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+7907, dvrr_stack+10879, dvrr_stack+2364);
 tmp = dvrr_stack + 7907;
 target_ptr = Libderiv->deriv2_classes[2][6][107];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+1555, dvrr_stack+9133, dvrr_stack+8294);
 tmp = dvrr_stack + 1555;
 target_ptr = Libderiv->deriv2_classes[1][3][106];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+15206, dvrr_stack+8342, dvrr_stack+8312);
 tmp = dvrr_stack + 15206;
 target_ptr = Libderiv->deriv2_classes[1][4][106];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+1294, dvrr_stack+8405, dvrr_stack+9133);
 tmp = dvrr_stack + 1294;
 target_ptr = Libderiv->deriv2_classes[1][5][106];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+2789, dvrr_stack+6530, dvrr_stack+8342);
 tmp = dvrr_stack + 2789;
 target_ptr = Libderiv->deriv2_classes[1][6][106];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+8075, dvrr_stack+3524, dvrr_stack+8489);
 tmp = dvrr_stack + 8075;
 target_ptr = Libderiv->deriv2_classes[2][3][106];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+8135, dvrr_stack+3674, dvrr_stack+3614);
 tmp = dvrr_stack + 8135;
 target_ptr = Libderiv->deriv2_classes[2][4][106];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+6317, dvrr_stack+13471, dvrr_stack+3524);
 tmp = dvrr_stack + 6317;
 target_ptr = Libderiv->deriv2_classes[2][5][106];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+4970, dvrr_stack+13639, dvrr_stack+3674);
 tmp = dvrr_stack + 4970;
 target_ptr = Libderiv->deriv2_classes[2][6][106];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+65, dvrr_stack+10114, dvrr_stack+8525);
 tmp = dvrr_stack + 65;
 target_ptr = Libderiv->deriv2_classes[1][3][105];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+4160, dvrr_stack+174, dvrr_stack+403);
 tmp = dvrr_stack + 4160;
 target_ptr = Libderiv->deriv2_classes[1][4][105];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+8225, dvrr_stack+1003, dvrr_stack+10114);
 tmp = dvrr_stack + 8225;
 target_ptr = Libderiv->deriv2_classes[1][5][105];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+6443, dvrr_stack+1860, dvrr_stack+174);
 tmp = dvrr_stack + 6443;
 target_ptr = Libderiv->deriv2_classes[1][6][105];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+5138, dvrr_stack+10507, dvrr_stack+3020);
 tmp = dvrr_stack + 5138;
 target_ptr = Libderiv->deriv2_classes[2][3][105];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+11095, dvrr_stack+4343, dvrr_stack+4712);
 tmp = dvrr_stack + 11095;
 target_ptr = Libderiv->deriv2_classes[2][4][105];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+11185, dvrr_stack+5774, dvrr_stack+10507);
 tmp = dvrr_stack + 11185;
 target_ptr = Libderiv->deriv2_classes[2][5][105];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+11311, dvrr_stack+7286, dvrr_stack+4343);
 tmp = dvrr_stack + 11311;
 target_ptr = Libderiv->deriv2_classes[2][6][105];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+22361, dvrr_stack+9448, dvrr_stack+10597);
 tmp = dvrr_stack + 22361;
 target_ptr = Libderiv->deriv2_classes[1][3][104];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+0, dvrr_stack+9493, dvrr_stack+237);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[1][4][104];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+11479, dvrr_stack+9556, dvrr_stack+9448);
 tmp = dvrr_stack + 11479;
 target_ptr = Libderiv->deriv2_classes[1][5][104];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_i(Data,3,1,dvrr_stack+11542, dvrr_stack+9640, dvrr_stack+9493);
 tmp = dvrr_stack + 11542;
 target_ptr = Libderiv->deriv2_classes[1][6][104];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+11626, dvrr_stack+9802, dvrr_stack+9766);
 tmp = dvrr_stack + 11626;
 target_ptr = Libderiv->deriv2_classes[2][3][104];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+11686, dvrr_stack+9952, dvrr_stack+9892);
 tmp = dvrr_stack + 11686;
 target_ptr = Libderiv->deriv2_classes[2][4][104];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+11776, dvrr_stack+7502, dvrr_stack+9802);
 tmp = dvrr_stack + 11776;
 target_ptr = Libderiv->deriv2_classes[2][5][104];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+11902, dvrr_stack+3056, dvrr_stack+9952);
 tmp = dvrr_stack + 11902;
 target_ptr = Libderiv->deriv2_classes[2][6][104];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+267, dvrr_stack+673, dvrr_stack+10168);
 tmp = dvrr_stack + 267;
 target_ptr = Libderiv->deriv2_classes[1][3][95];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+5198, dvrr_stack+10216, dvrr_stack+10186);
 tmp = dvrr_stack + 5198;
 target_ptr = Libderiv->deriv2_classes[1][4][95];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+12070, dvrr_stack+10279, dvrr_stack+673);
 tmp = dvrr_stack + 12070;
 target_ptr = Libderiv->deriv2_classes[1][5][95];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+12133, dvrr_stack+10363, dvrr_stack+10216);
 tmp = dvrr_stack + 12133;
 target_ptr = Libderiv->deriv2_classes[1][6][95];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+12217, dvrr_stack+5252, dvrr_stack+10615);
 tmp = dvrr_stack + 12217;
 target_ptr = Libderiv->deriv2_classes[2][3][95];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+12277, dvrr_stack+2364, dvrr_stack+10651);
 tmp = dvrr_stack + 12277;
 target_ptr = Libderiv->deriv2_classes[2][4][95];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+12367, dvrr_stack+10711, dvrr_stack+5252);
 tmp = dvrr_stack + 12367;
 target_ptr = Libderiv->deriv2_classes[2][5][95];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+12493, dvrr_stack+10879, dvrr_stack+2364);
 tmp = dvrr_stack + 12493;
 target_ptr = Libderiv->deriv2_classes[2][6][95];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+8927, dvrr_stack+9133, dvrr_stack+8294);
 tmp = dvrr_stack + 8927;
 target_ptr = Libderiv->deriv2_classes[1][3][94];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+12661, dvrr_stack+8342, dvrr_stack+8312);
 tmp = dvrr_stack + 12661;
 target_ptr = Libderiv->deriv2_classes[1][4][94];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+12706, dvrr_stack+8405, dvrr_stack+9133);
 tmp = dvrr_stack + 12706;
 target_ptr = Libderiv->deriv2_classes[1][5][94];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+12769, dvrr_stack+6530, dvrr_stack+8342);
 tmp = dvrr_stack + 12769;
 target_ptr = Libderiv->deriv2_classes[1][6][94];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+12853, dvrr_stack+3524, dvrr_stack+8489);
 tmp = dvrr_stack + 12853;
 target_ptr = Libderiv->deriv2_classes[2][3][94];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+12913, dvrr_stack+3674, dvrr_stack+3614);
 tmp = dvrr_stack + 12913;
 target_ptr = Libderiv->deriv2_classes[2][4][94];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+13003, dvrr_stack+13471, dvrr_stack+3524);
 tmp = dvrr_stack + 13003;
 target_ptr = Libderiv->deriv2_classes[2][5][94];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+13129, dvrr_stack+13639, dvrr_stack+3674);
 tmp = dvrr_stack + 13129;
 target_ptr = Libderiv->deriv2_classes[2][6][94];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+13297, dvrr_stack+10114, dvrr_stack+8525);
 tmp = dvrr_stack + 13297;
 target_ptr = Libderiv->deriv2_classes[1][3][93];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+13327, dvrr_stack+174, dvrr_stack+403);
 tmp = dvrr_stack + 13327;
 target_ptr = Libderiv->deriv2_classes[1][4][93];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+13372, dvrr_stack+1003, dvrr_stack+10114);
 tmp = dvrr_stack + 13372;
 target_ptr = Libderiv->deriv2_classes[1][5][93];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+3800, dvrr_stack+1860, dvrr_stack+174);
 tmp = dvrr_stack + 3800;
 target_ptr = Libderiv->deriv2_classes[1][6][93];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+3884, dvrr_stack+10507, dvrr_stack+3020);
 tmp = dvrr_stack + 3884;
 target_ptr = Libderiv->deriv2_classes[2][3][93];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+3944, dvrr_stack+4343, dvrr_stack+4712);
 tmp = dvrr_stack + 3944;
 target_ptr = Libderiv->deriv2_classes[2][4][93];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+24611, dvrr_stack+5774, dvrr_stack+10507);
 tmp = dvrr_stack + 24611;
 target_ptr = Libderiv->deriv2_classes[2][5][93];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+24737, dvrr_stack+7286, dvrr_stack+4343);
 tmp = dvrr_stack + 24737;
 target_ptr = Libderiv->deriv2_classes[2][6][93];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+13435, dvrr_stack+9448, dvrr_stack+10597);
 tmp = dvrr_stack + 13435;
 target_ptr = Libderiv->deriv2_classes[1][3][92];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+4034, dvrr_stack+9493, dvrr_stack+237);
 tmp = dvrr_stack + 4034;
 target_ptr = Libderiv->deriv2_classes[1][4][92];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+24905, dvrr_stack+9556, dvrr_stack+9448);
 tmp = dvrr_stack + 24905;
 target_ptr = Libderiv->deriv2_classes[1][5][92];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+24968, dvrr_stack+9640, dvrr_stack+9493);
 tmp = dvrr_stack + 24968;
 target_ptr = Libderiv->deriv2_classes[1][6][92];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+25052, dvrr_stack+9802, dvrr_stack+9766);
 tmp = dvrr_stack + 25052;
 target_ptr = Libderiv->deriv2_classes[2][3][92];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+25112, dvrr_stack+9952, dvrr_stack+9892);
 tmp = dvrr_stack + 25112;
 target_ptr = Libderiv->deriv2_classes[2][4][92];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+25202, dvrr_stack+7502, dvrr_stack+9802);
 tmp = dvrr_stack + 25202;
 target_ptr = Libderiv->deriv2_classes[2][5][92];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+25328, dvrr_stack+3056, dvrr_stack+9952);
 tmp = dvrr_stack + 25328;
 target_ptr = Libderiv->deriv2_classes[2][6][92];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+25496, dvrr_stack+7670, dvrr_stack+10078);
 tmp = dvrr_stack + 25496;
 target_ptr = Libderiv->deriv2_classes[1][3][91];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+25526, dvrr_stack+5942, dvrr_stack+7715);
 tmp = dvrr_stack + 25526;
 target_ptr = Libderiv->deriv2_classes[1][4][91];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+25571, dvrr_stack+6005, dvrr_stack+7670);
 tmp = dvrr_stack + 25571;
 target_ptr = Libderiv->deriv2_classes[1][5][91];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_i(Data,3,1,dvrr_stack+25634, dvrr_stack+4772, dvrr_stack+5942);
 tmp = dvrr_stack + 25634;
 target_ptr = Libderiv->deriv2_classes[1][6][91];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+25718, dvrr_stack+4880, dvrr_stack+7745);
 tmp = dvrr_stack + 25718;
 target_ptr = Libderiv->deriv2_classes[2][3][91];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+25778, dvrr_stack+1968, dvrr_stack+6089);
 tmp = dvrr_stack + 25778;
 target_ptr = Libderiv->deriv2_classes[2][4][91];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+25868, dvrr_stack+13855, dvrr_stack+4880);
 tmp = dvrr_stack + 25868;
 target_ptr = Libderiv->deriv2_classes[2][5][91];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+25994, dvrr_stack+14023, dvrr_stack+1968);
 tmp = dvrr_stack + 25994;
 target_ptr = Libderiv->deriv2_classes[2][6][91];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+26162, dvrr_stack+673, dvrr_stack+10168);
 tmp = dvrr_stack + 26162;
 target_ptr = Libderiv->deriv2_classes[1][3][83];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+26192, dvrr_stack+10216, dvrr_stack+10186);
 tmp = dvrr_stack + 26192;
 target_ptr = Libderiv->deriv2_classes[1][4][83];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+26237, dvrr_stack+10279, dvrr_stack+673);
 tmp = dvrr_stack + 26237;
 target_ptr = Libderiv->deriv2_classes[1][5][83];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+26300, dvrr_stack+10363, dvrr_stack+10216);
 tmp = dvrr_stack + 26300;
 target_ptr = Libderiv->deriv2_classes[1][6][83];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+10363, dvrr_stack+5252, dvrr_stack+10615);
 tmp = dvrr_stack + 10363;
 target_ptr = Libderiv->deriv2_classes[2][3][83];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+26384, dvrr_stack+2364, dvrr_stack+10651);
 tmp = dvrr_stack + 26384;
 target_ptr = Libderiv->deriv2_classes[2][4][83];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+26474, dvrr_stack+10711, dvrr_stack+5252);
 tmp = dvrr_stack + 26474;
 target_ptr = Libderiv->deriv2_classes[2][5][83];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+26600, dvrr_stack+10879, dvrr_stack+2364);
 tmp = dvrr_stack + 26600;
 target_ptr = Libderiv->deriv2_classes[2][6][83];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+10879, dvrr_stack+9133, dvrr_stack+8294);
 tmp = dvrr_stack + 10879;
 target_ptr = Libderiv->deriv2_classes[1][3][82];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+10909, dvrr_stack+8342, dvrr_stack+8312);
 tmp = dvrr_stack + 10909;
 target_ptr = Libderiv->deriv2_classes[1][4][82];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+10954, dvrr_stack+8405, dvrr_stack+9133);
 tmp = dvrr_stack + 10954;
 target_ptr = Libderiv->deriv2_classes[1][5][82];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+26768, dvrr_stack+6530, dvrr_stack+8342);
 tmp = dvrr_stack + 26768;
 target_ptr = Libderiv->deriv2_classes[1][6][82];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+11017, dvrr_stack+3524, dvrr_stack+8489);
 tmp = dvrr_stack + 11017;
 target_ptr = Libderiv->deriv2_classes[2][3][82];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+6527, dvrr_stack+3674, dvrr_stack+3614);
 tmp = dvrr_stack + 6527;
 target_ptr = Libderiv->deriv2_classes[2][4][82];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+26852, dvrr_stack+13471, dvrr_stack+3524);
 tmp = dvrr_stack + 26852;
 target_ptr = Libderiv->deriv2_classes[2][5][82];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+26978, dvrr_stack+13639, dvrr_stack+3674);
 tmp = dvrr_stack + 26978;
 target_ptr = Libderiv->deriv2_classes[2][6][82];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+13639, dvrr_stack+10114, dvrr_stack+8525);
 tmp = dvrr_stack + 13639;
 target_ptr = Libderiv->deriv2_classes[1][3][81];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+13669, dvrr_stack+174, dvrr_stack+403);
 tmp = dvrr_stack + 13669;
 target_ptr = Libderiv->deriv2_classes[1][4][81];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+13714, dvrr_stack+1003, dvrr_stack+10114);
 tmp = dvrr_stack + 13714;
 target_ptr = Libderiv->deriv2_classes[1][5][81];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+27146, dvrr_stack+1860, dvrr_stack+174);
 tmp = dvrr_stack + 27146;
 target_ptr = Libderiv->deriv2_classes[1][6][81];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+13777, dvrr_stack+10507, dvrr_stack+3020);
 tmp = dvrr_stack + 13777;
 target_ptr = Libderiv->deriv2_classes[2][3][81];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+1852, dvrr_stack+4343, dvrr_stack+4712);
 tmp = dvrr_stack + 1852;
 target_ptr = Libderiv->deriv2_classes[2][4][81];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+27230, dvrr_stack+5774, dvrr_stack+10507);
 tmp = dvrr_stack + 27230;
 target_ptr = Libderiv->deriv2_classes[2][5][81];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+27356, dvrr_stack+7286, dvrr_stack+4343);
 tmp = dvrr_stack + 27356;
 target_ptr = Libderiv->deriv2_classes[2][6][81];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+8489, dvrr_stack+9448, dvrr_stack+10597);
 tmp = dvrr_stack + 8489;
 target_ptr = Libderiv->deriv2_classes[1][3][80];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+3011, dvrr_stack+9493, dvrr_stack+237);
 tmp = dvrr_stack + 3011;
 target_ptr = Libderiv->deriv2_classes[1][4][80];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+7262, dvrr_stack+9556, dvrr_stack+9448);
 tmp = dvrr_stack + 7262;
 target_ptr = Libderiv->deriv2_classes[1][5][80];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+7325, dvrr_stack+9640, dvrr_stack+9493);
 tmp = dvrr_stack + 7325;
 target_ptr = Libderiv->deriv2_classes[1][6][80];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+9640, dvrr_stack+9802, dvrr_stack+9766);
 tmp = dvrr_stack + 9640;
 target_ptr = Libderiv->deriv2_classes[2][3][80];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+9700, dvrr_stack+9952, dvrr_stack+9892);
 tmp = dvrr_stack + 9700;
 target_ptr = Libderiv->deriv2_classes[2][4][80];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+27524, dvrr_stack+7502, dvrr_stack+9802);
 tmp = dvrr_stack + 27524;
 target_ptr = Libderiv->deriv2_classes[2][5][80];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+27650, dvrr_stack+3056, dvrr_stack+9952);
 tmp = dvrr_stack + 27650;
 target_ptr = Libderiv->deriv2_classes[2][6][80];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+3056, dvrr_stack+7670, dvrr_stack+10078);
 tmp = dvrr_stack + 3056;
 target_ptr = Libderiv->deriv2_classes[1][3][79];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+3086, dvrr_stack+5942, dvrr_stack+7715);
 tmp = dvrr_stack + 3086;
 target_ptr = Libderiv->deriv2_classes[1][4][79];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+3131, dvrr_stack+6005, dvrr_stack+7670);
 tmp = dvrr_stack + 3131;
 target_ptr = Libderiv->deriv2_classes[1][5][79];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+7409, dvrr_stack+4772, dvrr_stack+5942);
 tmp = dvrr_stack + 7409;
 target_ptr = Libderiv->deriv2_classes[1][6][79];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+4772, dvrr_stack+4880, dvrr_stack+7745);
 tmp = dvrr_stack + 4772;
 target_ptr = Libderiv->deriv2_classes[2][3][79];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+27818, dvrr_stack+1968, dvrr_stack+6089);
 tmp = dvrr_stack + 27818;
 target_ptr = Libderiv->deriv2_classes[2][4][79];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+27908, dvrr_stack+13855, dvrr_stack+4880);
 tmp = dvrr_stack + 27908;
 target_ptr = Libderiv->deriv2_classes[2][5][79];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+28034, dvrr_stack+14023, dvrr_stack+1968);
 tmp = dvrr_stack + 28034;
 target_ptr = Libderiv->deriv2_classes[2][6][79];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+14023, dvrr_stack+4469, dvrr_stack+2094);
 tmp = dvrr_stack + 14023;
 target_ptr = Libderiv->deriv2_classes[1][3][78];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+14053, dvrr_stack+1087, dvrr_stack+3272);
 tmp = dvrr_stack + 14053;
 target_ptr = Libderiv->deriv2_classes[1][4][78];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+14098, dvrr_stack+433, dvrr_stack+4469);
 tmp = dvrr_stack + 14098;
 target_ptr = Libderiv->deriv2_classes[1][5][78];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_i(Data,3,1,dvrr_stack+28202, dvrr_stack+14239, dvrr_stack+1087);
 tmp = dvrr_stack + 28202;
 target_ptr = Libderiv->deriv2_classes[1][6][78];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+14161, dvrr_stack+2921, dvrr_stack+2885);
 tmp = dvrr_stack + 14161;
 target_ptr = Libderiv->deriv2_classes[2][3][78];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+14221, dvrr_stack+14347, dvrr_stack+1381);
 tmp = dvrr_stack + 14221;
 target_ptr = Libderiv->deriv2_classes[2][4][78];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+28286, dvrr_stack+14473, dvrr_stack+2921);
 tmp = dvrr_stack + 28286;
 target_ptr = Libderiv->deriv2_classes[2][5][78];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+28412, dvrr_stack+14641, dvrr_stack+14347);
 tmp = dvrr_stack + 28412;
 target_ptr = Libderiv->deriv2_classes[2][6][78];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,10,dvrr_stack+14641, dvrr_stack+10651, dvrr_stack+8960);
 tmp = dvrr_stack + 14641;
 target_ptr = Libderiv->deriv2_classes[1][3][35];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,15,dvrr_stack+14671, dvrr_stack+5252, dvrr_stack+9253);
 tmp = dvrr_stack + 14671;
 target_ptr = Libderiv->deriv2_classes[1][4][35];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,21,dvrr_stack+14716, dvrr_stack+2364, dvrr_stack+517);
 tmp = dvrr_stack + 14716;
 target_ptr = Libderiv->deriv2_classes[1][5][35];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,28,dvrr_stack+28580, dvrr_stack+10711, dvrr_stack+9415);
 tmp = dvrr_stack + 28580;
 target_ptr = Libderiv->deriv2_classes[1][6][35];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,10,dvrr_stack+14779, dvrr_stack+15557, dvrr_stack+10186);
 tmp = dvrr_stack + 14779;
 target_ptr = Libderiv->deriv2_classes[2][3][35];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,15,dvrr_stack+28664, dvrr_stack+16506, dvrr_stack+673);
 tmp = dvrr_stack + 28664;
 target_ptr = Libderiv->deriv2_classes[2][4][35];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,21,dvrr_stack+28754, dvrr_stack+17818, dvrr_stack+10216);
 tmp = dvrr_stack + 28754;
 target_ptr = Libderiv->deriv2_classes[2][5][35];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,28,dvrr_stack+28880, dvrr_stack+6962, dvrr_stack+10279);
 tmp = dvrr_stack + 28880;
 target_ptr = Libderiv->deriv2_classes[2][6][35];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+14311, dvrr_stack+3614, dvrr_stack+7242);
 tmp = dvrr_stack + 14311;
 target_ptr = Libderiv->deriv2_classes[1][3][34];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+4832, dvrr_stack+3524, dvrr_stack+99);
 tmp = dvrr_stack + 4832;
 target_ptr = Libderiv->deriv2_classes[1][4][34];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+3194, dvrr_stack+3674, dvrr_stack+298);
 tmp = dvrr_stack + 3194;
 target_ptr = Libderiv->deriv2_classes[1][5][34];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+29048, dvrr_stack+13471, dvrr_stack+863);
 tmp = dvrr_stack + 29048;
 target_ptr = Libderiv->deriv2_classes[1][6][34];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+29132, dvrr_stack+1752, dvrr_stack+8312);
 tmp = dvrr_stack + 29132;
 target_ptr = Libderiv->deriv2_classes[2][3][34];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+29192, dvrr_stack+18028, dvrr_stack+9133);
 tmp = dvrr_stack + 29192;
 target_ptr = Libderiv->deriv2_classes[2][4][34];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+29282, dvrr_stack+19552, dvrr_stack+8342);
 tmp = dvrr_stack + 29282;
 target_ptr = Libderiv->deriv2_classes[2][5][34];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+29408, dvrr_stack+19762, dvrr_stack+8405);
 tmp = dvrr_stack + 29408;
 target_ptr = Libderiv->deriv2_classes[2][6][34];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+7745, dvrr_stack+4712, dvrr_stack+7252);
 tmp = dvrr_stack + 7745;
 target_ptr = Libderiv->deriv2_classes[1][3][33];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+10597, dvrr_stack+10507, dvrr_stack+9178);
 tmp = dvrr_stack + 10597;
 target_ptr = Libderiv->deriv2_classes[1][4][33];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+29576, dvrr_stack+4343, dvrr_stack+1441);
 tmp = dvrr_stack + 29576;
 target_ptr = Libderiv->deriv2_classes[1][5][33];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+29639, dvrr_stack+5774, dvrr_stack+9268);
 tmp = dvrr_stack + 29639;
 target_ptr = Libderiv->deriv2_classes[1][6][33];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+29723, dvrr_stack+9296, dvrr_stack+403);
 tmp = dvrr_stack + 29723;
 target_ptr = Libderiv->deriv2_classes[2][3][33];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+29783, dvrr_stack+15257, dvrr_stack+10114);
 tmp = dvrr_stack + 29783;
 target_ptr = Libderiv->deriv2_classes[2][4][33];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+29873, dvrr_stack+16056, dvrr_stack+174);
 tmp = dvrr_stack + 29873;
 target_ptr = Libderiv->deriv2_classes[2][5][33];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+29999, dvrr_stack+17188, dvrr_stack+1003);
 tmp = dvrr_stack + 29999;
 target_ptr = Libderiv->deriv2_classes[2][6][33];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+10423, dvrr_stack+9892, dvrr_stack+18712);
 tmp = dvrr_stack + 10423;
 target_ptr = Libderiv->deriv2_classes[1][3][32];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+2873, dvrr_stack+9802, dvrr_stack+18722);
 tmp = dvrr_stack + 2873;
 target_ptr = Libderiv->deriv2_classes[1][4][32];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+30167, dvrr_stack+9952, dvrr_stack+18737);
 tmp = dvrr_stack + 30167;
 target_ptr = Libderiv->deriv2_classes[1][5][32];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+30230, dvrr_stack+7502, dvrr_stack+18758);
 tmp = dvrr_stack + 30230;
 target_ptr = Libderiv->deriv2_classes[1][6][32];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+30314, dvrr_stack+18786, dvrr_stack+237);
 tmp = dvrr_stack + 30314;
 target_ptr = Libderiv->deriv2_classes[2][3][32];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+30374, dvrr_stack+15407, dvrr_stack+9448);
 tmp = dvrr_stack + 30374;
 target_ptr = Libderiv->deriv2_classes[2][4][32];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+30464, dvrr_stack+18886, dvrr_stack+9493);
 tmp = dvrr_stack + 30464;
 target_ptr = Libderiv->deriv2_classes[2][5][32];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+30590, dvrr_stack+19096, dvrr_stack+9556);
 tmp = dvrr_stack + 30590;
 target_ptr = Libderiv->deriv2_classes[2][6][32];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+30758, dvrr_stack+6089, dvrr_stack+19376);
 tmp = dvrr_stack + 30758;
 target_ptr = Libderiv->deriv2_classes[1][3][31];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+30788, dvrr_stack+4880, dvrr_stack+19386);
 tmp = dvrr_stack + 30788;
 target_ptr = Libderiv->deriv2_classes[1][4][31];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+30833, dvrr_stack+1968, dvrr_stack+19401);
 tmp = dvrr_stack + 30833;
 target_ptr = Libderiv->deriv2_classes[1][5][31];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+30896, dvrr_stack+13855, dvrr_stack+19422);
 tmp = dvrr_stack + 30896;
 target_ptr = Libderiv->deriv2_classes[1][6][31];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+30980, dvrr_stack+19450, dvrr_stack+7715);
 tmp = dvrr_stack + 30980;
 target_ptr = Libderiv->deriv2_classes[2][3][31];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+31040, dvrr_stack+17468, dvrr_stack+7670);
 tmp = dvrr_stack + 31040;
 target_ptr = Libderiv->deriv2_classes[2][4][31];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+31130, dvrr_stack+16266, dvrr_stack+5942);
 tmp = dvrr_stack + 31130;
 target_ptr = Libderiv->deriv2_classes[2][5][31];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+31256, dvrr_stack+20042, dvrr_stack+6005);
 tmp = dvrr_stack + 31256;
 target_ptr = Libderiv->deriv2_classes[2][6][31];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+31424, dvrr_stack+1381, dvrr_stack+10096);
 tmp = dvrr_stack + 31424;
 target_ptr = Libderiv->deriv2_classes[1][3][30];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+31454, dvrr_stack+2921, dvrr_stack+10471);
 tmp = dvrr_stack + 31454;
 target_ptr = Libderiv->deriv2_classes[1][4][30];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+31499, dvrr_stack+14347, dvrr_stack+10486);
 tmp = dvrr_stack + 31499;
 target_ptr = Libderiv->deriv2_classes[1][5][30];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+31562, dvrr_stack+14473, dvrr_stack+17618);
 tmp = dvrr_stack + 31562;
 target_ptr = Libderiv->deriv2_classes[1][6][30];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+31646, dvrr_stack+17646, dvrr_stack+3272);
 tmp = dvrr_stack + 31646;
 target_ptr = Libderiv->deriv2_classes[2][3][30];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+31706, dvrr_stack+18178, dvrr_stack+4469);
 tmp = dvrr_stack + 31706;
 target_ptr = Libderiv->deriv2_classes[2][4][30];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+31796, dvrr_stack+20322, dvrr_stack+1087);
 tmp = dvrr_stack + 31796;
 target_ptr = Libderiv->deriv2_classes[2][5][30];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+31922, dvrr_stack+20532, dvrr_stack+433);
 tmp = dvrr_stack + 31922;
 target_ptr = Libderiv->deriv2_classes[2][6][30];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+32090, dvrr_stack+17746, dvrr_stack+9193);
 tmp = dvrr_stack + 32090;
 target_ptr = Libderiv->deriv2_classes[1][3][26];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+32120, dvrr_stack+18328, dvrr_stack+9203);
 tmp = dvrr_stack + 32120;
 target_ptr = Libderiv->deriv2_classes[1][4][26];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+32165, dvrr_stack+18418, dvrr_stack+9218);
 tmp = dvrr_stack + 32165;
 target_ptr = Libderiv->deriv2_classes[1][5][26];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,28,dvrr_stack+32228, dvrr_stack+18544, dvrr_stack+16476);
 tmp = dvrr_stack + 32228;
 target_ptr = Libderiv->deriv2_classes[1][6][26];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+32312, dvrr_stack+21122, dvrr_stack+2759);
 tmp = dvrr_stack + 32312;
 target_ptr = Libderiv->deriv2_classes[2][3][26];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+32372, dvrr_stack+20812, dvrr_stack+1680);
 tmp = dvrr_stack + 32372;
 target_ptr = Libderiv->deriv2_classes[2][4][26];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+32462, dvrr_stack+21222, dvrr_stack+15017);
 tmp = dvrr_stack + 32462;
 target_ptr = Libderiv->deriv2_classes[2][5][26];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+32588, dvrr_stack+22401, dvrr_stack+15657);
 tmp = dvrr_stack + 32588;
 target_ptr = Libderiv->deriv2_classes[2][6][26];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,10,dvrr_stack+32756, dvrr_stack+10651, dvrr_stack+8960);
 tmp = dvrr_stack + 32756;
 target_ptr = Libderiv->deriv2_classes[1][3][23];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,15,dvrr_stack+32786, dvrr_stack+5252, dvrr_stack+9253);
 tmp = dvrr_stack + 32786;
 target_ptr = Libderiv->deriv2_classes[1][4][23];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,21,dvrr_stack+32831, dvrr_stack+2364, dvrr_stack+517);
 tmp = dvrr_stack + 32831;
 target_ptr = Libderiv->deriv2_classes[1][5][23];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,28,dvrr_stack+32894, dvrr_stack+10711, dvrr_stack+9415);
 tmp = dvrr_stack + 32894;
 target_ptr = Libderiv->deriv2_classes[1][6][23];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,10,dvrr_stack+32978, dvrr_stack+15557, dvrr_stack+10186);
 tmp = dvrr_stack + 32978;
 target_ptr = Libderiv->deriv2_classes[2][3][23];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,15,dvrr_stack+33038, dvrr_stack+16506, dvrr_stack+673);
 tmp = dvrr_stack + 33038;
 target_ptr = Libderiv->deriv2_classes[2][4][23];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,21,dvrr_stack+33128, dvrr_stack+17818, dvrr_stack+10216);
 tmp = dvrr_stack + 33128;
 target_ptr = Libderiv->deriv2_classes[2][5][23];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,28,dvrr_stack+33254, dvrr_stack+6962, dvrr_stack+10279);
 tmp = dvrr_stack + 33254;
 target_ptr = Libderiv->deriv2_classes[2][6][23];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+33422, dvrr_stack+3614, dvrr_stack+7242);
 tmp = dvrr_stack + 33422;
 target_ptr = Libderiv->deriv2_classes[1][3][22];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+33452, dvrr_stack+3524, dvrr_stack+99);
 tmp = dvrr_stack + 33452;
 target_ptr = Libderiv->deriv2_classes[1][4][22];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+33497, dvrr_stack+3674, dvrr_stack+298);
 tmp = dvrr_stack + 33497;
 target_ptr = Libderiv->deriv2_classes[1][5][22];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+33560, dvrr_stack+13471, dvrr_stack+863);
 tmp = dvrr_stack + 33560;
 target_ptr = Libderiv->deriv2_classes[1][6][22];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+33644, dvrr_stack+1752, dvrr_stack+8312);
 tmp = dvrr_stack + 33644;
 target_ptr = Libderiv->deriv2_classes[2][3][22];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+33704, dvrr_stack+18028, dvrr_stack+9133);
 tmp = dvrr_stack + 33704;
 target_ptr = Libderiv->deriv2_classes[2][4][22];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+33794, dvrr_stack+19552, dvrr_stack+8342);
 tmp = dvrr_stack + 33794;
 target_ptr = Libderiv->deriv2_classes[2][5][22];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+33920, dvrr_stack+19762, dvrr_stack+8405);
 tmp = dvrr_stack + 33920;
 target_ptr = Libderiv->deriv2_classes[2][6][22];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+34088, dvrr_stack+4712, dvrr_stack+7252);
 tmp = dvrr_stack + 34088;
 target_ptr = Libderiv->deriv2_classes[1][3][21];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+34118, dvrr_stack+10507, dvrr_stack+9178);
 tmp = dvrr_stack + 34118;
 target_ptr = Libderiv->deriv2_classes[1][4][21];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+34163, dvrr_stack+4343, dvrr_stack+1441);
 tmp = dvrr_stack + 34163;
 target_ptr = Libderiv->deriv2_classes[1][5][21];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+34226, dvrr_stack+5774, dvrr_stack+9268);
 tmp = dvrr_stack + 34226;
 target_ptr = Libderiv->deriv2_classes[1][6][21];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+34310, dvrr_stack+9296, dvrr_stack+403);
 tmp = dvrr_stack + 34310;
 target_ptr = Libderiv->deriv2_classes[2][3][21];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+34370, dvrr_stack+15257, dvrr_stack+10114);
 tmp = dvrr_stack + 34370;
 target_ptr = Libderiv->deriv2_classes[2][4][21];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+34460, dvrr_stack+16056, dvrr_stack+174);
 tmp = dvrr_stack + 34460;
 target_ptr = Libderiv->deriv2_classes[2][5][21];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+34586, dvrr_stack+17188, dvrr_stack+1003);
 tmp = dvrr_stack + 34586;
 target_ptr = Libderiv->deriv2_classes[2][6][21];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+34754, dvrr_stack+9892, dvrr_stack+18712);
 tmp = dvrr_stack + 34754;
 target_ptr = Libderiv->deriv2_classes[1][3][20];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+34784, dvrr_stack+9802, dvrr_stack+18722);
 tmp = dvrr_stack + 34784;
 target_ptr = Libderiv->deriv2_classes[1][4][20];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+34829, dvrr_stack+9952, dvrr_stack+18737);
 tmp = dvrr_stack + 34829;
 target_ptr = Libderiv->deriv2_classes[1][5][20];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+34892, dvrr_stack+7502, dvrr_stack+18758);
 tmp = dvrr_stack + 34892;
 target_ptr = Libderiv->deriv2_classes[1][6][20];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+34976, dvrr_stack+18786, dvrr_stack+237);
 tmp = dvrr_stack + 34976;
 target_ptr = Libderiv->deriv2_classes[2][3][20];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+35036, dvrr_stack+15407, dvrr_stack+9448);
 tmp = dvrr_stack + 35036;
 target_ptr = Libderiv->deriv2_classes[2][4][20];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+35126, dvrr_stack+18886, dvrr_stack+9493);
 tmp = dvrr_stack + 35126;
 target_ptr = Libderiv->deriv2_classes[2][5][20];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+35252, dvrr_stack+19096, dvrr_stack+9556);
 tmp = dvrr_stack + 35252;
 target_ptr = Libderiv->deriv2_classes[2][6][20];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+35420, dvrr_stack+6089, dvrr_stack+19376);
 tmp = dvrr_stack + 35420;
 target_ptr = Libderiv->deriv2_classes[1][3][19];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+35450, dvrr_stack+4880, dvrr_stack+19386);
 tmp = dvrr_stack + 35450;
 target_ptr = Libderiv->deriv2_classes[1][4][19];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+35495, dvrr_stack+1968, dvrr_stack+19401);
 tmp = dvrr_stack + 35495;
 target_ptr = Libderiv->deriv2_classes[1][5][19];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+35558, dvrr_stack+13855, dvrr_stack+19422);
 tmp = dvrr_stack + 35558;
 target_ptr = Libderiv->deriv2_classes[1][6][19];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+35642, dvrr_stack+19450, dvrr_stack+7715);
 tmp = dvrr_stack + 35642;
 target_ptr = Libderiv->deriv2_classes[2][3][19];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+35702, dvrr_stack+17468, dvrr_stack+7670);
 tmp = dvrr_stack + 35702;
 target_ptr = Libderiv->deriv2_classes[2][4][19];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+35792, dvrr_stack+16266, dvrr_stack+5942);
 tmp = dvrr_stack + 35792;
 target_ptr = Libderiv->deriv2_classes[2][5][19];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+35918, dvrr_stack+20042, dvrr_stack+6005);
 tmp = dvrr_stack + 35918;
 target_ptr = Libderiv->deriv2_classes[2][6][19];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+36086, dvrr_stack+1381, dvrr_stack+10096);
 tmp = dvrr_stack + 36086;
 target_ptr = Libderiv->deriv2_classes[1][3][18];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+36116, dvrr_stack+2921, dvrr_stack+10471);
 tmp = dvrr_stack + 36116;
 target_ptr = Libderiv->deriv2_classes[1][4][18];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+36161, dvrr_stack+14347, dvrr_stack+10486);
 tmp = dvrr_stack + 36161;
 target_ptr = Libderiv->deriv2_classes[1][5][18];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+36224, dvrr_stack+14473, dvrr_stack+17618);
 tmp = dvrr_stack + 36224;
 target_ptr = Libderiv->deriv2_classes[1][6][18];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+36308, dvrr_stack+17646, dvrr_stack+3272);
 tmp = dvrr_stack + 36308;
 target_ptr = Libderiv->deriv2_classes[2][3][18];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+36368, dvrr_stack+18178, dvrr_stack+4469);
 tmp = dvrr_stack + 36368;
 target_ptr = Libderiv->deriv2_classes[2][4][18];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+36458, dvrr_stack+20322, dvrr_stack+1087);
 tmp = dvrr_stack + 36458;
 target_ptr = Libderiv->deriv2_classes[2][5][18];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+36584, dvrr_stack+20532, dvrr_stack+433);
 tmp = dvrr_stack + 36584;
 target_ptr = Libderiv->deriv2_classes[2][6][18];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+36752, dvrr_stack+17746, dvrr_stack+9193);
 tmp = dvrr_stack + 36752;
 target_ptr = Libderiv->deriv2_classes[1][3][14];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+36782, dvrr_stack+18328, dvrr_stack+9203);
 tmp = dvrr_stack + 36782;
 target_ptr = Libderiv->deriv2_classes[1][4][14];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+36827, dvrr_stack+18418, dvrr_stack+9218);
 tmp = dvrr_stack + 36827;
 target_ptr = Libderiv->deriv2_classes[1][5][14];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+36890, dvrr_stack+18544, dvrr_stack+16476);
 tmp = dvrr_stack + 36890;
 target_ptr = Libderiv->deriv2_classes[1][6][14];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+36974, dvrr_stack+21122, dvrr_stack+2759);
 tmp = dvrr_stack + 36974;
 target_ptr = Libderiv->deriv2_classes[2][3][14];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+37034, dvrr_stack+20812, dvrr_stack+1680);
 tmp = dvrr_stack + 37034;
 target_ptr = Libderiv->deriv2_classes[2][4][14];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+37124, dvrr_stack+21222, dvrr_stack+15017);
 tmp = dvrr_stack + 37124;
 target_ptr = Libderiv->deriv2_classes[2][5][14];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+37250, dvrr_stack+22401, dvrr_stack+15657);
 tmp = dvrr_stack + 37250;
 target_ptr = Libderiv->deriv2_classes[2][6][14];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+37418, dvrr_stack+14857, dvrr_stack+20962);
 tmp = dvrr_stack + 37418;
 target_ptr = Libderiv->deriv2_classes[1][3][13];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+37448, dvrr_stack+16671, dvrr_stack+16656);
 tmp = dvrr_stack + 37448;
 target_ptr = Libderiv->deriv2_classes[1][4][13];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+37493, dvrr_stack+16782, dvrr_stack+16761);
 tmp = dvrr_stack + 37493;
 target_ptr = Libderiv->deriv2_classes[1][5][13];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,28,dvrr_stack+37556, dvrr_stack+22681, dvrr_stack+8970);
 tmp = dvrr_stack + 37556;
 target_ptr = Libderiv->deriv2_classes[1][6][13];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+37640, dvrr_stack+21687, dvrr_stack+21432);
 tmp = dvrr_stack + 37640;
 target_ptr = Libderiv->deriv2_classes[2][3][13];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+37700, dvrr_stack+21787, dvrr_stack+4208);
 tmp = dvrr_stack + 37700;
 target_ptr = Libderiv->deriv2_classes[2][4][13];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+37790, dvrr_stack+23269, dvrr_stack+5522);
 tmp = dvrr_stack + 37790;
 target_ptr = Libderiv->deriv2_classes[2][5][13];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+37916, dvrr_stack+23479, dvrr_stack+21937);
 tmp = dvrr_stack + 37916;
 target_ptr = Libderiv->deriv2_classes[2][6][13];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,10,dvrr_stack+38084, dvrr_stack+10651, dvrr_stack+8960);
 tmp = dvrr_stack + 38084;
 target_ptr = Libderiv->deriv2_classes[1][3][11];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,15,dvrr_stack+10642, dvrr_stack+5252, dvrr_stack+9253);
 tmp = dvrr_stack + 10642;
 target_ptr = Libderiv->deriv2_classes[1][4][11];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,21,dvrr_stack+5243, dvrr_stack+2364, dvrr_stack+517);
 tmp = dvrr_stack + 5243;
 target_ptr = Libderiv->deriv2_classes[1][5][11];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,28,dvrr_stack+2337, dvrr_stack+10711, dvrr_stack+9415);
 tmp = dvrr_stack + 2337;
 target_ptr = Libderiv->deriv2_classes[1][6][11];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,10,dvrr_stack+10687, dvrr_stack+15557, dvrr_stack+10186);
 tmp = dvrr_stack + 10687;
 target_ptr = Libderiv->deriv2_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,15,dvrr_stack+15557, dvrr_stack+16506, dvrr_stack+673);
 tmp = dvrr_stack + 15557;
 target_ptr = Libderiv->deriv2_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,21,dvrr_stack+10747, dvrr_stack+17818, dvrr_stack+10216);
 tmp = dvrr_stack + 10747;
 target_ptr = Libderiv->deriv2_classes[2][5][11];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,28,dvrr_stack+17806, dvrr_stack+6962, dvrr_stack+10279);
 tmp = dvrr_stack + 17806;
 target_ptr = Libderiv->deriv2_classes[2][6][11];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+6962, dvrr_stack+3614, dvrr_stack+7242);
 tmp = dvrr_stack + 6962;
 target_ptr = Libderiv->deriv2_classes[1][3][10];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3614, dvrr_stack+3524, dvrr_stack+99);
 tmp = dvrr_stack + 3614;
 target_ptr = Libderiv->deriv2_classes[1][4][10];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+6992, dvrr_stack+3674, dvrr_stack+298);
 tmp = dvrr_stack + 6992;
 target_ptr = Libderiv->deriv2_classes[1][5][10];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+3659, dvrr_stack+13471, dvrr_stack+863);
 tmp = dvrr_stack + 3659;
 target_ptr = Libderiv->deriv2_classes[1][6][10];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+7055, dvrr_stack+1752, dvrr_stack+8312);
 tmp = dvrr_stack + 7055;
 target_ptr = Libderiv->deriv2_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+7115, dvrr_stack+18028, dvrr_stack+9133);
 tmp = dvrr_stack + 7115;
 target_ptr = Libderiv->deriv2_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+1725, dvrr_stack+19552, dvrr_stack+8342);
 tmp = dvrr_stack + 1725;
 target_ptr = Libderiv->deriv2_classes[2][5][10];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+19550, dvrr_stack+19762, dvrr_stack+8405);
 tmp = dvrr_stack + 19550;
 target_ptr = Libderiv->deriv2_classes[2][6][10];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+3743, dvrr_stack+4712, dvrr_stack+7252);
 tmp = dvrr_stack + 3743;
 target_ptr = Libderiv->deriv2_classes[1][3][9];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+4712, dvrr_stack+10507, dvrr_stack+9178);
 tmp = dvrr_stack + 4712;
 target_ptr = Libderiv->deriv2_classes[1][4][9];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+10507, dvrr_stack+4343, dvrr_stack+1441);
 tmp = dvrr_stack + 10507;
 target_ptr = Libderiv->deriv2_classes[1][5][9];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+4316, dvrr_stack+5774, dvrr_stack+9268);
 tmp = dvrr_stack + 4316;
 target_ptr = Libderiv->deriv2_classes[1][6][9];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+2421, dvrr_stack+9296, dvrr_stack+403);
 tmp = dvrr_stack + 2421;
 target_ptr = Libderiv->deriv2_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+5766, dvrr_stack+15257, dvrr_stack+10114);
 tmp = dvrr_stack + 5766;
 target_ptr = Libderiv->deriv2_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+15251, dvrr_stack+16056, dvrr_stack+174);
 tmp = dvrr_stack + 15251;
 target_ptr = Libderiv->deriv2_classes[2][5][9];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+16056, dvrr_stack+17188, dvrr_stack+1003);
 tmp = dvrr_stack + 16056;
 target_ptr = Libderiv->deriv2_classes[2][6][9];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+15377, dvrr_stack+9892, dvrr_stack+18712);
 tmp = dvrr_stack + 15377;
 target_ptr = Libderiv->deriv2_classes[1][3][8];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+9892, dvrr_stack+9802, dvrr_stack+18722);
 tmp = dvrr_stack + 9892;
 target_ptr = Libderiv->deriv2_classes[1][4][8];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+17188, dvrr_stack+9952, dvrr_stack+18737);
 tmp = dvrr_stack + 17188;
 target_ptr = Libderiv->deriv2_classes[1][5][8];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+9937, dvrr_stack+7502, dvrr_stack+18758);
 tmp = dvrr_stack + 9937;
 target_ptr = Libderiv->deriv2_classes[1][6][8];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+18712, dvrr_stack+18786, dvrr_stack+237);
 tmp = dvrr_stack + 18712;
 target_ptr = Libderiv->deriv2_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+18772, dvrr_stack+15407, dvrr_stack+9448);
 tmp = dvrr_stack + 18772;
 target_ptr = Libderiv->deriv2_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+15407, dvrr_stack+18886, dvrr_stack+9493);
 tmp = dvrr_stack + 15407;
 target_ptr = Libderiv->deriv2_classes[2][5][8];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+18862, dvrr_stack+19096, dvrr_stack+9556);
 tmp = dvrr_stack + 18862;
 target_ptr = Libderiv->deriv2_classes[2][6][8];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+19030, dvrr_stack+6089, dvrr_stack+19376);
 tmp = dvrr_stack + 19030;
 target_ptr = Libderiv->deriv2_classes[1][3][7];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+6089, dvrr_stack+4880, dvrr_stack+19386);
 tmp = dvrr_stack + 6089;
 target_ptr = Libderiv->deriv2_classes[1][4][7];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+19060, dvrr_stack+1968, dvrr_stack+19401);
 tmp = dvrr_stack + 19060;
 target_ptr = Libderiv->deriv2_classes[1][5][7];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+19123, dvrr_stack+13855, dvrr_stack+19422);
 tmp = dvrr_stack + 19123;
 target_ptr = Libderiv->deriv2_classes[1][6][7];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+19207, dvrr_stack+19450, dvrr_stack+7715);
 tmp = dvrr_stack + 19207;
 target_ptr = Libderiv->deriv2_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+19267, dvrr_stack+17468, dvrr_stack+7670);
 tmp = dvrr_stack + 19267;
 target_ptr = Libderiv->deriv2_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+19357, dvrr_stack+16266, dvrr_stack+5942);
 tmp = dvrr_stack + 19357;
 target_ptr = Libderiv->deriv2_classes[2][5][7];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+17251, dvrr_stack+20042, dvrr_stack+6005);
 tmp = dvrr_stack + 17251;
 target_ptr = Libderiv->deriv2_classes[2][6][7];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+19483, dvrr_stack+1381, dvrr_stack+10096);
 tmp = dvrr_stack + 19483;
 target_ptr = Libderiv->deriv2_classes[1][3][6];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+10021, dvrr_stack+2921, dvrr_stack+10471);
 tmp = dvrr_stack + 10021;
 target_ptr = Libderiv->deriv2_classes[1][4][6];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+10066, dvrr_stack+14347, dvrr_stack+10486);
 tmp = dvrr_stack + 10066;
 target_ptr = Libderiv->deriv2_classes[1][5][6];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+10129, dvrr_stack+14473, dvrr_stack+17618);
 tmp = dvrr_stack + 10129;
 target_ptr = Libderiv->deriv2_classes[1][6][6];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+10213, dvrr_stack+17646, dvrr_stack+3272);
 tmp = dvrr_stack + 10213;
 target_ptr = Libderiv->deriv2_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+10273, dvrr_stack+18178, dvrr_stack+4469);
 tmp = dvrr_stack + 10273;
 target_ptr = Libderiv->deriv2_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+17419, dvrr_stack+20322, dvrr_stack+1087);
 tmp = dvrr_stack + 17419;
 target_ptr = Libderiv->deriv2_classes[2][5][6];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+17545, dvrr_stack+20532, dvrr_stack+433);
 tmp = dvrr_stack + 17545;
 target_ptr = Libderiv->deriv2_classes[2][6][6];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+19513, dvrr_stack+17746, dvrr_stack+9193);
 tmp = dvrr_stack + 19513;
 target_ptr = Libderiv->deriv2_classes[1][3][2];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+3257, dvrr_stack+18328, dvrr_stack+9203);
 tmp = dvrr_stack + 3257;
 target_ptr = Libderiv->deriv2_classes[1][4][2];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+17713, dvrr_stack+18418, dvrr_stack+9218);
 tmp = dvrr_stack + 17713;
 target_ptr = Libderiv->deriv2_classes[1][5][2];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+16224, dvrr_stack+18544, dvrr_stack+16476);
 tmp = dvrr_stack + 16224;
 target_ptr = Libderiv->deriv2_classes[1][6][2];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+16308, dvrr_stack+21122, dvrr_stack+2759);
 tmp = dvrr_stack + 16308;
 target_ptr = Libderiv->deriv2_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+16368, dvrr_stack+20812, dvrr_stack+1680);
 tmp = dvrr_stack + 16368;
 target_ptr = Libderiv->deriv2_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+16458, dvrr_stack+21222, dvrr_stack+15017);
 tmp = dvrr_stack + 16458;
 target_ptr = Libderiv->deriv2_classes[2][5][2];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+14341, dvrr_stack+22401, dvrr_stack+15657);
 tmp = dvrr_stack + 14341;
 target_ptr = Libderiv->deriv2_classes[2][6][2];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+17776, dvrr_stack+14857, dvrr_stack+20962);
 tmp = dvrr_stack + 17776;
 target_ptr = Libderiv->deriv2_classes[1][3][1];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+16584, dvrr_stack+16671, dvrr_stack+16656);
 tmp = dvrr_stack + 16584;
 target_ptr = Libderiv->deriv2_classes[1][4][1];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+16629, dvrr_stack+16782, dvrr_stack+16761);
 tmp = dvrr_stack + 16629;
 target_ptr = Libderiv->deriv2_classes[1][5][1];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+16692, dvrr_stack+22681, dvrr_stack+8970);
 tmp = dvrr_stack + 16692;
 target_ptr = Libderiv->deriv2_classes[1][6][1];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+16776, dvrr_stack+21687, dvrr_stack+21432);
 tmp = dvrr_stack + 16776;
 target_ptr = Libderiv->deriv2_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+15647, dvrr_stack+21787, dvrr_stack+4208);
 tmp = dvrr_stack + 15647;
 target_ptr = Libderiv->deriv2_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+14509, dvrr_stack+23269, dvrr_stack+5522);
 tmp = dvrr_stack + 14509;
 target_ptr = Libderiv->deriv2_classes[2][5][1];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+13837, dvrr_stack+23479, dvrr_stack+21937);
 tmp = dvrr_stack + 13837;
 target_ptr = Libderiv->deriv2_classes[2][6][1];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+16836, dvrr_stack+22021, dvrr_stack+5585);
 tmp = dvrr_stack + 16836;
 target_ptr = Libderiv->deriv2_classes[1][3][0];
 for(i=0;i<30;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+7205, dvrr_stack+14917, dvrr_stack+50);
 tmp = dvrr_stack + 7205;
 target_ptr = Libderiv->deriv2_classes[1][4][0];
 for(i=0;i<45;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+14839, dvrr_stack+15080, dvrr_stack+129);
 tmp = dvrr_stack + 14839;
 target_ptr = Libderiv->deriv2_classes[1][5][0];
 for(i=0;i<63;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,28,dvrr_stack+14902, dvrr_stack+15741, dvrr_stack+340);
 tmp = dvrr_stack + 14902;
 target_ptr = Libderiv->deriv2_classes[1][6][0];
 for(i=0;i<84;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+15737, dvrr_stack+16908, dvrr_stack+919);
 tmp = dvrr_stack + 15737;
 target_ptr = Libderiv->deriv2_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+16866, dvrr_stack+17008, dvrr_stack+20972);
 tmp = dvrr_stack + 16866;
 target_ptr = Libderiv->deriv2_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+16956, dvrr_stack+21462, dvrr_stack+4253);
 tmp = dvrr_stack + 16956;
 target_ptr = Libderiv->deriv2_classes[2][5][0];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+14986, dvrr_stack+22081, dvrr_stack+21017);
 tmp = dvrr_stack + 14986;
 target_ptr = Libderiv->deriv2_classes[2][6][0];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];


}

