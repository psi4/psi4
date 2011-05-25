#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (dd|fp) integrals */

void d12vrr_order_ddfp(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+0, dvrr_stack+3, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+12, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+6, NULL, NULL, dvrr_stack+0);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+39, dvrr_stack+15, dvrr_stack+6, dvrr_stack+12, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+49, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+52, dvrr_stack+49, dvrr_stack+12, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+58, dvrr_stack+52, dvrr_stack+15, dvrr_stack+49, dvrr_stack+12, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+68, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+71, dvrr_stack+3, dvrr_stack+68, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+77, dvrr_stack+6, dvrr_stack+71, dvrr_stack+0, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+87, dvrr_stack+39, dvrr_stack+77, NULL, NULL, dvrr_stack+6);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+117, dvrr_stack+58, dvrr_stack+39, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+147, dvrr_stack+117, dvrr_stack+87, dvrr_stack+58, dvrr_stack+39, dvrr_stack+21);
 tmp = dvrr_stack + 147;
 target_ptr = Libderiv->dvrr_classes[2][3];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+207, dvrr_stack+39, dvrr_stack+77, dvrr_stack+15, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+222, dvrr_stack+58, dvrr_stack+39, dvrr_stack+52, dvrr_stack+15, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+237, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+240, dvrr_stack+68, dvrr_stack+237, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+246, dvrr_stack+71, dvrr_stack+240, dvrr_stack+3, dvrr_stack+68, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+256, dvrr_stack+77, dvrr_stack+246, dvrr_stack+6, dvrr_stack+71, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+271, dvrr_stack+207, dvrr_stack+256, NULL, NULL, dvrr_stack+77);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+316, dvrr_stack+222, dvrr_stack+207, NULL, NULL, dvrr_stack+39);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+361, dvrr_stack+316, dvrr_stack+271, dvrr_stack+222, dvrr_stack+207, dvrr_stack+87);
 tmp = dvrr_stack + 361;
 target_ptr = Libderiv->dvrr_classes[2][4];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+451,dvrr_stack+361,dvrr_stack+147,6);


 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+631, dvrr_stack+207, dvrr_stack+256, dvrr_stack+39, dvrr_stack+77, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+652, dvrr_stack+222, dvrr_stack+207, dvrr_stack+58, dvrr_stack+39, NULL);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+673, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+676, dvrr_stack+237, dvrr_stack+673, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+682, dvrr_stack+240, dvrr_stack+676, dvrr_stack+68, dvrr_stack+237, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+692, dvrr_stack+246, dvrr_stack+682, dvrr_stack+71, dvrr_stack+240, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+707, dvrr_stack+256, dvrr_stack+692, dvrr_stack+77, dvrr_stack+246, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+728, dvrr_stack+631, dvrr_stack+707, NULL, NULL, dvrr_stack+256);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+791, dvrr_stack+652, dvrr_stack+631, NULL, NULL, dvrr_stack+207);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+854, dvrr_stack+791, dvrr_stack+728, dvrr_stack+652, dvrr_stack+631, dvrr_stack+271);

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+980,dvrr_stack+854,dvrr_stack+361,6);


 /* compute (2 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+1250,dvrr_stack+980,dvrr_stack+451,6);


 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,60,dvrr_stack+1610, dvrr_stack+1250, dvrr_stack+147);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1790, dvrr_stack+631, dvrr_stack+707, dvrr_stack+207, dvrr_stack+256, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1818, dvrr_stack+652, dvrr_stack+631, dvrr_stack+222, dvrr_stack+207, NULL);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+652, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+655, dvrr_stack+673, dvrr_stack+652, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+661, dvrr_stack+676, dvrr_stack+655, dvrr_stack+237, dvrr_stack+673, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1846, dvrr_stack+682, dvrr_stack+661, dvrr_stack+240, dvrr_stack+676, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1861, dvrr_stack+692, dvrr_stack+1846, dvrr_stack+246, dvrr_stack+682, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1882, dvrr_stack+707, dvrr_stack+1861, dvrr_stack+256, dvrr_stack+692, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1910, dvrr_stack+1790, dvrr_stack+1882, NULL, NULL, dvrr_stack+707);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1994, dvrr_stack+1818, dvrr_stack+1790, NULL, NULL, dvrr_stack+631);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2078, dvrr_stack+1994, dvrr_stack+1910, dvrr_stack+1818, dvrr_stack+1790, dvrr_stack+728);

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+2246,dvrr_stack+2078,dvrr_stack+854,6);


 /* compute (2 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+2624,dvrr_stack+2246,dvrr_stack+980,6);


 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,90,dvrr_stack+3164, dvrr_stack+2624, dvrr_stack+361);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1818, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+3);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+1827, dvrr_stack+6, dvrr_stack+71, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3434, dvrr_stack+21, dvrr_stack+1827, dvrr_stack+15, dvrr_stack+6, dvrr_stack+1818);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+3470, dvrr_stack+77, dvrr_stack+246, NULL, NULL, dvrr_stack+71);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+3500, dvrr_stack+87, dvrr_stack+3470, dvrr_stack+39, dvrr_stack+77, dvrr_stack+1827);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+3560, dvrr_stack+147, dvrr_stack+3500, dvrr_stack+117, dvrr_stack+87, dvrr_stack+3434);
 tmp = dvrr_stack + 3560;
 target_ptr = Libderiv->dvrr_classes[3][3];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3660, dvrr_stack+256, dvrr_stack+692, NULL, NULL, dvrr_stack+246);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3705, dvrr_stack+271, dvrr_stack+3660, dvrr_stack+207, dvrr_stack+256, dvrr_stack+3470);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3795, dvrr_stack+361, dvrr_stack+3705, dvrr_stack+316, dvrr_stack+271, dvrr_stack+3500);
 tmp = dvrr_stack + 3795;
 target_ptr = Libderiv->dvrr_classes[3][4];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+3945,dvrr_stack+3795,dvrr_stack+3560,10);


 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4245, dvrr_stack+707, dvrr_stack+1861, NULL, NULL, dvrr_stack+692);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4308, dvrr_stack+728, dvrr_stack+4245, dvrr_stack+631, dvrr_stack+707, dvrr_stack+3660);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4434, dvrr_stack+854, dvrr_stack+4308, dvrr_stack+791, dvrr_stack+728, dvrr_stack+3705);

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+4644,dvrr_stack+4434,dvrr_stack+3795,10);


 /* compute (3 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+5094,dvrr_stack+4644,dvrr_stack+3945,10);


 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,100,dvrr_stack+5694, dvrr_stack+5094, dvrr_stack+3560);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+631, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+634, dvrr_stack+652, dvrr_stack+631, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+39, dvrr_stack+655, dvrr_stack+634, dvrr_stack+673, dvrr_stack+652, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+207, dvrr_stack+661, dvrr_stack+39, dvrr_stack+676, dvrr_stack+655, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+5994, dvrr_stack+1846, dvrr_stack+207, dvrr_stack+682, dvrr_stack+661, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6015, dvrr_stack+1861, dvrr_stack+5994, dvrr_stack+692, dvrr_stack+1846, NULL);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6043, dvrr_stack+1882, dvrr_stack+6015, NULL, NULL, dvrr_stack+1861);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6127, dvrr_stack+1910, dvrr_stack+6043, dvrr_stack+1790, dvrr_stack+1882, dvrr_stack+4245);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6295, dvrr_stack+2078, dvrr_stack+6127, dvrr_stack+1994, dvrr_stack+1910, dvrr_stack+4308);

 /* compute (3 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+6575,dvrr_stack+6295,dvrr_stack+4434,10);


 /* compute (3 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+7205,dvrr_stack+6575,dvrr_stack+4644,10);


 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+8105, dvrr_stack+7205, dvrr_stack+3795);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+1994, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1997, dvrr_stack+3, dvrr_stack+68, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+2006, dvrr_stack+1818, dvrr_stack+1997, dvrr_stack+0, dvrr_stack+3, dvrr_stack+1994);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2024, dvrr_stack+71, dvrr_stack+240, NULL, NULL, dvrr_stack+68);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+2042, dvrr_stack+1827, dvrr_stack+2024, dvrr_stack+6, dvrr_stack+71, dvrr_stack+1997);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+8555, dvrr_stack+3434, dvrr_stack+2042, dvrr_stack+21, dvrr_stack+1827, dvrr_stack+2006);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+8615, dvrr_stack+246, dvrr_stack+682, NULL, NULL, dvrr_stack+240);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+8645, dvrr_stack+3470, dvrr_stack+8615, dvrr_stack+77, dvrr_stack+246, dvrr_stack+2024);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+8705, dvrr_stack+3500, dvrr_stack+8645, dvrr_stack+87, dvrr_stack+3470, dvrr_stack+2042);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+8805, dvrr_stack+3560, dvrr_stack+8705, dvrr_stack+147, dvrr_stack+3500, dvrr_stack+8555);
 tmp = dvrr_stack + 8805;
 target_ptr = Libderiv->dvrr_classes[4][3];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8955, dvrr_stack+692, dvrr_stack+1846, NULL, NULL, dvrr_stack+682);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+9000, dvrr_stack+3660, dvrr_stack+8955, dvrr_stack+256, dvrr_stack+692, dvrr_stack+8615);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+9090, dvrr_stack+3705, dvrr_stack+9000, dvrr_stack+271, dvrr_stack+3660, dvrr_stack+8645);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+9240, dvrr_stack+3795, dvrr_stack+9090, dvrr_stack+361, dvrr_stack+3705, dvrr_stack+8705);

 /* compute (4 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+9465,dvrr_stack+9240,dvrr_stack+8805,15);


 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+9915, dvrr_stack+1861, dvrr_stack+5994, NULL, NULL, dvrr_stack+1846);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+9978, dvrr_stack+4245, dvrr_stack+9915, dvrr_stack+707, dvrr_stack+1861, dvrr_stack+8955);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+10104, dvrr_stack+4308, dvrr_stack+9978, dvrr_stack+728, dvrr_stack+4245, dvrr_stack+9000);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+10314, dvrr_stack+4434, dvrr_stack+10104, dvrr_stack+854, dvrr_stack+4308, dvrr_stack+9090);

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+10629,dvrr_stack+10314,dvrr_stack+9240,15);


 /* compute (4 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+11304,dvrr_stack+10629,dvrr_stack+9465,15);


 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+12204, dvrr_stack+11304, dvrr_stack+8805);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+707, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+631, dvrr_stack+707, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+707, dvrr_stack+634, dvrr_stack+6, dvrr_stack+652, dvrr_stack+631, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+640, dvrr_stack+39, dvrr_stack+707, dvrr_stack+655, dvrr_stack+634, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+707, dvrr_stack+207, dvrr_stack+640, dvrr_stack+661, dvrr_stack+39, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1790, dvrr_stack+5994, dvrr_stack+707, dvrr_stack+1846, dvrr_stack+207, NULL);

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+12654, dvrr_stack+6015, dvrr_stack+1790, NULL, NULL, dvrr_stack+5994);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+12738, dvrr_stack+6043, dvrr_stack+12654, dvrr_stack+1882, dvrr_stack+6015, dvrr_stack+9915);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+12906, dvrr_stack+6127, dvrr_stack+12738, dvrr_stack+1910, dvrr_stack+6043, dvrr_stack+9978);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+13186, dvrr_stack+6295, dvrr_stack+12906, dvrr_stack+2078, dvrr_stack+6127, dvrr_stack+10104);

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+13606,dvrr_stack+13186,dvrr_stack+10314,15);


 /* compute (4 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+14551,dvrr_stack+13606,dvrr_stack+10629,15);


 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,225,dvrr_stack+15901, dvrr_stack+14551, dvrr_stack+9240);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,60,dvrr_stack+6015, dvrr_stack+1250, dvrr_stack+147);

 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,90,dvrr_stack+12654, dvrr_stack+2624, dvrr_stack+361);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,100,dvrr_stack+16576, dvrr_stack+5094, dvrr_stack+3560);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+16876, dvrr_stack+7205, dvrr_stack+3795);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+17326, dvrr_stack+11304, dvrr_stack+8805);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,225,dvrr_stack+17776, dvrr_stack+14551, dvrr_stack+9240);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,60,dvrr_stack+12924, dvrr_stack+1250, dvrr_stack+147);

 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,90,dvrr_stack+1250, dvrr_stack+2624, dvrr_stack+361);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,100,dvrr_stack+2624, dvrr_stack+5094, dvrr_stack+3560);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+5094, dvrr_stack+7205, dvrr_stack+3795);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+7205, dvrr_stack+11304, dvrr_stack+8805);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,225,dvrr_stack+11304, dvrr_stack+14551, dvrr_stack+9240);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+631, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+14551, dvrr_stack+52, dvrr_stack+15, NULL, NULL, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+14569, dvrr_stack+14551, dvrr_stack+21, dvrr_stack+52, dvrr_stack+15, dvrr_stack+631);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+14605,dvrr_stack+147,dvrr_stack+14569,6);


 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+14713, dvrr_stack+14605, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+1520, dvrr_stack+980, NULL);
 tmp = dvrr_stack + 1520;
 target_ptr = Libderiv->deriv_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+256, dvrr_stack+451, NULL);
 tmp = dvrr_stack + 256;
 target_ptr = Libderiv->deriv_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+14749, dvrr_stack+2246, NULL);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+15, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+14875, dvrr_stack+631, dvrr_stack+1818, dvrr_stack+12, dvrr_stack+0, dvrr_stack+15);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+14893, dvrr_stack+14569, dvrr_stack+3434, dvrr_stack+14551, dvrr_stack+21, dvrr_stack+14875);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+14953,dvrr_stack+3560,dvrr_stack+14893,10);


 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+15133, dvrr_stack+14953, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+5544, dvrr_stack+4644, NULL);
 tmp = dvrr_stack + 5544;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+6195, dvrr_stack+3945, NULL);
 tmp = dvrr_stack + 6195;
 target_ptr = Libderiv->deriv_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+15193, dvrr_stack+6575, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+52, dvrr_stack+15, dvrr_stack+1994, Data->F+2, Data->F+3, NULL);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+15403, dvrr_stack+14875, dvrr_stack+2006, dvrr_stack+631, dvrr_stack+1818, dvrr_stack+52);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+15433, dvrr_stack+14893, dvrr_stack+8555, dvrr_stack+14569, dvrr_stack+3434, dvrr_stack+15403);

 /* compute (4 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+15523,dvrr_stack+8805,dvrr_stack+15433,15);


 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+15793, dvrr_stack+15523, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+11979, dvrr_stack+10629, NULL);
 tmp = dvrr_stack + 11979;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+7655, dvrr_stack+9465, NULL);
 tmp = dvrr_stack + 7655;
 target_ptr = Libderiv->deriv_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+18451, dvrr_stack+13606, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+7805, dvrr_stack+14605, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+7841, dvrr_stack+980, NULL);
 tmp = dvrr_stack + 7841;
 target_ptr = Libderiv->deriv_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+7931, dvrr_stack+451, NULL);
 tmp = dvrr_stack + 7931;
 target_ptr = Libderiv->deriv_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+2924, dvrr_stack+2246, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+7991, dvrr_stack+14953, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+18766, dvrr_stack+4644, NULL);
 tmp = dvrr_stack + 18766;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+3050, dvrr_stack+3945, NULL);
 tmp = dvrr_stack + 3050;
 target_ptr = Libderiv->deriv_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+18916, dvrr_stack+6575, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+1882, dvrr_stack+15523, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+19126, dvrr_stack+10629, NULL);
 tmp = dvrr_stack + 19126;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+19351, dvrr_stack+9465, NULL);
 tmp = dvrr_stack + 19351;
 target_ptr = Libderiv->deriv_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+19501, dvrr_stack+13606, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+8051, dvrr_stack+14605, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+14605, dvrr_stack+980, NULL);
 tmp = dvrr_stack + 14605;
 target_ptr = Libderiv->deriv_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+980, dvrr_stack+451, NULL);
 tmp = dvrr_stack + 980;
 target_ptr = Libderiv->deriv_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+451, dvrr_stack+2246, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+2246, dvrr_stack+14953, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+14953, dvrr_stack+4644, NULL);
 tmp = dvrr_stack + 14953;
 target_ptr = Libderiv->deriv_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+4644, dvrr_stack+3945, NULL);
 tmp = dvrr_stack + 4644;
 target_ptr = Libderiv->deriv_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+1040, dvrr_stack+6575, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+6575, dvrr_stack+15523, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+15523, dvrr_stack+10629, NULL);
 tmp = dvrr_stack + 15523;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+10629, dvrr_stack+9465, NULL);
 tmp = dvrr_stack + 10629;
 target_ptr = Libderiv->deriv_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+9465, dvrr_stack+13606, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+13606, dvrr_stack+49, dvrr_stack+12, NULL, NULL, Data->F+1);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+14695, dvrr_stack+13606, dvrr_stack+631, dvrr_stack+49, dvrr_stack+12, dvrr_stack+0);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+13615, dvrr_stack+147, dvrr_stack+14695);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+13651, dvrr_stack+854, dvrr_stack+147);
 tmp = dvrr_stack + 13651;
 target_ptr = Libderiv->deriv_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+13741, dvrr_stack+361, dvrr_stack+14569);
 tmp = dvrr_stack + 13741;
 target_ptr = Libderiv->deriv_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+13801, dvrr_stack+2078, dvrr_stack+361);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+13927, dvrr_stack+0, dvrr_stack+15, Data->F+1, Data->F+2, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+15103, dvrr_stack+14695, dvrr_stack+14875, dvrr_stack+13606, dvrr_stack+631, dvrr_stack+13927);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+13933, dvrr_stack+3560, dvrr_stack+15103);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+13993, dvrr_stack+4434, dvrr_stack+3560);
 tmp = dvrr_stack + 13993;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+14143, dvrr_stack+3795, dvrr_stack+14893);
 tmp = dvrr_stack + 14143;
 target_ptr = Libderiv->deriv_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+14243, dvrr_stack+6295, dvrr_stack+3795);

 /* compute (3 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+14453, dvrr_stack+13927, dvrr_stack+52, dvrr_stack+0, dvrr_stack+15, NULL);

 /* compute (4 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+15748, dvrr_stack+15103, dvrr_stack+15403, dvrr_stack+14695, dvrr_stack+14875, dvrr_stack+14453);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,15,1,dvrr_stack+14453, dvrr_stack+8805, dvrr_stack+15748);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+10779, dvrr_stack+10314, dvrr_stack+8805);
 tmp = dvrr_stack + 10779;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+11004, dvrr_stack+9240, dvrr_stack+15433);
 tmp = dvrr_stack + 11004;
 target_ptr = Libderiv->deriv_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+6665, dvrr_stack+13186, dvrr_stack+9240);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+9780, dvrr_stack+147, dvrr_stack+14695);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+9816, dvrr_stack+854, dvrr_stack+147);
 tmp = dvrr_stack + 9816;
 target_ptr = Libderiv->deriv_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+11154, dvrr_stack+361, dvrr_stack+14569);
 tmp = dvrr_stack + 11154;
 target_ptr = Libderiv->deriv_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+6980, dvrr_stack+2078, dvrr_stack+361);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+11214, dvrr_stack+3560, dvrr_stack+15103);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+3945, dvrr_stack+4434, dvrr_stack+3560);
 tmp = dvrr_stack + 3945;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+4095, dvrr_stack+3795, dvrr_stack+14893);
 tmp = dvrr_stack + 4095;
 target_ptr = Libderiv->deriv_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+4744, dvrr_stack+6295, dvrr_stack+3795);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,15,1,dvrr_stack+7106, dvrr_stack+8805, dvrr_stack+15748);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+2306, dvrr_stack+10314, dvrr_stack+8805);
 tmp = dvrr_stack + 2306;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+19816, dvrr_stack+9240, dvrr_stack+15433);
 tmp = dvrr_stack + 19816;
 target_ptr = Libderiv->deriv_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+19966, dvrr_stack+13186, dvrr_stack+9240);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+4195, dvrr_stack+147, dvrr_stack+14695);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+4954, dvrr_stack+854, dvrr_stack+147);
 tmp = dvrr_stack + 4954;
 target_ptr = Libderiv->deriv_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+854, dvrr_stack+361, dvrr_stack+14569);
 tmp = dvrr_stack + 854;
 target_ptr = Libderiv->deriv_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+20281, dvrr_stack+2078, dvrr_stack+361);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+2078, dvrr_stack+3560, dvrr_stack+15103);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+20407, dvrr_stack+4434, dvrr_stack+3560);
 tmp = dvrr_stack + 20407;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+2138, dvrr_stack+3795, dvrr_stack+14893);
 tmp = dvrr_stack + 2138;
 target_ptr = Libderiv->deriv_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+20557, dvrr_stack+6295, dvrr_stack+3795);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,15,1,dvrr_stack+6295, dvrr_stack+8805, dvrr_stack+15748);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+20767, dvrr_stack+10314, dvrr_stack+8805);
 tmp = dvrr_stack + 20767;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+6385, dvrr_stack+9240, dvrr_stack+15433);
 tmp = dvrr_stack + 6385;
 target_ptr = Libderiv->deriv_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+20992, dvrr_stack+13186, dvrr_stack+9240);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+2531,dvrr_stack+316,dvrr_stack+117,3);


 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+15103, dvrr_stack+2531, NULL);

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+13104,dvrr_stack+791,dvrr_stack+316,3);


 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+15748, dvrr_stack+13104, NULL);

 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+13927, dvrr_stack+1994, dvrr_stack+0, Data->F+3, Data->F+4, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+9906, dvrr_stack+68, dvrr_stack+237, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+14695, dvrr_stack+1997, dvrr_stack+9906, dvrr_stack+3, dvrr_stack+68, dvrr_stack+0);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+11274, dvrr_stack+2006, dvrr_stack+14695, dvrr_stack+1818, dvrr_stack+1997, dvrr_stack+13927);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+15883, dvrr_stack+240, dvrr_stack+676, NULL, NULL, dvrr_stack+237);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+14569, dvrr_stack+2024, dvrr_stack+15883, dvrr_stack+71, dvrr_stack+240, dvrr_stack+9906);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+914, dvrr_stack+2042, dvrr_stack+14569, dvrr_stack+1827, dvrr_stack+2024, dvrr_stack+14695);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+13239, dvrr_stack+8555, dvrr_stack+914, dvrr_stack+3434, dvrr_stack+2042, dvrr_stack+11274);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+3434, dvrr_stack+682, dvrr_stack+661, NULL, NULL, dvrr_stack+676);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+577, dvrr_stack+8615, dvrr_stack+3434, dvrr_stack+246, dvrr_stack+682, dvrr_stack+15883);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+13329, dvrr_stack+8645, dvrr_stack+577, dvrr_stack+3470, dvrr_stack+8615, dvrr_stack+14569);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+13429, dvrr_stack+8705, dvrr_stack+13329, dvrr_stack+3500, dvrr_stack+8645, dvrr_stack+914);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+21307, dvrr_stack+8805, dvrr_stack+13429, dvrr_stack+3560, dvrr_stack+8705, dvrr_stack+13239);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3464, dvrr_stack+1846, dvrr_stack+207, NULL, NULL, dvrr_stack+661);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21517, dvrr_stack+8955, dvrr_stack+3464, dvrr_stack+692, dvrr_stack+1846, dvrr_stack+3434);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21607, dvrr_stack+9000, dvrr_stack+21517, dvrr_stack+3660, dvrr_stack+8955, dvrr_stack+577);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21757, dvrr_stack+9090, dvrr_stack+21607, dvrr_stack+3705, dvrr_stack+9000, dvrr_stack+13329);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21982, dvrr_stack+9240, dvrr_stack+21757, dvrr_stack+3795, dvrr_stack+9090, dvrr_stack+13429);

 /* compute (5 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+22297,dvrr_stack+21982,dvrr_stack+21307,21);


 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+22927, dvrr_stack+22297, NULL);

 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+728, dvrr_stack+5994, dvrr_stack+707, NULL, NULL, dvrr_stack+207);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3660, dvrr_stack+9915, dvrr_stack+728, dvrr_stack+1861, dvrr_stack+5994, dvrr_stack+3464);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+23137, dvrr_stack+9978, dvrr_stack+3660, dvrr_stack+4245, dvrr_stack+9915, dvrr_stack+21517);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+23347, dvrr_stack+10104, dvrr_stack+23137, dvrr_stack+4308, dvrr_stack+9978, dvrr_stack+21607);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+23662, dvrr_stack+10314, dvrr_stack+23347, dvrr_stack+4434, dvrr_stack+10104, dvrr_stack+21757);

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+24103,dvrr_stack+23662,dvrr_stack+21982,21);


 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+23137, dvrr_stack+24103, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+23452, dvrr_stack+2531, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+23482, dvrr_stack+13104, NULL);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+9915, dvrr_stack+22297, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+10125, dvrr_stack+24103, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+23527, dvrr_stack+2531, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+2531, dvrr_stack+13104, NULL);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+4231, dvrr_stack+22297, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+22297, dvrr_stack+24103, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+24103, dvrr_stack+316, dvrr_stack+14551);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+24133, dvrr_stack+791, dvrr_stack+117);

 /* compute (3 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+246, dvrr_stack+52, dvrr_stack+13927, dvrr_stack+15, dvrr_stack+1994, NULL);

 /* compute (4 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+24178, dvrr_stack+15403, dvrr_stack+11274, dvrr_stack+14875, dvrr_stack+2006, dvrr_stack+246);

 /* compute (5 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+24223, dvrr_stack+15433, dvrr_stack+13239, dvrr_stack+14893, dvrr_stack+8555, dvrr_stack+24178);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+24349, dvrr_stack+21982, dvrr_stack+24223);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+22612, dvrr_stack+23662, dvrr_stack+21307);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+24178, dvrr_stack+316, dvrr_stack+14551);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+14875, dvrr_stack+791, dvrr_stack+117);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+24559, dvrr_stack+21982, dvrr_stack+24223);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+24769, dvrr_stack+23662, dvrr_stack+21307);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+14920, dvrr_stack+316, dvrr_stack+14551);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+15403, dvrr_stack+791, dvrr_stack+117);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+25084, dvrr_stack+21982, dvrr_stack+24223);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+25294, dvrr_stack+23662, dvrr_stack+21307);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+24208, dvrr_stack+147, dvrr_stack+58);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+24238, dvrr_stack+8805, dvrr_stack+147);
 tmp = dvrr_stack + 24238;
 target_ptr = Libderiv->deriv_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+15448, dvrr_stack+361, dvrr_stack+222);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+25609, dvrr_stack+9240, dvrr_stack+361);
 tmp = dvrr_stack + 25609;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+25759, dvrr_stack+3560, dvrr_stack+117);
 tmp = dvrr_stack + 25759;
 target_ptr = Libderiv->deriv_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+25819, dvrr_stack+21307, dvrr_stack+3560);
 tmp = dvrr_stack + 25819;
 target_ptr = Libderiv->deriv_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+13104, dvrr_stack+3795, dvrr_stack+316);
 tmp = dvrr_stack + 13104;
 target_ptr = Libderiv->deriv_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+23557, dvrr_stack+21982, dvrr_stack+3795);
 tmp = dvrr_stack + 23557;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+14950, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+974, dvrr_stack+0, dvrr_stack+14950, Data->F+4, Data->F+5, NULL);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+246, dvrr_stack+13927, dvrr_stack+974, dvrr_stack+1994, dvrr_stack+0, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+49, dvrr_stack+237, dvrr_stack+673, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+8087, dvrr_stack+9906, dvrr_stack+49, dvrr_stack+68, dvrr_stack+237, dvrr_stack+14950);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+15493, dvrr_stack+14695, dvrr_stack+8087, dvrr_stack+1997, dvrr_stack+9906, dvrr_stack+974);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+13194, dvrr_stack+11274, dvrr_stack+15493, dvrr_stack+2006, dvrr_stack+14695, dvrr_stack+246);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+14695, dvrr_stack+676, dvrr_stack+655, NULL, NULL, dvrr_stack+673);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+13579, dvrr_stack+15883, dvrr_stack+14695, dvrr_stack+240, dvrr_stack+676, dvrr_stack+49);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+23782, dvrr_stack+14569, dvrr_stack+13579, dvrr_stack+2024, dvrr_stack+15883, dvrr_stack+8087);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+23842, dvrr_stack+914, dvrr_stack+23782, dvrr_stack+2042, dvrr_stack+14569, dvrr_stack+15493);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+23932, dvrr_stack+13239, dvrr_stack+23842, dvrr_stack+8555, dvrr_stack+914, dvrr_stack+13194);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+15493, dvrr_stack+661, dvrr_stack+39, NULL, NULL, dvrr_stack+655);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+8555, dvrr_stack+3434, dvrr_stack+15493, dvrr_stack+682, dvrr_stack+661, dvrr_stack+14695);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+655, dvrr_stack+577, dvrr_stack+8555, dvrr_stack+8615, dvrr_stack+3434, dvrr_stack+13579);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10440, dvrr_stack+13329, dvrr_stack+655, dvrr_stack+8645, dvrr_stack+577, dvrr_stack+23782);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+25969, dvrr_stack+13429, dvrr_stack+10440, dvrr_stack+8705, dvrr_stack+13329, dvrr_stack+23842);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+26179, dvrr_stack+21307, dvrr_stack+25969, dvrr_stack+8805, dvrr_stack+13429, dvrr_stack+23932);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+23782, dvrr_stack+26179, dvrr_stack+8805);

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+23992, dvrr_stack+207, dvrr_stack+640, NULL, NULL, dvrr_stack+39);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8615, dvrr_stack+3464, dvrr_stack+23992, dvrr_stack+1846, dvrr_stack+207, dvrr_stack+15493);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+13194, dvrr_stack+21517, dvrr_stack+8615, dvrr_stack+8955, dvrr_stack+3464, dvrr_stack+8555);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8555, dvrr_stack+21607, dvrr_stack+13194, dvrr_stack+9000, dvrr_stack+21517, dvrr_stack+655);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+13194, dvrr_stack+21757, dvrr_stack+8555, dvrr_stack+9090, dvrr_stack+21607, dvrr_stack+10440);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+26459, dvrr_stack+21982, dvrr_stack+13194, dvrr_stack+9240, dvrr_stack+21757, dvrr_stack+25969);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+13194, dvrr_stack+26459, dvrr_stack+9240);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+15493, dvrr_stack+147, dvrr_stack+58);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+25969, dvrr_stack+8805, dvrr_stack+147);
 tmp = dvrr_stack + 25969;
 target_ptr = Libderiv->deriv_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+26069, dvrr_stack+361, dvrr_stack+222);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+10440, dvrr_stack+9240, dvrr_stack+361);
 tmp = dvrr_stack + 10440;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+26114, dvrr_stack+3560, dvrr_stack+117);
 tmp = dvrr_stack + 26114;
 target_ptr = Libderiv->deriv_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+8555, dvrr_stack+21307, dvrr_stack+3560);
 tmp = dvrr_stack + 8555;
 target_ptr = Libderiv->deriv_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+13509, dvrr_stack+3795, dvrr_stack+316);
 tmp = dvrr_stack + 13509;
 target_ptr = Libderiv->deriv_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+21517, dvrr_stack+21982, dvrr_stack+3795);
 tmp = dvrr_stack + 21517;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+21742, dvrr_stack+26179, dvrr_stack+8805);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+26879, dvrr_stack+26459, dvrr_stack+9240);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+21952, dvrr_stack+147, dvrr_stack+58);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+8705, dvrr_stack+8805, dvrr_stack+147);
 tmp = dvrr_stack + 8705;
 target_ptr = Libderiv->deriv_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+147, dvrr_stack+361, dvrr_stack+222);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+8955, dvrr_stack+9240, dvrr_stack+361);
 tmp = dvrr_stack + 8955;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+361, dvrr_stack+3560, dvrr_stack+117);
 tmp = dvrr_stack + 361;
 target_ptr = Libderiv->deriv_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+577, dvrr_stack+21307, dvrr_stack+3560);
 tmp = dvrr_stack + 577;
 target_ptr = Libderiv->deriv_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+21307, dvrr_stack+3795, dvrr_stack+316);
 tmp = dvrr_stack + 21307;
 target_ptr = Libderiv->deriv_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+3434, dvrr_stack+21982, dvrr_stack+3795);
 tmp = dvrr_stack + 3434;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+21982, dvrr_stack+26179, dvrr_stack+8805);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+27194, dvrr_stack+26459, dvrr_stack+9240);

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,60,dvrr_stack+8805, dvrr_stack+1610, NULL);
 tmp = dvrr_stack + 8805;
 target_ptr = Libderiv->deriv2_classes[2][3][143];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,90,dvrr_stack+8865, dvrr_stack+3164, NULL);
 tmp = dvrr_stack + 8865;
 target_ptr = Libderiv->deriv2_classes[2][4][143];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,100,dvrr_stack+22192, dvrr_stack+5694, NULL);
 tmp = dvrr_stack + 22192;
 target_ptr = Libderiv->deriv2_classes[3][3][143];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+26174, dvrr_stack+8105, NULL);
 tmp = dvrr_stack + 26174;
 target_ptr = Libderiv->deriv2_classes[3][4][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+26324, dvrr_stack+12204, NULL);
 tmp = dvrr_stack + 26324;
 target_ptr = Libderiv->deriv2_classes[4][3][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,225,dvrr_stack+26474, dvrr_stack+15901, NULL);
 tmp = dvrr_stack + 26474;
 target_ptr = Libderiv->deriv2_classes[4][4][143];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,60,dvrr_stack+21397, dvrr_stack+1610, NULL);
 tmp = dvrr_stack + 21397;
 target_ptr = Libderiv->deriv2_classes[2][3][131];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,90,dvrr_stack+26699, dvrr_stack+3164, NULL);
 tmp = dvrr_stack + 26699;
 target_ptr = Libderiv->deriv2_classes[2][4][131];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,100,dvrr_stack+9105, dvrr_stack+5694, NULL);
 tmp = dvrr_stack + 9105;
 target_ptr = Libderiv->deriv2_classes[3][3][131];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+9205, dvrr_stack+8105, NULL);
 tmp = dvrr_stack + 9205;
 target_ptr = Libderiv->deriv2_classes[3][4][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+3659, dvrr_stack+12204, NULL);
 tmp = dvrr_stack + 3659;
 target_ptr = Libderiv->deriv2_classes[4][3][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,225,dvrr_stack+27509, dvrr_stack+15901, NULL);
 tmp = dvrr_stack + 27509;
 target_ptr = Libderiv->deriv2_classes[4][4][131];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+21457, dvrr_stack+6015, NULL);
 tmp = dvrr_stack + 21457;
 target_ptr = Libderiv->deriv2_classes[2][3][130];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+26789, dvrr_stack+12654, NULL);
 tmp = dvrr_stack + 26789;
 target_ptr = Libderiv->deriv2_classes[2][4][130];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+9355, dvrr_stack+16576, NULL);
 tmp = dvrr_stack + 9355;
 target_ptr = Libderiv->deriv2_classes[3][3][130];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+27734, dvrr_stack+16876, NULL);
 tmp = dvrr_stack + 27734;
 target_ptr = Libderiv->deriv2_classes[3][4][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+4441, dvrr_stack+17326, NULL);
 tmp = dvrr_stack + 4441;
 target_ptr = Libderiv->deriv2_classes[4][3][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+27884, dvrr_stack+17776, NULL);
 tmp = dvrr_stack + 27884;
 target_ptr = Libderiv->deriv2_classes[4][4][130];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,60,dvrr_stack+192, dvrr_stack+1610, NULL);
 tmp = dvrr_stack + 192;
 target_ptr = Libderiv->deriv2_classes[2][3][119];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,90,dvrr_stack+1610, dvrr_stack+3164, NULL);
 tmp = dvrr_stack + 1610;
 target_ptr = Libderiv->deriv2_classes[2][4][119];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,100,dvrr_stack+1700, dvrr_stack+5694, NULL);
 tmp = dvrr_stack + 1700;
 target_ptr = Libderiv->deriv2_classes[3][3][119];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+5694, dvrr_stack+8105, NULL);
 tmp = dvrr_stack + 5694;
 target_ptr = Libderiv->deriv2_classes[3][4][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+5844, dvrr_stack+12204, NULL);
 tmp = dvrr_stack + 5844;
 target_ptr = Libderiv->deriv2_classes[4][3][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,225,dvrr_stack+12204, dvrr_stack+15901, NULL);
 tmp = dvrr_stack + 12204;
 target_ptr = Libderiv->deriv2_classes[4][4][119];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+12429, dvrr_stack+6015, NULL);
 tmp = dvrr_stack + 12429;
 target_ptr = Libderiv->deriv2_classes[2][3][118];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+12489, dvrr_stack+12654, NULL);
 tmp = dvrr_stack + 12489;
 target_ptr = Libderiv->deriv2_classes[2][4][118];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+12579, dvrr_stack+16576, NULL);
 tmp = dvrr_stack + 12579;
 target_ptr = Libderiv->deriv2_classes[3][3][118];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+12679, dvrr_stack+16876, NULL);
 tmp = dvrr_stack + 12679;
 target_ptr = Libderiv->deriv2_classes[3][4][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+5994, dvrr_stack+17326, NULL);
 tmp = dvrr_stack + 5994;
 target_ptr = Libderiv->deriv2_classes[4][3][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+8087, dvrr_stack+17776, NULL);
 tmp = dvrr_stack + 8087;
 target_ptr = Libderiv->deriv2_classes[4][4][118];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+12829, dvrr_stack+12924, NULL);
 tmp = dvrr_stack + 12829;
 target_ptr = Libderiv->deriv2_classes[2][3][117];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+12889, dvrr_stack+1250, NULL);
 tmp = dvrr_stack + 12889;
 target_ptr = Libderiv->deriv2_classes[2][4][117];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+1250, dvrr_stack+2624, NULL);
 tmp = dvrr_stack + 1250;
 target_ptr = Libderiv->deriv2_classes[3][3][117];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+1350, dvrr_stack+5094, NULL);
 tmp = dvrr_stack + 1350;
 target_ptr = Libderiv->deriv2_classes[3][4][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+8312, dvrr_stack+7205, NULL);
 tmp = dvrr_stack + 8312;
 target_ptr = Libderiv->deriv2_classes[4][3][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+15883, dvrr_stack+11304, NULL);
 tmp = dvrr_stack + 15883;
 target_ptr = Libderiv->deriv2_classes[4][4][117];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+12979, dvrr_stack+1520, dvrr_stack+14713);
 tmp = dvrr_stack + 12979;
 target_ptr = Libderiv->deriv2_classes[2][3][107];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+23992, dvrr_stack+14749, dvrr_stack+256);
 tmp = dvrr_stack + 23992;
 target_ptr = Libderiv->deriv2_classes[2][4][107];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+727, dvrr_stack+5544, dvrr_stack+15133);
 tmp = dvrr_stack + 727;
 target_ptr = Libderiv->deriv2_classes[3][3][107];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+16108, dvrr_stack+15193, dvrr_stack+6195);
 tmp = dvrr_stack + 16108;
 target_ptr = Libderiv->deriv2_classes[3][4][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+16258, dvrr_stack+11979, dvrr_stack+15793);
 tmp = dvrr_stack + 16258;
 target_ptr = Libderiv->deriv2_classes[4][3][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+16408, dvrr_stack+18451, dvrr_stack+7655);
 tmp = dvrr_stack + 16408;
 target_ptr = Libderiv->deriv2_classes[4][4][107];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+13039, dvrr_stack+7841, dvrr_stack+7805);
 tmp = dvrr_stack + 13039;
 target_ptr = Libderiv->deriv2_classes[2][3][106];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+3809, dvrr_stack+2924, dvrr_stack+7931);
 tmp = dvrr_stack + 3809;
 target_ptr = Libderiv->deriv2_classes[2][4][106];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+16633, dvrr_stack+18766, dvrr_stack+7991);
 tmp = dvrr_stack + 16633;
 target_ptr = Libderiv->deriv2_classes[3][3][106];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+16733, dvrr_stack+18916, dvrr_stack+3050);
 tmp = dvrr_stack + 16733;
 target_ptr = Libderiv->deriv2_classes[3][4][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+16883, dvrr_stack+19126, dvrr_stack+1882);
 tmp = dvrr_stack + 16883;
 target_ptr = Libderiv->deriv2_classes[4][3][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+17033, dvrr_stack+19501, dvrr_stack+19351);
 tmp = dvrr_stack + 17033;
 target_ptr = Libderiv->deriv2_classes[4][4][106];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+1800, dvrr_stack+14605, dvrr_stack+8051);
 tmp = dvrr_stack + 1800;
 target_ptr = Libderiv->deriv2_classes[2][3][105];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+8462, dvrr_stack+451, dvrr_stack+980);
 tmp = dvrr_stack + 8462;
 target_ptr = Libderiv->deriv2_classes[2][4][105];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+17258, dvrr_stack+14953, dvrr_stack+2246);
 tmp = dvrr_stack + 17258;
 target_ptr = Libderiv->deriv2_classes[3][3][105];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+17358, dvrr_stack+1040, dvrr_stack+4644);
 tmp = dvrr_stack + 17358;
 target_ptr = Libderiv->deriv2_classes[3][4][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+17508, dvrr_stack+15523, dvrr_stack+6575);
 tmp = dvrr_stack + 17508;
 target_ptr = Libderiv->deriv2_classes[4][3][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+17658, dvrr_stack+9465, dvrr_stack+10629);
 tmp = dvrr_stack + 17658;
 target_ptr = Libderiv->deriv2_classes[4][4][105];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+914, dvrr_stack+13651, dvrr_stack+13615);
 tmp = dvrr_stack + 914;
 target_ptr = Libderiv->deriv2_classes[2][3][104];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+17883, dvrr_stack+13801, dvrr_stack+13741);
 tmp = dvrr_stack + 17883;
 target_ptr = Libderiv->deriv2_classes[2][4][104];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+17973, dvrr_stack+13993, dvrr_stack+13933);
 tmp = dvrr_stack + 17973;
 target_ptr = Libderiv->deriv2_classes[3][3][104];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+18073, dvrr_stack+14243, dvrr_stack+14143);
 tmp = dvrr_stack + 18073;
 target_ptr = Libderiv->deriv2_classes[3][4][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+18223, dvrr_stack+10779, dvrr_stack+14453);
 tmp = dvrr_stack + 18223;
 target_ptr = Libderiv->deriv2_classes[4][3][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+11274, dvrr_stack+6665, dvrr_stack+11004);
 tmp = dvrr_stack + 11274;
 target_ptr = Libderiv->deriv2_classes[4][4][104];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+18373, dvrr_stack+1520, dvrr_stack+14713);
 tmp = dvrr_stack + 18373;
 target_ptr = Libderiv->deriv2_classes[2][3][95];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+11499, dvrr_stack+14749, dvrr_stack+256);
 tmp = dvrr_stack + 11499;
 target_ptr = Libderiv->deriv2_classes[2][4][95];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+11589, dvrr_stack+5544, dvrr_stack+15133);
 tmp = dvrr_stack + 11589;
 target_ptr = Libderiv->deriv2_classes[3][3][95];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+11689, dvrr_stack+15193, dvrr_stack+6195);
 tmp = dvrr_stack + 11689;
 target_ptr = Libderiv->deriv2_classes[3][4][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+28109, dvrr_stack+11979, dvrr_stack+15793);
 tmp = dvrr_stack + 28109;
 target_ptr = Libderiv->deriv2_classes[4][3][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+28259, dvrr_stack+18451, dvrr_stack+7655);
 tmp = dvrr_stack + 28259;
 target_ptr = Libderiv->deriv2_classes[4][4][95];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+11839, dvrr_stack+7841, dvrr_stack+7805);
 tmp = dvrr_stack + 11839;
 target_ptr = Libderiv->deriv2_classes[2][3][94];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+0, dvrr_stack+2924, dvrr_stack+7931);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[2][4][94];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+28484, dvrr_stack+18766, dvrr_stack+7991);
 tmp = dvrr_stack + 28484;
 target_ptr = Libderiv->deriv2_classes[3][3][94];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+28584, dvrr_stack+18916, dvrr_stack+3050);
 tmp = dvrr_stack + 28584;
 target_ptr = Libderiv->deriv2_classes[3][4][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+28734, dvrr_stack+19126, dvrr_stack+1882);
 tmp = dvrr_stack + 28734;
 target_ptr = Libderiv->deriv2_classes[4][3][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+2576, dvrr_stack+19501, dvrr_stack+19351);
 tmp = dvrr_stack + 2576;
 target_ptr = Libderiv->deriv2_classes[4][4][94];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+11899, dvrr_stack+14605, dvrr_stack+8051);
 tmp = dvrr_stack + 11899;
 target_ptr = Libderiv->deriv2_classes[2][3][93];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+28884, dvrr_stack+451, dvrr_stack+980);
 tmp = dvrr_stack + 28884;
 target_ptr = Libderiv->deriv2_classes[2][4][93];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+2801, dvrr_stack+14953, dvrr_stack+2246);
 tmp = dvrr_stack + 2801;
 target_ptr = Libderiv->deriv2_classes[3][3][93];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+7196, dvrr_stack+1040, dvrr_stack+4644);
 tmp = dvrr_stack + 7196;
 target_ptr = Libderiv->deriv2_classes[3][4][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+7346, dvrr_stack+15523, dvrr_stack+6575);
 tmp = dvrr_stack + 7346;
 target_ptr = Libderiv->deriv2_classes[4][3][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+5044, dvrr_stack+9465, dvrr_stack+10629);
 tmp = dvrr_stack + 5044;
 target_ptr = Libderiv->deriv2_classes[4][4][93];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+14543, dvrr_stack+13651, dvrr_stack+13615);
 tmp = dvrr_stack + 14543;
 target_ptr = Libderiv->deriv2_classes[2][3][92];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+7496, dvrr_stack+13801, dvrr_stack+13741);
 tmp = dvrr_stack + 7496;
 target_ptr = Libderiv->deriv2_classes[2][4][92];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+5269, dvrr_stack+13993, dvrr_stack+13933);
 tmp = dvrr_stack + 5269;
 target_ptr = Libderiv->deriv2_classes[3][3][92];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+5369, dvrr_stack+14243, dvrr_stack+14143);
 tmp = dvrr_stack + 5369;
 target_ptr = Libderiv->deriv2_classes[3][4][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+3150, dvrr_stack+10779, dvrr_stack+14453);
 tmp = dvrr_stack + 3150;
 target_ptr = Libderiv->deriv2_classes[4][3][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+28974, dvrr_stack+6665, dvrr_stack+11004);
 tmp = dvrr_stack + 28974;
 target_ptr = Libderiv->deriv2_classes[4][4][92];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+29199, dvrr_stack+9816, dvrr_stack+9780);
 tmp = dvrr_stack + 29199;
 target_ptr = Libderiv->deriv2_classes[2][3][91];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+29259, dvrr_stack+6980, dvrr_stack+11154);
 tmp = dvrr_stack + 29259;
 target_ptr = Libderiv->deriv2_classes[2][4][91];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+29349, dvrr_stack+3945, dvrr_stack+11214);
 tmp = dvrr_stack + 29349;
 target_ptr = Libderiv->deriv2_classes[3][3][91];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+29449, dvrr_stack+4744, dvrr_stack+4095);
 tmp = dvrr_stack + 29449;
 target_ptr = Libderiv->deriv2_classes[3][4][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+29599, dvrr_stack+2306, dvrr_stack+7106);
 tmp = dvrr_stack + 29599;
 target_ptr = Libderiv->deriv2_classes[4][3][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+29749, dvrr_stack+19966, dvrr_stack+19816);
 tmp = dvrr_stack + 29749;
 target_ptr = Libderiv->deriv2_classes[4][4][91];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+7586, dvrr_stack+1520, dvrr_stack+14713);
 tmp = dvrr_stack + 7586;
 target_ptr = Libderiv->deriv2_classes[2][3][83];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+3300, dvrr_stack+14749, dvrr_stack+256);
 tmp = dvrr_stack + 3300;
 target_ptr = Libderiv->deriv2_classes[2][4][83];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+14695, dvrr_stack+5544, dvrr_stack+15133);
 tmp = dvrr_stack + 14695;
 target_ptr = Libderiv->deriv2_classes[3][3][83];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+29974, dvrr_stack+15193, dvrr_stack+6195);
 tmp = dvrr_stack + 29974;
 target_ptr = Libderiv->deriv2_classes[3][4][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+15133, dvrr_stack+11979, dvrr_stack+15793);
 tmp = dvrr_stack + 15133;
 target_ptr = Libderiv->deriv2_classes[4][3][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+30124, dvrr_stack+18451, dvrr_stack+7655);
 tmp = dvrr_stack + 30124;
 target_ptr = Libderiv->deriv2_classes[4][4][83];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+15793, dvrr_stack+7841, dvrr_stack+7805);
 tmp = dvrr_stack + 15793;
 target_ptr = Libderiv->deriv2_classes[2][3][82];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+15283, dvrr_stack+2924, dvrr_stack+7931);
 tmp = dvrr_stack + 15283;
 target_ptr = Libderiv->deriv2_classes[2][4][82];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+18433, dvrr_stack+18766, dvrr_stack+7991);
 tmp = dvrr_stack + 18433;
 target_ptr = Libderiv->deriv2_classes[3][3][82];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+18533, dvrr_stack+18916, dvrr_stack+3050);
 tmp = dvrr_stack + 18533;
 target_ptr = Libderiv->deriv2_classes[3][4][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+18916, dvrr_stack+19126, dvrr_stack+1882);
 tmp = dvrr_stack + 18916;
 target_ptr = Libderiv->deriv2_classes[4][3][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+30349, dvrr_stack+19501, dvrr_stack+19351);
 tmp = dvrr_stack + 30349;
 target_ptr = Libderiv->deriv2_classes[4][4][82];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+19066, dvrr_stack+14605, dvrr_stack+8051);
 tmp = dvrr_stack + 19066;
 target_ptr = Libderiv->deriv2_classes[2][3][81];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+19501, dvrr_stack+451, dvrr_stack+980);
 tmp = dvrr_stack + 19501;
 target_ptr = Libderiv->deriv2_classes[2][4][81];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+19591, dvrr_stack+14953, dvrr_stack+2246);
 tmp = dvrr_stack + 19591;
 target_ptr = Libderiv->deriv2_classes[3][3][81];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+1860, dvrr_stack+1040, dvrr_stack+4644);
 tmp = dvrr_stack + 1860;
 target_ptr = Libderiv->deriv2_classes[3][4][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+1040, dvrr_stack+15523, dvrr_stack+6575);
 tmp = dvrr_stack + 1040;
 target_ptr = Libderiv->deriv2_classes[4][3][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+30574, dvrr_stack+9465, dvrr_stack+10629);
 tmp = dvrr_stack + 30574;
 target_ptr = Libderiv->deriv2_classes[4][4][81];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+1190, dvrr_stack+13651, dvrr_stack+13615);
 tmp = dvrr_stack + 1190;
 target_ptr = Libderiv->deriv2_classes[2][3][80];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+19691, dvrr_stack+13801, dvrr_stack+13741);
 tmp = dvrr_stack + 19691;
 target_ptr = Libderiv->deriv2_classes[2][4][80];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+13801, dvrr_stack+13993, dvrr_stack+13933);
 tmp = dvrr_stack + 13801;
 target_ptr = Libderiv->deriv2_classes[3][3][80];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+421, dvrr_stack+14243, dvrr_stack+14143);
 tmp = dvrr_stack + 421;
 target_ptr = Libderiv->deriv2_classes[3][4][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+14243, dvrr_stack+10779, dvrr_stack+14453);
 tmp = dvrr_stack + 14243;
 target_ptr = Libderiv->deriv2_classes[4][3][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+9455, dvrr_stack+6665, dvrr_stack+11004);
 tmp = dvrr_stack + 9455;
 target_ptr = Libderiv->deriv2_classes[4][4][80];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+14393, dvrr_stack+9816, dvrr_stack+9780);
 tmp = dvrr_stack + 14393;
 target_ptr = Libderiv->deriv2_classes[2][3][79];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+14453, dvrr_stack+6980, dvrr_stack+11154);
 tmp = dvrr_stack + 14453;
 target_ptr = Libderiv->deriv2_classes[2][4][79];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+9680, dvrr_stack+3945, dvrr_stack+11214);
 tmp = dvrr_stack + 9680;
 target_ptr = Libderiv->deriv2_classes[3][3][79];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+30799, dvrr_stack+4744, dvrr_stack+4095);
 tmp = dvrr_stack + 30799;
 target_ptr = Libderiv->deriv2_classes[3][4][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+4744, dvrr_stack+2306, dvrr_stack+7106);
 tmp = dvrr_stack + 4744;
 target_ptr = Libderiv->deriv2_classes[4][3][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+6535, dvrr_stack+19966, dvrr_stack+19816);
 tmp = dvrr_stack + 6535;
 target_ptr = Libderiv->deriv2_classes[4][4][79];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+4894, dvrr_stack+4954, dvrr_stack+4195);
 tmp = dvrr_stack + 4894;
 target_ptr = Libderiv->deriv2_classes[2][3][78];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+19966, dvrr_stack+20281, dvrr_stack+854);
 tmp = dvrr_stack + 19966;
 target_ptr = Libderiv->deriv2_classes[2][4][78];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+20056, dvrr_stack+20407, dvrr_stack+2078);
 tmp = dvrr_stack + 20056;
 target_ptr = Libderiv->deriv2_classes[3][3][78];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+20156, dvrr_stack+20557, dvrr_stack+2138);
 tmp = dvrr_stack + 20156;
 target_ptr = Libderiv->deriv2_classes[3][4][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+20557, dvrr_stack+20767, dvrr_stack+6295);
 tmp = dvrr_stack + 20557;
 target_ptr = Libderiv->deriv2_classes[4][3][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+6760, dvrr_stack+20992, dvrr_stack+6385);
 tmp = dvrr_stack + 6760;
 target_ptr = Libderiv->deriv2_classes[4][4][78];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,10,dvrr_stack+20707, dvrr_stack+6195, dvrr_stack+15103);
 tmp = dvrr_stack + 20707;
 target_ptr = Libderiv->deriv2_classes[2][3][35];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,15,dvrr_stack+6295, dvrr_stack+5544, dvrr_stack+15748);
 tmp = dvrr_stack + 6295;
 target_ptr = Libderiv->deriv2_classes[2][4][35];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,10,dvrr_stack+20992, dvrr_stack+7655, dvrr_stack+256);
 tmp = dvrr_stack + 20992;
 target_ptr = Libderiv->deriv2_classes[3][3][35];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,15,dvrr_stack+21092, dvrr_stack+11979, dvrr_stack+1520);
 tmp = dvrr_stack + 21092;
 target_ptr = Libderiv->deriv2_classes[3][4][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,10,dvrr_stack+6985, dvrr_stack+22927, dvrr_stack+6195);
 tmp = dvrr_stack + 6985;
 target_ptr = Libderiv->deriv2_classes[4][3][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,15,dvrr_stack+30949, dvrr_stack+23137, dvrr_stack+5544);
 tmp = dvrr_stack + 30949;
 target_ptr = Libderiv->deriv2_classes[4][4][35];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+11214, dvrr_stack+3050, dvrr_stack+23452);
 tmp = dvrr_stack + 11214;
 target_ptr = Libderiv->deriv2_classes[2][3][34];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+20306, dvrr_stack+18766, dvrr_stack+23482);
 tmp = dvrr_stack + 20306;
 target_ptr = Libderiv->deriv2_classes[2][4][34];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+2010, dvrr_stack+19351, dvrr_stack+7931);
 tmp = dvrr_stack + 2010;
 target_ptr = Libderiv->deriv2_classes[3][3][34];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+31174, dvrr_stack+19126, dvrr_stack+7841);
 tmp = dvrr_stack + 31174;
 target_ptr = Libderiv->deriv2_classes[3][4][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+31324, dvrr_stack+9915, dvrr_stack+3050);
 tmp = dvrr_stack + 31324;
 target_ptr = Libderiv->deriv2_classes[4][3][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+31474, dvrr_stack+10125, dvrr_stack+18766);
 tmp = dvrr_stack + 31474;
 target_ptr = Libderiv->deriv2_classes[4][4][34];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+21242, dvrr_stack+4644, dvrr_stack+23527);
 tmp = dvrr_stack + 21242;
 target_ptr = Libderiv->deriv2_classes[2][3][33];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+13901, dvrr_stack+14953, dvrr_stack+2531);
 tmp = dvrr_stack + 13901;
 target_ptr = Libderiv->deriv2_classes[2][4][33];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+31699, dvrr_stack+10629, dvrr_stack+980);
 tmp = dvrr_stack + 31699;
 target_ptr = Libderiv->deriv2_classes[3][3][33];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+31799, dvrr_stack+15523, dvrr_stack+14605);
 tmp = dvrr_stack + 31799;
 target_ptr = Libderiv->deriv2_classes[3][4][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+31949, dvrr_stack+4231, dvrr_stack+4644);
 tmp = dvrr_stack + 31949;
 target_ptr = Libderiv->deriv2_classes[4][3][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+32099, dvrr_stack+22297, dvrr_stack+14953);
 tmp = dvrr_stack + 32099;
 target_ptr = Libderiv->deriv2_classes[4][4][33];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+7991, dvrr_stack+14143, dvrr_stack+24103);
 tmp = dvrr_stack + 7991;
 target_ptr = Libderiv->deriv2_classes[2][3][32];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+32324, dvrr_stack+13993, dvrr_stack+24133);
 tmp = dvrr_stack + 32324;
 target_ptr = Libderiv->deriv2_classes[2][4][32];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+32414, dvrr_stack+11004, dvrr_stack+13741);
 tmp = dvrr_stack + 32414;
 target_ptr = Libderiv->deriv2_classes[3][3][32];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+32514, dvrr_stack+10779, dvrr_stack+13651);
 tmp = dvrr_stack + 32514;
 target_ptr = Libderiv->deriv2_classes[3][4][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+32664, dvrr_stack+24349, dvrr_stack+14143);
 tmp = dvrr_stack + 32664;
 target_ptr = Libderiv->deriv2_classes[4][3][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+32814, dvrr_stack+22612, dvrr_stack+13993);
 tmp = dvrr_stack + 32814;
 target_ptr = Libderiv->deriv2_classes[4][4][32];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+14795, dvrr_stack+4095, dvrr_stack+24178);
 tmp = dvrr_stack + 14795;
 target_ptr = Libderiv->deriv2_classes[2][3][31];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+33039, dvrr_stack+3945, dvrr_stack+14875);
 tmp = dvrr_stack + 33039;
 target_ptr = Libderiv->deriv2_classes[2][4][31];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+33129, dvrr_stack+19816, dvrr_stack+11154);
 tmp = dvrr_stack + 33129;
 target_ptr = Libderiv->deriv2_classes[3][3][31];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+33229, dvrr_stack+2306, dvrr_stack+9816);
 tmp = dvrr_stack + 33229;
 target_ptr = Libderiv->deriv2_classes[3][4][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+33379, dvrr_stack+24559, dvrr_stack+4095);
 tmp = dvrr_stack + 33379;
 target_ptr = Libderiv->deriv2_classes[4][3][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+33529, dvrr_stack+24769, dvrr_stack+3945);
 tmp = dvrr_stack + 33529;
 target_ptr = Libderiv->deriv2_classes[4][4][31];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+18683, dvrr_stack+2138, dvrr_stack+14920);
 tmp = dvrr_stack + 18683;
 target_ptr = Libderiv->deriv2_classes[2][3][30];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+33754, dvrr_stack+20407, dvrr_stack+15403);
 tmp = dvrr_stack + 33754;
 target_ptr = Libderiv->deriv2_classes[2][4][30];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+33844, dvrr_stack+6385, dvrr_stack+854);
 tmp = dvrr_stack + 33844;
 target_ptr = Libderiv->deriv2_classes[3][3][30];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+33944, dvrr_stack+20767, dvrr_stack+4954);
 tmp = dvrr_stack + 33944;
 target_ptr = Libderiv->deriv2_classes[3][4][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+34094, dvrr_stack+25084, dvrr_stack+2138);
 tmp = dvrr_stack + 34094;
 target_ptr = Libderiv->deriv2_classes[4][3][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+34244, dvrr_stack+25294, dvrr_stack+20407);
 tmp = dvrr_stack + 34244;
 target_ptr = Libderiv->deriv2_classes[4][4][30];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+34469, dvrr_stack+24238, dvrr_stack+24208);
 tmp = dvrr_stack + 34469;
 target_ptr = Libderiv->deriv2_classes[2][3][26];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+34529, dvrr_stack+25609, dvrr_stack+15448);
 tmp = dvrr_stack + 34529;
 target_ptr = Libderiv->deriv2_classes[2][4][26];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+34619, dvrr_stack+25819, dvrr_stack+25759);
 tmp = dvrr_stack + 34619;
 target_ptr = Libderiv->deriv2_classes[3][3][26];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+34719, dvrr_stack+23557, dvrr_stack+13104);
 tmp = dvrr_stack + 34719;
 target_ptr = Libderiv->deriv2_classes[3][4][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+34869, dvrr_stack+23782, dvrr_stack+24238);
 tmp = dvrr_stack + 34869;
 target_ptr = Libderiv->deriv2_classes[4][3][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+35019, dvrr_stack+13194, dvrr_stack+25609);
 tmp = dvrr_stack + 35019;
 target_ptr = Libderiv->deriv2_classes[4][4][26];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,10,dvrr_stack+35244, dvrr_stack+6195, dvrr_stack+15103);
 tmp = dvrr_stack + 35244;
 target_ptr = Libderiv->deriv2_classes[2][3][23];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,15,dvrr_stack+35304, dvrr_stack+5544, dvrr_stack+15748);
 tmp = dvrr_stack + 35304;
 target_ptr = Libderiv->deriv2_classes[2][4][23];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,10,dvrr_stack+35394, dvrr_stack+7655, dvrr_stack+256);
 tmp = dvrr_stack + 35394;
 target_ptr = Libderiv->deriv2_classes[3][3][23];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,15,dvrr_stack+35494, dvrr_stack+11979, dvrr_stack+1520);
 tmp = dvrr_stack + 35494;
 target_ptr = Libderiv->deriv2_classes[3][4][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,10,dvrr_stack+35644, dvrr_stack+22927, dvrr_stack+6195);
 tmp = dvrr_stack + 35644;
 target_ptr = Libderiv->deriv2_classes[4][3][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,15,dvrr_stack+35794, dvrr_stack+23137, dvrr_stack+5544);
 tmp = dvrr_stack + 35794;
 target_ptr = Libderiv->deriv2_classes[4][4][23];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+36019, dvrr_stack+3050, dvrr_stack+23452);
 tmp = dvrr_stack + 36019;
 target_ptr = Libderiv->deriv2_classes[2][3][22];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+36079, dvrr_stack+18766, dvrr_stack+23482);
 tmp = dvrr_stack + 36079;
 target_ptr = Libderiv->deriv2_classes[2][4][22];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+36169, dvrr_stack+19351, dvrr_stack+7931);
 tmp = dvrr_stack + 36169;
 target_ptr = Libderiv->deriv2_classes[3][3][22];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+36269, dvrr_stack+19126, dvrr_stack+7841);
 tmp = dvrr_stack + 36269;
 target_ptr = Libderiv->deriv2_classes[3][4][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+36419, dvrr_stack+9915, dvrr_stack+3050);
 tmp = dvrr_stack + 36419;
 target_ptr = Libderiv->deriv2_classes[4][3][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+36569, dvrr_stack+10125, dvrr_stack+18766);
 tmp = dvrr_stack + 36569;
 target_ptr = Libderiv->deriv2_classes[4][4][22];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+36794, dvrr_stack+4644, dvrr_stack+23527);
 tmp = dvrr_stack + 36794;
 target_ptr = Libderiv->deriv2_classes[2][3][21];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+36854, dvrr_stack+14953, dvrr_stack+2531);
 tmp = dvrr_stack + 36854;
 target_ptr = Libderiv->deriv2_classes[2][4][21];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+2901, dvrr_stack+10629, dvrr_stack+980);
 tmp = dvrr_stack + 2901;
 target_ptr = Libderiv->deriv2_classes[3][3][21];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+36944, dvrr_stack+15523, dvrr_stack+14605);
 tmp = dvrr_stack + 36944;
 target_ptr = Libderiv->deriv2_classes[3][4][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+37094, dvrr_stack+4231, dvrr_stack+4644);
 tmp = dvrr_stack + 37094;
 target_ptr = Libderiv->deriv2_classes[4][3][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+37244, dvrr_stack+22297, dvrr_stack+14953);
 tmp = dvrr_stack + 37244;
 target_ptr = Libderiv->deriv2_classes[4][4][21];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+37469, dvrr_stack+14143, dvrr_stack+24103);
 tmp = dvrr_stack + 37469;
 target_ptr = Libderiv->deriv2_classes[2][3][20];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+37529, dvrr_stack+13993, dvrr_stack+24133);
 tmp = dvrr_stack + 37529;
 target_ptr = Libderiv->deriv2_classes[2][4][20];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+37619, dvrr_stack+11004, dvrr_stack+13741);
 tmp = dvrr_stack + 37619;
 target_ptr = Libderiv->deriv2_classes[3][3][20];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+37719, dvrr_stack+10779, dvrr_stack+13651);
 tmp = dvrr_stack + 37719;
 target_ptr = Libderiv->deriv2_classes[3][4][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+37869, dvrr_stack+24349, dvrr_stack+14143);
 tmp = dvrr_stack + 37869;
 target_ptr = Libderiv->deriv2_classes[4][3][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+38019, dvrr_stack+22612, dvrr_stack+13993);
 tmp = dvrr_stack + 38019;
 target_ptr = Libderiv->deriv2_classes[4][4][20];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+38244, dvrr_stack+4095, dvrr_stack+24178);
 tmp = dvrr_stack + 38244;
 target_ptr = Libderiv->deriv2_classes[2][3][19];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+38304, dvrr_stack+3945, dvrr_stack+14875);
 tmp = dvrr_stack + 38304;
 target_ptr = Libderiv->deriv2_classes[2][4][19];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+38394, dvrr_stack+19816, dvrr_stack+11154);
 tmp = dvrr_stack + 38394;
 target_ptr = Libderiv->deriv2_classes[3][3][19];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+38494, dvrr_stack+2306, dvrr_stack+9816);
 tmp = dvrr_stack + 38494;
 target_ptr = Libderiv->deriv2_classes[3][4][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+38644, dvrr_stack+24559, dvrr_stack+4095);
 tmp = dvrr_stack + 38644;
 target_ptr = Libderiv->deriv2_classes[4][3][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+38794, dvrr_stack+24769, dvrr_stack+3945);
 tmp = dvrr_stack + 38794;
 target_ptr = Libderiv->deriv2_classes[4][4][19];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+39019, dvrr_stack+2138, dvrr_stack+14920);
 tmp = dvrr_stack + 39019;
 target_ptr = Libderiv->deriv2_classes[2][3][18];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+39079, dvrr_stack+20407, dvrr_stack+15403);
 tmp = dvrr_stack + 39079;
 target_ptr = Libderiv->deriv2_classes[2][4][18];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+39169, dvrr_stack+6385, dvrr_stack+854);
 tmp = dvrr_stack + 39169;
 target_ptr = Libderiv->deriv2_classes[3][3][18];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+39269, dvrr_stack+20767, dvrr_stack+4954);
 tmp = dvrr_stack + 39269;
 target_ptr = Libderiv->deriv2_classes[3][4][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+39419, dvrr_stack+25084, dvrr_stack+2138);
 tmp = dvrr_stack + 39419;
 target_ptr = Libderiv->deriv2_classes[4][3][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+39569, dvrr_stack+25294, dvrr_stack+20407);
 tmp = dvrr_stack + 39569;
 target_ptr = Libderiv->deriv2_classes[4][4][18];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+39794, dvrr_stack+24238, dvrr_stack+24208);
 tmp = dvrr_stack + 39794;
 target_ptr = Libderiv->deriv2_classes[2][3][14];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+39854, dvrr_stack+25609, dvrr_stack+15448);
 tmp = dvrr_stack + 39854;
 target_ptr = Libderiv->deriv2_classes[2][4][14];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+39944, dvrr_stack+25819, dvrr_stack+25759);
 tmp = dvrr_stack + 39944;
 target_ptr = Libderiv->deriv2_classes[3][3][14];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+40044, dvrr_stack+23557, dvrr_stack+13104);
 tmp = dvrr_stack + 40044;
 target_ptr = Libderiv->deriv2_classes[3][4][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+40194, dvrr_stack+23782, dvrr_stack+24238);
 tmp = dvrr_stack + 40194;
 target_ptr = Libderiv->deriv2_classes[4][3][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+40344, dvrr_stack+13194, dvrr_stack+25609);
 tmp = dvrr_stack + 40344;
 target_ptr = Libderiv->deriv2_classes[4][4][14];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+40569, dvrr_stack+25969, dvrr_stack+15493);
 tmp = dvrr_stack + 40569;
 target_ptr = Libderiv->deriv2_classes[2][3][13];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+40629, dvrr_stack+10440, dvrr_stack+26069);
 tmp = dvrr_stack + 40629;
 target_ptr = Libderiv->deriv2_classes[2][4][13];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+40719, dvrr_stack+8555, dvrr_stack+26114);
 tmp = dvrr_stack + 40719;
 target_ptr = Libderiv->deriv2_classes[3][3][13];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+40819, dvrr_stack+21517, dvrr_stack+13509);
 tmp = dvrr_stack + 40819;
 target_ptr = Libderiv->deriv2_classes[3][4][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+40969, dvrr_stack+21742, dvrr_stack+25969);
 tmp = dvrr_stack + 40969;
 target_ptr = Libderiv->deriv2_classes[4][3][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+41119, dvrr_stack+26879, dvrr_stack+10440);
 tmp = dvrr_stack + 41119;
 target_ptr = Libderiv->deriv2_classes[4][4][13];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,10,dvrr_stack+41344, dvrr_stack+6195, dvrr_stack+15103);
 tmp = dvrr_stack + 41344;
 target_ptr = Libderiv->deriv2_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,15,dvrr_stack+41404, dvrr_stack+5544, dvrr_stack+15748);
 tmp = dvrr_stack + 41404;
 target_ptr = Libderiv->deriv2_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,10,dvrr_stack+41494, dvrr_stack+7655, dvrr_stack+256);
 tmp = dvrr_stack + 41494;
 target_ptr = Libderiv->deriv2_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,15,dvrr_stack+41594, dvrr_stack+11979, dvrr_stack+1520);
 tmp = dvrr_stack + 41594;
 target_ptr = Libderiv->deriv2_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,10,dvrr_stack+11959, dvrr_stack+22927, dvrr_stack+6195);
 tmp = dvrr_stack + 11959;
 target_ptr = Libderiv->deriv2_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,15,dvrr_stack+41744, dvrr_stack+23137, dvrr_stack+5544);
 tmp = dvrr_stack + 41744;
 target_ptr = Libderiv->deriv2_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+22927, dvrr_stack+3050, dvrr_stack+23452);
 tmp = dvrr_stack + 22927;
 target_ptr = Libderiv->deriv2_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+22987, dvrr_stack+18766, dvrr_stack+23482);
 tmp = dvrr_stack + 22987;
 target_ptr = Libderiv->deriv2_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+23077, dvrr_stack+19351, dvrr_stack+7931);
 tmp = dvrr_stack + 23077;
 target_ptr = Libderiv->deriv2_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+19351, dvrr_stack+19126, dvrr_stack+7841);
 tmp = dvrr_stack + 19351;
 target_ptr = Libderiv->deriv2_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+19126, dvrr_stack+9915, dvrr_stack+3050);
 tmp = dvrr_stack + 19126;
 target_ptr = Libderiv->deriv2_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+23177, dvrr_stack+10125, dvrr_stack+18766);
 tmp = dvrr_stack + 23177;
 target_ptr = Libderiv->deriv2_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+19276, dvrr_stack+4644, dvrr_stack+23527);
 tmp = dvrr_stack + 19276;
 target_ptr = Libderiv->deriv2_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+23402, dvrr_stack+14953, dvrr_stack+2531);
 tmp = dvrr_stack + 23402;
 target_ptr = Libderiv->deriv2_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+1500, dvrr_stack+10629, dvrr_stack+980);
 tmp = dvrr_stack + 1500;
 target_ptr = Libderiv->deriv2_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+6144, dvrr_stack+15523, dvrr_stack+14605);
 tmp = dvrr_stack + 6144;
 target_ptr = Libderiv->deriv2_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+15523, dvrr_stack+4231, dvrr_stack+4644);
 tmp = dvrr_stack + 15523;
 target_ptr = Libderiv->deriv2_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+4195, dvrr_stack+22297, dvrr_stack+14953);
 tmp = dvrr_stack + 4195;
 target_ptr = Libderiv->deriv2_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+15673, dvrr_stack+14143, dvrr_stack+24103);
 tmp = dvrr_stack + 15673;
 target_ptr = Libderiv->deriv2_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+22292, dvrr_stack+13993, dvrr_stack+24133);
 tmp = dvrr_stack + 22292;
 target_ptr = Libderiv->deriv2_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+22382, dvrr_stack+11004, dvrr_stack+13741);
 tmp = dvrr_stack + 22382;
 target_ptr = Libderiv->deriv2_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+11004, dvrr_stack+10779, dvrr_stack+13651);
 tmp = dvrr_stack + 11004;
 target_ptr = Libderiv->deriv2_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+13599, dvrr_stack+24349, dvrr_stack+14143);
 tmp = dvrr_stack + 13599;
 target_ptr = Libderiv->deriv2_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+10590, dvrr_stack+22612, dvrr_stack+13993);
 tmp = dvrr_stack + 10590;
 target_ptr = Libderiv->deriv2_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+15733, dvrr_stack+4095, dvrr_stack+24178);
 tmp = dvrr_stack + 15733;
 target_ptr = Libderiv->deriv2_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+13991, dvrr_stack+3945, dvrr_stack+14875);
 tmp = dvrr_stack + 13991;
 target_ptr = Libderiv->deriv2_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+14081, dvrr_stack+19816, dvrr_stack+11154);
 tmp = dvrr_stack + 14081;
 target_ptr = Libderiv->deriv2_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+19781, dvrr_stack+2306, dvrr_stack+9816);
 tmp = dvrr_stack + 19781;
 target_ptr = Libderiv->deriv2_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+22482, dvrr_stack+24559, dvrr_stack+4095);
 tmp = dvrr_stack + 22482;
 target_ptr = Libderiv->deriv2_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+22632, dvrr_stack+24769, dvrr_stack+3945);
 tmp = dvrr_stack + 22632;
 target_ptr = Libderiv->deriv2_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+11154, dvrr_stack+2138, dvrr_stack+14920);
 tmp = dvrr_stack + 11154;
 target_ptr = Libderiv->deriv2_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+252, dvrr_stack+20407, dvrr_stack+15403);
 tmp = dvrr_stack + 252;
 target_ptr = Libderiv->deriv2_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+10815, dvrr_stack+6385, dvrr_stack+854);
 tmp = dvrr_stack + 10815;
 target_ptr = Libderiv->deriv2_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+6385, dvrr_stack+20767, dvrr_stack+4954);
 tmp = dvrr_stack + 6385;
 target_ptr = Libderiv->deriv2_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+20767, dvrr_stack+25084, dvrr_stack+2138);
 tmp = dvrr_stack + 20767;
 target_ptr = Libderiv->deriv2_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+2110, dvrr_stack+25294, dvrr_stack+20407);
 tmp = dvrr_stack + 2110;
 target_ptr = Libderiv->deriv2_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+4954, dvrr_stack+24238, dvrr_stack+24208);
 tmp = dvrr_stack + 4954;
 target_ptr = Libderiv->deriv2_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+20396, dvrr_stack+25609, dvrr_stack+15448);
 tmp = dvrr_stack + 20396;
 target_ptr = Libderiv->deriv2_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+15373, dvrr_stack+25819, dvrr_stack+25759);
 tmp = dvrr_stack + 15373;
 target_ptr = Libderiv->deriv2_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+25759, dvrr_stack+23557, dvrr_stack+13104);
 tmp = dvrr_stack + 25759;
 target_ptr = Libderiv->deriv2_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+23492, dvrr_stack+23782, dvrr_stack+24238);
 tmp = dvrr_stack + 23492;
 target_ptr = Libderiv->deriv2_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+23642, dvrr_stack+13194, dvrr_stack+25609);
 tmp = dvrr_stack + 23642;
 target_ptr = Libderiv->deriv2_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+25909, dvrr_stack+25969, dvrr_stack+15493);
 tmp = dvrr_stack + 25909;
 target_ptr = Libderiv->deriv2_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+23867, dvrr_stack+10440, dvrr_stack+26069);
 tmp = dvrr_stack + 23867;
 target_ptr = Libderiv->deriv2_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+13099, dvrr_stack+8555, dvrr_stack+26114);
 tmp = dvrr_stack + 13099;
 target_ptr = Libderiv->deriv2_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+13199, dvrr_stack+21517, dvrr_stack+13509);
 tmp = dvrr_stack + 13199;
 target_ptr = Libderiv->deriv2_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+21517, dvrr_stack+21742, dvrr_stack+25969);
 tmp = dvrr_stack + 21517;
 target_ptr = Libderiv->deriv2_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+21667, dvrr_stack+26879, dvrr_stack+10440);
 tmp = dvrr_stack + 21667;
 target_ptr = Libderiv->deriv2_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+21892, dvrr_stack+8705, dvrr_stack+21952);
 tmp = dvrr_stack + 21892;
 target_ptr = Libderiv->deriv2_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+26879, dvrr_stack+8955, dvrr_stack+147);
 tmp = dvrr_stack + 26879;
 target_ptr = Libderiv->deriv2_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+26969, dvrr_stack+577, dvrr_stack+361);
 tmp = dvrr_stack + 26969;
 target_ptr = Libderiv->deriv2_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+25969, dvrr_stack+3434, dvrr_stack+21307);
 tmp = dvrr_stack + 25969;
 target_ptr = Libderiv->deriv2_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+13349, dvrr_stack+21982, dvrr_stack+8705);
 tmp = dvrr_stack + 13349;
 target_ptr = Libderiv->deriv2_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+21952, dvrr_stack+27194, dvrr_stack+8955);
 tmp = dvrr_stack + 21952;
 target_ptr = Libderiv->deriv2_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];


}

