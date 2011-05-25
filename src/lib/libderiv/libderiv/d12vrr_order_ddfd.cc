#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (dd|fd) integrals */

void d12vrr_order_ddfd(Libderiv_t *Libderiv, prim_data *Data)
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
 tmp = dvrr_stack + 854;
 target_ptr = Libderiv->dvrr_classes[2][5];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

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
 _BUILD_00p0(Data,dvrr_stack+1846, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1849, dvrr_stack+673, dvrr_stack+1846, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1855, dvrr_stack+676, dvrr_stack+1849, dvrr_stack+237, dvrr_stack+673, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1865, dvrr_stack+682, dvrr_stack+1855, dvrr_stack+240, dvrr_stack+676, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1880, dvrr_stack+692, dvrr_stack+1865, dvrr_stack+246, dvrr_stack+682, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1901, dvrr_stack+707, dvrr_stack+1880, dvrr_stack+256, dvrr_stack+692, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1929, dvrr_stack+1790, dvrr_stack+1901, NULL, NULL, dvrr_stack+707);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2013, dvrr_stack+1818, dvrr_stack+1790, NULL, NULL, dvrr_stack+631);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2097, dvrr_stack+2013, dvrr_stack+1929, dvrr_stack+1818, dvrr_stack+1790, dvrr_stack+728);

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+2265,dvrr_stack+2097,dvrr_stack+854,6);


 /* compute (2 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+2643,dvrr_stack+2265,dvrr_stack+980,6);


 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,90,dvrr_stack+3183, dvrr_stack+2643, dvrr_stack+361);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3453, dvrr_stack+1790, dvrr_stack+1901, dvrr_stack+631, dvrr_stack+707, NULL);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3489, dvrr_stack+1818, dvrr_stack+1790, dvrr_stack+652, dvrr_stack+631, NULL);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1818, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1821, dvrr_stack+1846, dvrr_stack+1818, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1827, dvrr_stack+1849, dvrr_stack+1821, dvrr_stack+673, dvrr_stack+1846, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3525, dvrr_stack+1855, dvrr_stack+1827, dvrr_stack+676, dvrr_stack+1849, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3540, dvrr_stack+1865, dvrr_stack+3525, dvrr_stack+682, dvrr_stack+1855, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3561, dvrr_stack+1880, dvrr_stack+3540, dvrr_stack+692, dvrr_stack+1865, NULL);

 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3589, dvrr_stack+1901, dvrr_stack+3561, dvrr_stack+707, dvrr_stack+1880, NULL);

 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3625, dvrr_stack+3453, dvrr_stack+3589, NULL, NULL, dvrr_stack+1901);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3733, dvrr_stack+3489, dvrr_stack+3453, NULL, NULL, dvrr_stack+1790);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3841, dvrr_stack+3733, dvrr_stack+3625, dvrr_stack+3489, dvrr_stack+3453, dvrr_stack+1929);

 /* compute (2 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+4057,dvrr_stack+3841,dvrr_stack+2097,6);


 /* compute (2 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+4561,dvrr_stack+4057,dvrr_stack+2265,6);


 /* compute (2 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,126,dvrr_stack+5317, dvrr_stack+4561, dvrr_stack+854);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+1837, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+3);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+3489, dvrr_stack+6, dvrr_stack+71, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+5695, dvrr_stack+21, dvrr_stack+3489, dvrr_stack+15, dvrr_stack+6, dvrr_stack+1837);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+5731, dvrr_stack+77, dvrr_stack+246, NULL, NULL, dvrr_stack+71);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+5761, dvrr_stack+87, dvrr_stack+5731, dvrr_stack+39, dvrr_stack+77, dvrr_stack+3489);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+5821, dvrr_stack+147, dvrr_stack+5761, dvrr_stack+117, dvrr_stack+87, dvrr_stack+5695);
 tmp = dvrr_stack + 5821;
 target_ptr = Libderiv->dvrr_classes[3][3];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5921, dvrr_stack+256, dvrr_stack+692, NULL, NULL, dvrr_stack+246);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5966, dvrr_stack+271, dvrr_stack+5921, dvrr_stack+207, dvrr_stack+256, dvrr_stack+5731);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6056, dvrr_stack+361, dvrr_stack+5966, dvrr_stack+316, dvrr_stack+271, dvrr_stack+5761);
 tmp = dvrr_stack + 6056;
 target_ptr = Libderiv->dvrr_classes[3][4];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+6206,dvrr_stack+6056,dvrr_stack+5821,10);


 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6506, dvrr_stack+707, dvrr_stack+1880, NULL, NULL, dvrr_stack+692);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6569, dvrr_stack+728, dvrr_stack+6506, dvrr_stack+631, dvrr_stack+707, dvrr_stack+5921);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6695, dvrr_stack+854, dvrr_stack+6569, dvrr_stack+791, dvrr_stack+728, dvrr_stack+5966);
 tmp = dvrr_stack + 6695;
 target_ptr = Libderiv->dvrr_classes[3][5];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+6905,dvrr_stack+6695,dvrr_stack+6056,10);


 /* compute (3 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+7355,dvrr_stack+6905,dvrr_stack+6206,10);


 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,100,dvrr_stack+7955, dvrr_stack+7355, dvrr_stack+5821);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+8255, dvrr_stack+1901, dvrr_stack+3561, NULL, NULL, dvrr_stack+1880);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+8339, dvrr_stack+1929, dvrr_stack+8255, dvrr_stack+1790, dvrr_stack+1901, dvrr_stack+6506);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+8507, dvrr_stack+2097, dvrr_stack+8339, dvrr_stack+2013, dvrr_stack+1929, dvrr_stack+6569);

 /* compute (3 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+8787,dvrr_stack+8507,dvrr_stack+6695,10);


 /* compute (3 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+9417,dvrr_stack+8787,dvrr_stack+6905,10);


 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+10317, dvrr_stack+9417, dvrr_stack+6056);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1790, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1793, dvrr_stack+1818, dvrr_stack+1790, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+39, dvrr_stack+1821, dvrr_stack+1793, dvrr_stack+1846, dvrr_stack+1818, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+207, dvrr_stack+1827, dvrr_stack+39, dvrr_stack+1849, dvrr_stack+1821, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+631, dvrr_stack+3525, dvrr_stack+207, dvrr_stack+1855, dvrr_stack+1827, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+10767, dvrr_stack+3540, dvrr_stack+631, dvrr_stack+1865, dvrr_stack+3525, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+10795, dvrr_stack+3561, dvrr_stack+10767, dvrr_stack+1880, dvrr_stack+3540, NULL);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+10831, dvrr_stack+3589, dvrr_stack+10795, NULL, NULL, dvrr_stack+3561);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+10939, dvrr_stack+3625, dvrr_stack+10831, dvrr_stack+3453, dvrr_stack+3589, dvrr_stack+8255);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+11155, dvrr_stack+3841, dvrr_stack+10939, dvrr_stack+3733, dvrr_stack+3625, dvrr_stack+8339);

 /* compute (3 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+11515,dvrr_stack+11155,dvrr_stack+8507,10);


 /* compute (3 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+12355,dvrr_stack+11515,dvrr_stack+8787,10);


 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,210,dvrr_stack+13615, dvrr_stack+12355, dvrr_stack+6695);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+1846, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3733, dvrr_stack+3, dvrr_stack+68, NULL, NULL, Data->F+4);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+3507, dvrr_stack+1837, dvrr_stack+3733, dvrr_stack+0, dvrr_stack+3, dvrr_stack+1846);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+3742, dvrr_stack+71, dvrr_stack+240, NULL, NULL, dvrr_stack+68);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3453, dvrr_stack+3489, dvrr_stack+3742, dvrr_stack+6, dvrr_stack+71, dvrr_stack+3733);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+3760, dvrr_stack+5695, dvrr_stack+3453, dvrr_stack+21, dvrr_stack+3489, dvrr_stack+3507);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+14245, dvrr_stack+246, dvrr_stack+682, NULL, NULL, dvrr_stack+240);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+14275, dvrr_stack+5731, dvrr_stack+14245, dvrr_stack+77, dvrr_stack+246, dvrr_stack+3742);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+14335, dvrr_stack+5761, dvrr_stack+14275, dvrr_stack+87, dvrr_stack+5731, dvrr_stack+3453);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+14435, dvrr_stack+5821, dvrr_stack+14335, dvrr_stack+147, dvrr_stack+5761, dvrr_stack+3760);
 tmp = dvrr_stack + 14435;
 target_ptr = Libderiv->dvrr_classes[4][3];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14585, dvrr_stack+692, dvrr_stack+1865, NULL, NULL, dvrr_stack+682);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14630, dvrr_stack+5921, dvrr_stack+14585, dvrr_stack+256, dvrr_stack+692, dvrr_stack+14245);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14720, dvrr_stack+5966, dvrr_stack+14630, dvrr_stack+271, dvrr_stack+5921, dvrr_stack+14275);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14870, dvrr_stack+6056, dvrr_stack+14720, dvrr_stack+361, dvrr_stack+5966, dvrr_stack+14335);
 tmp = dvrr_stack + 14870;
 target_ptr = Libderiv->dvrr_classes[4][4];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+15095,dvrr_stack+14870,dvrr_stack+14435,15);


 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15545, dvrr_stack+1880, dvrr_stack+3540, NULL, NULL, dvrr_stack+1865);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15608, dvrr_stack+6506, dvrr_stack+15545, dvrr_stack+707, dvrr_stack+1880, dvrr_stack+14585);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15734, dvrr_stack+6569, dvrr_stack+15608, dvrr_stack+728, dvrr_stack+6506, dvrr_stack+14630);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15944, dvrr_stack+6695, dvrr_stack+15734, dvrr_stack+854, dvrr_stack+6569, dvrr_stack+14720);

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+16259,dvrr_stack+15944,dvrr_stack+14870,15);


 /* compute (4 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+16934,dvrr_stack+16259,dvrr_stack+15095,15);


 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,150,dvrr_stack+17834, dvrr_stack+16934, dvrr_stack+14435);

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+707, dvrr_stack+3561, dvrr_stack+10767, NULL, NULL, dvrr_stack+3540);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+18284, dvrr_stack+8255, dvrr_stack+707, dvrr_stack+1901, dvrr_stack+3561, dvrr_stack+15545);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+18452, dvrr_stack+8339, dvrr_stack+18284, dvrr_stack+1929, dvrr_stack+8255, dvrr_stack+15608);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+18732, dvrr_stack+8507, dvrr_stack+18452, dvrr_stack+2097, dvrr_stack+8339, dvrr_stack+15734);

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+19152,dvrr_stack+18732,dvrr_stack+15944,15);


 /* compute (4 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+20097,dvrr_stack+19152,dvrr_stack+16259,15);


 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,225,dvrr_stack+21447, dvrr_stack+20097, dvrr_stack+14870);

 /* compute (0 0 | 1 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1901, Data->F+10, Data->F+11, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+1790, dvrr_stack+1901, Data->F+9, Data->F+10, NULL);

 /* compute (0 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1901, dvrr_stack+1793, dvrr_stack+6, dvrr_stack+1818, dvrr_stack+1790, NULL);

 /* compute (0 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1911, dvrr_stack+39, dvrr_stack+1901, dvrr_stack+1821, dvrr_stack+1793, NULL);

 /* compute (0 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3820, dvrr_stack+207, dvrr_stack+1911, dvrr_stack+1827, dvrr_stack+39, NULL);

 /* compute (0 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1901, dvrr_stack+631, dvrr_stack+3820, dvrr_stack+3525, dvrr_stack+207, NULL);

 /* compute (0 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1929, dvrr_stack+10767, dvrr_stack+1901, dvrr_stack+3540, dvrr_stack+631, NULL);

 /* compute (1 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+22122, dvrr_stack+10795, dvrr_stack+1929, NULL, NULL, dvrr_stack+10767);

 /* compute (2 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+22230, dvrr_stack+10831, dvrr_stack+22122, dvrr_stack+3589, dvrr_stack+10795, dvrr_stack+707);

 /* compute (3 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+22446, dvrr_stack+10939, dvrr_stack+22230, dvrr_stack+3625, dvrr_stack+10831, dvrr_stack+18284);

 /* compute (4 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+22806, dvrr_stack+11155, dvrr_stack+22446, dvrr_stack+3841, dvrr_stack+10939, dvrr_stack+18452);

 /* compute (4 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+23346,dvrr_stack+22806,dvrr_stack+18732,15);


 /* compute (4 0 | 5 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hd(Libderiv->CD,dvrr_stack+24606,dvrr_stack+23346,dvrr_stack+19152,15);


 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,315,dvrr_stack+26496, dvrr_stack+24606, dvrr_stack+15944);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,60,dvrr_stack+10795, dvrr_stack+1250, dvrr_stack+147);

 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,90,dvrr_stack+22122, dvrr_stack+2643, dvrr_stack+361);

 /* compute (2 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,126,dvrr_stack+22392, dvrr_stack+4561, dvrr_stack+854);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,100,dvrr_stack+27441, dvrr_stack+7355, dvrr_stack+5821);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+27741, dvrr_stack+9417, dvrr_stack+6056);

 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,210,dvrr_stack+28191, dvrr_stack+12355, dvrr_stack+6695);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,150,dvrr_stack+28821, dvrr_stack+16934, dvrr_stack+14435);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,225,dvrr_stack+29271, dvrr_stack+20097, dvrr_stack+14870);

 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,315,dvrr_stack+29946, dvrr_stack+24606, dvrr_stack+15944);

 /* compute (2 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,60,dvrr_stack+10975, dvrr_stack+1250, dvrr_stack+147);

 /* compute (2 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,90,dvrr_stack+1250, dvrr_stack+2643, dvrr_stack+361);

 /* compute (2 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,126,dvrr_stack+2643, dvrr_stack+4561, dvrr_stack+854);

 /* compute (3 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,100,dvrr_stack+4561, dvrr_stack+7355, dvrr_stack+5821);

 /* compute (3 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+7355, dvrr_stack+9417, dvrr_stack+6056);

 /* compute (3 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,210,dvrr_stack+9417, dvrr_stack+12355, dvrr_stack+6695);

 /* compute (4 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,150,dvrr_stack+12355, dvrr_stack+16934, dvrr_stack+14435);

 /* compute (4 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,225,dvrr_stack+16934, dvrr_stack+20097, dvrr_stack+14870);

 /* compute (4 0 | 5 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,315,dvrr_stack+20097, dvrr_stack+24606, dvrr_stack+15944);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+24606, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+2);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+24615, dvrr_stack+52, dvrr_stack+15, NULL, NULL, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+22770, dvrr_stack+24615, dvrr_stack+21, dvrr_stack+52, dvrr_stack+15, dvrr_stack+24606);

 /* compute (2 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+24633,dvrr_stack+147,dvrr_stack+22770,6);


 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,36,dvrr_stack+24741, dvrr_stack+24633, NULL);

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
 deriv_build_DZ_0(Data,126,dvrr_stack+24777, dvrr_stack+2265, NULL);
 tmp = dvrr_stack + 24777;
 target_ptr = Libderiv->deriv_classes[2][5][11];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,168,dvrr_stack+24903, dvrr_stack+4057, NULL);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+15, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+25071, dvrr_stack+24606, dvrr_stack+1837, dvrr_stack+12, dvrr_stack+0, dvrr_stack+15);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+25089, dvrr_stack+22770, dvrr_stack+5695, dvrr_stack+24615, dvrr_stack+21, dvrr_stack+25071);

 /* compute (3 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+25149,dvrr_stack+5821,dvrr_stack+25089,10);


 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,60,dvrr_stack+25329, dvrr_stack+25149, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+7805, dvrr_stack+6905, NULL);
 tmp = dvrr_stack + 7805;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,100,dvrr_stack+25389, dvrr_stack+6206, NULL);
 tmp = dvrr_stack + 25389;
 target_ptr = Libderiv->deriv_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+25489, dvrr_stack+8787, NULL);
 tmp = dvrr_stack + 25489;
 target_ptr = Libderiv->deriv_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,280,dvrr_stack+25699, dvrr_stack+11515, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+52, dvrr_stack+15, dvrr_stack+1846, Data->F+2, Data->F+3, NULL);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+18, dvrr_stack+25071, dvrr_stack+3507, dvrr_stack+24606, dvrr_stack+1837, dvrr_stack+52);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+25979, dvrr_stack+25089, dvrr_stack+3760, dvrr_stack+22770, dvrr_stack+5695, dvrr_stack+18);

 /* compute (4 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+10047,dvrr_stack+14435,dvrr_stack+25979,15);


 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+26069, dvrr_stack+10047, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+17609, dvrr_stack+16259, NULL);
 tmp = dvrr_stack + 17609;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+26159, dvrr_stack+15095, NULL);
 tmp = dvrr_stack + 26159;
 target_ptr = Libderiv->deriv_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+21042, dvrr_stack+19152, NULL);
 tmp = dvrr_stack + 21042;
 target_ptr = Libderiv->deriv_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+12805, dvrr_stack+23346, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,36,dvrr_stack+26309, dvrr_stack+24633, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+21357, dvrr_stack+980, NULL);
 tmp = dvrr_stack + 21357;
 target_ptr = Libderiv->deriv_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+26345, dvrr_stack+451, NULL);
 tmp = dvrr_stack + 26345;
 target_ptr = Libderiv->deriv_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+13225, dvrr_stack+2265, NULL);
 tmp = dvrr_stack + 13225;
 target_ptr = Libderiv->deriv_classes[2][5][10];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,168,dvrr_stack+13351, dvrr_stack+4057, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+26405, dvrr_stack+25149, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+4861, dvrr_stack+6905, NULL);
 tmp = dvrr_stack + 4861;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+5011, dvrr_stack+6206, NULL);
 tmp = dvrr_stack + 5011;
 target_ptr = Libderiv->deriv_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+30891, dvrr_stack+8787, NULL);
 tmp = dvrr_stack + 30891;
 target_ptr = Libderiv->deriv_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+31101, dvrr_stack+11515, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+13519, dvrr_stack+10047, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+31381, dvrr_stack+16259, NULL);
 tmp = dvrr_stack + 31381;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+5111, dvrr_stack+15095, NULL);
 tmp = dvrr_stack + 5111;
 target_ptr = Libderiv->deriv_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+31606, dvrr_stack+19152, NULL);
 tmp = dvrr_stack + 31606;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+31921, dvrr_stack+23346, NULL);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,36,dvrr_stack+5261, dvrr_stack+24633, NULL);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+24633, dvrr_stack+980, NULL);
 tmp = dvrr_stack + 24633;
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
 deriv_build_DX_0(Data,126,dvrr_stack+451, dvrr_stack+2265, NULL);
 tmp = dvrr_stack + 451;
 target_ptr = Libderiv->deriv_classes[2][5][9];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,168,dvrr_stack+2265, dvrr_stack+4057, NULL);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+4057, dvrr_stack+25149, NULL);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+25149, dvrr_stack+6905, NULL);
 tmp = dvrr_stack + 25149;
 target_ptr = Libderiv->deriv_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+6905, dvrr_stack+6206, NULL);
 tmp = dvrr_stack + 6905;
 target_ptr = Libderiv->deriv_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+2433, dvrr_stack+8787, NULL);
 tmp = dvrr_stack + 2433;
 target_ptr = Libderiv->deriv_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+8787, dvrr_stack+11515, NULL);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+11515, dvrr_stack+10047, NULL);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+10047, dvrr_stack+16259, NULL);
 tmp = dvrr_stack + 10047;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+16259, dvrr_stack+15095, NULL);
 tmp = dvrr_stack + 16259;
 target_ptr = Libderiv->deriv_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+15095, dvrr_stack+19152, NULL);
 tmp = dvrr_stack + 15095;
 target_ptr = Libderiv->deriv_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+19152, dvrr_stack+23346, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+23346, dvrr_stack+49, dvrr_stack+12, NULL, NULL, Data->F+1);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+24723, dvrr_stack+23346, dvrr_stack+24606, dvrr_stack+49, dvrr_stack+12, dvrr_stack+0);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,6,1,dvrr_stack+23355, dvrr_stack+147, dvrr_stack+24723);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+23391, dvrr_stack+854, dvrr_stack+147);
 tmp = dvrr_stack + 23391;
 target_ptr = Libderiv->deriv_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+23481, dvrr_stack+361, dvrr_stack+22770);
 tmp = dvrr_stack + 23481;
 target_ptr = Libderiv->deriv_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+23541, dvrr_stack+2097, dvrr_stack+361);
 tmp = dvrr_stack + 23541;
 target_ptr = Libderiv->deriv_classes[2][5][8];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+23667, dvrr_stack+3841, dvrr_stack+854);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+13609, dvrr_stack+0, dvrr_stack+15, Data->F+1, Data->F+2, NULL);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+25299, dvrr_stack+24723, dvrr_stack+25071, dvrr_stack+23346, dvrr_stack+24606, dvrr_stack+13609);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,10,1,dvrr_stack+23835, dvrr_stack+5821, dvrr_stack+25299);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+23895, dvrr_stack+6695, dvrr_stack+5821);
 tmp = dvrr_stack + 23895;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+24045, dvrr_stack+6056, dvrr_stack+25089);
 tmp = dvrr_stack + 24045;
 target_ptr = Libderiv->deriv_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+1040, dvrr_stack+8507, dvrr_stack+6056);
 tmp = dvrr_stack + 1040;
 target_ptr = Libderiv->deriv_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+24145, dvrr_stack+11155, dvrr_stack+6695);

 /* compute (3 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+24425, dvrr_stack+13609, dvrr_stack+52, dvrr_stack+0, dvrr_stack+15, NULL);

 /* compute (4 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+10272, dvrr_stack+25299, dvrr_stack+18, dvrr_stack+24723, dvrr_stack+25071, dvrr_stack+24425);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,15,1,dvrr_stack+24425, dvrr_stack+14435, dvrr_stack+10272);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+19572, dvrr_stack+15944, dvrr_stack+14435);
 tmp = dvrr_stack + 19572;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+19797, dvrr_stack+14870, dvrr_stack+25979);
 tmp = dvrr_stack + 19797;
 target_ptr = Libderiv->deriv_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+16409, dvrr_stack+18732, dvrr_stack+14870);
 tmp = dvrr_stack + 16409;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+11605, dvrr_stack+22806, dvrr_stack+15944);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,6,1,dvrr_stack+24515, dvrr_stack+147, dvrr_stack+24723);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+19947, dvrr_stack+854, dvrr_stack+147);
 tmp = dvrr_stack + 19947;
 target_ptr = Libderiv->deriv_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+20037, dvrr_stack+361, dvrr_stack+22770);
 tmp = dvrr_stack + 20037;
 target_ptr = Libderiv->deriv_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+15410, dvrr_stack+2097, dvrr_stack+361);
 tmp = dvrr_stack + 15410;
 target_ptr = Libderiv->deriv_classes[2][5][7];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+16724, dvrr_stack+3841, dvrr_stack+854);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,10,1,dvrr_stack+24551, dvrr_stack+5821, dvrr_stack+25299);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+12025, dvrr_stack+6695, dvrr_stack+5821);
 tmp = dvrr_stack + 12025;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+12175, dvrr_stack+6056, dvrr_stack+25089);
 tmp = dvrr_stack + 12175;
 target_ptr = Libderiv->deriv_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+9067, dvrr_stack+8507, dvrr_stack+6056);
 tmp = dvrr_stack + 9067;
 target_ptr = Libderiv->deriv_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+6206, dvrr_stack+11155, dvrr_stack+6695);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,15,1,dvrr_stack+9277, dvrr_stack+14435, dvrr_stack+10272);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+7005, dvrr_stack+15944, dvrr_stack+14435);
 tmp = dvrr_stack + 7005;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+4117, dvrr_stack+14870, dvrr_stack+25979);
 tmp = dvrr_stack + 4117;
 target_ptr = Libderiv->deriv_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+32341, dvrr_stack+18732, dvrr_stack+14870);
 tmp = dvrr_stack + 32341;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+32656, dvrr_stack+22806, dvrr_stack+15944);

 /* compute (2 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,6,1,dvrr_stack+16892, dvrr_stack+147, dvrr_stack+24723);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+7230, dvrr_stack+854, dvrr_stack+147);
 tmp = dvrr_stack + 7230;
 target_ptr = Libderiv->deriv_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+12275, dvrr_stack+361, dvrr_stack+22770);
 tmp = dvrr_stack + 12275;
 target_ptr = Libderiv->deriv_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+4267, dvrr_stack+2097, dvrr_stack+361);
 tmp = dvrr_stack + 4267;
 target_ptr = Libderiv->deriv_classes[2][5][6];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+2097, dvrr_stack+3841, dvrr_stack+854);

 /* compute (3 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,10,1,dvrr_stack+3841, dvrr_stack+5821, dvrr_stack+25299);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+3901, dvrr_stack+6695, dvrr_stack+5821);
 tmp = dvrr_stack + 3901;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+4393, dvrr_stack+6056, dvrr_stack+25089);
 tmp = dvrr_stack + 4393;
 target_ptr = Libderiv->deriv_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+33076, dvrr_stack+8507, dvrr_stack+6056);
 tmp = dvrr_stack + 33076;
 target_ptr = Libderiv->deriv_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+33286, dvrr_stack+11155, dvrr_stack+6695);

 /* compute (4 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,15,1,dvrr_stack+11155, dvrr_stack+14435, dvrr_stack+10272);

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+11245, dvrr_stack+15944, dvrr_stack+14435);
 tmp = dvrr_stack + 11245;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+3021, dvrr_stack+14870, dvrr_stack+25979);
 tmp = dvrr_stack + 3021;
 target_ptr = Libderiv->deriv_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+33566, dvrr_stack+18732, dvrr_stack+14870);
 tmp = dvrr_stack + 33566;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+33881, dvrr_stack+22806, dvrr_stack+15944);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+22770,dvrr_stack+316,dvrr_stack+117,3);


 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+25299, dvrr_stack+22770, NULL);

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+22860,dvrr_stack+791,dvrr_stack+316,3);


 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+10272, dvrr_stack+22860, NULL);

 /* compute (1 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+22995,dvrr_stack+2013,dvrr_stack+791,3);


 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,63,dvrr_stack+23184, dvrr_stack+22995, NULL);

 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+4051, dvrr_stack+1846, dvrr_stack+0, Data->F+3, Data->F+4, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+15536, dvrr_stack+68, dvrr_stack+237, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+24723, dvrr_stack+3733, dvrr_stack+15536, dvrr_stack+3, dvrr_stack+68, dvrr_stack+0);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+11470, dvrr_stack+3507, dvrr_stack+24723, dvrr_stack+1837, dvrr_stack+3733, dvrr_stack+4051);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+23247, dvrr_stack+240, dvrr_stack+676, NULL, NULL, dvrr_stack+237);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+23265, dvrr_stack+3742, dvrr_stack+23247, dvrr_stack+71, dvrr_stack+240, dvrr_stack+15536);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+4493, dvrr_stack+3453, dvrr_stack+23265, dvrr_stack+3489, dvrr_stack+3742, dvrr_stack+24723);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+3589, dvrr_stack+3760, dvrr_stack+4493, dvrr_stack+5695, dvrr_stack+3453, dvrr_stack+11470);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+5695, dvrr_stack+682, dvrr_stack+1855, NULL, NULL, dvrr_stack+676);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1929, dvrr_stack+14245, dvrr_stack+5695, dvrr_stack+246, dvrr_stack+682, dvrr_stack+23247);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+34301, dvrr_stack+14275, dvrr_stack+1929, dvrr_stack+5731, dvrr_stack+14245, dvrr_stack+23265);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+34401, dvrr_stack+14335, dvrr_stack+34301, dvrr_stack+5761, dvrr_stack+14275, dvrr_stack+4493);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+34551, dvrr_stack+14435, dvrr_stack+34401, dvrr_stack+5821, dvrr_stack+14335, dvrr_stack+3589);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5725, dvrr_stack+1865, dvrr_stack+3525, NULL, NULL, dvrr_stack+1855);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+34761, dvrr_stack+14585, dvrr_stack+5725, dvrr_stack+692, dvrr_stack+1865, dvrr_stack+5695);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+34851, dvrr_stack+14630, dvrr_stack+34761, dvrr_stack+5921, dvrr_stack+14585, dvrr_stack+1929);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+35001, dvrr_stack+14720, dvrr_stack+34851, dvrr_stack+5966, dvrr_stack+14630, dvrr_stack+34301);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+35226, dvrr_stack+14870, dvrr_stack+35001, dvrr_stack+6056, dvrr_stack+14720, dvrr_stack+34401);

 /* compute (5 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+35541,dvrr_stack+35226,dvrr_stack+34551,21);


 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+36171, dvrr_stack+35541, NULL);

 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+5921, dvrr_stack+3540, dvrr_stack+631, NULL, NULL, dvrr_stack+3525);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+36381, dvrr_stack+15545, dvrr_stack+5921, dvrr_stack+1880, dvrr_stack+3540, dvrr_stack+5725);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+36507, dvrr_stack+15608, dvrr_stack+36381, dvrr_stack+6506, dvrr_stack+15545, dvrr_stack+34761);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+36717, dvrr_stack+15734, dvrr_stack+36507, dvrr_stack+6569, dvrr_stack+15608, dvrr_stack+34851);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+37032, dvrr_stack+15944, dvrr_stack+36717, dvrr_stack+6695, dvrr_stack+15734, dvrr_stack+35001);

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+37473,dvrr_stack+37032,dvrr_stack+35226,21);


 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+38418, dvrr_stack+37473, NULL);

 /* compute (1 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6486, dvrr_stack+10767, dvrr_stack+1901, NULL, NULL, dvrr_stack+631);

 /* compute (2 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+38733, dvrr_stack+707, dvrr_stack+6486, dvrr_stack+3561, dvrr_stack+10767, dvrr_stack+5921);

 /* compute (3 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+38901, dvrr_stack+18284, dvrr_stack+38733, dvrr_stack+8255, dvrr_stack+707, dvrr_stack+36381);

 /* compute (4 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+39181, dvrr_stack+18452, dvrr_stack+38901, dvrr_stack+8339, dvrr_stack+18284, dvrr_stack+36507);

 /* compute (5 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+39601, dvrr_stack+18732, dvrr_stack+39181, dvrr_stack+8507, dvrr_stack+18452, dvrr_stack+36717);

 /* compute (5 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+40189,dvrr_stack+39601,dvrr_stack+37032,21);


 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,441,dvrr_stack+18284, dvrr_stack+40189, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+18725, dvrr_stack+22770, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+18755, dvrr_stack+22860, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,63,dvrr_stack+18800, dvrr_stack+22995, NULL);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+18863, dvrr_stack+35541, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+8255, dvrr_stack+37473, NULL);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+38733, dvrr_stack+40189, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+19073, dvrr_stack+22770, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+22770, dvrr_stack+22860, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,63,dvrr_stack+22815, dvrr_stack+22995, NULL);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+22878, dvrr_stack+35541, NULL);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+35541, dvrr_stack+37473, NULL);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+37473, dvrr_stack+40189, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+40189, dvrr_stack+316, dvrr_stack+24615);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+40219, dvrr_stack+791, dvrr_stack+117);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,3,1,dvrr_stack+40264, dvrr_stack+2013, dvrr_stack+316);

 /* compute (3 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+246, dvrr_stack+52, dvrr_stack+4051, dvrr_stack+15, dvrr_stack+1846, NULL);

 /* compute (4 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+40327, dvrr_stack+18, dvrr_stack+11470, dvrr_stack+25071, dvrr_stack+3507, dvrr_stack+246);

 /* compute (5 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+40372, dvrr_stack+25979, dvrr_stack+3589, dvrr_stack+25089, dvrr_stack+3760, dvrr_stack+40327);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+40498, dvrr_stack+35226, dvrr_stack+40372);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+35856, dvrr_stack+37032, dvrr_stack+34551);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+40708, dvrr_stack+39601, dvrr_stack+35226);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+40327, dvrr_stack+316, dvrr_stack+24615);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+25979, dvrr_stack+791, dvrr_stack+117);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,3,1,dvrr_stack+25071, dvrr_stack+2013, dvrr_stack+316);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+41149, dvrr_stack+35226, dvrr_stack+40372);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+41359, dvrr_stack+37032, dvrr_stack+34551);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+41674, dvrr_stack+39601, dvrr_stack+35226);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+26024, dvrr_stack+316, dvrr_stack+24615);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+42115, dvrr_stack+791, dvrr_stack+117);

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,3,1,dvrr_stack+42160, dvrr_stack+2013, dvrr_stack+316);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+42223, dvrr_stack+35226, dvrr_stack+40372);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+37914, dvrr_stack+37032, dvrr_stack+34551);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+42433, dvrr_stack+39601, dvrr_stack+35226);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+40357, dvrr_stack+147, dvrr_stack+58);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+40387, dvrr_stack+14435, dvrr_stack+147);
 tmp = dvrr_stack + 40387;
 target_ptr = Libderiv->deriv_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+42874, dvrr_stack+361, dvrr_stack+222);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+42919, dvrr_stack+14870, dvrr_stack+361);
 tmp = dvrr_stack + 42919;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,21,dvrr_stack+43069, dvrr_stack+854, dvrr_stack+652);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+43132, dvrr_stack+15944, dvrr_stack+854);
 tmp = dvrr_stack + 43132;
 target_ptr = Libderiv->deriv_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+43342, dvrr_stack+5821, dvrr_stack+117);
 tmp = dvrr_stack + 43342;
 target_ptr = Libderiv->deriv_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+38229, dvrr_stack+34551, dvrr_stack+5821);
 tmp = dvrr_stack + 38229;
 target_ptr = Libderiv->deriv_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+43402, dvrr_stack+6056, dvrr_stack+316);
 tmp = dvrr_stack + 43402;
 target_ptr = Libderiv->deriv_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+39174, dvrr_stack+35226, dvrr_stack+6056);
 tmp = dvrr_stack + 39174;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+8570, dvrr_stack+6695, dvrr_stack+791);
 tmp = dvrr_stack + 8570;
 target_ptr = Libderiv->deriv_classes[2][5][2];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+39399, dvrr_stack+37032, dvrr_stack+6695);
 tmp = dvrr_stack + 39399;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+40487, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+13609, dvrr_stack+0, dvrr_stack+40487, Data->F+4, Data->F+5, NULL);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+246, dvrr_stack+4051, dvrr_stack+13609, dvrr_stack+1846, dvrr_stack+0, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+43492, dvrr_stack+237, dvrr_stack+673, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+3489, dvrr_stack+15536, dvrr_stack+43492, dvrr_stack+68, dvrr_stack+237, dvrr_stack+40487);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+68, dvrr_stack+24723, dvrr_stack+3489, dvrr_stack+3733, dvrr_stack+15536, dvrr_stack+13609);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+11470, dvrr_stack+68, dvrr_stack+3507, dvrr_stack+24723, dvrr_stack+246);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+24723, dvrr_stack+676, dvrr_stack+1849, NULL, NULL, dvrr_stack+673);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+11470, dvrr_stack+23247, dvrr_stack+24723, dvrr_stack+240, dvrr_stack+676, dvrr_stack+43492);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+43492, dvrr_stack+23265, dvrr_stack+11470, dvrr_stack+3742, dvrr_stack+23247, dvrr_stack+3489);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+43552, dvrr_stack+4493, dvrr_stack+43492, dvrr_stack+3453, dvrr_stack+23265, dvrr_stack+68);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+43642, dvrr_stack+3589, dvrr_stack+43552, dvrr_stack+3760, dvrr_stack+4493, dvrr_stack+0);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+0, dvrr_stack+1855, dvrr_stack+1827, NULL, NULL, dvrr_stack+1849);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4493, dvrr_stack+5695, dvrr_stack+0, dvrr_stack+682, dvrr_stack+1855, dvrr_stack+24723);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+23247, dvrr_stack+1929, dvrr_stack+4493, dvrr_stack+14245, dvrr_stack+5695, dvrr_stack+11470);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+43768, dvrr_stack+34301, dvrr_stack+23247, dvrr_stack+14275, dvrr_stack+1929, dvrr_stack+43492);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+39714, dvrr_stack+34401, dvrr_stack+43768, dvrr_stack+14335, dvrr_stack+34301, dvrr_stack+43552);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+43918, dvrr_stack+34551, dvrr_stack+39714, dvrr_stack+14435, dvrr_stack+34401, dvrr_stack+43642);

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+34301, dvrr_stack+43918, dvrr_stack+14435);

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11470, dvrr_stack+3525, dvrr_stack+207, NULL, NULL, dvrr_stack+1827);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+43492, dvrr_stack+5725, dvrr_stack+11470, dvrr_stack+1865, dvrr_stack+3525, dvrr_stack+0);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+43582, dvrr_stack+34761, dvrr_stack+43492, dvrr_stack+14585, dvrr_stack+5725, dvrr_stack+4493);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+44198, dvrr_stack+34851, dvrr_stack+43582, dvrr_stack+14630, dvrr_stack+34761, dvrr_stack+23247);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+44423, dvrr_stack+35001, dvrr_stack+44198, dvrr_stack+14720, dvrr_stack+34851, dvrr_stack+43768);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+44738, dvrr_stack+35226, dvrr_stack+44423, dvrr_stack+14870, dvrr_stack+35001, dvrr_stack+39714);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+39714, dvrr_stack+44738, dvrr_stack+14870);

 /* compute (1 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+40029, dvrr_stack+631, dvrr_stack+3820, NULL, NULL, dvrr_stack+207);

 /* compute (2 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+5695, dvrr_stack+5921, dvrr_stack+40029, dvrr_stack+3540, dvrr_stack+631, dvrr_stack+11470);

 /* compute (3 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+34761, dvrr_stack+36381, dvrr_stack+5695, dvrr_stack+15545, dvrr_stack+5921, dvrr_stack+43492);

 /* compute (4 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3453, dvrr_stack+36507, dvrr_stack+34761, dvrr_stack+15608, dvrr_stack+36381, dvrr_stack+43582);

 /* compute (5 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+34761, dvrr_stack+36717, dvrr_stack+3453, dvrr_stack+15734, dvrr_stack+36507, dvrr_stack+44198);

 /* compute (6 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+45158, dvrr_stack+37032, dvrr_stack+34761, dvrr_stack+15944, dvrr_stack+36717, dvrr_stack+44423);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+34761, dvrr_stack+45158, dvrr_stack+15944);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+44198, dvrr_stack+147, dvrr_stack+58);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+44228, dvrr_stack+14435, dvrr_stack+147);
 tmp = dvrr_stack + 44228;
 target_ptr = Libderiv->deriv_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+11470, dvrr_stack+361, dvrr_stack+222);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+44328, dvrr_stack+14870, dvrr_stack+361);
 tmp = dvrr_stack + 44328;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,21,dvrr_stack+44478, dvrr_stack+854, dvrr_stack+652);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+3453, dvrr_stack+15944, dvrr_stack+854);
 tmp = dvrr_stack + 3453;
 target_ptr = Libderiv->deriv_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+44541, dvrr_stack+5821, dvrr_stack+117);
 tmp = dvrr_stack + 44541;
 target_ptr = Libderiv->deriv_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+3663, dvrr_stack+34551, dvrr_stack+5821);
 tmp = dvrr_stack + 3663;
 target_ptr = Libderiv->deriv_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+44601, dvrr_stack+6056, dvrr_stack+316);
 tmp = dvrr_stack + 44601;
 target_ptr = Libderiv->deriv_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+36381, dvrr_stack+35226, dvrr_stack+6056);
 tmp = dvrr_stack + 36381;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+5695, dvrr_stack+6695, dvrr_stack+791);
 tmp = dvrr_stack + 5695;
 target_ptr = Libderiv->deriv_classes[2][5][1];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+36606, dvrr_stack+37032, dvrr_stack+6695);
 tmp = dvrr_stack + 36606;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+43492, dvrr_stack+43918, dvrr_stack+14435);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+15536, dvrr_stack+44738, dvrr_stack+14870);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+45746, dvrr_stack+45158, dvrr_stack+15944);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+44691, dvrr_stack+147, dvrr_stack+58);

 /* compute (3 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+36921, dvrr_stack+14435, dvrr_stack+147);
 tmp = dvrr_stack + 36921;
 target_ptr = Libderiv->deriv_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+147, dvrr_stack+361, dvrr_stack+222);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+43702, dvrr_stack+14870, dvrr_stack+361);
 tmp = dvrr_stack + 43702;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,21,dvrr_stack+361, dvrr_stack+854, dvrr_stack+652);

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+14585, dvrr_stack+15944, dvrr_stack+854);
 tmp = dvrr_stack + 14585;
 target_ptr = Libderiv->deriv_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+854, dvrr_stack+5821, dvrr_stack+117);
 tmp = dvrr_stack + 854;
 target_ptr = Libderiv->deriv_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+40029, dvrr_stack+34551, dvrr_stack+5821);
 tmp = dvrr_stack + 40029;
 target_ptr = Libderiv->deriv_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+5821, dvrr_stack+6056, dvrr_stack+316);
 tmp = dvrr_stack + 5821;
 target_ptr = Libderiv->deriv_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+34511, dvrr_stack+35226, dvrr_stack+6056);
 tmp = dvrr_stack + 34511;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+5911, dvrr_stack+6695, dvrr_stack+791);
 tmp = dvrr_stack + 5911;
 target_ptr = Libderiv->deriv_classes[2][5][0];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+35202, dvrr_stack+37032, dvrr_stack+6695);
 tmp = dvrr_stack + 35202;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+37021, dvrr_stack+43918, dvrr_stack+14435);

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+43852, dvrr_stack+44738, dvrr_stack+14870);

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+46187, dvrr_stack+45158, dvrr_stack+15944);

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,60,dvrr_stack+6037, dvrr_stack+1610, NULL);
 tmp = dvrr_stack + 6037;
 target_ptr = Libderiv->deriv2_classes[2][3][143];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,90,dvrr_stack+6097, dvrr_stack+3183, NULL);
 tmp = dvrr_stack + 6097;
 target_ptr = Libderiv->deriv2_classes[2][4][143];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,126,dvrr_stack+44721, dvrr_stack+5317, NULL);
 tmp = dvrr_stack + 44721;
 target_ptr = Libderiv->deriv2_classes[2][5][143];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,100,dvrr_stack+44847, dvrr_stack+7955, NULL);
 tmp = dvrr_stack + 44847;
 target_ptr = Libderiv->deriv2_classes[3][3][143];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+44947, dvrr_stack+10317, NULL);
 tmp = dvrr_stack + 44947;
 target_ptr = Libderiv->deriv2_classes[3][4][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,210,dvrr_stack+45097, dvrr_stack+13615, NULL);
 tmp = dvrr_stack + 45097;
 target_ptr = Libderiv->deriv2_classes[3][5][143];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,150,dvrr_stack+45307, dvrr_stack+17834, NULL);
 tmp = dvrr_stack + 45307;
 target_ptr = Libderiv->deriv2_classes[4][3][143];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,225,dvrr_stack+45457, dvrr_stack+21447, NULL);
 tmp = dvrr_stack + 45457;
 target_ptr = Libderiv->deriv2_classes[4][4][143];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,315,dvrr_stack+14245, dvrr_stack+26496, NULL);
 tmp = dvrr_stack + 14245;
 target_ptr = Libderiv->deriv2_classes[4][5][143];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,60,dvrr_stack+914, dvrr_stack+1610, NULL);
 tmp = dvrr_stack + 914;
 target_ptr = Libderiv->deriv2_classes[2][3][131];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,90,dvrr_stack+37231, dvrr_stack+3183, NULL);
 tmp = dvrr_stack + 37231;
 target_ptr = Libderiv->deriv2_classes[2][4][131];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,126,dvrr_stack+37321, dvrr_stack+5317, NULL);
 tmp = dvrr_stack + 37321;
 target_ptr = Libderiv->deriv2_classes[2][5][131];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,100,dvrr_stack+23247, dvrr_stack+7955, NULL);
 tmp = dvrr_stack + 23247;
 target_ptr = Libderiv->deriv2_classes[3][3][131];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+14795, dvrr_stack+10317, NULL);
 tmp = dvrr_stack + 14795;
 target_ptr = Libderiv->deriv2_classes[3][4][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,210,dvrr_stack+46628, dvrr_stack+13615, NULL);
 tmp = dvrr_stack + 46628;
 target_ptr = Libderiv->deriv2_classes[3][5][131];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,150,dvrr_stack+14945, dvrr_stack+17834, NULL);
 tmp = dvrr_stack + 14945;
 target_ptr = Libderiv->deriv2_classes[4][3][131];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,225,dvrr_stack+15851, dvrr_stack+21447, NULL);
 tmp = dvrr_stack + 15851;
 target_ptr = Libderiv->deriv2_classes[4][4][131];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,315,dvrr_stack+6486, dvrr_stack+26496, NULL);
 tmp = dvrr_stack + 6486;
 target_ptr = Libderiv->deriv2_classes[4][5][131];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,60,dvrr_stack+192, dvrr_stack+10795, NULL);
 tmp = dvrr_stack + 192;
 target_ptr = Libderiv->deriv2_classes[2][3][130];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+0, dvrr_stack+22122, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[2][4][130];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+46838, dvrr_stack+22392, NULL);
 tmp = dvrr_stack + 46838;
 target_ptr = Libderiv->deriv2_classes[2][5][130];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,100,dvrr_stack+16076, dvrr_stack+27441, NULL);
 tmp = dvrr_stack + 16076;
 target_ptr = Libderiv->deriv2_classes[3][3][130];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+577, dvrr_stack+27741, NULL);
 tmp = dvrr_stack + 577;
 target_ptr = Libderiv->deriv2_classes[3][4][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+1790, dvrr_stack+28191, NULL);
 tmp = dvrr_stack + 1790;
 target_ptr = Libderiv->deriv2_classes[3][5][130];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+46964, dvrr_stack+28821, NULL);
 tmp = dvrr_stack + 46964;
 target_ptr = Libderiv->deriv2_classes[4][3][130];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+47114, dvrr_stack+29271, NULL);
 tmp = dvrr_stack + 47114;
 target_ptr = Libderiv->deriv2_classes[4][4][130];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+47339, dvrr_stack+29946, NULL);
 tmp = dvrr_stack + 47339;
 target_ptr = Libderiv->deriv2_classes[4][5][130];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,60,dvrr_stack+45682, dvrr_stack+1610, NULL);
 tmp = dvrr_stack + 45682;
 target_ptr = Libderiv->deriv2_classes[2][3][119];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,90,dvrr_stack+1610, dvrr_stack+3183, NULL);
 tmp = dvrr_stack + 1610;
 target_ptr = Libderiv->deriv2_classes[2][4][119];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,126,dvrr_stack+47654, dvrr_stack+5317, NULL);
 tmp = dvrr_stack + 47654;
 target_ptr = Libderiv->deriv2_classes[2][5][119];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,100,dvrr_stack+47780, dvrr_stack+7955, NULL);
 tmp = dvrr_stack + 47780;
 target_ptr = Libderiv->deriv2_classes[3][3][119];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+7955, dvrr_stack+10317, NULL);
 tmp = dvrr_stack + 7955;
 target_ptr = Libderiv->deriv2_classes[3][4][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,210,dvrr_stack+10317, dvrr_stack+13615, NULL);
 tmp = dvrr_stack + 10317;
 target_ptr = Libderiv->deriv2_classes[3][5][119];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,150,dvrr_stack+8105, dvrr_stack+17834, NULL);
 tmp = dvrr_stack + 8105;
 target_ptr = Libderiv->deriv2_classes[4][3][119];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,225,dvrr_stack+17834, dvrr_stack+21447, NULL);
 tmp = dvrr_stack + 17834;
 target_ptr = Libderiv->deriv2_classes[4][4][119];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,315,dvrr_stack+21447, dvrr_stack+26496, NULL);
 tmp = dvrr_stack + 21447;
 target_ptr = Libderiv->deriv2_classes[4][5][119];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+21762, dvrr_stack+10795, NULL);
 tmp = dvrr_stack + 21762;
 target_ptr = Libderiv->deriv2_classes[2][3][118];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+1700, dvrr_stack+22122, NULL);
 tmp = dvrr_stack + 1700;
 target_ptr = Libderiv->deriv2_classes[2][4][118];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+21822, dvrr_stack+22392, NULL);
 tmp = dvrr_stack + 21822;
 target_ptr = Libderiv->deriv2_classes[2][5][118];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+21948, dvrr_stack+27441, NULL);
 tmp = dvrr_stack + 21948;
 target_ptr = Libderiv->deriv2_classes[3][3][118];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+22048, dvrr_stack+27741, NULL);
 tmp = dvrr_stack + 22048;
 target_ptr = Libderiv->deriv2_classes[3][4][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+22198, dvrr_stack+28191, NULL);
 tmp = dvrr_stack + 22198;
 target_ptr = Libderiv->deriv2_classes[3][5][118];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+22408, dvrr_stack+28821, NULL);
 tmp = dvrr_stack + 22408;
 target_ptr = Libderiv->deriv2_classes[4][3][118];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+18059, dvrr_stack+29271, NULL);
 tmp = dvrr_stack + 18059;
 target_ptr = Libderiv->deriv2_classes[4][4][118];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+10527, dvrr_stack+29946, NULL);
 tmp = dvrr_stack + 10527;
 target_ptr = Libderiv->deriv2_classes[4][5][118];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,60,dvrr_stack+22558, dvrr_stack+10975, NULL);
 tmp = dvrr_stack + 22558;
 target_ptr = Libderiv->deriv2_classes[2][3][117];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+22618, dvrr_stack+1250, NULL);
 tmp = dvrr_stack + 22618;
 target_ptr = Libderiv->deriv2_classes[2][4][117];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+1250, dvrr_stack+2643, NULL);
 tmp = dvrr_stack + 1250;
 target_ptr = Libderiv->deriv2_classes[2][5][117];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,100,dvrr_stack+2643, dvrr_stack+4561, NULL);
 tmp = dvrr_stack + 2643;
 target_ptr = Libderiv->deriv2_classes[3][3][117];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+2743, dvrr_stack+7355, NULL);
 tmp = dvrr_stack + 2743;
 target_ptr = Libderiv->deriv2_classes[3][4][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+10842, dvrr_stack+9417, NULL);
 tmp = dvrr_stack + 10842;
 target_ptr = Libderiv->deriv2_classes[3][5][117];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+4493, dvrr_stack+12355, NULL);
 tmp = dvrr_stack + 4493;
 target_ptr = Libderiv->deriv2_classes[4][3][117];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+13609, dvrr_stack+16934, NULL);
 tmp = dvrr_stack + 13609;
 target_ptr = Libderiv->deriv2_classes[4][4][117];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+13834, dvrr_stack+20097, NULL);
 tmp = dvrr_stack + 13834;
 target_ptr = Libderiv->deriv2_classes[4][5][117];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+20097, dvrr_stack+1520, dvrr_stack+24741);
 tmp = dvrr_stack + 20097;
 target_ptr = Libderiv->deriv2_classes[2][3][107];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+20157, dvrr_stack+24777, dvrr_stack+256);
 tmp = dvrr_stack + 20157;
 target_ptr = Libderiv->deriv2_classes[2][4][107];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+20247, dvrr_stack+24903, dvrr_stack+1520);
 tmp = dvrr_stack + 20247;
 target_ptr = Libderiv->deriv2_classes[2][5][107];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+20373, dvrr_stack+7805, dvrr_stack+25329);
 tmp = dvrr_stack + 20373;
 target_ptr = Libderiv->deriv2_classes[3][3][107];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+20473, dvrr_stack+25489, dvrr_stack+25389);
 tmp = dvrr_stack + 20473;
 target_ptr = Libderiv->deriv2_classes[3][4][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+20623, dvrr_stack+25699, dvrr_stack+7805);
 tmp = dvrr_stack + 20623;
 target_ptr = Libderiv->deriv2_classes[3][5][107];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+20833, dvrr_stack+17609, dvrr_stack+26069);
 tmp = dvrr_stack + 20833;
 target_ptr = Libderiv->deriv2_classes[4][3][107];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+16928, dvrr_stack+21042, dvrr_stack+26159);
 tmp = dvrr_stack + 16928;
 target_ptr = Libderiv->deriv2_classes[4][4][107];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+17153, dvrr_stack+12805, dvrr_stack+17609);
 tmp = dvrr_stack + 17153;
 target_ptr = Libderiv->deriv2_classes[4][5][107];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+2893, dvrr_stack+21357, dvrr_stack+26309);
 tmp = dvrr_stack + 2893;
 target_ptr = Libderiv->deriv2_classes[2][3][106];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+1376, dvrr_stack+13225, dvrr_stack+26345);
 tmp = dvrr_stack + 1376;
 target_ptr = Libderiv->deriv2_classes[2][4][106];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+4643, dvrr_stack+13351, dvrr_stack+21357);
 tmp = dvrr_stack + 4643;
 target_ptr = Libderiv->deriv2_classes[2][5][106];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+11052, dvrr_stack+4861, dvrr_stack+26405);
 tmp = dvrr_stack + 11052;
 target_ptr = Libderiv->deriv2_classes[3][3][106];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+12335, dvrr_stack+30891, dvrr_stack+5011);
 tmp = dvrr_stack + 12335;
 target_ptr = Libderiv->deriv2_classes[3][4][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+12485, dvrr_stack+31101, dvrr_stack+4861);
 tmp = dvrr_stack + 12485;
 target_ptr = Libderiv->deriv2_classes[3][5][106];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+9367, dvrr_stack+31381, dvrr_stack+13519);
 tmp = dvrr_stack + 9367;
 target_ptr = Libderiv->deriv2_classes[4][3][106];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+9517, dvrr_stack+31606, dvrr_stack+5111);
 tmp = dvrr_stack + 9517;
 target_ptr = Libderiv->deriv2_classes[4][4][106];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+7320, dvrr_stack+31921, dvrr_stack+31381);
 tmp = dvrr_stack + 7320;
 target_ptr = Libderiv->deriv2_classes[4][5][106];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+2953, dvrr_stack+24633, dvrr_stack+5261);
 tmp = dvrr_stack + 2953;
 target_ptr = Libderiv->deriv2_classes[2][3][105];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+4769, dvrr_stack+451, dvrr_stack+980);
 tmp = dvrr_stack + 4769;
 target_ptr = Libderiv->deriv2_classes[2][4][105];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+17468, dvrr_stack+2265, dvrr_stack+24633);
 tmp = dvrr_stack + 17468;
 target_ptr = Libderiv->deriv2_classes[2][5][105];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+47880, dvrr_stack+25149, dvrr_stack+4057);
 tmp = dvrr_stack + 47880;
 target_ptr = Libderiv->deriv2_classes[3][3][105];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+9742, dvrr_stack+2433, dvrr_stack+6905);
 tmp = dvrr_stack + 9742;
 target_ptr = Libderiv->deriv2_classes[3][4][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+26465, dvrr_stack+8787, dvrr_stack+25149);
 tmp = dvrr_stack + 26465;
 target_ptr = Libderiv->deriv2_classes[3][5][105];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+9892, dvrr_stack+10047, dvrr_stack+11515);
 tmp = dvrr_stack + 9892;
 target_ptr = Libderiv->deriv2_classes[4][3][105];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+26675, dvrr_stack+15095, dvrr_stack+16259);
 tmp = dvrr_stack + 26675;
 target_ptr = Libderiv->deriv2_classes[4][4][105];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+26900, dvrr_stack+19152, dvrr_stack+10047);
 tmp = dvrr_stack + 26900;
 target_ptr = Libderiv->deriv2_classes[4][5][105];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,6,1,dvrr_stack+22708, dvrr_stack+23391, dvrr_stack+23355);
 tmp = dvrr_stack + 22708;
 target_ptr = Libderiv->deriv2_classes[2][3][104];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+14149, dvrr_stack+23541, dvrr_stack+23481);
 tmp = dvrr_stack + 14149;
 target_ptr = Libderiv->deriv2_classes[2][4][104];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+7635, dvrr_stack+23667, dvrr_stack+23391);
 tmp = dvrr_stack + 7635;
 target_ptr = Libderiv->deriv2_classes[2][5][104];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,10,1,dvrr_stack+6801, dvrr_stack+23895, dvrr_stack+23835);
 tmp = dvrr_stack + 6801;
 target_ptr = Libderiv->deriv2_classes[3][3][104];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+27215, dvrr_stack+1040, dvrr_stack+24045);
 tmp = dvrr_stack + 27215;
 target_ptr = Libderiv->deriv2_classes[3][4][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+27365, dvrr_stack+24145, dvrr_stack+23895);
 tmp = dvrr_stack + 27365;
 target_ptr = Libderiv->deriv2_classes[3][5][104];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+27575, dvrr_stack+19572, dvrr_stack+24425);
 tmp = dvrr_stack + 27575;
 target_ptr = Libderiv->deriv2_classes[4][3][104];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+27725, dvrr_stack+16409, dvrr_stack+19797);
 tmp = dvrr_stack + 27725;
 target_ptr = Libderiv->deriv2_classes[4][4][104];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+27950, dvrr_stack+11605, dvrr_stack+19572);
 tmp = dvrr_stack + 27950;
 target_ptr = Libderiv->deriv2_classes[4][5][104];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+16176, dvrr_stack+1520, dvrr_stack+24741);
 tmp = dvrr_stack + 16176;
 target_ptr = Libderiv->deriv2_classes[2][3][95];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+23088, dvrr_stack+24777, dvrr_stack+256);
 tmp = dvrr_stack + 23088;
 target_ptr = Libderiv->deriv2_classes[2][4][95];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+727, dvrr_stack+24903, dvrr_stack+1520);
 tmp = dvrr_stack + 727;
 target_ptr = Libderiv->deriv2_classes[2][5][95];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+12695, dvrr_stack+7805, dvrr_stack+25329);
 tmp = dvrr_stack + 12695;
 target_ptr = Libderiv->deriv2_classes[3][3][95];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+28265, dvrr_stack+25489, dvrr_stack+25389);
 tmp = dvrr_stack + 28265;
 target_ptr = Libderiv->deriv2_classes[3][4][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+28415, dvrr_stack+25699, dvrr_stack+7805);
 tmp = dvrr_stack + 28415;
 target_ptr = Libderiv->deriv2_classes[3][5][95];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+28625, dvrr_stack+17609, dvrr_stack+26069);
 tmp = dvrr_stack + 28625;
 target_ptr = Libderiv->deriv2_classes[4][3][95];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+28775, dvrr_stack+21042, dvrr_stack+26159);
 tmp = dvrr_stack + 28775;
 target_ptr = Libderiv->deriv2_classes[4][4][95];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+29000, dvrr_stack+12805, dvrr_stack+17609);
 tmp = dvrr_stack + 29000;
 target_ptr = Libderiv->deriv2_classes[4][5][95];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+8696, dvrr_stack+21357, dvrr_stack+26309);
 tmp = dvrr_stack + 8696;
 target_ptr = Libderiv->deriv2_classes[2][3][94];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+29315, dvrr_stack+13225, dvrr_stack+26345);
 tmp = dvrr_stack + 29315;
 target_ptr = Libderiv->deriv2_classes[2][4][94];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+29405, dvrr_stack+13351, dvrr_stack+21357);
 tmp = dvrr_stack + 29405;
 target_ptr = Libderiv->deriv2_classes[2][5][94];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+29531, dvrr_stack+4861, dvrr_stack+26405);
 tmp = dvrr_stack + 29531;
 target_ptr = Libderiv->deriv2_classes[3][3][94];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+29631, dvrr_stack+30891, dvrr_stack+5011);
 tmp = dvrr_stack + 29631;
 target_ptr = Libderiv->deriv2_classes[3][4][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+29781, dvrr_stack+31101, dvrr_stack+4861);
 tmp = dvrr_stack + 29781;
 target_ptr = Libderiv->deriv2_classes[3][5][94];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+29991, dvrr_stack+31381, dvrr_stack+13519);
 tmp = dvrr_stack + 29991;
 target_ptr = Libderiv->deriv2_classes[4][3][94];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+30141, dvrr_stack+31606, dvrr_stack+5111);
 tmp = dvrr_stack + 30141;
 target_ptr = Libderiv->deriv2_classes[4][4][94];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+30366, dvrr_stack+31921, dvrr_stack+31381);
 tmp = dvrr_stack + 30366;
 target_ptr = Libderiv->deriv2_classes[4][5][94];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+30681, dvrr_stack+24633, dvrr_stack+5261);
 tmp = dvrr_stack + 30681;
 target_ptr = Libderiv->deriv2_classes[2][3][93];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+30741, dvrr_stack+451, dvrr_stack+980);
 tmp = dvrr_stack + 30741;
 target_ptr = Libderiv->deriv2_classes[2][4][93];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+5297, dvrr_stack+2265, dvrr_stack+24633);
 tmp = dvrr_stack + 5297;
 target_ptr = Libderiv->deriv2_classes[2][5][93];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+5423, dvrr_stack+25149, dvrr_stack+4057);
 tmp = dvrr_stack + 5423;
 target_ptr = Libderiv->deriv2_classes[3][3][93];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+5523, dvrr_stack+2433, dvrr_stack+6905);
 tmp = dvrr_stack + 5523;
 target_ptr = Libderiv->deriv2_classes[3][4][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+3171, dvrr_stack+8787, dvrr_stack+25149);
 tmp = dvrr_stack + 3171;
 target_ptr = Libderiv->deriv2_classes[3][5][93];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+47980, dvrr_stack+10047, dvrr_stack+11515);
 tmp = dvrr_stack + 47980;
 target_ptr = Libderiv->deriv2_classes[4][3][93];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+48130, dvrr_stack+15095, dvrr_stack+16259);
 tmp = dvrr_stack + 48130;
 target_ptr = Libderiv->deriv2_classes[4][4][93];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+48355, dvrr_stack+19152, dvrr_stack+10047);
 tmp = dvrr_stack + 48355;
 target_ptr = Libderiv->deriv2_classes[4][5][93];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+30831, dvrr_stack+23391, dvrr_stack+23355);
 tmp = dvrr_stack + 30831;
 target_ptr = Libderiv->deriv2_classes[2][3][92];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+48670, dvrr_stack+23541, dvrr_stack+23481);
 tmp = dvrr_stack + 48670;
 target_ptr = Libderiv->deriv2_classes[2][4][92];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+48760, dvrr_stack+23667, dvrr_stack+23391);
 tmp = dvrr_stack + 48760;
 target_ptr = Libderiv->deriv2_classes[2][5][92];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+48886, dvrr_stack+23895, dvrr_stack+23835);
 tmp = dvrr_stack + 48886;
 target_ptr = Libderiv->deriv2_classes[3][3][92];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+48986, dvrr_stack+1040, dvrr_stack+24045);
 tmp = dvrr_stack + 48986;
 target_ptr = Libderiv->deriv2_classes[3][4][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+49136, dvrr_stack+24145, dvrr_stack+23895);
 tmp = dvrr_stack + 49136;
 target_ptr = Libderiv->deriv2_classes[3][5][92];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+49346, dvrr_stack+19572, dvrr_stack+24425);
 tmp = dvrr_stack + 49346;
 target_ptr = Libderiv->deriv2_classes[4][3][92];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+49496, dvrr_stack+16409, dvrr_stack+19797);
 tmp = dvrr_stack + 49496;
 target_ptr = Libderiv->deriv2_classes[4][4][92];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+49721, dvrr_stack+11605, dvrr_stack+19572);
 tmp = dvrr_stack + 49721;
 target_ptr = Libderiv->deriv2_classes[4][5][92];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,6,1,dvrr_stack+50036, dvrr_stack+19947, dvrr_stack+24515);
 tmp = dvrr_stack + 50036;
 target_ptr = Libderiv->deriv2_classes[2][3][91];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+50096, dvrr_stack+15410, dvrr_stack+20037);
 tmp = dvrr_stack + 50096;
 target_ptr = Libderiv->deriv2_classes[2][4][91];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+50186, dvrr_stack+16724, dvrr_stack+19947);
 tmp = dvrr_stack + 50186;
 target_ptr = Libderiv->deriv2_classes[2][5][91];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,10,1,dvrr_stack+50312, dvrr_stack+12025, dvrr_stack+24551);
 tmp = dvrr_stack + 50312;
 target_ptr = Libderiv->deriv2_classes[3][3][91];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+50412, dvrr_stack+9067, dvrr_stack+12175);
 tmp = dvrr_stack + 50412;
 target_ptr = Libderiv->deriv2_classes[3][4][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+50562, dvrr_stack+6206, dvrr_stack+12025);
 tmp = dvrr_stack + 50562;
 target_ptr = Libderiv->deriv2_classes[3][5][91];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+50772, dvrr_stack+7005, dvrr_stack+9277);
 tmp = dvrr_stack + 50772;
 target_ptr = Libderiv->deriv2_classes[4][3][91];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+50922, dvrr_stack+32341, dvrr_stack+4117);
 tmp = dvrr_stack + 50922;
 target_ptr = Libderiv->deriv2_classes[4][4][91];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+51147, dvrr_stack+32656, dvrr_stack+7005);
 tmp = dvrr_stack + 51147;
 target_ptr = Libderiv->deriv2_classes[4][5][91];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+51462, dvrr_stack+1520, dvrr_stack+24741);
 tmp = dvrr_stack + 51462;
 target_ptr = Libderiv->deriv2_classes[2][3][83];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+51522, dvrr_stack+24777, dvrr_stack+256);
 tmp = dvrr_stack + 51522;
 target_ptr = Libderiv->deriv2_classes[2][4][83];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+51612, dvrr_stack+24903, dvrr_stack+1520);
 tmp = dvrr_stack + 51612;
 target_ptr = Libderiv->deriv2_classes[2][5][83];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+24903, dvrr_stack+7805, dvrr_stack+25329);
 tmp = dvrr_stack + 24903;
 target_ptr = Libderiv->deriv2_classes[3][3][83];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+51738, dvrr_stack+25489, dvrr_stack+25389);
 tmp = dvrr_stack + 51738;
 target_ptr = Libderiv->deriv2_classes[3][4][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+51888, dvrr_stack+25699, dvrr_stack+7805);
 tmp = dvrr_stack + 51888;
 target_ptr = Libderiv->deriv2_classes[3][5][83];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+25699, dvrr_stack+17609, dvrr_stack+26069);
 tmp = dvrr_stack + 25699;
 target_ptr = Libderiv->deriv2_classes[4][3][83];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+52098, dvrr_stack+21042, dvrr_stack+26159);
 tmp = dvrr_stack + 52098;
 target_ptr = Libderiv->deriv2_classes[4][4][83];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+52323, dvrr_stack+12805, dvrr_stack+17609);
 tmp = dvrr_stack + 52323;
 target_ptr = Libderiv->deriv2_classes[4][5][83];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+25329, dvrr_stack+21357, dvrr_stack+26309);
 tmp = dvrr_stack + 25329;
 target_ptr = Libderiv->deriv2_classes[2][3][82];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+25849, dvrr_stack+13225, dvrr_stack+26345);
 tmp = dvrr_stack + 25849;
 target_ptr = Libderiv->deriv2_classes[2][4][82];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+52638, dvrr_stack+13351, dvrr_stack+21357);
 tmp = dvrr_stack + 52638;
 target_ptr = Libderiv->deriv2_classes[2][5][82];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+13351, dvrr_stack+4861, dvrr_stack+26405);
 tmp = dvrr_stack + 13351;
 target_ptr = Libderiv->deriv2_classes[3][3][82];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+52764, dvrr_stack+30891, dvrr_stack+5011);
 tmp = dvrr_stack + 52764;
 target_ptr = Libderiv->deriv2_classes[3][4][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+12795, dvrr_stack+31101, dvrr_stack+4861);
 tmp = dvrr_stack + 12795;
 target_ptr = Libderiv->deriv2_classes[3][5][82];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+31101, dvrr_stack+31381, dvrr_stack+13519);
 tmp = dvrr_stack + 31101;
 target_ptr = Libderiv->deriv2_classes[4][3][82];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+52914, dvrr_stack+31606, dvrr_stack+5111);
 tmp = dvrr_stack + 52914;
 target_ptr = Libderiv->deriv2_classes[4][4][82];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+53139, dvrr_stack+31921, dvrr_stack+31381);
 tmp = dvrr_stack + 53139;
 target_ptr = Libderiv->deriv2_classes[4][5][82];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+26405, dvrr_stack+24633, dvrr_stack+5261);
 tmp = dvrr_stack + 26405;
 target_ptr = Libderiv->deriv2_classes[2][3][81];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+31921, dvrr_stack+451, dvrr_stack+980);
 tmp = dvrr_stack + 31921;
 target_ptr = Libderiv->deriv2_classes[2][4][81];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+32011, dvrr_stack+2265, dvrr_stack+24633);
 tmp = dvrr_stack + 32011;
 target_ptr = Libderiv->deriv2_classes[2][5][81];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+2265, dvrr_stack+25149, dvrr_stack+4057);
 tmp = dvrr_stack + 2265;
 target_ptr = Libderiv->deriv2_classes[3][3][81];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+32137, dvrr_stack+2433, dvrr_stack+6905);
 tmp = dvrr_stack + 32137;
 target_ptr = Libderiv->deriv2_classes[3][4][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+53454, dvrr_stack+8787, dvrr_stack+25149);
 tmp = dvrr_stack + 53454;
 target_ptr = Libderiv->deriv2_classes[3][5][81];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+13451, dvrr_stack+10047, dvrr_stack+11515);
 tmp = dvrr_stack + 13451;
 target_ptr = Libderiv->deriv2_classes[4][3][81];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+53664, dvrr_stack+15095, dvrr_stack+16259);
 tmp = dvrr_stack + 53664;
 target_ptr = Libderiv->deriv2_classes[4][4][81];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+53889, dvrr_stack+19152, dvrr_stack+10047);
 tmp = dvrr_stack + 53889;
 target_ptr = Libderiv->deriv2_classes[4][5][81];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+11515, dvrr_stack+23391, dvrr_stack+23355);
 tmp = dvrr_stack + 11515;
 target_ptr = Libderiv->deriv2_classes[2][3][80];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+31251, dvrr_stack+23541, dvrr_stack+23481);
 tmp = dvrr_stack + 31251;
 target_ptr = Libderiv->deriv2_classes[2][4][80];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+54204, dvrr_stack+23667, dvrr_stack+23391);
 tmp = dvrr_stack + 54204;
 target_ptr = Libderiv->deriv2_classes[2][5][80];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+23667, dvrr_stack+23895, dvrr_stack+23835);
 tmp = dvrr_stack + 23667;
 target_ptr = Libderiv->deriv2_classes[3][3][80];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+54330, dvrr_stack+1040, dvrr_stack+24045);
 tmp = dvrr_stack + 54330;
 target_ptr = Libderiv->deriv2_classes[3][4][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+54480, dvrr_stack+24145, dvrr_stack+23895);
 tmp = dvrr_stack + 54480;
 target_ptr = Libderiv->deriv2_classes[3][5][80];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+24145, dvrr_stack+19572, dvrr_stack+24425);
 tmp = dvrr_stack + 24145;
 target_ptr = Libderiv->deriv2_classes[4][3][80];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+54690, dvrr_stack+16409, dvrr_stack+19797);
 tmp = dvrr_stack + 54690;
 target_ptr = Libderiv->deriv2_classes[4][4][80];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+19103, dvrr_stack+11605, dvrr_stack+19572);
 tmp = dvrr_stack + 19103;
 target_ptr = Libderiv->deriv2_classes[4][5][80];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+24295, dvrr_stack+19947, dvrr_stack+24515);
 tmp = dvrr_stack + 24295;
 target_ptr = Libderiv->deriv2_classes[2][3][79];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+24355, dvrr_stack+15410, dvrr_stack+20037);
 tmp = dvrr_stack + 24355;
 target_ptr = Libderiv->deriv2_classes[2][4][79];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+23767, dvrr_stack+16724, dvrr_stack+19947);
 tmp = dvrr_stack + 23767;
 target_ptr = Libderiv->deriv2_classes[2][5][79];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+16724, dvrr_stack+12025, dvrr_stack+24551);
 tmp = dvrr_stack + 16724;
 target_ptr = Libderiv->deriv2_classes[3][3][79];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+24445, dvrr_stack+9067, dvrr_stack+12175);
 tmp = dvrr_stack + 24445;
 target_ptr = Libderiv->deriv2_classes[3][4][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+11575, dvrr_stack+6206, dvrr_stack+12025);
 tmp = dvrr_stack + 11575;
 target_ptr = Libderiv->deriv2_classes[3][5][79];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+11785, dvrr_stack+7005, dvrr_stack+9277);
 tmp = dvrr_stack + 11785;
 target_ptr = Libderiv->deriv2_classes[4][3][79];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+6187, dvrr_stack+32341, dvrr_stack+4117);
 tmp = dvrr_stack + 6187;
 target_ptr = Libderiv->deriv2_classes[4][4][79];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+54915, dvrr_stack+32656, dvrr_stack+7005);
 tmp = dvrr_stack + 54915;
 target_ptr = Libderiv->deriv2_classes[4][5][79];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,6,1,dvrr_stack+32656, dvrr_stack+7230, dvrr_stack+16892);
 tmp = dvrr_stack + 32656;
 target_ptr = Libderiv->deriv2_classes[2][3][78];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+9277, dvrr_stack+4267, dvrr_stack+12275);
 tmp = dvrr_stack + 9277;
 target_ptr = Libderiv->deriv2_classes[2][4][78];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+32716, dvrr_stack+2097, dvrr_stack+7230);
 tmp = dvrr_stack + 32716;
 target_ptr = Libderiv->deriv2_classes[2][5][78];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,10,1,dvrr_stack+32842, dvrr_stack+3901, dvrr_stack+3841);
 tmp = dvrr_stack + 32842;
 target_ptr = Libderiv->deriv2_classes[3][3][78];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+55230, dvrr_stack+33076, dvrr_stack+4393);
 tmp = dvrr_stack + 55230;
 target_ptr = Libderiv->deriv2_classes[3][4][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+55380, dvrr_stack+33286, dvrr_stack+3901);
 tmp = dvrr_stack + 55380;
 target_ptr = Libderiv->deriv2_classes[3][5][78];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+33286, dvrr_stack+11245, dvrr_stack+11155);
 tmp = dvrr_stack + 33286;
 target_ptr = Libderiv->deriv2_classes[4][3][78];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+55590, dvrr_stack+33566, dvrr_stack+3021);
 tmp = dvrr_stack + 55590;
 target_ptr = Libderiv->deriv2_classes[4][4][78];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+55815, dvrr_stack+33881, dvrr_stack+11245);
 tmp = dvrr_stack + 55815;
 target_ptr = Libderiv->deriv2_classes[4][5][78];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,10,dvrr_stack+33881, dvrr_stack+25389, dvrr_stack+25299);
 tmp = dvrr_stack + 33881;
 target_ptr = Libderiv->deriv2_classes[2][3][35];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,15,dvrr_stack+11935, dvrr_stack+7805, dvrr_stack+10272);
 tmp = dvrr_stack + 11935;
 target_ptr = Libderiv->deriv2_classes[2][4][35];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_d(Data,21,dvrr_stack+33941, dvrr_stack+25489, dvrr_stack+23184);
 tmp = dvrr_stack + 33941;
 target_ptr = Libderiv->deriv2_classes[2][5][35];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,10,dvrr_stack+34067, dvrr_stack+26159, dvrr_stack+256);
 tmp = dvrr_stack + 34067;
 target_ptr = Libderiv->deriv2_classes[3][3][35];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,15,dvrr_stack+56130, dvrr_stack+17609, dvrr_stack+1520);
 tmp = dvrr_stack + 56130;
 target_ptr = Libderiv->deriv2_classes[3][4][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_f(Data,21,dvrr_stack+56280, dvrr_stack+21042, dvrr_stack+24777);
 tmp = dvrr_stack + 56280;
 target_ptr = Libderiv->deriv2_classes[3][5][35];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,10,dvrr_stack+56490, dvrr_stack+36171, dvrr_stack+25389);
 tmp = dvrr_stack + 56490;
 target_ptr = Libderiv->deriv2_classes[4][3][35];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,15,dvrr_stack+56640, dvrr_stack+38418, dvrr_stack+7805);
 tmp = dvrr_stack + 56640;
 target_ptr = Libderiv->deriv2_classes[4][4][35];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_g(Data,21,dvrr_stack+56865, dvrr_stack+18284, dvrr_stack+25489);
 tmp = dvrr_stack + 56865;
 target_ptr = Libderiv->deriv2_classes[4][5][35];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+34167, dvrr_stack+5011, dvrr_stack+18725);
 tmp = dvrr_stack + 34167;
 target_ptr = Libderiv->deriv2_classes[2][3][34];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+33436, dvrr_stack+4861, dvrr_stack+18755);
 tmp = dvrr_stack + 33436;
 target_ptr = Libderiv->deriv2_classes[2][4][34];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+32942, dvrr_stack+30891, dvrr_stack+18800);
 tmp = dvrr_stack + 32942;
 target_ptr = Libderiv->deriv2_classes[2][5][34];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+16824, dvrr_stack+5111, dvrr_stack+26345);
 tmp = dvrr_stack + 16824;
 target_ptr = Libderiv->deriv2_classes[3][3][34];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+57180, dvrr_stack+31381, dvrr_stack+21357);
 tmp = dvrr_stack + 57180;
 target_ptr = Libderiv->deriv2_classes[3][4][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+57330, dvrr_stack+31606, dvrr_stack+13225);
 tmp = dvrr_stack + 57330;
 target_ptr = Libderiv->deriv2_classes[3][5][34];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+57540, dvrr_stack+18863, dvrr_stack+5011);
 tmp = dvrr_stack + 57540;
 target_ptr = Libderiv->deriv2_classes[4][3][34];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+57690, dvrr_stack+8255, dvrr_stack+4861);
 tmp = dvrr_stack + 57690;
 target_ptr = Libderiv->deriv2_classes[4][4][34];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+57915, dvrr_stack+38733, dvrr_stack+30891);
 tmp = dvrr_stack + 57915;
 target_ptr = Libderiv->deriv2_classes[4][5][34];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+34227, dvrr_stack+6905, dvrr_stack+19073);
 tmp = dvrr_stack + 34227;
 target_ptr = Libderiv->deriv2_classes[2][3][33];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+11152, dvrr_stack+25149, dvrr_stack+22770);
 tmp = dvrr_stack + 11152;
 target_ptr = Libderiv->deriv2_classes[2][4][33];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+58230, dvrr_stack+2433, dvrr_stack+22815);
 tmp = dvrr_stack + 58230;
 target_ptr = Libderiv->deriv2_classes[2][5][33];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+58356, dvrr_stack+16259, dvrr_stack+980);
 tmp = dvrr_stack + 58356;
 target_ptr = Libderiv->deriv2_classes[3][3][33];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+58456, dvrr_stack+10047, dvrr_stack+24633);
 tmp = dvrr_stack + 58456;
 target_ptr = Libderiv->deriv2_classes[3][4][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+58606, dvrr_stack+15095, dvrr_stack+451);
 tmp = dvrr_stack + 58606;
 target_ptr = Libderiv->deriv2_classes[3][5][33];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+58816, dvrr_stack+22878, dvrr_stack+6905);
 tmp = dvrr_stack + 58816;
 target_ptr = Libderiv->deriv2_classes[4][3][33];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+8756, dvrr_stack+35541, dvrr_stack+25149);
 tmp = dvrr_stack + 8756;
 target_ptr = Libderiv->deriv2_classes[4][4][33];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+58966, dvrr_stack+37473, dvrr_stack+2433);
 tmp = dvrr_stack + 58966;
 target_ptr = Libderiv->deriv2_classes[4][5][33];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+2365, dvrr_stack+24045, dvrr_stack+40189);
 tmp = dvrr_stack + 2365;
 target_ptr = Libderiv->deriv2_classes[2][3][32];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+59281, dvrr_stack+23895, dvrr_stack+40219);
 tmp = dvrr_stack + 59281;
 target_ptr = Libderiv->deriv2_classes[2][4][32];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+59371, dvrr_stack+1040, dvrr_stack+40264);
 tmp = dvrr_stack + 59371;
 target_ptr = Libderiv->deriv2_classes[2][5][32];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+59497, dvrr_stack+19797, dvrr_stack+23481);
 tmp = dvrr_stack + 59497;
 target_ptr = Libderiv->deriv2_classes[3][3][32];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+59597, dvrr_stack+19572, dvrr_stack+23391);
 tmp = dvrr_stack + 59597;
 target_ptr = Libderiv->deriv2_classes[3][4][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+59747, dvrr_stack+16409, dvrr_stack+23541);
 tmp = dvrr_stack + 59747;
 target_ptr = Libderiv->deriv2_classes[3][5][32];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+19418, dvrr_stack+40498, dvrr_stack+24045);
 tmp = dvrr_stack + 19418;
 target_ptr = Libderiv->deriv2_classes[4][3][32];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+2000, dvrr_stack+35856, dvrr_stack+23895);
 tmp = dvrr_stack + 2000;
 target_ptr = Libderiv->deriv2_classes[4][4][32];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+59957, dvrr_stack+40708, dvrr_stack+1040);
 tmp = dvrr_stack + 59957;
 target_ptr = Libderiv->deriv2_classes[4][5][32];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+25003, dvrr_stack+12175, dvrr_stack+40327);
 tmp = dvrr_stack + 25003;
 target_ptr = Libderiv->deriv2_classes[2][3][31];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+60272, dvrr_stack+12025, dvrr_stack+25979);
 tmp = dvrr_stack + 60272;
 target_ptr = Libderiv->deriv2_classes[2][4][31];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+60362, dvrr_stack+9067, dvrr_stack+25071);
 tmp = dvrr_stack + 60362;
 target_ptr = Libderiv->deriv2_classes[2][5][31];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+60488, dvrr_stack+4117, dvrr_stack+20037);
 tmp = dvrr_stack + 60488;
 target_ptr = Libderiv->deriv2_classes[3][3][31];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+60588, dvrr_stack+7005, dvrr_stack+19947);
 tmp = dvrr_stack + 60588;
 target_ptr = Libderiv->deriv2_classes[3][4][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+60738, dvrr_stack+32341, dvrr_stack+15410);
 tmp = dvrr_stack + 60738;
 target_ptr = Libderiv->deriv2_classes[3][5][31];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+13005, dvrr_stack+41149, dvrr_stack+12175);
 tmp = dvrr_stack + 13005;
 target_ptr = Libderiv->deriv2_classes[4][3][31];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+60948, dvrr_stack+41359, dvrr_stack+12025);
 tmp = dvrr_stack + 60948;
 target_ptr = Libderiv->deriv2_classes[4][4][31];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+61173, dvrr_stack+41674, dvrr_stack+9067);
 tmp = dvrr_stack + 61173;
 target_ptr = Libderiv->deriv2_classes[4][5][31];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+6412, dvrr_stack+4393, dvrr_stack+26024);
 tmp = dvrr_stack + 6412;
 target_ptr = Libderiv->deriv2_classes[2][3][30];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+61488, dvrr_stack+3901, dvrr_stack+42115);
 tmp = dvrr_stack + 61488;
 target_ptr = Libderiv->deriv2_classes[2][4][30];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+61578, dvrr_stack+33076, dvrr_stack+42160);
 tmp = dvrr_stack + 61578;
 target_ptr = Libderiv->deriv2_classes[2][5][30];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+61704, dvrr_stack+3021, dvrr_stack+12275);
 tmp = dvrr_stack + 61704;
 target_ptr = Libderiv->deriv2_classes[3][3][30];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+61804, dvrr_stack+11245, dvrr_stack+7230);
 tmp = dvrr_stack + 61804;
 target_ptr = Libderiv->deriv2_classes[3][4][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+61954, dvrr_stack+33566, dvrr_stack+4267);
 tmp = dvrr_stack + 61954;
 target_ptr = Libderiv->deriv2_classes[3][5][30];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+62164, dvrr_stack+42223, dvrr_stack+4393);
 tmp = dvrr_stack + 62164;
 target_ptr = Libderiv->deriv2_classes[4][3][30];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+62314, dvrr_stack+37914, dvrr_stack+3901);
 tmp = dvrr_stack + 62314;
 target_ptr = Libderiv->deriv2_classes[4][4][30];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+62539, dvrr_stack+42433, dvrr_stack+33076);
 tmp = dvrr_stack + 62539;
 target_ptr = Libderiv->deriv2_classes[4][5][30];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,10,dvrr_stack+3813, dvrr_stack+40387, dvrr_stack+40357);
 tmp = dvrr_stack + 3813;
 target_ptr = Libderiv->deriv2_classes[2][3][26];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+62854, dvrr_stack+42919, dvrr_stack+42874);
 tmp = dvrr_stack + 62854;
 target_ptr = Libderiv->deriv2_classes[2][4][26];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+62944, dvrr_stack+43132, dvrr_stack+43069);
 tmp = dvrr_stack + 62944;
 target_ptr = Libderiv->deriv2_classes[2][5][26];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,10,dvrr_stack+63070, dvrr_stack+38229, dvrr_stack+43342);
 tmp = dvrr_stack + 63070;
 target_ptr = Libderiv->deriv2_classes[3][3][26];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+63170, dvrr_stack+39174, dvrr_stack+43402);
 tmp = dvrr_stack + 63170;
 target_ptr = Libderiv->deriv2_classes[3][4][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+63320, dvrr_stack+39399, dvrr_stack+8570);
 tmp = dvrr_stack + 63320;
 target_ptr = Libderiv->deriv2_classes[3][5][26];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+63530, dvrr_stack+34301, dvrr_stack+40387);
 tmp = dvrr_stack + 63530;
 target_ptr = Libderiv->deriv2_classes[4][3][26];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+63680, dvrr_stack+39714, dvrr_stack+42919);
 tmp = dvrr_stack + 63680;
 target_ptr = Libderiv->deriv2_classes[4][4][26];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+63905, dvrr_stack+34761, dvrr_stack+43132);
 tmp = dvrr_stack + 63905;
 target_ptr = Libderiv->deriv2_classes[4][5][26];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,10,dvrr_stack+64220, dvrr_stack+25389, dvrr_stack+25299);
 tmp = dvrr_stack + 64220;
 target_ptr = Libderiv->deriv2_classes[2][3][23];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,15,dvrr_stack+64280, dvrr_stack+7805, dvrr_stack+10272);
 tmp = dvrr_stack + 64280;
 target_ptr = Libderiv->deriv2_classes[2][4][23];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_d(Data,21,dvrr_stack+64370, dvrr_stack+25489, dvrr_stack+23184);
 tmp = dvrr_stack + 64370;
 target_ptr = Libderiv->deriv2_classes[2][5][23];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,10,dvrr_stack+64496, dvrr_stack+26159, dvrr_stack+256);
 tmp = dvrr_stack + 64496;
 target_ptr = Libderiv->deriv2_classes[3][3][23];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,15,dvrr_stack+64596, dvrr_stack+17609, dvrr_stack+1520);
 tmp = dvrr_stack + 64596;
 target_ptr = Libderiv->deriv2_classes[3][4][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_f(Data,21,dvrr_stack+64746, dvrr_stack+21042, dvrr_stack+24777);
 tmp = dvrr_stack + 64746;
 target_ptr = Libderiv->deriv2_classes[3][5][23];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,10,dvrr_stack+64956, dvrr_stack+36171, dvrr_stack+25389);
 tmp = dvrr_stack + 64956;
 target_ptr = Libderiv->deriv2_classes[4][3][23];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,15,dvrr_stack+65106, dvrr_stack+38418, dvrr_stack+7805);
 tmp = dvrr_stack + 65106;
 target_ptr = Libderiv->deriv2_classes[4][4][23];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_g(Data,21,dvrr_stack+65331, dvrr_stack+18284, dvrr_stack+25489);
 tmp = dvrr_stack + 65331;
 target_ptr = Libderiv->deriv2_classes[4][5][23];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+65646, dvrr_stack+5011, dvrr_stack+18725);
 tmp = dvrr_stack + 65646;
 target_ptr = Libderiv->deriv2_classes[2][3][22];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+65706, dvrr_stack+4861, dvrr_stack+18755);
 tmp = dvrr_stack + 65706;
 target_ptr = Libderiv->deriv2_classes[2][4][22];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+65796, dvrr_stack+30891, dvrr_stack+18800);
 tmp = dvrr_stack + 65796;
 target_ptr = Libderiv->deriv2_classes[2][5][22];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+26054, dvrr_stack+5111, dvrr_stack+26345);
 tmp = dvrr_stack + 26054;
 target_ptr = Libderiv->deriv2_classes[3][3][22];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+65922, dvrr_stack+31381, dvrr_stack+21357);
 tmp = dvrr_stack + 65922;
 target_ptr = Libderiv->deriv2_classes[3][4][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+66072, dvrr_stack+31606, dvrr_stack+13225);
 tmp = dvrr_stack + 66072;
 target_ptr = Libderiv->deriv2_classes[3][5][22];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+66282, dvrr_stack+18863, dvrr_stack+5011);
 tmp = dvrr_stack + 66282;
 target_ptr = Libderiv->deriv2_classes[4][3][22];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+66432, dvrr_stack+8255, dvrr_stack+4861);
 tmp = dvrr_stack + 66432;
 target_ptr = Libderiv->deriv2_classes[4][4][22];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+66657, dvrr_stack+38733, dvrr_stack+30891);
 tmp = dvrr_stack + 66657;
 target_ptr = Libderiv->deriv2_classes[4][5][22];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+4051, dvrr_stack+6905, dvrr_stack+19073);
 tmp = dvrr_stack + 4051;
 target_ptr = Libderiv->deriv2_classes[2][3][21];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+66972, dvrr_stack+25149, dvrr_stack+22770);
 tmp = dvrr_stack + 66972;
 target_ptr = Libderiv->deriv2_classes[2][4][21];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+67062, dvrr_stack+2433, dvrr_stack+22815);
 tmp = dvrr_stack + 67062;
 target_ptr = Libderiv->deriv2_classes[2][5][21];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+67188, dvrr_stack+16259, dvrr_stack+980);
 tmp = dvrr_stack + 67188;
 target_ptr = Libderiv->deriv2_classes[3][3][21];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+67288, dvrr_stack+10047, dvrr_stack+24633);
 tmp = dvrr_stack + 67288;
 target_ptr = Libderiv->deriv2_classes[3][4][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+67438, dvrr_stack+15095, dvrr_stack+451);
 tmp = dvrr_stack + 67438;
 target_ptr = Libderiv->deriv2_classes[3][5][21];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+67648, dvrr_stack+22878, dvrr_stack+6905);
 tmp = dvrr_stack + 67648;
 target_ptr = Libderiv->deriv2_classes[4][3][21];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+67798, dvrr_stack+35541, dvrr_stack+25149);
 tmp = dvrr_stack + 67798;
 target_ptr = Libderiv->deriv2_classes[4][4][21];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+68023, dvrr_stack+37473, dvrr_stack+2433);
 tmp = dvrr_stack + 68023;
 target_ptr = Libderiv->deriv2_classes[4][5][21];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+68338, dvrr_stack+24045, dvrr_stack+40189);
 tmp = dvrr_stack + 68338;
 target_ptr = Libderiv->deriv2_classes[2][3][20];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+68398, dvrr_stack+23895, dvrr_stack+40219);
 tmp = dvrr_stack + 68398;
 target_ptr = Libderiv->deriv2_classes[2][4][20];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+68488, dvrr_stack+1040, dvrr_stack+40264);
 tmp = dvrr_stack + 68488;
 target_ptr = Libderiv->deriv2_classes[2][5][20];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+68614, dvrr_stack+19797, dvrr_stack+23481);
 tmp = dvrr_stack + 68614;
 target_ptr = Libderiv->deriv2_classes[3][3][20];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+68714, dvrr_stack+19572, dvrr_stack+23391);
 tmp = dvrr_stack + 68714;
 target_ptr = Libderiv->deriv2_classes[3][4][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+68864, dvrr_stack+16409, dvrr_stack+23541);
 tmp = dvrr_stack + 68864;
 target_ptr = Libderiv->deriv2_classes[3][5][20];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+69074, dvrr_stack+40498, dvrr_stack+24045);
 tmp = dvrr_stack + 69074;
 target_ptr = Libderiv->deriv2_classes[4][3][20];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+69224, dvrr_stack+35856, dvrr_stack+23895);
 tmp = dvrr_stack + 69224;
 target_ptr = Libderiv->deriv2_classes[4][4][20];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+69449, dvrr_stack+40708, dvrr_stack+1040);
 tmp = dvrr_stack + 69449;
 target_ptr = Libderiv->deriv2_classes[4][5][20];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+69764, dvrr_stack+12175, dvrr_stack+40327);
 tmp = dvrr_stack + 69764;
 target_ptr = Libderiv->deriv2_classes[2][3][19];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+69824, dvrr_stack+12025, dvrr_stack+25979);
 tmp = dvrr_stack + 69824;
 target_ptr = Libderiv->deriv2_classes[2][4][19];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+69914, dvrr_stack+9067, dvrr_stack+25071);
 tmp = dvrr_stack + 69914;
 target_ptr = Libderiv->deriv2_classes[2][5][19];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+70040, dvrr_stack+4117, dvrr_stack+20037);
 tmp = dvrr_stack + 70040;
 target_ptr = Libderiv->deriv2_classes[3][3][19];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+70140, dvrr_stack+7005, dvrr_stack+19947);
 tmp = dvrr_stack + 70140;
 target_ptr = Libderiv->deriv2_classes[3][4][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+70290, dvrr_stack+32341, dvrr_stack+15410);
 tmp = dvrr_stack + 70290;
 target_ptr = Libderiv->deriv2_classes[3][5][19];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+70500, dvrr_stack+41149, dvrr_stack+12175);
 tmp = dvrr_stack + 70500;
 target_ptr = Libderiv->deriv2_classes[4][3][19];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+70650, dvrr_stack+41359, dvrr_stack+12025);
 tmp = dvrr_stack + 70650;
 target_ptr = Libderiv->deriv2_classes[4][4][19];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+70875, dvrr_stack+41674, dvrr_stack+9067);
 tmp = dvrr_stack + 70875;
 target_ptr = Libderiv->deriv2_classes[4][5][19];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+71190, dvrr_stack+4393, dvrr_stack+26024);
 tmp = dvrr_stack + 71190;
 target_ptr = Libderiv->deriv2_classes[2][3][18];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+71250, dvrr_stack+3901, dvrr_stack+42115);
 tmp = dvrr_stack + 71250;
 target_ptr = Libderiv->deriv2_classes[2][4][18];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+71340, dvrr_stack+33076, dvrr_stack+42160);
 tmp = dvrr_stack + 71340;
 target_ptr = Libderiv->deriv2_classes[2][5][18];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+71466, dvrr_stack+3021, dvrr_stack+12275);
 tmp = dvrr_stack + 71466;
 target_ptr = Libderiv->deriv2_classes[3][3][18];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+71566, dvrr_stack+11245, dvrr_stack+7230);
 tmp = dvrr_stack + 71566;
 target_ptr = Libderiv->deriv2_classes[3][4][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+71716, dvrr_stack+33566, dvrr_stack+4267);
 tmp = dvrr_stack + 71716;
 target_ptr = Libderiv->deriv2_classes[3][5][18];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+71926, dvrr_stack+42223, dvrr_stack+4393);
 tmp = dvrr_stack + 71926;
 target_ptr = Libderiv->deriv2_classes[4][3][18];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+72076, dvrr_stack+37914, dvrr_stack+3901);
 tmp = dvrr_stack + 72076;
 target_ptr = Libderiv->deriv2_classes[4][4][18];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+72301, dvrr_stack+42433, dvrr_stack+33076);
 tmp = dvrr_stack + 72301;
 target_ptr = Libderiv->deriv2_classes[4][5][18];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+72616, dvrr_stack+40387, dvrr_stack+40357);
 tmp = dvrr_stack + 72616;
 target_ptr = Libderiv->deriv2_classes[2][3][14];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+72676, dvrr_stack+42919, dvrr_stack+42874);
 tmp = dvrr_stack + 72676;
 target_ptr = Libderiv->deriv2_classes[2][4][14];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+72766, dvrr_stack+43132, dvrr_stack+43069);
 tmp = dvrr_stack + 72766;
 target_ptr = Libderiv->deriv2_classes[2][5][14];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+72892, dvrr_stack+38229, dvrr_stack+43342);
 tmp = dvrr_stack + 72892;
 target_ptr = Libderiv->deriv2_classes[3][3][14];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+72992, dvrr_stack+39174, dvrr_stack+43402);
 tmp = dvrr_stack + 72992;
 target_ptr = Libderiv->deriv2_classes[3][4][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+73142, dvrr_stack+39399, dvrr_stack+8570);
 tmp = dvrr_stack + 73142;
 target_ptr = Libderiv->deriv2_classes[3][5][14];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+73352, dvrr_stack+34301, dvrr_stack+40387);
 tmp = dvrr_stack + 73352;
 target_ptr = Libderiv->deriv2_classes[4][3][14];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+73502, dvrr_stack+39714, dvrr_stack+42919);
 tmp = dvrr_stack + 73502;
 target_ptr = Libderiv->deriv2_classes[4][4][14];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+73727, dvrr_stack+34761, dvrr_stack+43132);
 tmp = dvrr_stack + 73727;
 target_ptr = Libderiv->deriv2_classes[4][5][14];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,10,dvrr_stack+74042, dvrr_stack+44228, dvrr_stack+44198);
 tmp = dvrr_stack + 74042;
 target_ptr = Libderiv->deriv2_classes[2][3][13];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+74102, dvrr_stack+44328, dvrr_stack+11470);
 tmp = dvrr_stack + 74102;
 target_ptr = Libderiv->deriv2_classes[2][4][13];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+74192, dvrr_stack+3453, dvrr_stack+44478);
 tmp = dvrr_stack + 74192;
 target_ptr = Libderiv->deriv2_classes[2][5][13];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,10,dvrr_stack+74318, dvrr_stack+3663, dvrr_stack+44541);
 tmp = dvrr_stack + 74318;
 target_ptr = Libderiv->deriv2_classes[3][3][13];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+74418, dvrr_stack+36381, dvrr_stack+44601);
 tmp = dvrr_stack + 74418;
 target_ptr = Libderiv->deriv2_classes[3][4][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+74568, dvrr_stack+36606, dvrr_stack+5695);
 tmp = dvrr_stack + 74568;
 target_ptr = Libderiv->deriv2_classes[3][5][13];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+74778, dvrr_stack+43492, dvrr_stack+44228);
 tmp = dvrr_stack + 74778;
 target_ptr = Libderiv->deriv2_classes[4][3][13];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+74928, dvrr_stack+15536, dvrr_stack+44328);
 tmp = dvrr_stack + 74928;
 target_ptr = Libderiv->deriv2_classes[4][4][13];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+75153, dvrr_stack+45746, dvrr_stack+3453);
 tmp = dvrr_stack + 75153;
 target_ptr = Libderiv->deriv2_classes[4][5][13];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,10,dvrr_stack+75468, dvrr_stack+25389, dvrr_stack+25299);
 tmp = dvrr_stack + 75468;
 target_ptr = Libderiv->deriv2_classes[2][3][11];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,15,dvrr_stack+75528, dvrr_stack+7805, dvrr_stack+10272);
 tmp = dvrr_stack + 75528;
 target_ptr = Libderiv->deriv2_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_d(Data,21,dvrr_stack+75618, dvrr_stack+25489, dvrr_stack+23184);
 tmp = dvrr_stack + 75618;
 target_ptr = Libderiv->deriv2_classes[2][5][11];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,10,dvrr_stack+75744, dvrr_stack+26159, dvrr_stack+256);
 tmp = dvrr_stack + 75744;
 target_ptr = Libderiv->deriv2_classes[3][3][11];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,15,dvrr_stack+75844, dvrr_stack+17609, dvrr_stack+1520);
 tmp = dvrr_stack + 75844;
 target_ptr = Libderiv->deriv2_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_f(Data,21,dvrr_stack+17594, dvrr_stack+21042, dvrr_stack+24777);
 tmp = dvrr_stack + 17594;
 target_ptr = Libderiv->deriv2_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,10,dvrr_stack+20983, dvrr_stack+36171, dvrr_stack+25389);
 tmp = dvrr_stack + 20983;
 target_ptr = Libderiv->deriv2_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,15,dvrr_stack+75994, dvrr_stack+38418, dvrr_stack+7805);
 tmp = dvrr_stack + 75994;
 target_ptr = Libderiv->deriv2_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_g(Data,21,dvrr_stack+76219, dvrr_stack+18284, dvrr_stack+25489);
 tmp = dvrr_stack + 76219;
 target_ptr = Libderiv->deriv2_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+18284, dvrr_stack+5011, dvrr_stack+18725);
 tmp = dvrr_stack + 18284;
 target_ptr = Libderiv->deriv2_classes[2][3][10];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+18344, dvrr_stack+4861, dvrr_stack+18755);
 tmp = dvrr_stack + 18344;
 target_ptr = Libderiv->deriv2_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+18434, dvrr_stack+30891, dvrr_stack+18800);
 tmp = dvrr_stack + 18434;
 target_ptr = Libderiv->deriv2_classes[2][5][10];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+18560, dvrr_stack+5111, dvrr_stack+26345);
 tmp = dvrr_stack + 18560;
 target_ptr = Libderiv->deriv2_classes[3][3][10];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+5111, dvrr_stack+31381, dvrr_stack+21357);
 tmp = dvrr_stack + 5111;
 target_ptr = Libderiv->deriv2_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+36171, dvrr_stack+31606, dvrr_stack+13225);
 tmp = dvrr_stack + 36171;
 target_ptr = Libderiv->deriv2_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+18660, dvrr_stack+18863, dvrr_stack+5011);
 tmp = dvrr_stack + 18660;
 target_ptr = Libderiv->deriv2_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+18810, dvrr_stack+8255, dvrr_stack+4861);
 tmp = dvrr_stack + 18810;
 target_ptr = Libderiv->deriv2_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+8255, dvrr_stack+38733, dvrr_stack+30891);
 tmp = dvrr_stack + 8255;
 target_ptr = Libderiv->deriv2_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+30891, dvrr_stack+6905, dvrr_stack+19073);
 tmp = dvrr_stack + 30891;
 target_ptr = Libderiv->deriv2_classes[2][3][9];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+30951, dvrr_stack+25149, dvrr_stack+22770);
 tmp = dvrr_stack + 30951;
 target_ptr = Libderiv->deriv2_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+25389, dvrr_stack+2433, dvrr_stack+22815);
 tmp = dvrr_stack + 25389;
 target_ptr = Libderiv->deriv2_classes[2][5][9];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+25515, dvrr_stack+16259, dvrr_stack+980);
 tmp = dvrr_stack + 25515;
 target_ptr = Libderiv->deriv2_classes[3][3][9];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+31341, dvrr_stack+10047, dvrr_stack+24633);
 tmp = dvrr_stack + 31341;
 target_ptr = Libderiv->deriv2_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+24595, dvrr_stack+15095, dvrr_stack+451);
 tmp = dvrr_stack + 24595;
 target_ptr = Libderiv->deriv2_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+15095, dvrr_stack+22878, dvrr_stack+6905);
 tmp = dvrr_stack + 15095;
 target_ptr = Libderiv->deriv2_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+31491, dvrr_stack+35541, dvrr_stack+25149);
 tmp = dvrr_stack + 31491;
 target_ptr = Libderiv->deriv2_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+22768, dvrr_stack+37473, dvrr_stack+2433);
 tmp = dvrr_stack + 22768;
 target_ptr = Libderiv->deriv2_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+31041, dvrr_stack+24045, dvrr_stack+40189);
 tmp = dvrr_stack + 31041;
 target_ptr = Libderiv->deriv2_classes[2][3][8];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+15245, dvrr_stack+23895, dvrr_stack+40219);
 tmp = dvrr_stack + 15245;
 target_ptr = Libderiv->deriv2_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+2425, dvrr_stack+1040, dvrr_stack+40264);
 tmp = dvrr_stack + 2425;
 target_ptr = Libderiv->deriv2_classes[2][5][8];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+31716, dvrr_stack+19797, dvrr_stack+23481);
 tmp = dvrr_stack + 31716;
 target_ptr = Libderiv->deriv2_classes[3][3][8];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+19797, dvrr_stack+19572, dvrr_stack+23391);
 tmp = dvrr_stack + 19797;
 target_ptr = Libderiv->deriv2_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+21133, dvrr_stack+16409, dvrr_stack+23541);
 tmp = dvrr_stack + 21133;
 target_ptr = Libderiv->deriv2_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+424, dvrr_stack+40498, dvrr_stack+24045);
 tmp = dvrr_stack + 424;
 target_ptr = Libderiv->deriv2_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+35517, dvrr_stack+35856, dvrr_stack+23895);
 tmp = dvrr_stack + 35517;
 target_ptr = Libderiv->deriv2_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+35742, dvrr_stack+40708, dvrr_stack+1040);
 tmp = dvrr_stack + 35742;
 target_ptr = Libderiv->deriv2_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+15335, dvrr_stack+12175, dvrr_stack+40327);
 tmp = dvrr_stack + 15335;
 target_ptr = Libderiv->deriv2_classes[2][3][7];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+24805, dvrr_stack+12025, dvrr_stack+25979);
 tmp = dvrr_stack + 24805;
 target_ptr = Libderiv->deriv2_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+23893, dvrr_stack+9067, dvrr_stack+25071);
 tmp = dvrr_stack + 23893;
 target_ptr = Libderiv->deriv2_classes[2][5][7];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+24019, dvrr_stack+4117, dvrr_stack+20037);
 tmp = dvrr_stack + 24019;
 target_ptr = Libderiv->deriv2_classes[3][3][7];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+25063, dvrr_stack+7005, dvrr_stack+19947);
 tmp = dvrr_stack + 25063;
 target_ptr = Libderiv->deriv2_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+974, dvrr_stack+32341, dvrr_stack+15410);
 tmp = dvrr_stack + 974;
 target_ptr = Libderiv->deriv2_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+19947, dvrr_stack+41149, dvrr_stack+12175);
 tmp = dvrr_stack + 19947;
 target_ptr = Libderiv->deriv2_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+32287, dvrr_stack+41359, dvrr_stack+12025);
 tmp = dvrr_stack + 32287;
 target_ptr = Libderiv->deriv2_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+37447, dvrr_stack+41674, dvrr_stack+9067);
 tmp = dvrr_stack + 37447;
 target_ptr = Libderiv->deriv2_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+12025, dvrr_stack+4393, dvrr_stack+26024);
 tmp = dvrr_stack + 12025;
 target_ptr = Libderiv->deriv2_classes[2][3][6];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+12085, dvrr_stack+3901, dvrr_stack+42115);
 tmp = dvrr_stack + 12085;
 target_ptr = Libderiv->deriv2_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+15395, dvrr_stack+33076, dvrr_stack+42160);
 tmp = dvrr_stack + 15395;
 target_ptr = Libderiv->deriv2_classes[2][5][6];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+12175, dvrr_stack+3021, dvrr_stack+12275);
 tmp = dvrr_stack + 12175;
 target_ptr = Libderiv->deriv2_classes[3][3][6];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+3013, dvrr_stack+11245, dvrr_stack+7230);
 tmp = dvrr_stack + 3013;
 target_ptr = Libderiv->deriv2_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+11242, dvrr_stack+33566, dvrr_stack+4267);
 tmp = dvrr_stack + 11242;
 target_ptr = Libderiv->deriv2_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+33526, dvrr_stack+42223, dvrr_stack+4393);
 tmp = dvrr_stack + 33526;
 target_ptr = Libderiv->deriv2_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+23347, dvrr_stack+37914, dvrr_stack+3901);
 tmp = dvrr_stack + 23347;
 target_ptr = Libderiv->deriv2_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+37762, dvrr_stack+42433, dvrr_stack+33076);
 tmp = dvrr_stack + 37762;
 target_ptr = Libderiv->deriv2_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+12275, dvrr_stack+40387, dvrr_stack+40357);
 tmp = dvrr_stack + 12275;
 target_ptr = Libderiv->deriv2_classes[2][3][2];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+33676, dvrr_stack+42919, dvrr_stack+42874);
 tmp = dvrr_stack + 33676;
 target_ptr = Libderiv->deriv2_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+33068, dvrr_stack+43132, dvrr_stack+43069);
 tmp = dvrr_stack + 33068;
 target_ptr = Libderiv->deriv2_classes[2][5][2];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+33766, dvrr_stack+38229, dvrr_stack+43342);
 tmp = dvrr_stack + 33766;
 target_ptr = Libderiv->deriv2_classes[3][3][2];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+3873, dvrr_stack+39174, dvrr_stack+43402);
 tmp = dvrr_stack + 3873;
 target_ptr = Libderiv->deriv2_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+38077, dvrr_stack+39399, dvrr_stack+8570);
 tmp = dvrr_stack + 38077;
 target_ptr = Libderiv->deriv2_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+43342, dvrr_stack+34301, dvrr_stack+40387);
 tmp = dvrr_stack + 43342;
 target_ptr = Libderiv->deriv2_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+38287, dvrr_stack+39714, dvrr_stack+42919);
 tmp = dvrr_stack + 38287;
 target_ptr = Libderiv->deriv2_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+38512, dvrr_stack+34761, dvrr_stack+43132);
 tmp = dvrr_stack + 38512;
 target_ptr = Libderiv->deriv2_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+8570, dvrr_stack+44228, dvrr_stack+44198);
 tmp = dvrr_stack + 8570;
 target_ptr = Libderiv->deriv2_classes[2][3][1];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+34287, dvrr_stack+44328, dvrr_stack+11470);
 tmp = dvrr_stack + 34287;
 target_ptr = Libderiv->deriv2_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+34377, dvrr_stack+3453, dvrr_stack+44478);
 tmp = dvrr_stack + 34377;
 target_ptr = Libderiv->deriv2_classes[2][5][1];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+32512, dvrr_stack+3663, dvrr_stack+44541);
 tmp = dvrr_stack + 32512;
 target_ptr = Libderiv->deriv2_classes[3][3][1];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+3663, dvrr_stack+36381, dvrr_stack+44601);
 tmp = dvrr_stack + 3663;
 target_ptr = Libderiv->deriv2_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+36381, dvrr_stack+36606, dvrr_stack+5695);
 tmp = dvrr_stack + 36381;
 target_ptr = Libderiv->deriv2_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+36591, dvrr_stack+43492, dvrr_stack+44228);
 tmp = dvrr_stack + 36591;
 target_ptr = Libderiv->deriv2_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+38827, dvrr_stack+15536, dvrr_stack+44328);
 tmp = dvrr_stack + 38827;
 target_ptr = Libderiv->deriv2_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+15521, dvrr_stack+45746, dvrr_stack+3453);
 tmp = dvrr_stack + 15521;
 target_ptr = Libderiv->deriv2_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,10,dvrr_stack+43492, dvrr_stack+36921, dvrr_stack+44691);
 tmp = dvrr_stack + 43492;
 target_ptr = Libderiv->deriv2_classes[2][3][0];
 for(i=0;i<60;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+43552, dvrr_stack+43702, dvrr_stack+147);
 tmp = dvrr_stack + 43552;
 target_ptr = Libderiv->deriv2_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+36741, dvrr_stack+14585, dvrr_stack+361);
 tmp = dvrr_stack + 36741;
 target_ptr = Libderiv->deriv2_classes[2][5][0];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,10,dvrr_stack+31816, dvrr_stack+40029, dvrr_stack+854);
 tmp = dvrr_stack + 31816;
 target_ptr = Libderiv->deriv2_classes[3][3][0];
 for(i=0;i<100;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+252, dvrr_stack+34511, dvrr_stack+5821);
 tmp = dvrr_stack + 252;
 target_ptr = Libderiv->deriv2_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+34503, dvrr_stack+35202, dvrr_stack+5911);
 tmp = dvrr_stack + 34503;
 target_ptr = Libderiv->deriv2_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+34713, dvrr_stack+37021, dvrr_stack+36921);
 tmp = dvrr_stack + 34713;
 target_ptr = Libderiv->deriv2_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+36867, dvrr_stack+43852, dvrr_stack+43702);
 tmp = dvrr_stack + 36867;
 target_ptr = Libderiv->deriv2_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+43642, dvrr_stack+46187, dvrr_stack+14585);
 tmp = dvrr_stack + 43642;
 target_ptr = Libderiv->deriv2_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];


}

