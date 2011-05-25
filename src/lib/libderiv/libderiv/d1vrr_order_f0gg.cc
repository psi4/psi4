#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (f0|gg) integrals */

void d1vrr_order_f0gg(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+0, dvrr_stack+3, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+12, dvrr_stack+0, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+6, NULL, NULL, dvrr_stack+0);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+39, dvrr_stack+15, dvrr_stack+6, dvrr_stack+12, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+49, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+52, dvrr_stack+49, dvrr_stack+12, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+58, dvrr_stack+52, dvrr_stack+15, dvrr_stack+49, dvrr_stack+12, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+68, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+71, dvrr_stack+3, dvrr_stack+68, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+77, dvrr_stack+6, dvrr_stack+71, dvrr_stack+0, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+87, dvrr_stack+39, dvrr_stack+77, NULL, NULL, dvrr_stack+6);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+117, dvrr_stack+58, dvrr_stack+39, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+147, dvrr_stack+117, dvrr_stack+87, dvrr_stack+58, dvrr_stack+39, dvrr_stack+21);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+207, dvrr_stack+39, dvrr_stack+77, dvrr_stack+15, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+222, dvrr_stack+58, dvrr_stack+39, dvrr_stack+52, dvrr_stack+15, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+237, dvrr_stack+222, dvrr_stack+207, NULL, NULL, dvrr_stack+39);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+282, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+285, dvrr_stack+282, dvrr_stack+49, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+291, dvrr_stack+285, dvrr_stack+52, dvrr_stack+282, dvrr_stack+49, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+301, dvrr_stack+291, dvrr_stack+58, dvrr_stack+285, dvrr_stack+52, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+316, dvrr_stack+301, dvrr_stack+222, NULL, NULL, dvrr_stack+58);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+49, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+282, dvrr_stack+68, dvrr_stack+49, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+361, dvrr_stack+71, dvrr_stack+282, dvrr_stack+3, dvrr_stack+68, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+371, dvrr_stack+77, dvrr_stack+361, dvrr_stack+6, dvrr_stack+71, NULL);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+386, dvrr_stack+207, dvrr_stack+371, NULL, NULL, dvrr_stack+77);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+431, dvrr_stack+237, dvrr_stack+386, dvrr_stack+222, dvrr_stack+207, dvrr_stack+87);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+521, dvrr_stack+316, dvrr_stack+237, dvrr_stack+301, dvrr_stack+222, dvrr_stack+117);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+611, dvrr_stack+521, dvrr_stack+431, dvrr_stack+316, dvrr_stack+237, dvrr_stack+147);
 tmp = dvrr_stack + 611;
 target_ptr = Libderiv->dvrr_classes[3][4];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+316, dvrr_stack+207, dvrr_stack+371, dvrr_stack+39, dvrr_stack+77, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+337, dvrr_stack+222, dvrr_stack+207, dvrr_stack+58, dvrr_stack+39, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+761, dvrr_stack+337, dvrr_stack+316, NULL, NULL, dvrr_stack+207);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+824, dvrr_stack+301, dvrr_stack+222, dvrr_stack+291, dvrr_stack+58, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+845, dvrr_stack+824, dvrr_stack+337, NULL, NULL, dvrr_stack+222);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+358, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+908, dvrr_stack+49, dvrr_stack+358, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+914, dvrr_stack+282, dvrr_stack+908, dvrr_stack+68, dvrr_stack+49, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+924, dvrr_stack+361, dvrr_stack+914, dvrr_stack+71, dvrr_stack+282, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+939, dvrr_stack+371, dvrr_stack+924, dvrr_stack+77, dvrr_stack+361, NULL);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+960, dvrr_stack+316, dvrr_stack+939, NULL, NULL, dvrr_stack+371);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1023, dvrr_stack+761, dvrr_stack+960, dvrr_stack+337, dvrr_stack+316, dvrr_stack+386);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1149, dvrr_stack+845, dvrr_stack+761, dvrr_stack+824, dvrr_stack+337, dvrr_stack+237);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1275, dvrr_stack+1149, dvrr_stack+1023, dvrr_stack+845, dvrr_stack+761, dvrr_stack+431);
 tmp = dvrr_stack + 1275;
 target_ptr = Libderiv->dvrr_classes[3][5];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+1485,dvrr_stack+1275,dvrr_stack+611,10);


 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+845, dvrr_stack+316, dvrr_stack+939, dvrr_stack+207, dvrr_stack+371, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+873, dvrr_stack+337, dvrr_stack+316, dvrr_stack+222, dvrr_stack+207, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1935, dvrr_stack+873, dvrr_stack+845, NULL, NULL, dvrr_stack+316);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2019, dvrr_stack+824, dvrr_stack+337, dvrr_stack+301, dvrr_stack+222, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2047, dvrr_stack+2019, dvrr_stack+873, NULL, NULL, dvrr_stack+337);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+68, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+222, dvrr_stack+358, dvrr_stack+68, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+301, dvrr_stack+908, dvrr_stack+222, dvrr_stack+49, dvrr_stack+358, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2131, dvrr_stack+914, dvrr_stack+301, dvrr_stack+282, dvrr_stack+908, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2146, dvrr_stack+924, dvrr_stack+2131, dvrr_stack+361, dvrr_stack+914, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2167, dvrr_stack+939, dvrr_stack+2146, dvrr_stack+371, dvrr_stack+924, NULL);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2195, dvrr_stack+845, dvrr_stack+2167, NULL, NULL, dvrr_stack+939);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2279, dvrr_stack+1935, dvrr_stack+2195, dvrr_stack+873, dvrr_stack+845, dvrr_stack+960);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2447, dvrr_stack+2047, dvrr_stack+1935, dvrr_stack+2019, dvrr_stack+873, dvrr_stack+761);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2615, dvrr_stack+2447, dvrr_stack+2279, dvrr_stack+2047, dvrr_stack+1935, dvrr_stack+1023);
 tmp = dvrr_stack + 2615;
 target_ptr = Libderiv->dvrr_classes[3][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+2895,dvrr_stack+2615,dvrr_stack+1275,10);


 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2047, dvrr_stack+845, dvrr_stack+2167, dvrr_stack+316, dvrr_stack+939, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2083, dvrr_stack+873, dvrr_stack+845, dvrr_stack+337, dvrr_stack+316, NULL);

 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3525, dvrr_stack+2083, dvrr_stack+2047, NULL, NULL, dvrr_stack+845);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3633, dvrr_stack+2019, dvrr_stack+873, dvrr_stack+824, dvrr_stack+337, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3669, dvrr_stack+3633, dvrr_stack+2083, NULL, NULL, dvrr_stack+873);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+49, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+337, dvrr_stack+68, dvrr_stack+49, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+343, dvrr_stack+222, dvrr_stack+337, dvrr_stack+358, dvrr_stack+68, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+824, dvrr_stack+301, dvrr_stack+343, dvrr_stack+908, dvrr_stack+222, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3777, dvrr_stack+2131, dvrr_stack+824, dvrr_stack+914, dvrr_stack+301, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3798, dvrr_stack+2146, dvrr_stack+3777, dvrr_stack+924, dvrr_stack+2131, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3826, dvrr_stack+2167, dvrr_stack+3798, dvrr_stack+939, dvrr_stack+2146, NULL);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3862, dvrr_stack+2047, dvrr_stack+3826, NULL, NULL, dvrr_stack+2167);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3970, dvrr_stack+3525, dvrr_stack+3862, dvrr_stack+2083, dvrr_stack+2047, dvrr_stack+2195);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4186, dvrr_stack+3669, dvrr_stack+3525, dvrr_stack+3633, dvrr_stack+2083, dvrr_stack+1935);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4402, dvrr_stack+4186, dvrr_stack+3970, dvrr_stack+3669, dvrr_stack+3525, dvrr_stack+2279);
 tmp = dvrr_stack + 4402;
 target_ptr = Libderiv->dvrr_classes[3][7];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+4762,dvrr_stack+4402,dvrr_stack+2615,10);


 /* compute (0 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3669, dvrr_stack+2047, dvrr_stack+3826, dvrr_stack+845, dvrr_stack+2167, NULL);

 /* compute (0 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3714, dvrr_stack+2083, dvrr_stack+2047, dvrr_stack+873, dvrr_stack+845, NULL);

 /* compute (1 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+5602, dvrr_stack+3714, dvrr_stack+3669, NULL, NULL, dvrr_stack+2047);

 /* compute (0 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+5737, dvrr_stack+3633, dvrr_stack+2083, dvrr_stack+2019, dvrr_stack+873, NULL);

 /* compute (1 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+5782, dvrr_stack+5737, dvrr_stack+3714, NULL, NULL, dvrr_stack+2083);

 /* compute (0 0 | 1 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+873, Data->F+10, Data->F+11, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+839, dvrr_stack+49, dvrr_stack+873, Data->F+9, Data->F+10, NULL);

 /* compute (0 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+876, dvrr_stack+337, dvrr_stack+839, dvrr_stack+68, dvrr_stack+49, NULL);

 /* compute (0 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+886, dvrr_stack+343, dvrr_stack+876, dvrr_stack+222, dvrr_stack+337, NULL);

 /* compute (0 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+901, dvrr_stack+824, dvrr_stack+886, dvrr_stack+301, dvrr_stack+343, NULL);

 /* compute (0 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2019, dvrr_stack+3777, dvrr_stack+901, dvrr_stack+2131, dvrr_stack+824, NULL);

 /* compute (0 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+5917, dvrr_stack+3798, dvrr_stack+2019, dvrr_stack+2146, dvrr_stack+3777, NULL);

 /* compute (0 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+5953, dvrr_stack+3826, dvrr_stack+5917, dvrr_stack+2167, dvrr_stack+3798, NULL);

 /* compute (1 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+5998, dvrr_stack+3669, dvrr_stack+5953, NULL, NULL, dvrr_stack+3826);

 /* compute (2 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6133, dvrr_stack+5602, dvrr_stack+5998, dvrr_stack+3714, dvrr_stack+3669, dvrr_stack+3862);

 /* compute (2 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6403, dvrr_stack+5782, dvrr_stack+5602, dvrr_stack+5737, dvrr_stack+3714, dvrr_stack+3525);

 /* compute (3 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+6673, dvrr_stack+6403, dvrr_stack+6133, dvrr_stack+5782, dvrr_stack+5602, dvrr_stack+3970);

 /* compute (3 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+7123,dvrr_stack+6673,dvrr_stack+4402,10);


 /* compute (0 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+5782, dvrr_stack+3669, dvrr_stack+5953, dvrr_stack+2047, dvrr_stack+3826, NULL);

 /* compute (0 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+5837, dvrr_stack+3714, dvrr_stack+3669, dvrr_stack+2083, dvrr_stack+2047, NULL);

 /* compute (1 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+8203, dvrr_stack+5837, dvrr_stack+5782, NULL, NULL, dvrr_stack+3669);

 /* compute (0 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+8368, dvrr_stack+5737, dvrr_stack+3714, dvrr_stack+3633, dvrr_stack+2083, NULL);

 /* compute (1 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+8423, dvrr_stack+8368, dvrr_stack+5837, NULL, NULL, dvrr_stack+3714);

 /* compute (0 0 | 1 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+68, Data->F+11, Data->F+12, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+3714, dvrr_stack+873, dvrr_stack+68, Data->F+10, Data->F+11, NULL);

 /* compute (0 0 | 3 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+3720, dvrr_stack+839, dvrr_stack+3714, dvrr_stack+49, dvrr_stack+873, NULL);

 /* compute (0 0 | 4 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+301, dvrr_stack+876, dvrr_stack+3720, dvrr_stack+337, dvrr_stack+839, NULL);

 /* compute (0 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3714, dvrr_stack+886, dvrr_stack+301, dvrr_stack+343, dvrr_stack+876, NULL);

 /* compute (0 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3735, dvrr_stack+901, dvrr_stack+3714, dvrr_stack+824, dvrr_stack+886, NULL);

 /* compute (0 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3633, dvrr_stack+2019, dvrr_stack+3735, dvrr_stack+3777, dvrr_stack+901, NULL);

 /* compute (0 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+5737, dvrr_stack+5917, dvrr_stack+3633, dvrr_stack+3798, dvrr_stack+2019, NULL);

 /* compute (0 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+3714, dvrr_stack+5953, dvrr_stack+5737, dvrr_stack+3826, dvrr_stack+5917, NULL);

 /* compute (1 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+8588, dvrr_stack+5782, dvrr_stack+3714, NULL, NULL, dvrr_stack+5953);

 /* compute (2 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+8753, dvrr_stack+8203, dvrr_stack+8588, dvrr_stack+5837, dvrr_stack+5782, dvrr_stack+5998);

 /* compute (2 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+9083, dvrr_stack+8423, dvrr_stack+8203, dvrr_stack+8368, dvrr_stack+5837, dvrr_stack+5602);

 /* compute (3 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+9413, dvrr_stack+9083, dvrr_stack+8753, dvrr_stack+8423, dvrr_stack+8203, dvrr_stack+6133);

 /* compute (3 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+9963,dvrr_stack+9413,dvrr_stack+6673,10);


 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+282, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+3);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+8203, dvrr_stack+52, dvrr_stack+15, NULL, NULL, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3633, dvrr_stack+8203, dvrr_stack+21, dvrr_stack+52, dvrr_stack+15, dvrr_stack+282);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+8221, dvrr_stack+291, dvrr_stack+58, NULL, NULL, dvrr_stack+52);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+8251, dvrr_stack+8221, dvrr_stack+117, dvrr_stack+291, dvrr_stack+58, dvrr_stack+8203);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+8311, dvrr_stack+8251, dvrr_stack+147, dvrr_stack+8221, dvrr_stack+117, dvrr_stack+3633);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3633, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+4);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+3642, dvrr_stack+6, dvrr_stack+71, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+8203, dvrr_stack+21, dvrr_stack+3642, dvrr_stack+15, dvrr_stack+6, dvrr_stack+3633);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+0, dvrr_stack+77, dvrr_stack+361, NULL, NULL, dvrr_stack+71);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+8239, dvrr_stack+87, dvrr_stack+0, dvrr_stack+39, dvrr_stack+77, dvrr_stack+3642);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+8411, dvrr_stack+147, dvrr_stack+8239, dvrr_stack+117, dvrr_stack+87, dvrr_stack+8203);

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+30, dvrr_stack+371, dvrr_stack+924, NULL, NULL, dvrr_stack+361);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+75, dvrr_stack+386, dvrr_stack+30, dvrr_stack+207, dvrr_stack+371, dvrr_stack+0);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8511, dvrr_stack+431, dvrr_stack+75, dvrr_stack+237, dvrr_stack+386, dvrr_stack+8239);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+8661, dvrr_stack+611, dvrr_stack+8511, dvrr_stack+521, dvrr_stack+431, dvrr_stack+8411);

 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2083, dvrr_stack+939, dvrr_stack+2146, NULL, NULL, dvrr_stack+924);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+165, dvrr_stack+960, dvrr_stack+2083, dvrr_stack+316, dvrr_stack+939, dvrr_stack+30);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+291, dvrr_stack+1023, dvrr_stack+165, dvrr_stack+761, dvrr_stack+960, dvrr_stack+75);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+8886, dvrr_stack+1275, dvrr_stack+291, dvrr_stack+1149, dvrr_stack+1023, dvrr_stack+8511);

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+761, dvrr_stack+2167, dvrr_stack+3798, NULL, NULL, dvrr_stack+2146);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+8411, dvrr_stack+2195, dvrr_stack+761, dvrr_stack+845, dvrr_stack+2167, dvrr_stack+2083);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+845, dvrr_stack+2279, dvrr_stack+8411, dvrr_stack+1935, dvrr_stack+2195, dvrr_stack+165);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+11313, dvrr_stack+2615, dvrr_stack+845, dvrr_stack+2447, dvrr_stack+2279, dvrr_stack+291);

 /* compute (1 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+8203, dvrr_stack+3826, dvrr_stack+5917, NULL, NULL, dvrr_stack+3798);

 /* compute (2 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2083, dvrr_stack+3862, dvrr_stack+8203, dvrr_stack+2047, dvrr_stack+3826, dvrr_stack+761);

 /* compute (3 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+3970, dvrr_stack+2083, dvrr_stack+3525, dvrr_stack+3862, dvrr_stack+8411);

 /* compute (4 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+11733, dvrr_stack+4402, dvrr_stack+0, dvrr_stack+4186, dvrr_stack+3970, dvrr_stack+845);

 /* compute (1 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+5782, dvrr_stack+5953, dvrr_stack+5737, NULL, NULL, dvrr_stack+5917);

 /* compute (2 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+761, dvrr_stack+5998, dvrr_stack+5782, dvrr_stack+3669, dvrr_stack+5953, dvrr_stack+8203);

 /* compute (3 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3525, dvrr_stack+6133, dvrr_stack+761, dvrr_stack+5602, dvrr_stack+5998, dvrr_stack+2083);

 /* compute (4 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+12273, dvrr_stack+6673, dvrr_stack+3525, dvrr_stack+6403, dvrr_stack+6133, dvrr_stack+0);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+0, dvrr_stack+1485, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+150, dvrr_stack+2895, NULL);
 tmp = dvrr_stack + 150;
 target_ptr = Libderiv->deriv_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,280,dvrr_stack+3525, dvrr_stack+4762, NULL);
 tmp = dvrr_stack + 3525;
 target_ptr = Libderiv->deriv_classes[3][6][11];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,360,dvrr_stack+3805, dvrr_stack+7123, NULL);
 tmp = dvrr_stack + 3805;
 target_ptr = Libderiv->deriv_classes[3][7][11];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,450,dvrr_stack+5602, dvrr_stack+9963, NULL);
 tmp = dvrr_stack + 5602;
 target_ptr = Libderiv->deriv_classes[3][8][11];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+360, dvrr_stack+1485, NULL);
 tmp = dvrr_stack + 360;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+6052, dvrr_stack+2895, NULL);
 tmp = dvrr_stack + 6052;
 target_ptr = Libderiv->deriv_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+761, dvrr_stack+4762, NULL);
 tmp = dvrr_stack + 761;
 target_ptr = Libderiv->deriv_classes[3][6][10];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,360,dvrr_stack+1935, dvrr_stack+7123, NULL);
 tmp = dvrr_stack + 1935;
 target_ptr = Libderiv->deriv_classes[3][7][10];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,450,dvrr_stack+12948, dvrr_stack+9963, NULL);
 tmp = dvrr_stack + 12948;
 target_ptr = Libderiv->deriv_classes[3][8][10];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+8411, dvrr_stack+1485, NULL);
 tmp = dvrr_stack + 8411;
 target_ptr = Libderiv->deriv_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+1485, dvrr_stack+2895, NULL);
 tmp = dvrr_stack + 1485;
 target_ptr = Libderiv->deriv_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+2895, dvrr_stack+4762, NULL);
 tmp = dvrr_stack + 2895;
 target_ptr = Libderiv->deriv_classes[3][6][9];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,360,dvrr_stack+4762, dvrr_stack+7123, NULL);
 tmp = dvrr_stack + 4762;
 target_ptr = Libderiv->deriv_classes[3][7][9];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,450,dvrr_stack+7123, dvrr_stack+9963, NULL);
 tmp = dvrr_stack + 7123;
 target_ptr = Libderiv->deriv_classes[3][8][9];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+9963, dvrr_stack+1275, dvrr_stack+8311);
 tmp = dvrr_stack + 9963;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+10113, dvrr_stack+2615, dvrr_stack+611);
 tmp = dvrr_stack + 10113;
 target_ptr = Libderiv->deriv_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+10323, dvrr_stack+4402, dvrr_stack+1275);
 tmp = dvrr_stack + 10323;
 target_ptr = Libderiv->deriv_classes[3][6][8];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,10,1,dvrr_stack+10603, dvrr_stack+6673, dvrr_stack+2615);
 tmp = dvrr_stack + 10603;
 target_ptr = Libderiv->deriv_classes[3][7][8];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,10,1,dvrr_stack+7573, dvrr_stack+9413, dvrr_stack+4402);
 tmp = dvrr_stack + 7573;
 target_ptr = Libderiv->deriv_classes[3][8][8];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+10963, dvrr_stack+1275, dvrr_stack+8311);
 tmp = dvrr_stack + 10963;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+8023, dvrr_stack+2615, dvrr_stack+611);
 tmp = dvrr_stack + 8023;
 target_ptr = Libderiv->deriv_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+5122, dvrr_stack+4402, dvrr_stack+1275);
 tmp = dvrr_stack + 5122;
 target_ptr = Libderiv->deriv_classes[3][6][7];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,10,1,dvrr_stack+13398, dvrr_stack+6673, dvrr_stack+2615);
 tmp = dvrr_stack + 13398;
 target_ptr = Libderiv->deriv_classes[3][7][7];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,10,1,dvrr_stack+13758, dvrr_stack+9413, dvrr_stack+4402);
 tmp = dvrr_stack + 13758;
 target_ptr = Libderiv->deriv_classes[3][8][7];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+11113, dvrr_stack+1275, dvrr_stack+8311);
 tmp = dvrr_stack + 11113;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+3175, dvrr_stack+2615, dvrr_stack+611);
 tmp = dvrr_stack + 3175;
 target_ptr = Libderiv->deriv_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+14208, dvrr_stack+4402, dvrr_stack+1275);
 tmp = dvrr_stack + 14208;
 target_ptr = Libderiv->deriv_classes[3][6][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,10,1,dvrr_stack+14488, dvrr_stack+6673, dvrr_stack+2615);
 tmp = dvrr_stack + 14488;
 target_ptr = Libderiv->deriv_classes[3][7][6];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,10,1,dvrr_stack+6673, dvrr_stack+9413, dvrr_stack+4402);
 tmp = dvrr_stack + 6673;
 target_ptr = Libderiv->deriv_classes[3][8][6];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+611, dvrr_stack+8661, dvrr_stack+521);
 tmp = dvrr_stack + 611;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+1275, dvrr_stack+8886, dvrr_stack+1149);
 tmp = dvrr_stack + 1275;
 target_ptr = Libderiv->deriv_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+2615, dvrr_stack+11313, dvrr_stack+2447);
 tmp = dvrr_stack + 2615;
 target_ptr = Libderiv->deriv_classes[3][6][2];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,36,dvrr_stack+4402, dvrr_stack+11733, dvrr_stack+4186);
 tmp = dvrr_stack + 4402;
 target_ptr = Libderiv->deriv_classes[3][7][2];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,45,dvrr_stack+9201, dvrr_stack+12273, dvrr_stack+6403);
 tmp = dvrr_stack + 9201;
 target_ptr = Libderiv->deriv_classes[3][8][2];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+8233, dvrr_stack+8661, dvrr_stack+521);
 tmp = dvrr_stack + 8233;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+1695, dvrr_stack+8886, dvrr_stack+1149);
 tmp = dvrr_stack + 1695;
 target_ptr = Libderiv->deriv_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+9651, dvrr_stack+11313, dvrr_stack+2447);
 tmp = dvrr_stack + 9651;
 target_ptr = Libderiv->deriv_classes[3][6][1];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,36,dvrr_stack+14848, dvrr_stack+11733, dvrr_stack+4186);
 tmp = dvrr_stack + 14848;
 target_ptr = Libderiv->deriv_classes[3][7][1];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,45,dvrr_stack+15208, dvrr_stack+12273, dvrr_stack+6403);
 tmp = dvrr_stack + 15208;
 target_ptr = Libderiv->deriv_classes[3][8][1];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+5402, dvrr_stack+8661, dvrr_stack+521);
 tmp = dvrr_stack + 5402;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+8561, dvrr_stack+8886, dvrr_stack+1149);
 tmp = dvrr_stack + 8561;
 target_ptr = Libderiv->deriv_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+8771, dvrr_stack+11313, dvrr_stack+2447);
 tmp = dvrr_stack + 8771;
 target_ptr = Libderiv->deriv_classes[3][6][0];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,36,dvrr_stack+11263, dvrr_stack+11733, dvrr_stack+4186);
 tmp = dvrr_stack + 11263;
 target_ptr = Libderiv->deriv_classes[3][7][0];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,45,dvrr_stack+11623, dvrr_stack+12273, dvrr_stack+6403);
 tmp = dvrr_stack + 11623;
 target_ptr = Libderiv->deriv_classes[3][8][0];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];


}

