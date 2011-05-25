#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (ff|gd) integrals */

void d1vrr_order_ffgd(Libderiv_t *Libderiv, prim_data *Data)
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
 _BUILD_00p0(Data,dvrr_stack+288, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+222, dvrr_stack+358, dvrr_stack+288, Data->F+7, Data->F+8, NULL);

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
 _BUILD_00p0(Data,dvrr_stack+873, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+876, dvrr_stack+288, dvrr_stack+873, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+882, dvrr_stack+222, dvrr_stack+876, dvrr_stack+358, dvrr_stack+288, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+892, dvrr_stack+301, dvrr_stack+882, dvrr_stack+908, dvrr_stack+222, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+337, dvrr_stack+2131, dvrr_stack+892, dvrr_stack+914, dvrr_stack+301, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2019, dvrr_stack+2146, dvrr_stack+337, dvrr_stack+924, dvrr_stack+2131, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3777, dvrr_stack+2167, dvrr_stack+2019, dvrr_stack+939, dvrr_stack+2146, NULL);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3813, dvrr_stack+2047, dvrr_stack+3777, NULL, NULL, dvrr_stack+2167);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3921, dvrr_stack+3525, dvrr_stack+3813, dvrr_stack+2083, dvrr_stack+2047, dvrr_stack+2195);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4137, dvrr_stack+3669, dvrr_stack+3525, dvrr_stack+3633, dvrr_stack+2083, dvrr_stack+1935);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4353, dvrr_stack+4137, dvrr_stack+3921, dvrr_stack+3669, dvrr_stack+3525, dvrr_stack+2279);

 /* compute (3 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+4713,dvrr_stack+4353,dvrr_stack+2615,10);


 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+228, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+4);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+2083, dvrr_stack+6, dvrr_stack+71, NULL, NULL, dvrr_stack+3);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3633, dvrr_stack+21, dvrr_stack+2083, dvrr_stack+15, dvrr_stack+6, dvrr_stack+228);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+2101, dvrr_stack+77, dvrr_stack+361, NULL, NULL, dvrr_stack+71);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+3669, dvrr_stack+87, dvrr_stack+2101, dvrr_stack+39, dvrr_stack+77, dvrr_stack+2083);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+5553, dvrr_stack+147, dvrr_stack+3669, dvrr_stack+117, dvrr_stack+87, dvrr_stack+3633);

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3729, dvrr_stack+371, dvrr_stack+924, NULL, NULL, dvrr_stack+361);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5653, dvrr_stack+386, dvrr_stack+3729, dvrr_stack+207, dvrr_stack+371, dvrr_stack+2101);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5743, dvrr_stack+431, dvrr_stack+5653, dvrr_stack+237, dvrr_stack+386, dvrr_stack+3669);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+5893, dvrr_stack+611, dvrr_stack+5743, dvrr_stack+521, dvrr_stack+431, dvrr_stack+5553);
 tmp = dvrr_stack + 5893;
 target_ptr = Libderiv->dvrr_classes[4][4];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6118, dvrr_stack+939, dvrr_stack+2146, NULL, NULL, dvrr_stack+924);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6181, dvrr_stack+960, dvrr_stack+6118, dvrr_stack+316, dvrr_stack+939, dvrr_stack+3729);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6307, dvrr_stack+1023, dvrr_stack+6181, dvrr_stack+761, dvrr_stack+960, dvrr_stack+5653);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6517, dvrr_stack+1275, dvrr_stack+6307, dvrr_stack+1149, dvrr_stack+1023, dvrr_stack+5743);
 tmp = dvrr_stack + 6517;
 target_ptr = Libderiv->dvrr_classes[4][5];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+6832,dvrr_stack+6517,dvrr_stack+5893,15);


 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+761, dvrr_stack+2167, dvrr_stack+2019, NULL, NULL, dvrr_stack+2146);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+7507, dvrr_stack+2195, dvrr_stack+761, dvrr_stack+845, dvrr_stack+2167, dvrr_stack+6118);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+7675, dvrr_stack+2279, dvrr_stack+7507, dvrr_stack+1935, dvrr_stack+2195, dvrr_stack+6181);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+7955, dvrr_stack+2615, dvrr_stack+7675, dvrr_stack+2447, dvrr_stack+2279, dvrr_stack+6307);
 tmp = dvrr_stack + 7955;
 target_ptr = Libderiv->dvrr_classes[4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+8375,dvrr_stack+7955,dvrr_stack+6517,15);


 /* compute (0 0 | 1 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3774, Data->F+10, Data->F+11, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1935, dvrr_stack+873, dvrr_stack+3774, Data->F+9, Data->F+10, NULL);

 /* compute (0 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+39, dvrr_stack+876, dvrr_stack+1935, dvrr_stack+288, dvrr_stack+873, NULL);

 /* compute (0 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+207, dvrr_stack+882, dvrr_stack+39, dvrr_stack+222, dvrr_stack+876, NULL);

 /* compute (0 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1941, dvrr_stack+892, dvrr_stack+207, dvrr_stack+301, dvrr_stack+882, NULL);

 /* compute (0 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+845, dvrr_stack+337, dvrr_stack+1941, dvrr_stack+2131, dvrr_stack+892, NULL);

 /* compute (0 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1962, dvrr_stack+2019, dvrr_stack+845, dvrr_stack+2146, dvrr_stack+337, NULL);

 /* compute (1 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+9320, dvrr_stack+3777, dvrr_stack+1962, NULL, NULL, dvrr_stack+2019);

 /* compute (2 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+9428, dvrr_stack+3813, dvrr_stack+9320, dvrr_stack+2047, dvrr_stack+3777, dvrr_stack+761);

 /* compute (3 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+9644, dvrr_stack+3921, dvrr_stack+9428, dvrr_stack+3525, dvrr_stack+3813, dvrr_stack+7507);

 /* compute (4 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+10004, dvrr_stack+4353, dvrr_stack+9644, dvrr_stack+4137, dvrr_stack+3921, dvrr_stack+7675);

 /* compute (4 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+10544,dvrr_stack+10004,dvrr_stack+7955,15);


 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+288, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+4137, dvrr_stack+3, dvrr_stack+68, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+4146, dvrr_stack+228, dvrr_stack+4137, dvrr_stack+0, dvrr_stack+3, dvrr_stack+288);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+4164, dvrr_stack+71, dvrr_stack+282, NULL, NULL, dvrr_stack+68);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+2047, dvrr_stack+2083, dvrr_stack+4164, dvrr_stack+6, dvrr_stack+71, dvrr_stack+4137);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+4182, dvrr_stack+3633, dvrr_stack+2047, dvrr_stack+21, dvrr_stack+2083, dvrr_stack+4146);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+4242, dvrr_stack+361, dvrr_stack+914, NULL, NULL, dvrr_stack+282);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4272, dvrr_stack+2101, dvrr_stack+4242, dvrr_stack+77, dvrr_stack+361, dvrr_stack+4164);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+3525, dvrr_stack+3669, dvrr_stack+4272, dvrr_stack+87, dvrr_stack+2101, dvrr_stack+2047);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+11804, dvrr_stack+5553, dvrr_stack+3525, dvrr_stack+147, dvrr_stack+3669, dvrr_stack+4182);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+237, dvrr_stack+924, dvrr_stack+2131, NULL, NULL, dvrr_stack+914);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11954, dvrr_stack+3729, dvrr_stack+237, dvrr_stack+371, dvrr_stack+924, dvrr_stack+4242);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+12044, dvrr_stack+5653, dvrr_stack+11954, dvrr_stack+386, dvrr_stack+3729, dvrr_stack+4272);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+12194, dvrr_stack+5743, dvrr_stack+12044, dvrr_stack+431, dvrr_stack+5653, dvrr_stack+3525);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+12419, dvrr_stack+5893, dvrr_stack+12194, dvrr_stack+611, dvrr_stack+5743, dvrr_stack+11804);
 tmp = dvrr_stack + 12419;
 target_ptr = Libderiv->dvrr_classes[5][4];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+371, dvrr_stack+2146, dvrr_stack+337, NULL, NULL, dvrr_stack+2131);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+12734, dvrr_stack+6118, dvrr_stack+371, dvrr_stack+939, dvrr_stack+2146, dvrr_stack+237);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+12860, dvrr_stack+6181, dvrr_stack+12734, dvrr_stack+960, dvrr_stack+6118, dvrr_stack+11954);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+13070, dvrr_stack+6307, dvrr_stack+12860, dvrr_stack+1023, dvrr_stack+6181, dvrr_stack+12044);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+13385, dvrr_stack+6517, dvrr_stack+13070, dvrr_stack+1275, dvrr_stack+6307, dvrr_stack+12194);
 tmp = dvrr_stack + 13385;
 target_ptr = Libderiv->dvrr_classes[5][5];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+13826,dvrr_stack+13385,dvrr_stack+12419,21);


 /* compute (1 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+939, dvrr_stack+2019, dvrr_stack+845, NULL, NULL, dvrr_stack+337);

 /* compute (2 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14771, dvrr_stack+761, dvrr_stack+939, dvrr_stack+2167, dvrr_stack+2019, dvrr_stack+371);

 /* compute (3 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14939, dvrr_stack+7507, dvrr_stack+14771, dvrr_stack+2195, dvrr_stack+761, dvrr_stack+12734);

 /* compute (4 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+15219, dvrr_stack+7675, dvrr_stack+14939, dvrr_stack+2279, dvrr_stack+7507, dvrr_stack+12860);

 /* compute (5 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+15639, dvrr_stack+7955, dvrr_stack+15219, dvrr_stack+2615, dvrr_stack+7675, dvrr_stack+13070);
 tmp = dvrr_stack + 15639;
 target_ptr = Libderiv->dvrr_classes[5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+16227,dvrr_stack+15639,dvrr_stack+13385,21);


 /* compute (0 0 | 1 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+2167, Data->F+11, Data->F+12, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3774, dvrr_stack+2167, Data->F+10, Data->F+11, NULL);

 /* compute (0 0 | 3 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+2170, dvrr_stack+1935, dvrr_stack+6, dvrr_stack+873, dvrr_stack+3774, NULL);

 /* compute (0 0 | 4 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2180, dvrr_stack+39, dvrr_stack+2170, dvrr_stack+876, dvrr_stack+1935, NULL);

 /* compute (0 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4332, dvrr_stack+207, dvrr_stack+2180, dvrr_stack+882, dvrr_stack+39, NULL);

 /* compute (0 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2195, dvrr_stack+1941, dvrr_stack+4332, dvrr_stack+892, dvrr_stack+207, NULL);

 /* compute (0 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2223, dvrr_stack+845, dvrr_stack+2195, dvrr_stack+337, dvrr_stack+1941, NULL);

 /* compute (1 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2259, dvrr_stack+1962, dvrr_stack+2223, NULL, NULL, dvrr_stack+845);

 /* compute (2 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+17550, dvrr_stack+9320, dvrr_stack+2259, dvrr_stack+3777, dvrr_stack+1962, dvrr_stack+939);

 /* compute (3 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+17766, dvrr_stack+9428, dvrr_stack+17550, dvrr_stack+3813, dvrr_stack+9320, dvrr_stack+14771);

 /* compute (4 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18126, dvrr_stack+9644, dvrr_stack+17766, dvrr_stack+3921, dvrr_stack+9428, dvrr_stack+14939);

 /* compute (5 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18666, dvrr_stack+10004, dvrr_stack+18126, dvrr_stack+4353, dvrr_stack+9644, dvrr_stack+15219);

 /* compute (5 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+19422,dvrr_stack+18666,dvrr_stack+15639,21);


 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3777, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+3780, dvrr_stack+288, dvrr_stack+3777, Data->F+4, Data->F+5, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+873, dvrr_stack+68, dvrr_stack+49, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+3786, dvrr_stack+4137, dvrr_stack+873, dvrr_stack+3, dvrr_stack+68, dvrr_stack+3777);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+3804, dvrr_stack+4146, dvrr_stack+3786, dvrr_stack+228, dvrr_stack+4137, dvrr_stack+3780);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+3834, dvrr_stack+282, dvrr_stack+908, NULL, NULL, dvrr_stack+49);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+3852, dvrr_stack+4164, dvrr_stack+3834, dvrr_stack+71, dvrr_stack+282, dvrr_stack+873);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+3888, dvrr_stack+2047, dvrr_stack+3852, dvrr_stack+2083, dvrr_stack+4164, dvrr_stack+3786);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+3948, dvrr_stack+4182, dvrr_stack+3888, dvrr_stack+3633, dvrr_stack+2047, dvrr_stack+3804);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+71, dvrr_stack+914, dvrr_stack+301, NULL, NULL, dvrr_stack+908);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4038, dvrr_stack+4242, dvrr_stack+71, dvrr_stack+361, dvrr_stack+914, dvrr_stack+3834);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+1023, dvrr_stack+4272, dvrr_stack+4038, dvrr_stack+2101, dvrr_stack+4242, dvrr_stack+3852);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+21186, dvrr_stack+3525, dvrr_stack+1023, dvrr_stack+3669, dvrr_stack+4272, dvrr_stack+3888);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+21336, dvrr_stack+11804, dvrr_stack+21186, dvrr_stack+5553, dvrr_stack+3525, dvrr_stack+3948);

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3669, dvrr_stack+2131, dvrr_stack+892, NULL, NULL, dvrr_stack+301);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21546, dvrr_stack+237, dvrr_stack+3669, dvrr_stack+924, dvrr_stack+2131, dvrr_stack+71);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21636, dvrr_stack+11954, dvrr_stack+21546, dvrr_stack+3729, dvrr_stack+237, dvrr_stack+4038);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21786, dvrr_stack+12044, dvrr_stack+21636, dvrr_stack+5653, dvrr_stack+11954, dvrr_stack+1023);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+22011, dvrr_stack+12194, dvrr_stack+21786, dvrr_stack+5743, dvrr_stack+12044, dvrr_stack+21186);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+22326, dvrr_stack+12419, dvrr_stack+22011, dvrr_stack+5893, dvrr_stack+12194, dvrr_stack+21336);
 tmp = dvrr_stack + 22326;
 target_ptr = Libderiv->dvrr_classes[6][4];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+5653, dvrr_stack+337, dvrr_stack+1941, NULL, NULL, dvrr_stack+892);

 /* compute (2 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+5716, dvrr_stack+371, dvrr_stack+5653, dvrr_stack+2146, dvrr_stack+337, dvrr_stack+3669);

 /* compute (3 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+22746, dvrr_stack+12734, dvrr_stack+5716, dvrr_stack+6118, dvrr_stack+371, dvrr_stack+21546);

 /* compute (4 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+22956, dvrr_stack+12860, dvrr_stack+22746, dvrr_stack+6181, dvrr_stack+12734, dvrr_stack+21636);

 /* compute (5 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+23271, dvrr_stack+13070, dvrr_stack+22956, dvrr_stack+6307, dvrr_stack+12860, dvrr_stack+21786);

 /* compute (6 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+23712, dvrr_stack+13385, dvrr_stack+23271, dvrr_stack+6517, dvrr_stack+13070, dvrr_stack+22011);
 tmp = dvrr_stack + 23712;
 target_ptr = Libderiv->dvrr_classes[6][5];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+24300,dvrr_stack+23712,dvrr_stack+22326,28);


 /* compute (1 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6118, dvrr_stack+845, dvrr_stack+2195, NULL, NULL, dvrr_stack+1941);

 /* compute (2 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+6202, dvrr_stack+939, dvrr_stack+6118, dvrr_stack+2019, dvrr_stack+845, dvrr_stack+5653);

 /* compute (3 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+25560, dvrr_stack+14771, dvrr_stack+6202, dvrr_stack+761, dvrr_stack+939, dvrr_stack+5716);

 /* compute (4 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+25840, dvrr_stack+14939, dvrr_stack+25560, dvrr_stack+7507, dvrr_stack+14771, dvrr_stack+22746);

 /* compute (5 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+26260, dvrr_stack+15219, dvrr_stack+25840, dvrr_stack+7675, dvrr_stack+14939, dvrr_stack+22956);

 /* compute (6 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+26848, dvrr_stack+15639, dvrr_stack+26260, dvrr_stack+7955, dvrr_stack+15219, dvrr_stack+23271);

 /* compute (6 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+27632,dvrr_stack+26848,dvrr_stack+23712,28);


 /* compute (0 0 | 1 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+12, Data->F+13, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+7507, dvrr_stack+2167, dvrr_stack+3, Data->F+11, Data->F+12, NULL);

 /* compute (0 0 | 3 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+361, dvrr_stack+6, dvrr_stack+7507, dvrr_stack+3774, dvrr_stack+2167, NULL);

 /* compute (0 0 | 4 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+924, dvrr_stack+2170, dvrr_stack+361, dvrr_stack+1935, dvrr_stack+6, NULL);

 /* compute (0 0 | 5 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+7507, dvrr_stack+2180, dvrr_stack+924, dvrr_stack+39, dvrr_stack+2170, NULL);

 /* compute (0 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+7528, dvrr_stack+4332, dvrr_stack+7507, dvrr_stack+207, dvrr_stack+2180, NULL);

 /* compute (0 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7556, dvrr_stack+2195, dvrr_stack+7528, dvrr_stack+1941, dvrr_stack+4332, NULL);

 /* compute (1 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7592, dvrr_stack+2223, dvrr_stack+7556, NULL, NULL, dvrr_stack+2195);

 /* compute (2 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7700, dvrr_stack+2259, dvrr_stack+7592, dvrr_stack+1962, dvrr_stack+2223, dvrr_stack+6118);

 /* compute (3 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+29396, dvrr_stack+17550, dvrr_stack+7700, dvrr_stack+9320, dvrr_stack+2259, dvrr_stack+6202);

 /* compute (4 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+29756, dvrr_stack+17766, dvrr_stack+29396, dvrr_stack+9428, dvrr_stack+17550, dvrr_stack+25560);

 /* compute (5 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+30296, dvrr_stack+18126, dvrr_stack+29756, dvrr_stack+9644, dvrr_stack+17766, dvrr_stack+25840);

 /* compute (6 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+31052, dvrr_stack+18666, dvrr_stack+30296, dvrr_stack+10004, dvrr_stack+18126, dvrr_stack+26260);

 /* compute (6 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+32060,dvrr_stack+31052,dvrr_stack+26848,28);


 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3, dvrr_stack+12, dvrr_stack+0, NULL, NULL, Data->F+3);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+17550, dvrr_stack+52, dvrr_stack+15, NULL, NULL, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+17568, dvrr_stack+17550, dvrr_stack+21, dvrr_stack+52, dvrr_stack+15, dvrr_stack+3);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+17604, dvrr_stack+291, dvrr_stack+58, NULL, NULL, dvrr_stack+52);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+17634, dvrr_stack+17604, dvrr_stack+117, dvrr_stack+291, dvrr_stack+58, dvrr_stack+17550);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+17694, dvrr_stack+17634, dvrr_stack+147, dvrr_stack+17604, dvrr_stack+117, dvrr_stack+17568);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+17604, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+17607, dvrr_stack+3, dvrr_stack+228, dvrr_stack+12, dvrr_stack+0, dvrr_stack+17604);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+17794, dvrr_stack+17568, dvrr_stack+3633, dvrr_stack+17550, dvrr_stack+21, dvrr_stack+17607);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+17854, dvrr_stack+17694, dvrr_stack+5553, dvrr_stack+17634, dvrr_stack+147, dvrr_stack+17794);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+1935, dvrr_stack+17604, dvrr_stack+288, Data->F+3, Data->F+4, NULL);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+12, dvrr_stack+17607, dvrr_stack+4146, dvrr_stack+3, dvrr_stack+228, dvrr_stack+1935);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+18004, dvrr_stack+17794, dvrr_stack+4182, dvrr_stack+17568, dvrr_stack+3633, dvrr_stack+12);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+18094, dvrr_stack+17854, dvrr_stack+11804, dvrr_stack+17694, dvrr_stack+5553, dvrr_stack+18004);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+291, dvrr_stack+1935, dvrr_stack+3780, dvrr_stack+17604, dvrr_stack+288, NULL);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+5553, dvrr_stack+12, dvrr_stack+3804, dvrr_stack+17607, dvrr_stack+4146, dvrr_stack+291);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+17550, dvrr_stack+18004, dvrr_stack+3948, dvrr_stack+17794, dvrr_stack+4182, dvrr_stack+5553);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+18304, dvrr_stack+18094, dvrr_stack+21336, dvrr_stack+17854, dvrr_stack+11804, dvrr_stack+17550);

 /* compute (1 0 | 0 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+17550, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+1935, dvrr_stack+3777, dvrr_stack+17550, Data->F+5, Data->F+6, NULL);

 /* compute (3 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+291, dvrr_stack+3780, dvrr_stack+1935, dvrr_stack+288, dvrr_stack+3777, NULL);

 /* compute (1 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+228, dvrr_stack+49, dvrr_stack+358, NULL, NULL, Data->F+7);

 /* compute (2 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+17553, dvrr_stack+873, dvrr_stack+228, dvrr_stack+68, dvrr_stack+49, dvrr_stack+17550);

 /* compute (3 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+17571, dvrr_stack+3786, dvrr_stack+17553, dvrr_stack+4137, dvrr_stack+873, dvrr_stack+1935);

 /* compute (4 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+17601, dvrr_stack+3804, dvrr_stack+17571, dvrr_stack+4146, dvrr_stack+3786, dvrr_stack+291);

 /* compute (1 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+17646, dvrr_stack+908, dvrr_stack+222, NULL, NULL, dvrr_stack+358);

 /* compute (2 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+5553, dvrr_stack+3834, dvrr_stack+17646, dvrr_stack+282, dvrr_stack+908, dvrr_stack+228);

 /* compute (3 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+17794, dvrr_stack+3852, dvrr_stack+5553, dvrr_stack+4164, dvrr_stack+3834, dvrr_stack+17553);

 /* compute (4 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+18004, dvrr_stack+3888, dvrr_stack+17794, dvrr_stack+2047, dvrr_stack+3852, dvrr_stack+17571);

 /* compute (5 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+29396, dvrr_stack+3948, dvrr_stack+18004, dvrr_stack+4182, dvrr_stack+3888, dvrr_stack+17601);

 /* compute (1 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+17664, dvrr_stack+301, dvrr_stack+882, NULL, NULL, dvrr_stack+222);

 /* compute (2 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+17550, dvrr_stack+71, dvrr_stack+17664, dvrr_stack+914, dvrr_stack+301, dvrr_stack+17646);

 /* compute (3 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+29522, dvrr_stack+4038, dvrr_stack+17550, dvrr_stack+4242, dvrr_stack+71, dvrr_stack+5553);

 /* compute (4 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+1023, dvrr_stack+29522, dvrr_stack+4272, dvrr_stack+4038, dvrr_stack+17794);

 /* compute (5 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+29622, dvrr_stack+21186, dvrr_stack+0, dvrr_stack+3525, dvrr_stack+1023, dvrr_stack+18004);

 /* compute (6 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+29832, dvrr_stack+21336, dvrr_stack+29622, dvrr_stack+11804, dvrr_stack+21186, dvrr_stack+29396);

 /* compute (1 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+29396, dvrr_stack+892, dvrr_stack+207, NULL, NULL, dvrr_stack+882);

 /* compute (2 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+18004, dvrr_stack+3669, dvrr_stack+29396, dvrr_stack+2131, dvrr_stack+892, dvrr_stack+17664);

 /* compute (3 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11804, dvrr_stack+21546, dvrr_stack+18004, dvrr_stack+237, dvrr_stack+3669, dvrr_stack+17550);

 /* compute (4 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21186, dvrr_stack+21636, dvrr_stack+11804, dvrr_stack+11954, dvrr_stack+21546, dvrr_stack+29522);

 /* compute (5 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3525, dvrr_stack+21786, dvrr_stack+21186, dvrr_stack+12044, dvrr_stack+21636, dvrr_stack+0);

 /* compute (6 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+3840, dvrr_stack+22011, dvrr_stack+3525, dvrr_stack+12194, dvrr_stack+21786, dvrr_stack+29622);

 /* compute (7 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+21411, dvrr_stack+22326, dvrr_stack+3840, dvrr_stack+12419, dvrr_stack+22011, dvrr_stack+29832);

 /* compute (1 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+1941, dvrr_stack+4332, NULL, NULL, dvrr_stack+207);

 /* compute (2 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1023, dvrr_stack+5653, dvrr_stack+0, dvrr_stack+337, dvrr_stack+1941, dvrr_stack+29396);

 /* compute (3 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+29396, dvrr_stack+5716, dvrr_stack+1023, dvrr_stack+371, dvrr_stack+5653, dvrr_stack+18004);

 /* compute (4 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+29606, dvrr_stack+22746, dvrr_stack+29396, dvrr_stack+12734, dvrr_stack+5716, dvrr_stack+11804);

 /* compute (5 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+11804, dvrr_stack+22956, dvrr_stack+29606, dvrr_stack+12860, dvrr_stack+22746, dvrr_stack+21186);

 /* compute (6 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+29921, dvrr_stack+23271, dvrr_stack+11804, dvrr_stack+13070, dvrr_stack+22956, dvrr_stack+3525);

 /* compute (7 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+34412, dvrr_stack+23712, dvrr_stack+29921, dvrr_stack+13385, dvrr_stack+23271, dvrr_stack+3840);

 /* compute (1 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+761, dvrr_stack+2195, dvrr_stack+7528, NULL, NULL, dvrr_stack+4332);

 /* compute (2 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3525, dvrr_stack+6118, dvrr_stack+761, dvrr_stack+845, dvrr_stack+2195, dvrr_stack+0);

 /* compute (3 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+6202, dvrr_stack+3525, dvrr_stack+939, dvrr_stack+6118, dvrr_stack+1023);

 /* compute (4 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3525, dvrr_stack+25560, dvrr_stack+0, dvrr_stack+14771, dvrr_stack+6202, dvrr_stack+29396);

 /* compute (5 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+22746, dvrr_stack+25840, dvrr_stack+3525, dvrr_stack+14939, dvrr_stack+25560, dvrr_stack+29606);

 /* compute (6 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3525, dvrr_stack+26260, dvrr_stack+22746, dvrr_stack+15219, dvrr_stack+25840, dvrr_stack+11804);

 /* compute (7 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 7;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+35168, dvrr_stack+26848, dvrr_stack+3525, dvrr_stack+15639, dvrr_stack+26260, dvrr_stack+29921);

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+3525, dvrr_stack+1485, NULL);
 tmp = dvrr_stack + 3525;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+3675, dvrr_stack+2895, NULL);
 tmp = dvrr_stack + 3675;
 target_ptr = Libderiv->deriv_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,280,dvrr_stack+3885, dvrr_stack+4713, NULL);
 tmp = dvrr_stack + 3885;
 target_ptr = Libderiv->deriv_classes[3][6][11];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+21186, dvrr_stack+6832, NULL);
 tmp = dvrr_stack + 21186;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+11804, dvrr_stack+8375, NULL);
 tmp = dvrr_stack + 11804;
 target_ptr = Libderiv->deriv_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+22746, dvrr_stack+10544, NULL);
 tmp = dvrr_stack + 22746;
 target_ptr = Libderiv->deriv_classes[4][6][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+23166, dvrr_stack+13826, NULL);
 tmp = dvrr_stack + 23166;
 target_ptr = Libderiv->deriv_classes[5][4][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,441,dvrr_stack+25560, dvrr_stack+16227, NULL);
 tmp = dvrr_stack + 25560;
 target_ptr = Libderiv->deriv_classes[5][5][11];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+26001, dvrr_stack+19422, NULL);
 tmp = dvrr_stack + 26001;
 target_ptr = Libderiv->deriv_classes[5][6][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+29396, dvrr_stack+24300, NULL);
 tmp = dvrr_stack + 29396;
 target_ptr = Libderiv->deriv_classes[6][4][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+29816, dvrr_stack+27632, NULL);
 tmp = dvrr_stack + 29816;
 target_ptr = Libderiv->deriv_classes[6][5][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,784,dvrr_stack+14771, dvrr_stack+32060, NULL);
 tmp = dvrr_stack + 14771;
 target_ptr = Libderiv->deriv_classes[6][6][11];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+4165, dvrr_stack+1485, NULL);
 tmp = dvrr_stack + 4165;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+12119, dvrr_stack+2895, NULL);
 tmp = dvrr_stack + 12119;
 target_ptr = Libderiv->deriv_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+30404, dvrr_stack+4713, NULL);
 tmp = dvrr_stack + 30404;
 target_ptr = Libderiv->deriv_classes[3][6][10];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+23481, dvrr_stack+6832, NULL);
 tmp = dvrr_stack + 23481;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+30684, dvrr_stack+8375, NULL);
 tmp = dvrr_stack + 30684;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+0, dvrr_stack+10544, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[4][6][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+6118, dvrr_stack+13826, NULL);
 tmp = dvrr_stack + 6118;
 target_ptr = Libderiv->deriv_classes[5][4][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+12734, dvrr_stack+16227, NULL);
 tmp = dvrr_stack + 12734;
 target_ptr = Libderiv->deriv_classes[5][5][10];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+9320, dvrr_stack+19422, NULL);
 tmp = dvrr_stack + 9320;
 target_ptr = Libderiv->deriv_classes[5][6][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+1935, dvrr_stack+24300, NULL);
 tmp = dvrr_stack + 1935;
 target_ptr = Libderiv->deriv_classes[6][4][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+36176, dvrr_stack+27632, NULL);
 tmp = dvrr_stack + 36176;
 target_ptr = Libderiv->deriv_classes[6][5][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,784,dvrr_stack+36764, dvrr_stack+32060, NULL);
 tmp = dvrr_stack + 36764;
 target_ptr = Libderiv->deriv_classes[6][6][10];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+26589, dvrr_stack+1485, NULL);
 tmp = dvrr_stack + 26589;
 target_ptr = Libderiv->deriv_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+13175, dvrr_stack+2895, NULL);
 tmp = dvrr_stack + 13175;
 target_ptr = Libderiv->deriv_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+2895, dvrr_stack+4713, NULL);
 tmp = dvrr_stack + 2895;
 target_ptr = Libderiv->deriv_classes[3][6][9];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+4713, dvrr_stack+6832, NULL);
 tmp = dvrr_stack + 4713;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+6832, dvrr_stack+8375, NULL);
 tmp = dvrr_stack + 6832;
 target_ptr = Libderiv->deriv_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+8375, dvrr_stack+10544, NULL);
 tmp = dvrr_stack + 8375;
 target_ptr = Libderiv->deriv_classes[4][6][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+10544, dvrr_stack+13826, NULL);
 tmp = dvrr_stack + 10544;
 target_ptr = Libderiv->deriv_classes[5][4][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+13826, dvrr_stack+16227, NULL);
 tmp = dvrr_stack + 13826;
 target_ptr = Libderiv->deriv_classes[5][5][9];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+16227, dvrr_stack+19422, NULL);
 tmp = dvrr_stack + 16227;
 target_ptr = Libderiv->deriv_classes[5][6][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+19422, dvrr_stack+24300, NULL);
 tmp = dvrr_stack + 19422;
 target_ptr = Libderiv->deriv_classes[6][4][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+24300, dvrr_stack+27632, NULL);
 tmp = dvrr_stack + 24300;
 target_ptr = Libderiv->deriv_classes[6][5][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,784,dvrr_stack+27632, dvrr_stack+32060, NULL);
 tmp = dvrr_stack + 27632;
 target_ptr = Libderiv->deriv_classes[6][6][9];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+32060, dvrr_stack+1275, dvrr_stack+17694);
 tmp = dvrr_stack + 32060;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+32210, dvrr_stack+2615, dvrr_stack+611);
 tmp = dvrr_stack + 32210;
 target_ptr = Libderiv->deriv_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+32420, dvrr_stack+4353, dvrr_stack+1275);
 tmp = dvrr_stack + 32420;
 target_ptr = Libderiv->deriv_classes[3][6][8];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+32700, dvrr_stack+6517, dvrr_stack+17854);
 tmp = dvrr_stack + 32700;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+32925, dvrr_stack+7955, dvrr_stack+5893);
 tmp = dvrr_stack + 32925;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+33240, dvrr_stack+10004, dvrr_stack+6517);
 tmp = dvrr_stack + 33240;
 target_ptr = Libderiv->deriv_classes[4][6][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+33660, dvrr_stack+13385, dvrr_stack+18094);
 tmp = dvrr_stack + 33660;
 target_ptr = Libderiv->deriv_classes[5][4][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+28416, dvrr_stack+15639, dvrr_stack+12419);
 tmp = dvrr_stack + 28416;
 target_ptr = Libderiv->deriv_classes[5][5][8];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+24888, dvrr_stack+18666, dvrr_stack+13385);
 tmp = dvrr_stack + 24888;
 target_ptr = Libderiv->deriv_classes[5][6][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,28,1,dvrr_stack+33975, dvrr_stack+23712, dvrr_stack+18304);
 tmp = dvrr_stack + 33975;
 target_ptr = Libderiv->deriv_classes[6][4][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,28,1,dvrr_stack+19842, dvrr_stack+26848, dvrr_stack+22326);
 tmp = dvrr_stack + 19842;
 target_ptr = Libderiv->deriv_classes[6][5][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,28,1,dvrr_stack+16815, dvrr_stack+31052, dvrr_stack+23712);
 tmp = dvrr_stack + 16815;
 target_ptr = Libderiv->deriv_classes[6][6][8];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+28857, dvrr_stack+1275, dvrr_stack+17694);
 tmp = dvrr_stack + 28857;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+29007, dvrr_stack+2615, dvrr_stack+611);
 tmp = dvrr_stack + 29007;
 target_ptr = Libderiv->deriv_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+20430, dvrr_stack+4353, dvrr_stack+1275);
 tmp = dvrr_stack + 20430;
 target_ptr = Libderiv->deriv_classes[3][6][7];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+20710, dvrr_stack+6517, dvrr_stack+17854);
 tmp = dvrr_stack + 20710;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+14267, dvrr_stack+7955, dvrr_stack+5893);
 tmp = dvrr_stack + 14267;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+10859, dvrr_stack+10004, dvrr_stack+6517);
 tmp = dvrr_stack + 10859;
 target_ptr = Libderiv->deriv_classes[4][6][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+11279, dvrr_stack+13385, dvrr_stack+18094);
 tmp = dvrr_stack + 11279;
 target_ptr = Libderiv->deriv_classes[5][4][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+8795, dvrr_stack+15639, dvrr_stack+12419);
 tmp = dvrr_stack + 8795;
 target_ptr = Libderiv->deriv_classes[5][5][7];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+7147, dvrr_stack+18666, dvrr_stack+13385);
 tmp = dvrr_stack + 7147;
 target_ptr = Libderiv->deriv_classes[5][6][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,28,1,dvrr_stack+4938, dvrr_stack+23712, dvrr_stack+18304);
 tmp = dvrr_stack + 4938;
 target_ptr = Libderiv->deriv_classes[6][4][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,28,1,dvrr_stack+37548, dvrr_stack+26848, dvrr_stack+22326);
 tmp = dvrr_stack + 37548;
 target_ptr = Libderiv->deriv_classes[6][5][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,28,1,dvrr_stack+38136, dvrr_stack+31052, dvrr_stack+23712);
 tmp = dvrr_stack + 38136;
 target_ptr = Libderiv->deriv_classes[6][6][7];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+29217, dvrr_stack+1275, dvrr_stack+17694);
 tmp = dvrr_stack + 29217;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+11594, dvrr_stack+2615, dvrr_stack+611);
 tmp = dvrr_stack + 11594;
 target_ptr = Libderiv->deriv_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+5358, dvrr_stack+4353, dvrr_stack+1275);
 tmp = dvrr_stack + 5358;
 target_ptr = Libderiv->deriv_classes[3][6][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+20935, dvrr_stack+6517, dvrr_stack+17854);
 tmp = dvrr_stack + 20935;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+17599, dvrr_stack+7955, dvrr_stack+5893);
 tmp = dvrr_stack + 17599;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+1485, dvrr_stack+10004, dvrr_stack+6517);
 tmp = dvrr_stack + 1485;
 target_ptr = Libderiv->deriv_classes[4][6][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+3175, dvrr_stack+13385, dvrr_stack+18094);
 tmp = dvrr_stack + 3175;
 target_ptr = Libderiv->deriv_classes[5][4][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+9908, dvrr_stack+15639, dvrr_stack+12419);
 tmp = dvrr_stack + 9908;
 target_ptr = Libderiv->deriv_classes[5][5][6];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+38920, dvrr_stack+18666, dvrr_stack+13385);
 tmp = dvrr_stack + 38920;
 target_ptr = Libderiv->deriv_classes[5][6][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,28,1,dvrr_stack+18584, dvrr_stack+23712, dvrr_stack+18304);
 tmp = dvrr_stack + 18584;
 target_ptr = Libderiv->deriv_classes[6][4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,28,1,dvrr_stack+17914, dvrr_stack+26848, dvrr_stack+22326);
 tmp = dvrr_stack + 17914;
 target_ptr = Libderiv->deriv_classes[6][5][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,28,1,dvrr_stack+39508, dvrr_stack+31052, dvrr_stack+23712);
 tmp = dvrr_stack + 39508;
 target_ptr = Libderiv->deriv_classes[6][6][6];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+14582, dvrr_stack+5893, dvrr_stack+521);
 tmp = dvrr_stack + 14582;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+7735, dvrr_stack+6517, dvrr_stack+1149);
 tmp = dvrr_stack + 7735;
 target_ptr = Libderiv->deriv_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+4315, dvrr_stack+7955, dvrr_stack+2447);
 tmp = dvrr_stack + 4315;
 target_ptr = Libderiv->deriv_classes[3][6][2];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+5638, dvrr_stack+12419, dvrr_stack+611);
 tmp = dvrr_stack + 5638;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+30999, dvrr_stack+13385, dvrr_stack+1275);
 tmp = dvrr_stack + 30999;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+31314, dvrr_stack+15639, dvrr_stack+2615);
 tmp = dvrr_stack + 31314;
 target_ptr = Libderiv->deriv_classes[4][6][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+31734, dvrr_stack+22326, dvrr_stack+5893);
 tmp = dvrr_stack + 31734;
 target_ptr = Libderiv->deriv_classes[5][4][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+40292, dvrr_stack+23712, dvrr_stack+6517);
 tmp = dvrr_stack + 40292;
 target_ptr = Libderiv->deriv_classes[5][5][2];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+40733, dvrr_stack+26848, dvrr_stack+7955);
 tmp = dvrr_stack + 40733;
 target_ptr = Libderiv->deriv_classes[5][6][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,15,dvrr_stack+41321, dvrr_stack+21411, dvrr_stack+12419);
 tmp = dvrr_stack + 41321;
 target_ptr = Libderiv->deriv_classes[6][4][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,21,dvrr_stack+41741, dvrr_stack+34412, dvrr_stack+13385);
 tmp = dvrr_stack + 41741;
 target_ptr = Libderiv->deriv_classes[6][5][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_i(Data,28,dvrr_stack+42329, dvrr_stack+35168, dvrr_stack+15639);
 tmp = dvrr_stack + 42329;
 target_ptr = Libderiv->deriv_classes[6][6][2];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+761, dvrr_stack+5893, dvrr_stack+521);
 tmp = dvrr_stack + 761;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+911, dvrr_stack+6517, dvrr_stack+1149);
 tmp = dvrr_stack + 911;
 target_ptr = Libderiv->deriv_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+21951, dvrr_stack+7955, dvrr_stack+2447);
 tmp = dvrr_stack + 21951;
 target_ptr = Libderiv->deriv_classes[3][6][1];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+19004, dvrr_stack+12419, dvrr_stack+611);
 tmp = dvrr_stack + 19004;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+43113, dvrr_stack+13385, dvrr_stack+1275);
 tmp = dvrr_stack + 43113;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+43428, dvrr_stack+15639, dvrr_stack+2615);
 tmp = dvrr_stack + 43428;
 target_ptr = Libderiv->deriv_classes[4][6][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+43848, dvrr_stack+22326, dvrr_stack+5893);
 tmp = dvrr_stack + 43848;
 target_ptr = Libderiv->deriv_classes[5][4][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+44163, dvrr_stack+23712, dvrr_stack+6517);
 tmp = dvrr_stack + 44163;
 target_ptr = Libderiv->deriv_classes[5][5][1];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+44604, dvrr_stack+26848, dvrr_stack+7955);
 tmp = dvrr_stack + 44604;
 target_ptr = Libderiv->deriv_classes[5][6][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,15,dvrr_stack+45192, dvrr_stack+21411, dvrr_stack+12419);
 tmp = dvrr_stack + 45192;
 target_ptr = Libderiv->deriv_classes[6][4][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,21,dvrr_stack+45612, dvrr_stack+34412, dvrr_stack+13385);
 tmp = dvrr_stack + 45612;
 target_ptr = Libderiv->deriv_classes[6][5][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_i(Data,28,dvrr_stack+46200, dvrr_stack+35168, dvrr_stack+15639);
 tmp = dvrr_stack + 46200;
 target_ptr = Libderiv->deriv_classes[6][6][1];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+19229, dvrr_stack+5893, dvrr_stack+521);
 tmp = dvrr_stack + 19229;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+46984, dvrr_stack+6517, dvrr_stack+1149);
 tmp = dvrr_stack + 46984;
 target_ptr = Libderiv->deriv_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+47194, dvrr_stack+7955, dvrr_stack+2447);
 tmp = dvrr_stack + 47194;
 target_ptr = Libderiv->deriv_classes[3][6][0];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+2355, dvrr_stack+12419, dvrr_stack+611);
 tmp = dvrr_stack + 2355;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+420, dvrr_stack+13385, dvrr_stack+1275);
 tmp = dvrr_stack + 420;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+47474, dvrr_stack+15639, dvrr_stack+2615);
 tmp = dvrr_stack + 47474;
 target_ptr = Libderiv->deriv_classes[4][6][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+2580, dvrr_stack+22326, dvrr_stack+5893);
 tmp = dvrr_stack + 2580;
 target_ptr = Libderiv->deriv_classes[5][4][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+22231, dvrr_stack+23712, dvrr_stack+6517);
 tmp = dvrr_stack + 22231;
 target_ptr = Libderiv->deriv_classes[5][5][0];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+23706, dvrr_stack+26848, dvrr_stack+7955);
 tmp = dvrr_stack + 23706;
 target_ptr = Libderiv->deriv_classes[5][6][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,15,dvrr_stack+7945, dvrr_stack+21411, dvrr_stack+12419);
 tmp = dvrr_stack + 7945;
 target_ptr = Libderiv->deriv_classes[6][4][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,21,dvrr_stack+26739, dvrr_stack+34412, dvrr_stack+13385);
 tmp = dvrr_stack + 26739;
 target_ptr = Libderiv->deriv_classes[6][5][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (6 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_i(Data,28,dvrr_stack+47894, dvrr_stack+35168, dvrr_stack+15639);
 tmp = dvrr_stack + 47894;
 target_ptr = Libderiv->deriv_classes[6][6][0];
 for(i=0;i<784;i++)
   target_ptr[i] += tmp[i];


}

