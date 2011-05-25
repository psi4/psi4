#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (dd|gg) integrals */

void d1vrr_order_ddgg(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+0, dvrr_stack+12, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+6, dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+34, dvrr_stack+31, dvrr_stack+3, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+34, dvrr_stack+6, dvrr_stack+31, dvrr_stack+3, NULL);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+50, dvrr_stack+40, dvrr_stack+21, NULL, NULL, dvrr_stack+6);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+80, dvrr_stack+40, dvrr_stack+21, dvrr_stack+34, dvrr_stack+6, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+95, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+98, dvrr_stack+95, dvrr_stack+31, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+104, dvrr_stack+98, dvrr_stack+34, dvrr_stack+95, dvrr_stack+31, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+114, dvrr_stack+104, dvrr_stack+40, dvrr_stack+98, dvrr_stack+34, NULL);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+95, dvrr_stack+12, dvrr_stack+31, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+129, dvrr_stack+15, dvrr_stack+95, dvrr_stack+0, dvrr_stack+12, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+139, dvrr_stack+21, dvrr_stack+129, dvrr_stack+6, dvrr_stack+15, NULL);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+154, dvrr_stack+80, dvrr_stack+139, NULL, NULL, dvrr_stack+21);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+199, dvrr_stack+114, dvrr_stack+80, NULL, NULL, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+244, dvrr_stack+199, dvrr_stack+154, dvrr_stack+114, dvrr_stack+80, dvrr_stack+50);
 tmp = dvrr_stack + 244;
 target_ptr = Libderiv->dvrr_classes[2][4];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+334, dvrr_stack+80, dvrr_stack+139, dvrr_stack+40, dvrr_stack+21, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+355, dvrr_stack+114, dvrr_stack+80, dvrr_stack+104, dvrr_stack+40, NULL);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+101, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+376, dvrr_stack+31, dvrr_stack+101, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+382, dvrr_stack+95, dvrr_stack+376, dvrr_stack+12, dvrr_stack+31, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+392, dvrr_stack+129, dvrr_stack+382, dvrr_stack+15, dvrr_stack+95, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+407, dvrr_stack+139, dvrr_stack+392, dvrr_stack+21, dvrr_stack+129, NULL);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+428, dvrr_stack+334, dvrr_stack+407, NULL, NULL, dvrr_stack+139);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+491, dvrr_stack+355, dvrr_stack+334, NULL, NULL, dvrr_stack+80);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+554, dvrr_stack+491, dvrr_stack+428, dvrr_stack+355, dvrr_stack+334, dvrr_stack+154);
 tmp = dvrr_stack + 554;
 target_ptr = Libderiv->dvrr_classes[2][5];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+680,dvrr_stack+554,dvrr_stack+244,6);


 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+950, dvrr_stack+334, dvrr_stack+407, dvrr_stack+80, dvrr_stack+139, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+978, dvrr_stack+355, dvrr_stack+334, dvrr_stack+114, dvrr_stack+80, NULL);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+114, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+117, dvrr_stack+101, dvrr_stack+114, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1006, dvrr_stack+376, dvrr_stack+117, dvrr_stack+31, dvrr_stack+101, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1016, dvrr_stack+382, dvrr_stack+1006, dvrr_stack+95, dvrr_stack+376, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1031, dvrr_stack+392, dvrr_stack+1016, dvrr_stack+129, dvrr_stack+382, NULL);

 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1052, dvrr_stack+407, dvrr_stack+1031, dvrr_stack+139, dvrr_stack+392, NULL);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1080, dvrr_stack+950, dvrr_stack+1052, NULL, NULL, dvrr_stack+407);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1164, dvrr_stack+978, dvrr_stack+950, NULL, NULL, dvrr_stack+334);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1248, dvrr_stack+1164, dvrr_stack+1080, dvrr_stack+978, dvrr_stack+950, dvrr_stack+428);
 tmp = dvrr_stack + 1248;
 target_ptr = Libderiv->dvrr_classes[2][6];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+1416,dvrr_stack+1248,dvrr_stack+554,6);


 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1794, dvrr_stack+950, dvrr_stack+1052, dvrr_stack+334, dvrr_stack+407, NULL);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1830, dvrr_stack+978, dvrr_stack+950, dvrr_stack+355, dvrr_stack+334, NULL);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+355, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+123, dvrr_stack+114, dvrr_stack+355, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+358, dvrr_stack+117, dvrr_stack+123, dvrr_stack+101, dvrr_stack+114, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1866, dvrr_stack+1006, dvrr_stack+358, dvrr_stack+376, dvrr_stack+117, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1881, dvrr_stack+1016, dvrr_stack+1866, dvrr_stack+382, dvrr_stack+1006, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1902, dvrr_stack+1031, dvrr_stack+1881, dvrr_stack+392, dvrr_stack+1016, NULL);

 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1930, dvrr_stack+1052, dvrr_stack+1902, dvrr_stack+407, dvrr_stack+1031, NULL);

 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1966, dvrr_stack+1794, dvrr_stack+1930, NULL, NULL, dvrr_stack+1052);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2074, dvrr_stack+1830, dvrr_stack+1794, NULL, NULL, dvrr_stack+950);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+2182, dvrr_stack+2074, dvrr_stack+1966, dvrr_stack+1830, dvrr_stack+1794, dvrr_stack+1080);
 tmp = dvrr_stack + 2182;
 target_ptr = Libderiv->dvrr_classes[2][7];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+2398,dvrr_stack+2182,dvrr_stack+1248,6);


 /* compute (0 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+2902, dvrr_stack+1794, dvrr_stack+1930, dvrr_stack+950, dvrr_stack+1052, NULL);

 /* compute (0 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+2947, dvrr_stack+1830, dvrr_stack+1794, dvrr_stack+978, dvrr_stack+950, NULL);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+101, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+978, dvrr_stack+355, dvrr_stack+101, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+984, dvrr_stack+123, dvrr_stack+978, dvrr_stack+114, dvrr_stack+355, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2992, dvrr_stack+358, dvrr_stack+984, dvrr_stack+117, dvrr_stack+123, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3007, dvrr_stack+1866, dvrr_stack+2992, dvrr_stack+1006, dvrr_stack+358, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3028, dvrr_stack+1881, dvrr_stack+3007, dvrr_stack+1016, dvrr_stack+1866, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3056, dvrr_stack+1902, dvrr_stack+3028, dvrr_stack+1031, dvrr_stack+1881, NULL);

 /* compute (0 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3092, dvrr_stack+1930, dvrr_stack+3056, dvrr_stack+1052, dvrr_stack+1902, NULL);

 /* compute (1 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3137, dvrr_stack+2902, dvrr_stack+3092, NULL, NULL, dvrr_stack+1930);

 /* compute (1 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3272, dvrr_stack+2947, dvrr_stack+2902, NULL, NULL, dvrr_stack+1794);

 /* compute (2 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3407, dvrr_stack+3272, dvrr_stack+3137, dvrr_stack+2947, dvrr_stack+2902, dvrr_stack+1966);
 tmp = dvrr_stack + 3407;
 target_ptr = Libderiv->dvrr_classes[2][8];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+3677,dvrr_stack+3407,dvrr_stack+2182,6);


 /* compute (0 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+4325, dvrr_stack+2902, dvrr_stack+3092, dvrr_stack+1794, dvrr_stack+1930, NULL);

 /* compute (0 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+4380, dvrr_stack+2947, dvrr_stack+2902, dvrr_stack+1830, dvrr_stack+1794, NULL);

 /* compute (0 0 | 1 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+1830, Data->F+10, Data->F+11, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+1833, dvrr_stack+101, dvrr_stack+1830, Data->F+9, Data->F+10, NULL);

 /* compute (0 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1839, dvrr_stack+978, dvrr_stack+1833, dvrr_stack+355, dvrr_stack+101, NULL);

 /* compute (0 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1849, dvrr_stack+984, dvrr_stack+1839, dvrr_stack+123, dvrr_stack+978, NULL);

 /* compute (0 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2947, dvrr_stack+2992, dvrr_stack+1849, dvrr_stack+358, dvrr_stack+984, NULL);

 /* compute (0 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4435, dvrr_stack+3007, dvrr_stack+2947, dvrr_stack+1866, dvrr_stack+2992, NULL);

 /* compute (0 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4463, dvrr_stack+3028, dvrr_stack+4435, dvrr_stack+1881, dvrr_stack+3007, NULL);

 /* compute (0 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+4499, dvrr_stack+3056, dvrr_stack+4463, dvrr_stack+1902, dvrr_stack+3028, NULL);

 /* compute (0 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+4544, dvrr_stack+3092, dvrr_stack+4499, dvrr_stack+1930, dvrr_stack+3056, NULL);

 /* compute (1 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+4599, dvrr_stack+4325, dvrr_stack+4544, NULL, NULL, dvrr_stack+3092);

 /* compute (1 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+4764, dvrr_stack+4380, dvrr_stack+4325, NULL, NULL, dvrr_stack+2902);

 /* compute (2 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+4929, dvrr_stack+4764, dvrr_stack+4599, dvrr_stack+4380, dvrr_stack+4325, dvrr_stack+3137);

 /* compute (2 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+5259,dvrr_stack+4929,dvrr_stack+3407,6);


 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+4380, dvrr_stack+6, dvrr_stack+15, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+4398, dvrr_stack+21, dvrr_stack+129, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+6069, dvrr_stack+50, dvrr_stack+4398, dvrr_stack+40, dvrr_stack+21, dvrr_stack+4380);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6129, dvrr_stack+139, dvrr_stack+392, NULL, NULL, dvrr_stack+129);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6174, dvrr_stack+154, dvrr_stack+6129, dvrr_stack+80, dvrr_stack+139, dvrr_stack+4398);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+6264, dvrr_stack+244, dvrr_stack+6174, dvrr_stack+199, dvrr_stack+154, dvrr_stack+6069);
 tmp = dvrr_stack + 6264;
 target_ptr = Libderiv->dvrr_classes[3][4];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6414, dvrr_stack+407, dvrr_stack+1031, NULL, NULL, dvrr_stack+392);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6477, dvrr_stack+428, dvrr_stack+6414, dvrr_stack+334, dvrr_stack+407, dvrr_stack+6129);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+6603, dvrr_stack+554, dvrr_stack+6477, dvrr_stack+491, dvrr_stack+428, dvrr_stack+6174);
 tmp = dvrr_stack + 6603;
 target_ptr = Libderiv->dvrr_classes[3][5];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+6813,dvrr_stack+6603,dvrr_stack+6264,10);


 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+7263, dvrr_stack+1052, dvrr_stack+1902, NULL, NULL, dvrr_stack+1031);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+7347, dvrr_stack+1080, dvrr_stack+7263, dvrr_stack+950, dvrr_stack+1052, dvrr_stack+6414);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+7515, dvrr_stack+1248, dvrr_stack+7347, dvrr_stack+1164, dvrr_stack+1080, dvrr_stack+6477);
 tmp = dvrr_stack + 7515;
 target_ptr = Libderiv->dvrr_classes[3][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+7795,dvrr_stack+7515,dvrr_stack+6603,10);


 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+8425, dvrr_stack+1930, dvrr_stack+3056, NULL, NULL, dvrr_stack+1902);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+8533, dvrr_stack+1966, dvrr_stack+8425, dvrr_stack+1794, dvrr_stack+1930, dvrr_stack+7263);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+8749, dvrr_stack+2182, dvrr_stack+8533, dvrr_stack+2074, dvrr_stack+1966, dvrr_stack+7347);
 tmp = dvrr_stack + 8749;
 target_ptr = Libderiv->dvrr_classes[3][7];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+9109,dvrr_stack+8749,dvrr_stack+7515,10);


 /* compute (1 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+9949, dvrr_stack+3092, dvrr_stack+4499, NULL, NULL, dvrr_stack+3056);

 /* compute (2 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10084, dvrr_stack+3137, dvrr_stack+9949, dvrr_stack+2902, dvrr_stack+3092, dvrr_stack+8425);

 /* compute (3 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+10354, dvrr_stack+3407, dvrr_stack+10084, dvrr_stack+3272, dvrr_stack+3137, dvrr_stack+8533);
 tmp = dvrr_stack + 10354;
 target_ptr = Libderiv->dvrr_classes[3][8];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+10804,dvrr_stack+10354,dvrr_stack+8749,10);


 /* compute (0 0 | 1 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+2902, Data->F+11, Data->F+12, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+2905, dvrr_stack+1830, dvrr_stack+2902, Data->F+10, Data->F+11, NULL);

 /* compute (0 0 | 3 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+2911, dvrr_stack+1833, dvrr_stack+2905, dvrr_stack+101, dvrr_stack+1830, NULL);

 /* compute (0 0 | 4 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+80, dvrr_stack+1839, dvrr_stack+2911, dvrr_stack+978, dvrr_stack+1833, NULL);

 /* compute (0 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2921, dvrr_stack+1849, dvrr_stack+80, dvrr_stack+984, dvrr_stack+1839, NULL);

 /* compute (0 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1794, dvrr_stack+2947, dvrr_stack+2921, dvrr_stack+2992, dvrr_stack+1849, NULL);

 /* compute (0 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+950, dvrr_stack+4435, dvrr_stack+1794, dvrr_stack+3007, dvrr_stack+2947, NULL);

 /* compute (0 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+11884, dvrr_stack+4463, dvrr_stack+950, dvrr_stack+3028, dvrr_stack+4435, NULL);

 /* compute (0 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+11929, dvrr_stack+4499, dvrr_stack+11884, dvrr_stack+3056, dvrr_stack+4463, NULL);

 /* compute (1 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+11984, dvrr_stack+4544, dvrr_stack+11929, NULL, NULL, dvrr_stack+4499);

 /* compute (2 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+12149, dvrr_stack+4599, dvrr_stack+11984, dvrr_stack+4325, dvrr_stack+4544, dvrr_stack+9949);

 /* compute (3 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+12479, dvrr_stack+4929, dvrr_stack+12149, dvrr_stack+4764, dvrr_stack+4599, dvrr_stack+10084);

 /* compute (3 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+13029,dvrr_stack+12479,dvrr_stack+10354,10);


 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+4764, dvrr_stack+0, dvrr_stack+12, NULL, NULL, Data->F+4);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+4773, dvrr_stack+15, dvrr_stack+95, NULL, NULL, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+4791, dvrr_stack+4380, dvrr_stack+4773, dvrr_stack+6, dvrr_stack+15, dvrr_stack+4764);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+4827, dvrr_stack+129, dvrr_stack+382, NULL, NULL, dvrr_stack+95);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+4857, dvrr_stack+4398, dvrr_stack+4827, dvrr_stack+21, dvrr_stack+129, dvrr_stack+4773);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+14379, dvrr_stack+6069, dvrr_stack+4857, dvrr_stack+50, dvrr_stack+4398, dvrr_stack+4791);

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4325, dvrr_stack+392, dvrr_stack+1016, NULL, NULL, dvrr_stack+382);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14479, dvrr_stack+6129, dvrr_stack+4325, dvrr_stack+139, dvrr_stack+392, dvrr_stack+4827);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14569, dvrr_stack+6174, dvrr_stack+14479, dvrr_stack+154, dvrr_stack+6129, dvrr_stack+4857);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14719, dvrr_stack+6264, dvrr_stack+14569, dvrr_stack+244, dvrr_stack+6174, dvrr_stack+14379);
 tmp = dvrr_stack + 14719;
 target_ptr = Libderiv->dvrr_classes[4][4];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+14944, dvrr_stack+1031, dvrr_stack+1881, NULL, NULL, dvrr_stack+1016);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15007, dvrr_stack+6414, dvrr_stack+14944, dvrr_stack+407, dvrr_stack+1031, dvrr_stack+4325);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15133, dvrr_stack+6477, dvrr_stack+15007, dvrr_stack+428, dvrr_stack+6414, dvrr_stack+14479);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+15343, dvrr_stack+6603, dvrr_stack+15133, dvrr_stack+554, dvrr_stack+6477, dvrr_stack+14569);
 tmp = dvrr_stack + 15343;
 target_ptr = Libderiv->dvrr_classes[4][5];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+15658,dvrr_stack+15343,dvrr_stack+14719,15);


 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+407, dvrr_stack+1902, dvrr_stack+3028, NULL, NULL, dvrr_stack+1881);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16333, dvrr_stack+7263, dvrr_stack+407, dvrr_stack+1052, dvrr_stack+1902, dvrr_stack+14944);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16501, dvrr_stack+7347, dvrr_stack+16333, dvrr_stack+1080, dvrr_stack+7263, dvrr_stack+15007);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+16781, dvrr_stack+7515, dvrr_stack+16501, dvrr_stack+1248, dvrr_stack+7347, dvrr_stack+15133);
 tmp = dvrr_stack + 16781;
 target_ptr = Libderiv->dvrr_classes[4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+17201,dvrr_stack+16781,dvrr_stack+15343,15);


 /* compute (1 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+1052, dvrr_stack+3056, dvrr_stack+4463, NULL, NULL, dvrr_stack+3028);

 /* compute (2 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18146, dvrr_stack+8425, dvrr_stack+1052, dvrr_stack+1930, dvrr_stack+3056, dvrr_stack+407);

 /* compute (3 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18362, dvrr_stack+8533, dvrr_stack+18146, dvrr_stack+1966, dvrr_stack+8425, dvrr_stack+16333);

 /* compute (4 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18722, dvrr_stack+8749, dvrr_stack+18362, dvrr_stack+2182, dvrr_stack+8533, dvrr_stack+16501);
 tmp = dvrr_stack + 18722;
 target_ptr = Libderiv->dvrr_classes[4][7];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+19262,dvrr_stack+18722,dvrr_stack+16781,15);


 /* compute (1 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+1930, dvrr_stack+4499, dvrr_stack+11884, NULL, NULL, dvrr_stack+4463);

 /* compute (2 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+20522, dvrr_stack+9949, dvrr_stack+1930, dvrr_stack+3092, dvrr_stack+4499, dvrr_stack+1052);

 /* compute (3 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+20792, dvrr_stack+10084, dvrr_stack+20522, dvrr_stack+3137, dvrr_stack+9949, dvrr_stack+18146);

 /* compute (4 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+21242, dvrr_stack+10354, dvrr_stack+20792, dvrr_stack+3407, dvrr_stack+10084, dvrr_stack+18362);

 /* compute (4 0 | 7 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_kp(Libderiv->CD,dvrr_stack+21917,dvrr_stack+21242,dvrr_stack+18722,15);


 /* compute (0 0 | 1 0) m=12 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+101, Data->F+12, Data->F+13, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+3092, dvrr_stack+2902, dvrr_stack+101, Data->F+11, Data->F+12, NULL);

 /* compute (0 0 | 3 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+2905, dvrr_stack+3092, dvrr_stack+1830, dvrr_stack+2902, NULL);

 /* compute (0 0 | 4 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+114, dvrr_stack+2911, dvrr_stack+21, dvrr_stack+1833, dvrr_stack+2905, NULL);

 /* compute (0 0 | 5 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3092, dvrr_stack+80, dvrr_stack+114, dvrr_stack+1839, dvrr_stack+2911, NULL);

 /* compute (0 0 | 6 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+3113, dvrr_stack+2921, dvrr_stack+3092, dvrr_stack+1849, dvrr_stack+80, NULL);

 /* compute (0 0 | 7 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+3141, dvrr_stack+1794, dvrr_stack+3113, dvrr_stack+2947, dvrr_stack+2921, NULL);

 /* compute (0 0 | 8 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3092, dvrr_stack+950, dvrr_stack+3141, dvrr_stack+4435, dvrr_stack+1794, NULL);

 /* compute (0 0 | 9 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+1794, dvrr_stack+11884, dvrr_stack+3092, dvrr_stack+4463, dvrr_stack+950, NULL);

 /* compute (1 0 | 9 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+23537, dvrr_stack+11929, dvrr_stack+1794, NULL, NULL, dvrr_stack+11884);

 /* compute (2 0 | 9 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+23702, dvrr_stack+11984, dvrr_stack+23537, dvrr_stack+4544, dvrr_stack+11929, dvrr_stack+1930);

 /* compute (3 0 | 9 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+24032, dvrr_stack+12149, dvrr_stack+23702, dvrr_stack+4599, dvrr_stack+11984, dvrr_stack+20522);

 /* compute (4 0 | 9 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 9;
 vrr_build_xxxx(am,Data,dvrr_stack+24582, dvrr_stack+12479, dvrr_stack+24032, dvrr_stack+4929, dvrr_stack+12149, dvrr_stack+20792);

 /* compute (4 0 | 8 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_lp(Libderiv->CD,dvrr_stack+25407,dvrr_stack+24582,dvrr_stack+21242,15);


 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+11929, dvrr_stack+34, dvrr_stack+6, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+11947, dvrr_stack+104, dvrr_stack+40, NULL, NULL, dvrr_stack+34);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+139, dvrr_stack+11947, dvrr_stack+50, dvrr_stack+104, dvrr_stack+40, dvrr_stack+11929);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+2065, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+3);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+11977, dvrr_stack+11929, dvrr_stack+4380, dvrr_stack+34, dvrr_stack+6, dvrr_stack+2065);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+12013, dvrr_stack+139, dvrr_stack+6069, dvrr_stack+11947, dvrr_stack+50, dvrr_stack+11977);

 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+11947, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+11950, dvrr_stack+2065, dvrr_stack+4764, dvrr_stack+3, dvrr_stack+0, dvrr_stack+11947);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+34, dvrr_stack+11977, dvrr_stack+4791, dvrr_stack+11929, dvrr_stack+4380, dvrr_stack+11950);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+12113, dvrr_stack+12013, dvrr_stack+14379, dvrr_stack+139, dvrr_stack+6069, dvrr_stack+34);

 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+34, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+3, dvrr_stack+12, dvrr_stack+31, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+37, dvrr_stack+4764, dvrr_stack+3, dvrr_stack+0, dvrr_stack+12, dvrr_stack+34);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+55, dvrr_stack+95, dvrr_stack+376, NULL, NULL, dvrr_stack+31);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+11929, dvrr_stack+4773, dvrr_stack+55, dvrr_stack+15, dvrr_stack+95, dvrr_stack+3);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+12263, dvrr_stack+4791, dvrr_stack+11929, dvrr_stack+4380, dvrr_stack+4773, dvrr_stack+37);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+0, dvrr_stack+382, dvrr_stack+1006, NULL, NULL, dvrr_stack+376);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+12323, dvrr_stack+4827, dvrr_stack+0, dvrr_stack+129, dvrr_stack+382, dvrr_stack+55);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+30, dvrr_stack+4857, dvrr_stack+12323, dvrr_stack+4398, dvrr_stack+4827, dvrr_stack+11929);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+4544, dvrr_stack+14379, dvrr_stack+30, dvrr_stack+6069, dvrr_stack+4857, dvrr_stack+12263);

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+12263, dvrr_stack+1016, dvrr_stack+1866, NULL, NULL, dvrr_stack+1006);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+14379, dvrr_stack+4325, dvrr_stack+12263, dvrr_stack+392, dvrr_stack+1016, dvrr_stack+0);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4694, dvrr_stack+14479, dvrr_stack+14379, dvrr_stack+6129, dvrr_stack+4325, dvrr_stack+12323);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+23537, dvrr_stack+14569, dvrr_stack+4694, dvrr_stack+6174, dvrr_stack+14479, dvrr_stack+30);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+23762, dvrr_stack+14719, dvrr_stack+23537, dvrr_stack+6264, dvrr_stack+14569, dvrr_stack+4544);

 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4544, dvrr_stack+1881, dvrr_stack+3007, NULL, NULL, dvrr_stack+1866);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+14944, dvrr_stack+4544, dvrr_stack+1031, dvrr_stack+1881, dvrr_stack+12263);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+12263, dvrr_stack+15007, dvrr_stack+0, dvrr_stack+6414, dvrr_stack+14944, dvrr_stack+14379);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+14379, dvrr_stack+15133, dvrr_stack+12263, dvrr_stack+6477, dvrr_stack+15007, dvrr_stack+4694);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+24077, dvrr_stack+15343, dvrr_stack+14379, dvrr_stack+6603, dvrr_stack+15133, dvrr_stack+23537);

 /* compute (1 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+11929, dvrr_stack+3028, dvrr_stack+4435, NULL, NULL, dvrr_stack+3007);

 /* compute (2 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+23537, dvrr_stack+407, dvrr_stack+11929, dvrr_stack+1902, dvrr_stack+3028, dvrr_stack+4544);

 /* compute (3 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4544, dvrr_stack+16333, dvrr_stack+23537, dvrr_stack+7263, dvrr_stack+407, dvrr_stack+0);

 /* compute (4 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+27432, dvrr_stack+16501, dvrr_stack+4544, dvrr_stack+7347, dvrr_stack+16333, dvrr_stack+12263);

 /* compute (5 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+27852, dvrr_stack+16781, dvrr_stack+27432, dvrr_stack+7515, dvrr_stack+16501, dvrr_stack+14379);

 /* compute (1 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+14379, dvrr_stack+4463, dvrr_stack+950, NULL, NULL, dvrr_stack+4435);

 /* compute (2 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+12263, dvrr_stack+1052, dvrr_stack+14379, dvrr_stack+3056, dvrr_stack+4463, dvrr_stack+11929);

 /* compute (3 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+16333, dvrr_stack+18146, dvrr_stack+12263, dvrr_stack+8425, dvrr_stack+1052, dvrr_stack+23537);

 /* compute (4 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+28440, dvrr_stack+18362, dvrr_stack+16333, dvrr_stack+8533, dvrr_stack+18146, dvrr_stack+4544);

 /* compute (5 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+28980, dvrr_stack+18722, dvrr_stack+28440, dvrr_stack+8749, dvrr_stack+18362, dvrr_stack+27432);

 /* compute (1 0 | 8 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+3137, dvrr_stack+11884, dvrr_stack+3092, NULL, NULL, dvrr_stack+950);

 /* compute (2 0 | 8 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+27432, dvrr_stack+1930, dvrr_stack+3137, dvrr_stack+4499, dvrr_stack+11884, dvrr_stack+14379);

 /* compute (3 0 | 8 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+18146, dvrr_stack+20522, dvrr_stack+27432, dvrr_stack+9949, dvrr_stack+1930, dvrr_stack+12263);

 /* compute (4 0 | 8 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+29736, dvrr_stack+20792, dvrr_stack+18146, dvrr_stack+10084, dvrr_stack+20522, dvrr_stack+16333);

 /* compute (5 0 | 8 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 8;
 vrr_build_xxxx(am,Data,dvrr_stack+30411, dvrr_stack+21242, dvrr_stack+29736, dvrr_stack+10354, dvrr_stack+20792, dvrr_stack+28440);

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,90,dvrr_stack+28440, dvrr_stack+680, NULL);
 tmp = dvrr_stack + 28440;
 target_ptr = Libderiv->deriv_classes[2][4][11];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,126,dvrr_stack+28530, dvrr_stack+1416, NULL);
 tmp = dvrr_stack + 28530;
 target_ptr = Libderiv->deriv_classes[2][5][11];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,168,dvrr_stack+28656, dvrr_stack+2398, NULL);
 tmp = dvrr_stack + 28656;
 target_ptr = Libderiv->deriv_classes[2][6][11];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,216,dvrr_stack+12263, dvrr_stack+3677, NULL);
 tmp = dvrr_stack + 12263;
 target_ptr = Libderiv->deriv_classes[2][7][11];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,270,dvrr_stack+29736, dvrr_stack+5259, NULL);
 tmp = dvrr_stack + 29736;
 target_ptr = Libderiv->deriv_classes[2][8][11];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+28824, dvrr_stack+6813, NULL);
 tmp = dvrr_stack + 28824;
 target_ptr = Libderiv->deriv_classes[3][4][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+30006, dvrr_stack+7795, NULL);
 tmp = dvrr_stack + 30006;
 target_ptr = Libderiv->deriv_classes[3][5][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,280,dvrr_stack+1794, dvrr_stack+9109, NULL);
 tmp = dvrr_stack + 1794;
 target_ptr = Libderiv->deriv_classes[3][6][11];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,360,dvrr_stack+16333, dvrr_stack+10804, NULL);
 tmp = dvrr_stack + 16333;
 target_ptr = Libderiv->deriv_classes[3][7][11];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,450,dvrr_stack+20522, dvrr_stack+13029, NULL);
 tmp = dvrr_stack + 20522;
 target_ptr = Libderiv->deriv_classes[3][8][11];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+23537, dvrr_stack+15658, NULL);
 tmp = dvrr_stack + 23537;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+18146, dvrr_stack+17201, NULL);
 tmp = dvrr_stack + 18146;
 target_ptr = Libderiv->deriv_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+27432, dvrr_stack+19262, NULL);
 tmp = dvrr_stack + 27432;
 target_ptr = Libderiv->deriv_classes[4][6][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,540,dvrr_stack+4325, dvrr_stack+21917, NULL);
 tmp = dvrr_stack + 4325;
 target_ptr = Libderiv->deriv_classes[4][7][11];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,675,dvrr_stack+31356, dvrr_stack+25407, NULL);
 tmp = dvrr_stack + 31356;
 target_ptr = Libderiv->deriv_classes[4][8][11];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,90,dvrr_stack+30216, dvrr_stack+680, NULL);
 tmp = dvrr_stack + 30216;
 target_ptr = Libderiv->deriv_classes[2][4][10];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,126,dvrr_stack+20972, dvrr_stack+1416, NULL);
 tmp = dvrr_stack + 20972;
 target_ptr = Libderiv->deriv_classes[2][5][10];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,168,dvrr_stack+18461, dvrr_stack+2398, NULL);
 tmp = dvrr_stack + 18461;
 target_ptr = Libderiv->deriv_classes[2][6][10];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,216,dvrr_stack+9949, dvrr_stack+3677, NULL);
 tmp = dvrr_stack + 9949;
 target_ptr = Libderiv->deriv_classes[2][7][10];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,270,dvrr_stack+14379, dvrr_stack+5259, NULL);
 tmp = dvrr_stack + 14379;
 target_ptr = Libderiv->deriv_classes[2][8][10];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+10165, dvrr_stack+6813, NULL);
 tmp = dvrr_stack + 10165;
 target_ptr = Libderiv->deriv_classes[3][4][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+950, dvrr_stack+7795, NULL);
 tmp = dvrr_stack + 950;
 target_ptr = Libderiv->deriv_classes[3][5][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,280,dvrr_stack+8425, dvrr_stack+9109, NULL);
 tmp = dvrr_stack + 8425;
 target_ptr = Libderiv->deriv_classes[3][6][10];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,360,dvrr_stack+14944, dvrr_stack+10804, NULL);
 tmp = dvrr_stack + 14944;
 target_ptr = Libderiv->deriv_classes[3][7][10];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,450,dvrr_stack+32031, dvrr_stack+13029, NULL);
 tmp = dvrr_stack + 32031;
 target_ptr = Libderiv->deriv_classes[3][8][10];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+7263, dvrr_stack+15658, NULL);
 tmp = dvrr_stack + 7263;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+2902, dvrr_stack+17201, NULL);
 tmp = dvrr_stack + 2902;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+32481, dvrr_stack+19262, NULL);
 tmp = dvrr_stack + 32481;
 target_ptr = Libderiv->deriv_classes[4][6][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,540,dvrr_stack+32901, dvrr_stack+21917, NULL);
 tmp = dvrr_stack + 32901;
 target_ptr = Libderiv->deriv_classes[4][7][10];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,675,dvrr_stack+33441, dvrr_stack+25407, NULL);
 tmp = dvrr_stack + 33441;
 target_ptr = Libderiv->deriv_classes[4][8][10];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,90,dvrr_stack+30306, dvrr_stack+680, NULL);
 tmp = dvrr_stack + 30306;
 target_ptr = Libderiv->deriv_classes[2][4][9];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,126,dvrr_stack+680, dvrr_stack+1416, NULL);
 tmp = dvrr_stack + 680;
 target_ptr = Libderiv->deriv_classes[2][5][9];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,168,dvrr_stack+1416, dvrr_stack+2398, NULL);
 tmp = dvrr_stack + 1416;
 target_ptr = Libderiv->deriv_classes[2][6][9];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,216,dvrr_stack+2398, dvrr_stack+3677, NULL);
 tmp = dvrr_stack + 2398;
 target_ptr = Libderiv->deriv_classes[2][7][9];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,270,dvrr_stack+3677, dvrr_stack+5259, NULL);
 tmp = dvrr_stack + 3677;
 target_ptr = Libderiv->deriv_classes[2][8][9];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+5259, dvrr_stack+6813, NULL);
 tmp = dvrr_stack + 5259;
 target_ptr = Libderiv->deriv_classes[3][4][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+1584, dvrr_stack+7795, NULL);
 tmp = dvrr_stack + 1584;
 target_ptr = Libderiv->deriv_classes[3][5][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,280,dvrr_stack+7795, dvrr_stack+9109, NULL);
 tmp = dvrr_stack + 7795;
 target_ptr = Libderiv->deriv_classes[3][6][9];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,360,dvrr_stack+9109, dvrr_stack+10804, NULL);
 tmp = dvrr_stack + 9109;
 target_ptr = Libderiv->deriv_classes[3][7][9];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,450,dvrr_stack+6813, dvrr_stack+13029, NULL);
 tmp = dvrr_stack + 6813;
 target_ptr = Libderiv->deriv_classes[3][8][9];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+13029, dvrr_stack+15658, NULL);
 tmp = dvrr_stack + 13029;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+15658, dvrr_stack+17201, NULL);
 tmp = dvrr_stack + 15658;
 target_ptr = Libderiv->deriv_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+17201, dvrr_stack+19262, NULL);
 tmp = dvrr_stack + 17201;
 target_ptr = Libderiv->deriv_classes[4][6][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,540,dvrr_stack+19262, dvrr_stack+21917, NULL);
 tmp = dvrr_stack + 19262;
 target_ptr = Libderiv->deriv_classes[4][7][9];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,675,dvrr_stack+21917, dvrr_stack+25407, NULL);
 tmp = dvrr_stack + 21917;
 target_ptr = Libderiv->deriv_classes[4][8][9];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,6,1,dvrr_stack+25407, dvrr_stack+554, dvrr_stack+139);
 tmp = dvrr_stack + 25407;
 target_ptr = Libderiv->deriv_classes[2][4][8];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,6,1,dvrr_stack+25497, dvrr_stack+1248, dvrr_stack+244);
 tmp = dvrr_stack + 25497;
 target_ptr = Libderiv->deriv_classes[2][5][8];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,6,1,dvrr_stack+25623, dvrr_stack+2182, dvrr_stack+554);
 tmp = dvrr_stack + 25623;
 target_ptr = Libderiv->deriv_classes[2][6][8];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,6,1,dvrr_stack+25791, dvrr_stack+3407, dvrr_stack+1248);
 tmp = dvrr_stack + 25791;
 target_ptr = Libderiv->deriv_classes[2][7][8];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,6,1,dvrr_stack+26007, dvrr_stack+4929, dvrr_stack+2182);
 tmp = dvrr_stack + 26007;
 target_ptr = Libderiv->deriv_classes[2][8][8];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,10,1,dvrr_stack+26277, dvrr_stack+6603, dvrr_stack+12013);
 tmp = dvrr_stack + 26277;
 target_ptr = Libderiv->deriv_classes[3][4][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,10,1,dvrr_stack+26427, dvrr_stack+7515, dvrr_stack+6264);
 tmp = dvrr_stack + 26427;
 target_ptr = Libderiv->deriv_classes[3][5][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,10,1,dvrr_stack+26637, dvrr_stack+8749, dvrr_stack+6603);
 tmp = dvrr_stack + 26637;
 target_ptr = Libderiv->deriv_classes[3][6][8];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,10,1,dvrr_stack+15973, dvrr_stack+10354, dvrr_stack+7515);
 tmp = dvrr_stack + 15973;
 target_ptr = Libderiv->deriv_classes[3][7][8];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,10,1,dvrr_stack+26917, dvrr_stack+12479, dvrr_stack+8749);
 tmp = dvrr_stack + 26917;
 target_ptr = Libderiv->deriv_classes[3][8][8];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+22592, dvrr_stack+15343, dvrr_stack+12113);
 tmp = dvrr_stack + 22592;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+22817, dvrr_stack+16781, dvrr_stack+14719);
 tmp = dvrr_stack + 22817;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+19802, dvrr_stack+18722, dvrr_stack+15343);
 tmp = dvrr_stack + 19802;
 target_ptr = Libderiv->deriv_classes[4][6][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_k(Data,15,1,dvrr_stack+13254, dvrr_stack+21242, dvrr_stack+16781);
 tmp = dvrr_stack + 13254;
 target_ptr = Libderiv->deriv_classes[4][7][8];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_l(Data,15,1,dvrr_stack+10804, dvrr_stack+24582, dvrr_stack+18722);
 tmp = dvrr_stack + 10804;
 target_ptr = Libderiv->deriv_classes[4][8][8];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,6,1,dvrr_stack+23132, dvrr_stack+554, dvrr_stack+139);
 tmp = dvrr_stack + 23132;
 target_ptr = Libderiv->deriv_classes[2][4][7];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,6,1,dvrr_stack+23222, dvrr_stack+1248, dvrr_stack+244);
 tmp = dvrr_stack + 23222;
 target_ptr = Libderiv->deriv_classes[2][5][7];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,6,1,dvrr_stack+23348, dvrr_stack+2182, dvrr_stack+554);
 tmp = dvrr_stack + 23348;
 target_ptr = Libderiv->deriv_classes[2][6][7];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,6,1,dvrr_stack+20222, dvrr_stack+3407, dvrr_stack+1248);
 tmp = dvrr_stack + 20222;
 target_ptr = Libderiv->deriv_classes[2][7][7];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,6,1,dvrr_stack+17621, dvrr_stack+4929, dvrr_stack+2182);
 tmp = dvrr_stack + 17621;
 target_ptr = Libderiv->deriv_classes[2][8][7];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,10,1,dvrr_stack+17891, dvrr_stack+6603, dvrr_stack+12013);
 tmp = dvrr_stack + 17891;
 target_ptr = Libderiv->deriv_classes[3][4][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,10,1,dvrr_stack+13794, dvrr_stack+7515, dvrr_stack+6264);
 tmp = dvrr_stack + 13794;
 target_ptr = Libderiv->deriv_classes[3][5][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,10,1,dvrr_stack+14004, dvrr_stack+8749, dvrr_stack+6603);
 tmp = dvrr_stack + 14004;
 target_ptr = Libderiv->deriv_classes[3][6][7];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,10,1,dvrr_stack+11479, dvrr_stack+10354, dvrr_stack+7515);
 tmp = dvrr_stack + 11479;
 target_ptr = Libderiv->deriv_classes[3][7][7];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,10,1,dvrr_stack+9469, dvrr_stack+12479, dvrr_stack+8749);
 tmp = dvrr_stack + 9469;
 target_ptr = Libderiv->deriv_classes[3][8][7];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+8075, dvrr_stack+15343, dvrr_stack+12113);
 tmp = dvrr_stack + 8075;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+5409, dvrr_stack+16781, dvrr_stack+14719);
 tmp = dvrr_stack + 5409;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+5724, dvrr_stack+18722, dvrr_stack+15343);
 tmp = dvrr_stack + 5724;
 target_ptr = Libderiv->deriv_classes[4][6][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_k(Data,15,1,dvrr_stack+34116, dvrr_stack+21242, dvrr_stack+16781);
 tmp = dvrr_stack + 34116;
 target_ptr = Libderiv->deriv_classes[4][7][7];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_l(Data,15,1,dvrr_stack+34656, dvrr_stack+24582, dvrr_stack+18722);
 tmp = dvrr_stack + 34656;
 target_ptr = Libderiv->deriv_classes[4][8][7];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,6,1,dvrr_stack+18041, dvrr_stack+554, dvrr_stack+139);
 tmp = dvrr_stack + 18041;
 target_ptr = Libderiv->deriv_classes[2][4][6];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,6,1,dvrr_stack+11839, dvrr_stack+1248, dvrr_stack+244);
 tmp = dvrr_stack + 11839;
 target_ptr = Libderiv->deriv_classes[2][5][6];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,6,1,dvrr_stack+3947, dvrr_stack+2182, dvrr_stack+554);
 tmp = dvrr_stack + 3947;
 target_ptr = Libderiv->deriv_classes[2][6][6];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,6,1,dvrr_stack+2614, dvrr_stack+3407, dvrr_stack+1248);
 tmp = dvrr_stack + 2614;
 target_ptr = Libderiv->deriv_classes[2][7][6];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,6,1,dvrr_stack+35331, dvrr_stack+4929, dvrr_stack+2182);
 tmp = dvrr_stack + 35331;
 target_ptr = Libderiv->deriv_classes[2][8][6];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,10,1,dvrr_stack+4115, dvrr_stack+6603, dvrr_stack+12013);
 tmp = dvrr_stack + 4115;
 target_ptr = Libderiv->deriv_classes[3][4][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,10,1,dvrr_stack+4865, dvrr_stack+7515, dvrr_stack+6264);
 tmp = dvrr_stack + 4865;
 target_ptr = Libderiv->deriv_classes[3][5][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,10,1,dvrr_stack+35601, dvrr_stack+8749, dvrr_stack+6603);
 tmp = dvrr_stack + 35601;
 target_ptr = Libderiv->deriv_classes[3][6][6];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,10,1,dvrr_stack+35881, dvrr_stack+10354, dvrr_stack+7515);
 tmp = dvrr_stack + 35881;
 target_ptr = Libderiv->deriv_classes[3][7][6];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,10,1,dvrr_stack+36241, dvrr_stack+12479, dvrr_stack+8749);
 tmp = dvrr_stack + 36241;
 target_ptr = Libderiv->deriv_classes[3][8][6];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+12479, dvrr_stack+15343, dvrr_stack+12113);
 tmp = dvrr_stack + 12479;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+12704, dvrr_stack+16781, dvrr_stack+14719);
 tmp = dvrr_stack + 12704;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+36691, dvrr_stack+18722, dvrr_stack+15343);
 tmp = dvrr_stack + 36691;
 target_ptr = Libderiv->deriv_classes[4][6][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_k(Data,15,1,dvrr_stack+37111, dvrr_stack+21242, dvrr_stack+16781);
 tmp = dvrr_stack + 37111;
 target_ptr = Libderiv->deriv_classes[4][7][6];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_l(Data,15,1,dvrr_stack+37651, dvrr_stack+24582, dvrr_stack+18722);
 tmp = dvrr_stack + 37651;
 target_ptr = Libderiv->deriv_classes[4][8][6];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,15,dvrr_stack+14284, dvrr_stack+6264, dvrr_stack+199);
 tmp = dvrr_stack + 14284;
 target_ptr = Libderiv->deriv_classes[2][4][2];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,21,dvrr_stack+11965, dvrr_stack+6603, dvrr_stack+491);
 tmp = dvrr_stack + 11965;
 target_ptr = Libderiv->deriv_classes[2][5][2];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,28,dvrr_stack+12091, dvrr_stack+7515, dvrr_stack+1164);
 tmp = dvrr_stack + 12091;
 target_ptr = Libderiv->deriv_classes[2][6][2];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,36,dvrr_stack+24518, dvrr_stack+8749, dvrr_stack+2074);
 tmp = dvrr_stack + 24518;
 target_ptr = Libderiv->deriv_classes[2][7][2];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,45,dvrr_stack+24734, dvrr_stack+10354, dvrr_stack+3272);
 tmp = dvrr_stack + 24734;
 target_ptr = Libderiv->deriv_classes[2][8][2];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,15,dvrr_stack+0, dvrr_stack+14719, dvrr_stack+244);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[3][4][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,21,dvrr_stack+25004, dvrr_stack+15343, dvrr_stack+554);
 tmp = dvrr_stack + 25004;
 target_ptr = Libderiv->deriv_classes[3][5][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,28,dvrr_stack+38326, dvrr_stack+16781, dvrr_stack+1248);
 tmp = dvrr_stack + 38326;
 target_ptr = Libderiv->deriv_classes[3][6][2];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,36,dvrr_stack+38606, dvrr_stack+18722, dvrr_stack+2182);
 tmp = dvrr_stack + 38606;
 target_ptr = Libderiv->deriv_classes[3][7][2];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_f(Data,45,dvrr_stack+38966, dvrr_stack+21242, dvrr_stack+3407);
 tmp = dvrr_stack + 38966;
 target_ptr = Libderiv->deriv_classes[3][8][2];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+39416, dvrr_stack+23762, dvrr_stack+6264);
 tmp = dvrr_stack + 39416;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+39641, dvrr_stack+24077, dvrr_stack+6603);
 tmp = dvrr_stack + 39641;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+39956, dvrr_stack+27852, dvrr_stack+7515);
 tmp = dvrr_stack + 39956;
 target_ptr = Libderiv->deriv_classes[4][6][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,36,dvrr_stack+40376, dvrr_stack+28980, dvrr_stack+8749);
 tmp = dvrr_stack + 40376;
 target_ptr = Libderiv->deriv_classes[4][7][2];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,45,dvrr_stack+40916, dvrr_stack+30411, dvrr_stack+10354);
 tmp = dvrr_stack + 40916;
 target_ptr = Libderiv->deriv_classes[4][8][2];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,15,dvrr_stack+8300, dvrr_stack+6264, dvrr_stack+199);
 tmp = dvrr_stack + 8300;
 target_ptr = Libderiv->deriv_classes[2][4][1];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,21,dvrr_stack+806, dvrr_stack+6603, dvrr_stack+491);
 tmp = dvrr_stack + 806;
 target_ptr = Libderiv->deriv_classes[2][5][1];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,28,dvrr_stack+6414, dvrr_stack+7515, dvrr_stack+1164);
 tmp = dvrr_stack + 6414;
 target_ptr = Libderiv->deriv_classes[2][6][1];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,36,dvrr_stack+41591, dvrr_stack+8749, dvrr_stack+2074);
 tmp = dvrr_stack + 41591;
 target_ptr = Libderiv->deriv_classes[2][7][1];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,45,dvrr_stack+41807, dvrr_stack+10354, dvrr_stack+3272);
 tmp = dvrr_stack + 41807;
 target_ptr = Libderiv->deriv_classes[2][8][1];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,15,dvrr_stack+5075, dvrr_stack+14719, dvrr_stack+244);
 tmp = dvrr_stack + 5075;
 target_ptr = Libderiv->deriv_classes[3][4][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,21,dvrr_stack+42077, dvrr_stack+15343, dvrr_stack+554);
 tmp = dvrr_stack + 42077;
 target_ptr = Libderiv->deriv_classes[3][5][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,28,dvrr_stack+42287, dvrr_stack+16781, dvrr_stack+1248);
 tmp = dvrr_stack + 42287;
 target_ptr = Libderiv->deriv_classes[3][6][1];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,36,dvrr_stack+42567, dvrr_stack+18722, dvrr_stack+2182);
 tmp = dvrr_stack + 42567;
 target_ptr = Libderiv->deriv_classes[3][7][1];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_f(Data,45,dvrr_stack+42927, dvrr_stack+21242, dvrr_stack+3407);
 tmp = dvrr_stack + 42927;
 target_ptr = Libderiv->deriv_classes[3][8][1];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+43377, dvrr_stack+23762, dvrr_stack+6264);
 tmp = dvrr_stack + 43377;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+43602, dvrr_stack+24077, dvrr_stack+6603);
 tmp = dvrr_stack + 43602;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+43917, dvrr_stack+27852, dvrr_stack+7515);
 tmp = dvrr_stack + 43917;
 target_ptr = Libderiv->deriv_classes[4][6][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,36,dvrr_stack+44337, dvrr_stack+28980, dvrr_stack+8749);
 tmp = dvrr_stack + 44337;
 target_ptr = Libderiv->deriv_classes[4][7][1];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,45,dvrr_stack+44877, dvrr_stack+30411, dvrr_stack+10354);
 tmp = dvrr_stack + 44877;
 target_ptr = Libderiv->deriv_classes[4][8][1];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,15,dvrr_stack+6144, dvrr_stack+6264, dvrr_stack+199);
 tmp = dvrr_stack + 6144;
 target_ptr = Libderiv->deriv_classes[2][4][0];
 for(i=0;i<90;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,21,dvrr_stack+21098, dvrr_stack+6603, dvrr_stack+491);
 tmp = dvrr_stack + 21098;
 target_ptr = Libderiv->deriv_classes[2][5][0];
 for(i=0;i<126;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,28,dvrr_stack+25214, dvrr_stack+7515, dvrr_stack+1164);
 tmp = dvrr_stack + 25214;
 target_ptr = Libderiv->deriv_classes[2][6][0];
 for(i=0;i<168;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,36,dvrr_stack+334, dvrr_stack+8749, dvrr_stack+2074);
 tmp = dvrr_stack + 334;
 target_ptr = Libderiv->deriv_classes[2][7][0];
 for(i=0;i<216;i++)
   target_ptr[i] += tmp[i];

 /* compute (2 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,45,dvrr_stack+45552, dvrr_stack+10354, dvrr_stack+3272);
 tmp = dvrr_stack + 45552;
 target_ptr = Libderiv->deriv_classes[2][8][0];
 for(i=0;i<270;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,15,dvrr_stack+3217, dvrr_stack+14719, dvrr_stack+244);
 tmp = dvrr_stack + 3217;
 target_ptr = Libderiv->deriv_classes[3][4][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,21,dvrr_stack+14649, dvrr_stack+15343, dvrr_stack+554);
 tmp = dvrr_stack + 14649;
 target_ptr = Libderiv->deriv_classes[3][5][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,28,dvrr_stack+15304, dvrr_stack+16781, dvrr_stack+1248);
 tmp = dvrr_stack + 15304;
 target_ptr = Libderiv->deriv_classes[3][6][0];
 for(i=0;i<280;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,36,dvrr_stack+16693, dvrr_stack+18722, dvrr_stack+2182);
 tmp = dvrr_stack + 16693;
 target_ptr = Libderiv->deriv_classes[3][7][0];
 for(i=0;i<360;i++)
   target_ptr[i] += tmp[i];

 /* compute (3 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_f(Data,45,dvrr_stack+18629, dvrr_stack+21242, dvrr_stack+3407);
 tmp = dvrr_stack + 18629;
 target_ptr = Libderiv->deriv_classes[3][8][0];
 for(i=0;i<450;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+2074, dvrr_stack+23762, dvrr_stack+6264);
 tmp = dvrr_stack + 2074;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+23762, dvrr_stack+24077, dvrr_stack+6603);
 tmp = dvrr_stack + 23762;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+24077, dvrr_stack+27852, dvrr_stack+7515);
 tmp = dvrr_stack + 24077;
 target_ptr = Libderiv->deriv_classes[4][6][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 7 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,36,dvrr_stack+27852, dvrr_stack+28980, dvrr_stack+8749);
 tmp = dvrr_stack + 27852;
 target_ptr = Libderiv->deriv_classes[4][7][0];
 for(i=0;i<540;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 8 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,45,dvrr_stack+28974, dvrr_stack+30411, dvrr_stack+10354);
 tmp = dvrr_stack + 28974;
 target_ptr = Libderiv->deriv_classes[4][8][0];
 for(i=0;i<675;i++)
   target_ptr[i] += tmp[i];


}

