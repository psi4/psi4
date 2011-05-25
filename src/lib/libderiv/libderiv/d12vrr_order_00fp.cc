#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|fp) integrals */

void d12vrr_order_00fp(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+9, dvrr_stack+0, dvrr_stack+6, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+9, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+31, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+34, dvrr_stack+6, dvrr_stack+31, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+9, dvrr_stack+34, dvrr_stack+0, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+50, dvrr_stack+21, dvrr_stack+40, dvrr_stack+15, dvrr_stack+9, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+65, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+68, dvrr_stack+65, dvrr_stack+3, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+74, dvrr_stack+68, dvrr_stack+15, dvrr_stack+65, dvrr_stack+3, NULL);
 tmp = dvrr_stack + 74;
 target_ptr = Libderiv->dvrr_classes[0][3];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+84, dvrr_stack+74, dvrr_stack+21, dvrr_stack+68, dvrr_stack+15, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+99, dvrr_stack+84, dvrr_stack+50, NULL, NULL, dvrr_stack+21);

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+144,dvrr_stack+84,dvrr_stack+74,1);


 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+174, dvrr_stack+84, dvrr_stack+50, dvrr_stack+74, dvrr_stack+21, NULL);

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+195,dvrr_stack+174,dvrr_stack+84,1);


 /* compute (0 0 | 3 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fd(Libderiv->CD,dvrr_stack+240,dvrr_stack+195,dvrr_stack+144,1);


 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,10,dvrr_stack+300, dvrr_stack+240, dvrr_stack+74);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+330, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+333, dvrr_stack+31, dvrr_stack+330, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+339, dvrr_stack+34, dvrr_stack+333, dvrr_stack+6, dvrr_stack+31, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+349, dvrr_stack+40, dvrr_stack+339, dvrr_stack+9, dvrr_stack+34, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+364, dvrr_stack+50, dvrr_stack+349, dvrr_stack+21, dvrr_stack+40, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+385, dvrr_stack+174, dvrr_stack+364, dvrr_stack+84, dvrr_stack+50, NULL);

 /* compute (0 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+413,dvrr_stack+385,dvrr_stack+174,1);


 /* compute (0 0 | 4 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gd(Libderiv->CD,dvrr_stack+476,dvrr_stack+413,dvrr_stack+195,1);


 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,15,dvrr_stack+566, dvrr_stack+476, dvrr_stack+84);

 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,10,dvrr_stack+611, dvrr_stack+240, dvrr_stack+74);

 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,15,dvrr_stack+641, dvrr_stack+476, dvrr_stack+84);

 /* compute (0 0 | 3 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,10,dvrr_stack+686, dvrr_stack+240, dvrr_stack+74);

 /* compute (0 0 | 4 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,15,dvrr_stack+240, dvrr_stack+476, dvrr_stack+84);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+476,dvrr_stack+74,dvrr_stack+68,1);


 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+494, dvrr_stack+476, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+285, dvrr_stack+195, NULL);
 tmp = dvrr_stack + 285;
 target_ptr = Libderiv->deriv_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+500, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 500;
 target_ptr = Libderiv->deriv_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,21,dvrr_stack+510, dvrr_stack+413, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+531, dvrr_stack+476, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+537, dvrr_stack+195, NULL);
 tmp = dvrr_stack + 537;
 target_ptr = Libderiv->deriv_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+552, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 552;
 target_ptr = Libderiv->deriv_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+716, dvrr_stack+413, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+31, dvrr_stack+476, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+476, dvrr_stack+195, NULL);
 tmp = dvrr_stack + 476;
 target_ptr = Libderiv->deriv_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+195, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 195;
 target_ptr = Libderiv->deriv_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+144, dvrr_stack+413, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+413, dvrr_stack+74, dvrr_stack+65);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+419, dvrr_stack+174, dvrr_stack+74);
 tmp = dvrr_stack + 419;
 target_ptr = Libderiv->deriv_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+434, dvrr_stack+84, dvrr_stack+68);
 tmp = dvrr_stack + 434;
 target_ptr = Libderiv->deriv_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+444, dvrr_stack+385, dvrr_stack+84);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+465, dvrr_stack+74, dvrr_stack+65);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+205, dvrr_stack+174, dvrr_stack+74);
 tmp = dvrr_stack + 205;
 target_ptr = Libderiv->deriv_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+220, dvrr_stack+84, dvrr_stack+68);
 tmp = dvrr_stack + 220;
 target_ptr = Libderiv->deriv_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+737, dvrr_stack+385, dvrr_stack+84);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+165, dvrr_stack+74, dvrr_stack+65);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+330, dvrr_stack+174, dvrr_stack+74);
 tmp = dvrr_stack + 330;
 target_ptr = Libderiv->deriv_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+230, dvrr_stack+84, dvrr_stack+68);
 tmp = dvrr_stack + 230;
 target_ptr = Libderiv->deriv_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+758, dvrr_stack+385, dvrr_stack+84);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+779, dvrr_stack+74, dvrr_stack+21, NULL, NULL, dvrr_stack+15);

 /* compute (1 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+809,dvrr_stack+99,dvrr_stack+779,3);


 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,30,dvrr_stack+899, dvrr_stack+809, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+929, dvrr_stack+174, dvrr_stack+364, NULL, NULL, dvrr_stack+50);

 /* compute (1 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+992,dvrr_stack+929,dvrr_stack+99,3);


 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,45,dvrr_stack+364, dvrr_stack+992, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,30,dvrr_stack+1127, dvrr_stack+809, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,45,dvrr_stack+1157, dvrr_stack+992, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,30,dvrr_stack+1202, dvrr_stack+809, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,45,dvrr_stack+809, dvrr_stack+992, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+992, dvrr_stack+68, dvrr_stack+15, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,3,1,dvrr_stack+1010, dvrr_stack+99, dvrr_stack+992);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,3,1,dvrr_stack+854, dvrr_stack+929, dvrr_stack+779);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,3,1,dvrr_stack+1040, dvrr_stack+99, dvrr_stack+992);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,3,1,dvrr_stack+1070, dvrr_stack+929, dvrr_stack+779);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,3,1,dvrr_stack+1232, dvrr_stack+99, dvrr_stack+992);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,3,1,dvrr_stack+1262, dvrr_stack+929, dvrr_stack+779);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+929, dvrr_stack+15, dvrr_stack+9, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+947, dvrr_stack+21, dvrr_stack+40, NULL, NULL, dvrr_stack+9);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+1307, dvrr_stack+779, dvrr_stack+947, dvrr_stack+74, dvrr_stack+21, dvrr_stack+929);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,10,dvrr_stack+0, dvrr_stack+1307, dvrr_stack+74);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1367, dvrr_stack+50, dvrr_stack+349, NULL, NULL, dvrr_stack+40);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1412, dvrr_stack+99, dvrr_stack+1367, dvrr_stack+84, dvrr_stack+50, dvrr_stack+947);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,15,dvrr_stack+1367, dvrr_stack+1412, dvrr_stack+84);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,10,dvrr_stack+929, dvrr_stack+1307, dvrr_stack+74);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,15,dvrr_stack+959, dvrr_stack+1412, dvrr_stack+84);

 /* compute (1 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,10,dvrr_stack+37, dvrr_stack+1307, dvrr_stack+74);

 /* compute (1 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,15,dvrr_stack+1307, dvrr_stack+1412, dvrr_stack+84);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1352, dvrr_stack+99, NULL);
 tmp = dvrr_stack + 1352;
 target_ptr = Libderiv->deriv_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+1412, dvrr_stack+99, NULL);
 tmp = dvrr_stack + 1412;
 target_ptr = Libderiv->deriv_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1427, dvrr_stack+99, NULL);
 tmp = dvrr_stack + 1427;
 target_ptr = Libderiv->deriv_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,10,dvrr_stack+1442, dvrr_stack+300, NULL);
 tmp = dvrr_stack + 1442;
 target_ptr = Libderiv->deriv2_classes[0][3][143];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,15,dvrr_stack+1452, dvrr_stack+566, NULL);
 tmp = dvrr_stack + 1452;
 target_ptr = Libderiv->deriv2_classes[0][4][143];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,10,dvrr_stack+1467, dvrr_stack+300, NULL);
 tmp = dvrr_stack + 1467;
 target_ptr = Libderiv->deriv2_classes[0][3][131];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,15,dvrr_stack+1477, dvrr_stack+566, NULL);
 tmp = dvrr_stack + 1477;
 target_ptr = Libderiv->deriv2_classes[0][4][131];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+1492, dvrr_stack+611, NULL);
 tmp = dvrr_stack + 1492;
 target_ptr = Libderiv->deriv2_classes[0][3][130];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+1502, dvrr_stack+641, NULL);
 tmp = dvrr_stack + 1502;
 target_ptr = Libderiv->deriv2_classes[0][4][130];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,10,dvrr_stack+1517, dvrr_stack+300, NULL);
 tmp = dvrr_stack + 1517;
 target_ptr = Libderiv->deriv2_classes[0][3][119];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,15,dvrr_stack+300, dvrr_stack+566, NULL);
 tmp = dvrr_stack + 300;
 target_ptr = Libderiv->deriv2_classes[0][4][119];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+315, dvrr_stack+611, NULL);
 tmp = dvrr_stack + 315;
 target_ptr = Libderiv->deriv2_classes[0][3][118];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+1527, dvrr_stack+641, NULL);
 tmp = dvrr_stack + 1527;
 target_ptr = Libderiv->deriv2_classes[0][4][118];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+1542, dvrr_stack+686, NULL);
 tmp = dvrr_stack + 1542;
 target_ptr = Libderiv->deriv2_classes[0][3][117];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+1552, dvrr_stack+240, NULL);
 tmp = dvrr_stack + 1552;
 target_ptr = Libderiv->deriv2_classes[0][4][117];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+240, dvrr_stack+285, dvrr_stack+494);
 tmp = dvrr_stack + 240;
 target_ptr = Libderiv->deriv2_classes[0][3][107];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+250, dvrr_stack+510, dvrr_stack+500);
 tmp = dvrr_stack + 250;
 target_ptr = Libderiv->deriv2_classes[0][4][107];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+265, dvrr_stack+537, dvrr_stack+531);
 tmp = dvrr_stack + 265;
 target_ptr = Libderiv->deriv2_classes[0][3][106];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+1567, dvrr_stack+716, dvrr_stack+552);
 tmp = dvrr_stack + 1567;
 target_ptr = Libderiv->deriv2_classes[0][4][106];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+275, dvrr_stack+476, dvrr_stack+31);
 tmp = dvrr_stack + 275;
 target_ptr = Libderiv->deriv2_classes[0][3][105];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+1582, dvrr_stack+144, dvrr_stack+195);
 tmp = dvrr_stack + 1582;
 target_ptr = Libderiv->deriv2_classes[0][4][105];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+1597, dvrr_stack+419, dvrr_stack+413);
 tmp = dvrr_stack + 1597;
 target_ptr = Libderiv->deriv2_classes[0][3][104];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+1607, dvrr_stack+444, dvrr_stack+434);
 tmp = dvrr_stack + 1607;
 target_ptr = Libderiv->deriv2_classes[0][4][104];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+1622, dvrr_stack+285, dvrr_stack+494);
 tmp = dvrr_stack + 1622;
 target_ptr = Libderiv->deriv2_classes[0][3][95];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+1632, dvrr_stack+510, dvrr_stack+500);
 tmp = dvrr_stack + 1632;
 target_ptr = Libderiv->deriv2_classes[0][4][95];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+1647, dvrr_stack+537, dvrr_stack+531);
 tmp = dvrr_stack + 1647;
 target_ptr = Libderiv->deriv2_classes[0][3][94];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+1657, dvrr_stack+716, dvrr_stack+552);
 tmp = dvrr_stack + 1657;
 target_ptr = Libderiv->deriv2_classes[0][4][94];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+1672, dvrr_stack+476, dvrr_stack+31);
 tmp = dvrr_stack + 1672;
 target_ptr = Libderiv->deriv2_classes[0][3][93];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+1682, dvrr_stack+144, dvrr_stack+195);
 tmp = dvrr_stack + 1682;
 target_ptr = Libderiv->deriv2_classes[0][4][93];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+1697, dvrr_stack+419, dvrr_stack+413);
 tmp = dvrr_stack + 1697;
 target_ptr = Libderiv->deriv2_classes[0][3][92];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+1707, dvrr_stack+444, dvrr_stack+434);
 tmp = dvrr_stack + 1707;
 target_ptr = Libderiv->deriv2_classes[0][4][92];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+1722, dvrr_stack+205, dvrr_stack+465);
 tmp = dvrr_stack + 1722;
 target_ptr = Libderiv->deriv2_classes[0][3][91];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+1732, dvrr_stack+737, dvrr_stack+220);
 tmp = dvrr_stack + 1732;
 target_ptr = Libderiv->deriv2_classes[0][4][91];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1747, dvrr_stack+285, dvrr_stack+494);
 tmp = dvrr_stack + 1747;
 target_ptr = Libderiv->deriv2_classes[0][3][83];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+285, dvrr_stack+510, dvrr_stack+500);
 tmp = dvrr_stack + 285;
 target_ptr = Libderiv->deriv2_classes[0][4][83];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1757, dvrr_stack+537, dvrr_stack+531);
 tmp = dvrr_stack + 1757;
 target_ptr = Libderiv->deriv2_classes[0][3][82];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+1767, dvrr_stack+716, dvrr_stack+552);
 tmp = dvrr_stack + 1767;
 target_ptr = Libderiv->deriv2_classes[0][4][82];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1782, dvrr_stack+476, dvrr_stack+31);
 tmp = dvrr_stack + 1782;
 target_ptr = Libderiv->deriv2_classes[0][3][81];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+1792, dvrr_stack+144, dvrr_stack+195);
 tmp = dvrr_stack + 1792;
 target_ptr = Libderiv->deriv2_classes[0][4][81];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1807, dvrr_stack+419, dvrr_stack+413);
 tmp = dvrr_stack + 1807;
 target_ptr = Libderiv->deriv2_classes[0][3][80];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+1817, dvrr_stack+444, dvrr_stack+434);
 tmp = dvrr_stack + 1817;
 target_ptr = Libderiv->deriv2_classes[0][4][80];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1832, dvrr_stack+205, dvrr_stack+465);
 tmp = dvrr_stack + 1832;
 target_ptr = Libderiv->deriv2_classes[0][3][79];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+1842, dvrr_stack+737, dvrr_stack+220);
 tmp = dvrr_stack + 1842;
 target_ptr = Libderiv->deriv2_classes[0][4][79];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+1857, dvrr_stack+330, dvrr_stack+165);
 tmp = dvrr_stack + 1857;
 target_ptr = Libderiv->deriv2_classes[0][3][78];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+325, dvrr_stack+758, dvrr_stack+230);
 tmp = dvrr_stack + 325;
 target_ptr = Libderiv->deriv2_classes[0][4][78];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,10,dvrr_stack+340, dvrr_stack+899, NULL);
 tmp = dvrr_stack + 340;
 target_ptr = Libderiv->deriv2_classes[0][3][35];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1867, dvrr_stack+364, NULL);
 tmp = dvrr_stack + 1867;
 target_ptr = Libderiv->deriv2_classes[0][4][35];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+350, dvrr_stack+1127, NULL);
 tmp = dvrr_stack + 350;
 target_ptr = Libderiv->deriv2_classes[0][3][34];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1882, dvrr_stack+1157, NULL);
 tmp = dvrr_stack + 1882;
 target_ptr = Libderiv->deriv2_classes[0][4][34];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1897, dvrr_stack+1202, NULL);
 tmp = dvrr_stack + 1897;
 target_ptr = Libderiv->deriv2_classes[0][3][33];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1907, dvrr_stack+809, NULL);
 tmp = dvrr_stack + 1907;
 target_ptr = Libderiv->deriv2_classes[0][4][33];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1922, dvrr_stack+1010, NULL);
 tmp = dvrr_stack + 1922;
 target_ptr = Libderiv->deriv2_classes[0][3][32];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1932, dvrr_stack+854, NULL);
 tmp = dvrr_stack + 1932;
 target_ptr = Libderiv->deriv2_classes[0][4][32];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1947, dvrr_stack+1040, NULL);
 tmp = dvrr_stack + 1947;
 target_ptr = Libderiv->deriv2_classes[0][3][31];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+1957, dvrr_stack+1070, NULL);
 tmp = dvrr_stack + 1957;
 target_ptr = Libderiv->deriv2_classes[0][4][31];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1972, dvrr_stack+779, NULL);
 tmp = dvrr_stack + 1972;
 target_ptr = Libderiv->deriv_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1982, dvrr_stack+1232, NULL);
 tmp = dvrr_stack + 1982;
 target_ptr = Libderiv->deriv2_classes[0][3][30];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+409, dvrr_stack+1262, NULL);
 tmp = dvrr_stack + 409;
 target_ptr = Libderiv->deriv2_classes[0][4][30];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+1115, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 1115;
 target_ptr = Libderiv->deriv2_classes[0][3][26];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+424, dvrr_stack+1367, NULL);
 tmp = dvrr_stack + 424;
 target_ptr = Libderiv->deriv2_classes[0][4][26];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,10,dvrr_stack+439, dvrr_stack+899, NULL);
 tmp = dvrr_stack + 439;
 target_ptr = Libderiv->deriv2_classes[0][3][23];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,15,dvrr_stack+449, dvrr_stack+364, NULL);
 tmp = dvrr_stack + 449;
 target_ptr = Libderiv->deriv2_classes[0][4][23];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+464, dvrr_stack+1127, NULL);
 tmp = dvrr_stack + 464;
 target_ptr = Libderiv->deriv2_classes[0][3][22];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+474, dvrr_stack+1157, NULL);
 tmp = dvrr_stack + 474;
 target_ptr = Libderiv->deriv2_classes[0][4][22];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+489, dvrr_stack+1202, NULL);
 tmp = dvrr_stack + 489;
 target_ptr = Libderiv->deriv2_classes[0][3][21];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+499, dvrr_stack+809, NULL);
 tmp = dvrr_stack + 499;
 target_ptr = Libderiv->deriv2_classes[0][4][21];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+514, dvrr_stack+1010, NULL);
 tmp = dvrr_stack + 514;
 target_ptr = Libderiv->deriv2_classes[0][3][20];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+524, dvrr_stack+854, NULL);
 tmp = dvrr_stack + 524;
 target_ptr = Libderiv->deriv2_classes[0][4][20];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+539, dvrr_stack+1040, NULL);
 tmp = dvrr_stack + 539;
 target_ptr = Libderiv->deriv2_classes[0][3][19];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+549, dvrr_stack+1070, NULL);
 tmp = dvrr_stack + 549;
 target_ptr = Libderiv->deriv2_classes[0][4][19];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+564, dvrr_stack+779, NULL);
 tmp = dvrr_stack + 564;
 target_ptr = Libderiv->deriv_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+574, dvrr_stack+1232, NULL);
 tmp = dvrr_stack + 574;
 target_ptr = Libderiv->deriv2_classes[0][3][18];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+584, dvrr_stack+1262, NULL);
 tmp = dvrr_stack + 584;
 target_ptr = Libderiv->deriv2_classes[0][4][18];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+599, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 599;
 target_ptr = Libderiv->deriv2_classes[0][3][14];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+609, dvrr_stack+1367, NULL);
 tmp = dvrr_stack + 609;
 target_ptr = Libderiv->deriv2_classes[0][4][14];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+624, dvrr_stack+929, NULL);
 tmp = dvrr_stack + 624;
 target_ptr = Libderiv->deriv2_classes[0][3][13];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+634, dvrr_stack+959, NULL);
 tmp = dvrr_stack + 634;
 target_ptr = Libderiv->deriv2_classes[0][4][13];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,10,dvrr_stack+649, dvrr_stack+899, NULL);
 tmp = dvrr_stack + 649;
 target_ptr = Libderiv->deriv2_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,15,dvrr_stack+899, dvrr_stack+364, NULL);
 tmp = dvrr_stack + 899;
 target_ptr = Libderiv->deriv2_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+914, dvrr_stack+1127, NULL);
 tmp = dvrr_stack + 914;
 target_ptr = Libderiv->deriv2_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+360, dvrr_stack+1157, NULL);
 tmp = dvrr_stack + 360;
 target_ptr = Libderiv->deriv2_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+375, dvrr_stack+1202, NULL);
 tmp = dvrr_stack + 375;
 target_ptr = Libderiv->deriv2_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+385, dvrr_stack+809, NULL);
 tmp = dvrr_stack + 385;
 target_ptr = Libderiv->deriv2_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+809, dvrr_stack+1010, NULL);
 tmp = dvrr_stack + 809;
 target_ptr = Libderiv->deriv2_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+819, dvrr_stack+854, NULL);
 tmp = dvrr_stack + 819;
 target_ptr = Libderiv->deriv2_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+834, dvrr_stack+1040, NULL);
 tmp = dvrr_stack + 834;
 target_ptr = Libderiv->deriv2_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+844, dvrr_stack+1070, NULL);
 tmp = dvrr_stack + 844;
 target_ptr = Libderiv->deriv2_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+859, dvrr_stack+779, NULL);
 tmp = dvrr_stack + 859;
 target_ptr = Libderiv->deriv_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+869, dvrr_stack+1232, NULL);
 tmp = dvrr_stack + 869;
 target_ptr = Libderiv->deriv2_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+879, dvrr_stack+1262, NULL);
 tmp = dvrr_stack + 879;
 target_ptr = Libderiv->deriv2_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1004, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 1004;
 target_ptr = Libderiv->deriv2_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+0, dvrr_stack+1367, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv2_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1367, dvrr_stack+929, NULL);
 tmp = dvrr_stack + 1367;
 target_ptr = Libderiv->deriv2_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+1377, dvrr_stack+959, NULL);
 tmp = dvrr_stack + 1377;
 target_ptr = Libderiv->deriv2_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+1392, dvrr_stack+37, NULL);
 tmp = dvrr_stack + 1392;
 target_ptr = Libderiv->deriv2_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+15, dvrr_stack+1307, NULL);
 tmp = dvrr_stack + 15;
 target_ptr = Libderiv->deriv2_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];


}

