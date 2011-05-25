#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|fd) integrals */

void d1vrr_order_00fd(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+9, dvrr_stack+0, dvrr_stack+6, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+3, dvrr_stack+0, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+21, dvrr_stack+15, dvrr_stack+9, dvrr_stack+3, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->dvrr_classes[0][3];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+31, dvrr_stack+6, dvrr_stack+3, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+37, dvrr_stack+9, dvrr_stack+31, dvrr_stack+0, dvrr_stack+6, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+47, dvrr_stack+21, dvrr_stack+37, dvrr_stack+15, dvrr_stack+9, NULL);
 tmp = dvrr_stack + 47;
 target_ptr = Libderiv->dvrr_classes[0][4];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+62,dvrr_stack+47,dvrr_stack+21,1);


 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+92, dvrr_stack+3, dvrr_stack+0, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+98, dvrr_stack+31, dvrr_stack+92, dvrr_stack+6, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+108, dvrr_stack+37, dvrr_stack+98, dvrr_stack+9, dvrr_stack+31, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+123, dvrr_stack+47, dvrr_stack+108, dvrr_stack+21, dvrr_stack+37, NULL);

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+144,dvrr_stack+123,dvrr_stack+47,1);


 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+189, dvrr_stack+0, dvrr_stack+6, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+195, dvrr_stack+92, dvrr_stack+189, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+205, dvrr_stack+98, dvrr_stack+195, dvrr_stack+31, dvrr_stack+92, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+220, dvrr_stack+108, dvrr_stack+205, dvrr_stack+37, dvrr_stack+98, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+189, dvrr_stack+123, dvrr_stack+220, dvrr_stack+47, dvrr_stack+108, NULL);

 /* compute (0 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+241,dvrr_stack+189,dvrr_stack+123,1);


 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+304, dvrr_stack+21, dvrr_stack+37, NULL, NULL, dvrr_stack+9);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+334, dvrr_stack+47, dvrr_stack+108, NULL, NULL, dvrr_stack+37);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+379, dvrr_stack+123, dvrr_stack+220, NULL, NULL, dvrr_stack+108);

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,10,dvrr_stack+92, dvrr_stack+62, NULL);
 tmp = dvrr_stack + 92;
 target_ptr = Libderiv->deriv_classes[0][3][11];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+0, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,21,dvrr_stack+102, dvrr_stack+241, NULL);
 tmp = dvrr_stack + 102;
 target_ptr = Libderiv->deriv_classes[0][5][11];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,10,dvrr_stack+31, dvrr_stack+62, NULL);
 tmp = dvrr_stack + 31;
 target_ptr = Libderiv->deriv_classes[0][3][10];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+217, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 217;
 target_ptr = Libderiv->deriv_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+442, dvrr_stack+241, NULL);
 tmp = dvrr_stack + 442;
 target_ptr = Libderiv->deriv_classes[0][5][10];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,10,dvrr_stack+463, dvrr_stack+62, NULL);
 tmp = dvrr_stack + 463;
 target_ptr = Libderiv->deriv_classes[0][3][9];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+62, dvrr_stack+144, NULL);
 tmp = dvrr_stack + 62;
 target_ptr = Libderiv->deriv_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+144, dvrr_stack+241, NULL);
 tmp = dvrr_stack + 144;
 target_ptr = Libderiv->deriv_classes[0][5][9];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,1,1,dvrr_stack+165, dvrr_stack+47, dvrr_stack+15);
 tmp = dvrr_stack + 165;
 target_ptr = Libderiv->deriv_classes[0][3][8];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+77, dvrr_stack+123, dvrr_stack+21);
 tmp = dvrr_stack + 77;
 target_ptr = Libderiv->deriv_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+232, dvrr_stack+189, dvrr_stack+47);
 tmp = dvrr_stack + 232;
 target_ptr = Libderiv->deriv_classes[0][5][8];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,1,1,dvrr_stack+175, dvrr_stack+47, dvrr_stack+15);
 tmp = dvrr_stack + 175;
 target_ptr = Libderiv->deriv_classes[0][3][7];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+253, dvrr_stack+123, dvrr_stack+21);
 tmp = dvrr_stack + 253;
 target_ptr = Libderiv->deriv_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+268, dvrr_stack+189, dvrr_stack+47);
 tmp = dvrr_stack + 268;
 target_ptr = Libderiv->deriv_classes[0][5][7];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,1,1,dvrr_stack+289, dvrr_stack+47, dvrr_stack+15);
 tmp = dvrr_stack + 289;
 target_ptr = Libderiv->deriv_classes[0][3][6];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+473, dvrr_stack+123, dvrr_stack+21);
 tmp = dvrr_stack + 473;
 target_ptr = Libderiv->deriv_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+123, dvrr_stack+189, dvrr_stack+47);
 tmp = dvrr_stack + 123;
 target_ptr = Libderiv->deriv_classes[0][5][6];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,10,dvrr_stack+15, dvrr_stack+304, NULL);
 tmp = dvrr_stack + 15;
 target_ptr = Libderiv->deriv_classes[0][3][2];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+185, dvrr_stack+334, NULL);
 tmp = dvrr_stack + 185;
 target_ptr = Libderiv->deriv_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+41, dvrr_stack+379, NULL);
 tmp = dvrr_stack + 41;
 target_ptr = Libderiv->deriv_classes[0][5][2];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,10,dvrr_stack+200, dvrr_stack+304, NULL);
 tmp = dvrr_stack + 200;
 target_ptr = Libderiv->deriv_classes[0][3][1];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+488, dvrr_stack+334, NULL);
 tmp = dvrr_stack + 488;
 target_ptr = Libderiv->deriv_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+503, dvrr_stack+379, NULL);
 tmp = dvrr_stack + 503;
 target_ptr = Libderiv->deriv_classes[0][5][1];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,10,dvrr_stack+524, dvrr_stack+304, NULL);
 tmp = dvrr_stack + 524;
 target_ptr = Libderiv->deriv_classes[0][3][0];
 for(i=0;i<10;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+299, dvrr_stack+334, NULL);
 tmp = dvrr_stack + 299;
 target_ptr = Libderiv->deriv_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+314, dvrr_stack+379, NULL);
 tmp = dvrr_stack + 314;
 target_ptr = Libderiv->deriv_classes[0][5][0];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];


}

