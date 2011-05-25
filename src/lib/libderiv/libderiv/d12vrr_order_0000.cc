#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|00) integrals */

void d12vrr_order_0000(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 0 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0p(Libderiv->CD,dvrr_stack+3,dvrr_stack+0,Data->F,1);


 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+6, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+9, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+12, dvrr_stack+0, dvrr_stack+9, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+18,dvrr_stack+12,dvrr_stack+0,1);


 /* compute (0 0 | 0 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0d(Libderiv->CD,dvrr_stack+27,dvrr_stack+18,dvrr_stack+3,1);


 /* compute (0 0 | 0 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,1,dvrr_stack+33, dvrr_stack+27, Data->F+0);

 /* compute (0 0 | 0 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,1,dvrr_stack+36, dvrr_stack+27, Data->F+0);

 /* compute (0 0 | 0 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,1,dvrr_stack+39, dvrr_stack+27, Data->F+0);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+27, dvrr_stack+18, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+30, dvrr_stack+18, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+42, dvrr_stack+18, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+18, dvrr_stack+12, Data->F+0);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+21, dvrr_stack+12, Data->F+0);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+24, dvrr_stack+12, Data->F+0);

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+45, dvrr_stack+0, dvrr_stack+9, NULL, NULL, Data->F+1);

 /* compute (1 0 | 0 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0p(Libderiv->CD,dvrr_stack+9,dvrr_stack+45,dvrr_stack+6,3);


 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+54, dvrr_stack+9, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+57, dvrr_stack+9, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+60, dvrr_stack+9, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_0(Data,3,1,dvrr_stack+9, dvrr_stack+45, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_0(Data,3,1,dvrr_stack+12, dvrr_stack+45, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_0(Data,3,1,dvrr_stack+15, dvrr_stack+45, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+45, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+48, dvrr_stack+6, dvrr_stack+45, Data->F+0, Data->F+1, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,1,dvrr_stack+45, dvrr_stack+48, Data->F+0);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,1,dvrr_stack+63, dvrr_stack+48, Data->F+0);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,1,dvrr_stack+66, dvrr_stack+48, Data->F+0);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1,dvrr_stack+48, dvrr_stack+3, NULL);
 tmp = dvrr_stack + 48;
 target_ptr = Libderiv->deriv_classes[0][0][11];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1,dvrr_stack+49, dvrr_stack+3, NULL);
 tmp = dvrr_stack + 49;
 target_ptr = Libderiv->deriv_classes[0][0][10];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1,dvrr_stack+50, dvrr_stack+3, NULL);
 tmp = dvrr_stack + 50;
 target_ptr = Libderiv->deriv_classes[0][0][9];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_0(Data,1,1,dvrr_stack+3, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 3;
 target_ptr = Libderiv->deriv_classes[0][0][8];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+4, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 4;
 target_ptr = Libderiv->deriv_classes[0][0][7];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+5, dvrr_stack+0, NULL);
 tmp = dvrr_stack + 5;
 target_ptr = Libderiv->deriv_classes[0][0][6];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,1,dvrr_stack+0, dvrr_stack+6, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[0][0][2];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+1, dvrr_stack+6, NULL);
 tmp = dvrr_stack + 1;
 target_ptr = Libderiv->deriv_classes[0][0][1];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+2, dvrr_stack+6, NULL);
 tmp = dvrr_stack + 2;
 target_ptr = Libderiv->deriv_classes[0][0][0];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,1,dvrr_stack+6, dvrr_stack+33, NULL);
 tmp = dvrr_stack + 6;
 target_ptr = Libderiv->deriv2_classes[0][0][143];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,1,dvrr_stack+7, dvrr_stack+33, NULL);
 tmp = dvrr_stack + 7;
 target_ptr = Libderiv->deriv2_classes[0][0][131];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,1,dvrr_stack+8, dvrr_stack+36, NULL);
 tmp = dvrr_stack + 8;
 target_ptr = Libderiv->deriv2_classes[0][0][130];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,1,dvrr_stack+51, dvrr_stack+33, NULL);
 tmp = dvrr_stack + 51;
 target_ptr = Libderiv->deriv2_classes[0][0][119];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,1,dvrr_stack+33, dvrr_stack+36, NULL);
 tmp = dvrr_stack + 33;
 target_ptr = Libderiv->deriv2_classes[0][0][118];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,1,dvrr_stack+34, dvrr_stack+39, NULL);
 tmp = dvrr_stack + 34;
 target_ptr = Libderiv->deriv2_classes[0][0][117];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_0(Data,1,1,dvrr_stack+35, dvrr_stack+27, NULL);
 tmp = dvrr_stack + 35;
 target_ptr = Libderiv->deriv2_classes[0][0][107];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_0(Data,1,1,dvrr_stack+36, dvrr_stack+30, NULL);
 tmp = dvrr_stack + 36;
 target_ptr = Libderiv->deriv2_classes[0][0][106];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_0(Data,1,1,dvrr_stack+37, dvrr_stack+42, NULL);
 tmp = dvrr_stack + 37;
 target_ptr = Libderiv->deriv2_classes[0][0][105];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_0(Data,1,1,dvrr_stack+38, dvrr_stack+18, NULL);
 tmp = dvrr_stack + 38;
 target_ptr = Libderiv->deriv2_classes[0][0][104];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+39, dvrr_stack+27, NULL);
 tmp = dvrr_stack + 39;
 target_ptr = Libderiv->deriv2_classes[0][0][95];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+40, dvrr_stack+30, NULL);
 tmp = dvrr_stack + 40;
 target_ptr = Libderiv->deriv2_classes[0][0][94];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+41, dvrr_stack+42, NULL);
 tmp = dvrr_stack + 41;
 target_ptr = Libderiv->deriv2_classes[0][0][93];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+52, dvrr_stack+18, NULL);
 tmp = dvrr_stack + 52;
 target_ptr = Libderiv->deriv2_classes[0][0][92];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+53, dvrr_stack+21, NULL);
 tmp = dvrr_stack + 53;
 target_ptr = Libderiv->deriv2_classes[0][0][91];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+69, dvrr_stack+27, NULL);
 tmp = dvrr_stack + 69;
 target_ptr = Libderiv->deriv2_classes[0][0][83];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+27, dvrr_stack+30, NULL);
 tmp = dvrr_stack + 27;
 target_ptr = Libderiv->deriv2_classes[0][0][82];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+28, dvrr_stack+42, NULL);
 tmp = dvrr_stack + 28;
 target_ptr = Libderiv->deriv2_classes[0][0][81];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+42, dvrr_stack+18, NULL);
 tmp = dvrr_stack + 42;
 target_ptr = Libderiv->deriv2_classes[0][0][80];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+18, dvrr_stack+21, NULL);
 tmp = dvrr_stack + 18;
 target_ptr = Libderiv->deriv2_classes[0][0][79];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+19, dvrr_stack+24, NULL);
 tmp = dvrr_stack + 19;
 target_ptr = Libderiv->deriv2_classes[0][0][78];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,1,dvrr_stack+20, dvrr_stack+54, NULL);
 tmp = dvrr_stack + 20;
 target_ptr = Libderiv->deriv2_classes[0][0][35];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,1,dvrr_stack+21, dvrr_stack+57, NULL);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv2_classes[0][0][34];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,1,dvrr_stack+22, dvrr_stack+60, NULL);
 tmp = dvrr_stack + 22;
 target_ptr = Libderiv->deriv2_classes[0][0][33];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,1,dvrr_stack+23, dvrr_stack+9, NULL);
 tmp = dvrr_stack + 23;
 target_ptr = Libderiv->deriv2_classes[0][0][32];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,1,dvrr_stack+24, dvrr_stack+12, NULL);
 tmp = dvrr_stack + 24;
 target_ptr = Libderiv->deriv2_classes[0][0][31];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,1,dvrr_stack+25, dvrr_stack+15, NULL);
 tmp = dvrr_stack + 25;
 target_ptr = Libderiv->deriv2_classes[0][0][30];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,1,dvrr_stack+26, dvrr_stack+45, NULL);
 tmp = dvrr_stack + 26;
 target_ptr = Libderiv->deriv2_classes[0][0][26];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,1,dvrr_stack+43, dvrr_stack+54, NULL);
 tmp = dvrr_stack + 43;
 target_ptr = Libderiv->deriv2_classes[0][0][23];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+44, dvrr_stack+57, NULL);
 tmp = dvrr_stack + 44;
 target_ptr = Libderiv->deriv2_classes[0][0][22];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+29, dvrr_stack+60, NULL);
 tmp = dvrr_stack + 29;
 target_ptr = Libderiv->deriv2_classes[0][0][21];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+30, dvrr_stack+9, NULL);
 tmp = dvrr_stack + 30;
 target_ptr = Libderiv->deriv2_classes[0][0][20];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+31, dvrr_stack+12, NULL);
 tmp = dvrr_stack + 31;
 target_ptr = Libderiv->deriv2_classes[0][0][19];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+32, dvrr_stack+15, NULL);
 tmp = dvrr_stack + 32;
 target_ptr = Libderiv->deriv2_classes[0][0][18];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+70, dvrr_stack+45, NULL);
 tmp = dvrr_stack + 70;
 target_ptr = Libderiv->deriv2_classes[0][0][14];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,1,dvrr_stack+71, dvrr_stack+63, NULL);
 tmp = dvrr_stack + 71;
 target_ptr = Libderiv->deriv2_classes[0][0][13];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,1,dvrr_stack+72, dvrr_stack+54, NULL);
 tmp = dvrr_stack + 72;
 target_ptr = Libderiv->deriv2_classes[0][0][11];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+54, dvrr_stack+57, NULL);
 tmp = dvrr_stack + 54;
 target_ptr = Libderiv->deriv2_classes[0][0][10];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+55, dvrr_stack+60, NULL);
 tmp = dvrr_stack + 55;
 target_ptr = Libderiv->deriv2_classes[0][0][9];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+56, dvrr_stack+9, NULL);
 tmp = dvrr_stack + 56;
 target_ptr = Libderiv->deriv2_classes[0][0][8];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+9, dvrr_stack+12, NULL);
 tmp = dvrr_stack + 9;
 target_ptr = Libderiv->deriv2_classes[0][0][7];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+10, dvrr_stack+15, NULL);
 tmp = dvrr_stack + 10;
 target_ptr = Libderiv->deriv2_classes[0][0][6];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+11, dvrr_stack+45, NULL);
 tmp = dvrr_stack + 11;
 target_ptr = Libderiv->deriv2_classes[0][0][2];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+45, dvrr_stack+63, NULL);
 tmp = dvrr_stack + 45;
 target_ptr = Libderiv->deriv2_classes[0][0][1];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 0 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,1,dvrr_stack+46, dvrr_stack+66, NULL);
 tmp = dvrr_stack + 46;
 target_ptr = Libderiv->deriv2_classes[0][0][0];
 for(i=0;i<1;i++)
   target_ptr[i] += tmp[i];


}

