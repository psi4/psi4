#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|p0) integrals */

void d12vrr_order_00p0(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+0, dvrr_stack+3, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+12,dvrr_stack+6,dvrr_stack+0,1);


 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+21, dvrr_stack+0, dvrr_stack+3, NULL, NULL, Data->F+1);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+30, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+33, dvrr_stack+3, dvrr_stack+30, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+39, dvrr_stack+6, dvrr_stack+33, dvrr_stack+0, dvrr_stack+3, NULL);

 /* compute (0 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+49,dvrr_stack+39,dvrr_stack+6,1);


 /* compute (0 0 | 1 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pd(Libderiv->CD,dvrr_stack+67,dvrr_stack+49,dvrr_stack+12,1);


 /* compute (0 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,3,dvrr_stack+85, dvrr_stack+67, dvrr_stack+0);

 /* compute (0 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,3,dvrr_stack+94, dvrr_stack+67, dvrr_stack+0);

 /* compute (0 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,3,dvrr_stack+103, dvrr_stack+67, dvrr_stack+0);

 /* compute (0 0 | 0 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0p(Libderiv->CD,dvrr_stack+67,dvrr_stack+0,Data->F,1);


 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,1,dvrr_stack+70, dvrr_stack+67, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,6,dvrr_stack+71, dvrr_stack+49, NULL);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,1,dvrr_stack+77, dvrr_stack+67, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,6,dvrr_stack+78, dvrr_stack+49, NULL);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,1,dvrr_stack+84, dvrr_stack+67, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,6,dvrr_stack+112, dvrr_stack+49, NULL);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_0(Data,1,1,dvrr_stack+49, dvrr_stack+0, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,1,1,dvrr_stack+50, dvrr_stack+39, dvrr_stack+0);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_0(Data,1,1,dvrr_stack+56, dvrr_stack+0, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,1,1,dvrr_stack+57, dvrr_stack+39, dvrr_stack+0);

 /* compute (0 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_0(Data,1,1,dvrr_stack+63, dvrr_stack+0, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,1,1,dvrr_stack+64, dvrr_stack+39, dvrr_stack+0);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+118, dvrr_stack+6, dvrr_stack+33, NULL, NULL, dvrr_stack+3);

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+136,dvrr_stack+118,dvrr_stack+21,3);


 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+33, dvrr_stack+136, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+163, dvrr_stack+136, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+172, dvrr_stack+136, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+136, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+139, dvrr_stack+118, dvrr_stack+136);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+148, dvrr_stack+118, dvrr_stack+136);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+181, dvrr_stack+118, dvrr_stack+136);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+118, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+121, dvrr_stack+3, dvrr_stack+30, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+190, dvrr_stack+21, dvrr_stack+121, dvrr_stack+0, dvrr_stack+3, dvrr_stack+118);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+118, dvrr_stack+190, dvrr_stack+0);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+127, dvrr_stack+190, dvrr_stack+0);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+208, dvrr_stack+190, dvrr_stack+0);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+136, dvrr_stack+12, NULL);
 tmp = dvrr_stack + 136;
 target_ptr = Libderiv->deriv_classes[0][1][11];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+30, dvrr_stack+12, NULL);
 tmp = dvrr_stack + 30;
 target_ptr = Libderiv->deriv_classes[0][1][10];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+0, dvrr_stack+12, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[0][1][9];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+3, dvrr_stack+6, Data->F+0);
 tmp = dvrr_stack + 3;
 target_ptr = Libderiv->deriv_classes[0][1][8];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+12, dvrr_stack+6, Data->F+0);
 tmp = dvrr_stack + 12;
 target_ptr = Libderiv->deriv_classes[0][1][7];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+15, dvrr_stack+6, Data->F+0);
 tmp = dvrr_stack + 15;
 target_ptr = Libderiv->deriv_classes[0][1][6];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+18, dvrr_stack+21, NULL);
 tmp = dvrr_stack + 18;
 target_ptr = Libderiv->deriv_classes[0][1][2];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+6, dvrr_stack+21, NULL);
 tmp = dvrr_stack + 6;
 target_ptr = Libderiv->deriv_classes[0][1][1];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+9, dvrr_stack+21, NULL);
 tmp = dvrr_stack + 9;
 target_ptr = Libderiv->deriv_classes[0][1][0];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,3,dvrr_stack+21, dvrr_stack+85, NULL);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv2_classes[0][1][143];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,3,dvrr_stack+24, dvrr_stack+85, NULL);
 tmp = dvrr_stack + 24;
 target_ptr = Libderiv->deriv2_classes[0][1][131];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+27, dvrr_stack+94, NULL);
 tmp = dvrr_stack + 27;
 target_ptr = Libderiv->deriv2_classes[0][1][130];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,3,dvrr_stack+190, dvrr_stack+85, NULL);
 tmp = dvrr_stack + 190;
 target_ptr = Libderiv->deriv2_classes[0][1][119];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+85, dvrr_stack+94, NULL);
 tmp = dvrr_stack + 85;
 target_ptr = Libderiv->deriv2_classes[0][1][118];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+88, dvrr_stack+103, NULL);
 tmp = dvrr_stack + 88;
 target_ptr = Libderiv->deriv2_classes[0][1][117];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+91, dvrr_stack+71, dvrr_stack+70);
 tmp = dvrr_stack + 91;
 target_ptr = Libderiv->deriv2_classes[0][1][107];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+94, dvrr_stack+78, dvrr_stack+77);
 tmp = dvrr_stack + 94;
 target_ptr = Libderiv->deriv2_classes[0][1][106];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+97, dvrr_stack+112, dvrr_stack+84);
 tmp = dvrr_stack + 97;
 target_ptr = Libderiv->deriv2_classes[0][1][105];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+100, dvrr_stack+50, dvrr_stack+49);
 tmp = dvrr_stack + 100;
 target_ptr = Libderiv->deriv2_classes[0][1][104];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+103, dvrr_stack+71, dvrr_stack+70);
 tmp = dvrr_stack + 103;
 target_ptr = Libderiv->deriv2_classes[0][1][95];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+106, dvrr_stack+78, dvrr_stack+77);
 tmp = dvrr_stack + 106;
 target_ptr = Libderiv->deriv2_classes[0][1][94];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+109, dvrr_stack+112, dvrr_stack+84);
 tmp = dvrr_stack + 109;
 target_ptr = Libderiv->deriv2_classes[0][1][93];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+193, dvrr_stack+50, dvrr_stack+49);
 tmp = dvrr_stack + 193;
 target_ptr = Libderiv->deriv2_classes[0][1][92];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+196, dvrr_stack+57, dvrr_stack+56);
 tmp = dvrr_stack + 196;
 target_ptr = Libderiv->deriv2_classes[0][1][91];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+199, dvrr_stack+71, dvrr_stack+70);
 tmp = dvrr_stack + 199;
 target_ptr = Libderiv->deriv2_classes[0][1][83];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+70, dvrr_stack+78, dvrr_stack+77);
 tmp = dvrr_stack + 70;
 target_ptr = Libderiv->deriv2_classes[0][1][82];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+73, dvrr_stack+112, dvrr_stack+84);
 tmp = dvrr_stack + 73;
 target_ptr = Libderiv->deriv2_classes[0][1][81];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+112, dvrr_stack+50, dvrr_stack+49);
 tmp = dvrr_stack + 112;
 target_ptr = Libderiv->deriv2_classes[0][1][80];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+115, dvrr_stack+57, dvrr_stack+56);
 tmp = dvrr_stack + 115;
 target_ptr = Libderiv->deriv2_classes[0][1][79];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+76, dvrr_stack+64, dvrr_stack+63);
 tmp = dvrr_stack + 76;
 target_ptr = Libderiv->deriv2_classes[0][1][78];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_0(Data,3,dvrr_stack+79, dvrr_stack+33, NULL);
 tmp = dvrr_stack + 79;
 target_ptr = Libderiv->deriv2_classes[0][1][35];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+82, dvrr_stack+163, NULL);
 tmp = dvrr_stack + 82;
 target_ptr = Libderiv->deriv2_classes[0][1][34];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+202, dvrr_stack+172, NULL);
 tmp = dvrr_stack + 202;
 target_ptr = Libderiv->deriv2_classes[0][1][33];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+205, dvrr_stack+139, NULL);
 tmp = dvrr_stack + 205;
 target_ptr = Libderiv->deriv2_classes[0][1][32];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+157, dvrr_stack+148, NULL);
 tmp = dvrr_stack + 157;
 target_ptr = Libderiv->deriv2_classes[0][1][31];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+160, dvrr_stack+181, NULL);
 tmp = dvrr_stack + 160;
 target_ptr = Libderiv->deriv2_classes[0][1][30];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+42, dvrr_stack+118, NULL);
 tmp = dvrr_stack + 42;
 target_ptr = Libderiv->deriv2_classes[0][1][26];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_0(Data,3,dvrr_stack+45, dvrr_stack+33, NULL);
 tmp = dvrr_stack + 45;
 target_ptr = Libderiv->deriv2_classes[0][1][23];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+48, dvrr_stack+163, NULL);
 tmp = dvrr_stack + 48;
 target_ptr = Libderiv->deriv2_classes[0][1][22];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+51, dvrr_stack+172, NULL);
 tmp = dvrr_stack + 51;
 target_ptr = Libderiv->deriv2_classes[0][1][21];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+54, dvrr_stack+139, NULL);
 tmp = dvrr_stack + 54;
 target_ptr = Libderiv->deriv2_classes[0][1][20];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+57, dvrr_stack+148, NULL);
 tmp = dvrr_stack + 57;
 target_ptr = Libderiv->deriv2_classes[0][1][19];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+60, dvrr_stack+181, NULL);
 tmp = dvrr_stack + 60;
 target_ptr = Libderiv->deriv2_classes[0][1][18];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+63, dvrr_stack+118, NULL);
 tmp = dvrr_stack + 63;
 target_ptr = Libderiv->deriv2_classes[0][1][14];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+66, dvrr_stack+127, NULL);
 tmp = dvrr_stack + 66;
 target_ptr = Libderiv->deriv2_classes[0][1][13];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_0(Data,3,dvrr_stack+217, dvrr_stack+33, NULL);
 tmp = dvrr_stack + 217;
 target_ptr = Libderiv->deriv2_classes[0][1][11];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+33, dvrr_stack+163, NULL);
 tmp = dvrr_stack + 33;
 target_ptr = Libderiv->deriv2_classes[0][1][10];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+163, dvrr_stack+172, NULL);
 tmp = dvrr_stack + 163;
 target_ptr = Libderiv->deriv2_classes[0][1][9];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+166, dvrr_stack+139, NULL);
 tmp = dvrr_stack + 166;
 target_ptr = Libderiv->deriv2_classes[0][1][8];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+139, dvrr_stack+148, NULL);
 tmp = dvrr_stack + 139;
 target_ptr = Libderiv->deriv2_classes[0][1][7];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+142, dvrr_stack+181, NULL);
 tmp = dvrr_stack + 142;
 target_ptr = Libderiv->deriv2_classes[0][1][6];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+145, dvrr_stack+118, NULL);
 tmp = dvrr_stack + 145;
 target_ptr = Libderiv->deriv2_classes[0][1][2];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+118, dvrr_stack+127, NULL);
 tmp = dvrr_stack + 118;
 target_ptr = Libderiv->deriv2_classes[0][1][1];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+121, dvrr_stack+208, NULL);
 tmp = dvrr_stack + 121;
 target_ptr = Libderiv->deriv2_classes[0][1][0];
 for(i=0;i<3;i++)
   target_ptr[i] += tmp[i];


}

