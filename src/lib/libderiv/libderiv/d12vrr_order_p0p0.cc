#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (p0|p0) integrals */

void d12vrr_order_p0p0(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (1 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, NULL, NULL, Data->F+1);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+15, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+18, dvrr_stack+0, dvrr_stack+15, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+24, dvrr_stack+3, dvrr_stack+0, Data->F+0, Data->F+1, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+30, dvrr_stack+24, dvrr_stack+18, NULL, NULL, dvrr_stack+0);

 /* compute (1 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+48,dvrr_stack+30,dvrr_stack+6,3);


 /* compute (1 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+75, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (1 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+78, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+81, dvrr_stack+0, dvrr_stack+15, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+90, dvrr_stack+6, dvrr_stack+81, dvrr_stack+3, dvrr_stack+0, dvrr_stack+78);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+108, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+111, dvrr_stack+15, dvrr_stack+108, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+117, dvrr_stack+18, dvrr_stack+111, dvrr_stack+0, dvrr_stack+15, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+127, dvrr_stack+24, dvrr_stack+18, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+137, dvrr_stack+127, dvrr_stack+117, NULL, NULL, dvrr_stack+18);

 /* compute (1 0 | 2 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_dp(Libderiv->CD,dvrr_stack+167,dvrr_stack+137,dvrr_stack+30,3);


 /* compute (1 0 | 1 2) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pd(Libderiv->CD,dvrr_stack+221,dvrr_stack+167,dvrr_stack+48,3);


 /* compute (1 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_p(Data,9,dvrr_stack+275, dvrr_stack+221, dvrr_stack+6);

 /* compute (1 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_p(Data,9,dvrr_stack+302, dvrr_stack+221, dvrr_stack+6);

 /* compute (1 0 | 1 1) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_p(Data,9,dvrr_stack+329, dvrr_stack+221, dvrr_stack+6);

 /* compute (1 0 | 0 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_0p(Libderiv->CD,dvrr_stack+221,dvrr_stack+6,dvrr_stack+75,3);


 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+230, dvrr_stack+221, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+233, dvrr_stack+167, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+251, dvrr_stack+221, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+254, dvrr_stack+167, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+272, dvrr_stack+221, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+117, dvrr_stack+167, NULL);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_0(Data,3,1,dvrr_stack+167, dvrr_stack+6, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_d(Data,3,1,dvrr_stack+170, dvrr_stack+137, dvrr_stack+6);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_0(Data,3,1,dvrr_stack+188, dvrr_stack+6, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_d(Data,3,1,dvrr_stack+191, dvrr_stack+137, dvrr_stack+6);

 /* compute (1 0 | 0 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_0(Data,3,1,dvrr_stack+209, dvrr_stack+6, NULL);

 /* compute (1 0 | 2 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_d(Data,3,1,dvrr_stack+212, dvrr_stack+137, dvrr_stack+6);

 /* compute (0 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+135,dvrr_stack+24,dvrr_stack+3,1);


 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,3,dvrr_stack+144, dvrr_stack+135, NULL);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+147, dvrr_stack+18, dvrr_stack+111, NULL, NULL, dvrr_stack+15);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+356, dvrr_stack+30, dvrr_stack+147, dvrr_stack+24, dvrr_stack+18, dvrr_stack+81);

 /* compute (2 0 | 1 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_pp(Libderiv->CD,dvrr_stack+392,dvrr_stack+356,dvrr_stack+90,6);


 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,18,dvrr_stack+147, dvrr_stack+392, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,3,dvrr_stack+18, dvrr_stack+135, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,18,dvrr_stack+446, dvrr_stack+392, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,3,dvrr_stack+21, dvrr_stack+135, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,18,dvrr_stack+464, dvrr_stack+392, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,1,1,dvrr_stack+392, dvrr_stack+24, Data->F+0);

 /* compute (2 0 | 0 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+111, dvrr_stack+75, dvrr_stack+78, Data->F+0, Data->F+1, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,6,1,dvrr_stack+395, dvrr_stack+356, dvrr_stack+111);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,1,1,dvrr_stack+413, dvrr_stack+24, Data->F+0);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,6,1,dvrr_stack+416, dvrr_stack+356, dvrr_stack+111);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,1,1,dvrr_stack+434, dvrr_stack+24, Data->F+0);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,6,1,dvrr_stack+482, dvrr_stack+356, dvrr_stack+111);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,3,dvrr_stack+111, dvrr_stack+6, NULL);

 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+114, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+24, dvrr_stack+78, dvrr_stack+114, Data->F+1, Data->F+2, NULL);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+437, dvrr_stack+15, dvrr_stack+108, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+356, dvrr_stack+81, dvrr_stack+437, dvrr_stack+0, dvrr_stack+15, dvrr_stack+114);

 /* compute (3 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+500, dvrr_stack+90, dvrr_stack+356, dvrr_stack+6, dvrr_stack+81, dvrr_stack+24);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_d(Data,3,dvrr_stack+356, dvrr_stack+500, dvrr_stack+6);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,3,dvrr_stack+114, dvrr_stack+6, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_d(Data,3,dvrr_stack+374, dvrr_stack+500, dvrr_stack+6);

 /* compute (0 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,3,dvrr_stack+15, dvrr_stack+6, NULL);

 /* compute (2 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_d(Data,3,dvrr_stack+530, dvrr_stack+500, dvrr_stack+6);

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,9,dvrr_stack+6, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 6;
 target_ptr = Libderiv->deriv_classes[1][1][11];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+437, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 437;
 target_ptr = Libderiv->deriv_classes[1][1][10];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+135, dvrr_stack+48, NULL);
 tmp = dvrr_stack + 135;
 target_ptr = Libderiv->deriv_classes[1][1][9];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+48, dvrr_stack+30, dvrr_stack+75);
 tmp = dvrr_stack + 48;
 target_ptr = Libderiv->deriv_classes[1][1][8];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+57, dvrr_stack+30, dvrr_stack+75);
 tmp = dvrr_stack + 57;
 target_ptr = Libderiv->deriv_classes[1][1][7];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+66, dvrr_stack+30, dvrr_stack+75);
 tmp = dvrr_stack + 66;
 target_ptr = Libderiv->deriv_classes[1][1][6];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+75, dvrr_stack+90, dvrr_stack+3);
 tmp = dvrr_stack + 75;
 target_ptr = Libderiv->deriv_classes[1][1][2];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+500, dvrr_stack+90, dvrr_stack+3);
 tmp = dvrr_stack + 500;
 target_ptr = Libderiv->deriv_classes[1][1][1];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+509, dvrr_stack+90, dvrr_stack+3);
 tmp = dvrr_stack + 509;
 target_ptr = Libderiv->deriv_classes[1][1][0];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 2 */
 deriv_build_DZ_0(Data,9,dvrr_stack+84, dvrr_stack+275, NULL);
 tmp = dvrr_stack + 84;
 target_ptr = Libderiv->deriv2_classes[1][1][143];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 1 */
 deriv_build_DY_0(Data,9,dvrr_stack+93, dvrr_stack+275, NULL);
 tmp = dvrr_stack + 93;
 target_ptr = Libderiv->deriv2_classes[1][1][131];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 2 0 */
 deriv_build_DY_0(Data,9,dvrr_stack+102, dvrr_stack+302, NULL);
 tmp = dvrr_stack + 102;
 target_ptr = Libderiv->deriv2_classes[1][1][130];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 1 */
 deriv_build_DX_0(Data,9,dvrr_stack+518, dvrr_stack+275, NULL);
 tmp = dvrr_stack + 518;
 target_ptr = Libderiv->deriv2_classes[1][1][119];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 1 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+275, dvrr_stack+302, NULL);
 tmp = dvrr_stack + 275;
 target_ptr = Libderiv->deriv2_classes[1][1][118];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  2 0 0 */
 deriv_build_DX_0(Data,9,dvrr_stack+284, dvrr_stack+329, NULL);
 tmp = dvrr_stack + 284;
 target_ptr = Libderiv->deriv2_classes[1][1][117];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 1 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+293, dvrr_stack+233, dvrr_stack+230);
 tmp = dvrr_stack + 293;
 target_ptr = Libderiv->deriv2_classes[1][1][107];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 1 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+302, dvrr_stack+254, dvrr_stack+251);
 tmp = dvrr_stack + 302;
 target_ptr = Libderiv->deriv2_classes[1][1][106];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  1 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+311, dvrr_stack+117, dvrr_stack+272);
 tmp = dvrr_stack + 311;
 target_ptr = Libderiv->deriv2_classes[1][1][105];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 2  0 0 0 */
 deriv_build_CZ_p(Data,3,1,dvrr_stack+320, dvrr_stack+170, dvrr_stack+167);
 tmp = dvrr_stack + 320;
 target_ptr = Libderiv->deriv2_classes[1][1][104];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 1 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+329, dvrr_stack+233, dvrr_stack+230);
 tmp = dvrr_stack + 329;
 target_ptr = Libderiv->deriv2_classes[1][1][95];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 1 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+338, dvrr_stack+254, dvrr_stack+251);
 tmp = dvrr_stack + 338;
 target_ptr = Libderiv->deriv2_classes[1][1][94];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  1 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+347, dvrr_stack+117, dvrr_stack+272);
 tmp = dvrr_stack + 347;
 target_ptr = Libderiv->deriv2_classes[1][1][93];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 1  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+24, dvrr_stack+170, dvrr_stack+167);
 tmp = dvrr_stack + 24;
 target_ptr = Libderiv->deriv2_classes[1][1][92];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  0 2 0  0 0 0 */
 deriv_build_CY_p(Data,3,1,dvrr_stack+33, dvrr_stack+191, dvrr_stack+188);
 tmp = dvrr_stack + 33;
 target_ptr = Libderiv->deriv2_classes[1][1][91];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 1 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+548, dvrr_stack+233, dvrr_stack+230);
 tmp = dvrr_stack + 548;
 target_ptr = Libderiv->deriv2_classes[1][1][83];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 1 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+230, dvrr_stack+254, dvrr_stack+251);
 tmp = dvrr_stack + 230;
 target_ptr = Libderiv->deriv2_classes[1][1][82];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  1 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+239, dvrr_stack+117, dvrr_stack+272);
 tmp = dvrr_stack + 239;
 target_ptr = Libderiv->deriv2_classes[1][1][81];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 1  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+117, dvrr_stack+170, dvrr_stack+167);
 tmp = dvrr_stack + 117;
 target_ptr = Libderiv->deriv2_classes[1][1][80];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  1 1 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+126, dvrr_stack+191, dvrr_stack+188);
 tmp = dvrr_stack + 126;
 target_ptr = Libderiv->deriv2_classes[1][1][79];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 0  0 0 0  2 0 0  0 0 0 */
 deriv_build_CX_p(Data,3,1,dvrr_stack+248, dvrr_stack+212, dvrr_stack+209);
 tmp = dvrr_stack + 248;
 target_ptr = Libderiv->deriv2_classes[1][1][78];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 1 */
 deriv_build_AZ_p(Data,3,dvrr_stack+257, dvrr_stack+147, dvrr_stack+144);
 tmp = dvrr_stack + 257;
 target_ptr = Libderiv->deriv2_classes[1][1][35];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 1 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+266, dvrr_stack+446, dvrr_stack+18);
 tmp = dvrr_stack + 266;
 target_ptr = Libderiv->deriv2_classes[1][1][34];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  1 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+165, dvrr_stack+464, dvrr_stack+21);
 tmp = dvrr_stack + 165;
 target_ptr = Libderiv->deriv2_classes[1][1][33];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 1  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+174, dvrr_stack+395, dvrr_stack+392);
 tmp = dvrr_stack + 174;
 target_ptr = Libderiv->deriv2_classes[1][1][32];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  0 1 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+183, dvrr_stack+416, dvrr_stack+413);
 tmp = dvrr_stack + 183;
 target_ptr = Libderiv->deriv2_classes[1][1][31];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 1  0 0 0  1 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+192, dvrr_stack+482, dvrr_stack+434);
 tmp = dvrr_stack + 192;
 target_ptr = Libderiv->deriv2_classes[1][1][30];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 0 2  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_p(Data,3,dvrr_stack+201, dvrr_stack+356, dvrr_stack+111);
 tmp = dvrr_stack + 201;
 target_ptr = Libderiv->deriv2_classes[1][1][26];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AY_p(Data,3,dvrr_stack+210, dvrr_stack+147, dvrr_stack+144);
 tmp = dvrr_stack + 210;
 target_ptr = Libderiv->deriv2_classes[1][1][23];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+219, dvrr_stack+446, dvrr_stack+18);
 tmp = dvrr_stack + 219;
 target_ptr = Libderiv->deriv2_classes[1][1][22];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+557, dvrr_stack+464, dvrr_stack+21);
 tmp = dvrr_stack + 557;
 target_ptr = Libderiv->deriv2_classes[1][1][21];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+566, dvrr_stack+395, dvrr_stack+392);
 tmp = dvrr_stack + 566;
 target_ptr = Libderiv->deriv2_classes[1][1][20];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+575, dvrr_stack+416, dvrr_stack+413);
 tmp = dvrr_stack + 575;
 target_ptr = Libderiv->deriv2_classes[1][1][19];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+584, dvrr_stack+482, dvrr_stack+434);
 tmp = dvrr_stack + 584;
 target_ptr = Libderiv->deriv2_classes[1][1][18];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 1 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+593, dvrr_stack+356, dvrr_stack+111);
 tmp = dvrr_stack + 593;
 target_ptr = Libderiv->deriv2_classes[1][1][14];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 0 2 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_p(Data,3,dvrr_stack+602, dvrr_stack+374, dvrr_stack+114);
 tmp = dvrr_stack + 602;
 target_ptr = Libderiv->deriv2_classes[1][1][13];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_AX_p(Data,3,dvrr_stack+611, dvrr_stack+147, dvrr_stack+144);
 tmp = dvrr_stack + 611;
 target_ptr = Libderiv->deriv2_classes[1][1][11];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+144, dvrr_stack+446, dvrr_stack+18);
 tmp = dvrr_stack + 144;
 target_ptr = Libderiv->deriv2_classes[1][1][10];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+446, dvrr_stack+464, dvrr_stack+21);
 tmp = dvrr_stack + 446;
 target_ptr = Libderiv->deriv2_classes[1][1][9];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+455, dvrr_stack+395, dvrr_stack+392);
 tmp = dvrr_stack + 455;
 target_ptr = Libderiv->deriv2_classes[1][1][8];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+392, dvrr_stack+416, dvrr_stack+413);
 tmp = dvrr_stack + 392;
 target_ptr = Libderiv->deriv2_classes[1][1][7];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+401, dvrr_stack+482, dvrr_stack+434);
 tmp = dvrr_stack + 401;
 target_ptr = Libderiv->deriv2_classes[1][1][6];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+410, dvrr_stack+356, dvrr_stack+111);
 tmp = dvrr_stack + 410;
 target_ptr = Libderiv->deriv2_classes[1][1][2];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 1 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+356, dvrr_stack+374, dvrr_stack+114);
 tmp = dvrr_stack + 356;
 target_ptr = Libderiv->deriv2_classes[1][1][1];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 1 0) m=0 deriv level 2 */
 /* deriv_ind: 2 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_p(Data,3,dvrr_stack+365, dvrr_stack+530, dvrr_stack+15);
 tmp = dvrr_stack + 365;
 target_ptr = Libderiv->deriv2_classes[1][1][0];
 for(i=0;i<9;i++)
   target_ptr[i] += tmp[i];


}

