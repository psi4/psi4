#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (00|gd) integrals */

void d1vrr_order_00gd(Libderiv_t *Libderiv, prim_data *Data)
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

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, Data->F+1, Data->F+2, NULL);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+12, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+15, dvrr_stack+12, dvrr_stack+3, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+21, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+24, dvrr_stack+0, dvrr_stack+21, Data->F+2, Data->F+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+30, dvrr_stack+6, dvrr_stack+24, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+40, dvrr_stack+15, dvrr_stack+6, dvrr_stack+12, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+50, dvrr_stack+40, dvrr_stack+30, dvrr_stack+15, dvrr_stack+6, NULL);
 tmp = dvrr_stack + 50;
 target_ptr = Libderiv->dvrr_classes[0][4];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+12, dvrr_stack+21, dvrr_stack+3, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+65, dvrr_stack+24, dvrr_stack+12, dvrr_stack+0, dvrr_stack+21, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+75, dvrr_stack+30, dvrr_stack+65, dvrr_stack+6, dvrr_stack+24, NULL);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+90, dvrr_stack+50, dvrr_stack+75, dvrr_stack+40, dvrr_stack+30, NULL);
 tmp = dvrr_stack + 90;
 target_ptr = Libderiv->dvrr_classes[0][5];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+111,dvrr_stack+90,dvrr_stack+50,1);


 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+0, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+6, dvrr_stack+3, dvrr_stack+0, Data->F+4, Data->F+5, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+156, dvrr_stack+12, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, NULL);

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+166, dvrr_stack+65, dvrr_stack+156, dvrr_stack+24, dvrr_stack+12, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+181, dvrr_stack+75, dvrr_stack+166, dvrr_stack+30, dvrr_stack+65, NULL);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+202, dvrr_stack+90, dvrr_stack+181, dvrr_stack+50, dvrr_stack+75, NULL);

 /* compute (0 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+230,dvrr_stack+202,dvrr_stack+90,1);


 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+18, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+21, dvrr_stack+0, dvrr_stack+18, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+293, dvrr_stack+6, dvrr_stack+21, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+303, dvrr_stack+156, dvrr_stack+293, dvrr_stack+12, dvrr_stack+6, NULL);

 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+0, dvrr_stack+166, dvrr_stack+303, dvrr_stack+65, dvrr_stack+156, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+293, dvrr_stack+181, dvrr_stack+0, dvrr_stack+75, dvrr_stack+166, NULL);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+321, dvrr_stack+202, dvrr_stack+293, dvrr_stack+90, dvrr_stack+181, NULL);

 /* compute (0 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+357,dvrr_stack+321,dvrr_stack+202,1);


 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+441, dvrr_stack+50, dvrr_stack+75, NULL, NULL, dvrr_stack+30);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+486, dvrr_stack+90, dvrr_stack+181, NULL, NULL, dvrr_stack+75);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+549, dvrr_stack+202, dvrr_stack+293, NULL, NULL, dvrr_stack+181);

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,15,dvrr_stack+293, dvrr_stack+111, NULL);
 tmp = dvrr_stack + 293;
 target_ptr = Libderiv->deriv_classes[0][4][11];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,21,dvrr_stack+0, dvrr_stack+230, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[0][5][11];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,28,dvrr_stack+156, dvrr_stack+357, NULL);
 tmp = dvrr_stack + 156;
 target_ptr = Libderiv->deriv_classes[0][6][11];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,15,dvrr_stack+21, dvrr_stack+111, NULL);
 tmp = dvrr_stack + 21;
 target_ptr = Libderiv->deriv_classes[0][4][10];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,21,dvrr_stack+65, dvrr_stack+230, NULL);
 tmp = dvrr_stack + 65;
 target_ptr = Libderiv->deriv_classes[0][5][10];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,28,dvrr_stack+633, dvrr_stack+357, NULL);
 tmp = dvrr_stack + 633;
 target_ptr = Libderiv->deriv_classes[0][6][10];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,15,dvrr_stack+184, dvrr_stack+111, NULL);
 tmp = dvrr_stack + 184;
 target_ptr = Libderiv->deriv_classes[0][4][9];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,21,dvrr_stack+111, dvrr_stack+230, NULL);
 tmp = dvrr_stack + 111;
 target_ptr = Libderiv->deriv_classes[0][5][9];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,28,dvrr_stack+230, dvrr_stack+357, NULL);
 tmp = dvrr_stack + 230;
 target_ptr = Libderiv->deriv_classes[0][6][9];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,1,1,dvrr_stack+357, dvrr_stack+90, dvrr_stack+40);
 tmp = dvrr_stack + 357;
 target_ptr = Libderiv->deriv_classes[0][4][8];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,1,1,dvrr_stack+372, dvrr_stack+202, dvrr_stack+50);
 tmp = dvrr_stack + 372;
 target_ptr = Libderiv->deriv_classes[0][5][8];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,1,1,dvrr_stack+393, dvrr_stack+321, dvrr_stack+90);
 tmp = dvrr_stack + 393;
 target_ptr = Libderiv->deriv_classes[0][6][8];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,1,1,dvrr_stack+421, dvrr_stack+90, dvrr_stack+40);
 tmp = dvrr_stack + 421;
 target_ptr = Libderiv->deriv_classes[0][4][7];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,1,1,dvrr_stack+258, dvrr_stack+202, dvrr_stack+50);
 tmp = dvrr_stack + 258;
 target_ptr = Libderiv->deriv_classes[0][5][7];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,1,1,dvrr_stack+661, dvrr_stack+321, dvrr_stack+90);
 tmp = dvrr_stack + 661;
 target_ptr = Libderiv->deriv_classes[0][6][7];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,1,1,dvrr_stack+132, dvrr_stack+90, dvrr_stack+40);
 tmp = dvrr_stack + 132;
 target_ptr = Libderiv->deriv_classes[0][4][6];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,1,1,dvrr_stack+689, dvrr_stack+202, dvrr_stack+50);
 tmp = dvrr_stack + 689;
 target_ptr = Libderiv->deriv_classes[0][5][6];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,1,1,dvrr_stack+36, dvrr_stack+321, dvrr_stack+90);
 tmp = dvrr_stack + 36;
 target_ptr = Libderiv->deriv_classes[0][6][6];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,15,dvrr_stack+308, dvrr_stack+441, NULL);
 tmp = dvrr_stack + 308;
 target_ptr = Libderiv->deriv_classes[0][4][2];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,21,dvrr_stack+323, dvrr_stack+486, NULL);
 tmp = dvrr_stack + 323;
 target_ptr = Libderiv->deriv_classes[0][5][2];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_0(Data,28,dvrr_stack+199, dvrr_stack+549, NULL);
 tmp = dvrr_stack + 199;
 target_ptr = Libderiv->deriv_classes[0][6][2];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,15,dvrr_stack+86, dvrr_stack+441, NULL);
 tmp = dvrr_stack + 86;
 target_ptr = Libderiv->deriv_classes[0][4][1];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,21,dvrr_stack+710, dvrr_stack+486, NULL);
 tmp = dvrr_stack + 710;
 target_ptr = Libderiv->deriv_classes[0][5][1];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_0(Data,28,dvrr_stack+731, dvrr_stack+549, NULL);
 tmp = dvrr_stack + 731;
 target_ptr = Libderiv->deriv_classes[0][6][1];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,15,dvrr_stack+759, dvrr_stack+441, NULL);
 tmp = dvrr_stack + 759;
 target_ptr = Libderiv->deriv_classes[0][4][0];
 for(i=0;i<15;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,21,dvrr_stack+436, dvrr_stack+486, NULL);
 tmp = dvrr_stack + 436;
 target_ptr = Libderiv->deriv_classes[0][5][0];
 for(i=0;i<21;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_0(Data,28,dvrr_stack+457, dvrr_stack+549, NULL);
 tmp = dvrr_stack + 457;
 target_ptr = Libderiv->deriv_classes[0][6][0];
 for(i=0;i<28;i++)
   target_ptr[i] += tmp[i];


}

