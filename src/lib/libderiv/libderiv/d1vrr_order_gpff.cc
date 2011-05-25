#include <stdio.h>
#include <libint/libint.h>
#include "libderiv.h"
#include <libint/vrr_header.h>
#include <libint/hrr_header.h>
#include "deriv_header.h"

  /* Computes quartets necessary to compute derivatives of (gp|ff) integrals */

void d1vrr_order_gpff(Libderiv_t *Libderiv, prim_data *Data)
{
 double *dvrr_stack = Libderiv->dvrr_stack;
 double *tmp, *target_ptr;
 int i, am[2];
 /* compute (1 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+6, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (0 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+9, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (1 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+12, dvrr_stack+3, dvrr_stack+9, NULL, NULL, Data->F+4);

 /* compute (1 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+21, dvrr_stack+6, dvrr_stack+3, NULL, NULL, Data->F+3);

 /* compute (2 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+30, dvrr_stack+21, dvrr_stack+12, dvrr_stack+6, dvrr_stack+3, dvrr_stack+0);

 /* compute (0 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+48, dvrr_stack+3, dvrr_stack+9, Data->F+3, Data->F+4, NULL);

 /* compute (0 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+54, dvrr_stack+6, dvrr_stack+3, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+60, dvrr_stack+54, dvrr_stack+48, NULL, NULL, dvrr_stack+3);

 /* compute (0 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+78, Data->F+1, Data->F+2, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+81, dvrr_stack+78, dvrr_stack+6, Data->F+1, Data->F+2, NULL);

 /* compute (1 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+87, dvrr_stack+81, dvrr_stack+54, NULL, NULL, dvrr_stack+6);

 /* compute (0 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+105, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+108, dvrr_stack+9, dvrr_stack+105, Data->F+4, Data->F+5, NULL);

 /* compute (1 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+114, dvrr_stack+48, dvrr_stack+108, NULL, NULL, dvrr_stack+9);

 /* compute (2 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+132, dvrr_stack+60, dvrr_stack+114, dvrr_stack+54, dvrr_stack+48, dvrr_stack+12);

 /* compute (2 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+168, dvrr_stack+87, dvrr_stack+60, dvrr_stack+81, dvrr_stack+54, dvrr_stack+21);

 /* compute (3 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+204, dvrr_stack+168, dvrr_stack+132, dvrr_stack+87, dvrr_stack+60, dvrr_stack+30);

 /* compute (0 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+264, dvrr_stack+54, dvrr_stack+48, dvrr_stack+6, dvrr_stack+3, NULL);

 /* compute (0 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+274, dvrr_stack+81, dvrr_stack+54, dvrr_stack+78, dvrr_stack+6, NULL);

 /* compute (0 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+284, dvrr_stack+48, dvrr_stack+108, dvrr_stack+3, dvrr_stack+9, NULL);

 /* compute (1 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+294, dvrr_stack+264, dvrr_stack+284, NULL, NULL, dvrr_stack+48);

 /* compute (1 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+324, dvrr_stack+274, dvrr_stack+264, NULL, NULL, dvrr_stack+54);

 /* compute (2 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+354, dvrr_stack+324, dvrr_stack+294, dvrr_stack+274, dvrr_stack+264, dvrr_stack+60);

 /* compute (0 0 | 1 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+414, Data->F+0, Data->F+1, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+417, dvrr_stack+414, dvrr_stack+78, Data->F+0, Data->F+1, NULL);

 /* compute (0 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+423, dvrr_stack+417, dvrr_stack+81, dvrr_stack+414, dvrr_stack+78, NULL);

 /* compute (1 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+433, dvrr_stack+423, dvrr_stack+274, NULL, NULL, dvrr_stack+81);

 /* compute (2 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+463, dvrr_stack+433, dvrr_stack+324, dvrr_stack+423, dvrr_stack+274, dvrr_stack+87);

 /* compute (0 0 | 1 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+414, Data->F+6, Data->F+7, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+523, dvrr_stack+105, dvrr_stack+414, Data->F+5, Data->F+6, NULL);

 /* compute (0 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+529, dvrr_stack+108, dvrr_stack+523, dvrr_stack+9, dvrr_stack+105, NULL);

 /* compute (1 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+539, dvrr_stack+284, dvrr_stack+529, NULL, NULL, dvrr_stack+108);

 /* compute (2 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+569, dvrr_stack+294, dvrr_stack+539, dvrr_stack+264, dvrr_stack+284, dvrr_stack+114);

 /* compute (3 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+629, dvrr_stack+354, dvrr_stack+569, dvrr_stack+324, dvrr_stack+294, dvrr_stack+132);

 /* compute (3 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+729, dvrr_stack+463, dvrr_stack+354, dvrr_stack+433, dvrr_stack+324, dvrr_stack+168);

 /* compute (4 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+829, dvrr_stack+729, dvrr_stack+629, dvrr_stack+463, dvrr_stack+354, dvrr_stack+204);
 tmp = dvrr_stack + 829;
 target_ptr = Libderiv->dvrr_classes[4][3];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (0 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+433, dvrr_stack+264, dvrr_stack+284, dvrr_stack+54, dvrr_stack+48, NULL);

 /* compute (0 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+448, dvrr_stack+274, dvrr_stack+264, dvrr_stack+81, dvrr_stack+54, NULL);

 /* compute (0 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+463, dvrr_stack+284, dvrr_stack+529, dvrr_stack+48, dvrr_stack+108, NULL);

 /* compute (1 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+478, dvrr_stack+433, dvrr_stack+463, NULL, NULL, dvrr_stack+284);

 /* compute (1 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+979, dvrr_stack+448, dvrr_stack+433, NULL, NULL, dvrr_stack+264);

 /* compute (2 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1024, dvrr_stack+979, dvrr_stack+478, dvrr_stack+448, dvrr_stack+433, dvrr_stack+294);

 /* compute (0 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1114, dvrr_stack+423, dvrr_stack+274, dvrr_stack+417, dvrr_stack+81, NULL);

 /* compute (1 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1129, dvrr_stack+1114, dvrr_stack+448, NULL, NULL, dvrr_stack+274);

 /* compute (2 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1174, dvrr_stack+1129, dvrr_stack+979, dvrr_stack+1114, dvrr_stack+448, dvrr_stack+324);

 /* compute (0 0 | 1 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+324, Data->F+7, Data->F+8, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+54, dvrr_stack+414, dvrr_stack+324, Data->F+6, Data->F+7, NULL);

 /* compute (0 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+327, dvrr_stack+523, dvrr_stack+54, dvrr_stack+105, dvrr_stack+414, NULL);

 /* compute (0 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+337, dvrr_stack+529, dvrr_stack+327, dvrr_stack+108, dvrr_stack+523, NULL);

 /* compute (1 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1264, dvrr_stack+463, dvrr_stack+337, NULL, NULL, dvrr_stack+529);

 /* compute (2 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1309, dvrr_stack+478, dvrr_stack+1264, dvrr_stack+433, dvrr_stack+463, dvrr_stack+539);

 /* compute (3 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1399, dvrr_stack+1024, dvrr_stack+1309, dvrr_stack+979, dvrr_stack+478, dvrr_stack+569);

 /* compute (3 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1549, dvrr_stack+1174, dvrr_stack+1024, dvrr_stack+1129, dvrr_stack+979, dvrr_stack+354);

 /* compute (4 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1699, dvrr_stack+1549, dvrr_stack+1399, dvrr_stack+1174, dvrr_stack+1024, dvrr_stack+629);
 tmp = dvrr_stack + 1699;
 target_ptr = Libderiv->dvrr_classes[4][4];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+1924,dvrr_stack+1699,dvrr_stack+829,15);


 /* compute (0 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1129, dvrr_stack+433, dvrr_stack+463, dvrr_stack+264, dvrr_stack+284, NULL);

 /* compute (0 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1150, dvrr_stack+448, dvrr_stack+433, dvrr_stack+274, dvrr_stack+264, NULL);

 /* compute (0 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1171, dvrr_stack+463, dvrr_stack+337, dvrr_stack+284, dvrr_stack+529, NULL);

 /* compute (1 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1192, dvrr_stack+1129, dvrr_stack+1171, NULL, NULL, dvrr_stack+463);

 /* compute (1 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2374, dvrr_stack+1150, dvrr_stack+1129, NULL, NULL, dvrr_stack+433);

 /* compute (2 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2437, dvrr_stack+2374, dvrr_stack+1192, dvrr_stack+1150, dvrr_stack+1129, dvrr_stack+478);

 /* compute (0 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2563, dvrr_stack+1114, dvrr_stack+448, dvrr_stack+423, dvrr_stack+274, NULL);

 /* compute (1 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2584, dvrr_stack+2563, dvrr_stack+1150, NULL, NULL, dvrr_stack+448);

 /* compute (2 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2647, dvrr_stack+2584, dvrr_stack+2374, dvrr_stack+2563, dvrr_stack+1150, dvrr_stack+979);

 /* compute (0 0 | 1 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+979, Data->F+8, Data->F+9, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+982, dvrr_stack+324, dvrr_stack+979, Data->F+7, Data->F+8, NULL);

 /* compute (0 0 | 3 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+423, dvrr_stack+54, dvrr_stack+982, dvrr_stack+414, dvrr_stack+324, NULL);

 /* compute (0 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+988, dvrr_stack+327, dvrr_stack+423, dvrr_stack+523, dvrr_stack+54, NULL);

 /* compute (0 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1003, dvrr_stack+337, dvrr_stack+988, dvrr_stack+529, dvrr_stack+327, NULL);

 /* compute (1 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2773, dvrr_stack+1171, dvrr_stack+1003, NULL, NULL, dvrr_stack+337);

 /* compute (2 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2836, dvrr_stack+1192, dvrr_stack+2773, dvrr_stack+1129, dvrr_stack+1171, dvrr_stack+1264);

 /* compute (3 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2962, dvrr_stack+2437, dvrr_stack+2836, dvrr_stack+2374, dvrr_stack+1192, dvrr_stack+1309);

 /* compute (3 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3172, dvrr_stack+2647, dvrr_stack+2437, dvrr_stack+2584, dvrr_stack+2374, dvrr_stack+1024);

 /* compute (4 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+3382, dvrr_stack+3172, dvrr_stack+2962, dvrr_stack+2647, dvrr_stack+2437, dvrr_stack+1399);
 tmp = dvrr_stack + 3382;
 target_ptr = Libderiv->dvrr_classes[4][5];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+3697,dvrr_stack+3382,dvrr_stack+1699,15);


 /* compute (0 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2584, dvrr_stack+1129, dvrr_stack+1171, dvrr_stack+433, dvrr_stack+463, NULL);

 /* compute (0 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2612, dvrr_stack+1150, dvrr_stack+1129, dvrr_stack+448, dvrr_stack+433, NULL);

 /* compute (0 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2640, dvrr_stack+1171, dvrr_stack+1003, dvrr_stack+463, dvrr_stack+337, NULL);

 /* compute (1 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2668, dvrr_stack+2584, dvrr_stack+2640, NULL, NULL, dvrr_stack+1171);

 /* compute (1 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4372, dvrr_stack+2612, dvrr_stack+2584, NULL, NULL, dvrr_stack+1129);

 /* compute (2 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4456, dvrr_stack+4372, dvrr_stack+2668, dvrr_stack+2612, dvrr_stack+2584, dvrr_stack+1192);

 /* compute (0 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4624, dvrr_stack+2563, dvrr_stack+1150, dvrr_stack+1114, dvrr_stack+448, NULL);

 /* compute (1 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4652, dvrr_stack+4624, dvrr_stack+2612, NULL, NULL, dvrr_stack+1150);

 /* compute (2 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4736, dvrr_stack+4652, dvrr_stack+4372, dvrr_stack+4624, dvrr_stack+2612, dvrr_stack+2374);

 /* compute (0 0 | 1 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+2374, Data->F+9, Data->F+10, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+2377, dvrr_stack+979, dvrr_stack+2374, Data->F+8, Data->F+9, NULL);

 /* compute (0 0 | 3 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+2383, dvrr_stack+982, dvrr_stack+2377, dvrr_stack+324, dvrr_stack+979, NULL);

 /* compute (0 0 | 4 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1114, dvrr_stack+423, dvrr_stack+2383, dvrr_stack+54, dvrr_stack+982, NULL);

 /* compute (0 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+2752, dvrr_stack+988, dvrr_stack+1114, dvrr_stack+327, dvrr_stack+423, NULL);

 /* compute (0 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2393, dvrr_stack+1003, dvrr_stack+2752, dvrr_stack+337, dvrr_stack+988, NULL);

 /* compute (1 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4904, dvrr_stack+2640, dvrr_stack+2393, NULL, NULL, dvrr_stack+1003);

 /* compute (2 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4988, dvrr_stack+2668, dvrr_stack+4904, dvrr_stack+2584, dvrr_stack+2640, dvrr_stack+2773);

 /* compute (3 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+5156, dvrr_stack+4456, dvrr_stack+4988, dvrr_stack+4372, dvrr_stack+2668, dvrr_stack+2836);

 /* compute (3 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+5436, dvrr_stack+4736, dvrr_stack+4456, dvrr_stack+4652, dvrr_stack+4372, dvrr_stack+2437);

 /* compute (4 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+5716, dvrr_stack+5436, dvrr_stack+5156, dvrr_stack+4736, dvrr_stack+4456, dvrr_stack+2962);
 tmp = dvrr_stack + 5716;
 target_ptr = Libderiv->dvrr_classes[4][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+6136,dvrr_stack+5716,dvrr_stack+3382,15);


 /* compute (0 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4652, dvrr_stack+2584, dvrr_stack+2640, dvrr_stack+1129, dvrr_stack+1171, NULL);

 /* compute (0 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4688, dvrr_stack+2612, dvrr_stack+2584, dvrr_stack+1150, dvrr_stack+1129, NULL);

 /* compute (0 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4724, dvrr_stack+2640, dvrr_stack+2393, dvrr_stack+1171, dvrr_stack+1003, NULL);

 /* compute (1 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4760, dvrr_stack+4652, dvrr_stack+4724, NULL, NULL, dvrr_stack+2640);

 /* compute (1 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7081, dvrr_stack+4688, dvrr_stack+4652, NULL, NULL, dvrr_stack+2584);

 /* compute (2 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7189, dvrr_stack+7081, dvrr_stack+4760, dvrr_stack+4688, dvrr_stack+4652, dvrr_stack+2668);

 /* compute (0 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4868, dvrr_stack+4624, dvrr_stack+2612, dvrr_stack+2563, dvrr_stack+1150, NULL);

 /* compute (1 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7405, dvrr_stack+4868, dvrr_stack+4688, NULL, NULL, dvrr_stack+2612);

 /* compute (2 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7513, dvrr_stack+7405, dvrr_stack+7081, dvrr_stack+4868, dvrr_stack+4688, dvrr_stack+4372);

 /* compute (0 0 | 1 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+324, Data->F+10, Data->F+11, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+4372, dvrr_stack+2374, dvrr_stack+324, Data->F+9, Data->F+10, NULL);

 /* compute (0 0 | 3 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+4378, dvrr_stack+2377, dvrr_stack+4372, dvrr_stack+979, dvrr_stack+2374, NULL);

 /* compute (0 0 | 4 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+4388, dvrr_stack+2383, dvrr_stack+4378, dvrr_stack+982, dvrr_stack+2377, NULL);

 /* compute (0 0 | 5 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4403, dvrr_stack+1114, dvrr_stack+4388, dvrr_stack+423, dvrr_stack+2383, NULL);

 /* compute (0 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4624, dvrr_stack+2752, dvrr_stack+4403, dvrr_stack+988, dvrr_stack+1114, NULL);

 /* compute (0 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4688, dvrr_stack+2393, dvrr_stack+4624, dvrr_stack+1003, dvrr_stack+2752, NULL);

 /* compute (1 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7729, dvrr_stack+4724, dvrr_stack+4688, NULL, NULL, dvrr_stack+2393);

 /* compute (2 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+7837, dvrr_stack+4760, dvrr_stack+7729, dvrr_stack+4652, dvrr_stack+4724, dvrr_stack+4904);

 /* compute (3 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+8053, dvrr_stack+7189, dvrr_stack+7837, dvrr_stack+7081, dvrr_stack+4760, dvrr_stack+4988);

 /* compute (3 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+8413, dvrr_stack+7513, dvrr_stack+7189, dvrr_stack+7405, dvrr_stack+7081, dvrr_stack+4456);

 /* compute (4 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+8773, dvrr_stack+8413, dvrr_stack+8053, dvrr_stack+7513, dvrr_stack+7189, dvrr_stack+5156);

 /* compute (4 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+9313,dvrr_stack+8773,dvrr_stack+5716,15);


 /* compute (1 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+7081, Data->F+4, Data->F+5, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+7084, dvrr_stack+0, dvrr_stack+7081, Data->F+3, Data->F+4, NULL);

 /* compute (1 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+979, dvrr_stack+9, dvrr_stack+105, NULL, NULL, Data->F+5);

 /* compute (2 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+7090, dvrr_stack+12, dvrr_stack+979, dvrr_stack+3, dvrr_stack+9, dvrr_stack+7081);

 /* compute (3 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+433, dvrr_stack+30, dvrr_stack+7090, dvrr_stack+21, dvrr_stack+12, dvrr_stack+7084);

 /* compute (1 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+7108, dvrr_stack+108, dvrr_stack+523, NULL, NULL, dvrr_stack+105);

 /* compute (2 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+4652, dvrr_stack+114, dvrr_stack+7108, dvrr_stack+48, dvrr_stack+108, dvrr_stack+979);

 /* compute (3 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+7126, dvrr_stack+132, dvrr_stack+4652, dvrr_stack+60, dvrr_stack+114, dvrr_stack+7090);

 /* compute (4 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+7405, dvrr_stack+204, dvrr_stack+7126, dvrr_stack+168, dvrr_stack+132, dvrr_stack+433);

 /* compute (1 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+7495, dvrr_stack+529, dvrr_stack+327, NULL, NULL, dvrr_stack+523);

 /* compute (2 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+7525, dvrr_stack+539, dvrr_stack+7495, dvrr_stack+284, dvrr_stack+529, dvrr_stack+7108);

 /* compute (3 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+7585, dvrr_stack+569, dvrr_stack+7525, dvrr_stack+294, dvrr_stack+539, dvrr_stack+4652);

 /* compute (4 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10573, dvrr_stack+629, dvrr_stack+7585, dvrr_stack+354, dvrr_stack+569, dvrr_stack+7126);

 /* compute (5 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+10723, dvrr_stack+829, dvrr_stack+10573, dvrr_stack+729, dvrr_stack+629, dvrr_stack+7405);
 tmp = dvrr_stack + 10723;
 target_ptr = Libderiv->dvrr_classes[5][3];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (1 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+2563, dvrr_stack+337, dvrr_stack+988, NULL, NULL, dvrr_stack+327);

 /* compute (2 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10933, dvrr_stack+1264, dvrr_stack+2563, dvrr_stack+463, dvrr_stack+337, dvrr_stack+7495);

 /* compute (3 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11023, dvrr_stack+1309, dvrr_stack+10933, dvrr_stack+478, dvrr_stack+1264, dvrr_stack+7525);

 /* compute (4 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11173, dvrr_stack+1399, dvrr_stack+11023, dvrr_stack+1024, dvrr_stack+1309, dvrr_stack+7585);

 /* compute (5 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+11398, dvrr_stack+1699, dvrr_stack+11173, dvrr_stack+1549, dvrr_stack+1399, dvrr_stack+10573);
 tmp = dvrr_stack + 11398;
 target_ptr = Libderiv->dvrr_classes[5][4];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_fp(Libderiv->CD,dvrr_stack+11713,dvrr_stack+11398,dvrr_stack+10723,21);


 /* compute (1 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+1024, dvrr_stack+1003, dvrr_stack+2752, NULL, NULL, dvrr_stack+988);

 /* compute (2 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+12343, dvrr_stack+2773, dvrr_stack+1024, dvrr_stack+1171, dvrr_stack+1003, dvrr_stack+2563);

 /* compute (3 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+12469, dvrr_stack+2836, dvrr_stack+12343, dvrr_stack+1192, dvrr_stack+2773, dvrr_stack+10933);

 /* compute (4 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+12679, dvrr_stack+2962, dvrr_stack+12469, dvrr_stack+2437, dvrr_stack+2836, dvrr_stack+11023);

 /* compute (5 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+12994, dvrr_stack+3382, dvrr_stack+12679, dvrr_stack+3172, dvrr_stack+2962, dvrr_stack+11173);
 tmp = dvrr_stack + 12994;
 target_ptr = Libderiv->dvrr_classes[5][5];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_gp(Libderiv->CD,dvrr_stack+13435,dvrr_stack+12994,dvrr_stack+11398,21);


 /* compute (1 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+1129, dvrr_stack+2393, dvrr_stack+4624, NULL, NULL, dvrr_stack+2752);

 /* compute (2 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14380, dvrr_stack+4904, dvrr_stack+1129, dvrr_stack+2640, dvrr_stack+2393, dvrr_stack+1024);

 /* compute (3 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14548, dvrr_stack+4988, dvrr_stack+14380, dvrr_stack+2668, dvrr_stack+4904, dvrr_stack+12343);

 /* compute (4 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+14828, dvrr_stack+5156, dvrr_stack+14548, dvrr_stack+4456, dvrr_stack+4988, dvrr_stack+12469);

 /* compute (5 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+15248, dvrr_stack+5716, dvrr_stack+14828, dvrr_stack+5436, dvrr_stack+5156, dvrr_stack+12679);

 /* compute (5 0 | 5 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_hp(Libderiv->CD,dvrr_stack+15836,dvrr_stack+15248,dvrr_stack+12994,21);


 /* compute (0 0 | 1 0) m=11 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00p0(Data,dvrr_stack+3, Data->F+11, Data->F+12, NULL, NULL, NULL);

 /* compute (0 0 | 2 0) m=10 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00d0(Data,dvrr_stack+48, dvrr_stack+324, dvrr_stack+3, Data->F+10, Data->F+11, NULL);

 /* compute (0 0 | 3 0) m=9 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_00f0(Data,dvrr_stack+1087, dvrr_stack+4372, dvrr_stack+48, dvrr_stack+2374, dvrr_stack+324, NULL);

 /* compute (0 0 | 4 0) m=8 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+1097, dvrr_stack+4378, dvrr_stack+1087, dvrr_stack+2377, dvrr_stack+4372, NULL);

 /* compute (0 0 | 5 0) m=7 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+463, dvrr_stack+4388, dvrr_stack+1097, dvrr_stack+2383, dvrr_stack+4378, NULL);

 /* compute (0 0 | 6 0) m=6 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+484, dvrr_stack+4403, dvrr_stack+463, dvrr_stack+1114, dvrr_stack+4388, NULL);

 /* compute (0 0 | 7 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4868, dvrr_stack+4624, dvrr_stack+484, dvrr_stack+2752, dvrr_stack+4403, NULL);

 /* compute (1 0 | 7 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+4424, dvrr_stack+4688, dvrr_stack+4868, NULL, NULL, dvrr_stack+4624);

 /* compute (2 0 | 7 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+17159, dvrr_stack+7729, dvrr_stack+4424, dvrr_stack+4724, dvrr_stack+4688, dvrr_stack+1129);

 /* compute (3 0 | 7 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+17375, dvrr_stack+7837, dvrr_stack+17159, dvrr_stack+4760, dvrr_stack+7729, dvrr_stack+14380);

 /* compute (4 0 | 7 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+17735, dvrr_stack+8053, dvrr_stack+17375, dvrr_stack+7189, dvrr_stack+7837, dvrr_stack+14548);

 /* compute (5 0 | 7 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 7;
 vrr_build_xxxx(am,Data,dvrr_stack+18275, dvrr_stack+8773, dvrr_stack+17735, dvrr_stack+8413, dvrr_stack+8053, dvrr_stack+14828);

 /* compute (5 0 | 6 1) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 hrr3_build_ip(Libderiv->CD,dvrr_stack+19031,dvrr_stack+18275,dvrr_stack+15248,21);


 /* compute (1 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+48, dvrr_stack+3, dvrr_stack+0, Data->F+2, Data->F+3, NULL);

 /* compute (1 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+17159, dvrr_stack+78, dvrr_stack+6, NULL, NULL, Data->F+2);

 /* compute (2 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+60, dvrr_stack+17159, dvrr_stack+21, dvrr_stack+78, dvrr_stack+6, dvrr_stack+3);

 /* compute (3 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+17168, dvrr_stack+60, dvrr_stack+30, dvrr_stack+17159, dvrr_stack+21, dvrr_stack+48);

 /* compute (1 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+17198, dvrr_stack+417, dvrr_stack+81, NULL, NULL, dvrr_stack+78);

 /* compute (2 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+17216, dvrr_stack+17198, dvrr_stack+87, dvrr_stack+417, dvrr_stack+81, dvrr_stack+17159);

 /* compute (3 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+17252, dvrr_stack+17216, dvrr_stack+168, dvrr_stack+17198, dvrr_stack+87, dvrr_stack+60);

 /* compute (4 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+17312, dvrr_stack+17252, dvrr_stack+204, dvrr_stack+17216, dvrr_stack+168, dvrr_stack+17168);

 /* compute (3 0 | 0 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+168, dvrr_stack+48, dvrr_stack+7084, dvrr_stack+3, dvrr_stack+0, NULL);

 /* compute (4 0 | 1 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+17198, dvrr_stack+17168, dvrr_stack+433, dvrr_stack+60, dvrr_stack+30, dvrr_stack+168);

 /* compute (5 0 | 2 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+17402, dvrr_stack+17312, dvrr_stack+7405, dvrr_stack+17252, dvrr_stack+204, dvrr_stack+17198);

 /* compute (1 0 | 0 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p000(Data,dvrr_stack+168, Data->F+5, Data->F+6, NULL, NULL, NULL);

 /* compute (2 0 | 0 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d000(Data,dvrr_stack+3, dvrr_stack+7081, dvrr_stack+168, Data->F+4, Data->F+5, NULL);

 /* compute (3 0 | 0 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f000(Data,dvrr_stack+171, dvrr_stack+7084, dvrr_stack+3, dvrr_stack+0, dvrr_stack+7081, NULL);

 /* compute (1 0 | 1 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0p0(Data,dvrr_stack+7081, dvrr_stack+105, dvrr_stack+414, NULL, NULL, Data->F+6);

 /* compute (2 0 | 1 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0p0(Data,dvrr_stack+181, dvrr_stack+979, dvrr_stack+7081, dvrr_stack+9, dvrr_stack+105, dvrr_stack+168);

 /* compute (3 0 | 1 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0p0(Data,dvrr_stack+60, dvrr_stack+7090, dvrr_stack+181, dvrr_stack+12, dvrr_stack+979, dvrr_stack+3);

 /* compute (4 0 | 1 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 1;
 vrr_build_xxxx(am,Data,dvrr_stack+17159, dvrr_stack+433, dvrr_stack+60, dvrr_stack+30, dvrr_stack+7090, dvrr_stack+171);

 /* compute (1 0 | 2 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0d0(Data,dvrr_stack+7090, dvrr_stack+523, dvrr_stack+54, NULL, NULL, dvrr_stack+414);

 /* compute (2 0 | 2 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0d0(Data,dvrr_stack+433, dvrr_stack+7108, dvrr_stack+7090, dvrr_stack+108, dvrr_stack+523, dvrr_stack+7081);

 /* compute (3 0 | 2 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0d0(Data,dvrr_stack+17204, dvrr_stack+4652, dvrr_stack+433, dvrr_stack+114, dvrr_stack+7108, dvrr_stack+181);

 /* compute (4 0 | 2 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+17528, dvrr_stack+7126, dvrr_stack+17204, dvrr_stack+132, dvrr_stack+4652, dvrr_stack+60);

 /* compute (5 0 | 2 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 2;
 vrr_build_xxxx(am,Data,dvrr_stack+60, dvrr_stack+7405, dvrr_stack+17528, dvrr_stack+204, dvrr_stack+7126, dvrr_stack+17159);

 /* compute (1 0 | 3 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_p0f0(Data,dvrr_stack+17159, dvrr_stack+327, dvrr_stack+423, NULL, NULL, dvrr_stack+54);

 /* compute (2 0 | 3 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_d0f0(Data,dvrr_stack+0, dvrr_stack+7495, dvrr_stack+17159, dvrr_stack+529, dvrr_stack+327, dvrr_stack+7090);

 /* compute (3 0 | 3 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 _BUILD_f0f0(Data,dvrr_stack+186, dvrr_stack+7525, dvrr_stack+0, dvrr_stack+539, dvrr_stack+7495, dvrr_stack+433);

 /* compute (4 0 | 3 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+4652, dvrr_stack+7585, dvrr_stack+186, dvrr_stack+569, dvrr_stack+7525, dvrr_stack+17204);

 /* compute (5 0 | 3 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+7081, dvrr_stack+10573, dvrr_stack+4652, dvrr_stack+629, dvrr_stack+7585, dvrr_stack+17528);

 /* compute (6 0 | 3 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 3;
 vrr_build_xxxx(am,Data,dvrr_stack+17528, dvrr_stack+10723, dvrr_stack+7081, dvrr_stack+829, dvrr_stack+10573, dvrr_stack+60);

 /* compute (1 0 | 4 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+60, dvrr_stack+988, dvrr_stack+1114, NULL, NULL, dvrr_stack+423);

 /* compute (2 0 | 4 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+10573, dvrr_stack+2563, dvrr_stack+60, dvrr_stack+337, dvrr_stack+988, dvrr_stack+17159);

 /* compute (3 0 | 4 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+17159, dvrr_stack+10933, dvrr_stack+10573, dvrr_stack+1264, dvrr_stack+2563, dvrr_stack+0);

 /* compute (4 0 | 4 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+17808, dvrr_stack+11023, dvrr_stack+17159, dvrr_stack+1309, dvrr_stack+10933, dvrr_stack+186);

 /* compute (5 0 | 4 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+105, dvrr_stack+11173, dvrr_stack+17808, dvrr_stack+1399, dvrr_stack+11023, dvrr_stack+4652);

 /* compute (6 0 | 4 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 4;
 vrr_build_xxxx(am,Data,dvrr_stack+7291, dvrr_stack+11398, dvrr_stack+105, dvrr_stack+1699, dvrr_stack+11173, dvrr_stack+7081);

 /* compute (1 0 | 5 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+7081, dvrr_stack+2752, dvrr_stack+4403, NULL, NULL, dvrr_stack+1114);

 /* compute (2 0 | 5 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+7144, dvrr_stack+1024, dvrr_stack+7081, dvrr_stack+1003, dvrr_stack+2752, dvrr_stack+60);

 /* compute (3 0 | 5 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+4652, dvrr_stack+12343, dvrr_stack+7144, dvrr_stack+2773, dvrr_stack+1024, dvrr_stack+10573);

 /* compute (4 0 | 5 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+10933, dvrr_stack+12469, dvrr_stack+4652, dvrr_stack+2836, dvrr_stack+12343, dvrr_stack+17159);

 /* compute (5 0 | 5 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+7711, dvrr_stack+12679, dvrr_stack+10933, dvrr_stack+2962, dvrr_stack+12469, dvrr_stack+17808);

 /* compute (6 0 | 5 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 5;
 vrr_build_xxxx(am,Data,dvrr_stack+8152, dvrr_stack+12994, dvrr_stack+7711, dvrr_stack+3382, dvrr_stack+12679, dvrr_stack+105);

 /* compute (1 0 | 6 0) m=5 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 1;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+17808, dvrr_stack+4624, dvrr_stack+484, NULL, NULL, dvrr_stack+4403);

 /* compute (2 0 | 6 0) m=4 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 2;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+17892, dvrr_stack+1129, dvrr_stack+17808, dvrr_stack+2393, dvrr_stack+4624, dvrr_stack+7081);

 /* compute (3 0 | 6 0) m=3 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 3;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4372, dvrr_stack+14380, dvrr_stack+17892, dvrr_stack+4904, dvrr_stack+1129, dvrr_stack+7144);

 /* compute (4 0 | 6 0) m=2 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 4;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+17808, dvrr_stack+14548, dvrr_stack+4372, dvrr_stack+4988, dvrr_stack+14380, dvrr_stack+4652);

 /* compute (5 0 | 6 0) m=1 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 5;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+4372, dvrr_stack+14828, dvrr_stack+17808, dvrr_stack+5156, dvrr_stack+14548, dvrr_stack+10933);

 /* compute (6 0 | 6 0) m=0 deriv level 0 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 0 */
 am[0] = 6;  am[1] = 6;
 vrr_build_xxxx(am,Data,dvrr_stack+2374, dvrr_stack+15248, dvrr_stack+4372, dvrr_stack+5716, dvrr_stack+14828, dvrr_stack+7711);

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,150,dvrr_stack+10573, dvrr_stack+1924, NULL);
 tmp = dvrr_stack + 10573;
 target_ptr = Libderiv->deriv_classes[4][3][11];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,225,dvrr_stack+7711, dvrr_stack+3697, NULL);
 tmp = dvrr_stack + 7711;
 target_ptr = Libderiv->deriv_classes[4][4][11];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+4372, dvrr_stack+6136, NULL);
 tmp = dvrr_stack + 4372;
 target_ptr = Libderiv->deriv_classes[4][5][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,420,dvrr_stack+4687, dvrr_stack+9313, NULL);
 tmp = dvrr_stack + 4687;
 target_ptr = Libderiv->deriv_classes[4][6][11];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,210,dvrr_stack+7081, dvrr_stack+11713, NULL);
 tmp = dvrr_stack + 7081;
 target_ptr = Libderiv->deriv_classes[5][3][11];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,315,dvrr_stack+5107, dvrr_stack+13435, NULL);
 tmp = dvrr_stack + 5107;
 target_ptr = Libderiv->deriv_classes[5][4][11];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,441,dvrr_stack+10933, dvrr_stack+15836, NULL);
 tmp = dvrr_stack + 10933;
 target_ptr = Libderiv->deriv_classes[5][5][11];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 0 1 */
 deriv_build_DZ_0(Data,588,dvrr_stack+14380, dvrr_stack+19031, NULL);
 tmp = dvrr_stack + 14380;
 target_ptr = Libderiv->deriv_classes[5][6][11];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,150,dvrr_stack+7936, dvrr_stack+1924, NULL);
 tmp = dvrr_stack + 7936;
 target_ptr = Libderiv->deriv_classes[4][3][10];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,225,dvrr_stack+17808, dvrr_stack+3697, NULL);
 tmp = dvrr_stack + 17808;
 target_ptr = Libderiv->deriv_classes[4][4][10];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+12343, dvrr_stack+6136, NULL);
 tmp = dvrr_stack + 12343;
 target_ptr = Libderiv->deriv_classes[4][5][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,420,dvrr_stack+0, dvrr_stack+9313, NULL);
 tmp = dvrr_stack + 0;
 target_ptr = Libderiv->deriv_classes[4][6][10];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,210,dvrr_stack+18033, dvrr_stack+11713, NULL);
 tmp = dvrr_stack + 18033;
 target_ptr = Libderiv->deriv_classes[5][3][10];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,315,dvrr_stack+12658, dvrr_stack+13435, NULL);
 tmp = dvrr_stack + 12658;
 target_ptr = Libderiv->deriv_classes[5][4][10];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,441,dvrr_stack+979, dvrr_stack+15836, NULL);
 tmp = dvrr_stack + 979;
 target_ptr = Libderiv->deriv_classes[5][5][10];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  0 1 0 */
 deriv_build_DY_0(Data,588,dvrr_stack+20795, dvrr_stack+19031, NULL);
 tmp = dvrr_stack + 20795;
 target_ptr = Libderiv->deriv_classes[5][6][10];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,150,dvrr_stack+14968, dvrr_stack+1924, NULL);
 tmp = dvrr_stack + 14968;
 target_ptr = Libderiv->deriv_classes[4][3][9];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,225,dvrr_stack+1924, dvrr_stack+3697, NULL);
 tmp = dvrr_stack + 1924;
 target_ptr = Libderiv->deriv_classes[4][4][9];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+3697, dvrr_stack+6136, NULL);
 tmp = dvrr_stack + 3697;
 target_ptr = Libderiv->deriv_classes[4][5][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,420,dvrr_stack+6136, dvrr_stack+9313, NULL);
 tmp = dvrr_stack + 6136;
 target_ptr = Libderiv->deriv_classes[4][6][9];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,210,dvrr_stack+9313, dvrr_stack+11713, NULL);
 tmp = dvrr_stack + 9313;
 target_ptr = Libderiv->deriv_classes[5][3][9];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,315,dvrr_stack+11713, dvrr_stack+13435, NULL);
 tmp = dvrr_stack + 11713;
 target_ptr = Libderiv->deriv_classes[5][4][9];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,441,dvrr_stack+13435, dvrr_stack+15836, NULL);
 tmp = dvrr_stack + 13435;
 target_ptr = Libderiv->deriv_classes[5][5][9];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 0  1 0 0 */
 deriv_build_DX_0(Data,588,dvrr_stack+15836, dvrr_stack+19031, NULL);
 tmp = dvrr_stack + 15836;
 target_ptr = Libderiv->deriv_classes[5][6][9];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,15,1,dvrr_stack+19031, dvrr_stack+1699, dvrr_stack+17312);
 tmp = dvrr_stack + 19031;
 target_ptr = Libderiv->deriv_classes[4][3][8];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,15,1,dvrr_stack+2149, dvrr_stack+3382, dvrr_stack+829);
 tmp = dvrr_stack + 2149;
 target_ptr = Libderiv->deriv_classes[4][4][8];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,15,1,dvrr_stack+12028, dvrr_stack+5716, dvrr_stack+1699);
 tmp = dvrr_stack + 12028;
 target_ptr = Libderiv->deriv_classes[4][5][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,15,1,dvrr_stack+19181, dvrr_stack+8773, dvrr_stack+3382);
 tmp = dvrr_stack + 19181;
 target_ptr = Libderiv->deriv_classes[4][6][8];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_f(Data,21,1,dvrr_stack+19601, dvrr_stack+11398, dvrr_stack+17402);
 tmp = dvrr_stack + 19601;
 target_ptr = Libderiv->deriv_classes[5][3][8];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_g(Data,21,1,dvrr_stack+19811, dvrr_stack+12994, dvrr_stack+10723);
 tmp = dvrr_stack + 19811;
 target_ptr = Libderiv->deriv_classes[5][4][8];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_h(Data,21,1,dvrr_stack+20126, dvrr_stack+15248, dvrr_stack+11398);
 tmp = dvrr_stack + 20126;
 target_ptr = Libderiv->deriv_classes[5][5][8];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 0 1  0 0 0 */
 deriv_build_CZ_i(Data,21,1,dvrr_stack+16424, dvrr_stack+18275, dvrr_stack+12994);
 tmp = dvrr_stack + 16424;
 target_ptr = Libderiv->deriv_classes[5][6][8];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,15,1,dvrr_stack+20567, dvrr_stack+1699, dvrr_stack+17312);
 tmp = dvrr_stack + 20567;
 target_ptr = Libderiv->deriv_classes[4][3][7];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,15,1,dvrr_stack+17012, dvrr_stack+3382, dvrr_stack+829);
 tmp = dvrr_stack + 17012;
 target_ptr = Libderiv->deriv_classes[4][4][7];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,15,1,dvrr_stack+13876, dvrr_stack+5716, dvrr_stack+1699);
 tmp = dvrr_stack + 13876;
 target_ptr = Libderiv->deriv_classes[4][5][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,15,1,dvrr_stack+9523, dvrr_stack+8773, dvrr_stack+3382);
 tmp = dvrr_stack + 9523;
 target_ptr = Libderiv->deriv_classes[4][6][7];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_f(Data,21,1,dvrr_stack+9943, dvrr_stack+11398, dvrr_stack+17402);
 tmp = dvrr_stack + 9943;
 target_ptr = Libderiv->deriv_classes[5][3][7];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_g(Data,21,1,dvrr_stack+10153, dvrr_stack+12994, dvrr_stack+10723);
 tmp = dvrr_stack + 10153;
 target_ptr = Libderiv->deriv_classes[5][4][7];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_h(Data,21,1,dvrr_stack+6556, dvrr_stack+15248, dvrr_stack+11398);
 tmp = dvrr_stack + 6556;
 target_ptr = Libderiv->deriv_classes[5][5][7];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  0 1 0  0 0 0 */
 deriv_build_CY_i(Data,21,1,dvrr_stack+21383, dvrr_stack+18275, dvrr_stack+12994);
 tmp = dvrr_stack + 21383;
 target_ptr = Libderiv->deriv_classes[5][6][7];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,15,1,dvrr_stack+14191, dvrr_stack+1699, dvrr_stack+17312);
 tmp = dvrr_stack + 14191;
 target_ptr = Libderiv->deriv_classes[4][3][6];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,15,1,dvrr_stack+4012, dvrr_stack+3382, dvrr_stack+829);
 tmp = dvrr_stack + 4012;
 target_ptr = Libderiv->deriv_classes[4][4][6];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,15,1,dvrr_stack+21971, dvrr_stack+5716, dvrr_stack+1699);
 tmp = dvrr_stack + 21971;
 target_ptr = Libderiv->deriv_classes[4][5][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,15,1,dvrr_stack+22286, dvrr_stack+8773, dvrr_stack+3382);
 tmp = dvrr_stack + 22286;
 target_ptr = Libderiv->deriv_classes[4][6][6];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_f(Data,21,1,dvrr_stack+420, dvrr_stack+11398, dvrr_stack+17402);
 tmp = dvrr_stack + 420;
 target_ptr = Libderiv->deriv_classes[5][3][6];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_g(Data,21,1,dvrr_stack+8740, dvrr_stack+12994, dvrr_stack+10723);
 tmp = dvrr_stack + 8740;
 target_ptr = Libderiv->deriv_classes[5][4][6];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_h(Data,21,1,dvrr_stack+22706, dvrr_stack+15248, dvrr_stack+11398);
 tmp = dvrr_stack + 22706;
 target_ptr = Libderiv->deriv_classes[5][5][6];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 0  0 0 0  1 0 0  0 0 0 */
 deriv_build_CX_i(Data,21,1,dvrr_stack+23147, dvrr_stack+18275, dvrr_stack+12994);
 tmp = dvrr_stack + 23147;
 target_ptr = Libderiv->deriv_classes[5][6][6];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,10,dvrr_stack+17237, dvrr_stack+10723, dvrr_stack+729);
 tmp = dvrr_stack + 17237;
 target_ptr = Libderiv->deriv_classes[4][3][2];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,15,dvrr_stack+18243, dvrr_stack+11398, dvrr_stack+1549);
 tmp = dvrr_stack + 18243;
 target_ptr = Libderiv->deriv_classes[4][4][2];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,21,dvrr_stack+18468, dvrr_stack+12994, dvrr_stack+3172);
 tmp = dvrr_stack + 18468;
 target_ptr = Libderiv->deriv_classes[4][5][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_g(Data,28,dvrr_stack+23735, dvrr_stack+15248, dvrr_stack+5436);
 tmp = dvrr_stack + 23735;
 target_ptr = Libderiv->deriv_classes[4][6][2];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,10,dvrr_stack+18783, dvrr_stack+17528, dvrr_stack+829);
 tmp = dvrr_stack + 18783;
 target_ptr = Libderiv->deriv_classes[5][3][2];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,15,dvrr_stack+24155, dvrr_stack+7291, dvrr_stack+1699);
 tmp = dvrr_stack + 24155;
 target_ptr = Libderiv->deriv_classes[5][4][2];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,21,dvrr_stack+24470, dvrr_stack+8152, dvrr_stack+3382);
 tmp = dvrr_stack + 24470;
 target_ptr = Libderiv->deriv_classes[5][5][2];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 0 1  0 0 0  0 0 0  0 0 0 */
 deriv_build_AZ_h(Data,28,dvrr_stack+24911, dvrr_stack+2374, dvrr_stack+5716);
 tmp = dvrr_stack + 24911;
 target_ptr = Libderiv->deriv_classes[5][6][2];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,10,dvrr_stack+9055, dvrr_stack+10723, dvrr_stack+729);
 tmp = dvrr_stack + 9055;
 target_ptr = Libderiv->deriv_classes[4][3][1];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,15,dvrr_stack+25499, dvrr_stack+11398, dvrr_stack+1549);
 tmp = dvrr_stack + 25499;
 target_ptr = Libderiv->deriv_classes[4][4][1];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,21,dvrr_stack+25724, dvrr_stack+12994, dvrr_stack+3172);
 tmp = dvrr_stack + 25724;
 target_ptr = Libderiv->deriv_classes[4][5][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_g(Data,28,dvrr_stack+26039, dvrr_stack+15248, dvrr_stack+5436);
 tmp = dvrr_stack + 26039;
 target_ptr = Libderiv->deriv_classes[4][6][1];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,10,dvrr_stack+26459, dvrr_stack+17528, dvrr_stack+829);
 tmp = dvrr_stack + 26459;
 target_ptr = Libderiv->deriv_classes[5][3][1];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,15,dvrr_stack+26669, dvrr_stack+7291, dvrr_stack+1699);
 tmp = dvrr_stack + 26669;
 target_ptr = Libderiv->deriv_classes[5][4][1];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,21,dvrr_stack+26984, dvrr_stack+8152, dvrr_stack+3382);
 tmp = dvrr_stack + 26984;
 target_ptr = Libderiv->deriv_classes[5][5][1];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 0 1 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AY_h(Data,28,dvrr_stack+27425, dvrr_stack+2374, dvrr_stack+5716);
 tmp = dvrr_stack + 27425;
 target_ptr = Libderiv->deriv_classes[5][6][1];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,10,dvrr_stack+28013, dvrr_stack+10723, dvrr_stack+729);
 tmp = dvrr_stack + 28013;
 target_ptr = Libderiv->deriv_classes[4][3][0];
 for(i=0;i<150;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,15,dvrr_stack+28163, dvrr_stack+11398, dvrr_stack+1549);
 tmp = dvrr_stack + 28163;
 target_ptr = Libderiv->deriv_classes[4][4][0];
 for(i=0;i<225;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,21,dvrr_stack+11374, dvrr_stack+12994, dvrr_stack+3172);
 tmp = dvrr_stack + 11374;
 target_ptr = Libderiv->deriv_classes[4][5][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (4 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_g(Data,28,dvrr_stack+12973, dvrr_stack+15248, dvrr_stack+5436);
 tmp = dvrr_stack + 12973;
 target_ptr = Libderiv->deriv_classes[4][6][0];
 for(i=0;i<420;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 3 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,10,dvrr_stack+10723, dvrr_stack+17528, dvrr_stack+829);
 tmp = dvrr_stack + 10723;
 target_ptr = Libderiv->deriv_classes[5][3][0];
 for(i=0;i<210;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 4 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,15,dvrr_stack+17387, dvrr_stack+7291, dvrr_stack+1699);
 tmp = dvrr_stack + 17387;
 target_ptr = Libderiv->deriv_classes[5][4][0];
 for(i=0;i<315;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 5 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,21,dvrr_stack+15118, dvrr_stack+8152, dvrr_stack+3382);
 tmp = dvrr_stack + 15118;
 target_ptr = Libderiv->deriv_classes[5][5][0];
 for(i=0;i<441;i++)
   target_ptr[i] += tmp[i];

 /* compute (5 0 | 6 0) m=0 deriv level 1 */
 /* deriv_ind: 1 0 0  0 0 0  0 0 0  0 0 0 */
 deriv_build_AX_h(Data,28,dvrr_stack+8086, dvrr_stack+2374, dvrr_stack+5716);
 tmp = dvrr_stack + 8086;
 target_ptr = Libderiv->deriv_classes[5][6][0];
 for(i=0;i<588;i++)
   target_ptr[i] += tmp[i];


}

