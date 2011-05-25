#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (dd|dd) integrals */

void vrr_order_dddd(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_p000(Data,vrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_p000(Data,vrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_d000(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, Data->F+2, Data->F+3, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+15, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+18, vrr_stack+15, vrr_stack+12, NULL, NULL, Data->F+3);
 _BUILD_00p0(Data,vrr_stack+27, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+30, vrr_stack+27, vrr_stack+15, NULL, NULL, Data->F+2);
 _BUILD_00p0(Data,vrr_stack+39, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+42, vrr_stack+12, vrr_stack+39, NULL, NULL, Data->F+4);
 _BUILD_d0p0(Data,vrr_stack+51, vrr_stack+18, vrr_stack+42, vrr_stack+15, vrr_stack+12, vrr_stack+0);
 _BUILD_d0p0(Data,vrr_stack+69, vrr_stack+30, vrr_stack+18, vrr_stack+27, vrr_stack+15, vrr_stack+3);
 _BUILD_f0p0(Data,vrr_stack+87, vrr_stack+69, vrr_stack+51, vrr_stack+30, vrr_stack+18, vrr_stack+6);
 _BUILD_00d0(Data,vrr_stack+0, vrr_stack+15, vrr_stack+12, Data->F+2, Data->F+3, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+27, vrr_stack+15, Data->F+1, Data->F+2, NULL);
 _BUILD_00d0(Data,vrr_stack+117, vrr_stack+12, vrr_stack+39, Data->F+3, Data->F+4, NULL);
 _BUILD_p0d0(Data,vrr_stack+123, vrr_stack+0, vrr_stack+117, NULL, NULL, vrr_stack+12);
 _BUILD_p0d0(Data,vrr_stack+141, vrr_stack+6, vrr_stack+0, NULL, NULL, vrr_stack+15);
 _BUILD_d0d0(Data,vrr_stack+159, vrr_stack+141, vrr_stack+123, vrr_stack+6, vrr_stack+0, vrr_stack+18);
 _BUILD_00p0(Data,vrr_stack+18, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+21, vrr_stack+18, vrr_stack+27, Data->F+0, Data->F+1, NULL);
 _BUILD_p0d0(Data,vrr_stack+195, vrr_stack+21, vrr_stack+6, NULL, NULL, vrr_stack+27);
 _BUILD_d0d0(Data,vrr_stack+213, vrr_stack+195, vrr_stack+141, vrr_stack+21, vrr_stack+6, vrr_stack+30);
   tmp = vrr_stack + 213;
   target_ptr = Libint->vrr_classes[2][2];
   for(i=0;i<36;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+30, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+33, vrr_stack+39, vrr_stack+30, Data->F+4, Data->F+5, NULL);
 _BUILD_p0d0(Data,vrr_stack+249, vrr_stack+117, vrr_stack+33, NULL, NULL, vrr_stack+39);
 _BUILD_d0d0(Data,vrr_stack+267, vrr_stack+123, vrr_stack+249, vrr_stack+0, vrr_stack+117, vrr_stack+42);
 _BUILD_f0d0(Data,vrr_stack+303, vrr_stack+159, vrr_stack+267, vrr_stack+141, vrr_stack+123, vrr_stack+51);
 _BUILD_f0d0(Data,vrr_stack+363, vrr_stack+213, vrr_stack+159, vrr_stack+195, vrr_stack+141, vrr_stack+69);
   tmp = vrr_stack + 363;
   target_ptr = Libint->vrr_classes[3][2];
   for(i=0;i<60;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00f0(Data,vrr_stack+195, vrr_stack+0, vrr_stack+117, vrr_stack+15, vrr_stack+12, NULL);
 _BUILD_00f0(Data,vrr_stack+42, vrr_stack+6, vrr_stack+0, vrr_stack+27, vrr_stack+15, NULL);
 _BUILD_00f0(Data,vrr_stack+52, vrr_stack+117, vrr_stack+33, vrr_stack+12, vrr_stack+39, NULL);
 _BUILD_p0f0(Data,vrr_stack+423, vrr_stack+195, vrr_stack+52, NULL, NULL, vrr_stack+117);
 _BUILD_p0f0(Data,vrr_stack+453, vrr_stack+42, vrr_stack+195, NULL, NULL, vrr_stack+0);
 _BUILD_d0f0(Data,vrr_stack+483, vrr_stack+453, vrr_stack+423, vrr_stack+42, vrr_stack+195, vrr_stack+123);
 _BUILD_00f0(Data,vrr_stack+123, vrr_stack+21, vrr_stack+6, vrr_stack+18, vrr_stack+27, NULL);
 _BUILD_p0f0(Data,vrr_stack+543, vrr_stack+123, vrr_stack+42, NULL, NULL, vrr_stack+6);
 _BUILD_d0f0(Data,vrr_stack+573, vrr_stack+543, vrr_stack+453, vrr_stack+123, vrr_stack+42, vrr_stack+141);
   tmp = vrr_stack + 573;
   target_ptr = Libint->vrr_classes[2][3];
   for(i=0;i<60;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+27, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+133, vrr_stack+30, vrr_stack+27, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+139, vrr_stack+33, vrr_stack+133, vrr_stack+39, vrr_stack+30, NULL);
 _BUILD_p0f0(Data,vrr_stack+633, vrr_stack+52, vrr_stack+139, NULL, NULL, vrr_stack+33);
 _BUILD_d0f0(Data,vrr_stack+663, vrr_stack+423, vrr_stack+633, vrr_stack+195, vrr_stack+52, vrr_stack+249);
 _BUILD_f0f0(Data,vrr_stack+723, vrr_stack+483, vrr_stack+663, vrr_stack+453, vrr_stack+423, vrr_stack+267);
 _BUILD_f0f0(Data,vrr_stack+823, vrr_stack+573, vrr_stack+483, vrr_stack+543, vrr_stack+453, vrr_stack+159);
   tmp = vrr_stack + 823;
   target_ptr = Libint->vrr_classes[3][3];
   for(i=0;i<100;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00g0(Data,vrr_stack+543, vrr_stack+195, vrr_stack+52, vrr_stack+0, vrr_stack+117, NULL);
 _BUILD_00g0(Data,vrr_stack+558, vrr_stack+42, vrr_stack+195, vrr_stack+6, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+249, vrr_stack+52, vrr_stack+139, vrr_stack+117, vrr_stack+33, NULL);
 _BUILD_p0g0(Data,vrr_stack+923, vrr_stack+543, vrr_stack+249, NULL, NULL, vrr_stack+52);
 _BUILD_p0g0(Data,vrr_stack+968, vrr_stack+558, vrr_stack+543, NULL, NULL, vrr_stack+195);
 _BUILD_d0g0(Data,vrr_stack+1013, vrr_stack+968, vrr_stack+923, vrr_stack+558, vrr_stack+543, vrr_stack+423);
 _BUILD_00g0(Data,vrr_stack+423, vrr_stack+123, vrr_stack+42, vrr_stack+21, vrr_stack+6, NULL);
 _BUILD_p0g0(Data,vrr_stack+1103, vrr_stack+423, vrr_stack+558, NULL, NULL, vrr_stack+42);
 _BUILD_d0g0(Data,vrr_stack+1148, vrr_stack+1103, vrr_stack+968, vrr_stack+423, vrr_stack+558, vrr_stack+453);
   tmp = vrr_stack + 1148;
   target_ptr = Libint->vrr_classes[2][4];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+558, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+561, vrr_stack+27, vrr_stack+558, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+149, vrr_stack+133, vrr_stack+561, vrr_stack+30, vrr_stack+27, NULL);
 _BUILD_00g0(Data,vrr_stack+558, vrr_stack+139, vrr_stack+149, vrr_stack+33, vrr_stack+133, NULL);
 _BUILD_p0g0(Data,vrr_stack+423, vrr_stack+249, vrr_stack+558, NULL, NULL, vrr_stack+139);
 _BUILD_d0g0(Data,vrr_stack+1238, vrr_stack+923, vrr_stack+423, vrr_stack+543, vrr_stack+249, vrr_stack+633);
 _BUILD_f0g0(Data,vrr_stack+1328, vrr_stack+1013, vrr_stack+1238, vrr_stack+968, vrr_stack+923, vrr_stack+663);
 _BUILD_f0g0(Data,vrr_stack+1478, vrr_stack+1148, vrr_stack+1013, vrr_stack+1103, vrr_stack+968, vrr_stack+483);
   tmp = vrr_stack + 1478;
   target_ptr = Libint->vrr_classes[3][4];
   for(i=0;i<150;i++)
     target_ptr[i] += tmp[i];
 _BUILD_g0d0(Data,vrr_stack+923, vrr_stack+363, vrr_stack+303, vrr_stack+213, vrr_stack+159, vrr_stack+87);
   tmp = vrr_stack + 923;
   target_ptr = Libint->vrr_classes[4][2];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];
 _BUILD_g0f0(Data,vrr_stack+0, vrr_stack+823, vrr_stack+723, vrr_stack+573, vrr_stack+483, vrr_stack+303);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[4][3];
   for(i=0;i<150;i++)
     target_ptr[i] += tmp[i];
 _BUILD_g0g0(Data,vrr_stack+150, vrr_stack+1478, vrr_stack+1328, vrr_stack+1148, vrr_stack+1013, vrr_stack+723);
   tmp = vrr_stack + 150;
   target_ptr = Libint->vrr_classes[4][4];
   for(i=0;i<225;i++)
     target_ptr[i] += tmp[i];

}

