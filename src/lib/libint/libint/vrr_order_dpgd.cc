#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (dp|gd) integrals */

void vrr_order_dpgd(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+0, vrr_stack+3, Data->F+3, Data->F+4, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+0, Data->F+2, Data->F+3, NULL);
 _BUILD_p0d0(Data,vrr_stack+21, vrr_stack+15, vrr_stack+6, NULL, NULL, vrr_stack+0);
 _BUILD_00f0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+6, vrr_stack+12, vrr_stack+0, NULL);
 _BUILD_00p0(Data,vrr_stack+49, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+52, vrr_stack+49, vrr_stack+12, Data->F+1, Data->F+2, NULL);
 _BUILD_00f0(Data,vrr_stack+58, vrr_stack+52, vrr_stack+15, vrr_stack+49, vrr_stack+12, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+68, vrr_stack+3, vrr_stack+12, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+74, vrr_stack+6, vrr_stack+68, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+84, vrr_stack+39, vrr_stack+74, NULL, NULL, vrr_stack+6);
 _BUILD_p0f0(Data,vrr_stack+114, vrr_stack+58, vrr_stack+39, NULL, NULL, vrr_stack+15);
 _BUILD_d0f0(Data,vrr_stack+144, vrr_stack+114, vrr_stack+84, vrr_stack+58, vrr_stack+39, vrr_stack+21);
 _BUILD_00g0(Data,vrr_stack+21, vrr_stack+39, vrr_stack+74, vrr_stack+15, vrr_stack+6, NULL);
 _BUILD_00g0(Data,vrr_stack+204, vrr_stack+58, vrr_stack+39, vrr_stack+52, vrr_stack+15, NULL);
 _BUILD_p0g0(Data,vrr_stack+219, vrr_stack+204, vrr_stack+21, NULL, NULL, vrr_stack+39);
 _BUILD_00p0(Data,vrr_stack+36, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+36, vrr_stack+49, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+264, vrr_stack+15, vrr_stack+52, vrr_stack+36, vrr_stack+49, NULL);
 _BUILD_00g0(Data,vrr_stack+274, vrr_stack+264, vrr_stack+58, vrr_stack+15, vrr_stack+52, NULL);
 _BUILD_p0g0(Data,vrr_stack+289, vrr_stack+274, vrr_stack+204, NULL, NULL, vrr_stack+58);
 _BUILD_00p0(Data,vrr_stack+36, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+36, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+334, vrr_stack+68, vrr_stack+15, vrr_stack+3, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+344, vrr_stack+74, vrr_stack+334, vrr_stack+6, vrr_stack+68, NULL);
 _BUILD_p0g0(Data,vrr_stack+359, vrr_stack+21, vrr_stack+344, NULL, NULL, vrr_stack+74);
 _BUILD_d0g0(Data,vrr_stack+404, vrr_stack+219, vrr_stack+359, vrr_stack+204, vrr_stack+21, vrr_stack+84);
 _BUILD_d0g0(Data,vrr_stack+494, vrr_stack+289, vrr_stack+219, vrr_stack+274, vrr_stack+204, vrr_stack+114);
   tmp = vrr_stack + 494;
   target_ptr = Libint->vrr_classes[2][4];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00h0(Data,vrr_stack+84, vrr_stack+21, vrr_stack+344, vrr_stack+39, vrr_stack+74, NULL);
 _BUILD_00h0(Data,vrr_stack+105, vrr_stack+204, vrr_stack+21, vrr_stack+58, vrr_stack+39, NULL);
 _BUILD_p0h0(Data,vrr_stack+584, vrr_stack+105, vrr_stack+84, NULL, NULL, vrr_stack+21);
 _BUILD_00h0(Data,vrr_stack+647, vrr_stack+274, vrr_stack+204, vrr_stack+264, vrr_stack+58, NULL);
 _BUILD_p0h0(Data,vrr_stack+668, vrr_stack+647, vrr_stack+105, NULL, NULL, vrr_stack+204);
 _BUILD_00p0(Data,vrr_stack+264, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+267, vrr_stack+36, vrr_stack+264, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+267, vrr_stack+12, vrr_stack+36, NULL);
 _BUILD_00g0(Data,vrr_stack+0, vrr_stack+334, vrr_stack+39, vrr_stack+68, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+49, vrr_stack+344, vrr_stack+0, vrr_stack+74, vrr_stack+334, NULL);
 _BUILD_p0h0(Data,vrr_stack+731, vrr_stack+84, vrr_stack+49, NULL, NULL, vrr_stack+344);
 _BUILD_d0h0(Data,vrr_stack+794, vrr_stack+584, vrr_stack+731, vrr_stack+105, vrr_stack+84, vrr_stack+359);
 _BUILD_d0h0(Data,vrr_stack+920, vrr_stack+668, vrr_stack+584, vrr_stack+647, vrr_stack+105, vrr_stack+219);
   tmp = vrr_stack + 920;
   target_ptr = Libint->vrr_classes[2][5];
   for(i=0;i<126;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00i0(Data,vrr_stack+359, vrr_stack+84, vrr_stack+49, vrr_stack+21, vrr_stack+344, NULL);
 _BUILD_00i0(Data,vrr_stack+1046, vrr_stack+105, vrr_stack+84, vrr_stack+204, vrr_stack+21, NULL);
 _BUILD_p0i0(Data,vrr_stack+1074, vrr_stack+1046, vrr_stack+359, NULL, NULL, vrr_stack+84);
 _BUILD_00i0(Data,vrr_stack+70, vrr_stack+647, vrr_stack+105, vrr_stack+274, vrr_stack+204, NULL);
 _BUILD_p0i0(Data,vrr_stack+1158, vrr_stack+70, vrr_stack+1046, NULL, NULL, vrr_stack+105);
 _BUILD_00p0(Data,vrr_stack+204, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+207, vrr_stack+264, vrr_stack+204, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+647, vrr_stack+267, vrr_stack+207, vrr_stack+36, vrr_stack+264, NULL);
 _BUILD_00g0(Data,vrr_stack+204, vrr_stack+39, vrr_stack+647, vrr_stack+15, vrr_stack+267, NULL);
 _BUILD_00h0(Data,vrr_stack+647, vrr_stack+0, vrr_stack+204, vrr_stack+334, vrr_stack+39, NULL);
 _BUILD_00i0(Data,vrr_stack+15, vrr_stack+49, vrr_stack+647, vrr_stack+344, vrr_stack+0, NULL);
 _BUILD_p0i0(Data,vrr_stack+1242, vrr_stack+359, vrr_stack+15, NULL, NULL, vrr_stack+49);
 _BUILD_d0i0(Data,vrr_stack+1326, vrr_stack+1074, vrr_stack+1242, vrr_stack+1046, vrr_stack+359, vrr_stack+731);
 _BUILD_d0i0(Data,vrr_stack+1494, vrr_stack+1158, vrr_stack+1074, vrr_stack+70, vrr_stack+1046, vrr_stack+584);
   tmp = vrr_stack + 1494;
   target_ptr = Libint->vrr_classes[2][6];
   for(i=0;i<168;i++)
     target_ptr[i] += tmp[i];
 _BUILD_f0g0(Data,vrr_stack+1662, vrr_stack+494, vrr_stack+404, vrr_stack+289, vrr_stack+219, vrr_stack+144);
   tmp = vrr_stack + 1662;
   target_ptr = Libint->vrr_classes[3][4];
   for(i=0;i<150;i++)
     target_ptr[i] += tmp[i];
 _BUILD_f0h0(Data,vrr_stack+0, vrr_stack+920, vrr_stack+794, vrr_stack+668, vrr_stack+584, vrr_stack+404);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[3][5];
   for(i=0;i<210;i++)
     target_ptr[i] += tmp[i];
 _BUILD_f0i0(Data,vrr_stack+210, vrr_stack+1494, vrr_stack+1326, vrr_stack+1158, vrr_stack+1074, vrr_stack+794);
   tmp = vrr_stack + 210;
   target_ptr = Libint->vrr_classes[3][6];
   for(i=0;i<280;i++)
     target_ptr[i] += tmp[i];

}

