#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (f0|g0) integrals */

void vrr_order_f0g0(Libint_t * Libint, prim_data *Data)
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
 _BUILD_00p0(Data,vrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+0, vrr_stack+49, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+36, vrr_stack+15, vrr_stack+52, vrr_stack+0, vrr_stack+49, NULL);
 _BUILD_00g0(Data,vrr_stack+264, vrr_stack+36, vrr_stack+58, vrr_stack+15, vrr_stack+52, NULL);
 _BUILD_p0g0(Data,vrr_stack+279, vrr_stack+264, vrr_stack+204, NULL, NULL, vrr_stack+58);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+0, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+36, vrr_stack+68, vrr_stack+15, vrr_stack+3, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+46, vrr_stack+74, vrr_stack+36, vrr_stack+6, vrr_stack+68, NULL);
 _BUILD_p0g0(Data,vrr_stack+324, vrr_stack+21, vrr_stack+46, NULL, NULL, vrr_stack+74);
 _BUILD_d0g0(Data,vrr_stack+369, vrr_stack+219, vrr_stack+324, vrr_stack+204, vrr_stack+21, vrr_stack+84);
 _BUILD_d0g0(Data,vrr_stack+0, vrr_stack+279, vrr_stack+219, vrr_stack+264, vrr_stack+204, vrr_stack+114);
 _BUILD_f0g0(Data,vrr_stack+459, vrr_stack+0, vrr_stack+369, vrr_stack+279, vrr_stack+219, vrr_stack+144);
   tmp = vrr_stack + 459;
   target_ptr = Libint->vrr_classes[3][4];
   for(i=0;i<150;i++)
     target_ptr[i] += tmp[i];

}

