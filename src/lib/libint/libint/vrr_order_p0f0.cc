#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (p0|f0) integrals */

void vrr_order_p0f0(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, Data->F+1, Data->F+2, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+0, vrr_stack+12, Data->F+2, Data->F+3, NULL);
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+6, vrr_stack+15, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+12, vrr_stack+0, vrr_stack+3, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+31, vrr_stack+12, vrr_stack+6, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+41, vrr_stack+31, vrr_stack+21, NULL, NULL, vrr_stack+6);
   tmp = vrr_stack + 41;
   target_ptr = Libint->vrr_classes[1][3];
   for(i=0;i<30;i++)
     target_ptr[i] += tmp[i];

}

