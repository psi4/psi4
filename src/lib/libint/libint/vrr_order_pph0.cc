#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (pp|h0) integrals */

void vrr_order_pph0(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+6, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+9, vrr_stack+0, vrr_stack+6, Data->F+3, Data->F+4, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+3, vrr_stack+0, Data->F+2, Data->F+3, NULL);
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+15, vrr_stack+9, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00p0(Data,vrr_stack+31, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+34, vrr_stack+6, vrr_stack+31, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+40, vrr_stack+9, vrr_stack+34, vrr_stack+0, vrr_stack+6, NULL);
 _BUILD_00g0(Data,vrr_stack+50, vrr_stack+21, vrr_stack+40, vrr_stack+15, vrr_stack+9, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+65, vrr_stack+0, vrr_stack+3, Data->F+1, Data->F+2, NULL);
 _BUILD_00f0(Data,vrr_stack+71, vrr_stack+65, vrr_stack+15, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+81, vrr_stack+71, vrr_stack+21, vrr_stack+65, vrr_stack+15, NULL);
 _BUILD_p0g0(Data,vrr_stack+96, vrr_stack+81, vrr_stack+50, NULL, NULL, vrr_stack+21);
 _BUILD_00h0(Data,vrr_stack+141, vrr_stack+81, vrr_stack+50, vrr_stack+71, vrr_stack+21, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+3, vrr_stack+0, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+162, vrr_stack+15, vrr_stack+65, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+172, vrr_stack+162, vrr_stack+71, vrr_stack+15, vrr_stack+65, NULL);
 _BUILD_00h0(Data,vrr_stack+187, vrr_stack+172, vrr_stack+81, vrr_stack+162, vrr_stack+71, NULL);
 _BUILD_00p0(Data,vrr_stack+162, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+31, vrr_stack+162, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+162, vrr_stack+34, vrr_stack+15, vrr_stack+6, vrr_stack+31, NULL);
 _BUILD_00g0(Data,vrr_stack+172, vrr_stack+40, vrr_stack+162, vrr_stack+9, vrr_stack+34, NULL);
 _BUILD_00h0(Data,vrr_stack+0, vrr_stack+50, vrr_stack+172, vrr_stack+21, vrr_stack+40, NULL);
 _BUILD_p0h0(Data,vrr_stack+208, vrr_stack+141, vrr_stack+0, NULL, NULL, vrr_stack+50);
 _BUILD_p0h0(Data,vrr_stack+0, vrr_stack+187, vrr_stack+141, NULL, NULL, vrr_stack+81);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[1][5];
   for(i=0;i<63;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0h0(Data,vrr_stack+271, vrr_stack+0, vrr_stack+208, vrr_stack+187, vrr_stack+141, vrr_stack+96);
   tmp = vrr_stack + 271;
   target_ptr = Libint->vrr_classes[2][5];
   for(i=0;i<126;i++)
     target_ptr[i] += tmp[i];

}

