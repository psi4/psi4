#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (00|gd) integrals */

void vrr_order_00gd(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, Data->F+2, Data->F+3, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+3, Data->F+1, Data->F+2, NULL);
 _BUILD_00p0(Data,vrr_stack+21, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+24, vrr_stack+0, vrr_stack+21, Data->F+3, Data->F+4, NULL);
 _BUILD_00f0(Data,vrr_stack+30, vrr_stack+6, vrr_stack+24, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00f0(Data,vrr_stack+40, vrr_stack+15, vrr_stack+6, vrr_stack+12, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+50, vrr_stack+40, vrr_stack+30, vrr_stack+15, vrr_stack+6, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+65, vrr_stack+3, vrr_stack+12, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+71, vrr_stack+65, vrr_stack+15, vrr_stack+3, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+81, vrr_stack+71, vrr_stack+40, vrr_stack+65, vrr_stack+15, NULL);
   tmp = vrr_stack + 81;
   target_ptr = Libint->vrr_classes[0][4];
   for(i=0;i<15;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+3, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+65, vrr_stack+21, vrr_stack+3, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+96, vrr_stack+24, vrr_stack+65, vrr_stack+0, vrr_stack+21, NULL);
 _BUILD_00g0(Data,vrr_stack+106, vrr_stack+30, vrr_stack+96, vrr_stack+6, vrr_stack+24, NULL);
 _BUILD_00h0(Data,vrr_stack+0, vrr_stack+50, vrr_stack+106, vrr_stack+40, vrr_stack+30, NULL);
 _BUILD_00h0(Data,vrr_stack+96, vrr_stack+81, vrr_stack+50, vrr_stack+71, vrr_stack+40, NULL);
   tmp = vrr_stack + 96;
   target_ptr = Libint->vrr_classes[0][5];
   for(i=0;i<21;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00i0(Data,vrr_stack+117, vrr_stack+96, vrr_stack+0, vrr_stack+81, vrr_stack+50, NULL);
   tmp = vrr_stack + 117;
   target_ptr = Libint->vrr_classes[0][6];
   for(i=0;i<28;i++)
     target_ptr[i] += tmp[i];

}

