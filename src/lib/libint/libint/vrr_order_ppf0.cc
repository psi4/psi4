#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (pp|f0) integrals */

void vrr_order_ppf0(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+0, vrr_stack+3, Data->F+2, Data->F+3, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+0, Data->F+1, Data->F+2, NULL);
 _BUILD_p0d0(Data,vrr_stack+21, vrr_stack+15, vrr_stack+6, NULL, NULL, vrr_stack+0);
 _BUILD_00f0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+6, vrr_stack+12, vrr_stack+0, NULL);
 _BUILD_00p0(Data,vrr_stack+49, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+52, vrr_stack+49, vrr_stack+12, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+58, vrr_stack+52, vrr_stack+15, vrr_stack+49, vrr_stack+12, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+49, vrr_stack+3, vrr_stack+12, Data->F+3, Data->F+4, NULL);
 _BUILD_00f0(Data,vrr_stack+68, vrr_stack+6, vrr_stack+49, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+78, vrr_stack+39, vrr_stack+68, NULL, NULL, vrr_stack+6);
 _BUILD_p0f0(Data,vrr_stack+108, vrr_stack+58, vrr_stack+39, NULL, NULL, vrr_stack+15);
   tmp = vrr_stack + 108;
   target_ptr = Libint->vrr_classes[1][3];
   for(i=0;i<30;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0f0(Data,vrr_stack+138, vrr_stack+108, vrr_stack+78, vrr_stack+58, vrr_stack+39, vrr_stack+21);
   tmp = vrr_stack + 138;
   target_ptr = Libint->vrr_classes[2][3];
   for(i=0;i<60;i++)
     target_ptr[i] += tmp[i];

}

