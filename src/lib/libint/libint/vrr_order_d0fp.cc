#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (d0|fp) integrals */

void vrr_order_d0fp(Libint_t * Libint, prim_data *Data)
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
 _BUILD_00d0(Data,vrr_stack+68, vrr_stack+3, vrr_stack+12, Data->F+3, Data->F+4, NULL);
 _BUILD_00f0(Data,vrr_stack+74, vrr_stack+6, vrr_stack+68, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+84, vrr_stack+39, vrr_stack+74, NULL, NULL, vrr_stack+6);
 _BUILD_p0f0(Data,vrr_stack+114, vrr_stack+58, vrr_stack+39, NULL, NULL, vrr_stack+15);
 _BUILD_00g0(Data,vrr_stack+144, vrr_stack+39, vrr_stack+74, vrr_stack+15, vrr_stack+6, NULL);
 _BUILD_00g0(Data,vrr_stack+159, vrr_stack+58, vrr_stack+39, vrr_stack+52, vrr_stack+15, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+0, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+174, vrr_stack+68, vrr_stack+15, vrr_stack+3, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+184, vrr_stack+74, vrr_stack+174, vrr_stack+6, vrr_stack+68, NULL);
 _BUILD_p0g0(Data,vrr_stack+199, vrr_stack+144, vrr_stack+184, NULL, NULL, vrr_stack+74);
 _BUILD_p0g0(Data,vrr_stack+244, vrr_stack+159, vrr_stack+144, NULL, NULL, vrr_stack+39);
 _BUILD_d0f0(Data,vrr_stack+289, vrr_stack+114, vrr_stack+84, vrr_stack+58, vrr_stack+39, vrr_stack+21);
   tmp = vrr_stack + 289;
   target_ptr = Libint->vrr_classes[2][3];
   for(i=0;i<60;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0g0(Data,vrr_stack+349, vrr_stack+244, vrr_stack+199, vrr_stack+159, vrr_stack+144, vrr_stack+84);
   tmp = vrr_stack + 349;
   target_ptr = Libint->vrr_classes[2][4];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];

}

