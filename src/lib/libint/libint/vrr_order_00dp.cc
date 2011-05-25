#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (00|dp) integrals */

void vrr_order_00dp(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+6, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+9, vrr_stack+0, vrr_stack+6, Data->F+1, Data->F+2, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+3, vrr_stack+0, Data->F+0, Data->F+1, NULL);
   tmp = vrr_stack + 15;
   target_ptr = Libint->vrr_classes[0][2];
   for(i=0;i<6;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+15, vrr_stack+9, vrr_stack+3, vrr_stack+0, NULL);
   tmp = vrr_stack + 21;
   target_ptr = Libint->vrr_classes[0][3];
   for(i=0;i<10;i++)
     target_ptr[i] += tmp[i];

}

