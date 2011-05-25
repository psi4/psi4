#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (p0|p0) integrals */

void vrr_order_p0p0(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, NULL, NULL, Data->F+1);
   tmp = vrr_stack + 6;
   target_ptr = Libint->vrr_classes[1][1];
   for(i=0;i<9;i++)
     target_ptr[i] += tmp[i];

}

