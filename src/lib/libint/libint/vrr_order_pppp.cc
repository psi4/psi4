#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (pp|pp) integrals */

void vrr_order_pppp(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_p000(Data,vrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+6, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+9, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+12, vrr_stack+3, vrr_stack+9, NULL, NULL, Data->F+2);
 _BUILD_p0p0(Data,vrr_stack+21, vrr_stack+6, vrr_stack+3, NULL, NULL, Data->F+1);
   tmp = vrr_stack + 21;
   target_ptr = Libint->vrr_classes[1][1];
   for(i=0;i<9;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00d0(Data,vrr_stack+30, vrr_stack+3, vrr_stack+9, Data->F+1, Data->F+2, NULL);
 _BUILD_00d0(Data,vrr_stack+36, vrr_stack+6, vrr_stack+3, Data->F+0, Data->F+1, NULL);
 _BUILD_00p0(Data,vrr_stack+42, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+45, vrr_stack+9, vrr_stack+42, Data->F+2, Data->F+3, NULL);
 _BUILD_p0d0(Data,vrr_stack+51, vrr_stack+30, vrr_stack+45, NULL, NULL, vrr_stack+9);
 _BUILD_p0d0(Data,vrr_stack+69, vrr_stack+36, vrr_stack+30, NULL, NULL, vrr_stack+3);
   tmp = vrr_stack + 69;
   target_ptr = Libint->vrr_classes[1][2];
   for(i=0;i<18;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0p0(Data,vrr_stack+87, vrr_stack+21, vrr_stack+12, vrr_stack+6, vrr_stack+3, vrr_stack+0);
   tmp = vrr_stack + 87;
   target_ptr = Libint->vrr_classes[2][1];
   for(i=0;i<18;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0d0(Data,vrr_stack+105, vrr_stack+69, vrr_stack+51, vrr_stack+36, vrr_stack+30, vrr_stack+12);
   tmp = vrr_stack + 105;
   target_ptr = Libint->vrr_classes[2][2];
   for(i=0;i<36;i++)
     target_ptr[i] += tmp[i];

}

