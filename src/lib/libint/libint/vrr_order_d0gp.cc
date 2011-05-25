#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (d0|gp) integrals */

void vrr_order_d0gp(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, Data->F+2, Data->F+3, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+0, vrr_stack+12, Data->F+3, Data->F+4, NULL);
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+6, vrr_stack+15, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00p0(Data,vrr_stack+31, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+34, vrr_stack+31, vrr_stack+3, Data->F+1, Data->F+2, NULL);
 _BUILD_00f0(Data,vrr_stack+40, vrr_stack+34, vrr_stack+6, vrr_stack+31, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+50, vrr_stack+40, vrr_stack+21, NULL, NULL, vrr_stack+6);
 _BUILD_00g0(Data,vrr_stack+80, vrr_stack+40, vrr_stack+21, vrr_stack+34, vrr_stack+6, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+95, vrr_stack+3, vrr_stack+31, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+101, vrr_stack+95, vrr_stack+34, vrr_stack+3, vrr_stack+31, NULL);
 _BUILD_00g0(Data,vrr_stack+111, vrr_stack+101, vrr_stack+40, vrr_stack+95, vrr_stack+34, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+95, vrr_stack+12, vrr_stack+3, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+126, vrr_stack+15, vrr_stack+95, vrr_stack+0, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+136, vrr_stack+21, vrr_stack+126, vrr_stack+6, vrr_stack+15, NULL);
 _BUILD_p0g0(Data,vrr_stack+151, vrr_stack+80, vrr_stack+136, NULL, NULL, vrr_stack+21);
 _BUILD_p0g0(Data,vrr_stack+196, vrr_stack+111, vrr_stack+80, NULL, NULL, vrr_stack+40);
 _BUILD_00h0(Data,vrr_stack+241, vrr_stack+80, vrr_stack+136, vrr_stack+40, vrr_stack+21, NULL);
 _BUILD_00h0(Data,vrr_stack+262, vrr_stack+111, vrr_stack+80, vrr_stack+101, vrr_stack+40, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+101, vrr_stack+95, vrr_stack+6, vrr_stack+12, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+0, vrr_stack+126, vrr_stack+101, vrr_stack+15, vrr_stack+95, NULL);
 _BUILD_00h0(Data,vrr_stack+283, vrr_stack+136, vrr_stack+0, vrr_stack+21, vrr_stack+126, NULL);
 _BUILD_p0h0(Data,vrr_stack+304, vrr_stack+241, vrr_stack+283, NULL, NULL, vrr_stack+136);
 _BUILD_p0h0(Data,vrr_stack+367, vrr_stack+262, vrr_stack+241, NULL, NULL, vrr_stack+80);
 _BUILD_d0g0(Data,vrr_stack+430, vrr_stack+196, vrr_stack+151, vrr_stack+111, vrr_stack+80, vrr_stack+50);
   tmp = vrr_stack + 430;
   target_ptr = Libint->vrr_classes[2][4];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0h0(Data,vrr_stack+0, vrr_stack+367, vrr_stack+304, vrr_stack+262, vrr_stack+241, vrr_stack+151);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[2][5];
   for(i=0;i<126;i++)
     target_ptr[i] += tmp[i];

}

