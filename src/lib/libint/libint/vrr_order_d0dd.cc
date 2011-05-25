#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (d0|dd) integrals */

void vrr_order_d0dd(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, NULL, NULL, Data->F+2);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+3, vrr_stack+0, Data->F+1, Data->F+2, NULL);
 _BUILD_00p0(Data,vrr_stack+21, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+24, vrr_stack+21, vrr_stack+3, Data->F+0, Data->F+1, NULL);
 _BUILD_00p0(Data,vrr_stack+30, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+33, vrr_stack+0, vrr_stack+30, Data->F+2, Data->F+3, NULL);
 _BUILD_p0d0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+33, NULL, NULL, vrr_stack+0);
 _BUILD_p0d0(Data,vrr_stack+57, vrr_stack+24, vrr_stack+15, NULL, NULL, vrr_stack+3);
 _BUILD_00f0(Data,vrr_stack+75, vrr_stack+15, vrr_stack+33, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00f0(Data,vrr_stack+85, vrr_stack+24, vrr_stack+15, vrr_stack+21, vrr_stack+3, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+95, vrr_stack+30, vrr_stack+3, Data->F+3, Data->F+4, NULL);
 _BUILD_00f0(Data,vrr_stack+101, vrr_stack+33, vrr_stack+95, vrr_stack+0, vrr_stack+30, NULL);
 _BUILD_p0f0(Data,vrr_stack+111, vrr_stack+75, vrr_stack+101, NULL, NULL, vrr_stack+33);
 _BUILD_p0f0(Data,vrr_stack+141, vrr_stack+85, vrr_stack+75, NULL, NULL, vrr_stack+15);
 _BUILD_00g0(Data,vrr_stack+171, vrr_stack+75, vrr_stack+101, vrr_stack+15, vrr_stack+33, NULL);
 _BUILD_00g0(Data,vrr_stack+186, vrr_stack+85, vrr_stack+75, vrr_stack+24, vrr_stack+15, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+201, vrr_stack+3, vrr_stack+0, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+207, vrr_stack+95, vrr_stack+201, vrr_stack+30, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+217, vrr_stack+101, vrr_stack+207, vrr_stack+33, vrr_stack+95, NULL);
 _BUILD_p0g0(Data,vrr_stack+232, vrr_stack+171, vrr_stack+217, NULL, NULL, vrr_stack+101);
 _BUILD_p0g0(Data,vrr_stack+277, vrr_stack+186, vrr_stack+171, NULL, NULL, vrr_stack+75);
 _BUILD_d0d0(Data,vrr_stack+322, vrr_stack+57, vrr_stack+39, vrr_stack+24, vrr_stack+15, vrr_stack+6);
   tmp = vrr_stack + 322;
   target_ptr = Libint->vrr_classes[2][2];
   for(i=0;i<36;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0f0(Data,vrr_stack+358, vrr_stack+141, vrr_stack+111, vrr_stack+85, vrr_stack+75, vrr_stack+39);
   tmp = vrr_stack + 358;
   target_ptr = Libint->vrr_classes[2][3];
   for(i=0;i<60;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0g0(Data,vrr_stack+0, vrr_stack+277, vrr_stack+232, vrr_stack+186, vrr_stack+171, vrr_stack+111);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[2][4];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];

}

