#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (f0|dp) integrals */

void vrr_order_f0dp(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_p000(Data,vrr_stack+0, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+6, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+9, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+12, vrr_stack+3, vrr_stack+9, NULL, NULL, Data->F+3);
 _BUILD_p0p0(Data,vrr_stack+21, vrr_stack+6, vrr_stack+3, NULL, NULL, Data->F+2);
 _BUILD_d0p0(Data,vrr_stack+30, vrr_stack+21, vrr_stack+12, vrr_stack+6, vrr_stack+3, vrr_stack+0);
 _BUILD_00d0(Data,vrr_stack+48, vrr_stack+3, vrr_stack+9, Data->F+2, Data->F+3, NULL);
 _BUILD_00d0(Data,vrr_stack+54, vrr_stack+6, vrr_stack+3, Data->F+1, Data->F+2, NULL);
 _BUILD_p0d0(Data,vrr_stack+60, vrr_stack+54, vrr_stack+48, NULL, NULL, vrr_stack+3);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+78, vrr_stack+0, vrr_stack+6, Data->F+0, Data->F+1, NULL);
 _BUILD_p0d0(Data,vrr_stack+84, vrr_stack+78, vrr_stack+54, NULL, NULL, vrr_stack+6);
 _BUILD_00p0(Data,vrr_stack+102, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+105, vrr_stack+9, vrr_stack+102, Data->F+3, Data->F+4, NULL);
 _BUILD_p0d0(Data,vrr_stack+111, vrr_stack+48, vrr_stack+105, NULL, NULL, vrr_stack+9);
 _BUILD_d0d0(Data,vrr_stack+129, vrr_stack+60, vrr_stack+111, vrr_stack+54, vrr_stack+48, vrr_stack+12);
 _BUILD_d0d0(Data,vrr_stack+165, vrr_stack+84, vrr_stack+60, vrr_stack+78, vrr_stack+54, vrr_stack+21);
 _BUILD_00f0(Data,vrr_stack+12, vrr_stack+48, vrr_stack+105, vrr_stack+3, vrr_stack+9, NULL);
 _BUILD_00f0(Data,vrr_stack+201, vrr_stack+54, vrr_stack+48, vrr_stack+6, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+211, vrr_stack+201, vrr_stack+12, NULL, NULL, vrr_stack+48);
 _BUILD_00f0(Data,vrr_stack+241, vrr_stack+78, vrr_stack+54, vrr_stack+0, vrr_stack+6, NULL);
 _BUILD_p0f0(Data,vrr_stack+251, vrr_stack+241, vrr_stack+201, NULL, NULL, vrr_stack+54);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+3, vrr_stack+102, vrr_stack+0, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+48, vrr_stack+105, vrr_stack+3, vrr_stack+9, vrr_stack+102, NULL);
 _BUILD_p0f0(Data,vrr_stack+281, vrr_stack+12, vrr_stack+48, NULL, NULL, vrr_stack+105);
 _BUILD_d0f0(Data,vrr_stack+311, vrr_stack+211, vrr_stack+281, vrr_stack+201, vrr_stack+12, vrr_stack+111);
 _BUILD_d0f0(Data,vrr_stack+371, vrr_stack+251, vrr_stack+211, vrr_stack+241, vrr_stack+201, vrr_stack+60);
 _BUILD_f0d0(Data,vrr_stack+431, vrr_stack+165, vrr_stack+129, vrr_stack+84, vrr_stack+60, vrr_stack+30);
   tmp = vrr_stack + 431;
   target_ptr = Libint->vrr_classes[3][2];
   for(i=0;i<60;i++)
     target_ptr[i] += tmp[i];
 _BUILD_f0f0(Data,vrr_stack+0, vrr_stack+371, vrr_stack+311, vrr_stack+251, vrr_stack+211, vrr_stack+129);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[3][3];
   for(i=0;i<100;i++)
     target_ptr[i] += tmp[i];

}

