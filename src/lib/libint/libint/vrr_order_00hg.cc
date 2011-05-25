#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (00|hg) integrals */

void vrr_order_00hg(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+6, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+9, vrr_stack+0, vrr_stack+6, Data->F+4, Data->F+5, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+3, vrr_stack+0, Data->F+3, Data->F+4, NULL);
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+15, vrr_stack+9, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00p0(Data,vrr_stack+31, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+34, vrr_stack+31, vrr_stack+3, Data->F+2, Data->F+3, NULL);
 _BUILD_00f0(Data,vrr_stack+40, vrr_stack+34, vrr_stack+15, vrr_stack+31, vrr_stack+3, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+50, vrr_stack+6, vrr_stack+3, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+56, vrr_stack+9, vrr_stack+50, vrr_stack+0, vrr_stack+6, NULL);
 _BUILD_00g0(Data,vrr_stack+66, vrr_stack+21, vrr_stack+56, vrr_stack+15, vrr_stack+9, NULL);
 _BUILD_00g0(Data,vrr_stack+81, vrr_stack+40, vrr_stack+21, vrr_stack+34, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+96, vrr_stack+81, vrr_stack+66, vrr_stack+40, vrr_stack+21, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+0, vrr_stack+31, Data->F+1, Data->F+2, NULL);
 _BUILD_00f0(Data,vrr_stack+117, vrr_stack+15, vrr_stack+34, vrr_stack+0, vrr_stack+31, NULL);
 _BUILD_00g0(Data,vrr_stack+127, vrr_stack+117, vrr_stack+40, vrr_stack+15, vrr_stack+34, NULL);
 _BUILD_00h0(Data,vrr_stack+142, vrr_stack+127, vrr_stack+81, vrr_stack+117, vrr_stack+40, NULL);
 _BUILD_00p0(Data,vrr_stack+31, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+34, vrr_stack+3, vrr_stack+31, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+40, vrr_stack+50, vrr_stack+34, vrr_stack+6, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+163, vrr_stack+56, vrr_stack+40, vrr_stack+9, vrr_stack+50, NULL);
 _BUILD_00h0(Data,vrr_stack+178, vrr_stack+66, vrr_stack+163, vrr_stack+21, vrr_stack+56, NULL);
 _BUILD_00i0(Data,vrr_stack+199, vrr_stack+96, vrr_stack+178, vrr_stack+81, vrr_stack+66, NULL);
 _BUILD_00i0(Data,vrr_stack+227, vrr_stack+142, vrr_stack+96, vrr_stack+127, vrr_stack+81, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+255, vrr_stack+227, vrr_stack+199, vrr_stack+142, vrr_stack+96, NULL);
 _BUILD_00p0(Data,vrr_stack+81, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+84, vrr_stack+81, vrr_stack+0, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+84, vrr_stack+15, vrr_stack+81, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+291, vrr_stack+21, vrr_stack+117, vrr_stack+84, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+306, vrr_stack+291, vrr_stack+127, vrr_stack+21, vrr_stack+117, NULL);
   tmp = vrr_stack + 306;
   target_ptr = Libint->vrr_classes[0][5];
   for(i=0;i<21;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00i0(Data,vrr_stack+327, vrr_stack+306, vrr_stack+142, vrr_stack+291, vrr_stack+127, NULL);
   tmp = vrr_stack + 327;
   target_ptr = Libint->vrr_classes[0][6];
   for(i=0;i<28;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+355, vrr_stack+327, vrr_stack+227, vrr_stack+306, vrr_stack+142, NULL);
   tmp = vrr_stack + 355;
   target_ptr = Libint->vrr_classes[0][7];
   for(i=0;i<36;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+291, vrr_stack+31, vrr_stack+0, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+297, vrr_stack+34, vrr_stack+291, vrr_stack+3, vrr_stack+31, NULL);
 _BUILD_00g0(Data,vrr_stack+81, vrr_stack+40, vrr_stack+297, vrr_stack+50, vrr_stack+34, NULL);
 _BUILD_00h0(Data,vrr_stack+291, vrr_stack+163, vrr_stack+81, vrr_stack+56, vrr_stack+40, NULL);
 _BUILD_00i0(Data,vrr_stack+0, vrr_stack+178, vrr_stack+291, vrr_stack+66, vrr_stack+163, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+291, vrr_stack+199, vrr_stack+0, vrr_stack+96, vrr_stack+178, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+0, vrr_stack+255, vrr_stack+291, vrr_stack+227, vrr_stack+199, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+45, vrr_stack+355, vrr_stack+255, vrr_stack+327, vrr_stack+227, NULL);
   tmp = vrr_stack + 45;
   target_ptr = Libint->vrr_classes[0][8];
   for(i=0;i<45;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+291, vrr_stack+45, vrr_stack+0, vrr_stack+355, vrr_stack+255, NULL);
   tmp = vrr_stack + 291;
   target_ptr = Libint->vrr_classes[0][9];
   for(i=0;i<55;i++)
     target_ptr[i] += tmp[i];

}

