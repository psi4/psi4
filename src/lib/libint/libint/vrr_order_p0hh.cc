#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (p0|hh) integrals */

void vrr_order_p0hh(Libint_t * Libint, prim_data *Data)
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
 _BUILD_00p0(Data,vrr_stack+3, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+65, vrr_stack+21, vrr_stack+3, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+71, vrr_stack+24, vrr_stack+65, vrr_stack+0, vrr_stack+21, NULL);
 _BUILD_00g0(Data,vrr_stack+81, vrr_stack+30, vrr_stack+71, vrr_stack+6, vrr_stack+24, NULL);
 _BUILD_00h0(Data,vrr_stack+96, vrr_stack+50, vrr_stack+81, vrr_stack+40, vrr_stack+30, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+0, vrr_stack+12, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+117, vrr_stack+6, vrr_stack+15, vrr_stack+0, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+127, vrr_stack+117, vrr_stack+40, vrr_stack+6, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+142, vrr_stack+127, vrr_stack+50, vrr_stack+117, vrr_stack+40, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+40, vrr_stack+3, vrr_stack+0, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+117, vrr_stack+65, vrr_stack+40, vrr_stack+21, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+6, vrr_stack+71, vrr_stack+117, vrr_stack+24, vrr_stack+65, NULL);
 _BUILD_00h0(Data,vrr_stack+163, vrr_stack+81, vrr_stack+6, vrr_stack+30, vrr_stack+71, NULL);
 _BUILD_00i0(Data,vrr_stack+184, vrr_stack+96, vrr_stack+163, vrr_stack+50, vrr_stack+81, NULL);
 _BUILD_00i0(Data,vrr_stack+212, vrr_stack+142, vrr_stack+96, vrr_stack+127, vrr_stack+50, NULL);
 _BUILD_00p0(Data,vrr_stack+127, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+130, vrr_stack+0, vrr_stack+127, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+40, vrr_stack+130, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+240, vrr_stack+117, vrr_stack+21, vrr_stack+65, vrr_stack+40, NULL);
 _BUILD_00h0(Data,vrr_stack+255, vrr_stack+6, vrr_stack+240, vrr_stack+71, vrr_stack+117, NULL);
 _BUILD_00i0(Data,vrr_stack+276, vrr_stack+163, vrr_stack+255, vrr_stack+81, vrr_stack+6, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+304, vrr_stack+184, vrr_stack+276, vrr_stack+96, vrr_stack+163, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+340, vrr_stack+212, vrr_stack+184, vrr_stack+142, vrr_stack+96, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+136, vrr_stack+127, vrr_stack+3, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+65, vrr_stack+130, vrr_stack+136, vrr_stack+0, vrr_stack+127, NULL);
 _BUILD_00g0(Data,vrr_stack+75, vrr_stack+21, vrr_stack+65, vrr_stack+40, vrr_stack+130, NULL);
 _BUILD_00h0(Data,vrr_stack+376, vrr_stack+240, vrr_stack+75, vrr_stack+117, vrr_stack+21, NULL);
 _BUILD_00i0(Data,vrr_stack+397, vrr_stack+255, vrr_stack+376, vrr_stack+6, vrr_stack+240, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+425, vrr_stack+276, vrr_stack+397, vrr_stack+163, vrr_stack+255, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+461, vrr_stack+304, vrr_stack+425, vrr_stack+184, vrr_stack+276, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+506, vrr_stack+340, vrr_stack+304, vrr_stack+212, vrr_stack+184, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+9, Data->F+10, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+90, vrr_stack+3, vrr_stack+0, Data->F+8, Data->F+9, NULL);
 _BUILD_00f0(Data,vrr_stack+117, vrr_stack+136, vrr_stack+90, vrr_stack+127, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+6, vrr_stack+65, vrr_stack+117, vrr_stack+130, vrr_stack+136, NULL);
 _BUILD_00h0(Data,vrr_stack+163, vrr_stack+75, vrr_stack+6, vrr_stack+21, vrr_stack+65, NULL);
 _BUILD_00i0(Data,vrr_stack+21, vrr_stack+376, vrr_stack+163, vrr_stack+240, vrr_stack+75, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+551, vrr_stack+397, vrr_stack+21, vrr_stack+255, vrr_stack+376, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+587, vrr_stack+425, vrr_stack+551, vrr_stack+276, vrr_stack+397, NULL);
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+240, vrr_stack+461, vrr_stack+587, vrr_stack+304, vrr_stack+425, NULL);
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+632, vrr_stack+506, vrr_stack+461, vrr_stack+340, vrr_stack+304, NULL);
 _BUILD_00p0(Data,vrr_stack+295, Data->F+10, Data->F+11, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+298, vrr_stack+0, vrr_stack+295, Data->F+9, Data->F+10, NULL);
 _BUILD_00f0(Data,vrr_stack+687, vrr_stack+90, vrr_stack+298, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+697, vrr_stack+117, vrr_stack+687, vrr_stack+136, vrr_stack+90, NULL);
 _BUILD_00h0(Data,vrr_stack+712, vrr_stack+6, vrr_stack+697, vrr_stack+65, vrr_stack+117, NULL);
 _BUILD_00i0(Data,vrr_stack+733, vrr_stack+163, vrr_stack+712, vrr_stack+75, vrr_stack+6, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+687, vrr_stack+21, vrr_stack+733, vrr_stack+376, vrr_stack+163, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+723, vrr_stack+551, vrr_stack+687, vrr_stack+397, vrr_stack+21, NULL);
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+768, vrr_stack+587, vrr_stack+723, vrr_stack+425, vrr_stack+551, NULL);
 am[0] = 0;  am[1] = 10;
 vrr_build_xxxx(am,Data,vrr_stack+687, vrr_stack+240, vrr_stack+768, vrr_stack+461, vrr_stack+587, NULL);
 am[0] = 0;  am[1] = 10;
 vrr_build_xxxx(am,Data,vrr_stack+551, vrr_stack+632, vrr_stack+240, vrr_stack+506, vrr_stack+461, NULL);
 _BUILD_p0h0(Data,vrr_stack+753, vrr_stack+142, vrr_stack+96, NULL, NULL, vrr_stack+50);
   tmp = vrr_stack + 753;
   target_ptr = Libint->vrr_classes[1][5];
   for(i=0;i<63;i++)
     target_ptr[i] += tmp[i];
 _BUILD_p0i0(Data,vrr_stack+816, vrr_stack+212, vrr_stack+184, NULL, NULL, vrr_stack+96);
   tmp = vrr_stack + 816;
   target_ptr = Libint->vrr_classes[1][6];
   for(i=0;i<84;i++)
     target_ptr[i] += tmp[i];
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+0, vrr_stack+340, vrr_stack+304, NULL, NULL, vrr_stack+184);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[1][7];
   for(i=0;i<108;i++)
     target_ptr[i] += tmp[i];
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+900, vrr_stack+506, vrr_stack+461, NULL, NULL, vrr_stack+304);
   tmp = vrr_stack + 900;
   target_ptr = Libint->vrr_classes[1][8];
   for(i=0;i<135;i++)
     target_ptr[i] += tmp[i];
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+1035, vrr_stack+632, vrr_stack+240, NULL, NULL, vrr_stack+461);
   tmp = vrr_stack + 1035;
   target_ptr = Libint->vrr_classes[1][9];
   for(i=0;i<165;i++)
     target_ptr[i] += tmp[i];
 am[0] = 1;  am[1] = 10;
 vrr_build_xxxx(am,Data,vrr_stack+1200, vrr_stack+551, vrr_stack+687, NULL, NULL, vrr_stack+240);
   tmp = vrr_stack + 1200;
   target_ptr = Libint->vrr_classes[1][10];
   for(i=0;i<198;i++)
     target_ptr[i] += tmp[i];

}

