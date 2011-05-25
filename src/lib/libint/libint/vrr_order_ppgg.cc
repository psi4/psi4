#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (pp|gg) integrals */

void vrr_order_ppgg(Libint_t * Libint, prim_data *Data)
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
   tmp = vrr_stack + 196;
   target_ptr = Libint->vrr_classes[1][4];
   for(i=0;i<45;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00h0(Data,vrr_stack+241, vrr_stack+80, vrr_stack+136, vrr_stack+40, vrr_stack+21, NULL);
 _BUILD_00h0(Data,vrr_stack+262, vrr_stack+111, vrr_stack+80, vrr_stack+101, vrr_stack+40, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+101, vrr_stack+95, vrr_stack+6, vrr_stack+12, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+31, vrr_stack+126, vrr_stack+101, vrr_stack+15, vrr_stack+95, NULL);
 _BUILD_00h0(Data,vrr_stack+283, vrr_stack+136, vrr_stack+31, vrr_stack+21, vrr_stack+126, NULL);
 _BUILD_p0h0(Data,vrr_stack+304, vrr_stack+241, vrr_stack+283, NULL, NULL, vrr_stack+136);
 _BUILD_p0h0(Data,vrr_stack+367, vrr_stack+262, vrr_stack+241, NULL, NULL, vrr_stack+80);
   tmp = vrr_stack + 367;
   target_ptr = Libint->vrr_classes[1][5];
   for(i=0;i<63;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00i0(Data,vrr_stack+430, vrr_stack+241, vrr_stack+283, vrr_stack+80, vrr_stack+136, NULL);
 _BUILD_00i0(Data,vrr_stack+458, vrr_stack+262, vrr_stack+241, vrr_stack+111, vrr_stack+80, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+0, vrr_stack+12, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+6, vrr_stack+15, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+486, vrr_stack+101, vrr_stack+21, vrr_stack+95, vrr_stack+6, NULL);
 _BUILD_00h0(Data,vrr_stack+501, vrr_stack+31, vrr_stack+486, vrr_stack+126, vrr_stack+101, NULL);
 _BUILD_00i0(Data,vrr_stack+522, vrr_stack+283, vrr_stack+501, vrr_stack+136, vrr_stack+31, NULL);
 _BUILD_p0i0(Data,vrr_stack+550, vrr_stack+430, vrr_stack+522, NULL, NULL, vrr_stack+283);
 _BUILD_p0i0(Data,vrr_stack+634, vrr_stack+458, vrr_stack+430, NULL, NULL, vrr_stack+241);
   tmp = vrr_stack + 634;
   target_ptr = Libint->vrr_classes[1][6];
   for(i=0;i<84;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+718, vrr_stack+430, vrr_stack+522, vrr_stack+241, vrr_stack+283, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+754, vrr_stack+458, vrr_stack+430, vrr_stack+262, vrr_stack+241, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+95, vrr_stack+12, vrr_stack+3, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+126, vrr_stack+15, vrr_stack+95, vrr_stack+0, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+136, vrr_stack+21, vrr_stack+126, vrr_stack+6, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+790, vrr_stack+486, vrr_stack+136, vrr_stack+101, vrr_stack+21, NULL);
 _BUILD_00i0(Data,vrr_stack+811, vrr_stack+501, vrr_stack+790, vrr_stack+31, vrr_stack+486, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+839, vrr_stack+522, vrr_stack+811, vrr_stack+283, vrr_stack+501, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+875, vrr_stack+718, vrr_stack+839, NULL, NULL, vrr_stack+522);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+983, vrr_stack+754, vrr_stack+718, NULL, NULL, vrr_stack+430);
   tmp = vrr_stack + 983;
   target_ptr = Libint->vrr_classes[1][7];
   for(i=0;i<108;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1091, vrr_stack+718, vrr_stack+839, vrr_stack+430, vrr_stack+522, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1136, vrr_stack+754, vrr_stack+718, vrr_stack+458, vrr_stack+430, NULL);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+9, Data->F+10, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, Data->F+8, Data->F+9, NULL);
 _BUILD_00f0(Data,vrr_stack+101, vrr_stack+95, vrr_stack+6, vrr_stack+12, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+0, vrr_stack+126, vrr_stack+101, vrr_stack+15, vrr_stack+95, NULL);
 _BUILD_00h0(Data,vrr_stack+283, vrr_stack+136, vrr_stack+0, vrr_stack+21, vrr_stack+126, NULL);
 _BUILD_00i0(Data,vrr_stack+0, vrr_stack+790, vrr_stack+283, vrr_stack+486, vrr_stack+136, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1181, vrr_stack+811, vrr_stack+0, vrr_stack+501, vrr_stack+790, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+0, vrr_stack+839, vrr_stack+1181, vrr_stack+522, vrr_stack+811, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1181, vrr_stack+1091, vrr_stack+0, NULL, NULL, vrr_stack+839);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1316, vrr_stack+1136, vrr_stack+1091, NULL, NULL, vrr_stack+718);
   tmp = vrr_stack + 1316;
   target_ptr = Libint->vrr_classes[1][8];
   for(i=0;i<135;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0g0(Data,vrr_stack+1451, vrr_stack+196, vrr_stack+151, vrr_stack+111, vrr_stack+80, vrr_stack+50);
   tmp = vrr_stack + 1451;
   target_ptr = Libint->vrr_classes[2][4];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0h0(Data,vrr_stack+0, vrr_stack+367, vrr_stack+304, vrr_stack+262, vrr_stack+241, vrr_stack+151);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[2][5];
   for(i=0;i<126;i++)
     target_ptr[i] += tmp[i];
 _BUILD_d0i0(Data,vrr_stack+126, vrr_stack+634, vrr_stack+550, vrr_stack+458, vrr_stack+430, vrr_stack+304);
   tmp = vrr_stack + 126;
   target_ptr = Libint->vrr_classes[2][6];
   for(i=0;i<168;i++)
     target_ptr[i] += tmp[i];
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+294, vrr_stack+983, vrr_stack+875, vrr_stack+754, vrr_stack+718, vrr_stack+550);
   tmp = vrr_stack + 294;
   target_ptr = Libint->vrr_classes[2][7];
   for(i=0;i<216;i++)
     target_ptr[i] += tmp[i];
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+510, vrr_stack+1316, vrr_stack+1181, vrr_stack+1136, vrr_stack+1091, vrr_stack+875);
   tmp = vrr_stack + 510;
   target_ptr = Libint->vrr_classes[2][8];
   for(i=0;i<270;i++)
     target_ptr[i] += tmp[i];

}

