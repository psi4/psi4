#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (dp|gg) integrals */

void vrr_order_dpgg(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+0, vrr_stack+3, Data->F+3, Data->F+4, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+0, Data->F+2, Data->F+3, NULL);
 _BUILD_p0d0(Data,vrr_stack+21, vrr_stack+15, vrr_stack+6, NULL, NULL, vrr_stack+0);
 _BUILD_00f0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+6, vrr_stack+12, vrr_stack+0, NULL);
 _BUILD_00p0(Data,vrr_stack+49, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+52, vrr_stack+49, vrr_stack+12, Data->F+1, Data->F+2, NULL);
 _BUILD_00f0(Data,vrr_stack+58, vrr_stack+52, vrr_stack+15, vrr_stack+49, vrr_stack+12, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+68, vrr_stack+3, vrr_stack+12, Data->F+4, Data->F+5, NULL);
 _BUILD_00f0(Data,vrr_stack+74, vrr_stack+6, vrr_stack+68, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+84, vrr_stack+39, vrr_stack+74, NULL, NULL, vrr_stack+6);
 _BUILD_p0f0(Data,vrr_stack+114, vrr_stack+58, vrr_stack+39, NULL, NULL, vrr_stack+15);
 _BUILD_d0f0(Data,vrr_stack+144, vrr_stack+114, vrr_stack+84, vrr_stack+58, vrr_stack+39, vrr_stack+21);
 _BUILD_00g0(Data,vrr_stack+21, vrr_stack+39, vrr_stack+74, vrr_stack+15, vrr_stack+6, NULL);
 _BUILD_00g0(Data,vrr_stack+204, vrr_stack+58, vrr_stack+39, vrr_stack+52, vrr_stack+15, NULL);
 _BUILD_p0g0(Data,vrr_stack+219, vrr_stack+204, vrr_stack+21, NULL, NULL, vrr_stack+39);
 _BUILD_00p0(Data,vrr_stack+36, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+36, vrr_stack+49, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+264, vrr_stack+15, vrr_stack+52, vrr_stack+36, vrr_stack+49, NULL);
 _BUILD_00g0(Data,vrr_stack+274, vrr_stack+264, vrr_stack+58, vrr_stack+15, vrr_stack+52, NULL);
 _BUILD_p0g0(Data,vrr_stack+289, vrr_stack+274, vrr_stack+204, NULL, NULL, vrr_stack+58);
 _BUILD_00p0(Data,vrr_stack+36, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+36, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+334, vrr_stack+68, vrr_stack+15, vrr_stack+3, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+344, vrr_stack+74, vrr_stack+334, vrr_stack+6, vrr_stack+68, NULL);
 _BUILD_p0g0(Data,vrr_stack+359, vrr_stack+21, vrr_stack+344, NULL, NULL, vrr_stack+74);
 _BUILD_d0g0(Data,vrr_stack+404, vrr_stack+219, vrr_stack+359, vrr_stack+204, vrr_stack+21, vrr_stack+84);
 _BUILD_d0g0(Data,vrr_stack+494, vrr_stack+289, vrr_stack+219, vrr_stack+274, vrr_stack+204, vrr_stack+114);
   tmp = vrr_stack + 494;
   target_ptr = Libint->vrr_classes[2][4];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00h0(Data,vrr_stack+84, vrr_stack+21, vrr_stack+344, vrr_stack+39, vrr_stack+74, NULL);
 _BUILD_00h0(Data,vrr_stack+105, vrr_stack+204, vrr_stack+21, vrr_stack+58, vrr_stack+39, NULL);
 _BUILD_p0h0(Data,vrr_stack+584, vrr_stack+105, vrr_stack+84, NULL, NULL, vrr_stack+21);
 _BUILD_00h0(Data,vrr_stack+647, vrr_stack+274, vrr_stack+204, vrr_stack+264, vrr_stack+58, NULL);
 _BUILD_p0h0(Data,vrr_stack+668, vrr_stack+647, vrr_stack+105, NULL, NULL, vrr_stack+204);
 _BUILD_00p0(Data,vrr_stack+264, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+267, vrr_stack+36, vrr_stack+264, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+267, vrr_stack+12, vrr_stack+36, NULL);
 _BUILD_00g0(Data,vrr_stack+0, vrr_stack+334, vrr_stack+39, vrr_stack+68, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+49, vrr_stack+344, vrr_stack+0, vrr_stack+74, vrr_stack+334, NULL);
 _BUILD_p0h0(Data,vrr_stack+731, vrr_stack+84, vrr_stack+49, NULL, NULL, vrr_stack+344);
 _BUILD_d0h0(Data,vrr_stack+794, vrr_stack+584, vrr_stack+731, vrr_stack+105, vrr_stack+84, vrr_stack+359);
 _BUILD_d0h0(Data,vrr_stack+920, vrr_stack+668, vrr_stack+584, vrr_stack+647, vrr_stack+105, vrr_stack+219);
   tmp = vrr_stack + 920;
   target_ptr = Libint->vrr_classes[2][5];
   for(i=0;i<126;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00i0(Data,vrr_stack+359, vrr_stack+84, vrr_stack+49, vrr_stack+21, vrr_stack+344, NULL);
 _BUILD_00i0(Data,vrr_stack+1046, vrr_stack+105, vrr_stack+84, vrr_stack+204, vrr_stack+21, NULL);
 _BUILD_p0i0(Data,vrr_stack+1074, vrr_stack+1046, vrr_stack+359, NULL, NULL, vrr_stack+84);
 _BUILD_00i0(Data,vrr_stack+1158, vrr_stack+647, vrr_stack+105, vrr_stack+274, vrr_stack+204, NULL);
 _BUILD_p0i0(Data,vrr_stack+1186, vrr_stack+1158, vrr_stack+1046, NULL, NULL, vrr_stack+105);
 _BUILD_00p0(Data,vrr_stack+204, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+207, vrr_stack+264, vrr_stack+204, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+21, vrr_stack+267, vrr_stack+207, vrr_stack+36, vrr_stack+264, NULL);
 _BUILD_00g0(Data,vrr_stack+387, vrr_stack+39, vrr_stack+21, vrr_stack+15, vrr_stack+267, NULL);
 _BUILD_00h0(Data,vrr_stack+1270, vrr_stack+0, vrr_stack+387, vrr_stack+334, vrr_stack+39, NULL);
 _BUILD_00i0(Data,vrr_stack+1291, vrr_stack+49, vrr_stack+1270, vrr_stack+344, vrr_stack+0, NULL);
 _BUILD_p0i0(Data,vrr_stack+1319, vrr_stack+359, vrr_stack+1291, NULL, NULL, vrr_stack+49);
 _BUILD_d0i0(Data,vrr_stack+1403, vrr_stack+1074, vrr_stack+1319, vrr_stack+1046, vrr_stack+359, vrr_stack+731);
 _BUILD_d0i0(Data,vrr_stack+1571, vrr_stack+1186, vrr_stack+1074, vrr_stack+1158, vrr_stack+1046, vrr_stack+584);
   tmp = vrr_stack + 1571;
   target_ptr = Libint->vrr_classes[2][6];
   for(i=0;i<168;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+731, vrr_stack+359, vrr_stack+1291, vrr_stack+84, vrr_stack+49, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1739, vrr_stack+1046, vrr_stack+359, vrr_stack+105, vrr_stack+84, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1775, vrr_stack+1739, vrr_stack+731, NULL, NULL, vrr_stack+359);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1883, vrr_stack+1158, vrr_stack+1046, vrr_stack+647, vrr_stack+105, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1919, vrr_stack+1883, vrr_stack+1739, NULL, NULL, vrr_stack+1046);
 _BUILD_00p0(Data,vrr_stack+647, Data->F+9, Data->F+10, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+204, vrr_stack+647, Data->F+8, Data->F+9, NULL);
 _BUILD_00f0(Data,vrr_stack+650, vrr_stack+207, vrr_stack+15, vrr_stack+264, vrr_stack+204, NULL);
 _BUILD_00g0(Data,vrr_stack+767, vrr_stack+21, vrr_stack+650, vrr_stack+267, vrr_stack+207, NULL);
 _BUILD_00h0(Data,vrr_stack+264, vrr_stack+387, vrr_stack+767, vrr_stack+39, vrr_stack+21, NULL);
 _BUILD_00i0(Data,vrr_stack+70, vrr_stack+1270, vrr_stack+264, vrr_stack+0, vrr_stack+387, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+98, vrr_stack+1291, vrr_stack+70, vrr_stack+49, vrr_stack+1270, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+2027, vrr_stack+731, vrr_stack+98, NULL, NULL, vrr_stack+1291);
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+2135, vrr_stack+1775, vrr_stack+2027, vrr_stack+1739, vrr_stack+731, vrr_stack+1319);
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+2351, vrr_stack+1919, vrr_stack+1775, vrr_stack+1883, vrr_stack+1739, vrr_stack+1074);
   tmp = vrr_stack + 2351;
   target_ptr = Libint->vrr_classes[2][7];
   for(i=0;i<216;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1319, vrr_stack+731, vrr_stack+98, vrr_stack+359, vrr_stack+1291, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+2567, vrr_stack+1739, vrr_stack+731, vrr_stack+1046, vrr_stack+359, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+2612, vrr_stack+2567, vrr_stack+1319, NULL, NULL, vrr_stack+731);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+334, vrr_stack+1883, vrr_stack+1739, vrr_stack+1158, vrr_stack+1046, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+2747, vrr_stack+334, vrr_stack+2567, NULL, NULL, vrr_stack+1739);
 _BUILD_00p0(Data,vrr_stack+1739, Data->F+10, Data->F+11, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+213, vrr_stack+647, vrr_stack+1739, Data->F+9, Data->F+10, NULL);
 _BUILD_00f0(Data,vrr_stack+134, vrr_stack+15, vrr_stack+213, vrr_stack+204, vrr_stack+647, NULL);
 _BUILD_00g0(Data,vrr_stack+0, vrr_stack+650, vrr_stack+134, vrr_stack+207, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+1739, vrr_stack+767, vrr_stack+0, vrr_stack+21, vrr_stack+650, NULL);
 _BUILD_00i0(Data,vrr_stack+1046, vrr_stack+264, vrr_stack+1739, vrr_stack+387, vrr_stack+767, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1739, vrr_stack+70, vrr_stack+1046, vrr_stack+1270, vrr_stack+264, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+0, vrr_stack+98, vrr_stack+1739, vrr_stack+1291, vrr_stack+70, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+2882, vrr_stack+1319, vrr_stack+0, NULL, NULL, vrr_stack+98);
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+3017, vrr_stack+2612, vrr_stack+2882, vrr_stack+2567, vrr_stack+1319, vrr_stack+2027);
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+3287, vrr_stack+2747, vrr_stack+2612, vrr_stack+334, vrr_stack+2567, vrr_stack+1775);
   tmp = vrr_stack + 3287;
   target_ptr = Libint->vrr_classes[2][8];
   for(i=0;i<270;i++)
     target_ptr[i] += tmp[i];
 _BUILD_f0g0(Data,vrr_stack+3557, vrr_stack+494, vrr_stack+404, vrr_stack+289, vrr_stack+219, vrr_stack+144);
   tmp = vrr_stack + 3557;
   target_ptr = Libint->vrr_classes[3][4];
   for(i=0;i<150;i++)
     target_ptr[i] += tmp[i];
 _BUILD_f0h0(Data,vrr_stack+0, vrr_stack+920, vrr_stack+794, vrr_stack+668, vrr_stack+584, vrr_stack+404);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[3][5];
   for(i=0;i<210;i++)
     target_ptr[i] += tmp[i];
 _BUILD_f0i0(Data,vrr_stack+210, vrr_stack+1571, vrr_stack+1403, vrr_stack+1186, vrr_stack+1074, vrr_stack+794);
   tmp = vrr_stack + 210;
   target_ptr = Libint->vrr_classes[3][6];
   for(i=0;i<280;i++)
     target_ptr[i] += tmp[i];
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+490, vrr_stack+2351, vrr_stack+2135, vrr_stack+1919, vrr_stack+1775, vrr_stack+1403);
   tmp = vrr_stack + 490;
   target_ptr = Libint->vrr_classes[3][7];
   for(i=0;i<360;i++)
     target_ptr[i] += tmp[i];
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+850, vrr_stack+3287, vrr_stack+3017, vrr_stack+2747, vrr_stack+2612, vrr_stack+2135);
   tmp = vrr_stack + 850;
   target_ptr = Libint->vrr_classes[3][8];
   for(i=0;i<450;i++)
     target_ptr[i] += tmp[i];

}

