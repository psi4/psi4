#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (dd|gg) integrals */

void vrr_order_ddgg(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, NULL, NULL, Data->F+4);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+3, vrr_stack+0, Data->F+3, Data->F+4, NULL);
 _BUILD_00p0(Data,vrr_stack+21, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+24, vrr_stack+21, vrr_stack+3, Data->F+2, Data->F+3, NULL);
 _BUILD_00p0(Data,vrr_stack+30, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+33, vrr_stack+0, vrr_stack+30, Data->F+4, Data->F+5, NULL);
 _BUILD_p0d0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+33, NULL, NULL, vrr_stack+0);
 _BUILD_p0d0(Data,vrr_stack+57, vrr_stack+24, vrr_stack+15, NULL, NULL, vrr_stack+3);
 _BUILD_d0d0(Data,vrr_stack+75, vrr_stack+57, vrr_stack+39, vrr_stack+24, vrr_stack+15, vrr_stack+6);
 _BUILD_00f0(Data,vrr_stack+111, vrr_stack+15, vrr_stack+33, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00f0(Data,vrr_stack+121, vrr_stack+24, vrr_stack+15, vrr_stack+21, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+131, vrr_stack+121, vrr_stack+111, NULL, NULL, vrr_stack+15);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+21, Data->F+1, Data->F+2, NULL);
 _BUILD_00f0(Data,vrr_stack+161, vrr_stack+6, vrr_stack+24, vrr_stack+3, vrr_stack+21, NULL);
 _BUILD_p0f0(Data,vrr_stack+171, vrr_stack+161, vrr_stack+121, NULL, NULL, vrr_stack+24);
 _BUILD_00p0(Data,vrr_stack+21, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+201, vrr_stack+30, vrr_stack+21, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+207, vrr_stack+33, vrr_stack+201, vrr_stack+0, vrr_stack+30, NULL);
 _BUILD_p0f0(Data,vrr_stack+217, vrr_stack+111, vrr_stack+207, NULL, NULL, vrr_stack+33);
 _BUILD_d0f0(Data,vrr_stack+247, vrr_stack+131, vrr_stack+217, vrr_stack+121, vrr_stack+111, vrr_stack+39);
 _BUILD_d0f0(Data,vrr_stack+307, vrr_stack+171, vrr_stack+131, vrr_stack+161, vrr_stack+121, vrr_stack+57);
 _BUILD_f0f0(Data,vrr_stack+367, vrr_stack+307, vrr_stack+247, vrr_stack+171, vrr_stack+131, vrr_stack+75);
 _BUILD_00g0(Data,vrr_stack+39, vrr_stack+121, vrr_stack+111, vrr_stack+24, vrr_stack+15, NULL);
 _BUILD_00g0(Data,vrr_stack+54, vrr_stack+161, vrr_stack+121, vrr_stack+6, vrr_stack+24, NULL);
 _BUILD_00g0(Data,vrr_stack+69, vrr_stack+111, vrr_stack+207, vrr_stack+15, vrr_stack+33, NULL);
 _BUILD_p0g0(Data,vrr_stack+467, vrr_stack+39, vrr_stack+69, NULL, NULL, vrr_stack+111);
 _BUILD_p0g0(Data,vrr_stack+512, vrr_stack+54, vrr_stack+39, NULL, NULL, vrr_stack+121);
 _BUILD_d0g0(Data,vrr_stack+557, vrr_stack+512, vrr_stack+467, vrr_stack+54, vrr_stack+39, vrr_stack+131);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+24, vrr_stack+0, vrr_stack+3, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+131, vrr_stack+24, vrr_stack+6, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+141, vrr_stack+131, vrr_stack+161, vrr_stack+24, vrr_stack+6, NULL);
 _BUILD_p0g0(Data,vrr_stack+647, vrr_stack+141, vrr_stack+54, NULL, NULL, vrr_stack+161);
 _BUILD_d0g0(Data,vrr_stack+692, vrr_stack+647, vrr_stack+512, vrr_stack+141, vrr_stack+54, vrr_stack+171);
   tmp = vrr_stack + 692;
   target_ptr = Libint->vrr_classes[2][4];
   for(i=0;i<90;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+171, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+24, vrr_stack+21, vrr_stack+171, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+174, vrr_stack+201, vrr_stack+24, vrr_stack+30, vrr_stack+21, NULL);
 _BUILD_00g0(Data,vrr_stack+184, vrr_stack+207, vrr_stack+174, vrr_stack+33, vrr_stack+201, NULL);
 _BUILD_p0g0(Data,vrr_stack+782, vrr_stack+69, vrr_stack+184, NULL, NULL, vrr_stack+207);
 _BUILD_d0g0(Data,vrr_stack+827, vrr_stack+467, vrr_stack+782, vrr_stack+39, vrr_stack+69, vrr_stack+217);
 _BUILD_f0g0(Data,vrr_stack+917, vrr_stack+557, vrr_stack+827, vrr_stack+512, vrr_stack+467, vrr_stack+247);
 _BUILD_f0g0(Data,vrr_stack+1067, vrr_stack+692, vrr_stack+557, vrr_stack+647, vrr_stack+512, vrr_stack+307);
   tmp = vrr_stack + 1067;
   target_ptr = Libint->vrr_classes[3][4];
   for(i=0;i<150;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00h0(Data,vrr_stack+0, vrr_stack+39, vrr_stack+69, vrr_stack+121, vrr_stack+111, NULL);
 _BUILD_00h0(Data,vrr_stack+647, vrr_stack+54, vrr_stack+39, vrr_stack+161, vrr_stack+121, NULL);
 _BUILD_00h0(Data,vrr_stack+668, vrr_stack+69, vrr_stack+184, vrr_stack+111, vrr_stack+207, NULL);
 _BUILD_p0h0(Data,vrr_stack+217, vrr_stack+0, vrr_stack+668, NULL, NULL, vrr_stack+69);
 _BUILD_p0h0(Data,vrr_stack+280, vrr_stack+647, vrr_stack+0, NULL, NULL, vrr_stack+39);
 _BUILD_d0h0(Data,vrr_stack+1217, vrr_stack+280, vrr_stack+217, vrr_stack+647, vrr_stack+0, vrr_stack+467);
 _BUILD_00h0(Data,vrr_stack+467, vrr_stack+141, vrr_stack+54, vrr_stack+131, vrr_stack+161, NULL);
 _BUILD_p0h0(Data,vrr_stack+1343, vrr_stack+467, vrr_stack+647, NULL, NULL, vrr_stack+54);
 _BUILD_d0h0(Data,vrr_stack+1406, vrr_stack+1343, vrr_stack+280, vrr_stack+467, vrr_stack+647, vrr_stack+512);
   tmp = vrr_stack + 1406;
   target_ptr = Libint->vrr_classes[2][5];
   for(i=0;i<126;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+689, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+488, vrr_stack+171, vrr_stack+689, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+494, vrr_stack+24, vrr_stack+488, vrr_stack+21, vrr_stack+171, NULL);
 _BUILD_00g0(Data,vrr_stack+156, vrr_stack+174, vrr_stack+494, vrr_stack+201, vrr_stack+24, NULL);
 _BUILD_00h0(Data,vrr_stack+504, vrr_stack+184, vrr_stack+156, vrr_stack+207, vrr_stack+174, NULL);
 _BUILD_p0h0(Data,vrr_stack+1532, vrr_stack+668, vrr_stack+504, NULL, NULL, vrr_stack+184);
 _BUILD_d0h0(Data,vrr_stack+1595, vrr_stack+217, vrr_stack+1532, vrr_stack+0, vrr_stack+668, vrr_stack+782);
 _BUILD_f0h0(Data,vrr_stack+1721, vrr_stack+1217, vrr_stack+1595, vrr_stack+280, vrr_stack+217, vrr_stack+827);
 _BUILD_f0h0(Data,vrr_stack+1931, vrr_stack+1406, vrr_stack+1217, vrr_stack+1343, vrr_stack+280, vrr_stack+557);
   tmp = vrr_stack + 1931;
   target_ptr = Libint->vrr_classes[3][5];
   for(i=0;i<210;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00i0(Data,vrr_stack+1343, vrr_stack+0, vrr_stack+668, vrr_stack+39, vrr_stack+69, NULL);
 _BUILD_00i0(Data,vrr_stack+1371, vrr_stack+647, vrr_stack+0, vrr_stack+54, vrr_stack+39, NULL);
 _BUILD_00i0(Data,vrr_stack+782, vrr_stack+668, vrr_stack+504, vrr_stack+69, vrr_stack+184, NULL);
 _BUILD_p0i0(Data,vrr_stack+810, vrr_stack+1343, vrr_stack+782, NULL, NULL, vrr_stack+668);
 _BUILD_p0i0(Data,vrr_stack+2141, vrr_stack+1371, vrr_stack+1343, NULL, NULL, vrr_stack+0);
 _BUILD_d0i0(Data,vrr_stack+2225, vrr_stack+2141, vrr_stack+810, vrr_stack+1371, vrr_stack+1343, vrr_stack+217);
 _BUILD_00i0(Data,vrr_stack+69, vrr_stack+467, vrr_stack+647, vrr_stack+141, vrr_stack+54, NULL);
 _BUILD_p0i0(Data,vrr_stack+2393, vrr_stack+69, vrr_stack+1371, NULL, NULL, vrr_stack+647);
 _BUILD_d0i0(Data,vrr_stack+2477, vrr_stack+2393, vrr_stack+2141, vrr_stack+69, vrr_stack+1371, vrr_stack+280);
   tmp = vrr_stack + 2477;
   target_ptr = Libint->vrr_classes[2][6];
   for(i=0;i<168;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+21, Data->F+9, Data->F+10, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+97, vrr_stack+689, vrr_stack+21, Data->F+8, Data->F+9, NULL);
 _BUILD_00f0(Data,vrr_stack+103, vrr_stack+488, vrr_stack+97, vrr_stack+171, vrr_stack+689, NULL);
 _BUILD_00g0(Data,vrr_stack+113, vrr_stack+494, vrr_stack+103, vrr_stack+24, vrr_stack+488, NULL);
 _BUILD_00h0(Data,vrr_stack+24, vrr_stack+156, vrr_stack+113, vrr_stack+174, vrr_stack+494, NULL);
 _BUILD_00i0(Data,vrr_stack+128, vrr_stack+504, vrr_stack+24, vrr_stack+184, vrr_stack+156, NULL);
 _BUILD_p0i0(Data,vrr_stack+171, vrr_stack+782, vrr_stack+128, NULL, NULL, vrr_stack+504);
 _BUILD_d0i0(Data,vrr_stack+2645, vrr_stack+810, vrr_stack+171, vrr_stack+1343, vrr_stack+782, vrr_stack+1532);
 _BUILD_f0i0(Data,vrr_stack+2813, vrr_stack+2225, vrr_stack+2645, vrr_stack+2141, vrr_stack+810, vrr_stack+1595);
 _BUILD_f0i0(Data,vrr_stack+3093, vrr_stack+2477, vrr_stack+2225, vrr_stack+2393, vrr_stack+2141, vrr_stack+1217);
   tmp = vrr_stack + 3093;
   target_ptr = Libint->vrr_classes[3][6];
   for(i=0;i<280;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+2393, vrr_stack+1343, vrr_stack+782, vrr_stack+0, vrr_stack+668, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+2429, vrr_stack+1371, vrr_stack+1343, vrr_stack+647, vrr_stack+0, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1532, vrr_stack+782, vrr_stack+128, vrr_stack+668, vrr_stack+504, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1568, vrr_stack+2393, vrr_stack+1532, NULL, NULL, vrr_stack+782);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+255, vrr_stack+2429, vrr_stack+2393, NULL, NULL, vrr_stack+1343);
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3373, vrr_stack+255, vrr_stack+1568, vrr_stack+2429, vrr_stack+2393, vrr_stack+810);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+810, vrr_stack+69, vrr_stack+1371, vrr_stack+467, vrr_stack+647, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3589, vrr_stack+810, vrr_stack+2429, NULL, NULL, vrr_stack+1371);
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3697, vrr_stack+3589, vrr_stack+255, vrr_stack+810, vrr_stack+2429, vrr_stack+2141);
   tmp = vrr_stack + 3697;
   target_ptr = Libint->vrr_classes[2][7];
   for(i=0;i<216;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+2141, Data->F+10, Data->F+11, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+2144, vrr_stack+21, vrr_stack+2141, Data->F+9, Data->F+10, NULL);
 _BUILD_00f0(Data,vrr_stack+2150, vrr_stack+97, vrr_stack+2144, vrr_stack+689, vrr_stack+21, NULL);
 _BUILD_00g0(Data,vrr_stack+2160, vrr_stack+103, vrr_stack+2150, vrr_stack+488, vrr_stack+97, NULL);
 _BUILD_00h0(Data,vrr_stack+0, vrr_stack+113, vrr_stack+2160, vrr_stack+494, vrr_stack+103, NULL);
 _BUILD_00i0(Data,vrr_stack+2175, vrr_stack+24, vrr_stack+0, vrr_stack+156, vrr_stack+113, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+647, vrr_stack+128, vrr_stack+2175, vrr_stack+504, vrr_stack+24, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3913, vrr_stack+1532, vrr_stack+647, NULL, NULL, vrr_stack+128);
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+4021, vrr_stack+1568, vrr_stack+3913, vrr_stack+2393, vrr_stack+1532, vrr_stack+171);
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+4237, vrr_stack+3373, vrr_stack+4021, vrr_stack+255, vrr_stack+1568, vrr_stack+2645);
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+4597, vrr_stack+3697, vrr_stack+3373, vrr_stack+3589, vrr_stack+255, vrr_stack+2225);
   tmp = vrr_stack + 4597;
   target_ptr = Libint->vrr_classes[3][7];
   for(i=0;i<360;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1676, vrr_stack+2393, vrr_stack+1532, vrr_stack+1343, vrr_stack+782, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+3589, vrr_stack+2429, vrr_stack+2393, vrr_stack+1371, vrr_stack+1343, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+3634, vrr_stack+1532, vrr_stack+647, vrr_stack+782, vrr_stack+128, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+2645, vrr_stack+1676, vrr_stack+3634, NULL, NULL, vrr_stack+1532);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+4957, vrr_stack+3589, vrr_stack+1676, NULL, NULL, vrr_stack+2393);
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+5092, vrr_stack+4957, vrr_stack+2645, vrr_stack+3589, vrr_stack+1676, vrr_stack+1568);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1532, vrr_stack+810, vrr_stack+2429, vrr_stack+69, vrr_stack+1371, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+782, vrr_stack+1532, vrr_stack+3589, NULL, NULL, vrr_stack+2429);
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+5362, vrr_stack+782, vrr_stack+4957, vrr_stack+1532, vrr_stack+3589, vrr_stack+255);
   tmp = vrr_stack + 5362;
   target_ptr = Libint->vrr_classes[2][8];
   for(i=0;i<270;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+3589, Data->F+11, Data->F+12, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+3592, vrr_stack+2141, vrr_stack+3589, Data->F+10, Data->F+11, NULL);
 _BUILD_00f0(Data,vrr_stack+3598, vrr_stack+2144, vrr_stack+3592, vrr_stack+21, vrr_stack+2141, NULL);
 _BUILD_00g0(Data,vrr_stack+3608, vrr_stack+2150, vrr_stack+3598, vrr_stack+97, vrr_stack+2144, NULL);
 _BUILD_00h0(Data,vrr_stack+1532, vrr_stack+2160, vrr_stack+3608, vrr_stack+103, vrr_stack+2150, NULL);
 _BUILD_00i0(Data,vrr_stack+3589, vrr_stack+0, vrr_stack+1532, vrr_stack+113, vrr_stack+2160, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1532, vrr_stack+2175, vrr_stack+3589, vrr_stack+24, vrr_stack+0, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+3589, vrr_stack+647, vrr_stack+1532, vrr_stack+128, vrr_stack+2175, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1532, vrr_stack+3634, vrr_stack+3589, NULL, NULL, vrr_stack+647);
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+0, vrr_stack+2645, vrr_stack+1532, vrr_stack+1676, vrr_stack+3634, vrr_stack+3913);
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+5632, vrr_stack+5092, vrr_stack+0, vrr_stack+4957, vrr_stack+2645, vrr_stack+4021);
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+6082, vrr_stack+5362, vrr_stack+5092, vrr_stack+782, vrr_stack+4957, vrr_stack+3373);
   tmp = vrr_stack + 6082;
   target_ptr = Libint->vrr_classes[3][8];
   for(i=0;i<450;i++)
     target_ptr[i] += tmp[i];
 _BUILD_g0g0(Data,vrr_stack+0, vrr_stack+1067, vrr_stack+917, vrr_stack+692, vrr_stack+557, vrr_stack+367);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[4][4];
   for(i=0;i<225;i++)
     target_ptr[i] += tmp[i];
 _BUILD_g0h0(Data,vrr_stack+225, vrr_stack+1931, vrr_stack+1721, vrr_stack+1406, vrr_stack+1217, vrr_stack+917);
   tmp = vrr_stack + 225;
   target_ptr = Libint->vrr_classes[4][5];
   for(i=0;i<315;i++)
     target_ptr[i] += tmp[i];
 _BUILD_g0i0(Data,vrr_stack+540, vrr_stack+3093, vrr_stack+2813, vrr_stack+2477, vrr_stack+2225, vrr_stack+1721);
   tmp = vrr_stack + 540;
   target_ptr = Libint->vrr_classes[4][6];
   for(i=0;i<420;i++)
     target_ptr[i] += tmp[i];
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+960, vrr_stack+4597, vrr_stack+4237, vrr_stack+3697, vrr_stack+3373, vrr_stack+2813);
   tmp = vrr_stack + 960;
   target_ptr = Libint->vrr_classes[4][7];
   for(i=0;i<540;i++)
     target_ptr[i] += tmp[i];
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1500, vrr_stack+6082, vrr_stack+5632, vrr_stack+5362, vrr_stack+5092, vrr_stack+4237);
   tmp = vrr_stack + 1500;
   target_ptr = Libint->vrr_classes[4][8];
   for(i=0;i<675;i++)
     target_ptr[i] += tmp[i];

}

