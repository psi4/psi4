#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (dd|hg) integrals */

void vrr_order_ddhg(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+0, vrr_stack+3, Data->F+4, Data->F+5, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+12, vrr_stack+0, Data->F+3, Data->F+4, NULL);
 _BUILD_p0d0(Data,vrr_stack+21, vrr_stack+15, vrr_stack+6, NULL, NULL, vrr_stack+0);
 _BUILD_00f0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+6, vrr_stack+12, vrr_stack+0, NULL);
 _BUILD_00p0(Data,vrr_stack+49, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+52, vrr_stack+49, vrr_stack+12, Data->F+2, Data->F+3, NULL);
 _BUILD_00f0(Data,vrr_stack+58, vrr_stack+52, vrr_stack+15, vrr_stack+49, vrr_stack+12, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+68, vrr_stack+3, vrr_stack+12, Data->F+5, Data->F+6, NULL);
 _BUILD_00f0(Data,vrr_stack+74, vrr_stack+6, vrr_stack+68, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+84, vrr_stack+39, vrr_stack+74, NULL, NULL, vrr_stack+6);
 _BUILD_p0f0(Data,vrr_stack+114, vrr_stack+58, vrr_stack+39, NULL, NULL, vrr_stack+15);
 _BUILD_d0f0(Data,vrr_stack+144, vrr_stack+114, vrr_stack+84, vrr_stack+58, vrr_stack+39, vrr_stack+21);
 _BUILD_00g0(Data,vrr_stack+21, vrr_stack+39, vrr_stack+74, vrr_stack+15, vrr_stack+6, NULL);
 _BUILD_00g0(Data,vrr_stack+204, vrr_stack+58, vrr_stack+39, vrr_stack+52, vrr_stack+15, NULL);
 _BUILD_p0g0(Data,vrr_stack+219, vrr_stack+204, vrr_stack+21, NULL, NULL, vrr_stack+39);
 _BUILD_00p0(Data,vrr_stack+36, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+36, vrr_stack+49, Data->F+1, Data->F+2, NULL);
 _BUILD_00f0(Data,vrr_stack+264, vrr_stack+15, vrr_stack+52, vrr_stack+36, vrr_stack+49, NULL);
 _BUILD_00g0(Data,vrr_stack+274, vrr_stack+264, vrr_stack+58, vrr_stack+15, vrr_stack+52, NULL);
 _BUILD_p0g0(Data,vrr_stack+289, vrr_stack+274, vrr_stack+204, NULL, NULL, vrr_stack+58);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+49, vrr_stack+12, vrr_stack+0, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+334, vrr_stack+68, vrr_stack+49, vrr_stack+3, vrr_stack+12, NULL);
 _BUILD_00g0(Data,vrr_stack+344, vrr_stack+74, vrr_stack+334, vrr_stack+6, vrr_stack+68, NULL);
 _BUILD_p0g0(Data,vrr_stack+359, vrr_stack+21, vrr_stack+344, NULL, NULL, vrr_stack+74);
 _BUILD_d0g0(Data,vrr_stack+404, vrr_stack+219, vrr_stack+359, vrr_stack+204, vrr_stack+21, vrr_stack+84);
 _BUILD_d0g0(Data,vrr_stack+494, vrr_stack+289, vrr_stack+219, vrr_stack+274, vrr_stack+204, vrr_stack+114);
 _BUILD_f0g0(Data,vrr_stack+584, vrr_stack+494, vrr_stack+404, vrr_stack+289, vrr_stack+219, vrr_stack+144);
 _BUILD_00h0(Data,vrr_stack+84, vrr_stack+204, vrr_stack+21, vrr_stack+58, vrr_stack+39, NULL);
 _BUILD_00h0(Data,vrr_stack+105, vrr_stack+274, vrr_stack+204, vrr_stack+264, vrr_stack+58, NULL);
 _BUILD_00h0(Data,vrr_stack+126, vrr_stack+21, vrr_stack+344, vrr_stack+39, vrr_stack+74, NULL);
 _BUILD_p0h0(Data,vrr_stack+734, vrr_stack+84, vrr_stack+126, NULL, NULL, vrr_stack+21);
 _BUILD_p0h0(Data,vrr_stack+797, vrr_stack+105, vrr_stack+84, NULL, NULL, vrr_stack+204);
 _BUILD_d0h0(Data,vrr_stack+860, vrr_stack+797, vrr_stack+734, vrr_stack+105, vrr_stack+84, vrr_stack+219);
 _BUILD_00p0(Data,vrr_stack+219, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+222, vrr_stack+219, vrr_stack+36, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+39, vrr_stack+222, vrr_stack+15, vrr_stack+219, vrr_stack+36, NULL);
 _BUILD_00g0(Data,vrr_stack+986, vrr_stack+39, vrr_stack+264, vrr_stack+222, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+219, vrr_stack+986, vrr_stack+274, vrr_stack+39, vrr_stack+264, NULL);
 _BUILD_p0h0(Data,vrr_stack+1001, vrr_stack+219, vrr_stack+105, NULL, NULL, vrr_stack+274);
 _BUILD_d0h0(Data,vrr_stack+1064, vrr_stack+1001, vrr_stack+797, vrr_stack+219, vrr_stack+105, vrr_stack+289);
   tmp = vrr_stack + 1064;
   target_ptr = Libint->vrr_classes[2][5];
   for(i=0;i<126;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+289, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+0, vrr_stack+289, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+292, vrr_stack+49, vrr_stack+15, vrr_stack+12, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+302, vrr_stack+334, vrr_stack+292, vrr_stack+68, vrr_stack+49, NULL);
 _BUILD_00h0(Data,vrr_stack+240, vrr_stack+344, vrr_stack+302, vrr_stack+74, vrr_stack+334, NULL);
 _BUILD_p0h0(Data,vrr_stack+1190, vrr_stack+126, vrr_stack+240, NULL, NULL, vrr_stack+344);
 _BUILD_d0h0(Data,vrr_stack+1253, vrr_stack+734, vrr_stack+1190, vrr_stack+84, vrr_stack+126, vrr_stack+359);
 _BUILD_f0h0(Data,vrr_stack+1379, vrr_stack+860, vrr_stack+1253, vrr_stack+797, vrr_stack+734, vrr_stack+404);
 _BUILD_f0h0(Data,vrr_stack+1589, vrr_stack+1064, vrr_stack+860, vrr_stack+1001, vrr_stack+797, vrr_stack+494);
   tmp = vrr_stack + 1589;
   target_ptr = Libint->vrr_classes[3][5];
   for(i=0;i<210;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00i0(Data,vrr_stack+1001, vrr_stack+84, vrr_stack+126, vrr_stack+204, vrr_stack+21, NULL);
 _BUILD_00i0(Data,vrr_stack+1029, vrr_stack+105, vrr_stack+84, vrr_stack+274, vrr_stack+204, NULL);
 _BUILD_00i0(Data,vrr_stack+359, vrr_stack+126, vrr_stack+240, vrr_stack+21, vrr_stack+344, NULL);
 _BUILD_p0i0(Data,vrr_stack+387, vrr_stack+1001, vrr_stack+359, NULL, NULL, vrr_stack+126);
 _BUILD_p0i0(Data,vrr_stack+471, vrr_stack+1029, vrr_stack+1001, NULL, NULL, vrr_stack+84);
 _BUILD_d0i0(Data,vrr_stack+1799, vrr_stack+471, vrr_stack+387, vrr_stack+1029, vrr_stack+1001, vrr_stack+734);
 _BUILD_00i0(Data,vrr_stack+21, vrr_stack+219, vrr_stack+105, vrr_stack+986, vrr_stack+274, NULL);
 _BUILD_p0i0(Data,vrr_stack+1967, vrr_stack+21, vrr_stack+1029, NULL, NULL, vrr_stack+105);
 _BUILD_d0i0(Data,vrr_stack+2051, vrr_stack+1967, vrr_stack+471, vrr_stack+21, vrr_stack+1029, vrr_stack+797);
   tmp = vrr_stack + 2051;
   target_ptr = Libint->vrr_classes[2][6];
   for(i=0;i<168;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+986, Data->F+9, Data->F+10, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+989, vrr_stack+289, vrr_stack+986, Data->F+8, Data->F+9, NULL);
 _BUILD_00f0(Data,vrr_stack+734, vrr_stack+15, vrr_stack+989, vrr_stack+0, vrr_stack+289, NULL);
 _BUILD_00g0(Data,vrr_stack+0, vrr_stack+292, vrr_stack+734, vrr_stack+49, vrr_stack+15, NULL);
 _BUILD_00h0(Data,vrr_stack+49, vrr_stack+302, vrr_stack+0, vrr_stack+334, vrr_stack+292, NULL);
 _BUILD_00i0(Data,vrr_stack+261, vrr_stack+240, vrr_stack+49, vrr_stack+344, vrr_stack+302, NULL);
 _BUILD_p0i0(Data,vrr_stack+744, vrr_stack+359, vrr_stack+261, NULL, NULL, vrr_stack+240);
 _BUILD_d0i0(Data,vrr_stack+2219, vrr_stack+387, vrr_stack+744, vrr_stack+1001, vrr_stack+359, vrr_stack+1190);
 _BUILD_f0i0(Data,vrr_stack+2387, vrr_stack+1799, vrr_stack+2219, vrr_stack+471, vrr_stack+387, vrr_stack+1253);
 _BUILD_f0i0(Data,vrr_stack+2667, vrr_stack+2051, vrr_stack+1799, vrr_stack+1967, vrr_stack+471, vrr_stack+860);
   tmp = vrr_stack + 2667;
   target_ptr = Libint->vrr_classes[3][6];
   for(i=0;i<280;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1967, vrr_stack+1001, vrr_stack+359, vrr_stack+84, vrr_stack+126, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+2003, vrr_stack+1029, vrr_stack+1001, vrr_stack+105, vrr_stack+84, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1190, vrr_stack+359, vrr_stack+261, vrr_stack+126, vrr_stack+240, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+1226, vrr_stack+1967, vrr_stack+1190, NULL, NULL, vrr_stack+359);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+2947, vrr_stack+2003, vrr_stack+1967, NULL, NULL, vrr_stack+1001);
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3055, vrr_stack+2947, vrr_stack+1226, vrr_stack+2003, vrr_stack+1967, vrr_stack+387);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+387, vrr_stack+21, vrr_stack+1029, vrr_stack+219, vrr_stack+105, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+70, vrr_stack+387, vrr_stack+2003, NULL, NULL, vrr_stack+1029);
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3271, vrr_stack+70, vrr_stack+2947, vrr_stack+387, vrr_stack+2003, vrr_stack+471);
   tmp = vrr_stack + 3271;
   target_ptr = Libint->vrr_classes[2][7];
   for(i=0;i<216;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+423, Data->F+10, Data->F+11, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+995, vrr_stack+986, vrr_stack+423, Data->F+9, Data->F+10, NULL);
 _BUILD_00f0(Data,vrr_stack+426, vrr_stack+989, vrr_stack+995, vrr_stack+289, vrr_stack+986, NULL);
 _BUILD_00g0(Data,vrr_stack+436, vrr_stack+734, vrr_stack+426, vrr_stack+15, vrr_stack+989, NULL);
 _BUILD_00h0(Data,vrr_stack+451, vrr_stack+0, vrr_stack+436, vrr_stack+292, vrr_stack+734, NULL);
 _BUILD_00i0(Data,vrr_stack+472, vrr_stack+49, vrr_stack+451, vrr_stack+302, vrr_stack+0, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+289, vrr_stack+261, vrr_stack+472, vrr_stack+240, vrr_stack+49, NULL);
 am[0] = 1;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3487, vrr_stack+1190, vrr_stack+289, NULL, NULL, vrr_stack+261);
 am[0] = 2;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3595, vrr_stack+1226, vrr_stack+3487, vrr_stack+1967, vrr_stack+1190, vrr_stack+744);
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+3811, vrr_stack+3055, vrr_stack+3595, vrr_stack+2947, vrr_stack+1226, vrr_stack+2219);
 am[0] = 3;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+4171, vrr_stack+3271, vrr_stack+3055, vrr_stack+70, vrr_stack+2947, vrr_stack+1799);
   tmp = vrr_stack + 4171;
   target_ptr = Libint->vrr_classes[3][7];
   for(i=0;i<360;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1334, vrr_stack+1967, vrr_stack+1190, vrr_stack+1001, vrr_stack+359, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+70, vrr_stack+2003, vrr_stack+1967, vrr_stack+1029, vrr_stack+1001, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+115, vrr_stack+1190, vrr_stack+289, vrr_stack+359, vrr_stack+261, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+2219, vrr_stack+1334, vrr_stack+115, NULL, NULL, vrr_stack+1190);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+4531, vrr_stack+70, vrr_stack+1334, NULL, NULL, vrr_stack+1967);
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+4666, vrr_stack+4531, vrr_stack+2219, vrr_stack+70, vrr_stack+1334, vrr_stack+1226);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+1226, vrr_stack+387, vrr_stack+2003, vrr_stack+21, vrr_stack+1029, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+4936, vrr_stack+1226, vrr_stack+70, NULL, NULL, vrr_stack+2003);
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+5071, vrr_stack+4936, vrr_stack+4531, vrr_stack+1226, vrr_stack+70, vrr_stack+2947);
   tmp = vrr_stack + 5071;
   target_ptr = Libint->vrr_classes[2][8];
   for(i=0;i<270;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+2947, Data->F+11, Data->F+12, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+2950, vrr_stack+423, vrr_stack+2947, Data->F+10, Data->F+11, NULL);
 _BUILD_00f0(Data,vrr_stack+2956, vrr_stack+995, vrr_stack+2950, vrr_stack+986, vrr_stack+423, NULL);
 _BUILD_00g0(Data,vrr_stack+2966, vrr_stack+426, vrr_stack+2956, vrr_stack+989, vrr_stack+995, NULL);
 _BUILD_00h0(Data,vrr_stack+2981, vrr_stack+436, vrr_stack+2966, vrr_stack+734, vrr_stack+426, NULL);
 _BUILD_00i0(Data,vrr_stack+734, vrr_stack+451, vrr_stack+2981, vrr_stack+0, vrr_stack+436, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+0, vrr_stack+472, vrr_stack+734, vrr_stack+49, vrr_stack+451, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+762, vrr_stack+289, vrr_stack+0, vrr_stack+261, vrr_stack+472, NULL);
 am[0] = 1;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+5341, vrr_stack+115, vrr_stack+762, NULL, NULL, vrr_stack+289);
 am[0] = 2;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+5476, vrr_stack+2219, vrr_stack+5341, vrr_stack+1334, vrr_stack+115, vrr_stack+3487);
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+5746, vrr_stack+4666, vrr_stack+5476, vrr_stack+4531, vrr_stack+2219, vrr_stack+3595);
 am[0] = 3;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+6196, vrr_stack+5071, vrr_stack+4666, vrr_stack+4936, vrr_stack+4531, vrr_stack+3055);
   tmp = vrr_stack + 6196;
   target_ptr = Libint->vrr_classes[3][8];
   for(i=0;i<450;i++)
     target_ptr[i] += tmp[i];
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+4936, vrr_stack+1334, vrr_stack+115, vrr_stack+1967, vrr_stack+1190, NULL);
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+4991, vrr_stack+70, vrr_stack+1334, vrr_stack+2003, vrr_stack+1967, NULL);
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+3487, vrr_stack+115, vrr_stack+762, vrr_stack+1190, vrr_stack+289, NULL);
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+3542, vrr_stack+4936, vrr_stack+3487, NULL, NULL, vrr_stack+115);
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+115, vrr_stack+4991, vrr_stack+4936, NULL, NULL, vrr_stack+1334);
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+6646, vrr_stack+115, vrr_stack+3542, vrr_stack+4991, vrr_stack+4936, vrr_stack+2219);
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+2219, vrr_stack+1226, vrr_stack+70, vrr_stack+387, vrr_stack+2003, NULL);
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+1190, vrr_stack+2219, vrr_stack+4991, NULL, NULL, vrr_stack+70);
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+6976, vrr_stack+1190, vrr_stack+115, vrr_stack+2219, vrr_stack+4991, vrr_stack+4531);
   tmp = vrr_stack + 6976;
   target_ptr = Libint->vrr_classes[2][9];
   for(i=0;i<330;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00p0(Data,vrr_stack+4531, Data->F+12, Data->F+13, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+4534, vrr_stack+2947, vrr_stack+4531, Data->F+11, Data->F+12, NULL);
 _BUILD_00f0(Data,vrr_stack+4540, vrr_stack+2950, vrr_stack+4534, vrr_stack+423, vrr_stack+2947, NULL);
 _BUILD_00g0(Data,vrr_stack+4550, vrr_stack+2956, vrr_stack+4540, vrr_stack+995, vrr_stack+2950, NULL);
 _BUILD_00h0(Data,vrr_stack+4565, vrr_stack+2966, vrr_stack+4550, vrr_stack+426, vrr_stack+2956, NULL);
 _BUILD_00i0(Data,vrr_stack+4531, vrr_stack+2981, vrr_stack+4565, vrr_stack+436, vrr_stack+2966, NULL);
 am[0] = 0;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+4559, vrr_stack+734, vrr_stack+4531, vrr_stack+451, vrr_stack+2981, NULL);
 am[0] = 0;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+2947, vrr_stack+0, vrr_stack+4559, vrr_stack+472, vrr_stack+734, NULL);
 am[0] = 0;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+4531, vrr_stack+762, vrr_stack+2947, vrr_stack+289, vrr_stack+0, NULL);
 am[0] = 1;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+2219, vrr_stack+3487, vrr_stack+4531, NULL, NULL, vrr_stack+762);
 am[0] = 2;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+7306, vrr_stack+3542, vrr_stack+2219, vrr_stack+4936, vrr_stack+3487, vrr_stack+5341);
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+7636, vrr_stack+6646, vrr_stack+7306, vrr_stack+115, vrr_stack+3542, vrr_stack+5476);
 am[0] = 3;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+8186, vrr_stack+6976, vrr_stack+6646, vrr_stack+1190, vrr_stack+115, vrr_stack+4666);
   tmp = vrr_stack + 8186;
   target_ptr = Libint->vrr_classes[3][9];
   for(i=0;i<550;i++)
     target_ptr[i] += tmp[i];
 _BUILD_g0h0(Data,vrr_stack+7306, vrr_stack+1589, vrr_stack+1379, vrr_stack+1064, vrr_stack+860, vrr_stack+584);
   tmp = vrr_stack + 7306;
   target_ptr = Libint->vrr_classes[4][5];
   for(i=0;i<315;i++)
     target_ptr[i] += tmp[i];
 _BUILD_g0i0(Data,vrr_stack+0, vrr_stack+2667, vrr_stack+2387, vrr_stack+2051, vrr_stack+1799, vrr_stack+1379);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[4][6];
   for(i=0;i<420;i++)
     target_ptr[i] += tmp[i];
 am[0] = 4;  am[1] = 7;
 vrr_build_xxxx(am,Data,vrr_stack+420, vrr_stack+4171, vrr_stack+3811, vrr_stack+3271, vrr_stack+3055, vrr_stack+2387);
   tmp = vrr_stack + 420;
   target_ptr = Libint->vrr_classes[4][7];
   for(i=0;i<540;i++)
     target_ptr[i] += tmp[i];
 am[0] = 4;  am[1] = 8;
 vrr_build_xxxx(am,Data,vrr_stack+960, vrr_stack+6196, vrr_stack+5746, vrr_stack+5071, vrr_stack+4666, vrr_stack+3811);
   tmp = vrr_stack + 960;
   target_ptr = Libint->vrr_classes[4][8];
   for(i=0;i<675;i++)
     target_ptr[i] += tmp[i];
 am[0] = 4;  am[1] = 9;
 vrr_build_xxxx(am,Data,vrr_stack+1635, vrr_stack+8186, vrr_stack+7636, vrr_stack+6976, vrr_stack+6646, vrr_stack+5746);
   tmp = vrr_stack + 1635;
   target_ptr = Libint->vrr_classes[4][9];
   for(i=0;i<825;i++)
     target_ptr[i] += tmp[i];

}

