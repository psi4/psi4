#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (gp|hp) integrals */

void vrr_order_gphp(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_00p0(Data,vrr_stack+0, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, NULL, NULL, Data->F+5);
 _BUILD_00d0(Data,vrr_stack+15, vrr_stack+3, vrr_stack+0, Data->F+4, Data->F+5, NULL);
 _BUILD_00p0(Data,vrr_stack+21, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+24, vrr_stack+21, vrr_stack+3, Data->F+3, Data->F+4, NULL);
 _BUILD_00p0(Data,vrr_stack+30, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+33, vrr_stack+0, vrr_stack+30, Data->F+5, Data->F+6, NULL);
 _BUILD_p0d0(Data,vrr_stack+39, vrr_stack+15, vrr_stack+33, NULL, NULL, vrr_stack+0);
 _BUILD_p0d0(Data,vrr_stack+57, vrr_stack+24, vrr_stack+15, NULL, NULL, vrr_stack+3);
 _BUILD_d0d0(Data,vrr_stack+75, vrr_stack+57, vrr_stack+39, vrr_stack+24, vrr_stack+15, vrr_stack+6);
 _BUILD_00f0(Data,vrr_stack+111, vrr_stack+15, vrr_stack+33, vrr_stack+3, vrr_stack+0, NULL);
 _BUILD_00f0(Data,vrr_stack+121, vrr_stack+24, vrr_stack+15, vrr_stack+21, vrr_stack+3, NULL);
 _BUILD_p0f0(Data,vrr_stack+131, vrr_stack+121, vrr_stack+111, NULL, NULL, vrr_stack+15);
 _BUILD_00p0(Data,vrr_stack+3, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+3, vrr_stack+21, Data->F+2, Data->F+3, NULL);
 _BUILD_00f0(Data,vrr_stack+161, vrr_stack+6, vrr_stack+24, vrr_stack+3, vrr_stack+21, NULL);
 _BUILD_p0f0(Data,vrr_stack+171, vrr_stack+161, vrr_stack+121, NULL, NULL, vrr_stack+24);
 _BUILD_00p0(Data,vrr_stack+21, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+201, vrr_stack+30, vrr_stack+21, Data->F+6, Data->F+7, NULL);
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
 _BUILD_00p0(Data,vrr_stack+0, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+24, vrr_stack+0, vrr_stack+3, Data->F+1, Data->F+2, NULL);
 _BUILD_00f0(Data,vrr_stack+131, vrr_stack+24, vrr_stack+6, vrr_stack+0, vrr_stack+3, NULL);
 _BUILD_00g0(Data,vrr_stack+141, vrr_stack+131, vrr_stack+161, vrr_stack+24, vrr_stack+6, NULL);
 _BUILD_p0g0(Data,vrr_stack+647, vrr_stack+141, vrr_stack+54, NULL, NULL, vrr_stack+161);
 _BUILD_d0g0(Data,vrr_stack+692, vrr_stack+647, vrr_stack+512, vrr_stack+141, vrr_stack+54, vrr_stack+171);
 _BUILD_00p0(Data,vrr_stack+171, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+174, vrr_stack+21, vrr_stack+171, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+180, vrr_stack+201, vrr_stack+174, vrr_stack+30, vrr_stack+21, NULL);
 _BUILD_00g0(Data,vrr_stack+3, vrr_stack+207, vrr_stack+180, vrr_stack+33, vrr_stack+201, NULL);
 _BUILD_p0g0(Data,vrr_stack+782, vrr_stack+69, vrr_stack+3, NULL, NULL, vrr_stack+207);
 _BUILD_d0g0(Data,vrr_stack+827, vrr_stack+467, vrr_stack+782, vrr_stack+39, vrr_stack+69, vrr_stack+217);
 _BUILD_f0g0(Data,vrr_stack+917, vrr_stack+557, vrr_stack+827, vrr_stack+512, vrr_stack+467, vrr_stack+247);
 _BUILD_f0g0(Data,vrr_stack+1067, vrr_stack+692, vrr_stack+557, vrr_stack+647, vrr_stack+512, vrr_stack+307);
 _BUILD_g0g0(Data,vrr_stack+1217, vrr_stack+1067, vrr_stack+917, vrr_stack+692, vrr_stack+557, vrr_stack+367);
 _BUILD_00h0(Data,vrr_stack+217, vrr_stack+39, vrr_stack+69, vrr_stack+121, vrr_stack+111, NULL);
 _BUILD_00h0(Data,vrr_stack+238, vrr_stack+54, vrr_stack+39, vrr_stack+161, vrr_stack+121, NULL);
 _BUILD_p0h0(Data,vrr_stack+259, vrr_stack+238, vrr_stack+217, NULL, NULL, vrr_stack+39);
 _BUILD_00h0(Data,vrr_stack+322, vrr_stack+141, vrr_stack+54, vrr_stack+131, vrr_stack+161, NULL);
 _BUILD_p0h0(Data,vrr_stack+343, vrr_stack+322, vrr_stack+238, NULL, NULL, vrr_stack+54);
 _BUILD_00h0(Data,vrr_stack+406, vrr_stack+69, vrr_stack+3, vrr_stack+111, vrr_stack+207, NULL);
 _BUILD_p0h0(Data,vrr_stack+1442, vrr_stack+217, vrr_stack+406, NULL, NULL, vrr_stack+69);
 _BUILD_d0h0(Data,vrr_stack+1505, vrr_stack+259, vrr_stack+1442, vrr_stack+238, vrr_stack+217, vrr_stack+467);
 _BUILD_d0h0(Data,vrr_stack+1631, vrr_stack+343, vrr_stack+259, vrr_stack+322, vrr_stack+238, vrr_stack+512);
 _BUILD_f0h0(Data,vrr_stack+1757, vrr_stack+1631, vrr_stack+1505, vrr_stack+343, vrr_stack+259, vrr_stack+557);
 _BUILD_00p0(Data,vrr_stack+18, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+427, vrr_stack+18, vrr_stack+0, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+433, vrr_stack+427, vrr_stack+24, vrr_stack+18, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+156, vrr_stack+433, vrr_stack+131, vrr_stack+427, vrr_stack+24, NULL);
 _BUILD_00h0(Data,vrr_stack+443, vrr_stack+156, vrr_stack+141, vrr_stack+433, vrr_stack+131, NULL);
 _BUILD_p0h0(Data,vrr_stack+464, vrr_stack+443, vrr_stack+322, NULL, NULL, vrr_stack+141);
 _BUILD_d0h0(Data,vrr_stack+1967, vrr_stack+464, vrr_stack+343, vrr_stack+443, vrr_stack+322, vrr_stack+647);
 _BUILD_f0h0(Data,vrr_stack+2093, vrr_stack+1967, vrr_stack+1631, vrr_stack+464, vrr_stack+343, vrr_stack+692);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+9, Data->F+10, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+464, vrr_stack+171, vrr_stack+0, Data->F+8, Data->F+9, NULL);
 _BUILD_00f0(Data,vrr_stack+470, vrr_stack+174, vrr_stack+464, vrr_stack+21, vrr_stack+171, NULL);
 _BUILD_00g0(Data,vrr_stack+480, vrr_stack+180, vrr_stack+470, vrr_stack+201, vrr_stack+174, NULL);
 _BUILD_00h0(Data,vrr_stack+18, vrr_stack+3, vrr_stack+480, vrr_stack+207, vrr_stack+180, NULL);
 _BUILD_p0h0(Data,vrr_stack+495, vrr_stack+406, vrr_stack+18, NULL, NULL, vrr_stack+3);
 _BUILD_d0h0(Data,vrr_stack+558, vrr_stack+1442, vrr_stack+495, vrr_stack+217, vrr_stack+406, vrr_stack+782);
 _BUILD_f0h0(Data,vrr_stack+2303, vrr_stack+1505, vrr_stack+558, vrr_stack+259, vrr_stack+1442, vrr_stack+827);
 _BUILD_g0h0(Data,vrr_stack+2513, vrr_stack+1757, vrr_stack+2303, vrr_stack+1631, vrr_stack+1505, vrr_stack+917);
 _BUILD_g0h0(Data,vrr_stack+684, vrr_stack+2093, vrr_stack+1757, vrr_stack+1967, vrr_stack+1631, vrr_stack+1067);
   tmp = vrr_stack + 684;
   target_ptr = Libint->vrr_classes[4][5];
   for(i=0;i<315;i++)
     target_ptr[i] += tmp[i];
 _BUILD_00i0(Data,vrr_stack+1967, vrr_stack+217, vrr_stack+406, vrr_stack+39, vrr_stack+69, NULL);
 _BUILD_00i0(Data,vrr_stack+1995, vrr_stack+238, vrr_stack+217, vrr_stack+54, vrr_stack+39, NULL);
 _BUILD_p0i0(Data,vrr_stack+999, vrr_stack+1995, vrr_stack+1967, NULL, NULL, vrr_stack+217);
 _BUILD_00i0(Data,vrr_stack+2023, vrr_stack+322, vrr_stack+238, vrr_stack+141, vrr_stack+54, NULL);
 _BUILD_p0i0(Data,vrr_stack+1083, vrr_stack+2023, vrr_stack+1995, NULL, NULL, vrr_stack+238);
 _BUILD_00i0(Data,vrr_stack+39, vrr_stack+406, vrr_stack+18, vrr_stack+69, vrr_stack+3, NULL);
 _BUILD_p0i0(Data,vrr_stack+2828, vrr_stack+1967, vrr_stack+39, NULL, NULL, vrr_stack+406);
 _BUILD_d0i0(Data,vrr_stack+2912, vrr_stack+999, vrr_stack+2828, vrr_stack+1995, vrr_stack+1967, vrr_stack+1442);
 _BUILD_d0i0(Data,vrr_stack+3080, vrr_stack+1083, vrr_stack+999, vrr_stack+2023, vrr_stack+1995, vrr_stack+259);
 _BUILD_f0i0(Data,vrr_stack+3248, vrr_stack+3080, vrr_stack+2912, vrr_stack+1083, vrr_stack+999, vrr_stack+1505);
 _BUILD_00i0(Data,vrr_stack+1995, vrr_stack+443, vrr_stack+322, vrr_stack+156, vrr_stack+141, NULL);
 _BUILD_p0i0(Data,vrr_stack+1442, vrr_stack+1995, vrr_stack+2023, NULL, NULL, vrr_stack+322);
 _BUILD_d0i0(Data,vrr_stack+3528, vrr_stack+1442, vrr_stack+1083, vrr_stack+1995, vrr_stack+2023, vrr_stack+343);
 _BUILD_f0i0(Data,vrr_stack+3696, vrr_stack+3528, vrr_stack+3080, vrr_stack+1442, vrr_stack+1083, vrr_stack+1631);
 _BUILD_00p0(Data,vrr_stack+1083, Data->F+10, Data->F+11, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+1086, vrr_stack+0, vrr_stack+1083, Data->F+9, Data->F+10, NULL);
 _BUILD_00f0(Data,vrr_stack+1092, vrr_stack+464, vrr_stack+1086, vrr_stack+171, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+1102, vrr_stack+470, vrr_stack+1092, vrr_stack+174, vrr_stack+464, NULL);
 _BUILD_00h0(Data,vrr_stack+1117, vrr_stack+480, vrr_stack+1102, vrr_stack+180, vrr_stack+470, NULL);
 _BUILD_00i0(Data,vrr_stack+1083, vrr_stack+18, vrr_stack+1117, vrr_stack+3, vrr_stack+480, NULL);
 _BUILD_p0i0(Data,vrr_stack+1111, vrr_stack+39, vrr_stack+1083, NULL, NULL, vrr_stack+18);
 _BUILD_d0i0(Data,vrr_stack+1442, vrr_stack+2828, vrr_stack+1111, vrr_stack+1967, vrr_stack+39, vrr_stack+495);
 _BUILD_f0i0(Data,vrr_stack+0, vrr_stack+2912, vrr_stack+1442, vrr_stack+999, vrr_stack+2828, vrr_stack+558);
 _BUILD_g0i0(Data,vrr_stack+3976, vrr_stack+3248, vrr_stack+0, vrr_stack+3080, vrr_stack+2912, vrr_stack+2303);
 _BUILD_g0i0(Data,vrr_stack+0, vrr_stack+3696, vrr_stack+3248, vrr_stack+3528, vrr_stack+3080, vrr_stack+1757);
   tmp = vrr_stack + 0;
   target_ptr = Libint->vrr_classes[4][6];
   for(i=0;i<420;i++)
     target_ptr[i] += tmp[i];
 _BUILD_h0h0(Data,vrr_stack+4396, vrr_stack+684, vrr_stack+2513, vrr_stack+2093, vrr_stack+1757, vrr_stack+1217);
   tmp = vrr_stack + 4396;
   target_ptr = Libint->vrr_classes[5][5];
   for(i=0;i<441;i++)
     target_ptr[i] += tmp[i];
 _BUILD_h0i0(Data,vrr_stack+420, vrr_stack+0, vrr_stack+3976, vrr_stack+3696, vrr_stack+3248, vrr_stack+2513);
   tmp = vrr_stack + 420;
   target_ptr = Libint->vrr_classes[5][6];
   for(i=0;i<588;i++)
     target_ptr[i] += tmp[i];

}

