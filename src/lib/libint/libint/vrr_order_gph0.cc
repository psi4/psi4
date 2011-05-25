#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (gp|h0) integrals */

void vrr_order_gph0(Libint_t * Libint, prim_data *Data)
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
 _BUILD_00d0(Data,vrr_stack+238, vrr_stack+18, vrr_stack+0, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+244, vrr_stack+238, vrr_stack+24, vrr_stack+18, vrr_stack+0, NULL);
 _BUILD_00g0(Data,vrr_stack+156, vrr_stack+244, vrr_stack+131, vrr_stack+238, vrr_stack+24, NULL);
 _BUILD_00h0(Data,vrr_stack+24, vrr_stack+156, vrr_stack+141, vrr_stack+244, vrr_stack+131, NULL);
 _BUILD_p0h0(Data,vrr_stack+45, vrr_stack+24, vrr_stack+322, NULL, NULL, vrr_stack+141);
 _BUILD_d0h0(Data,vrr_stack+427, vrr_stack+45, vrr_stack+343, vrr_stack+24, vrr_stack+322, vrr_stack+647);
 _BUILD_f0h0(Data,vrr_stack+1967, vrr_stack+427, vrr_stack+1631, vrr_stack+45, vrr_stack+343, vrr_stack+692);
 _BUILD_00p0(Data,vrr_stack+0, Data->F+9, Data->F+10, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+322, vrr_stack+171, vrr_stack+0, Data->F+8, Data->F+9, NULL);
 _BUILD_00f0(Data,vrr_stack+328, vrr_stack+174, vrr_stack+322, vrr_stack+21, vrr_stack+171, NULL);
 _BUILD_00g0(Data,vrr_stack+338, vrr_stack+180, vrr_stack+328, vrr_stack+201, vrr_stack+174, NULL);
 _BUILD_00h0(Data,vrr_stack+238, vrr_stack+3, vrr_stack+338, vrr_stack+207, vrr_stack+180, NULL);
 _BUILD_p0h0(Data,vrr_stack+322, vrr_stack+406, vrr_stack+238, NULL, NULL, vrr_stack+3);
 _BUILD_d0h0(Data,vrr_stack+0, vrr_stack+1442, vrr_stack+322, vrr_stack+217, vrr_stack+406, vrr_stack+782);
 _BUILD_f0h0(Data,vrr_stack+553, vrr_stack+1505, vrr_stack+0, vrr_stack+259, vrr_stack+1442, vrr_stack+827);
 _BUILD_g0h0(Data,vrr_stack+0, vrr_stack+1757, vrr_stack+553, vrr_stack+1631, vrr_stack+1505, vrr_stack+917);
 _BUILD_g0h0(Data,vrr_stack+553, vrr_stack+1967, vrr_stack+1757, vrr_stack+427, vrr_stack+1631, vrr_stack+1067);
   tmp = vrr_stack + 553;
   target_ptr = Libint->vrr_classes[4][5];
   for(i=0;i<315;i++)
     target_ptr[i] += tmp[i];
 _BUILD_h0h0(Data,vrr_stack+2177, vrr_stack+553, vrr_stack+0, vrr_stack+1967, vrr_stack+1757, vrr_stack+1217);
   tmp = vrr_stack + 2177;
   target_ptr = Libint->vrr_classes[5][5];
   for(i=0;i<441;i++)
     target_ptr[i] += tmp[i];

}

