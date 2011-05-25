#include <stdio.h>
#include "libint.h"
#include "vrr_header.h"

  /* Computes quartets necessary to compute (h0|fd) integrals */

void vrr_order_h0fd(Libint_t * Libint, prim_data *Data)
{
 REALTYPE *vrr_stack = Libint->vrr_stack;
 REALTYPE *tmp, *target_ptr;
 int i, am[2];
 _BUILD_p000(Data,vrr_stack+0, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_p000(Data,vrr_stack+3, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_d000(Data,vrr_stack+6, vrr_stack+3, vrr_stack+0, Data->F+3, Data->F+4, NULL);
 _BUILD_00p0(Data,vrr_stack+12, Data->F+4, Data->F+5, NULL, NULL, NULL);
 _BUILD_00p0(Data,vrr_stack+15, Data->F+3, Data->F+4, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+18, vrr_stack+15, vrr_stack+12, NULL, NULL, Data->F+4);
 _BUILD_00p0(Data,vrr_stack+27, Data->F+2, Data->F+3, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+30, vrr_stack+27, vrr_stack+15, NULL, NULL, Data->F+3);
 _BUILD_00p0(Data,vrr_stack+39, Data->F+5, Data->F+6, NULL, NULL, NULL);
 _BUILD_p0p0(Data,vrr_stack+42, vrr_stack+12, vrr_stack+39, NULL, NULL, Data->F+5);
 _BUILD_d0p0(Data,vrr_stack+51, vrr_stack+18, vrr_stack+42, vrr_stack+15, vrr_stack+12, vrr_stack+0);
 _BUILD_d0p0(Data,vrr_stack+69, vrr_stack+30, vrr_stack+18, vrr_stack+27, vrr_stack+15, vrr_stack+3);
 _BUILD_f0p0(Data,vrr_stack+87, vrr_stack+69, vrr_stack+51, vrr_stack+30, vrr_stack+18, vrr_stack+6);
 _BUILD_00d0(Data,vrr_stack+0, vrr_stack+15, vrr_stack+12, Data->F+3, Data->F+4, NULL);
 _BUILD_00d0(Data,vrr_stack+6, vrr_stack+27, vrr_stack+15, Data->F+2, Data->F+3, NULL);
 _BUILD_00d0(Data,vrr_stack+117, vrr_stack+12, vrr_stack+39, Data->F+4, Data->F+5, NULL);
 _BUILD_p0d0(Data,vrr_stack+123, vrr_stack+0, vrr_stack+117, NULL, NULL, vrr_stack+12);
 _BUILD_p0d0(Data,vrr_stack+141, vrr_stack+6, vrr_stack+0, NULL, NULL, vrr_stack+15);
 _BUILD_d0d0(Data,vrr_stack+159, vrr_stack+141, vrr_stack+123, vrr_stack+6, vrr_stack+0, vrr_stack+18);
 _BUILD_00p0(Data,vrr_stack+18, Data->F+1, Data->F+2, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+21, vrr_stack+18, vrr_stack+27, Data->F+1, Data->F+2, NULL);
 _BUILD_p0d0(Data,vrr_stack+195, vrr_stack+21, vrr_stack+6, NULL, NULL, vrr_stack+27);
 _BUILD_d0d0(Data,vrr_stack+213, vrr_stack+195, vrr_stack+141, vrr_stack+21, vrr_stack+6, vrr_stack+30);
 _BUILD_00p0(Data,vrr_stack+30, Data->F+6, Data->F+7, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+33, vrr_stack+39, vrr_stack+30, Data->F+5, Data->F+6, NULL);
 _BUILD_p0d0(Data,vrr_stack+249, vrr_stack+117, vrr_stack+33, NULL, NULL, vrr_stack+39);
 _BUILD_d0d0(Data,vrr_stack+267, vrr_stack+123, vrr_stack+249, vrr_stack+0, vrr_stack+117, vrr_stack+42);
 _BUILD_f0d0(Data,vrr_stack+303, vrr_stack+159, vrr_stack+267, vrr_stack+141, vrr_stack+123, vrr_stack+51);
 _BUILD_f0d0(Data,vrr_stack+363, vrr_stack+213, vrr_stack+159, vrr_stack+195, vrr_stack+141, vrr_stack+69);
 _BUILD_g0d0(Data,vrr_stack+423, vrr_stack+363, vrr_stack+303, vrr_stack+213, vrr_stack+159, vrr_stack+87);
 _BUILD_00f0(Data,vrr_stack+42, vrr_stack+0, vrr_stack+117, vrr_stack+15, vrr_stack+12, NULL);
 _BUILD_00f0(Data,vrr_stack+52, vrr_stack+6, vrr_stack+0, vrr_stack+27, vrr_stack+15, NULL);
 _BUILD_p0f0(Data,vrr_stack+62, vrr_stack+52, vrr_stack+42, NULL, NULL, vrr_stack+0);
 _BUILD_00f0(Data,vrr_stack+92, vrr_stack+21, vrr_stack+6, vrr_stack+18, vrr_stack+27, NULL);
 _BUILD_p0f0(Data,vrr_stack+513, vrr_stack+92, vrr_stack+52, NULL, NULL, vrr_stack+6);
 _BUILD_00f0(Data,vrr_stack+102, vrr_stack+117, vrr_stack+33, vrr_stack+12, vrr_stack+39, NULL);
 _BUILD_p0f0(Data,vrr_stack+543, vrr_stack+42, vrr_stack+102, NULL, NULL, vrr_stack+117);
 _BUILD_d0f0(Data,vrr_stack+573, vrr_stack+62, vrr_stack+543, vrr_stack+52, vrr_stack+42, vrr_stack+123);
 _BUILD_d0f0(Data,vrr_stack+633, vrr_stack+513, vrr_stack+62, vrr_stack+92, vrr_stack+52, vrr_stack+141);
 _BUILD_f0f0(Data,vrr_stack+693, vrr_stack+633, vrr_stack+573, vrr_stack+513, vrr_stack+62, vrr_stack+159);
 _BUILD_00p0(Data,vrr_stack+27, Data->F+0, Data->F+1, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+12, vrr_stack+27, vrr_stack+18, Data->F+0, Data->F+1, NULL);
 _BUILD_00f0(Data,vrr_stack+123, vrr_stack+12, vrr_stack+21, vrr_stack+27, vrr_stack+18, NULL);
 _BUILD_p0f0(Data,vrr_stack+133, vrr_stack+123, vrr_stack+92, NULL, NULL, vrr_stack+21);
 _BUILD_d0f0(Data,vrr_stack+793, vrr_stack+133, vrr_stack+513, vrr_stack+123, vrr_stack+92, vrr_stack+195);
 _BUILD_f0f0(Data,vrr_stack+853, vrr_stack+793, vrr_stack+633, vrr_stack+133, vrr_stack+513, vrr_stack+213);
 _BUILD_00p0(Data,vrr_stack+18, Data->F+7, Data->F+8, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+133, vrr_stack+30, vrr_stack+18, Data->F+6, Data->F+7, NULL);
 _BUILD_00f0(Data,vrr_stack+139, vrr_stack+33, vrr_stack+133, vrr_stack+39, vrr_stack+30, NULL);
 _BUILD_p0f0(Data,vrr_stack+149, vrr_stack+102, vrr_stack+139, NULL, NULL, vrr_stack+33);
 _BUILD_d0f0(Data,vrr_stack+179, vrr_stack+543, vrr_stack+149, vrr_stack+42, vrr_stack+102, vrr_stack+249);
 _BUILD_f0f0(Data,vrr_stack+953, vrr_stack+573, vrr_stack+179, vrr_stack+62, vrr_stack+543, vrr_stack+267);
 _BUILD_g0f0(Data,vrr_stack+1053, vrr_stack+693, vrr_stack+953, vrr_stack+633, vrr_stack+573, vrr_stack+303);
 _BUILD_g0f0(Data,vrr_stack+1203, vrr_stack+853, vrr_stack+693, vrr_stack+793, vrr_stack+633, vrr_stack+363);
 _BUILD_00g0(Data,vrr_stack+793, vrr_stack+42, vrr_stack+102, vrr_stack+0, vrr_stack+117, NULL);
 _BUILD_00g0(Data,vrr_stack+808, vrr_stack+52, vrr_stack+42, vrr_stack+6, vrr_stack+0, NULL);
 _BUILD_p0g0(Data,vrr_stack+239, vrr_stack+808, vrr_stack+793, NULL, NULL, vrr_stack+42);
 _BUILD_00g0(Data,vrr_stack+823, vrr_stack+92, vrr_stack+52, vrr_stack+21, vrr_stack+6, NULL);
 _BUILD_p0g0(Data,vrr_stack+284, vrr_stack+823, vrr_stack+808, NULL, NULL, vrr_stack+52);
 _BUILD_00g0(Data,vrr_stack+838, vrr_stack+102, vrr_stack+139, vrr_stack+117, vrr_stack+33, NULL);
 _BUILD_p0g0(Data,vrr_stack+329, vrr_stack+793, vrr_stack+838, NULL, NULL, vrr_stack+102);
 _BUILD_d0g0(Data,vrr_stack+1353, vrr_stack+239, vrr_stack+329, vrr_stack+808, vrr_stack+793, vrr_stack+543);
 _BUILD_d0g0(Data,vrr_stack+1443, vrr_stack+284, vrr_stack+239, vrr_stack+823, vrr_stack+808, vrr_stack+62);
 _BUILD_f0g0(Data,vrr_stack+1533, vrr_stack+1443, vrr_stack+1353, vrr_stack+284, vrr_stack+239, vrr_stack+573);
 _BUILD_00g0(Data,vrr_stack+62, vrr_stack+123, vrr_stack+92, vrr_stack+12, vrr_stack+21, NULL);
 _BUILD_p0g0(Data,vrr_stack+543, vrr_stack+62, vrr_stack+823, NULL, NULL, vrr_stack+92);
 _BUILD_d0g0(Data,vrr_stack+1683, vrr_stack+543, vrr_stack+284, vrr_stack+62, vrr_stack+823, vrr_stack+513);
 _BUILD_f0g0(Data,vrr_stack+1773, vrr_stack+1683, vrr_stack+1443, vrr_stack+543, vrr_stack+284, vrr_stack+633);
 _BUILD_00p0(Data,vrr_stack+39, Data->F+8, Data->F+9, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+513, vrr_stack+18, vrr_stack+39, Data->F+7, Data->F+8, NULL);
 _BUILD_00f0(Data,vrr_stack+519, vrr_stack+133, vrr_stack+513, vrr_stack+30, vrr_stack+18, NULL);
 _BUILD_00g0(Data,vrr_stack+77, vrr_stack+139, vrr_stack+519, vrr_stack+33, vrr_stack+133, NULL);
 _BUILD_p0g0(Data,vrr_stack+529, vrr_stack+838, vrr_stack+77, NULL, NULL, vrr_stack+139);
 _BUILD_d0g0(Data,vrr_stack+574, vrr_stack+329, vrr_stack+529, vrr_stack+793, vrr_stack+838, vrr_stack+149);
 _BUILD_f0g0(Data,vrr_stack+1923, vrr_stack+1353, vrr_stack+574, vrr_stack+239, vrr_stack+329, vrr_stack+179);
 _BUILD_g0g0(Data,vrr_stack+2073, vrr_stack+1533, vrr_stack+1923, vrr_stack+1443, vrr_stack+1353, vrr_stack+953);
 _BUILD_g0g0(Data,vrr_stack+2298, vrr_stack+1773, vrr_stack+1533, vrr_stack+1683, vrr_stack+1443, vrr_stack+693);
 _BUILD_00h0(Data,vrr_stack+1683, vrr_stack+793, vrr_stack+838, vrr_stack+42, vrr_stack+102, NULL);
 _BUILD_00h0(Data,vrr_stack+1704, vrr_stack+808, vrr_stack+793, vrr_stack+52, vrr_stack+42, NULL);
 _BUILD_p0h0(Data,vrr_stack+953, vrr_stack+1704, vrr_stack+1683, NULL, NULL, vrr_stack+793);
 _BUILD_00h0(Data,vrr_stack+1725, vrr_stack+823, vrr_stack+808, vrr_stack+92, vrr_stack+52, NULL);
 _BUILD_p0h0(Data,vrr_stack+149, vrr_stack+1725, vrr_stack+1704, NULL, NULL, vrr_stack+808);
 _BUILD_00h0(Data,vrr_stack+793, vrr_stack+838, vrr_stack+77, vrr_stack+102, vrr_stack+139, NULL);
 _BUILD_p0h0(Data,vrr_stack+2523, vrr_stack+1683, vrr_stack+793, NULL, NULL, vrr_stack+838);
 _BUILD_d0h0(Data,vrr_stack+2586, vrr_stack+953, vrr_stack+2523, vrr_stack+1704, vrr_stack+1683, vrr_stack+329);
 _BUILD_d0h0(Data,vrr_stack+2712, vrr_stack+149, vrr_stack+953, vrr_stack+1725, vrr_stack+1704, vrr_stack+239);
 _BUILD_f0h0(Data,vrr_stack+2838, vrr_stack+2712, vrr_stack+2586, vrr_stack+149, vrr_stack+953, vrr_stack+1353);
 _BUILD_00h0(Data,vrr_stack+1704, vrr_stack+62, vrr_stack+823, vrr_stack+123, vrr_stack+92, NULL);
 _BUILD_p0h0(Data,vrr_stack+1353, vrr_stack+1704, vrr_stack+1725, NULL, NULL, vrr_stack+823);
 _BUILD_d0h0(Data,vrr_stack+3048, vrr_stack+1353, vrr_stack+149, vrr_stack+1704, vrr_stack+1725, vrr_stack+284);
 _BUILD_f0h0(Data,vrr_stack+212, vrr_stack+3048, vrr_stack+2712, vrr_stack+1353, vrr_stack+149, vrr_stack+1443);
 _BUILD_00p0(Data,vrr_stack+149, Data->F+9, Data->F+10, NULL, NULL, NULL);
 _BUILD_00d0(Data,vrr_stack+152, vrr_stack+39, vrr_stack+149, Data->F+8, Data->F+9, NULL);
 _BUILD_00f0(Data,vrr_stack+158, vrr_stack+513, vrr_stack+152, vrr_stack+18, vrr_stack+39, NULL);
 _BUILD_00g0(Data,vrr_stack+168, vrr_stack+519, vrr_stack+158, vrr_stack+133, vrr_stack+513, NULL);
 _BUILD_00h0(Data,vrr_stack+183, vrr_stack+77, vrr_stack+168, vrr_stack+139, vrr_stack+519, NULL);
 _BUILD_p0h0(Data,vrr_stack+1353, vrr_stack+793, vrr_stack+183, NULL, NULL, vrr_stack+77);
 _BUILD_d0h0(Data,vrr_stack+0, vrr_stack+2523, vrr_stack+1353, vrr_stack+1683, vrr_stack+793, vrr_stack+529);
 _BUILD_f0h0(Data,vrr_stack+3174, vrr_stack+2586, vrr_stack+0, vrr_stack+953, vrr_stack+2523, vrr_stack+574);
 _BUILD_g0h0(Data,vrr_stack+3384, vrr_stack+2838, vrr_stack+3174, vrr_stack+2712, vrr_stack+2586, vrr_stack+1923);
 _BUILD_g0h0(Data,vrr_stack+3699, vrr_stack+212, vrr_stack+2838, vrr_stack+3048, vrr_stack+2712, vrr_stack+1533);
 _BUILD_h0f0(Data,vrr_stack+3048, vrr_stack+1203, vrr_stack+1053, vrr_stack+853, vrr_stack+693, vrr_stack+423);
   tmp = vrr_stack + 3048;
   target_ptr = Libint->vrr_classes[5][3];
   for(i=0;i<210;i++)
     target_ptr[i] += tmp[i];
 _BUILD_h0g0(Data,vrr_stack+2523, vrr_stack+2298, vrr_stack+2073, vrr_stack+1773, vrr_stack+1533, vrr_stack+1053);
   tmp = vrr_stack + 2523;
   target_ptr = Libint->vrr_classes[5][4];
   for(i=0;i<315;i++)
     target_ptr[i] += tmp[i];
 _BUILD_h0h0(Data,vrr_stack+422, vrr_stack+3699, vrr_stack+3384, vrr_stack+212, vrr_stack+2838, vrr_stack+2073);
   tmp = vrr_stack + 422;
   target_ptr = Libint->vrr_classes[5][5];
   for(i=0;i<441;i++)
     target_ptr[i] += tmp[i];

}

