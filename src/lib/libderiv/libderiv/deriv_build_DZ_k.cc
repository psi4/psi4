#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DZ_k(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[2] ;
    *(vp++) = twotzeta*I0[4] ;
    *(vp++) = twotzeta*I0[5] - 1.000000*I1[0];
    *(vp++) = twotzeta*I0[7] ;
    *(vp++) = twotzeta*I0[8] - 1.000000*I1[1];
    *(vp++) = twotzeta*I0[9] - 2.000000*I1[2];
    *(vp++) = twotzeta*I0[11] ;
    *(vp++) = twotzeta*I0[12] - 1.000000*I1[3];
    *(vp++) = twotzeta*I0[13] - 2.000000*I1[4];
    *(vp++) = twotzeta*I0[14] - 3.000000*I1[5];
    *(vp++) = twotzeta*I0[16] ;
    *(vp++) = twotzeta*I0[17] - 1.000000*I1[6];
    *(vp++) = twotzeta*I0[18] - 2.000000*I1[7];
    *(vp++) = twotzeta*I0[19] - 3.000000*I1[8];
    *(vp++) = twotzeta*I0[20] - 4.000000*I1[9];
    *(vp++) = twotzeta*I0[22] ;
    *(vp++) = twotzeta*I0[23] - 1.000000*I1[10];
    *(vp++) = twotzeta*I0[24] - 2.000000*I1[11];
    *(vp++) = twotzeta*I0[25] - 3.000000*I1[12];
    *(vp++) = twotzeta*I0[26] - 4.000000*I1[13];
    *(vp++) = twotzeta*I0[27] - 5.000000*I1[14];
    *(vp++) = twotzeta*I0[29] ;
    *(vp++) = twotzeta*I0[30] - 1.000000*I1[15];
    *(vp++) = twotzeta*I0[31] - 2.000000*I1[16];
    *(vp++) = twotzeta*I0[32] - 3.000000*I1[17];
    *(vp++) = twotzeta*I0[33] - 4.000000*I1[18];
    *(vp++) = twotzeta*I0[34] - 5.000000*I1[19];
    *(vp++) = twotzeta*I0[35] - 6.000000*I1[20];
    *(vp++) = twotzeta*I0[37] ;
    *(vp++) = twotzeta*I0[38] - 1.000000*I1[21];
    *(vp++) = twotzeta*I0[39] - 2.000000*I1[22];
    *(vp++) = twotzeta*I0[40] - 3.000000*I1[23];
    *(vp++) = twotzeta*I0[41] - 4.000000*I1[24];
    *(vp++) = twotzeta*I0[42] - 5.000000*I1[25];
    *(vp++) = twotzeta*I0[43] - 6.000000*I1[26];
    *(vp++) = twotzeta*I0[44] - 7.000000*I1[27];
  I0 += 45;  I1 += 28;
  }
}
