#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DX_k(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[0] - 7.000000*I1[0];
    *(vp++) = twotzeta*I0[1] - 6.000000*I1[1];
    *(vp++) = twotzeta*I0[2] - 6.000000*I1[2];
    *(vp++) = twotzeta*I0[3] - 5.000000*I1[3];
    *(vp++) = twotzeta*I0[4] - 5.000000*I1[4];
    *(vp++) = twotzeta*I0[5] - 5.000000*I1[5];
    *(vp++) = twotzeta*I0[6] - 4.000000*I1[6];
    *(vp++) = twotzeta*I0[7] - 4.000000*I1[7];
    *(vp++) = twotzeta*I0[8] - 4.000000*I1[8];
    *(vp++) = twotzeta*I0[9] - 4.000000*I1[9];
    *(vp++) = twotzeta*I0[10] - 3.000000*I1[10];
    *(vp++) = twotzeta*I0[11] - 3.000000*I1[11];
    *(vp++) = twotzeta*I0[12] - 3.000000*I1[12];
    *(vp++) = twotzeta*I0[13] - 3.000000*I1[13];
    *(vp++) = twotzeta*I0[14] - 3.000000*I1[14];
    *(vp++) = twotzeta*I0[15] - 2.000000*I1[15];
    *(vp++) = twotzeta*I0[16] - 2.000000*I1[16];
    *(vp++) = twotzeta*I0[17] - 2.000000*I1[17];
    *(vp++) = twotzeta*I0[18] - 2.000000*I1[18];
    *(vp++) = twotzeta*I0[19] - 2.000000*I1[19];
    *(vp++) = twotzeta*I0[20] - 2.000000*I1[20];
    *(vp++) = twotzeta*I0[21] - 1.000000*I1[21];
    *(vp++) = twotzeta*I0[22] - 1.000000*I1[22];
    *(vp++) = twotzeta*I0[23] - 1.000000*I1[23];
    *(vp++) = twotzeta*I0[24] - 1.000000*I1[24];
    *(vp++) = twotzeta*I0[25] - 1.000000*I1[25];
    *(vp++) = twotzeta*I0[26] - 1.000000*I1[26];
    *(vp++) = twotzeta*I0[27] - 1.000000*I1[27];
    *(vp++) = twotzeta*I0[28] ;
    *(vp++) = twotzeta*I0[29] ;
    *(vp++) = twotzeta*I0[30] ;
    *(vp++) = twotzeta*I0[31] ;
    *(vp++) = twotzeta*I0[32] ;
    *(vp++) = twotzeta*I0[33] ;
    *(vp++) = twotzeta*I0[34] ;
    *(vp++) = twotzeta*I0[35] ;
  I0 += 45;  I1 += 28;
  }
}
