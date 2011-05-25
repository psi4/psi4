#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DX_l(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[0] - 8.000000*I1[0];
    *(vp++) = twotzeta*I0[1] - 7.000000*I1[1];
    *(vp++) = twotzeta*I0[2] - 7.000000*I1[2];
    *(vp++) = twotzeta*I0[3] - 6.000000*I1[3];
    *(vp++) = twotzeta*I0[4] - 6.000000*I1[4];
    *(vp++) = twotzeta*I0[5] - 6.000000*I1[5];
    *(vp++) = twotzeta*I0[6] - 5.000000*I1[6];
    *(vp++) = twotzeta*I0[7] - 5.000000*I1[7];
    *(vp++) = twotzeta*I0[8] - 5.000000*I1[8];
    *(vp++) = twotzeta*I0[9] - 5.000000*I1[9];
    *(vp++) = twotzeta*I0[10] - 4.000000*I1[10];
    *(vp++) = twotzeta*I0[11] - 4.000000*I1[11];
    *(vp++) = twotzeta*I0[12] - 4.000000*I1[12];
    *(vp++) = twotzeta*I0[13] - 4.000000*I1[13];
    *(vp++) = twotzeta*I0[14] - 4.000000*I1[14];
    *(vp++) = twotzeta*I0[15] - 3.000000*I1[15];
    *(vp++) = twotzeta*I0[16] - 3.000000*I1[16];
    *(vp++) = twotzeta*I0[17] - 3.000000*I1[17];
    *(vp++) = twotzeta*I0[18] - 3.000000*I1[18];
    *(vp++) = twotzeta*I0[19] - 3.000000*I1[19];
    *(vp++) = twotzeta*I0[20] - 3.000000*I1[20];
    *(vp++) = twotzeta*I0[21] - 2.000000*I1[21];
    *(vp++) = twotzeta*I0[22] - 2.000000*I1[22];
    *(vp++) = twotzeta*I0[23] - 2.000000*I1[23];
    *(vp++) = twotzeta*I0[24] - 2.000000*I1[24];
    *(vp++) = twotzeta*I0[25] - 2.000000*I1[25];
    *(vp++) = twotzeta*I0[26] - 2.000000*I1[26];
    *(vp++) = twotzeta*I0[27] - 2.000000*I1[27];
    *(vp++) = twotzeta*I0[28] - 1.000000*I1[28];
    *(vp++) = twotzeta*I0[29] - 1.000000*I1[29];
    *(vp++) = twotzeta*I0[30] - 1.000000*I1[30];
    *(vp++) = twotzeta*I0[31] - 1.000000*I1[31];
    *(vp++) = twotzeta*I0[32] - 1.000000*I1[32];
    *(vp++) = twotzeta*I0[33] - 1.000000*I1[33];
    *(vp++) = twotzeta*I0[34] - 1.000000*I1[34];
    *(vp++) = twotzeta*I0[35] - 1.000000*I1[35];
    *(vp++) = twotzeta*I0[36] ;
    *(vp++) = twotzeta*I0[37] ;
    *(vp++) = twotzeta*I0[38] ;
    *(vp++) = twotzeta*I0[39] ;
    *(vp++) = twotzeta*I0[40] ;
    *(vp++) = twotzeta*I0[41] ;
    *(vp++) = twotzeta*I0[42] ;
    *(vp++) = twotzeta*I0[43] ;
    *(vp++) = twotzeta*I0[44] ;
  I0 += 55;  I1 += 36;
  }
}
