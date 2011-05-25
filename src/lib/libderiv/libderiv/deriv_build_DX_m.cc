#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DX_m(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[0] - 9.000000*I1[0];
    *(vp++) = twotzeta*I0[1] - 8.000000*I1[1];
    *(vp++) = twotzeta*I0[2] - 8.000000*I1[2];
    *(vp++) = twotzeta*I0[3] - 7.000000*I1[3];
    *(vp++) = twotzeta*I0[4] - 7.000000*I1[4];
    *(vp++) = twotzeta*I0[5] - 7.000000*I1[5];
    *(vp++) = twotzeta*I0[6] - 6.000000*I1[6];
    *(vp++) = twotzeta*I0[7] - 6.000000*I1[7];
    *(vp++) = twotzeta*I0[8] - 6.000000*I1[8];
    *(vp++) = twotzeta*I0[9] - 6.000000*I1[9];
    *(vp++) = twotzeta*I0[10] - 5.000000*I1[10];
    *(vp++) = twotzeta*I0[11] - 5.000000*I1[11];
    *(vp++) = twotzeta*I0[12] - 5.000000*I1[12];
    *(vp++) = twotzeta*I0[13] - 5.000000*I1[13];
    *(vp++) = twotzeta*I0[14] - 5.000000*I1[14];
    *(vp++) = twotzeta*I0[15] - 4.000000*I1[15];
    *(vp++) = twotzeta*I0[16] - 4.000000*I1[16];
    *(vp++) = twotzeta*I0[17] - 4.000000*I1[17];
    *(vp++) = twotzeta*I0[18] - 4.000000*I1[18];
    *(vp++) = twotzeta*I0[19] - 4.000000*I1[19];
    *(vp++) = twotzeta*I0[20] - 4.000000*I1[20];
    *(vp++) = twotzeta*I0[21] - 3.000000*I1[21];
    *(vp++) = twotzeta*I0[22] - 3.000000*I1[22];
    *(vp++) = twotzeta*I0[23] - 3.000000*I1[23];
    *(vp++) = twotzeta*I0[24] - 3.000000*I1[24];
    *(vp++) = twotzeta*I0[25] - 3.000000*I1[25];
    *(vp++) = twotzeta*I0[26] - 3.000000*I1[26];
    *(vp++) = twotzeta*I0[27] - 3.000000*I1[27];
    *(vp++) = twotzeta*I0[28] - 2.000000*I1[28];
    *(vp++) = twotzeta*I0[29] - 2.000000*I1[29];
    *(vp++) = twotzeta*I0[30] - 2.000000*I1[30];
    *(vp++) = twotzeta*I0[31] - 2.000000*I1[31];
    *(vp++) = twotzeta*I0[32] - 2.000000*I1[32];
    *(vp++) = twotzeta*I0[33] - 2.000000*I1[33];
    *(vp++) = twotzeta*I0[34] - 2.000000*I1[34];
    *(vp++) = twotzeta*I0[35] - 2.000000*I1[35];
    *(vp++) = twotzeta*I0[36] - 1.000000*I1[36];
    *(vp++) = twotzeta*I0[37] - 1.000000*I1[37];
    *(vp++) = twotzeta*I0[38] - 1.000000*I1[38];
    *(vp++) = twotzeta*I0[39] - 1.000000*I1[39];
    *(vp++) = twotzeta*I0[40] - 1.000000*I1[40];
    *(vp++) = twotzeta*I0[41] - 1.000000*I1[41];
    *(vp++) = twotzeta*I0[42] - 1.000000*I1[42];
    *(vp++) = twotzeta*I0[43] - 1.000000*I1[43];
    *(vp++) = twotzeta*I0[44] - 1.000000*I1[44];
    *(vp++) = twotzeta*I0[45] ;
    *(vp++) = twotzeta*I0[46] ;
    *(vp++) = twotzeta*I0[47] ;
    *(vp++) = twotzeta*I0[48] ;
    *(vp++) = twotzeta*I0[49] ;
    *(vp++) = twotzeta*I0[50] ;
    *(vp++) = twotzeta*I0[51] ;
    *(vp++) = twotzeta*I0[52] ;
    *(vp++) = twotzeta*I0[53] ;
    *(vp++) = twotzeta*I0[54] ;
  I0 += 66;  I1 += 45;
  }
}
