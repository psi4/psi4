#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DY_m(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[1] ;
    *(vp++) = twotzeta*I0[3] - 1.000000*I1[0];
    *(vp++) = twotzeta*I0[4] ;
    *(vp++) = twotzeta*I0[6] - 2.000000*I1[1];
    *(vp++) = twotzeta*I0[7] - 1.000000*I1[2];
    *(vp++) = twotzeta*I0[8] ;
    *(vp++) = twotzeta*I0[10] - 3.000000*I1[3];
    *(vp++) = twotzeta*I0[11] - 2.000000*I1[4];
    *(vp++) = twotzeta*I0[12] - 1.000000*I1[5];
    *(vp++) = twotzeta*I0[13] ;
    *(vp++) = twotzeta*I0[15] - 4.000000*I1[6];
    *(vp++) = twotzeta*I0[16] - 3.000000*I1[7];
    *(vp++) = twotzeta*I0[17] - 2.000000*I1[8];
    *(vp++) = twotzeta*I0[18] - 1.000000*I1[9];
    *(vp++) = twotzeta*I0[19] ;
    *(vp++) = twotzeta*I0[21] - 5.000000*I1[10];
    *(vp++) = twotzeta*I0[22] - 4.000000*I1[11];
    *(vp++) = twotzeta*I0[23] - 3.000000*I1[12];
    *(vp++) = twotzeta*I0[24] - 2.000000*I1[13];
    *(vp++) = twotzeta*I0[25] - 1.000000*I1[14];
    *(vp++) = twotzeta*I0[26] ;
    *(vp++) = twotzeta*I0[28] - 6.000000*I1[15];
    *(vp++) = twotzeta*I0[29] - 5.000000*I1[16];
    *(vp++) = twotzeta*I0[30] - 4.000000*I1[17];
    *(vp++) = twotzeta*I0[31] - 3.000000*I1[18];
    *(vp++) = twotzeta*I0[32] - 2.000000*I1[19];
    *(vp++) = twotzeta*I0[33] - 1.000000*I1[20];
    *(vp++) = twotzeta*I0[34] ;
    *(vp++) = twotzeta*I0[36] - 7.000000*I1[21];
    *(vp++) = twotzeta*I0[37] - 6.000000*I1[22];
    *(vp++) = twotzeta*I0[38] - 5.000000*I1[23];
    *(vp++) = twotzeta*I0[39] - 4.000000*I1[24];
    *(vp++) = twotzeta*I0[40] - 3.000000*I1[25];
    *(vp++) = twotzeta*I0[41] - 2.000000*I1[26];
    *(vp++) = twotzeta*I0[42] - 1.000000*I1[27];
    *(vp++) = twotzeta*I0[43] ;
    *(vp++) = twotzeta*I0[45] - 8.000000*I1[28];
    *(vp++) = twotzeta*I0[46] - 7.000000*I1[29];
    *(vp++) = twotzeta*I0[47] - 6.000000*I1[30];
    *(vp++) = twotzeta*I0[48] - 5.000000*I1[31];
    *(vp++) = twotzeta*I0[49] - 4.000000*I1[32];
    *(vp++) = twotzeta*I0[50] - 3.000000*I1[33];
    *(vp++) = twotzeta*I0[51] - 2.000000*I1[34];
    *(vp++) = twotzeta*I0[52] - 1.000000*I1[35];
    *(vp++) = twotzeta*I0[53] ;
    *(vp++) = twotzeta*I0[55] - 9.000000*I1[36];
    *(vp++) = twotzeta*I0[56] - 8.000000*I1[37];
    *(vp++) = twotzeta*I0[57] - 7.000000*I1[38];
    *(vp++) = twotzeta*I0[58] - 6.000000*I1[39];
    *(vp++) = twotzeta*I0[59] - 5.000000*I1[40];
    *(vp++) = twotzeta*I0[60] - 4.000000*I1[41];
    *(vp++) = twotzeta*I0[61] - 3.000000*I1[42];
    *(vp++) = twotzeta*I0[62] - 2.000000*I1[43];
    *(vp++) = twotzeta*I0[63] - 1.000000*I1[44];
    *(vp++) = twotzeta*I0[64] ;
  I0 += 66;  I1 += 45;
  }
}
