#include <libint/libint.h>
#include "libderiv.h"

void deriv_build_DX_h(prim_data *Data, const int abc_num, double *vp, const double *I0, const double *I1)
{
  const double twotzeta = Data->twozeta_d;
  const double *i0, *i1;
  int abc;

  for(abc=0;abc<abc_num;abc++) {
    *(vp++) = twotzeta*I0[0] - 5.000000*I1[0];
    *(vp++) = twotzeta*I0[1] - 4.000000*I1[1];
    *(vp++) = twotzeta*I0[2] - 4.000000*I1[2];
    *(vp++) = twotzeta*I0[3] - 3.000000*I1[3];
    *(vp++) = twotzeta*I0[4] - 3.000000*I1[4];
    *(vp++) = twotzeta*I0[5] - 3.000000*I1[5];
    *(vp++) = twotzeta*I0[6] - 2.000000*I1[6];
    *(vp++) = twotzeta*I0[7] - 2.000000*I1[7];
    *(vp++) = twotzeta*I0[8] - 2.000000*I1[8];
    *(vp++) = twotzeta*I0[9] - 2.000000*I1[9];
    *(vp++) = twotzeta*I0[10] - 1.000000*I1[10];
    *(vp++) = twotzeta*I0[11] - 1.000000*I1[11];
    *(vp++) = twotzeta*I0[12] - 1.000000*I1[12];
    *(vp++) = twotzeta*I0[13] - 1.000000*I1[13];
    *(vp++) = twotzeta*I0[14] - 1.000000*I1[14];
    *(vp++) = twotzeta*I0[15] ;
    *(vp++) = twotzeta*I0[16] ;
    *(vp++) = twotzeta*I0[17] ;
    *(vp++) = twotzeta*I0[18] ;
    *(vp++) = twotzeta*I0[19] ;
    *(vp++) = twotzeta*I0[20] ;
  I0 += 28;  I1 += 15;
  }
}
