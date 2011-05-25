  /* This machine-generated function computes a quartet of |pd) first derivative ERIs */

void d1hrr3_build_pd(const double *CD, double *vp, const double *I0, const double *I1,
        double c2, const double *I2, double c3, const double *I3, double c4, const double *I4,
        double c5, const double *I5, double c6, const double *I6, double c7, const double *I7, int ab_num)
{
  int ab;
  const double CD0 = CD[0];
  const double CD1 = CD[1];
  const double CD2 = CD[2];
  for(ab=0;ab<ab_num;ab++) {
    *(vp++) = I0[0] + CD0*I1[0] + c2*I2[0] - c5*I5[0];
    *(vp++) = I0[1] + CD0*I1[1] + c2*I2[1] - c5*I5[1];
    *(vp++) = I0[2] + CD0*I1[2] + c2*I2[2] - c5*I5[2];
    *(vp++) = I0[4] + CD1*I1[1] + c3*I3[1] - c6*I6[1];
    *(vp++) = I0[5] + CD1*I1[2] + c3*I3[2] - c6*I6[2];
    *(vp++) = I0[8] + CD2*I1[2] + c4*I4[2] - c7*I7[2];
    *(vp++) = I0[3] + CD0*I1[3] + c2*I2[3] - c5*I5[3];
    *(vp++) = I0[4] + CD0*I1[4] + c2*I2[4] - c5*I5[4];
    *(vp++) = I0[5] + CD0*I1[5] + c2*I2[5] - c5*I5[5];
    *(vp++) = I0[10] + CD1*I1[4] + c3*I3[4] - c6*I6[4];
    *(vp++) = I0[11] + CD1*I1[5] + c3*I3[5] - c6*I6[5];
    *(vp++) = I0[14] + CD2*I1[5] + c4*I4[5] - c7*I7[5];
    *(vp++) = I0[6] + CD0*I1[6] + c2*I2[6] - c5*I5[6];
    *(vp++) = I0[7] + CD0*I1[7] + c2*I2[7] - c5*I5[7];
    *(vp++) = I0[8] + CD0*I1[8] + c2*I2[8] - c5*I5[8];
    *(vp++) = I0[13] + CD1*I1[7] + c3*I3[7] - c6*I6[7];
    *(vp++) = I0[14] + CD1*I1[8] + c3*I3[8] - c6*I6[8];
    *(vp++) = I0[17] + CD2*I1[8] + c4*I4[8] - c7*I7[8];
    I0 += 18;
    I1 += 9;
    I2 += 9;
    I3 += 9;
    I4 += 9;
    I5 += 9;
    I6 += 9;
    I7 += 9;
  }
}
