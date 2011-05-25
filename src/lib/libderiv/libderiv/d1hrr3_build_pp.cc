  /* This machine-generated function computes a quartet of |pp) first derivative ERIs */

void d1hrr3_build_pp(const double *CD, double *vp, const double *I0, const double *I1,
        double c2, const double *I2, double c3, const double *I3, double c4, const double *I4,
        double c5, const double *I5, double c6, const double *I6, double c7, const double *I7, int ab_num)
{
  int ab;
  const double CD0 = CD[0];
  const double CD1 = CD[1];
  const double CD2 = CD[2];
  for(ab=0;ab<ab_num;ab++) {
    *(vp++) = I0[0] + CD0*I1[0] + c2*I2[0] - c5*I5[0];
    *(vp++) = I0[1] + CD1*I1[0] + c3*I3[0] - c6*I6[0];
    *(vp++) = I0[2] + CD2*I1[0] + c4*I4[0] - c7*I7[0];
    *(vp++) = I0[1] + CD0*I1[1] + c2*I2[1] - c5*I5[1];
    *(vp++) = I0[3] + CD1*I1[1] + c3*I3[1] - c6*I6[1];
    *(vp++) = I0[4] + CD2*I1[1] + c4*I4[1] - c7*I7[1];
    *(vp++) = I0[2] + CD0*I1[2] + c2*I2[2] - c5*I5[2];
    *(vp++) = I0[4] + CD1*I1[2] + c3*I3[2] - c6*I6[2];
    *(vp++) = I0[5] + CD2*I1[2] + c4*I4[2] - c7*I7[2];
    I0 += 6;
    I1 += 3;
    I2 += 3;
    I3 += 3;
    I4 += 3;
    I5 += 3;
    I6 += 3;
    I7 += 3;
  }
}
