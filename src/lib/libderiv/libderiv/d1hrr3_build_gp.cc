  /* This machine-generated function computes a quartet of |gp) first derivative ERIs */

void d1hrr3_build_gp(const double *CD, double *vp, const double *I0, const double *I1,
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
    *(vp++) = I0[3] + CD0*I1[3] + c2*I2[3] - c5*I5[3];
    *(vp++) = I0[6] + CD1*I1[3] + c3*I3[3] - c6*I6[3];
    *(vp++) = I0[7] + CD2*I1[3] + c4*I4[3] - c7*I7[3];
    *(vp++) = I0[4] + CD0*I1[4] + c2*I2[4] - c5*I5[4];
    *(vp++) = I0[7] + CD1*I1[4] + c3*I3[4] - c6*I6[4];
    *(vp++) = I0[8] + CD2*I1[4] + c4*I4[4] - c7*I7[4];
    *(vp++) = I0[5] + CD0*I1[5] + c2*I2[5] - c5*I5[5];
    *(vp++) = I0[8] + CD1*I1[5] + c3*I3[5] - c6*I6[5];
    *(vp++) = I0[9] + CD2*I1[5] + c4*I4[5] - c7*I7[5];
    *(vp++) = I0[6] + CD0*I1[6] + c2*I2[6] - c5*I5[6];
    *(vp++) = I0[10] + CD1*I1[6] + c3*I3[6] - c6*I6[6];
    *(vp++) = I0[11] + CD2*I1[6] + c4*I4[6] - c7*I7[6];
    *(vp++) = I0[7] + CD0*I1[7] + c2*I2[7] - c5*I5[7];
    *(vp++) = I0[11] + CD1*I1[7] + c3*I3[7] - c6*I6[7];
    *(vp++) = I0[12] + CD2*I1[7] + c4*I4[7] - c7*I7[7];
    *(vp++) = I0[8] + CD0*I1[8] + c2*I2[8] - c5*I5[8];
    *(vp++) = I0[12] + CD1*I1[8] + c3*I3[8] - c6*I6[8];
    *(vp++) = I0[13] + CD2*I1[8] + c4*I4[8] - c7*I7[8];
    *(vp++) = I0[9] + CD0*I1[9] + c2*I2[9] - c5*I5[9];
    *(vp++) = I0[13] + CD1*I1[9] + c3*I3[9] - c6*I6[9];
    *(vp++) = I0[14] + CD2*I1[9] + c4*I4[9] - c7*I7[9];
    *(vp++) = I0[10] + CD0*I1[10] + c2*I2[10] - c5*I5[10];
    *(vp++) = I0[15] + CD1*I1[10] + c3*I3[10] - c6*I6[10];
    *(vp++) = I0[16] + CD2*I1[10] + c4*I4[10] - c7*I7[10];
    *(vp++) = I0[11] + CD0*I1[11] + c2*I2[11] - c5*I5[11];
    *(vp++) = I0[16] + CD1*I1[11] + c3*I3[11] - c6*I6[11];
    *(vp++) = I0[17] + CD2*I1[11] + c4*I4[11] - c7*I7[11];
    *(vp++) = I0[12] + CD0*I1[12] + c2*I2[12] - c5*I5[12];
    *(vp++) = I0[17] + CD1*I1[12] + c3*I3[12] - c6*I6[12];
    *(vp++) = I0[18] + CD2*I1[12] + c4*I4[12] - c7*I7[12];
    *(vp++) = I0[13] + CD0*I1[13] + c2*I2[13] - c5*I5[13];
    *(vp++) = I0[18] + CD1*I1[13] + c3*I3[13] - c6*I6[13];
    *(vp++) = I0[19] + CD2*I1[13] + c4*I4[13] - c7*I7[13];
    *(vp++) = I0[14] + CD0*I1[14] + c2*I2[14] - c5*I5[14];
    *(vp++) = I0[19] + CD1*I1[14] + c3*I3[14] - c6*I6[14];
    *(vp++) = I0[20] + CD2*I1[14] + c4*I4[14] - c7*I7[14];
    I0 += 21;
    I1 += 15;
    I2 += 15;
    I3 += 15;
    I4 += 15;
    I5 += 15;
    I6 += 15;
    I7 += 15;
  }
}
