/*! \file 
    \ingroup QT
    \brief simple functions for 3-vectors
*/
#include <stdio.h>
#include <stdlib.h>
#include <libqt/qt.h>
#include <math.h>

namespace psi {

double dot_prod(double *v1, double *v2) {
  return v1[0]*v2[0]+ v1[1]*v2[1]+ v1[2]*v2[2];
}

void cross_prod(double *v1, double *v2, double *out) {
  out[0] = v1[1]*v2[2]-v1[2]*v2[1];
  out[1] = -v1[0]*v2[2]+v1[2]*v2[0];
  out[2] = v1[0]*v2[1]-v1[1]*v2[0];
  return;
}

void unit_vec(double *B, double *A, double *AB) {
  double norm = 0.0;
  int i;
  
  for (i=0; i<3; i++)
    norm += (A[i]-B[i])*(A[i]-B[i]);
  norm = sqrt(norm);
  for (i=0; i<3; i++)
    AB[i] = (B[i] - A[i]) / norm;
  return;
}

}

