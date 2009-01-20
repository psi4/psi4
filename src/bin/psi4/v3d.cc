/*! \file v3d.cc
    \ingroup OPT08
    \brief Functions for manipulating 3d space vectors
*/

#include "v3d.h"

namespace psi { namespace v3d {

// Compute a*X + Y -> Z
void v3d_axpy(const double a, const double *X, const double *Y, double *Z) {
  Z[0] = a*X[0] + Y[0];
  Z[1] = a*X[1] + Y[1];
  Z[2] = a*X[2] + Y[2];
}

// return angle A-B-C (between vector B->A and vector B->C)
double v3d_angle(const double *A, const double *B, const double *C) {
  int i;
  double eBA[3],eBC[3],dotprod;

  v3d_axpy(-1,B,A,eBA);
  v3d_axpy(-1,B,C,eBC);

  v3d_normalize(eBA);
  v3d_normalize(eBC);

  dotprod = v3d_dot(eBA,eBC);

  if (dotprod >= 1.0) return 0.0;
  else if (dotprod <= -1.0) return 180.0;
  else return (acos(dotprod)*180.0/acos(-1.0));
}

}}

