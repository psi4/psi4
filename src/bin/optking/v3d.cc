/*! \file v3d.cc
    \ingroup OPT08
    \brief Functions for manipulating 3d space vectors
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <psifiles.h>
#include <physconst.h>

#include "v3d.h"

namespace psi { namespace v3d {

extern "C" { extern FILE *outfile;}

// returns distance between 2 points
double v3d_dist(const double *A, const double *B) {
  double tval;
  tval = sqr(B[0]-A[0]) + sqr(B[1]-A[1]) + sqr(B[2]-A[2]);
  if (tval < 0.0) {
    fprintf(outfile,"Warning: distance in v3d_dist() is imaginary!\n");
    return 0.0;
  }
  else
    return sqrt(tval);
}

// returns norm of a 3d vector
double v3d_norm(const double *A) {
  double tval; 
  tval = sqr(A[0]) + sqr(A[1]) + sqr(A[2]);
  if (tval < 0.0) {
    fprintf(stderr,"Warning: distance in v3d_norm() is imaginary!\n");
    return 0.0;
  }
  else
    return sqrt(tval);
}

// normalize a 3d vector
void v3d_normalize(double *A) {
  double tval = v3d_norm(A);
  if (tval < 0) {
    fprintf(stderr,"Warning: could not normalize vector in v3d_normalize()!\n");
    std::exit(PSI_RETURN_FAILURE);
  }
  v3d_scm(1.0/tval,A);
}

// return angle A-B-C (between vector B->A and vector B->C)
double v3d_angle(const double *A, const double *B, const double *C) {
  int i;
  double eBA[3],eBC[3],dotprod;

  for (i=0;i<3;++i) {
    eBA[i] = A[i] - B[i];
    eBC[i] = C[i] - B[i];
  }

  v3d_normalize(eBA);
  v3d_normalize(eBC);

  dotprod = v3d_dot(eBA,eBC);

  if (dotprod > 1.0) return 0.0;
  else if (dotprod < -1.0) return 180.0;
  else return (acos(dotprod)*180.0/_pi);
}

}}

