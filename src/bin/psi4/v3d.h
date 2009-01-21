/*! \file misc-3d.cc
    \ingroup OPT08
    \brief Functions for manipulating 3d space vectors
*/

#ifndef _psi4_src_bin_opt09_v3d_h_
#define _psi4_src_bin_opt09_v3d_h_

#include <cmath>
#define V3D_SQR(x) ((x)*(x))

namespace psi { namespace v3d {

// take dot product of 3d vectors
inline double v3d_dot(const double *A, const double *B) {
  return (A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
}

// multiply a 3d vector by a scalar
inline void v3d_scm(const double a, double *A) {
  A[0] *= a; A[1] *= a; A[2] *= a;
}

// returns distance between 2 points
inline double v3d_dist(const double *A, const double *B) {
  return sqrt(V3D_SQR(B[0]-A[0])+V3D_SQR(B[1]-A[1])+V3D_SQR(B[2]-A[2]));
}

// returns norm of a 3d vector
inline double v3d_norm(const double *A) {
  return sqrt(V3D_SQR(A[0])+V3D_SQR(A[1])+V3D_SQR(A[2]));
}

// normalize a 3d vector
inline void v3d_normalize(double *A) {
  double tval = v3d_norm(A);
  v3d_scm(1.0/tval,A);
}

// Compute a*X + Y -> Z
void v3d_axpy(const double a, const double *X, const double *Y, double *Z);

// compute angle A-B-C, i.e., between vector B->A and B->C
double v3d_angle(const double *A, const double *B, const double *C);

}}

#endif

