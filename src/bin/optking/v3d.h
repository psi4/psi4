/*! \file misc-3d.cc
    \ingroup OPT08
    \brief Functions for manipulating 3d space vectors
*/

#ifndef _psi3_src_bin_optking_v3d_h_
#define _psi3_src_bin_optking_v3d_h_

namespace psi { namespace v3d {

inline double sqr(const double val) { return val * val; }

// take dot product of 3d vectors
inline double v3d_dot(const double *A, const double *B) {
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

// multiply a 3d vector by a scalar
inline void v3d_scm(const double a, double *A) {
  A[0] *= a; A[1] *= a; A[2] *= a;
}

// returns distance between 2 points
double v3d_dist(const double *A, const double *B);

// returns norm of a 3d vector
double v3d_norm(const double *A);

// normalize a 3d vector
void v3d_normalize(double *A);

// compute angle A-B-C, i.e., between vector B->A and B->C
double v3d_angle(const double *A, const double *B, const double *C);

}}

#endif

