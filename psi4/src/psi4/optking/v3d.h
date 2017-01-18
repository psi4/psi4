/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file v3d.h
    \ingroup optking
    \brief Functions for real-space vectors of length 3
*/

#ifndef _opt_v3d_h_
#define _opt_v3d_h_

#include "package.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

#define V3D_SQR(x) ((x)*(x))

#define PARALLEL_LIMIT (1.0e-10)

namespace opt { namespace v3d {

// scalar multiply a vector
inline void v3d_scm(const double a, double *A) {
  A[0] *= a; A[1] *= a; A[2] *= a;
}

// compute norm of a vector
inline double v3d_norm(const double *A) {
  return sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
}

// take dot product of two vectors
inline double v3d_dot(const double *A, const double *B) {
  return (A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
}

// returns distance between 2 points
inline double v3d_dist(const double *A, const double *B) {
  return sqrt(V3D_SQR(B[0]-A[0]) + V3D_SQR(B[1]-A[1]) + V3D_SQR(B[2]-A[2]));
}

// Compute a*X + Y -> Z
inline void v3d_axpy(const double a, const double *X, const double *Y, double *Z) {
  Z[0] = a*X[0] + Y[0];
  Z[1] = a*X[1] + Y[1];
  Z[2] = a*X[2] + Y[2];
}

// normalize vector.  Return "false and leave A unchanged if norm of A is
// less than min_norm or greater than max_norm
inline bool v3d_normalize(double *A, const double min_norm=1.0e-8, const double max_norm=1.0e8) {
  double tval = v3d_norm(A);
  if ( tval < min_norm || tval > max_norm)
    return false;
  else
    v3d_scm(1.0/tval, A);
  return true;
}

// Compute cross product of two vectors
inline void v3d_cross_product(const double *u, const double *v, double *X) {
  X[0] = u[1]*v[2]-u[2]*v[1];
  X[1] = -1.0*(u[0]*v[2]-u[2]*v[0]);
  X[2] = u[0]*v[1]-u[1]*v[0];
  return;
}

// check to see if two vectors are parallel, i.e., is |A.B| close to 1
inline bool v3d_is_parallel(const double *A, const double *B) {
  if ( fabs( fabs(v3d_dot(A,B)) - 1.0e0) > PARALLEL_LIMIT )
    return false;
  else
    return true;
}

// Compute vector A->B.  Normalize eAB.  Return "false" and do not normalize
// eAB if points are too close or distant
inline bool v3d_eAB(const double *A, const double *B, double *eAB,
const double min_norm=1.0e-8, const double max_norm=1.0e15) {
  v3d_axpy(-1, A, B, eAB);
  return ( v3d_normalize(eAB, min_norm, max_norm) );
}

// Computed angle in radians A-B-C (between vector B->A and vector B->C)
// if points are absurdly close or far apart, returns false
// tol is nearness of cos to 1/-1 to make angle 0/pi
bool v3d_angle(const double *A, const double *B, const double *C, double & phi, double tol=1e-14);

// Computed torsional angle in radians A-B-C-D
// Returns false if bond angles ABC or BCD are too close to 0 or 180
bool v3d_tors(const double *A, const double *B, const double *C, const double *D, double & tau);

// Computed out-of-plane angle in radians A-B-C-D for angle between A-B and C-B-D planes.
// Returns false if bond angles BCD is too close to 0 or 180
bool v3d_oofp(const double *A, const double *B, const double *C, const double *D,
  double & oop_angle);

}}

#endif