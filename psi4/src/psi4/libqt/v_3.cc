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

/*! \file
    \ingroup QT
    \brief simple functions for 3-vectors
*/
#include <stdio.h>
#include <stdlib.h>
#include "psi4/libqt/qt.h"
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
