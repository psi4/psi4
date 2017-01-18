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

/*!
** rotate_vecs(): Rotate a set of vectors around an arbitrary axis
**
** \brief Rotate a set of vectors around an arbitrary axis
** Vectors are rows of input matrix
**
** \param  w     double *  : axis to rotate around (wx, wy, wz) - gets normalized here
** \param  phi   double    : magnitude of rotation
** \param  v   double ** : points to rotate - column dim is 3; overwritten on exit
** \param  num_v  int       :
**
** Returns: none
**
** Rollin King, Feb. 2008
** \ingroup QT
*/

#include <stdio.h>
#include <math.h>
#include "psi4/libciomr/libciomr.h"

namespace psi {

void rotate_vecs(double *w, double phi, double **v, int num_v)
{
  int i, j;
  double **R, **v_new, wx, wy, wz, cp, norm;

  norm = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
  w[0] /= norm;
  w[1] /= norm;
  w[2] /= norm;

  wx = w[0]; wy = w[1]; wz = w[2];
  cp = 1.0 - cos(phi);

  R = block_matrix(3,3);

  R[0][0] =     cos(phi) + wx*wx*cp;
  R[0][1] = -wz*sin(phi) + wx*wy*cp;
  R[0][2] =  wy*sin(phi) + wx*wz*cp;
  R[1][0] =  wz*sin(phi) + wx*wy*cp;
  R[1][1] =     cos(phi) + wy*wy*cp;
  R[1][2] = -wx*sin(phi) + wy*wz*cp;
  R[2][0] = -wy*sin(phi) + wx*wz*cp;
  R[2][1] =  wx*sin(phi) + wy*wz*cp;
  R[2][2] =     cos(phi) + wz*wz*cp;

  v_new = block_matrix(num_v,3);
  mmult(R, 0, v, 1, v_new, 1, 3, 3, num_v, 0);

  for (i=0; i<num_v; ++i)
    for (j=0; j<3; ++j)
      v[i][j] = v_new[i][j];

  free_block(v_new);
  free_block(R);
}

}
