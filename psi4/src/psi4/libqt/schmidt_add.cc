/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
  \file
  \brief Gram-Schmidt orthogonalize vector and add to set
  \ingroup QT
*/

#include <cstdio>
#include <cmath>

#include "psi4/pragma.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

namespace psi {

#define NORM_TOL 1.0E-5

/*!
** SCHMIDT_ADD(): Assume A is a orthogonal matrix.  This function Gram-Schmidt
** orthogonalizes a new vector v and adds it to matrix A.  A must contain
** a free row pointer for a new row.  Don't add orthogonalized v' if
** norm(v') < NORM_TOL.
**
** David Sherrill, Feb 1994
**
** \param A    = matrix to add new vector to
** \param rows = current number of rows in A
**               (A must have ptr for 'rows+1' row.)
** \param cols = columns in A
** \parm v     = vector to add to A after it has been made orthogonal
**               to rest of A
**
** Returns: 1 if a vector is added to A, 0 otherwise
** \ingroup QT
*/
PSI_API int schmidt_add(double **A, int rows, int cols, double *v) {
    double dotval, normval;
    int i, I;

    for (i = 0; i < rows; i++) {
        // dot_arr(A[i], v, cols, &dotval) ;
        dotval = C_DDOT(cols, A[i], 1, v, 1);
        for (I = 0; I < cols; I++) v[I] -= dotval * A[i][I];
    }

    // dot_arr(v, v, cols, &normval) ;
    normval = C_DDOT(cols, v, 1, v, 1);
    normval = sqrt(normval);

    if (normval < NORM_TOL)
        return (0);
    else {
        if (A[rows] == nullptr) A[rows] = init_array(cols);
        for (I = 0; I < cols; I++) A[rows][I] = v[I] / normval;
        return (1);
    }
}
}
