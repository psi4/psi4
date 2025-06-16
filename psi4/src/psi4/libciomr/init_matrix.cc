/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
** \file
** \brief Initialize a matrix of doubles
** \ingroup CIOMR
*/

#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace psi {

/**
 *  WARNING: Psi 3 init/free_matrix routines deprecated
 *  by Robert Parrish, robparrish@gmail.com
 *
 *  block_matrix() replaces this routine
 *
 *  the signature of this method remains the same
 *
 *  June 22, 2010
 **/
/*!
** init_matrix(): Initialize an nxm matrix of doubles and return a pointer to
** the first row.  Note that this does not form a matrix which is
** necessarily contiguous in memory.  Use block_matrix() for that.
**
** \param n = number of rows (size_t to allow large matrices)
** \param m = number of columns (size_t to allow large matrices)
**
** Returns: pointer to first row
**
** \ingroup CIOMR
*/
double **init_matrix(size_t n, size_t m) {
    double **A = nullptr;
    double *B = nullptr;
    size_t i;

    if (!m || !n) return (static_cast<double **>(nullptr));

    //  if ((A = (double **) malloc(n * (size_t)sizeof(double *)))==nullptr) {
    if ((A = new double *[n]) == nullptr) {
        std::ostringstream oss;
        oss << "block_matrix: trouble allocating memory, n = " << n << "\n";
        throw std::runtime_error(oss.str());
    }

    //  if ((B = (double *) malloc(m*n * (size_t)sizeof(double)))==nullptr) {
    if ((B = new double[n * m]) == nullptr) {
        std::ostringstream oss;
        oss << "block_matrix: trouble allocating memory, n = " << n << " m = " << m << "\n";
        throw std::runtime_error(oss.str());
    }

    // bzero is not in the C standard, use memset instead.
    // bzero(B, m*n*(size_t)sizeof(double));
    memset(static_cast<void *>(B), 0, m * n * sizeof(double));

    for (i = 0; i < n; i++) {
        A[i] = &(B[i * m]);
    }

    return (A);
}

/**
 *  WARNING: Psi 3 init/free_matrix routines deprecated
 *  by Robert Parrish, robparrish@gmail.com
 *
 *  use block_matrix allocation/free calls instead
 *
 *  the signature of this method remains the same
 *
 *  June 22, 2010
 **/
/*!
** free_matrix(): Free a 2D matrix allocated with init_matrix().
**
** \param array = matrix to free
**
** Returns: none
**
** \ingroup CIOMR
*/
void free_matrix(double **array, size_t /*size*/) {
    if (array == nullptr) return;
    delete[] array[0];
    delete[] array;
    // <<<<<<<<<<<<<<<<<<<<<
    // BEGIN DEPRECATED CODE
    // <<<<<<<<<<<<<<<<<<<<<

    /**
    size_t i;

    for (i=0; i < size ; i++) {
      free(array[i]);
    }

    free(array);
  **/
}
}  // namespace psi
