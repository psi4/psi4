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
** \file
** \brief Diagonalize a symmetric matrix in packed (lower triangular) form
** \ingroup CIOMR
*/

#include "libciomr.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

#include <cstdlib>

namespace psi {
/*!
** rsp(): diagonalize a symmetric matrix in packed (row-major lower triangular) form
** in 'array'. For square symmetric matrices, see sq_rsp().
**
** \param nm     = rows of matrix (unused, present for historic reasons)
** \param n      = columns of matrix
** \param nv     = number of elements in lower triangle (n*(n+1)/2)
** \param array  = matrix to diagonalize (packed as linear array)
** \param e_vals = array to hold eigenvalues
** \param matz   = 0 (no eigenvectors, eigenvals in ascending order)
**               = 1 (eigenvectors and eigenvalues in ascending order)
** \param e_vecs = matrix of eigenvectors (one column for each eigvector)
** \param toler  = tolerance for eigenvalues?  Often 1.0E-14. (unused, present for historic reasons)
**
** Returns: none
**
** \ingroup CIOMR
*/
void rsp(int /*nm*/, const int n, const int nv, const double * const array, double *e_vals, const int matz,
             double * const * const e_vecs, double /*toler*/){
    if ((matz > 1) || (matz < 0)){
        outfile->Printf("matz values other than 0 and 1 are no longer supported by rsp(...)");
        exit(PSI_RETURN_FAILURE);
    };
    // Do you want eigenvectors?
    bool eigenvectors = (matz == 1 || matz == 3) ? true : false;
    // Ascending or Descending?
    bool ascending = true;

    // LAPACK expects column-major-packed arrays, so allocate a temporary and rearrange.
    // This also takes care of users not expecting rsp to destroy the matrix they have passed in.
    double* tmp_array = init_array(nv);
    // To do the rearangement we run a column-major loop and calculate the packed row-major index (simpler)
    // If anyone ever needs it, here is the index calculation for packed column-major arrays that start at zero:
    // k = i-j + j*N - (j*(j-1))/2
    for(int j=0, ij=0; j<n; j++){
        for(int i=j; i<n; i++,ij++){
            const int f = j + (i*(i+1))/2;
            tmp_array[ij] = array[f];
        }
    }
    // LAPACK needs more memory for temporary values
    double* tmp_work = init_array(3*n);
    if (eigenvectors){
        // LAPACK needs a 1D array with N*N elements to put the eigenvectors in and cannot work with double**
        double* tmp_eigvecs = init_array(n*n);
        const auto info = C_DSPEV('V', 'L', n, tmp_array, e_vals, tmp_eigvecs, n, tmp_work);
        if (info != 0){
            outfile->Printf("DSPEV failed in a call of rsp(...), WITH eigenvectors");
            exit(PSI_RETURN_FAILURE);
        }
        // The eigenvectors are now in tmp_eigvecs in columns, flattened into 1D column major
        // Copy them to the 2D row major array as columns
        for(int j=0, ij=0; j<n; j++){
            for(int i=0; i<n; i++,ij++){
                e_vecs[i][j] = tmp_eigvecs[ij];
            }
        }
        free(tmp_eigvecs);
    } else {
        // LAPACK promises to not reference the eigenvector array if we request no eigenvectors, so don't allocate
        const auto info = C_DSPEV('N', 'L', n, tmp_array, e_vals, nullptr, n, tmp_work);
        if (info != 0){
            outfile->Printf("DSPEV failed in a call of rsp(...), WITHOUT eigenvectors");
            exit(PSI_RETURN_FAILURE);
        }
    }
    free(tmp_work);
    free(tmp_array);
}
}
