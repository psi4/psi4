/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
\brief Allocate a blocked (memory-contiguous) 2D matrix of doubles
\ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#ifdef _POSIX_MEMLOCK
#include <sys/mman.h>
#endif

#include "psi4/psi4-dec.h"

namespace psi {

/*!
** \brief block_matrix() allocates a 2D array of doubles using contiguous memory.
**
** \details Allocates a contiguous block of memory for an array of
** doubles, allocates an array of pointers to the beginning of each row and
** returns the pointer to the first row pointer.  This allows transparent
** 2d-array style access, but keeps memory together such that the matrix
** could be used in conjunction with FORTRAN matrix routines.
**
** Allocates memory for an n x m matrix and returns a pointer to the
** first row.
**
** \param n = number of rows (size_t to allow large matrices)
** \param m = number of columns (size_t to allow large matrices)
** \param memlock = optional bool indicating whether to lock memory
**   into physical RAM or not, and available only where _POSIX_MEMLOCK
**   is defined. Defaults to false if not specified.
**
** \return double star pointer to newly allocated matrix
**
** \author T. Daniel Crawford
** \date Sometime in 1994
**
** \remark Based on init_matrix() from libciomr
** \ingroup CIOMR
*/
PSI_API [[nodiscard]] double **block_matrix(size_t n, size_t m, bool memlock) {
    double **A = nullptr;
    double *B = nullptr;
    size_t i;

    if (!m || !n) return (static_cast<double **>(nullptr));

    A = new double *[n];
    if (A == nullptr) {
        outfile->Printf("block_matrix: trouble allocating memory \n");
        outfile->Printf("n = %ld\n", n);
        throw PSIEXCEPTION("Could not allocate memory in block_matrix! Tried to allocate " +
                           std::to_string(n * sizeof(double *)) + " bytes.");
    }

    B = new double[n * m];
    if (B == nullptr) {
        outfile->Printf("block_matrix: trouble allocating memory \n");
        outfile->Printf("m = %ld\n", m);
        throw PSIEXCEPTION("Could not allocate memory in block_matrix! Tried to allocate " +
                           std::to_string(n * m * sizeof(double)) + " bytes.");
    }
    memset(static_cast<void *>(B), 0, m * n * sizeof(double));

    for (i = 0; i < n; i++) {
        A[i] = &(B[i * m]);
    }

#ifdef _POSIX_MEMLOCK
    if (memlock) {
        char *addr = (char *)B;
        size_t size = m * n * (size_t)sizeof(double);
        size_t page_offset, page_size;

        page_size = sysconf(_SC_PAGESIZE);
        page_offset = (size_t)addr % page_size;

        addr -= page_offset; /* Adjust addr to page boundary */
        size += page_offset; /* Adjust size with page_offset */

        if (mlock(addr, size)) { /* Lock the memory */
            outfile->Printf("block_matrix: trouble locking memory \n");
            fflush(stderr);
            exit(PSI_RETURN_FAILURE);
        }

        addr = (char *)A;
        size = n * (size_t)sizeof(double *);

        page_offset = (size_t)addr % page_size;

        addr -= page_offset; /* Adjust addr to page boundary */
        size += page_offset; /* Adjust size with page_offset */

        if (mlock(addr, size)) { /* Lock the memory */
            outfile->Printf("block_matrix: trouble locking memory \n");
            fflush(stderr);
            exit(PSI_RETURN_FAILURE);
        }
    }
#endif

    return (A);
}

/*!
** free_block(): Free a block matrix
**
** \param array = pointer to matrix to be freed
**
** Returns: none
**
** \ingroup CIOMR
*/
void PSI_API free_block(double **array) {
    if (array == nullptr) return;
    delete[] array[0];
    delete[] array;
}
}  // namespace psi
