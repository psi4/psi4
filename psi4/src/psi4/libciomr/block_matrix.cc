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
\file
\brief Allocate a blocked (memory-contiguous) 2D matrix of doubles
\ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <strings.h>
#include "psi4/psifiles.h"
#include <unistd.h>
#ifdef _POSIX_MEMLOCK
#include <sys/mman.h>
#endif

#include "psi4/psi4-dec.h"

namespace psi {

/*!
** block_matrix(): Allocate a 2D array of doubles using contiguous memory
**
** Allocates a contiguous block of memory for an array of
** doubles, allocates an array of pointers to the beginning of each row and
** returns the pointer to the first row pointer.  This allows transparent
** 2d-array style access, but keeps memory together such that the matrix
** could be used in conjunction with FORTRAN matrix routines.
**
** Allocates memory for an n x m matrix and returns a pointer to the
** first row.
**
** \param n = number of rows (unsigned long to allow large matrices)
** \param m = number of columns (unsigned long to allow large matrices)
** \param memlock = optional bool indicating whether to lock memory
**   into physical RAM or not, and available only where _POSIX_MEMLOCK
**   is defined. Defaults to false if not specified.
**
** Returns: double star pointer to newly allocated matrix
**
** T. Daniel Crawford
** Sometime in 1994
**
** Based on init_matrix() from libciomr
** \ingroup CIOMR
*/

double ** block_matrix(unsigned long int n, unsigned long int m, bool memlock)
{
    double **A=NULL;
    double *B=NULL;
    unsigned long int i;

    if(!m || !n) return(static_cast<double **>(0));

    A = new double*[n];
    if (A==NULL) {
        outfile->Printf("block_matrix: trouble allocating memory \n");
        outfile->Printf("n = %ld\n",n);
        exit(PSI_RETURN_FAILURE);
    }

    B = new double[n*m];
    if (B == NULL) {
        outfile->Printf("block_matrix: trouble allocating memory \n");
        outfile->Printf("m = %ld\n",m);
        exit(PSI_RETURN_FAILURE);
    }
    memset(static_cast<void*>(B), 0, m*n*sizeof(double));

    for (i = 0; i < n; i++) {
        A[i] = &(B[i*m]);
    }

#ifdef _POSIX_MEMLOCK
    if (memlock) {

        char* addr = (char*) B;
        unsigned long size = m*n*(unsigned long)sizeof(double);
        unsigned long page_offset, page_size;

        page_size = sysconf(_SC_PAGESIZE);
        page_offset = (unsigned long) addr % page_size;

        addr -= page_offset;  /* Adjust addr to page boundary */
        size += page_offset;  /* Adjust size with page_offset */

        if ( mlock(addr, size) ) {  /* Lock the memory */
            outfile->Printf("block_matrix: trouble locking memory \n");
            fflush(stderr);
            exit(PSI_RETURN_FAILURE);
        }

        addr = (char*) A;
        size = n*(unsigned long)sizeof(double*);

        page_offset = (unsigned long) addr % page_size;

        addr -= page_offset;  /* Adjust addr to page boundary */
        size += page_offset;  /* Adjust size with page_offset */

        if ( mlock(addr, size) ) {  /* Lock the memory */
            outfile->Printf("block_matrix: trouble locking memory \n");
            fflush(stderr);
            exit(PSI_RETURN_FAILURE);
        }
    }
#endif

    return(A);
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
void free_block(double **array)
{
    if(array == NULL) return;
    delete [] array[0];
    delete [] array;
}

}
