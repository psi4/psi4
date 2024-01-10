/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
** \brief This file includes the long integer versions of several psi routines
** for handling arrays and matrices of doubles
**
** David Sherrill, 1996
**
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

/*!
** init_long_int_array(): Allocates memory for one-D array of long ints of
** dimension  'size' and returns pointer to 1st element.  Zeroes all elements.
**
** Just modified the init_int_array() routine to do long int's instead.
**
** Returns: pointer to new array
**
** C. David Sherrill
** \ingroup CIOMR
*/
long int *init_long_int_array(int size) {
    long int *array;

    if ((array = (long int *)malloc(sizeof(long int) * size)) == nullptr) {
        outfile->Printf("init_array:  trouble allocating memory \n");
        outfile->Printf("size = %d\n", size);
        exit(PSI_RETURN_FAILURE);
    }
    memset(array, 0, sizeof(long int) * size);
    return (array);
}

/*!
** init_size_t_array(): Allocates memory for one-D array of size_t of
**  dimension
** 'size' and returns pointer to 1st element.  Zeroes all elements.
**
** \param size = length of array to allocate
**
** Returns: pointer to new array
**
** Added by RAK, 2020
** \ingroup CIOMR
*/
PSI_API size_t *init_size_t_array(int size) {
    size_t *array;

    if ((array = (size_t *)malloc(sizeof(size_t) * size)) == nullptr) {
        outfile->Printf("init_size_t_array:  trouble allocating memory \n");
        outfile->Printf("size = %d\n", size);
        exit(PSI_RETURN_FAILURE);
    }
    memset(array, 0, sizeof(size_t) * size);
    return (array);
}
}
