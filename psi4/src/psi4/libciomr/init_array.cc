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

/*
** \file
** \brief Initialize an array of doubles
** \ingroup CIOMR
*/

#include "psi4/psifiles.h"
#include <cstdio>
#include <cstdlib>
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
namespace psi {

/*!
** init_array(): This function initializes an array of doubles of
** length 'size' and returns a pointer to the first element
**
** \param size = length of array (size_t to allow large arrays)
**
** Returns: pointer to new array
**
** \ingroup CIOMR
*/
double *init_array(size_t size) {
    double *array;

    if ((array = (double *)malloc(size * (size_t)sizeof(double))) == nullptr) {
        std::ostringstream oss;
        oss << "init_array: trouble allocating memory, size = " << size << "\n";
        throw std::runtime_error(oss.str());
    }
    memset(array, 0, size * (size_t)sizeof(double));
    return (array);
}
}  // namespace psi
