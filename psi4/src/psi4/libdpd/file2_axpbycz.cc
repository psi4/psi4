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
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

/* file2_axpbycz(): Evaluates the standard operation aX + bY -> cZ for
** dpdfile2's.
**
** Arguments:
**   dpdfile2 *FileA: A pointer to the leftmost dpdfile2.
**   dpdfile2 *FileB: A pointer to the rightmost summand dpdfile2.
**   dpdfile2 *FileC: A pointer to the target dpdfile2.
**   double a, b, c, scalar prefactors
*/

int DPD::file2_axpbycz(dpdfile2 *FileA, dpdfile2 *FileB, dpdfile2 *FileC,
                       double a, double b, double c)
{
    file2_scm(FileC, c);

    file2_axpy(FileB, FileC, b, 0);

    file2_axpy(FileA, FileC, a, 0);
    return 0;
}


}
