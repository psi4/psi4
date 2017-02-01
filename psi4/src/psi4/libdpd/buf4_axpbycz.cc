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

/* buf4_axpbycz(): Evaluates the standard operation aX + bY -> cZ for
** dpdbuf4's.
**
** Arguments:
**   dpdbuf4 *FileA: A pointer to the leftmost dpdbuf4.
**   dpdbuf4 *FileB: A pointer to the rightmost summand dpdbuf4.
**   dpdbuf4 *FileC: A pointer to the target dpdbuf4.
**   double a, b, c, scalar prefactors
*/

int DPD::buf4_axpbycz(dpdbuf4 *FileA, dpdbuf4 *FileB, dpdbuf4 *FileC,
                      double a, double b, double c)
{
    buf4_scm(FileC, c);

    buf4_axpy(FileB, FileC, b);

    buf4_axpy(FileA, FileC, a);
    return 0;
}

}
