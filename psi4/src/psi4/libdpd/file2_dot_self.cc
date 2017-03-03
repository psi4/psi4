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
#include "dpd.h"

namespace psi {

/* dpd_file2_dot_self(): Evaluates the sum of the squares of the elements of a
** given dpdfile2.
*/

double DPD::file2_dot_self(dpdfile2 *BufX)
{
    int h, nirreps, my_irrep;
    int row, col;
    double alpha=0.0;

    nirreps = BufX->params->nirreps;
    my_irrep = BufX->my_irrep;

    file2_mat_init(BufX);
    file2_mat_rd(BufX);

    for(h=0; h < nirreps; h++) {

        for(row=0; row < BufX->params->rowtot[h]; row++)
            for(col=0; col < BufX->params->coltot[h^my_irrep]; col++)
                alpha += BufX->matrix[h][row][col] * BufX->matrix[h][row][col];

    }

    file2_mat_close(BufX);

    return alpha;

}

}
