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

double DPD::file2_dot(dpdfile2 *FileA, dpdfile2 *FileB)
{
    int h, nirreps, my_irrep;
    double dot;

    nirreps = FileA->params->nirreps;
    my_irrep = FileA->my_irrep;

    dot = 0.0;

    file2_mat_init(FileA);
    file2_mat_init(FileB);
    file2_mat_rd(FileA);
    file2_mat_rd(FileB);

    for(h=0; h < nirreps; h++) {

        dot += dot_block(FileA->matrix[h], FileB->matrix[h],
                         FileA->params->rowtot[h],
                         FileA->params->coltot[h^my_irrep], 1.0);

    }

    file2_mat_close(FileA);
    file2_mat_close(FileB);

    return dot;

}

}
