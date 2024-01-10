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

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "dpd.h"

namespace psi {

int DPD::trans4_init(dpdtrans4 *Trans, dpdbuf4 *Buf) {
    int nirreps;

    nirreps = Buf->params->nirreps;

    /* Assign the input dpdbuf */
    Trans->buf = *Buf;

    Trans->matrix = (double ***)malloc(nirreps * sizeof(double **));

    /* Set up shifted matrix info */
    Trans->shift.shift_type = 0;
    Trans->shift.rowtot = init_int_matrix(nirreps, nirreps);
    Trans->shift.coltot = init_int_matrix(nirreps, nirreps);
    Trans->shift.matrix = (double ****)malloc(nirreps * sizeof(double ***));

    return 0;
}

}  // namespace psi
