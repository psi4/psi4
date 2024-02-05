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

#include "dpd.h"

namespace psi {

double DPD::buf4_trace(dpdbuf4 *Buf) {
    double trace = 0.0;
    for (int h = 0; h < Buf->params->nirreps; h++) {
        if (Buf->params->rowtot[h] == Buf->params->coltot[h]) {
            buf4_mat_irrep_init(Buf, h);
            buf4_mat_irrep_rd(Buf, h);
            for (int row = 0; row < Buf->params->rowtot[h]; row++) trace += Buf->matrix[h][row][row];
            buf4_mat_irrep_close(Buf, h);
        }
    }

    return trace;
}

}  // namespace psi
