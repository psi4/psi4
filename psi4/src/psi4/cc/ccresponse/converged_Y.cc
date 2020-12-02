/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

double converged_Y(const char *pert, int irrep, double omega) {
    dpdfile2 Y1, Y1new;
    dpdbuf4 Y2, Y2new;
    double rms = 0.0, value;
    int row, col, h, nirreps;
    char lbl[32];

    nirreps = moinfo.nirreps;

    sprintf(lbl, "New Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_mat_init(&Y1new);
    global_dpd_->file2_mat_rd(&Y1new);
    sprintf(lbl, "Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_mat_init(&Y1);
    global_dpd_->file2_mat_rd(&Y1);

    for (h = 0; h < nirreps; h++)
        for (row = 0; row < Y1.params->rowtot[h]; row++)
            for (col = 0; col < Y1.params->coltot[h ^ irrep]; col++) {
                value = Y1new.matrix[h][row][col] - Y1.matrix[h][row][col];
                rms += value * value;
            }
    global_dpd_->file2_mat_close(&Y1new);
    global_dpd_->file2_close(&Y1new);
    global_dpd_->file2_mat_close(&Y1);
    global_dpd_->file2_close(&Y1);

    sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    for (h = 0; h < nirreps; h++) {
        global_dpd_->buf4_mat_irrep_init(&Y2new, h);
        global_dpd_->buf4_mat_irrep_rd(&Y2new, h);
        global_dpd_->buf4_mat_irrep_init(&Y2, h);
        global_dpd_->buf4_mat_irrep_rd(&Y2, h);

        for (row = 0; row < Y2.params->rowtot[h]; row++)
            for (col = 0; col < Y2.params->coltot[h ^ irrep]; col++) {
                value = Y2new.matrix[h][row][col] - Y2.matrix[h][row][col];
                rms += value * value;
            }

        global_dpd_->buf4_mat_irrep_close(&Y2new, h);
        global_dpd_->buf4_mat_irrep_close(&Y2, h);
    }
    global_dpd_->buf4_close(&Y2new);
    global_dpd_->buf4_close(&Y2);

    return sqrt(rms);
}

}  // namespace ccresponse
}  // namespace psi
