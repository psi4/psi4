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

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

double CCEnergyWavefunction::diagnostic() {
    dpdfile2 T1A, T1B;

    auto nirreps = moinfo_.nirreps;
    auto clsdpi = moinfo_.clsdpi;
    auto uoccpi = moinfo_.uoccpi;
    auto openpi = moinfo_.openpi;

    Dimension occpi, virtpi;
    int *occ_sym, *vir_sym;
    if (params_.ref != 2) {
        occpi = moinfo_.occpi;
        virtpi = moinfo_.virtpi;
        occ_sym = moinfo_.occ_sym;
        vir_sym = moinfo_.vir_sym;
    }

    /* Compute the number of electrons */
    auto num_elec_a = 0;
    auto num_elec_b = 0;
    for (int h = 0; h < nirreps; h++) {
        num_elec_a += clsdpi[h] + openpi[h];
        num_elec_b += clsdpi[h];
    }
    auto num_elec = num_elec_a + num_elec_b;

    auto t1diag = 0.0;

    if (params_.ref == 0) { /** RHF **/

        global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
        t1diag = global_dpd_->file2_dot_self(&T1A);
        global_dpd_->file2_close(&T1A);

        t1diag /= num_elec;
        t1diag = sqrt(t1diag);

    } else if (params_.ref == 1) { /** ROHF **/

        global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_mat_init(&T1A);
        global_dpd_->file2_mat_rd(&T1A);
        global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
        global_dpd_->file2_mat_init(&T1B);
        global_dpd_->file2_mat_rd(&T1B);

        for (int h = 0; h < nirreps; h++) {
            for (int i = 0; i < (occpi[h] - openpi[h]); i++) {
                for (int a = 0; a < (virtpi[h] - openpi[h]); a++) {
                    t1diag += (T1A.matrix[h][i][a] + T1B.matrix[h][i][a]) * (T1A.matrix[h][i][a] + T1B.matrix[h][i][a]);
                }
            }

            for (int i = 0; i < (occpi[h] - openpi[h]); i++) {
                for (int a = 0; a < openpi[h]; a++) {
                    auto A = a + uoccpi[h];

                    t1diag += 2 * T1B.matrix[h][i][A] * T1B.matrix[h][i][A];
                }
            }

            for (int i = 0; i < openpi[h]; i++) {
                auto I = i + clsdpi[h];
                for (int a = 0; a < (virtpi[h] - openpi[h]); a++) {
                    t1diag += 2 * T1A.matrix[h][I][a] * T1A.matrix[h][I][a];
                }
            }
        }

        t1diag /= num_elec;
        t1diag = sqrt(t1diag);
        t1diag *= 0.5;

        global_dpd_->file2_mat_close(&T1A);
        global_dpd_->file2_close(&T1A);
        global_dpd_->file2_mat_close(&T1B);
        global_dpd_->file2_close(&T1B);

    } else if (params_.ref == 2) { /** UHF **/

        global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_mat_init(&T1A);
        global_dpd_->file2_mat_rd(&T1A);
        global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
        global_dpd_->file2_mat_init(&T1B);
        global_dpd_->file2_mat_rd(&T1B);

        auto t1diag_a = 0.0;
        auto t1diag_b = 0.0;
        for (int h = 0; h < nirreps; h++) {
            for (int row = 0; row < T1A.params->rowtot[h]; row++)
                for (int col = 0; col < T1A.params->coltot[h]; col++)
                    t1diag_a += (T1A.matrix[h][row][col] * T1A.matrix[h][row][col]);

            for (int row = 0; row < T1B.params->rowtot[h]; row++)
                for (int col = 0; col < T1B.params->coltot[h]; col++)
                    t1diag_b += (T1B.matrix[h][row][col] * T1B.matrix[h][row][col]);
        }

        t1diag = sqrt((t1diag_a + t1diag_b) / (num_elec_a + num_elec_b));

        global_dpd_->file2_mat_close(&T1A);
        global_dpd_->file2_mat_close(&T1B);
        global_dpd_->file2_close(&T1A);
        global_dpd_->file2_close(&T1B);
    }

    return t1diag;
}
}  // namespace ccenergy
}  // namespace psi
