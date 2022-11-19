/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace ccenergy {

/* Computes the D2 diagnostic as defined in:
 *
 * I.M.B. Nielsen and C.L. Janssen, Chem. Phys. Lett. 310, 568 (1999).
 *
 * */

double CCEnergyWavefunction::d2diag_rhf() {
    double *Eo;
    double *Ev;
    dpdbuf4 Tikab, Tjkab;
    dpdbuf4 Tijac, Tijbc;
    dpdfile2 To, Tv;

    auto nirreps = moinfo_.nirreps;
    auto max = 0.0;

    global_dpd_->buf4_init(&Tikab, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&Tjkab, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&To, PSIF_CC_TMP0, 0, 0, 0, "To");
    global_dpd_->buf4_init(&Tijac, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&Tijbc, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&Tv, PSIF_CC_TMP0, 0, 1, 1, "Tv");

    // Build diagnostic matrices To and Tv //
    global_dpd_->contract442(&Tikab, &Tjkab, &To, 0, 0, 1, 0);
    global_dpd_->contract442(&Tijac, &Tijbc, &Tv, 3, 3, 1, 0);

    global_dpd_->buf4_close(&Tikab);
    global_dpd_->buf4_close(&Tjkab);
    global_dpd_->file2_close(&To);
    global_dpd_->buf4_close(&Tijac);
    global_dpd_->buf4_close(&Tijbc);
    global_dpd_->file2_close(&Tv);

    global_dpd_->file2_init(&To, PSIF_CC_TMP0, 0, 0, 0, "To");
    global_dpd_->file2_mat_init(&To);
    global_dpd_->file2_mat_rd(&To);

    global_dpd_->file2_init(&Tv, PSIF_CC_TMP0, 0, 1, 1, "Tv");
    global_dpd_->file2_mat_init(&Tv);
    global_dpd_->file2_mat_rd(&Tv);

    for (int h = 0; h < nirreps; h++) {
        if (To.params->rowtot[h]) {
            // Diagonalize To //
            Eo = init_array(To.params->rowtot[h]);
            if (DSYEV_ascending(To.params->rowtot[h], To.matrix[h], Eo) != 0){
                throw PSIEXCEPTION("DSYEV diagonalizer failed in D2 diagnostic!");
            }
            // Find maximum To eigenvalue //
            for (int i = 0; i < To.params->rowtot[h]; i++) {
                if (Eo[i] > max) max = Eo[i];
            }
            free(Eo);
        }

        if (Tv.params->rowtot[h]) {
            // Diagonalize Tv //
            Ev = init_array(Tv.params->rowtot[h]);
            if (DSYEV_ascending(Tv.params->rowtot[h], Tv.matrix[h], Ev) != 0){
                throw PSIEXCEPTION("DSYEV diagonalizer failed in D2 diagnostic!");
            }
            // Find maximum Tv eigenvalue //
            for (int i = 0; i < Tv.params->rowtot[h]; i++) {
                if (Ev[i] > max) max = Ev[i];
            }
            free(Ev);
        }
    }

    global_dpd_->file2_mat_close(&To);
    global_dpd_->file2_mat_close(&Tv);
    global_dpd_->file2_close(&To);
    global_dpd_->file2_close(&Tv);

    /*
    // Original algorithm
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&T2, CC_TMP0, qrsp, 10, 11, "tjAbI");
    dpd_buf4_close(&T2);

  //  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  //  dpd_buf4_sort(&T2, CC_TMP0, rpqs, 11, 10, "tAIjb");
  //  dpd_buf4_sort(&T2, CC_TMP0, pqsr, 0, 5, "tIjbA");
  //  dpd_buf4_close(&T2);

    dpd_buf4_init(&Tijab, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&Tjabi, CC_TMP0, 0, 10, 11, 10, 11, 0, "tjAbI");

    for(int h = 0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&Tijab, h);
      dpd_buf4_mat_irrep_rd(&Tijab, h);
  //    dpd_buf4_mat_irrep_shift13(&Tijab, h);

      dpd_buf4_mat_irrep_init(&Tjabi, h);
      dpd_buf4_mat_irrep_rd(&Tjabi, h);
      dpd_bufg_mat_irrep_shift31(&Tjabi, h);

      nrows = moinfo.occpi[h];
      ncols = moinfo.occpi[h] * moinfo.virtpi[h] * moinfo.virtpi[h];
      To = block_matrix(nrows, nrows);
      C_DGEMM('n', 'n', nrows, nrows, ncols, 1.0, Tijab.shift.matrix[h][h][0],
              ncols, Tjabi.shift.matrix[h][h][0], nrows, 0.0, To[0], nrows);

      Eo = init_array(nrows);
      Co = block_matrix(nrows, nrows);
      sq_rsp(nrows, nrows, To, Eo, 0, Co, 1e-12);

      for(int i = 0; i < nrows; i++) {
        if(Eo[i] > max) max = Eo[i];
      }

      dpd_buf4_mat_irrep_close(&Tijab, h);
      dpd_buf4_mat_irrep_close(&Tjabi, h);

      free_block(To);
      free_block(Co);
      free(Eo);
    }

    dpd_buf4_close(&Tijab);
    dpd_buf4_close(&Tjabi);
    // END original algorithm
    */

    max = sqrt(max);

    return max;
}

double CCEnergyWavefunction::d2diag() {
    double norm = 0.0;

    if (params_.ref == 0) { /** RHF **/
        norm = d2diag_rhf();
    }
    return norm;

    // There's no open shell definitions, but I've set it up just incase -TJM
}
}  // namespace ccenergy
}  // namespace psi
