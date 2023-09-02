/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "dct.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/psifiles.h"
#include "psi4/psifiles.h"

#include <cmath>

namespace psi {
namespace dct {

/**
 * Computes the residual for the lambda equations
 * R = G + F
 * @return RMS residual
 */
double DCTSolver::compute_cumulant_residual_RHF() {
    dct_timer_on("DCTSolver::compute_lambda_residual()");

    dpdbuf4 R, G, F;
    double sumSQ = 0.0;
    size_t nElements = 0;

    /*
     * R_ijab = G_ijab + F_ijab
     */

    // R_IjAb = G_IjAb
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "G <OO|VV>");                       // G <Oo|Vv>
    global_dpd_->buf4_copy(&G, PSIF_DCT_DPD, "R SF <OO|VV>");  // R <Oo|Vv>
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "R SF <OO|VV>");  // R <Oo|Vv>

    // R_IjAb += F_IjAb
    global_dpd_->buf4_init(&F, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "F <OO|VV>");  // F <Oo|Vv>
    dpd_buf4_add(&R, &F, 1.0);
    global_dpd_->buf4_close(&F);
    for (int h = 0; h < nirrep_; ++h) nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += global_dpd_->buf4_dot_self(&R);
    global_dpd_->buf4_close(&R);

    dct_timer_off("DCTSolver::compute_lambda_residual()");

    if (nElements > 0)
        return sqrt(sumSQ / nElements);
    else
        return 0.0;
}

/**
 * Builds the new lambda tensor from the intermediates
 */
void DCTSolver::update_cumulant_jacobi_RHF() {
    dct_timer_on("DCTSolver::update_lambda_from_residual()");

    dpdbuf4 L, D, R;

    /*
     * Amplitude_ijab += R_ijab / D_ijab
     */

    // L_IjAb += R_IjAb / D_IjAb
    global_dpd_->buf4_init(&D, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0,
                           "D <OO|VV>");  // D <Oo|Vv>
    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "R SF <OO|VV>");  // R <Oo|Vv>
    global_dpd_->buf4_dirprd(&D, &R);
    global_dpd_->buf4_close(&D);
    // Update new cumulant
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>
    dpd_buf4_add(&L, &R, 1.0);
    global_dpd_->buf4_close(&L);

    global_dpd_->buf4_close(&R);

    /* update lambda <OO|VV> for tau and G intermediates */
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           "Amplitude SF <OO|VV>");
    global_dpd_->buf4_copy(&L, PSIF_DCT_DPD, "Amplitude <OO|VV>");
    global_dpd_->buf4_copy(&L, PSIF_DCT_DPD, "Amplitude <oo|vv>");
    global_dpd_->buf4_close(&L);

    dct_timer_off("DCTSolver::update_lambda_from_residual()");
}

/**
 * Compute R_OOVV and R_oovv from R_OoVv, used as DIIS error vectors
 * this is an unnecessary step, but can reduce # Iterations by around 5
 */
void DCTSolver::compute_R_AA_and_BB() {
    dct_timer_on("DCTSolver::compute_R_AA_and_BB");
    dpdbuf4 R;

    /*
     * R_IJAB = R_IjAb - R_JiAb
     * Copy R_IJAB to R_ijab
     */
    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1, "R SF <OO|VV>");
    global_dpd_->buf4_copy(&R, PSIF_DCT_DPD, "R <OO|VV>");
    global_dpd_->buf4_copy(&R, PSIF_DCT_DPD, "R <oo|vv>");
    global_dpd_->buf4_close(&R);

    dct_timer_off("DCTSolver::compute_R_AA_and_BB");
}

}  // namespace dct
}  // namespace psi
