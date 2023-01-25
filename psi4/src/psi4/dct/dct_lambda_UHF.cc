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
#include "psi4/psifiles.h"

#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/psifiles.h"

#include <cmath>

namespace psi {
namespace dct {

/**
 * Computes the residual for the lambda equations
 * R = G + F
 * @return RMS residual
 */
double DCTSolver::compute_cumulant_residual() {
    dct_timer_on("DCTSolver::compute_lambda_residual()");

    dpdbuf4 R, G, F, I, V, W;
    double sumSQ = 0.0;
    size_t nElements = 0;

    /*
     * R_ijab = G_ijab + F_ijab + gbar_ijab
     */

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // OOVV

    // R_IJAB = G_IJAB
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    global_dpd_->buf4_copy(&G, PSIF_DCT_DPD, "R <OO|VV>");
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");

    // R_IJAB += F_IJAB
    global_dpd_->buf4_init(&F, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
    dpd_buf4_add(&R, &F, 1.0);
    global_dpd_->buf4_close(&F);

    // R_IJAB += gbar_IJAB
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    dpd_buf4_add(&R, &I, 1.0);
    global_dpd_->buf4_close(&I);

    // Add third-order N-representability terms if needed
    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        // R_IJAB += V_IJAB
        global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "V <OO|VV>");
        dpd_buf4_add(&R, &V, 1.0);
        global_dpd_->buf4_close(&V);

        // R_IJAB += W_IJAB
        global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "W <OO|VV>");
        dpd_buf4_add(&R, &W, 1.0);
        global_dpd_->buf4_close(&W);
    }

    for (int h = 0; h < nirrep_; ++h) nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += global_dpd_->buf4_dot_self(&R);
    global_dpd_->buf4_close(&R);

    // OoVv

    // R_IjAb = G_IjAb
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    global_dpd_->buf4_copy(&G, PSIF_DCT_DPD, "R <Oo|Vv>");
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");

    // R_IjAb += F_IjAb
    global_dpd_->buf4_init(&F, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    dpd_buf4_add(&R, &F, 1.0);
    global_dpd_->buf4_close(&F);

    // R_IjAb += gbar_IjAb
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "MO Ints <Oo|Vv>");
    dpd_buf4_add(&R, &I, 1.0);
    global_dpd_->buf4_close(&I);

    // Add third-order N-representability terms if needed
    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        // R_IjAb += V_IjAb
        global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
        dpd_buf4_add(&R, &V, 1.0);
        global_dpd_->buf4_close(&V);

        // R_IJAB += W_IJAB
        global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");
        dpd_buf4_add(&R, &W, 1.0);
        global_dpd_->buf4_close(&W);
    }

    for (int h = 0; h < nirrep_; ++h) nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += global_dpd_->buf4_dot_self(&R);
    global_dpd_->buf4_close(&R);

    // oovv

    // R_ijab = G_ijab
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    global_dpd_->buf4_copy(&G, PSIF_DCT_DPD, "R <oo|vv>");
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");

    // R_ijab += F_ijab
    global_dpd_->buf4_init(&F, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
    dpd_buf4_add(&R, &F, 1.0);
    global_dpd_->buf4_close(&F);

    // R_ijab += gbar_ijab
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 1,
                           "MO Ints <oo|vv>");
    dpd_buf4_add(&R, &I, 1.0);
    global_dpd_->buf4_close(&I);

    // Add third-order N-representability terms if needed
    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        // R_ijab += V_ijab
        global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "V <oo|vv>");
        dpd_buf4_add(&R, &V, 1.0);
        global_dpd_->buf4_close(&V);

        // R_IJAB += W_IJAB
        global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "W <oo|vv>");
        dpd_buf4_add(&R, &W, 1.0);
        global_dpd_->buf4_close(&W);
    }

    for (int h = 0; h < nirrep_; ++h) nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += global_dpd_->buf4_dot_self(&R);
    global_dpd_->buf4_close(&R);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dct_timer_off("DCTSolver::compute_lambda_residual()");

    if (nElements > 0)
        return sqrt(sumSQ / nElements);
    else
        return 0.0;
}

/**
 * Builds the new lambda tensor from the intermediates
 */
void DCTSolver::update_cumulant_jacobi() {
    dct_timer_on("DCTSolver::update_lambda_from_residual()");

    dpdbuf4 L, D, R;

    /*
     * Amplitude_ijab += R_ijab / D_ijab
     */

    // L_IJAB += R_IJAB / D_IJAB
    global_dpd_->buf4_init(&D, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
    global_dpd_->buf4_dirprd(&D, &R);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    dpd_buf4_add(&L, &R, 1.0);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&L);

    // L_IjAb += R_IjAb / D_IjAb
    global_dpd_->buf4_init(&D, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
    global_dpd_->buf4_dirprd(&D, &R);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    dpd_buf4_add(&L, &R, 1.0);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&L);

    // L_IJAB += R_ijab / D_ijab
    global_dpd_->buf4_init(&D, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>=o]+"), ID("[v>=v]+"), 0, "D <oo|vv>");
    global_dpd_->buf4_init(&R, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
    global_dpd_->buf4_dirprd(&D, &R);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    dpd_buf4_add(&L, &R, 1.0);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&L);

    dct_timer_off("DCTSolver::update_lambda_from_residual()");
}

}  // namespace dct
}  // namespace psi
