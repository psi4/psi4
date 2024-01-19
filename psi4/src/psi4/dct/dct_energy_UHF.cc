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

#include "dct.h"
#include "psi4/psifiles.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"

namespace psi {
namespace dct {

/**
 * Uses the intermediates to compute the energy
 */
void DCTSolver::compute_dct_energy() {
    dct_timer_on("DCTSolver::compute_dct_energy()");

    dpdbuf4 L, G, I;
    double eGaa, eGab, eGbb, eIaa, eIab, eIbb;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // E += 1/4 L_IJAB G_IJAB
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    eGaa = global_dpd_->buf4_dot(&G, &L);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&L);

    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

        // E += gbar_IJAB Gamma_IJAB
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "Gamma <OO|VV>");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 1,
                               "MO Ints <OO|VV>");
        eIaa = 4.0 * global_dpd_->buf4_dot(&I, &G);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&G);

    } else {
        // E += 1/2 gbar_IJAB L_IJAB
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "Amplitude <OO|VV>");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 1,
                               "MO Ints <OO|VV>");
        eIaa = 2.0 * global_dpd_->buf4_dot(&I, &L);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
    }

    // E += L_IjAb G_IjAb
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    eGab = global_dpd_->buf4_dot(&G, &L);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&L);

    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        // E += 4.0 * gbar_IjAb Gamma_IjAb
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Gamma <Oo|Vv>");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "MO Ints <Oo|Vv>");
        eIab = 4.0 * global_dpd_->buf4_dot(&I, &G);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&G);
    } else {
        // E += 2.0 gbar_IjAb L_IjAb
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Amplitude <Oo|Vv>");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "MO Ints <Oo|Vv>");
        eIab = 2.0 * global_dpd_->buf4_dot(&I, &L);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
    }

    // E += 1/4 L_ijab G_ijab
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    eGbb = global_dpd_->buf4_dot(&G, &L);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&L);

    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        // E += gbar_ijab Gamma_ijab
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "Gamma <oo|vv>");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 1,
                               "MO Ints <oo|vv>");
        eIbb = 4.0 * global_dpd_->buf4_dot(&I, &G);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&G);

        psio_->close(PSIF_DCT_DENSITY, 1);
    } else {
        // E += 1/2 gbar_ijab L_ijab
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "Amplitude <oo|vv>");
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 1,
                               "MO Ints <oo|vv>");
        eIbb = 2.0 * global_dpd_->buf4_dot(&I, &L);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
    }
    psio_->close(PSIF_LIBTRANS_DPD, 1);

#if PRINT_ENERGY_COMPONENTS
    outfile->Printf("\tAA G Energy = %20.12f\n", eGaa);
    outfile->Printf("\tAB G Energy = %20.12f\n", eGab);
    outfile->Printf("\tBB G Energy = %20.12f\n", eGbb);
    outfile->Printf("\tAA I Energy = %20.12f\n", eIaa);
    outfile->Printf("\tAB I Energy = %20.12f\n", eIab);
    outfile->Printf("\tBB I Energy = %20.12f\n", eIbb);
    outfile->Printf("\tTotal G Energy = %20.12f\n", eGaa + eGab + eGbb);
    outfile->Printf("\tTotal I Energy = %20.12f\n", eIaa + eIab + eIbb);
#endif

    lambda_energy_ = eGaa + eGab + eGbb + eIaa + eIab + eIbb;

    dct_timer_off("DCTSolver::compute_dct_energy()");
}

}  // namespace dct
}  // namespace psi
