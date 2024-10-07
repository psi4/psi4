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

#include "psi4/psifiles.h"
#include "dct.h"

#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/liboptions/liboptions.h"

namespace psi {
namespace dct {

/**
 * Builds the intermediate tensors
 */
void DCTSolver::build_cumulant_intermediates() {
    dct_timer_on("DCTSolver::build_intermediates()");

    compute_G_intermediate();

    if (exact_tau_) {
        form_density_weighted_fock();
    }

    compute_F_intermediate();

    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        // Add third-order N-representability terms
        compute_V_intermediate();

        // Fourth-order one-electron terms
        compute_W_intermediate();
    }

    dct_timer_off("DCTSolver::build_intermediates()");
}

void DCTSolver::compute_G_intermediate() {
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 I, L, G, T, Taa, Tab, Tbb, Laa, Lab, Lbb;

    /*
     * G_ijab = 1/2 Sum_cd gbar_cdab lambda_ijcd
     */
    dct_timer_on("DCTSolver::gbar_cdab lambda_ijcd");
    if (options_.get_str("AO_BASIS") == "NONE" && options_.get_str("DCT_TYPE") == "CONV") {
        // G_IJAB += 1/2 Sum_CD gbar_CDAB lambda_IJCD
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>V]-"), ID("[V>V]-"), ID("[V,V]"), ID("[V,V]"), 1,
                               "MO Ints <VV|VV>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "Amplitude <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "G <OO|VV>");
        global_dpd_->contract444(&L, &I, &G, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        // G_IjAb += Sum_Cd g_CdAb lambda_IjCd
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), 0,
                               "MO Ints <Vv|Vv>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Amplitude <Oo|Vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
        global_dpd_->contract444(&L, &I, &G, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        // G_ijab += 1/2 Sum_cd gbar_cdab lambda_ijcd
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v>v]-"), ID("[v>v]-"), ID("[v,v]"), ID("[v,v]"), 1,
                               "MO Ints <vv|vv>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "Amplitude <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "G <oo|vv>");
        global_dpd_->contract444(&L, &I, &G, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);
    } else {
        /***********AA***********/
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                               "tau(temp) <OO|VV>");
        global_dpd_->buf4_copy(&L, PSIF_DCT_DPD, "G <OO|VV>");
        global_dpd_->buf4_close(&L);

        /***********BB***********/
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                               "tau(temp) <oo|vv>");
        global_dpd_->buf4_copy(&L, PSIF_DCT_DPD, "G <oo|vv>");
        global_dpd_->buf4_close(&L);

        /***********AB***********/
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "tau(temp) <Oo|Vv>");
        global_dpd_->buf4_copy(&L, PSIF_DCT_DPD, "G <Oo|Vv>");
        global_dpd_->buf4_close(&L);
    }
    dct_timer_off("DCTSolver::gbar_cdab lambda_ijcd");

    /*
     * G_ijab += 1/2 Sum_kl gbar_ijkl lambda_klab
     */
    // G_IJAB += 1/2 Sum_KL gbar_IJKL lambda_KLAB
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O,O]"), ID("[O,O]"), 1,
                           "MO Ints <OO|OO>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    global_dpd_->contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // G_IjAb += Sum_Kl g_IjKl lambda_KlAb
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0,
                           "MO Ints <Oo|Oo>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    global_dpd_->contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // G_ijab += 1/2 Sum_kl gbar_ijkl lambda_klab
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o,o]"), ID("[o,o]"), 1,
                           "MO Ints <oo|oo>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    global_dpd_->contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    /*
     * G_ijab -= P(ij)P(ab) Sum_kc gbar_jckb lambda_ikac
     */
    global_dpd_->buf4_init(&Taa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Temp (OV|OV)");
    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           "Temp (Ov|oV)");
    global_dpd_->buf4_init(&Tbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Temp (ov|ov)");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_sort(&Laa, PSIF_DCT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Amplitude (OV|OV)");
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "Amplitude (Ov|oV)");
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           "Amplitude (Ov|oV)");

    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_sort(&Lbb, PSIF_DCT_DPD, prqs, ID("[o,v]"), ID("[o,v]"), "Amplitude (ov|ov)");
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Amplitude (ov|ov)");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0,
                           "MO Ints <oV|oV>");

    // T_IbjA = -Sum_kC lambda_IbkC g_jAkC
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, -1.0, 0.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0,
                           "MO Ints <Ov|Ov>");

    // T_IbjA -= Sum_Kc g_IbKc lambda_KcjA
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);

    // T_IbjA -> T_IAjb
    global_dpd_->buf4_sort(&Tab, PSIF_DCT_DPD, psrq, ID("[O,V]"), ID("[o,v]"), "Temp (OV|ov)");
    global_dpd_->buf4_close(&Tab);
    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Temp (OV|ov)");

    // Amplitude_IbkC -> Amplitude_ICkb
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DPD, psrq, ID("[O,V]"), ID("[o,v]"), "Amplitude (OV|ov)");
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "MO Ints (OV|ov)");

    // T_IAJB = Sum_kc lambda_IAkc g_JBkc
    global_dpd_->contract444(&Lab, &I, &Taa, 0, 0, 1.0, 0.0);

    // T_iajb = Sum_KC lambda_KCia g_KCjb
    global_dpd_->contract444(&Lab, &I, &Tbb, 1, 1, 1.0, 0.0);

    // T_IAjb += Sum_kc g_IAkc lambda_jbkc
    global_dpd_->contract444(&I, &Lbb, &Tab, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC lambda_IAKC g_KCjb
    global_dpd_->contract444(&Laa, &I, &Tab, 0, 1, 1.0, 1.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints <OV|OV>");

    // T_IAJB -= Sum_KC lambda_IAKC g_JBKC
    global_dpd_->contract444(&Laa, &I, &Taa, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC g_IAKC lambda_KCjb
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints (OV|OV)");

    // T_IAJB += Sum_KC lambda_IAKC (JB|KC)
    global_dpd_->contract444(&Laa, &I, &Taa, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC (JB|KC) lambda_KCjb
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "MO Ints <ov|ov>");

    // T_iajb -= Sum_kc lambda_iakc g_jbkc
    global_dpd_->contract444(&Lbb, &I, &Tbb, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC lambda_IAkc g_jbkc
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "MO Ints (ov|ov)");

    // T_iajb += Sum_kc lambda_iakc (JB|KC)
    global_dpd_->contract444(&Lbb, &I, &Tbb, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC lambda_IAkc (kc|jb)
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);

    // T_IAJB -> T_IJAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    global_dpd_->buf4_close(&Taa);
    // G_IJAB += T_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Taa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // G_IJAB -= T_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_IJBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB -= T_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_JIBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB += T_JIBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Taa);

    // T_IAjb -> T_IjAb
    global_dpd_->buf4_sort(&Tab, PSIF_DCT_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&Tab);
    // G_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_iajb -> T_ijab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    global_dpd_->buf4_close(&Tbb);
    // G_ijab += T_ijab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");

    // T_ijab -> T_jiab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_ijba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_jiba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab += T_jiba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Tbb);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void DCTSolver::compute_F_intermediate() {
    dpdfile2 F_OO, F_oo, F_VV, F_vv;
    dpdbuf4 F, T, Laa, Lab, Lbb;

    /*
     * F_ijab += P(ab) F_ca lambda_ijcb - P(ij) F_ki lambda_jkab
     */
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    // Temp_IJAB = lambda_IJCB F_AC
    global_dpd_->file2_init(&F_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->contract244(&F_VV, &Laa, &T, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_VV);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&T);
    // F_IJAB = Temp_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");
    global_dpd_->buf4_copy(&T, PSIF_DCT_DPD, "F <OO|VV>");
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&F, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);
    // F_IJAB -= Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    // Temp_IJAB = -lambda_KJAB F_IK
    global_dpd_->file2_init(&F_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract244(&F_OO, &Laa, &T, 1, 0, 0, -1.0, 0.0);

    global_dpd_->file2_close(&F_OO);
    global_dpd_->buf4_close(&Laa);
    // F_IJAB += Temp_IJAB
    dpd_buf4_add(&F, &T, 1.0);
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);
    // F_IJAB -= Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    // F_IjAb += lambda_IjCb F_AC
    global_dpd_->file2_init(&F_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->contract244(&F_VV, &Lab, &F, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_VV);
    // F_IjAb += lambda_IjAc F_bc
    global_dpd_->file2_init(&F_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    global_dpd_->contract424(&Lab, &F_vv, &F, 3, 1, 0, 1.0, 1.0);
    global_dpd_->file2_close(&F_vv);
    // F_IjAb -= lambda_KjAb F_IK
    global_dpd_->file2_init(&F_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract244(&F_OO, &Lab, &F, 1, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&F_OO);
    // F_IjAb -= lambda_IkAb F_jk
    global_dpd_->file2_init(&F_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->contract424(&Lab, &F_oo, &F, 1, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    // Temp_ijab = lambda_ijcb F_ac
    global_dpd_->file2_init(&F_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    global_dpd_->contract244(&F_vv, &Lbb, &T, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_vv);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");
    // F_ijab = Temp_ijab
    global_dpd_->buf4_copy(&T, PSIF_DCT_DPD, "F <oo|vv>");
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&F, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);
    // F_ijab -= Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    // Temp_ijab = -lambda_kjab F_ik
    global_dpd_->file2_init(&F_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->contract244(&F_oo, &Lbb, &T, 1, 0, 0, -1.0, 0.0);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->buf4_close(&Lbb);
    // F_ijab += Temp_ijab
    dpd_buf4_add(&F, &T, 1.0);
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);
    // F_ijab -= Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&F);
}

void DCTSolver::form_density_weighted_fock() {
    dpdfile2 T_OO, T_oo, T_VV, T_vv, F_OO, F_oo, F_VV, F_vv;

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    // Copy Tau in MO basis from the DPD library
    auto a_tau_mo = Matrix("Alpha Tau in the MO basis", nmopi_, nmopi_);
    auto b_tau_mo = Matrix("Beta Tau in the MO basis", nmopi_, nmopi_);

    a_tau_mo.set_block(slices_.at("ACTIVE_OCC_A"), Matrix(&T_OO));
    a_tau_mo.set_block(slices_.at("ACTIVE_VIR_A"), Matrix(&T_VV));
    b_tau_mo.set_block(slices_.at("ACTIVE_OCC_B"), Matrix(&T_oo));
    b_tau_mo.set_block(slices_.at("ACTIVE_VIR_B"), Matrix(&T_vv));

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    auto a_evecs = std::make_shared<Matrix>("Tau Eigenvectors (Alpha)", nmopi_, nmopi_);
    auto b_evecs = std::make_shared<Matrix>("Tau Eigenvectors (Beta)", nmopi_, nmopi_);
    auto a_evals = std::make_shared<Vector>("Tau Eigenvalues (Alpha)",  nmopi_);
    auto b_evals = std::make_shared<Vector>("Tau Eigenvalues (Beta)", nmopi_);

    // Diagonalize Tau
    a_tau_mo.diagonalize(a_evecs, a_evals);
    a_tau_mo.zero();
    a_tau_mo.set_diagonal(a_evals);
    b_tau_mo.diagonalize(b_evecs, b_evals);
    b_tau_mo.zero();
    b_tau_mo.set_diagonal(b_evals);

    // Transform the Fock matrix to NSO basis
    auto nso_Fa = std::make_shared<Matrix>("Alpha Fock in the NSO basis", nirrep_, nmopi_, nmopi_);
    auto nso_Fb = std::make_shared<Matrix>("Beta Fock in the NSO basis", nirrep_, nmopi_, nmopi_);

    // Alpha spin
    nso_Fa->transform(moFa_, a_evecs);
    // Beta spin
    nso_Fb->transform(moFb_, b_evecs);

    // Form density-weighted Fock matrix in the NSO basis

    for (int h = 0; h < nirrep_; ++h) {
        if (nsopi_[h] == 0) continue;

        // Alpha occupied
        for (int ip = 0; ip < naoccpi_[h]; ++ip) {
            for (int kp = 0; kp < naoccpi_[h]; ++kp) {
                double value = nso_Fa->get(h, ip, kp) / (1.0 + a_evals->get(h, ip) + a_evals->get(h, kp));
                nso_Fa->set(h, ip, kp, value);
            }
        }

        // Beta occupied
        for (int ip = 0; ip < nboccpi_[h]; ++ip) {
            for (int kp = 0; kp < nboccpi_[h]; ++kp) {
                double value = nso_Fb->get(h, ip, kp) / (1.0 + b_evals->get(h, ip) + b_evals->get(h, kp));
                nso_Fb->set(h, ip, kp, value);
            }
        }

        // Alpha virtual
        for (int ap = naoccpi_[h]; ap < nmopi_[h]; ++ap) {
            for (int dp = naoccpi_[h]; dp < nmopi_[h]; ++dp) {
                double value = nso_Fa->get(h, dp, ap) / (1.0 - a_evals->get(h, dp) - a_evals->get(h, ap));
                nso_Fa->set(h, dp, ap, value);
            }
        }

        // Beta virtual
        for (int ap = nboccpi_[h]; ap < nmopi_[h]; ++ap) {
            for (int dp = nboccpi_[h]; dp < nmopi_[h]; ++dp) {
                double value = nso_Fb->get(h, dp, ap) / (1.0 - b_evals->get(h, dp) - b_evals->get(h, ap));
                nso_Fb->set(h, dp, ap, value);
            }
        }
    }

    // Transform density-weighted density matrix back to the original basis

    // Alpha spin
    nso_Fa->back_transform(a_evecs);
    // Beta spin
    nso_Fb->back_transform(b_evecs);

    Ftilde_a_->copy(nso_Fa);
    Ftilde_b_->copy(nso_Fb);

    global_dpd_->file2_init(&F_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->file2_init(&F_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->file2_init(&F_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->file2_init(&F_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "F <v|v>");

    global_dpd_->file2_mat_init(&F_OO);
    global_dpd_->file2_mat_init(&F_oo);
    global_dpd_->file2_mat_init(&F_VV);
    global_dpd_->file2_mat_init(&F_vv);

    for (int h = 0; h < nirrep_; ++h) {
        if (nsopi_[h] == 0) continue;

        // Alpha occupied
        for (int ip = 0; ip < naoccpi_[h]; ++ip) {
            for (int kp = 0; kp < naoccpi_[h]; ++kp) {
                F_OO.matrix[h][ip][kp] = nso_Fa->get(h, ip, kp);
            }
        }

        // Beta occupied
        for (int ip = 0; ip < nboccpi_[h]; ++ip) {
            for (int kp = 0; kp < nboccpi_[h]; ++kp) {
                F_oo.matrix[h][ip][kp] = nso_Fb->get(h, ip, kp);
            }
        }

        // Alpha virtual
        for (int ap = 0; ap < navirpi_[h]; ++ap) {
            for (int dp = 0; dp < navirpi_[h]; ++dp) {
                F_VV.matrix[h][ap][dp] = nso_Fa->get(h, ap + naoccpi_[h], dp + naoccpi_[h]);
                if (ap == dp) F_VV.matrix[h][ap][dp] += energy_level_shift_;
            }
        }

        // Beta virtual
        for (int ap = 0; ap < nbvirpi_[h]; ++ap) {
            for (int dp = 0; dp < nbvirpi_[h]; ++dp) {
                F_vv.matrix[h][ap][dp] = nso_Fb->get(h, ap + nboccpi_[h], dp + nboccpi_[h]);
                if (ap == dp) F_vv.matrix[h][ap][dp] += energy_level_shift_;
            }
        }
    }

    global_dpd_->file2_mat_wrt(&F_OO);
    global_dpd_->file2_mat_wrt(&F_oo);
    global_dpd_->file2_mat_wrt(&F_VV);
    global_dpd_->file2_mat_wrt(&F_vv);

    global_dpd_->file2_close(&F_OO);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->file2_close(&F_VV);
    global_dpd_->file2_close(&F_vv);
}

void DCTSolver::compute_V_intermediate() {
    dpdbuf4 I, L, V, T, II, J, Taa, Tab, Tbb, Kaa, Kab, Kbb, Laa, Lab, Lbb;
    dpdfile2 T_OO, T_oo, T_VV, T_vv, H_OO, H_oo, H_VV, H_vv;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "T <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "T <v|v>");

    /*
     * V_ijab = 1/3 P_(ij) gbar_abik T_kj
     */

    // OOVV

    // V_IJAB = 1/3 gbar_IKAB * T_KJ
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    global_dpd_->contract424(&I, &T_OO, &T, 1, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    global_dpd_->buf4_scm(&V, 0.0);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_IJAB -= 1/3 gbar_JKAB * T_KI
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // OoVv

    // V_IjAb += 1/3 gbar_IkAb * T_kj
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "MO Ints <Oo|Vv>");
    global_dpd_->contract424(&I, &T_oo, &V, 1, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // V_IjAb += 1/3 T_IK * gbar_KjAb
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "MO Ints <Oo|Vv>");
    global_dpd_->contract244(&T_OO, &I, &V, 1, 0, 0, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // oovv

    // V_ijab = 1/3 gbar_ikab * T_kj
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 1,
                           "MO Ints <oo|vv>");
    global_dpd_->contract424(&I, &T_oo, &T, 1, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    global_dpd_->buf4_scm(&V, 0.0);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_ijab -= 1/3 gbar_jkab * T_ki
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    /*
     * V_ijab -= 1/3 P_(ab) gbar_ijac T_cb
     */

    // OOVV

    // V_IJAB -= 1/3 gbar_IJAC * T_CB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    global_dpd_->contract424(&I, &T_VV, &T, 3, 0, 0, -1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_IJAB += 1/3 gbar_IJBC * T_CA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // OoVv

    // V_IjAb -= 1/3 gbar_IjAc * T_cb
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "MO Ints <Oo|Vv>");
    global_dpd_->contract424(&I, &T_vv, &V, 3, 0, 0, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // V_IjAb -= 1/3 gbar_IjCb * T_CA
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "MO Ints <Oo|Vv>");
    global_dpd_->contract244(&T_VV, &I, &V, 1, 2, 1, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // oovv

    // V_ijab -= 1/3 gbar_ijac * T_cb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 1,
                           "MO Ints <oo|vv>");
    global_dpd_->contract424(&I, &T_vv, &T, 3, 0, 0, -1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_ijab += 1/3 gbar_ijbc * T_ca
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    /*
     * V_ijab += 1/3 P_(ij) lambda_abik * H_kj
     */

    compute_H_intermediate();

    global_dpd_->file2_init(&H_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&H_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_init(&H_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&H_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "H <v|v>");

    // OOVV

    // V_IJAB += 1/3 lambda_IKAB * H_KJ
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract424(&L, &H_OO, &T, 1, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_IJAB -= 1/3 lambda_JKAB * H_KI
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // OoVv

    // V_IjAb += 1/3 lambda_IkAb * H_kj
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract424(&L, &H_oo, &V, 1, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // V_IjAb += 1/3 H_IK * lambda_KjAb
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract244(&H_OO, &L, &V, 1, 0, 0, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // oovv

    // V_ijab += 1/3 lambda_ikab * H_kj
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract424(&L, &H_oo, &T, 1, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_ijab -= 1/3 lambda_jkab * H_ki
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    /*
     * V_ijab -= 1/3 P_(ab) lambda_ijac H_cb
     */

    // OOVV

    // V_IJAB -= 1/3 lambda_IJAC * H_CB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract424(&L, &H_VV, &T, 3, 0, 0, -1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_IJAB += 1/3 lambda_IJBC * H_CA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // OoVv

    // V_IjAb -= 1/3 lambda_IjAc * H_cb
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract424(&L, &H_vv, &V, 3, 0, 0, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // V_IjAb -= 1/3 lambda_IjCb * H_CA
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract244(&H_VV, &L, &V, 1, 2, 1, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // oovv

    // V_ijab -= 1/3 lambda_ijac * H_cb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract424(&L, &H_vv, &T, 3, 0, 0, -1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_ijab += 1/3 lambda_ijbc * H_ca
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->file2_close(&H_OO);
    global_dpd_->file2_close(&H_oo);
    global_dpd_->file2_close(&H_VV);
    global_dpd_->file2_close(&H_vv);

    /*
     * V_ijab += 1/6 gbar_abkl I_klij
     */

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // V_IJAB += 1/6 * gbar_ABKL * I_KLIJ
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           "I <OO|OO>");
    global_dpd_->contract444(&II, &I, &V, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // V_IjAb += 1/3 * gbar_AbKl * I_KlIj
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "MO Ints <Oo|Vv>");
    global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    global_dpd_->contract444(&II, &I, &V, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // V_ijab += 1/6 * gbar_abkl * I_klij
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 1,
                           "MO Ints <oo|vv>");
    global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           "I <oo|oo>");
    global_dpd_->contract444(&II, &I, &V, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    /*
     * V_ijab += 1/6 lambda_abkl J_klij
     */

    compute_J_intermediate();

    // V_IJAB += 1/6 * lambda_ABKL * J_KLIJ
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&J, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0, "J <OO|OO>");
    global_dpd_->contract444(&J, &L, &V, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&J);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // V_IjAb += 1/3 * lambda_AbKl * J_KlIj
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&J, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "J <Oo|Oo>");
    global_dpd_->contract444(&J, &L, &V, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&J);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // V_ijab += 1/6 * lambda_abkl * J_klij
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&J, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0, "J <oo|oo>");
    global_dpd_->contract444(&J, &L, &V, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&J);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    /*
     * V_ijab += 2/3 P_(ij) P_(ab) gbar_acik K_kbjc
     */

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // OOVV
    global_dpd_->buf4_init(&Taa, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "Temp (VO|OV)");

    global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");

    // T_AIJB = 2/3 (AI|KC) K_(KC|JB)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "MO Ints (VO|OV)");
    global_dpd_->contract444(&I, &Kaa, &Taa, 0, 0, 2.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);

    // T_AIJB -= 2/3 <AI|KC> K_(KC|JB)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "MO Ints <VO|OV>");
    global_dpd_->contract444(&I, &Kaa, &Taa, 0, 0, -2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);

    // T_AIJB += 2/3 (AI|kc) K_(kc|JB)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0,
                           "MO Ints (VO|ov)");
    global_dpd_->contract444(&I, &Kab, &Taa, 0, 0, 2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Kaa);
    global_dpd_->buf4_close(&Kab);

    // T_AIJB -> T_IJAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qrps, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    global_dpd_->buf4_close(&Taa);

    // V_IJAB += T_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Taa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB -= T_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_IJBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB -= T_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_JIBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB += T_JIBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Taa);

    // OoVv
    global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");
    global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");

    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0,
                           "Temp (VO|ov)");

    // T_AIjb = 2/3 (AI|kc) K_(kc|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0,
                           "MO Ints (VO|ov)");
    global_dpd_->contract444(&I, &Kbb, &Tab, 0, 1, 2.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);

    // T_AIjb += 2/3 (AI|KC) K_(KC|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "MO Ints (VO|OV)");
    global_dpd_->contract444(&I, &Kab, &Tab, 0, 1, 2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);

    // T_AIjb -= 2/3 <AI|KC> K_(KC|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "MO Ints <VO|OV>");
    global_dpd_->contract444(&I, &Kab, &Tab, 0, 1, -2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Tab);

    // T_AIjb -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0, "Temp (VO|ov)");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qrps, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), 0,
                           "Temp (OV|vo)");

    // T_IAbj = 2/3 (bj|KC) K_(KC|IA)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), 0,
                           "MO Ints (vo|OV)");
    global_dpd_->contract444(&Kaa, &I, &Tab, 0, 0, 2.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);

    // T_IAbj += 2/3 (bj|kc) K_(kc|IA)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "MO Ints (vo|ov)");
    global_dpd_->contract444(&Kab, &I, &Tab, 0, 0, 2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);

    // T_IAbj -= 2/3 <bj|kc> K_(kc|IA)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "MO Ints <vo|ov>");
    global_dpd_->contract444(&Kab, &I, &Tab, 0, 0, -2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Tab);

    // T_IAbj -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), 0, "Temp (OV|vo)");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, psqr, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Kaa);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Kbb);

    // T_IbAj = 2/3 K_(Ib|Kc) <Kc|Aj>
    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "Temp (Ov|Vo)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "MO Ints <Ov|Vo>");
    global_dpd_->contract444(&Kab, &I, &Tab, 0, 1, 2.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Tab);

    // T_IbAj += 2/3 <Ib|Ck> K_(Ck|Aj)
    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "Temp (Ov|Vo)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), 0, "K <Vo|Vo>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "MO Ints <Ov|Vo>");
    global_dpd_->contract444(&I, &Kab, &Tab, 0, 1, 2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Tab);

    // T_IbAj -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0, "Temp (Ov|Vo)");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, psrq, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // oovv
    global_dpd_->buf4_init(&Tbb, PSIF_DCT_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "Temp (vo|ov)");

    global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");

    // T_aijb = 2/3 (ai|kc) K_(kc|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "MO Ints (vo|ov)");
    global_dpd_->contract444(&I, &Kbb, &Tbb, 0, 0, 2.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&I);

    // T_aijb -= 2/3 <ai|kc> K_(kc|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "MO Ints <vo|ov>");
    global_dpd_->contract444(&I, &Kbb, &Tbb, 0, 0, -2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);

    // T_aijb += 2/3 (ai|KC) K_(KC|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), 0,
                           "MO Ints (vo|OV)");
    global_dpd_->contract444(&I, &Kab, &Tbb, 0, 1, 2.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Kbb);

    // T_aijb -> T_ijab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qrps, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    global_dpd_->buf4_close(&Tbb);

    // V_ijab += T_ijab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");

    // T_ijab -> T_jiab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab -= T_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_ijba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab -= T_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_jiba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab += T_jiba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Tbb);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    /*
     * V_ijab += 1/3 P_(ij) P_(ab) lambda_acik L_kbjc
     */

    compute_L_intermediate();

    // Compute the residual contribution
    // OOVV
    global_dpd_->buf4_init(&Taa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), 0,
                           "Temp (OV|VO)");

    // T_IABJ = 1/3 Amplitude_(IA|KC) L_(KC|BJ)
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0, "L (VO|O'V')");
    global_dpd_->contract444(&L, &Laa, &Taa, 0, 0, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&L);

    // T_IABJ += 1/3 Amplitude_(IA|kc) L_(BJ|kc)
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0, "L (VO|o'v')");
    global_dpd_->contract444(&L, &Lab, &Taa, 0, 0, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&L);

    // T_IABJ -> T_IJAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, psqr, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    global_dpd_->buf4_close(&Taa);

    // V_IJAB += T_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Taa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB -= T_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_IJBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB -= T_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_JIBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB += T_JIBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Taa);

    // OoVv

    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), 0,
                           "Temp (OV|vo)");

    // T_IAbj = 1/3 Amplitude_(IA|KC) L_(KC|bj)
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), 0, "L (O'V'|vo)");
    global_dpd_->contract444(&L, &Lab, &Tab, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&L);

    // T_IAbj += 1/3 Amplitude_(IA|kc) L_(bj|kc)
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0, "L (vo|o'v')");
    global_dpd_->contract444(&L, &Lbb, &Tab, 0, 0, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&L);

    global_dpd_->buf4_close(&Tab);

    // T_IAbj -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), 0, "Temp (OV|vo)");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, psqr, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), 0,
                           "Temp (ov|VO)");

    // T_jbAI = 1/3 Amplitude_(KC|jb) L_(AI|KC)
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0, "L (VO|O'V')");
    global_dpd_->contract444(&L, &Laa, &Tab, 1, 0, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&L);

    // T_jbAI += 1/3 Amplitude_(kc|jb) L_(AI|kc)
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Amplitude (ov|ov)");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0, "L (VO|o'v')");
    global_dpd_->contract444(&L, &Lab, &Tab, 1, 0, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&L);

    global_dpd_->buf4_close(&Tab);

    // T_jbAI -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), 0, "Temp (ov|VO)");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, sprq, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "Temp (Ov|Vo)");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "Amplitude (Ov|Vo)");

    // T_IbAj = 1/3 L_<Ib|Kc> Amplitude_(Kc|Aj)
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "L <Ov|O'v'>");
    global_dpd_->contract444(&Lab, &L, &Tab, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&Lab);

    // T_IbAj += 1/3 Amplitude_(Ib|Ck) L_<Aj|Ck>
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), 0, "L <Vo|V'o'>");
    global_dpd_->contract444(&L, &Lab, &Tab, 0, 0, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&Lab);

    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Tab);

    // T_IbAj -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0, "Temp (Ov|Vo)");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, psrq, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // oovv
    global_dpd_->buf4_init(&Tbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), 0,
                           "Temp (ov|vo)");

    // T_iabj = 1/3 Amplitude_(ia|kc) L_(kc|bj)
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Amplitude (ov|ov)");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0, "L (vo|o'v')");
    global_dpd_->contract444(&L, &Lbb, &Tbb, 0, 0, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&L);

    // T_iabj += 1/3 Amplitude_(KC|ia) L_(KC|bj)
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), 0, "L (O'V'|vo)");
    global_dpd_->contract444(&L, &Lab, &Tbb, 1, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&L);

    // T_iabj -> T_ijab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, psqr, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    global_dpd_->buf4_close(&Tbb);

    // V_ijab += T_ijab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");

    // T_ijab -> T_jiab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab -= T_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_ijba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab -= T_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_jiba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab += T_jiba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Tbb);
}

void DCTSolver::compute_H_intermediate() {
    dpdbuf4 L, I;
    dpdfile2 H_OO, H_oo, H_VV, H_vv, HH_OO, HH_oo, HH_VV, HH_vv;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->file2_init(&H_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&H_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_init(&H_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&H_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "H <v|v>");

    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    /*
     * H_IJ = -1/2 Amplitude_IKAB gbar_JKAB
     */
    global_dpd_->contract442(&L, &I, &H_OO, 0, 0, -1.0, 0.0);
    /*
     * H_AB = +1/2 Amplitude_IJAC gbar_IJBC
     */
    global_dpd_->contract442(&L, &I, &H_VV, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 1,
                           "MO Ints <oo|vv>");
    /*
     * H_ij = -1/2 Amplitude_ikab gbar_jkab
     */
    global_dpd_->contract442(&L, &I, &H_oo, 0, 0, -1.0, 0.0);
    /*
     * H_ab = +1/2 Amplitude_ijac gbar_ijbc
     */
    global_dpd_->contract442(&L, &I, &H_vv, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "MO Ints <Oo|Vv>");
    /*
     * H_IJ -= Amplitude_IkAb gbar_JkAb - Amplitude_IkaB gbar_JkaB
     */
    global_dpd_->contract442(&L, &I, &H_OO, 0, 0, -2.0, 1.0);
    /*
     * H_ij -= Amplitude_KiAb gbar_KjAb - Amplitude_KiaB gbar_KjaB
     */
    global_dpd_->contract442(&L, &I, &H_oo, 1, 1, -2.0, 1.0);
    /*
     * H_AB += Amplitude_IjAc gbar_IjBc + Amplitude_iJAc gbar_iJBc
     */
    global_dpd_->contract442(&L, &I, &H_VV, 2, 2, 2.0, 1.0);
    /*
     * H_ab += Amplitude_IjCa gbar_IjCb + Amplitude_iJCa gbar_iJCb
     */
    global_dpd_->contract442(&L, &I, &H_vv, 3, 3, 2.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&I);

    global_dpd_->file2_close(&H_OO);
    global_dpd_->file2_close(&H_oo);
    global_dpd_->file2_close(&H_VV);
    global_dpd_->file2_close(&H_vv);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // Symmetrize H
    global_dpd_->file2_init(&H_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&H_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_init(&H_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&H_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "H <v|v>");
    global_dpd_->file2_init(&HH_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&HH_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_init(&HH_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&HH_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "H <v|v>");

    global_dpd_->file2_axpy(&HH_OO, &H_OO, 1.0, 1);
    global_dpd_->file2_axpy(&HH_oo, &H_oo, 1.0, 1);
    global_dpd_->file2_axpy(&HH_VV, &H_VV, 1.0, 1);
    global_dpd_->file2_axpy(&HH_vv, &H_vv, 1.0, 1);

    global_dpd_->file2_close(&HH_OO);
    global_dpd_->file2_close(&HH_oo);
    global_dpd_->file2_close(&HH_VV);
    global_dpd_->file2_close(&HH_vv);
    global_dpd_->file2_close(&H_OO);
    global_dpd_->file2_close(&H_oo);
    global_dpd_->file2_close(&H_VV);
    global_dpd_->file2_close(&H_vv);
}

void DCTSolver::compute_I_intermediate() {
    dpdbuf4 LLaa, LLab, LLbb, Laa, Lab, Lbb, Iaa, Iab, Ibb;

    // I_ijkl = Amplitude_ijab * Amplitude_klab
    global_dpd_->buf4_init(&Iaa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           "I <OO|OO>");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&LLaa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract444(&Laa, &LLaa, &Iaa, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&LLaa);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Iaa);

    global_dpd_->buf4_init(&Iab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&LLab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract444(&Lab, &LLab, &Iab, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Iab);

    global_dpd_->buf4_init(&Ibb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           "I <oo|oo>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&LLbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract444(&Lbb, &LLbb, &Ibb, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&LLbb);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&Ibb);

    // Sort the I intermediate to chemist's notation
    global_dpd_->buf4_init(&Iaa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>O]-"), ID("[O>O]-"), 0, "I <OO|OO>");
    global_dpd_->buf4_sort(&Iaa, PSIF_DCT_DPD, prqs, ID("[O,O]"), ID("[O,O]"), "I (OO|OO)");
    global_dpd_->buf4_close(&Iaa);

    global_dpd_->buf4_init(&Iab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    global_dpd_->buf4_sort(&Iab, PSIF_DCT_DPD, prqs, ID("[O,O]"), ID("[o,o]"), "I (OO|oo)");
    global_dpd_->buf4_close(&Iab);

    global_dpd_->buf4_init(&Iab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[o,o]"), ID("[O,O]"), ID("[o,o]"), 0, "I (OO|oo)");
    global_dpd_->buf4_sort(&Iab, PSIF_DCT_DPD, rspq, ID("[o,o]"), ID("[O,O]"), "I (oo|OO)");
    global_dpd_->buf4_close(&Iab);

    global_dpd_->buf4_init(&Ibb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[o,o]"), ID("[o>o]-"), ID("[o>o]-"), 0, "I <oo|oo>");
    global_dpd_->buf4_sort(&Ibb, PSIF_DCT_DPD, prqs, ID("[o,o]"), ID("[o,o]"), "I (oo|oo)");
    global_dpd_->buf4_close(&Ibb);
}

void DCTSolver::compute_J_intermediate() {
    dpdbuf4 Iaa, Iab, Ibb, Laa, Lab, Lbb, Jaa, Jab, Jbb;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // J_ijkl = Amplitude_ijab * gbar_klab + gbar_ijab * lambda_klab
    global_dpd_->buf4_init(&Jaa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           "J <OO|OO>");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&Iaa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    global_dpd_->contract444(&Laa, &Iaa, &Jaa, 0, 0, 4.0, 0.0);
    global_dpd_->buf4_close(&Iaa);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Jaa);

    global_dpd_->buf4_init(&Jab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "J <Oo|Oo>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&Iab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "MO Ints <Oo|Vv>");
    global_dpd_->contract444(&Lab, &Iab, &Jab, 0, 0, 4.0, 0.0);
    global_dpd_->buf4_close(&Iab);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Jab);

    global_dpd_->buf4_init(&Jbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           "J <oo|oo>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&Ibb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 1,
                           "MO Ints <oo|vv>");
    global_dpd_->contract444(&Lbb, &Ibb, &Jbb, 0, 0, 4.0, 0.0);
    global_dpd_->buf4_close(&Ibb);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&Jbb);

    // Bra-ket symmetrize J
    global_dpd_->buf4_init(&Jaa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           "J <OO|OO>");
    global_dpd_->buf4_symm(&Jaa);
    global_dpd_->buf4_close(&Jaa);

    global_dpd_->buf4_init(&Jab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "J <Oo|Oo>");
    global_dpd_->buf4_symm(&Jab);
    global_dpd_->buf4_close(&Jab);

    global_dpd_->buf4_init(&Jbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           "J <oo|oo>");
    global_dpd_->buf4_symm(&Jbb);
    global_dpd_->buf4_close(&Jbb);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void DCTSolver::compute_K_intermediate() {
    dpdbuf4 LLaa, LLab, LLbb, Laa, Lab, Lbb, Kaa, Kab, Kba, Kbb, K;

    // There are five unique spin cases: K<IAJB>, K<iajb>, K<IaJb>, K<iAjB>, K<IajB>

    // Sort the cumulant

    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_sort(&Laa, PSIF_DCT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Amplitude (OV|OV)");
    global_dpd_->buf4_close(&Laa);

    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "Amplitude (Ov|oV)");
    global_dpd_->buf4_close(&Lab);

    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           "Amplitude (Ov|oV)");
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DPD, psrq, ID("[O,V]"), ID("[o,v]"), "Amplitude (OV|ov)");
    global_dpd_->buf4_close(&Lab);

    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_sort(&Lbb, PSIF_DCT_DPD, prqs, ID("[o,v]"), ID("[o,v]"), "Amplitude (ov|ov)");
    global_dpd_->buf4_close(&Lbb);

    // K<IAJB> spin case

    global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->buf4_init(&LLaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->contract444(&Laa, &LLaa, &Kaa, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&LLaa);
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&LLab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->contract444(&Lab, &LLab, &Kaa, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Kaa);

    // Resort K(OV|OV) to the K<OV|OV>
    global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_sort(&Kaa, PSIF_DCT_DPD, psrq, ID("[O,V]"), ID("[O,V]"), "K <OV|OV>");
    global_dpd_->buf4_close(&Kaa);

    // Resort K<OV|OV> to the K(OO|VV)
    global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K <OV|OV>");
    global_dpd_->buf4_sort(&Kaa, PSIF_DCT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "K (OO|VV)");
    global_dpd_->buf4_close(&Kaa);

    // K<IaJb> and K<iAjB> spin cases:

    // Although we denote K <Ov|Ov> and K <Vo|Vo> as in physist's notation, it is actually in chemist's notation.
    // However, it is convenient to store that intermediate in this form to avoid unnecessary tensor sorts for the
    // density computation
    // E.g. K_(Kb|Ic) <- Amplitude_(Kb|mE) Amplitude_(Ic|mE) = Amplitude_<Km|Eb> Amplitude_<Im|Ec>
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           "Amplitude (Ov|oV)");
    global_dpd_->buf4_init(&LLab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           "Amplitude (Ov|oV)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->contract444(&Lab, &LLab, &Kab, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_init(&Kba, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->contract444(&Lab, &LLab, &Kba, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Kba);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&LLab);

    // Resort K<Ov|Ov> to the K(OO|vv)
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->buf4_sort(&Kab, PSIF_DCT_DPD, prqs, ID("[O,O]"), ID("[v,v]"), "K (OO|vv)");
    global_dpd_->buf4_close(&Kab);

    // Resort K<oV|oV> to the K(oo|VV)
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->buf4_sort(&Kab, PSIF_DCT_DPD, prqs, ID("[o,o]"), ID("[V,V]"), "K (oo|VV)");
    global_dpd_->buf4_close(&Kab);

    // K_kCjA -> K_CkAj
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->buf4_sort(&Kab, PSIF_DCT_DPD, qpsr, ID("[V,o]"), ID("[V,o]"), "K <Vo|Vo>");
    global_dpd_->buf4_close(&Kab);

    // K<IajB> spin case:

    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->contract444(&Laa, &Lab, &Kab, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Amplitude (ov|ov)");
    global_dpd_->contract444(&Lab, &Lbb, &Kab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");
    global_dpd_->buf4_sort(&Kab, PSIF_DCT_DPD, psrq, ID("[O,v]"), ID("[o,V]"), "K <Ov|oV>");
    // Resort to get the K_oVOv. Used for the MO Lagrangian
    global_dpd_->buf4_sort(&Kab, PSIF_DCT_DPD, rqps, ID("[o,V]"), ID("[O,v]"), "K <oV|Ov>");
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Lab);

    // K<iajb> spin case:

    global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Amplitude (ov|ov)");
    global_dpd_->buf4_init(&LLbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Amplitude (ov|ov)");
    global_dpd_->contract444(&Lbb, &LLbb, &Kbb, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&LLbb);
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&LLab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->contract444(&Lab, &LLab, &Kbb, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Kbb);

    // Resort K(ov|ov) to the K<ov|ov>
    global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");
    global_dpd_->buf4_sort(&Kbb, PSIF_DCT_DPD, psrq, ID("[o,v]"), ID("[o,v]"), "K <ov|ov>");
    global_dpd_->buf4_close(&Kbb);

    // Resort K<ov|ov> to the K(oo|vv)
    global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K <ov|ov>");
    global_dpd_->buf4_sort(&Kbb, PSIF_DCT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "K (oo|vv)");
    global_dpd_->buf4_close(&Kbb);

    // Resort all chemist's notation K intermediates to VVOO format (needed for fourth-order Tau terms)
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K (OO|VV)");
    global_dpd_->buf4_sort(&K, PSIF_DCT_DPD, rspq, ID("[V,V]"), ID("[O,O]"), "K (VV|OO)");
    global_dpd_->buf4_close(&K);

    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[v,v]"), ID("[O,O]"), ID("[v,v]"), 0, "K (OO|vv)");
    global_dpd_->buf4_sort(&K, PSIF_DCT_DPD, rspq, ID("[v,v]"), ID("[O,O]"), "K (vv|OO)");
    global_dpd_->buf4_close(&K);

    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[V,V]"), ID("[o,o]"), ID("[V,V]"), 0, "K (oo|VV)");
    global_dpd_->buf4_sort(&K, PSIF_DCT_DPD, rspq, ID("[V,V]"), ID("[o,o]"), "K (VV|oo)");
    global_dpd_->buf4_close(&K);

    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "K (oo|vv)");
    global_dpd_->buf4_sort(&K, PSIF_DCT_DPD, rspq, ID("[v,v]"), ID("[o,o]"), "K (vv|oo)");
    global_dpd_->buf4_close(&K);
}

void DCTSolver::compute_L_intermediate() {
    dpdbuf4 I, Iaa, Iab, Ibb, Amplitude_aa, Amplitude_ab, Amplitude_bb, Laa, Lab, Lba, Lbb;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // There are six unique spin cases: L<IAJB>, L<iajb>, L<IaJb>, L<iAjB>, L<IajB>, L<iAJb>
    // Primed indices denote that they are located on ERI tensor

    // L<IAJB> spin case
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), 0, "L (OV|V'O')");

    global_dpd_->buf4_init(&Amplitude_aa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");

    // L_(IB|AJ) = Amplitude_(IB|KC) * (AJ|KC)
    global_dpd_->buf4_init(&Iaa, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "MO Ints (VO|OV)");
    global_dpd_->contract444(&Amplitude_aa, &Iaa, &Laa, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Iaa);

    // L_(IB|AJ) -= Amplitude_(IB|KC) * <AJ|KC>
    global_dpd_->buf4_init(&Iaa, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "MO Ints <VO|OV>");
    global_dpd_->contract444(&Amplitude_aa, &Iaa, &Laa, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&Iaa);

    global_dpd_->buf4_close(&Amplitude_aa);

    // L_(IB|AJ) += Amplitude_(IB|kc) * (AJ|kc)
    global_dpd_->buf4_init(&Amplitude_ab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&Iab, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0,
                           "MO Ints (VO|ov)");
    global_dpd_->contract444(&Amplitude_ab, &Iab, &Laa, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Iab);
    global_dpd_->buf4_close(&Amplitude_ab);

    global_dpd_->buf4_close(&Laa);

    // L<IaJb> and L<iAjB> spin cases:

    // Although we denote L <Ov|Ov> and L <Vo|Vo> as in physist's notation, it is actually in chemist's notation.
    // However, it is convenient to store that intermediate in this form to avoid unnecessary tensor sorts
    global_dpd_->buf4_init(&Amplitude_ab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           "Amplitude (Ov|oV)");
    global_dpd_->buf4_sort(&Amplitude_ab, PSIF_DCT_DPD, pqsr, ID("[O,v]"), ID("[V,o]"), "Amplitude (Ov|Vo)");
    global_dpd_->buf4_close(&Amplitude_ab);

    global_dpd_->buf4_init(&Amplitude_ab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "Amplitude (Ov|Vo)");
    global_dpd_->buf4_init(&Iab, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "MO Ints <Ov|Vo>");

    // L_<Ib|Ja> = Amplitude_(Ib|Ck) * <Ja|Ck>
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "L <Ov|O'v'>");
    global_dpd_->contract444(&Amplitude_ab, &Iab, &Lab, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lab);

    // L_<Bi|Aj> = Amplitude_(Kc|Bi) * <Kc|Aj>
    global_dpd_->buf4_init(&Lba, PSIF_DCT_DPD, 0, ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), 0, "L <Vo|V'o'>");
    global_dpd_->contract444(&Amplitude_ab, &Iab, &Lba, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Lba);

    global_dpd_->buf4_close(&Iab);
    global_dpd_->buf4_close(&Amplitude_ab);

    // L<IajB> spin case:
    global_dpd_->buf4_init(&Amplitude_aa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->buf4_init(&Amplitude_ab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&Amplitude_bb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Amplitude (ov|ov)");

    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), 0, "L (OV|v'o')");

    // L_(IB|aj) = Amplitude_(IB|kc) * (aj|kc)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "MO Ints (vo|ov)");
    global_dpd_->contract444(&Amplitude_ab, &I, &Lab, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);

    // L_(IB|aj) -= Amplitude_(IB|kc) * <aj|kc>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "MO Ints <vo|ov>");
    global_dpd_->contract444(&Amplitude_ab, &I, &Lab, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // L_(IB|aj) += Amplitude_(IB|KC) * (aj|KC)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), 0,
                           "MO Ints (vo|OV)");
    global_dpd_->contract444(&Amplitude_aa, &I, &Lab, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Lab);

    // L<iAJb> spin case:
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0, "L (V'O'|ov)");

    // L_(AJ|ib) += (AJ|KC) * Amplitude_(KC|ib)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "MO Ints (VO|OV)");
    global_dpd_->contract444(&I, &Amplitude_ab, &Lab, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&I);

    // L_(AJ|ib) -= <AJ|KC> * Amplitude_(KC|ib)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), 0,
                           "MO Ints <VO|OV>");
    global_dpd_->contract444(&I, &Amplitude_ab, &Lab, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // L_(AJ|ib) += (AJ|kc) * Amplitude_(ib|kc)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0,
                           "MO Ints (VO|ov)");
    global_dpd_->contract444(&I, &Amplitude_bb, &Lab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Lab);

    global_dpd_->buf4_close(&Amplitude_bb);
    global_dpd_->buf4_close(&Amplitude_ab);
    global_dpd_->buf4_close(&Amplitude_aa);

    // L<iajb> spin case
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), 0, "L (ov|v'o')");

    global_dpd_->buf4_init(&Amplitude_bb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Amplitude (ov|ov)");

    // L_(ib|aj) = Amplitude_(ib|kc) * (aj|kc)
    global_dpd_->buf4_init(&Ibb, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "MO Ints (vo|ov)");
    global_dpd_->contract444(&Amplitude_bb, &Ibb, &Lbb, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Ibb);

    // L_(ib|aj) -= Amplitude_(ib|kc) * <aj|kc>
    global_dpd_->buf4_init(&Ibb, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), 0,
                           "MO Ints <vo|ov>");
    global_dpd_->contract444(&Amplitude_bb, &Ibb, &Lbb, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&Ibb);

    global_dpd_->buf4_close(&Amplitude_bb);

    // L_(ib|aj) += Amplitude_(ib|KC) * (aj|KC)
    global_dpd_->buf4_init(&Amplitude_ab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Amplitude (OV|ov)");
    global_dpd_->buf4_init(&Iab, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), 0,
                           "MO Ints (vo|OV)");
    global_dpd_->contract444(&Amplitude_ab, &Iab, &Lbb, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Iab);
    global_dpd_->buf4_close(&Amplitude_ab);

    global_dpd_->buf4_close(&Lbb);

    // Perform necessary sorts of L intermediate
    // L_(JB|CK) -> L_(BJ|KC)
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[V,O]"), ID("[O,V]"), ID("[V,O]"), 0, "L (OV|V'O')");
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DPD, qpsr, ID("[V,O]"), ID("[O,V]"), "L (VO|O'V')");
    global_dpd_->buf4_close(&Lab);

    // L_(JB|ck) -> L_(BJ|kc)
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[v,o]"), ID("[O,V]"), ID("[v,o]"), 0, "L (OV|v'o')");
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DPD, qpsr, ID("[V,O]"), ID("[o,v]"), "L (VO|o'v')");
    global_dpd_->buf4_close(&Lab);

    // L_(JB|ck) -> L_(BJ|kc)
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[V,O]"), ID("[o,v]"), ID("[V,O]"), ID("[o,v]"), 0, "L (V'O'|ov)");
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DPD, qpsr, ID("[O,V]"), ID("[v,o]"), "L (O'V'|vo)");
    global_dpd_->buf4_close(&Lab);

    // L_(jb|ck) -> L_(bj|kc)
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[v,o]"), ID("[o,v]"), ID("[v,o]"), 0, "L (ov|v'o')");
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DPD, qpsr, ID("[v,o]"), ID("[o,v]"), "L (vo|o'v')");
    global_dpd_->buf4_close(&Lab);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void DCTSolver::compute_O_intermediate() {
    dpdbuf4 O, L, I;

    /*
     * O_ijab = Amplitude_abkl I_klij
     */

    // O_IJAB = Amplitude_ABKL * I_KLIJ
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "O <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0, "I <OO|OO>");
    global_dpd_->contract444(&I, &L, &O, 0, 1, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&O);

    // O_IjAb = Amplitude_AbKl * I_KlIj
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "O <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    global_dpd_->contract444(&I, &L, &O, 0, 1, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&O);

    // O_ijab = Amplitude_abkl * I_klij
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "O <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0, "I <oo|oo>");
    global_dpd_->contract444(&I, &L, &O, 0, 1, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&O);
}

void DCTSolver::compute_W_intermediate() {
    dpdbuf4 I, Laa, Lab, Lbb, W, T, TT, N, M, K, O, L, Kaa, Kab, Kbb;
    dpdfile2 T_OO, T_oo, T_VV, T_vv, F_OO, F_oo, F_VV, F_vv, Temp_OO, Temp_oo, Temp_VV, Temp_vv;

    // Compute M and N intermediates
    compute_M_intermediate();
    compute_N_intermediate();

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "T <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "T <v|v>");

    global_dpd_->file2_init(&F_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->file2_init(&F_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->file2_init(&F_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->file2_init(&F_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "F <v|v>");

    // Precompute Temp_im = 4 (Ft_il T_lm + T_il Ft_lm) - I_ilmk Ft_kl + 2 K_idmc Ft_cd
    global_dpd_->file2_init(&Temp_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    global_dpd_->file2_init(&Temp_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");

    // Temp_IM = 4.0 * (Ft_IL T_LM + T_IL F_LM)
    global_dpd_->contract222(&F_OO, &T_OO, &Temp_OO, 0, 1, 4.0, 0.0);
    global_dpd_->contract222(&T_OO, &F_OO, &Temp_OO, 0, 1, 4.0, 1.0);

    // Temp_im = 4.0 * (Ft_il T_lm + T_il Ft_lm)
    global_dpd_->contract222(&F_oo, &T_oo, &Temp_oo, 0, 1, 4.0, 0.0);
    global_dpd_->contract222(&T_oo, &F_oo, &Temp_oo, 0, 1, 4.0, 1.0);

    // Temp_IM -= I_(IM|LK) F_LK
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0, "I (OO|OO)");
    global_dpd_->contract422(&I, &F_OO, &Temp_OO, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // Temp_IM -= I_(IM|lk) F_lk
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[o,o]"), ID("[O,O]"), ID("[o,o]"), 0, "I (OO|oo)");
    global_dpd_->contract422(&I, &F_oo, &Temp_OO, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // Temp_im -= I_(im|LK) F_LK
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[O,O]"), ID("[o,o]"), ID("[O,O]"), 0, "I (oo|OO)");
    global_dpd_->contract422(&I, &F_OO, &Temp_oo, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // Temp_im -= I_(im|lk) F_lk
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[o,o]"), ID("[o,o]"), ID("[o,o]"), 0, "I (oo|oo)");
    global_dpd_->contract422(&I, &F_oo, &Temp_oo, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // Temp_IM += 2.0 * K_(IM|CD) * F_CD
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K (OO|VV)");
    global_dpd_->contract422(&K, &F_VV, &Temp_OO, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&K);

    // Temp_IM += 2.0 * K_(IM|cd) * F_cd
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[v,v]"), ID("[O,O]"), ID("[v,v]"), 0, "K (OO|vv)");
    global_dpd_->contract422(&K, &F_vv, &Temp_OO, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&K);

    // Temp_im += 2.0 * K_(im|CD) * F_CD
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[V,V]"), ID("[o,o]"), ID("[V,V]"), 0, "K (oo|VV)");
    global_dpd_->contract422(&K, &F_VV, &Temp_oo, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&K);

    // Temp_im += 2.0 * K_(im|cd) * F_cd
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "K (oo|vv)");
    global_dpd_->contract422(&K, &F_vv, &Temp_oo, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&K);

    global_dpd_->file2_close(&Temp_oo);
    global_dpd_->file2_close(&Temp_OO);

    // Precompute Temp_ae = 4 (Ft_da T_ed + T_dc Ft_ed) + M_klca Amplitude_cekl - 2 Ft_lk K_kela
    global_dpd_->file2_init(&Temp_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    global_dpd_->file2_init(&Temp_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");

    // Temp_AC = 4.0 * (Ft_AD T_DC + T_AD Ft_DC)
    global_dpd_->contract222(&F_VV, &T_VV, &Temp_VV, 0, 1, 4.0, 0.0);
    global_dpd_->contract222(&T_VV, &F_VV, &Temp_VV, 0, 1, 4.0, 1.0);

    // Temp_ac = 4.0 * (Ft_ad T_dc + T_ad Ft_dc)
    global_dpd_->contract222(&F_vv, &T_vv, &Temp_vv, 0, 1, 4.0, 0.0);
    global_dpd_->contract222(&T_vv, &F_vv, &Temp_vv, 0, 1, 4.0, 1.0);

    // Temp_AD += M_KLCA Amplitude_KLCD
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "M <OO|V'V>");
    global_dpd_->contract442(&M, &Laa, &Temp_VV, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&Laa);

    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");

    // Temp_AD += 2.0 * M_LkAc Amplitude_LkDc (M_OoVv')
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv'>");
    global_dpd_->contract442(&M, &Lab, &Temp_VV, 2, 2, 2.0, 1.0);
    global_dpd_->buf4_close(&M);

    // Temp_ad += 2.0 * M_KlCa Amplitude_KlCd (M_OoV'v)
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|V'v>");
    global_dpd_->contract442(&M, &Lab, &Temp_vv, 3, 3, 2.0, 1.0);
    global_dpd_->buf4_close(&M);

    global_dpd_->buf4_close(&Lab);

    // Temp_ad += M_klca Amplitude_klcd
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "M <oo|v'v>");
    global_dpd_->contract442(&M, &Lbb, &Temp_vv, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&Lbb);

    // Temp_AD -= 2.0 * K_(AD|KL) * F_KL
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), 0, "K (VV|OO)");
    global_dpd_->contract422(&K, &F_OO, &Temp_VV, 0, 0, -2.0, 1.0);
    global_dpd_->buf4_close(&K);

    // Temp_AD -= 2.0 * K_(AD|kl) * F_kl
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[V,V]"), ID("[o,o]"), ID("[V,V]"), ID("[o,o]"), 0, "K (VV|oo)");
    global_dpd_->contract422(&K, &F_oo, &Temp_VV, 0, 0, -2.0, 1.0);
    global_dpd_->buf4_close(&K);

    // Temp_ad -= 2.0 * K_(ad|KL) * F_KL
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[v,v]"), ID("[O,O]"), ID("[v,v]"), ID("[O,O]"), 0, "K (vv|OO)");
    global_dpd_->contract422(&K, &F_OO, &Temp_vv, 0, 0, -2.0, 1.0);
    global_dpd_->buf4_close(&K);

    // Temp_ad -= 2.0 * K_(ad|kl) * F_kl
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), 0, "K (vv|oo)");
    global_dpd_->contract422(&K, &F_oo, &Temp_vv, 0, 0, -2.0, 1.0);
    global_dpd_->buf4_close(&K);

    global_dpd_->file2_close(&Temp_vv);
    global_dpd_->file2_close(&Temp_VV);

    /*
     * 1. W_ijab = -1/6 P_(ij) [Temp_ik Amplitude_kjab]
     */

    /*
     * 2. W_ijab -= 1/6 P_(ab) [Temp_ca Amplitude_ijcb]
     */

    // OOVV
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    // Temp_IJAB = -1/6 lambda_IJCB Temp_AC
    global_dpd_->file2_init(&Temp_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    global_dpd_->contract244(&Temp_VV, &Laa, &T, 1, 2, 1, -1.0 / 6.0, 0.0);
    global_dpd_->file2_close(&Temp_VV);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&T);
    // W_IJAB = Temp_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");
    global_dpd_->buf4_copy(&T, PSIF_DCT_DPD, "W <OO|VV>");
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);
    // W_IJAB -= Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    // Temp_IJAB = -1/6 lambda_KJAB Temp_IK
    global_dpd_->file2_init(&Temp_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    global_dpd_->contract244(&Temp_OO, &Laa, &T, 1, 0, 0, -1.0 / 6.0, 0.0);

    global_dpd_->file2_close(&Temp_OO);
    global_dpd_->buf4_close(&Laa);
    // W_IJAB += Temp_IJAB
    dpd_buf4_add(&W, &T, 1.0);
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);
    // W_IJAB -= Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    // OoVv
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    // W_IjAb += -1/6 lambda_IjCb Temp_AC
    global_dpd_->file2_init(&Temp_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    global_dpd_->contract244(&Temp_VV, &Lab, &W, 1, 2, 1, -1.0 / 6.0, 0.0);
    global_dpd_->file2_close(&Temp_VV);
    // W_IjAb += -1/6 lambda_IjAc Temp_bc
    global_dpd_->file2_init(&Temp_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");
    global_dpd_->contract424(&Lab, &Temp_vv, &W, 3, 1, 0, -1.0 / 6.0, 1.0);
    global_dpd_->file2_close(&Temp_vv);
    // W_IjAb -= -1/6 lambda_KjAb Temp_IK
    global_dpd_->file2_init(&Temp_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    global_dpd_->contract244(&Temp_OO, &Lab, &W, 1, 0, 0, -1.0 / 6.0, 1.0);
    global_dpd_->file2_close(&Temp_OO);
    // W_IjAb -= lambda_IkAb Temp_jk
    global_dpd_->file2_init(&Temp_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    global_dpd_->contract424(&Lab, &Temp_oo, &W, 1, 1, 1, -1.0 / 6.0, 1.0);
    global_dpd_->file2_close(&Temp_oo);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&W);

    // oovv
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    // Temp_ijab = -1/6 lambda_ijcb Temp_ac
    global_dpd_->file2_init(&Temp_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");
    global_dpd_->contract244(&Temp_vv, &Lbb, &T, 1, 2, 1, -1.0 / 6.0, 0.0);
    global_dpd_->file2_close(&Temp_vv);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");
    // W_ijab = Temp_ijab
    global_dpd_->buf4_copy(&T, PSIF_DCT_DPD, "W <oo|vv>");
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);
    // W_ijab -= Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    // Temp_ijab = -1/6 lambda_kjab Temp_ik
    global_dpd_->file2_init(&Temp_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    global_dpd_->contract244(&Temp_oo, &Lbb, &T, 1, 0, 0, -1.0 / 6.0, 0.0);
    global_dpd_->file2_close(&Temp_oo);
    global_dpd_->buf4_close(&Lbb);
    // W_ijab += Temp_ijab
    dpd_buf4_add(&W, &T, 1.0);
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);
    // W_ijab -= Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    /*
     * 3. W_ijab += 1/6 P_(ij) [(P_(ab) [M_ikab] - N_ikab) T_jk]
     */

    // OOVV
    // Precompute Temp_IKAB = P_(AB) [M_IKAB] - N_IKAB
    // Temp_IJAB = M_IJAB
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "M <OO|V'V>");
    global_dpd_->buf4_copy(&M, PSIF_DCT_DPD, "Temp <OO|VV>");
    global_dpd_->buf4_close(&M);

    // Temp_IJAB -> Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -= Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&T, &TT, -1.0);
    global_dpd_->buf4_close(&TT);

    // Temp_IJAB -= N_IJAB
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "N <O'O|VV>");
    dpd_buf4_add(&T, &N, -1.0);
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_close(&T);

    // W_IJAB += 1/6 P_(IJ) [Temp_IKAB T_JK]
    // Temp2_IJAB = 1/6 Temp_IKAB T_JK
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp2 <OO|VV>");
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->contract424(&T, &T_OO, &TT, 1, 1, 1, 1.0 / 6.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&TT);

    // W_IJAB += Temp2_IJAB
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp2 <OO|VV>");
    dpd_buf4_add(&W, &TT, 1.0);
    global_dpd_->buf4_close(&TT);

    // Temp2_IJAB -> Temp2_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp2 <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    // W_IJAB -= Temp2_JIAB
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&W, &TT, -1.0);
    global_dpd_->buf4_close(&TT);
    global_dpd_->buf4_close(&W);

    // OoVv
    // Temp_IkAb = M_IkAb - M_IkbA = M_IkAb + M_kIbA = (M_OoV'v + M_OoVv')
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|V'v>");
    global_dpd_->buf4_copy(&M, PSIF_DCT_DPD, "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&M);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv'>");
    dpd_buf4_add(&T, &N, 1.0);
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    // W_IjAb += 1/6 Temp_IkAb T_jk
    global_dpd_->contract424(&T, &T_oo, &W, 1, 1, 1, 1.0 / 6.0, 1.0);
    // W_IjAb += 1/6 T_KI Temp_KjAb
    global_dpd_->contract244(&T_OO, &T, &W, 1, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&T);

    // W_IjAb += -1/6 N_IkAb T_jk
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <O'o|Vv>");
    global_dpd_->contract424(&N, &T_oo, &W, 1, 1, 1, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&N);

    // W_IjAb += -1/6 T_IK N_KjAb
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <Oo'|Vv>");
    global_dpd_->contract244(&T_OO, &N, &W, 1, 0, 0, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&N);

    global_dpd_->buf4_close(&W);

    // oovv
    // Precompute Temp_ikab = P_(ab) [M_ikab] - N_ikab
    // Temp_ijab = M_ijab
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "M <oo|v'v>");
    global_dpd_->buf4_copy(&M, PSIF_DCT_DPD, "Temp <oo|vv>");
    global_dpd_->buf4_close(&M);

    // Temp_ijab -> Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    // Temp_ijab -= Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&T, &TT, -1.0);
    global_dpd_->buf4_close(&TT);

    // Temp_ijab -= N_ijab
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "N <o'o|vv>");
    dpd_buf4_add(&T, &N, -1.0);
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_close(&T);

    // W_ijab += 1/6 P_(ij) [Temp_ikab T_jk]
    // Temp2_ijab = 1/6 Temp_ikab T_jk
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp2 <oo|vv>");
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->contract424(&T, &T_oo, &TT, 1, 1, 1, 1.0 / 6.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&TT);

    // W_ijab += Temp2_ijab
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp2 <oo|vv>");
    dpd_buf4_add(&W, &TT, 1.0);
    global_dpd_->buf4_close(&TT);

    // Temp2_ijab -> Temp2_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp2 <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    // W_ijab -= Temp2_jiab
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&W, &TT, -1.0);
    global_dpd_->buf4_close(&TT);
    global_dpd_->buf4_close(&W);

    /*
     * 4. W_ijab += 1/6 P_(ab) [(P_(ij) [N_ijac] - M_ijac) T_cb]
     */

    // OOVV
    // Precompute Temp_IJAC = P_(IJ) [N_IJAC] - M_IJAC
    // Temp_IJAB = N_IJAB
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "N <O'O|VV>");
    global_dpd_->buf4_copy(&N, PSIF_DCT_DPD, "Temp <OO|VV>");
    global_dpd_->buf4_close(&N);

    // Temp_IJAB -> Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -= Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&T, &TT, -1.0);
    global_dpd_->buf4_close(&TT);

    // Temp_IJAB -= M_IJAB
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "M <OO|V'V>");
    dpd_buf4_add(&T, &M, -1.0);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&T);

    // W_IJAB += 1/6 P_(AB) [Temp_IJAC T_CB]
    // Temp2_IJAB = 1/6 Temp_IJAC T_CB
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp2 <OO|VV>");
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->contract424(&T, &T_VV, &TT, 3, 1, 0, 1.0 / 6.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&TT);

    // W_IJAB += Temp2_IJAB
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp2 <OO|VV>");
    dpd_buf4_add(&W, &TT, 1.0);
    global_dpd_->buf4_close(&TT);

    // Temp2_IJAB -> Temp2_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp2 <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    // W_IJAB -= Temp2_IJBA
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&W, &TT, -1.0);
    global_dpd_->buf4_close(&TT);
    global_dpd_->buf4_close(&W);

    // OoVv
    // Temp_IjAc = N_IjAc - N_jIAc = N_IjAc + N_jIcA = (N_O'oVv + N_Oo'Vv)
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <O'o|Vv>");
    global_dpd_->buf4_copy(&N, PSIF_DCT_DPD, "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&N);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <Oo'|Vv>");
    dpd_buf4_add(&T, &N, 1.0);
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    // W_IjAb += 1/6 Temp_IjCb T_AC
    global_dpd_->contract244(&T_VV, &T, &W, 1, 2, 1, 1.0 / 6.0, 1.0);
    // W_IjAb += 1/6 Temp_IjAc T_bc
    global_dpd_->contract424(&T, &T_vv, &W, 3, 1, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&T);

    // W_IjAb += -1/6 T_AC M_IjCb
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv'>");
    global_dpd_->contract244(&T_VV, &M, &W, 1, 2, 1, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&M);

    // W_IjAb += -1/6 M_IjAc T_cb
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|V'v>");
    global_dpd_->contract424(&M, &T_vv, &W, 3, 1, 0, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&M);

    global_dpd_->buf4_close(&W);

    // oovv
    // Precompute Temp_ijac = P_(ij) [N_ijac] - M_ijac
    // Temp_ijab = N_ijab
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "N <o'o|vv>");
    global_dpd_->buf4_copy(&N, PSIF_DCT_DPD, "Temp <oo|vv>");
    global_dpd_->buf4_close(&N);

    // Temp_ijab -> Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    // Temp_ijab -= Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&T, &TT, -1.0);
    global_dpd_->buf4_close(&TT);

    // Temp_ijab -= M_ijab
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "M <oo|v'v>");
    dpd_buf4_add(&T, &M, -1.0);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&T);

    // W_ijab += 1/6 P_(ab) [Temp_ijac T_cb]
    // Temp2_ijab = 1/6 Temp_ijac T_cb
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp2 <oo|vv>");
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->contract424(&T, &T_vv, &TT, 3, 1, 0, 1.0 / 6.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&TT);

    // W_ijab += Temp2_ijab
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp2 <oo|vv>");
    dpd_buf4_add(&W, &TT, 1.0);
    global_dpd_->buf4_close(&TT);

    // Temp2_ijab -> Temp2_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp2 <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    // W_ijab -= Temp2_ijba
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&W, &TT, -1.0);
    global_dpd_->buf4_close(&TT);
    global_dpd_->buf4_close(&W);

    /*
     * 5. W_ijab += 1/12 P_(ab) [F_ca O_ijcb]
     */

    /*
     * 6. W_ijab -= 1/12 P_(ij) [F_ik O_kjab]
     */

    // OOVV
    // Temp_IJAB = 1/12 F_AC O_IJCB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "O <OO|VV>");
    global_dpd_->contract244(&F_VV, &O, &T, 1, 2, 1, 1.0 / 12.0, 0.0);
    global_dpd_->buf4_close(&O);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");

    // W_IJAB += Temp_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_add(&W, &T, 1.0);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    // W_IJAB -= Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB = -1/12 F_IK O_KJAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "O <OO|VV>");
    global_dpd_->contract244(&F_OO, &O, &T, 1, 0, 0, -1.0 / 12.0, 0.0);
    global_dpd_->buf4_close(&O);

    // W_IJAB += Temp_IJAB
    dpd_buf4_add(&W, &T, 1.0);

    // Temp_IJAB -> Temp_JIAB
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    // W_IJAB -= Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    // OoVv
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "O <Oo|Vv>");
    // W_IjAb += 1/12 F_AC O_IjCb
    global_dpd_->contract244(&F_VV, &O, &W, 1, 2, 1, 1.0 / 12.0, 1.0);
    // W_IjAb += 1/12 O_IjAc F_bc
    global_dpd_->contract424(&O, &F_vv, &W, 3, 1, 0, 1.0 / 12.0, 1.0);
    // W_IjAb -= 1/12 F_IK O_KjAb
    global_dpd_->contract244(&F_OO, &O, &W, 1, 0, 0, -1.0 / 12.0, 1.0);
    // W_IjAb -= 1/12 O_IkAb F_jk
    global_dpd_->contract424(&O, &F_oo, &W, 1, 1, 1, -1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&O);
    global_dpd_->buf4_close(&W);

    // oovv
    // Temp_ijab = 1/12 F_ac O_ijcb
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "O <oo|vv>");
    global_dpd_->contract244(&F_vv, &O, &T, 1, 2, 1, 1.0 / 12.0, 0.0);
    global_dpd_->buf4_close(&O);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");

    // W_ijab += Temp_ijab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_add(&W, &T, 1.0);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    // W_ijab -= Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&T);

    // Temp_ijab = -1/12 F_ik O_kjab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "O <oo|vv>");
    global_dpd_->contract244(&F_oo, &O, &T, 1, 0, 0, -1.0 / 12.0, 0.0);
    global_dpd_->buf4_close(&O);

    // W_ijab += Temp_ijab
    dpd_buf4_add(&W, &T, 1.0);

    // Temp_ijab -> Temp_jiab
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    // W_ijab -= Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&F_OO);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->file2_close(&F_VV);
    global_dpd_->file2_close(&F_vv);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    /*
     * 7. W_ijab -= 1/6 I_ijkl N_klab
     */

    // W_IJAB -= 1/6 I_IJKL N_KLAB
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>O]-"), ID("[O>O]-"), 0, "I <OO|OO>");
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "N <O'O|VV>");
    global_dpd_->contract444(&I, &N, &W, 0, 1, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");

    // W_IjAb -= 1/6 I_IjKl N_KlAb (O'oVv)
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <O'o|Vv>");
    global_dpd_->contract444(&I, &N, &W, 0, 1, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&N);

    // W_IjAb -= 1/6 I_IjLk N_LkAb (Oo'Vv)
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <Oo'|Vv>");
    global_dpd_->contract444(&I, &N, &W, 0, 1, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&N);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&W);

    // W_ijab -= 1/6 I_ijkl N_klab
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[o,o]"), ID("[o>o]-"), ID("[o>o]-"), 0, "I <oo|oo>");
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "N <o'o|vv>");
    global_dpd_->contract444(&I, &N, &W, 0, 1, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&W);

    /*
     * 8. W_ijab += 1/6 M_ijcd Amplitude_cdkl Amplitude_klab
     */

    // Precompute Temp_ijkl = M_ijcd Amplitude_cdkl
    // Temp_IJKL = M_IJCD Amplitude_CDKL
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           "Temp <OO|OO>");
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "M <OO|V'V>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract444(&M, &L, &T, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "Temp <Oo|Oo>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");

    // Temp_IjKl = M_IjCd Amplitude_CdKl (M_OoV'v)
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|V'v>");
    global_dpd_->contract444(&M, &L, &T, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&M);

    // Temp_IjKl += M_IjcD Amplitude_cDKl = M_IjDc Amplitude_DcKl (M_OoVv')
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv'>");
    global_dpd_->contract444(&M, &L, &T, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&M);

    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_ijkl = M_ijcd Amplitude_cdkl
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[o,o]"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           "Temp <oo|oo>");
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "M <oo|v'v>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract444(&M, &L, &T, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&T);

    // W_IJAB += 1/6 Temp_IJKL Amplitude_KLAB
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           "Temp <OO|OO>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract444(&T, &L, &W, 0, 1, 2.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    // W_IjAb += 2/6 Temp_IjKl Amplitude_KlAb
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "Temp <Oo|Oo>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract444(&T, &L, &W, 0, 1, 2.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    // W_ijab += 1/6 Temp_ijkl Amplitude_klab
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           "Temp <oo|oo>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract444(&T, &L, &W, 0, 1, 2.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    /*
     * 9. W_ijab += 1/3 P_(ij) P_(ab) [(M_ikac + M_kica - N_ikac - N_kica) K_jckb]
     */

    // OOVV presort
    // Temp_IKAC = M_IKAC'
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "M <OO|V'V>");
    global_dpd_->buf4_sort(&M, PSIF_DCT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    global_dpd_->buf4_close(&M);

    // N_I'KAC -> N_IK'AC
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "N <O'O|VV>");
    global_dpd_->buf4_sort(&N, PSIF_DCT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "Temp2 <OO|VV>");
    global_dpd_->buf4_close(&N);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");

    // Temp_IKAC += M_IKA'C
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "M <OO|V'V>");
    dpd_buf4_add(&T, &M, 1.0);
    global_dpd_->buf4_close(&M);

    // Temp_IKAC -= N_I'KAC
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "N <O'O|VV>");
    dpd_buf4_add(&T, &N, -1.0);
    global_dpd_->buf4_close(&N);

    // Temp_IKAC -= N_IK'AC
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp2 <OO|VV>");
    dpd_buf4_add(&T, &N, -1.0);
    global_dpd_->buf4_close(&N);

    global_dpd_->buf4_close(&T);

    // Temp_IKAC -> Temp_(IA|KC)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Temp (OV|OV)");
    global_dpd_->buf4_close(&T);

    // OoVv presort
    // Temp_IkAc = M_IkA'c
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|V'v>");
    global_dpd_->buf4_copy(&M, PSIF_DCT_DPD, "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&M);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");

    // Temp_IkAc += M_IkAc'
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv'>");
    dpd_buf4_add(&T, &M, 1.0);
    global_dpd_->buf4_close(&M);

    // Temp_IkAc -= N_I'kAc
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <O'o|Vv>");
    dpd_buf4_add(&T, &N, -1.0);
    global_dpd_->buf4_close(&N);

    // Temp_IkAc -= N_Ik'Ac
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <Oo'|Vv>");
    dpd_buf4_add(&T, &N, -1.0);
    global_dpd_->buf4_close(&N);

    global_dpd_->buf4_close(&T);

    // Temp_IkAc -> Temp_(IA|kc)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, prqs, ID("[O,V]"), ID("[o,v]"), "Temp (OV|ov)");
    global_dpd_->buf4_close(&T);

    // oovv presort
    // Temp_ikac = M_ikac'
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "M <oo|v'v>");
    global_dpd_->buf4_sort(&M, PSIF_DCT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    global_dpd_->buf4_close(&M);

    // N_i'kac -> N_ik'ac
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "N <o'o|vv>");
    global_dpd_->buf4_sort(&N, PSIF_DCT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "Temp2 <oo|vv>");
    global_dpd_->buf4_close(&N);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");

    // Temp_ikac += M_ika'c
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "M <oo|v'v>");
    dpd_buf4_add(&T, &M, 1.0);
    global_dpd_->buf4_close(&M);

    // Temp_ikac -= N_i'kac
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "N <o'o|vv>");
    dpd_buf4_add(&T, &N, -1.0);
    global_dpd_->buf4_close(&N);

    // Temp_ikac -= N_ik'ac
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp2 <oo|vv>");
    dpd_buf4_add(&T, &N, -1.0);
    global_dpd_->buf4_close(&N);

    global_dpd_->buf4_close(&T);

    // Temp_ikac -> Temp_(ia|kc)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, prqs, ID("[o,v]"), ID("[o,v]"), "Temp (ov|ov)");
    global_dpd_->buf4_close(&T);

    // Compute W contributions
    // OOVV
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Temp2 (OV|OV)");
    global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");

    // Temp2_(IA|JB) += 1/3 Temp_(IA|KC) K_(KC|JB)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
    global_dpd_->contract444(&T, &Kaa, &TT, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&T);

    // Temp2_(IA|JB) += 1/3 Temp_(IA|kc) K_(kc|JB)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    global_dpd_->contract444(&T, &Kab, &TT, 0, 0, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Kaa);

    // Temp2_(IA|JB) -> Temp2_<IJ|AB>
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp2 <OO|VV>");
    global_dpd_->buf4_close(&TT);

    // W_IJAB += Temp2_IJAB
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp2 <OO|VV>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    dpd_buf4_add(&W, &TT, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&TT);

    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp2 <OO|VV>");

    // Temp2_IJAB -> Temp2_IJBA
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // W_IJAB -= Temp2_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T);

    // Temp2_IJAB -> Temp2_JIAB
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // W_IJAB -= Temp2_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T);

    // Temp2_IJAB -> Temp2_JIBA
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // W_IJAB += Temp2_JIBA
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                           "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0, "W <OO|VV>");
    dpd_buf4_add(&W, &T, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&TT);

    // OoVv
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                           "Temp2 (OV|ov)");

    global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");
    global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");

    // Temp2_(IA|jb) = 1/3 Temp_(IA|KC) K_(KC|jb)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
    global_dpd_->contract444(&T, &Kab, &TT, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    // Temp2_(IA|jb) += 1/3 Temp_(IA|kc) K_(kc|jb)
    global_dpd_->contract444(&T, &Kbb, &TT, 0, 1, 1.0 / 3.0, 1.0);
    // Temp2_(IA|jb) += 1/3 K_(IA|KC) Temp_(KC|jb)
    global_dpd_->contract444(&Kaa, &T, &TT, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&T);

    // Temp2_(IA|jb) += 1/3 K_(IA|kc) Temp_(kc|jb)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "Temp (ov|ov)");
    global_dpd_->contract444(&Kab, &T, &TT, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Kbb);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Kaa);

    // Temp2_(IA|jb) -> Temp2_<Ij|Ab>
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "Temp2 <Oo|Vv>");
    global_dpd_->buf4_close(&TT);

    // W_IjAb += Temp2_IjAb
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Temp2 <Oo|Vv>");
    dpd_buf4_add(&W, &TT, 1.0);
    global_dpd_->buf4_close(&TT);
    global_dpd_->buf4_close(&W);

    // Temp_(IC|kb) -> Temp_(Ib|Ck)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, psqr, ID("[O,v]"), ID("[V,o]"), "Temp (Ov|Vo)");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                           "Temp2 (Ov|Vo)");
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0, "Temp (Ov|Vo)");

    // Temp2_(Ib|Aj) = 1/3 Temp_(Ib|Ck) K_<Ck|Aj>
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), 0, "K <Vo|Vo>");
    global_dpd_->contract444(&T, &Kab, &TT, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&Kab);

    // Temp2_(Ib|Aj) += 1/3 K_<Kc|Ib> Temp_(Kc|Aj)
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->contract444(&Kab, &T, &TT, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&Kab);

    global_dpd_->buf4_close(&T);

    // Temp2_(Ib|Aj) -> Temp2_<Ij|Ab>
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, psrq, ID("[O,o]"), ID("[V,v]"), "Temp2 <Oo|Vv>");
    global_dpd_->buf4_close(&TT);

    // W_IjAb += Temp2_<Ij|Ab>
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "W <Oo|Vv>");
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Temp2 <Oo|Vv>");
    dpd_buf4_add(&W, &TT, 1.0);
    global_dpd_->buf4_close(&TT);
    global_dpd_->buf4_close(&W);

    // oovv
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Temp2 (ov|ov)");
    global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");

    // Temp2_(ia|jb) += 1/3 Temp_(ia|kc) K_(kc|jb)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "Temp (ov|ov)");
    global_dpd_->contract444(&T, &Kbb, &TT, 0, 1, 1.0 / 3.0, 0.0);
    global_dpd_->buf4_close(&T);

    // Temp2_(ia|jb) += 1/3 Temp_(ia|KC) K_(KC|JB)
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    global_dpd_->contract444(&T, &Kab, &TT, 1, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Kbb);

    // Temp2_(ia|jb) -> Temp2_<ij|ab>
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "Temp2 <oo|vv>");
    global_dpd_->buf4_close(&TT);

    // W_ijab += Temp2_ijab
    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp2 <oo|vv>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    dpd_buf4_add(&W, &TT, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&TT);

    global_dpd_->buf4_init(&TT, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp2 <oo|vv>");

    // Temp2_ijab -> Temp2_ijba
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // W_ijab -= Temp2_ijba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T);

    // Temp2_ijab -> Temp2_jiab
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // W_ijab -= Temp2_jiab
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    dpd_buf4_add(&W, &T, -1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T);

    // Temp2_ijab -> Temp2_jiba
    global_dpd_->buf4_sort(&TT, PSIF_DCT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // W_ijab += Temp2_jiba
    global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                           "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0, "W <oo|vv>");
    dpd_buf4_add(&W, &T, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&TT);
}

void DCTSolver::compute_M_intermediate() {
    dpdbuf4 M, Laa, Lab, Lbb;
    dpdfile2 F_VV, F_vv;

    // M_IJAB = F_AC Amplitude_IJCB
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "M <OO|V'V>");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->file2_init(&F_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->contract244(&F_VV, &Laa, &M, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_VV);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&M);

    // M_IjAb = F_AC Amplitude_IjCb
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|V'v>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->file2_init(&F_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->contract244(&F_VV, &Lab, &M, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_VV);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&M);

    // M_JiBa = Amplitude_JiBc F_ac
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv'>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->file2_init(&F_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    global_dpd_->contract424(&Lab, &F_vv, &M, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F_vv);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&M);

    // M_ijab = F_ac lambda_ijcb
    global_dpd_->buf4_init(&M, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "M <oo|v'v>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->file2_init(&F_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    global_dpd_->contract244(&F_vv, &Lbb, &M, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_vv);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&M);
}

void DCTSolver::compute_N_intermediate() {
    dpdbuf4 N, Laa, Lab, Lbb;
    dpdfile2 F_OO, F_oo;

    // N_IJAB = F_IK lambda_KJAB
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "N <O'O|VV>");
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->file2_init(&F_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract244(&F_OO, &Laa, &N, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F_OO);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&N);

    // N_IjAb = F_IK lambda_KjAb
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <O'o|Vv>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->file2_init(&F_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract244(&F_OO, &Lab, &N, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F_OO);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&N);

    // N_JiBa = lambda_JkBa F_ik
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "N <Oo'|Vv>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->file2_init(&F_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->contract424(&Lab, &F_oo, &N, 1, 1, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&N);

    // N_ijab = F_ik lambda_kjab
    global_dpd_->buf4_init(&N, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "N <o'o|vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->file2_init(&F_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->contract244(&F_oo, &Lbb, &N, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&N);
}
}  // namespace dct
}  // namespace psi
