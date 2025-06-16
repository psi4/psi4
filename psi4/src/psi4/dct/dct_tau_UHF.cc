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

#include "dct.h"
#include "psi4/psifiles.h"

#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/molecule.h"
#include "psi4/psifiles.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

#include <algorithm>
#include <functional>

namespace psi {
namespace dct {

void DCTSolver::compute_SO_tau_U() {
    build_d_U();
    if (exact_tau_) {
        build_tau_U();
    }
    transform_tau_U();
}

/**
 * Forms Tau in the MO basis from the Amplitude tensors.
 * Several variables used in this function are labeled Tau. At this time, they are NOT Tau, but d.
 */
void DCTSolver::build_d_U() {
    dct_timer_on("DCTSolver::build_d()");
    dpdbuf4 L1, L2;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    global_dpd_->buf4_init(&L1, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&L2, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    /*
     * d_IJ = -1/2 Amplitude_IKAB Amplitude_JKAB
     */
    global_dpd_->contract442(&L1, &L2, &T_OO, 0, 0, -0.5, 0.0);
    /*
     * d_AB = +1/2 Amplitude_IJAC Amplitude_IJBC
     */
    global_dpd_->contract442(&L1, &L2, &T_VV, 2, 2, 0.5, 0.0);
    global_dpd_->buf4_close(&L1);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L1, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&L2, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    /*
     * d_ij = -1/2 Amplitude_ikab Amplitude_jkab
     */
    global_dpd_->contract442(&L1, &L2, &T_oo, 0, 0, -0.5, 0.0);
    /*
     * d_ab = +1/2 Amplitude_ijac Amplitude_ijbc
     */
    global_dpd_->contract442(&L1, &L2, &T_vv, 2, 2, 0.5, 0.0);
    global_dpd_->buf4_close(&L1);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L1, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&L2, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    /*
     * d_IJ -= 1/2 Amplitude_IkAb Amplitude_JkAb - 1/2 Amplitude_IkaB Amplitude_JkaB
     */
    global_dpd_->contract442(&L1, &L2, &T_OO, 0, 0, -1.0, 1.0);
    /*
     * d_ij -= 1/2 Amplitude_KiAb Amplitude_KjAb - 1/2 Amplitude_KiaB Amplitude_KjaB
     */
    global_dpd_->contract442(&L1, &L2, &T_oo, 1, 1, -1.0, 1.0);
    /*
     * d_AB += 1/2 Amplitude_IjAc Amplitude_IjBc + 1/2 Amplitude_iJAc Amplitude_iJBc
     */
    global_dpd_->contract442(&L1, &L2, &T_VV, 2, 2, 1.0, 1.0);
    /*
     * d_ab += 1/2 Amplitude_IjCa Amplitude_IjCb + 1/2 Amplitude_iJCa Amplitude_iJCb
     */
    global_dpd_->contract442(&L1, &L2, &T_vv, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    // Compute fourth-order Tau terms if requested
    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        build_d_fourth_order_U();
    }

    // Read MO-basis Tau from disk into the memory
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    aocc_tau_ = Matrix(&T_OO);
    avir_tau_ = Matrix(&T_VV);
    bocc_tau_ = Matrix(&T_oo);
    bvir_tau_ = Matrix(&T_vv);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    dct_timer_off("DCTSolver::build_d()");
}

void DCTSolver::build_d_fourth_order_U() {
    dpdbuf4 I, K, II, KK, Temp, L, O;
    dpdfile2 T_OO, T_oo, T_VV, T_vv, TT_OO, TT_oo, TT_VV, TT_vv, Tau_OO, Tau_oo, Tau_VV, Tau_vv;

    // Prepare T intermediates
    // Copy Tau
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    global_dpd_->file2_copy(&T_OO, PSIF_DCT_DPD, "T <O|O>");
    global_dpd_->file2_copy(&T_oo, PSIF_DCT_DPD, "T <o|o>");
    global_dpd_->file2_copy(&T_VV, PSIF_DCT_DPD, "T <V|V>");
    global_dpd_->file2_copy(&T_vv, PSIF_DCT_DPD, "T <v|v>");

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    // Multiply by 2.0
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "T <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "T <v|v>");

    global_dpd_->file2_scm(&T_OO, 2.0);
    global_dpd_->file2_scm(&T_oo, 2.0);
    global_dpd_->file2_scm(&T_VV, 2.0);
    global_dpd_->file2_scm(&T_vv, 2.0);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    // Compute I, K, and O intermediates needed for the fourth-order d terms
    compute_I_intermediate();

    compute_K_intermediate();

    compute_O_intermediate();

    // Fourth-order Tau terms

    /*
     * d_ij += 1/3 T_ik T_kj
     */

    // Open Tau_OO
    global_dpd_->file2_init(&Tau_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&Tau_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");

    // Open T
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "T <o|o>");

    // 1. d_IJ += 1/3 T_IK T_KJ
    global_dpd_->file2_init(&TT_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
    global_dpd_->contract222(&T_OO, &TT_OO, &Tau_OO, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->file2_close(&TT_OO);

    // 1. d_ij += 1/3 T_ik T_kj
    global_dpd_->file2_init(&TT_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "T <o|o>");
    global_dpd_->contract222(&T_oo, &TT_oo, &Tau_oo, 0, 1, 1.0 / 3.0, 1.0);
    global_dpd_->file2_close(&TT_oo);

    /*
     * d_ij -= 1/12 I_ikjl T_lk
     */

    // 2. d_IJ -= 1/12 I_(IJ|KL) T_KL
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0, "I (OO|OO)");
    global_dpd_->contract422(&I, &T_OO, &Tau_OO, 0, 0, -1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&I);

    // 2. d_IJ -= 1/12 I_(IJ|kl) T_kl
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[o,o]"), ID("[O,O]"), ID("[o,o]"), 0, "I (OO|oo)");
    global_dpd_->contract422(&I, &T_oo, &Tau_OO, 0, 0, -1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&I);

    // 2. d_ij -= 1/12 I_(ij|KL) T_KL
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[O,O]"), ID("[o,o]"), ID("[O,O]"), 0, "I (oo|OO)");
    global_dpd_->contract422(&I, &T_OO, &Tau_oo, 0, 0, -1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&I);

    // 2. d_ij -= 1/12 I_(ij|kl) T_kl
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[o,o]"), ID("[o,o]"), ID("[o,o]"), 0, "I (oo|oo)");
    global_dpd_->contract422(&I, &T_oo, &Tau_oo, 0, 0, -1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&I);

    // Close T
    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);

    /*
     * d_ij += 1/6 K_icjd T_dc
     */

    // 3. d_IJ += 1/6 K_(IJ|CD) T_CD
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K (OO|VV)");
    global_dpd_->contract422(&K, &T_VV, &Tau_OO, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->file2_close(&T_VV);

    // 3. d_IJ += 1/6 K_(IJ|cd) T_cd
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "T <v|v>");
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[v,v]"), ID("[O,O]"), ID("[v,v]"), 0, "K (OO|vv)");
    global_dpd_->contract422(&K, &T_vv, &Tau_OO, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->file2_close(&T_vv);

    // 3. d_ij += 1/6 K_(ij|CD) T_CD
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[V,V]"), ID("[o,o]"), ID("[V,V]"), 0, "K (oo|VV)");
    global_dpd_->contract422(&K, &T_VV, &Tau_oo, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->file2_close(&T_VV);

    // 3. d_ij += 1/6 K_(ij|cd) T_cd
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "T <v|v>");
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0, "K (oo|vv)");
    global_dpd_->contract422(&K, &T_vv, &Tau_oo, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->file2_close(&T_vv);

    /*
     * d_ij -= 1/24 I_iklm I_lmjk
     */

    // 4. d_IJ -= 1/24 I<IK|LM> I<JK|LM>
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>O]-"), ID("[O>O]-"), 0, "I <OO|OO>");
    global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>O]-"), ID("[O>O]-"), 0, "I <OO|OO>");
    global_dpd_->contract442(&I, &II, &Tau_OO, 0, 0, -1.0 / 24.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    // 4. d_IJ -= 1/12 I<Ik|Lm> I<Jk|Lm>
    global_dpd_->contract442(&I, &II, &Tau_OO, 0, 0, -1.0 / 12.0, 1.0);
    // 4. d_ij -= 1/12 I<Ki|Lm> I<Kj|Lm>
    global_dpd_->contract442(&I, &II, &Tau_oo, 1, 1, -1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);

    // 4. d_ij -= 1/24 I<ik|lm> I<jk|lm>
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[o,o]"), ID("[o>o]-"), ID("[o>o]-"), 0, "I <oo|oo>");
    global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[o,o]"), ID("[o>o]-"), ID("[o>o]-"), 0, "I <oo|oo>");
    global_dpd_->contract442(&I, &II, &Tau_oo, 0, 0, -1.0 / 24.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);

    /*
     * d_ij -= 1/3 K_ickd I_kdjc
     */

    // 5. d_IJ -= 1/3 K<IC|KD> K<JC|KD>
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K <OV|OV>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K <OV|OV>");
    global_dpd_->contract442(&K, &KK, &Tau_OO, 0, 0, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    // 5. d_IJ -= 1/3 K<Ic|Kd> K<Jc|Kd>
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->contract442(&K, &KK, &Tau_OO, 0, 0, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0, "K <Ov|oV>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0, "K <Ov|oV>");
    // 5. d_IJ -= 1/3 K<Ic|kD> K<Jc|kD>
    global_dpd_->contract442(&K, &KK, &Tau_OO, 0, 0, -1.0 / 3.0, 1.0);
    // 5. d_ij -= 1/3 K<Kc|iD> K<Kc|jD>
    global_dpd_->contract442(&K, &KK, &Tau_oo, 2, 2, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    // 5. d_ij -= 1/3 K<iC|kD> K<jC|kD>
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->contract442(&K, &KK, &Tau_oo, 0, 0, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    // 5. d_ij -= 1/3 K<ic|kd> K<jc|kd>
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K <ov|ov>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K <ov|ov>");
    global_dpd_->contract442(&K, &KK, &Tau_oo, 0, 0, -1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    // Close Tau_OO
    global_dpd_->file2_close(&Tau_oo);
    global_dpd_->file2_close(&Tau_OO);

    /*
     * d_ab -= 1/3 T_ac T_cb
     */

    // Open Tau_VV
    global_dpd_->file2_init(&Tau_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&Tau_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    // Open T
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "T <v|v>");

    // 1. d_AB -= 1/3 T_AC T_CB
    global_dpd_->file2_init(&TT_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
    global_dpd_->contract222(&T_VV, &TT_VV, &Tau_VV, 0, 1, -1.0 / 3.0, 1.0);
    global_dpd_->file2_close(&TT_VV);

    // 1. d_ab -= 1/3 T_ac T_cb
    global_dpd_->file2_init(&TT_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "T <v|v>");
    global_dpd_->contract222(&T_vv, &TT_vv, &Tau_vv, 0, 1, -1.0 / 3.0, 1.0);
    global_dpd_->file2_close(&TT_vv);

    /*
     * d_ab -= 1/12 Amplitude_ackl Amplitude_klbd T_dc
     */

    // Precompute Temp_klbc = Amplitude_klbd T_dc

    // OOVV

    // Temp_IJAB = Amplitude_IJAC * T_CB

    global_dpd_->buf4_init(&Temp, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract424(&L, &T_VV, &Temp, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Temp);

    // OoVv

    // Temp_IjAb = Amplitude_IjAc * T_cb
    global_dpd_->buf4_init(&Temp, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract424(&L, &T_vv, &Temp, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Temp);

    // Temp2_JiBa = T_CB * Amplitude_JiCa
    global_dpd_->buf4_init(&Temp, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Temp2 <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract244(&T_VV, &L, &Temp, 1, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Temp);

    // oovv

    // Temp_ijab = Amplitude_ijac * T_cb
    global_dpd_->buf4_init(&Temp, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract424(&L, &T_vv, &Temp, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Temp);

    // Close T
    global_dpd_->file2_close(&T_vv);
    global_dpd_->file2_close(&T_VV);

    // Compute Tau contribution

    // 2. d_AB -= 1/12 Amplitude<KL|AC> Temp_<KL|BC>
    global_dpd_->buf4_init(&Temp, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Temp <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract442(&L, &Temp, &Tau_VV, 2, 2, -1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Temp);

    // 2. d_AB -= 1/6 Amplitude<Kl|Ac> Temp_<Kl|Bc>
    global_dpd_->buf4_init(&Temp, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract442(&L, &Temp, &Tau_VV, 2, 2, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Temp);

    // 2. d_ab -= 1/6 Amplitude<Kl|Ca> Temp2_<Kl|Cb>
    global_dpd_->buf4_init(&Temp, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Temp2 <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->contract442(&L, &Temp, &Tau_vv, 3, 3, -1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Temp);

    // 2. d_ab -= 1/12 Amplitude<kl|ac> Temp_<kl|bc>
    global_dpd_->buf4_init(&Temp, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                           "Temp <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract442(&L, &Temp, &Tau_vv, 2, 2, -1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&Temp);

    /*
     * d_ab += 1/6 K_kalb T_lk
     */

    // 3. d_AB += 1/6 K_(AB|KL) T_KL
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), 0, "K (VV|OO)");
    global_dpd_->contract422(&K, &T_OO, &Tau_VV, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->file2_close(&T_OO);

    // 3. d_AB += 1/6 K_(AB|kl) T_kl
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "T <o|o>");
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[V,V]"), ID("[o,o]"), ID("[V,V]"), ID("[o,o]"), 0, "K (VV|oo)");
    global_dpd_->contract422(&K, &T_oo, &Tau_VV, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->file2_close(&T_oo);

    // 3. d_ab += 1/6 K_(ab|KL) T_KL
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[v,v]"), ID("[O,O]"), ID("[v,v]"), ID("[O,O]"), 0, "K (vv|OO)");
    global_dpd_->contract422(&K, &T_OO, &Tau_vv, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->file2_close(&T_OO);

    // 3. d_ab += 1/6 K_(ab|kl) T_kl
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "T <o|o>");
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), 0, "K (vv|oo)");
    global_dpd_->contract422(&K, &T_oo, &Tau_vv, 0, 0, 1.0 / 6.0, 1.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->file2_close(&T_oo);

    /*
     * d_ab += 1/24 Amplitude_ackl O_klbc
     */

    // 4. d_AB += 1/24 Amplitude<KL|AC> O<KL|BC>
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0, "O <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract442(&L, &O, &Tau_VV, 2, 2, 1.0 / 24.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&O);

    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0, "O <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    // 4. d_AB += 1/12 Amplitude<Kl|Ac> O<Kl|Bc>
    global_dpd_->contract442(&L, &O, &Tau_VV, 2, 2, 1.0 / 12.0, 1.0);
    // 4. d_AB += 1/24 Amplitude<Lk|Ca> O<Lk|Cb>
    global_dpd_->contract442(&L, &O, &Tau_vv, 3, 3, 1.0 / 12.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&O);

    // 4. d_ab += 1/24 Amplitude<kl|ac> O<kl|bc>
    global_dpd_->buf4_init(&O, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0, "O <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract442(&L, &O, &Tau_vv, 2, 2, 1.0 / 24.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&O);

    /*
     * d_ab += 1/3 K_kalc K_lckb
     */

    // 5. d_AB += 1/3 K<KA|LC> K<KB|LC>
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K <OV|OV>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K <OV|OV>");
    global_dpd_->contract442(&K, &KK, &Tau_VV, 1, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    // 5. d_AB += 1/3 K<kA|lC> K<kB|lC>
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->contract442(&K, &KK, &Tau_VV, 1, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0, "K <Ov|oV>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0, "K <Ov|oV>");
    // 5. d_AB += 1/3 K<Lc|kA> K<Lc|kB>
    global_dpd_->contract442(&K, &KK, &Tau_VV, 3, 3, 1.0 / 3.0, 1.0);
    // 5. d_ab += 1/3 K<Ka|lC> K<Kb|lC>
    global_dpd_->contract442(&K, &KK, &Tau_vv, 1, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    // 5. d_ab += 1/3 K<Ka|Lc> K<Kb|Lc>
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->contract442(&K, &KK, &Tau_vv, 1, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    // 5. d_ab += 1/3 K<ka|lc> K<kb|lc>
    global_dpd_->buf4_init(&K, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K <ov|ov>");
    global_dpd_->buf4_init(&KK, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K <ov|ov>");
    global_dpd_->contract442(&K, &KK, &Tau_vv, 1, 1, 1.0 / 3.0, 1.0);
    global_dpd_->buf4_close(&KK);
    global_dpd_->buf4_close(&K);

    // Close Tau_VV
    global_dpd_->file2_close(&Tau_vv);
    global_dpd_->file2_close(&Tau_VV);
}

void DCTSolver::transform_tau_U() {
    dct_timer_on("DCTSolver::transform_tau()");

    dpdfile2 T_OO, T_oo, T_VV, T_vv;

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    // Zero SO tau arrays before computing it in the MO basis
    tau_so_a_->zero();
    tau_so_b_->zero();

    tau_so_a_->add(linalg::triplet(*aocc_c_, Matrix(&T_OO), *aocc_c_, false, false, true));
    tau_so_a_->add(linalg::triplet(*avir_c_, Matrix(&T_VV), *avir_c_, false, false, true));
    tau_so_b_->add(linalg::triplet(*bocc_c_, Matrix(&T_oo), *bocc_c_, false, false, true));
    tau_so_b_->add(linalg::triplet(*bvir_c_, Matrix(&T_vv), *bvir_c_, false, false, true));

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    dct_timer_off("DCTSolver::transform_tau()");
}

/**
 * Prints the occupation numbers from the OPDM
 */
void DCTSolver::print_opdm() {
    // Obtain natural occupancies
    // Form one-particle density matrix
    auto a_opdm = std::make_shared<Matrix>("MO basis OPDM (Alpha)", nirrep_, nmopi_, nmopi_);
    auto b_opdm = std::make_shared<Matrix>("MO basis OPDM (Beta)", nirrep_, nmopi_, nmopi_);

    // Alpha spin
    for (int h = 0; h < nirrep_; ++h) {
        // O-O
        for (int i = 0; i < naoccpi_[h]; ++i) {
            for (int j = 0; j <= i; ++j) {
                a_opdm->set(h, i, j, (aocc_tau_.get(h, i, j) + kappa_mo_a_->get(h, i, j)));
                if (i != j) a_opdm->set(h, j, i, (aocc_tau_.get(h, i, j) + kappa_mo_a_->get(h, i, j)));
            }
        }
        // V-V
        for (int a = 0; a < navirpi_[h]; ++a) {
            for (int b = 0; b <= a; ++b) {
                a_opdm->set(h, a + naoccpi_[h], b + naoccpi_[h], avir_tau_.get(h, a, b));
                if (a != b) a_opdm->set(h, b + naoccpi_[h], a + naoccpi_[h], avir_tau_.get(h, a, b));
            }
        }
    }

    // Beta spin
    for (int h = 0; h < nirrep_; ++h) {
        // O-O
        for (int i = 0; i < nboccpi_[h]; ++i) {
            for (int j = 0; j <= i; ++j) {
                b_opdm->set(h, i, j, (bocc_tau_.get(h, i, j) + kappa_mo_b_->get(h, i, j)));
                if (i != j) b_opdm->set(h, j, i, (bocc_tau_.get(h, i, j) + kappa_mo_b_->get(h, i, j)));
            }
        }
        // V-V
        for (int a = 0; a < nbvirpi_[h]; ++a) {
            for (int b = 0; b <= a; ++b) {
                b_opdm->set(h, a + nboccpi_[h], b + nboccpi_[h], bvir_tau_.get(h, a, b));
                if (a != b) b_opdm->set(h, b + nboccpi_[h], a + nboccpi_[h], bvir_tau_.get(h, a, b));
            }
        }
    }

    // Diagonalize OPDM to obtain NOs
    auto aevecs = std::make_shared<Matrix>("Eigenvectors (Alpha)", nirrep_, nmopi_, nmopi_);
    auto bevecs = std::make_shared<Matrix>("Eigenvectors (Beta)", nirrep_, nmopi_, nmopi_);
    auto aevals = std::make_shared<Vector>("Eigenvalues (Alpha)", nmopi_);
    auto bevals = std::make_shared<Vector>("Eigenvalues (Beta)", nmopi_);

    a_opdm->diagonalize(aevecs, aevals, descending);
    b_opdm->diagonalize(bevecs, bevals, descending);

    std::vector<std::pair<double, int> > aPairs;
    std::vector<std::pair<double, int> > bPairs;

    for (int h = 0; h < nirrep_; ++h) {
        for (int p = 0; p < nmopi_[h]; ++p) {
            aPairs.push_back(std::make_pair(aevals->get(h, p), h));
            bPairs.push_back(std::make_pair(bevals->get(h, p), h));
        }
    }

    sort(aPairs.begin(), aPairs.end(), std::greater<std::pair<double, int> >());
    sort(bPairs.begin(), bPairs.end(), std::greater<std::pair<double, int> >());

    auto aIrrepCount = std::vector<int>(nirrep_, 0);
    auto bIrrepCount = std::vector<int>(nirrep_, 0);
    std::vector<std::string> irrepLabels = molecule_->irrep_labels();

    outfile->Printf("\n\tOrbital occupations:\n\t\tAlpha occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < nalpha_; ++i, ++count) {
        int irrep = aPairs[i].second;
        outfile->Printf("%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep].c_str(), aPairs[i].first);
        if (count % 4 == 3 && i != nalpha_) outfile->Printf("\n\t\t");
    }
    outfile->Printf("\n\n\t\tBeta occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < nbeta_; ++i, ++count) {
        int irrep = bPairs[i].second;
        outfile->Printf("%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep].c_str(), bPairs[i].first);
        if (count % 4 == 3 && i != nbeta_) outfile->Printf("\n\t\t");
    }
    outfile->Printf("\n\n\t\tAlpha virtual orbitals\n\t\t");
    for (int i = nalpha_, count = 0; i < nmo_; ++i, ++count) {
        int irrep = aPairs[i].second;
        outfile->Printf("%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep].c_str(), aPairs[i].first);
        if (count % 4 == 3 && i != nmo_) outfile->Printf("\n\t\t");
    }
    outfile->Printf("\n\n\t\tBeta virtual orbitals\n\t\t");
    for (int i = nbeta_, count = 0; i < nmo_; ++i, ++count) {
        int irrep = bPairs[i].second;
        outfile->Printf("%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep].c_str(), bPairs[i].first);
        if (count % 4 == 3 && i != nmo_) outfile->Printf("\n\t\t");
    }
    outfile->Printf("\n\n");
}

void DCTSolver::build_tau_U() {
    // Read MO-basis Tau from disk into the memory
    dpdfile2 T_OO, T_oo, T_VV, T_vv;

    // Iteratively compute the exact Tau

    auto aocc_tau_old = std::make_shared<Matrix>("MO basis Tau (Alpha Occupied, old)", nirrep_, naoccpi_, naoccpi_);
    auto bocc_tau_old = std::make_shared<Matrix>("MO basis Tau (Beta Occupied, old)", nirrep_, nboccpi_, nboccpi_);
    auto avir_tau_old = std::make_shared<Matrix>("MO basis Tau (Alpha Virtual, old)", nirrep_, navirpi_, navirpi_);
    auto bvir_tau_old = std::make_shared<Matrix>("MO basis Tau (Beta Virtual, old)", nirrep_, nbvirpi_, nbvirpi_);
    auto aocc_d =
        std::make_shared<Matrix>("Non-idempotency of OPDM (Alpha Occupied, old)", nirrep_, naoccpi_, naoccpi_);
    auto bocc_d = std::make_shared<Matrix>("Non-idempotency of OPDM (Beta Occupied, old)", nirrep_, nboccpi_, nboccpi_);
    auto avir_d = std::make_shared<Matrix>("Non-idempotency of OPDM (Alpha Virtual, old)", nirrep_, navirpi_, navirpi_);
    auto bvir_d = std::make_shared<Matrix>("Non-idempotency of OPDM (Beta Virtual, old)", nirrep_, nbvirpi_, nbvirpi_);

    bool converged = false;
    bool failed = false;
    int cycle = 0;

    // Copy approximate Tau as the non-idempotency of OPDM
    aocc_d->copy(aocc_tau_);
    avir_d->copy(avir_tau_);
    bocc_d->copy(bocc_tau_);
    bvir_d->copy(bvir_tau_);

    DIISManager diisManager(maxdiis_, "DCT DIIS Tau", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);
    if ((nalpha_ + nbeta_) > 1) {
        diisManager.set_error_vector_size(&aocc_tau_, &bocc_tau_, &avir_tau_, &bvir_tau_);
        diisManager.set_vector_size(&aocc_tau_, &bocc_tau_, &avir_tau_, &bvir_tau_);
    }

    auto aocc_r = std::make_shared<Matrix>("Residual (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_);
    auto bocc_r = std::make_shared<Matrix>("Residual (Beta Occupied)", nirrep_, nboccpi_, nboccpi_);
    auto avir_r = std::make_shared<Matrix>("Residual (Alpha Virtual)", nirrep_, navirpi_, navirpi_);
    auto bvir_r = std::make_shared<Matrix>("Residual (Beta Virtual)", nirrep_, nbvirpi_, nbvirpi_);

    while (!converged && !failed) {
        std::string diisString;

        // Save old tau from previous iteration
        aocc_tau_old->copy(aocc_tau_);
        avir_tau_old->copy(avir_tau_);
        bocc_tau_old->copy(bocc_tau_);
        bvir_tau_old->copy(bvir_tau_);

        // R_ij = d_ij
        // R_ab = -d_ab
        aocc_r->copy(aocc_d);
        bocc_r->copy(bocc_d);
        avir_r->copy(avir_d);
        bvir_r->copy(bvir_d);

        // R_ij -= Tau_ik * Tau_kj
        // R_ab += Tau_ac * Tau_cb
        aocc_r->gemm(false, false, -1.0, aocc_tau_old, aocc_tau_old, 1.0);
        bocc_r->gemm(false, false, -1.0, bocc_tau_old, bocc_tau_old, 1.0);
        avir_r->gemm(false, false, 1.0, avir_tau_old, avir_tau_old, 1.0);
        bvir_r->gemm(false, false, 1.0, bvir_tau_old, bvir_tau_old, 1.0);

        // R_ij -= Tau_ij
        // R_ab -= Tau_ab
        aocc_r->subtract(aocc_tau_old);
        bocc_r->subtract(bocc_tau_old);
        avir_r->subtract(avir_tau_old);
        bvir_r->subtract(bvir_tau_old);

        // Compute RMS
        double rms = aocc_r->rms();
        rms += bocc_r->rms();
        rms += avir_r->rms();
        rms += bvir_r->rms();

        converged = (rms < cumulant_threshold_);
        failed = (++cycle == maxiter_);

        aocc_tau_.add(aocc_r);
        bocc_tau_.add(bocc_r);
        avir_tau_.add(avir_r);
        bvir_tau_.add(bvir_r);

        if (rms < diis_start_thresh_) {
            // Store the DIIS vectors
            if (diisManager.add_entry(aocc_r.get(), bocc_r.get(), avir_r.get(), bvir_r.get(), &aocc_tau_, &bocc_tau_, &avir_tau_, &bvir_tau_)) {
                diisString += "S";
            }

            if (diisManager.subspace_size() > mindiisvecs_) {
                diisString += "/E";
                diisManager.extrapolate(&aocc_tau_, &bocc_tau_, &avir_tau_, &bvir_tau_);
            }
        }

        if (print_ > 1) outfile->Printf("\t Exact Tau Iterations: %-3d %20.12f %-3s\n", cycle, rms, diisString.c_str());

    }  // end of macroiterations

    // If exact tau iterations failed, throw a message about it and compute it non-iteratively
    if (failed) {
        outfile->Printf("\t Exact Tau didn't converge. Evaluating it non-iteratively\n");
        // Set old tau matrices to identity
        aocc_tau_old->identity();
        bocc_tau_old->identity();
        avir_tau_old->identity();
        bvir_tau_old->identity();
        // Scale the non-idempotency elements
        aocc_d->scale(4.0);
        bocc_d->scale(4.0);
        avir_d->scale(-4.0);
        bvir_d->scale(-4.0);
        // Add them to the old tau
        aocc_tau_old->add(aocc_d);
        bocc_tau_old->add(bocc_d);
        avir_tau_old->add(avir_d);
        bvir_tau_old->add(bvir_d);
        // Zero out new tau
        aocc_tau_.zero();
        avir_tau_.zero();
        bocc_tau_.zero();
        bvir_tau_.zero();
        // Diagonalize and take a square root
        auto aocc_evecs = std::make_shared<Matrix>("Eigenvectors (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_);
        auto bocc_evecs = std::make_shared<Matrix>("Eigenvectors (Beta Occupied)", nirrep_, nboccpi_, nboccpi_);
        auto avir_evecs = std::make_shared<Matrix>("Eigenvectors (Alpha Virtual)", nirrep_, navirpi_, navirpi_);
        auto bvir_evecs = std::make_shared<Matrix>("Eigenvectors (Beta Virtual)", nirrep_, nbvirpi_, nbvirpi_);
        auto aocc_evals = std::make_shared<Vector>("Eigenvalues (Alpha Occupied)", naoccpi_);
        auto bocc_evals = std::make_shared<Vector>("Eigenvalues (Beta Occupied)", nboccpi_);
        auto avir_evals = std::make_shared<Vector>("Eigenvalues (Alpha Virtual)", navirpi_);
        auto bvir_evals = std::make_shared<Vector>("Eigenvalues (Beta Virtual)", nbvirpi_);
        aocc_tau_old->diagonalize(aocc_evecs, aocc_evals);
        bocc_tau_old->diagonalize(bocc_evecs, bocc_evals);
        avir_tau_old->diagonalize(avir_evecs, avir_evals);
        bvir_tau_old->diagonalize(bvir_evecs, bvir_evals);

        for (int h = 0; h < nirrep_; ++h) {
            if (nsopi_[h] == 0) continue;

            // Alpha occupied
            for (int p = 0; p < naoccpi_[h]; ++p) aocc_tau_.set(h, p, p, (-1.0 + sqrt(aocc_evals->get(h, p))) / 2.0);

            // Beta occupied
            for (int p = 0; p < nboccpi_[h]; ++p) bocc_tau_.set(h, p, p, (-1.0 + sqrt(bocc_evals->get(h, p))) / 2.0);

            // Alpha virtual
            for (int p = 0; p < navirpi_[h]; ++p) avir_tau_.set(h, p, p, (1.0 - sqrt(avir_evals->get(h, p))) / 2.0);

            // Beta virtual
            for (int p = 0; p < nbvirpi_[h]; ++p) bvir_tau_.set(h, p, p, (1.0 - sqrt(bvir_evals->get(h, p))) / 2.0);
        }

        // Back-transform the diagonal Tau to the original basis
        aocc_tau_.back_transform(aocc_evecs);
        bocc_tau_.back_transform(bocc_evecs);
        avir_tau_.back_transform(avir_evecs);
        bvir_tau_.back_transform(bvir_evecs);
    }

    // Write the exact tau back to disk

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    aocc_tau_.write_to_dpdfile2(&T_OO);
    avir_tau_.write_to_dpdfile2(&T_VV);
    bocc_tau_.write_to_dpdfile2(&T_oo);
    bvir_tau_.write_to_dpdfile2(&T_vv);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);
}

}  // namespace dct
}  // namespace psi
