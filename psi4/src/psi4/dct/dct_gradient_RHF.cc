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

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.h"

namespace psi {
namespace dct {

void DCTSolver::compute_gradient_RHF() {
    bool is_df = options_.get_str("DCT_TYPE") == "DF";

    // Compute the VVVV block of the relaxed TPDM
    // TODO: The VVVV density requires V^4 memory, which we'd rather avoid. Implement a lower memory algorithm.
    // The obvious one is to assemble aVVV "slices" for all a in compute_lagrangian_VV, which requires V^3
    // memory. That will be slower, but may be worth it in some cases.
    compute_unrelaxed_density_VVVV_RHF(is_df);
    outfile->Printf("\t Computing energy-weighted density matrix from one- and two-particle densities...\n");
    if (is_df) {
        three_idx_separable_density();
        three_idx_cumulant_density_RHF();
        construct_metric_density("Reference");
        construct_metric_density("Correlation");
    }
    // Compute the OO block of MO Lagrangian
    compute_lagrangian_OO_RHF(is_df);
    // Compute the VV block of MO Lagrangian
    compute_lagrangian_VV_RHF(is_df);
    // Compute the energy-weighted density matrix
    compute_ewdm_odc_RHF();
}

void DCTSolver::compute_lagrangian_OO_RHF(bool separate_gbargamma) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, pT;

    // X_OO: One-electron contributions

    // X_KI = H_JK Tau_JI
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&pT, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_mat_init(&X);
    global_dpd_->file2_mat_init(&H);
    global_dpd_->file2_mat_init(&pT);
    global_dpd_->file2_mat_rd(&H);
    global_dpd_->file2_mat_rd(&pT);
    for (int h = 0; h < nirrep_; ++h) {
#pragma omp parallel for
        for (int i = 0; i < naoccpi_[h]; ++i) {
            for (int k = 0; k < naoccpi_[h]; ++k) {
                double value = 0.0;
                for (int j = 0; j < naoccpi_[h]; ++j) {
                    value += H.matrix[h][j][k] * (pT.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                }
                X.matrix[h][k][i] = value;
            }
        }
    }
    global_dpd_->file2_mat_wrt(&X);
    global_dpd_->file2_mat_close(&pT);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_mat_close(&X);
    global_dpd_->file2_close(&pT);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_close(&X);

    // X_OO: Two-electron contributions

    // If we have gbar_gamma separate, just like when computing the OO gradient, we compute these terms special.
    if (separate_gbargamma) {
        auto zero = Dimension(nirrep_);
        auto gbar_alpha_block = mo_gbarGamma_A_.get_block(Slice(zero, nalphapi_), Slice(zero, nsopi_));
        auto gamma_alpha_block = mo_gammaA_.get_block(Slice(zero, nsopi_), Slice(zero, nalphapi_));
        auto alpha_jk = linalg::doublet(gbar_alpha_block, gamma_alpha_block, false, false);

        global_dpd_->file2_init(&H, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X JK <O|O>");
        alpha_jk->write_to_dpdfile2(&H);
        global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
        global_dpd_->file2_axpy(&H, &X, 1.0, 0);
        global_dpd_->file2_close(&X);
        global_dpd_->file2_close(&H);
    }

    const std::string density_variable = separate_gbargamma ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string &x) { return (density_variable + x); };

    //
    // 2 * <OO||OO> Г_OOOO
    //

    // X_MI += 2 * <MJ||KL> Г_IJKL
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 1,
                           "MO Ints <OO|OO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                           varname("<OO|OO>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_MI += 4 * <Mj|Kl> Г_IjKl
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                           "MO Ints <OO|OO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                           varname("SF <OO|OO>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <OO||VV> Г_OOVV
    //

    // X_MI += <JM||BC> Г_JIBC
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           varname("<OO|VV>"));

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_MI += <Mj|Bc> Г_IjBc
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           varname("SF <OO|VV>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <OV||OV> Г_OVOV
    //

    // X_MI += <JB||MC> Г_JBIC
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints <OV|OV> - <OV|VO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("<OV|OV>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_MI += <Jb|Mc> Г_JbIc
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints <OV|OV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("SF <OV|OV>:<Ov|Ov>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_MI -= <Ma|jB> Г_IajB
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints SF <OV|OV>:<Ov|oV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("SF <OV|OV>:<Ov|oV>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    psio_->close(PSIF_DCT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void DCTSolver::compute_lagrangian_VV_RHF(bool separate_gbargamma) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, pT;

    // X_VV: One-electron contributions

    // X_CA = H_CB Tau_BA
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&pT, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

    global_dpd_->contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&pT);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_close(&X);

    // X_VV: Two-electron contributions

    if (separate_gbargamma) {
        auto zero = Dimension(nirrep_);
        auto gbar_alpha_block = mo_gbarGamma_A_.get_block(Slice(nalphapi_, nsopi_), Slice(zero, nsopi_));
        auto gamma_alpha_block = mo_gammaA_.get_block(Slice(zero, nsopi_), Slice(nalphapi_, nsopi_));
        auto alpha_jk = linalg::doublet(gbar_alpha_block, gamma_alpha_block, false, false);

        global_dpd_->file2_init(&H, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X JK <V|V>");
        alpha_jk->write_to_dpdfile2(&H);
        global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
        global_dpd_->file2_axpy(&H, &X, 1.0, 0);
        global_dpd_->file2_close(&X);
        global_dpd_->file2_close(&H);
    }

    const std::string density_variable = separate_gbargamma ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string &x) { return (density_variable + x); };

    //
    // 2 * <VV||VV> Г_VVVV
    //

    // X_EA += 2 * <EB||CD> Г_ABCD
    dct_timer_on("DCTSolver::2 * g_EBCD Gamma_ABCD");
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 1,
                           "MO Ints <VV|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                           varname("<VV|VV>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);
    dct_timer_off("DCTSolver::2 * g_EBCD Gamma_ABCD");

    // X_EA += 4 * <Eb|Cd> Г_AbCd
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                           "MO Ints <VV|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                           varname("SF <VV|VV>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <OO||VV> Г_OOVV
    //

    // X_CA += <BC||JK> Г_BAJK
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           varname("<OO|VV>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_CA += 2 * <Jk|Cb> Г_JkAb
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           varname("SF <OO|VV>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <OO||OV> Г_OVOV
    //

    // X_CA += <JB||KC> Г_JBKA
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints <OV|OV> - <OV|VO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("<OV|OV>"));

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_CA += <kB|jC> Г_kBjA
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints <OV|OV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("SF <OV|OV>:<Ov|Ov>"));  // Г<oV|oV> = Г <Ov|Ov>

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_CA -= <Kb|jC> Г_KbjA
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints SF <OV|OV>:<Ov|oV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("SF <OV|OV>:<Ov|oV>"));

    global_dpd_->contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    psio_->close(PSIF_DCT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void DCTSolver::compute_ewdm_odc_RHF() {
    dpdfile2 X_OV, X_VO, X_OO, X_VV;

    Matrix aW("Energy-weighted density matrix (Alpha)", nmopi_, nmopi_);

    auto a_opdm = Matrix("MO basis OPDM (Alpha)", nmopi_, nmopi_);

    const int *alpha_corr_to_pitzer = _ints->alpha_corr_to_pitzer();
    auto *alpha_pitzer_to_corr = new int[nmo_];
    ::memset(alpha_pitzer_to_corr, '\0', nmo_ * sizeof(int));

    for (int n = 0; n < nmo_; ++n) {
        alpha_pitzer_to_corr[alpha_corr_to_pitzer[n]] = n;
    }

    const int *beta_corr_to_pitzer = _ints->beta_corr_to_pitzer();
    auto *beta_pitzer_to_corr = new int[nmo_];
    ::memset(beta_pitzer_to_corr, '\0', nmo_ * sizeof(int));

    for (int n = 0; n < nmo_; ++n) {
        beta_pitzer_to_corr[beta_corr_to_pitzer[n]] = n;
    }

    // Alpha spin
    global_dpd_->file2_init(&X_OV, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->file2_init(&X_VO, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->file2_init(&X_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    global_dpd_->file2_init(&X_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "X <V|V>");

    global_dpd_->file2_mat_init(&X_OV);
    global_dpd_->file2_mat_init(&X_VO);
    global_dpd_->file2_mat_init(&X_OO);
    global_dpd_->file2_mat_init(&X_VV);

    global_dpd_->file2_mat_rd(&X_OV);
    global_dpd_->file2_mat_rd(&X_VO);
    global_dpd_->file2_mat_rd(&X_OO);
    global_dpd_->file2_mat_rd(&X_VV);

    for (int h = 0; h < nirrep_; ++h) {
// O-V and V-O
#pragma omp parallel for
        for (int i = 0; i < naoccpi_[h]; ++i) {
            for (int a = 0; a < navirpi_[h]; ++a) {
                double value = -0.5 * (X_OV.matrix[h][i][a] + X_VO.matrix[h][a][i]);
                aW.set(h, i, a + naoccpi_[h], value);
                aW.set(h, a + naoccpi_[h], i, value);
            }
        }
// O-O
#pragma omp parallel for
        for (int i = 0; i < naoccpi_[h]; ++i) {
            for (int j = 0; j <= i; ++j) {
                double value = -0.5 * (X_OO.matrix[h][i][j] + X_OO.matrix[h][j][i]);
                aW.set(h, i, j, value);
                aW.set(h, j, i, value);
                a_opdm.set(h, i, j, (aocc_tau_.get(h, i, j) + kappa_mo_a_->get(h, i, j)));
                if (i != j) a_opdm.set(h, j, i, (aocc_tau_.get(h, i, j) + kappa_mo_a_->get(h, i, j)));
            }
        }
// V-V
#pragma omp parallel for
        for (int a = 0; a < navirpi_[h]; ++a) {
            for (int b = 0; b <= a; ++b) {
                double value = -0.5 * (X_VV.matrix[h][a][b] + X_VV.matrix[h][b][a]);
                aW.set(h, a + naoccpi_[h], b + naoccpi_[h], value);
                aW.set(h, b + naoccpi_[h], a + naoccpi_[h], value);
                a_opdm.set(h, a + naoccpi_[h], b + naoccpi_[h], avir_tau_.get(h, a, b));
                if (a != b) a_opdm.set(h, b + naoccpi_[h], a + naoccpi_[h], avir_tau_.get(h, a, b));
            }
        }
    }
    global_dpd_->file2_mat_close(&X_VO);
    global_dpd_->file2_mat_close(&X_OV);
    global_dpd_->file2_mat_close(&X_OO);
    global_dpd_->file2_mat_close(&X_VV);
    global_dpd_->file2_close(&X_VO);
    global_dpd_->file2_close(&X_OV);
    global_dpd_->file2_close(&X_OO);
    global_dpd_->file2_close(&X_VV);

    /* Save the SO EWDM to the wavefunction as Lagrangian_.
     * All other EWDM processing operations are then redundant, but it's what deriv.cc:compute expects...
     */
    Lagrangian_ = std::make_shared<Matrix>("Lagrangian matrix", nirrep_, nsopi_, nsopi_);
    Lagrangian_->back_transform(aW, *Ca_);
    Lagrangian_->scale(-2.0);  // Doubling is equivalent to adding beta spin case

    if (options_.get_str("DCT_TYPE") == "DF") {
        return;
    }

    auto *aocc_qt = new int[nalpha_];
    auto *bocc_qt = new int[nbeta_];
    auto *avir_qt = new int[navir_];
    auto *bvir_qt = new int[nbvir_];

    int aocc_count = 0;
    int bocc_count = 0;
    int avir_count = 0;
    int bvir_count = 0;
    int offset = 0;

    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < naoccpi_[h]; ++i) {
            int pitzer = offset + i;
            aocc_qt[aocc_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = 0; i < nboccpi_[h]; ++i) {
            int pitzer = offset + i;
            bocc_qt[bocc_count++] = beta_pitzer_to_corr[pitzer];
        }
        for (int i = naoccpi_[h]; i < (nmopi_[h] - frzvpi_[h]); ++i) {
            int pitzer = offset + i;
            avir_qt[avir_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = nboccpi_[h]; i < (nmopi_[h] - frzvpi_[h]); ++i) {
            int pitzer = offset + i;
            bvir_qt[bvir_count++] = beta_pitzer_to_corr[pitzer];
        }
        offset += nmopi_[h];
    }

    dpdbuf4 G;

    struct iwlbuf AA, AB;

    // Dump the density to IWL

    // Initialize IWL buffer for AB case
    iwl_buf_init(&AB, PSIF_MO_TPDM, 1.0E-15, 0, 0);

    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);
    // VVVV
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                           "Gamma SF <VV|VV>");  // Gamma <Vv|Vv>
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for (size_t ab = 0; ab < G.params->rowtot[h]; ++ab) {
            int a = G.params->roworb[h][ab][0];
            int b = G.params->roworb[h][ab][1];
            int A = avir_qt[a];
            int B = bvir_qt[b];
            for (size_t cd = 0; cd < G.params->coltot[h]; ++cd) {
                int c = G.params->colorb[h][cd][0];
                int d = G.params->colorb[h][cd][1];
                int C = avir_qt[c];
                int D = bvir_qt[d];
                double value = 4.0 * G.matrix[h][ab][cd];
                iwl_buf_wrt_val(&AB, A, C, B, D, value, 0, "outfile", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // OOOO
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                           "Gamma SF <OO|OO>");

    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for (size_t ij = 0; ij < G.params->rowtot[h]; ++ij) {
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for (size_t kl = 0; kl < G.params->coltot[h]; ++kl) {
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
                int K = aocc_qt[k];
                int L = bocc_qt[l];
                double value = 4.0 * G.matrix[h][ij][kl];
                iwl_buf_wrt_val(&AB, I, K, J, L, value, 0, "outfile", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // OOVV
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Gamma SF <OO|VV>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for (size_t ij = 0; ij < G.params->rowtot[h]; ++ij) {
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for (size_t ab = 0; ab < G.params->coltot[h]; ++ab) {
                int a = G.params->colorb[h][ab][0];
                int b = G.params->colorb[h][ab][1];
                int A = avir_qt[a];
                int B = bvir_qt[b];
                double value = 4.0 * G.matrix[h][ij][ab];
                iwl_buf_wrt_val(&AB, I, A, J, B, value, 0, "outfile", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // OVOV
    // Г<OV|OV> must be antisymmetrized before contracting it with TEI derivatives:
    // Г<OV|OV> contribution

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Gamma SF <OV|OV>:<Ov|Ov>");
    global_dpd_->buf4_dump(&G, &AB, aocc_qt, bvir_qt, aocc_qt, bvir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Gamma SF <OV|OV>:<Ov|Ov>");  // Gamma <oV|oV>
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for (size_t ia = 0; ia < G.params->rowtot[h]; ++ia) {
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = avir_qt[a];
            for (size_t jb = 0; jb < G.params->coltot[h]; ++jb) {
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = avir_qt[b];
                double value = 1.0 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AB, A, B, I, J, value, 0, "outfile", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // Resort Г<Ia|jB> -> (-1.0) * Г(IB|ja) and dump it
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Gamma SF <OV|OV>:<Ov|oV>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for (size_t ia = 0; ia < G.params->rowtot[h]; ++ia) {
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = bvir_qt[a];
            for (size_t jb = 0; jb < G.params->coltot[h]; ++jb) {
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = avir_qt[b];
                double value = (-2.0) * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AB, I, B, J, A, value, 0, "outfile", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    iwl_buf_flush(&AB, 1);
    iwl_buf_close(&AB, 1);

    // Presort TPDM (MO, AB case) for Deriv use
    presort_mo_tpdm_AB();

    // Initialize IWL buffer for AA case
    iwl_buf_init(&AA, PSIF_MO_TPDM, 1.0E-16, 0, 0);
    // VVVV
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                           "Gamma <VV|VV>");
    global_dpd_->buf4_scm(&G, 2.0);  // Prefactor was 1.0. Use 2.0 because BB case = AA case
    global_dpd_->buf4_dump(&G, &AA, avir_qt, avir_qt, avir_qt, avir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // OOOO
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                           "Gamma <OO|OO>");
    global_dpd_->buf4_scm(&G, 2.0);  // Prefactor was 1.0. Use 2.0 because BB case = AA case
    global_dpd_->buf4_dump(&G, &AA, aocc_qt, aocc_qt, aocc_qt, aocc_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // OOVV
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Gamma <OO|VV>");
    global_dpd_->buf4_scm(&G, 2.0);  // Prefactor was 1.0. Use 2.0 because BB case = AA case
    global_dpd_->buf4_dump(&G, &AA, aocc_qt, aocc_qt, avir_qt, avir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // OVOV
    // Г<OV|OV> must be antisymmetrized before contracting it with TEI derivatives:
    // Г<OV|OV> contribution

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Gamma <OV|OV>");
    global_dpd_->buf4_scm(&G, 1.0);
    global_dpd_->buf4_dump(&G, &AA, aocc_qt, avir_qt, aocc_qt, avir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // -Г<ov|vo> contribution:
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Gamma <OV|OV>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for (size_t ia = 0; ia < G.params->rowtot[h]; ++ia) {
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = avir_qt[a];
            for (size_t jb = 0; jb < G.params->coltot[h]; ++jb) {
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = aocc_qt[j];
                int B = avir_qt[b];
                double value = -1.0 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AA, I, B, A, J, value, 0, "outfile", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    iwl_buf_flush(&AA, 1);
    iwl_buf_close(&AA, 1);

    // Presort TPDM (MO) for Deriv use
    presort_mo_tpdm_AA();

    delete[] alpha_pitzer_to_corr;
    delete[] beta_pitzer_to_corr;
    delete[] aocc_qt;
    delete[] bocc_qt;
    delete[] avir_qt;
    delete[] bvir_qt;

    psio_->close(PSIF_DCT_DENSITY, 1);
}

}  // namespace dct
}  // namespace psi
