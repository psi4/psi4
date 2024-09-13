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

#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/molecule.h"
#include "psi4/psifiles.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <algorithm>
#include <functional>

namespace psi {
namespace dct {

void DCTSolver::compute_SO_tau_R() {
    build_d_R();
    if (exact_tau_) {
        build_tau_R();
    }
    transform_tau_R();
}

/**
 * Forms d in the MO basis from the Amplitude tensors.
 * Several variables used in this function are labeled Tau. At this time, they are NOT Tau, but d.
 */
void DCTSolver::build_d_R() {
    dct_timer_on("DCTSolver::build_d()");
    dpdbuf4 L1, L2;
    dpdfile2 T_OO, T_VV;

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

    global_dpd_->buf4_init(&L1, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&L2, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
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

    global_dpd_->buf4_init(&L1, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>
    global_dpd_->buf4_init(&L2, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>

    /*
     * d_IJ -= 1/2 Amplitude_IkAb Amplitude_JkAb - 1/2 Amplitude_IkaB Amplitude_JkaB
     */
    global_dpd_->contract442(&L1, &L2, &T_OO, 0, 0, -1.0, 1.0);
    /*
     * d_AB += 1/2 Amplitude_IjAc Amplitude_IjBc + 1/2 Amplitude_iJAc Amplitude_iJBc
     */
    global_dpd_->contract442(&L1, &L2, &T_VV, 2, 2, 1.0, 1.0);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->buf4_close(&L1);
    global_dpd_->buf4_close(&L2);

    // Read MO-basis Tau from disk into the memory
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

    aocc_tau_ = Matrix(&T_OO);
    avir_tau_ = Matrix(&T_VV);

    bocc_tau_.copy(aocc_tau_);
    bvir_tau_.copy(avir_tau_);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_VV);

    dct_timer_off("DCTSolver::build_d()");
}

void DCTSolver::build_tau_R() {
    dct_timer_on("DCTSolver::build_tau()");

    dpdfile2 T_OO, T_VV;

    // Iteratively compute the exact Tau

    auto aocc_tau_old = Matrix("MO basis Tau (Alpha Occupied, old)", naoccpi_, naoccpi_);
    auto avir_tau_old = Matrix("MO basis Tau (Alpha Virtual, old)", navirpi_, navirpi_);
    auto aocc_d = Matrix("Non-idempotency of OPDM (Alpha Occupied, old)", naoccpi_, naoccpi_);
    auto avir_d = Matrix("Non-idempotency of OPDM (Alpha Virtual, old)", navirpi_, navirpi_);

    bool converged = false;
    bool failed = false;
    int cycle = 0;

    // Copy approximate Tau as the non-idempotency of OPDM
    aocc_d.copy(aocc_tau_);
    avir_d.copy(avir_tau_);

    while (!converged && !failed) {
        // Save old tau from previous iteration
        aocc_tau_old.copy(aocc_tau_);
        avir_tau_old.copy(avir_tau_);

        // Tau_ij = d_ij
        // Tau_ab = -d_ab
        aocc_tau_.copy(aocc_d);
        avir_tau_.copy(avir_d);

        // Tau_ij -= Tau_ik * Tau_kj
        // Tau_ab += Tau_ac * Tau_cb
        aocc_tau_.gemm(false, false, -1.0, aocc_tau_old, aocc_tau_old, 1.0);
        avir_tau_.gemm(false, false, 1.0, avir_tau_old, avir_tau_old, 1.0);

        // Compute RMS
        aocc_tau_old.subtract(aocc_tau_);
        avir_tau_old.subtract(avir_tau_);

        double rms = aocc_tau_old.rms();
        rms += avir_tau_old.rms();
        rms *= 2.0;

        converged = (rms < cumulant_threshold_);
        failed = (++cycle == maxiter_);

        if (print_ > 2) outfile->Printf("\t Exact Tau Iterations: %-3d %20.12f\n", cycle, rms);
        //            if (print_ > 0) outfile->Printf( "\t Exact Tau Iterations: %-3d %20.12f\n", cycle, rms);

    }  // end of macroiterations

    // If exact tau iterations failed, throw a message about it and compute it non-iteratively
    if (failed) {
        outfile->Printf("\t Exact Tau didn't converge. Evaluating it non-iteratively\n");
        // Set old tau matrices to identity
        aocc_tau_old.identity();
        avir_tau_old.identity();
        // Scale the non-idempotency elements
        aocc_d.scale(4.0);
        avir_d.scale(-4.0);
        // Add them to the old tau
        aocc_tau_old.add(aocc_d);
        avir_tau_old.add(avir_d);
        // Zero out new tau
        aocc_tau_.zero();
        avir_tau_.zero();
        // Diagonalize and take a square root
        auto aocc_evecs = std::make_shared<Matrix>("Eigenvectors (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_);
        auto avir_evecs = std::make_shared<Matrix>("Eigenvectors (Alpha Virtual)", nirrep_, navirpi_, navirpi_);
        auto aocc_evals = std::make_shared<Vector>("Eigenvalues (Alpha Occupied)", naoccpi_);
        auto avir_evals = std::make_shared<Vector>("Eigenvalues (Alpha Virtual)", navirpi_);
        aocc_tau_old.diagonalize(aocc_evecs, aocc_evals);
        avir_tau_old.diagonalize(avir_evecs, avir_evals);

        for (int h = 0; h < nirrep_; ++h) {
            if (nsopi_[h] == 0) continue;

            // Alpha occupied
            for (int p = 0; p < naoccpi_[h]; ++p) aocc_tau_.set(h, p, p, (-1.0 + sqrt(aocc_evals->get(h, p))) / 2.0);

            // Alpha virtual
            for (int p = 0; p < navirpi_[h]; ++p) avir_tau_.set(h, p, p, (1.0 - sqrt(avir_evals->get(h, p))) / 2.0);
        }

        // Back-transform the diagonal Tau to the original basis
        aocc_tau_.back_transform(aocc_evecs);
        avir_tau_.back_transform(avir_evecs);
    }

    // Copy Tau_alpha to Tau_beta
    bocc_tau_.copy(aocc_tau_);
    bvir_tau_.copy(avir_tau_);

    // Write the exact tau back to disk

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

    aocc_tau_.write_to_dpdfile2(&T_OO);
    avir_tau_.write_to_dpdfile2(&T_VV);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_VV);

    dct_timer_off("DCTSolver::build_tau()");
}

void DCTSolver::transform_tau_R() {
    dct_timer_on("DCTSolver::transform_tau()");

    dpdfile2 T_OO, T_VV;

    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

    // Zero SO tau arrays before computing it in the MO basis
    tau_so_a_->zero();
    tau_so_a_->add(linalg::triplet(*aocc_c_, Matrix(&T_OO), *aocc_c_, false, false, true));
    tau_so_a_->add(linalg::triplet(*avir_c_, Matrix(&T_VV), *avir_c_, false, false, true));

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_VV);

    // Copy tau_so_alpha to tau_so_beta
    tau_so_b_->copy(tau_so_a_);

    dct_timer_off("DCTSolver::transform_tau()");
}

/**
 * Prints the occupation numbers from the OPDM
 */
void DCTSolver::print_opdm_RHF() {
    dpdbuf4 L1, L2;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

    global_dpd_->file2_mat_init(&T_OO);
    global_dpd_->file2_mat_init(&T_VV);
    global_dpd_->file2_mat_rd(&T_OO);
    global_dpd_->file2_mat_rd(&T_VV);

    std::vector<std::pair<double, int> > aPairs;

    for (int h = 0; h < nirrep_; ++h) {
        for (int row = 0; row < T_OO.params->coltot[h]; ++row)
            aPairs.push_back(std::make_pair(1.0 + T_OO.matrix[h][row][row], h));
        for (int row = 0; row < T_VV.params->coltot[h]; ++row)
            aPairs.push_back(std::make_pair(T_VV.matrix[h][row][row], h));
    }

    global_dpd_->file2_mat_close(&T_OO);
    global_dpd_->file2_mat_close(&T_VV);
    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_VV);

    sort(aPairs.begin(), aPairs.end(), std::greater<std::pair<double, int> >());

    auto aIrrepCount = std::vector<int>(nirrep_, 0);
    std::vector<std::string> irrepLabels = molecule_->irrep_labels();

    outfile->Printf("\n\tOrbital occupations:\n\t\tDoubly occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < nalpha_; ++i, ++count) {
        int irrep = aPairs[i].second;
        outfile->Printf("%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep].c_str(), 2 * aPairs[i].first);
        if (count % 4 == 3 && i != nalpha_) outfile->Printf("\n\t\t");
    }
    outfile->Printf("\n\n\t\tVirtual orbitals\n\t\t");
    for (int i = nalpha_, count = 0; i < nmo_; ++i, ++count) {
        int irrep = aPairs[i].second;
        outfile->Printf("%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep].c_str(), 2 * aPairs[i].first);
        if (count % 4 == 3 && i != nmo_) outfile->Printf("\n\t\t");
    }
    outfile->Printf("\n\n");
}
}  // namespace dct
}  // namespace psi
