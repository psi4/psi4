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

#include "dlpno.h"
#include "sparse.h"

#include "psi4/lib3index/3index.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/orthog.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libqt/qt.h"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dlpno {

DLPNOMP2::DLPNOMP2(SharedWavefunction ref_wfn, Options& options) : DLPNO(ref_wfn, options) {}
DLPNOMP2::~DLPNOMP2() {}

void DLPNOMP2::compute_pno_overlaps() {
    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    S_pno_ij_kj_.resize(n_lmo_pairs);
    S_pno_ij_ik_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        if (n_pno_[ij] == 0) continue;
        if (i < j) continue;

        S_pno_ij_kj_[ij] = std::vector<SharedMatrix>(naocc);
        S_pno_ij_ik_[ij] = std::vector<SharedMatrix>(naocc);

        for (int k = 0; k < naocc; ++k) {
            int kj = i_j_to_ij_[k][j];
            if (kj != -1 && i != k && fabs(F_lmo_->get(i, k)) > options_.get_double("F_CUT") && n_pno_[kj] > 0) {
                S_pno_ij_kj_[ij][k] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[kj]);
                S_pno_ij_kj_[ij][k] = linalg::triplet(X_pno_[ij], S_pno_ij_kj_[ij][k], X_pno_[kj], true, false, false);
            }

            int ik = i_j_to_ij_[i][k];
            if (ik != -1 && j != k && fabs(F_lmo_->get(k, j)) > options_.get_double("F_CUT") && n_pno_[ik] > 0) {
                S_pno_ij_ik_[ij][k] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ik]);
                S_pno_ij_ik_[ij][k] = linalg::triplet(X_pno_[ij], S_pno_ij_ik_[ij][k], X_pno_[ik], true, false, false);
            }
        }

        if (i > j) {
            S_pno_ij_kj_[ji] = S_pno_ij_ik_[ij];
            S_pno_ij_ik_[ji] = S_pno_ij_kj_[ij];
        }
    }
}

void DLPNOMP2::lmp2_iterations() {
    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    outfile->Printf("\n  ==> Local MP2 <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                     Corr. Energy    Delta E     Max R\n");

    std::vector<SharedMatrix> R_iajb(n_lmo_pairs);

    int iteration = 0, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, r_curr = 0.0;
    bool e_converged = false, r_converged = false;
    DIISManager diis(options_.get_int("DIIS_MAX_VECS"), "LMP2 DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    while (!(e_converged && r_converged)) {
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);

        // Calculate residuals from current amplitudes
#pragma omp parallel for schedule(static, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];

            R_iajb[ij] = std::make_shared<Matrix>("Residual", n_pno_[ij], n_pno_[ij]);

            if (n_pno_[ij] == 0) continue;

            for (int a = 0; a < n_pno_[ij]; ++a) {
                for (int b = 0; b < n_pno_[ij]; ++b) {
                    R_iajb[ij]->set(a, b,
                                    K_iajb_[ij]->get(a, b) +
                                        (e_pno_[ij]->get(a) + e_pno_[ij]->get(b) - F_lmo_->get(i, i) - F_lmo_->get(j, j)) *
                                            T_iajb_[ij]->get(a, b));
                }
            }

            for (int k = 0; k < naocc; ++k) {
                int kj = i_j_to_ij_[k][j];
                int ik = i_j_to_ij_[i][k];

                if (kj != -1 && i != k && fabs(F_lmo_->get(i, k)) > options_.get_double("F_CUT") && n_pno_[kj] > 0) {
                    auto temp =
                        linalg::triplet(S_pno_ij_kj_[ij][k], T_iajb_[kj], S_pno_ij_kj_[ij][k], false, false, true);
                    temp->scale(-1.0 * F_lmo_->get(i, k));
                    R_iajb[ij]->add(temp);
                }
                if (ik != -1 && j != k && fabs(F_lmo_->get(k, j)) > options_.get_double("F_CUT") && n_pno_[ik] > 0) {
                    auto temp =
                        linalg::triplet(S_pno_ij_ik_[ij][k], T_iajb_[ik], S_pno_ij_ik_[ij][k], false, false, true);
                    temp->scale(-1.0 * F_lmo_->get(k, j));
                    R_iajb[ij]->add(temp);
                }
            }

            R_iajb_rms[ij] = R_iajb[ij]->rms();
        }

        // evaluate convergence using current amplitudes and residuals
        e_prev = e_curr;
        e_curr = compute_iteration_energy(R_iajb);
        r_curr = *max_element(R_iajb_rms.begin(), R_iajb_rms.end());

        r_converged = (fabs(r_curr) < options_.get_double("R_CONVERGENCE"));
        e_converged = (fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE"));

// use residuals to get next amplitudes
#pragma omp parallel for schedule(static, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            for (int a = 0; a < n_pno_[ij]; ++a) {
                for (int b = 0; b < n_pno_[ij]; ++b) {
                    T_iajb_[ij]->add(a, b, -R_iajb[ij]->get(a, b) / ((e_pno_[ij]->get(a) + e_pno_[ij]->get(b)) -
                                                                    (F_lmo_->get(i, i) + F_lmo_->get(j, j))));
                }
            }
        }

        // DIIS extrapolation
        auto T_iajb_flat = flatten_mats(T_iajb_);
        auto R_iajb_flat = flatten_mats(R_iajb);

        if (iteration == 0) {
            diis.set_error_vector_size(R_iajb_flat.get());
            diis.set_vector_size(T_iajb_flat.get());
        }

        diis.add_entry(R_iajb_flat.get(), T_iajb_flat.get());
        diis.extrapolate(T_iajb_flat.get());

        copy_flat_mats(T_iajb_flat, T_iajb_);

#pragma omp parallel for schedule(static, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            Tt_iajb_[ij]->copy(T_iajb_[ij]);
            Tt_iajb_[ij]->scale(2.0);
            Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());
        }

        outfile->Printf("  @LMP2 iter %3d: %16.12f %10.3e %10.3e\n", iteration, e_curr, e_curr - e_prev, r_curr);

        iteration++;

        if(iteration > max_iteration) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    e_lmp2_ = e_curr;

    e_lmp2_os_ = 0.0;
    for (int ij = 0; ij < T_iajb_.size(); ++ij) {
        e_lmp2_os_ += K_iajb_[ij]->vector_dot(T_iajb_[ij]);
    }
    e_lmp2_ss_ = e_curr - e_lmp2_os_;

}

double DLPNOMP2::compute_iteration_energy(const std::vector<SharedMatrix> &R_iajb) {
    double e_mp2 = 0.0;
    for (int ij = 0; ij < Tt_iajb_.size(); ++ij) {
        e_mp2 += K_iajb_[ij]->vector_dot(Tt_iajb_[ij]);
        e_mp2 += R_iajb[ij]->vector_dot(Tt_iajb_[ij]);
    }
    return e_mp2;
}

double DLPNOMP2::compute_energy() {
    timer_on("DLPNO-MP2");

    print_header();

    timer_on("Setup Orbitals");
    setup_orbitals();
    timer_off("Setup Orbitals");

    timer_on("Overlap Ints");
    compute_overlap_ints();
    timer_off("Overlap Ints");

    timer_on("Dipole Ints");
    compute_dipole_ints();
    timer_off("Dipole Ints");

    timer_on("Sparsity");
    prep_sparsity(true, false);
    timer_off("Sparsity");

    timer_on("DF Ints");
    print_integral_sparsity();
    compute_metric();
    compute_qia();
    timer_off("DF Ints");

    timer_on("PNO Transform");
    pno_transform();
    timer_off("PNO Transform");

    timer_on("PNO Overlaps");
    compute_pno_overlaps();
    timer_off("PNO Overlaps");

    timer_on("LMP2");
    lmp2_iterations();
    timer_off("LMP2");

    print_results();

    timer_off("DLPNO-MP2");

    double e_scf = reference_wavefunction_->energy();
    double e_mp2_corr = e_lmp2_ + de_dipole_ + de_pno_total_;
    double e_mp2_corr_os = e_lmp2_os_ + 0.5 * de_dipole_ + de_pno_total_os_;
    double e_mp2_corr_ss = e_lmp2_ss_ + 0.5 * de_dipole_ + de_pno_total_ss_;
    double e_mp2_total = e_scf + e_mp2_corr;

    double sss = options_.get_double("MP2_SS_SCALE");
    double oss = options_.get_double("MP2_OS_SCALE");

    set_scalar_variable("MP2 CORRELATION ENERGY", e_mp2_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_mp2_corr);
    set_scalar_variable("MP2 TOTAL ENERGY", e_mp2_total);
    set_scalar_variable("CURRENT ENERGY", e_mp2_total);

    set_scalar_variable("MP2 SINGLES ENERGY", 0.0);
    set_scalar_variable("MP2 DOUBLES ENERGY", e_mp2_corr);
    set_scalar_variable("MP2 SAME-SPIN CORRELATION ENERGY", e_mp2_corr_ss);
    set_scalar_variable("MP2 OPPOSITE-SPIN CORRELATION ENERGY", e_mp2_corr_os);
    set_scalar_variable("SCS-MP2 CORRELATION ENERGY", 6.0/5.0 * e_mp2_corr_os +
                                                      1.0/3.0 * e_mp2_corr_ss);
    set_scalar_variable("SCS-MP2 TOTAL ENERGY", e_scf +
                                                6.0/5.0 * e_mp2_corr_os +
                                                1.0/3.0 * e_mp2_corr_ss);
    set_scalar_variable("CUSTOM SCS-MP2 CORRELATION ENERGY", oss * e_mp2_corr_os +
                                                             sss * e_mp2_corr_ss);
    set_scalar_variable("CUSTOM SCS-MP2 TOTAL ENERGY", e_scf +
                                                       oss * e_mp2_corr_os +
                                                       sss * e_mp2_corr_ss);

    return e_mp2_total;
}

void DLPNOMP2::print_header() {
    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                     DLPNO-MP2                 \n");
    outfile->Printf("                   by Zach Glick               \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Detailed DLPNO thresholds and cutoffs:\n");
    outfile->Printf("    T_CUT_DO     = %6.3e \n", T_CUT_DO_);
    outfile->Printf("    T_CUT_PNO    = %6.3e \n", T_CUT_PNO_);
    outfile->Printf("    T_CUT_DO_IJ  = %6.3e \n", options_.get_double("T_CUT_DO_IJ"));
    outfile->Printf("    T_CUT_PRE    = %6.3e \n", options_.get_double("T_CUT_PRE"));
    outfile->Printf("    T_CUT_DO_PRE = %6.3e \n", options_.get_double("T_CUT_DO_PRE"));
    outfile->Printf("    T_CUT_MKN    = %6.3e \n", options_.get_double("T_CUT_MKN"));
    outfile->Printf("    T_CUT_CLMO   = %6.3e \n", options_.get_double("T_CUT_CLMO"));
    outfile->Printf("    T_CUT_CPAO   = %6.3e \n", options_.get_double("T_CUT_CPAO"));
    outfile->Printf("    S_CUT        = %6.3e \n", options_.get_double("S_CUT"));
    outfile->Printf("    F_CUT        = %6.3e \n", options_.get_double("F_CUT"));
    outfile->Printf("\n");
}

void DLPNOMP2::print_results() {
    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-MP2 Correlation Energy: %16.12f \n", e_lmp2_ + de_pno_total_ + de_dipole_);
    outfile->Printf("    MP2 Correlation Energy:           %16.12f \n", e_lmp2_);
    outfile->Printf("    LMO Truncation Correction:        %16.12f \n", de_dipole_);
    outfile->Printf("    PNO Truncation Correction:        %16.12f \n", de_pno_total_);
}

void DLPNOMP2::print_integral_sparsity() {
    // statistics for number of (MN|K) shell triplets we need to compute

    int nbf = basisset_->nbf();
    int nshell = basisset_->nshell();
    int naux = ribasis_->nbf();
    int naocc = nalpha_ - nfrzc();

    size_t triplets = 0;       // computed (MN|K) triplets with no screening
    size_t triplets_lmo = 0;   // computed (MN|K) triplets with only LMO screening
    size_t triplets_pao = 0;   // computed (MN|K) triplets with only PAO screening
    size_t triplets_both = 0;  // computed (MN|K) triplets with LMO and PAO screening

    for (size_t atom = 0; atom < riatom_to_shells1_.size(); atom++) {
        size_t nshellri_atom = atom_to_rishell_[atom].size();
        triplets += nshell * nshell * nshellri_atom;
        triplets_lmo += riatom_to_shells1_[atom].size() * nshell * nshellri_atom;
        triplets_pao += nshell * riatom_to_shells2_[atom].size() * nshellri_atom;
        triplets_both += riatom_to_shells1_[atom].size() * riatom_to_shells2_[atom].size() * nshellri_atom;
    }
    size_t screened_total = triplets - triplets_both;
    size_t screened_lmo = triplets - triplets_lmo;
    size_t screened_pao = triplets - triplets_pao;

    // statistics for the number of (iu|Q) integrals we're left with after the transformation

    size_t total_integrals = (size_t)naocc * nbf * naux;
    size_t actual_integrals = 0;

    for (size_t atom = 0; atom < riatom_to_shells1_.size(); atom++) {
        actual_integrals +=
            riatom_to_lmos_ext_[atom].size() * riatom_to_paos_ext_[atom].size() * atom_to_ribf_[atom].size();
    }

    // number of doubles * (2^3 bytes / double) * (1 GiB / 2^30 bytes)
    double total_memory = total_integrals * pow(2.0, -27);
    double actual_memory = actual_integrals * pow(2.0, -27);
    double screened_memory = total_memory - actual_memory;

    outfile->Printf("\n");
    outfile->Printf("    Coefficient sparsity in AO -> LMO transform: %6.2f %% \n", screened_lmo * 100.0 / triplets);
    outfile->Printf("    Coefficient sparsity in AO -> PAO transform: %6.2f %% \n", screened_pao * 100.0 / triplets);
    outfile->Printf("    Coefficient sparsity in combined transforms: %6.2f %% \n", screened_total * 100.0 / triplets);
    outfile->Printf("\n");
    outfile->Printf("    Storing transformed LMO/PAO integrals in sparse format.\n");
    outfile->Printf("    Required memory: %.3f GiB (%.2f %% reduction from dense format) \n", actual_memory,
                    screened_memory * 100.0 / total_memory);
}

}  // namespace dlpno
}  // namespace psi
