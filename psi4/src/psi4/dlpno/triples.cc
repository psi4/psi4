/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include <ctime>
#include <algorithm>
#include <array>
#include <unordered_map>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dlpno {

#ifdef USING_Einsums
using einsums::All;
using einsums::Indices;
using einsums::Tensor;
using einsums::TensorView;
using einsums::tensor_algebra::einsum;
using einsums::tensor_algebra::permute;
namespace index = einsums::index;
namespace linear_algebra = einsums::linear_algebra;
#endif

namespace {

double iterative_tno_cutoff(Options& options, const char* absolute_option,
                            const char* legacy_scale_option) {
    if (options[absolute_option].has_changed()) {
        return options.get_double(absolute_option);
    }
    return options.get_double("T_CUT_TNO") * options.get_double(legacy_scale_option);
}

}  // namespace

DLPNOCCSD_T::DLPNOCCSD_T(SharedWavefunction ref_wfn, Options &options) : DLPNOCCSD(ref_wfn, options) {}
DLPNOCCSD_T::~DLPNOCCSD_T() {}

void DLPNOCCSD_T::print_header() {
    const bool full_triples_follow = algorithm_ != DLPNOMethod::CCSD_T;
    const bool t0_only = !full_triples_follow && options_.get_bool("T0_APPROXIMATION");
    std::string triples_algorithm = (t0_only) ? "SEMICANONICAL (T0)" : "ITERATIVE (T)";
    double t_cut_tno = options_.get_double("T_CUT_TNO");
    const double t_cut_tno_full = options_.get_double("T_CUT_TNO_FULL");
    const double t_cut_tno_strong =
        full_triples_follow
            ? t_cut_tno_full
            : iterative_tno_cutoff(options_, "T_CUT_TNO_STRONG", "T_CUT_TNO_STRONG_SCALE");
    const double t_cut_tno_weak =
        full_triples_follow
            ? t_cut_tno_full
            : iterative_tno_cutoff(options_, "T_CUT_TNO_WEAK", "T_CUT_TNO_WEAK_SCALE");

    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                    DLPNO-CCSD(T)              \n");
    outfile->Printf("                    by Andy Jiang              \n");
    outfile->Printf("               DOI: 10.1063/5.0219963          \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Detailed DLPNO thresholds and cutoffs:\n");
    outfile->Printf("    ALGORITHM    = %6s   \n", triples_algorithm.c_str());
    outfile->Printf("    T_CUT_TNO (T0)             = %6.3e \n", t_cut_tno);
    outfile->Printf("    T_CUT_DO_TRIPLES (T0)      = %6.3e \n", options_.get_double("T_CUT_DO_TRIPLES"));
    outfile->Printf("    T_CUT_MKN_TRIPLES (T0)     = %6.3e \n", options_.get_double("T_CUT_MKN_TRIPLES"));
    outfile->Printf("    T_CUT_TRIPLES_WEAK (T0)    = %6.3e \n", options_.get_double("T_CUT_TRIPLES_WEAK"));
    outfile->Printf("    T_CUT_TNO_PRE (T0)         = %6.3e \n", options_.get_double("T_CUT_TNO_PRE"));
    outfile->Printf("    T_CUT_DO_TRIPLES_PRE (T0)  = %6.3e \n", options_.get_double("T_CUT_DO_TRIPLES_PRE"));
    outfile->Printf("    T_CUT_MKN_TRIPLES_PRE (T0) = %6.3e \n", options_.get_double("T_CUT_MKN_TRIPLES_PRE"));
    if (!t0_only) {
        outfile->Printf("    T_CUT_TNO_STRONG (T)       = %6.3e \n", t_cut_tno_strong);
        outfile->Printf("    T_CUT_TNO_WEAK (T)         = %6.3e \n", t_cut_tno_weak);
        outfile->Printf("    F_CUT_T (T)                = %6.3e \n", options_.get_double("F_CUT_T"));
        outfile->Printf("    T_CUT_ITER (T)             = %6.3e \n", options_.get_double("T_CUT_ITER"));
    }
    outfile->Printf("    MIN_TNOS                   = %6d   \n", options_.get_int("MIN_TNOS"));
    outfile->Printf("    TRIPLES_MAX_WEAK_PAIRS     = %6d   \n", options_.get_int("TRIPLES_MAX_WEAK_PAIRS"));
    if (full_triples_follow &&
        (options_["T_CUT_TNO_STRONG"].has_changed() || options_["T_CUT_TNO_WEAK"].has_changed() ||
         options_["T_CUT_TNO_STRONG_SCALE"].has_changed() || options_["T_CUT_TNO_WEAK_SCALE"].has_changed())) {
        outfile->Printf(
            "\n    WARNING: Iterative-(T) strong/weak TNO cutoff options apply only when DLPNO-CCSD(T)\n"
            "             is the final method. Full triples were requested, so both iterative-(T)\n"
            "             cutoffs are overridden by T_CUT_TNO_FULL = %6.3e.\n",
            t_cut_tno_full);
    }
    outfile->Printf("\n\n");
}

SharedMatrix DLPNOCCSD_T::matmul_3d(SharedMatrix A, SharedMatrix X, int dim_old, int dim_new) {
    /*
    Performs the operation A'[i,j,k] = A[I,J,K] * X[i,I] * X[j,J] * X[k,K] for cube 3d tensors
    */

    SharedMatrix A_new = linalg::doublet(X, A, false, false);
    A_new->reshape(dim_new * dim_old, dim_old);
    A_new = linalg::doublet(A_new, X, false, true);

    SharedMatrix A_T = std::make_shared<Matrix>(dim_new * dim_new, dim_old);
    for (int ind = 0; ind < dim_new * dim_new * dim_old; ++ind) {
        int a = ind / (dim_new * dim_old), b = (ind / dim_old) % dim_new, c = ind % dim_old;
        (*A_T)(a *dim_new + b, c) = (*A_new)(a * dim_old + c, b);
    }
    A_T = linalg::doublet(A_T, X, false, true);

    A_new = std::make_shared<Matrix>(dim_new, dim_new * dim_new);

    for (int ind = 0; ind < dim_new * dim_new * dim_new; ++ind) {
        int a = ind / (dim_new * dim_new), b = (ind / dim_new) % dim_new, c = ind % dim_new;
        (*A_new)(a, b *dim_new + c) = (*A_T)(a * dim_new + c, b);
    }

    return A_new;
}

void DLPNOCCSD_T::triples_sparsity(bool prescreening) {
    /* 
    In the prescreening step, this generates the initial list of triplets from
    strong and weak pairs ijk from strong and weak pairs ij, jk, and ik 
    (weak pairs set by TRIPLES_MAX_WEAK_PAIRS).

    In the second prescreening step, screened triplets are determined and removed
    from consideration from the rest of the computation with their energy accounted
    for by de_lccsd_t_screened_
    */

    timer_on("Triples Sparsity");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();
    int npao = C_pao_->ncol();

    int MAX_WEAK_PAIRS = options_.get_int("TRIPLES_MAX_WEAK_PAIRS");

    if (prescreening) {
        int ijk = 0;
        for (int ij = 0; ij < n_lmo_pairs; ij++) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            if (i > j) continue;
            for (int k : lmopair_to_lmos_[ij]) {
                if (i > k || j > k) continue;
                if (i == j && j == k) continue;
                int ij_weak = i_j_to_ij_weak_[i][j], ik_weak = i_j_to_ij_weak_[i][k], kj_weak = i_j_to_ij_weak_[k][j];

                int weak_pair_count = 0;
                if (ij_weak != -1) weak_pair_count += 1;
                if (ik_weak != -1) weak_pair_count += 1;
                if (kj_weak != -1) weak_pair_count += 1;

                if (weak_pair_count > MAX_WEAK_PAIRS) continue;

                ijk_to_i_j_k_.push_back(std::make_tuple(i, j, k));
                i_j_k_to_ijk_[triplet_key(i, j, k, naocc)] = ijk;
                i_j_k_to_ijk_[triplet_key(i, k, j, naocc)] = ijk;
                i_j_k_to_ijk_[triplet_key(j, i, k, naocc)] = ijk;
                i_j_k_to_ijk_[triplet_key(j, k, i, naocc)] = ijk;
                i_j_k_to_ijk_[triplet_key(k, i, j, naocc)] = ijk;
                i_j_k_to_ijk_[triplet_key(k, j, i, naocc)] = ijk;
                ++ijk;
            }
        }
    } else {
        std::unordered_map<size_t, int> i_j_k_to_ijk_new;
        std::vector<std::tuple<int, int, int>> ijk_to_i_j_k_new;

        double t_cut_triples_weak = options_.get_double("T_CUT_TRIPLES_WEAK");
        de_lccsd_t_screened_ = 0.0;

        int ijk_new = 0;
        for (int ijk = 0; ijk < ijk_to_i_j_k_.size(); ++ijk) {
            int i, j, k;
            std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

            if (std::fabs(e_ijk_[ijk]) >= t_cut_triples_weak) {
                ijk_to_i_j_k_new.push_back(std::make_tuple(i, j, k));
                i_j_k_to_ijk_new[triplet_key(i, j, k, naocc)] = ijk_new;
                i_j_k_to_ijk_new[triplet_key(i, k, j, naocc)] = ijk_new;
                i_j_k_to_ijk_new[triplet_key(j, i, k, naocc)] = ijk_new;
                i_j_k_to_ijk_new[triplet_key(j, k, i, naocc)] = ijk_new;
                i_j_k_to_ijk_new[triplet_key(k, i, j, naocc)] = ijk_new;
                i_j_k_to_ijk_new[triplet_key(k, j, i, naocc)] = ijk_new;
                ++ijk_new;
            } else {
                de_lccsd_t_screened_ += e_ijk_[ijk];
            }
        }
        i_j_k_to_ijk_ = i_j_k_to_ijk_new;
        ijk_to_i_j_k_ = ijk_to_i_j_k_new;
    }

    int n_lmo_triplets = ijk_to_i_j_k_.size();
    int natom = molecule_->natom();
    int nbf = basisset_->nbf();

    // => Local density fitting domains <= //

    SparseMap lmo_to_ribfs(naocc);
    SparseMap lmo_to_riatoms(naocc);

    double t_cut_mkn_triples = (prescreening) ? options_.get_double("T_CUT_MKN_TRIPLES_PRE") : options_.get_double("T_CUT_MKN_TRIPLES");

    for (size_t i = 0; i < naocc; ++i) {
        // atomic mulliken populations for this orbital
        std::vector<double> mkn_pop(natom, 0.0);

        auto P_i = reference_wavefunction_->S()->clone();

        for (size_t u = 0; u < nbf; u++) {
            P_i->scale_row(0, u, C_lmo_->get(u, i));
            P_i->scale_column(0, u, C_lmo_->get(u, i));
        }

        for (size_t u = 0; u < nbf; u++) {
            int centerU = basisset_->function_to_center(u);
            double p_uu = P_i->get(u, u);

            for (size_t v = 0; v < nbf; v++) {
                int centerV = basisset_->function_to_center(v);
                double p_vv = P_i->get(v, v);

                // off-diag pops (p_uv) split between u and v prop to diag pops
                double p_uv = P_i->get(u, v);
                double p_diag_sum = p_uu + p_vv;
                // Symmetry-exact orbitals can have C_ui = C_vi = 0. In that
                // case p_uv is also zero and the limiting contribution is
                // zero; avoid contaminating the atomic population with 0/0.
                if (p_diag_sum == 0.0) continue;
                mkn_pop[centerU] += p_uv * (p_uu / p_diag_sum);
                mkn_pop[centerV] += p_uv * (p_vv / p_diag_sum);
            }
        }

        // if non-zero mulliken pop on atom, include atom in the LMO's fitting domain
        for (size_t a = 0; a < natom; a++) {
            if (fabs(mkn_pop[a]) > t_cut_mkn_triples) {
                lmo_to_riatoms[i].push_back(a);

                // each atom's aux orbitals are all-or-nothing for each LMO
                for (int u : atom_to_ribf_[a]) {
                    lmo_to_ribfs[i].push_back(u);
                }
            }
        }
    }

    // => PAO domains <= //

    SparseMap lmo_to_paos(naocc);

    double t_cut_do_triples = (prescreening) ? options_.get_double("T_CUT_DO_TRIPLES_PRE") : options_.get_double("T_CUT_DO_TRIPLES");

    for (size_t i = 0; i < naocc; ++i) {
        // PAO domains determined by differential overlap integral
        std::vector<int> lmo_to_paos_temp;
        for (size_t u = 0; u < nbf; ++u) {
            if (fabs(DOI_iu_->get(i, u)) > t_cut_do_triples) {
                lmo_to_paos_temp.push_back(u);
            }
        }

        // if any PAO on an atom is in the list, we take all of the PAOs on that atom
        lmo_to_paos[i] = contract_lists(lmo_to_paos_temp, atom_to_bf_);
    }

    if (!prescreening) {
        lmotriplet_to_ribfs_.clear();
        lmotriplet_to_lmos_.clear();
        lmotriplet_to_paos_.clear();
    }

    lmotriplet_to_ribfs_.resize(n_lmo_triplets);
    lmotriplet_to_lmos_.resize(n_lmo_triplets);
    lmotriplet_to_paos_.resize(n_lmo_triplets);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];
        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];

        lmotriplet_to_ribfs_[ijk] = merge_lists(lmo_to_ribfs[i], merge_lists(lmo_to_ribfs[j], lmo_to_ribfs[k]));
        for (int l = 0; l < naocc; ++l) {
            int il = i_j_to_ij_[i][l], jl = i_j_to_ij_[j][l], kl = i_j_to_ij_[k][l];
            if (il != -1 && jl != -1 && kl != -1) lmotriplet_to_lmos_[ijk].push_back(l);
        }
        lmotriplet_to_paos_[ijk] = merge_lists(lmo_to_paos[i], merge_lists(lmo_to_paos[j], lmo_to_paos[k]));
    }


    timer_off("Triples Sparsity");
}

void DLPNOCCSD_T::sort_triplets(double e_total) {
    /* Sorting triplets by energy values to determine strong and weak
       triplets for the iterative (T) portion of the computation */

    timer_on("Sort Triplets");

    outfile->Printf("  ==> Sorting Triplets <== \n\n");

    int n_lmo_triplets = ijk_to_i_j_k_.size();
    std::vector<std::pair<int, double>> ijk_e_pairs(n_lmo_triplets);

#pragma omp parallel for
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        ijk_e_pairs[ijk] = std::make_pair(ijk, e_ijk_[ijk]);
    }

    std::sort(ijk_e_pairs.begin(), ijk_e_pairs.end(), [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return (std::fabs(a.second) > std::fabs(b.second));
    });

    double e_curr = 0.0;
    const bool full_triples_follow = algorithm_ != DLPNOMethod::CCSD_T;
    const double t_cut_tno_full = options_.get_double("T_CUT_TNO_FULL");
    const double strong_cutoff =
        full_triples_follow
            ? t_cut_tno_full
            : iterative_tno_cutoff(options_, "T_CUT_TNO_STRONG", "T_CUT_TNO_STRONG_SCALE");
    const double weak_cutoff =
        full_triples_follow
            ? t_cut_tno_full
            : iterative_tno_cutoff(options_, "T_CUT_TNO_WEAK", "T_CUT_TNO_WEAK_SCALE");
    is_strong_triplet_.resize(n_lmo_triplets, false);
    tno_cutoff_.assign(n_lmo_triplets, weak_cutoff);

    int strong_count = 0;
    for (int idx = 0; idx < n_lmo_triplets; ++idx) {
        is_strong_triplet_[ijk_e_pairs[idx].first] = true;
        tno_cutoff_[ijk_e_pairs[idx].first] = strong_cutoff;
        e_curr += ijk_e_pairs[idx].second;
        ++strong_count;
        if (e_curr / e_total > 0.9) break;
    }

    outfile->Printf("    Number of Strong Triplets: %6d, Total Triplets: %6d, Ratio: %.4f\n\n", strong_count, n_lmo_triplets, 
                            (double) strong_count / n_lmo_triplets);

    timer_off("Sort Triplets");
}

void DLPNOCCSD_T::tno_transform(double t_cut_tno, bool use_tuple_cutoffs) {
    // Computes TNOs from converged LCCSD densities of triplets ij, jk, and ik
    
    timer_on("TNO transform");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();
    int n_lmo_triplets = ijk_to_i_j_k_.size();
    const int MIN_TNOS = options_.get_int("MIN_TNOS");

    if (use_tuple_cutoffs && tno_cutoff_.size() != n_lmo_triplets) {
        throw PSIEXCEPTION("Triplet-specific TNO cutoffs were not initialized before the iterative (T) transform.");
    }

    X_tno_.clear();
    e_tno_.clear();
    n_tno_.clear();

    X_tno_.resize(n_lmo_triplets);
    e_tno_.resize(n_lmo_triplets);
    n_tno_.resize(n_lmo_triplets);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];
        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];

        // number of PAOs in the triplet domain (before removing linear dependencies)
        int npao_ijk = lmotriplet_to_paos_[ijk].size();

        // number of auxiliary basis in the domain
        int naux_ijk = lmotriplet_to_ribfs_[ijk].size();

        //                                          //
        // ==> Canonicalize PAOs of triplet ijk <== //
        //                                          //

        auto S_pao_ijk = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk], lmotriplet_to_paos_[ijk]);
        auto F_pao_ijk = submatrix_rows_and_cols(*F_pao_, lmotriplet_to_paos_[ijk], lmotriplet_to_paos_[ijk]);

        SharedMatrix X_pao_ijk;
        SharedVector e_pao_ijk;
        std::tie(X_pao_ijk, e_pao_ijk) = orthocanonicalizer(S_pao_ijk, F_pao_ijk);

        F_pao_ijk = linalg::triplet(X_pao_ijk, F_pao_ijk, X_pao_ijk, true, false, false);

        // number of PAOs in the domain after removing linear dependencies
        int npao_can_ijk = X_pao_ijk->ncol();

        // S_ijk partially transformed overlap matrix
        std::vector<int> triples_ext_domain = merge_lists(lmo_to_paos_[i], merge_lists(lmo_to_paos_[j], lmo_to_paos_[k]));
        auto S_ijk = submatrix_rows_and_cols(*S_pao_, triples_ext_domain, lmotriplet_to_paos_[ijk]);
        S_ijk = linalg::doublet(S_ijk, X_pao_ijk, false, false);
        

        //                                           //
        // ==> Canonical PAOs  to Canonical TNOs <== //
        //                                           //

        size_t nvir_ijk = F_pao_ijk->nrow();

        // Construct pair densities from amplitudes
        auto D_ij = linalg::doublet(Tt_iajb_[ij], T_iajb_[ij], false, true);
        D_ij->add(linalg::doublet(Tt_iajb_[ij], T_iajb_[ij], true, false));
        if (i == j) D_ij->scale(0.5);

        auto D_jk = linalg::doublet(Tt_iajb_[jk], T_iajb_[jk], false, true);
        D_jk->add(linalg::doublet(Tt_iajb_[jk], T_iajb_[jk], true, false));
        if (j == k) D_jk->scale(0.5);

        auto D_ik = linalg::doublet(Tt_iajb_[ik], T_iajb_[ik], false, true);
        D_ik->add(linalg::doublet(Tt_iajb_[ik], T_iajb_[ik], true, false));
        if (i == k) D_ik->scale(0.5);

        // Project pair densities into combined PAO space of ijk
        std::vector<int> ij_index = index_list(triples_ext_domain, lmopair_to_paos_[ij]);
        auto S_ij = linalg::doublet(X_pno_[ij], submatrix_rows(*S_ijk, ij_index), true, false);
        D_ij = linalg::triplet(S_ij, D_ij, S_ij, true, false, false);

        std::vector<int> jk_index = index_list(triples_ext_domain, lmopair_to_paos_[jk]);
        auto S_jk = linalg::doublet(X_pno_[jk], submatrix_rows(*S_ijk, jk_index), true, false);
        D_jk = linalg::triplet(S_jk, D_jk, S_jk, true, false, false);

        std::vector<int> ik_index = index_list(triples_ext_domain, lmopair_to_paos_[ik]);
        auto S_ik = linalg::doublet(X_pno_[ik], submatrix_rows(*S_ijk, ik_index), true, false);
        D_ik = linalg::triplet(S_ik, D_ik, S_ik, true, false, false);

        // Construct triplet density from pair densities
        auto D_ijk = D_ij->clone();
        D_ijk->add(D_jk);
        D_ijk->add(D_ik);
        D_ijk->scale(1.0 / 3.0);

        // Diagonalization of triplet density gives TNOs (in basis of LMO's virtual domain)
        // as well as TNO occ numbers
        auto X_tno_ijk = std::make_shared<Matrix>("eigenvectors", nvir_ijk, nvir_ijk);
        Vector tno_occ("eigenvalues", nvir_ijk);
        D_ijk->diagonalize(*X_tno_ijk, tno_occ, descending);

        const double tno_cutoff = use_tuple_cutoffs ? tno_cutoff_[ijk] : t_cut_tno;

        int nvir_ijk_final = 0;
        for (size_t a = 0; a < nvir_ijk; ++a) {
            if (fabs(tno_occ.get(a)) >= tno_cutoff || a < MIN_TNOS) {
                nvir_ijk_final++;
            }
        }

        nvir_ijk_final = std::max(1, nvir_ijk_final);

        Dimension zero(1);
        Dimension dim_final(1);
        dim_final.fill(nvir_ijk_final);

        // This transformation gives orbitals that are orthonormal but not canonical
        X_tno_ijk = X_tno_ijk->get_block({zero, X_tno_ijk->rowspi()}, {zero, dim_final});
        tno_occ = tno_occ.get_block({zero, dim_final});

        SharedMatrix tno_canon;
        SharedVector e_tno_ijk;
        std::tie(tno_canon, e_tno_ijk) = canonicalizer(X_tno_ijk, F_pao_ijk);

        X_tno_ijk = linalg::doublet(X_tno_ijk, tno_canon, false, false);
        X_tno_ijk = linalg::doublet(X_pao_ijk, X_tno_ijk, false, false);

        X_tno_[ijk] = X_tno_ijk;
        e_tno_[ijk] = e_tno_ijk;
        n_tno_[ijk] = X_tno_ijk->ncol();
    }

    int tno_count_total = 0, tno_count_min = C_pao_->ncol(), tno_count_max = 0;
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        tno_count_total += n_tno_[ijk];
        tno_count_min = std::min(tno_count_min, n_tno_[ijk]);
        tno_count_max = std::max(tno_count_max, n_tno_[ijk]);
    }

    // (n + 3 choose 3) minus triples from the same orbital (iii)
    size_t n_total_possible = (naocc + 2) * (naocc + 1) * (naocc) / 6 - naocc;

    outfile->Printf("  \n");
    outfile->Printf("    Number of (Unique) Local MO triplets: %d\n", n_lmo_triplets);
    outfile->Printf("    Max Number of Possible (Unique) LMO Triplets: %d (Ratio: %.4f)\n", n_total_possible,
                    (double)n_lmo_triplets / n_total_possible);
    outfile->Printf("    Natural Orbitals per Local MO triplet:\n");
    outfile->Printf("      Avg: %3d NOs \n", tno_count_total / n_lmo_triplets);
    outfile->Printf("      Min: %3d NOs \n", tno_count_min);
    outfile->Printf("      Max: %3d NOs \n", tno_count_max);
    outfile->Printf("  \n");

    timer_off("TNO transform");
}

void DLPNOCCSD_T::estimate_memory() {
    outfile->Printf("\n  ==> DLPNO-(T) Memory Requirements <== \n\n");

    int n_lmo_triplets = ijk_to_i_j_k_.size();

    size_t tno_total_memory = 0;
#pragma omp parallel for reduction(+ : tno_total_memory)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        tno_total_memory += n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk];
    }

    // Depending on which intermediates are written to disk, the RAM usage varies
    double W3_memory = tno_total_memory;
    double V3_memory = tno_total_memory;
    double T3_memory = tno_total_memory;

    // Write W and V intermediates to disk (this is more efficient than writing amplitudes)
    write_intermediates_ = options_.get_bool("WRITE_TRIPLES_INTERMEDIATES");
    if (write_intermediates_) W3_memory = 0.0, V3_memory = 0.0;

    // Write T3 amplitudes to disk (this should only be turned on as a last resort)
    write_amplitudes_ = options_.get_bool("WRITE_TRIPLES_AMPLITUDES");
    if (write_amplitudes_) T3_memory = 0.0;

    size_t total_memory = qij_memory_ + qia_memory_ + qab_memory_ + W3_memory + V3_memory + T3_memory;

    // 1 GB = 1000^3 = 10^9 Bytes
    const double DOUBLES_TO_GB = pow(10.0, -9) * sizeof(double);
    const double WORDS_TO_GB = pow(10.0, -9);

    outfile->Printf("    (q | i j) integrals    : %.3f [GB]\n", qij_memory_ * DOUBLES_TO_GB);
    outfile->Printf("    (q | i a) integrals    : %.3f [GB]\n", qia_memory_ * DOUBLES_TO_GB);
    outfile->Printf("    (q | a b) integrals    : %.3f [GB]\n", qab_memory_ * DOUBLES_TO_GB);
    outfile->Printf("    W_{ijk}^{abc}          : %.3f [GB]\n", W3_memory * DOUBLES_TO_GB);
    outfile->Printf("    V_{ijk}^{abc}          : %.3f [GB]\n", V3_memory * DOUBLES_TO_GB);
    outfile->Printf("    T_{ijk}^{abc}          : %.3f [GB]\n", T3_memory * DOUBLES_TO_GB);
    outfile->Printf("    Total Memory Required  : %.3f [GB]\n", total_memory * DOUBLES_TO_GB);
    outfile->Printf("    Total Memory Given     : %.3f [GB]\n\n", memory_ * WORDS_TO_GB);
    

    // Memory checks!!!
    bool memory_changed = false;

    if (toggle_memory_ && !write_intermediates_ && total_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is more than 90%% of Available Memory!\n");
        outfile->Printf("    Attempting to switch to disk IO for W and V intermediates...\n");

        total_memory -= (W3_memory + V3_memory);
        W3_memory = 0.0, V3_memory = 0.0;
        write_intermediates_ = true;
        memory_changed = true;
        outfile->Printf("    Required Memory Reduced to %.3f [GB]\n\n", total_memory * DOUBLES_TO_GB);
    }

    if (toggle_memory_ && !write_amplitudes_ && total_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is (still) more than 90%% of Available Memory!\n");
        outfile->Printf("    Attempting to switch to disk IO for T3 amplitudes...\n");

        total_memory -= T3_memory;
        T3_memory = 0.0;
        write_amplitudes_ = true;
        memory_changed = true;
        outfile->Printf("    Required Memory Reduced to %.3f [GB]\n\n", total_memory * DOUBLES_TO_GB);
    }

    // This will likely never be executed, barring a pathological case 
    // (as the memory here is less than what is needed for CCSD)
    if (toggle_memory_ && total_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is (still) more than 90%% of Available Memory!\n");
        throw PSIEXCEPTION("   Too little memory given for DLPNO-(T) Algorithm!");
    }

    if (memory_changed) {
        outfile->Printf("\n  ==> (Updated) DLPNO-(T) Memory Requirements <== \n\n");
        outfile->Printf("    (q | i j) integrals    : %.3f [GB]\n", qij_memory_ * DOUBLES_TO_GB);
        outfile->Printf("    (q | i a) integrals    : %.3f [GB]\n", qia_memory_ * DOUBLES_TO_GB);
        outfile->Printf("    (q | a b) integrals    : %.3f [GB]\n", qab_memory_ * DOUBLES_TO_GB);
        outfile->Printf("    W_{ijk}^{abc}          : %.3f [GB]\n", W3_memory * DOUBLES_TO_GB);
        outfile->Printf("    V_{ijk}^{abc}          : %.3f [GB]\n", T3_memory * DOUBLES_TO_GB);
        outfile->Printf("    T_{ijk}^{abc}          : %.3f [GB]\n", V3_memory * DOUBLES_TO_GB);
        outfile->Printf("    Total Memory Required  : %.3f [GB]\n", total_memory * DOUBLES_TO_GB);
        outfile->Printf("    Total Memory Given     : %.3f [GB]\n\n", memory_ * WORDS_TO_GB);
        
    }

    if (write_intermediates_) {
        outfile->Printf("    Writing W_{ijk}^{abc} and W_{ijk}^{abc} to disk...\n\n");
    }

    if (write_amplitudes_) {
        outfile->Printf("    Writing T_{ijk}^{abc} to disk...\n\n");
    }
    
    if (!write_intermediates_ && !write_amplitudes_) {
        outfile->Printf("    Storing all X_{ijk}^{abc} quantities in RAM...\n\n");
    }
}

double DLPNOCCSD_T::compute_lccsd_t0(bool save_memory) {
    timer_on("LCCSD(T0)");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_triplets = ijk_to_i_j_k_.size();

    double E_T0 = 0.0;

    if (save_memory) {
        W_iajbkc_.resize(n_lmo_triplets);
        V_iajbkc_.resize(n_lmo_triplets);
        T_iajbkc_.resize(n_lmo_triplets);
    }

    e_ijk_.clear();
    e_ijk_.resize(n_lmo_triplets, 0.0);

    std::time_t time_start = std::time(nullptr);
    std::time_t time_lap = std::time(nullptr);

    // Sort Triplets by the approximate number of operations (for maximal parallel efficiency)
    std::vector<std::pair<int, size_t>> ijk_cost_tuple(n_lmo_triplets);
    
#pragma omp parallel for
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        const int npao_ijk = lmotriplet_to_paos_[ijk].size();
        const int naux_ijk = lmotriplet_to_ribfs_[ijk].size();

        // Cost of transforming q_vv from PAO to TNO basis
        size_t cost = naux_ijk * n_tno_[ijk] * npao_ijk * npao_ijk;
        cost += naux_ijk * n_tno_[ijk] * n_tno_[ijk] * npao_ijk;

        ijk_cost_tuple[ijk] = std::make_pair(ijk, cost);
    }
    
    std::sort(ijk_cost_tuple.begin(), ijk_cost_tuple.end(), [&](const std::pair<int, size_t>& a, const std::pair<int, size_t>& b) {
        return (a.second > b.second);
    });

    std::vector<int> ijk_sorted_by_cost(n_lmo_triplets);
    
#pragma omp parallel for
    for (int ijk_idx = 0; ijk_idx < n_lmo_triplets; ++ijk_idx) {
        ijk_sorted_by_cost[ijk_idx] = ijk_cost_tuple[ijk_idx].first;
    }

#pragma omp parallel for schedule(dynamic) reduction(+ : E_T0)
    for (int ijk_idx = 0; ijk_idx < n_lmo_triplets; ++ijk_idx) {
        // Triplets assigned to threads dynamically, sorted in descending order of cost
        // This maximizes parallel efficiency
        int ijk = ijk_sorted_by_cost[ijk_idx];

        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];
        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];

        int ntno_ijk = n_tno_[ijk];

        if (ntno_ijk == 0) continue;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        if (thread == 0) timer_on("LCCSD(T0): Setup Integrals");

        // => Step 1: Compute all necessary integrals

        // number of LMOs in the triplet domain
        const int nlmo_ijk = lmotriplet_to_lmos_[ijk].size();
        // number of PAOs in the triplet domain (before removing linear dependencies)
        const int npao_ijk = lmotriplet_to_paos_[ijk].size();
        // number of auxiliary functions in the triplet domain
        const int naux_ijk = lmotriplet_to_ribfs_[ijk].size();

        // number of PAOs in the pair domains of ij, jk, and ik
        const int npao_ij = lmopair_to_paos_[ij].size(), npao_jk = lmopair_to_paos_[jk].size(), npao_ik = lmopair_to_paos_[ik].size();

        /// => Build (i a_ijk | b_ijk d_jk) and (k c_ijk | j l) integrals <= ///

        auto q_iv = std::make_shared<Matrix>(naux_ijk, npao_ijk); // (Q_{ijk} | i u_{ijk})
        auto q_jv = std::make_shared<Matrix>(naux_ijk, npao_ijk); // (Q_{ijk} | j u_{ijk})
        auto q_kv = std::make_shared<Matrix>(naux_ijk, npao_ijk); // (Q_{ijk} | k u_{ijk})

        auto q_io = std::make_shared<Matrix>(naux_ijk, nlmo_ijk); // (Q_{ijk} | m_{ijk} i)
        auto q_jo = std::make_shared<Matrix>(naux_ijk, nlmo_ijk); // (Q_{ijk} | m_{ijk} j)
        auto q_ko = std::make_shared<Matrix>(naux_ijk, nlmo_ijk); // (Q_{ijk} | m_{ijk} k)

        auto q_vv = std::make_shared<Matrix>(naux_ijk, ntno_ijk * ntno_ijk); // (Q_{ijk} | a_{ijk} b_{ijk})

        for (int q_ijk = 0; q_ijk < naux_ijk; q_ijk++) {
            const int q = lmotriplet_to_ribfs_[ijk][q_ijk];
            const int centerq = ribasis_->function_to_center(q);

            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                (*q_io)(q_ijk, l_ijk) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][i], riatom_to_lmos_ext_dense_[centerq][l]);
                (*q_jo)(q_ijk, l_ijk) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][j], riatom_to_lmos_ext_dense_[centerq][l]);
                (*q_ko)(q_ijk, l_ijk) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][k], riatom_to_lmos_ext_dense_[centerq][l]);
            }

            for (int u_ijk = 0; u_ijk < npao_ijk; ++u_ijk) {
                int u = lmotriplet_to_paos_[ijk][u_ijk];
                (*q_iv)(q_ijk, u_ijk) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][i], riatom_to_paos_ext_dense_[centerq][u]);
                (*q_jv)(q_ijk, u_ijk) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][j], riatom_to_paos_ext_dense_[centerq][u]);
                (*q_kv)(q_ijk, u_ijk) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][k], riatom_to_paos_ext_dense_[centerq][u]);
            }

            auto q_vv_tmp = std::make_shared<Matrix>(npao_ijk, npao_ijk);
            q_vv_tmp->zero();

            for (int u_ijk = 0; u_ijk < npao_ijk; ++u_ijk) {
                int u = lmotriplet_to_paos_[ijk][u_ijk];
                for (int v_ijk = 0; v_ijk < npao_ijk; ++v_ijk) {
                    int v = lmotriplet_to_paos_[ijk][v_ijk];
                    int uv_idx = riatom_to_pao_pairs_dense_[centerq][u][v];
                    if (uv_idx == -1) continue;
                    (*q_vv_tmp)(u_ijk, v_ijk) = (*qab_[q])(uv_idx, 0);
                } // end v_ijk
            } // end u_ijk
            
            // naux_{ijk} * npao_{ijk}^{2} * ntno_{ijk} (this is the most expensive operation in this loop)
            q_vv_tmp = linalg::triplet(X_tno_[ijk], q_vv_tmp, X_tno_[ijk], true, false, false);
            ::memcpy(&(*q_vv)(q_ijk, 0), &(*q_vv_tmp)(0, 0), ntno_ijk * ntno_ijk * sizeof(double));

            // naux_ijk * npao_ijk^2 * ntno_{ijk}
        } // end q_ijk

        q_iv = linalg::doublet(q_iv, X_tno_[ijk]); // (Q_{ijk} | i u_{ijk}) -> (Q_{ijk} | i a_{ijk})
        q_jv = linalg::doublet(q_jv, X_tno_[ijk]); // (Q_{ijk} | j u_{ijk}) -> (Q_{ijk} | j a_{ijk})
        q_kv = linalg::doublet(q_kv, X_tno_[ijk]); // (Q_{ijk} | k u_{ijk}) -> (Q_{ijk} | k a_{ijk})
        
        auto q_iv_clone = q_iv->clone();
        auto q_jv_clone = q_jv->clone();
        auto q_kv_clone = q_kv->clone();

        auto A_solve = submatrix_rows_and_cols(*full_metric_, lmotriplet_to_ribfs_[ijk], lmotriplet_to_ribfs_[ijk]);

        /* These are cloned and inverted by the full coulomb metric (not to the half power)
            to make formation of (i a | b c)-type integrals more efficient later */

        C_DGESV_wrapper(A_solve->clone(), q_iv_clone);
        C_DGESV_wrapper(A_solve->clone(), q_jv_clone);
        C_DGESV_wrapper(A_solve->clone(), q_kv_clone);

        A_solve->power(-0.5, 1.0e-14);

        q_iv = linalg::doublet(A_solve, q_iv);
        q_jv = linalg::doublet(A_solve, q_jv);
        q_kv = linalg::doublet(A_solve, q_kv);
        q_io = linalg::doublet(A_solve, q_io);
        q_jo = linalg::doublet(A_solve, q_jo);
        q_ko = linalg::doublet(A_solve, q_ko);
        A_solve.reset();

        if (thread == 0) timer_off("LCCSD(T0): Setup Integrals");

        if (thread == 0) timer_on("LCCSD(T0): Contract Integrals");

        // W integrals
        auto K_ivvv = linalg::doublet(q_iv_clone, q_vv, true, false); // (i a_{ijk} | b_{ijk} d_{ijk})
        auto K_jvvv = linalg::doublet(q_jv_clone, q_vv, true, false); // (j b_{ijk} | c_{ijk} d_{ijk})
        auto K_kvvv = linalg::doublet(q_kv_clone, q_vv, true, false); // (k c_{ijk} | a_{ijk} d_{ijk})

        auto K_iojv = linalg::doublet(q_io, q_jv, true, false); // (i l_{ijk} | j b_{ijk})
        auto K_joiv = linalg::doublet(q_jo, q_iv, true, false); // (j l_{ijk} | i a_{ijk})
        auto K_kojv = linalg::doublet(q_ko, q_jv, true, false); // (k l_{ijk} | j b_{ijk})
        auto K_jokv = linalg::doublet(q_jo, q_kv, true, false); // (j l_{ijk} | k c_{ijk})
        auto K_iokv = linalg::doublet(q_io, q_kv, true, false); // (i l_{ijk} | k c_{ijk})
        auto K_koiv = linalg::doublet(q_ko, q_iv, true, false); // (k l_{ijk} | i a_{ijk})

        // V integrals
        auto K_jk = linalg::doublet(q_jv, q_kv, true, false); // (j b_{ijk} | k c_{ijk})
        auto K_ik = linalg::doublet(q_iv, q_kv, true, false); // (i a_{ijk} | k c_{ijk})
        auto K_ij = linalg::doublet(q_iv, q_jv, true, false); // (i a_{ijk} | j b_{ijk})

        // S integrals (semi-direct algorithm)
        std::vector<int> triples_ext_domain = merge_lists(lmo_to_paos_[i], merge_lists(lmo_to_paos_[j], lmo_to_paos_[k]));
        for (int l_ijk = 0; l_ijk < lmotriplet_to_lmos_[ijk].size(); ++l_ijk) {
            int l = lmotriplet_to_lmos_[ijk][l_ijk];
            triples_ext_domain = merge_lists(triples_ext_domain, lmo_to_paos_[l]);
        }
        auto S_ijk = submatrix_rows_and_cols(*S_pao_, triples_ext_domain, lmotriplet_to_paos_[ijk]);
        S_ijk = linalg::doublet(S_ijk, X_tno_[ijk], false, false);

        // => Step 1: Compute W_ijk (Jiang Eq. 109) <= //
        // W_{ijk}^{abc} = P_{ijk}^{abc}[(ia|bd)t_{kj}^{cd} - t_{il}^{ab}(jl|kc)]
        // P_{ijk}^{abc} is explicitly applied through the perms

        std::stringstream w_name;
        w_name << "W " << (ijk);
        auto W_ijk = std::make_shared<Matrix>(w_name.str(), ntno_ijk, ntno_ijk * ntno_ijk);
        W_ijk->zero();

        std::vector<std::tuple<int, int, int>> perms = {std::make_tuple(i, j, k), std::make_tuple(i, k, j),
                                                        std::make_tuple(j, i, k), std::make_tuple(j, k, i),
                                                        std::make_tuple(k, i, j), std::make_tuple(k, j, i)};
        std::vector<SharedMatrix> Wperms(perms.size());

        std::vector<SharedMatrix> K_ovvv_list = {K_ivvv, K_ivvv, K_jvvv, K_jvvv, K_kvvv, K_kvvv};
        std::vector<SharedMatrix> K_ooov_list = {K_jokv, K_kojv, K_iokv, K_koiv, K_iojv, K_joiv};

        if (thread == 0) timer_off("LCCSD(T0): Contract Integrals");

        if (thread == 0) timer_on("LCCSD(T0): Form W");

        for (int idx = 0; idx < perms.size(); ++idx) {
            int i, j, k;
            std::tie(i, j, k) = perms[idx];

            int ii = i_j_to_ij_[i][i];
            int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];
            int kj = ij_to_ji_[jk];

            Wperms[idx] = std::make_shared<Matrix>(ntno_ijk, ntno_ijk * ntno_ijk);
            Wperms[idx]->zero();
            
            // Compute overlap between TNOs of triplet ijk and PNOs of pair kj
            std::vector<int> kj_idx_list = index_list(triples_ext_domain, lmopair_to_paos_[kj]);
            auto S_kj_ijk = linalg::doublet(X_pno_[kj], submatrix_rows(*S_ijk, kj_idx_list), true, false);
            // (c_{kj}, d_{kj}) -> (c_{ijk}, d_{ijk})
            auto T_kj = linalg::triplet(S_kj_ijk, T_iajb_[kj], S_kj_ijk, true, false, false); 

            auto K_ovvv = K_ovvv_list[idx]->clone(); // (i a | b d) stored as: (a, b * d)

            // Jiang Eq. 109a
            // W_{ijk}^{abc} += (ia|bd)t_{kj}^{cd}
            K_ovvv->reshape(ntno_ijk * ntno_ijk, ntno_ijk);  // (a, b * d) -> (a * b, d)
            K_ovvv = linalg::doublet(K_ovvv, T_kj, false, true); // (a * b, d) (c, d) -> (a * b, c)
            K_ovvv->reshape(ntno_ijk, ntno_ijk * ntno_ijk); // (a * b, c) -> (a, b * c)
            Wperms[idx]->add(K_ovvv);

            for (int l_ijk = 0; l_ijk < lmotriplet_to_lmos_[ijk].size(); ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                int il = i_j_to_ij_[i][l];

                // Compute overlap between TNOs of triplet ijk and PNOs of pair il
                std::vector<int> il_idx_list = index_list(triples_ext_domain, lmopair_to_paos_[il]);
                auto S_il_ijk = linalg::doublet(X_pno_[il], submatrix_rows(*S_ijk, il_idx_list), true, false);
                // (a_{il}, b_{il}) -> (a_{ijk}, b_{ijk})
                auto T_il = linalg::triplet(S_il_ijk, T_iajb_[il], S_il_ijk, true, false, false);

                // Jiang Eq. 109b
                // W_{ijk}^{abc} -= t_{il}^{ab}(jl|kc)
                for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
                    for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                        for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                            (*Wperms[idx])(a_ijk, b_ijk * ntno_ijk + c_ijk) -=
                                (*T_il)(a_ijk, b_ijk) * (*K_ooov_list[idx])(l_ijk, c_ijk); // (a, b) * (c) -> (a, b * c)
                        }
                    }
                }  // end a_ijk
            }      // end l_ijk
        }

        // Encapsulates the P_{ijk}^{abc} permutation
        // Reminder: P_{ijk}^{abc}X_{ijk}^{abc} =>
        // X_{ijk}^{abc} + X_{ikj}^{acb} + X_{jik}^{bac} + X_{jki}^{bca} + X_{kij}^{cab} + X_{kji}^{cba}
        for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
            for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                    (*W_ijk)(a_ijk, b_ijk *ntno_ijk + c_ijk) =
                        (*Wperms[0])(a_ijk, b_ijk * ntno_ijk + c_ijk) + (*Wperms[1])(a_ijk, c_ijk * ntno_ijk + b_ijk) +
                        (*Wperms[2])(b_ijk, a_ijk * ntno_ijk + c_ijk) + (*Wperms[3])(b_ijk, c_ijk * ntno_ijk + a_ijk) +
                        (*Wperms[4])(c_ijk, a_ijk * ntno_ijk + b_ijk) + (*Wperms[5])(c_ijk, b_ijk * ntno_ijk + a_ijk);
                }
            }
        }

        if (thread == 0) timer_off("LCCSD(T0): Form W");

        if (thread == 0) timer_on("LCCSD(T0): Form V");

        // => Step 2: Compute V_ijk (Jiang Eq. 110) <= //
        // V_{ijk}^{abc} = W_{ijk}^{abc} + T_{i}^{a}(jb|kc) + T_{j}^{b}(ia|kc) + T_{k}^{c}(ia|jb)

        auto V_ijk = W_ijk->clone();
        std::stringstream v_name;
        v_name << "V " << (ijk);
        V_ijk->set_name(v_name.str());

        // Compute overlap between TNOs of triplet ijk and PNOs of pair ii, jj, and kk
        int ii = i_j_to_ij_[i][i];
        std::vector<int> ii_idx_list = index_list(triples_ext_domain, lmopair_to_paos_[ii]);
        auto S_ii_ijk = linalg::doublet(X_pno_[ii], submatrix_rows(*S_ijk, ii_idx_list), true, false);

        int jj = i_j_to_ij_[j][j];
        std::vector<int> jj_idx_list = index_list(triples_ext_domain, lmopair_to_paos_[jj]);
        auto S_jj_ijk = linalg::doublet(X_pno_[jj], submatrix_rows(*S_ijk, jj_idx_list), true, false);

        int kk = i_j_to_ij_[k][k];
        std::vector<int> kk_idx_list = index_list(triples_ext_domain, lmopair_to_paos_[kk]);
        auto S_kk_ijk = linalg::doublet(X_pno_[kk], submatrix_rows(*S_ijk, kk_idx_list), true, false);

        // Transform singles amplitude to TNO space
        auto T_i = linalg::doublet(S_ii_ijk, T_ia_[i], true, false); // (i, a_{ii}) -> (i, a_{ijk})
        auto T_j = linalg::doublet(S_jj_ijk, T_ia_[j], true, false); // (j, b_{ii}) -> (j, b_{ijk})
        auto T_k = linalg::doublet(S_kk_ijk, T_ia_[k], true, false); // (k, c_{ii}) -> (k, c_{ijk})

        for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
            for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                    (*V_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) += (*T_i)(a_ijk, 0) * (*K_jk)(b_ijk, c_ijk) +
                        (*T_j)(b_ijk, 0) * (*K_ik)(a_ijk, c_ijk) + (*T_k)(c_ijk, 0) * (*K_ij)(a_ijk, b_ijk);
                } // end c_ijk
            } // end b_ijk
        } // end a_ijk

        if (thread == 0) timer_off("LCCSD(T0): Form V");

        // Step 3: Compute T0 energy through amplitudes (Jiang Eq. 53)

        // T_{ijk}^{abc} = W_{ijk}^{abc} (\eps_{ijk}^{abc})^{-1} 
        // (initial semicanonical T3 amplitudes)
        auto T_ijk = W_ijk->clone();
        std::stringstream t_name;
        t_name << "T " << (ijk);
        T_ijk->set_name(t_name.str());

        for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
            for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                    (*T_ijk)(a_ijk, b_ijk *ntno_ijk + c_ijk) =
                        -(*T_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) /
                        (e_tno_[ijk]->get(a_ijk) + e_tno_[ijk]->get(b_ijk) + e_tno_[ijk]->get(c_ijk) - (*F_lmo_)(i, i) -
                         (*F_lmo_)(j, j) - (*F_lmo_)(k, k));
                }
            }
        }

        /* E_{(T)} = prefactor * T_{ijk}^{abc} * (8 V_{ijk}^{abc} - 4 V_{ijk}^{bac} - 4 V_{ijk}^{acb}
                        - 4 V_{ijk}^{cab} + 2 V_{ijk}^{bca} + 2 V_{ijk}^{cab}) (Jiang Eq. 53) */

        double prefactor = 1.0;
        if (i == j && j == k) {
            prefactor /= 6.0;
        } else if (i == j || j == k || i == k) {
            prefactor /= 2.0;
        }

        e_ijk_[ijk] += 8.0 * prefactor * V_ijk->vector_dot(T_ijk);
        e_ijk_[ijk] -= 4.0 * prefactor * triples_permuter(V_ijk, k, j, i)->vector_dot(T_ijk);
        e_ijk_[ijk] -= 4.0 * prefactor * triples_permuter(V_ijk, i, k, j)->vector_dot(T_ijk);
        e_ijk_[ijk] -= 4.0 * prefactor * triples_permuter(V_ijk, j, i, k)->vector_dot(T_ijk);
        e_ijk_[ijk] += 2.0 * prefactor * triples_permuter(V_ijk, j, k, i)->vector_dot(T_ijk);
        e_ijk_[ijk] += 2.0 * prefactor * triples_permuter(V_ijk, k, i, j)->vector_dot(T_ijk);

        E_T0 += e_ijk_[ijk];

        // Step 4: Save Matrices (if doing full (T))

        if (save_memory && !write_intermediates_) {
            W_iajbkc_[ijk] = W_ijk;
            V_iajbkc_[ijk] = V_ijk;
        } else if (save_memory && write_intermediates_) {
#pragma omp critical
            W_ijk->save(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
#pragma omp critical
            V_ijk->save(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
        }

        if (save_memory && !write_amplitudes_) {
            T_iajbkc_[ijk] = T_ijk;
        } else if (save_memory && write_amplitudes_) {
#pragma omp critical
            T_ijk->save(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
        }

        if (thread == 0) {
            std::time_t time_curr = std::time(nullptr);
            int time_elapsed = (int) time_curr - (int) time_lap;
            if (time_elapsed > 60) {
                outfile->Printf("  Time Elapsed from last checkpoint %4d (s), Progress %2d %%, Amplitudes for (%6d / %6d) Triplets Computed\n", time_elapsed, 
                                    (100 * ijk_idx) / n_lmo_triplets, ijk_idx, n_lmo_triplets);
                time_lap = std::time(nullptr);
            }
        }
    }

    timer_off("LCCSD(T0)");

    std::time_t time_stop = std::time(nullptr);
    int time_elapsed = (int) time_stop - (int) time_start;
    outfile->Printf("    (Relavent) Semicanonical LCCSD(T0) Computation Complete!!! Time Elapsed: %4d seconds\n\n", time_elapsed);

    return E_T0;
}

double DLPNOCCSD_T::compute_t_iteration_energy() {
    /* Jiang Eq. 53 */
    /* E_{(T)} = prefactor * T_{ijk}^{abc} *(8 V_{ijk}^{abc} - 4 V_{ijk}^{bac} - 4 V_{ijk}^{acb}
                    - 4 V_{ijk}^{cab} + 2 V_{ijk}^{bca} + 2 V_{ijk}^{cab}) */

    timer_on("Compute (T) Energy");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_triplets = ijk_to_i_j_k_.size();

    double E_T = 0.0;

#pragma omp parallel for schedule(dynamic) reduction(+ : E_T)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

        int ntno_ijk = n_tno_[ijk];
        if (ntno_ijk == 0) continue;

        int kji = i_j_k_to_ijk_[triplet_key(k, j, i, naocc)];
        int ikj = i_j_k_to_ijk_[triplet_key(i, k, j, naocc)];
        int jik = i_j_k_to_ijk_[triplet_key(j, i, k, naocc)];
        int jki = i_j_k_to_ijk_[triplet_key(j, k, i, naocc)];
        int kij = i_j_k_to_ijk_[triplet_key(k, i, j, naocc)];

        double prefactor = 1.0;
        if (i == j && j == k) {
            prefactor /= 6.0;
        } else if (i == j || j == k || i == k) {
            prefactor /= 2.0;
        }

        SharedMatrix V_ijk;
        SharedMatrix T_ijk;

        // Grab V3 and T3 as needed
        if (write_intermediates_) {
            std::stringstream v_name;
            v_name << "V " << (ijk);
            V_ijk = std::make_shared<Matrix>(v_name.str(), ntno_ijk, ntno_ijk * ntno_ijk);
#pragma omp critical
            V_ijk->load(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
        } else {
            V_ijk = V_iajbkc_[ijk];
        }

        if (write_amplitudes_) {
            std::stringstream t_name;
            t_name << "T " << (ijk);
            T_ijk = std::make_shared<Matrix>(t_name.str(), ntno_ijk, ntno_ijk * ntno_ijk);
#pragma omp critical
            T_ijk->load(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
        } else {
            T_ijk = T_iajbkc_[ijk];
        }

        e_ijk_[ijk] = 8.0 * prefactor * V_ijk->vector_dot(T_ijk);
        e_ijk_[ijk] -= 4.0 * prefactor * triples_permuter(V_ijk, k, j, i)->vector_dot(T_ijk);
        e_ijk_[ijk] -= 4.0 * prefactor * triples_permuter(V_ijk, i, k, j)->vector_dot(T_ijk);
        e_ijk_[ijk] -= 4.0 * prefactor * triples_permuter(V_ijk, j, i, k)->vector_dot(T_ijk);
        e_ijk_[ijk] += 2.0 * prefactor * triples_permuter(V_ijk, j, k, i)->vector_dot(T_ijk);
        e_ijk_[ijk] += 2.0 * prefactor * triples_permuter(V_ijk, k, i, j)->vector_dot(T_ijk);

        E_T += e_ijk_[ijk];
    }

    timer_off("Compute (T) Energy");

    return E_T;
}

SharedMatrix DLPNOCCSD_T::triples_permuter(const SharedMatrix &X, int i, int j, int k, bool reverse) {
    /* A helper function that permutes the indices of a triples amplitude based on ordering i, j, k.
        This is helpful for getting the true permutation of the virtual indices for triples that are
        not in i <= j <= k ordering */

    SharedMatrix Xperm = X->clone();
    int ntno_ijk = X->nrow();

    int perm_idx;
    if (i <= j && j <= k && i <= k) {
        perm_idx = 0;
    } else if (i <= k && k <= j && i <= j) {
        perm_idx = 1;
    } else if (j <= i && i <= k && j <= k) {
        perm_idx = 2;
    } else if (j <= k && k <= i && j <= i) {
        perm_idx = 3;
    } else if (k <= i && i <= j && k <= j) {
        perm_idx = 4;
    } else {
        perm_idx = 5;
    }

    for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
        for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
            for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                if (perm_idx == 0)
                    (*Xperm)(a_ijk, b_ijk * ntno_ijk + c_ijk) = (*X)(a_ijk, b_ijk * ntno_ijk + c_ijk);
                else if (perm_idx == 1)
                    (*Xperm)(a_ijk, b_ijk * ntno_ijk + c_ijk) = (*X)(a_ijk, c_ijk * ntno_ijk + b_ijk);
                else if (perm_idx == 2)
                    (*Xperm)(a_ijk, b_ijk * ntno_ijk + c_ijk) = (*X)(b_ijk, a_ijk * ntno_ijk + c_ijk);
                else if ((perm_idx == 3 && !reverse) || (perm_idx == 4 && reverse))
                    (*Xperm)(a_ijk, b_ijk * ntno_ijk + c_ijk) = (*X)(b_ijk, c_ijk * ntno_ijk + a_ijk);
                else if ((perm_idx == 4 && !reverse) || (perm_idx == 3 && reverse))
                    (*Xperm)(a_ijk, b_ijk * ntno_ijk + c_ijk) = (*X)(c_ijk, a_ijk * ntno_ijk + b_ijk);
                else
                    (*Xperm)(a_ijk, b_ijk * ntno_ijk + c_ijk) = (*X)(c_ijk, b_ijk * ntno_ijk + a_ijk);
            }
        }
    }

    return Xperm;
}

double DLPNOCCSD_T::lccsd_t_iterations() {
    timer_on("LCCSD(T) Iterations");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_triplets = ijk_to_i_j_k_.size();

    outfile->Printf("\n  ==> Local CCSD(T) <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                         Corr. Energy    Delta E     Max R     Time (s)\n");

    int iteration = 1, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, r_curr = 0.0;
    bool e_converged = false, r_converged = false;

    double F_CUT = options_.get_double("F_CUT_T");
    double T_CUT_ITER = options_.get_double("T_CUT_ITER");

    std::vector<double> e_ijk_old(n_lmo_triplets, 0.0);

    // Sort Triplets by the approximate number of operations (for maximal parallel efficiency)
    std::vector<std::pair<int, size_t>> ijk_cost_tuple(n_lmo_triplets);
    
#pragma omp parallel for
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

        size_t cost = 0;

        // Compute the cost of the TNO projections per triple
        for (int l = 0; l < naocc; ++l) {
            const size_t ijl_dense = triplet_key(i, j, l, naocc);
            const size_t ilk_dense = triplet_key(i, l, k, naocc);
            const size_t ljk_dense = triplet_key(l, j, k, naocc);
            
            if (l != k && i_j_k_to_ijk_.count(ijl_dense) && std::fabs((*F_lmo_)(l, k)) >= F_CUT) {
                int ijl = i_j_k_to_ijk_[ijl_dense];
                cost += n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijl];
                cost += n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijl] * n_tno_[ijl];
                cost += n_tno_[ijk] * n_tno_[ijl] * n_tno_[ijl] * n_tno_[ijl];
            }
            
            if (l != j && i_j_k_to_ijk_.count(ilk_dense) && std::fabs((*F_lmo_)(l, j)) >= F_CUT) {
                int ilk = i_j_k_to_ijk_[ilk_dense];
                cost += n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * n_tno_[ilk];
                cost += n_tno_[ijk] * n_tno_[ijk] * n_tno_[ilk] * n_tno_[ilk];
                cost += n_tno_[ijk] * n_tno_[ilk] * n_tno_[ilk] * n_tno_[ilk];
            }
            
            if (l != i && i_j_k_to_ijk_.count(ljk_dense) && std::fabs((*F_lmo_)(l, i)) >= F_CUT) {
                int ljk = i_j_k_to_ijk_[ljk_dense];
                cost += n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * n_tno_[ljk];
                cost += n_tno_[ijk] * n_tno_[ijk] * n_tno_[ljk] * n_tno_[ljk];
                cost += n_tno_[ijk] * n_tno_[ljk] * n_tno_[ljk] * n_tno_[ljk];
            }
        }

        ijk_cost_tuple[ijk] = std::make_pair(ijk, cost);
    }
    
    std::sort(ijk_cost_tuple.begin(), ijk_cost_tuple.end(), [&](const std::pair<int, size_t>& a, const std::pair<int, size_t>& b) {
        return (a.second > b.second);
    });

    std::vector<int> ijk_sorted_by_cost(n_lmo_triplets);
    
#pragma omp parallel for
    for (int ijk_idx = 0; ijk_idx < n_lmo_triplets; ++ijk_idx) {
        ijk_sorted_by_cost[ijk_idx] = ijk_cost_tuple[ijk_idx].first;
    }

    while (!(e_converged && r_converged)) {
        // RMS of residual per single LMO, for assesing convergence
        std::vector<double> R_iajbkc_rms(n_lmo_triplets, 0.0);

        std::time_t time_start = std::time(nullptr);

#pragma omp parallel for schedule(dynamic, 1)
        for (int ijk_idx = 0; ijk_idx < n_lmo_triplets; ++ijk_idx) {
            // Triplets assigned to threads dynamically, sorted in descending order of cost
            // This maximizes parallel efficiency
            int ijk = ijk_sorted_by_cost[ijk_idx];

            int i, j, k;
            std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

            int ntno_ijk = n_tno_[ijk];

            if (std::fabs(e_ijk_[ijk] - e_ijk_old[ijk]) < std::fabs(e_ijk_old[ijk] * T_CUT_ITER)) continue;

            // S integrals
            std::vector<int> triplet_ext_domain;
            for (int l = 0; l < naocc; ++l) {
                const size_t ijl_dense = triplet_key(i, j, l, naocc);
                const size_t ilk_dense = triplet_key(i, l, k, naocc);
                const size_t ljk_dense = triplet_key(l, j, k, naocc);
                
                if (l != k && i_j_k_to_ijk_.count(ijl_dense) && std::fabs((*F_lmo_)(l, k)) >= F_CUT) {
                    int ijl = i_j_k_to_ijk_[ijl_dense];
                    triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[ijl]);
                }
                
                if (l != j && i_j_k_to_ijk_.count(ilk_dense) && std::fabs((*F_lmo_)(l, j)) >= F_CUT) {
                    int ilk = i_j_k_to_ijk_[ilk_dense];
                    triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[ilk]);
                }
                
                if (l != i && i_j_k_to_ijk_.count(ljk_dense) && std::fabs((*F_lmo_)(l, i)) >= F_CUT) {
                    int ljk = i_j_k_to_ijk_[ljk_dense];
                    triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[ljk]);
                    
                }
            }

            // Overlap integrals are formed using semi-direct algorithm
            auto S_ijk = submatrix_rows_and_cols(*S_pao_, triplet_ext_domain, lmotriplet_to_paos_[ijk]);
            S_ijk = linalg::doublet(S_ijk, X_tno_[ijk], false, false);

            auto R_ijk = std::make_shared<Matrix>("R_ijk", ntno_ijk, ntno_ijk * ntno_ijk);
            
            SharedMatrix W_ijk;
            SharedMatrix T_ijk;

            // Grab W3 and T3 as needed
            if (write_intermediates_) {
                std::stringstream w_name;
                w_name << "W " << (ijk);
                W_ijk = std::make_shared<Matrix>(w_name.str(), ntno_ijk, ntno_ijk * ntno_ijk);
#pragma omp critical
                W_ijk->load(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
            } else {
                W_ijk = W_iajbkc_[ijk];
            }

            if (write_amplitudes_) {
                std::stringstream t_name;
                t_name << "T " << (ijk);
                T_ijk = std::make_shared<Matrix>(t_name.str(), ntno_ijk, ntno_ijk * ntno_ijk);
#pragma omp critical
                T_ijk->load(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
            } else {
                T_ijk = T_iajbkc_[ijk];
            }

            // => Jiang Eq. 111 <= //
            // R_{ijk}^{abc} = W_{ijk}^{abc} + (typo in original work, plus, not minus)
            // T_{ijk}^{abc} (e_{a} + e_{b} + e_{c} - f_{ii} - f_{jj} - f_{kk})
            // - f_{il} t_{ljk}^{abc} - f_{jl} t_{ilk}^{abc} - f_{kl} t_{ijl}^{abc}
            R_ijk->copy(W_ijk);

            for (int a_ijk = 0; a_ijk < ntno_ijk; ++a_ijk) {
                for (int b_ijk = 0; b_ijk < ntno_ijk; ++b_ijk) {
                    for (int c_ijk = 0; c_ijk < ntno_ijk; ++c_ijk) {
                        (*R_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) += (*T_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) *
                            ((*e_tno_[ijk])(a_ijk) + (*e_tno_[ijk])(b_ijk) + (*e_tno_[ijk])(c_ijk) 
                                - (*F_lmo_)(i, i) - (*F_lmo_)(j, j) - (*F_lmo_)(k, k));
                    }
                }
            }

            // This inner loop performs the operation
            // R_{ijk}^{abc} += - f_{il} t_{ljk}^{abc} - f_{jl} t_{ilk}^{abc} - f_{kl} t_{ijl}^{abc}
            for (int l = 0; l < naocc; l++) {
                const size_t ijl_dense = triplet_key(i, j, l, naocc);
                if (l != k && i_j_k_to_ijk_.count(ijl_dense) && std::fabs((*F_lmo_)(l, k)) >= F_CUT) {
                    int ijl = i_j_k_to_ijk_[ijl_dense];

                    // S(a_{ijk}, a_{ijl})
                    std::vector<int> ijl_idx_list = index_list(triplet_ext_domain, lmotriplet_to_paos_[ijl]);
                    auto S_ijk_ijl = linalg::doublet(submatrix_rows(*S_ijk, ijl_idx_list), X_tno_[ijl], true, false);

                    SharedMatrix T_ijl;
                    if (write_amplitudes_) {
                        std::stringstream t_name;
                        t_name << "T " << (ijl);
                        T_ijl = std::make_shared<Matrix>(t_name.str(), n_tno_[ijl], n_tno_[ijl] * n_tno_[ijl]);
#pragma omp critical
                        T_ijl->load(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
                    } else {
                        T_ijl = T_iajbkc_[ijl];
                    }

                    // (a_{ijl}, b_{ijl}, c_{ijl}) -> (a_{ijk}, b_{ijk}, c_{ijk})
                    auto T_temp1 =
                        matmul_3d(triples_permuter(T_ijl, i, j, l), S_ijk_ijl, n_tno_[ijl], n_tno_[ijk]);
                    C_DAXPY(ntno_ijk * ntno_ijk * ntno_ijk, -(*F_lmo_)(l, k), &(*T_temp1)(0, 0), 1,
                            &(*R_ijk)(0, 0), 1);
                }

                const size_t ilk_dense = triplet_key(i, l, k, naocc);
                if (l != j && i_j_k_to_ijk_.count(ilk_dense) && std::fabs((*F_lmo_)(l, j)) >= F_CUT) {
                    int ilk = i_j_k_to_ijk_[ilk_dense];

                    // S(a_{ijk}, a_{ilk})
                    std::vector<int> ilk_idx_list = index_list(triplet_ext_domain, lmotriplet_to_paos_[ilk]);
                    auto S_ijk_ilk = linalg::doublet(submatrix_rows(*S_ijk, ilk_idx_list), X_tno_[ilk], true, false);

                    SharedMatrix T_ilk;
                    if (write_amplitudes_) {
                        std::stringstream t_name;
                        t_name << "T " << (ilk);
                        T_ilk = std::make_shared<Matrix>(t_name.str(), n_tno_[ilk], n_tno_[ilk] * n_tno_[ilk]);
#pragma omp critical
                        T_ilk->load(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
                    } else {
                        T_ilk = T_iajbkc_[ilk];
                    }

                    // (a_{ilk}, b_{ilk}, c_{ilk}) -> (a_{ijk}, b_{ijk}, c_{ijk})
                    auto T_temp1 =
                        matmul_3d(triples_permuter(T_ilk, i, l, k), S_ijk_ilk, n_tno_[ilk], n_tno_[ijk]);
                    C_DAXPY(ntno_ijk * ntno_ijk * ntno_ijk, -(*F_lmo_)(l, j), &(*T_temp1)(0, 0), 1,
                            &(*R_ijk)(0, 0), 1);
                }

                const size_t ljk_dense = triplet_key(l, j, k, naocc);
                if (l != i && i_j_k_to_ijk_.count(ljk_dense) && std::fabs((*F_lmo_)(l, i)) >= F_CUT) {
                    int ljk = i_j_k_to_ijk_[ljk_dense];

                    // S(a_{ijk}, a_{ljk})
                    std::vector<int> ljk_idx_list = index_list(triplet_ext_domain, lmotriplet_to_paos_[ljk]);
                    auto S_ijk_ljk = linalg::doublet(submatrix_rows(*S_ijk, ljk_idx_list), X_tno_[ljk], true, false);

                    SharedMatrix T_ljk;
                    if (write_amplitudes_) {
                        std::stringstream t_name;
                        t_name << "T " << (ljk);
                        T_ljk = std::make_shared<Matrix>(t_name.str(), n_tno_[ljk], n_tno_[ljk] * n_tno_[ljk]);
#pragma omp critical
                        T_ljk->load(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
                    } else {
                        T_ljk = T_iajbkc_[ljk];
                    }

                    // (a_{ljk}, b_{ljk}, c_{ljk}) -> (a_{ijk}, b_{ijk}, c_{ijk})
                    auto T_temp1 =
                        matmul_3d(triples_permuter(T_ljk, l, j, k), S_ijk_ljk, n_tno_[ljk], n_tno_[ijk]);
                    C_DAXPY(ntno_ijk * ntno_ijk * ntno_ijk, -(*F_lmo_)(l, i), &(*T_temp1)(0, 0), 1,
                            &(*R_ijk)(0, 0), 1);
                }
            }

            // => Update T3 Amplitudes (Jiang Eq. 112) <= //
            for (int a_ijk = 0; a_ijk < ntno_ijk; ++a_ijk) {
                for (int b_ijk = 0; b_ijk < ntno_ijk; ++b_ijk) {
                    for (int c_ijk = 0; c_ijk < ntno_ijk; ++c_ijk) {
                        (*T_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) -= (*R_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) /
                            ((*e_tno_[ijk])(a_ijk) + (*e_tno_[ijk])(b_ijk) + (*e_tno_[ijk])(c_ijk) 
                                - (*F_lmo_)(i, i) - (*F_lmo_)(j, j) - (*F_lmo_)(k, k));
                    }
                }
            }

            if (write_amplitudes_) {
#pragma omp critical
                T_ijk->save(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
            }
            
            R_iajbkc_rms[ijk] = R_ijk->rms();
        }

        // evaluate convergence
        e_prev = e_curr;
        e_ijk_old = e_ijk_;
        // Compute LCCSD(T) energy
        e_curr = compute_t_iteration_energy();

        double r_curr = *max_element(R_iajbkc_rms.begin(), R_iajbkc_rms.end());

        r_converged = fabs(r_curr) < options_.get_double("R_CONVERGENCE");
        e_converged = fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE");

        std::time_t time_stop = std::time(nullptr);

        outfile->Printf("  @LCCSD(T) iter %3d: %16.12f %10.3e %10.3e %8d\n", iteration, e_curr, e_curr - e_prev, r_curr, (int)time_stop - (int)time_start);

        iteration++;

        if (iteration > max_iteration + 1) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    timer_off("LCCSD(T) Iterations");

    return e_curr;
}

double DLPNOCCSD_T::compute_energy() {
    timer_on("DLPNO-CCSD(T)");

    // Run DLPNO-CCSD
    double e_dlpno_ccsd = DLPNOCCSD::compute_energy();

    psio_->open(PSIF_DLPNO_TRIPLES, PSIO_OPEN_NEW);

    // Clear CCSD integrals if CCSD(T) is last step
    if (algorithm_ == DLPNOMethod::CCSD_T) {
        K_mibj_.clear();
        J_ijmb_.clear();
        L_mibj_.clear();
        L_iajb_.clear();
        J_ikac_non_proj_.clear();
        K_iakc_non_proj_.clear();
        K_ivvv_.clear();
        Qma_ij_.clear();
        Qab_ij_.clear();
        i_Qk_ij_.clear();
        i_Qa_ij_.clear();
        i_Qk_t1_.clear();
        i_Qa_t1_.clear();
        S_pno_ij_kj_.clear();
        S_pno_ij_nn_.clear();
        S_pno_ij_mn_.clear();
        Fkc_.clear();
        Fai_.clear();
        Fab_.clear();
        T_n_ij_.clear();
    }

    print_header();

    double t_cut_tno_pre = options_.get_double("T_CUT_TNO_PRE");
    double t_cut_tno = options_.get_double("T_CUT_TNO");

    // Step 1: Perform the prescreening
    outfile->Printf("\n   Starting Triplet Prescreening...\n");
    outfile->Printf("     T_CUT_TNO set to %6.3e \n", t_cut_tno_pre);
    outfile->Printf("     T_CUT_DO  set to %6.3e \n", options_.get_double("T_CUT_DO_TRIPLES_PRE"));
    outfile->Printf("     T_CUT_MKN set to %6.3e \n\n", options_.get_double("T_CUT_MKN_TRIPLES_PRE"));

    triples_sparsity(true);
    tno_transform(t_cut_tno_pre);
    double E_T0_pre = compute_lccsd_t0();

    // Step 2: Compute DLPNO-CCSD(T0) energy with surviving triplets
    outfile->Printf("\n   Continuing computation with surviving triplets...\n");
    outfile->Printf("     Eliminated all triples with energy less than %6.3e Eh... \n\n", options_.get_double("T_CUT_TRIPLES_WEAK"));
    triples_sparsity(false);
    outfile->Printf("    * Energy Contribution From Screened Triplets: %.12f \n\n", de_lccsd_t_screened_);

    outfile->Printf("     T_CUT_TNO (re)set to %6.3e \n", t_cut_tno);
    outfile->Printf("     T_CUT_DO  (re)set to %6.3e \n", options_.get_double("T_CUT_DO_TRIPLES"));
    outfile->Printf("     T_CUT_MKN (re)set to %6.3e \n\n", options_.get_double("T_CUT_MKN_TRIPLES"));
    
    tno_transform(t_cut_tno);
    double E_T0 = compute_lccsd_t0();
    e_lccsd_t_ = e_lccsd_ + E_T0 + de_lccsd_t_screened_;

    outfile->Printf("    DLPNO-CCSD(T0) Correlation Energy: %16.12f \n", e_lccsd_ + E_T0 + de_lccsd_t_screened_);
    outfile->Printf("    * DLPNO-CCSD Contribution:         %16.12f \n", e_lccsd_);
    outfile->Printf("    * DLPNO-(T0) Contribution:         %16.12f \n", E_T0);
    outfile->Printf("    * Screened Triplets Contribution:  %16.12f \n\n", de_lccsd_t_screened_);


    // Step 3: Compute full DLPNO-CCSD(T) energy if NOT using T0 approximation

    // Full triples and higher methods always form the iterative-(T) amplitudes
    // used as their starting point; T0_APPROXIMATION applies only when
    // DLPNO-CCSD(T) is the final requested method.
    const bool t0_only = algorithm_ == DLPNOMethod::CCSD_T && options_.get_bool("T0_APPROXIMATION");
    if (!t0_only) {
        outfile->Printf("\n\n  ==> Computing Full Iterative (T) <==\n\n");

        sort_triplets(E_T0);

        const bool full_triples_follow = algorithm_ != DLPNOMethod::CCSD_T;
        const double t_cut_tno_full = options_.get_double("T_CUT_TNO_FULL");
        const double t_cut_tno_strong =
            full_triples_follow
                ? t_cut_tno_full
                : iterative_tno_cutoff(options_, "T_CUT_TNO_STRONG", "T_CUT_TNO_STRONG_SCALE");
        const double t_cut_tno_weak =
            full_triples_follow
                ? t_cut_tno_full
                : iterative_tno_cutoff(options_, "T_CUT_TNO_WEAK", "T_CUT_TNO_WEAK_SCALE");
        outfile->Printf("     T_CUT_TNO set to %6.3e for strong triplets \n", t_cut_tno_strong);
        outfile->Printf("     T_CUT_TNO set to %6.3e for weak triplets   \n\n", t_cut_tno_weak);

        tno_transform(t_cut_tno, true);
        estimate_memory();

        double E_T0_crude = compute_lccsd_t0(true);
        E_T_ = lccsd_t_iterations();
        double dE_T = E_T_ - E_T0_crude;

        outfile->Printf("\n");
        outfile->Printf("    DLPNO-CCSD(T0) energy at looser tolerance: %16.12f\n", E_T0_crude);
        outfile->Printf("    DLPNO-CCSD(T)  energy at looser tolerance: %16.12f\n", E_T_);
        outfile->Printf("    * Net Iterative (T) contribution:          %16.12f\n\n", dE_T);

        e_lccsd_t_ += dE_T;
    }

    double e_scf = reference_wavefunction_->energy();
    double e_ccsd_t_corr = e_lccsd_t_ + de_weak_ + de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
    double e_ccsd_t_total = e_scf + e_ccsd_t_corr;

    const std::string triples_energy_label = t0_only ? "CCSD(T0)" : "CCSD(T)";
    const std::string triples_correction_label = t0_only ? "(T0)" : "(T)";

    set_scalar_variable(triples_energy_label + " CORRELATION ENERGY", e_ccsd_t_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsd_t_corr);
    set_scalar_variable(triples_energy_label + " TOTAL ENERGY", e_ccsd_t_total);
    set_scalar_variable("CURRENT ENERGY", e_ccsd_t_total);

    // psivars for the selected perturbative-triples energy components
    set_scalar_variable(triples_correction_label + " CORRECTION ENERGY", e_lccsd_t_ - e_lccsd_);
    if (t0_only) {
        // Preserve the long-standing generic aliases used by CBS assembly,
        // scripts, and external harvesters while also publishing T0-specific
        // labels for callers that need to distinguish the approximation.
        set_scalar_variable("CCSD(T) CORRELATION ENERGY", e_ccsd_t_corr);
        set_scalar_variable("CCSD(T) TOTAL ENERGY", e_ccsd_t_total);
        set_scalar_variable("(T) CORRECTION ENERGY", e_lccsd_t_ - e_lccsd_);
    }
    set_scalar_variable("DLPNO SEMICANONICAL (T0) ENERGY", E_T0 + de_lccsd_t_screened_);
    set_scalar_variable("DLPNO SCREENED TRIPLETS ENERGY", de_lccsd_t_screened_);

    // Iterative CCSDT needs the converged T3 amplitudes in memory. If the
    // preceding (T) stage used disk storage, hydrate that final generation
    // before deleting its PSIO file.
    const bool ccsdt_follows = algorithm_ != DLPNOMethod::CCSD_T;
    if (ccsdt_follows && write_amplitudes_) {
        const int n_lmo_triplets = ijk_to_i_j_k_.size();
#pragma omp parallel for schedule(dynamic, 1)
        for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
            const int ntno_ijk = n_tno_[ijk];
            std::stringstream t_name;
            t_name << "T " << ijk;
            auto T_ijk = std::make_shared<Matrix>(t_name.str(), ntno_ijk, ntno_ijk * ntno_ijk);
#pragma omp critical
            T_ijk->load(psio_, PSIF_DLPNO_TRIPLES, Matrix::SubBlocks);
            T_iajbkc_[ijk] = std::move(T_ijk);
        }
    }

    print_results();
    
    psio_->close(PSIF_DLPNO_TRIPLES, 0);

    timer_off("DLPNO-CCSD(T)");

    return e_ccsd_t_total;
}

void DLPNOCCSD_T::print_results() {
    double e_dlpno_ccsd = e_lccsd_ + de_weak_ + de_lmp2_eliminated_ + de_pno_total_ + de_dipole_;
    double e_total = e_lccsd_t_ + de_weak_ + de_lmp2_eliminated_ + de_pno_total_ + de_dipole_;
    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-CCSD(T) Correlation Energy: %16.12f \n", e_total);
    outfile->Printf("    DLPNO-CCSD Contribution:              %16.12f \n", e_dlpno_ccsd);
    outfile->Printf("    DLPNO-(T) Contribution:               %16.12f \n", e_lccsd_t_ - e_lccsd_ - de_lccsd_t_screened_);
    outfile->Printf("    Screened Triplets Contribution:       %16.12f \n", de_lccsd_t_screened_);
    outfile->Printf("\n\n  @Total DLPNO-CCSD(T) Energy: %16.12f \n",
                    variables_["SCF TOTAL ENERGY"] + de_weak_ + de_lmp2_eliminated_ + e_lccsd_t_ + de_pno_total_ + de_dipole_);
    outfile->Printf("    *** Andy Jiang... FOR THREEEEEEEEEEE!!!\n\n");
}

#ifdef USING_Einsums

DLPNOCCSDT::DLPNOCCSDT(SharedWavefunction ref_wfn, Options& options)
    : DLPNOCCSD_T(ref_wfn, options), disk_ints_(options.get_bool("DLPNO_CCSDT_DISK_INTS")) {}
DLPNOCCSDT::~DLPNOCCSDT() {}

Tensor<double, 3> DLPNOCCSDT::matmul_3d_einsums(const Tensor<double, 3> &A, const SharedMatrix &X, int dim_old, int dim_new) {
    /* Transform all three indices of a rank-3 tensor: A'[i,j,k] = A[I,J,K] X[i,I] X[j,J] X[k,K]. */

    // TODO: Change this into a TensorView
    Tensor<double, 2> Xview("Xview", dim_new, dim_old);
    ::memcpy(Xview.data(), X->get_pointer(), dim_new * dim_old * sizeof(double));

    // Move each next contracted old index to the contiguous final axis.  The
    // three contractions are therefore GEMMs, and nested lifetimes keep no
    // more than the current and next rank-3 work buffers alive.
    Tensor<double, 3> work("rank3_transform_1", dim_old, dim_old, dim_new);
    einsum(0.0, Indices{index::I, index::J, index::k}, &work, 1.0,
           Indices{index::I, index::J, index::K}, A,
           Indices{index::k, index::K}, Xview);

    {
        Tensor<double, 3> contiguous("rank3_transform_perm_2", dim_old,
                                      dim_new, dim_old);
        permute(Indices{index::I, index::k, index::J}, &contiguous,
                Indices{index::I, index::J, index::k}, work);
        work = Tensor<double, 3>("released", 0, 0, 0);
        Tensor<double, 3> next("rank3_transform_2", dim_old, dim_new, dim_new);
        einsum(0.0, Indices{index::I, index::k, index::j}, &next, 1.0,
               Indices{index::I, index::k, index::J}, contiguous,
               Indices{index::j, index::J}, Xview);
        work = std::move(next);
    }

    Tensor<double, 3> transformed("rank3_transform_3", dim_new, dim_new,
                                   dim_new);
    einsum(0.0, Indices{index::i, index::k, index::j}, &transformed, 1.0,
           Indices{index::i, index::I}, Xview,
           Indices{index::I, index::k, index::j}, work);
    work = Tensor<double, 3>("released", 0, 0, 0);

    Tensor<double, 3> result("rank3_transform_result", dim_new, dim_new,
                              dim_new);
    permute(Indices{index::i, index::j, index::k}, &result,
            Indices{index::i, index::k, index::j}, transformed);
    return result;
}

Tensor<double, 3> DLPNOCCSDT::matmul_3d_index(const Tensor<double, 3> &A, const SharedMatrix &X, int index) {
    /* Transform one selected index of a rank-3 tensor: A'[i,j,k] = A[I,j,k] X[i,I]. */

    // TODO: Change this into a TensorView
    int dim_new = X->rowspi(0);
    int dim_old = X->colspi(0);
    Tensor<double, 2> Xview("Xview", dim_new, dim_old);
    ::memcpy(Xview.data(), X->get_pointer(), dim_new * dim_old * sizeof(double));

    Tensor<double, 3> A_new;

    if (index == 0) {
        A_new = Tensor<double, 3>("A_new", dim_new, A.dim(1), A.dim(2));
        einsum(0.0, Indices{index::i, index::J, index::K}, &A_new, 1.0, Indices{index::I, index::J, index::K}, A, Indices{index::i, index::I}, Xview);
    } else if (index == 1) {
        Tensor<double, 3> A_temp("A_temp", A.dim(1), A.dim(0), A.dim(2));
        permute(Indices{index::J, index::I, index::K}, &A_temp, Indices{index::I, index::J, index::K}, A);

        Tensor<double, 3> A_temp2("A_temp2", dim_new, A.dim(0), A.dim(2));
        einsum(0.0, Indices{index::j, index::I, index::K}, &A_temp2, 1.0, Indices{index::J, index::I, index::K}, A_temp, Indices{index::j, index::J}, Xview);

        A_new = Tensor<double, 3>("A_new", A.dim(0), dim_new, A.dim(2));
        permute(Indices{index::I, index::j, index::K}, &A_new, Indices{index::j, index::I, index::K}, A_temp2);
    } else if (index == 2) {
        A_new = Tensor<double, 3>("A_new", A.dim(0), A.dim(1), dim_new);
        einsum(0.0, Indices{index::I, index::J, index::k}, &A_new, 1.0, Indices{index::I, index::J, index::K}, A, Indices{index::k, index::K}, Xview);
    } else {
        throw PSIEXCEPTION("Index out of bounds for 3D Tensor");
    }

    return A_new;
}

Tensor<double, 3> DLPNOCCSDT::triples_permuter_einsums(const Tensor<double, 3> &X, int i, int j, int k, bool reverse) {
    // Generate an occupied/virtual-column permutation from the stored i <= j <= k tensor.
    // The direct permutation operator is defined by manuscript Eqs. (81)-(86); reverse=true
    // applies the conjugate operator of Eqs. (91)-(96).
    Tensor<double, 3> Xperm = X;

    if (!reverse) {
        if (i <= k && k <= j && i <= j) {
            permute(Indices{index::a, index::b, index::c}, &Xperm, Indices{index::a, index::c, index::b}, X);
        } else if (j <= i && i <= k && j <= k) {
            permute(Indices{index::a, index::b, index::c}, &Xperm, Indices{index::b, index::a, index::c}, X);
        } else if (j <= k && k <= i && j <= i) {
            permute(Indices{index::a, index::b, index::c}, &Xperm, Indices{index::b, index::c, index::a}, X);
        } else if (k <= i && i <= j && k <= j) {
            permute(Indices{index::a, index::b, index::c}, &Xperm, Indices{index::c, index::a, index::b}, X);
        } else if (k <= j && j <= i && k <= i) {
            permute(Indices{index::a, index::b, index::c}, &Xperm, Indices{index::c, index::b, index::a}, X);
        }
    } else {
        if (i <= k && k <= j && i <= j) {
            permute(Indices{index::a, index::c, index::b}, &Xperm, Indices{index::a, index::b, index::c}, X);
        } else if (j <= i && i <= k && j <= k) {
            permute(Indices{index::b, index::a, index::c}, &Xperm, Indices{index::a, index::b, index::c}, X);
        } else if (j <= k && k <= i && j <= i) {
            permute(Indices{index::b, index::c, index::a}, &Xperm, Indices{index::a, index::b, index::c}, X);
        } else if (k <= i && i <= j && k <= j) {
            permute(Indices{index::c, index::a, index::b}, &Xperm, Indices{index::a, index::b, index::c}, X);
        } else if (k <= j && j <= i && k <= i) {
            permute(Indices{index::c, index::b, index::a}, &Xperm, Indices{index::a, index::b, index::c}, X);
        }
    }

    return Xperm;
}

Tensor<double, 3> DLPNOCCSDT::triples_spin_summation(const Tensor<double, 3> &X) {

    // Successive spin summation of an orbital triples tensor. This is the six-permutation
    // specialization of the procedure described by Matthews and Stanton immediately before
    // their Eqs. (23)-(25), J. Chem. Phys. 142, 064108 (2015), DOI: 10.1063/1.4907278.

    Tensor<double, 3> Xnew = X;
    Xnew *= 4.0;
    Tensor<double, 3> Xtemp = triples_permuter_einsums(X, 0, 2, 1);
    Xtemp += triples_permuter_einsums(X, 2, 1, 0);
    Xtemp += triples_permuter_einsums(X, 1, 0, 2);
    Xtemp *= 2.0;
    Xnew -= Xtemp;
    Xnew += triples_permuter_einsums(X, 2, 0, 1);
    Xnew += triples_permuter_einsums(X, 1, 2, 0);

    return Xnew;
}

Tensor<double, 3> DLPNOCCSDT::triples_spin_desummation(const Tensor<double, 3> &X) {

    // Matthews and Stanton Eq. (27), using their fully spin-summed three-body form:
    // t = (1/4) X - (1/12) P_(jki) X - (1/12) P_(kij) X. The inverse is not unique
    // because the spin-summation transformation is singular; this sparse pseudoinverse
    // recovers an equivalent orbital tensor without retaining a second amplitude copy.

    Tensor<double, 3> Xnew = X;
    Xnew *= 3.0;
    Xnew -= triples_permuter_einsums(X, 2, 0, 1);
    Xnew -= triples_permuter_einsums(X, 1, 2, 0);
    Xnew *= 1.0 / 12.0;

    return Xnew;
}

void DLPNOCCSDT::print_header() {
    double t_cut_tno = options_.get_double("T_CUT_TNO");
    const double t_cut_tno_full = options_.get_double("T_CUT_TNO_FULL");

    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                    DLPNO-CCSDT                \n");
    outfile->Printf("                   by Andy Jiang               \n");
    outfile->Printf("              DOI: 10.1021/acs.jctc.4c01716   \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Inherited CCSD parameters are reported in the preceding DLPNO-CCSD header.\n");
    outfile->Printf("  Full-triples and inherited triples parameters:\n");
    outfile->Printf("    T_CUT_TNO                  = %6.3e \n", t_cut_tno);
    outfile->Printf("    T_CUT_TNO_FULL             = %6.3e \n", t_cut_tno_full);
    outfile->Printf("    T_CUT_TNO_STRONG           = %6.3e \n", t_cut_tno_full);
    outfile->Printf("    T_CUT_TNO_WEAK             = %6.3e \n", t_cut_tno_full);
    outfile->Printf("    T_CUT_TNO_PRE              = %6.3e \n", options_.get_double("T_CUT_TNO_PRE"));
    outfile->Printf("    T_CUT_TRIPLES_WEAK         = %6.3e \n", options_.get_double("T_CUT_TRIPLES_WEAK"));
    outfile->Printf("    T_CUT_DO_TRIPLES           = %6.3e \n", options_.get_double("T_CUT_DO_TRIPLES"));
    outfile->Printf("    T_CUT_DO_TRIPLES_PRE       = %6.3e \n", options_.get_double("T_CUT_DO_TRIPLES_PRE"));
    outfile->Printf("    T_CUT_MKN_TRIPLES          = %6.3e \n", options_.get_double("T_CUT_MKN_TRIPLES"));
    outfile->Printf("    T_CUT_MKN_TRIPLES_PRE      = %6.3e \n", options_.get_double("T_CUT_MKN_TRIPLES_PRE"));
    outfile->Printf("    F_CUT_T                    = %6.3e \n", options_.get_double("F_CUT_T"));
    outfile->Printf("    T_CUT_ITER                 = %6.3e \n", options_.get_double("T_CUT_ITER"));
    outfile->Printf("    MIN_TNOS                   = %6d   \n", options_.get_int("MIN_TNOS"));
    outfile->Printf("    TRIPLES_MAX_WEAK_PAIRS     = %6d   \n", options_.get_int("TRIPLES_MAX_WEAK_PAIRS"));
    outfile->Printf("    TRIPLES_MICROITERATIONS    = %6d   \n", options_.get_int("DLPNO_TRIPLES_MICROITERATIONS"));
    outfile->Printf("    TRIPLES_DAMPING_RATIO      = %6.3e \n", options_.get_double("TRIPLES_DAMPING_RATIO"));
    outfile->Printf("    E_CONVERGENCE              = %6.3e \n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE              = %6.3e \n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("    DLPNO_MAXITER              = %6d   \n", options_.get_int("DLPNO_MAXITER"));
    outfile->Printf("    DIIS_MAX_VECS              = %6d   \n", options_.get_int("DIIS_MAX_VECS"));
    outfile->Printf("    DLPNO_TOGGLE_MEMORY        = %6s   \n",
                    options_.get_bool("DLPNO_TOGGLE_MEMORY") ? "TRUE" : "FALSE");
    outfile->Printf("    DLPNO_CCSDT_DISK_INTS      = %6s   \n", disk_ints_ ? "TRUE" : "FALSE");
    outfile->Printf("    WRITE_TRIPLES_INTERMEDIATES= %6s   \n", write_intermediates_ ? "TRUE" : "FALSE");
    outfile->Printf("    WRITE_TRIPLES_AMPLITUDES   = %6s   \n\n", write_amplitudes_ ? "TRUE" : "FALSE");
}

void DLPNOCCSDT::estimate_memory() {
    const size_t n_lmo_triplets = ijk_to_i_j_k_.size();
    const size_t n_lmo_pairs = ij_to_i_j_.size();
    const size_t nthreads = static_cast<size_t>(nthread_);

    size_t tno_basis_memory = 0;
    size_t triples_amplitude_memory = 0;
    size_t perturbative_triples_memory = 0;
    size_t inexpensive_df_memory = 0;
    size_t k_ovov_cache_memory = 0;
    size_t disk_eligible_df_memory = 0;
    size_t triples_iteration_memory = 0;
    size_t triples_iteration_resident_memory = 0;
    size_t integral_workspace_per_thread = 0;
    size_t contraction_workspace_per_thread = 0;

    // The per-thread expressions below are conservative upper bounds built from the tensors
    // simultaneously scoped in compute_integrals() and the R1/R2/R3 residual builders. They
    // intentionally favor a safe estimate over relying only on the leading asymptotic term.
    for (size_t ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        const size_t naux = lmotriplet_to_ribfs_[ijk].size();
        const size_t nlmo = lmotriplet_to_lmos_[ijk].size();
        const size_t npao = lmotriplet_to_paos_[ijk].size();
        const size_t ntno = n_tno_[ijk];
        const size_t ntno2 = ntno * ntno;
        const size_t ntno3 = ntno2 * ntno;
        const size_t ntno4 = ntno3 * ntno;
        const size_t nlmo2 = nlmo * nlmo;

        // X_tno, e_tno, and the triples amplitude produced by the preceding (T) calculation.
        tno_basis_memory += npao * ntno + ntno;
        triples_amplitude_memory += ntno3;

        // W and V remain resident after the initial (T) calculation unless disk-backed.
        if (!write_intermediates_) perturbative_triples_memory += 2 * ntno3;

        // q_io/q_jo/q_ko and q_iv/q_jv/q_kv are always retained in core. q_ov and q_vv
        // are the two tensors controlled by DLPNO_CCSDT_DISK_INTS.
        inexpensive_df_memory += 3 * naux * (nlmo + ntno);
        k_ovov_cache_memory += 3 * ntno2;
        disk_eligible_df_memory += naux * nlmo * ntno + naux * ntno2;

        // T_n_ijk, the Einsums T3 clone, contravariant T3, and the R3 matrix.
        triples_iteration_memory += nlmo * ntno + 3 * ntno3;
        // R3 is iteration-local; the other three quantities persist into a following
        // CCSDT(Q)/CCSDTQ stage and therefore belong in the inherited resident estimate.
        triples_iteration_resident_memory += nlmo * ntno + 2 * ntno3;

        const size_t integral_workspace =
            3 * naux * nlmo + 3 * naux * npao + naux * nlmo * ntno + naux * ntno2 +
            2 * naux * naux + nlmo * ntno + npao * npao +
            std::max(naux * nlmo * ntno, naux * ntno2);
        integral_workspace_per_thread = std::max(integral_workspace_per_thread, integral_workspace);

        // One contiguous U3 column permutation is formed at a time so the
        // subsequent contraction is always a GEMV.
        const size_t r1_workspace = ntno3 + 3 * ntno2 + 2 * ntno;
        const size_t r2_workspace =
            ntno3 + 14 * nlmo * ntno + 14 * ntno2 +
            6 * nlmo * ntno2;

        // R3 includes the long-permutation W tensors, short-permutation V tensors,
        // rho/chi intermediates, T_lm/K_ldme blocks, and T1-dressed DF work arrays.
        const size_t r3_workspace =
            3 * ntno4 + 5 * nlmo2 * ntno2 + 24 * ntno3 + 14 * nlmo * ntno2 +
            6 * nlmo2 * ntno + 3 * nlmo2 + 3 * naux * ntno2 +
            4 * naux * nlmo * ntno + 2 * naux * nlmo2 + 6 * naux * (nlmo + ntno) +
            static_cast<size_t>(basisset_->nbf()) * ntno;
        contraction_workspace_per_thread =
            std::max(contraction_workspace_per_thread, std::max(r1_workspace, std::max(r2_workspace, r3_workspace)));
    }

    size_t singles_words = 0;
    size_t doubles_words = 0;
    size_t projected_pair_singles_memory = 0;
    for (size_t i = 0; i < i_j_to_ij_.size(); ++i) {
        const int ii = i_j_to_ij_[i][i];
        singles_words += n_pno_[ii];
    }
    for (size_t ij = 0; ij < n_lmo_pairs; ++ij) {
        const size_t npno = n_pno_[ij];
        doubles_words += npno * npno;
        projected_pair_singles_memory += lmopair_to_lmos_[ij].size() * npno;
    }

    // T1, R1, R2, Rn2, T_n_ij, and the thread-private R1/R2 accumulation buffers.
    const size_t lower_rank_iteration_memory =
        2 * singles_words + 2 * doubles_words + projected_pair_singles_memory +
        nthreads * (singles_words + doubles_words);

    // The in-core DIIS manager retains a solution and error vector for every slot. Include
    // one additional solution/error pair for the flattened vectors used during extrapolation.
    const size_t iteration_vector_words = singles_words + doubles_words + triples_amplitude_memory;
    const size_t diis_memory =
        2 * (static_cast<size_t>(options_.get_int("DIIS_MAX_VECS")) + 1) * iteration_vector_words;

    // rho_dbck_contribution() returns one diagonal-PNO cube per occupied orbital.
    // It packs B_Q^{me} once and streams one npno^2 integral slice, rather than
    // retaining the former nlmo^2*npno^2 g_menf tensor.
    size_t rho_intermediate_memory = 0;
    size_t rho_workspace_per_thread = 0;
    for (size_t i = 0; i < i_j_to_ij_.size(); ++i) {
        const int ii = i_j_to_ij_[i][i];
        const size_t npno = n_pno_[ii];
        const size_t nlmo = lmopair_to_lmos_[ii].size();
        const size_t naux = lmopair_to_ribfs_[ii].size();
        const size_t npno2 = npno * npno;
        const size_t npno3 = npno2 * npno;
        rho_intermediate_memory += npno3;
        const size_t packed_df = naux * nlmo * npno;
        const size_t packing_peak = 2 * packed_df;
        const size_t contraction_peak = packed_df + 3 * npno3 + npno2;
        rho_workspace_per_thread = std::max(
            rho_workspace_per_thread, std::max(packing_peak, contraction_peak));
    }
    contraction_workspace_per_thread = std::max(contraction_workspace_per_thread, rho_workspace_per_thread);

    size_t in_core_disk_eligible_df_memory = disk_ints_ ? 0 : disk_eligible_df_memory;

    auto memory_peaks = [&]() {
        const size_t common_memory = ccsd_baseline_memory_doubles_ + tno_basis_memory +
                                     triples_amplitude_memory + perturbative_triples_memory +
                                     inexpensive_df_memory + k_ovov_cache_memory +
                                     in_core_disk_eligible_df_memory;
        const size_t integral_peak = common_memory + nthreads * integral_workspace_per_thread;
        const size_t iteration_resident = common_memory + triples_iteration_memory +
                                          lower_rank_iteration_memory + diis_memory;
        const size_t ccsd_residual_peak = iteration_resident + ccsd_iteration_workspace_doubles_;
        const size_t triples_residual_peak = iteration_resident + rho_intermediate_memory +
                                             nthreads * contraction_workspace_per_thread;
        const size_t iteration_peak = std::max(ccsd_residual_peak, triples_residual_peak);
        return std::make_pair(integral_peak, iteration_peak);
    };

    const double DOUBLES_TO_GB = 1.0e-9 * sizeof(double);
    const double BYTES_TO_GB = 1.0e-9;

    auto print_estimate = [&](const char* title) {
        const auto [integral_peak, iteration_peak] = memory_peaks();
        auto print_memory_line = [&](const std::string& label, size_t words) {
            outfile->Printf("    %-46s : %8.3f [GB]\n", label.c_str(), words * DOUBLES_TO_GB);
        };
        outfile->Printf("\n  ==> %s <==\n\n", title);
        print_memory_line("Inherited DLPNO-CCSD baseline estimate", ccsd_baseline_memory_doubles_);
        print_memory_line("Inherited CCSD iteration workspace", ccsd_iteration_workspace_doubles_);
        print_memory_line("Inherited DLPNO-CCSD peak estimate", ccsd_peak_memory_doubles_);
        print_memory_line("TNO transforms and orbital energies", tno_basis_memory);
        print_memory_line("(T) triples amplitudes", triples_amplitude_memory);
        print_memory_line("Retained (T) W/V intermediates", perturbative_triples_memory);
        print_memory_line("Always-in-core CCSDT DF integrals", inexpensive_df_memory);
        print_memory_line("Shared T3-to-R1/R2 exchange cache", k_ovov_cache_memory);
        print_memory_line("In-core disk-eligible CCSDT integrals", in_core_disk_eligible_df_memory);
        print_memory_line("CCSDT amplitudes/intermediates", triples_iteration_memory);
        print_memory_line("Lower-rank amplitudes/residual buffers", lower_rank_iteration_memory);
        print_memory_line("In-core DIIS upper bound", diis_memory);
        print_memory_line("rho_dbck cross-domain intermediates", rho_intermediate_memory);
        print_memory_line("Integral workspace per thread (" + std::to_string(nthreads) + ")",
                          integral_workspace_per_thread);
        print_memory_line("Contraction workspace per thread (" + std::to_string(nthreads) + ")",
                          contraction_workspace_per_thread);
        print_memory_line("Estimated integral-build peak", integral_peak);
        print_memory_line("Estimated CCSDT-iteration peak", iteration_peak);
        outfile->Printf("    %-46s : %8.3f [GB]\n\n", "Total memory given", memory_ * BYTES_TO_GB);
    };

    print_estimate("DLPNO-CCSDT Memory Requirements");

    auto [integral_peak, iteration_peak] = memory_peaks();
    size_t required_memory = std::max(integral_peak, iteration_peak);
    if (toggle_memory_ && !disk_ints_ && required_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total required memory is more than 90%% of available memory.\n");
        outfile->Printf("    Switching (Q_{ijk}|m_{ijk} a_{ijk}) and (Q_{ijk}|a_{ijk} b_{ijk}) to disk...\n");
        disk_ints_ = true;
        in_core_disk_eligible_df_memory = 0;
        std::tie(integral_peak, iteration_peak) = memory_peaks();
        required_memory = std::max(integral_peak, iteration_peak);
        print_estimate("Updated DLPNO-CCSDT Memory Requirements");
    }

    if (toggle_memory_ && required_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total required memory remains more than 90%% of available memory after all safe toggles.\n");
        throw PSIEXCEPTION("Too little memory given for the DLPNO-CCSDT algorithm.");
    }

    ccsdt_resident_memory_doubles_ =
        ccsd_baseline_memory_doubles_ + tno_basis_memory + triples_amplitude_memory +
        perturbative_triples_memory + inexpensive_df_memory + k_ovov_cache_memory +
        in_core_disk_eligible_df_memory +
        triples_iteration_resident_memory;
    ccsdt_iteration_workspace_doubles_ =
        std::max(ccsd_iteration_workspace_doubles_,
                 rho_intermediate_memory + nthreads * contraction_workspace_per_thread);

    if (disk_ints_) {
        outfile->Printf("    Writing the largest TNO-basis DF integrals to disk.\n\n");
    } else {
        outfile->Printf("    Keeping all TNO-basis DF integrals in core.\n\n");
    }
}

void DLPNOCCSDT::compute_integrals() {

    size_t n_lmo_triplets = ijk_to_i_j_k_.size();

    // Three-center quantities
    q_io_.resize(n_lmo_triplets);
    q_jo_.resize(n_lmo_triplets);
    q_ko_.resize(n_lmo_triplets);

    q_iv_.resize(n_lmo_triplets);
    q_jv_.resize(n_lmo_triplets);
    q_kv_.resize(n_lmo_triplets);

    q_ov_.resize(n_lmo_triplets);
    q_vv_.resize(n_lmo_triplets);
    K_ovov_tno_cache_.resize(n_lmo_triplets);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];
        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];

        int ntno_ijk = n_tno_[ijk];

        if (ntno_ijk == 0) continue;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // => Compute all necessary integrals <= //

        // number of auxiliary functions in the triplet domain
        const int naux_ijk = lmotriplet_to_ribfs_[ijk].size();
        // number of LMOs in the triplet domain
        const int nlmo_ijk = lmotriplet_to_lmos_[ijk].size();
        // number of PAOs in the triplet domain (before removing linear dependencies)
        const int npao_ijk = lmotriplet_to_paos_[ijk].size();

        auto q_io = std::make_shared<Matrix>("(Q_ijk | m i)", naux_ijk, nlmo_ijk);
        auto q_jo = std::make_shared<Matrix>("(Q_ijk | m j)", naux_ijk, nlmo_ijk);
        auto q_ko = std::make_shared<Matrix>("(Q_ijk | m k)", naux_ijk, nlmo_ijk);

        auto q_iv = std::make_shared<Matrix>("(Q_ijk | i a)", naux_ijk, npao_ijk);
        auto q_jv = std::make_shared<Matrix>("(Q_ijk | j b)", naux_ijk, npao_ijk);
        auto q_kv = std::make_shared<Matrix>("(Q_ijk | k c)", naux_ijk, npao_ijk);

        std::stringstream q_ov_name;
        q_ov_name << "(Q_ijk | m a) " << (ijk);
        std::stringstream q_vv_name;
        q_vv_name << "(Q_ijk | a b) " << (ijk);

        auto q_ov = std::make_shared<Matrix>(q_ov_name.str(), naux_ijk, nlmo_ijk * ntno_ijk);
        auto q_vv = std::make_shared<Matrix>(q_vv_name.str(), naux_ijk, ntno_ijk * ntno_ijk);

        for (int q_ijk = 0; q_ijk < naux_ijk; q_ijk++) {
            const int q = lmotriplet_to_ribfs_[ijk][q_ijk];
            const int centerq = ribasis_->function_to_center(q);

            // Cheaper Integrals
            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                (*q_io)(q_ijk, l_ijk) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][i], riatom_to_lmos_ext_dense_[centerq][l]);
                (*q_jo)(q_ijk, l_ijk) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][j], riatom_to_lmos_ext_dense_[centerq][l]);
                (*q_ko)(q_ijk, l_ijk) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][k], riatom_to_lmos_ext_dense_[centerq][l]);
            }


            for (int u_ijk = 0; u_ijk < npao_ijk; ++u_ijk) {
                int u = lmotriplet_to_paos_[ijk][u_ijk];
                (*q_iv)(q_ijk, u_ijk) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][i], riatom_to_paos_ext_dense_[centerq][u]);
                (*q_jv)(q_ijk, u_ijk) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][j], riatom_to_paos_ext_dense_[centerq][u]);
                (*q_kv)(q_ijk, u_ijk) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][k], riatom_to_paos_ext_dense_[centerq][u]);
            }

            // More expensive integrals
            auto q_ov_tmp = std::make_shared<Matrix>(nlmo_ijk, npao_ijk);

            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                for (int u_ijk = 0; u_ijk < npao_ijk; ++u_ijk) {
                    int u = lmotriplet_to_paos_[ijk][u_ijk];
                    (*q_ov_tmp)(l_ijk, u_ijk) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][l], riatom_to_paos_ext_dense_[centerq][u]);
                }
            }
            q_ov_tmp = linalg::doublet(q_ov_tmp, X_tno_[ijk], false, false);
            ::memcpy(&(*q_ov)(q_ijk, 0), &(*q_ov_tmp)(0, 0), nlmo_ijk * ntno_ijk * sizeof(double));

            auto q_vv_tmp = std::make_shared<Matrix>(npao_ijk, npao_ijk);

            for (int u_ijk = 0; u_ijk < npao_ijk; ++u_ijk) {
                int u = lmotriplet_to_paos_[ijk][u_ijk];
                for (int v_ijk = 0; v_ijk < npao_ijk; ++v_ijk) {
                    int v = lmotriplet_to_paos_[ijk][v_ijk];
                    int uv_idx = riatom_to_pao_pairs_dense_[centerq][u][v];
                    if (uv_idx == -1) continue;
                    (*q_vv_tmp)(u_ijk, v_ijk) = (*qab_[q])(uv_idx, 0);
                } // end v_ijk
            } // end u_ijk
            q_vv_tmp = linalg::triplet(X_tno_[ijk], q_vv_tmp, X_tno_[ijk], true, false, false);
            ::memcpy(&(*q_vv)(q_ijk, 0), &(*q_vv_tmp)(0, 0), ntno_ijk * ntno_ijk * sizeof(double));

        } // end q_ijk

        // Contract Intermediates
        q_iv = linalg::doublet(q_iv, X_tno_[ijk]);
        q_jv = linalg::doublet(q_jv, X_tno_[ijk]);
        q_kv = linalg::doublet(q_kv, X_tno_[ijk]);

        // Multiply by (P|Q)^{-1/2}
        auto metric_inverse_sqrt =
            submatrix_rows_and_cols(*full_metric_, lmotriplet_to_ribfs_[ijk], lmotriplet_to_ribfs_[ijk]);
        metric_inverse_sqrt->power(-0.5, 1.0e-14);

        q_io = linalg::doublet(metric_inverse_sqrt, q_io);
        q_jo = linalg::doublet(metric_inverse_sqrt, q_jo);
        q_ko = linalg::doublet(metric_inverse_sqrt, q_ko);
        q_iv = linalg::doublet(metric_inverse_sqrt, q_iv);
        q_jv = linalg::doublet(metric_inverse_sqrt, q_jv);
        q_kv = linalg::doublet(metric_inverse_sqrt, q_kv);
        q_ov = linalg::doublet(metric_inverse_sqrt, q_ov);
        q_vv = linalg::doublet(metric_inverse_sqrt, q_vv);

        // linalg::doublet() names its result "T". Restore the unique names used
        // as PSIO keys before the disk-backed integral tensors are saved.
        q_ov->set_name(q_ov_name.str());
        q_vv->set_name(q_vv_name.str());
        metric_inverse_sqrt.reset();

        q_io_[ijk] = Tensor<double, 2>("(Q_ijk | m i)", naux_ijk, nlmo_ijk);
        q_jo_[ijk] = Tensor<double, 2>("(Q_ijk | m j)", naux_ijk, nlmo_ijk);
        q_ko_[ijk] = Tensor<double, 2>("(Q_ijk | m k)", naux_ijk, nlmo_ijk);

        q_iv_[ijk] = Tensor<double, 2>("(Q_ijk | i a)", naux_ijk, ntno_ijk);
        q_jv_[ijk] = Tensor<double, 2>("(Q_ijk | j b)", naux_ijk, ntno_ijk);
        q_kv_[ijk] = Tensor<double, 2>("(Q_ijk | k c)", naux_ijk, ntno_ijk);

        q_ov_[ijk] = Tensor<double, 3>(q_ov->name(), naux_ijk, nlmo_ijk, ntno_ijk);
        q_vv_[ijk] = Tensor<double, 3>(q_vv->name(), naux_ijk, ntno_ijk, ntno_ijk);

        ::memcpy(q_io_[ijk].data(), q_io->get_pointer(), naux_ijk * nlmo_ijk * sizeof(double));
        ::memcpy(q_jo_[ijk].data(), q_jo->get_pointer(), naux_ijk * nlmo_ijk * sizeof(double));
        ::memcpy(q_ko_[ijk].data(), q_ko->get_pointer(), naux_ijk * nlmo_ijk * sizeof(double));
        ::memcpy(q_iv_[ijk].data(), q_iv->get_pointer(), naux_ijk * ntno_ijk * sizeof(double));
        ::memcpy(q_jv_[ijk].data(), q_jv->get_pointer(), naux_ijk * ntno_ijk * sizeof(double));
        ::memcpy(q_kv_[ijk].data(), q_kv->get_pointer(), naux_ijk * ntno_ijk * sizeof(double));
        ::memcpy(q_ov_[ijk].data(), q_ov->get_pointer(), naux_ijk * nlmo_ijk * ntno_ijk * sizeof(double));
        ::memcpy(q_vv_[ijk].data(), q_vv->get_pointer(), naux_ijk * ntno_ijk * ntno_ijk * sizeof(double));

        // These three exchange matrices are used independently by the T3 -> R1
        // and T3 -> R2 passes, but neither the DF factors nor the TNO basis changes
        // during an iteration. Build them once here and share the cache.
        auto& K_cache = K_ovov_tno_cache_[ijk];
        K_cache[0] = Tensor<double, 2>("K_ivjv", ntno_ijk, ntno_ijk);
        K_cache[1] = Tensor<double, 2>("K_ivkv", ntno_ijk, ntno_ijk);
        K_cache[2] = Tensor<double, 2>("K_jvkv", ntno_ijk, ntno_ijk);
        einsum(0.0, Indices{index::a, index::b}, &K_cache[0], 1.0,
               Indices{index::Q, index::a}, q_iv_[ijk],
               Indices{index::Q, index::b}, q_jv_[ijk]);
        einsum(0.0, Indices{index::a, index::b}, &K_cache[1], 1.0,
               Indices{index::Q, index::a}, q_iv_[ijk],
               Indices{index::Q, index::b}, q_kv_[ijk]);
        einsum(0.0, Indices{index::a, index::b}, &K_cache[2], 1.0,
               Indices{index::Q, index::a}, q_jv_[ijk],
               Indices{index::Q, index::b}, q_kv_[ijk]);

        if (disk_ints_) {
#pragma omp critical
            q_ov->save(psio_.get(), PSIF_DLPNO_QIA_TNO, Matrix::SubBlocks);

#pragma omp critical
            q_vv->save(psio_.get(), PSIF_DLPNO_QAB_TNO, Matrix::ThreeIndexLowerTriangle);

            q_ov_[ijk] = Tensor<double, 3>(q_ov->name(), 0, 0, 0);
            q_vv_[ijk] = Tensor<double, 3>(q_vv->name(), 0, 0, 0);
        }
    } // end ijk
}

void DLPNOCCSDT::load_qia_tno(int ijk) {
    if (disk_ints_) {
        int naux_ijk = lmotriplet_to_ribfs_[ijk].size();
        int nlmo_ijk = lmotriplet_to_lmos_[ijk].size();
        int ntno_ijk = n_tno_[ijk];

        std::stringstream q_ov_name;
        q_ov_name << "(Q_ijk | m a) " << (ijk);

        auto q_ov = std::make_shared<Matrix>(q_ov_name.str(), naux_ijk, nlmo_ijk * ntno_ijk);
#pragma omp critical
        q_ov->load(psio_.get(), PSIF_DLPNO_QIA_TNO, Matrix::SubBlocks);

        q_ov_[ijk] = Tensor<double, 3>(q_ov->name(), naux_ijk, nlmo_ijk, ntno_ijk);
        ::memcpy(q_ov_[ijk].data(), q_ov->get_pointer(), naux_ijk * nlmo_ijk * ntno_ijk * sizeof(double));
    }
}

void DLPNOCCSDT::load_qab_tno(int ijk) {
    if (disk_ints_) {
        int naux_ijk = lmotriplet_to_ribfs_[ijk].size();
        int ntno_ijk = n_tno_[ijk];

        std::stringstream q_vv_name;
        q_vv_name << "(Q_ijk | a b) " << (ijk);

        auto q_vv = std::make_shared<Matrix>(q_vv_name.str(), naux_ijk, ntno_ijk * ntno_ijk);
#pragma omp critical
        q_vv->load(psio_.get(), PSIF_DLPNO_QAB_TNO, Matrix::ThreeIndexLowerTriangle);

        q_vv_[ijk] = Tensor<double, 3>(q_vv->name(), naux_ijk, ntno_ijk, ntno_ijk);
        ::memcpy(q_vv_[ijk].data(), q_vv->get_pointer(), naux_ijk * ntno_ijk * ntno_ijk * sizeof(double));
    }
}

void DLPNOCCSDT::compute_R_ia_triples(std::vector<SharedMatrix>& R_ia, std::vector<std::vector<SharedMatrix>>& R_ia_buffer) {

    size_t naocc = i_j_to_ij_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();

    // Start from the CCSD part, then add the canonical T3 term of manuscript Eq. (32)
    // in its domain-projected form, Eq. (80) and Algorithm 1.
    DLPNOCCSD::compute_R_ia(R_ia, R_ia_buffer);

    // Clean buffers
#pragma omp parallel for
    for (int thread = 0; thread < nthread_; ++thread) {
        for (int i = 0; i < naocc; ++i) {
            R_ia_buffer[thread][i]->zero();
        } // end thread
    } // end i

    // Loop over unique triplets as prescribed by DLPNO-CCSDT Algorithm 1.
#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        auto &[i, j, k] = ijk_to_i_j_k_[ijk];

        int ntno_ijk = n_tno_[ijk];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        double prefactor = 1.0;
        if (i == j || j == k) {
            prefactor = 0.5;
        }

        std::vector<std::tuple<int, int, int>> P_S = {std::make_tuple(i, j, k), std::make_tuple(j, i, k), std::make_tuple(k, i, j)};
        const auto& K_cache = K_ovov_tno_cache_[ijk];
        std::array<const Tensor<double, 2>*, 3> K_ovov_list = {
            &K_cache[2], &K_cache[1], &K_cache[0]};

        // Each diagonal-pair overlap is independent of the T3 column
        // permutation. Cache the three transforms locally instead of rebuilding
        // them in every contribution.
        std::array<SharedMatrix, 3> S_ijk_diag;
        const std::array<int, 3> output_lmos = {i, j, k};
        for (int output_idx = 0; output_idx < 3; ++output_idx) {
            const int ii = i_j_to_ij_[output_lmos[output_idx]][output_lmos[output_idx]];
            auto S = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk],
                                             lmopair_to_paos_[ii]);
            S_ijk_diag[output_idx] =
                linalg::triplet(X_tno_[ijk], S, X_pno_[ii], true, false, false);
        }

        for (int perm_idx = 0; perm_idx < P_S.size(); ++perm_idx) {
            auto &[i, j, k] = P_S[perm_idx];
            const Tensor<double, 3>* U_ijk = &U_iajbkc_[ijk];
            Tensor<double, 3> U_ijk_permuted;
            if (perm_idx == 1) {
                U_ijk_permuted = Tensor<double, 3>(
                    "U_ijk R1 permutation", ntno_ijk, ntno_ijk, ntno_ijk);
                permute(Indices{index::b, index::a, index::c},
                        &U_ijk_permuted,
                        Indices{index::a, index::b, index::c},
                        U_iajbkc_[ijk]);
                U_ijk = &U_ijk_permuted;
            } else if (perm_idx == 2) {
                U_ijk_permuted = Tensor<double, 3>(
                    "U_ijk R1 permutation", ntno_ijk, ntno_ijk, ntno_ijk);
                permute(Indices{index::c, index::a, index::b},
                        &U_ijk_permuted,
                        Indices{index::a, index::b, index::c},
                        U_iajbkc_[ijk]);
                U_ijk = &U_ijk_permuted;
            }

            Tensor<double, 1> R_ia_cont("R_ia_cont", n_tno_[ijk]);
            einsum(0.0, Indices{index::a}, &R_ia_cont, prefactor,
                   Indices{index::a, index::b, index::c}, *U_ijk,
                   Indices{index::b, index::c}, *K_ovov_list[perm_idx]);

            auto R_ia_psi = std::make_shared<Matrix>(n_tno_[ijk], 1);
            ::memcpy(R_ia_psi->get_pointer(), &(R_ia_cont)(0), n_tno_[ijk] * sizeof(double));

            R_ia_buffer[thread][i]->add(
                linalg::doublet(S_ijk_diag[perm_idx], R_ia_psi, true, false));
        }

    } // end ijk

    // Flush buffers
#pragma omp parallel for
    for (int i = 0; i < naocc; ++i) {
        for (int thread = 0; thread < nthread_; ++thread) {
            R_ia[i]->add(R_ia_buffer[thread][i]);
        } // end thread
    } // end i
}

void DLPNOCCSDT::compute_R_iajb_triples(std::vector<SharedMatrix>& R_iajb, std::vector<SharedMatrix>& Rn_iajb,
                                        std::vector<std::vector<SharedMatrix>>& R_iajb_buffer) {

    size_t n_lmo_pairs = ij_to_i_j_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();

    // Start from the CCSD part, then add the canonical T3 terms of manuscript Eq. (33)
    // and apply Eq. (34) through the local contractions in Eqs. (87)-(90), Algorithm 2.
    DLPNOCCSD::compute_R_iajb(R_iajb, Rn_iajb);

    // Clean buffers
#pragma omp parallel for
    for (int thread = 0; thread < nthread_; ++thread) {
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            R_iajb_buffer[thread][ij]->zero();
        } // end thread
    } // end i

    // Loop over unique triplets as prescribed by DLPNO-CCSDT Algorithm 2.
#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        auto &[i, j, k] = ijk_to_i_j_k_[ijk];
        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];

        int nlmo_ijk = lmotriplet_to_lmos_[ijk].size();
        int naux_ijk = lmotriplet_to_ribfs_[ijk].size();
        int ntno_ijk = n_tno_[ijk];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // Buffers
        Tensor<double, 2> R_iajb_cont_a("R_iajb_cont_a", n_tno_[ijk], n_tno_[ijk]);
        Tensor<double, 3> R_iajb_cont_b("R_iajb_cont_b", nlmo_ijk, n_tno_[ijk], n_tno_[ijk]);
        SharedMatrix psi_buffer = std::make_shared<Matrix>(n_tno_[ijk], n_tno_[ijk]);

        double prefactor = 1.0;
        if (i == j || j == k) {
            prefactor = 0.5;
        }

        // Read integrals when disk-backed storage is enabled.
        if (disk_ints_) {
            load_qia_tno(ijk);
            load_qab_tno(ijk);
        }

        // => (i l_{ijk} | j b_{ijk}) integrals <= //

        Tensor<double, 2> K_iojv = Tensor<double, 2>("K_iojv", nlmo_ijk, ntno_ijk);
        Tensor<double, 2> K_joiv = Tensor<double, 2>("K_joiv", nlmo_ijk, ntno_ijk);
        Tensor<double, 2> K_jokv = Tensor<double, 2>("K_jokv", nlmo_ijk, ntno_ijk);
        Tensor<double, 2> K_kojv = Tensor<double, 2>("K_kojv", nlmo_ijk, ntno_ijk);
        Tensor<double, 2> K_iokv = Tensor<double, 2>("K_iokv", nlmo_ijk, ntno_ijk);
        Tensor<double, 2> K_koiv = Tensor<double, 2>("K_koiv", nlmo_ijk, ntno_ijk);

        einsum(0.0, Indices{index::l, index::a}, &K_iojv, 1.0, Indices{index::Q, index::l}, q_io_[ijk], Indices{index::Q, index::a}, q_jv_[ijk]);
        einsum(0.0, Indices{index::l, index::a}, &K_joiv, 1.0, Indices{index::Q, index::l}, q_jo_[ijk], Indices{index::Q, index::a}, q_iv_[ijk]);
        einsum(0.0, Indices{index::l, index::a}, &K_jokv, 1.0, Indices{index::Q, index::l}, q_jo_[ijk], Indices{index::Q, index::a}, q_kv_[ijk]);
        einsum(0.0, Indices{index::l, index::a}, &K_kojv, 1.0, Indices{index::Q, index::l}, q_ko_[ijk], Indices{index::Q, index::a}, q_jv_[ijk]);
        einsum(0.0, Indices{index::l, index::a}, &K_iokv, 1.0, Indices{index::Q, index::l}, q_io_[ijk], Indices{index::Q, index::a}, q_kv_[ijk]);
        einsum(0.0, Indices{index::l, index::a}, &K_koiv, 1.0, Indices{index::Q, index::l}, q_ko_[ijk], Indices{index::Q, index::a}, q_iv_[ijk]);

        // => (i a_{ijk} | j b_{ijk}) integrals <= //
        // Reuse the iteration-invariant cache shared with the T3 -> R1 pass.
        // Reversed orientations are expressed by index labels below, avoiding
        // three full ntno x ntno transpose copies.
        const auto& K_cache = K_ovov_tno_cache_[ijk];

        // => (i a_{ijk} | l_{ijk} d_{ijk}) integrals <= //

        Tensor<double, 3> K_ivov = Tensor<double, 3>("K_ivov", nlmo_ijk, ntno_ijk, ntno_ijk);
        Tensor<double, 3> K_jvov = Tensor<double, 3>("K_jvov", nlmo_ijk, ntno_ijk, ntno_ijk);
        Tensor<double, 3> K_kvov = Tensor<double, 3>("K_kvov", nlmo_ijk, ntno_ijk, ntno_ijk);

        einsum(0.0, Indices{index::m, index::a, index::b}, &K_ivov, 1.0, Indices{index::Q, index::m, index::a}, q_ov_[ijk], 
                    Indices{index::Q, index::b}, q_iv_[ijk]);
        einsum(0.0, Indices{index::m, index::a, index::b}, &K_jvov, 1.0, Indices{index::Q, index::m, index::a}, q_ov_[ijk],
                    Indices{index::Q, index::b}, q_jv_[ijk]);
        einsum(0.0, Indices{index::m, index::a, index::b}, &K_kvov, 1.0, Indices{index::Q, index::m, index::a}, q_ov_[ijk],
                    Indices{index::Q, index::b}, q_kv_[ijk]);

        // The (ia|bc), (ja|bc), and (ka|bc) contributions are consumed in
        // streamed DF form below.  Avoid three resident ntno^3 integral cubes.

        // Form the extended domain used for TNO overlaps.
        std::vector<int> triplet_ext_domain = merge_lists(merge_lists(lmopair_to_paos_[ij], lmopair_to_paos_[jk]), lmopair_to_paos_[ik]);

        for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
            int l = lmotriplet_to_lmos_[ijk][l_ijk];
            int il = i_j_to_ij_[i][l], jl = i_j_to_ij_[j][l], kl = i_j_to_ij_[k][l];
            triplet_ext_domain = merge_lists(triplet_ext_domain, merge_lists(merge_lists(lmopair_to_paos_[il], 
                                    lmopair_to_paos_[jl]), lmopair_to_paos_[kl]));

        } // end l

        // Semi-direct TNO overlap algorithm
        auto S_ijk = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk], triplet_ext_domain);
        S_ijk = linalg::doublet(X_tno_[ijk], S_ijk, true, false);

        std::unordered_map<int, SharedMatrix> tno_to_pno_overlap;
        auto overlap_for_pair = [&](int pair) -> SharedMatrix {
            auto found = tno_to_pno_overlap.find(pair);
            if (found != tno_to_pno_overlap.end()) return found->second;
            auto S_pair = submatrix_cols(
                *S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[pair]));
            S_pair = linalg::doublet(S_pair, X_pno_[pair], false, false);
            tno_to_pno_overlap.emplace(pair, S_pair);
            return S_pair;
        };

        // => F_kc contribution: DLPNO Eq. (87), from canonical Eq. (33) <= //

        std::vector<std::tuple<int, int, int>> P_S = {std::make_tuple(i, j, k), std::make_tuple(i, k, j), std::make_tuple(j, k, i)};
        std::array<const Tensor<double, 3>*, 3> K_kvov_fock = {
            &K_kvov, &K_jvov, &K_ivov};

        for (int idx = 0; idx < P_S.size(); ++idx) {
            auto &[i, j, k] = P_S[idx];
            int ij = i_j_to_ij_[i][j];

            const Tensor<double, 3>* U_ijk = &U_iajbkc_[ijk];
            Tensor<double, 3> U_ijk_permuted;
            if (idx != 0) {
                U_ijk_permuted = Tensor<double, 3>(
                    "U_ijk R2 Fock permutation", ntno_ijk, ntno_ijk,
                    ntno_ijk);
                // The destination labels encode the physical occupied-column
                // order stored in this fresh contiguous tensor. The following
                // contraction then assigns its physical axes the canonical
                // (a,b,c) labels. This direction matters for the three-cycle:
                // P(j,k,i) and its inverse are distinct.
                if (idx == 1) {
                    permute(Indices{index::a, index::c, index::b},
                            &U_ijk_permuted,
                            Indices{index::a, index::b, index::c},
                            U_iajbkc_[ijk]);
                } else {
                    permute(Indices{index::b, index::c, index::a},
                            &U_ijk_permuted,
                            Indices{index::a, index::b, index::c},
                            U_iajbkc_[ijk]);
                }
                U_ijk = &U_ijk_permuted;
            }

            // (T1-dressed Fock Matrix) F_kc = [2.0 * (kc|ld) - (kd|lc)] t_{ld}
            Tensor<double, 1> Fkc("Fkc", ntno_ijk);
            Tensor<double, 3> K_lckd("K_lckd", nlmo_ijk, ntno_ijk,
                                      ntno_ijk);
            permute(Indices{index::l, index::c, index::d}, &K_lckd,
                    Indices{index::l, index::d, index::c},
                    *K_kvov_fock[idx]);
            einsum(0.0, Indices{index::c}, &Fkc, 2.0, Indices{index::l, index::d, index::c}, *K_kvov_fock[idx], Indices{index::l, index::d}, T_n_ijk_[ijk]);
            einsum(1.0, Indices{index::c}, &Fkc, -1.0,
                   Indices{index::l, index::d, index::c}, K_lckd,
                   Indices{index::l, index::d}, T_n_ijk_[ijk]);

            einsum(0.0, Indices{index::a, index::b}, &R_iajb_cont_a,
                   prefactor, Indices{index::a, index::b, index::c}, *U_ijk,
                   Indices{index::c}, Fkc);
            ::memcpy(psi_buffer->get_pointer(), R_iajb_cont_a.data(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));

            auto S_ijk_ij = overlap_for_pair(ij);
            R_iajb_buffer[thread][ij]->add(linalg::triplet(S_ijk_ij, psi_buffer, S_ijk_ij, true, false, false));
        }

        // => (db|kc) and (jl|kc) contributions: DLPNO Eqs. (88)-(89), canonical Eq. (33) <= //
        std::vector<std::tuple<int, int, int>> perms = {std::make_tuple(i, j, k), std::make_tuple(i, k, j),
                                                        std::make_tuple(j, i, k), std::make_tuple(j, k, i),
                                                        std::make_tuple(k, i, j), std::make_tuple(k, j, i)};

        std::array<const Tensor<double, 2>*, 6> K_jokv_list = {
            &K_jokv, &K_kojv, &K_iokv, &K_koiv, &K_iojv, &K_joiv};
        constexpr std::array<int, 6> K_cache_component = {2, 2, 1, 1, 0, 0};
        constexpr std::array<bool, 6> K_cache_transposed = {false, true, false, true, false, true};
        std::array<const Tensor<double, 2>*, 6> q_kv_list = {
            &q_kv_[ijk], &q_jv_[ijk], &q_kv_[ijk],
            &q_iv_[ijk], &q_jv_[ijk], &q_iv_[ijk]};
        std::array<const Tensor<double, 3>*, 6> K_kvov_list = {
            &K_kvov, &K_jvov, &K_kvov, &K_ivov, &K_jvov, &K_ivov};

        for (int idx = 0; idx < perms.size(); ++idx) {
            auto &[i, j, k] = perms[idx];
            int ij = i_j_to_ij_[i][j];

            const Tensor<double, 3>* U_ijk = &U_iajbkc_[ijk];
            Tensor<double, 3> U_ijk_permuted;
            if (idx != 0) {
                U_ijk_permuted = Tensor<double, 3>(
                    "U_ijk R2 permutation", ntno_ijk, ntno_ijk, ntno_ijk);
                if (idx == 1) {
                    permute(Indices{index::a, index::c, index::b},
                            &U_ijk_permuted,
                            Indices{index::a, index::b, index::c},
                            U_iajbkc_[ijk]);
                } else if (idx == 2) {
                    permute(Indices{index::b, index::a, index::c},
                            &U_ijk_permuted,
                            Indices{index::a, index::b, index::c},
                            U_iajbkc_[ijk]);
                } else if (idx == 3) {
                    permute(Indices{index::b, index::c, index::a},
                            &U_ijk_permuted,
                            Indices{index::a, index::b, index::c},
                            U_iajbkc_[ijk]);
                } else if (idx == 4) {
                    permute(Indices{index::c, index::a, index::b},
                            &U_ijk_permuted,
                            Indices{index::a, index::b, index::c},
                            U_iajbkc_[ijk]);
                } else {
                    permute(Indices{index::c, index::b, index::a},
                            &U_ijk_permuted,
                            Indices{index::a, index::b, index::c},
                            U_iajbkc_[ijk]);
                }
                U_ijk = &U_ijk_permuted;
            }

            // (T1-dressed integral g_jlkc)
            // (jl|kc)_t1 = (jl|kc) + (jd|kc)t_{l}^{d}
            Tensor<double, 2> g_jlkc = *K_jokv_list[idx];
            const auto& K_ovov = K_cache[K_cache_component[idx]];
            if (K_cache_transposed[idx]) {
                einsum(1.0, Indices{index::l, index::c}, &g_jlkc, 1.0,
                       Indices{index::c, index::d}, K_ovov,
                       Indices{index::l, index::d}, T_n_ijk_[ijk]);
            } else {
                einsum(1.0, Indices{index::l, index::c}, &g_jlkc, 1.0,
                       Indices{index::d, index::c}, K_ovov,
                       Indices{index::l, index::d}, T_n_ijk_[ijk]);
            }

            auto contract_u3_terms = [&]() {
                // The (jl|kc) contraction is a GEMM after explicitly packing
                // U3 in the requested occupied-column order.
                einsum(0.0, Indices{index::l, index::a, index::b},
                       &R_iajb_cont_b, prefactor,
                       Indices{index::a, index::b, index::c}, *U_ijk,
                       Indices{index::l, index::c}, g_jlkc);

                // Direct DF form of U_3^{abc}(db|kc): for each auxiliary
                // function, contract c first (GEMV), then b (GEMM).  This
                // replaces both a resident v^3 integral and the equally large
                // T1-dressed g_dbkc work tensor with one v^2 matrix.
                R_iajb_cont_a.zero();
                for (int q = 0; q < naux_ijk; ++q) {
                    Tensor<double, 1> q_kc("q_kc", ntno_ijk);
                    q_kc = (*q_kv_list[idx])(q, All);
                    Tensor<double, 2> U_ab("U_ab", ntno_ijk, ntno_ijk);
                    einsum(0.0, Indices{index::a, index::b}, &U_ab, 1.0,
                           Indices{index::a, index::b, index::c}, *U_ijk,
                           Indices{index::c}, q_kc);
                    auto q_db = q_vv_[ijk](q, All, All);
                    einsum(1.0, Indices{index::a, index::d},
                           &R_iajb_cont_a, prefactor,
                           Indices{index::a, index::b}, U_ab,
                           Indices{index::d, index::b}, q_db);
                }

                // T1 dressing of the first virtual integral index, as in
                // manuscript Eq. (26) and DLPNO Algorithm 2:
                // (db|kc)_t1 = (db|kc) - (lb|kc)t_l^d. This is intentionally
                // different from the second-index q_vv_t1 convention used
                // later in the R3 implementation. Contract (b,c) first and
                // finish with one GEMM over l.
                Tensor<double, 2> U_al("U_al", ntno_ijk, nlmo_ijk);
                einsum(0.0, Indices{index::a, index::l}, &U_al, 1.0,
                       Indices{index::a, index::b, index::c}, *U_ijk,
                       Indices{index::l, index::b, index::c},
                       *K_kvov_list[idx]);
                einsum(1.0, Indices{index::a, index::d},
                       &R_iajb_cont_a, -prefactor,
                       Indices{index::a, index::l}, U_al,
                       Indices{index::l, index::d}, T_n_ijk_[ijk]);
            };
            contract_u3_terms();

            // Copy into a Psi4 Matrix for the subsequent Matrix-based basis transformation.
            ::memcpy(psi_buffer->get_pointer(), R_iajb_cont_a.data(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));

            auto S_ijk_ij = overlap_for_pair(ij);
            R_iajb_buffer[thread][ij]->add(linalg::triplet(S_ijk_ij, psi_buffer, S_ijk_ij, true, false, false));

            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                int il = i_j_to_ij_[i][l];
                // Flush il buffer
                ::memcpy(psi_buffer->get_pointer(), &(R_iajb_cont_b)(l_ijk, 0, 0), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));

                auto S_ijk_il = overlap_for_pair(il);
                R_iajb_buffer[thread][il]->subtract(linalg::triplet(S_ijk_il, psi_buffer, S_ijk_il, true, false, false));
            } // end l_ijk
        } // end perms

        if (disk_ints_) {
            q_ov_[ijk] = Tensor<double, 3>(q_ov_[ijk].name(), 0, 0, 0);
            q_vv_[ijk] = Tensor<double, 3>(q_vv_[ijk].name(), 0, 0, 0);
        }

    } // end ijk

    // Apply the spin-adapted pair permutation of DLPNO Eq. (90), corresponding to
    // the canonical permutation operator in Eq. (34), while flushing the buffers.
#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int ji = ij_to_ji_[ij];
        for (int thread = 0; thread < nthread_; ++thread) {
            auto cont_a = R_iajb_buffer[thread][ij]->clone();
            cont_a->scale(2.0 / 3.0);
            auto cont_b = R_iajb_buffer[thread][ji]->transpose();
            cont_b->scale(2.0 / 3.0);
            auto cont_c = R_iajb_buffer[thread][ij]->transpose();
            cont_c->scale(1.0 / 3.0);
            auto cont_d = R_iajb_buffer[thread][ji]->clone();
            cont_d->scale(1.0 / 3.0);

            // Add all contributions
            R_iajb[ij]->add(cont_a);
            R_iajb[ij]->add(cont_b);
            R_iajb[ij]->add(cont_c);
            R_iajb[ij]->add(cont_d);
        } // end thread
    } // end ij
}

std::vector<Tensor<double, 3>> DLPNOCCSDT::rho_dbck_contribution() {

    // This helper extracts the most expensive cross-domain term implied by rho_k(d,b,c),
    // manuscript Eq. (101) (canonical counterpart: Eqs. (41)-(42)). In SI Algorithm S3,
    // it is the triples-amplitude contraction in the final update of rho_k. Forming that
    // contribution once in the diagonal PNO domain of k and projecting it into each TNO
    // domain avoids repeating the same contraction for every target triplet. This is an
    // algebraic factorization of the published term, not a separate approximation.

    int naocc = i_j_to_ij_.size();

    std::vector<Tensor<double, 3>> rho_dbck_cont(naocc);

#pragma omp parallel for
    for (int k = 0; k < naocc; ++k) {
        int kk = i_j_to_ij_[k][k];
        int nlmo_k = lmopair_to_lmos_[kk].size();
        int naux_k = lmopair_to_ribfs_[kk].size();
        int npno_k = n_pno_[kk];

        rho_dbck_cont[k] = Tensor<double, 3>("rho_dbck_cont", npno_k, npno_k, npno_k);
        rho_dbck_cont[k].zero();

        std::vector<SharedMatrix> q_ov_k = QIA_PNO(kk);

        // Pack B_Q^{me} once with (Q,e) contiguous for every m.  The old
        // implementation materialized every (me|nf) block at once, requiring
        // nlmo^2*npno^2 doubles.  Each (m,l) slice below is instead one GEMM
        // B_m^T B_l and is consumed immediately.
        Tensor<double, 3> q_mqe("q_mqe", nlmo_k, naux_k, npno_k);
        for (int q_k = 0; q_k < naux_k; ++q_k) {
            for (int m_k = 0; m_k < nlmo_k; ++m_k) {
                ::memcpy(&q_mqe(m_k, q_k, 0), q_ov_k[q_k]->pointer()[m_k],
                         static_cast<size_t>(npno_k) * sizeof(double));
            }
        }
        q_ov_k.clear();

        for (int m_k = 0; m_k < nlmo_k; ++m_k) {
            int m = lmopair_to_lmos_[kk][m_k];
            for (int l_k = 0; l_k < nlmo_k; ++l_k) {
                int l = lmopair_to_lmos_[kk][l_k];

                const size_t mlk_dense = triplet_key(m, l, k, naocc);
                if (!i_j_k_to_ijk_.count(mlk_dense)) continue;
                int mlk = i_j_k_to_ijk_[mlk_dense];

                // (2 t_{mlk}^{ebc} - t_{mlk}^{cbe} - t_{mlk}^{bec})
                auto S_mlk_k = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[mlk], lmopair_to_paos_[kk]);
                S_mlk_k = linalg::triplet(X_tno_[mlk], S_mlk_k, X_pno_[kk], true, false, false);

                Tensor<double, 3> T_mlk = matmul_3d_einsums(triples_permuter_einsums(
                                            T_iajbkc_clone_[mlk], m, l, k), S_mlk_k->transpose(), n_tno_[mlk], n_pno_[kk]);
                Tensor<double, 3> Z_mlk = T_mlk;
                Z_mlk *= 2.0;
                Z_mlk -= triples_permuter_einsums(T_mlk, 1, 0, 2);
                Z_mlk -= triples_permuter_einsums(T_mlk, 2, 1, 0);

                Tensor<double, 2> g_ed("(me|ld)", npno_k, npno_k);
                TensorView<double, 2> q_m = q_mqe(m_k, All, All);
                TensorView<double, 2> q_l = q_mqe(l_k, All, All);
                einsum(0.0, Indices{index::e, index::d}, &g_ed, 1.0,
                       Indices{index::Q, index::e}, q_m,
                       Indices{index::Q, index::d}, q_l);
                einsum(1.0, Indices{index::d, index::b, index::c},
                       &rho_dbck_cont[k], -1.0,
                       Indices{index::e, index::d}, g_ed,
                       Indices{index::e, index::b, index::c}, Z_mlk);
            } // end m_ijk
        } // end l_ijk
    } // end k

    return rho_dbck_cont;
}

void DLPNOCCSDT::compute_R_iajbkc(std::vector<SharedMatrix>& R_iajbkc) {

    size_t naocc = i_j_to_ij_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();

    // Compute the reusable part of rho_k(d,b,c), manuscript Eq. (101) / SI Algorithm S3.
    std::vector<Tensor<double, 3>> rho_dbck_cont = rho_dbck_contribution();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        auto &[i, j, k] = ijk_to_i_j_k_[ijk];
        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];

        int ntno_ijk = n_tno_[ijk];
        int naux_ijk = lmotriplet_to_ribfs_[ijk].size();
        int nlmo_ijk = lmotriplet_to_lmos_[ijk].size();

        auto R_ijk = R_iajbkc[ijk];
        auto T_ijk = T_iajbkc_[ijk];

        R_ijk->zero();

        // Read integrals when disk-backed storage is enabled.
        if (disk_ints_) {
            load_qia_tno(ijk);
            load_qab_tno(ijk);
        }

        // Form extended domain for overlap integrals
        std::vector<int> triplet_ext_domain = merge_lists(merge_lists(lmopair_to_paos_[ij], lmopair_to_paos_[jk]), lmopair_to_paos_[ik]);

        for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
            int m = lmotriplet_to_lmos_[ijk][m_ijk];
            int mi = i_j_to_ij_[m][i], mj = i_j_to_ij_[m][j], mk = i_j_to_ij_[m][k];

            triplet_ext_domain = merge_lists(triplet_ext_domain, lmopair_to_paos_[mi]);
            triplet_ext_domain = merge_lists(triplet_ext_domain, lmopair_to_paos_[mj]);
            triplet_ext_domain = merge_lists(triplet_ext_domain, lmopair_to_paos_[mk]);

            const size_t mjk_dense = triplet_key(m, j, k, naocc);
            if (i_j_k_to_ijk_.count(mjk_dense)) {
                int mjk = i_j_k_to_ijk_[mjk_dense];
                triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[mjk]);
            } // end if

            const size_t imk_dense = triplet_key(i, m, k, naocc);
            if (i_j_k_to_ijk_.count(imk_dense)) {
                int imk = i_j_k_to_ijk_[imk_dense];
                triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[imk]);
            } // end if

            const size_t ijm_dense = triplet_key(i, j, m, naocc);
            if (i_j_k_to_ijk_.count(ijm_dense)) {
                int ijm = i_j_k_to_ijk_[ijm_dense];
                triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[ijm]);
            } // end if

            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                int li = i_j_to_ij_[l][i], lj = i_j_to_ij_[l][j], lk = i_j_to_ij_[l][k], ml = i_j_to_ij_[m][l];
                if (ml == -1) continue;

                triplet_ext_domain = merge_lists(triplet_ext_domain, lmopair_to_paos_[li]);
                triplet_ext_domain = merge_lists(triplet_ext_domain, lmopair_to_paos_[lj]);
                triplet_ext_domain = merge_lists(triplet_ext_domain, lmopair_to_paos_[lk]);
                triplet_ext_domain = merge_lists(triplet_ext_domain, lmopair_to_paos_[ml]);

                const size_t mli_dense = triplet_key(m, l, i, naocc);
                if (i_j_k_to_ijk_.count(mli_dense)) {
                    int mli = i_j_k_to_ijk_[mli_dense];
                    triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[mli]);
                } // end if

                const size_t mlj_dense = triplet_key(m, l, j, naocc);
                if (i_j_k_to_ijk_.count(mlj_dense)) {
                    int mlj = i_j_k_to_ijk_[mlj_dense];
                    triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[mlj]);
                } // end if

                const size_t mlk_dense = triplet_key(m, l, k, naocc);
                if (i_j_k_to_ijk_.count(mlk_dense)) {
                    int mlk = i_j_k_to_ijk_[mlk_dense];
                    triplet_ext_domain = merge_lists(triplet_ext_domain, lmotriplet_to_paos_[mlk]);
                }
            } // end l_ijk
        } // end m_ijk

        // Semi-direct computation of TNOs
        auto S_ijk = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk], triplet_ext_domain);
        S_ijk = linalg::doublet(X_tno_[ijk], S_ijk, true, false);

        std::vector<int> ijk_idx = {i, j, k};

        std::vector<Tensor<double, 2>> q_io_list = {q_io_[ijk], q_jo_[ijk], q_ko_[ijk]};
        std::vector<Tensor<double, 2>> q_iv_list = {q_iv_[ijk], q_jv_[ijk], q_kv_[ijk]};

        std::vector<Tensor<double, 1>> T_i_list(ijk_idx.size());
        std::vector<Tensor<double, 2>> q_io_t1_list(ijk_idx.size());
        std::vector<Tensor<double, 2>> q_iv_t1_list(ijk_idx.size());

        // => DF factorization and T1-dressed three-index integrals: manuscript Eqs. (20)-(26),
        //    SI Eq. (S1) and Algorithm S1 <= //

        for (int i_idx = 0; i_idx < ijk_idx.size(); ++i_idx) {
            int i = ijk_idx[i_idx];
            auto i_ijk = std::find(lmotriplet_to_lmos_[ijk].begin(), lmotriplet_to_lmos_[ijk].end(), i) - lmotriplet_to_lmos_[ijk].begin();

            // T_{i}^{a_{ijk}}
            Tensor<double, 1> T_i("T_i", ntno_ijk);
            ::memcpy(T_i.data(), &(T_n_ijk_[ijk])(i_ijk, 0), ntno_ijk * sizeof(double));
            T_i_list[i_idx] = T_i;

            // T1-dressing of q_iv
            q_iv_t1_list[i_idx] = q_iv_list[i_idx];
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], -1.0, Indices{index::Q, index::l}, q_io_list[i_idx], Indices{index::l, index::a}, T_n_ijk_[ijk]);
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], 1.0, Indices{index::Q, index::a, index::b}, q_vv_[ijk], Indices{index::b}, T_i);

            Tensor<double, 2> q_iv_t1_temp("q_iv_t1_temp", naux_ijk, nlmo_ijk);
            einsum(0.0, Indices{index::Q, index::l}, &q_iv_t1_temp, 1.0, Indices{index::Q, index::l, index::b}, q_ov_[ijk], Indices{index::b}, T_i);
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], -1.0, Indices{index::Q, index::l}, q_iv_t1_temp, Indices{index::l, index::a}, T_n_ijk_[ijk]);

            // T1-dressing of q_io
            q_io_t1_list[i_idx] = q_io_list[i_idx];
            einsum(1.0, Indices{index::Q, index::l}, &q_io_t1_list[i_idx], 1.0, Indices{index::Q, index::l, index::a}, q_ov_[ijk], Indices{index::a}, T_i);
        } // end i_idx

        // This one is special... the second index is dressed instead of the first (per convention), to increase computational efficiency
        Tensor<double, 3> q_vv_t1 = q_vv_[ijk];
        Tensor<double, 3> q_vo("q_vo", naux_ijk, ntno_ijk, nlmo_ijk);
        permute(Indices{index::Q, index::a, index::l}, &q_vo, Indices{index::Q, index::l, index::a}, q_ov_[ijk]);
        einsum(1.0, Indices{index::Q, index::a, index::b}, &q_vv_t1, -1.0, Indices{index::Q, index::a, index::l}, q_vo, Indices{index::l, index::b}, T_n_ijk_[ijk]);

        // => Form dressed Fock intermediates: manuscript Eqs. (27)-(31), SI Eqs. (S2)-(S4)
        //    and Algorithm S2 <= //

        // Gamma_Q is used universally for J-like contractions
        Tensor<double, 1> gamma_Q("gamma_Q", naux_ijk);
        einsum(0.0, Indices{index::Q}, &gamma_Q, 1.0, Indices{index::Q, index::m, index::e}, q_ov_[ijk], Indices{index::m, index::e}, T_n_ijk_[ijk]);

        // => F_ld (this is scoped to ensure that the intermediate tensors are not persistent in memory <= //
        Tensor<double, 2> F_ld("F_ld", nlmo_ijk, ntno_ijk); {
            // J contractions
            einsum(0.0, Indices{index::l, index::d}, &F_ld, 2.0, Indices{index::Q, index::l, index::d}, q_ov_[ijk], Indices{index::Q}, gamma_Q);

            // K contractions
            Tensor<double, 3> F_ld_K_temp("F_ld_K_temp", naux_ijk, nlmo_ijk, nlmo_ijk);
            einsum(0.0, Indices{index::Q, index::l, index::m}, &F_ld_K_temp, 1.0, Indices{index::Q, index::l, index::e}, q_ov_[ijk], Indices{index::m, index::e}, T_n_ijk_[ijk]);
            Tensor<double, 3> F_ld_K_temp2("F_ld_K_temp2", naux_ijk, nlmo_ijk, nlmo_ijk);
            permute(Indices{index::Q, index::m, index::l}, &F_ld_K_temp2, Indices{index::Q, index::l, index::m}, F_ld_K_temp);
            einsum(1.0, Indices{index::l, index::d}, &F_ld, -1.0, Indices{index::Q, index::m, index::l}, F_ld_K_temp2, Indices{index::Q, index::m, index::d}, q_ov_[ijk]);
        }

        // => F_li <= //
        std::vector<Tensor<double, 1>> F_li_list(ijk_idx.size()); // F_li, F_lj, F_lk

        for (int idx = 0; idx < ijk_idx.size(); ++idx) {
            int i = ijk_idx[idx];

            // F_lmo (non-T1 contribution)
            F_li_list[idx] = Tensor<double, 1>("F_li", nlmo_ijk);
            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                F_li_list[idx](l_ijk) = (*F_lmo_)(l, i);
            }

            // J Contractions
            einsum(1.0, Indices{index::l}, &F_li_list[idx], 2.0, Indices{index::Q, index::l}, q_io_list[idx], Indices{index::Q}, gamma_Q);

            // K contractions
            Tensor<double, 2> F_li_K_temp("F_li_K_temp", naux_ijk, ntno_ijk);
            einsum(0.0, Indices{index::Q, index::a}, &F_li_K_temp, 1.0, Indices{index::Q, index::m}, q_io_list[idx], Indices{index::m, index::a}, T_n_ijk_[ijk]);
            einsum(1.0, Indices{index::l}, &F_li_list[idx], -1.0, Indices{index::Q, index::a, index::l}, q_vo, Indices{index::Q, index::a}, F_li_K_temp);

            // Add F_ld contribution
            einsum(1.0, Indices{index::l}, &F_li_list[idx], 1.0, Indices{index::l, index::d}, F_ld, Indices{index::d}, T_i_list[idx]);
        }

        // => F_ad <= //
        Tensor<double, 2> F_ad("F_ad", ntno_ijk, ntno_ijk); {
            F_ad.zero();
            // e_tno (non-t1 contribution)
            for (int a_ijk = 0; a_ijk < ntno_ijk; ++a_ijk) {
                F_ad(a_ijk, a_ijk) = (*e_tno_[ijk])(a_ijk);
            }

            // J contribution
            einsum(1.0, Indices{index::a, index::d}, &F_ad, 2.0, Indices{index::Q, index::a, index::d}, q_vv_[ijk], Indices{index::Q}, gamma_Q);
            // K contribution
            Tensor<double, 3> F_ad_K_temp("F_ad_K_temp", naux_ijk, ntno_ijk, nlmo_ijk);
            einsum(0.0, Indices{index::Q, index::a, index::m}, &F_ad_K_temp, 1.0, Indices{index::Q, index::a, index::e}, q_vv_[ijk], Indices{index::m, index::e}, T_n_ijk_[ijk]);
            Tensor<double, 3> F_ad_K_temp2("F_ad_K_temp2", naux_ijk, nlmo_ijk, ntno_ijk);
            permute(Indices{index::Q, index::m, index::a}, &F_ad_K_temp2, Indices{index::Q, index::a, index::m}, F_ad_K_temp);
            einsum(1.0, Indices{index::a, index::d}, &F_ad, -1.0, Indices{index::Q, index::m, index::a}, F_ad_K_temp2, Indices{index::Q, index::m, index::d}, q_ov_[ijk]);

            // Add the F_ld contribution to F_ad
            einsum(1.0, Indices{index::a, index::d}, &F_ad, -1.0, Indices{index::m, index::a}, T_n_ijk_[ijk], Indices{index::m, index::d}, F_ld);
        }

        // rho_i/rho_j/rho_k(d,b,c): manuscript Eq. (101), canonical Eqs. (41)-(42),
        // SI Eq. (S5) and Algorithm S3.
        std::vector<Tensor<double, 3>> rho_dbck_list(ijk_idx.size());

        // T_lm amplitudes
        Tensor<double, 4> T_lm("T_lm", nlmo_ijk, nlmo_ijk, ntno_ijk, ntno_ijk);
        T_lm.zero();
        for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
            int l = lmotriplet_to_lmos_[ijk][l_ijk];

            for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
                int m = lmotriplet_to_lmos_[ijk][m_ijk];
                int lm = i_j_to_ij_[l][m];
                if (lm == -1) continue;
                int lm_idx = (l_ijk > m_ijk) ? m_ijk * nlmo_ijk + l_ijk : l_ijk * nlmo_ijk + m_ijk;

                auto S_ijk_lm = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[lm]));
                S_ijk_lm = linalg::doublet(S_ijk_lm, X_pno_[lm], false, false);

                auto T_lm_ijk = linalg::triplet(S_ijk_lm, T_iajb_[lm], S_ijk_lm, false, false, true);
                ::memcpy(&T_lm(l_ijk, m_ijk, 0, 0), T_lm_ijk->get_pointer(), ntno_ijk * ntno_ijk * sizeof(double));
            }
        }

        // Erase overlap matrices from RAM

        // Factorize the contractions below to avoid storing four-index K_dble intermediates.

        // Reuse contraction buffers to avoid allocations inside the innermost loops.
        Tensor<double, 3> contraction_buffer_a("contraction_buffer_a", ntno_ijk, ntno_ijk, ntno_ijk);
        Tensor<double, 3> contraction_buffer_b("contraction_buffer_b", ntno_ijk, ntno_ijk, ntno_ijk);

        // Build the natural (l,d,m,e) density-fitting product first: grouping
        // (l,d) and (m,e) makes this call one GEMM.  Then pack the downstream
        // (l,m,d,e) order explicitly.  The temporary is scoped so only two
        // rank-four K tensors coexist, rather than the legacy three-bank set.
        Tensor<double, 4> K_lmde("K_lmde", nlmo_ijk, nlmo_ijk, ntno_ijk,
                                  ntno_ijk);
        {
            Tensor<double, 4> K_ldme_gemm("K_ldme GEMM", nlmo_ijk, ntno_ijk,
                                           nlmo_ijk, ntno_ijk);
            einsum(0.0, Indices{index::l, index::d, index::m, index::e},
                   &K_ldme_gemm, 1.0,
                   Indices{index::Q, index::l, index::d}, q_ov_[ijk],
                   Indices{index::Q, index::m, index::e}, q_ov_[ijk]);
            permute(Indices{index::l, index::m, index::d, index::e}, &K_lmde,
                    Indices{index::l, index::d, index::m, index::e},
                    K_ldme_gemm);
        }

        // Load overlap integrals from disk.

        for (int idx = 0; idx < ijk_idx.size(); ++idx) {
            int i = ijk_idx[idx];
            int ii = i_j_to_ij_[i][i];

            rho_dbck_list[idx] = Tensor<double, 3>("rho_dbck_list", ntno_ijk, ntno_ijk, ntno_ijk);
            einsum(0.0, Indices{index::d, index::b, index::c}, &rho_dbck_list[idx], 1.0, Indices{index::Q, index::d, index::b}, q_vv_t1, Indices{index::Q, index::c}, q_iv_t1_list[idx]);

            Tensor<double, 3> T_li("T_li", nlmo_ijk, ntno_ijk, ntno_ijk);
            Tensor<double, 3> U_li("U_li", nlmo_ijk, ntno_ijk, ntno_ijk);

            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                int li = i_j_to_ij_[l][i];

                auto S_ijk_li = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[li]));
                S_ijk_li = linalg::doublet(S_ijk_li, X_pno_[li], false, false);

                auto T_li_ijk = linalg::triplet(S_ijk_li, T_iajb_[li], S_ijk_li, false, false, true);
                auto U_li_ijk = linalg::triplet(S_ijk_li, Tt_iajb_[li], S_ijk_li, false, false, true);

                ::memcpy(&T_li(l_ijk, 0, 0), T_li_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
                ::memcpy(&U_li(l_ijk, 0, 0), U_li_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
            } // end l_ijk

            einsum(1.0, Indices{index::d, index::b, index::c}, &rho_dbck_list[idx], -1.0, Indices{index::l, index::d}, F_ld, Indices{index::l, index::b, index::c}, T_li);

            Tensor<double, 3> K_limd("K_limd", nlmo_ijk, nlmo_ijk, ntno_ijk);
            einsum(0.0, Indices{index::l, index::m, index::d}, &K_limd, 1.0, Indices{index::Q, index::l}, q_io_t1_list[idx], Indices{index::Q, index::m, index::d}, q_ov_[ijk]);
            Tensor<double, 3> K_limd_T("K_limd_T", nlmo_ijk, nlmo_ijk, ntno_ijk);
            permute(Indices{index::m, index::l, index::d}, &K_limd_T, Indices{index::l, index::m, index::d}, K_limd);
            einsum(1.0, Indices{index::d, index::b, index::c}, &rho_dbck_list[idx], 1.0, Indices{index::m, index::l, index::d}, K_limd_T, Indices{index::m, index::l, index::b, index::c}, T_lm);

            // Factorized evaluation of the direct term:
            // (d, b, c) <- (Q, d, b) (q_vv_t1) * (Q, l, e) (q_ov) * (l, e, c) (U_li)
            // (Q, c) (q_c_intermediate) <- (Q, l, e) (q_ov) * (l, e, c) (U_li)
            Tensor<double, 2> q_c_intermediate("q_c_intermediate", naux_ijk, ntno_ijk);
            einsum(0.0, Indices{index::Q, index::c}, &q_c_intermediate, 1.0, Indices{index::Q, index::l, index::e}, q_ov_[ijk], Indices{index::l, index::e, index::c}, U_li);
            // (d, b, c) <- (Q, d, b) (q_vv_t1) * (Q, c) (q_c_intermediate)
            einsum(1.0, Indices{index::d, index::b, index::c}, &rho_dbck_list[idx], 1.0, Indices{index::Q, index::d, index::b}, q_vv_t1, Indices{index::Q, index::c}, q_c_intermediate);

            // Factorized evaluation of the exchange-like term over auxiliary-function slices:
            // (d, b, c) <- -1.0 (Q, e, b) (q_vv_t1) * (Q, l, d) (q_ov_[ijk]) * (l, e, c) (T_li)
            for (int q_ijk = 0; q_ijk < naux_ijk; ++q_ijk) {
                TensorView<double, 2> q_vv_t1_slice = q_vv_t1(q_ijk, All, All); // (e, b)
                TensorView<double, 2> q_ov_slice = q_ov_[ijk](q_ijk, All, All); // (l, d)

                // For fixed Q: (d, b, c) <- -1.0 (e, b) (q_vv_t1_slice) * (l, d) (q_ov_slice) * (l, e, c) (T_li)
                // First step: (d, e, c) (contraction_buffer_a) <- (l, d) (q_ov_slice) * (l, e, c) (T_li)
                einsum(0.0, Indices{index::d, index::e, index::c}, &contraction_buffer_a, 1.0, Indices{index::l, index::d}, q_ov_slice, Indices{index::l, index::e, index::c}, T_li);
                // Second step: (e, d, c) (contraction_buffer_b) <- (d, e, c) contraction_buffer_a
                permute(Indices{index::e, index::d, index::c}, &contraction_buffer_b, Indices{index::d, index::e, index::c}, contraction_buffer_a);
                // Third step: (b, d, c) (contraction_buffer_a) <- -1.0 * (e, d, c) (contraction_buffer_b) * (e, b) (q_vv_t1_slice)
                einsum(0.0, Indices{index::b, index::d, index::c}, &contraction_buffer_a, -1.0, Indices{index::e, index::d, index::c}, contraction_buffer_b, Indices{index::e, index::b}, q_vv_t1_slice);
                // Final step (d, b, c) contraction_buffer_b <- (b, d, c) (contraction_buffer_a)
                permute(Indices{index::d, index::b, index::c}, &contraction_buffer_b, Indices{index::b, index::d, index::c}, contraction_buffer_a);

                rho_dbck_list[idx] += contraction_buffer_b;

                // (d, b, c) (rho_dbck_list) <- -1.0 * (Q, e, c) (q_vv_t1) * (Q, l, d) (q_ov) * (l, b, e) (T_li)
                // For fixed Q: (d, b, c) (rho_dbck_list) <- -1.0 * (e, c) (q_vv_slice) * (l, d) (q_ov_slice) * (l, b, e) (T_li)
                // (d, b, e) (contraction_buffer_a) <- (l, d) (q_ov_slice) * (l, b, e) (T_li) 
                einsum(0.0, Indices{index::d, index::b, index::e}, &contraction_buffer_a, 1.0, Indices{index::l, index::b, index::e}, T_li, Indices{index::l, index::d}, q_ov_slice);
                // (d, b, c) (rho_dbck_list) <- -1.0 * (d, b, e) * (e, c) (q_vv_t1)
                einsum(1.0, Indices{index::d, index::b, index::c}, &rho_dbck_list[idx], -1.0, Indices{index::d, index::b, index::e}, contraction_buffer_a, Indices{index::e, index::c}, q_vv_t1_slice);
            }

            // Remaining triples contributions
            auto S_ijk_ii = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[ii]));
            S_ijk_ii = linalg::doublet(S_ijk_ii, X_pno_[ii], false, false);
            rho_dbck_list[idx] += matmul_3d_einsums(rho_dbck_cont[i], S_ijk_ii, n_pno_[ii], n_tno_[ijk]);
        } // end idx

        std::vector<std::tuple<int, int, int>> long_perms = {std::make_tuple(i, j, k), std::make_tuple(i, k, j),
                                                                std::make_tuple(j, i, k), std::make_tuple(j, k, i),
                                                                std::make_tuple(k, i, j), std::make_tuple(k, j, i)};

        std::vector<std::tuple<int, int, int>> long_perms_idx = {std::make_tuple(0, 1, 2), std::make_tuple(0, 2, 1),
                                                                    std::make_tuple(1, 0, 2), std::make_tuple(1, 2, 0),
                                                                    std::make_tuple(2, 0, 1), std::make_tuple(2, 1, 0)};

        // rho_jk(l,c): manuscript Eq. (102), canonical Eqs. (39)-(40),
        // SI Eq. (S6) and Algorithm S4.
        std::vector<Tensor<double, 2>> rho_ljck_list(long_perms.size());

        for (int idx = 0; idx < long_perms.size(); ++idx) {
            auto &[i, j, k] = long_perms[idx];
            auto &[i_idx, j_idx, k_idx] = long_perms_idx[idx];
            int jk = i_j_to_ij_[j][k];

            rho_ljck_list[idx] = Tensor<double, 2>("rho_ljck", nlmo_ijk, ntno_ijk);
            einsum(0.0, Indices{index::l, index::c}, &rho_ljck_list[idx], 1.0, Indices{index::Q, index::l}, q_io_t1_list[j_idx], Indices{index::Q, index::c}, q_iv_t1_list[k_idx]);

            Tensor<double, 3> T_jm("T_jm", nlmo_ijk, ntno_ijk, ntno_ijk);
            Tensor<double, 3> T_mk("T_mk", nlmo_ijk, ntno_ijk, ntno_ijk);
            Tensor<double, 3> U_mk("U_mk", nlmo_ijk, ntno_ijk, ntno_ijk);

            for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
                int m = lmotriplet_to_lmos_[ijk][m_ijk];
                int mk = i_j_to_ij_[m][k], jm = i_j_to_ij_[j][m];

                // jm
                auto S_ijk_jm = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[jm]));
                S_ijk_jm = linalg::doublet(S_ijk_jm, X_pno_[jm], false, false);

                auto T_jm_ijk = linalg::triplet(S_ijk_jm, T_iajb_[jm], S_ijk_jm, false, false, true);
                ::memcpy(&T_jm(m_ijk, 0, 0), T_jm_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));

                // mk
                auto S_ijk_mk = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[mk]));
                S_ijk_mk = linalg::doublet(S_ijk_mk, X_pno_[mk], false, false);

                auto T_mk_ijk = linalg::triplet(S_ijk_mk, T_iajb_[mk], S_ijk_mk, false, false, true);
                auto U_mk_ijk = linalg::triplet(S_ijk_mk, Tt_iajb_[mk], S_ijk_mk, false, false, true);

                ::memcpy(&T_mk(m_ijk, 0, 0), T_mk_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
                ::memcpy(&U_mk(m_ijk, 0, 0), U_mk_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
            }

            Tensor<double, 3> K_mdlj("K_mdlj", nlmo_ijk, ntno_ijk, nlmo_ijk);
            einsum(0.0, Indices{index::m, index::d, index::l}, &K_mdlj, 1.0, Indices{index::Q, index::m, index::d}, q_ov_[ijk], Indices{index::Q, index::l}, q_io_t1_list[j_idx]);

            einsum(1.0, Indices{index::l, index::c}, &rho_ljck_list[idx], 1.0, Indices{index::m, index::d, index::l}, K_mdlj, Indices{index::m, index::d, index::c}, U_mk);

            Tensor<double, 3> K_mdlj_T("K_mdlj_T", nlmo_ijk, nlmo_ijk, ntno_ijk);
            permute(Indices{index::l, index::m, index::d}, &K_mdlj_T, Indices{index::l, index::d, index::m}, K_mdlj);
            einsum(1.0, Indices{index::l, index::c}, &rho_ljck_list[idx], -1.0, Indices{index::l, index::m, index::d}, K_mdlj_T, Indices{index::m, index::d, index::c}, T_mk);

            Tensor<double, 3> K_ldmk("K_ldmk", nlmo_ijk, ntno_ijk, nlmo_ijk);
            einsum(0.0, Indices{index::l, index::d, index::m}, &K_ldmk, 1.0, Indices{index::Q, index::l, index::d}, q_ov_[ijk], Indices{index::Q, index::m}, q_io_t1_list[k_idx]);

            Tensor<double, 3> K_ldmk_T("K_ldmk_T", nlmo_ijk, nlmo_ijk, ntno_ijk);
            permute(Indices{index::l, index::m, index::d}, &K_ldmk_T, Indices{index::l, index::d, index::m}, K_ldmk);

            einsum(1.0, Indices{index::l, index::c}, &rho_ljck_list[idx], -1.0, Indices{index::l, index::m, index::d}, K_ldmk_T, Indices{index::m, index::d, index::c}, T_jm);

            Tensor<double, 2> T_jk("T_jk", ntno_ijk, ntno_ijk);

            auto S_ijk_jk = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[jk]));
            S_ijk_jk = linalg::doublet(S_ijk_jk, X_pno_[jk], false, false);

            auto T_jk_ijk = linalg::triplet(S_ijk_jk, T_iajb_[jk], S_ijk_jk, false, false, true);
            ::memcpy(T_jk.data(), T_jk_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));

            // (l, c) [rho_ljck] <- (Q, d, c) [q_vv_t1] * (Q, l, e) [q_ov] * (e, d) [T_jk]
            for (int q_ijk = 0; q_ijk < naux_ijk; ++q_ijk) {
                TensorView<double, 2> q_vv_t1_slice = q_vv_t1(q_ijk, All, All); // (d, c)
                TensorView<double, 2> q_ov_slice = q_ov_[ijk](q_ijk, All, All); // (l, e)

                // (l, c) [rho_ljck] <- (d, c) [q_vv_t1_slice] * (l, e) [q_ov_slice] * (e, d) [T_jk]
                // (l, d) [ld_intermediate] <- (l, e) [q_ov_slice] * (e, d) [T_jk]
                Tensor<double, 2> ld_intermediate("ld_intermediate", nlmo_ijk, ntno_ijk);
                einsum(0.0, Indices{index::l, index::d}, &ld_intermediate, 1.0, Indices{index::l, index::e}, q_ov_slice, Indices{index::e, index::d}, T_jk);
                // (l, c) <- (l, d) [ld_intermediate] * (d, c) [q_vv_t1_slice]
                einsum(1.0, Indices{index::l, index::c}, &rho_ljck_list[idx], 1.0, Indices{index::l, index::d}, ld_intermediate, Indices{index::d, index::c}, q_vv_t1_slice);
            }

            for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
                int m = lmotriplet_to_lmos_[ijk][m_ijk];
                const size_t mjk_dense = triplet_key(m, j, k, naocc);
                if (!i_j_k_to_ijk_.count(mjk_dense)) continue;

                int mjk = i_j_k_to_ijk_[mjk_dense];
                auto S_ijk_mjk = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmotriplet_to_paos_[mjk]));
                S_ijk_mjk = linalg::doublet(S_ijk_mjk, X_tno_[mjk], false, false);

                // (2 t_{mjk}^{edc} - t_{mjk}^{cde} - t_{mjk}^{dec})
                Tensor<double, 3> T_mjk = matmul_3d_einsums(
                    triples_permuter_einsums(T_iajbkc_clone_[mjk], m, j, k), S_ijk_mjk, n_tno_[mjk], n_tno_[ijk]);

                Tensor<double, 3> T_kjm("T_kjm", ntno_ijk, ntno_ijk, ntno_ijk);
                permute(Indices{index::c, index::b, index::a}, &T_kjm, Indices{index::a, index::b, index::c}, T_mjk);

                Tensor<double, 3> T_jmk("T_jmk", ntno_ijk, ntno_ijk, ntno_ijk);
                permute(Indices{index::b, index::a, index::c}, &T_jmk, Indices{index::a, index::b, index::c}, T_mjk);

                T_mjk *= 2.0;
                T_mjk -= T_kjm;
                T_mjk -= T_jmk;

                // (l, c) [rho_ljck_list] <- (Q, l, d) [q_ov] * (Q, m, e) [q_ov] * (e, d, c) [T_mjk]
                TensorView<double, 3> K_mled_slice =
                    K_lmde(m_ijk, All, All, All);  // physical axes (l,e,d)
                einsum(1.0, Indices{index::l, index::c},
                       &rho_ljck_list[idx], 1.0,
                       Indices{index::l, index::e, index::d}, K_mled_slice,
                       Indices{index::e, index::d, index::c}, T_mjk);
            }
        }

        // Long-permutation W contribution: manuscript Eq. (100), canonical Eqs. (38)-(42),
        // and Algorithm 3. The sign is absorbed into the local rho contraction below.
        std::vector<Tensor<double, 3>> Wperms(long_perms.size());

        for (int idx = 0; idx < long_perms.size(); ++idx) {
            auto &[i, j, k] = long_perms[idx];
            auto &[i_idx, j_idx, k_idx] = long_perms_idx[idx];
            int ij = i_j_to_ij_[i][j];

            Wperms[idx] = Tensor("Wperm", n_tno_[ijk], n_tno_[ijk], n_tno_[ijk]);

            /// => rho_dbck contribution <= ///

            // Compute overlap between TNOs of triplet ijk and PNOs of pair ij
            auto S_ijk_ij = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[ij]));
            S_ijk_ij = linalg::doublet(S_ijk_ij, X_pno_[ij], false, false);

            auto T_ij = linalg::triplet(S_ijk_ij, T_iajb_[ij], S_ijk_ij, false, false, true);
            Tensor<double, 2> T_ij_einsums("T_ij", n_tno_[ijk], n_tno_[ijk]);
            ::memcpy(T_ij_einsums.data(), T_ij->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));

            einsum(0.0, Indices{index::a, index::b, index::c}, &Wperms[idx], 1.0, Indices{index::a, index::d}, T_ij_einsums, 
                    Indices{index::d, index::b, index::c}, rho_dbck_list[k_idx]);

            /// => rho_ljck contribution <= //
            Tensor<double, 3> T_il_einsums("T_il", nlmo_ijk, n_tno_[ijk], n_tno_[ijk]);

            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                int il = i_j_to_ij_[i][l];

                auto S_ijk_il = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[il]));
                S_ijk_il = linalg::doublet(S_ijk_il, X_pno_[il], false, false);

                auto T_il = linalg::triplet(S_ijk_il, T_iajb_[il], S_ijk_il, false, false, true);
                ::memcpy(&T_il_einsums(l_ijk, 0, 0), T_il->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
            } // end l_ijk

            einsum(1.0, Indices{index::a, index::b, index::c}, &Wperms[idx], -1.0, Indices{index::l, index::a, index::b}, T_il_einsums, 
                    Indices{index::l, index::c}, rho_ljck_list[idx]);
        }

        for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
            for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                    (*R_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) +=
                        (Wperms[0])(a_ijk, b_ijk, c_ijk) + (Wperms[1])(a_ijk, c_ijk, b_ijk) + (Wperms[2])(b_ijk, a_ijk, c_ijk) + 
                        (Wperms[3])(b_ijk, c_ijk, a_ijk) + (Wperms[4])(c_ijk, a_ijk, b_ijk) + (Wperms[5])(c_ijk, b_ijk, a_ijk);
                }
            }
        }

        // => Short-permutation V contribution: manuscript Eq. (103), canonical Eqs. (43)-(50),
        //    and Algorithm 4 <= //
        std::vector<std::tuple<int, int, int>> short_perms = {std::make_tuple(i, j, k), std::make_tuple(j, i, k), std::make_tuple(k, j, i)};
        std::vector<std::tuple<int, int, int>> short_perms_idx = {std::make_tuple(0, 1, 2), std::make_tuple(1, 0, 2), std::make_tuple(2, 1, 0)};

        // chi_i(l), chi_j(l), chi_k(l): manuscript Eq. (106), canonical Eq. (45),
        // SI Eq. (S12) and Algorithm S7.
        std::vector<Tensor<double, 1>> chi_li_list(ijk_idx.size());

        for (int idx = 0; idx < ijk_idx.size(); ++idx) {
            int i = ijk_idx[idx];

            chi_li_list[idx] = F_li_list[idx];

            Tensor<double, 3> U_mi("U_mi", nlmo_ijk, ntno_ijk, ntno_ijk);
            for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
                int m = lmotriplet_to_lmos_[ijk][m_ijk];
                int mi = i_j_to_ij_[m][i];

                auto S_ijk_mi = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[mi]));
                S_ijk_mi = linalg::doublet(S_ijk_mi, X_pno_[mi], false, false);

                auto U_mi_ijk = linalg::triplet(S_ijk_mi, Tt_iajb_[mi], S_ijk_mi, false, false, true);
                ::memcpy(&U_mi(m_ijk, 0, 0), U_mi_ijk->get_pointer(), ntno_ijk * ntno_ijk * sizeof(double));
            }

            Tensor<double, 2> chi_li_temp("chi_li_temp", naux_ijk, ntno_ijk);
            einsum(0.0, Indices{index::Q, index::d}, &chi_li_temp, 1.0, Indices{index::Q, index::m, index::e}, q_ov_[ijk], Indices{index::m, index::e, index::d}, U_mi);
            einsum(1.0, Indices{index::l}, &chi_li_list[idx], 1.0, Indices{index::Q, index::d, index::l}, q_vo, Indices{index::Q, index::d}, chi_li_temp);
        }

        // chi(a,d): manuscript Eq. (104), canonical Eq. (46), SI Eq. (S10) and Algorithm S5.
        Tensor<double, 2> chi_ad("chi_ad", ntno_ijk, ntno_ijk); {
            chi_ad = F_ad;

            Tensor<double, 4> U_lm = T_lm;
            Tensor<double, 4> T_lm_T("T_lm_T", nlmo_ijk, nlmo_ijk, ntno_ijk, ntno_ijk);
            permute(Indices{index::l, index::m, index::e, index::d}, &T_lm_T, Indices{index::l, index::m, index::d, index::e}, T_lm);
            U_lm *= 2.0;
            U_lm -= T_lm_T;
            einsum(1.0, Indices{index::a, index::d}, &chi_ad, -1.0,
                   Indices{index::m, index::l, index::e, index::a}, U_lm,
                   Indices{index::m, index::l, index::e, index::d}, K_lmde);
        }

        // chi_jk(l,m), chi_ik(l,m), chi_ji(l,m): manuscript Eq. (109), canonical Eq. (47),
        // SI Eq. (S15) and Algorithm S10.
        std::vector<Tensor<double, 2>> chi_jk_list(short_perms.size());

        for (int idx = 0; idx < short_perms.size(); ++idx) {
            auto &[i, j, k] = short_perms[idx];
            auto &[i_idx, j_idx, k_idx] = short_perms_idx[idx];
            int jk = i_j_to_ij_[j][k];

            chi_jk_list[idx] = Tensor<double, 2>("chi_jk", nlmo_ijk, nlmo_ijk);
            einsum(0.0, Indices{index::l, index::m}, &chi_jk_list[idx], 1.0, Indices{index::Q, index::l}, q_io_t1_list[j_idx], Indices{index::Q, index::m}, q_io_t1_list[k_idx]);

            auto S_ijk_jk = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[jk]));
            S_ijk_jk = linalg::doublet(S_ijk_jk, X_pno_[jk], false, false);

            auto T_jk_psi = linalg::triplet(S_ijk_jk, T_iajb_[jk], S_ijk_jk, false, false, true);

            Tensor<double, 2> T_jk("T_jk", ntno_ijk, ntno_ijk);
            ::memcpy(T_jk.data(), T_jk_psi->get_pointer(), ntno_ijk * ntno_ijk * sizeof(double));

            Tensor<double, 3> chi_jk_temp("chi_jk_temp", naux_ijk, nlmo_ijk, ntno_ijk);
            einsum(0.0, Indices{index::Q, index::l, index::e}, &chi_jk_temp, 1.0, Indices{index::Q, index::l, index::d}, q_ov_[ijk], Indices{index::d, index::e}, T_jk);

            Tensor<double, 3> chi_jk_temp_T("chi_jk_temp_T", naux_ijk, ntno_ijk, nlmo_ijk);
            permute(Indices{index::Q, index::e, index::l}, &chi_jk_temp_T, Indices{index::Q, index::l, index::e}, chi_jk_temp);

            einsum(1.0, Indices{index::l, index::m}, &chi_jk_list[idx], 1.0, Indices{index::Q, index::e, index::l}, chi_jk_temp_T, Indices{index::Q, index::e, index::m}, q_vo);
        }

        // chi(d,e,b,c): manuscript Eq. (105), canonical Eq. (48), SI Eq. (S11) and Algorithm S6.
        // The in-memory index order follows the T1-dressing convention used above.
        Tensor<double, 4> chi_dbec("chi_dbec", ntno_ijk, ntno_ijk, ntno_ijk, ntno_ijk); {
            einsum(0.0, Indices{index::d, index::b, index::e, index::c}, &chi_dbec, 1.0, Indices{index::Q, index::d, index::b}, q_vv_t1, Indices{index::Q, index::e, index::c}, q_vv_t1);

            Tensor<double, 4> chi_dbec_temp("chi_dbec_temp", ntno_ijk, ntno_ijk, ntno_ijk, ntno_ijk);
            einsum(0.0, Indices{index::d, index::e, index::b, index::c},
                   &chi_dbec_temp, 1.0,
                   Indices{index::l, index::m, index::d, index::e}, K_lmde,
                   Indices{index::l, index::m, index::b, index::c}, T_lm);

            Tensor<double, 4> chi_dbec_temp_T("chi_dbec_temp_T", ntno_ijk, ntno_ijk, ntno_ijk, ntno_ijk);
            permute(Indices{index::d, index::b, index::e, index::c}, &chi_dbec_temp_T, Indices{index::d, index::e, index::b, index::c}, chi_dbec_temp);
            chi_dbec += chi_dbec_temp_T;

            // Rearrange chi_dbec
            permute(Indices{index::d, index::e, index::b, index::c}, &chi_dbec_temp, Indices{index::d, index::b, index::e, index::c}, chi_dbec);
            chi_dbec = chi_dbec_temp;
        }

        // The remaining chi contractions need (l,e,m,d). Build that one
        // alternate order only for their shared lexical scope; at most two
        // rank-four K tensors coexist, and both contractions expose contiguous
        // GEMM dimensions.
        Tensor<double, 4> K_lemd("K_lemd", nlmo_ijk, ntno_ijk, nlmo_ijk,
                                  ntno_ijk);
        permute(Indices{index::l, index::e, index::m, index::d}, &K_lemd,
                Indices{index::l, index::m, index::d, index::e}, K_lmde);
        K_lmde = Tensor<double, 4>("released", 0, 0, 0, 0);
        Tensor<double, 4> K_ldme("K_ldme", nlmo_ijk, ntno_ijk, nlmo_ijk,
                                  ntno_ijk);
        permute(Indices{index::l, index::d, index::m, index::e}, &K_ldme,
                Indices{index::l, index::e, index::m, index::d}, K_lemd);

        // chi_i(l,d,a), chi_j(l,d,a), chi_k(l,d,a): manuscript Eq. (107), canonical Eq. (49),
        // SI Eq. (S13) and Algorithm S8.
        std::vector<Tensor<double, 3>> chi_lida_list(short_perms.size());

        for (int idx = 0; idx < ijk_idx.size(); ++idx) {
            int i = ijk_idx[idx];

            chi_lida_list[idx] = Tensor<double, 3>("chi_lida", nlmo_ijk, ntno_ijk, ntno_ijk);
            einsum(0.0, Indices{index::l, index::d, index::a}, &chi_lida_list[idx], 1.0, Indices{index::Q, index::l}, q_io_t1_list[idx], Indices{index::Q, index::d, index::a}, q_vv_t1);

            Tensor<double, 3> T_im("T_im", nlmo_ijk, ntno_ijk, ntno_ijk);
            for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
                int m = lmotriplet_to_lmos_[ijk][m_ijk];
                int im = i_j_to_ij_[i][m];

                auto S_ijk_im = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[im]));
                S_ijk_im = linalg::doublet(S_ijk_im, X_pno_[im], false, false);

                auto T_im_ijk = linalg::triplet(S_ijk_im, T_iajb_[im], S_ijk_im, false, false, true);
                ::memcpy(&T_im(m_ijk, 0, 0), T_im_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
            }

            einsum(1.0, Indices{index::l, index::d, index::a},
                   &chi_lida_list[idx], -1.0,
                   Indices{index::m, index::e, index::l, index::d}, K_lemd,
                   Indices{index::m, index::e, index::a}, T_im);
        }

        // tilde-chi_i(l,d,a), tilde-chi_j(l,d,a), tilde-chi_k(l,d,a): manuscript Eq. (108),
        // canonical Eq. (50), SI Eq. (S14) and Algorithm S9.
        std::vector<Tensor<double, 3>> chi_ldai_list(short_perms.size());
        for (int idx = 0; idx < ijk_idx.size(); ++idx) {
            int i = ijk_idx[idx];

            chi_ldai_list[idx] = Tensor<double, 3>("chi_ldai", nlmo_ijk, ntno_ijk, ntno_ijk);
            einsum(0.0, Indices{index::l, index::d, index::a}, &chi_ldai_list[idx], 1.0, Indices{index::Q, index::l, index::d}, q_ov_[ijk], Indices{index::Q, index::a}, q_iv_t1_list[idx]);

            Tensor<double, 3> T_mi("T_mi", nlmo_ijk, ntno_ijk, ntno_ijk);
            Tensor<double, 3> U_mi("U_mi", nlmo_ijk, ntno_ijk, ntno_ijk);
            for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
                int m = lmotriplet_to_lmos_[ijk][m_ijk];
                int mi = i_j_to_ij_[m][i];

                auto S_ijk_mi = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmopair_to_paos_[mi]));
                S_ijk_mi = linalg::doublet(S_ijk_mi, X_pno_[mi], false, false);

                auto T_mi_ijk = linalg::triplet(S_ijk_mi, T_iajb_[mi], S_ijk_mi, false, false, true);
                auto U_mi_ijk = linalg::triplet(S_ijk_mi, Tt_iajb_[mi], S_ijk_mi, false, false, true);
                ::memcpy(&T_mi(m_ijk, 0, 0), T_mi_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
                ::memcpy(&U_mi(m_ijk, 0, 0), U_mi_ijk->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
            }

            // (l, e, m, d) (m, a, e) => (a, d, l)
            einsum(1.0, Indices{index::l, index::d, index::a},
                   &chi_ldai_list[idx], -1.0,
                   Indices{index::l, index::d, index::m, index::e}, K_lemd,
                   Indices{index::m, index::e, index::a}, T_mi);
            // (l, d, m, e) (m, a, e) => (a, d, l)
            einsum(1.0, Indices{index::l, index::d, index::a},
                   &chi_ldai_list[idx], 1.0,
                   Indices{index::l, index::d, index::m, index::e}, K_ldme,
                   Indices{index::m, index::e, index::a}, U_mi);
        }

        std::vector<Tensor<double, 3>> Vperms(short_perms.size());

        for (int idx = 0; idx < short_perms.size(); ++idx) {
            auto &[i, j, k] = short_perms[idx];
            auto &[i_idx, j_idx, k_idx] = short_perms_idx[idx];
            int ii = i_j_to_ij_[i][i];

            Vperms[idx] = Tensor<double, 3>("V_iajbkc", ntno_ijk, ntno_ijk, ntno_ijk);

            // T_ijk contributions
            Tensor<double, 3> T_ijk = triples_permuter_einsums(T_iajbkc_clone_[ijk], i, j, k);
            einsum(0.0, Indices{index::a, index::b, index::c}, &Vperms[idx], 1.0, Indices{index::a, index::d}, chi_ad, Indices{index::d, index::b, index::c}, T_ijk);

            einsum(1.0, Indices{index::a, index::b, index::c}, &Vperms[idx], 1.0, Indices{index::a, index::d, index::e}, T_ijk, Indices{index::d, index::e, index::b, index::c}, chi_dbec);

            // T_ljk contributions
            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];
                const size_t ljk_dense = triplet_key(l, j, k, naocc);
                if (!i_j_k_to_ijk_.count(ljk_dense)) continue;
                int ljk = i_j_k_to_ijk_[ljk_dense];

                auto S_ijk_ljk = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmotriplet_to_paos_[ljk]));
                S_ijk_ljk = linalg::doublet(S_ijk_ljk, X_tno_[ljk], false, false);

                Tensor<double, 3> T_ljk = matmul_3d_einsums(
                    triples_permuter_einsums(T_iajbkc_clone_[ljk], l, j, k), S_ijk_ljk, n_tno_[ljk], n_tno_[ijk]);

                Tensor<double, 3> T_ljk_clone = T_ljk;
                T_ljk_clone *= -(chi_li_list[i_idx])(l_ijk);
                Vperms[idx] += T_ljk_clone;

                // chi_lida contributions
                Tensor<double, 2> chi_lida_temp = chi_lida_list[i_idx](l_ijk, All, All);
                einsum(1.0, Indices{index::a, index::b, index::c}, &Vperms[idx], -1.0, Indices{index::d, index::a}, chi_lida_temp, Indices{index::d, index::b, index::c}, T_ljk);
                Tensor<double, 3> Vperm_temp("Vperm_temp", ntno_ijk, ntno_ijk, ntno_ijk);
                permute(Indices{index::b, index::a, index::c}, &Vperm_temp, Indices{index::a, index::b, index::c}, T_ljk);
                Tensor<double, 3> Vperm_temp2("Vperm_temp2", ntno_ijk, ntno_ijk, ntno_ijk);
                einsum(0.0, Indices{index::b, index::a, index::c}, &Vperm_temp2, 1.0, Indices{index::d, index::b}, chi_lida_temp, Indices{index::d, index::a, index::c}, Vperm_temp);
                permute(Indices{index::a, index::b, index::c}, &Vperm_temp, Indices{index::b, index::a, index::c}, Vperm_temp2);
                Vperms[idx] -= Vperm_temp;
                einsum(1.0, Indices{index::a, index::b, index::c}, &Vperms[idx], -1.0, Indices{index::d, index::c}, chi_lida_temp, Indices{index::a, index::b, index::d}, T_ljk);

                // chi_ldai contributions
                Tensor<double, 3> T_jlk("T_jlk", ntno_ijk, ntno_ijk, ntno_ijk);
                permute(Indices{index::b, index::a, index::c}, &T_jlk, Indices{index::a, index::b, index::c}, T_ljk);

                Tensor<double, 3> T_kjl("T_kjl", ntno_ijk, ntno_ijk, ntno_ijk);
                permute(Indices{index::c, index::b, index::a}, &T_kjl, Indices{index::a, index::b, index::c}, T_ljk);

                T_ljk *= 2.0;
                T_ljk -= T_jlk;
                T_ljk -= T_kjl;

                Tensor<double, 2> chi_ldai_temp = chi_ldai_list[i_idx](l_ijk, All, All);
                einsum(1.0, Indices{index::a, index::b, index::c}, &Vperms[idx], 1.0, Indices{index::d, index::a}, chi_ldai_temp, Indices{index::d, index::b, index::c}, T_ljk);
            }

            // T_ilm contributions
            for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                int l = lmotriplet_to_lmos_[ijk][l_ijk];

                for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
                    int m = lmotriplet_to_lmos_[ijk][m_ijk];
                    const size_t ilm_dense = triplet_key(i, l, m, naocc);
                    if (!i_j_k_to_ijk_.count(ilm_dense)) continue;
                    int ilm = i_j_k_to_ijk_[ilm_dense];
                    int lm_idx = (l_ijk > m_ijk) ? m_ijk * nlmo_ijk + l_ijk : l_ijk * nlmo_ijk + m_ijk;

                    auto S_ijk_ilm = submatrix_cols(*S_ijk, index_list(triplet_ext_domain, lmotriplet_to_paos_[ilm]));
                    S_ijk_ilm = linalg::doublet(S_ijk_ilm, X_tno_[ilm], false, false);

                    Tensor<double, 3> T_ilm = matmul_3d_einsums(
                        triples_permuter_einsums(T_iajbkc_clone_[ilm], i, l, m), S_ijk_ilm, n_tno_[ilm], n_tno_[ijk]);

                    T_ilm *= (chi_jk_list[idx])(l_ijk, m_ijk);
                    Vperms[idx] += T_ilm;
                } // end m_ijk
            } // end l_ijk
        } // end idx

        if (disk_ints_) {
            q_ov_[ijk] = Tensor<double, 3>(q_ov_[ijk].name(), 0, 0, 0);
            q_vv_[ijk] = Tensor<double, 3>(q_vv_[ijk].name(), 0, 0, 0);
        }

        // Flush Vperms
        for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
            for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                    (*R_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) +=
                        (Vperms[0])(a_ijk, b_ijk, c_ijk) + (Vperms[1])(b_ijk, a_ijk, c_ijk) + (Vperms[2])(c_ijk, b_ijk, a_ijk);
                } // end c_ijk
            } // end b_ijk
        } // end a_ijk

    } // end ijk
}

void DLPNOCCSDT::lccsdt_iterations() {

    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();
    int n_lmo_triplets = ijk_to_i_j_k_.size();

    // Thread and OMP Parallel Info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // => Initialize Residuals and Amplitude <= //

    std::vector<SharedMatrix> R_ia(naocc);
    std::vector<SharedMatrix> Rn_iajb(n_lmo_pairs);
    std::vector<SharedMatrix> R_iajb(n_lmo_pairs);
    std::vector<SharedMatrix> R_iajbkc(n_lmo_triplets);

#pragma omp parallel for
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        R_ia[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
    }

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        R_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
        Rn_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
    }

#pragma omp parallel for
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        R_iajbkc[ijk] = std::make_shared<Matrix>(n_tno_[ijk], n_tno_[ijk] * n_tno_[ijk]);
    }

    std::vector<std::vector<SharedMatrix>> R_ia_buffer(nthreads);
    std::vector<std::vector<SharedMatrix>> R_iajb_buffer(nthreads);

#pragma omp parallel for
    for (int thread = 0; thread < nthreads; ++thread) {
        R_ia_buffer[thread].resize(naocc);
        R_iajb_buffer[thread].resize(n_lmo_pairs);

        for (int i = 0; i < naocc; ++i) {
            int ii = i_j_to_ij_[i][i];
            R_ia_buffer[thread][i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        }

        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            R_iajb_buffer[thread][ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
        }
    }

    // Amplitude intermediates
    T_n_ijk_.resize(n_lmo_triplets);
    T_iajbkc_clone_.resize(n_lmo_triplets);
    U_iajbkc_.resize(n_lmo_triplets);

    // LCCSDT iterations

    outfile->Printf("\n  ==> Local CCSDT <==\n\n");

    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                       Corr. Energy    Delta E     RMS R1     RMS R2     RMS R3     Time (s)\n");

    int iteration = 1, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, e_weak = 0.0, r_curr1 = 1.0, r_curr2 = 1.0, r_curr3 = 1.0;
    bool e_converged = false, r_converged = false;
    const int n_microiterations = options_.get_int("DLPNO_TRIPLES_MICROITERATIONS");

    DIISManager diis = DIISManager(options_.get_int("DIIS_MAX_VECS"), "LCCSDT DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    while (!(e_converged && r_converged)) {
        // RMS of residual per LMO orbital, for assessing convergence
        std::vector<double> R_ia_rms(naocc, 0.0);
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);
        // RMS of residual per LMO triplet, for assessing convergence
        std::vector<double> R_iajbkc_rms(n_lmo_triplets, 0.0);

        std::time_t time_start = std::time(nullptr);

        for (int miter = 0; miter < n_microiterations; ++miter) {

            // Create T_n_ij and T_n_ijk intermediates
    #pragma omp parallel for schedule(dynamic, 1)
            for (int ij = 0; ij < n_lmo_pairs; ++ij) {
                auto &[i, j] = ij_to_i_j_[ij];

                int nlmo_ij = lmopair_to_lmos_[ij].size();
                int npno_ij = n_pno_[ij];
                int ij_idx = (i <= j) ? ij : ij_to_ji_[ij];

                T_n_ij_[ij] = std::make_shared<Matrix>(nlmo_ij, npno_ij);

                for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                    int n = lmopair_to_lmos_[ij][n_ij];
                    int nn = i_j_to_ij_[n][n];
                    auto T_n_temp = linalg::doublet(S_pno_ij_nn_[ij_idx][n], T_ia_[n], false, false);

                    for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                        (*T_n_ij_[ij])(n_ij, a_ij) = (*T_n_temp)(a_ij, 0);
                    } // end a_ij
                } // end n_ij
            }

    #pragma omp parallel for schedule(dynamic, 1)
            for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                int ijk = sorted_triplets_[ijk_sorted];
                auto &[i, j, k] = ijk_to_i_j_k_[ijk];
                int nlmo_ijk = lmotriplet_to_lmos_[ijk].size();

                T_n_ijk_[ijk] = Tensor<double, 2>("T_n_ijk", nlmo_ijk, n_tno_[ijk]);

                for (int l_ijk = 0; l_ijk < nlmo_ijk; ++l_ijk) {
                    int l = lmotriplet_to_lmos_[ijk][l_ijk];
                    int ll = i_j_to_ij_[l][l];

                    auto S_ijk_ll = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk], lmopair_to_paos_[ll]);
                    S_ijk_ll = linalg::triplet(X_tno_[ijk], S_ijk_ll, X_pno_[ll], true, false, false);

                    auto T_l_temp = linalg::doublet(S_ijk_ll, T_ia_[l]);

                    for (int a_ijk = 0; a_ijk < n_tno_[ijk]; ++a_ijk) {
                        (T_n_ijk_[ijk])(l_ijk, a_ijk) = (*T_l_temp)(a_ijk, 0);
                    }
                } // end l_ijk
            } // end ijk

            // Create T_iajbkc_clone intermediate
    #pragma omp parallel for schedule(dynamic, 1)
            for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                int ijk = sorted_triplets_[ijk_sorted];
                auto &[i, j, k] = ijk_to_i_j_k_[ijk];

                T_iajbkc_clone_[ijk] = Tensor<double, 3>("T_ijk", n_tno_[ijk], n_tno_[ijk], n_tno_[ijk]);
                ::memcpy(T_iajbkc_clone_[ijk].data(), T_iajbkc_[ijk]->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));

                U_iajbkc_[ijk] = triples_spin_summation(T_iajbkc_clone_[ijk]);
            } // end ijk

            // Canonical t1 dressing and its DF factorization, manuscript Eqs. (20)-(31); the local forms are
            // detailed in SI Algorithms S1-S2.
            t1_ints();
            t1_fock();

            // Project the orbital triples tensor through the spin-summed representation
            // using Matthews and Stanton's spin-summation procedure and Eq. (27).
        #pragma omp parallel for schedule(dynamic, 1)
            for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                int ijk = sorted_triplets_[ijk_sorted];
                auto &[i, j, k] = ijk_to_i_j_k_[ijk];

                T_iajbkc_clone_[ijk] = triples_spin_desummation(triples_spin_summation(T_iajbkc_clone_[ijk]));
                U_iajbkc_[ijk] = triples_spin_summation(T_iajbkc_clone_[ijk]);
                ::memcpy(T_iajbkc_[ijk]->get_pointer(), T_iajbkc_clone_[ijk].data(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
            }

            // Triples residual: canonical Eqs. (37)-(50), DLPNO Eqs. (99)-(110).
            timer_on("DLPNO-CCSDT : R_iajbkc");
            if (miter == 0) {
                compute_R_iajbkc(R_iajbkc);

                // spin adapt and then de-adapt triples residual
        #pragma omp parallel for schedule(dynamic, 1)
                for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                    int ijk = sorted_triplets_[ijk_sorted];
                    auto &[i, j, k] = ijk_to_i_j_k_[ijk];

                    Tensor<double, 3> R3_spinad("R3_spinad", n_tno_[ijk], n_tno_[ijk], n_tno_[ijk]);
                    ::memcpy(R3_spinad.data(), R_iajbkc[ijk]->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
                    R3_spinad = triples_spin_desummation(triples_spin_summation(R3_spinad));
                    ::memcpy(R_iajbkc[ijk]->get_pointer(), R3_spinad.data(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
                }

                // Update triples amplitude
                r_curr3 = 0.0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+ : r_curr3)
                for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                    int ijk = sorted_triplets_[ijk_sorted];
                    auto &[i, j, k] = ijk_to_i_j_k_[ijk];
                    double alpha = (fabs(R_iajbkc[ijk]->rms()) > fabs(R_iajbkc_rms[ijk])) ? damping_ratio_ : 0.0;

                    for (int a_ijk = 0; a_ijk < n_tno_[ijk]; ++a_ijk) {
                        for (int b_ijk = 0; b_ijk < n_tno_[ijk]; ++b_ijk) {
                            for (int c_ijk = 0; c_ijk < n_tno_[ijk]; ++c_ijk) {
                                (*T_iajbkc_[ijk])(a_ijk, b_ijk * n_tno_[ijk] + c_ijk) -= (1.0 - alpha) * (*R_iajbkc[ijk])(a_ijk, b_ijk * n_tno_[ijk] + c_ijk) /
                                                    (e_tno_[ijk]->get(a_ijk) + e_tno_[ijk]->get(b_ijk) + e_tno_[ijk]->get(c_ijk) - F_lmo_->get(i,i) - F_lmo_->get(j,j) - F_lmo_->get(k,k));
                            }
                        }
                    }
                    // Copy information
                    ::memcpy(T_iajbkc_clone_[ijk].data(), T_iajbkc_[ijk]->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
                    U_iajbkc_[ijk] = triples_spin_summation(T_iajbkc_clone_[ijk]);

                    R_iajbkc_rms[ijk] = R_iajbkc[ijk]->rms();
                    r_curr3 += R_iajbkc_rms[ijk] * R_iajbkc_rms[ijk];
                }
                r_curr3 = std::sqrt(r_curr3 / n_lmo_triplets);
            } // end if
            timer_off("DLPNO-CCSDT : R_iajbkc");

            // Doubles residual: canonical Eqs. (33)-(34), DLPNO Eqs. (87)-(90).
            timer_on("DLPNO-CCSDT : R_iajb");
            compute_R_iajb_triples(R_iajb, Rn_iajb, R_iajb_buffer);

            // Update doubles amplitude
            r_curr2 = 0.0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : r_curr2)
            for (int ij = 0; ij < n_lmo_pairs; ++ij) {
                auto &[i, j] = ij_to_i_j_[ij];
                double alpha = (fabs(R_iajb[ij]->rms()) > fabs(R_iajb_rms[ij])) ? damping_ratio_ : 0.0;

                for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                    for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                        (*T_iajb_[ij])(a_ij, b_ij) -= (1.0 - alpha) * (*R_iajb[ij])(a_ij, b_ij) / 
                                        (e_pno_[ij]->get(a_ij) + e_pno_[ij]->get(b_ij) - F_lmo_->get(i,i) - F_lmo_->get(j,j));
                    }
                }
                Tt_iajb_[ij] = T_iajb_[ij]->clone();
                Tt_iajb_[ij]->scale(2.0);
                Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());

                R_iajb_rms[ij] = R_iajb[ij]->rms();
                r_curr2 += R_iajb_rms[ij] * R_iajb_rms[ij];
            }
            r_curr2 = std::sqrt(r_curr2 / n_lmo_pairs);
            timer_off("DLPNO-CCSDT : R_iajb");

            // Singles residual: canonical Eq. (32), DLPNO Eq. (80).
            timer_on("DLPNO-CCSDT : R_ia");
            compute_R_ia_triples(R_ia, R_ia_buffer);

            // Update singles amplitude
            r_curr1 = 0.0;
    #pragma omp parallel for reduction(+ : r_curr1)
            for (int i = 0; i < naocc; ++i) {
                int ii = i_j_to_ij_[i][i];
                double alpha = (fabs(R_ia[i]->rms()) > fabs(R_ia_rms[i])) ? damping_ratio_ : 0.0;

                for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                    (*T_ia_[i])(a_ii, 0) -= (1.0 - alpha) * (*R_ia[i])(a_ii, 0) / (e_pno_[ii]->get(a_ii) - F_lmo_->get(i,i));
                }
                R_ia_rms[i] = R_ia[i]->rms();
                r_curr1 += R_ia_rms[i] * R_ia_rms[i];
            }
            r_curr1 = std::sqrt(r_curr1 / naocc); 
            timer_off("DLPNO-CCSDT : R_ia");
        } // end miter

        // DIIS Extrapolation
        std::vector<SharedMatrix> T_vecs;
        T_vecs.reserve(T_ia_.size() + T_iajb_.size() + T_iajbkc_.size());
        T_vecs.insert(T_vecs.end(), T_ia_.begin(), T_ia_.end());
        T_vecs.insert(T_vecs.end(), T_iajb_.begin(), T_iajb_.end());
        T_vecs.insert(T_vecs.end(), T_iajbkc_.begin(), T_iajbkc_.end());

        std::vector<SharedMatrix> R_vecs;
        R_vecs.reserve(R_ia.size() + R_iajb.size() + R_iajbkc.size());
        R_vecs.insert(R_vecs.end(), R_ia.begin(), R_ia.end());
        R_vecs.insert(R_vecs.end(), R_iajb.begin(), R_iajb.end());
        R_vecs.insert(R_vecs.end(), R_iajbkc.begin(), R_iajbkc.end());

        auto T_vecs_flat = flatten_mats(T_vecs);
        auto R_vecs_flat = flatten_mats(R_vecs);

        if (iteration == 1) {
            diis.set_error_vector_size(R_vecs_flat);
            diis.set_vector_size(T_vecs_flat);
        }

        diis.add_entry(R_vecs_flat.get(), T_vecs_flat.get());
        diis.extrapolate(T_vecs_flat.get());

        copy_flat_mats(T_vecs_flat, T_vecs);

        // Evaluate the canonical CCSDT energy of Eq. (51) in the local representation;
        // the final DLPNO energy assembly is given by manuscript Eqs. (111)-(113).
        e_prev = e_curr;
        e_curr = 0.0;
        e_weak = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr, e_weak)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];
            int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j];
            int ij_idx = (i < j) ? ij : ij_to_ji_[ij];

            // Update anti-symmetrized amplitudes
            Tt_iajb_[ij] = T_iajb_[ij]->clone();
            Tt_iajb_[ij]->scale(2.0);
            Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());

            auto tau = T_iajb_[ij]->clone();
            auto S_ij_ii = S_pno_ij_nn_[ij_idx][i];
            auto S_ij_jj = S_pno_ij_nn_[ij_idx][j];
            auto tia_temp = linalg::doublet(S_ij_ii, T_ia_[i]);
            auto tjb_temp = linalg::doublet(S_ij_jj, T_ia_[j]);

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    (*tau)(a_ij, b_ij) += (*tia_temp)(a_ij, 0) * (*tjb_temp)(b_ij, 0);
                } // end b_ij
            } // end a_ij

            double e_ij = tau->vector_dot(L_iajb_[ij]);

            e_curr += e_ij;
            if (i_j_to_ij_strong_[i][j] == -1) e_weak += e_ij;
        }

        r_converged = fabs(r_curr1) < options_.get_double("R_CONVERGENCE");
        r_converged &= fabs(r_curr2) < options_.get_double("R_CONVERGENCE");
        r_converged &= fabs(r_curr3) < options_.get_double("R_CONVERGENCE");
        e_converged = fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE");

        e_lccsdt_ = e_curr - e_weak;
        de_weak_ = e_weak;

        std::time_t time_stop = std::time(nullptr);

        outfile->Printf("  @LCCSDT iter %3d: %16.12f %10.3e %10.3e %10.3e %10.3e %8d\n", iteration, e_curr, e_curr - e_prev, r_curr1, r_curr2, r_curr3, (int)time_stop - (int)time_start);

        ++iteration;

        if (iteration > max_iteration) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }
}

double DLPNOCCSDT::compute_energy() {

    // Run DLPNO-CCSD(T) as initial step
    DLPNOCCSD_T::compute_energy();

    timer_on("DLPNO-CCSDT");

    einsums::profile::initialize();

    print_header();

    if (write_qia_pno_) {
        psio_->open(PSIF_DLPNO_QIA_PNO, PSIO_OPEN_OLD);
    }

    if (write_qab_pno_) {
        psio_->open(PSIF_DLPNO_QAB_PNO, PSIO_OPEN_OLD);
    }

    damping_ratio_ = options_.get_double("TRIPLES_DAMPING_RATIO");

    int n_lmo_triplets = ijk_to_i_j_k_.size();

    // Sort list of triplets based on number of TNOs (for parallel efficiency)
    std::vector<std::pair<int, int>> ijk_tnos(n_lmo_triplets);

#pragma omp parallel for
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        ijk_tnos[ijk] = std::make_pair(ijk, n_tno_[ijk]);
    }

    std::sort(ijk_tnos.begin(), ijk_tnos.end(), [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return (a.second > b.second);
    });

    sorted_triplets_.resize(n_lmo_triplets);
#pragma omp parallel for
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        sorted_triplets_[ijk] = ijk_tnos[ijk].first;
    }

    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    timer_on("DLPNO-CCSDT : Estimate Memory");
    estimate_memory();
    timer_off("DLPNO-CCSDT : Estimate Memory");

    // estimate_memory() may activate disk-backed storage, so open these files only
    // after the automatic memory decision has been made.
    if (disk_ints_) {
        psio_->open(PSIF_DLPNO_QIA_TNO, PSIO_OPEN_NEW);
        psio_->open(PSIF_DLPNO_QAB_TNO, PSIO_OPEN_NEW);
    }

    timer_on("DLPNO-CCSDT : Compute Integrals");
    compute_integrals();
    timer_off("DLPNO-CCSDT : Compute Integrals");

    // Compute DLPNO-CCSDT energy
    timer_on("DLPNO-CCSDT : LCCSDT iterations");
    lccsdt_iterations();
    timer_off("DLPNO-CCSDT : LCCSDT iterations");

    if (write_qia_pno_) {
        // Integrals deleted if CCSDT is the last iterative CC method used in sequence
        if (algorithm_ != DLPNOMethod::CCSDTQ) {
            psio_->close(PSIF_DLPNO_QIA_PNO, 0);
        } else {
            psio_->close(PSIF_DLPNO_QIA_PNO, 1);
        }
    }

    if (write_qab_pno_) {
        // Integrals deleted if CCSDT is the last iterative CC method used in sequence
        if (algorithm_ != DLPNOMethod::CCSDTQ) {
            psio_->close(PSIF_DLPNO_QAB_PNO, 0);
        } else {
            psio_->close(PSIF_DLPNO_QAB_PNO, 1);
        }
    }

    if (disk_ints_) {
        // Integrals deleted if CCSDT is the last iterative CC method used in sequence
        if (algorithm_ != DLPNOMethod::CCSDTQ) {
            psio_->close(PSIF_DLPNO_QIA_TNO, 0);
            psio_->close(PSIF_DLPNO_QAB_TNO, 0);
        } else {
            psio_->close(PSIF_DLPNO_QIA_TNO, 1);
            psio_->close(PSIF_DLPNO_QAB_TNO, 1);
        }
    }

    einsums::profile::finalize();

    timer_off("DLPNO-CCSDT");

    // de_lccsd_t_screened_ is the screened-triplet component carried forward from the
    // preceding DLPNO-(T) calculation. It is not a separate TNO-rank correction here.
    double e_scf = variables_["SCF TOTAL ENERGY"];
    double e_ccsdt_corr = e_lccsdt_ + de_lccsd_t_screened_ + de_weak_ + de_lmp2_eliminated_ + de_dipole_ +
                          de_pno_total_;
    double e_ccsdt_total = e_scf + e_ccsdt_corr;

    set_scalar_variable("CCSDT CORRELATION ENERGY", e_ccsdt_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsdt_corr);
    set_scalar_variable("CCSDT TOTAL ENERGY", e_ccsdt_total);
    set_scalar_variable("CURRENT ENERGY", e_ccsdt_total);
    set_scalar_variable("T CORRECTION ENERGY",
                        e_ccsdt_total - variables_["CCSD TOTAL ENERGY"]);

    print_results();

    return e_ccsdt_total;
}

void DLPNOCCSDT::print_results() {

    int naocc = i_j_to_ij_.size();
    double t1diag = 0.0;
#pragma omp parallel for reduction(+ : t1diag)
    for (int i = 0; i < naocc; ++i) {
        t1diag += T_ia_[i]->vector_dot(T_ia_[i]);
    }
    t1diag = std::sqrt(t1diag / (2.0 * naocc));
    outfile->Printf("\n  T1 Diagnostic: %8.8f \n", t1diag);
    // Lee and Taylor's conventional 0.02 closed-shell threshold was defined from CCSD
    // singles amplitudes (Int. J. Quantum Chem. 36, 199, 1989; DOI: 10.1002/qua.560360824).
    // No distinct, generally accepted CCSDT-specific cutoff is available, so retain 0.02
    // as a conservative warning about the single-reference description.
    constexpr double closed_shell_t1_warning = 0.02;
    if (t1diag > closed_shell_t1_warning) {
        outfile->Printf("    WARNING: T1 Diagnostic is greater than 0.02; single-reference coupled-cluster "
                        "results may be unreliable.\n");
    }
    set_scalar_variable("CC T1 DIAGNOSTIC", t1diag);

    double e_total = e_lccsdt_ + de_lccsd_t_screened_ + de_weak_ + de_lmp2_eliminated_ + de_dipole_ +
                     de_pno_total_;
    const double e_ccsdt_total = variables_["SCF TOTAL ENERGY"] + e_total;
    const double ccsdt_minus_ccsd = e_ccsdt_total - variables_["CCSD TOTAL ENERGY"];
    const double ccsdt_minus_ccsd_t = e_ccsdt_total - variables_["CCSD(T) TOTAL ENERGY"];

    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-CCSDT Correlation Energy: %16.12f \n", e_total);
    outfile->Printf("    LCCSDT Correlation Energy:          %16.12f \n", e_lccsdt_);
    outfile->Printf("    Screened Triplets Contribution:     %16.12f \n", de_lccsd_t_screened_);
    outfile->Printf("    CCSDT - CCSD Energy:                %16.12f \n", ccsdt_minus_ccsd);
    outfile->Printf("    CCSDT - CCSD(T) Energy:             %16.12f \n", ccsdt_minus_ccsd_t);
    outfile->Printf("\n\n  @Total DLPNO-CCSDT Energy: %16.12f \n", e_ccsdt_total);
}

#endif  // USING_Einsums

}  // namespace dlpno
}  // namespace psi
