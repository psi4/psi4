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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dlpno {

DLPNOCCSDT_Q::DLPNOCCSDT_Q(SharedWavefunction ref_wfn, Options &options) : DLPNOCCSDT(ref_wfn, options) {}
DLPNOCCSDT_Q::~DLPNOCCSDT_Q() {}

Tensor<double, 4> DLPNOCCSDT_Q::matmul_4d(const Tensor<double, 4> &A, const SharedMatrix &X, int dim_old, int dim_new) {
    /* Performs the operation A'[i,j,k,l] = A[I, J, K, L] * X[i, I] * X[j, J] * X[k, K] * X[l, L] for tesseract 4d tensors */

    // TODO: Change this into a TensorView
    Tensor<double, 2> Xview("Xview", dim_new, dim_old);
    ::memcpy(Xview.data(), X->get_pointer(), dim_new * dim_old * sizeof(double));

    Tensor<double, 4> A_new1("A_new1", dim_old, dim_old, dim_old, dim_new);
    einsum(0.0, Indices{index::I, index::J, index::K, index::l}, &A_new1, 1.0, Indices{index::I, index::J, index::K, index::L}, A, Indices{index::l, index::L}, Xview);

    Tensor<double, 4> A_new2("A_new2", dim_new, dim_old, dim_old, dim_new);
    einsum(0.0, Indices{index::i, index::J, index::K, index::l}, &A_new2, 1.0, Indices{index::I, index::J, index::K, index::l}, A_new1, Indices{index::i, index::I}, Xview);

    Tensor<double, 4> A_new3("A_new3", dim_old, dim_new, dim_new, dim_old);
    permute(Indices{index::J, index::i, index::l, index::K}, &A_new3, Indices{index::i, index::J, index::K, index::l}, A_new2);

    Tensor<double, 4> A_new4("A_new4", dim_old, dim_new, dim_new, dim_new);
    einsum(0.0, Indices{index::J, index::i, index::l, index::k}, &A_new4, 1.0, Indices{index::J, index::i, index::l, index::K}, A_new3, Indices{index::k, index::K}, Xview);

    Tensor<double, 4> A_new5("A_new5", dim_new, dim_new, dim_new, dim_new);
    einsum(0.0, Indices{index::j, index::i, index::l, index::k}, &A_new5, 1.0, Indices{index::J, index::i, index::l, index::k}, A_new4, Indices{index::j, index::J}, Xview);

    Tensor<double, 4> A_new("A_new", dim_new, dim_new, dim_new, dim_new);
    permute(Indices{index::i, index::j, index::k, index::l}, &A_new, Indices{index::j, index::i, index::l, index::k}, A_new5);

    return A_new;
}

Tensor<double, 4> DLPNOCCSDT_Q::quadruples_permuter(const Tensor<double, 4>& X, int i, int j, int k, int l) {

    Tensor<double, 4> Xperm = X;

    if (i <= j && j <= l && l <= k) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::a, index::b, index::d, index::c}, X);
    } else if (i <= k && k <= j && j <= l) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::a, index::c, index::b, index::d}, X);
    } else if (i <= k && k <= l && l <= j) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::a, index::c, index::d, index::b}, X);
    } else if (i <= l && l <= j && j <= k) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::a, index::d, index::b, index::c}, X);
    } else if (i <= l && l <= k && k <= j) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::a, index::d, index::c, index::b}, X);
    } else if (j <= i && i <= k && k <= l) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::b, index::a, index::c, index::d}, X);
    } else if (j <= i && i <= l && l <= k) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::b, index::a, index::d, index::c}, X);
    } else if (j <= k && k <= i && i <= l) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::b, index::c, index::a, index::d}, X);
    } else if (j <= k && k <= l && l <= i) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::b, index::c, index::d, index::a}, X);
    } else if (j <= l && l <= i && i <= k) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::b, index::d, index::a, index::c}, X);
    } else if (j <= l && l <= k && k <= i) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::b, index::d, index::c, index::a}, X);
    } else if (k <= i && i <= j && j <= l) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::c, index::a, index::b, index::d}, X);
    } else if (k <= i && i <= l && l <= j) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::c, index::a, index::d, index::b}, X);
    } else if (k <= j && j <= i && i <= l) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::c, index::b, index::a, index::d}, X);
    } else if (k <= j && j <= l && l <= i) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::c, index::b, index::d, index::a}, X);
    } else if (k <= l && l <= i && i <= j) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::c, index::d, index::a, index::b}, X);
    } else if (k <= l && l <= j && j <= i) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::c, index::d, index::b, index::a}, X);
    } else if (l <= i && i <= j && j <= k) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::d, index::a, index::b, index::c}, X);
    } else if (l <= i && i <= k && k <= j) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::d, index::a, index::c, index::b}, X);
    } else if (l <= j && j <= i && i <= k) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::d, index::b, index::a, index::c}, X);
    } else if (l <= j && j <= k && k <= i) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::d, index::b, index::c, index::a}, X);
    } else if (l <= k && k <= i && i <= j) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::d, index::c, index::a, index::b}, X);
    } else if (l <= k && k <= j && j <= i) {
        permute(Indices{index::a, index::b, index::c, index::d}, &Xperm, Indices{index::d, index::c, index::b, index::a}, X);
    }

    return Xperm;
}

void DLPNOCCSDT_Q::quadruples_sparsity(bool prescreening) {
    timer_on("Quadruples Sparsity");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();
    int npao = C_pao_->colspi(0);

    int MAX_WEAK_PAIRS = options_.get_int("QUADS_MAX_WEAK_PAIRS");

    if (prescreening) {
        int ijkl = 0;
        // Every quadruplet contains at least three strong pairs
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

                for (int l : lmopair_to_lmos_[ij]) {
                    int kl = i_j_to_ij_[k][l];
                    if (kl == -1) continue;
                    if (i > l || j > l || k > l) continue;
                    int il_weak = i_j_to_ij_weak_[i][l], jl_weak = i_j_to_ij_weak_[j][l], kl_weak = i_j_to_ij_weak_[k][l];
                    if (i == j && j == l || j == k && k == l || i == k && k == l) continue;

                    if (il_weak != -1) weak_pair_count += 1;
                    if (jl_weak != -1) weak_pair_count += 1;
                    if (kl_weak != -1) weak_pair_count += 1;

                    if (weak_pair_count > MAX_WEAK_PAIRS) continue;

                    ijkl_to_i_j_k_l_.push_back(std::make_tuple(i, j, k, l));

                    std::vector<int> ijkl_list = {i, j, k, l};
                    for (const auto &perm : quad_perms_long) {
                        auto &[i_idx, j_idx, k_idx, l_idx] = perm;
                        int ip = ijkl_list[i_idx], jp = ijkl_list[j_idx], kp = ijkl_list[k_idx], lp = ijkl_list[l_idx];
                        i_j_k_l_to_ijkl_[ip * std::pow(naocc, 3) + jp * std::pow(naocc, 2) + kp * naocc + lp] = ijkl;
                    } // end for
                    ++ijkl;
                } // end l
            } // end k
        } // end ij
    } else {
        std::unordered_map<int, int> i_j_k_l_to_ijkl_new;
        std::vector<std::tuple<int, int, int, int>> ijkl_to_i_j_k_l_new;
        std::vector<double> e_ijkl_new;

        double t_cut_quadruples_weak = options_.get_double("T_CUT_QUADS_WEAK");
        de_lccsdt_q_screened_ = 0.0;

        int ijkl_new = 0;
        for (int ijkl = 0; ijkl < ijkl_to_i_j_k_l_.size(); ++ijkl) {
            auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

            if (std::fabs(e_ijkl_[ijkl]) >= t_cut_quadruples_weak) {
                ijkl_to_i_j_k_l_new.push_back(std::make_tuple(i, j, k, l));
                e_ijkl_new.push_back(e_ijkl_[ijkl]);
                std::vector<int> ijkl_list = {i, j, k, l};
                for (const auto &perm : quad_perms_long) {
                    auto &[i_idx, j_idx, k_idx, l_idx] = perm;
                    int ip = ijkl_list[i_idx], jp = ijkl_list[j_idx], kp = ijkl_list[k_idx], lp = ijkl_list[l_idx];
                    i_j_k_l_to_ijkl_new[ip * std::pow(naocc, 3) + jp * std::pow(naocc, 2) + kp * naocc + lp] = ijkl_new;
                } // end for
                ++ijkl_new;
            } else {
                de_lccsdt_q_screened_ += e_ijkl_[ijkl];
            }
        }
        i_j_k_l_to_ijkl_ = i_j_k_l_to_ijkl_new;
        ijkl_to_i_j_k_l_ = ijkl_to_i_j_k_l_new;
        e_ijkl_ = e_ijkl_new;
    }

    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    int natom = molecule_->natom();
    int nbf = basisset_->nbf();

    qno_scale_.clear();
    qno_scale_.resize(n_lmo_quadruplets, 1.0);

    // => Local density fitting domains <= //

    SparseMap lmo_to_ribfs(naocc);
    SparseMap lmo_to_riatoms(naocc);

    double t_cut_mkn_quadruples = (prescreening) ? options_.get_double("T_CUT_MKN_QUADS_PRE") : options_.get_double("T_CUT_MKN_QUADS");

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
                mkn_pop[centerU] += p_uv * ((p_uu) / (p_uu + p_vv));
                mkn_pop[centerV] += p_uv * ((p_vv) / (p_uu + p_vv));
            }
        }

        // if non-zero mulliken pop on atom, include atom in the LMO's fitting domain
        for (size_t a = 0; a < natom; a++) {
            if (fabs(mkn_pop[a]) > t_cut_mkn_quadruples) {
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

    double t_cut_do_quadruples = (prescreening) ? options_.get_double("T_CUT_DO_QUADS_PRE") : options_.get_double("T_CUT_DO_QUADS");

    for (size_t i = 0; i < naocc; ++i) {
        // PAO domains determined by differential overlap integral
        std::vector<int> lmo_to_paos_temp;
        for (size_t u = 0; u < nbf; ++u) {
            if (fabs(DOI_iu_->get(i, u)) > t_cut_do_quadruples) {
                lmo_to_paos_temp.push_back(u);
            }
        }

        // if any PAO on an atom is in the list, we take all of the PAOs on that atom
        lmo_to_paos[i] = contract_lists(lmo_to_paos_temp, atom_to_bf_);
    }

    if (!prescreening) {
        lmoquadruplet_to_ribfs_.clear();
        lmoquadruplet_to_lmos_.clear();
        lmoquadruplet_to_paos_.clear();
    }

    lmoquadruplet_to_ribfs_.resize(n_lmo_quadruplets);
    lmoquadruplet_to_lmos_.resize(n_lmo_quadruplets);
    lmoquadruplet_to_paos_.resize(n_lmo_quadruplets);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];
        int ij = i_j_to_ij_[i][j], ik = i_j_to_ij_[i][k], il = i_j_to_ij_[i][l], 
            jk = i_j_to_ij_[j][k], jl = i_j_to_ij_[j][l], kl = i_j_to_ij_[k][l];

        lmoquadruplet_to_ribfs_[ijkl] = merge_lists(lmo_to_ribfs[i], merge_lists(lmo_to_ribfs[j], merge_lists(lmo_to_ribfs[k], lmo_to_ribfs[l])));
        for (int m = 0; m < naocc; ++m) {
            int im = i_j_to_ij_[i][m], jm = i_j_to_ij_[j][m], km = i_j_to_ij_[k][m], lm = i_j_to_ij_[l][m];
            if (im != -1 && jm != -1 && km != -1 && lm != -1) lmoquadruplet_to_lmos_[ijkl].push_back(m);
        }
        lmoquadruplet_to_paos_[ijkl] = merge_lists(lmo_to_paos[i], merge_lists(lmo_to_paos[j], merge_lists(lmo_to_paos[k], lmo_to_paos[l])));
    }

    // => Make Full Lists <= //
    if (!prescreening) {
        ijkl_to_i_j_k_l_full_.clear();
    }

    for (int ij = 0; ij < n_lmo_pairs; ij++) {
        auto &[i, j] = ij_to_i_j_[ij];
        for (int k : lmopair_to_lmos_[ij]) {
            if (i == j && j == k) continue;

            for (int l : lmopair_to_lmos_[ij]) {
                int kl = i_j_to_ij_[k][l];
                if (kl == -1) continue;
                if (i == j && j == l || j == k && k == l || i == k && k == l) continue;

                // Is any permutation of this a quadruplet?
                int ijkl_dense = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + l;
                if (i_j_k_l_to_ijkl_.count(ijkl_dense)) {
                    int ijkl = i_j_k_l_to_ijkl_[ijkl_dense];
                    ijkl_to_i_j_k_l_full_.push_back(std::make_tuple(i, j, k, l));
                }
            } // end l
        } // end k
    } // ij

    timer_off("Quadruples Sparsity");
}

void DLPNOCCSDT_Q::sort_quadruplets(double e_total) {
    timer_on("Sort Quadruplets");

    outfile->Printf("  ==> Sorting Quadruplets <== \n\n");

    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    std::vector<std::pair<int, double>> ijkl_e_pairs(n_lmo_quadruplets);

#pragma omp parallel for
    for (int ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        ijkl_e_pairs[ijkl] = std::make_pair(ijkl, e_ijkl_[ijkl]);
    }

    std::sort(ijkl_e_pairs.begin(), ijkl_e_pairs.end(), [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return (std::fabs(a.second) > std::fabs(b.second));
    });

    double e_curr = 0.0;
    double strong_scale = options_.get_double("T_CUT_QNO_STRONG_SCALE");
    double weak_scale = options_.get_double("T_CUT_QNO_WEAK_SCALE");
    is_strong_quadruplet_.resize(n_lmo_quadruplets, false);
    qno_scale_.clear();
    qno_scale_.resize(n_lmo_quadruplets, weak_scale);

    int strong_count = 0;
    for (int idx = 0; idx < n_lmo_quadruplets; ++idx) {
        is_strong_quadruplet_[ijkl_e_pairs[idx].first] = true;
        qno_scale_[ijkl_e_pairs[idx].first] = strong_scale;
        e_curr += ijkl_e_pairs[idx].second;
        ++strong_count;
        if (e_curr / e_total > 0.9) break;
    }

    outfile->Printf("    Number of Strong Quadruplets: %6d, Total Quadruplets: %6d, Ratio: %.4f\n\n", strong_count, n_lmo_quadruplets, 
                            (double) strong_count / n_lmo_quadruplets);

    timer_off("Sort Quadruplets");
}

void DLPNOCCSDT_Q::qno_transform(double t_cut_qno) {
    timer_on("QNO transform");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();
    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    int min_qnos = options_.get_int("MIN_QNOS");

    X_qno_.clear();
    e_qno_.clear();
    n_qno_.clear();

    X_qno_.resize(n_lmo_quadruplets);
    e_qno_.resize(n_lmo_quadruplets);
    n_qno_.resize(n_lmo_quadruplets);

    ijkl_scale_.resize(n_lmo_quadruplets, 1.0);

    std::vector<SharedMatrix> D_ij_list(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        SharedMatrix D_ij = linalg::doublet(Tt_iajb_[ij], T_iajb_[ij], false, true);
        D_ij->add(linalg::doublet(Tt_iajb_[ij], T_iajb_[ij], true, false));
        if (i == j) D_ij->scale(0.5);

        D_ij_list[ij] = D_ij;
    }

    std::vector<double> occ_qno(n_lmo_quadruplets, 0.0);
    std::vector<double> trace_qno(n_lmo_quadruplets, 0.0);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];
        int ij = i_j_to_ij_[i][j], ik = i_j_to_ij_[i][k], il = i_j_to_ij_[i][l], 
            jk = i_j_to_ij_[j][k], jl = i_j_to_ij_[j][l], kl = i_j_to_ij_[k][l];

        // number of PAOs in the triplet domain (before removing linear dependencies)
        int npao_ijkl = lmoquadruplet_to_paos_[ijkl].size();

        // number of auxiliary basis in the domain
        int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();

        //                                          //
        // ==> Canonicalize PAOs of quadruplet ijkl <== //
        //                                          //

        auto S_pao_ijkl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmoquadruplet_to_paos_[ijkl]);
        auto F_pao_ijkl = submatrix_rows_and_cols(*F_pao_, lmoquadruplet_to_paos_[ijkl], lmoquadruplet_to_paos_[ijkl]);

        SharedMatrix X_pao_ijkl;
        SharedVector e_pao_ijkl;
        std::tie(X_pao_ijkl, e_pao_ijkl) = orthocanonicalizer(S_pao_ijkl, F_pao_ijkl);

        F_pao_ijkl = linalg::triplet(X_pao_ijkl, F_pao_ijkl, X_pao_ijkl, true, false, false);

        // number of PAOs in the domain after removing linear dependencies
        int npao_can_ijkl = X_pao_ijkl->colspi(0);

        // S_ijkl partially transformed overlap matrix
        std::vector<int> pair_ext_domain = merge_lists(lmo_to_paos_[i], merge_lists(lmo_to_paos_[j], merge_lists(lmo_to_paos_[k], lmo_to_paos_[l])));
        auto S_ijkl = submatrix_rows_and_cols(*S_pao_, pair_ext_domain, lmoquadruplet_to_paos_[ijkl]);
        S_ijkl = linalg::doublet(S_ijkl, X_pao_ijkl, false, false);

        //                                           //
        // ==> Canonical PAOs  to Canonical QNOs <== //
        //                                           //

        size_t nvir_ijkl = F_pao_ijkl->rowspi(0);

        // Project pair densities into combined PAO space of ijkl
        std::vector<int> ij_index = index_list(pair_ext_domain, lmopair_to_paos_[ij]);
        auto S_ij = linalg::doublet(X_pno_[ij], submatrix_rows(*S_ijkl, ij_index), true, false);
        auto D_ij = linalg::triplet(S_ij, D_ij_list[ij], S_ij, true, false, false);

        std::vector<int> ik_index = index_list(pair_ext_domain, lmopair_to_paos_[ik]);
        auto S_ik = linalg::doublet(X_pno_[ik], submatrix_rows(*S_ijkl, ik_index), true, false);
        auto D_ik = linalg::triplet(S_ik, D_ij_list[ik], S_ik, true, false, false);

        std::vector<int> il_index = index_list(pair_ext_domain, lmopair_to_paos_[il]);
        auto S_il = linalg::doublet(X_pno_[il], submatrix_rows(*S_ijkl, il_index), true, false);
        auto D_il = linalg::triplet(S_il, D_ij_list[il], S_il, true, false, false);

        std::vector<int> jk_index = index_list(pair_ext_domain, lmopair_to_paos_[jk]);
        auto S_jk = linalg::doublet(X_pno_[jk], submatrix_rows(*S_ijkl, jk_index), true, false);
        auto D_jk = linalg::triplet(S_jk, D_ij_list[jk], S_jk, true, false, false);

        std::vector<int> jl_index = index_list(pair_ext_domain, lmopair_to_paos_[jl]);
        auto S_jl = linalg::doublet(X_pno_[jl], submatrix_rows(*S_ijkl, jl_index), true, false);
        auto D_jl = linalg::triplet(S_jl, D_ij_list[jl], S_jl, true, false, false);

        std::vector<int> kl_index = index_list(pair_ext_domain, lmopair_to_paos_[kl]);
        auto S_kl = linalg::doublet(X_pno_[kl], submatrix_rows(*S_ijkl, kl_index), true, false);
        auto D_kl = linalg::triplet(S_kl, D_ij_list[kl], S_kl, true, false, false);

        // Construct quadruplet density from pair densities
        auto D_ijkl = D_ij->clone();
        D_ijkl->add(D_ik);
        D_ijkl->add(D_il);
        D_ijkl->add(D_jk);
        D_ijkl->add(D_jl);
        D_ijkl->add(D_kl);
        D_ijkl->scale(1.0 / 6.0);

        // Diagonalization of quadruplet density gives QNOs (in basis of LMO's virtual domain)
        // as well as QNO occ numbers
        auto X_qno_ijkl = std::make_shared<Matrix>("eigenvectors", nvir_ijkl, nvir_ijkl);
        Vector qno_occ("eigenvalues", nvir_ijkl);
        D_ijkl->diagonalize(*X_qno_ijkl, qno_occ, descending);

        // Compute trace sum
        double occ_total = 0.0;
        for (size_t a = 0; a < nvir_ijkl; ++a) {
            occ_total += qno_occ.get(a);
        }

        double qno_scale = qno_scale_[ijkl];

        int nvir_ijkl_final = 0;
        double occ_curr = 0.0;
        for (size_t a = 0; a < nvir_ijkl; ++a) {
            if (fabs(qno_occ.get(a)) >= qno_scale * t_cut_qno || a < min_qnos) {
                occ_curr += qno_occ.get(a);
                nvir_ijkl_final++;
            }
        }

        nvir_ijkl_final = std::max(1, nvir_ijkl_final);

        Dimension zero(1);
        Dimension dim_final(1);
        dim_final.fill(nvir_ijkl_final);

        // This transformation gives orbitals that are orthonormal but not canonical
        X_qno_ijkl = X_qno_ijkl->get_block({zero, X_qno_ijkl->rowspi()}, {zero, dim_final});
        qno_occ = qno_occ.get_block({zero, dim_final});

        SharedMatrix qno_canon;
        SharedVector e_qno_ijkl;
        std::tie(qno_canon, e_qno_ijkl) = canonicalizer(X_qno_ijkl, F_pao_ijkl);

        X_qno_ijkl = linalg::doublet(X_qno_ijkl, qno_canon, false, false);
        X_qno_ijkl = linalg::doublet(X_pao_ijkl, X_qno_ijkl, false, false);

        X_qno_[ijkl] = X_qno_ijkl;
        e_qno_[ijkl] = e_qno_ijkl;
        n_qno_[ijkl] = X_qno_ijkl->colspi(0);
        occ_qno[ijkl] = qno_occ.get(n_qno_[ijkl] - 1);
        trace_qno[ijkl] = occ_curr / occ_total;
    }

    int qno_count_total = 0, qno_count_min = C_pao_->colspi(0), qno_count_max = 0;
    double occ_number_total = 0.0, occ_number_min = 2.0, occ_number_max = 0.0;
    double trace_total = 0.0, trace_min = 1.0, trace_max = 0.0;
    for (int ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        qno_count_total += n_qno_[ijkl];
        qno_count_min = std::min(qno_count_min, n_qno_[ijkl]);
        qno_count_max = std::max(qno_count_max, n_qno_[ijkl]);
        occ_number_total += occ_qno[ijkl];
        occ_number_min = std::min(occ_number_min, occ_qno[ijkl]);
        occ_number_max = std::max(occ_number_max, occ_qno[ijkl]);
        trace_total += trace_qno[ijkl];
        trace_min = std::min(trace_min, trace_qno[ijkl]);
        trace_max = std::max(trace_max, trace_qno[ijkl]);
    }

    // From ChatGPT: N choose 4 (all distinct) + N choose 2 (two values appear twice) 
    // + 3 * N choose 3 (one value appears twice, each other value appears once)
    size_t n_total_possible = (naocc) * (naocc - 1) * (naocc - 2) * (naocc - 3) / 24 + (naocc) * (naocc - 1) / 2 
        + 3 * (naocc) * (naocc - 1) * (naocc - 2) / 6;

    outfile->Printf("  \n");
    outfile->Printf("    Number of (Unique) Local MO quadruplets: %d\n", n_lmo_quadruplets);
    outfile->Printf("    Max Number of Possible (Unique) LMO Quadruplets: %d (Ratio: %.4f)\n", n_total_possible,
                    (double)n_lmo_quadruplets / n_total_possible);
    outfile->Printf("    Natural Orbitals per Local MO quadruplet:\n");
    outfile->Printf("      Avg: %3d NOs \n", qno_count_total / n_lmo_quadruplets);
    outfile->Printf("      Min: %3d NOs \n", qno_count_min);
    outfile->Printf("      Max: %3d NOs \n", qno_count_max);
    outfile->Printf("      Avg Occ Number Tol: %.3e \n", occ_number_total / n_lmo_quadruplets);
    outfile->Printf("      Min Occ Number Tol: %.3e \n", occ_number_min);
    outfile->Printf("      Max Occ Number Tol: %.3e \n", occ_number_max);
    outfile->Printf("      Avg Trace Sum: %.6f \n", trace_total / n_lmo_quadruplets);
    outfile->Printf("      Min Trace Sum: %.6f \n", trace_min);
    outfile->Printf("      Max Trace Sum: %.6f \n", trace_max);
    outfile->Printf("  \n");

    // Sort list of quadruplets based on number of QNOs (for parallel efficiency)
    std::vector<std::pair<int, int>> ijkl_qnos(n_lmo_quadruplets);

#pragma omp parallel for
    for (int ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        ijkl_qnos[ijkl] = std::make_pair(ijkl, n_qno_[ijkl]);
    }
    
    std::sort(ijkl_qnos.begin(), ijkl_qnos.end(), [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return (a.second > b.second);
    });

    sorted_quadruplets_.resize(n_lmo_quadruplets);
#pragma omp parallel for
    for (int ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        sorted_quadruplets_[ijkl] = ijkl_qnos[ijkl].first;
    }

    timer_off("QNO transform");
}

void DLPNOCCSDT_Q::estimate_memory() {
    outfile->Printf("\n  ==> DLPNO-(Q) Memory Requirements <== \n\n");

    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    size_t K_iabe_memory = 0;
    size_t K_iajm_memory = 0;
    size_t K_iajb_memory = 0;
    size_t qno_total_memory = 0;

#pragma omp parallel for reduction(+ : K_iabe_memory, K_iajm_memory, K_iajb_memory, qno_total_memory)
    for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
        int ijkl = sorted_quadruplets_[ijkl_sorted];
        int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();

        K_iabe_memory += 4 * n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl];
        K_iajm_memory += 16 * n_qno_[ijkl] * nlmo_ijkl;
        K_iajb_memory += 16 * n_qno_[ijkl] * n_qno_[ijkl];
        qno_total_memory += n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl];
    }

    size_t total_memory = qij_memory_ + qia_memory_ + qab_memory_ 
        + K_iabe_memory + K_iajm_memory + 2 * K_iajb_memory + 2 * qno_total_memory;

    outfile->Printf("    (q | i j) integrals    : %.3f [GiB]\n", qij_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (q | i a) integrals    : %.3f [GiB]\n", qia_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (q | a b) integrals    : %.3f [GiB]\n", qab_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (i a | b e)            : %.3f [GiB]\n", K_iabe_memory * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (i a | j m)            : %.3f [GiB]\n", K_iajm_memory * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (i a | j b)            : %.3f [GiB]\n", 2 * K_iajb_memory * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    T_{ijkl}^{abcd}        : %.3f [GiB]\n", qno_total_memory * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    G_{ijkl}^{abcd}        : %.3f [GiB]\n", qno_total_memory * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    Total Memory Given     : %.3f [GiB]\n", memory_ * pow(2.0, -30));
    outfile->Printf("    Total Memory Required  : %.3f [GiB]\n\n", total_memory * pow(2.0, -30) * sizeof(double));

}

double DLPNOCCSDT_Q::compute_gamma_ijkl(bool store_amplitudes) {
    timer_on("gamma ijkl");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    auto einsum_indices = std::make_tuple(Indices{a, b, c, d}, Indices{a, b, d, c}, Indices{a, c, b, d}, Indices{a, c, d, b}, 
        Indices{a, d, b, c}, Indices{a, d, c, b}, Indices{b, a, c, d}, Indices{b, a, d, c}, Indices{b, c, a, d}, Indices{b, c, d, a}, 
        Indices{b, d, a, c}, Indices{b, d, c, a}, Indices{c, a, b, d}, Indices{c, a, d, b}, Indices{c, b, a, d}, Indices{c, b, d, a}, 
        Indices{c, d, a, b}, Indices{c, d, b, a}, Indices{d, a, b, c}, Indices{d, a, c, b}, Indices{d, b, a, c}, Indices{d, b, c, a}, 
        Indices{d, c, a, b}, Indices{d, c, b, a});

    // Clear space in intermediates
    K_iabe_list_.clear();
    K_iajm_list_.clear();
    K_iajb_list_.clear();
    U_iajb_list_.clear();
    e_ijkl_.clear();

    // Allocate space for intermediates
    K_iabe_list_.resize(n_lmo_quadruplets);
    K_iajm_list_.resize(n_lmo_quadruplets);
    K_iajb_list_.resize(n_lmo_quadruplets);
    U_iajb_list_.resize(n_lmo_quadruplets);
    e_ijkl_.resize(n_lmo_quadruplets);

    if (store_amplitudes) {
        gamma_ijkl_.clear();
        T_iajbkcld_.clear();

        gamma_ijkl_.resize(n_lmo_quadruplets);
        T_iajbkcld_.resize(n_lmo_quadruplets);
    }

    std::time_t time_start = std::time(nullptr);
    std::time_t time_lap = std::time(nullptr);

    double E_Q0 = 0.0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+ : E_Q0)
    for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
        int ijkl = sorted_quadruplets_[ijkl_sorted];
        auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];
        int ij = i_j_to_ij_[i][j], ik = i_j_to_ij_[i][k], il = i_j_to_ij_[i][l], 
            jk = i_j_to_ij_[j][k], jl = i_j_to_ij_[j][l], kl = i_j_to_ij_[k][l];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // => Step 1: Compute all necessary integrals

        // number of LMOs in the quadruplet domain
        const int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
        // number of PAOs in the quadruplet domain (before removing linear dependencies)
        const int npao_ijkl = lmoquadruplet_to_paos_[ijkl].size();
        // number of auxiliary functions in the quadruplet domain
        const int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        // number of quadruplet natural orbitals in quadruplet domain
        const int nqno_ijkl = n_qno_[ijkl];

        /// => Necessary integrals <= ///

        // q_io integrals
        auto q_io = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl);
        auto q_jo = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl);
        auto q_ko = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl);
        auto q_lo = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl);

        // q_iv integrals
        auto q_iv = std::make_shared<Matrix>(naux_ijkl, npao_ijkl);
        auto q_jv = std::make_shared<Matrix>(naux_ijkl, npao_ijkl);
        auto q_kv = std::make_shared<Matrix>(naux_ijkl, npao_ijkl);
        auto q_lv = std::make_shared<Matrix>(naux_ijkl, npao_ijkl);

        // more expensive integrals
        auto q_ov = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl * nqno_ijkl);
        auto q_vv = std::make_shared<Matrix>(naux_ijkl, nqno_ijkl * nqno_ijkl);

        for (int q_ijkl = 0; q_ijkl < naux_ijkl; q_ijkl++) {
            const int q = lmoquadruplet_to_ribfs_[ijkl][q_ijkl];
            const int centerq = ribasis_->function_to_center(q);

            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                (*q_io)(q_ijkl, m_ijkl) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][i], riatom_to_lmos_ext_dense_[centerq][m]);
                (*q_jo)(q_ijkl, m_ijkl) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][j], riatom_to_lmos_ext_dense_[centerq][m]);
                (*q_ko)(q_ijkl, m_ijkl) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][k], riatom_to_lmos_ext_dense_[centerq][m]);
                (*q_lo)(q_ijkl, m_ijkl) = (*qij_[q])(riatom_to_lmos_ext_dense_[centerq][l], riatom_to_lmos_ext_dense_[centerq][m]);
            }

            for (int u_ijkl = 0; u_ijkl < npao_ijkl; ++u_ijkl) {
                int u = lmoquadruplet_to_paos_[ijkl][u_ijkl];
                (*q_iv)(q_ijkl, u_ijkl) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][i], riatom_to_paos_ext_dense_[centerq][u]);
                (*q_jv)(q_ijkl, u_ijkl) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][j], riatom_to_paos_ext_dense_[centerq][u]);
                (*q_kv)(q_ijkl, u_ijkl) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][k], riatom_to_paos_ext_dense_[centerq][u]);
                (*q_lv)(q_ijkl, u_ijkl) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][l], riatom_to_paos_ext_dense_[centerq][u]);
            }

            // More expensive integrals
            auto q_ov_tmp = std::make_shared<Matrix>(nlmo_ijkl, npao_ijkl);

            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                for (int u_ijkl = 0; u_ijkl < npao_ijkl; ++u_ijkl) {
                    int u = lmoquadruplet_to_paos_[ijkl][u_ijkl];
                    (*q_ov_tmp)(m_ijkl, u_ijkl) = (*qia_[q])(riatom_to_lmos_ext_dense_[centerq][m], riatom_to_paos_ext_dense_[centerq][u]);
                }
            }
            q_ov_tmp = linalg::doublet(q_ov_tmp, X_qno_[ijkl], false, false);
            ::memcpy(&(*q_ov)(q_ijkl, 0), &(*q_ov_tmp)(0, 0), nlmo_ijkl * nqno_ijkl * sizeof(double));

            auto q_vv_tmp = std::make_shared<Matrix>(npao_ijkl, npao_ijkl);

            for (int u_ijkl = 0; u_ijkl < npao_ijkl; ++u_ijkl) {
                int u = lmoquadruplet_to_paos_[ijkl][u_ijkl];
                for (int v_ijkl = 0; v_ijkl < npao_ijkl; ++v_ijkl) {
                    int v = lmoquadruplet_to_paos_[ijkl][v_ijkl];
                    int uv_idx = riatom_to_pao_pairs_dense_[centerq][u][v];
                    if (uv_idx == -1) continue;
                    (*q_vv_tmp)(u_ijkl, v_ijkl) = (*qab_[q])(uv_idx, 0);
                } // end v_ijk
            } // end u_ijk
            q_vv_tmp = linalg::triplet(X_qno_[ijkl], q_vv_tmp, X_qno_[ijkl], true, false, false);
            ::memcpy(&(*q_vv)(q_ijkl, 0), &(*q_vv_tmp)(0, 0), nqno_ijkl * nqno_ijkl * sizeof(double));
        }

        auto A_solve = submatrix_rows_and_cols(*full_metric_, lmoquadruplet_to_ribfs_[ijkl], lmoquadruplet_to_ribfs_[ijkl]);
        A_solve->power(0.5, 1.0e-14);

        C_DGESV_wrapper(A_solve->clone(), q_io);
        C_DGESV_wrapper(A_solve->clone(), q_jo);
        C_DGESV_wrapper(A_solve->clone(), q_ko);
        C_DGESV_wrapper(A_solve->clone(), q_lo);
        
        q_iv = linalg::doublet(q_iv, X_qno_[ijkl]);
        q_jv = linalg::doublet(q_jv, X_qno_[ijkl]);
        q_kv = linalg::doublet(q_kv, X_qno_[ijkl]);
        q_lv = linalg::doublet(q_lv, X_qno_[ijkl]);

        C_DGESV_wrapper(A_solve->clone(), q_iv);
        C_DGESV_wrapper(A_solve->clone(), q_jv);
        C_DGESV_wrapper(A_solve->clone(), q_kv);
        C_DGESV_wrapper(A_solve->clone(), q_lv);

        C_DGESV_wrapper(A_solve->clone(), q_ov);
        C_DGESV_wrapper(A_solve->clone(), q_vv);

        Tensor<double, 2> q_io_ein("(Q_ijkl | m i)", naux_ijkl, nlmo_ijkl);
        Tensor<double, 2> q_jo_ein("(Q_ijkl | m j)", naux_ijkl, nlmo_ijkl);
        Tensor<double, 2> q_ko_ein("(Q_ijkl | m k)", naux_ijkl, nlmo_ijkl);
        Tensor<double, 2> q_lo_ein("(Q_ijkl | m l)", naux_ijkl, nlmo_ijkl);
        ::memcpy(q_io_ein.data(), q_io->get_pointer(), naux_ijkl * nlmo_ijkl * sizeof(double));
        ::memcpy(q_jo_ein.data(), q_jo->get_pointer(), naux_ijkl * nlmo_ijkl * sizeof(double));
        ::memcpy(q_ko_ein.data(), q_ko->get_pointer(), naux_ijkl * nlmo_ijkl * sizeof(double));
        ::memcpy(q_lo_ein.data(), q_lo->get_pointer(), naux_ijkl * nlmo_ijkl * sizeof(double));

        Tensor<double, 2> q_iv_ein("(Q_ijkl | i a)", naux_ijkl, nqno_ijkl);
        Tensor<double, 2> q_jv_ein("(Q_ijkl | j b)", naux_ijkl, nqno_ijkl);
        Tensor<double, 2> q_kv_ein("(Q_ijkl | k c)", naux_ijkl, nqno_ijkl);
        Tensor<double, 2> q_lv_ein("(Q_ijkl | l d)", naux_ijkl, nqno_ijkl);
        ::memcpy(q_iv_ein.data(), q_iv->get_pointer(), naux_ijkl * nqno_ijkl * sizeof(double));
        ::memcpy(q_jv_ein.data(), q_jv->get_pointer(), naux_ijkl * nqno_ijkl * sizeof(double));
        ::memcpy(q_kv_ein.data(), q_kv->get_pointer(), naux_ijkl * nqno_ijkl * sizeof(double));
        ::memcpy(q_lv_ein.data(), q_lv->get_pointer(), naux_ijkl * nqno_ijkl * sizeof(double));

        Tensor<double, 3> q_ov_ein("(Q_ijkl | m a)", naux_ijkl, nlmo_ijkl, nqno_ijkl);
        Tensor<double, 3> q_vv_ein("(Q_ijkl | a b)", naux_ijkl, nqno_ijkl, nqno_ijkl);
        ::memcpy(q_ov_ein.data(), q_ov->get_pointer(), naux_ijkl * nlmo_ijkl * nqno_ijkl * sizeof(double));
        ::memcpy(q_vv_ein.data(), q_vv->get_pointer(), naux_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));

        // List of intermediates
        std::vector<int> ijkl_list = {i, j, k, l};
        std::vector<Tensor<double, 2>> q_io_list = {q_io_ein, q_jo_ein, q_ko_ein, q_lo_ein};
        std::vector<Tensor<double, 2>> q_iv_list = {q_iv_ein, q_jv_ein, q_kv_ein, q_lv_ein};

        const int FOUR = ijkl_list.size();

        // Term 1 intermediates
        for (int idx = 0; idx < ijkl_list.size(); ++idx) {
            int i = ijkl_list[idx];
            K_iabe_list_[ijkl][idx] = Tensor<double, 3>("(ia|be)", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::a, index::b, index::e}, &K_iabe_list_[ijkl][idx], 1.0, Indices{index::Q, index::a}, q_iv_list[idx], Indices{index::Q, index::b, index::e}, q_vv_ein);
        }

        // Term 2 intermediates
        std::array<Tensor<double, 4>, 16> T_mkl_list;
        for (int i_idx = 0; i_idx < ijkl_list.size(); ++i_idx) {
            int i = ijkl_list[i_idx], k = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < ijkl_list.size(); ++j_idx) {
                int j = ijkl_list[j_idx], l = ijkl_list[j_idx];
                
                // Form K_iajm intermediate
                K_iajm_list_[ijkl][i_idx * FOUR + j_idx] = Tensor<double, 2>("(ia|jm)", nqno_ijkl, nlmo_ijkl);
                einsum(0.0, Indices{index::a, index::m}, &K_iajm_list_[ijkl][i_idx * FOUR + j_idx], 1.0, Indices{index::Q, index::a}, q_iv_list[i_idx], Indices{index::Q, index::m}, q_io_list[j_idx]);

                // Form T_mkl intermediate
                T_mkl_list[i_idx * FOUR + j_idx] = Tensor<double, 4>("T_mkl", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                T_mkl_list[i_idx * FOUR + j_idx].zero();

                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    int mkl_dense = m * naocc * naocc + k * naocc + l;
                    if (i_j_k_to_ijk_.count(mkl_dense)) {
                        int mkl = i_j_k_to_ijk_[mkl_dense];
                        auto S_ijkl_mkl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[mkl]);
                        S_ijkl_mkl = linalg::triplet(X_qno_[ijkl], S_ijkl_mkl, X_tno_[mkl], true, false, false);
                        auto T_mkl = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[mkl], m, k, l), 
                                                        S_ijkl_mkl, n_tno_[mkl], n_qno_[ijkl]);

                        ::memcpy(&T_mkl_list[i_idx * FOUR + j_idx](m_ijkl, 0, 0, 0), T_mkl.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                    } // end if
                } // end m_ijkl
            } // end j_idx
        } // end i_idx

        // Term 3 intermediates
        std::array<Tensor<double, 2>, 16> K_minj_list;
        for (int i_idx = 0; i_idx < ijkl_list.size(); ++i_idx) {
            for (int j_idx = 0; j_idx < ijkl_list.size(); ++j_idx) {

                K_minj_list[i_idx * FOUR + j_idx] = Tensor<double, 2>("(mi|nj)", nlmo_ijkl, nlmo_ijkl);
                einsum(0.0, Indices{index::m, index::n}, &K_minj_list[i_idx * FOUR + j_idx], 1.0, Indices{index::Q, index::m}, q_io_list[i_idx], Indices{index::Q, index::n}, q_io_list[j_idx]);
            }
        }

        std::array<Tensor<double, 3>, 4> T_mkac_list;
        for (int idx = 0; idx < ijkl_list.size(); ++idx) {
            int i = ijkl_list[idx];

            T_mkac_list[idx] = Tensor<double, 3>("T_mkac", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            T_mkac_list[idx].zero();

            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                int mi = i_j_to_ij_[m][i];

                auto S_ijkl_mi = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[mi]);
                S_ijkl_mi = linalg::triplet(X_qno_[ijkl], S_ijkl_mi, X_pno_[mi], true, false, false);

                auto T_mi = linalg::triplet(S_ijkl_mi, T_iajb_[mi], S_ijkl_mi, false, false, true);
                ::memcpy(&T_mkac_list[idx](m_ijkl, 0, 0), T_mi->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
            } // end m_ijkl
        } // end idx

        // Term 4 intermediates
        std::array<Tensor<double, 3>, 4> K_iame_list;
        for (int idx = 0; idx < ijkl_list.size(); ++idx) {
            int i = ijkl_list[idx];
            K_iame_list[idx] = Tensor<double, 3>("(ia|me)", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::a, index::m, index::e}, &K_iame_list[idx], 1.0, Indices{index::Q, index::a}, q_iv_list[idx], Indices{index::Q, index::m, index::e}, q_ov_ein);
        }

        std::array<Tensor<double, 2>, 16> T_ijab_list;
        for (int i_idx = 0; i_idx < ijkl_list.size(); ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < ijkl_list.size(); ++j_idx) {
                int j = ijkl_list[j_idx];
                int ij = i_j_to_ij_[i][j];
                // if (i_idx > j_idx) continue;

                T_ijab_list[i_idx * FOUR + j_idx] = Tensor<double, 2>("T_ijab", nqno_ijkl, nqno_ijkl);

                auto S_ijkl_ij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij, X_pno_[ij], true, false, false);
                auto T_ij = linalg::triplet(S_ijkl_ij, T_iajb_[ij], S_ijkl_ij, false, false, true);
                ::memcpy(T_ijab_list[i_idx * FOUR + j_idx].data(), T_ij->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
            } // end j_idx
        } // end i_idx

        // Term 5 intermediates
        std::array<Tensor<double, 3>, 4> K_mibe_list;
        for (int idx = 0; idx < ijkl_list.size(); ++idx) {
            int i = ijkl_list[idx];
            K_mibe_list[idx] = Tensor<double, 3>("(mi|be)", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::m, index::b, index::e}, &K_mibe_list[idx], 1.0, Indices{index::Q, index::m}, q_io_list[idx], Indices{index::Q, index::b, index::e}, q_vv_ein);
        }

        // Term 6 intermediates
        std::array<Tensor<double, 3>, 16> theta_Qab;
        for (int i_idx = 0; i_idx < ijkl_list.size(); ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < ijkl_list.size(); ++j_idx) {
                int j = ijkl_list[j_idx];
                int ij = i_j_to_ij_[i][j], ij_idx = i_idx * FOUR + j_idx;

                theta_Qab[ij_idx] = Tensor<double, 3>("theta_Qab", naux_ijkl, nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::Q, index::a, index::b}, &theta_Qab[ij_idx], 1.0, Indices{index::Q, index::a, index::e}, q_vv_ein, Indices{index::e, index::b}, T_ijab_list[ij_idx]);
            }
        }

        // Energy intermediates
        for (int i_idx = 0; i_idx < ijkl_list.size(); ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < ijkl_list.size(); ++j_idx) {
                int j = ijkl_list[j_idx];
                int ij = i_j_to_ij_[i][j], ij_idx = i_idx * FOUR + j_idx;

                // K_iajb
                K_iajb_list_[ijkl][i_idx * FOUR + j_idx] = Tensor<double, 2>("K_iajb", nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::a, index::b}, &K_iajb_list_[ijkl][i_idx * FOUR + j_idx], 1.0, Indices{index::Q, index::a}, q_iv_list[i_idx], Indices{index::Q, index::b}, q_iv_list[j_idx]);

                // U_iajb
                U_iajb_list_[ijkl][i_idx * FOUR + j_idx] = Tensor<double, 2>("U_iajb", nqno_ijkl, nqno_ijkl);
                auto S_ijkl_ij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij, X_pno_[ij], true, false, false);
                auto U_ij_psi = linalg::triplet(S_ijkl_ij, Tt_iajb_[ij], S_ijkl_ij, false, false, true);
                ::memcpy(U_iajb_list_[ijkl][i_idx * FOUR + j_idx].data(), U_ij_psi->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
            }
        }

        // => Form all possible gamma_ijkl's over unique (i, j, k, l) permutations in quadruplet ijkl
        Tensor<double, 4> gamma_ijkl("gamma_ijkl", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
        gamma_ijkl.zero();

        std::unordered_map<int, Tensor<double, 4>> gamma_ijkl_list;
        
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quad_perms_long[perm_idx];
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
            int ijkl_idx = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + l;

            // This reduces the operation count if there is any symmetry in the quadruplet (i.e iikl, ikil, etc.)
            if (!gamma_ijkl_list.count(ijkl_idx)) {
                const int ij_idx = i_idx * FOUR + j_idx, kl_idx = k_idx * FOUR + l_idx, kj_idx = k_idx * FOUR + j_idx;

                Tensor<double, 4> gamma_ijkl_perm("gamma_ijkl_perm", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                gamma_ijkl_perm.zero();

                Tensor<double, 4> gamma_ijkl_buff_a("gamma_ijkl_buff_a", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                Tensor<double, 4> gamma_ijkl_buff_b("gamma_ijkl_buff_b", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);

                // Term 1
                int jkl_dense = j * naocc * naocc + k * naocc + l;
                if (i_j_k_to_ijk_.count(jkl_dense)) {
                    int jkl = i_j_k_to_ijk_[jkl_dense];
                    auto S_ijkl_jkl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[jkl]);
                    S_ijkl_jkl = linalg::triplet(X_qno_[ijkl], S_ijkl_jkl, X_tno_[jkl], true, false, false);
                    auto T_jkl = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[jkl], j, k, l), 
                                                    S_ijkl_jkl, n_tno_[jkl], n_qno_[ijkl]);

                    einsum(0.0, Indices{index::a, index::b, index::c, index::d}, &gamma_ijkl_perm, 1.0, Indices{index::a, index::b, index::e}, K_iabe_list_[ijkl][i_idx], Indices{index::e, index::c, index::d}, T_jkl);
                }

                // Term 2
                einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &gamma_ijkl_perm, -1.0, Indices{index::a, index::m}, K_iajm_list_[ijkl][ij_idx], Indices{index::m, index::b, index::c, index::d}, T_mkl_list[kl_idx]);

                // Term 3
                Tensor<double, 3> gamma_3_temp("gamma_3_temp", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::n, index::a, index::c}, &gamma_3_temp, 1.0, Indices{index::m, index::n}, K_minj_list[ij_idx], Indices{index::m, index::a, index::c}, T_mkac_list[k_idx]);
                einsum(0.0, Indices{index::a, index::c, index::b, index::d}, &gamma_ijkl_buff_a, 1.0, Indices{index::n, index::a, index::c}, gamma_3_temp, Indices{index::n, index::b, index::d}, T_mkac_list[l_idx]);
                permute(Indices{index::a, index::b, index::c, index::d}, &gamma_ijkl_buff_b, Indices{index::a, index::c, index::b, index::d}, gamma_ijkl_buff_a);
                gamma_ijkl_perm += gamma_ijkl_buff_b;

                // Term 4
                Tensor<double, 3> gamma_4_temp("gamma_4_temp", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::a, index::m, index::b}, &gamma_4_temp, 1.0, Indices{index::a, index::m, index::e}, K_iame_list[i_idx], Indices{index::e, index::b}, T_ijab_list[kj_idx]);
                Tensor<double, 3> gamma_4_temp_b("gamma_4_temp_b", nqno_ijkl, nqno_ijkl, nlmo_ijkl);
                permute(Indices{index::a, index::b, index::m}, &gamma_4_temp_b, Indices{index::a, index::m, index::b}, gamma_4_temp);
                einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &gamma_ijkl_perm, -2.0, Indices{index::a, index::b, index::m}, gamma_4_temp_b, Indices{index::m, index::c, index::d}, T_mkac_list[l_idx]);

                // Term 5
                Tensor<double, 3> gamma_5_temp("gamma_5_temp", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::m, index::b, index::c}, &gamma_5_temp, 1.0, Indices{index::m, index::b, index::e}, K_mibe_list[i_idx], Indices{index::c, index::e}, T_ijab_list[kj_idx]);
                einsum(0.0, Indices{index::a, index::d, index::b, index::c}, &gamma_ijkl_buff_a, 1.0, Indices{index::m, index::a, index::d}, T_mkac_list[l_idx], Indices{index::m, index::b, index::c}, gamma_5_temp);
                permute(Indices{index::a, index::b, index::c, index::d}, &gamma_ijkl_buff_b, Indices{index::a, index::d, index::b, index::c}, gamma_ijkl_buff_a);
                gamma_ijkl_buff_b *= 2.0;
                gamma_ijkl_perm -= gamma_ijkl_buff_b;

                // Term 6
                einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &gamma_ijkl_perm, 1.0, Indices{index::Q, index::a, index::b}, theta_Qab[ij_idx], Indices{index::Q, index::c, index::d}, theta_Qab[kl_idx]);

                gamma_ijkl_list[ijkl_idx] = gamma_ijkl_perm;
                permute(Indices{index::a, index::b, index::c, index::d}, &gamma_ijkl_buff_a, std::get<perm_idx>(einsum_indices), gamma_ijkl_perm);
                gamma_ijkl += gamma_ijkl_buff_a;
            } else {
                Tensor<double, 4> gamma_ijkl_buff_a("gamma_ijkl_buff_a", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                permute(Indices{index::a, index::b, index::c, index::d}, &gamma_ijkl_buff_a, std::get<perm_idx>(einsum_indices), gamma_ijkl_list[ijkl_idx]);
                gamma_ijkl += gamma_ijkl_buff_a;
            }
        });

        gamma_ijkl *= 0.5;

        // Form T4 amplitudes from gamma_ijkl
        Tensor<double, 4> T_ijkl("T_ijkl", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
        T_ijkl.zero();

        for (int a_ijkl = 0; a_ijkl < nqno_ijkl; ++a_ijkl) {
            for (int b_ijkl = 0; b_ijkl < nqno_ijkl; ++b_ijkl) {
                for (int c_ijkl = 0; c_ijkl < nqno_ijkl; ++c_ijkl) {
                    for (int d_ijkl = 0; d_ijkl < nqno_ijkl; ++d_ijkl) {
                        (T_ijkl)(a_ijkl, b_ijkl, c_ijkl, d_ijkl) = (gamma_ijkl)(a_ijkl, b_ijkl, c_ijkl, d_ijkl) / 
                            ((*F_lmo_)(i,i) + (*F_lmo_)(j,j) + (*F_lmo_)(k,k) + (*F_lmo_)(l,l) - (*e_qno_[ijkl])(a_ijkl) 
                            - (*e_qno_[ijkl])(b_ijkl) - (*e_qno_[ijkl])(c_ijkl) - (*e_qno_[ijkl])(d_ijkl));
                    } // end d_ijkl
                } // end c_ijkl
            } // end b_ijkl
        } // end a_ijkl

        // Compute energy contribution
        double e_quad = compute_quadruplet_energy(ijkl, T_ijkl);
        e_ijkl_[ijkl] = e_quad;
        E_Q0 += e_quad;

        if (store_amplitudes) {
            gamma_ijkl_[ijkl] = gamma_ijkl;
            T_iajbkcld_[ijkl] = T_ijkl;
        } else {
            einsums::for_sequence<4UL>([&](auto idx) {
                K_iabe_list_[ijkl][idx] = Tensor<double, 3>("null", 0, 0, 0);
            });
            einsums::for_sequence<16UL>([&](auto idx) {
                K_iajm_list_[ijkl][idx] = Tensor<double, 2>("null", 0, 0);
                K_iajb_list_[ijkl][idx] = Tensor<double, 2>("null", 0, 0);
                U_iajb_list_[ijkl][idx] = Tensor<double, 2>("null", 0, 0);
            });
        }

        if (thread == 0) {
            std::time_t time_curr = std::time(nullptr);
            int time_elapsed = (int) time_curr - (int) time_lap;
            if (time_elapsed > 60) {
                outfile->Printf("  Time Elapsed from last checkpoint %4d (s), Progress %2d %%, Amplitudes for (%6d / %6d) Quadruplets Computed\n", time_elapsed, 
                                    (100 * ijkl_sorted) / n_lmo_quadruplets, ijkl_sorted, n_lmo_quadruplets);
                time_lap = std::time(nullptr);
            }
        } // end thread
    } // end [Q0]

    std::time_t time_stop = std::time(nullptr);
    int time_elapsed = (int) time_stop - (int) time_start;
    outfile->Printf("    Computation of T4 amplitudes complete!!! Time Elapsed: %4d seconds\n\n", time_elapsed);

    timer_off("gamma ijkl");

    return E_Q0;
}

double DLPNOCCSDT_Q::compute_quadruplet_energy(int ijkl, const Tensor<double, 4>& T4) {

    int naocc = i_j_to_ij_.size();

    auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];
    std::vector<int> ijkl_list = {i, j, k, l};
    const int FOUR = ijkl_list.size();

    // number of LMOs in the quadruplet domain
    const int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
    // number of PAOs in the quadruplet domain (before removing linear dependencies)
    const int npao_ijkl = lmoquadruplet_to_paos_[ijkl].size();
    // number of auxiliary functions in the quadruplet domain
    const int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
    // number of quadruplet natural orbitals in quadruplet domain
    const int nqno_ijkl = n_qno_[ijkl];

    double quadruplet_energy = 0.0;
    std::unordered_map<int, double> e_perm_energy;

    std::array<Tensor<double, 4>, 16> T_ijm_list;
    for (int i_idx = 0; i_idx < ijkl_list.size(); ++i_idx) {
        int i = ijkl_list[i_idx];
        for (int j_idx = 0; j_idx < ijkl_list.size(); ++j_idx) {
            int j = ijkl_list[j_idx];

            T_ijm_list[i_idx * FOUR + j_idx] = Tensor<double, 4>("T_ijm", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            T_ijm_list[i_idx * FOUR + j_idx].zero();

            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                int ijm_dense = i * naocc * naocc + j * naocc + m;
                if (i_j_k_to_ijk_.count(ijm_dense)) {
                    int ijm = i_j_k_to_ijk_[ijm_dense];
                    auto S_ijkl_ijm = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[ijm]);
                    S_ijkl_ijm = linalg::triplet(X_qno_[ijkl], S_ijkl_ijm, X_tno_[ijm], true, false, false);
                    auto T_ijm = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[ijm], i, j, m), 
                                                    S_ijkl_ijm, n_tno_[ijm], n_qno_[ijkl]);

                    ::memcpy(&T_ijm_list[i_idx * FOUR + j_idx](m_ijkl, 0, 0, 0), T_ijm.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                } // end if
            } // end m_ijkl
        } // end j_idx
    } // end i_idx
        
    einsums::for_sequence<24UL>([&](auto perm_idx) {
        auto &[i_idx, j_idx, k_idx, l_idx] = quad_perms_long[perm_idx];
        int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
        int ijkl_idx = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + l;

        if (!e_perm_energy.count(ijkl_idx)) {
            // Set up e_perm_energy
            e_perm_energy[ijkl_idx] = 0.0;

            // Get quadruples amplitude
            Tensor<double, 4> T_ijkl = quadruples_permuter(T4, i, j, k, l);

            // [Q] intermediates
            // u_{kl}^{ab}K_{ij}_{cd} - 2u_{kl}^{bd}L_{ij}^{ac} + u_{kl}^{cd}L_{ij}^{ab}
            Tensor<double, 2> U_kl = U_iajb_list_[ijkl][4 * k_idx + l_idx];
            Tensor<double, 2> K_ij = K_iajb_list_[ijkl][4 * i_idx + j_idx];
            Tensor<double, 2> L_ij = K_iajb_list_[ijkl][4 * i_idx + j_idx];
            L_ij *= 2.0;
            L_ij -= K_iajb_list_[ijkl][4 * j_idx + i_idx];

            for (int a_ijkl = 0; a_ijkl < nqno_ijkl; ++a_ijkl) {
                for (int b_ijkl = 0; b_ijkl < nqno_ijkl; ++b_ijkl) {
                    for (int c_ijkl = 0; c_ijkl < nqno_ijkl; ++c_ijkl) {
                        for (int d_ijkl = 0; d_ijkl < nqno_ijkl; ++d_ijkl) {
                            e_perm_energy[ijkl_idx] += T_ijkl(a_ijkl, b_ijkl, c_ijkl, d_ijkl) * (U_kl(a_ijkl, b_ijkl) * K_ij(c_ijkl, d_ijkl)
                                - 2.0 * U_kl(b_ijkl, d_ijkl) * L_ij(a_ijkl, c_ijkl) + U_kl(c_ijkl, d_ijkl) * L_ij(a_ijkl, b_ijkl));
                        } // end d_ijkl
                    } // end c_ijkl
                } // end b_ijkl
            } // end a_ijkl

            // Make a buffer
            Tensor<double, 4> T_buffer("T_buffer", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);

            // bar{t}_{ijkl}^{abcd} = -2t_{ijkl}^{abcd} - t_{ijkl}^{cdab} + t_{ijkl}^{bacd}
            Tensor<double, 4> T_bar = T_ijkl;
            T_bar *= -2.0;
            permute(Indices{index::a, index::b, index::c, index::d}, &T_buffer, Indices{index::c, index::d, index::a, index::b}, T_ijkl);
            T_bar -= T_buffer;
            permute(Indices{index::a, index::b, index::c, index::d}, &T_buffer, Indices{index::b, index::a, index::c, index::d}, T_ijkl);
            T_bar += T_buffer;

            // tilde{t}_{ijkl}^{abcd} = (1 + P_{kl}^{cd})[2t_{ijkl}^{dbac} - t_{ijkl}^{bdac}] =>
            // [2t_{ijkl}^{dbac} - t_{ijkl}^{bdac} + 2t_{ijlk}^{cbad} - t_{ijlk}^{bcad}] =>
            // 2t_{ijkl}^{dbac} - t_{ijkl}^{bdac} + 2t_{ijkl}^{cbda} - t_{ijkl}^{bcda}
            permute(Indices{index::a, index::b, index::c, index::d}, &T_buffer, Indices{index::d, index::b, index::a, index::c}, T_ijkl);
            Tensor<double, 4> T_tilde = T_buffer;
            permute(Indices{index::a, index::b, index::c, index::d}, &T_buffer, Indices{index::c, index::b, index::d, index::a}, T_ijkl);
            T_tilde += T_buffer;
            T_tilde *= 2.0;
            permute(Indices{index::a, index::b, index::c, index::d}, &T_buffer, Indices{index::b, index::d, index::a, index::c}, T_ijkl);
            T_tilde -= T_buffer;
            permute(Indices{index::a, index::b, index::c, index::d}, &T_buffer, Indices{index::b, index::c, index::d, index::a}, T_ijkl);
            T_tilde -= T_buffer;

            // => alpha and beta contributions <= //

            // => 2 - P_{cd} contributions

            Tensor<double, 4> T_bar_dc("T_bar_dc", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::a, index::b, index::d, index::c}, &T_bar_dc, Indices{index::a, index::b, index::c, index::d}, T_bar);
            Tensor<double, 4> T_tilde_dc("T_tilde_dc", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::a, index::b, index::d, index::c}, &T_tilde_dc, Indices{index::a, index::b, index::c, index::d}, T_tilde);

            Tensor<double, 2> alpha_ijm_buffer("alpha_ijm_buffer", nlmo_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::m, index::d}, &alpha_ijm_buffer, 2.0, Indices{index::m, index::a, index::b, index::c}, T_ijm_list[i_idx * FOUR + j_idx],
                    Indices{index::a, index::b, index::c, index::d}, T_bar);
            einsum(1.0, Indices{index::m, index::d}, &alpha_ijm_buffer, -1.0, Indices{index::m, index::a, index::b, index::c}, T_ijm_list[i_idx * FOUR + j_idx],
                    Indices{index::a, index::b, index::c, index::d}, T_bar_dc);

            Tensor<double, 2> beta_ijm_buffer("beta_ijm_buffer", nlmo_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::m, index::d}, &beta_ijm_buffer, 2.0, Indices{index::m, index::a, index::b, index::c}, T_ijm_list[i_idx * FOUR + j_idx],
                    Indices{index::a, index::b, index::c, index::d}, T_tilde);
            einsum(1.0, Indices{index::m, index::d}, &beta_ijm_buffer, -1.0, Indices{index::m, index::a, index::b, index::c}, T_ijm_list[i_idx * FOUR + j_idx],
                    Indices{index::a, index::b, index::c, index::d}, T_tilde_dc);

            Tensor<double, 2> K_lk_T("K_lk_T", nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::d}, &K_lk_T, Indices{index::d, index::m}, K_iajm_list_[ijkl][l_idx * 4 + k_idx]);
            Tensor<double, 2> K_kl_T("K_kl_T", nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::d}, &K_kl_T, Indices{index::d, index::m}, K_iajm_list_[ijkl][k_idx * 4 + l_idx]);

            e_perm_energy[ijkl_idx] += 2.0 * (linear_algebra::dot(alpha_ijm_buffer, K_lk_T) + linear_algebra::dot(beta_ijm_buffer, K_kl_T));

            // 2 - P_{kl} contributions
            int k_ijkl = std::find(lmoquadruplet_to_lmos_[ijkl].begin(), lmoquadruplet_to_lmos_[ijkl].end(), k) - lmoquadruplet_to_lmos_[ijkl].begin();
            Tensor<double, 3> T_ijk = T_ijm_list[i_idx * FOUR + j_idx](k_ijkl, All, All, All);
            int l_ijkl = std::find(lmoquadruplet_to_lmos_[ijkl].begin(), lmoquadruplet_to_lmos_[ijkl].end(), l) - lmoquadruplet_to_lmos_[ijkl].begin();
            Tensor<double, 3> T_ijl = T_ijm_list[i_idx * FOUR + j_idx](l_ijkl, All, All, All);

            Tensor<double, 3> T_ijk_bar("T_ijk_bar", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::c, index::d}, &T_ijk_bar, 1.0, Indices{index::a, index::b, index::c, index::d}, T_bar,
                    Indices{index::a, index::b, index::e}, T_ijk);
            Tensor<double, 3> T_ijl_bar("T_ijl_bar", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::c, index::d}, &T_ijl_bar, 1.0, Indices{index::a, index::b, index::c, index::d}, T_bar,
                    Indices{index::a, index::b, index::e}, T_ijl);
            Tensor<double, 3> T_ijk_tilde("T_ijk_tilde", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::c, index::d}, &T_ijk_tilde, 1.0, Indices{index::a, index::b, index::c, index::d}, T_tilde,
                    Indices{index::a, index::b, index::e}, T_ijk);
            Tensor<double, 3> T_ijl_tilde("T_ijl_tilde", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::c, index::d}, &T_ijl_tilde, 1.0, Indices{index::a, index::b, index::c, index::d}, T_tilde,
                    Indices{index::a, index::b, index::e}, T_ijl);

            Tensor<double, 3> g_kdce_T("g_kdce_T", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::e, index::c, index::d}, &g_kdce_T, Indices{index::d, index::c, index::e}, K_iabe_list_[ijkl][k_idx]);
            Tensor<double, 3> g_ldce_T("g_ldce_T", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::e, index::c, index::d}, &g_ldce_T, Indices{index::d, index::c, index::e}, K_iabe_list_[ijkl][l_idx]);

            e_perm_energy[ijkl_idx] += -4.0 * linear_algebra::dot(T_ijk_bar, g_ldce_T) + 2.0 * linear_algebra::dot(T_ijl_bar, g_kdce_T);
            e_perm_energy[ijkl_idx] += -4.0 * linear_algebra::dot(T_ijl_tilde, g_kdce_T) + 2.0 * linear_algebra::dot(T_ijk_tilde, g_ldce_T);

            quadruplet_energy += e_perm_energy[ijkl_idx];
        }
    });
        
    return quadruplet_energy;
}

double DLPNOCCSDT_Q::lccsdt_q_iterations() {
    timer_on("LCCSDT(Q) Iterations");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    outfile->Printf("\n  ==> Local CCSDT(Q) <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                         Corr. Energy    Delta E     Max R     Time (s)\n");

    int iteration = 1, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, r_curr = 0.0;
    bool e_converged = false, r_converged = false;

    double F_CUT = options_.get_double("F_CUT_Q");
    double T_CUT_ITER = options_.get_double("T_CUT_ITER_Q");

    std::vector<double> e_ijkl_old(n_lmo_quadruplets, 0.0);

    while (!(e_converged && r_converged)) {
        // RMS of residual per single LMO, for assesing convergence
        std::vector<double> R_iajbkcld_rms(n_lmo_quadruplets, 0.0);

        std::time_t time_start = std::time(nullptr);

#pragma omp parallel for schedule(dynamic, 1)
        for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
            int ijkl = sorted_quadruplets_[ijkl_sorted];
            auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

            int nqno_ijkl = n_qno_[ijkl];
            if (std::fabs(e_ijkl_[ijkl] - e_ijkl_old[ijkl]) < std::fabs(e_ijkl_old[ijkl] * T_CUT_ITER)) continue;

            // S integrals
            std::vector<int> quadruplet_ext_domain;
            for (int m = 0; m < naocc; ++m) {
                int ijkm_dense = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + m;
                int ijml_dense = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + m * naocc + l;
                int imkl_dense = i * std::pow(naocc, 3) + m * std::pow(naocc, 2) + k * naocc + l;
                int mjkl_dense = m * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + l;

                if (l != m && i_j_k_l_to_ijkl_.count(ijkm_dense) && std::fabs((*F_lmo_)(l, m)) >= F_CUT) {
                    int ijkm = i_j_k_l_to_ijkl_[ijkm_dense];
                    quadruplet_ext_domain = merge_lists(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijkm]);
                }

                if (k != m && i_j_k_l_to_ijkl_.count(ijml_dense) && std::fabs((*F_lmo_)(k, m)) >= F_CUT) {
                    int ijml = i_j_k_l_to_ijkl_[ijml_dense];
                    quadruplet_ext_domain = merge_lists(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijml]);
                }

                if (j != m && i_j_k_l_to_ijkl_.count(imkl_dense) && std::fabs((*F_lmo_)(j, m)) >= F_CUT) {
                    int imkl = i_j_k_l_to_ijkl_[imkl_dense];
                    quadruplet_ext_domain = merge_lists(quadruplet_ext_domain, lmoquadruplet_to_paos_[imkl]);
                }

                if (i != m && i_j_k_l_to_ijkl_.count(mjkl_dense) && std::fabs((*F_lmo_)(i, m)) >= F_CUT) {
                    int mjkl = i_j_k_l_to_ijkl_[mjkl_dense];
                    quadruplet_ext_domain = merge_lists(quadruplet_ext_domain, lmoquadruplet_to_paos_[mjkl]);
                }
            }
            auto S_ijkl = submatrix_rows_and_cols(*S_pao_, quadruplet_ext_domain, lmoquadruplet_to_paos_[ijkl]);
            S_ijkl = linalg::doublet(S_ijkl, X_qno_[ijkl], false, false);

            Tensor<double, 4> R_ijkl("R_ijkl", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            R_ijkl = gamma_ijkl_[ijkl];

            for (int a_ijkl = 0; a_ijkl < nqno_ijkl; ++a_ijkl) {
                for (int b_ijkl = 0; b_ijkl < nqno_ijkl; ++b_ijkl) {
                    for (int c_ijkl = 0; c_ijkl < nqno_ijkl; ++c_ijkl) {
                        for (int d_ijkl = 0; d_ijkl < nqno_ijkl; ++d_ijkl) {
                            (R_ijkl)(a_ijkl, b_ijkl, c_ijkl, d_ijkl) += (T_iajbkcld_[ijkl])(a_ijkl, b_ijkl, c_ijkl, d_ijkl) *
                                ((*e_qno_[ijkl])(a_ijkl) + (*e_qno_[ijkl])(b_ijkl) + (*e_qno_[ijkl])(c_ijkl) + (*e_qno_[ijkl])(d_ijkl) 
                                    - (*F_lmo_)(i, i) - (*F_lmo_)(j, j) - (*F_lmo_)(k, k) - (*F_lmo_)(l, l));
                        } // end d_ijkl
                    } // end c_ijkl
                } // end b_ijkl
            } // end a_ijkl

            for (int m = 0; m < naocc; ++m) {
                int ijkm_dense = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + m;
                int ijml_dense = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + m * naocc + l;
                int imkl_dense = i * std::pow(naocc, 3) + m * std::pow(naocc, 2) + k * naocc + l;
                int mjkl_dense = m * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + l;

                if (l != m && i_j_k_l_to_ijkl_.count(ijkm_dense) && std::fabs((*F_lmo_)(l, m)) >= F_CUT) {
                    int ijkm = i_j_k_l_to_ijkl_[ijkm_dense];
                    std::vector<int> ijkm_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijkm]);
                    auto S_ijkl_ijkm = linalg::doublet(submatrix_rows(*S_ijkl, ijkm_idx_list), X_qno_[ijkm], true, false);
                    auto T_temp = matmul_4d(quadruples_permuter(T_iajbkcld_[ijkm], i, j, k, m), S_ijkl_ijkm, n_qno_[ijkm], n_qno_[ijkl]);
                    T_temp *= (*F_lmo_)(l, m);
                    R_ijkl -= T_temp;
                }

                if (k != m && i_j_k_l_to_ijkl_.count(ijml_dense) && std::fabs((*F_lmo_)(k, m)) >= F_CUT) {
                    int ijml = i_j_k_l_to_ijkl_[ijml_dense];
                    std::vector<int> ijml_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijml]);
                    auto S_ijkl_ijml = linalg::doublet(submatrix_rows(*S_ijkl, ijml_idx_list), X_qno_[ijml], true, false);
                    auto T_temp = matmul_4d(quadruples_permuter(T_iajbkcld_[ijml], i, j, m, l), S_ijkl_ijml, n_qno_[ijml], n_qno_[ijkl]);
                    T_temp *= (*F_lmo_)(k, m);
                    R_ijkl -= T_temp;
                }

                if (j != m && i_j_k_l_to_ijkl_.count(imkl_dense) && std::fabs((*F_lmo_)(j, m)) >= F_CUT) {
                    int imkl = i_j_k_l_to_ijkl_[imkl_dense];
                    std::vector<int> imkl_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[imkl]);
                    auto S_ijkl_imkl = linalg::doublet(submatrix_rows(*S_ijkl, imkl_idx_list), X_qno_[imkl], true, false);
                    auto T_temp = matmul_4d(quadruples_permuter(T_iajbkcld_[imkl], i, m, k, l), S_ijkl_imkl, n_qno_[imkl], n_qno_[ijkl]);
                    T_temp *= (*F_lmo_)(j, m);
                    R_ijkl -= T_temp;
                }

                if (i != m && i_j_k_l_to_ijkl_.count(mjkl_dense) && std::fabs((*F_lmo_)(i, m)) >= F_CUT) {
                    int mjkl = i_j_k_l_to_ijkl_[mjkl_dense];
                    std::vector<int> mjkl_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[mjkl]);
                    auto S_ijkl_mjkl = linalg::doublet(submatrix_rows(*S_ijkl, mjkl_idx_list), X_qno_[mjkl], true, false);
                    auto T_temp = matmul_4d(quadruples_permuter(T_iajbkcld_[mjkl], m, j, k, l), S_ijkl_mjkl, n_qno_[mjkl], n_qno_[ijkl]);
                    T_temp *= (*F_lmo_)(i, m);
                    R_ijkl -= T_temp;
                }
            }

            // => Update T4 Amplitudes <= //
            for (int a_ijkl = 0; a_ijkl < nqno_ijkl; ++a_ijkl) {
                for (int b_ijkl = 0; b_ijkl < nqno_ijkl; ++b_ijkl) {
                    for (int c_ijkl = 0; c_ijkl < nqno_ijkl; ++c_ijkl) {
                        for (int d_ijkl = 0; d_ijkl < nqno_ijkl; ++d_ijkl) {
                            (T_iajbkcld_[ijkl])(a_ijkl, b_ijkl, c_ijkl, d_ijkl) -= (R_ijkl)(a_ijkl, b_ijkl, c_ijkl, d_ijkl) /
                                ((*e_qno_[ijkl])(a_ijkl) + (*e_qno_[ijkl])(b_ijkl) + (*e_qno_[ijkl])(c_ijkl) + (*e_qno_[ijkl])(d_ijkl) 
                                    - (*F_lmo_)(i, i) - (*F_lmo_)(j, j) - (*F_lmo_)(k, k) - (*F_lmo_)(l, l));
                        } // end d_ijkl
                    } // end c_ijkl
                } // end b_ijkl
            } // end a_ijkl
            
            R_iajbkcld_rms[ijkl] = std::sqrt(linear_algebra::dot(R_ijkl, R_ijkl)) / (nqno_ijkl * nqno_ijkl);
        }

        // evaluate convergence
        e_prev = e_curr;
        e_ijkl_old = e_ijkl_;

        // Compute LCCSDT(Q) energy
        timer_on("Compute (Q) Energy");
        e_curr = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr)
        for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
            int ijkl = sorted_quadruplets_[ijkl_sorted];
            double e_ijkl = compute_quadruplet_energy(ijkl, T_iajbkcld_[ijkl]);
            e_ijkl_[ijkl] = e_ijkl;
            e_curr += e_ijkl;
        }
        timer_off("Compute (Q) Energy");

        double r_curr = *max_element(R_iajbkcld_rms.begin(), R_iajbkcld_rms.end());

        r_converged = fabs(r_curr) < options_.get_double("R_CONVERGENCE");
        e_converged = fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE");

        std::time_t time_stop = std::time(nullptr);

        outfile->Printf("  @LCCSDT(Q) iter %3d: %16.12f %10.3e %10.3e %8d\n", iteration, e_curr, e_curr - e_prev, r_curr, (int)time_stop - (int)time_start);

        iteration++;

        if (iteration > max_iteration) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    timer_off("LCCSDT(Q) Iterations");

    return e_curr;
}

double DLPNOCCSDT_Q::compute_energy() {
    timer_on("DLPNO-CCSDT(Q)");

    // Run DLPNO-CCSDT
    double e_dlpno_ccsdt = DLPNOCCSDT::compute_energy();

    einsums::profile::initialize();

    // Clear CCSD integrals
    K_mnij_.clear();
    K_bar_.clear();
    K_bar_chem_.clear();
    L_bar_.clear();
    J_ijab_.clear();
    L_iajb_.clear();
    M_iajb_.clear();
    J_ij_kj_.clear();
    K_ij_kj_.clear();
    K_tilde_chem_.clear();
    K_tilde_phys_.clear();
    L_tilde_.clear();
    Qma_ij_.clear();
    Qab_ij_.clear();
    i_Qk_ij_.clear();
    i_Qa_ij_.clear();
    i_Qa_ij_.clear();
    i_Qa_t1_.clear();
    S_pno_ij_kj_.clear();
    S_pno_ij_nn_.clear();
    S_pno_ij_mn_.clear();

    // Clear CCSDT integrals
    S_ijk_ii_.clear();
    S_ijk_jj_.clear();
    S_ijk_kk_.clear();
    S_ijk_ll_.clear();
    S_ijk_ij_.clear();
    S_ijk_jk_.clear();
    S_ijk_ik_.clear();
    S_ijk_il_.clear();
    S_ijk_jl_.clear();
    S_ijk_kl_.clear();
    S_ijk_lm_.clear();
    S_ijk_ljk_.clear();
    S_ijk_ilk_.clear();
    S_ijk_ijl_.clear();
    S_ijk_mli_.clear();
    S_ijk_mlj_.clear();
    S_ijk_mlk_.clear();
    K_iojv_.clear();
    K_joiv_.clear();
    K_jokv_.clear();
    K_kojv_.clear();
    K_iokv_.clear();
    K_koiv_.clear();
    K_ivjv_.clear();
    K_jvkv_.clear();
    K_ivkv_.clear();
    K_ivov_.clear();
    K_jvov_.clear();
    K_kvov_.clear();
    K_ivvv_.clear();
    K_jvvv_.clear();
    K_kvvv_.clear();
    q_io_.clear();
    q_jo_.clear();
    q_ko_.clear();
    q_iv_.clear();
    q_jv_.clear();
    q_kv_.clear();
    q_ov_.clear();
    q_vv_.clear();

    // Re-create T_iajbkc_clone intermediate
    int n_lmo_triplets = ijk_to_i_j_k_.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        auto &[i, j, k] = ijk_to_i_j_k_[ijk];

        T_iajbkc_clone_[ijk] = Tensor<double, 3>("T_ijk", n_tno_[ijk], n_tno_[ijk], n_tno_[ijk]);
        ::memcpy(T_iajbkc_clone_[ijk].data(), T_iajbkc_[ijk]->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
    }

    double t_cut_qno_pre = options_.get_double("T_CUT_QNO_PRE");
    double t_cut_qno = options_.get_double("T_CUT_QNO");

    // Step 1: Perform the prescreening
    outfile->Printf("\n   Starting Quadruplet Prescreening...\n");
    outfile->Printf("     T_CUT_QNO set to %6.3e \n", t_cut_qno_pre);
    outfile->Printf("     T_CUT_DO  set to %6.3e \n", options_.get_double("T_CUT_DO_QUADS_PRE"));
    outfile->Printf("     T_CUT_MKN set to %6.3e \n\n", options_.get_double("T_CUT_MKN_QUADS_PRE"));

    quadruples_sparsity(true);
    qno_transform(t_cut_qno_pre);
    double E_Q0_pre = compute_gamma_ijkl(false);
    outfile->Printf("    (Initial) DLPNO-(Q0) Correlation Energy: %16.12f\n\n", E_Q0_pre);

    // Step 2: Compute DLPNO-CCSDT(Q0) energy with surviving quadruplets
    outfile->Printf("\n   Continuing computation with surviving quadruplets...\n");
    outfile->Printf("     Eliminated all quadruplets with energy less than %6.3e Eh... \n\n", options_.get_double("T_CUT_QUADS_WEAK"));
    quadruples_sparsity(false);
    outfile->Printf("    * Energy Contribution From Screened Quadruplets: %.12f \n\n", de_lccsdt_q_screened_);
    
    outfile->Printf("     T_CUT_QNO (re)set to %6.3e \n", options_.get_double("T_CUT_QNO"));
    outfile->Printf("     T_CUT_DO  (re)set to %6.3e \n", options_.get_double("T_CUT_DO_QUADS"));
    outfile->Printf("     T_CUT_MKN (re)set to %6.3e \n\n", options_.get_double("T_CUT_MKN_QUADS"));

    qno_transform(t_cut_qno);
    double E_Q0 = compute_gamma_ijkl(false);
    outfile->Printf("    (Total) DLPNO-(Q0) Correlation Energy:      %16.12f\n", E_Q0 + de_lccsdt_q_screened_);
    outfile->Printf("    * Screened Quadruplets Contribution:        %16.12f\n", de_lccsdt_q_screened_);

    double e_scf = variables_["SCF TOTAL ENERGY"];
    double e_ccsdt_q_corr = E_Q0 + de_lccsdt_q_screened_ + e_lccsdt_ + dE_T_rank_ + de_weak_ + de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
    double e_ccsdt_q_total = e_scf + e_ccsdt_q_corr;

    outfile->Printf("\n\n  @Total DLPNO-CCSDT(Q0) Energy: %16.12f\n", e_ccsdt_q_total);

    double e_total = e_ccsdt_q_total;

    if (!options_.get_bool("Q0_ONLY")) {
        // STEP 3: Iterative (Q) computations
        outfile->Printf("\n\n  ==> Computing Full Iterative (Q) <==\n\n");

        double t_cut_qno_strong_scale = options_.get_double("T_CUT_QNO_STRONG_SCALE");
        double t_cut_qno_weak_scale = options_.get_double("T_CUT_QNO_WEAK_SCALE");
        outfile->Printf("     T_CUT_QNO (re)set to %6.3e for strong triples \n", t_cut_qno * t_cut_qno_strong_scale);
        outfile->Printf("     T_CUT_QNO (re)set to %6.3e for weak triples   \n\n", t_cut_qno * t_cut_qno_weak_scale);

        // Sort quadruplets into "strong" and "weak" quadruplets
        sort_quadruplets(E_Q0);
        qno_transform(t_cut_qno);
        estimate_memory();

        double E_Q0_crude = compute_gamma_ijkl(true);
        double E_Q = lccsdt_q_iterations();
        double dE_Q = E_Q - E_Q0_crude;

        outfile->Printf("\n");
        outfile->Printf("    DLPNO-(Q0) energy at looser tolerance:      %16.12f\n", E_Q0_crude);
        outfile->Printf("    DLPNO-(Q)  energy at looser tolerance:      %16.12f\n", E_Q);
        outfile->Printf("    * Net Iterative (Q) contribution:           %16.12f\n\n", dE_Q);

        outfile->Printf("    (Total) DLPNO-(Q) Correlation Energy:       %16.12f\n", E_Q0 + dE_Q + de_lccsdt_q_screened_);
        outfile->Printf("    * DLPNO-(Q0) Contribution:                  %16.12f\n", E_Q0);
        outfile->Printf("    * DLPNO-(Q) Contribution:                   %16.12f\n", dE_Q);
        outfile->Printf("    * Screened Quadruplets Contribution:        %16.12f\n", de_lccsdt_q_screened_);

        e_ccsdt_q_corr = E_Q0 + dE_Q + de_lccsdt_q_screened_ + e_lccsdt_ + dE_T_rank_ + de_weak_ + de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
        e_ccsdt_q_total = e_scf + e_ccsdt_q_corr;
        e_total = e_ccsdt_q_total;

        outfile->Printf("\n\n  @Total DLPNO-CCSDT(Q) Energy: %16.12f\n", e_ccsdt_q_total);
    }
    outfile->Printf("*** I can't believe this job finished! A thousand hallelujahs!\n\n");

    set_scalar_variable("CCSDT(Q) CORRELATION ENERGY", e_ccsdt_q_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsdt_q_corr);
    set_scalar_variable("CCSDT(Q) TOTAL ENERGY", e_ccsdt_q_total);
    set_scalar_variable("CURRENT ENERGY", e_ccsdt_q_total);

    einsums::profile::finalize();

    timer_off("DLPNO-CCSDT(Q)");

    return e_total;
}

}
}